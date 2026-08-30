"""Cryogenic thermal fuel management for the LNG circuit.

LNG arrives from the tank at around 110 K and must reach the injector warm enough to
flash usefully. That makes fuel temperature a design variable rather than a boundary
condition, and it is pulled in four directions at once:

- **enough superheat** at the nozzle for the intended flash regime, since flashing is
  what atomizes the fuel without needing high injection pressure;
- **enough subcooling margin upstream**, or the feed line boils and the pump sees
  two-phase fuel;
- **enough autoignition margin**, since warmer fuel means a warmer premixed mixture;
- **within the heat actually available** from engine and aircraft loads.

The window where all four hold is the answer to "how much should I heat the LNG", and
finding it turns a guess into a number.

The module also screens the circuit that is *not* running. A dual-fuel nozzle always has
one circuit hot and stagnant: Jet-A sitting in a hot nozzle cokes, and LNG sitting in one
boils and vapour-locks. Either can veto an otherwise attractive packaging.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from scipy.optimize import brentq

from fuelnozzle.crn.autoignition import AutoignitionMargin, AutoignitionVerdict
from fuelnozzle.crn.status import GateStatus
from fuelnozzle.feed import LNGFeedLineResult
from fuelnozzle.models import ModelWarning, ThermodynamicState, WarningSeverity
from fuelnozzle.properties import CoolPropLNGProvider, PropertyCalculationError

#: Jet-A begins to form deposits above roughly this wall temperature. It is an input
#: rather than a constant because it depends on fuel additives and residence time.
DEFAULT_JET_A_COKING_LIMIT_K = 450.0

#: Minimum subcooling to keep upstream of the nozzle before vapour lock is a concern.
DEFAULT_MINIMUM_SUBCOOLING_K = 5.0


class SupercriticalFeedError(ValueError):
    """The feed pressure is above the critical point, where saturation does not exist."""


class InfeasibleThermalTargetError(ValueError):
    """A requested inverse thermal target lies outside its declared bracket or phase limits."""


@dataclass(frozen=True)
class HeatSinkBudget:
    """How much cooling the LNG stream can absorb, and where it goes.

    Splitting the duty into three parts matters because they are not interchangeable.
    Sensible heating of the liquid is free of consequence; latent heat means the fuel is
    boiling, which is wanted at the nozzle and unwanted in the line; superheating the
    vapour is what actually shortens ignition delay.
    """

    mass_flow_kg_s: float
    tank_temperature_k: float
    nozzle_temperature_k: float
    pressure_pa: float
    total_duty_w: float
    sensible_liquid_w: float
    latent_w: float
    superheat_w: float
    vapor_quality: float
    warnings: tuple[ModelWarning, ...]


@dataclass(frozen=True)
class FeedPathHeatBudget:
    """Energy and phase ledger along an already solved LNG feed-pressure path."""

    total_heat_w: float
    enthalpy_rise_w: float
    energy_residual_w: float
    pressure_drop_pa: float
    outlet_vapor_quality_mass: float
    first_two_phase_position_m: float | None


def feed_path_heat_budget(
    result: LNGFeedLineResult,
    mass_flow_kg_s: float,
) -> FeedPathHeatBudget:
    """Cross-check tank-to-injector duty on the actual pressure/enthalpy march."""
    if mass_flow_kg_s <= 0.0:
        raise ValueError("Feed-path heat budget requires positive mass flow")
    enthalpy_rise = mass_flow_kg_s * (
        result.outlet_state.enthalpy_j_kg - result.inlet_state.enthalpy_j_kg
    )
    return FeedPathHeatBudget(
        total_heat_w=result.total_heat_leak_w,
        enthalpy_rise_w=enthalpy_rise,
        energy_residual_w=enthalpy_rise - result.total_heat_leak_w,
        pressure_drop_pa=result.pressure_drop_pa,
        outlet_vapor_quality_mass=result.outlet_state.vapor_quality_mass or 0.0,
        first_two_phase_position_m=result.first_two_phase_position_m,
    )


def heat_sink_budget(
    provider: CoolPropLNGProvider,
    mass_flow_kg_s: float,
    tank_temperature_k: float,
    nozzle_temperature_k: float,
    pressure_pa: float,
) -> HeatSinkBudget:
    """Cooling duty available between tank and injector, decomposed by mechanism."""
    if mass_flow_kg_s <= 0.0:
        raise ValueError("Heat sink budget requires a positive mass flow")

    warnings: list[ModelWarning] = []
    try:
        bubble = provider.bubble_state_at_pressure(pressure_pa)
        dew = provider.dew_state_at_pressure(pressure_pa)
    except PropertyCalculationError as error:
        raise ValueError(f"CoolProp could not evaluate LNG saturation: {error}") from error

    latent_heat = dew.enthalpy_j_kg - bubble.enthalpy_j_kg
    saturation_temperature = bubble.temperature_k

    def state_at_temperature(temperature_k: float) -> ThermodynamicState:
        if abs(temperature_k - bubble.temperature_k) <= 1.0e-6:
            return bubble
        if abs(temperature_k - dew.temperature_k) <= 1.0e-6:
            return dew
        return provider.state_pt(pressure_pa, temperature_k)

    tank_state = state_at_temperature(tank_temperature_k)
    nozzle_state = state_at_temperature(nozzle_temperature_k)
    tank_enthalpy = tank_state.enthalpy_j_kg
    nozzle_enthalpy = nozzle_state.enthalpy_j_kg
    if nozzle_enthalpy < tank_enthalpy:
        raise ValueError("Nozzle state is colder in enthalpy than the tank state")

    sensible_enthalpy = max(
        0.0,
        min(nozzle_enthalpy, bubble.enthalpy_j_kg) - tank_enthalpy,
    )
    latent_enthalpy = max(
        0.0,
        min(nozzle_enthalpy, dew.enthalpy_j_kg)
        - max(tank_enthalpy, bubble.enthalpy_j_kg),
    )
    superheat_enthalpy = max(
        0.0,
        nozzle_enthalpy - max(tank_enthalpy, dew.enthalpy_j_kg),
    )
    sensible = mass_flow_kg_s * sensible_enthalpy
    latent = mass_flow_kg_s * latent_enthalpy
    superheat = mass_flow_kg_s * superheat_enthalpy
    total = mass_flow_kg_s * (nozzle_enthalpy - tank_enthalpy)
    quality = (
        nozzle_state.vapor_quality_mass
        if nozzle_state.vapor_quality_mass is not None
        else (
            0.0
            if nozzle_enthalpy <= bubble.enthalpy_j_kg
            else 1.0
            if nozzle_enthalpy >= dew.enthalpy_j_kg
            else latent_enthalpy / latent_heat
        )
    )

    if nozzle_enthalpy < bubble.enthalpy_j_kg:
        # Still liquid at the injector: no flashing, so no flash atomization.
        warnings.append(
            ModelWarning(
                code="LNG_SUBCOOLED_AT_NOZZLE",
                severity=WarningSeverity.INFO,
                message=(
                    f"LNG reaches the injector at {nozzle_temperature_k:.1f} K, below its "
                    f"{saturation_temperature:.1f} K saturation temperature. It will not "
                    "flash, so atomization must come from injection pressure alone."
                ),
            )
        )
    return HeatSinkBudget(
        mass_flow_kg_s=mass_flow_kg_s,
        tank_temperature_k=tank_temperature_k,
        nozzle_temperature_k=nozzle_temperature_k,
        pressure_pa=pressure_pa,
        total_duty_w=total,
        sensible_liquid_w=sensible,
        latent_w=latent,
        superheat_w=superheat,
        vapor_quality=quality,
        warnings=tuple(warnings),
    )


def fuel_temperature_for_target_superheat(
    provider: CoolPropLNGProvider,
    chamber_pressure_pa: float,
    feed_pressure_pa: float,
    target_superheat_k: float,
) -> float:
    """Nozzle-inlet temperature giving a target superheat over chamber saturation.

    This is the inverse of the usual question. Rather than asking what a given fuel
    temperature produces, it asks what temperature is needed to reach an intended flash
    condition, which is what a designer actually wants to know.
    """
    if target_superheat_k <= 0.0:
        raise ValueError("Target superheat must be positive")
    try:
        chamber_saturation = provider.bubble_state_at_pressure(
            chamber_pressure_pa
        ).temperature_k
        feed_saturation = provider.bubble_state_at_pressure(feed_pressure_pa).temperature_k
    except PropertyCalculationError as error:
        raise SupercriticalFeedError(
            "Saturation is undefined above the critical point, so a target superheat "
            f"cannot be inverted at these pressures. CoolProp reported: {error}"
        ) from error
    target = chamber_saturation + target_superheat_k
    if target >= feed_saturation:
        raise InfeasibleThermalTargetError(
            f"Target temperature {target:.2f} K reaches or exceeds the "
            f"{feed_saturation:.2f} K feed bubble point; the requested superheat would "
            "boil the fuel upstream."
        )
    return target


@dataclass(frozen=True)
class HeatSourceModel:
    """Declared off-design heat source used to bound recoverable LNG heating."""

    hot_inlet_temperature_k: float
    hot_mass_flow_kg_s: float
    hot_specific_heat_j_kg_k: float
    effectiveness: float
    minimum_pinch_k: float
    fuel_side_pressure_loss_pa: float = 0.0
    evidence_id: str | None = None

    def __post_init__(self) -> None:
        if min(
            self.hot_inlet_temperature_k,
            self.hot_mass_flow_kg_s,
            self.hot_specific_heat_j_kg_k,
        ) <= 0.0:
            raise ValueError(
                "Heat-source temperature, mass flow, and heat capacity must be positive"
            )
        if not 0.0 < self.effectiveness <= 1.0:
            raise ValueError("Heat-exchanger effectiveness must lie in (0, 1]")
        if self.minimum_pinch_k < 0.0 or self.fuel_side_pressure_loss_pa < 0.0:
            raise ValueError("Pinch and pressure loss cannot be negative")

    def available_heat_w(
        self,
        fuel_inlet_temperature_k: float,
        fuel_outlet_temperature_k: float,
    ) -> float:
        hot_capacity_w_k = self.hot_mass_flow_kg_s * self.hot_specific_heat_j_kg_k
        effectiveness_limit = (
            self.effectiveness
            * hot_capacity_w_k
            * max(self.hot_inlet_temperature_k - fuel_inlet_temperature_k, 0.0)
        )
        pinch_limit = hot_capacity_w_k * max(
            self.hot_inlet_temperature_k
            - fuel_outlet_temperature_k
            - self.minimum_pinch_k,
            0.0,
        )
        return min(effectiveness_limit, pinch_limit)


@dataclass(frozen=True)
class ThermalWindowPoint:
    """One candidate fuel temperature, with every constraint evaluated."""

    fuel_temperature_k: float
    superheat_at_chamber_k: float
    subcooling_in_feed_k: float
    has_sufficient_superheat: bool
    has_sufficient_subcooling: bool
    heat_duty_w: float
    within_available_heat: bool
    autoignition_status: GateStatus = GateStatus.UNKNOWN
    heat_source_status: GateStatus = GateStatus.UNKNOWN

    @property
    def is_feasible(self) -> bool:
        return (
            self.has_sufficient_superheat
            and self.has_sufficient_subcooling
            and self.within_available_heat
            and self.autoignition_status is not GateStatus.FAIL
        )

    @property
    def acceptance_status(self) -> GateStatus:
        if not self.is_feasible:
            return GateStatus.FAIL
        if (
            self.autoignition_status is GateStatus.UNKNOWN
            or self.heat_source_status is GateStatus.UNKNOWN
        ):
            return GateStatus.UNKNOWN
        return GateStatus.PASS


@dataclass(frozen=True)
class ThermalWindow:
    """The set of fuel temperatures that satisfy every constraint at once.

    The width of this window is governed by the gap between the saturation temperature
    at feed pressure and at chamber pressure. Superheat at the injector and subcooling in
    the line compete for that gap, so raising pump pressure is what buys thermal design
    freedom, and below some pump pressure no fuel temperature works at all.
    """

    points: tuple[ThermalWindowPoint, ...]
    feasible_temperatures_k: tuple[float, ...]
    chamber_saturation_k: float
    feed_saturation_k: float
    warnings: tuple[ModelWarning, ...]

    @property
    def saturation_gap_k(self) -> float:
        """Temperature room the superheat and subcooling constraints must share."""
        return self.feed_saturation_k - self.chamber_saturation_k

    @property
    def is_empty(self) -> bool:
        return not self.feasible_temperatures_k

    @property
    def bounds_k(self) -> tuple[float, float] | None:
        if self.is_empty:
            return None
        return (min(self.feasible_temperatures_k), max(self.feasible_temperatures_k))


def thermal_window(
    provider: CoolPropLNGProvider,
    temperatures_k: tuple[float, ...] | list[float],
    *,
    chamber_pressure_pa: float,
    feed_pressure_pa: float,
    mass_flow_kg_s: float,
    tank_temperature_k: float,
    minimum_superheat_k: float = 5.0,
    minimum_subcooling_k: float = DEFAULT_MINIMUM_SUBCOOLING_K,
    available_heat_w: float | None = None,
    heat_source: HeatSourceModel | None = None,
    autoignition_evaluator: Callable[[float], AutoignitionMargin] | None = None,
) -> ThermalWindow:
    """Evaluate every constraint across a range of candidate fuel temperatures.

    Reporting a window rather than a single optimum is deliberate: the constraints move
    with operating point and with hardware detail, so what a designer needs is the room
    available, not one number sitting at an unknown distance from a boundary.
    """
    try:
        chamber_saturation = provider.bubble_state_at_pressure(
            chamber_pressure_pa
        ).temperature_k
        feed_saturation = provider.bubble_state_at_pressure(feed_pressure_pa).temperature_k
    except PropertyCalculationError as error:
        raise SupercriticalFeedError(
            "Saturation is undefined at one of the requested pressures, which means it "
            "is above the critical point of this LNG composition. Supercritical injection "
            "is a legitimate design choice, but the subcooling and superheat constraints "
            "used here are defined only below the critical point and do not apply to it. "
            f"CoolProp reported: {error}"
        ) from error

    points: list[ThermalWindowPoint] = []
    for temperature in sorted(temperatures_k):
        superheat = temperature - chamber_saturation
        subcooling = feed_saturation - temperature
        try:
            # Heating occurs upstream at feed pressure. The chamber-pressure flash is a
            # separate isenthalpic nozzle calculation.
            budget = heat_sink_budget(
                provider, mass_flow_kg_s, tank_temperature_k, temperature,
                feed_pressure_pa,
            )
            duty = budget.total_duty_w
        except (ValueError, PropertyCalculationError):
            duty = float("nan")

        source_limit = available_heat_w
        source_status = GateStatus.UNKNOWN
        if heat_source is not None:
            source_limit = heat_source.available_heat_w(tank_temperature_k, temperature)
            source_status = (
                GateStatus.PASS if heat_source.evidence_id else GateStatus.UNKNOWN
            )
        elif available_heat_w is not None:
            source_status = GateStatus.UNKNOWN

        ignition_status = GateStatus.UNKNOWN
        if autoignition_evaluator is not None:
            verdict = autoignition_evaluator(temperature).verdict
            ignition_status = {
                AutoignitionVerdict.SAFE: GateStatus.PASS,
                AutoignitionVerdict.MARGINAL: GateStatus.FAIL,
                AutoignitionVerdict.UNSAFE: GateStatus.FAIL,
                AutoignitionVerdict.UNKNOWN: GateStatus.UNKNOWN,
                AutoignitionVerdict.NO_PREMIXER: GateStatus.PASS,
            }[verdict]

        points.append(
            ThermalWindowPoint(
                fuel_temperature_k=temperature,
                superheat_at_chamber_k=superheat,
                subcooling_in_feed_k=subcooling,
                has_sufficient_superheat=superheat >= minimum_superheat_k,
                has_sufficient_subcooling=subcooling >= minimum_subcooling_k,
                heat_duty_w=duty,
                within_available_heat=(
                    source_limit is None
                    or (duty == duty and duty <= source_limit)
                ),
                autoignition_status=ignition_status,
                heat_source_status=source_status,
            )
        )

    feasible = tuple(point.fuel_temperature_k for point in points if point.is_feasible)
    warnings: list[ModelWarning] = []
    if not feasible:
        warnings.append(
            ModelWarning(
                code="THERMAL_WINDOW_EMPTY",
                severity=WarningSeverity.ERROR,
                message=(
                    "No fuel temperature satisfies superheat, subcooling, and available "
                    "heat at once. The nozzle pressure drop, feed pressure, or chamber "
                    "pressure must change; heating the fuel differently cannot fix it."
                ),
            )
        )
    elif len(feasible) < 3:
        warnings.append(
            ModelWarning(
                code="THERMAL_WINDOW_NARROW",
                severity=WarningSeverity.WARNING,
                message=(
                    f"Only {len(feasible)} of {len(points)} candidate temperatures are "
                    "feasible. The design sits close to a constraint boundary."
                ),
            )
        )
    if not feasible and (feed_saturation - chamber_saturation) < (
        minimum_superheat_k + minimum_subcooling_k
    ):
        warnings.append(
            ModelWarning(
                code="SATURATION_GAP_TOO_NARROW",
                severity=WarningSeverity.ERROR,
                message=(
                    f"Saturation temperature rises only "
                    f"{feed_saturation - chamber_saturation:.1f} K between chamber and "
                    f"feed pressure, less than the {minimum_superheat_k:.1f} K superheat "
                    f"plus {minimum_subcooling_k:.1f} K subcooling required. Raising pump "
                    "pressure widens this gap; adjusting fuel temperature cannot."
                ),
            )
        )
    return ThermalWindow(
        tuple(points), feasible, chamber_saturation, feed_saturation, tuple(warnings)
    )


@dataclass(frozen=True)
class IdleCircuitScreen:
    """Whether the fuel circuit that is shut off survives the mission segment."""

    idle_fuel: str
    wall_temperature_k: float
    coking_limit_k: float
    jet_a_coking_safe: bool | None
    lng_vapor_locked: bool | None
    warnings: tuple[ModelWarning, ...]
    transient_status: GateStatus = GateStatus.UNKNOWN

    @property
    def is_safe(self) -> bool:
        return not any(
            warning.severity is WarningSeverity.ERROR for warning in self.warnings
        )

    @property
    def acceptance_status(self) -> GateStatus:
        if not self.is_safe:
            return GateStatus.FAIL
        return self.transient_status


def idle_circuit_screen(
    active_fuel: str,
    wall_temperature_k: float,
    *,
    provider: CoolPropLNGProvider | None = None,
    circuit_pressure_pa: float = 5.0e5,
    coking_limit_k: float = DEFAULT_JET_A_COKING_LIMIT_K,
    soak_duration_s: float | None = None,
    purge_mass_flow_kg_s: float | None = None,
    purge_duration_s: float | None = None,
    restart_demonstrated: bool = False,
) -> IdleCircuitScreen:
    """Screen the stagnant circuit for coking or vapour lock.

    A dual-fuel nozzle always has one circuit hot and not flowing. With no flow there is
    no convective cooling, so the stagnant fuel sits at close to wall temperature. This
    is a packaging constraint that can veto an otherwise attractive design, and it is
    cheaper to discover here than on a test stand.
    """
    warnings: list[ModelWarning] = []
    jet_a_safe: bool | None = None
    vapor_locked: bool | None = None

    if active_fuel == "lng":
        jet_a_safe = wall_temperature_k < coking_limit_k
        if not jet_a_safe:
            warnings.append(
                ModelWarning(
                    code="IDLE_JET_A_COKING_RISK",
                    severity=WarningSeverity.ERROR,
                    message=(
                        f"With LNG burning, the idle Jet-A circuit sits at "
                        f"{wall_temperature_k:.1f} K, above the {coking_limit_k:.1f} K "
                        "coking limit. Stagnant fuel in a hot passage forms deposits that "
                        "block the injector."
                    ),
                )
            )
    elif active_fuel == "jet_a":
        if provider is None:
            warnings.append(
                ModelWarning(
                    code="IDLE_LNG_SCREEN_UNAVAILABLE",
                    severity=WarningSeverity.WARNING,
                    message=(
                        "No LNG property provider was supplied, so the idle LNG circuit "
                        "was not screened for vapour lock."
                    ),
                )
            )
        else:
            saturation = provider.bubble_state_at_pressure(
                circuit_pressure_pa
            ).temperature_k
            vapor_locked = wall_temperature_k > saturation
            if vapor_locked:
                warnings.append(
                    ModelWarning(
                        code="IDLE_LNG_VAPOR_LOCK",
                        severity=WarningSeverity.ERROR,
                        message=(
                            f"With Jet-A burning, the idle LNG circuit sits at "
                            f"{wall_temperature_k:.1f} K, above its "
                            f"{saturation:.1f} K saturation temperature at "
                            f"{circuit_pressure_pa:.3g} Pa. The stagnant fuel boils, and "
                            "restart delivers vapour rather than liquid."
                        ),
                    )
                )
    else:
        raise ValueError(f"Unknown active fuel {active_fuel!r}")

    transient_status = GateStatus.UNKNOWN
    if soak_duration_s is not None:
        if soak_duration_s < 0.0:
            raise ValueError("Soak duration cannot be negative")
        purge_complete = (
            purge_mass_flow_kg_s is not None
            and purge_mass_flow_kg_s > 0.0
            and purge_duration_s is not None
            and purge_duration_s > 0.0
        )
        transient_status = (
            GateStatus.PASS
            if purge_complete and restart_demonstrated
            else GateStatus.FAIL
        )
        if transient_status is GateStatus.FAIL:
            warnings.append(
                ModelWarning(
                    code="IDLE_CIRCUIT_TRANSIENT_UNSAFE",
                    severity=WarningSeverity.ERROR,
                    message=(
                        "A thermal-soak interval was declared without both a positive purge "
                        "and demonstrated restart capability."
                    ),
                )
            )
    else:
        warnings.append(
            ModelWarning(
                code="IDLE_CIRCUIT_TRANSIENT_UNKNOWN",
                severity=WarningSeverity.WARNING,
                message=(
                    "Soak duration, purge flow/time, deposit evidence, and restart evidence "
                    "were not supplied; steady wall-temperature screening cannot establish "
                    "mission restart safety."
                ),
            )
        )

    return IdleCircuitScreen(
        idle_fuel="jet_a" if active_fuel == "lng" else "lng",
        wall_temperature_k=wall_temperature_k,
        coking_limit_k=coking_limit_k,
        jet_a_coking_safe=jet_a_safe,
        lng_vapor_locked=vapor_locked,
        warnings=tuple(warnings),
        transient_status=transient_status,
    )


def required_feed_heat_w(
    provider: CoolPropLNGProvider,
    mass_flow_kg_s: float,
    tank_temperature_k: float,
    target_temperature_k: float,
    pressure_pa: float,
) -> float:
    """Heat that must be added in the feed line to reach a target injector temperature."""
    budget = heat_sink_budget(
        provider, mass_flow_kg_s, tank_temperature_k, target_temperature_k, pressure_pa
    )
    return budget.total_duty_w


def saturation_temperature_k(
    provider: CoolPropLNGProvider, pressure_pa: float
) -> float:
    """Bubble-point temperature, exposed for constraint construction."""
    return provider.bubble_state_at_pressure(pressure_pa).temperature_k


def solve_temperature_for_duty(
    provider: CoolPropLNGProvider,
    mass_flow_kg_s: float,
    tank_temperature_k: float,
    pressure_pa: float,
    duty_w: float,
    bracket_k: tuple[float, float] = (100.0, 400.0),
) -> float:
    """Injector temperature reached when a known heat load is applied to the fuel.

    The forward direction of :func:`required_feed_heat_w`: given the heat the engine has
    to dump, how warm does the fuel arrive.
    """
    if duty_w <= 0.0:
        raise ValueError("Duty must be positive")

    def residual(temperature_k: float) -> float:
        return (
            required_feed_heat_w(
                provider, mass_flow_kg_s, tank_temperature_k, temperature_k, pressure_pa
            )
            - duty_w
        )

    low, high = max(bracket_k[0], tank_temperature_k), bracket_k[1]
    if low >= high:
        raise InfeasibleThermalTargetError(
            "Temperature bracket does not extend above the tank temperature"
        )
    if residual(low) > 0.0:
        raise InfeasibleThermalTargetError(
            "Requested duty lies below the duty at the lower temperature bracket"
        )
    if residual(high) < 0.0:
        raise InfeasibleThermalTargetError(
            "Requested duty lies above the duty at the upper temperature bracket"
        )
    return float(brentq(residual, low, high, xtol=1.0e-4))
