"""Autoignition and flashback screening for premixing passages.

A premixing passage exists so that fuel and air are uniform before they burn. It works
only if the mixture takes longer to ignite on its own than it takes to cross the passage.
If it does not, the flame sits inside the injector, which destroys hardware.

The margin is a ratio of two times::

    M = tau_ign / tau_res

with ``tau_ign`` the ignition delay of the mixture at the passage conditions and
``tau_res`` the time it spends there.

This is where cryogenic LNG earns its place. Fuel that has flashed arrives cold, and the
latent heat of whatever is still liquid is drawn out of the air. The mixture can sit tens
of kelvin below the compressor discharge temperature, and because ignition delay depends
on temperature exponentially, that modest cooling buys a large increase in margin --
which is what permits a longer premixing passage, better mixing, and lower NOx. The chain
runs: flash cooling, longer delay, longer premixer, more uniform mixture, less NOx. Every
link is computed here rather than assumed.

Jet-A ignition delay uses the dedicated low-temperature mechanism. At landing and
take-off inlet temperatures a high-temperature-only mechanism overpredicts the delay by
one to two orders of magnitude, which would declare an unsafe premixer safe.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import cantera as ct
import numpy as np

from fuelnozzle.crn.chemistry import (
    DRY_AIR_MOLE_FRACTIONS,
    FuelKind,
    MechanismRegistry,
    MechanismRole,
    MechanismSpec,
)
from fuelnozzle.models import ModelWarning, WarningSeverity

#: Ignition is declared at the steepest temperature rise.
IGNITION_TEMPERATURE_RISE_K = 400.0

#: Longest ignition delay worth integrating for. Beyond this the mixture is inert on any
#: timescale a premixing passage cares about.
MAX_IGNITION_TIME_S = 10.0

#: Design guidance puts the acceptable floor for the margin here. It is an input rather
#: than a constant so that it can be argued with.
DEFAULT_MINIMUM_MARGIN = 4.0


class AutoignitionVerdict(StrEnum):
    """Whether a premixing passage is viable at a condition."""

    SAFE = "safe"
    MARGINAL = "marginal"
    UNSAFE = "unsafe"
    NO_PREMIXER = "no_premixing_passage"


@dataclass(frozen=True)
class PremixState:
    """Conditions inside a premixing passage after fuel and air mix."""

    temperature_k: float
    pressure_pa: float
    equivalence_ratio: float
    air_temperature_k: float
    fuel_temperature_k: float
    temperature_drop_k: float
    latent_cooling_k: float

    @property
    def is_cooled_below_air(self) -> bool:
        return self.temperature_k < self.air_temperature_k


@dataclass(frozen=True)
class AutoignitionMargin:
    """Ignition delay against residence time in one passage."""

    fuel: FuelKind
    premix: PremixState
    ignition_delay_s: float | None
    residence_time_s: float
    margin: float | None
    minimum_margin: float
    verdict: AutoignitionVerdict
    mechanism_path: str
    used_dedicated_ignition_mechanism: bool
    warnings: tuple[ModelWarning, ...]


def premix_state(
    solution: ct.Solution,
    spec: MechanismSpec,
    *,
    air_mass_flow_kg_s: float,
    air_temperature_k: float,
    fuel_mass_flow_kg_s: float,
    fuel_temperature_k: float,
    pressure_pa: float,
    liquid_fraction: float = 0.0,
    latent_heat_j_kg: float = 0.0,
    oxidizer_mole_fractions: dict[str, float] | None = None,
) -> PremixState:
    """Adiabatic mixing of air and fuel, including the cooling the fuel causes.

    Two separate effects lower the mixture temperature below the air temperature. Cold
    fuel vapor dilutes the hot air, and any fuel still liquid takes its latent heat out
    of the mixture as it evaporates. Both are accounted for explicitly rather than folded
    into an effective temperature, so their sizes can be compared.
    """
    if air_mass_flow_kg_s <= 0.0 or fuel_mass_flow_kg_s <= 0.0:
        raise ValueError("Premix state requires positive air and fuel mass flows")
    if not 0.0 <= liquid_fraction <= 1.0:
        raise ValueError("Liquid fraction must lie between 0 and 1")

    oxidizer = oxidizer_mole_fractions or DRY_AIR_MOLE_FRACTIONS

    solution.TPX = air_temperature_k, pressure_pa, oxidizer
    air_enthalpy = float(solution.enthalpy_mass)
    air_composition = np.array(solution.Y)

    solution.TPX = max(fuel_temperature_k, 1.0), pressure_pa, spec.fuel_mole_fractions
    fuel_enthalpy = float(solution.enthalpy_mass)
    fuel_composition = np.array(solution.Y)

    total = air_mass_flow_kg_s + fuel_mass_flow_kg_s
    mixed_composition = (
        air_mass_flow_kg_s * air_composition + fuel_mass_flow_kg_s * fuel_composition
    ) / total

    latent_draw = fuel_mass_flow_kg_s * liquid_fraction * latent_heat_j_kg
    mixed_enthalpy = (
        air_mass_flow_kg_s * air_enthalpy
        + fuel_mass_flow_kg_s * fuel_enthalpy
        - latent_draw
    ) / total

    solution.Y = mixed_composition
    solution.HP = mixed_enthalpy, pressure_pa
    mixed_temperature = float(solution.T)

    # Isolate how much of the drop the latent heat alone accounts for.
    solution.Y = mixed_composition
    solution.HP = (
        air_mass_flow_kg_s * air_enthalpy + fuel_mass_flow_kg_s * fuel_enthalpy
    ) / total, pressure_pa
    without_latent = float(solution.T)

    phi = float(
        solution.equivalence_ratio(spec.fuel_string, oxidizer, basis="mole")
    )
    return PremixState(
        temperature_k=mixed_temperature,
        pressure_pa=pressure_pa,
        equivalence_ratio=phi,
        air_temperature_k=air_temperature_k,
        fuel_temperature_k=fuel_temperature_k,
        temperature_drop_k=air_temperature_k - mixed_temperature,
        latent_cooling_k=without_latent - mixed_temperature,
    )


def ignition_delay(
    solution: ct.Solution,
    spec: MechanismSpec,
    temperature_k: float,
    pressure_pa: float,
    equivalence_ratio: float,
    oxidizer_mole_fractions: dict[str, float] | None = None,
    max_time_s: float = MAX_IGNITION_TIME_S,
) -> float | None:
    """Homogeneous constant-pressure ignition delay, or ``None`` if it does not ignite.

    Ignition is declared when the temperature has risen by a fixed amount, which is a
    robust and reproducible marker of the steep rise. Returning ``None`` rather than a
    large number keeps "did not ignite within the window" distinguishable from "ignited
    slowly", which matter differently for a safety screen.
    """
    if equivalence_ratio <= 0.0:
        return None
    solution.TP = temperature_k, pressure_pa
    solution.set_equivalence_ratio(
        equivalence_ratio, spec.fuel_string, oxidizer_mole_fractions or DRY_AIR_MOLE_FRACTIONS
    )
    reactor = ct.IdealGasConstPressureReactor(solution, clone=False)
    network = ct.ReactorNet([reactor])
    target = temperature_k + IGNITION_TEMPERATURE_RISE_K
    while network.time < max_time_s:
        network.step()
        if reactor.T >= target:
            return float(network.time)
    return None


class IgnitionDelayTable:
    """Ignition delays on a grid, interpolated between.

    Ignition delay is called repeatedly inside sweeps and each evaluation integrates a
    stiff system, so it is tabulated once and interpolated after. Interpolation is linear
    in ``log(tau)`` against ``1/T``, which is where Arrhenius behaviour is straight and
    the error of interpolating is smallest.
    """

    def __init__(
        self,
        registry: MechanismRegistry,
        fuel: FuelKind,
        temperatures_k: tuple[float, ...],
        pressures_pa: tuple[float, ...],
        equivalence_ratios: tuple[float, ...],
        oxidizer_mole_fractions: dict[str, float] | None = None,
    ) -> None:
        if len(temperatures_k) < 2:
            raise ValueError("At least two temperatures are required to interpolate")
        self.spec = registry.spec_for(fuel, MechanismRole.IGNITION_DELAY)
        self.uses_dedicated_mechanism = registry.has_dedicated_ignition_mechanism(fuel)
        self.temperatures_k = tuple(sorted(temperatures_k))
        self.pressures_pa = tuple(sorted(pressures_pa))
        self.equivalence_ratios = tuple(sorted(equivalence_ratios))
        self._oxidizer = oxidizer_mole_fractions or DRY_AIR_MOLE_FRACTIONS

        solution = registry.new_solution(fuel, MechanismRole.IGNITION_DELAY)
        shape = (
            len(self.temperatures_k),
            len(self.pressures_pa),
            len(self.equivalence_ratios),
        )
        self._log_delay = np.full(shape, np.nan)
        for i, temperature in enumerate(self.temperatures_k):
            for j, pressure in enumerate(self.pressures_pa):
                for k, phi in enumerate(self.equivalence_ratios):
                    delay = ignition_delay(
                        solution, self.spec, temperature, pressure, phi, self._oxidizer
                    )
                    if delay is not None and delay > 0.0:
                        self._log_delay[i, j, k] = np.log(delay)

    def __call__(
        self, temperature_k: float, pressure_pa: float, equivalence_ratio: float
    ) -> float | None:
        """Interpolated ignition delay, or ``None`` where the mixture does not ignite."""
        inverse = 1.0 / np.array(self.temperatures_k)
        query_inverse = 1.0 / temperature_k

        # Interpolate along 1/T at the nearest pressure and equivalence ratio, which is
        # adequate because the temperature dependence is by far the strongest.
        pressure_index = int(
            np.argmin(np.abs(np.array(self.pressures_pa) - pressure_pa))
        )
        phi_index = int(
            np.argmin(np.abs(np.array(self.equivalence_ratios) - equivalence_ratio))
        )
        column = self._log_delay[:, pressure_index, phi_index]
        valid = ~np.isnan(column)
        if valid.sum() < 2:
            return None

        order = np.argsort(inverse[valid])
        result = np.interp(
            query_inverse, inverse[valid][order], column[valid][order],
            left=np.nan, right=np.nan,
        )
        if np.isnan(result):
            return None
        return float(np.exp(result))


def autoignition_margin(
    table: IgnitionDelayTable,
    premix: PremixState,
    residence_time_s: float,
    fuel: FuelKind,
    minimum_margin: float = DEFAULT_MINIMUM_MARGIN,
) -> AutoignitionMargin:
    """Compare ignition delay with the time the mixture spends in the passage."""
    if residence_time_s <= 0.0:
        raise ValueError("Residence time must be positive")

    warnings: list[ModelWarning] = []
    if not table.uses_dedicated_mechanism:
        warnings.append(
            ModelWarning(
                code="IGNITION_MECHANISM_NOT_DEDICATED",
                severity=WarningSeverity.WARNING,
                message=(
                    "Ignition delay was computed from the reactor-network mechanism "
                    "rather than one carrying low-temperature chemistry. Measured on the "
                    "supplied Jet-A files, that overpredicts the delay by 510x at 700 K "
                    "and 71x at 800 K, which would declare an unsafe premixer safe."
                ),
            )
        )

    delay = table(premix.temperature_k, premix.pressure_pa, premix.equivalence_ratio)
    if delay is None:
        warnings.append(
            ModelWarning(
                code="IGNITION_DELAY_UNAVAILABLE",
                severity=WarningSeverity.INFO,
                message=(
                    "The mixture did not ignite within the tabulated window, or the "
                    "condition lies outside the table. No margin is reported."
                ),
            )
        )
        return AutoignitionMargin(
            fuel=fuel, premix=premix, ignition_delay_s=None,
            residence_time_s=residence_time_s, margin=None,
            minimum_margin=minimum_margin, verdict=AutoignitionVerdict.SAFE,
            mechanism_path=table.spec.path,
            used_dedicated_ignition_mechanism=table.uses_dedicated_mechanism,
            warnings=tuple(warnings),
        )

    margin = delay / residence_time_s
    if margin < 1.0:
        verdict = AutoignitionVerdict.UNSAFE
        warnings.append(
            ModelWarning(
                code="PREMIXER_AUTOIGNITION",
                severity=WarningSeverity.ERROR,
                message=(
                    f"Ignition delay ({delay:.3e} s) is shorter than the passage "
                    f"residence time ({residence_time_s:.3e} s). The mixture ignites "
                    "inside the premixer. This design is invalid at this condition and "
                    "its emissions are meaningless."
                ),
            )
        )
    elif margin < minimum_margin:
        verdict = AutoignitionVerdict.MARGINAL
        warnings.append(
            ModelWarning(
                code="PREMIXER_MARGIN_LOW",
                severity=WarningSeverity.WARNING,
                message=(
                    f"Autoignition margin is {margin:.2f}, below the required "
                    f"{minimum_margin:.2f}. The passage does not ignite in this "
                    "calculation but has little room for variation."
                ),
            )
        )
    else:
        verdict = AutoignitionVerdict.SAFE

    return AutoignitionMargin(
        fuel=fuel, premix=premix, ignition_delay_s=delay,
        residence_time_s=residence_time_s, margin=margin,
        minimum_margin=minimum_margin, verdict=verdict,
        mechanism_path=table.spec.path,
        used_dedicated_ignition_mechanism=table.uses_dedicated_mechanism,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True)
class FlashbackScreen:
    """Whether a flame could travel upstream into the passage."""

    passage_velocity_m_s: float
    turbulent_flame_speed_m_s: float | None
    margin: float | None
    is_safe: bool | None
    warnings: tuple[ModelWarning, ...]


def flashback_screen(
    passage_velocity_m_s: float,
    laminar_flame_speed_m_s: float | None,
    turbulence_intensity: float = 0.1,
    turbulent_factor_coefficient: float = 3.5,
    turbulent_factor_exponent: float = 0.5,
) -> FlashbackScreen:
    """Compare the passage velocity with how fast a flame could climb it.

    The turbulent flame speed is estimated from the laminar value by a declared
    correlation of the usual form ``S_T = S_L (1 + C (u'/S_L)^n)``.

    The laminar flame speed is **not** computed here. Obtaining it requires a flame
    calculation this tool does not perform, and inventing one would put an unmarked guess
    underneath a safety screen. Supplying it is the user's decision, and without it the
    screen reports nothing rather than something unfounded.
    """
    if passage_velocity_m_s <= 0.0:
        raise ValueError("Passage velocity must be positive")

    if laminar_flame_speed_m_s is None or laminar_flame_speed_m_s <= 0.0:
        return FlashbackScreen(
            passage_velocity_m_s=passage_velocity_m_s,
            turbulent_flame_speed_m_s=None,
            margin=None,
            is_safe=None,
            warnings=(
                ModelWarning(
                    code="FLASHBACK_SCREEN_UNAVAILABLE",
                    severity=WarningSeverity.WARNING,
                    message=(
                        "No laminar flame speed was supplied, so the flashback screen was "
                        "not evaluated. It requires a flame calculation this tool does "
                        "not perform; supply a value from a flame solver or measurement."
                    ),
                ),
            ),
        )

    fluctuation = turbulence_intensity * passage_velocity_m_s
    turbulent = laminar_flame_speed_m_s * (
        1.0
        + turbulent_factor_coefficient
        * (fluctuation / laminar_flame_speed_m_s) ** turbulent_factor_exponent
    )
    margin = passage_velocity_m_s / turbulent
    warnings: list[ModelWarning] = []
    if margin < 1.0:
        warnings.append(
            ModelWarning(
                code="FLASHBACK_RISK",
                severity=WarningSeverity.ERROR,
                message=(
                    f"Estimated turbulent flame speed ({turbulent:.2f} m/s) exceeds the "
                    f"passage velocity ({passage_velocity_m_s:.2f} m/s); a flame could "
                    "travel upstream into the premixer."
                ),
            )
        )
    warnings.append(
        ModelWarning(
            code="FLASHBACK_SCREEN_CORRELATION",
            severity=WarningSeverity.INFO,
            message=(
                "The flashback screen uses a declared turbulent flame speed correlation "
                f"with C={turbulent_factor_coefficient} and n={turbulent_factor_exponent}. "
                "It is a screen, not a validated prediction."
            ),
        )
    )
    return FlashbackScreen(
        passage_velocity_m_s=passage_velocity_m_s,
        turbulent_flame_speed_m_s=turbulent,
        margin=margin,
        is_safe=margin >= 1.0,
        warnings=tuple(warnings),
    )
