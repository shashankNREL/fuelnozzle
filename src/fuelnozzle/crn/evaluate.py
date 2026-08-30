"""Physics-gated evaluation of one design across a mission set.

This is the inner loop of every sweep and every search, so it has to be both correct and
quick. Two choices make it quick.

Mechanism objects and ignition-delay tables are cached across evaluations, because
building them dominates the cost of a single point and none of them depend on the design.

An explicit spray-model callback can supply the nozzle result and liquid properties for
conservative gas/liquid coupling. Omitting it retains a fast prevaporized diagnostic, but the
acceptance gate then remains UNKNOWN.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import cantera as ct

from fuelnozzle.crn.autoignition import (
    AutoignitionMargin,
    AutoignitionVerdict,
    FlashbackScreen,
    IgnitionDelayTable,
    autoignition_margin,
    premix_state,
)
from fuelnozzle.crn.chemistry import (
    FuelKind,
    MechanismRegistry,
    MechanismRole,
    emission_index_g_per_kg,
    stoichiometric_air_fuel_ratio,
    validate_mechanism,
)
from fuelnozzle.crn.coupling import CoupledSolution, solve_coupled
from fuelnozzle.crn.design import DesignVector, MissionPoint
from fuelnozzle.crn.droplets import LiquidPropertyProvider
from fuelnozzle.crn.emissions import LTOMode, LTOResult, lto_dp_foo
from fuelnozzle.crn.hardware import AirflowMode, DualFuelHardware
from fuelnozzle.crn.mission import MissionProfile, mission_point_from_operating
from fuelnozzle.crn.network import CombustorNetwork, NetworkError
from fuelnozzle.crn.reactors import InletSpec, OutletSpec
from fuelnozzle.crn.spray_source import SprayBoundary
from fuelnozzle.crn.status import GateResult, GateStatus, aggregate_gate_status
from fuelnozzle.crn.templates import (
    Architecture,
    ArchitectureInputs,
    check_quench_residence_time,
    ldi_architecture,
    lpp_architecture,
    quench_residence_time_s,
    rql_architecture,
)
from fuelnozzle.crn.thermal import IdleCircuitScreen, idle_circuit_screen
from fuelnozzle.models import ModelWarning, WarningSeverity
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider

#: Temperatures at which ignition delay is tabulated once per fuel.
IGNITION_TABLE_TEMPERATURES_K = (650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 1000.0)

ARCHITECTURES = {
    "rql": rql_architecture,
    "ldi": ldi_architecture,
    "lpp": lpp_architecture,
}


@dataclass(frozen=True)
class PointResult:
    """Prototype screening result at one mission point."""

    point: MissionPoint
    architecture: str
    exit_temperature_k: float
    peak_temperature_k: float
    ei_nox_g_per_kg: float
    equivalence_ratio_spread: float
    near_field_equivalence_ratio: float
    exit_temperature_spread_k: float
    quench_residence_time_s: float
    autoignition: AutoignitionMargin | None
    warnings: tuple[ModelWarning, ...]
    computational_status: GateStatus
    acceptance_status: GateStatus
    gates: tuple[GateResult, ...]
    mechanism_path: str
    mechanism_provenance: str
    coupled_spray: CoupledSolution | None = None
    idle_circuit: IdleCircuitScreen | None = None
    flashback: FlashbackScreen | None = None

    @property
    def feasible(self) -> bool:
        """Compatibility view: only an accepted result is feasible."""
        return self.acceptance_status is GateStatus.PASS

    @property
    def is_extinguished(self) -> bool:
        return any(
            warning.code == "NETWORK_EXTINGUISHED" for warning in self.warnings
        )


@dataclass(frozen=True)
class DesignResult:
    """Prototype screening result across the whole mission."""

    design: DesignVector
    points: tuple[PointResult, ...]
    warnings: tuple[ModelWarning, ...]
    computational_status: GateStatus
    acceptance_status: GateStatus
    gates: tuple[GateResult, ...]
    rated_thrust_kn: float | None = None

    @property
    def feasible(self) -> bool:
        """Compatibility view: only an accepted design is feasible."""
        return self.acceptance_status is GateStatus.PASS

    @property
    def exit_gate_satisfied(self) -> bool:
        """Stage-5/7 credibility gate; missing evidence never becomes success."""
        return self.acceptance_status is GateStatus.PASS and all(
            point.acceptance_status is GateStatus.PASS for point in self.points
        )

    def by_fuel(self, fuel: FuelKind) -> tuple[PointResult, ...]:
        return tuple(result for result in self.points if result.point.fuel is fuel)

    def weighted_ei_nox(self, fuel: FuelKind) -> float:
        """Fuel-mass-weighted emission index for one fuel."""
        results = self.by_fuel(fuel)
        if not results:
            return 0.0
        weights = [
            max(result.point.duration_s, 0.0)
            * max(result.point.fuel_mass_flow_kg_s, 0.0)
            for result in results
        ]
        if sum(weights) <= 0.0:
            return sum(result.ei_nox_g_per_kg for result in results) / len(results)
        return sum(
            result.ei_nox_g_per_kg * weight for result, weight in zip(results, weights, strict=True)
        ) / sum(weights)

    def lto_emissions(self) -> LTOResult | None:
        """Assemble the four-mode Jet-A Dp/Foo result when rated thrust is known."""
        if self.rated_thrust_kn is None:
            return None
        modes = tuple(
            LTOMode(
                name=result.point.name,
                thrust_fraction=result.point.thrust_fraction,
                duration_s=result.point.duration_s,
                fuel_mass_flow_kg_s=result.point.fuel_mass_flow_kg_s,
                ei_nox_g_per_kg=result.ei_nox_g_per_kg,
            )
            for result in self.by_fuel(FuelKind.JET_A)
        )
        return lto_dp_foo(modes, self.rated_thrust_kn) if modes else None

    def lng_cruise_ei_by_point(self) -> dict[str, float]:
        """Keep named cruise results visible rather than hiding them in one average."""
        return {
            result.point.name: result.ei_nox_g_per_kg
            for result in self.by_fuel(FuelKind.LNG)
        }


SprayModel = Callable[
    [MissionPoint, DesignVector],
    tuple[SprayBoundary, LiquidPropertyProvider],
]
FlashbackModel = Callable[[MissionPoint, DesignVector], FlashbackScreen]


class DesignEvaluator:
    """Runs a gated gas-only or gas/liquid model over a fixed mission set."""

    def __init__(
        self,
        registry: MechanismRegistry,
        mission: (
            MissionProfile
            | tuple[MissionPoint | OperatingPoint, ...]
            | list[MissionPoint | OperatingPoint]
        ),
        architecture: str = "rql",
        lng_architecture: str | None = None,
        minimum_autoignition_margin: float = 4.0,
        hardware: DualFuelHardware | None = None,
        spray_model: SprayModel | None = None,
        rated_thrust_kn: float | None = None,
        lng_properties: CoolPropLNGProvider | None = None,
        flashback_model: FlashbackModel | None = None,
    ) -> None:
        if architecture not in ARCHITECTURES:
            raise ValueError(f"Unknown architecture {architecture!r}")
        if lng_architecture is not None and lng_architecture not in ARCHITECTURES:
            raise ValueError(f"Unknown architecture {lng_architecture!r}")
        source_points = mission.points if isinstance(mission, MissionProfile) else mission
        self.registry = registry
        self.mission = tuple(
            mission_point_from_operating(point, hardware.sector if hardware else None)
            if isinstance(point, OperatingPoint)
            else point
            for point in source_points
        )
        self.architecture = architecture
        # The two fuel paths are separate hardware, so they may use different
        # architectures while sharing one liner.
        self.lng_architecture = lng_architecture or architecture
        self.minimum_margin = minimum_autoignition_margin
        self.hardware = hardware
        self.spray_model = spray_model
        if rated_thrust_kn is not None and rated_thrust_kn <= 0.0:
            raise ValueError("Rated thrust must be positive")
        self.rated_thrust_kn = rated_thrust_kn
        self.lng_properties = lng_properties
        self.flashback_model = flashback_model
        self._tables: dict[tuple[FuelKind, float], IgnitionDelayTable] = {}
        self._afr: dict[FuelKind, float] = {}

    def architecture_for(self, fuel: FuelKind) -> str:
        return self.lng_architecture if fuel is FuelKind.LNG else self.architecture

    def stoichiometric_afr(self, fuel: FuelKind) -> float:
        if fuel not in self._afr:
            self._afr[fuel] = stoichiometric_air_fuel_ratio(
                self.registry.template(fuel, MechanismRole.NETWORK),
                self.registry.spec_for(fuel, MechanismRole.NETWORK),
            )
        return self._afr[fuel]

    def ignition_table(self, fuel: FuelKind, pressure_pa: float) -> IgnitionDelayTable:
        key = (fuel, round(pressure_pa, -3))
        if key not in self._tables:
            self._tables[key] = IgnitionDelayTable(
                self.registry, fuel, IGNITION_TABLE_TEMPERATURES_K, (pressure_pa,),
                (0.4, 0.6, 0.8, 1.2, 1.6),
            )
        return self._tables[key]

    def _build(self, design: DesignVector, point: MissionPoint) -> Architecture:
        split = design.air_split(point.fuel)
        zone_volumes = (
            design.quench_volume_m3,
            design.flame_volume_m3,
            design.post_volume_m3,
        )
        if self.hardware is not None:
            zone_volumes = self.hardware.liner.zone_volumes_m3
            split = split.model_copy(
                update={
                    "cooling_destination": self.hardware.liner.cooling_destination,
                    "jet_a_passage_share": self.hardware.jet_a_passage_share,
                }
            )
        if (
            self.hardware is not None
            and self.hardware.airflow_mode is AirflowMode.AREA_DERIVED
        ):
            stations = point.pressure_stations
            if stations is None or stations.liner_pressure_loss_pa <= 0.0:
                raise ValueError(
                    "Area-derived hardware requires a canonical operating point with "
                    "positive liner pressure loss"
                )
            split = self.hardware.liner.area_derived_split(
                stations.compressor_discharge_pa,
                stations.combustor_exit_pa,
                point.air_temperature_k,
                jet_a_passage_share=self.hardware.jet_a_passage_share,
                idle_passage_mixing_fraction=design.idle_passage_mixing_fraction,
            )
        inputs = ArchitectureInputs(
            fuel=point.fuel,
            fuel_mass_flow_kg_s=point.fuel_mass_flow_kg_s,
            total_air_mass_flow_kg_s=point.air_mass_flow_kg_s,
            air_temperature_k=point.air_temperature_k,
            air_split=split,
            stoichiometric_air_fuel_ratio=self.stoichiometric_afr(point.fuel),
            quench_volume_m3=zone_volumes[0],
            flame_volume_m3=zone_volumes[1],
            post_volume_m3=zone_volumes[2],
        )
        builder = ARCHITECTURES[self.architecture_for(point.fuel)]
        if builder is rql_architecture:
            return builder(inputs, quench_stages=design.quench_stages)
        return builder(inputs)

    def evaluate_point(self, design: DesignVector, point: MissionPoint) -> PointResult:
        architecture = self._build(design, point)
        warnings = list(architecture.warnings) + list(design.packaging_warnings)
        spec = self.registry.spec_for(point.fuel, MechanismRole.NETWORK)

        def mechanism() -> ct.Solution:
            return self.registry.new_solution(point.fuel, MechanismRole.NETWORK)

        warnings.extend(
            validate_mechanism(
                spec,
                mechanism(),
                require_nox=True,
                pressure_pa=point.pressure_pa,
                temperature_k=point.air_temperature_k,
            )
        )
        total = point.air_mass_flow_kg_s + point.fuel_mass_flow_kg_s
        coupled: CoupledSolution | None = None

        try:
            if self.spray_model is not None:
                spray, liquid_provider = self.spray_model(point, design)
                if spray.fuel is not point.fuel:
                    raise ValueError("Spray boundary fuel does not match the mission point")
                mismatch = abs(
                    spray.total_fuel_mass_flow_kg_s - point.fuel_mass_flow_kg_s
                )
                if mismatch > 1.0e-9 * point.fuel_mass_flow_kg_s:
                    raise ValueError(
                        "Spray boundary fuel flow does not match the mission point"
                    )
                coupled = solve_coupled(
                    architecture.reactors,
                    architecture.air_inlets,
                    architecture.outlet_reactor,
                    architecture.internal_flows,
                    mechanism,
                    point.pressure_pa,
                    spray,
                    liquid_provider,
                    architecture.spray_path,
                    spec.fuel_mole_fractions,
                    fixed_internal_flows=architecture.fixed_internal_flows,
                )
                solution = coupled.network
                warnings.extend(coupled.warnings)
            else:
                inlets = list(architecture.air_inlets) + [
                    InletSpec(
                        name="fuel_vapor",
                        target_reactor=architecture.spray_path[0],
                        mass_flow_kg_s=point.fuel_mass_flow_kg_s,
                        temperature_k=design.fuel_temperature_k(point.fuel),
                        mole_fractions=spec.fuel_mole_fractions,
                    )
                ]
                network = CombustorNetwork(
                    architecture.reactors,
                    inlets,
                    OutletSpec(
                        source_reactor=architecture.outlet_reactor,
                        mass_flow_kg_s=total,
                    ),
                    architecture.internal_flows,
                    fixed_internal_flows=architecture.fixed_internal_flows,
                )
                solution = network.solve(mechanism, point.pressure_pa)
        except (NetworkError, ValueError) as error:
            warnings.append(
                ModelWarning(
                    code="DESIGN_EVALUATION_FAILED",
                    severity=WarningSeverity.ERROR,
                    message=f"{point.name}: {error}",
                )
            )
            return PointResult(
                point=point, architecture=self.architecture_for(point.fuel),
                exit_temperature_k=float("nan"), peak_temperature_k=float("nan"),
                ei_nox_g_per_kg=float("inf"), equivalence_ratio_spread=float("nan"),
                near_field_equivalence_ratio=architecture.near_field_equivalence_ratio,
                exit_temperature_spread_k=float("nan"),
                quench_residence_time_s=float("nan"), autoignition=None,
                warnings=tuple(warnings),
                computational_status=GateStatus.FAIL,
                acceptance_status=GateStatus.FAIL,
                gates=(
                    GateResult(
                        name="network_solve",
                        status=GateStatus.FAIL,
                        reason=str(error),
                        evidence="CombustorNetwork.solve",
                    ),
                ),
                mechanism_path=spec.path,
                mechanism_provenance=spec.provenance,
            )

        warnings.extend(solution.warnings)
        warnings.extend(check_quench_residence_time(solution))

        gas = mechanism()
        outlet = solution.outlet
        gas.TPY = outlet.temperature_k, point.pressure_pa, outlet.mass_fractions
        exhaust_mass_flow = total - (
            coupled.liquid_carryover_kg_s if coupled is not None else 0.0
        )
        ei_nox = emission_index_g_per_kg(
            float(gas.Y[gas.species_index("NO")]) * 46.0055 / 30.0061
            + (
                float(gas.Y[gas.species_index("NO2")])
                if "NO2" in gas.species_names
                else 0.0
            ),
            exhaust_mass_flow, point.fuel_mass_flow_kg_s,
        )

        ignition = self._autoignition(design, point, architecture)
        if ignition is not None:
            warnings.extend(ignition.warnings)
        idle_circuit: IdleCircuitScreen | None = None
        if point.nozzle_wall_temperature_k is not None:
            idle_circuit = idle_circuit_screen(
                point.fuel.value,
                point.nozzle_wall_temperature_k,
                provider=self.lng_properties,
            )
            warnings.extend(idle_circuit.warnings)
        flashback = (
            self.flashback_model(point, design)
            if self.flashback_model is not None
            else None
        )
        if flashback is not None:
            warnings.extend(flashback.warnings)

        computational_status = (
            GateStatus.PASS
            if solution.converged and (coupled is None or coupled.converged)
            else GateStatus.FAIL
        )
        gates = [
            GateResult(
                name="network_solve",
                status=computational_status,
                reason=(
                    "network converged without an error"
                    if computational_status is GateStatus.PASS
                    else "network did not converge or emitted an error"
                ),
                evidence="CombustorNetwork.solve",
            ),
            GateResult(
                name="model_fidelity",
                status=GateStatus.UNKNOWN,
                reason=(
                    "no nozzle/spray callback was supplied; fuel is prevaporized"
                    if coupled is None
                    else "coupled spray physics lacks the required rig validation"
                ),
                evidence="docs/CRN_EQUATION_REGISTER.md",
            ),
            GateResult(
                name="mechanism_validation",
                status=GateStatus.UNKNOWN,
                reason=(
                    "mechanism structure and declared range were checked, but no held-out "
                    "high-pressure emissions validation is registered"
                ),
                evidence=f"{spec.path}: {spec.provenance}",
            ),
            GateResult(
                name="mechanism_applicability",
                status=(
                    GateStatus.UNKNOWN
                    if all(
                        bound is None
                        for bound in (
                            spec.min_pressure_pa,
                            spec.max_pressure_pa,
                            spec.min_temperature_k,
                            spec.max_temperature_k,
                        )
                    )
                    or (
                        (
                            spec.min_pressure_pa is not None
                            and point.pressure_pa < spec.min_pressure_pa
                        )
                        or (
                            spec.max_pressure_pa is not None
                            and point.pressure_pa > spec.max_pressure_pa
                        )
                        or (
                            spec.min_temperature_k is not None
                            and point.air_temperature_k < spec.min_temperature_k
                        )
                        or (
                            spec.max_temperature_k is not None
                            and point.air_temperature_k > spec.max_temperature_k
                        )
                    )
                    else GateStatus.PASS
                ),
                reason=(
                    "mechanism state is checked against its declared "
                    "pressure/temperature domain"
                ),
                evidence=f"{spec.path}: {spec.provenance}",
            ),
            GateResult(
                name="canonical_hardware",
                status=(
                    GateStatus.PASS
                    if self.hardware is not None and point.operating_point is not None
                    else GateStatus.UNKNOWN
                ),
                reason=(
                    "mission point and shared hardware are explicitly defined"
                    if self.hardware is not None and point.operating_point is not None
                    else "legacy mission input or shared hardware definition is missing"
                ),
                evidence="fuelnozzle.crn.mission and fuelnozzle.crn.hardware",
            ),
        ]
        if coupled is not None:
            gates.extend(
                [
                    GateResult(
                        name="spray_mass_closure",
                        status=(
                            GateStatus.PASS
                            if abs(coupled.fuel_mass_residual_kg_s)
                            <= 1.0e-9 * point.fuel_mass_flow_kg_s
                            else GateStatus.FAIL
                        ),
                        reason=(
                            "nozzle vapor, evaporated fuel, and carryover close the "
                            "fuel-mass ledger"
                        ),
                        evidence="fuelnozzle.crn.coupling.solve_coupled",
                    ),
                    GateResult(
                        name="liquid_carryover",
                        status=(
                            GateStatus.PASS
                            if coupled.liquid_carryover_kg_s
                            <= 1.0e-6 * point.fuel_mass_flow_kg_s
                            else GateStatus.FAIL
                        ),
                        reason=(
                            f"liquid carryover is {coupled.liquid_carryover_kg_s:.4g} kg/s"
                        ),
                        evidence="fuelnozzle.crn.coupling.solve_coupled",
                    ),
                    GateResult(
                        name="spray_validation",
                        status=GateStatus.UNKNOWN,
                        reason=(
                            "spray size, breakup, impingement, and multicomponent "
                            "evaporation do not have full-range calibration"
                        ),
                        evidence="docs/CRN_EQUATION_REGISTER.md",
                    ),
                ]
            )
        gates.append(
            GateResult(
                name="idle_circuit",
                status=(
                    idle_circuit.acceptance_status
                    if idle_circuit is not None
                    else GateStatus.UNKNOWN
                ),
                reason=(
                    "idle circuit steady and transient screen evaluated"
                    if idle_circuit is not None
                    else "no nozzle wall temperature or idle-circuit transient evidence"
                ),
                evidence="fuelnozzle.crn.thermal.idle_circuit_screen",
            )
        )
        if ignition is not None and ignition.verdict is not AutoignitionVerdict.NO_PREMIXER:
            ignition_status = {
                AutoignitionVerdict.SAFE: GateStatus.PASS,
                AutoignitionVerdict.MARGINAL: GateStatus.FAIL,
                AutoignitionVerdict.UNSAFE: GateStatus.FAIL,
                AutoignitionVerdict.UNKNOWN: GateStatus.UNKNOWN,
            }[ignition.verdict]
            gates.append(
                GateResult(
                    name="autoignition",
                    status=ignition_status,
                    reason=f"autoignition verdict is {ignition.verdict.value}",
                    evidence=ignition.mechanism_path,
                )
            )
            gates.append(
                GateResult(
                    name="ignition_mechanism_bracket",
                    status=GateStatus.UNKNOWN,
                    reason=(
                        "only one dedicated ignition mechanism is available; the lower "
                        "confidence bound across defensible Jet-A/LNG mechanisms is unknown"
                    ),
                    evidence="mech/README.md",
                )
            )
            gates.append(
                GateResult(
                    name="flashback",
                    status=(
                        GateStatus.UNKNOWN
                        if flashback is None or flashback.is_safe is None
                        else GateStatus.FAIL
                        if not flashback.is_safe
                        else GateStatus.PASS
                        if flashback.calibration_id
                        else GateStatus.UNKNOWN
                    ),
                    reason=(
                        "no passage flame-speed screen was supplied"
                        if flashback is None
                        else "flashback correlation is unsafe or lacks passage calibration"
                    ),
                    evidence=(
                        flashback.calibration_id if flashback is not None else None
                    ),
                )
            )
        gates.extend(
            GateResult(
                name=name,
                status=GateStatus.UNKNOWN,
                reason="external gate: no validated model or hardware-matched rig evidence",
                evidence="docs/CRN_IMPLEMENTATION_LOG.md",
            )
            for name in (
                "lean_blowout",
                "transient_ignition",
                "relight",
                "fuel_switching",
                "thermoacoustics",
            )
        )
        acceptance_status = aggregate_gate_status(gates)
        return PointResult(
            point=point,
            architecture=self.architecture_for(point.fuel),
            exit_temperature_k=outlet.temperature_k,
            peak_temperature_k=solution.peak_temperature_k,
            ei_nox_g_per_kg=ei_nox,
            equivalence_ratio_spread=float("nan"),
            near_field_equivalence_ratio=architecture.near_field_equivalence_ratio,
            exit_temperature_spread_k=float("nan"),
            quench_residence_time_s=quench_residence_time_s(solution),
            autoignition=ignition,
            warnings=tuple(warnings),
            computational_status=computational_status,
            acceptance_status=acceptance_status,
            gates=tuple(gates),
            mechanism_path=spec.path,
            mechanism_provenance=spec.provenance,
            coupled_spray=coupled,
            idle_circuit=idle_circuit,
            flashback=flashback,
        )

    def _autoignition(
        self, design: DesignVector, point: MissionPoint, architecture: Architecture
    ) -> AutoignitionMargin | None:
        """Screen the premixing passage, if this path has one.

        A path with no premixing passage has nowhere for autoignition to happen, which is
        the whole reason lean direct injection exists. That is reported as its own verdict
        rather than as a very large margin.
        """
        residence = design.premix_residence_s(point.fuel)
        spec = self.registry.spec_for(point.fuel, MechanismRole.NETWORK)
        state = premix_state(
            self.registry.new_solution(point.fuel, MechanismRole.NETWORK), spec,
            air_mass_flow_kg_s=point.air_mass_flow_kg_s * design.dome_air_fraction,
            air_temperature_k=point.air_temperature_k,
            fuel_mass_flow_kg_s=point.fuel_mass_flow_kg_s,
            fuel_temperature_k=design.fuel_temperature_k(point.fuel),
            pressure_pa=point.pressure_pa,
        )
        if residence <= 0.0:
            return AutoignitionMargin(
                fuel=point.fuel, premix=state, ignition_delay_s=None,
                residence_time_s=0.0, margin=None,
                minimum_margin=self.minimum_margin,
                verdict=AutoignitionVerdict.NO_PREMIXER,
                mechanism_path=spec.path, used_dedicated_ignition_mechanism=False,
                warnings=(),
            )
        table = self.ignition_table(point.fuel, point.pressure_pa)
        return autoignition_margin(
            table, state, residence, point.fuel, minimum_margin=self.minimum_margin
        )

    def evaluate(self, design: DesignVector) -> DesignResult:
        results = tuple(self.evaluate_point(design, point) for point in self.mission)
        warnings = tuple(
            warning for result in results for warning in result.warnings
        )
        computational_gates = tuple(
            GateResult(
                name=f"{result.point.name}:computation",
                status=result.computational_status,
                reason="mission-point computational status",
            )
            for result in results
        )
        acceptance_gates = tuple(
            GateResult(
                name=f"{result.point.name}:acceptance",
                status=result.acceptance_status,
                reason="mission-point acceptance status",
            )
            for result in results
        )
        return DesignResult(
            design=design,
            points=results,
            warnings=warnings,
            computational_status=aggregate_gate_status(computational_gates),
            acceptance_status=aggregate_gate_status(acceptance_gates),
            gates=acceptance_gates,
            rated_thrust_kn=self.rated_thrust_kn,
        )
