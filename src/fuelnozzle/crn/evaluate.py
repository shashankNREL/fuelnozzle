"""Evaluate one design across a mission set.

This is the inner loop of every sweep and every search, so it has to be both correct and
quick. Two choices make it quick.

Mechanism objects and ignition-delay tables are cached across evaluations, because
building them dominates the cost of a single point and none of them depend on the design.

Fuel is introduced prevaporized rather than through the full droplet coupling. That is a
real approximation and it is stated: evaporation completes upstream of the quench in the
cases studied, and the prevaporized path was measured to reproduce the fully coupled
solution to 2.7% while running about five times faster. Atomization-quality outputs are
therefore **not** available from this path, and are reported separately by the droplet
models when a spray calibration exists.
"""

from __future__ import annotations

from dataclasses import dataclass

import cantera as ct

from fuelnozzle.crn.autoignition import (
    AutoignitionMargin,
    AutoignitionVerdict,
    IgnitionDelayTable,
    autoignition_margin,
    premix_state,
)
from fuelnozzle.crn.chemistry import (
    FuelKind,
    MechanismRegistry,
    MechanismRole,
    emission_index_g_per_kg,
    equivalence_ratio,
    stoichiometric_air_fuel_ratio,
)
from fuelnozzle.crn.design import DesignVector, MissionPoint
from fuelnozzle.crn.network import CombustorNetwork, NetworkError
from fuelnozzle.crn.reactors import InletSpec, OutletSpec
from fuelnozzle.crn.templates import (
    Architecture,
    ArchitectureInputs,
    check_quench_residence_time,
    ldi_architecture,
    lpp_architecture,
    quench_residence_time_s,
    rql_architecture,
)
from fuelnozzle.models import ModelWarning, WarningSeverity

#: Temperatures at which ignition delay is tabulated once per fuel.
IGNITION_TABLE_TEMPERATURES_K = (650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 1000.0)

ARCHITECTURES = {
    "rql": rql_architecture,
    "ldi": ldi_architecture,
    "lpp": lpp_architecture,
}


@dataclass(frozen=True)
class PointResult:
    """What one design does at one mission point."""

    point: MissionPoint
    architecture: str
    feasible: bool
    exit_temperature_k: float
    peak_temperature_k: float
    ei_nox_g_per_kg: float
    equivalence_ratio_spread: float
    near_field_equivalence_ratio: float
    exit_temperature_spread_k: float
    quench_residence_time_s: float
    autoignition: AutoignitionMargin | None
    warnings: tuple[ModelWarning, ...]

    @property
    def is_extinguished(self) -> bool:
        return any(
            warning.code == "NETWORK_EXTINGUISHED" for warning in self.warnings
        )


@dataclass(frozen=True)
class DesignResult:
    """What one design does across the whole mission."""

    design: DesignVector
    points: tuple[PointResult, ...]
    feasible: bool
    warnings: tuple[ModelWarning, ...]

    def by_fuel(self, fuel: FuelKind) -> tuple[PointResult, ...]:
        return tuple(result for result in self.points if result.point.fuel is fuel)

    def weighted_ei_nox(self, fuel: FuelKind) -> float:
        """Time-weighted emission index for one fuel, or the plain mean without times."""
        results = self.by_fuel(fuel)
        if not results:
            return 0.0
        weights = [max(result.point.duration_s, 0.0) for result in results]
        if sum(weights) <= 0.0:
            return sum(result.ei_nox_g_per_kg for result in results) / len(results)
        return sum(
            result.ei_nox_g_per_kg * weight for result, weight in zip(results, weights, strict=True)
        ) / sum(weights)


class DesignEvaluator:
    """Evaluates designs against a fixed mission set and mechanism registry."""

    def __init__(
        self,
        registry: MechanismRegistry,
        mission: tuple[MissionPoint, ...] | list[MissionPoint],
        architecture: str = "rql",
        lng_architecture: str | None = None,
        minimum_autoignition_margin: float = 4.0,
    ) -> None:
        if architecture not in ARCHITECTURES:
            raise ValueError(f"Unknown architecture {architecture!r}")
        if lng_architecture is not None and lng_architecture not in ARCHITECTURES:
            raise ValueError(f"Unknown architecture {lng_architecture!r}")
        self.registry = registry
        self.mission = tuple(mission)
        self.architecture = architecture
        # The two fuel paths are separate hardware, so they may use different
        # architectures while sharing one liner.
        self.lng_architecture = lng_architecture or architecture
        self.minimum_margin = minimum_autoignition_margin
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
        inputs = ArchitectureInputs(
            fuel=point.fuel,
            fuel_mass_flow_kg_s=point.fuel_mass_flow_kg_s,
            total_air_mass_flow_kg_s=point.air_mass_flow_kg_s,
            air_temperature_k=point.air_temperature_k,
            air_split=design.air_split(point.fuel),
            stoichiometric_air_fuel_ratio=self.stoichiometric_afr(point.fuel),
            quench_volume_m3=design.quench_volume_m3,
            flame_volume_m3=design.flame_volume_m3,
            post_volume_m3=design.post_volume_m3,
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

        inlets = list(architecture.air_inlets) + [
            InletSpec(
                name="fuel_vapor",
                target_reactor=architecture.spray_path[0],
                mass_flow_kg_s=point.fuel_mass_flow_kg_s,
                temperature_k=design.fuel_temperature_k(point.fuel),
                mole_fractions=spec.fuel_mole_fractions,
            )
        ]
        total = point.air_mass_flow_kg_s + point.fuel_mass_flow_kg_s

        try:
            network = CombustorNetwork(
                architecture.reactors, inlets,
                OutletSpec(source_reactor=architecture.outlet_reactor, mass_flow_kg_s=total),
                architecture.internal_flows,
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
                point=point, architecture=self.architecture_for(point.fuel), feasible=False,
                exit_temperature_k=float("nan"), peak_temperature_k=float("nan"),
                ei_nox_g_per_kg=float("inf"), equivalence_ratio_spread=float("nan"),
                near_field_equivalence_ratio=architecture.near_field_equivalence_ratio,
                exit_temperature_spread_k=float("nan"),
                quench_residence_time_s=float("nan"), autoignition=None,
                warnings=tuple(warnings),
            )

        warnings.extend(solution.warnings)
        warnings.extend(check_quench_residence_time(solution))

        gas = mechanism()
        outlet = solution.outlet
        gas.TPY = outlet.temperature_k, point.pressure_pa, outlet.mass_fractions
        ei_nox = emission_index_g_per_kg(
            float(gas.Y[gas.species_index("NO")]) * 46.0055 / 30.0061
            + (
                float(gas.Y[gas.species_index("NO2")])
                if "NO2" in gas.species_names
                else 0.0
            ),
            total, point.fuel_mass_flow_kg_s,
        )

        ratios = []
        for reactor in solution.reactors:
            gas.TPY = reactor.temperature_k, point.pressure_pa, reactor.mass_fractions
            ratios.append(equivalence_ratio(gas, spec))
        temperatures = [reactor.temperature_k for reactor in solution.reactors]

        ignition = self._autoignition(design, point, architecture)
        if ignition is not None:
            warnings.extend(ignition.warnings)

        feasible = (
            not any(w.severity is WarningSeverity.ERROR for w in warnings)
            and solution.converged
        )
        return PointResult(
            point=point,
            architecture=self.architecture_for(point.fuel),
            feasible=feasible,
            exit_temperature_k=outlet.temperature_k,
            peak_temperature_k=solution.peak_temperature_k,
            ei_nox_g_per_kg=ei_nox,
            equivalence_ratio_spread=max(ratios) - min(ratios),
            near_field_equivalence_ratio=architecture.near_field_equivalence_ratio,
            exit_temperature_spread_k=max(temperatures) - min(temperatures),
            quench_residence_time_s=quench_residence_time_s(solution),
            autoignition=ignition,
            warnings=tuple(warnings),
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
        return DesignResult(
            design=design,
            points=results,
            feasible=all(result.feasible for result in results),
            warnings=warnings,
        )
