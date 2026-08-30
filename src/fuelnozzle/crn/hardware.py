"""Shared-liner geometry, pressure stations, and area-derived air admission."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import sqrt

from fuelnozzle.crn.streams import AirSplit, CoolingAirDestination

AIR_GAMMA = 1.4
AIR_GAS_CONSTANT_J_KG_K = 287.05


class AirflowMode(StrEnum):
    """Whether air fractions come from hardware or a documented calibration."""

    AREA_DERIVED = "area_derived"
    PRESCRIBED_CALIBRATION = "prescribed_calibration"


@dataclass(frozen=True)
class SectorDefinition:
    """How much of a complete engine the modeled sector represents."""

    cups_per_engine: int
    modeled_cups: int = 1

    def __post_init__(self) -> None:
        if self.cups_per_engine <= 0 or not 0 < self.modeled_cups <= self.cups_per_engine:
            raise ValueError("Modeled cups must lie between one and cups per engine")

    @property
    def engine_fraction(self) -> float:
        return self.modeled_cups / self.cups_per_engine

    def from_engine_total(self, value: float) -> float:
        return value * self.engine_fraction

    def to_engine_total(self, value: float) -> float:
        return value / self.engine_fraction


@dataclass(frozen=True)
class AirAdmission:
    """One set of liner holes represented by an effective flow area."""

    name: str
    effective_area_m2: float
    discharge_coefficient: float

    def __post_init__(self) -> None:
        if self.effective_area_m2 <= 0.0:
            raise ValueError("Air-admission area must be positive")
        if not 0.0 < self.discharge_coefficient <= 1.0:
            raise ValueError("Air-admission discharge coefficient must lie in (0, 1]")

    def mass_flow_kg_s(
        self,
        upstream_pressure_pa: float,
        downstream_pressure_pa: float,
        upstream_temperature_k: float,
    ) -> float:
        """Compressible ideal-gas orifice flow, including the choked limit."""
        if not 0.0 < downstream_pressure_pa < upstream_pressure_pa:
            raise ValueError("Air admission requires 0 < downstream pressure < upstream pressure")
        if upstream_temperature_k <= 0.0:
            raise ValueError("Upstream temperature must be positive")
        pressure_ratio = downstream_pressure_pa / upstream_pressure_pa
        critical_ratio = (2.0 / (AIR_GAMMA + 1.0)) ** (
            AIR_GAMMA / (AIR_GAMMA - 1.0)
        )
        scale = (
            self.discharge_coefficient
            * self.effective_area_m2
            * upstream_pressure_pa
            / sqrt(AIR_GAS_CONSTANT_J_KG_K * upstream_temperature_k)
        )
        if pressure_ratio <= critical_ratio:
            flux_factor = sqrt(AIR_GAMMA) * (
                2.0 / (AIR_GAMMA + 1.0)
            ) ** ((AIR_GAMMA + 1.0) / (2.0 * (AIR_GAMMA - 1.0)))
        else:
            flux_factor = sqrt(
                2.0
                * AIR_GAMMA
                / (AIR_GAMMA - 1.0)
                * (
                    pressure_ratio ** (2.0 / AIR_GAMMA)
                    - pressure_ratio ** ((AIR_GAMMA + 1.0) / AIR_GAMMA)
                )
            )
        return scale * flux_factor


@dataclass(frozen=True)
class FuelPassageGeometry:
    """Packaging and flow area for one fuel injector passage."""

    fuel: str
    air_effective_area_m2: float
    envelope_area_m2: float

    def __post_init__(self) -> None:
        if self.air_effective_area_m2 <= 0.0 or self.envelope_area_m2 <= 0.0:
            raise ValueError("Fuel-passage areas must be positive")


@dataclass(frozen=True)
class SharedLinerGeometry:
    """One immutable liner used by both fuel systems."""

    quench_volume_m3: float
    flame_volume_m3: float
    post_volume_m3: float
    dome: AirAdmission
    primary: AirAdmission
    quench: AirAdmission
    dilution: AirAdmission
    cooling: AirAdmission
    cooling_destination: CoolingAirDestination = CoolingAirDestination.DILUTION

    def __post_init__(self) -> None:
        if min(self.zone_volumes_m3) <= 0.0:
            raise ValueError("Every shared-liner zone volume must be positive")

    @property
    def zone_volumes_m3(self) -> tuple[float, float, float]:
        return (self.quench_volume_m3, self.flame_volume_m3, self.post_volume_m3)

    @property
    def volume_m3(self) -> float:
        return sum(self.zone_volumes_m3)

    def area_derived_split(
        self,
        upstream_pressure_pa: float,
        downstream_pressure_pa: float,
        temperature_k: float,
        *,
        jet_a_passage_share: float,
        idle_passage_mixing_fraction: float,
    ) -> AirSplit:
        """Normalize station flows computed from fixed effective areas."""
        admissions = (self.dome, self.primary, self.quench, self.dilution, self.cooling)
        flows = [
            admission.mass_flow_kg_s(
                upstream_pressure_pa, downstream_pressure_pa, temperature_k
            )
            for admission in admissions
        ]
        total = sum(flows)
        if total <= 0.0:
            raise ValueError("Shared liner admits no air")
        fractions = [flow / total for flow in flows]
        return AirSplit(
            dome=fractions[0],
            primary=fractions[1],
            quench=fractions[2],
            dilution=fractions[3],
            cooling=fractions[4],
            cooling_destination=self.cooling_destination,
            jet_a_passage_share=jet_a_passage_share,
            idle_passage_mixing_fraction=idle_passage_mixing_fraction,
        )


@dataclass(frozen=True)
class DualFuelHardware:
    """Separate fuel passages installed in one shared physical liner."""

    liner: SharedLinerGeometry
    jet_a_passage: FuelPassageGeometry
    lng_passage: FuelPassageGeometry
    sector: SectorDefinition
    airflow_mode: AirflowMode = AirflowMode.AREA_DERIVED
    calibration_id: str | None = None

    def __post_init__(self) -> None:
        if (
            self.airflow_mode is AirflowMode.PRESCRIBED_CALIBRATION
            and not self.calibration_id
        ):
            raise ValueError("Prescribed airflow requires a calibration identifier")

    @property
    def jet_a_passage_share(self) -> float:
        """Dome-flow share implied by the two fixed passage effective areas."""
        total = (
            self.jet_a_passage.air_effective_area_m2
            + self.lng_passage.air_effective_area_m2
        )
        return self.jet_a_passage.air_effective_area_m2 / total
