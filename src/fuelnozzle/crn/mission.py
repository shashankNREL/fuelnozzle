"""Adapters from canonical operating points to CRN mission inputs."""

from __future__ import annotations

from dataclasses import dataclass

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.design import MissionPoint
from fuelnozzle.crn.hardware import SectorDefinition
from fuelnozzle.operating import OperatingPoint, resolve_pressure_stations

ICAO_LTO_NAMES = ("takeoff", "climb_out", "approach", "idle")


def mission_point_from_operating(
    point: OperatingPoint, sector: SectorDefinition | None = None
) -> MissionPoint:
    """Build a CRN mission point without duplicating operating conditions by hand."""
    active = point.active_fuel
    if active is None:
        raise ValueError("A CRN operating point must have exactly one active fuel")
    if point.combustor_air_mass_flow_kg_s is None:
        raise ValueError("A CRN operating point requires combustor_air_mass_flow_kg_s")
    stations = resolve_pressure_stations(point)
    fuel = FuelKind.JET_A if active == "jet_a" else FuelKind.LNG
    fuel_flow = (
        point.jet_a_mass_flow_kg_s if fuel is FuelKind.JET_A else point.lng_mass_flow_kg_s
    )
    scale = sector.engine_fraction if sector is not None else 1.0
    return MissionPoint(
        name=point.name,
        fuel=fuel,
        fuel_mass_flow_kg_s=fuel_flow * scale,
        air_mass_flow_kg_s=point.combustor_air_mass_flow_kg_s * scale,
        air_temperature_k=point.t3_k,
        pressure_pa=stations.dome_pa,
        duration_s=point.duration_s or 0.0,
        thrust_fraction=point.thrust_fraction or 1.0,
        nozzle_wall_temperature_k=point.nozzle_wall_temperature_k,
        operating_point=point,
        pressure_stations=stations,
    )


@dataclass(frozen=True)
class MissionProfile:
    """A validated collection of canonical mission points."""

    points: tuple[MissionPoint, ...]

    @classmethod
    def from_icao_lto(
        cls,
        points: tuple[OperatingPoint, ...] | list[OperatingPoint],
        sector: SectorDefinition | None = None,
    ) -> MissionProfile:
        """Require one Jet-A point for every ICAO LTO mode."""
        by_name = {point.name: point for point in points}
        missing = [name for name in ICAO_LTO_NAMES if name not in by_name]
        extra = [name for name in by_name if name not in ICAO_LTO_NAMES]
        duplicates = sorted(
            {point.name for point in points if sum(p.name == point.name for p in points) > 1}
        )
        if missing or extra or duplicates:
            raise ValueError(
                f"ICAO LTO profile mismatch; missing={missing or 'none'}, "
                f"unexpected={extra or 'none'}, duplicates={duplicates or 'none'}"
            )
        resolved = tuple(
            mission_point_from_operating(by_name[name], sector) for name in ICAO_LTO_NAMES
        )
        if any(point.fuel is not FuelKind.JET_A for point in resolved):
            raise ValueError("Every ICAO LTO point must burn Jet-A")
        if any(point.duration_s is None for point in points):
            raise ValueError("Every ICAO LTO point requires an explicit duration")
        return cls(resolved)

    @classmethod
    def from_cruise(
        cls,
        points: tuple[OperatingPoint, ...] | list[OperatingPoint],
        sector: SectorDefinition | None = None,
    ) -> MissionProfile:
        """Build named LNG cruise points from an engine-cycle deck."""
        if not points:
            raise ValueError("At least one cruise point is required")
        resolved = tuple(mission_point_from_operating(point, sector) for point in points)
        if any(point.fuel is not FuelKind.LNG for point in resolved):
            raise ValueError("Every cruise point must burn LNG")
        return cls(resolved)
