"""Canonical mission, pressure-station, hardware, and sector tests."""

from __future__ import annotations

import pytest

from fuelnozzle.crn.hardware import (
    AirAdmission,
    AirflowMode,
    DualFuelHardware,
    FuelPassageGeometry,
    SectorDefinition,
    SharedLinerGeometry,
)
from fuelnozzle.crn.mission import ICAO_LTO_NAMES, MissionProfile, mission_point_from_operating
from fuelnozzle.operating import OperatingPoint, resolve_pressure_stations


def operating(name: str, *, jet_a: float = 0.03, lng: float = 0.0) -> OperatingPoint:
    return OperatingPoint(
        name=name,
        duration_s=42.0,
        p3_pa=2.0e6,
        t3_k=800.0,
        lng_mass_flow_kg_s=lng,
        jet_a_mass_flow_kg_s=jet_a,
        lng_pump_outlet_pressure_pa=3.0e6,
        lng_pump_outlet_temperature_k=120.0,
        lng_nozzle_pressure_drop_pa=0.5e6,
        jet_a_pump_outlet_pressure_pa=3.0e6,
        jet_a_nozzle_inlet_temperature_k=400.0,
        jet_a_nozzle_pressure_drop_pa=0.5e6,
        combustor_air_mass_flow_kg_s=1.2,
        liner_pressure_loss_fraction=0.05,
        thrust_fraction=1.0,
    )


def admission(name: str, area: float) -> AirAdmission:
    return AirAdmission(name, area, 0.7)


def liner() -> SharedLinerGeometry:
    return SharedLinerGeometry(
        quench_volume_m3=1.0e-3,
        flame_volume_m3=2.0e-3,
        post_volume_m3=3.0e-3,
        dome=admission("dome", 3.0e-4),
        primary=admission("primary", 1.0e-4),
        quench=admission("quench", 2.0e-4),
        dilution=admission("dilution", 3.0e-4),
        cooling=admission("cooling", 1.0e-4),
    )


def test_pressure_stations_close_liner_and_nozzle_budgets():
    stations = resolve_pressure_stations(operating("takeoff"))
    assert stations.liner_pressure_loss_pa == pytest.approx(0.1e6)
    assert stations.combustor_exit_pa == pytest.approx(1.9e6)
    assert stations.jet_a_nozzle_inlet_pa == pytest.approx(2.5e6)


def test_pressure_station_adapter_rejects_active_fuel_pump_deficit():
    point = operating("takeoff").model_copy(
        update={"jet_a_pump_outlet_pressure_pa": 2.4e6}
    )
    with pytest.raises(ValueError, match="pump pressure"):
        resolve_pressure_stations(point)


def test_sector_scaling_is_reversible():
    sector = SectorDefinition(cups_per_engine=20, modeled_cups=2)
    assert sector.from_engine_total(10.0) == pytest.approx(1.0)
    assert sector.to_engine_total(1.0) == pytest.approx(10.0)


def test_area_derived_split_follows_effective_area_for_common_pressure_ratio():
    split = liner().area_derived_split(
        2.0e6,
        1.9e6,
        800.0,
        jet_a_passage_share=0.4,
        idle_passage_mixing_fraction=0.25,
    )
    assert split.dome == pytest.approx(0.3)
    assert split.primary == pytest.approx(0.1)
    assert split.quench == pytest.approx(0.2)
    assert split.dilution == pytest.approx(0.3)
    assert split.cooling == pytest.approx(0.1)


def test_prescribed_airflow_requires_traceable_calibration():
    passage = FuelPassageGeometry("jet_a", 1.0e-4, 2.0e-4)
    with pytest.raises(ValueError, match="calibration identifier"):
        DualFuelHardware(
            liner=liner(),
            jet_a_passage=passage,
            lng_passage=FuelPassageGeometry("lng", 1.0e-4, 2.0e-4),
            sector=SectorDefinition(20),
            airflow_mode=AirflowMode.PRESCRIBED_CALIBRATION,
        )


def test_passage_share_comes_from_installed_effective_areas():
    hardware = DualFuelHardware(
        liner=liner(),
        jet_a_passage=FuelPassageGeometry("jet_a", 1.0e-4, 2.0e-4),
        lng_passage=FuelPassageGeometry("lng", 3.0e-4, 2.0e-4),
        sector=SectorDefinition(20),
    )
    assert hardware.jet_a_passage_share == pytest.approx(0.25)


def test_mission_adapter_applies_explicit_sector_scale():
    point = mission_point_from_operating(operating("takeoff"), SectorDefinition(20))
    assert point.fuel_mass_flow_kg_s == pytest.approx(0.03 / 20.0)
    assert point.air_mass_flow_kg_s == pytest.approx(1.2 / 20.0)
    assert point.operating_point is not None


def test_icao_adapter_requires_and_orders_all_four_modes():
    unordered = [operating(name) for name in reversed(ICAO_LTO_NAMES)]
    profile = MissionProfile.from_icao_lto(unordered)
    assert tuple(point.name for point in profile.points) == ICAO_LTO_NAMES
    with pytest.raises(ValueError, match="missing"):
        MissionProfile.from_icao_lto(unordered[:-1])
