import pytest

from fuelnozzle.lng import (
    FlashLocation,
    LNGNozzleGeometry,
    RelaxationFlowSettings,
    solve_lng_relaxation_flow,
)
from fuelnozzle.models import LNGComposition
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider


@pytest.fixture(scope="module")
def methane() -> CoolPropLNGProvider:
    return CoolPropLNGProvider(LNGComposition.pure_methane())


@pytest.fixture
def warm_point() -> OperatingPoint:
    return OperatingPoint(
        name="warm",
        p3_pa=1.0e5,
        t3_k=700.0,
        lng_mass_flow_kg_s=0.10,
        jet_a_mass_flow_kg_s=0.10,
        lng_pump_outlet_pressure_pa=1.2e6,
        lng_pump_outlet_temperature_k=120.0,
        lng_nozzle_pressure_drop_pa=9.0e5,
        jet_a_pump_outlet_pressure_pa=1.2e6,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=9.0e5,
    )


def test_fast_relaxation_produces_internal_vapor(
    warm_point: OperatingPoint,
    methane: CoolPropLNGProvider,
) -> None:
    result = solve_lng_relaxation_flow(
        warm_point,
        LNGNozzleGeometry(number_of_orifices=4, orifice_length_m=2.0e-3),
        methane,
        relaxation_settings=RelaxationFlowSettings(relaxation_time_s=1.0e-7),
    )

    assert result.actual_flash_location in {FlashLocation.INTERNAL, FlashLocation.EXIT}
    assert result.actual_exit_vapor_quality_mass > 0.0
    assert result.actual_flash_onset_pressure_pa is not None


def test_pressure_delay_moves_flash_outside_nozzle(
    warm_point: OperatingPoint,
    methane: CoolPropLNGProvider,
) -> None:
    result = solve_lng_relaxation_flow(
        warm_point,
        LNGNozzleGeometry(number_of_orifices=4, orifice_length_m=1.0e-3),
        methane,
        relaxation_settings=RelaxationFlowSettings(
            relaxation_time_s=1.0,
            nucleation_pressure_delay_pa=1.0e5,
        ),
    )

    assert result.actual_flash_location == FlashLocation.EXTERNAL
    assert result.actual_exit_vapor_quality_mass == 0.0
    assert result.tier1.hem_operating_mass_flux_kg_m2_s <= result.operating_mass_flux_kg_m2_s
    assert result.operating_mass_flux_kg_m2_s <= result.tier1.single_phase_mass_flux_kg_m2_s


def test_no_flash_case_remains_liquid(methane: CoolPropLNGProvider) -> None:
    point = OperatingPoint(
        name="cold",
        p3_pa=1.0e5,
        t3_k=700.0,
        lng_mass_flow_kg_s=0.10,
        jet_a_mass_flow_kg_s=0.10,
        lng_pump_outlet_pressure_pa=1.2e6,
        lng_pump_outlet_temperature_k=110.0,
        lng_nozzle_pressure_drop_pa=9.0e5,
        jet_a_pump_outlet_pressure_pa=1.2e6,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=9.0e5,
    )
    result = solve_lng_relaxation_flow(
        point,
        LNGNozzleGeometry(),
        methane,
    )

    assert result.actual_flash_location == FlashLocation.NONE
    assert result.actual_exit_vapor_quality_mass == 0.0
