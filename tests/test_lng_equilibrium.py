import pytest

from fuelnozzle.lng import (
    FlashLocation,
    LNGNozzleGeometry,
    screen_lng_flash,
    solve_lng_equilibrium_flow,
)
from fuelnozzle.models import LNGComposition
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider


def make_point(*, name: str, p3_pa: float, lng_temperature_k: float) -> OperatingPoint:
    nozzle_drop = 9.0e5
    nozzle_inlet = p3_pa + nozzle_drop
    return OperatingPoint(
        name=name,
        p3_pa=p3_pa,
        t3_k=700.0,
        lng_mass_flow_kg_s=0.10,
        jet_a_mass_flow_kg_s=0.10,
        lng_pump_outlet_pressure_pa=nozzle_inlet + 2.0e5,
        lng_pump_outlet_temperature_k=lng_temperature_k,
        lng_nozzle_pressure_drop_pa=nozzle_drop,
        jet_a_pump_outlet_pressure_pa=nozzle_inlet + 2.0e5,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=nozzle_drop,
    )


@pytest.fixture(scope="module")
def methane() -> CoolPropLNGProvider:
    return CoolPropLNGProvider(LNGComposition.pure_methane())


def test_tier0_no_flash_when_p3_exceeds_bubble_pressure(
    methane: CoolPropLNGProvider,
) -> None:
    result = screen_lng_flash(
        make_point(name="cold", p3_pa=1.0e5, lng_temperature_k=110.0),
        methane,
    )

    assert result.flash_location == FlashLocation.NONE
    assert result.pressure_subcooling_margin_pa > 0.0
    assert result.equilibrium_flash_fraction_at_p3 == 0.0


def test_tier0_detects_internal_equilibrium_crossing(
    methane: CoolPropLNGProvider,
) -> None:
    result = screen_lng_flash(
        make_point(name="warm", p3_pa=1.0e5, lng_temperature_k=120.0),
        methane,
    )

    assert result.flash_location == FlashLocation.INTERNAL
    assert result.equilibrium_flash_onset_pressure_pa == pytest.approx(191430.0, rel=0.02)
    assert result.equilibrium_flash_fraction_at_p3 > 0.0


def test_tier1_returns_flow_bounds_and_sizes_parallel_orifices(
    methane: CoolPropLNGProvider,
) -> None:
    point = make_point(name="warm", p3_pa=1.0e5, lng_temperature_k=120.0)
    result = solve_lng_equilibrium_flow(
        point,
        LNGNozzleGeometry(number_of_orifices=4, orifice_length_m=1.0e-3),
        methane,
    )

    assert result.single_phase_mass_flux_kg_m2_s > 0.0
    assert result.hem_operating_mass_flux_kg_m2_s > 0.0
    assert result.required_geometric_area_m2 > 0.0
    assert result.required_orifice_diameter_m > 0.0
    assert result.critical_pressure_pa >= point.p3_pa
    assert len(result.path) == 180


def test_pressure_budget_deficit_is_reported(methane: CoolPropLNGProvider) -> None:
    point = make_point(name="bad-budget", p3_pa=1.0e5, lng_temperature_k=110.0)
    point = point.model_copy(update={"lng_pump_outlet_pressure_pa": 5.0e5})
    result = screen_lng_flash(point, methane)

    assert any(warning.code == "PRESSURE_BUDGET_DEFICIT" for warning in result.warnings)
