import pytest

from fuelnozzle.feed import LNGFeedLine, solve_lng_feed_line
from fuelnozzle.models import LNGComposition
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider


@pytest.fixture(scope="module")
def methane() -> CoolPropLNGProvider:
    return CoolPropLNGProvider(LNGComposition.pure_methane())


@pytest.fixture
def point() -> OperatingPoint:
    return OperatingPoint(
        name="cruise",
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


def test_measured_heat_leak_raises_enthalpy_and_ch_edl_computes_loss(
    point: OperatingPoint,
    methane: CoolPropLNGProvider,
) -> None:
    line = LNGFeedLine(
        length_m=5.0,
        inner_diameter_m=0.012,
        measured_heat_leak_w_per_m=20.0,
    )
    result = solve_lng_feed_line(point, line, methane)

    assert result.total_heat_leak_w == pytest.approx(100.0)
    assert result.outlet_state.enthalpy_j_kg > result.inlet_state.enthalpy_j_kg
    assert result.pressure_drop_pa > 0.0
    assert result.first_two_phase_position_m is None


def test_insulation_stack_uses_ht_and_returns_positive_heat_ingress(
    point: OperatingPoint,
    methane: CoolPropLNGProvider,
) -> None:
    line = LNGFeedLine(length_m=2.0, inner_diameter_m=0.012)
    result = solve_lng_feed_line(point, line, methane)

    assert result.total_heat_leak_w > 0.0
    assert all(path_point.heat_leak_w_per_m > 0.0 for path_point in result.path)
