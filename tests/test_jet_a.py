import pytest

from fuelnozzle.jet_a import (
    JetAProperties,
    JetAPropertyTable,
    PressureSwirlGeometry,
    solve_jet_a_pressure_swirl,
)
from fuelnozzle.operating import OperatingPoint


@pytest.fixture
def point() -> OperatingPoint:
    return OperatingPoint(
        name="takeoff",
        p3_pa=1.0e6,
        t3_k=750.0,
        lng_mass_flow_kg_s=0.10,
        jet_a_mass_flow_kg_s=0.12,
        lng_pump_outlet_pressure_pa=2.5e6,
        lng_pump_outlet_temperature_k=120.0,
        lng_nozzle_pressure_drop_pa=1.0e6,
        jet_a_pump_outlet_pressure_pa=2.5e6,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=1.0e6,
    )


@pytest.fixture
def geometry() -> PressureSwirlGeometry:
    return PressureSwirlGeometry(
        number_of_inlet_ports=4,
        inlet_port_diameter_m=0.7e-3,
        inlet_tangency_radius_m=2.0e-3,
        swirl_chamber_radius_m=2.5e-3,
        swirl_chamber_length_m=5.0e-3,
        smd_calibration_coefficient=8.0,
    )


def test_pressure_swirl_sizes_exit_and_returns_hollow_cone(
    point: OperatingPoint,
    geometry: PressureSwirlGeometry,
) -> None:
    result = solve_jet_a_pressure_swirl(
        point,
        geometry,
        JetAProperties(
            density_kg_m3=800.0,
            viscosity_pa_s=1.5e-3,
            surface_tension_n_m=0.025,
            vapor_pressure_pa=1000.0,
            source="declared test properties",
        ),
    )

    assert result.required_exit_diameter_m > 0.0
    assert result.predicted_mass_flow_kg_s == pytest.approx(point.jet_a_mass_flow_kg_s)
    assert result.air_core_radius_m > 0.0
    assert result.liquid_film_thickness_m > 0.0
    assert 0.0 < result.full_cone_angle_deg < 180.0
    assert result.smd_range_m is not None
    assert result.cavitation_pressure_margin_pa is not None
    assert result.cavitation_pressure_margin_pa > 0.0


def test_property_table_interpolates_and_warns_outside_range() -> None:
    table = JetAPropertyTable(
        temperature_k=(280.0, 320.0),
        density_kg_m3=(820.0, 780.0),
        viscosity_pa_s=(2.0e-3, 1.0e-3),
        surface_tension_n_m=(0.027, 0.023),
        source="measured batch table",
    )
    properties, warnings = table.at_temperature(300.0)
    _, extrapolation_warnings = table.at_temperature(350.0)

    assert properties.density_kg_m3 == pytest.approx(800.0)
    assert not warnings
    assert any(w.code == "JET_A_PROPERTY_EXTRAPOLATION" for w in extrapolation_warnings)


def test_uncalibrated_geometry_suppresses_smd(
    point: OperatingPoint,
    geometry: PressureSwirlGeometry,
) -> None:
    result = solve_jet_a_pressure_swirl(
        point,
        geometry.model_copy(update={"smd_calibration_coefficient": None}),
        JetAProperties(
            density_kg_m3=800.0,
            viscosity_pa_s=1.5e-3,
            surface_tension_n_m=0.025,
            source="declared test properties",
        ),
    )
    assert result.smd_estimate_m is None
    assert any(w.code == "PRESSURE_SWIRL_SMD_CALIBRATION_REQUIRED" for w in result.warnings)
