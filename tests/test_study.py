import pytest

from fuelnozzle.jet_a import JetAProperties, PressureSwirlGeometry
from fuelnozzle.lng import LNGNozzleGeometry
from fuelnozzle.models import LNGComposition
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.study import MissionFuelMasses, run_nozzle_study


def geometry() -> PressureSwirlGeometry:
    return PressureSwirlGeometry(
        number_of_inlet_ports=4,
        inlet_port_diameter_m=0.7e-3,
        inlet_tangency_radius_m=2.0e-3,
        swirl_chamber_radius_m=2.5e-3,
        swirl_chamber_length_m=5.0e-3,
    )


def properties() -> JetAProperties:
    return JetAProperties(
        density_kg_m3=800.0,
        viscosity_pa_s=1.5e-3,
        surface_tension_n_m=0.025,
        source="declared test properties",
    )


def point(name: str, temperature_k: float) -> OperatingPoint:
    return OperatingPoint(
        name=name,
        duration_s=100.0,
        flow_multiplier=2.0,
        p3_pa=1.0e5,
        t3_k=700.0,
        lng_mass_flow_kg_s=0.10,
        jet_a_mass_flow_kg_s=0.12,
        lng_pump_outlet_pressure_pa=1.2e6,
        lng_pump_outlet_temperature_k=temperature_k,
        lng_nozzle_pressure_drop_pa=9.0e5,
        jet_a_pump_outlet_pressure_pa=1.2e6,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=9.0e5,
    )


def test_study_runs_both_circuits_and_integrates_stage_mass() -> None:
    result = run_nozzle_study(
        [point("cold", 110.0), point("warm", 120.0)],
        LNGComposition.pure_methane(),
        LNGNozzleGeometry(number_of_orifices=4),
        geometry(),
        properties(),
        mission_fuel_masses=MissionFuelMasses(
            jet_a_kg=48.0,
            lng_kg=40.0,
            source="mission fixture",
        ),
    )

    assert result.coolprop_version.startswith("8.")
    assert len(result.operating_points) == 2
    assert all(item.jet_a is not None and item.lng is not None for item in result.operating_points)
    assert result.integrated_fuel_masses is not None
    assert result.integrated_fuel_masses.jet_a_kg == pytest.approx(48.0)
    assert result.integrated_fuel_masses.lng_kg == pytest.approx(40.0)
    assert not result.warnings


def test_mission_mismatch_is_a_warning_not_an_input_adjustment() -> None:
    result = run_nozzle_study(
        [point("cruise", 110.0)],
        LNGComposition.pure_methane(),
        LNGNozzleGeometry(),
        geometry(),
        properties(),
        mission_fuel_masses=MissionFuelMasses(
            jet_a_kg=100.0,
            lng_kg=100.0,
            source="external mission model",
        ),
    )
    assert len(result.warnings) == 2


def test_compositional_lng_runs_end_to_end() -> None:
    result = run_nozzle_study(
        [point("mixture", 120.0).model_copy(update={"duration_s": None})],
        LNGComposition(
            mole_fractions={"Methane": 0.90, "Ethane": 0.07, "Nitrogen": 0.03}
        ),
        LNGNozzleGeometry(number_of_orifices=4),
        geometry(),
        properties(),
    )

    lng_result = result.operating_points[0].lng
    assert lng_result is not None
    assert lng_result.tier2.tier1.required_orifice_diameter_m > 0.0
    assert result.integrated_fuel_masses is None
