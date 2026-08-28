import pytest

from fuelnozzle.lng import LNGNozzleGeometry, RelaxationFlowSettings
from fuelnozzle.models import LNGComposition
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider
from fuelnozzle.spray import (
    FlashSprayCalibration,
    FlashSprayRegime,
    solve_lng_flash_spray,
)


@pytest.fixture(scope="module")
def methane() -> CoolPropLNGProvider:
    return CoolPropLNGProvider(LNGComposition.pure_methane())


def make_point(temperature_k: float) -> OperatingPoint:
    return OperatingPoint(
        name="spray",
        p3_pa=1.0e5,
        t3_k=700.0,
        lng_mass_flow_kg_s=0.10,
        jet_a_mass_flow_kg_s=0.10,
        lng_pump_outlet_pressure_pa=1.2e6,
        lng_pump_outlet_temperature_k=temperature_k,
        lng_nozzle_pressure_drop_pa=9.0e5,
        jet_a_pump_outlet_pressure_pa=1.2e6,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=9.0e5,
    )


def test_uncalibrated_tier3_returns_regime_but_not_smd(
    methane: CoolPropLNGProvider,
) -> None:
    result = solve_lng_flash_spray(
        make_point(120.0),
        LNGNozzleGeometry(number_of_orifices=4),
        methane,
        relaxation_settings=RelaxationFlowSettings(
            relaxation_time_s=1.0,
            nucleation_pressure_delay_pa=1.0e5,
        ),
    )

    assert result.regime == FlashSprayRegime.EXTERNAL_FLASH
    assert result.smd_estimate_m is None
    assert result.cfd_boundary.mass_flow_kg_s == pytest.approx(0.10)
    assert any(w.code == "FLASH_SPRAY_CALIBRATION_REQUIRED" for w in result.warnings)


def test_calibration_produces_bounded_smd_and_cone_angle(
    methane: CoolPropLNGProvider,
) -> None:
    calibration = FlashSprayCalibration(
        calibration_id="example-ln2-screening-only",
        reference_smd_m=30.0e-6,
        reference_pressure_ratio=2.0,
        reference_flash_quality=0.10,
        reference_full_cone_angle_deg=80.0,
    )
    result = solve_lng_flash_spray(
        make_point(140.0),
        LNGNozzleGeometry(number_of_orifices=4, orifice_length_m=2.0e-3),
        methane,
        relaxation_settings=RelaxationFlowSettings(relaxation_time_s=1.0e-7),
        calibration=calibration,
    )

    assert result.regime == FlashSprayRegime.FULLY_FLASHING
    assert result.smd_estimate_m is not None
    assert result.smd_range_m is not None
    assert result.smd_range_m[0] < result.smd_estimate_m < result.smd_range_m[1]
    assert result.full_cone_angle_range_deg is not None
    assert result.calibration_id == calibration.calibration_id


def test_cold_lng_is_mechanical_breakup(methane: CoolPropLNGProvider) -> None:
    result = solve_lng_flash_spray(
        make_point(110.0),
        LNGNozzleGeometry(),
        methane,
    )
    assert result.regime == FlashSprayRegime.MECHANICAL
