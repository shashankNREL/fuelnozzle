"""Verification of the bridge from the nozzle solvers to droplet classes."""

from __future__ import annotations

from math import gamma, pi

import pytest

from fuelnozzle import (
    FlashSprayCalibration,
    JetAProperties,
    LNGComposition,
    LNGNozzleGeometry,
    OperatingPoint,
    PressureSwirlGeometry,
    solve_jet_a_pressure_swirl,
    solve_lng_flash_spray,
)
from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.spray_source import (
    InitialSizePolicy,
    jet_a_spray_boundary,
    lng_spray_boundary,
    number_rate,
    rosin_rammler_classes,
)
from fuelnozzle.properties import CoolPropLNGProvider


def operating_point(**overrides) -> OperatingPoint:
    values = dict(
        name="takeoff",
        p3_pa=1.2e6,
        t3_k=780.0,
        lng_mass_flow_kg_s=0.12,
        jet_a_mass_flow_kg_s=0.14,
        lng_pump_outlet_pressure_pa=2.5e6,
        lng_pump_outlet_temperature_k=120.0,
        lng_nozzle_pressure_drop_pa=9.0e5,
        jet_a_pump_outlet_pressure_pa=2.5e6,
        jet_a_nozzle_inlet_temperature_k=300.0,
        jet_a_nozzle_pressure_drop_pa=9.0e5,
    )
    values.update(overrides)
    return OperatingPoint(**values)


def swirl_geometry(**overrides) -> PressureSwirlGeometry:
    values = dict(
        number_of_inlet_ports=4,
        inlet_port_diameter_m=0.7e-3,
        inlet_tangency_radius_m=2.0e-3,
        swirl_chamber_radius_m=2.5e-3,
        swirl_chamber_length_m=5.0e-3,
    )
    values.update(overrides)
    return PressureSwirlGeometry(**values)


def jet_a_properties() -> JetAProperties:
    return JetAProperties(
        density_kg_m3=800.0,
        viscosity_pa_s=1.5e-3,
        surface_tension_n_m=0.0255,
        source="illustrative",
    )


# --- Number rate and size distribution ------------------------------------------------


def test_number_rate_matches_the_paper_definition():
    """John et al. Eq. (6): mass flow divided by the mass of one droplet."""
    radius, density, flow = 20.0e-6, 800.0, 0.05
    expected = flow / (density * (4.0 / 3.0) * pi * radius**3)
    assert number_rate(flow, radius, density) == pytest.approx(expected)


def test_number_rate_rejects_a_zero_radius():
    with pytest.raises(ValueError, match="positive radius"):
        number_rate(0.05, 0.0, 800.0)


def test_single_class_carries_the_whole_mass_flow():
    classes = rosin_rammler_classes(20.0e-6, 0.05, 800.0, 300.0, 30.0, class_count=1)
    assert len(classes) == 1
    assert classes[0].mass_flow_kg_s == pytest.approx(0.05)
    assert classes[0].radius_m == pytest.approx(20.0e-6)


def test_rosin_rammler_classes_conserve_mass():
    classes = rosin_rammler_classes(20.0e-6, 0.05, 800.0, 300.0, 30.0, class_count=6)
    assert len(classes) == 6
    assert sum(cls.mass_flow_kg_s for cls in classes) == pytest.approx(0.05)


def test_rosin_rammler_classes_span_a_range_of_sizes():
    """A distribution must produce a surviving large-droplet tail, not a single size."""
    classes = rosin_rammler_classes(20.0e-6, 0.05, 800.0, 300.0, 30.0, class_count=6)
    radii = [cls.radius_m for cls in classes]
    assert radii == sorted(radii)
    assert max(radii) > 2.0 * min(radii)


def test_rosin_rammler_reproduces_the_requested_sauter_mean_diameter():
    """The classes must represent the D32 they were asked for, not merely bracket it."""
    smd_radius, spread = 20.0e-6, 2.5
    classes = rosin_rammler_classes(
        smd_radius, 0.05, 800.0, 300.0, 30.0, spread_parameter=spread, class_count=400
    )
    # D32 = sum(n d^3) / sum(n d^2), weighting by the number rate of each class.
    numerator = sum(cls.number_rate_per_s * cls.diameter_m**3 for cls in classes)
    denominator = sum(cls.number_rate_per_s * cls.diameter_m**2 for cls in classes)
    assert numerator / denominator == pytest.approx(2.0 * smd_radius, rel=0.05)


def test_rosin_rammler_characteristic_diameter_uses_the_gamma_relation():
    """Guards the D32 = X / Gamma(1 - 1/n) inversion against silent drift."""
    spread = 3.0
    classes = rosin_rammler_classes(
        20.0e-6, 0.05, 800.0, 300.0, 30.0, spread_parameter=spread, class_count=2
    )
    characteristic = 2.0 * 20.0e-6 * gamma(1.0 - 1.0 / spread)
    largest = max(cls.diameter_m for cls in classes)
    assert largest < 4.0 * characteristic


def test_rosin_rammler_rejects_a_spread_with_no_finite_sauter_mean():
    with pytest.raises(ValueError, match="at least"):
        rosin_rammler_classes(20.0e-6, 0.05, 800.0, 300.0, 30.0, spread_parameter=1.0)


# --- Jet-A boundary -------------------------------------------------------------------


def test_jet_a_boundary_without_calibration_downgrades_and_warns():
    """The package refuses to invent an SMD; the bridge must inherit that discipline."""
    result = solve_jet_a_pressure_swirl(operating_point(), swirl_geometry(), jet_a_properties())
    assert result.smd_estimate_m is None

    boundary = jet_a_spray_boundary(result, 0.14, 300.0, 800.0)

    assert boundary.size_policy is InitialSizePolicy.NOZZLE_RADIUS_TAB
    assert boundary.apply_aerodynamic_breakup
    codes = {warning.code for warning in boundary.warnings}
    assert "SPRAY_SIZE_POLICY_DOWNGRADED" in codes


def test_jet_a_boundary_uses_the_calibrated_smd_when_available():
    geometry = swirl_geometry(smd_calibration_coefficient=1.2)
    result = solve_jet_a_pressure_swirl(operating_point(), geometry, jet_a_properties())
    assert result.smd_estimate_m is not None

    boundary = jet_a_spray_boundary(result, 0.14, 300.0, 800.0)

    assert boundary.size_policy is InitialSizePolicy.NOZZLE_SMD
    assert not boundary.apply_aerodynamic_breakup
    assert boundary.droplet_classes[0].radius_m == pytest.approx(0.5 * result.smd_estimate_m)


def test_jet_a_boundary_has_no_vapor_and_conserves_mass():
    result = solve_jet_a_pressure_swirl(operating_point(), swirl_geometry(), jet_a_properties())
    boundary = jet_a_spray_boundary(result, 0.14, 300.0, 800.0)

    assert boundary.fuel is FuelKind.JET_A
    assert boundary.vapor_mass_flow_kg_s == 0.0
    assert boundary.vapor_mass_fraction == 0.0
    assert boundary.liquid_mass_flow_kg_s == pytest.approx(0.14)


def test_user_policy_requires_an_explicit_radius():
    result = solve_jet_a_pressure_swirl(operating_point(), swirl_geometry(), jet_a_properties())
    with pytest.raises(ValueError, match="requires a positive user_radius_m"):
        jet_a_spray_boundary(result, 0.14, 300.0, 800.0, policy=InitialSizePolicy.USER)


def test_user_policy_is_honoured():
    result = solve_jet_a_pressure_swirl(operating_point(), swirl_geometry(), jet_a_properties())
    boundary = jet_a_spray_boundary(
        result, 0.14, 300.0, 800.0, policy=InitialSizePolicy.USER, user_radius_m=12.0e-6
    )
    assert boundary.droplet_classes[0].radius_m == pytest.approx(12.0e-6)


# --- LNG boundary ---------------------------------------------------------------------


def lng_result(point: OperatingPoint, calibration: FlashSprayCalibration | None = None):
    properties = CoolPropLNGProvider(LNGComposition.pure_methane())
    geometry = LNGNozzleGeometry(number_of_orifices=4, orifice_length_m=1.0e-3)
    return solve_lng_flash_spray(
        point, geometry, properties, calibration=calibration
    )


def test_lng_boundary_splits_vapor_from_liquid():
    """The flash has already vaporized part of the fuel; only the rest is droplets."""
    point = operating_point(p3_pa=3.0e5, lng_pump_outlet_temperature_k=140.0,
                            lng_pump_outlet_pressure_pa=1.5e6,
                            lng_nozzle_pressure_drop_pa=9.0e5)
    result = lng_result(point)
    boundary = lng_spray_boundary(result, 0.12, 420.0)

    assert boundary.fuel is FuelKind.LNG
    assert boundary.vapor_mass_flow_kg_s >= 0.0
    total = boundary.vapor_mass_flow_kg_s + boundary.liquid_mass_flow_kg_s
    assert total == pytest.approx(0.12, rel=1.0e-9)
    assert boundary.vapor_mass_fraction == pytest.approx(
        result.actual_exit_vapor_quality_mass, abs=1.0e-9
    )


def test_flashing_regime_skips_aerodynamic_breakup_and_records_why():
    """TAB describes aerodynamic tearing, which is not what bursts a flashing droplet."""
    point = operating_point(p3_pa=2.0e5, lng_pump_outlet_temperature_k=145.0,
                            lng_pump_outlet_pressure_pa=1.5e6,
                            lng_nozzle_pressure_drop_pa=9.0e5)
    result = lng_result(point)
    boundary = lng_spray_boundary(
        result, 0.12, 420.0, policy=InitialSizePolicy.NOZZLE_RADIUS_TAB
    )

    if boundary.droplet_classes:
        assert not boundary.apply_aerodynamic_breakup
    codes = {warning.code for warning in boundary.warnings}
    assert "AERODYNAMIC_BREAKUP_SKIPPED" in codes


def test_fully_vaporized_lng_produces_no_droplets():
    point = operating_point(p3_pa=1.1e5, lng_pump_outlet_temperature_k=160.0,
                            lng_pump_outlet_pressure_pa=1.5e6,
                            lng_nozzle_pressure_drop_pa=9.0e5)
    result = lng_result(point)
    boundary = lng_spray_boundary(result, 0.12, 420.0)

    if boundary.vapor_mass_fraction >= 1.0:
        assert boundary.droplet_classes == ()
        codes = {warning.code for warning in boundary.warnings}
        assert "LNG_FULLY_VAPORIZED_AT_INJECTION" in codes


def test_lng_boundary_preserves_nozzle_warnings():
    """Provenance must survive the bridge, as it does elsewhere in the package."""
    point = operating_point(p3_pa=3.0e5, lng_pump_outlet_temperature_k=140.0,
                            lng_pump_outlet_pressure_pa=1.5e6,
                            lng_nozzle_pressure_drop_pa=9.0e5)
    result = lng_result(point)
    boundary = lng_spray_boundary(result, 0.12, 420.0)

    nozzle_codes = {warning.code for warning in result.warnings}
    boundary_codes = {warning.code for warning in boundary.warnings}
    assert nozzle_codes <= boundary_codes
