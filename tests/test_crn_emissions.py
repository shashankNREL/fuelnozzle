"""Verification of emission indices, corrected concentrations, and the LTO metric."""

from __future__ import annotations

import cantera as ct
import pytest

from fuelnozzle.crn.emissions import (
    ICAO_LTO_MODES,
    LTOMode,
    lto_dp_foo,
    standard_lto_modes,
    summarize_emissions,
)
from fuelnozzle.models import WarningSeverity


def burned_state() -> tuple[ct.Solution, dict[str, float]]:
    gas = ct.Solution("gri30.yaml")
    gas.TP = 800.0, 20.0 * ct.one_atm
    gas.set_equivalence_ratio(0.5, "CH4:1", "O2:0.21,N2:0.79")
    reactor = ct.IdealGasConstPressureReactor(gas, clone=False)
    ct.ReactorNet([reactor]).advance(0.05)
    return ct.Solution("gri30.yaml"), dict(
        zip(gas.species_names, gas.Y, strict=True)
    )


def test_summary_reports_nox_as_the_sum_of_no_and_no2():
    solution, mass_fractions = burned_state()
    summary = summarize_emissions(
        solution, 1900.0, 20.0 * ct.one_atm, mass_fractions, 1.035, 0.035
    )
    assert summary.nox_ppmv_dry_15pct_o2 == pytest.approx(
        summary.no_ppmv_dry_15pct_o2 + summary.no2_ppmv_dry_15pct_o2
    )
    assert summary.ei_nox_g_per_kg > 0.0


def test_co_is_reported_but_flagged_uncalibrated():
    solution, mass_fractions = burned_state()
    summary = summarize_emissions(
        solution, 1900.0, 20.0 * ct.one_atm, mass_fractions, 1.035, 0.035
    )
    assert summary.co_ppmv_dry_15pct_o2 is not None
    codes = {warning.code for warning in summary.warnings}
    assert "CO_UNCALIBRATED" in codes


def test_summary_requires_a_positive_fuel_flow():
    solution, mass_fractions = burned_state()
    with pytest.raises(ValueError, match="positive fuel mass flow"):
        summarize_emissions(
            solution, 1900.0, 20.0 * ct.one_atm, mass_fractions, 1.035, 0.0
        )


def test_lto_mode_mass_is_index_times_fuel_times_time():
    mode = LTOMode("takeoff", 1.0, 42.0, 0.5, 20.0)
    assert mode.fuel_mass_kg == pytest.approx(21.0)
    assert mode.nox_mass_g == pytest.approx(20.0 * 21.0)


def test_dp_foo_is_cycle_nox_over_rated_thrust():
    modes = standard_lto_modes(
        {"takeoff": 0.5, "climb_out": 0.4, "approach": 0.2, "idle": 0.05},
        {"takeoff": 30.0, "climb_out": 25.0, "approach": 12.0, "idle": 4.0},
    )
    result = lto_dp_foo(modes, rated_thrust_kn=100.0)
    expected = sum(mode.nox_mass_g for mode in modes) / 100.0
    assert result.dp_foo_g_per_kn == pytest.approx(expected)
    assert len(result.modes) == len(ICAO_LTO_MODES)


def test_idle_dominates_the_time_in_mode():
    """26 minutes at idle against 0.7 at takeoff; a partial cycle is badly misleading."""
    modes = standard_lto_modes(
        {"takeoff": 0.5, "climb_out": 0.4, "approach": 0.2, "idle": 0.05},
        {"takeoff": 30.0, "climb_out": 25.0, "approach": 12.0, "idle": 4.0},
    )
    by_name = {mode.name: mode for mode in modes}
    assert by_name["idle"].duration_s > 30.0 * by_name["takeoff"].duration_s


def test_partial_cycle_is_flagged():
    modes = standard_lto_modes({"takeoff": 0.5}, {"takeoff": 30.0})
    result = lto_dp_foo(modes, rated_thrust_kn=100.0)
    codes = {warning.code for warning in result.warnings}
    assert "LTO_CYCLE_INCOMPLETE" in codes


def test_dp_foo_always_warns_it_is_not_certification():
    modes = standard_lto_modes(
        {"takeoff": 0.5, "climb_out": 0.4, "approach": 0.2, "idle": 0.05},
        {"takeoff": 30.0, "climb_out": 25.0, "approach": 12.0, "idle": 4.0},
    )
    result = lto_dp_foo(modes, rated_thrust_kn=100.0)
    warning = next(
        w for w in result.warnings if w.code == "LTO_NOT_A_CERTIFICATION_RESULT"
    )
    assert warning.severity is WarningSeverity.WARNING


def test_lto_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="Rated thrust must be positive"):
        lto_dp_foo([LTOMode("takeoff", 1.0, 42.0, 0.5, 20.0)], rated_thrust_kn=0.0)
    with pytest.raises(ValueError, match="At least one LTO mode"):
        lto_dp_foo([], rated_thrust_kn=100.0)


def test_operating_point_reports_its_active_fuel():
    from fuelnozzle import OperatingPoint

    def point(lng, jet_a):
        return OperatingPoint(
            name="p", p3_pa=1.2e6, t3_k=780.0,
            lng_mass_flow_kg_s=lng, jet_a_mass_flow_kg_s=jet_a,
            lng_pump_outlet_pressure_pa=2.5e6, lng_pump_outlet_temperature_k=120.0,
            lng_nozzle_pressure_drop_pa=9.0e5, jet_a_pump_outlet_pressure_pa=2.5e6,
            jet_a_nozzle_inlet_temperature_k=300.0, jet_a_nozzle_pressure_drop_pa=9.0e5,
        )

    assert point(0.1, 0.0).active_fuel == "lng"
    assert point(0.0, 0.1).active_fuel == "jet_a"
    assert point(0.1, 0.1).active_fuel is None
    assert point(0.0, 0.0).active_fuel is None
