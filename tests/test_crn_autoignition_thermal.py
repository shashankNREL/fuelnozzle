"""Verification of autoignition screening and cryogenic thermal management."""

from __future__ import annotations

from pathlib import Path

import cantera as ct
import numpy as np
import pytest

from fuelnozzle.crn.autoignition import (
    AutoignitionVerdict,
    IgnitionDelayTable,
    autoignition_margin,
    flashback_screen,
    ignition_delay,
    premix_state,
)
from fuelnozzle.crn.chemistry import (
    FuelKind,
    MechanismRegistry,
    MechanismRole,
    MechanismSpec,
)
from fuelnozzle.crn.thermal import (
    SupercriticalFeedError,
    heat_sink_budget,
    idle_circuit_screen,
    saturation_temperature_k,
    solve_temperature_for_duty,
    thermal_window,
)
from fuelnozzle.models import LNGComposition, WarningSeverity
from fuelnozzle.properties import CoolPropLNGProvider

MECH_DIR = Path(__file__).resolve().parents[1] / "mech"
PRESSURE_PA = 20.0 * ct.one_atm


def registry() -> MechanismRegistry:
    return MechanismRegistry(
        [
            MechanismSpec(
                path="gri30.yaml", fuel=FuelKind.LNG, role=MechanismRole.NETWORK,
                fuel_mole_fractions={"CH4": 1.0}, provenance="GRI-Mech 3.0",
            ),
            MechanismSpec(
                path=str(MECH_DIR / "A2NOx_skeletal.yaml"), fuel=FuelKind.JET_A,
                role=MechanismRole.NETWORK, fuel_mole_fractions={"POSF10325": 1.0},
                provenance="HyChem A2 + Glarborg",
            ),
            MechanismSpec(
                path=str(MECH_DIR / "A2NTCfast_ske.yaml"), fuel=FuelKind.JET_A,
                role=MechanismRole.IGNITION_DELAY,
                fuel_mole_fractions={"POSF10325": 1.0},
                provenance="HyChem A2 fast-NTC",
            ),
        ]
    )


def lng_provider() -> CoolPropLNGProvider:
    return CoolPropLNGProvider(LNGComposition.pure_methane())


# --- Premix state and flash cooling ---------------------------------------------------


def test_cold_fuel_cools_the_premixed_mixture_below_the_air():
    reg = registry()
    spec = reg.spec_for(FuelKind.LNG, MechanismRole.NETWORK)
    state = premix_state(
        reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK), spec,
        air_mass_flow_kg_s=1.0, air_temperature_k=800.0,
        fuel_mass_flow_kg_s=0.035, fuel_temperature_k=150.0, pressure_pa=PRESSURE_PA,
    )
    assert state.is_cooled_below_air
    assert state.temperature_drop_k > 40.0


def test_lng_cools_the_mixture_more_than_jet_a():
    """The flash-cooling advantage, which is what buys LNG its premixing margin."""
    reg = registry()

    def drop(fuel, fuel_temperature_k):
        spec = reg.spec_for(fuel, MechanismRole.NETWORK)
        return premix_state(
            reg.new_solution(fuel, MechanismRole.NETWORK), spec,
            air_mass_flow_kg_s=1.0, air_temperature_k=800.0,
            fuel_mass_flow_kg_s=0.035, fuel_temperature_k=fuel_temperature_k,
            pressure_pa=PRESSURE_PA,
        ).temperature_drop_k

    assert drop(FuelKind.LNG, 150.0) > 1.5 * drop(FuelKind.JET_A, 470.0)


def test_unevaporated_liquid_adds_latent_cooling():
    reg = registry()
    spec = reg.spec_for(FuelKind.LNG, MechanismRole.NETWORK)

    def state(liquid_fraction):
        return premix_state(
            reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK), spec,
            air_mass_flow_kg_s=1.0, air_temperature_k=800.0,
            fuel_mass_flow_kg_s=0.035, fuel_temperature_k=150.0,
            pressure_pa=PRESSURE_PA, liquid_fraction=liquid_fraction,
            latent_heat_j_kg=5.1e5,
        )

    dry, wet = state(0.0), state(0.3)
    assert wet.temperature_k < dry.temperature_k
    assert wet.latent_cooling_k > 0.0
    assert dry.latent_cooling_k == pytest.approx(0.0, abs=1.0e-6)


def test_latent_heat_is_the_smaller_part_of_the_cooling():
    """Most of the benefit is cold vapour, not the phase change; worth not overclaiming."""
    reg = registry()
    spec = reg.spec_for(FuelKind.LNG, MechanismRole.NETWORK)
    state = premix_state(
        reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK), spec,
        air_mass_flow_kg_s=1.0, air_temperature_k=800.0, fuel_mass_flow_kg_s=0.035,
        fuel_temperature_k=150.0, pressure_pa=PRESSURE_PA, liquid_fraction=0.3,
        latent_heat_j_kg=5.1e5,
    )
    assert state.latent_cooling_k < 0.25 * state.temperature_drop_k


def test_premix_state_rejects_impossible_inputs():
    reg = registry()
    spec = reg.spec_for(FuelKind.LNG, MechanismRole.NETWORK)
    solution = reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK)
    with pytest.raises(ValueError, match="positive air and fuel"):
        premix_state(
            solution, spec, air_mass_flow_kg_s=0.0, air_temperature_k=800.0,
            fuel_mass_flow_kg_s=0.035, fuel_temperature_k=150.0, pressure_pa=PRESSURE_PA,
        )
    with pytest.raises(ValueError, match="between 0 and 1"):
        premix_state(
            solution, spec, air_mass_flow_kg_s=1.0, air_temperature_k=800.0,
            fuel_mass_flow_kg_s=0.035, fuel_temperature_k=150.0,
            pressure_pa=PRESSURE_PA, liquid_fraction=1.5,
        )


# --- Ignition delay -------------------------------------------------------------------


def test_ignition_delay_shortens_with_temperature():
    reg = registry()
    spec = reg.spec_for(FuelKind.JET_A, MechanismRole.IGNITION_DELAY)
    solution = reg.new_solution(FuelKind.JET_A, MechanismRole.IGNITION_DELAY)

    cool = ignition_delay(solution, spec, 750.0, PRESSURE_PA, 0.5)
    hot = ignition_delay(solution, spec, 950.0, PRESSURE_PA, 0.5)
    assert cool is not None and hot is not None
    assert hot < cool


def test_methane_does_not_ignite_at_premixer_temperatures():
    """Why LNG tolerates a long premixer: at these conditions it is effectively inert."""
    reg = registry()
    spec = reg.spec_for(FuelKind.LNG, MechanismRole.NETWORK)
    solution = reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK)
    assert ignition_delay(solution, spec, 750.0, PRESSURE_PA, 0.6, max_time_s=1.0) is None


def test_table_interpolates_between_tabulated_temperatures():
    reg = registry()
    table = IgnitionDelayTable(
        reg, FuelKind.JET_A, (700.0, 800.0, 900.0), (PRESSURE_PA,), (0.5,)
    )
    middle = table(850.0, PRESSURE_PA, 0.5)
    at_800 = table(800.0, PRESSURE_PA, 0.5)
    at_900 = table(900.0, PRESSURE_PA, 0.5)
    assert middle is not None
    assert at_900 < middle < at_800


def test_table_returns_none_outside_its_range():
    reg = registry()
    table = IgnitionDelayTable(
        reg, FuelKind.JET_A, (700.0, 800.0), (PRESSURE_PA,), (0.5,)
    )
    assert table(400.0, PRESSURE_PA, 0.5) is None


def test_table_requires_at_least_two_temperatures():
    with pytest.raises(ValueError, match="At least two temperatures"):
        IgnitionDelayTable(registry(), FuelKind.JET_A, (800.0,), (PRESSURE_PA,), (0.5,))


# --- Autoignition margin ---------------------------------------------------------------


def jet_a_margin(residence_time_s: float, minimum_margin: float = 4.0):
    reg = registry()
    spec = reg.spec_for(FuelKind.JET_A, MechanismRole.NETWORK)
    state = premix_state(
        reg.new_solution(FuelKind.JET_A, MechanismRole.NETWORK), spec,
        air_mass_flow_kg_s=1.0, air_temperature_k=800.0, fuel_mass_flow_kg_s=0.035,
        fuel_temperature_k=470.0, pressure_pa=PRESSURE_PA,
    )
    table = IgnitionDelayTable(
        reg, FuelKind.JET_A, (700.0, 750.0, 800.0, 850.0, 900.0), (PRESSURE_PA,), (0.5,)
    )
    return autoignition_margin(
        table, state, residence_time_s, FuelKind.JET_A, minimum_margin=minimum_margin
    )


def test_short_passage_is_safe_and_long_passage_is_not():
    """Jet-A at landing and take-off tolerates only a short premixing passage."""
    assert jet_a_margin(2.0e-4).verdict is AutoignitionVerdict.SAFE
    assert jet_a_margin(2.0e-2).verdict is AutoignitionVerdict.UNSAFE


def test_autoignition_in_the_premixer_is_an_error_not_a_warning():
    result = jet_a_margin(2.0e-2)
    assert result.margin is not None and result.margin < 1.0
    severities = {warning.severity for warning in result.warnings}
    assert WarningSeverity.ERROR in severities


def test_margin_is_the_ratio_of_the_two_times():
    result = jet_a_margin(1.0e-3)
    assert result.ignition_delay_s is not None
    assert result.margin == pytest.approx(result.ignition_delay_s / 1.0e-3)


def test_marginal_verdict_between_one_and_the_threshold():
    result = jet_a_margin(1.0e-3, minimum_margin=1.0e6)
    assert result.verdict is AutoignitionVerdict.MARGINAL


def test_fallback_ignition_mechanism_is_flagged():
    """Using the high-temperature mechanism here would understate the danger."""
    reg = MechanismRegistry(
        [
            MechanismSpec(
                path=str(MECH_DIR / "A2NOx_skeletal.yaml"), fuel=FuelKind.JET_A,
                role=MechanismRole.NETWORK, fuel_mole_fractions={"POSF10325": 1.0},
                provenance="network mechanism only",
            )
        ]
    )
    table = IgnitionDelayTable(
        reg, FuelKind.JET_A, (800.0, 900.0), (PRESSURE_PA,), (0.5,)
    )
    assert not table.uses_dedicated_mechanism

    spec = reg.spec_for(FuelKind.JET_A, MechanismRole.NETWORK)
    state = premix_state(
        reg.new_solution(FuelKind.JET_A, MechanismRole.NETWORK), spec,
        air_mass_flow_kg_s=1.0, air_temperature_k=850.0, fuel_mass_flow_kg_s=0.035,
        fuel_temperature_k=470.0, pressure_pa=PRESSURE_PA,
    )
    result = autoignition_margin(table, state, 1.0e-3, FuelKind.JET_A)
    codes = {warning.code for warning in result.warnings}
    assert "IGNITION_MECHANISM_NOT_DEDICATED" in codes


def test_margin_rejects_a_nonpositive_residence_time():
    with pytest.raises(ValueError, match="Residence time must be positive"):
        jet_a_margin(0.0)


def test_unavailable_ignition_delay_is_unknown_and_fail_closed():
    reg = registry()
    spec = reg.spec_for(FuelKind.JET_A, MechanismRole.NETWORK)
    state = premix_state(
        reg.new_solution(FuelKind.JET_A, MechanismRole.NETWORK),
        spec,
        air_mass_flow_kg_s=1.0,
        air_temperature_k=1200.0,
        fuel_mass_flow_kg_s=0.035,
        fuel_temperature_k=470.0,
        pressure_pa=PRESSURE_PA,
    )
    table = IgnitionDelayTable(
        reg, FuelKind.JET_A, (700.0, 800.0), (PRESSURE_PA,), (0.5,)
    )
    result = autoignition_margin(table, state, 1.0e-3, FuelKind.JET_A)
    assert result.verdict is AutoignitionVerdict.UNKNOWN
    assert any(warning.severity is WarningSeverity.ERROR for warning in result.warnings)


# --- Flashback -------------------------------------------------------------------------


def test_flashback_screen_is_suppressed_without_a_flame_speed():
    """No invented flame speed underneath a safety screen."""
    screen = flashback_screen(50.0, laminar_flame_speed_m_s=None)
    assert screen.is_safe is None
    assert screen.margin is None
    codes = {warning.code for warning in screen.warnings}
    assert "FLASHBACK_SCREEN_UNAVAILABLE" in codes


def test_fast_passage_resists_flashback():
    screen = flashback_screen(80.0, laminar_flame_speed_m_s=0.4)
    assert screen.is_safe
    assert screen.margin is not None and screen.margin > 1.0


def test_slow_passage_is_flagged_for_flashback():
    screen = flashback_screen(1.0, laminar_flame_speed_m_s=0.9)
    assert screen.is_safe is False
    codes = {warning.code for warning in screen.warnings}
    assert "FLASHBACK_RISK" in codes


# --- Heat sink and thermal window -------------------------------------------------------


def test_heat_sink_duty_grows_with_delivery_temperature():
    provider = lng_provider()
    cold = heat_sink_budget(provider, 0.08, 112.0, 150.0, 2.0e6)
    warm = heat_sink_budget(provider, 0.08, 112.0, 200.0, 2.0e6)
    assert warm.total_duty_w > cold.total_duty_w


def test_subcooled_delivery_is_reported_as_not_flashing():
    provider = lng_provider()
    saturation = saturation_temperature_k(provider, 2.0e6)
    budget = heat_sink_budget(provider, 0.08, 112.0, saturation - 15.0, 2.0e6)
    codes = {warning.code for warning in budget.warnings}
    assert "LNG_SUBCOOLED_AT_NOZZLE" in codes
    assert budget.latent_w == 0.0


def test_superheated_delivery_splits_duty_across_all_three_mechanisms():
    provider = lng_provider()
    budget = heat_sink_budget(provider, 0.08, 112.0, 220.0, 2.0e6)
    assert budget.sensible_liquid_w > 0.0
    assert budget.latent_w > 0.0
    assert budget.superheat_w > 0.0
    assert budget.total_duty_w == pytest.approx(
        budget.sensible_liquid_w + budget.latent_w + budget.superheat_w
    )


def test_thermal_window_widens_with_feed_pressure():
    """Pump pressure is what buys thermal design freedom.

    Superheat at the injector and subcooling in the line compete for the temperature gap
    between saturation at feed pressure and at chamber pressure. Raising pump pressure
    widens that gap; changing fuel temperature cannot.
    """
    provider = lng_provider()
    candidates = np.arange(150.0, 260.0, 0.5).tolist()

    def window(feed_pressure_pa):
        return thermal_window(
            provider, candidates, chamber_pressure_pa=2.0e6,
            feed_pressure_pa=feed_pressure_pa, mass_flow_kg_s=0.08,
            tank_temperature_k=112.0,
        )

    narrow, wide = window(3.0e6), window(4.0e6)
    assert wide.saturation_gap_k > narrow.saturation_gap_k
    assert narrow.bounds_k is not None and wide.bounds_k is not None
    assert (wide.bounds_k[1] - wide.bounds_k[0]) > (narrow.bounds_k[1] - narrow.bounds_k[0])


def test_insufficient_pump_pressure_leaves_no_feasible_temperature():
    """Below a certain feed pressure the design is infeasible at any fuel temperature."""
    provider = lng_provider()
    window = thermal_window(
        provider, np.arange(150.0, 260.0, 0.5).tolist(), chamber_pressure_pa=2.0e6,
        feed_pressure_pa=2.5e6, mass_flow_kg_s=0.08, tank_temperature_k=112.0,
    )
    assert window.is_empty
    codes = {warning.code for warning in window.warnings}
    assert "THERMAL_WINDOW_EMPTY" in codes
    assert "SATURATION_GAP_TOO_NARROW" in codes


def test_supercritical_feed_pressure_is_explained_not_crashed():
    """Above the critical point there is no saturation, so the constraints do not apply."""
    provider = lng_provider()
    with pytest.raises(SupercriticalFeedError, match="critical point"):
        thermal_window(
            provider, [200.0], chamber_pressure_pa=2.0e6, feed_pressure_pa=6.0e6,
            mass_flow_kg_s=0.08, tank_temperature_k=112.0,
        )


def test_available_heat_can_bound_the_window_from_above():
    provider = lng_provider()
    candidates = np.arange(168.0, 200.0, 0.5).tolist()
    generous = thermal_window(
        provider, candidates, chamber_pressure_pa=2.0e6, feed_pressure_pa=4.0e6,
        mass_flow_kg_s=0.08, tank_temperature_k=112.0,
    )
    limited = thermal_window(
        provider, candidates, chamber_pressure_pa=2.0e6, feed_pressure_pa=4.0e6,
        mass_flow_kg_s=0.08, tank_temperature_k=112.0, available_heat_w=4.6e4,
    )
    assert len(limited.feasible_temperatures_k) < len(generous.feasible_temperatures_k)


def test_duty_inverts_back_to_a_temperature():
    provider = lng_provider()
    target = 200.0
    duty = heat_sink_budget(provider, 0.08, 112.0, target, 2.0e6).total_duty_w
    recovered = solve_temperature_for_duty(provider, 0.08, 112.0, 2.0e6, duty)
    assert recovered == pytest.approx(target, abs=1.0)


# --- Idle circuit -----------------------------------------------------------------------


def test_hot_idle_jet_a_circuit_is_a_coking_error():
    screen = idle_circuit_screen("lng", 480.0)
    assert screen.idle_fuel == "jet_a"
    assert screen.jet_a_coking_safe is False
    assert not screen.is_safe


def test_cool_idle_jet_a_circuit_passes():
    screen = idle_circuit_screen("lng", 400.0)
    assert screen.jet_a_coking_safe is True
    assert screen.is_safe


def test_warm_idle_lng_circuit_vapor_locks():
    screen = idle_circuit_screen("jet_a", 300.0, provider=lng_provider())
    assert screen.idle_fuel == "lng"
    assert screen.lng_vapor_locked is True
    codes = {warning.code for warning in screen.warnings}
    assert "IDLE_LNG_VAPOR_LOCK" in codes


def test_idle_lng_screen_needs_a_property_provider():
    screen = idle_circuit_screen("jet_a", 300.0)
    codes = {warning.code for warning in screen.warnings}
    assert "IDLE_LNG_SCREEN_UNAVAILABLE" in codes


def test_unknown_active_fuel_rejected():
    with pytest.raises(ValueError, match="Unknown active fuel"):
        idle_circuit_screen("hydrogen", 400.0)
