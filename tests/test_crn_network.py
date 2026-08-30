"""Verification of reactor specifications, mass balancing, and the network solve."""

from __future__ import annotations

import cantera as ct
import numpy as np
import pytest

from fuelnozzle.crn.network import (
    CombustorNetwork,
    NetworkError,
    minimum_norm_mass_correction,
)
from fuelnozzle.crn.reactors import InletSpec, OutletSpec, ReactorKind, ReactorSpec

PRESSURE_PA = 10.0 * ct.one_atm


def gri30() -> ct.Solution:
    return ct.Solution("gri30.yaml")


# --- Reactor specifications -----------------------------------------------------------


def test_plug_flow_segment_volume_divides_the_zone():
    spec = ReactorSpec(
        name="post", kind=ReactorKind.PFR, volume_m3=4.0e-3, plug_flow_segments=8
    )
    assert spec.segment_volume_m3 == pytest.approx(0.5e-3)


def test_stirred_reactor_segment_volume_is_the_whole_zone():
    spec = ReactorSpec(name="flame", kind=ReactorKind.PSR, volume_m3=1.5e-3)
    assert spec.segment_volume_m3 == pytest.approx(1.5e-3)


def test_only_droplet_hosting_zones_accept_a_spray_path():
    ReactorSpec(
        name="evap", kind=ReactorKind.EVAPORATOR, volume_m3=1.0e-3,
        spray_path_length_m=0.02,
    )
    with pytest.raises(ValueError, match="Only evaporator and mixer"):
        ReactorSpec(
            name="flame", kind=ReactorKind.PSR, volume_m3=1.0e-3,
            spray_path_length_m=0.02,
        )


def test_kind_reports_whether_droplets_may_be_present():
    assert ReactorKind.EVAPORATOR.hosts_droplets
    assert ReactorKind.MIXER.hosts_droplets
    assert not ReactorKind.PSR.hosts_droplets
    assert not ReactorKind.PFR.hosts_droplets


# --- Minimum-norm mass correction ------------------------------------------------------


def test_correction_matches_the_hand_solved_minimum():
    """Three reactors in a loop, solved analytically.

    Constraints force z1 = z2 = t and z3 = t - 1, leaving one free parameter. Minimizing
    the squared change from (1.10, 0.95, 0.30) gives t = 3.35/3 = 1.116667.
    """
    flows = {("a", "b"): 1.10, ("b", "c"): 0.95, ("c", "a"): 0.30}
    report = minimum_norm_mass_correction(
        ["a", "b", "c"], flows, {"a": 1.0}, {"c": 1.0}
    )

    assert report.corrected_flows[("a", "b")] == pytest.approx(1.1166667, rel=1e-6)
    assert report.corrected_flows[("b", "c")] == pytest.approx(1.1166667, rel=1e-6)
    assert report.corrected_flows[("c", "a")] == pytest.approx(0.1166667, rel=1e-6)
    assert report.final_residual_kg_s < 1.0e-12


def test_correction_is_genuinely_minimum_norm():
    """Any other feasible solution must be a larger change than the one chosen."""
    flows = {("a", "b"): 1.10, ("b", "c"): 0.95, ("c", "a"): 0.30}
    report = minimum_norm_mass_correction(
        ["a", "b", "c"], flows, {"a": 1.0}, {"c": 1.0}
    )
    chosen = np.array([report.correction[edge] for edge in flows])

    for offset in (-0.2, -0.05, 0.05, 0.2):
        # Walk along the feasible direction (1, 1, 1), which preserves the balances.
        alternative = chosen + offset * np.array([1.0, 1.0, 1.0])
        assert np.linalg.norm(alternative) > np.linalg.norm(chosen)


def test_balanced_network_is_left_untouched():
    flows = {("a", "b"): 1.0, ("b", "c"): 1.0}
    report = minimum_norm_mass_correction(
        ["a", "b", "c"], flows, {"a": 1.0}, {"c": 1.0}
    )
    assert report.correction_norm_kg_s == pytest.approx(0.0, abs=1.0e-12)
    assert report.largest_relative_correction == pytest.approx(0.0, abs=1.0e-12)


def test_correction_closes_imbalances_of_the_magnitude_the_paper_reports():
    """John et al. Table 1 shows imbalances up to 9.3% of a reactor's inflow.

    Their connectivity matrix is not published, so their exact case cannot be rebuilt.
    What is checked is that imbalances of that size are closed, and that the correction
    stays small relative to the flows.
    """
    names = ["evaporator", "mixer", "flame1", "flame2", "post", "recirc1", "recirc2"]
    through = 0.1235
    # A consistent seven-reactor topology with two recirculation loops.
    balanced = {
        ("evaporator", "mixer"): through,
        ("mixer", "flame1"): through + 0.05,
        ("flame1", "flame2"): through + 0.09,
        ("flame2", "post"): through,
        ("flame1", "recirc1"): 0.05,
        ("recirc1", "mixer"): 0.05,
        ("flame2", "recirc2"): 0.09,
        ("recirc2", "flame1"): 0.09,
    }
    clean = minimum_norm_mass_correction(
        names, balanced, {"evaporator": through}, {"post": through}
    )
    assert clean.correction_norm_kg_s == pytest.approx(0.0, abs=1.0e-12)

    # Now perturb by up to 9.3%, the largest relative imbalance the paper reports.
    rng = np.random.default_rng(20260828)
    perturbations = rng.uniform(-0.093, 0.093, size=len(balanced))
    perturbed = {
        edge: flow * (1.0 + shift)
        for (edge, flow), shift in zip(balanced.items(), perturbations, strict=True)
    }
    report = minimum_norm_mass_correction(
        names, perturbed, {"evaporator": through}, {"post": through}
    )

    assert report.initial_residual_kg_s > 0.0
    assert report.final_residual_kg_s < 1.0e-12
    # The repair should be comparable to the damage, not far larger.
    assert report.largest_relative_correction < 0.2


def test_globally_inconsistent_boundaries_raise():
    """No redistribution of internal flows can fix a boundary that does not balance."""
    with pytest.raises(NetworkError, match="do not balance globally"):
        minimum_norm_mass_correction(
            ["a", "b"], {("a", "b"): 1.0}, {"a": 1.0}, {"b": 2.0}
        )


def test_large_correction_is_reported():
    flows = {("a", "b"): 1.0, ("b", "c"): 0.2}
    report = minimum_norm_mass_correction(
        ["a", "b", "c"], flows, {"a": 1.0}, {"c": 1.0}
    )
    codes = {warning.code for warning in report.warnings}
    assert "LARGE_MASS_CORRECTION" in codes


def test_correction_that_reverses_a_flow_is_reported():
    """A negative corrected flow means the topology, not the numbers, is wrong."""
    # Balance forces z_ab - z_ba = 1. Starting from flows summing to less than 1, the
    # minimum-norm solution is z_ba = (0.1 + 0.2 - 1)/2 = -0.35, i.e. reversed.
    flows = {("a", "b"): 0.1, ("b", "a"): 0.2}
    report = minimum_norm_mass_correction(["a", "b"], flows, {"a": 1.0}, {"b": 1.0})

    assert report.corrected_flows[("b", "a")] == pytest.approx(-0.35, rel=1.0e-6)
    codes = {warning.code for warning in report.warnings}
    assert "MASS_CORRECTION_REVERSED_FLOW" in codes


def test_unknown_reactor_in_a_flow_raises():
    with pytest.raises(NetworkError, match="unknown reactor"):
        minimum_norm_mass_correction(["a"], {("a", "ghost"): 1.0}, {"a": 1.0}, {"a": 1.0})


def test_duplicate_reactor_names_raise():
    with pytest.raises(NetworkError, match="unique"):
        minimum_norm_mass_correction(["a", "a"], {}, {"a": 1.0}, {"a": 1.0})


# --- Network assembly and solution -----------------------------------------------------


def build_network(
    fuel_flow_kg_s: float = 0.035,
    air_flow_kg_s: float = 1.0,
    air_temperature_k: float = 750.0,
    flame_volume_m3: float = 1.5e-3,
    heat_loss_w: float = 0.0,
) -> CombustorNetwork:
    total = air_flow_kg_s + fuel_flow_kg_s
    reactors = [
        ReactorSpec(
            name="flame", kind=ReactorKind.PSR, volume_m3=flame_volume_m3,
            heat_loss_w=heat_loss_w,
        ),
        ReactorSpec(name="recirc", kind=ReactorKind.PSR, volume_m3=1.0e-3),
        ReactorSpec(name="post", kind=ReactorKind.PSR, volume_m3=4.0e-3),
    ]
    inlets = [
        InletSpec(
            name="air", target_reactor="flame", mass_flow_kg_s=air_flow_kg_s,
            temperature_k=air_temperature_k, mole_fractions={"O2": 0.21, "N2": 0.79},
        ),
        InletSpec(
            name="fuel", target_reactor="flame", mass_flow_kg_s=fuel_flow_kg_s,
            temperature_k=300.0, mole_fractions={"CH4": 1.0},
        ),
    ]
    flows = {
        ("flame", "recirc"): 0.4,
        ("recirc", "flame"): 0.4,
        ("flame", "post"): total,
    }
    return CombustorNetwork(
        reactors, inlets, OutletSpec(source_reactor="post", mass_flow_kg_s=total), flows
    )


def test_network_solves_and_conserves_elements():
    """Atoms cannot be created; this checks the flows and the solver together.

    The bound reflects the time-marching convergence tolerance rather than machine
    precision: marching stops once temperatures move less than 0.05 K per chunk, which
    leaves a residual element imbalance of order 1e-6. That is physically negligible but
    is not zero, and pretending otherwise would hide a genuine loosening of the solve.
    """
    solution = build_network().solve(gri30, PRESSURE_PA)

    assert solution.converged
    assert solution.element_balance_error < 1.0e-5


def test_recirculation_does_not_prevent_a_steady_solution():
    """A reactor whose inlet depends on its own downstream still converges."""
    network = build_network()
    assert ("recirc", "flame") in network.flows
    solution = network.solve(gri30, PRESSURE_PA)
    assert solution.converged
    assert solution.by_name("recirc").temperature_k > 1500.0


def test_lean_methane_flame_reaches_a_plausible_temperature():
    """phi = 0.6 methane-air preheated to 750 K burns near 2000 K at 10 atm."""
    solution = build_network().solve(gri30, PRESSURE_PA)
    assert 1800.0 < solution.peak_temperature_k < 2200.0


def test_no_accumulates_along_the_flow_path():
    """NO formation is residence-time limited, so it must build downstream."""
    solution = build_network().solve(gri30, PRESSURE_PA)
    flame = solution.by_name("flame").mole_fractions["NO"]
    post = solution.by_name("post").mole_fractions["NO"]
    assert post > flame > 0.0


def test_residence_time_follows_density_volume_over_mass_flow():
    solution = build_network().solve(gri30, PRESSURE_PA)
    flame = solution.by_name("flame")
    expected = 1.5e-3 * flame.density_kg_m3 / flame.mass_flow_kg_s
    assert flame.residence_time_s == pytest.approx(expected, rel=1.0e-9)


def test_inflow_accounts_for_external_and_recirculated_streams():
    network = build_network()
    assert network.inflow_of("flame") == pytest.approx(1.0 + 0.035 + 0.4)
    assert network.inflow_of("recirc") == pytest.approx(0.4)


def test_hotter_inlet_air_raises_temperature_and_nox():
    """The Arrhenius sensitivity that makes NOx a temperature problem."""
    cool = build_network(air_temperature_k=700.0).solve(gri30, PRESSURE_PA)
    hot = build_network(air_temperature_k=850.0).solve(gri30, PRESSURE_PA)

    assert hot.peak_temperature_k > cool.peak_temperature_k
    assert hot.outlet.mole_fractions["NO"] > cool.outlet.mole_fractions["NO"]


def test_heat_loss_lowers_the_flame_temperature():
    adiabatic = build_network().solve(gri30, PRESSURE_PA)
    cooled = build_network(heat_loss_w=2.0e4).solve(gri30, PRESSURE_PA)
    assert cooled.by_name("flame").temperature_k < adiabatic.by_name("flame").temperature_k


def test_extinguished_network_is_reported_not_returned_silently():
    """An unlit solution solves the equations but is not a burning combustor."""
    solution = build_network(fuel_flow_kg_s=0.008).solve(gri30, PRESSURE_PA)

    assert solution.peak_temperature_k < 800.0, "case was expected not to sustain a flame"
    codes = {warning.code for warning in solution.warnings}
    assert "NETWORK_EXTINGUISHED" in codes


def test_inlet_targeting_an_unknown_reactor_raises():
    with pytest.raises(NetworkError, match="targets unknown reactor"):
        CombustorNetwork(
            [ReactorSpec(name="a", kind=ReactorKind.PSR, volume_m3=1.0e-3)],
            [
                InletSpec(
                    name="air", target_reactor="ghost", mass_flow_kg_s=1.0,
                    temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
                )
            ],
            OutletSpec(source_reactor="a", mass_flow_kg_s=1.0),
            {},
        )


def test_outlet_from_an_unknown_reactor_raises():
    with pytest.raises(NetworkError, match="Outlet leaves unknown reactor"):
        CombustorNetwork(
            [ReactorSpec(name="a", kind=ReactorKind.PSR, volume_m3=1.0e-3)],
            [
                InletSpec(
                    name="air", target_reactor="a", mass_flow_kg_s=1.0,
                    temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
                )
            ],
            OutletSpec(source_reactor="ghost", mass_flow_kg_s=1.0),
            {},
        )


def test_solution_lookup_by_name_rejects_unknown_reactors():
    solution = build_network().solve(gri30, PRESSURE_PA)
    with pytest.raises(KeyError):
        solution.by_name("ghost")
