"""Verification of architecture templates and the mechanism merge checker."""

from __future__ import annotations

from pathlib import Path

import cantera as ct
import pytest

from fuelnozzle.crn.chemistry import (
    KNOWN_ALIAS_PAIRS,
    FuelKind,
    check_mechanism_merge,
)
from fuelnozzle.crn.network import CombustorNetwork
from fuelnozzle.crn.reactors import OutletSpec, ReactorKind
from fuelnozzle.crn.streams import AirSplit
from fuelnozzle.crn.templates import (
    RQL_RICH_BAND,
    ArchitectureInputs,
    QuenchSchedule,
    ldi_architecture,
    lpp_architecture,
    rql_architecture,
)
from fuelnozzle.models import WarningSeverity

JET_A_MECH = str(Path(__file__).resolve().parents[1] / "mech" / "A2NOx_skeletal.yaml")
JET_A_AFR = 14.69


def split(**overrides) -> AirSplit:
    values = dict(
        dome=0.20, primary=0.10, quench=0.35, dilution=0.25, cooling=0.10,
        jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0,
    )
    values.update(overrides)
    return AirSplit(**values)


def inputs(**overrides) -> ArchitectureInputs:
    values = dict(
        fuel=FuelKind.JET_A,
        fuel_mass_flow_kg_s=0.035,
        total_air_mass_flow_kg_s=1.0,
        air_temperature_k=800.0,
        air_split=split(),
        stoichiometric_air_fuel_ratio=JET_A_AFR,
    )
    values.update(overrides)
    return ArchitectureInputs(**values)


def as_network(architecture) -> CombustorNetwork:
    air_total = sum(inlet.mass_flow_kg_s for inlet in architecture.air_inlets)
    return CombustorNetwork(
        architecture.reactors,
        architecture.air_inlets,
        OutletSpec(source_reactor=architecture.outlet_reactor, mass_flow_kg_s=air_total),
        architecture.internal_flows,
    )


# --- Quench schedules -----------------------------------------------------------------


@pytest.mark.parametrize("schedule", list(QuenchSchedule))
def test_quench_weights_sum_to_one(schedule):
    assert sum(schedule.weights(5)) == pytest.approx(1.0)


def test_uniform_schedule_is_flat():
    assert QuenchSchedule.UNIFORM.weights(4) == pytest.approx((0.25, 0.25, 0.25, 0.25))


def test_front_loaded_schedule_puts_more_air_first():
    weights = QuenchSchedule.FRONT_LOADED.weights(4)
    assert weights[0] > weights[-1]
    assert list(weights) == sorted(weights, reverse=True)


def test_rear_loaded_schedule_puts_more_air_last():
    weights = QuenchSchedule.REAR_LOADED.weights(4)
    assert weights[-1] > weights[0]


def test_zero_quench_stages_rejected():
    with pytest.raises(ValueError, match="least one quench stage"):
        QuenchSchedule.UNIFORM.weights(0)


# --- Architectures --------------------------------------------------------------------


@pytest.mark.parametrize(
    "builder", [rql_architecture, ldi_architecture, lpp_architecture]
)
def test_architectures_are_mass_consistent_before_correction(builder):
    """A template that needs correcting has a bookkeeping error, not a physical one."""
    network = as_network(builder(inputs()))
    assert network.balance.correction_norm_kg_s < 1.0e-12


@pytest.mark.parametrize(
    "builder", [rql_architecture, ldi_architecture, lpp_architecture]
)
def test_all_supplied_air_enters_the_network(builder):
    architecture = builder(inputs())
    total = sum(inlet.mass_flow_kg_s for inlet in architecture.air_inlets)
    assert total == pytest.approx(1.0, rel=1.0e-9)


@pytest.mark.parametrize(
    "builder", [rql_architecture, ldi_architecture, lpp_architecture]
)
def test_spray_path_reactors_host_droplets(builder):
    architecture = builder(inputs())
    by_name = {spec.name: spec for spec in architecture.reactors}
    for name in architecture.spray_path:
        assert by_name[name].kind.hosts_droplets
        assert by_name[name].spray_path_length_m is not None


def test_rql_quench_is_staged_not_a_single_point():
    """A single mixing point would skip the stoichiometric crossing that makes the NOx."""
    architecture = rql_architecture(inputs(), quench_stages=5)
    quench_zones = [name for name in architecture.reactor_names if name.startswith("quench_")]
    assert len(quench_zones) == 5


def test_rql_quench_air_is_distributed_across_stages():
    architecture = rql_architecture(inputs(), quench_stages=4)
    quench_air = [
        inlet.mass_flow_kg_s
        for inlet in architecture.air_inlets
        if inlet.name.startswith("quench_air_")
    ]
    assert len(quench_air) == 4
    assert sum(quench_air) == pytest.approx(0.35, rel=1.0e-9)


def test_rql_front_loaded_schedule_changes_the_air_distribution():
    uniform = rql_architecture(inputs(), quench_stages=4, schedule=QuenchSchedule.UNIFORM)
    front = rql_architecture(
        inputs(), quench_stages=4, schedule=QuenchSchedule.FRONT_LOADED
    )

    def first_stage(architecture):
        return next(
            inlet.mass_flow_kg_s
            for inlet in architecture.air_inlets
            if inlet.name == "quench_air_1"
        )

    assert first_stage(front) > first_stage(uniform)


def test_rql_flags_a_rich_zone_outside_the_usual_band():
    """The default split gives an absurdly rich dome; the template must say so."""
    architecture = rql_architecture(inputs())
    assert architecture.near_field_equivalence_ratio > RQL_RICH_BAND[1]
    codes = {warning.code for warning in architecture.warnings}
    assert "RQL_RICH_ZONE_OUT_OF_BAND" in codes


def test_rql_accepts_a_realistic_rich_zone_without_complaint():
    """Sized for phi_rich near 1.5, the template should be satisfied."""
    dome_air = 0.035 * JET_A_AFR / 1.5
    architecture = rql_architecture(
        inputs(
            air_split=split(
                dome=dome_air, primary=0.05, quench=0.35,
                dilution=1.0 - dome_air - 0.05 - 0.35 - 0.05, cooling=0.05,
                jet_a_passage_share=1.0,
            )
        )
    )
    assert RQL_RICH_BAND[0] <= architecture.near_field_equivalence_ratio <= RQL_RICH_BAND[1]
    codes = {warning.code for warning in architecture.warnings}
    assert "RQL_RICH_ZONE_OUT_OF_BAND" not in codes


def test_rql_requires_at_least_one_quench_stage():
    with pytest.raises(ValueError, match="at least one quench stage"):
        rql_architecture(inputs(), quench_stages=0)


def test_ldi_has_no_premixing_passage():
    """LDI's immunity to autoignition comes from having nowhere for it to happen."""
    architecture = ldi_architecture(inputs())
    assert "premixer" not in architecture.reactor_names


def test_lpp_premixer_is_a_reacting_zone_and_says_so():
    architecture = lpp_architecture(inputs())
    assert "premixer" in architecture.reactor_names
    codes = {warning.code for warning in architecture.warnings}
    assert "LPP_PREMIXER_REACTS" in codes


def test_every_architecture_ends_in_a_plug_flow_zone():
    for builder in (rql_architecture, ldi_architecture, lpp_architecture):
        architecture = builder(inputs())
        outlet = next(
            spec for spec in architecture.reactors if spec.name == architecture.outlet_reactor
        )
        assert outlet.kind is ReactorKind.PFR


# --- The idle-passage lever (plan Section 8.2.1) --------------------------------------


def test_passage_sizing_changes_the_near_field_equivalence_ratio():
    """Same liner, same total air: only the passage split differs."""
    rich = rql_architecture(inputs(air_split=split(jet_a_passage_share=0.25)))
    lean = rql_architecture(inputs(air_split=split(jet_a_passage_share=1.0)))
    assert rich.near_field_equivalence_ratio > lean.near_field_equivalence_ratio


def test_idle_passage_air_is_routed_and_the_assumption_recorded():
    architecture = rql_architecture(
        inputs(air_split=split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0))
    )
    names = {inlet.name for inlet in architecture.air_inlets}
    assert "idle_passage_bypass" in names
    codes = {warning.code for warning in architecture.warnings}
    assert "IDLE_PASSAGE_AIR_ROUTED" in codes


def test_full_idle_mixing_removes_the_lever():
    """At full mixing the passages stop being a design variable."""
    segregated = rql_architecture(
        inputs(air_split=split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0))
    )
    mixed = rql_architecture(
        inputs(air_split=split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=1.0))
    )
    assert mixed.near_field_equivalence_ratio < segregated.near_field_equivalence_ratio


# --- Mechanism merge checker ----------------------------------------------------------


def test_merging_a_mechanism_with_itself_reports_duplicates_only():
    first, second = ct.Solution(JET_A_MECH), ct.Solution(JET_A_MECH)
    conflicts = check_mechanism_merge(first, second)
    kinds = {conflict.kind for conflict in conflicts}
    assert kinds == {"duplicate_reactions"}


def test_jet_a_and_gri_cannot_be_merged_naively():
    """Real conflicts, which is why the tool checks rather than merges."""
    conflicts = check_mechanism_merge(ct.Solution(JET_A_MECH), ct.Solution("gri30.yaml"))
    kinds = {conflict.kind for conflict in conflicts}

    assert "duplicate_reactions" in kinds
    assert "species_alias" in kinds
    assert "thermo_inconsistency" in kinds
    assert any(
        conflict.severity is WarningSeverity.ERROR for conflict in conflicts
    ), "an alias collision is not survivable and must be an error"


def test_alias_pairs_cover_the_naming_collision_actually_present():
    assert ("CH2*", "CH2(S)") in KNOWN_ALIAS_PAIRS


def test_merge_check_does_not_mutate_its_inputs():
    """The registry hands out cached templates; corrupting one would be invisible."""
    first, second = ct.Solution(JET_A_MECH), ct.Solution("gri30.yaml")
    first.TPX = 700.0, 5.0e5, {"O2": 0.21, "N2": 0.79}
    before = (first.T, first.P)

    check_mechanism_merge(first, second)

    assert (first.T, first.P) == pytest.approx(before)
