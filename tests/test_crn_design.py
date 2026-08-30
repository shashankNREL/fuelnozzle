"""Verification of design variables, objectives, and the search machinery.

These tests exercise the parts that do not need a combustor solve, so they run fast and
can be dense. The expensive end-to-end behaviour is covered by a single marked test.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.design import (
    DEFAULT_BOUNDS,
    DesignBound,
    DesignVector,
    MissionPoint,
    VariableClass,
    baseline_design,
    bounds_by_class,
    from_unit_cube,
    perturb,
)
from fuelnozzle.crn.objectives import (
    ConstraintViolation,
    ObjectiveName,
    ObjectiveVector,
    dominates,
    pareto_front,
    rank_key,
)
from fuelnozzle.crn.optimize import DEFAULT_ORDER
from fuelnozzle.models import WarningSeverity

# --- Design variables ------------------------------------------------------------------


def test_only_the_shared_liner_forces_a_compromise():
    """Separate injector hardware and schedulable values do not."""
    assert VariableClass.SHARED_COMBUSTOR.is_compromised
    assert not VariableClass.PER_FUEL_INJECTOR.is_compromised
    assert not VariableClass.SCHEDULABLE.is_compromised


def test_dilution_takes_whatever_the_other_stations_leave():
    design = baseline_design()
    total = (
        design.dome_air_fraction
        + design.quench_air_fraction
        + design.primary_air_fraction
        + design.cooling_air_fraction
        + design.dilution_air_fraction
    )
    assert total == pytest.approx(1.0)


def test_air_fractions_summing_above_one_are_rejected():
    with pytest.raises(ValueError, match="leaving no room for dilution"):
        DesignVector(
            dome_air_fraction=0.70, quench_air_fraction=0.40,
            jet_a_passage_share=0.5,
        )


def test_with_values_revalidates_rather_than_copying():
    """pydantic's model_copy does not re-run validators.

    Copying without validation would let a sweep build a liner whose air fractions sum
    above one, and the failure would surface much later as a negative dilution fraction
    inside the network builder.
    """
    design = baseline_design()
    with pytest.raises(ValueError, match="leaving no room for dilution"):
        design.with_values(dome_air_fraction=0.70, quench_air_fraction=0.40)


def test_with_values_leaves_the_original_untouched():
    design = baseline_design()
    changed = design.with_values(dome_air_fraction=0.45)
    assert changed.dome_air_fraction == pytest.approx(0.45)
    assert design.dome_air_fraction == pytest.approx(0.38)


def test_air_split_carries_the_passage_lever_through():
    design = baseline_design().with_values(
        jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0
    )
    split = design.air_split(FuelKind.JET_A)
    assert split.jet_a_passage == pytest.approx(design.dome_air_fraction * 0.25)
    assert split.lng_passage == pytest.approx(design.dome_air_fraction * 0.75)


def test_per_fuel_values_differ_by_fuel():
    design = baseline_design().with_values(
        jet_a_premix_residence_s=1.0e-3, lng_premix_residence_s=4.0e-3,
        jet_a_fuel_temperature_k=440.0, lng_fuel_temperature_k=190.0,
    )
    assert design.premix_residence_s(FuelKind.JET_A) == pytest.approx(1.0e-3)
    assert design.premix_residence_s(FuelKind.LNG) == pytest.approx(4.0e-3)
    assert design.fuel_temperature_k(FuelKind.JET_A) == pytest.approx(440.0)
    assert design.fuel_temperature_k(FuelKind.LNG) == pytest.approx(190.0)


def test_injectors_that_do_not_fit_the_dome_are_an_error():
    design = baseline_design().with_values(dome_packaging_budget=1.4)
    codes = {warning.code for warning in design.packaging_warnings}
    assert "DOME_PACKAGING_EXCEEDED" in codes
    assert all(
        warning.severity is WarningSeverity.ERROR
        for warning in design.packaging_warnings
    )


def test_bounds_can_be_filtered_by_class():
    shared = bounds_by_class(VariableClass.SHARED_COMBUSTOR)
    injector = bounds_by_class(VariableClass.PER_FUEL_INJECTOR)
    assert {bound.name for bound in shared} >= {"dome_air_fraction", "quench_air_fraction"}
    assert {bound.name for bound in injector} >= {"jet_a_passage_share"}
    assert not set(shared) & set(injector)


def test_bound_maps_between_unit_space_and_real_space():
    bound = DesignBound("x", 0.2, 0.7, VariableClass.SHARED_COMBUSTOR)
    assert bound.denormalize(0.0) == pytest.approx(0.2)
    assert bound.denormalize(1.0) == pytest.approx(0.7)
    assert bound.normalize(bound.denormalize(0.4)) == pytest.approx(0.4)
    assert bound.clip(9.0) == pytest.approx(0.7)


def test_unit_cube_sample_repairs_an_impossible_air_split():
    """A sampler that drew a liner with no dilution would otherwise waste the budget."""
    design = from_unit_cube({"dome_air_fraction": 1.0, "quench_air_fraction": 1.0})
    assert design.dilution_air_fraction >= 0.0
    assert design.dome_air_fraction + design.quench_air_fraction < 1.0


def test_unit_cube_sample_covers_the_declared_range():
    low = from_unit_cube({name: 0.0 for name in (b.name for b in DEFAULT_BOUNDS)})
    high = from_unit_cube({"jet_a_passage_share": 1.0})
    assert low.jet_a_passage_share == pytest.approx(0.10)
    assert high.jet_a_passage_share == pytest.approx(0.90)


def test_perturbation_stays_inside_the_bounds():
    design = baseline_design()
    moved = perturb(design, "jet_a_passage_share", +10.0)
    assert moved.jet_a_passage_share <= 0.90


def test_perturbing_an_unswept_variable_raises():
    with pytest.raises(KeyError, match="not a swept design variable"):
        perturb(baseline_design(), "post_volume_m3", 0.1)


def test_mission_point_scales_fuel_flow():
    point = MissionPoint("takeoff", FuelKind.JET_A, 0.035, 1.0, 800.0, 2.0e6)
    assert point.scaled(2.0).fuel_mass_flow_kg_s == pytest.approx(0.070)
    assert point.scaled(2.0).air_mass_flow_kg_s == pytest.approx(1.0)


# --- Objectives and ranking -------------------------------------------------------------


def objective(jet_a: float, lng: float, violations=()) -> ObjectiveVector:
    return ObjectiveVector(
        values={
            ObjectiveName.JET_A_NOX: jet_a,
            ObjectiveName.LNG_NOX: lng,
            ObjectiveName.MIXING_NONUNIFORMITY: 0.0,
            ObjectiveName.EXIT_TEMPERATURE_SPREAD: 0.0,
        },
        violations=tuple(violations),
    )


def unsafe(amount: float = 0.5) -> ConstraintViolation:
    return ConstraintViolation("autoignition", amount, "ignites in the premixer")


def test_a_feasible_design_outranks_a_cleaner_infeasible_one():
    """The point of feasibility-first: no trading a safety constraint for emissions."""
    clean_but_unsafe = objective(1.0, 1.0, [unsafe()])
    dirty_but_safe = objective(500.0, 500.0)

    ranked = sorted(
        [clean_but_unsafe, dirty_but_safe],
        key=lambda item: rank_key(item, DEFAULT_ORDER),
    )
    assert ranked[0] is dirty_but_safe


def test_infeasible_designs_are_ordered_by_how_much_they_violate():
    mild = objective(100.0, 100.0, [unsafe(0.1)])
    severe = objective(1.0, 1.0, [unsafe(0.9)])
    ranked = sorted([severe, mild], key=lambda item: rank_key(item, DEFAULT_ORDER))
    assert ranked[0] is mild


def test_an_infeasible_design_never_dominates_a_feasible_one():
    assert not dominates(objective(1.0, 1.0, [unsafe()]), objective(9.0, 9.0), DEFAULT_ORDER)
    assert dominates(objective(9.0, 9.0), objective(1.0, 1.0, [unsafe()]), DEFAULT_ORDER)


def test_domination_requires_better_somewhere_and_worse_nowhere():
    assert dominates(objective(1.0, 1.0), objective(2.0, 2.0), DEFAULT_ORDER)
    assert not dominates(objective(1.0, 3.0), objective(2.0, 2.0), DEFAULT_ORDER)
    assert not dominates(objective(2.0, 2.0), objective(2.0, 2.0), DEFAULT_ORDER)


def test_pareto_front_keeps_the_tradeoffs_and_drops_the_dominated():
    objectives = [
        objective(1.0, 9.0),   # best on Jet-A
        objective(9.0, 1.0),   # best on LNG
        objective(4.0, 4.0),   # a genuine compromise
        objective(8.0, 8.0),   # dominated by the compromise
    ]
    front = set(pareto_front(objectives, DEFAULT_ORDER))
    assert front == {0, 1, 2}


def test_feasible_designs_dominate_the_whole_infeasible_set():
    objectives = [objective(5.0, 5.0), objective(0.1, 0.1, [unsafe()])]
    assert pareto_front(objectives, DEFAULT_ORDER) == (0,)


def test_total_violation_sums_the_constraint_amounts():
    combined = objective(1.0, 1.0, [unsafe(0.2), unsafe(0.3)])
    assert combined.total_violation == pytest.approx(0.5)
    assert not combined.is_feasible


# --- Cost of sharing a liner, and the guards that stop it being over-read ---------------


def fake_sample(jet_a: float, lng: float, dome: float, violations=()):
    """A Sample stand-in carrying only what the cost calculation reads."""
    from fuelnozzle.crn.optimize import Sample

    return Sample(
        design=baseline_design().with_values(dome_air_fraction=dome),
        result=None,
        objectives=objective(jet_a, lng, violations),
    )


def test_cost_of_shared_liner_measures_each_fuel_against_its_own_optimum():
    from fuelnozzle.crn.optimize import cost_of_shared_liner

    samples = (
        fake_sample(10.0, 40.0, 0.30),   # Jet-A optimum
        fake_sample(20.0, 20.0, 0.40),   # genuine compromise
        fake_sample(40.0, 10.0, 0.50),   # LNG optimum
    )
    cost = cost_of_shared_liner(samples)

    assert cost.jet_a_only.design.dome_air_fraction == pytest.approx(0.30)
    assert cost.lng_only.design.dome_air_fraction == pytest.approx(0.50)
    assert cost.compromise.design.dome_air_fraction == pytest.approx(0.40)
    assert cost.jet_a_penalty == pytest.approx(1.0)
    assert cost.lng_penalty == pytest.approx(1.0)
    assert not cost.is_degenerate


def test_a_compromise_that_is_just_one_endpoint_is_refused():
    """With no intermediate design in the sample, the penalties mean nothing."""
    from fuelnozzle.crn.optimize import cost_of_shared_liner

    samples = (fake_sample(10.0, 40.0, 0.30), fake_sample(40.0, 10.0, 0.50))
    cost = cost_of_shared_liner(samples)

    assert cost.is_degenerate
    codes = {warning.code for warning in cost.warnings}
    assert "SHARED_LINER_COST_DEGENERATE" in codes
    assert any(
        warning.severity is WarningSeverity.ERROR for warning in cost.warnings
    )


def test_an_optimum_at_the_edge_of_the_sweep_is_refused():
    """A monotonic curve means the optimum is outside the range, not at its end.

    This is the artifact that a degeneracy check alone does not catch: the compromise can
    be interior while a single-fuel optimum is still pinned to the boundary.
    """
    from fuelnozzle.crn.optimize import cost_of_shared_liner

    samples = (
        fake_sample(140.0, 71.0, 0.22),
        fake_sample(62.0, 9.0, 0.38),
        fake_sample(21.0, 2.0, 0.50),
        fake_sample(13.0, 2.3, 0.54),
        fake_sample(9.0, 3.3, 0.58),   # Jet-A still improving at the last point
    )
    cost = cost_of_shared_liner(samples)

    assert "Jet-A" in cost.unbracketed
    assert "LNG" not in cost.unbracketed
    codes = {warning.code for warning in cost.warnings}
    assert "SINGLE_FUEL_OPTIMUM_UNBRACKETED" in codes


def test_a_properly_bracketed_sweep_passes_both_guards():
    from fuelnozzle.crn.optimize import cost_of_shared_liner

    samples = (
        fake_sample(30.0, 60.0, 0.28),
        fake_sample(10.0, 40.0, 0.34),   # Jet-A optimum, interior
        fake_sample(20.0, 20.0, 0.42),   # compromise
        fake_sample(40.0, 10.0, 0.50),   # LNG optimum, interior
        fake_sample(60.0, 30.0, 0.58),
    )
    cost = cost_of_shared_liner(samples)

    assert cost.unbracketed == ()
    assert not cost.is_degenerate
    assert cost.warnings == ()


# --- The RQL optimum is conditional on quench air ---------------------------------------


def rql_nox_curve(quench_air: float, domes: tuple[float, ...]) -> list[float]:
    """Jet-A NOx against dome air at a fixed quench-air fraction."""
    import cantera as ct

    from fuelnozzle.crn.chemistry import (
        MechanismRegistry,
        MechanismRole,
        MechanismSpec,
    )
    from fuelnozzle.crn.design import DesignVector, MissionPoint
    from fuelnozzle.crn.evaluate import DesignEvaluator
    from fuelnozzle.crn.optimize import evaluate_design

    mech_dir = Path(__file__).resolve().parents[1] / "mech"
    registry = MechanismRegistry(
        [
            MechanismSpec(
                path=str(mech_dir / "A2NOx_skeletal.yaml"), fuel=FuelKind.JET_A,
                role=MechanismRole.NETWORK, fuel_mole_fractions={"POSF10325": 1.0},
                provenance="quench-air study",
            )
        ]
    )
    mission = (
        MissionPoint(
            "takeoff", FuelKind.JET_A, 0.035, 1.0, 800.0, 20.0 * ct.one_atm,
            duration_s=42.0,
        ),
    )
    evaluator = DesignEvaluator(registry, mission, architecture="rql")
    values = []
    for dome in domes:
        design = DesignVector(
            dome_air_fraction=dome, quench_air_fraction=quench_air,
            primary_air_fraction=0.03, cooling_air_fraction=0.03,
            jet_a_passage_share=0.5, idle_passage_mixing_fraction=1.0,
        )
        values.append(
            evaluate_design(evaluator, design).objectives.values[ObjectiveName.JET_A_NOX]
        )
    return values


def test_rql_optimum_exists_only_with_enough_quench_air():
    """Rich-quench-lean works by crossing stoichiometric quickly.

    That requires enough quench air as well as a short quench time. With too little, the
    mixture lingers near stoichiometric instead of crossing it, the classic rich-zone
    optimum disappears, and phi_rich near 1.4 becomes the *worst* case rather than the
    best. This is the condition under which the O-005 resolution holds.
    """
    domes = (0.28, 0.34, 0.38, 0.44)

    starved = rql_nox_curve(0.20, domes)
    generous = rql_nox_curve(0.38, domes)

    # Starved: no interior minimum -- the best point is at the rich end of the sweep.
    assert starved.index(min(starved)) == 0

    # Generous: a genuine interior minimum near the classical band.
    best = generous.index(min(generous))
    assert 0 < best < len(generous) - 1, "expected an interior optimum with ample quench air"
    assert 1.2 <= 0.035 * 14.691 / domes[best] <= 1.6

    # And more quench air is better outright at the optimum.
    assert min(generous) < min(starved)
