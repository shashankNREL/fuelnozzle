"""Design-space exploration: sensitivity, Pareto search, and the cost of sharing a liner.

The strategy is staged on purpose. A single global optimization over a partly-validated
reduced-order model produces a confident number that nobody should believe. Sweeping first
tells you which variables actually matter and how strongly, and that ranking is worth more
than any single optimum. Only then is a search over the reduced space meaningful.

The headline output is not a design but a *price*: how much worse a shared-liner
compromise is than a liner dedicated to either fuel alone. Because the two circuits are
separate hardware, that price is the whole cost of dual-fuel operation on the combustor
side, and nobody has it for this configuration.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.design import (
    DEFAULT_BOUNDS,
    DesignBound,
    DesignVector,
    baseline_design,
    from_unit_cube,
    perturb,
)
from fuelnozzle.crn.evaluate import DesignEvaluator, DesignResult
from fuelnozzle.crn.objectives import (
    ObjectiveName,
    ObjectiveVector,
    evaluate_objectives,
    pareto_front,
    rank_key,
)
from fuelnozzle.models import ModelWarning, WarningSeverity

DEFAULT_ORDER = (
    ObjectiveName.JET_A_LTO_DP_FOO,
    ObjectiveName.LNG_CRUISE_NOX,
)


@dataclass(frozen=True)
class Sample:
    """One evaluated design."""

    design: DesignVector
    result: DesignResult
    objectives: ObjectiveVector


@dataclass(frozen=True)
class SensitivityEntry:
    """How much one variable moves one objective."""

    variable: str
    objective: ObjectiveName
    elasticity: float
    low_value: float
    high_value: float

    @property
    def magnitude(self) -> float:
        return abs(self.elasticity)


def evaluate_design(evaluator: DesignEvaluator, design: DesignVector) -> Sample:
    result = evaluator.evaluate(design)
    return Sample(design=design, result=result, objectives=evaluate_objectives(result))


def one_at_a_time_sensitivity(
    evaluator: DesignEvaluator,
    base: DesignVector | None = None,
    bounds: tuple[DesignBound, ...] = DEFAULT_BOUNDS,
    step: float = 0.25,
    objectives: tuple[ObjectiveName, ...] = DEFAULT_ORDER,
) -> tuple[SensitivityEntry, ...]:
    """Move each variable up and down and record what each objective does.

    Reported as an elasticity, the fractional change in the objective divided by the
    fractional span of the variable, so variables with different units can be ranked
    against each other. One-at-a-time misses interactions by construction, which is why
    it is a screening step and not a conclusion.
    """
    base = base or baseline_design()
    reference = evaluate_design(evaluator, base)
    entries: list[SensitivityEntry] = []

    for bound in bounds:
        try:
            low = evaluate_design(evaluator, perturb(base, bound.name, -step, bounds))
            high = evaluate_design(evaluator, perturb(base, bound.name, +step, bounds))
        except (ValueError, KeyError):
            continue
        for objective in objectives:
            base_value = reference.objectives.values[objective]
            low_value = low.objectives.values[objective]
            high_value = high.objectives.values[objective]
            scale = abs(base_value) if abs(base_value) > 1.0e-12 else 1.0
            elasticity = (high_value - low_value) / (2.0 * step * scale)
            entries.append(
                SensitivityEntry(
                    variable=bound.name, objective=objective, elasticity=elasticity,
                    low_value=low_value, high_value=high_value,
                )
            )
    return tuple(entries)


def rank_variables(
    entries: tuple[SensitivityEntry, ...], objective: ObjectiveName
) -> tuple[tuple[str, float], ...]:
    """Variables ordered by how strongly they move one objective."""
    filtered = [entry for entry in entries if entry.objective is objective]
    filtered.sort(key=lambda entry: entry.magnitude, reverse=True)
    return tuple((entry.variable, entry.elasticity) for entry in filtered)


def sample_designs(
    evaluator: DesignEvaluator,
    count: int,
    bounds: tuple[DesignBound, ...] = DEFAULT_BOUNDS,
    base: DesignVector | None = None,
    seed: int = 20260829,
) -> tuple[Sample, ...]:
    """Evaluate a Latin-hypercube sample of the design space.

    Latin hypercube rather than uniform random: with the sample sizes a reactor-network
    model can afford, plain random sampling leaves visible gaps along individual
    variables, and stratifying each one removes that without costing anything.
    """
    if count < 1:
        raise ValueError("At least one sample is required")
    rng = np.random.default_rng(seed)
    dimensions = len(bounds)
    cut = np.arange(count) / count
    matrix = np.empty((count, dimensions))
    for column in range(dimensions):
        matrix[:, column] = rng.permutation(cut + rng.random(count) / count)

    samples: list[Sample] = []
    for row in matrix:
        design = from_unit_cube(
            {bound.name: value for bound, value in zip(bounds, row, strict=True)},
            base=base, bounds=bounds,
        )
        samples.append(evaluate_design(evaluator, design))
    return tuple(samples)


def best_design(
    samples: tuple[Sample, ...], order: tuple[ObjectiveName, ...] = DEFAULT_ORDER
) -> Sample:
    """The single best design under feasibility-first ranking."""
    if not samples:
        raise ValueError("No samples to choose from")
    return min(samples, key=lambda sample: rank_key(sample.objectives, order))


def best_for_objective(
    samples: tuple[Sample, ...], objective: ObjectiveName
) -> Sample:
    """The best feasible design for one objective, ignoring the others.

    This is what a liner dedicated to a single fuel would look like, and it is the
    reference the shared-liner compromise is measured against.
    """
    feasible = [sample for sample in samples if sample.objectives.is_feasible]
    pool = feasible or list(samples)
    return min(pool, key=lambda sample: sample.objectives.values[objective])


def pareto_samples(
    samples: tuple[Sample, ...], order: tuple[ObjectiveName, ...] = DEFAULT_ORDER
) -> tuple[Sample, ...]:
    indices = pareto_front([sample.objectives for sample in samples], order)
    return tuple(samples[index] for index in indices)


@dataclass(frozen=True)
class SharedLinerCost:
    """The price of making one liner serve both fuels.

    ``jet_a_penalty`` and ``lng_penalty`` are each computed **within one fuel and one
    mechanism**, comparing the compromise design against the best this same model found
    for that fuel alone. That is what makes them defensible even though absolute NOx is
    not validated and the two fuels use different mechanisms.
    """

    jet_a_only: Sample
    lng_only: Sample
    compromise: Sample
    jet_a_penalty: float
    lng_penalty: float
    is_degenerate: bool
    sample_count: int
    unbracketed: tuple[str, ...] = ()
    varied_coordinate: str | None = None

    @property
    def worst_penalty(self) -> float:
        return max(self.jet_a_penalty, self.lng_penalty)

    @property
    def warnings(self) -> tuple[ModelWarning, ...]:
        """Refuse to present a degenerate result as a finding.

        When the best compromise the sample contains is simply one of the two single-fuel
        optima, the sample has no genuine intermediate design in it. The penalties are
        then a statement about sampling density, not about the cost of sharing a liner,
        and reporting them as the latter would be wrong.
        """
        issues: list[ModelWarning] = []
        if self.is_degenerate:
            issues.append(
                ModelWarning(
                    code="SHARED_LINER_COST_DEGENERATE",
                    severity=WarningSeverity.ERROR,
                    message=(
                        f"The best compromise among {self.sample_count} designs is itself "
                        "one of the single-fuel optima, so the sample contains no genuine "
                        "compromise. These penalties measure sampling density, not the "
                        "cost of sharing a liner. Sample more densely, or reduce the "
                        "design space to the variables that matter, before reporting them."
                    ),
                )
            )
        if self.unbracketed:
            issues.append(
                ModelWarning(
                    code="SINGLE_FUEL_OPTIMUM_UNBRACKETED",
                    severity=WarningSeverity.ERROR,
                    message=(
                        f"The best design for {', '.join(self.unbracketed)} lies at the "
                        "edge of the swept range, so its optimum was never bracketed and "
                        "the true optimum is outside the sweep. The penalty measured "
                        "against it is a lower bound at best. Widen the range, or check "
                        "that the swept variable actually reaches the intended condition."
                    ),
                )
            )
        return tuple(issues)


def cost_of_shared_liner(
    samples: tuple[Sample, ...], order: tuple[ObjectiveName, ...] = DEFAULT_ORDER
) -> SharedLinerCost:
    """Compare a shared-liner compromise against per-fuel optima.

    The compromise is the feasible design minimizing the worse of the two fuels' relative
    penalties. Minimizing the sum instead would let a large loss on one fuel hide behind a
    small gain on the other, which is exactly the trade a dual-fuel design must not make
    blindly.
    """
    if not samples:
        raise ValueError("No samples to choose from")

    jet_a_only = best_for_objective(samples, ObjectiveName.JET_A_NOX)
    lng_only = best_for_objective(samples, ObjectiveName.LNG_NOX)
    jet_a_best = jet_a_only.objectives.values[ObjectiveName.JET_A_NOX]
    lng_best = lng_only.objectives.values[ObjectiveName.LNG_NOX]

    def penalties(sample: Sample) -> tuple[float, float]:
        jet_a = sample.objectives.values[ObjectiveName.JET_A_NOX]
        lng = sample.objectives.values[ObjectiveName.LNG_NOX]
        return (
            (jet_a - jet_a_best) / jet_a_best if jet_a_best > 0.0 else 0.0,
            (lng - lng_best) / lng_best if lng_best > 0.0 else 0.0,
        )

    feasible = [sample for sample in samples if sample.objectives.is_feasible]
    pool = feasible or list(samples)
    compromise = min(pool, key=lambda sample: max(penalties(sample)))
    jet_a_penalty, lng_penalty = penalties(compromise)

    degenerate = (
        compromise.design == jet_a_only.design or compromise.design == lng_only.design
    )

    # Identify a one-dimensional sweep from the actual varied design coordinate. Insertion
    # order is not physical ordering (Latin-hypercube samples are deliberately shuffled).
    dumped = [sample.design.model_dump() for sample in pool]
    varied = [
        name
        for name in dumped[0]
        if len({values[name] for values in dumped}) > 1
    ] if dumped else []
    varied_coordinate = varied[0] if len(varied) == 1 else None
    unbracketed: list[str] = []
    if len(pool) >= 3 and varied_coordinate is not None:
        ordered_pool = sorted(
            pool,
            key=lambda sample: getattr(sample.design, varied_coordinate),
        )
        for label, objective_name in (
            ("Jet-A", ObjectiveName.JET_A_NOX),
            ("LNG", ObjectiveName.LNG_NOX),
        ):
            values = [sample.objectives.values[objective_name] for sample in ordered_pool]
            best_index = values.index(min(values))
            if best_index in (0, len(values) - 1):
                unbracketed.append(label)

    return SharedLinerCost(
        jet_a_only=jet_a_only, lng_only=lng_only, compromise=compromise,
        jet_a_penalty=jet_a_penalty, lng_penalty=lng_penalty,
        is_degenerate=degenerate, sample_count=len(samples),
        unbracketed=tuple(unbracketed),
        varied_coordinate=varied_coordinate,
    )


def passage_sweep(
    evaluator: DesignEvaluator,
    base: DesignVector | None = None,
    shares: tuple[float, ...] = (0.15, 0.30, 0.50, 0.70, 0.85),
    mixing_fractions: tuple[float, ...] = (0.0, 0.5, 1.0),
) -> tuple[tuple[float, float, float, float], ...]:
    """Sweep the injector passage split and the idle-passage mixing assumption.

    Returns ``(share, mixing, jet_a_near_field_phi, lng_near_field_phi)``.

    The passage split is the one lever that lets fixed hardware present a different
    near-field mixture to each fuel. Its strength depends entirely on how much
    idle-passage air reaches the near field, which a reactor network cannot determine, so
    the mixing fraction is swept alongside rather than assumed.
    """
    base = base or baseline_design()
    rows: list[tuple[float, float, float, float]] = []
    for mixing in mixing_fractions:
        for share in shares:
            design = base.with_values(
                jet_a_passage_share=share, idle_passage_mixing_fraction=mixing
            )
            phis: dict[FuelKind, float] = {}
            for fuel in (FuelKind.JET_A, FuelKind.LNG):
                point = next(
                    (item for item in evaluator.mission if item.fuel is fuel), None
                )
                if point is None:
                    phis[fuel] = float("nan")
                    continue
                split = design.air_split(fuel)
                near_air = point.air_mass_flow_kg_s * split.near_field_air_fraction(fuel)
                phis[fuel] = (
                    point.fuel_mass_flow_kg_s
                    * evaluator.stoichiometric_afr(fuel)
                    / near_air
                    if near_air > 0.0
                    else float("inf")
                )
            rows.append((share, mixing, phis[FuelKind.JET_A], phis[FuelKind.LNG]))
    return tuple(rows)


def focused_sweep(
    evaluator: DesignEvaluator,
    variable: str,
    values: tuple[float, ...],
    base: DesignVector | None = None,
) -> tuple[Sample, ...]:
    """Sweep one variable finely, holding the rest fixed.

    This is what the sensitivity stage is for. Screening identifies the handful of
    variables that actually move the objectives, and the search is then run over those
    rather than over the full space, where the affordable sample density is far too
    sparse to resolve a compromise between two fuels.
    """
    base = base or baseline_design()
    samples: list[Sample] = []
    skipped: list[float] = []
    for value in values:
        try:
            design = base.with_values(**{variable: value})
        except ValueError:
            # The value makes the liner impossible, for instance by leaving no air for
            # dilution. Skipping is right; silently clipping would report a result for a
            # design the caller did not ask for.
            skipped.append(value)
            continue
        samples.append(evaluate_design(evaluator, design))
    if skipped:
        print(
            f"  note: {len(skipped)} value(s) of {variable} skipped as infeasible: "
            + ", ".join(f"{value:g}" for value in skipped)
        )
    return tuple(samples)
