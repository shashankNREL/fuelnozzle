"""Objectives and constraints for dual-fuel design search.

Constraints are handled by **feasibility-first ranking** rather than by adding penalties
to the objective. A feasible design always outranks an infeasible one, and infeasible
designs are ordered among themselves by how much they violate.

The distinction is not cosmetic. With a weighted penalty, a design that ignites inside its
own premixer can outrank a safe one by being sufficiently clean, and the optimizer will
find that trade because nothing stops it. Feasibility-first makes the trade unavailable.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from fuelnozzle.crn.autoignition import AutoignitionVerdict
from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.evaluate import DesignResult
from fuelnozzle.models import WarningSeverity


class ObjectiveName(StrEnum):
    """The quantities being minimized. All are phrased so that smaller is better."""

    JET_A_NOX = "jet_a_nox"
    LNG_NOX = "lng_nox"
    MIXING_NONUNIFORMITY = "mixing_nonuniformity"
    EXIT_TEMPERATURE_SPREAD = "exit_temperature_spread"


@dataclass(frozen=True)
class ConstraintViolation:
    """One reason a design is not usable."""

    name: str
    amount: float
    detail: str


@dataclass(frozen=True)
class ObjectiveVector:
    """Objectives and constraint state for one design."""

    values: dict[ObjectiveName, float]
    violations: tuple[ConstraintViolation, ...]

    @property
    def is_feasible(self) -> bool:
        return not self.violations

    @property
    def total_violation(self) -> float:
        return sum(violation.amount for violation in self.violations)

    def as_tuple(self, order: tuple[ObjectiveName, ...]) -> tuple[float, ...]:
        return tuple(self.values[name] for name in order)


def evaluate_objectives(result: DesignResult) -> ObjectiveVector:
    """Reduce a multi-point result to objectives and constraint violations."""
    violations: list[ConstraintViolation] = []

    for point in result.points:
        label = point.point.name
        if not point.feasible:
            errors = [
                warning.message
                for warning in point.warnings
                if warning.severity is WarningSeverity.ERROR
            ]
            violations.append(
                ConstraintViolation(
                    name=f"{label}:infeasible",
                    amount=1.0,
                    detail=errors[0] if errors else "point did not solve",
                )
            )
        if point.is_extinguished:
            violations.append(
                ConstraintViolation(
                    name=f"{label}:extinguished",
                    amount=1.0,
                    detail="the network converged to an unlit solution",
                )
            )
        margin = point.autoignition
        if margin is not None and margin.margin is not None:
            if margin.verdict is AutoignitionVerdict.UNSAFE:
                violations.append(
                    ConstraintViolation(
                        name=f"{label}:autoignition",
                        amount=max(0.0, 1.0 - margin.margin),
                        detail=(
                            f"the mixture ignites inside the premixer "
                            f"(margin {margin.margin:.2f})"
                        ),
                    )
                )
            elif margin.verdict is AutoignitionVerdict.MARGINAL:
                violations.append(
                    ConstraintViolation(
                        name=f"{label}:autoignition_margin",
                        amount=(margin.minimum_margin - margin.margin)
                        / margin.minimum_margin,
                        detail=(
                            f"autoignition margin {margin.margin:.2f} is below the "
                            f"required {margin.minimum_margin:.2f}"
                        ),
                    )
                )

    jet_a = result.weighted_ei_nox(FuelKind.JET_A)
    lng = result.weighted_ei_nox(FuelKind.LNG)
    spreads = [point.equivalence_ratio_spread for point in result.points]
    exit_spreads = [point.exit_temperature_spread_k for point in result.points]

    return ObjectiveVector(
        values={
            ObjectiveName.JET_A_NOX: jet_a,
            ObjectiveName.LNG_NOX: lng,
            # Uniform mixing is the goal, so the spread itself is minimized.
            ObjectiveName.MIXING_NONUNIFORMITY: (
                sum(spreads) / len(spreads) if spreads else 0.0
            ),
            ObjectiveName.EXIT_TEMPERATURE_SPREAD: (
                sum(exit_spreads) / len(exit_spreads) if exit_spreads else 0.0
            ),
        },
        violations=tuple(violations),
    )


def rank_key(objective: ObjectiveVector, order: tuple[ObjectiveName, ...]) -> tuple:
    """Sort key implementing feasibility-first ranking.

    Feasible designs sort before infeasible ones regardless of objective value; among
    infeasible designs, less violation sorts first.
    """
    return (
        0 if objective.is_feasible else 1,
        objective.total_violation,
        sum(objective.as_tuple(order)),
    )


def dominates(
    left: ObjectiveVector, right: ObjectiveVector, order: tuple[ObjectiveName, ...]
) -> bool:
    """Whether ``left`` is at least as good everywhere and better somewhere.

    Feasibility is checked first, so an infeasible design never dominates a feasible one
    however good its objectives look.
    """
    if left.is_feasible != right.is_feasible:
        return left.is_feasible
    if not left.is_feasible:
        return left.total_violation < right.total_violation
    a, b = left.as_tuple(order), right.as_tuple(order)
    return all(x <= y for x, y in zip(a, b, strict=True)) and any(
        x < y for x, y in zip(a, b, strict=True)
    )


def pareto_front(
    objectives: list[ObjectiveVector], order: tuple[ObjectiveName, ...]
) -> tuple[int, ...]:
    """Indices of the non-dominated designs.

    A plain quadratic sweep. The sample sizes a reduced-order combustor model can afford
    are in the hundreds, where the cost of anything cleverer is not repaid.
    """
    keep: list[int] = []
    for index, candidate in enumerate(objectives):
        if not any(
            dominates(other, candidate, order)
            for position, other in enumerate(objectives)
            if position != index
        ):
            keep.append(index)
    return tuple(keep)
