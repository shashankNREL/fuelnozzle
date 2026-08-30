"""Empirical uncertainty ensembles for fail-closed design ranking."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from fuelnozzle.crn.design import DesignVector
from fuelnozzle.crn.evaluate import DesignEvaluator, DesignResult
from fuelnozzle.crn.objectives import ObjectiveName, ObjectiveVector, evaluate_objectives
from fuelnozzle.crn.status import GateStatus


class UncertaintyCategory(StrEnum):
    INPUT = "input"
    CALIBRATION = "calibration"
    MECHANISM = "mechanism"
    MANUFACTURING = "manufacturing"
    NUMERICAL = "numerical"
    MODEL_FORM = "model_form"


REQUIRED_UNCERTAINTY_CATEGORIES = frozenset(UncertaintyCategory)


@dataclass(frozen=True)
class UncertaintyCase:
    """One explicitly constructed perturbation or model substitution."""

    name: str
    category: UncertaintyCategory
    evaluator: DesignEvaluator


@dataclass(frozen=True)
class ObjectiveInterval:
    """Empirical central interval across declared uncertainty cases."""

    lower: float
    median: float
    upper: float


@dataclass(frozen=True)
class UncertaintyCaseResult:
    case: UncertaintyCase
    result: DesignResult
    objectives: ObjectiveVector


@dataclass(frozen=True)
class RobustDesignResult:
    """Nominal design evaluated under every supplied uncertainty case."""

    design: DesignVector
    cases: tuple[UncertaintyCaseResult, ...]
    intervals: dict[ObjectiveName, ObjectiveInterval]
    covered_categories: frozenset[UncertaintyCategory]
    status: GateStatus

    @property
    def robustly_feasible(self) -> bool:
        return self.status is GateStatus.PASS


def evaluate_uncertainty_ensemble(
    design: DesignVector,
    cases: tuple[UncertaintyCase, ...] | list[UncertaintyCase],
    *,
    confidence: float = 0.95,
    objectives: tuple[ObjectiveName, ...] = (
        ObjectiveName.JET_A_LTO_DP_FOO,
        ObjectiveName.LNG_CRUISE_NOX,
    ),
) -> RobustDesignResult:
    """Evaluate a design and rank it by worst-case feasibility plus empirical intervals."""
    if not 0.0 < confidence < 1.0:
        raise ValueError("Confidence must lie in (0, 1)")
    if not cases:
        raise ValueError("At least one uncertainty case is required")

    evaluated = tuple(
        UncertaintyCaseResult(
            case=case,
            result=(result := case.evaluator.evaluate(design)),
            objectives=evaluate_objectives(result),
        )
        for case in cases
    )
    tail = 0.5 * (1.0 - confidence)
    intervals: dict[ObjectiveName, ObjectiveInterval] = {}
    for objective in objectives:
        values = np.asarray(
            [case.objectives.values[objective] for case in evaluated],
            dtype=float,
        )
        finite = values[np.isfinite(values)]
        if finite.size:
            intervals[objective] = ObjectiveInterval(
                lower=float(np.quantile(finite, tail)),
                median=float(np.quantile(finite, 0.5)),
                upper=float(np.quantile(finite, 1.0 - tail)),
            )

    covered = frozenset(case.category for case in cases)
    if any(not case.objectives.is_feasible for case in evaluated):
        status = GateStatus.FAIL
    elif not REQUIRED_UNCERTAINTY_CATEGORIES <= covered:
        status = GateStatus.UNKNOWN
    else:
        status = GateStatus.PASS
    return RobustDesignResult(
        design=design,
        cases=evaluated,
        intervals=intervals,
        covered_categories=covered,
        status=status,
    )
