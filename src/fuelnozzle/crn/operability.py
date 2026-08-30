"""Continuation-based CRN lean-limit screens, distinct from validated blowout claims."""

from __future__ import annotations

from dataclasses import dataclass

from fuelnozzle.crn.network import (
    CombustorNetwork,
    InitializationBranch,
    NetworkSolution,
    solve_continuation,
)
from fuelnozzle.models import ModelWarning, WarningSeverity


@dataclass(frozen=True)
class LeanLimitPoint:
    equivalence_ratio: float
    solution: NetworkSolution
    extinguished: bool
    co_mass_fraction: float | None


@dataclass(frozen=True)
class LeanLimitScreen:
    """Lit-branch continuation result; not a physical LBO prediction without calibration."""

    points: tuple[LeanLimitPoint, ...]
    last_lit_equivalence_ratio: float | None
    first_extinguished_equivalence_ratio: float | None
    extinction_bracket: tuple[float, float] | None
    calibration_id: str | None
    warnings: tuple[ModelWarning, ...]


def continuation_lean_limit_screen(
    networks: tuple[CombustorNetwork, ...] | list[CombustorNetwork],
    equivalence_ratios: tuple[float, ...] | list[float],
    solution_factory,
    pressure_pa: float,
    *,
    calibration_id: str | None = None,
    **solve_options,
) -> LeanLimitScreen:
    """Follow the hot branch from rich to lean and bracket numerical extinction."""
    phis = tuple(equivalence_ratios)
    if len(networks) != len(phis) or not networks:
        raise ValueError("One equivalence ratio is required for every non-empty network")
    if any(phi <= 0.0 for phi in phis):
        raise ValueError("Equivalence ratios must be positive")
    if any(right >= left for left, right in zip(phis, phis[1:], strict=False)):
        raise ValueError("Lean-limit continuation must be ordered from rich to lean")

    solutions = solve_continuation(
        networks,
        solution_factory,
        pressure_pa,
        initialization_branch=InitializationBranch.HOT,
        **solve_options,
    )
    points: list[LeanLimitPoint] = []
    for phi, solution in zip(phis, solutions, strict=True):
        extinguished = any(
            warning.code == "NETWORK_EXTINGUISHED" for warning in solution.warnings
        )
        co = solution.outlet.mass_fractions.get("CO")
        points.append(
            LeanLimitPoint(
                equivalence_ratio=phi,
                solution=solution,
                extinguished=extinguished,
                co_mass_fraction=co,
            )
        )

    first_extinguished_index = next(
        (index for index, point in enumerate(points) if point.extinguished),
        None,
    )
    last_lit = next(
        (
            point.equivalence_ratio
            for point in reversed(points)
            if not point.extinguished
        ),
        None,
    )
    first_extinguished = (
        points[first_extinguished_index].equivalence_ratio
        if first_extinguished_index is not None
        else None
    )
    bracket = (
        (points[first_extinguished_index - 1].equivalence_ratio, first_extinguished)
        if first_extinguished_index is not None and first_extinguished_index > 0
        else None
    )
    warnings = [
        ModelWarning(
            code="CRN_LEAN_LIMIT_NOT_LBO",
            severity=WarningSeverity.WARNING,
            message=(
                "The bracket marks loss of the continued CRN hot branch. It is not a "
                "physical lean-blowout limit and cannot establish operability without a "
                "pressure- and hardware-matched rig calibration."
            ),
        )
    ]
    if calibration_id is None:
        warnings.append(
            ModelWarning(
                code="LEAN_LIMIT_CALIBRATION_UNAVAILABLE",
                severity=WarningSeverity.WARNING,
                message="No high-pressure rig calibration identifier was supplied.",
            )
        )
    return LeanLimitScreen(
        points=tuple(points),
        last_lit_equivalence_ratio=last_lit,
        first_extinguished_equivalence_ratio=first_extinguished,
        extinction_bracket=bracket,
        calibration_id=calibration_id,
        warnings=tuple(warnings),
    )
