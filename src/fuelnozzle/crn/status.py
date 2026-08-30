"""Tri-state acceptance gates for physics and safety evidence."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class GateStatus(StrEnum):
    """Result of one required design-acceptance gate."""

    PASS = "pass"
    FAIL = "fail"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class GateResult:
    """One acceptance decision with its reason and evidence source."""

    name: str
    status: GateStatus
    reason: str
    evidence: str | None = None


def aggregate_gate_status(gates: tuple[GateResult, ...] | list[GateResult]) -> GateStatus:
    """Combine required gates without treating missing evidence as success."""
    statuses = {gate.status for gate in gates}
    if GateStatus.FAIL in statuses:
        return GateStatus.FAIL
    if GateStatus.UNKNOWN in statuses or not statuses:
        return GateStatus.UNKNOWN
    return GateStatus.PASS
