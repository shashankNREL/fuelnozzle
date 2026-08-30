"""Evidence grades and fail-closed conceptual-design release reports."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from fuelnozzle.crn.evaluate import DesignResult
from fuelnozzle.crn.hardware import DualFuelHardware
from fuelnozzle.crn.status import GateResult, GateStatus, aggregate_gate_status
from fuelnozzle.crn.uncertainty import RobustDesignResult


class EvidenceGrade(StrEnum):
    VERIFIED = "verified"
    CALIBRATED = "calibrated"
    VALIDATED = "validated"
    EXTRAPOLATED = "extrapolated"
    UNAVAILABLE = "unavailable"
    MISSING = "missing"


class ReleaseDisposition(StrEnum):
    GO = "go"
    NO_GO = "no_go"


@dataclass(frozen=True)
class EvidenceRecord:
    """Evidence attached to one design-driving quantity."""

    quantity: str
    grade: EvidenceGrade
    domain: str
    uncertainty: str
    reference: str | None = None
    note: str = ""

    @property
    def gate_status(self) -> GateStatus:
        # Equation verification and calibration are necessary but do not demonstrate
        # predictive accuracy on independent data.
        return (
            GateStatus.PASS
            if self.grade is EvidenceGrade.VALIDATED
            else GateStatus.FAIL
            if self.grade is EvidenceGrade.MISSING
            else GateStatus.UNKNOWN
        )


DEFAULT_EVIDENCE: tuple[EvidenceRecord, ...] = (
    EvidenceRecord(
        "network_conservation",
        EvidenceGrade.VERIFIED,
        "declared CRN topology and solver tolerances",
        "numerical refinement still required per design",
        "tests/test_crn_network.py",
    ),
    EvidenceRecord(
        "engine_cycle_conditions",
        EvidenceGrade.UNAVAILABLE,
        "Jet-A four-mode LTO and named LNG cruise/off-design points",
        "unknown",
        note="A traceable engine cycle deck is not in the repository.",
    ),
    EvidenceRecord(
        "nozzle_hydraulics_and_spray",
        EvidenceGrade.UNAVAILABLE,
        "installed Jet-A and LNG injector operating ranges",
        "unknown",
        note="Hardware-matched flow and spray data are required.",
    ),
    EvidenceRecord(
        "lng_thermal_schedule",
        EvidenceGrade.UNAVAILABLE,
        "tank-to-injector path, heat source, soak, purge, and restart",
        "unknown",
    ),
    EvidenceRecord(
        "chemistry_and_emissions",
        EvidenceGrade.UNAVAILABLE,
        "LTO/cruise pressure, temperature, equivalence ratio, and residence time",
        "unknown",
        note="No held-out high-pressure emissions dataset is registered.",
    ),
    EvidenceRecord(
        "flashback_lbo_relight_switching_thermoacoustics",
        EvidenceGrade.MISSING,
        "installed dual-fuel combustor envelope",
        "unknown",
        note="These remain external gates until validated transient/rig models exist.",
    ),
    EvidenceRecord(
        "exit_pattern_factor_and_liner_thermal_state",
        EvidenceGrade.MISSING,
        "full annular exit and liner",
        "unknown",
    ),
)


@dataclass(frozen=True)
class ConceptualDesignReport:
    """Structured output that never substitutes absent evidence with a number."""

    dimensions_m: dict[str, float | None]
    effective_areas_m2: dict[str, float | None]
    mission_flows_kg_s: dict[str, dict[str, float]]
    pressure_loss_pa: dict[str, float | None]
    volumes_m3: dict[str, float]
    residence_times_s: dict[str, dict[str, float]]
    injector_ranges: dict[str, tuple[float, float] | None]
    thermal_schedule_k: dict[str, float | None]
    constraint_gates: tuple[GateResult, ...]
    pareto_position: int | None
    evidence: tuple[EvidenceRecord, ...]
    uncertainty_status: GateStatus
    independent_review_reference: str | None
    release_status: GateStatus
    disposition: ReleaseDisposition


def build_conceptual_design_report(
    result: DesignResult,
    hardware: DualFuelHardware,
    *,
    robust_result: RobustDesignResult | None = None,
    evidence: tuple[EvidenceRecord, ...] = DEFAULT_EVIDENCE,
    independent_review_reference: str | None = None,
    pareto_position: int | None = None,
) -> ConceptualDesignReport:
    """Build the Stage-7 schema and refuse release while any required gate is open."""
    mission_flows = {
        point.point.name: {
            "air": point.point.air_mass_flow_kg_s,
            "fuel": point.point.fuel_mass_flow_kg_s,
        }
        for point in result.points
    }
    pressure_losses = {
        point.point.name: (
            point.point.pressure_stations.liner_pressure_loss_pa
            if point.point.pressure_stations is not None
            else None
        )
        for point in result.points
    }
    residence_times = {
        point.point.name: {
            "quench": point.quench_residence_time_s,
            "premixer": result.design.premix_residence_s(point.point.fuel),
        }
        for point in result.points
    }
    thermal_schedule = {
        point.point.name: (
            result.design.fuel_temperature_k(point.point.fuel)
            if point.point.fuel.value == "lng"
            else None
        )
        for point in result.points
    }
    review_gate = GateResult(
        "independent_technical_review",
        GateStatus.PASS if independent_review_reference else GateStatus.UNKNOWN,
        (
            "independent review reference supplied"
            if independent_review_reference
            else "independent technical review has not been supplied"
        ),
        independent_review_reference,
    )
    evidence_gates = tuple(
        GateResult(
            f"evidence:{record.quantity}",
            record.gate_status,
            f"evidence grade is {record.grade.value}",
            record.reference,
        )
        for record in evidence
    )
    uncertainty_status = (
        robust_result.status if robust_result is not None else GateStatus.UNKNOWN
    )
    release_gates = (
        *result.gates,
        *evidence_gates,
        GateResult(
            "uncertainty_ensemble",
            uncertainty_status,
            "all required uncertainty categories must be covered and feasible",
        ),
        review_gate,
    )
    release_status = aggregate_gate_status(release_gates)
    # UNKNOWN is operationally NO-GO: release requires affirmative evidence.
    disposition = (
        ReleaseDisposition.GO
        if release_status is GateStatus.PASS
        else ReleaseDisposition.NO_GO
    )
    liner = hardware.liner
    return ConceptualDesignReport(
        dimensions_m={
            "liner_length": None,
            "liner_diameter": None,
            "jet_a_injector_diameter": None,
            "lng_injector_diameter": None,
        },
        effective_areas_m2={
            "dome": liner.dome.effective_area_m2,
            "primary": liner.primary.effective_area_m2,
            "quench": liner.quench.effective_area_m2,
            "dilution": liner.dilution.effective_area_m2,
            "cooling": liner.cooling.effective_area_m2,
            "jet_a_passage": hardware.jet_a_passage.air_effective_area_m2,
            "lng_passage": hardware.lng_passage.air_effective_area_m2,
        },
        mission_flows_kg_s=mission_flows,
        pressure_loss_pa=pressure_losses,
        volumes_m3={
            "quench": liner.quench_volume_m3,
            "flame": liner.flame_volume_m3,
            "post": liner.post_volume_m3,
        },
        residence_times_s=residence_times,
        injector_ranges={"jet_a": None, "lng": None},
        thermal_schedule_k=thermal_schedule,
        constraint_gates=release_gates,
        pareto_position=pareto_position,
        evidence=evidence,
        uncertainty_status=uncertainty_status,
        independent_review_reference=independent_review_reference,
        release_status=release_status,
        disposition=disposition,
    )
