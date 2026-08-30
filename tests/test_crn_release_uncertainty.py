"""Fast checks for Stage 6 robust metrics and Stage 7 release evidence."""

from __future__ import annotations

import pytest

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.design import MissionPoint, baseline_design
from fuelnozzle.crn.evaluate import DesignResult, PointResult
from fuelnozzle.crn.evidence import ReleaseDisposition, build_conceptual_design_report
from fuelnozzle.crn.hardware import (
    AirAdmission,
    DualFuelHardware,
    FuelPassageGeometry,
    SectorDefinition,
    SharedLinerGeometry,
)
from fuelnozzle.crn.objectives import ObjectiveName, evaluate_objectives
from fuelnozzle.crn.status import GateResult, GateStatus
from fuelnozzle.crn.uncertainty import (
    UncertaintyCase,
    UncertaintyCategory,
    evaluate_uncertainty_ensemble,
)


def point(name: str, fuel: FuelKind, flow: float, duration: float, ei_nox: float) -> PointResult:
    mission = MissionPoint(name, fuel, flow, 1.0, 800.0, 2.0e6, duration_s=duration)
    return PointResult(
        point=mission,
        architecture="ldi",
        exit_temperature_k=1400.0,
        peak_temperature_k=1800.0,
        ei_nox_g_per_kg=ei_nox,
        equivalence_ratio_spread=float("nan"),
        near_field_equivalence_ratio=0.6,
        exit_temperature_spread_k=float("nan"),
        quench_residence_time_s=1.0e-3,
        autoignition=None,
        warnings=(),
        computational_status=GateStatus.PASS,
        acceptance_status=GateStatus.PASS,
        gates=(GateResult("test", GateStatus.PASS, "test evidence"),),
        mechanism_path="test.yaml",
        mechanism_provenance="test",
    )


def design_result(points=(), rated_thrust_kn=None) -> DesignResult:
    return DesignResult(
        design=baseline_design(),
        points=tuple(points),
        warnings=(),
        computational_status=GateStatus.PASS,
        acceptance_status=GateStatus.PASS,
        gates=(GateResult("test", GateStatus.PASS, "test evidence"),),
        rated_thrust_kn=rated_thrust_kn,
    )


def hardware() -> DualFuelHardware:
    def admission(name: str) -> AirAdmission:
        return AirAdmission(name, 1.0e-4, 0.7)

    return DualFuelHardware(
        liner=SharedLinerGeometry(
            1.0e-4,
            2.0e-4,
            3.0e-4,
            admission("dome"),
            admission("primary"),
            admission("quench"),
            admission("dilution"),
            admission("cooling"),
        ),
        jet_a_passage=FuelPassageGeometry("jet_a", 1.0e-4, 2.0e-4),
        lng_passage=FuelPassageGeometry("lng", 1.0e-4, 2.0e-4),
        sector=SectorDefinition(20),
    )


def test_lto_objective_is_fuel_flow_time_weighted_dp_foo():
    results = (
        point("takeoff", FuelKind.JET_A, 1.0, 42.0, 10.0),
        point("climb_out", FuelKind.JET_A, 0.8, 132.0, 20.0),
        point("approach", FuelKind.JET_A, 0.2, 240.0, 30.0),
        point("idle", FuelKind.JET_A, 0.05, 1560.0, 40.0),
    )
    result = design_result(results, rated_thrust_kn=100.0)
    expected = sum(
        item.ei_nox_g_per_kg
        * item.point.fuel_mass_flow_kg_s
        * item.point.duration_s
        for item in results
    ) / 100.0
    objectives = evaluate_objectives(result)
    assert objectives.values[ObjectiveName.JET_A_LTO_DP_FOO] == pytest.approx(expected)
    assert "jet_a:idle:ei_nox_g_per_kg" in objectives.named_metrics


class FakeEvaluator:
    def __init__(self, result: DesignResult):
        self.result = result

    def evaluate(self, _design):
        return self.result


def test_robust_feasibility_requires_every_uncertainty_category():
    result = design_result()
    one = evaluate_uncertainty_ensemble(
        baseline_design(),
        [UncertaintyCase("input", UncertaintyCategory.INPUT, FakeEvaluator(result))],
    )
    assert one.status is GateStatus.UNKNOWN

    complete = evaluate_uncertainty_ensemble(
        baseline_design(),
        [
            UncertaintyCase(category.value, category, FakeEvaluator(result))
            for category in UncertaintyCategory
        ],
    )
    assert complete.status is GateStatus.PASS


def test_default_release_report_is_no_go_without_external_evidence():
    report = build_conceptual_design_report(design_result(), hardware())
    assert report.disposition is ReleaseDisposition.NO_GO
    assert report.release_status in (GateStatus.FAIL, GateStatus.UNKNOWN)
    assert report.dimensions_m["liner_length"] is None
