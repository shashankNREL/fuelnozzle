"""Flight-envelope orchestration for the Jet-A and LNG nozzle models."""

from __future__ import annotations

from dataclasses import dataclass

from pydantic import BaseModel, ConfigDict, Field

from fuelnozzle.feed import LNGFeedLine
from fuelnozzle.jet_a import (
    JetAProperties,
    JetAPropertyTable,
    PressureSwirlGeometry,
    PressureSwirlResult,
    solve_jet_a_pressure_swirl,
)
from fuelnozzle.lng import (
    EquilibriumFlowSettings,
    LNGNozzleGeometry,
    RelaxationFlowSettings,
)
from fuelnozzle.models import LNGComposition, ModelWarning, WarningSeverity
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider
from fuelnozzle.spray import (
    FlashSprayCalibration,
    FlashSpraySettings,
    Tier3FlashSpray,
    solve_lng_flash_spray,
)


class MissionFuelMasses(BaseModel):
    """Mission totals imported or copied from an external flight model."""

    model_config = ConfigDict(frozen=True)

    jet_a_kg: float = Field(ge=0.0)
    lng_kg: float = Field(ge=0.0)
    source: str = Field(min_length=1)


class StudySettings(BaseModel):
    """Cross-operating-point settings."""

    model_config = ConfigDict(frozen=True)

    mission_mass_relative_tolerance: float = Field(default=0.05, gt=0.0, lt=1.0)


@dataclass(frozen=True)
class IntegratedFuelMasses:
    jet_a_kg: float
    lng_kg: float


@dataclass(frozen=True)
class OperatingPointStudyResult:
    operating_point: OperatingPoint
    jet_a: PressureSwirlResult | None
    lng: Tier3FlashSpray | None
    warnings: tuple[ModelWarning, ...]


@dataclass(frozen=True)
class NozzleStudyResult:
    composition: LNGComposition
    coolprop_version: str
    operating_points: tuple[OperatingPointStudyResult, ...]
    integrated_fuel_masses: IntegratedFuelMasses | None
    mission_fuel_masses: MissionFuelMasses | None
    warnings: tuple[ModelWarning, ...]


def run_nozzle_study(
    operating_points: tuple[OperatingPoint, ...] | list[OperatingPoint],
    composition: LNGComposition,
    lng_geometry: LNGNozzleGeometry,
    jet_a_geometry: PressureSwirlGeometry,
    jet_a_properties: JetAProperties | JetAPropertyTable,
    *,
    feed_line: LNGFeedLine | None = None,
    equilibrium_settings: EquilibriumFlowSettings | None = None,
    relaxation_settings: RelaxationFlowSettings | None = None,
    spray_settings: FlashSpraySettings | None = None,
    spray_calibration: FlashSprayCalibration | None = None,
    mission_fuel_masses: MissionFuelMasses | None = None,
    study_settings: StudySettings | None = None,
) -> NozzleStudyResult:
    """Evaluate both fuel circuits at every user-specified operating point."""

    equilibrium_settings = equilibrium_settings or EquilibriumFlowSettings()
    relaxation_settings = relaxation_settings or RelaxationFlowSettings()
    spray_settings = spray_settings or FlashSpraySettings()
    study_settings = study_settings or StudySettings()
    points = tuple(operating_points)
    if not points:
        raise ValueError("At least one operating point is required")
    properties = CoolPropLNGProvider(composition)
    point_results: list[OperatingPointStudyResult] = []

    for point in points:
        jet_a_result = (
            solve_jet_a_pressure_swirl(point, jet_a_geometry, jet_a_properties)
            if point.jet_a_mass_flow_kg_s > 0.0
            else None
        )
        lng_result = (
            solve_lng_flash_spray(
                point,
                lng_geometry,
                properties,
                equilibrium_settings,
                relaxation_settings,
                spray_settings,
                spray_calibration,
                feed_line,
            )
            if point.lng_mass_flow_kg_s > 0.0
            else None
        )
        warnings = _deduplicate_warnings(
            (
                *(jet_a_result.warnings if jet_a_result is not None else ()),
                *(lng_result.warnings if lng_result is not None else ()),
            )
        )
        point_results.append(
            OperatingPointStudyResult(
                operating_point=point,
                jet_a=jet_a_result,
                lng=lng_result,
                warnings=warnings,
            )
        )

    integrated = _integrate_fuel_masses(points)
    study_warnings: list[ModelWarning] = []
    if mission_fuel_masses is not None:
        if integrated is None:
            study_warnings.append(
                ModelWarning(
                    code="MISSION_MASS_CHECK_UNAVAILABLE",
                    severity=WarningSeverity.INFO,
                    message=(
                        "Mission totals were provided, but at least one operating point lacks "
                        "duration; no mass-flow integration comparison was made."
                    ),
                )
            )
        else:
            study_warnings.extend(
                _mission_mass_warnings(
                    integrated,
                    mission_fuel_masses,
                    study_settings.mission_mass_relative_tolerance,
                )
            )

    return NozzleStudyResult(
        composition=composition,
        coolprop_version=properties.coolprop_version,
        operating_points=tuple(point_results),
        integrated_fuel_masses=integrated,
        mission_fuel_masses=mission_fuel_masses,
        warnings=tuple(study_warnings),
    )


def _integrate_fuel_masses(
    operating_points: tuple[OperatingPoint, ...],
) -> IntegratedFuelMasses | None:
    if any(point.duration_s is None for point in operating_points):
        return None
    jet_a_kg = sum(
        point.jet_a_mass_flow_kg_s * point.duration_s * point.flow_multiplier
        for point in operating_points
        if point.duration_s is not None
    )
    lng_kg = sum(
        point.lng_mass_flow_kg_s * point.duration_s * point.flow_multiplier
        for point in operating_points
        if point.duration_s is not None
    )
    return IntegratedFuelMasses(jet_a_kg=jet_a_kg, lng_kg=lng_kg)


def _mission_mass_warnings(
    integrated: IntegratedFuelMasses,
    mission: MissionFuelMasses,
    tolerance: float,
) -> list[ModelWarning]:
    warnings: list[ModelWarning] = []
    for fuel_name, calculated, expected in (
        ("JET_A", integrated.jet_a_kg, mission.jet_a_kg),
        ("LNG", integrated.lng_kg, mission.lng_kg),
    ):
        denominator = max(expected, 1.0e-12)
        relative_error = abs(calculated - expected) / denominator
        if relative_error > tolerance:
            warnings.append(
                ModelWarning(
                    code=f"{fuel_name}_MISSION_MASS_MISMATCH",
                    severity=WarningSeverity.WARNING,
                    message=(
                        f"Integrated stage mass differs from {mission.source} by "
                        f"{relative_error:.1%}. Check duration, per-engine flow, and multiplier."
                    ),
                )
            )
    return warnings


def _deduplicate_warnings(warnings: tuple[ModelWarning, ...]) -> tuple[ModelWarning, ...]:
    unique: dict[tuple[str, str], ModelWarning] = {}
    for warning in warnings:
        unique[(warning.code, warning.message)] = warning
    return tuple(unique.values())
