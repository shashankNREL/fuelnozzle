"""Tier 3 LNG flash-spray screening and calibrated output models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import log

from pydantic import BaseModel, ConfigDict, Field

from fuelnozzle.feed import LNGFeedLine
from fuelnozzle.lng import (
    EquilibriumFlowSettings,
    FlashLocation,
    LNGNozzleGeometry,
    RelaxationFlowSettings,
    Tier2RelaxationFlow,
    solve_lng_relaxation_flow,
)
from fuelnozzle.models import ModelWarning, WarningSeverity
from fuelnozzle.operating import OperatingPoint
from fuelnozzle.properties import CoolPropLNGProvider, PropertyCalculationError


class FlashSprayRegime(StrEnum):
    """Reduced-order spray regime at the LNG orifice exit."""

    MECHANICAL = "mechanical_breakup"
    EXTERNAL_FLASH = "external_flash"
    TRANSITIONAL_FLASH = "transitional_flash"
    FULLY_FLASHING = "fully_flashing"
    UPSTREAM_TWO_PHASE = "upstream_two_phase"


class FlashSpraySettings(BaseModel):
    """Transparent thresholds used only for Tier 3 regime screening."""

    model_config = ConfigDict(frozen=True)

    fully_flashing_quality_threshold: float = Field(default=0.10, gt=0.0, lt=1.0)
    fully_flashing_pressure_ratio: float = Field(default=3.0, gt=1.0)
    fully_flashing_jakob_threshold: float = Field(default=0.10, gt=0.0)


class FlashSprayCalibration(BaseModel):
    """Hardware-specific flash-spray calibration supplied by the user.

    The reference values and exponents must come from a documented experiment
    or higher-fidelity calculation with comparable fluid and nozzle geometry.
    """

    model_config = ConfigDict(frozen=True)

    calibration_id: str = Field(min_length=1)
    reference_smd_m: float = Field(gt=0.0)
    reference_pressure_ratio: float = Field(gt=1.0)
    reference_flash_quality: float = Field(gt=0.0, lt=1.0)
    pressure_ratio_exponent: float = Field(default=0.5, ge=0.0)
    flash_quality_exponent: float = Field(default=0.5, ge=0.0)
    smd_relative_uncertainty: float = Field(default=0.40, gt=0.0, lt=1.0)
    reference_full_cone_angle_deg: float = Field(gt=0.0, lt=180.0)
    cone_angle_gain_deg: float = 10.0
    cone_angle_uncertainty_deg: float = Field(default=10.0, ge=0.0)


@dataclass(frozen=True)
class CFDSprayBoundary:
    """Boundary values suitable for a later VOF or Euler-Lagrange export."""

    mass_flow_kg_s: float
    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    velocity_m_s: float
    vapor_quality_mass: float
    effective_orifice_diameter_m: float
    number_of_orifices: int
    spray_regime: FlashSprayRegime


@dataclass(frozen=True)
class Tier3FlashSpray:
    """Flash-spray regime, calibrated ranges, and CFD boundary data."""

    tier2: Tier2RelaxationFlow
    settings: FlashSpraySettings
    regime: FlashSprayRegime
    saturation_to_chamber_pressure_ratio: float
    superheat_at_p3_k: float
    equilibrium_flash_fraction_mass: float
    actual_exit_vapor_quality_mass: float
    jakob_number: float | None
    smd_estimate_m: float | None
    smd_range_m: tuple[float, float] | None
    full_cone_angle_estimate_deg: float | None
    full_cone_angle_range_deg: tuple[float, float] | None
    calibration_id: str | None
    cfd_boundary: CFDSprayBoundary
    warnings: tuple[ModelWarning, ...]


def solve_lng_flash_spray(
    operating_point: OperatingPoint,
    geometry: LNGNozzleGeometry,
    properties: CoolPropLNGProvider,
    equilibrium_settings: EquilibriumFlowSettings | None = None,
    relaxation_settings: RelaxationFlowSettings | None = None,
    spray_settings: FlashSpraySettings | None = None,
    calibration: FlashSprayCalibration | None = None,
    feed_line: LNGFeedLine | None = None,
) -> Tier3FlashSpray:
    """Tier 3: classify flash breakup and apply optional hardware calibration."""

    equilibrium_settings = equilibrium_settings or EquilibriumFlowSettings()
    relaxation_settings = relaxation_settings or RelaxationFlowSettings()
    spray_settings = spray_settings or FlashSpraySettings()
    tier2 = solve_lng_relaxation_flow(
        operating_point,
        geometry,
        properties,
        equilibrium_settings,
        relaxation_settings,
        feed_line,
    )
    tier0 = tier2.tier1.tier0
    pressure_ratio = (
        tier0.bubble_pressure_at_inlet_temperature_pa / operating_point.p3_pa
    )
    equilibrium_fraction = tier0.equilibrium_flash_fraction_at_p3
    superheat_k = 0.0
    jakob_number: float | None = None
    try:
        saturated_liquid = properties.bubble_state_at_pressure(operating_point.p3_pa)
        saturated_vapor = properties.dew_state_at_pressure(operating_point.p3_pa)
        superheat_k = max(
            0.0,
            tier0.nozzle_inlet_state.temperature_k - saturated_liquid.temperature_k,
        )
        latent_heat = saturated_vapor.enthalpy_j_kg - saturated_liquid.enthalpy_j_kg
        if tier0.nozzle_inlet_state.cp_j_kg_k is not None and latent_heat > 0.0:
            jakob_number = (
                tier0.nozzle_inlet_state.cp_j_kg_k * superheat_k / latent_heat
            )
    except PropertyCalculationError:
        pass

    regime = _classify_regime(tier2, pressure_ratio, jakob_number, spray_settings)
    warnings = list(tier2.warnings)
    smd_estimate: float | None = None
    smd_range: tuple[float, float] | None = None
    angle_estimate: float | None = None
    angle_range: tuple[float, float] | None = None
    calibration_id: str | None = None

    if calibration is None:
        warnings.append(
            ModelWarning(
                code="FLASH_SPRAY_CALIBRATION_REQUIRED",
                severity=WarningSeverity.WARNING,
                message=(
                    "Tier 3 does not report LNG SMD or cone angle without a traceable "
                    "cryogenic spray calibration for comparable geometry and conditions."
                ),
            )
        )
    else:
        calibration_id = calibration.calibration_id
        driving_quality = max(
            1.0e-6,
            tier2.actual_exit_vapor_quality_mass,
            equilibrium_fraction,
        )
        ratio_factor = max(1.0, pressure_ratio) / calibration.reference_pressure_ratio
        quality_factor = driving_quality / calibration.reference_flash_quality
        smd_estimate = calibration.reference_smd_m * (
            ratio_factor ** (-calibration.pressure_ratio_exponent)
        ) * (quality_factor ** (-calibration.flash_quality_exponent))
        smd_range = (
            smd_estimate * (1.0 - calibration.smd_relative_uncertainty),
            smd_estimate * (1.0 + calibration.smd_relative_uncertainty),
        )
        angle_estimate = calibration.reference_full_cone_angle_deg + (
            calibration.cone_angle_gain_deg * log(max(1.0e-9, ratio_factor))
        )
        angle_estimate = min(180.0, max(0.0, angle_estimate))
        angle_range = (
            max(0.0, angle_estimate - calibration.cone_angle_uncertainty_deg),
            min(180.0, angle_estimate + calibration.cone_angle_uncertainty_deg),
        )

    exit_point = tier2.path[-1]
    exit_temperature_k = tier2.tier1.path[-1].temperature_k
    effective_diameter = geometry.orifice_diameter_m or tier2.required_orifice_diameter_m
    cfd_boundary = CFDSprayBoundary(
        mass_flow_kg_s=operating_point.lng_mass_flow_kg_s,
        pressure_pa=exit_point.pressure_pa,
        temperature_k=exit_temperature_k,
        density_kg_m3=exit_point.actual_density_kg_m3,
        velocity_m_s=exit_point.velocity_m_s,
        vapor_quality_mass=exit_point.actual_vapor_quality_mass,
        effective_orifice_diameter_m=effective_diameter,
        number_of_orifices=geometry.number_of_orifices,
        spray_regime=regime,
    )

    return Tier3FlashSpray(
        tier2=tier2,
        settings=spray_settings,
        regime=regime,
        saturation_to_chamber_pressure_ratio=pressure_ratio,
        superheat_at_p3_k=superheat_k,
        equilibrium_flash_fraction_mass=equilibrium_fraction,
        actual_exit_vapor_quality_mass=tier2.actual_exit_vapor_quality_mass,
        jakob_number=jakob_number,
        smd_estimate_m=smd_estimate,
        smd_range_m=smd_range,
        full_cone_angle_estimate_deg=angle_estimate,
        full_cone_angle_range_deg=angle_range,
        calibration_id=calibration_id,
        cfd_boundary=cfd_boundary,
        warnings=tuple(warnings),
    )


def _classify_regime(
    tier2: Tier2RelaxationFlow,
    pressure_ratio: float,
    jakob_number: float | None,
    settings: FlashSpraySettings,
) -> FlashSprayRegime:
    if tier2.actual_flash_location == FlashLocation.UPSTREAM:
        return FlashSprayRegime.UPSTREAM_TWO_PHASE
    if tier2.tier1.tier0.flash_location == FlashLocation.NONE:
        return FlashSprayRegime.MECHANICAL
    if tier2.actual_flash_location == FlashLocation.EXTERNAL:
        return FlashSprayRegime.EXTERNAL_FLASH
    if (
        tier2.actual_exit_vapor_quality_mass
        >= settings.fully_flashing_quality_threshold
        or pressure_ratio >= settings.fully_flashing_pressure_ratio
        or (
            jakob_number is not None
            and jakob_number >= settings.fully_flashing_jakob_threshold
        )
    ):
        return FlashSprayRegime.FULLY_FLASHING
    return FlashSprayRegime.TRANSITIONAL_FLASH
