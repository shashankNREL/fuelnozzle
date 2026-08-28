"""Jet-A pressure-swirl injector sizing and spray screening."""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, degrees, pi, sqrt

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator

from fuelnozzle.models import ModelWarning, WarningSeverity
from fuelnozzle.operating import OperatingPoint, PressureBudget, fuel_pressure_budget

AIR_GAS_CONSTANT_J_KG_K = 287.05


class JetAProperties(BaseModel):
    """Measured or declared Jet-A properties at the nozzle inlet temperature."""

    model_config = ConfigDict(frozen=True)

    density_kg_m3: float = Field(gt=0.0)
    viscosity_pa_s: float = Field(gt=0.0)
    surface_tension_n_m: float = Field(gt=0.0)
    vapor_pressure_pa: float | None = Field(default=None, ge=0.0)
    source: str = Field(min_length=1)


class JetAPropertyTable(BaseModel):
    """Temperature-indexed measured Jet-A properties with linear interpolation."""

    model_config = ConfigDict(frozen=True)

    temperature_k: tuple[float, ...]
    density_kg_m3: tuple[float, ...]
    viscosity_pa_s: tuple[float, ...]
    surface_tension_n_m: tuple[float, ...]
    vapor_pressure_pa: tuple[float, ...] | None = None
    source: str = Field(min_length=1)

    @model_validator(mode="after")
    def validate_table(self) -> JetAPropertyTable:
        count = len(self.temperature_k)
        if count < 2:
            raise ValueError("Jet-A property table requires at least two temperatures")
        arrays = (self.density_kg_m3, self.viscosity_pa_s, self.surface_tension_n_m)
        if any(len(values) != count for values in arrays):
            raise ValueError("All Jet-A property arrays must have the same length")
        if self.vapor_pressure_pa is not None and len(self.vapor_pressure_pa) != count:
            raise ValueError("Jet-A vapor-pressure array must match the temperature array")
        if any(
            right <= left
            for left, right in zip(self.temperature_k, self.temperature_k[1:], strict=False)
        ):
            raise ValueError("Jet-A property temperatures must be strictly increasing")
        if any(value <= 0.0 for values in arrays for value in values):
            raise ValueError("Jet-A density, viscosity, and surface tension must be positive")
        return self

    def at_temperature(
        self, temperature_k: float
    ) -> tuple[JetAProperties, tuple[ModelWarning, ...]]:
        warnings: list[ModelWarning] = []
        if temperature_k < self.temperature_k[0] or temperature_k > self.temperature_k[-1]:
            warnings.append(
                ModelWarning(
                    code="JET_A_PROPERTY_EXTRAPOLATION",
                    severity=WarningSeverity.WARNING,
                    message=(
                        "Jet-A temperature is outside the supplied property table; endpoint "
                        "values are used."
                    ),
                )
            )

        def interpolate(values: tuple[float, ...]) -> float:
            return float(np.interp(temperature_k, self.temperature_k, values))

        vapor_pressure = (
            interpolate(self.vapor_pressure_pa)
            if self.vapor_pressure_pa is not None
            else None
        )
        return (
            JetAProperties(
                density_kg_m3=interpolate(self.density_kg_m3),
                viscosity_pa_s=interpolate(self.viscosity_pa_s),
                surface_tension_n_m=interpolate(self.surface_tension_n_m),
                vapor_pressure_pa=vapor_pressure,
                source=self.source,
            ),
            tuple(warnings),
        )


class PressureSwirlGeometry(BaseModel):
    """Simplex pressure-swirl geometry and transparent closure coefficients."""

    model_config = ConfigDict(frozen=True)

    number_of_inlet_ports: int = Field(default=4, ge=1)
    inlet_port_diameter_m: float = Field(gt=0.0)
    inlet_tangency_radius_m: float = Field(gt=0.0)
    swirl_chamber_radius_m: float = Field(gt=0.0)
    swirl_chamber_length_m: float = Field(gt=0.0)
    exit_orifice_diameter_m: float | None = Field(default=None, gt=0.0)
    exit_orifice_length_m: float = Field(default=1.0e-3, gt=0.0)
    design_discharge_coefficient: float = Field(default=0.65, gt=0.0, le=1.0)
    velocity_coefficient: float = Field(default=0.95, gt=0.0, le=1.0)
    angular_momentum_efficiency: float = Field(default=0.75, gt=0.0, le=1.0)
    smd_calibration_coefficient: float | None = Field(default=None, gt=0.0)
    smd_relative_uncertainty: float = Field(default=0.50, gt=0.0, lt=1.0)


@dataclass(frozen=True)
class PressureSwirlResult:
    """Hydraulic and reduced-order spray result for a Jet-A simplex atomizer."""

    pressure_budget: PressureBudget
    properties: JetAProperties
    required_exit_area_m2: float
    required_exit_diameter_m: float
    modeled_exit_diameter_m: float
    predicted_mass_flow_kg_s: float
    effective_discharge_coefficient: float
    inlet_port_velocity_m_s: float
    axial_exit_velocity_m_s: float
    tangential_exit_velocity_m_s: float
    air_core_radius_m: float
    liquid_film_thickness_m: float
    full_cone_angle_deg: float
    minimum_internal_pressure_pa: float
    cavitation_pressure_margin_pa: float | None
    liquid_reynolds_number: float
    sheet_weber_number: float
    sheet_ohnesorge_number: float
    ambient_air_density_kg_m3: float
    smd_estimate_m: float | None
    smd_range_m: tuple[float, float] | None
    warnings: tuple[ModelWarning, ...]


def solve_jet_a_pressure_swirl(
    operating_point: OperatingPoint,
    geometry: PressureSwirlGeometry,
    properties: JetAProperties | JetAPropertyTable,
) -> PressureSwirlResult:
    """Size and evaluate a Jet-A pressure-swirl injector at one operating point."""

    if operating_point.jet_a_mass_flow_kg_s <= 0.0:
        raise ValueError("Pressure-swirl calculation requires positive Jet-A mass flow")
    if isinstance(properties, JetAPropertyTable):
        local_properties, property_warnings = properties.at_temperature(
            operating_point.jet_a_nozzle_inlet_temperature_k
        )
    else:
        local_properties = properties
        property_warnings = ()

    budget = fuel_pressure_budget(
        operating_point.jet_a_pump_outlet_pressure_pa,
        operating_point.jet_a_nozzle_pressure_drop_pa,
        operating_point.p3_pa,
    )
    warnings = list(budget.warnings) + list(property_warnings)
    density = local_properties.density_kg_m3
    pressure_drop = operating_point.jet_a_nozzle_pressure_drop_pa
    ideal_velocity = sqrt(2.0 * pressure_drop / density)
    required_area = operating_point.jet_a_mass_flow_kg_s / (
        geometry.design_discharge_coefficient * density * ideal_velocity
    )
    required_diameter = sqrt(4.0 * required_area / pi)
    modeled_diameter = geometry.exit_orifice_diameter_m or required_diameter
    exit_radius = modeled_diameter / 2.0
    exit_area = pi * exit_radius**2
    predicted_mass_flow = (
        geometry.design_discharge_coefficient * density * exit_area * ideal_velocity
    )
    effective_discharge_coefficient = operating_point.jet_a_mass_flow_kg_s / (
        density * exit_area * ideal_velocity
    )

    inlet_area = (
        geometry.number_of_inlet_ports * pi * geometry.inlet_port_diameter_m**2 / 4.0
    )
    inlet_port_velocity = operating_point.jet_a_mass_flow_kg_s / (density * inlet_area)
    total_exit_velocity = geometry.velocity_coefficient * ideal_velocity
    unconstrained_tangential_velocity = (
        geometry.angular_momentum_efficiency
        * inlet_port_velocity
        * geometry.inlet_tangency_radius_m
        / exit_radius
    )
    minimum_axial_velocity = min(
        0.99 * total_exit_velocity,
        1.02 * geometry.design_discharge_coefficient * ideal_velocity,
    )
    annulus_limited_tangential_velocity = sqrt(
        max(0.0, total_exit_velocity**2 - minimum_axial_velocity**2)
    )
    tangential_velocity = min(
        0.98 * total_exit_velocity,
        annulus_limited_tangential_velocity,
        unconstrained_tangential_velocity,
    )
    if tangential_velocity < unconstrained_tangential_velocity:
        warnings.append(
            ModelWarning(
                code="SWIRL_ENERGY_LIMIT",
                severity=WarningSeverity.WARNING,
                message=(
                    "Port angular momentum exceeds the available nozzle pressure energy; "
                    "tangential velocity was capped at 98% of total exit velocity."
                ),
            )
        )
    axial_velocity = sqrt(max(1.0e-12, total_exit_velocity**2 - tangential_velocity**2))
    liquid_exit_area = operating_point.jet_a_mass_flow_kg_s / (density * axial_velocity)
    if liquid_exit_area >= exit_area:
        air_core_radius = 0.0
        film_thickness = exit_radius
        warnings.append(
            ModelWarning(
                code="NO_STABLE_AIR_CORE",
                severity=WarningSeverity.ERROR,
                message=(
                    "Required liquid annulus area equals or exceeds the exit area; this "
                    "geometry cannot support the predicted hollow-cone state."
                ),
            )
        )
    else:
        air_core_radius = sqrt(exit_radius**2 - liquid_exit_area / pi)
        film_thickness = exit_radius - air_core_radius

    full_cone_angle = 2.0 * degrees(atan2(tangential_velocity, axial_velocity))
    minimum_internal_pressure = (
        budget.nozzle_inlet_pressure_pa - 0.5 * density * tangential_velocity**2
    )
    cavitation_margin = (
        minimum_internal_pressure - local_properties.vapor_pressure_pa
        if local_properties.vapor_pressure_pa is not None
        else None
    )
    if cavitation_margin is not None and cavitation_margin <= 0.0:
        warnings.append(
            ModelWarning(
                code="JET_A_CAVITATION_RISK",
                severity=WarningSeverity.ERROR,
                message="Estimated minimum internal pressure is below Jet-A vapor pressure.",
            )
        )
    elif local_properties.vapor_pressure_pa is None:
        warnings.append(
            ModelWarning(
                code="JET_A_VAPOR_PRESSURE_NOT_PROVIDED",
                severity=WarningSeverity.INFO,
                message="Cavitation margin was not evaluated because vapor pressure is absent.",
            )
        )

    reynolds = density * total_exit_velocity * max(film_thickness, 1.0e-12) / (
        local_properties.viscosity_pa_s
    )
    weber = (
        density
        * total_exit_velocity**2
        * max(film_thickness, 1.0e-12)
        / local_properties.surface_tension_n_m
    )
    ohnesorge = local_properties.viscosity_pa_s / sqrt(
        density
        * local_properties.surface_tension_n_m
        * max(film_thickness, 1.0e-12)
    )
    ambient_air_density = operating_point.p3_pa / (
        AIR_GAS_CONSTANT_J_KG_K * operating_point.t3_k
    )

    smd_estimate: float | None = None
    smd_range: tuple[float, float] | None = None
    if geometry.smd_calibration_coefficient is None:
        warnings.append(
            ModelWarning(
                code="PRESSURE_SWIRL_SMD_CALIBRATION_REQUIRED",
                severity=WarningSeverity.WARNING,
                message=(
                    "Jet-A SMD is suppressed until a hardware-specific sheet-breakup "
                    "calibration coefficient is supplied."
                ),
            )
        )
    else:
        capillary_sheet_scale = sqrt(
            local_properties.surface_tension_n_m
            * max(film_thickness, 1.0e-12)
            / (density * total_exit_velocity**2)
        )
        smd_estimate = (
            geometry.smd_calibration_coefficient
            * capillary_sheet_scale
            * (1.0 + 3.0 * ohnesorge)
        )
        smd_range = (
            smd_estimate * (1.0 - geometry.smd_relative_uncertainty),
            smd_estimate * (1.0 + geometry.smd_relative_uncertainty),
        )

    return PressureSwirlResult(
        pressure_budget=budget,
        properties=local_properties,
        required_exit_area_m2=required_area,
        required_exit_diameter_m=required_diameter,
        modeled_exit_diameter_m=modeled_diameter,
        predicted_mass_flow_kg_s=predicted_mass_flow,
        effective_discharge_coefficient=effective_discharge_coefficient,
        inlet_port_velocity_m_s=inlet_port_velocity,
        axial_exit_velocity_m_s=axial_velocity,
        tangential_exit_velocity_m_s=tangential_velocity,
        air_core_radius_m=air_core_radius,
        liquid_film_thickness_m=film_thickness,
        full_cone_angle_deg=full_cone_angle,
        minimum_internal_pressure_pa=minimum_internal_pressure,
        cavitation_pressure_margin_pa=cavitation_margin,
        liquid_reynolds_number=reynolds,
        sheet_weber_number=weber,
        sheet_ohnesorge_number=ohnesorge,
        ambient_air_density_kg_m3=ambient_air_density,
        smd_estimate_m=smd_estimate,
        smd_range_m=smd_range,
        warnings=tuple(warnings),
    )
