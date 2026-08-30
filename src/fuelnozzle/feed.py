"""Cryogenic LNG feed-line hydraulic and thermal model."""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, degrees, pi, radians, sin

from fluids.friction import friction_factor
from fluids.two_phase import two_phase_dP
from ht.conduction import cylindrical_heat_transfer
from ht.conv_internal import Nu_conv_internal
from pydantic import BaseModel, ConfigDict, Field

from fuelnozzle.models import ModelWarning, ThermodynamicState, WarningSeverity
from fuelnozzle.operating import OperatingPoint, fuel_pressure_budget
from fuelnozzle.properties import CoolPropLNGProvider

GRAVITY_M_S2 = 9.80665


class LNGFeedLine(BaseModel):
    """Uniform feed line between the LNG pump and injector manifold."""

    model_config = ConfigDict(frozen=True)

    length_m: float = Field(gt=0.0)
    inner_diameter_m: float = Field(gt=0.0)
    roughness_m: float = Field(default=1.5e-6, ge=0.0)
    elevation_change_m: float = 0.0
    minor_loss_coefficient: float = Field(default=0.0, ge=0.0)
    segments: int = Field(default=80, ge=4, le=2000)
    measured_heat_leak_w_per_m: float | None = Field(default=None, ge=0.0)
    ambient_temperature_k: float = Field(default=300.0, gt=0.0)
    external_heat_transfer_coefficient_w_m2_k: float = Field(default=10.0, gt=0.0)
    wall_thickness_m: float = Field(default=1.0e-3, gt=0.0)
    wall_conductivity_w_m_k: float = Field(default=15.0, gt=0.0)
    insulation_thickness_m: float = Field(default=20.0e-3, ge=0.0)
    insulation_conductivity_w_m_k: float = Field(default=0.025, gt=0.0)


@dataclass(frozen=True)
class FeedLinePoint:
    position_m: float
    pressure_pa: float
    temperature_k: float
    enthalpy_j_kg: float
    vapor_quality_mass: float
    heat_leak_w_per_m: float
    liquid_mole_fractions: dict[str, float] | None = None
    vapor_mole_fractions: dict[str, float] | None = None


@dataclass(frozen=True)
class LNGFeedLineResult:
    inlet_state: ThermodynamicState
    outlet_state: ThermodynamicState
    required_nozzle_inlet_pressure_pa: float
    pressure_drop_pa: float
    total_heat_leak_w: float
    first_two_phase_position_m: float | None
    path: tuple[FeedLinePoint, ...]
    warnings: tuple[ModelWarning, ...]


def solve_lng_feed_line(
    operating_point: OperatingPoint,
    line: LNGFeedLine,
    properties: CoolPropLNGProvider,
) -> LNGFeedLineResult:
    """March pressure and enthalpy through a heated cryogenic feed line."""

    if operating_point.lng_mass_flow_kg_s <= 0.0:
        raise ValueError("LNG feed-line calculation requires positive LNG mass flow")

    budget = fuel_pressure_budget(
        operating_point.lng_pump_outlet_pressure_pa,
        operating_point.lng_nozzle_pressure_drop_pa,
        operating_point.p3_pa,
    )
    inlet = properties.state_pt(
        operating_point.lng_pump_outlet_pressure_pa,
        operating_point.lng_pump_outlet_temperature_k,
    )
    pressure_pa = inlet.pressure_pa
    enthalpy_j_kg = inlet.enthalpy_j_kg
    temperature_hint_k = inlet.temperature_k
    segment_length = line.length_m / line.segments
    segment_elevation = line.elevation_change_m / line.segments
    angle_degrees = 0.0
    if line.elevation_change_m != 0.0:
        angle_degrees = (
            90.0
            if abs(line.elevation_change_m) >= line.length_m
            else degrees(asin(line.elevation_change_m / line.length_m))
        )
    area_m2 = pi * line.inner_diameter_m**2 / 4.0
    mass_flux = operating_point.lng_mass_flow_kg_s / area_m2
    total_heat_leak_w = 0.0
    first_two_phase_position: float | None = None
    warnings = list(budget.warnings)
    path: list[FeedLinePoint] = []

    for index in range(line.segments + 1):
        position_m = index * segment_length
        state = properties.state_ph(
            pressure_pa,
            enthalpy_j_kg,
            temperature_hint_k=temperature_hint_k,
        )
        temperature_hint_k = state.temperature_k
        quality = state.vapor_quality_mass or 0.0
        heat_leak_w_per_m = _heat_leak_per_length(line, state, mass_flux)
        path.append(
            FeedLinePoint(
                position_m=position_m,
                pressure_pa=pressure_pa,
                temperature_k=state.temperature_k,
                enthalpy_j_kg=enthalpy_j_kg,
                vapor_quality_mass=quality,
                heat_leak_w_per_m=heat_leak_w_per_m,
                liquid_mole_fractions=state.liquid_mole_fractions,
                vapor_mole_fractions=state.vapor_mole_fractions,
            )
        )
        if quality > 0.0 and first_two_phase_position is None:
            first_two_phase_position = position_m
        if index == line.segments:
            break

        if quality > 0.0:
            pressure_drop = _two_phase_segment_pressure_drop(
                operating_point.lng_mass_flow_kg_s,
                quality,
                pressure_pa,
                segment_length,
                angle_degrees,
                line,
                properties,
            )
        else:
            pressure_drop = _liquid_segment_pressure_drop(
                mass_flux,
                segment_length,
                segment_elevation,
                line,
                state,
            )

        pressure_pa -= pressure_drop
        if pressure_pa <= 0.0:
            raise ValueError("Calculated feed-line pressure became non-positive")
        enthalpy_j_kg += heat_leak_w_per_m * segment_length / operating_point.lng_mass_flow_kg_s
        total_heat_leak_w += heat_leak_w_per_m * segment_length

    outlet = properties.state_ph(
        pressure_pa,
        enthalpy_j_kg,
        temperature_hint_k=temperature_hint_k,
    )
    pressure_drop_pa = inlet.pressure_pa - outlet.pressure_pa

    pressure_error = outlet.pressure_pa - budget.nozzle_inlet_pressure_pa
    tolerance_pa = max(10_000.0, 0.05 * max(budget.available_feed_pressure_drop_pa, 0.0))
    if pressure_error < -tolerance_pa:
        warnings.append(
            ModelWarning(
                code="FEED_PRESSURE_DROP_EXCEEDS_BUDGET",
                severity=WarningSeverity.ERROR,
                message=(
                    "Predicted feed-line outlet pressure is below the nozzle inlet pressure "
                    "required by P3 and the specified nozzle pressure drop."
                ),
            )
        )
    elif abs(pressure_error) > tolerance_pa:
        warnings.append(
            ModelWarning(
                code="FEED_PRESSURE_BUDGET_MISMATCH",
                severity=WarningSeverity.WARNING,
                message=(
                    "Predicted feed-line pressure loss does not use the available pump-to-nozzle "
                    "pressure allowance within tolerance."
                ),
            )
        )
    if first_two_phase_position is not None:
        warnings.append(
            ModelWarning(
                code="UPSTREAM_TWO_PHASE_FLOW",
                severity=WarningSeverity.ERROR,
                message="Equilibrium vapor appears in the LNG feed line before the nozzle.",
            )
        )
    if properties.transport_fallback_used:
        saturation_components = ", ".join(
            sorted(properties.transport_saturation_fallback_components)
        )
        omitted_components = ", ".join(sorted(properties.transport_omitted_components))
        details = ""
        if saturation_components:
            details += (
                " Same-temperature saturated-liquid transport was used for: "
                f"{saturation_components}."
            )
        if omitted_components:
            details += (
                " Unavailable components were omitted and weights renormalized: "
                f"{omitted_components}."
            )
        warnings.append(
            ModelWarning(
                code="LNG_TRANSPORT_MIXING_FALLBACK",
                severity=WarningSeverity.WARNING,
                message=(
                    "CoolProp did not supply native LNG-mixture transport properties. "
                    "Component log-mole viscosity and mole-linear conductivity estimates "
                    f"were used; validate these against mixture data.{details}"
                ),
            )
        )

    return LNGFeedLineResult(
        inlet_state=inlet,
        outlet_state=outlet,
        required_nozzle_inlet_pressure_pa=budget.nozzle_inlet_pressure_pa,
        pressure_drop_pa=pressure_drop_pa,
        total_heat_leak_w=total_heat_leak_w,
        first_two_phase_position_m=first_two_phase_position,
        path=tuple(path),
        warnings=tuple(warnings),
    )


def _heat_leak_per_length(
    line: LNGFeedLine,
    state: ThermodynamicState,
    mass_flux_kg_m2_s: float,
) -> float:
    if line.measured_heat_leak_w_per_m is not None:
        return line.measured_heat_leak_w_per_m
    if (
        state.viscosity_pa_s is None
        or state.conductivity_w_m_k is None
        or state.cp_j_kg_k is None
    ):
        raise ValueError("CoolProp transport properties are required for calculated heat leak")

    reynolds = mass_flux_kg_m2_s * line.inner_diameter_m / state.viscosity_pa_s
    prandtl = state.cp_j_kg_k * state.viscosity_pa_s / state.conductivity_w_m_k
    nusselt = Nu_conv_internal(
        Re=reynolds,
        Pr=prandtl,
        eD=line.roughness_m / line.inner_diameter_m,
        Di=line.inner_diameter_m,
    )
    internal_h = nusselt * state.conductivity_w_m_k / line.inner_diameter_m
    thicknesses = [line.wall_thickness_m]
    conductivities = [line.wall_conductivity_w_m_k]
    if line.insulation_thickness_m > 0.0:
        thicknesses.append(line.insulation_thickness_m)
        conductivities.append(line.insulation_conductivity_w_m_k)
    result = cylindrical_heat_transfer(
        Ti=state.temperature_k,
        To=line.ambient_temperature_k,
        hi=internal_h,
        ho=line.external_heat_transfer_coefficient_w_m2_k,
        Di=line.inner_diameter_m,
        ts=thicknesses,
        ks=conductivities,
    )
    return max(0.0, -float(result["Q"]))


def _liquid_segment_pressure_drop(
    mass_flux_kg_m2_s: float,
    segment_length_m: float,
    segment_elevation_m: float,
    line: LNGFeedLine,
    state: ThermodynamicState,
) -> float:
    if state.viscosity_pa_s is None:
        raise ValueError(
            "LNG viscosity is unavailable for feed-line pressure loss at "
            f"P={state.pressure_pa:g} Pa, T={state.temperature_k:g} K "
            f"({state.transport_model})"
        )
    reynolds = mass_flux_kg_m2_s * line.inner_diameter_m / state.viscosity_pa_s
    darcy_factor = friction_factor(
        Re=reynolds,
        eD=line.roughness_m / line.inner_diameter_m,
    )
    segment_minor_k = line.minor_loss_coefficient / line.segments
    dynamic_pressure = mass_flux_kg_m2_s**2 / (2.0 * state.density_kg_m3)
    friction_and_minor = (
        darcy_factor * segment_length_m / line.inner_diameter_m + segment_minor_k
    ) * dynamic_pressure
    static_head = state.density_kg_m3 * GRAVITY_M_S2 * segment_elevation_m
    return friction_and_minor + static_head


def _two_phase_segment_pressure_drop(
    mass_flow_kg_s: float,
    quality: float,
    pressure_pa: float,
    segment_length_m: float,
    angle_degrees: float,
    line: LNGFeedLine,
    properties: CoolPropLNGProvider,
) -> float:
    liquid = properties.state_pq(pressure_pa, 0.0)
    vapor = properties.state_pq(pressure_pa, 1.0)
    if (
        liquid.viscosity_pa_s is None
        or vapor.viscosity_pa_s is None
        or liquid.surface_tension_n_m is None
    ):
        raise ValueError("Saturated transport properties are required for two-phase line loss")
    method = "Kim_Mudawar" if line.inner_diameter_m <= 6.22e-3 else "Chisholm"
    friction_drop = two_phase_dP(
        m=mass_flow_kg_s,
        x=min(1.0 - 1.0e-8, max(1.0e-8, quality)),
        rhol=liquid.density_kg_m3,
        rhog=vapor.density_kg_m3,
        mul=liquid.viscosity_pa_s,
        mug=vapor.viscosity_pa_s,
        sigma=liquid.surface_tension_n_m,
        D=line.inner_diameter_m,
        L=segment_length_m,
        roughness=line.roughness_m,
        Method=method,
    )
    mixture_density = 1.0 / (
        (1.0 - quality) / liquid.density_kg_m3 + quality / vapor.density_kg_m3
    )
    static_head = (
        mixture_density * GRAVITY_M_S2 * segment_length_m * sin(radians(angle_degrees))
    )
    return friction_drop + static_head
