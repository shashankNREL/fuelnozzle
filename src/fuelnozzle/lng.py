"""Reduced-order liquid-LNG nozzle screening and critical-flow models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import exp, pi, sqrt

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from fuelnozzle.feed import LNGFeedLine, LNGFeedLineResult, solve_lng_feed_line
from fuelnozzle.models import ModelWarning, ThermodynamicState, WarningSeverity
from fuelnozzle.operating import OperatingPoint, PressureBudget, fuel_pressure_budget
from fuelnozzle.properties import CoolPropLNGProvider


class FlashLocation(StrEnum):
    """Equilibrium location where LNG first permits a vapor phase."""

    NONE = "none"
    UPSTREAM = "upstream"
    INTERNAL = "internal"
    EXIT = "exit"
    EXTERNAL = "external"


class LNGNozzleGeometry(BaseModel):
    """Equivalent parallel short-hole geometry for the LNG circuit."""

    model_config = ConfigDict(frozen=True)

    number_of_orifices: int = Field(default=1, ge=1)
    orifice_diameter_m: float | None = Field(default=None, gt=0.0)
    orifice_length_m: float = Field(default=1.0e-3, gt=0.0)
    discharge_coefficient: float = Field(default=0.80, gt=0.0, le=1.0)

    @property
    def geometric_area_m2(self) -> float | None:
        if self.orifice_diameter_m is None:
            return None
        return self.number_of_orifices * pi * self.orifice_diameter_m**2 / 4.0


class EquilibriumFlowSettings(BaseModel):
    """Numerical settings for the Tier 1 pressure-path search."""

    model_config = ConfigDict(frozen=True)

    pressure_steps: int = Field(default=180, ge=40, le=2000)
    choking_pressure_tolerance: float = Field(default=0.005, gt=0.0, lt=0.1)


class RelaxationFlowSettings(BaseModel):
    """Tier 2 finite-rate phase-change assumptions.

    Both parameters are calibration quantities. A larger relaxation time or
    nucleation delay preserves a metastable liquid state farther downstream.
    """

    model_config = ConfigDict(frozen=True)

    relaxation_time_s: float = Field(default=50.0e-6, gt=0.0)
    nucleation_pressure_delay_pa: float = Field(default=0.0, ge=0.0)
    quality_onset_threshold: float = Field(default=1.0e-6, gt=0.0, lt=0.1)


@dataclass(frozen=True)
class Tier0FlashScreen:
    """Thermodynamic flash-risk screen before a nozzle-flow model is selected."""

    pressure_budget: PressureBudget
    feed_line: LNGFeedLineResult | None
    pump_state: ThermodynamicState
    nozzle_inlet_state: ThermodynamicState
    bubble_temperature_at_inlet_pressure_k: float
    bubble_pressure_at_inlet_temperature_pa: float
    temperature_subcooling_margin_k: float
    pressure_subcooling_margin_pa: float
    equilibrium_flash_fraction_at_p3: float
    flash_location: FlashLocation
    equilibrium_flash_onset_pressure_pa: float | None
    equilibrium_flash_onset_fraction_of_length: float | None
    warnings: tuple[ModelWarning, ...]


@dataclass(frozen=True)
class EquilibriumPathPoint:
    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    velocity_m_s: float
    mass_flux_kg_m2_s: float
    vapor_quality_mass: float


@dataclass(frozen=True)
class Tier1EquilibriumFlow:
    """Single-phase and homogeneous-equilibrium mass-flow bounds."""

    tier0: Tier0FlashScreen
    single_phase_mass_flux_kg_m2_s: float
    hem_operating_mass_flux_kg_m2_s: float
    critical_mass_flux_kg_m2_s: float
    critical_pressure_pa: float
    critical_vapor_quality_mass: float
    is_choked: bool
    required_geometric_area_m2: float
    required_orifice_diameter_m: float
    predicted_mass_flow_kg_s: float | None
    path: tuple[EquilibriumPathPoint, ...]
    warnings: tuple[ModelWarning, ...]


@dataclass(frozen=True)
class RelaxationPathPoint:
    position_m: float
    pressure_pa: float
    velocity_m_s: float
    equilibrium_vapor_quality_mass: float
    actual_vapor_quality_mass: float
    actual_density_kg_m3: float
    mass_flux_kg_m2_s: float
    elapsed_time_s: float


@dataclass(frozen=True)
class Tier2RelaxationFlow:
    """Finite-rate flashing result bounded by Tier 1 SPI and HEM solutions."""

    tier1: Tier1EquilibriumFlow
    settings: RelaxationFlowSettings
    nucleation_threshold_pressure_pa: float
    actual_flash_location: FlashLocation
    actual_flash_onset_pressure_pa: float | None
    actual_flash_onset_fraction_of_length: float | None
    residence_time_s: float
    equilibrium_exit_vapor_quality_mass: float
    actual_exit_vapor_quality_mass: float
    operating_mass_flux_kg_m2_s: float
    critical_mass_flux_kg_m2_s: float
    critical_pressure_pa: float
    is_choked: bool
    required_geometric_area_m2: float
    required_orifice_diameter_m: float
    predicted_mass_flow_kg_s: float | None
    path: tuple[RelaxationPathPoint, ...]
    warnings: tuple[ModelWarning, ...]


def screen_lng_flash(
    operating_point: OperatingPoint,
    properties: CoolPropLNGProvider,
    feed_line: LNGFeedLine | None = None,
) -> Tier0FlashScreen:
    """Tier 0: determine subcooling and equilibrium saturation-crossing location."""

    budget = fuel_pressure_budget(
        operating_point.lng_pump_outlet_pressure_pa,
        operating_point.lng_nozzle_pressure_drop_pa,
        operating_point.p3_pa,
    )
    line_result = (
        solve_lng_feed_line(operating_point, feed_line, properties)
        if feed_line is not None
        else None
    )
    pump_state = (
        line_result.inlet_state
        if line_result is not None
        else properties.state_pt(
            operating_point.lng_pump_outlet_pressure_pa,
            operating_point.lng_pump_outlet_temperature_k,
        )
    )
    inlet_enthalpy = (
        line_result.outlet_state.enthalpy_j_kg
        if line_result is not None
        else pump_state.enthalpy_j_kg
    )
    inlet_state = properties.state_ph(
        budget.nozzle_inlet_pressure_pa,
        inlet_enthalpy,
        temperature_hint_k=(
            line_result.outlet_state.temperature_k
            if line_result is not None
            else pump_state.temperature_k
        ),
    )
    bubble_at_inlet = properties.bubble_state_at_pressure(budget.nozzle_inlet_pressure_pa)
    bubble_pressure = properties.bubble_pressure_at_temperature(inlet_state.temperature_k)
    temperature_margin = bubble_at_inlet.temperature_k - inlet_state.temperature_k
    pressure_margin = budget.nozzle_inlet_pressure_pa - bubble_pressure

    outlet_equilibrium = properties.state_ph(
        operating_point.p3_pa,
        inlet_enthalpy,
        temperature_hint_k=inlet_state.temperature_k,
    )
    equilibrium_flash_fraction = outlet_equilibrium.vapor_quality_mass or 0.0

    warnings = list(budget.warnings)
    if line_result is not None:
        warnings.extend(line_result.warnings)
    onset_pressure: float | None = None
    onset_fraction: float | None = None
    if inlet_state.phase == "two_phase" or pressure_margin <= 0.0:
        flash_location = FlashLocation.UPSTREAM
        onset_pressure = budget.nozzle_inlet_pressure_pa
        onset_fraction = 0.0
        warnings.append(
            ModelWarning(
                code="TWO_PHASE_AT_NOZZLE_INLET",
                severity=WarningSeverity.ERROR,
                message="The equilibrium LNG state is two-phase at or before the nozzle inlet.",
            )
        )
    elif bubble_pressure <= operating_point.p3_pa:
        flash_location = FlashLocation.NONE
    else:
        onset_pressure = bubble_pressure
        pressure_span = budget.nozzle_inlet_pressure_pa - operating_point.p3_pa
        onset_fraction = (budget.nozzle_inlet_pressure_pa - onset_pressure) / pressure_span
        if onset_fraction >= 0.98:
            flash_location = FlashLocation.EXIT
        else:
            flash_location = FlashLocation.INTERNAL
        warnings.append(
            ModelWarning(
                code="EQUILIBRIUM_FLASH_POSSIBLE",
                severity=WarningSeverity.WARNING,
                message=(
                    "Local pressure crosses the LNG bubble pressure. Tier 2 determines "
                    "whether nucleation delay moves actual flashing downstream."
                ),
            )
        )

    return Tier0FlashScreen(
        pressure_budget=budget,
        feed_line=line_result,
        pump_state=pump_state,
        nozzle_inlet_state=inlet_state,
        bubble_temperature_at_inlet_pressure_k=bubble_at_inlet.temperature_k,
        bubble_pressure_at_inlet_temperature_pa=bubble_pressure,
        temperature_subcooling_margin_k=temperature_margin,
        pressure_subcooling_margin_pa=pressure_margin,
        equilibrium_flash_fraction_at_p3=equilibrium_flash_fraction,
        flash_location=flash_location,
        equilibrium_flash_onset_pressure_pa=onset_pressure,
        equilibrium_flash_onset_fraction_of_length=onset_fraction,
        warnings=tuple(warnings),
    )


def solve_lng_equilibrium_flow(
    operating_point: OperatingPoint,
    geometry: LNGNozzleGeometry,
    properties: CoolPropLNGProvider,
    settings: EquilibriumFlowSettings | None = None,
    feed_line: LNGFeedLine | None = None,
) -> Tier1EquilibriumFlow:
    """Tier 1: size a nozzle with SPI and isentropic HEM mass-flux bounds."""

    settings = settings or EquilibriumFlowSettings()
    tier0 = screen_lng_flash(operating_point, properties, feed_line)
    inlet = tier0.nozzle_inlet_state
    inlet_pressure = tier0.pressure_budget.nozzle_inlet_pressure_pa
    outlet_pressure = operating_point.p3_pa
    pressures = np.linspace(inlet_pressure, outlet_pressure, settings.pressure_steps)

    path: list[EquilibriumPathPoint] = []
    temperature_hint = inlet.temperature_k
    for pressure_pa in pressures:
        state = properties.state_ps(
            float(pressure_pa),
            inlet.entropy_j_kg_k,
            temperature_hint_k=temperature_hint,
        )
        temperature_hint = state.temperature_k
        kinetic_energy_j_kg = max(0.0, inlet.enthalpy_j_kg - state.enthalpy_j_kg)
        velocity = sqrt(2.0 * kinetic_energy_j_kg)
        mass_flux = state.density_kg_m3 * velocity
        path.append(
            EquilibriumPathPoint(
                pressure_pa=float(pressure_pa),
                temperature_k=state.temperature_k,
                density_kg_m3=state.density_kg_m3,
                velocity_m_s=velocity,
                mass_flux_kg_m2_s=mass_flux,
                vapor_quality_mass=state.vapor_quality_mass or 0.0,
            )
        )

    fluxes = np.asarray([point.mass_flux_kg_m2_s for point in path])
    critical_index = int(np.argmax(fluxes))
    critical = path[critical_index]
    pressure_tolerance = settings.choking_pressure_tolerance * inlet_pressure
    is_choked = critical.pressure_pa > outlet_pressure + pressure_tolerance
    operating_flux = critical.mass_flux_kg_m2_s if is_choked else path[-1].mass_flux_kg_m2_s

    single_phase_flux = sqrt(
        2.0
        * inlet.density_kg_m3
        * tier0.pressure_budget.nozzle_pressure_drop_pa
    )
    if operating_flux <= 0.0:
        raise ValueError("HEM operating mass flux is zero; verify the pressure boundary conditions")

    required_area = operating_point.lng_mass_flow_kg_s / (
        geometry.discharge_coefficient * operating_flux
    )
    required_diameter = sqrt(
        4.0 * required_area / (pi * geometry.number_of_orifices)
    )
    geometric_area = geometry.geometric_area_m2
    predicted_mass_flow = (
        geometry.discharge_coefficient * geometric_area * operating_flux
        if geometric_area is not None
        else None
    )

    warnings = list(tier0.warnings)
    if is_choked:
        warnings.append(
            ModelWarning(
                code="TWO_PHASE_CHOKING",
                severity=WarningSeverity.WARNING,
                message=(
                    "The HEM mass flux reaches a maximum above P3; mass flow is insensitive "
                    "to further downstream-pressure reduction in the equilibrium limit."
                ),
            )
        )
    if tier0.flash_location != FlashLocation.NONE:
        warnings.append(
            ModelWarning(
                code="HEM_IS_EQUILIBRIUM_BOUND",
                severity=WarningSeverity.INFO,
                message=(
                    "HEM assumes instantaneous phase equilibrium and normally provides the "
                    "early-flashing/lower-mass-flux side of the design bracket."
                ),
            )
        )

    return Tier1EquilibriumFlow(
        tier0=tier0,
        single_phase_mass_flux_kg_m2_s=single_phase_flux,
        hem_operating_mass_flux_kg_m2_s=operating_flux,
        critical_mass_flux_kg_m2_s=critical.mass_flux_kg_m2_s,
        critical_pressure_pa=critical.pressure_pa,
        critical_vapor_quality_mass=critical.vapor_quality_mass,
        is_choked=is_choked,
        required_geometric_area_m2=required_area,
        required_orifice_diameter_m=required_diameter,
        predicted_mass_flow_kg_s=predicted_mass_flow,
        path=tuple(path),
        warnings=tuple(warnings),
    )


def solve_lng_relaxation_flow(
    operating_point: OperatingPoint,
    geometry: LNGNozzleGeometry,
    properties: CoolPropLNGProvider,
    equilibrium_settings: EquilibriumFlowSettings | None = None,
    relaxation_settings: RelaxationFlowSettings | None = None,
    feed_line: LNGFeedLine | None = None,
) -> Tier2RelaxationFlow:
    """Tier 2: integrate delayed vapor generation along the short nozzle."""

    equilibrium_settings = equilibrium_settings or EquilibriumFlowSettings()
    relaxation_settings = relaxation_settings or RelaxationFlowSettings()
    tier1 = solve_lng_equilibrium_flow(
        operating_point,
        geometry,
        properties,
        equilibrium_settings,
        feed_line,
    )
    tier0 = tier1.tier0
    inlet_pressure = tier0.pressure_budget.nozzle_inlet_pressure_pa
    outlet_pressure = operating_point.p3_pa
    nucleation_threshold = max(
        0.0,
        tier0.bubble_pressure_at_inlet_temperature_pa
        - relaxation_settings.nucleation_pressure_delay_pa,
    )
    number_of_points = len(tier1.path)
    positions = np.linspace(0.0, geometry.orifice_length_m, number_of_points)

    actual_quality = 0.0
    elapsed_time = 0.0
    actual_onset_pressure: float | None = None
    actual_onset_fraction: float | None = None
    path: list[RelaxationPathPoint] = []
    previous_position = 0.0
    previous_velocity = 0.0

    for index, (position_m, equilibrium) in enumerate(
        zip(positions, tier1.path, strict=True)
    ):
        if index > 0:
            distance = float(position_m) - previous_position
            average_velocity = max(0.5, 0.5 * (previous_velocity + equilibrium.velocity_m_s))
            time_step = distance / average_velocity
            elapsed_time += time_step
            if equilibrium.pressure_pa <= nucleation_threshold:
                relaxation_fraction = 1.0 - exp(
                    -time_step / relaxation_settings.relaxation_time_s
                )
                actual_quality += relaxation_fraction * (
                    equilibrium.vapor_quality_mass - actual_quality
                )
                actual_quality = min(
                    equilibrium.vapor_quality_mass,
                    max(0.0, actual_quality),
                )

        if (
            actual_onset_pressure is None
            and actual_quality >= relaxation_settings.quality_onset_threshold
        ):
            actual_onset_pressure = equilibrium.pressure_pa
            actual_onset_fraction = float(position_m) / geometry.orifice_length_m

        actual_density = _relaxed_mixture_density(
            equilibrium.pressure_pa,
            actual_quality,
            tier0.nozzle_inlet_state.density_kg_m3,
            properties,
        )
        unconstrained_flux = actual_density * equilibrium.velocity_m_s
        local_spi_flux = sqrt(
            max(
                0.0,
                2.0
                * tier0.nozzle_inlet_state.density_kg_m3
                * (inlet_pressure - equilibrium.pressure_pa),
            )
        )
        bounded_flux = min(
            local_spi_flux,
            max(equilibrium.mass_flux_kg_m2_s, unconstrained_flux),
        )
        path.append(
            RelaxationPathPoint(
                position_m=float(position_m),
                pressure_pa=equilibrium.pressure_pa,
                velocity_m_s=equilibrium.velocity_m_s,
                equilibrium_vapor_quality_mass=equilibrium.vapor_quality_mass,
                actual_vapor_quality_mass=actual_quality,
                actual_density_kg_m3=actual_density,
                mass_flux_kg_m2_s=bounded_flux,
                elapsed_time_s=elapsed_time,
            )
        )
        previous_position = float(position_m)
        previous_velocity = equilibrium.velocity_m_s

    fluxes = np.asarray([point.mass_flux_kg_m2_s for point in path])
    critical_index = int(np.argmax(fluxes))
    critical = path[critical_index]
    pressure_tolerance = equilibrium_settings.choking_pressure_tolerance * inlet_pressure
    is_choked = critical.pressure_pa > outlet_pressure + pressure_tolerance
    operating_flux = critical.mass_flux_kg_m2_s if is_choked else path[-1].mass_flux_kg_m2_s
    required_area = operating_point.lng_mass_flow_kg_s / (
        geometry.discharge_coefficient * operating_flux
    )
    required_diameter = sqrt(
        4.0 * required_area / (pi * geometry.number_of_orifices)
    )
    geometric_area = geometry.geometric_area_m2
    predicted_mass_flow = (
        geometry.discharge_coefficient * geometric_area * operating_flux
        if geometric_area is not None
        else None
    )

    if tier0.flash_location == FlashLocation.UPSTREAM:
        actual_location = FlashLocation.UPSTREAM
    elif actual_onset_fraction is None:
        actual_location = (
            FlashLocation.EXTERNAL
            if tier0.equilibrium_flash_fraction_at_p3 > 0.0
            else FlashLocation.NONE
        )
    elif actual_onset_fraction >= 0.98:
        actual_location = FlashLocation.EXIT
    else:
        actual_location = FlashLocation.INTERNAL

    warnings = list(tier1.warnings)
    warnings.append(
        ModelWarning(
            code="RELAXATION_PARAMETERS_REQUIRE_CALIBRATION",
            severity=WarningSeverity.WARNING,
            message=(
                "Tier 2 nucleation delay and relaxation time are empirical. Propagate their "
                "uncertainty or calibrate them to a geometrically similar cryogenic injector."
            ),
        )
    )
    if actual_location == FlashLocation.EXTERNAL:
        warnings.append(
            ModelWarning(
                code="METASTABLE_LIQUID_AT_EXIT",
                severity=WarningSeverity.WARNING,
                message=(
                    "The finite-rate model retains liquid through the nozzle while an "
                    "equilibrium flash is available after discharge."
                ),
            )
        )

    return Tier2RelaxationFlow(
        tier1=tier1,
        settings=relaxation_settings,
        nucleation_threshold_pressure_pa=nucleation_threshold,
        actual_flash_location=actual_location,
        actual_flash_onset_pressure_pa=actual_onset_pressure,
        actual_flash_onset_fraction_of_length=actual_onset_fraction,
        residence_time_s=elapsed_time,
        equilibrium_exit_vapor_quality_mass=path[-1].equilibrium_vapor_quality_mass,
        actual_exit_vapor_quality_mass=path[-1].actual_vapor_quality_mass,
        operating_mass_flux_kg_m2_s=operating_flux,
        critical_mass_flux_kg_m2_s=critical.mass_flux_kg_m2_s,
        critical_pressure_pa=critical.pressure_pa,
        is_choked=is_choked,
        required_geometric_area_m2=required_area,
        required_orifice_diameter_m=required_diameter,
        predicted_mass_flow_kg_s=predicted_mass_flow,
        path=tuple(path),
        warnings=tuple(warnings),
    )


def _relaxed_mixture_density(
    pressure_pa: float,
    vapor_quality_mass: float,
    compressed_liquid_density_kg_m3: float,
    properties: CoolPropLNGProvider,
) -> float:
    if vapor_quality_mass <= 0.0:
        return compressed_liquid_density_kg_m3
    saturated_liquid = properties.state_pq(pressure_pa, 0.0)
    saturated_vapor = properties.state_pq(pressure_pa, 1.0)
    return 1.0 / (
        (1.0 - vapor_quality_mass) / saturated_liquid.density_kg_m3
        + vapor_quality_mass / saturated_vapor.density_kg_m3
    )
