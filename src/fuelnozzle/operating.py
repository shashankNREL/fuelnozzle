"""Flight-envelope operating points and fuel pressure budgets."""

from __future__ import annotations

from dataclasses import dataclass

from pydantic import BaseModel, ConfigDict, Field

from fuelnozzle.models import ModelWarning, WarningSeverity


class OperatingPoint(BaseModel):
    """User-specified engine and fuel conditions for one flight stage.

    Mass flows are totals through the modeled injector circuit. Geometry
    objects specify how that flow is divided among individual orifices.
    """

    model_config = ConfigDict(frozen=True)

    name: str = Field(min_length=1)
    duration_s: float | None = Field(default=None, gt=0.0)
    flow_multiplier: float = Field(
        default=1.0,
        gt=0.0,
        description="Number of identical engine/nozzle circuits represented by the mass flows",
    )
    p3_pa: float = Field(
        gt=0.0,
        description="Combustor inlet/chamber pressure used by nozzle models",
    )
    t3_k: float = Field(gt=0.0, description="Combustor inlet air temperature")
    lng_mass_flow_kg_s: float = Field(ge=0.0)
    jet_a_mass_flow_kg_s: float = Field(ge=0.0)
    lng_pump_outlet_pressure_pa: float = Field(gt=0.0)
    lng_pump_outlet_temperature_k: float = Field(gt=0.0)
    lng_nozzle_pressure_drop_pa: float = Field(gt=0.0)
    jet_a_pump_outlet_pressure_pa: float = Field(gt=0.0)
    jet_a_nozzle_inlet_temperature_k: float = Field(gt=0.0)
    jet_a_nozzle_pressure_drop_pa: float = Field(gt=0.0)

    # Combustor-side inputs, used only by the reactor-network models. All optional, so
    # existing nozzle-only studies are unaffected.
    combustor_air_mass_flow_kg_s: float | None = Field(default=None, gt=0.0)
    overall_equivalence_ratio: float | None = Field(default=None, gt=0.0)
    liner_pressure_loss_fraction: float | None = Field(default=None, ge=0.0, lt=1.0)
    rated_thrust_kn: float | None = Field(default=None, gt=0.0)
    thrust_fraction: float | None = Field(default=None, gt=0.0, le=1.0)
    nozzle_wall_temperature_k: float | None = Field(
        default=None,
        gt=0.0,
        description="Wall temperature seen by the fuel circuit that is shut off, used "
        "for the idle-circuit coking and vapour-lock screen.",
    )


    @property
    def active_fuel(self) -> str | None:
        """Which fuel this point burns, or ``None`` if both or neither flow.

        One fuel is active per operating point by design. A point with both flowing is
        not rejected here, because the nozzle models predate that rule and may legitimately
        size both circuits, but the reactor-network models require a single answer.
        """
        lng = self.lng_mass_flow_kg_s > 0.0
        jet_a = self.jet_a_mass_flow_kg_s > 0.0
        if lng and not jet_a:
            return "lng"
        if jet_a and not lng:
            return "jet_a"
        return None


@dataclass(frozen=True)
class PressureBudget:
    """Pressure allocation between a pump, feed system, nozzle, and combustor."""

    pump_outlet_pressure_pa: float
    nozzle_inlet_pressure_pa: float
    chamber_pressure_pa: float
    nozzle_pressure_drop_pa: float
    available_feed_pressure_drop_pa: float
    warnings: tuple[ModelWarning, ...]


@dataclass(frozen=True)
class ResolvedPressureStations:
    """Pressures used by the compressor, liner, combustor, and fuel nozzles."""

    compressor_discharge_pa: float
    dome_pa: float
    combustor_exit_pa: float
    lng_nozzle_inlet_pa: float
    jet_a_nozzle_inlet_pa: float
    liner_pressure_loss_pa: float


def resolve_pressure_stations(point: OperatingPoint) -> ResolvedPressureStations:
    """Resolve one explicit pressure path without changing legacy nozzle semantics."""
    loss_fraction = point.liner_pressure_loss_fraction or 0.0
    liner_loss = point.p3_pa * loss_fraction
    combustor_exit = point.p3_pa - liner_loss
    if combustor_exit <= 0.0:
        raise ValueError("Liner pressure loss leaves a non-positive combustor pressure")
    active_pump = (
        point.jet_a_pump_outlet_pressure_pa
        if point.active_fuel == "jet_a"
        else point.lng_pump_outlet_pressure_pa
    )
    active_nozzle_inlet = point.p3_pa + (
        point.jet_a_nozzle_pressure_drop_pa
        if point.active_fuel == "jet_a"
        else point.lng_nozzle_pressure_drop_pa
    )
    if point.active_fuel is not None and active_pump < active_nozzle_inlet:
        raise ValueError("Active-fuel pump pressure cannot supply the requested nozzle pressure")
    return ResolvedPressureStations(
        compressor_discharge_pa=point.p3_pa,
        dome_pa=point.p3_pa,
        combustor_exit_pa=combustor_exit,
        lng_nozzle_inlet_pa=point.p3_pa + point.lng_nozzle_pressure_drop_pa,
        jet_a_nozzle_inlet_pa=point.p3_pa + point.jet_a_nozzle_pressure_drop_pa,
        liner_pressure_loss_pa=liner_loss,
    )


def fuel_pressure_budget(
    pump_outlet_pressure_pa: float,
    nozzle_pressure_drop_pa: float,
    chamber_pressure_pa: float,
) -> PressureBudget:
    """Build a physically closed pressure budget without modifying user inputs."""

    nozzle_inlet_pressure_pa = chamber_pressure_pa + nozzle_pressure_drop_pa
    available_feed_pressure_drop_pa = pump_outlet_pressure_pa - nozzle_inlet_pressure_pa
    warnings: list[ModelWarning] = []
    if available_feed_pressure_drop_pa < 0.0:
        warnings.append(
            ModelWarning(
                code="PRESSURE_BUDGET_DEFICIT",
                severity=WarningSeverity.ERROR,
                message=(
                    "Pump outlet pressure is below P3 plus the specified nozzle pressure "
                    "drop; the requested operating point is hydraulically infeasible."
                ),
            )
        )
    elif available_feed_pressure_drop_pa < 0.05 * nozzle_pressure_drop_pa:
        warnings.append(
            ModelWarning(
                code="LOW_FEED_PRESSURE_MARGIN",
                severity=WarningSeverity.WARNING,
                message="Less than 5% of nozzle pressure drop remains for feed-system losses.",
            )
        )

    return PressureBudget(
        pump_outlet_pressure_pa=pump_outlet_pressure_pa,
        nozzle_inlet_pressure_pa=nozzle_inlet_pressure_pa,
        chamber_pressure_pa=chamber_pressure_pa,
        nozzle_pressure_drop_pa=nozzle_pressure_drop_pa,
        available_feed_pressure_drop_pa=available_feed_pressure_drop_pa,
        warnings=tuple(warnings),
    )
