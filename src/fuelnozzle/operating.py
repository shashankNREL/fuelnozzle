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


@dataclass(frozen=True)
class PressureBudget:
    """Pressure allocation between a pump, feed system, nozzle, and combustor."""

    pump_outlet_pressure_pa: float
    nozzle_inlet_pressure_pa: float
    chamber_pressure_pa: float
    nozzle_pressure_drop_pa: float
    available_feed_pressure_drop_pa: float
    warnings: tuple[ModelWarning, ...]


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
