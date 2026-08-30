"""Design variables for the dual-fuel combustor.

The variables fall into three classes, and keeping them apart is the most important
modelling decision in the whole study.

**Class A, the shared combustor.** One liner serves both fuels, so its air split and zone
volumes are chosen once and both fuels live with the consequences. This is where the
compromise lives.

**Class A2, per-fuel injector hardware.** The two circuits are separate hardware. Each is
optimized for its own fuel with no cross-fuel compromise, and the only things tying them
together are that their air passages share the dome air and that both must fit the dome.

**Class B, schedulable.** Values that may differ between mission segments because they are
set by the fuel system rather than by geometry. LNG temperature is the important one: the
feed line can deliver whatever the heat budget allows, so it is a free variable at each
segment rather than a property of the hardware.

Confusing Class A with Class A2 would invent a compromise that does not exist; confusing
Class A with Class B would hide one that does.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, model_validator

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.streams import AirSplit
from fuelnozzle.models import ModelWarning, WarningSeverity
from fuelnozzle.operating import OperatingPoint, ResolvedPressureStations


class VariableClass(StrEnum):
    """Which of the three classes a variable belongs to."""

    SHARED_COMBUSTOR = "shared_combustor"
    PER_FUEL_INJECTOR = "per_fuel_injector"
    SCHEDULABLE = "schedulable"

    @property
    def is_compromised(self) -> bool:
        """Whether both fuels must accept one value.

        Only the shared liner forces a compromise. Separate injector hardware does not,
        and schedulable values can differ by segment.
        """
        return self is VariableClass.SHARED_COMBUSTOR


@dataclass(frozen=True)
class DesignBound:
    """One variable's range and class."""

    name: str
    low: float
    high: float
    variable_class: VariableClass

    def clip(self, value: float) -> float:
        return min(self.high, max(self.low, value))

    def denormalize(self, unit_value: float) -> float:
        """Map a value in [0, 1] onto the range, for samplers that work in unit space."""
        return self.low + (self.high - self.low) * min(1.0, max(0.0, unit_value))

    def normalize(self, value: float) -> float:
        if self.high == self.low:
            return 0.0
        return (value - self.low) / (self.high - self.low)


class DesignVector(BaseModel):
    """One candidate design.

    Air fractions are stored so that the liner is always internally consistent: the dome,
    quench, and cooling fractions are chosen and dilution takes the remainder, which
    cannot then be forgotten or double-counted.
    """

    model_config = ConfigDict(frozen=True)

    # --- Class A: shared combustor ---
    dome_air_fraction: float = Field(ge=0.05, le=0.80)
    quench_air_fraction: float = Field(ge=0.0, le=0.60)
    primary_air_fraction: float = Field(default=0.05, ge=0.0, le=0.30)
    cooling_air_fraction: float = Field(default=0.05, ge=0.0, le=0.30)
    quench_volume_m3: float = Field(default=2.0e-4, gt=0.0)
    quench_stages: int = Field(default=12, ge=1, le=40)
    flame_volume_m3: float = Field(default=1.5e-3, gt=0.0)
    post_volume_m3: float = Field(default=4.0e-3, gt=0.0)

    # --- Class A2: per-fuel injector hardware ---
    jet_a_passage_share: float = Field(ge=0.0, le=1.0)
    jet_a_premix_residence_s: float = Field(default=0.0, ge=0.0)
    lng_premix_residence_s: float = Field(default=0.0, ge=0.0)
    dome_packaging_budget: float = Field(
        default=1.0,
        gt=0.0,
        description="Sum of the two injectors' packaging demand relative to the dome "
        "envelope. Above 1.0 they do not both fit.",
    )

    # --- Class B: schedulable per segment ---
    lng_fuel_temperature_k: float = Field(default=175.0, gt=0.0)
    jet_a_fuel_temperature_k: float = Field(default=440.0, gt=0.0)

    # --- assumption the results must be reported across ---
    idle_passage_mixing_fraction: float = Field(default=0.0, ge=0.0, le=1.0)

    @model_validator(mode="after")
    def validate_air(self) -> DesignVector:
        fixed = (
            self.dome_air_fraction
            + self.quench_air_fraction
            + self.primary_air_fraction
            + self.cooling_air_fraction
        )
        if fixed > 1.0:
            raise ValueError(
                f"Dome, quench, primary, and cooling air sum to {fixed:.3f}, leaving no "
                "room for dilution. The liner cannot deliver this split."
            )
        return self

    @property
    def dilution_air_fraction(self) -> float:
        """Whatever the other stations do not take."""
        return 1.0 - (
            self.dome_air_fraction
            + self.quench_air_fraction
            + self.primary_air_fraction
            + self.cooling_air_fraction
        )

    def air_split(self, fuel: FuelKind) -> AirSplit:
        return AirSplit(
            dome=self.dome_air_fraction,
            primary=self.primary_air_fraction,
            quench=self.quench_air_fraction,
            dilution=self.dilution_air_fraction,
            cooling=self.cooling_air_fraction,
            jet_a_passage_share=self.jet_a_passage_share,
            idle_passage_mixing_fraction=self.idle_passage_mixing_fraction,
        )

    def premix_residence_s(self, fuel: FuelKind) -> float:
        return (
            self.jet_a_premix_residence_s
            if fuel is FuelKind.JET_A
            else self.lng_premix_residence_s
        )

    def fuel_temperature_k(self, fuel: FuelKind) -> float:
        return (
            self.jet_a_fuel_temperature_k
            if fuel is FuelKind.JET_A
            else self.lng_fuel_temperature_k
        )

    def with_values(self, **updates) -> DesignVector:
        """A copy with some variables changed, for sweeps and perturbations.

        Rebuilt rather than copied. ``model_copy(update=...)`` does **not** re-run
        validators, so it would happily produce a design whose air fractions sum above
        one, and the failure would surface much later as a negative dilution fraction
        deep inside the network builder.
        """
        return DesignVector(**{**self.model_dump(), **updates})

    @property
    def packaging_warnings(self) -> tuple[ModelWarning, ...]:
        if self.dome_packaging_budget > 1.0:
            return (
                ModelWarning(
                    code="DOME_PACKAGING_EXCEEDED",
                    severity=WarningSeverity.ERROR,
                    message=(
                        f"The two injectors demand {self.dome_packaging_budget:.2f} of the "
                        "dome envelope. Separate fuel circuits still have to fit in one "
                        "dome; this design is not buildable."
                    ),
                ),
            )
        return ()


#: The design space actually swept. Ranges are deliberately wide, because the point of
#: the sensitivity stage is to find which of them matter, not to pre-judge it.
DEFAULT_BOUNDS: tuple[DesignBound, ...] = (
    DesignBound("dome_air_fraction", 0.20, 0.70, VariableClass.SHARED_COMBUSTOR),
    DesignBound("quench_air_fraction", 0.05, 0.40, VariableClass.SHARED_COMBUSTOR),
    DesignBound("quench_volume_m3", 1.0e-4, 1.0e-3, VariableClass.SHARED_COMBUSTOR),
    DesignBound("flame_volume_m3", 0.8e-3, 3.0e-3, VariableClass.SHARED_COMBUSTOR),
    DesignBound("jet_a_passage_share", 0.10, 0.90, VariableClass.PER_FUEL_INJECTOR),
    DesignBound("lng_premix_residence_s", 0.0, 5.0e-3, VariableClass.PER_FUEL_INJECTOR),
    DesignBound("jet_a_premix_residence_s", 0.0, 2.0e-3, VariableClass.PER_FUEL_INJECTOR),
    DesignBound("lng_fuel_temperature_k", 160.0, 230.0, VariableClass.SCHEDULABLE),
    DesignBound("jet_a_fuel_temperature_k", 350.0, 460.0, VariableClass.SCHEDULABLE),
)


def bounds_by_class(
    variable_class: VariableClass, bounds: tuple[DesignBound, ...] = DEFAULT_BOUNDS
) -> tuple[DesignBound, ...]:
    return tuple(bound for bound in bounds if bound.variable_class is variable_class)


def baseline_design() -> DesignVector:
    """A plausible starting point, not an optimum."""
    return DesignVector(
        dome_air_fraction=0.38,
        quench_air_fraction=0.30,
        jet_a_passage_share=0.50,
    )


def from_unit_cube(
    values: dict[str, float],
    base: DesignVector | None = None,
    bounds: tuple[DesignBound, ...] = DEFAULT_BOUNDS,
) -> DesignVector:
    """Build a design from values in [0, 1], the form samplers produce.

    Air fractions are repaired rather than rejected when a sample would leave no room
    for dilution: a sampler that trips this on most draws would waste the budget, and
    clipping is a defensible projection back into the feasible set. The repair is
    reported so it can be seen in a sweep.
    """
    base = base or baseline_design()
    lookup = {bound.name: bound for bound in bounds}
    updates = {
        name: lookup[name].denormalize(value)
        for name, value in values.items()
        if name in lookup
    }
    candidate = {**base.model_dump(), **updates}

    fixed = (
        candidate["dome_air_fraction"]
        + candidate["quench_air_fraction"]
        + candidate["primary_air_fraction"]
        + candidate["cooling_air_fraction"]
    )
    if fixed >= 0.98:
        # Leave at least 2% for dilution by scaling the two swept fractions down.
        scale = (0.98 - candidate["primary_air_fraction"] - candidate["cooling_air_fraction"]) / (
            candidate["dome_air_fraction"] + candidate["quench_air_fraction"]
        )
        candidate["dome_air_fraction"] *= scale
        candidate["quench_air_fraction"] *= scale

    return DesignVector(**candidate)


def perturb(
    design: DesignVector, name: str, delta: float,
    bounds: tuple[DesignBound, ...] = DEFAULT_BOUNDS,
) -> DesignVector:
    """Move one variable by a fraction of its range, for one-at-a-time sensitivity."""
    lookup = {bound.name: bound for bound in bounds}
    if name not in lookup:
        raise KeyError(f"{name!r} is not a swept design variable")
    bound = lookup[name]
    current = getattr(design, name)
    moved = bound.clip(current + delta * (bound.high - bound.low))
    return design.with_values(**{name: moved})


@dataclass(frozen=True)
class MissionPoint:
    """One steady operating point, burning exactly one fuel."""

    name: str
    fuel: FuelKind
    fuel_mass_flow_kg_s: float
    air_mass_flow_kg_s: float
    air_temperature_k: float
    pressure_pa: float
    duration_s: float = 0.0
    thrust_fraction: float = 1.0
    nozzle_wall_temperature_k: float | None = None
    operating_point: OperatingPoint | None = None
    pressure_stations: ResolvedPressureStations | None = None

    def scaled(self, fuel_scale: float) -> MissionPoint:
        return replace(self, fuel_mass_flow_kg_s=self.fuel_mass_flow_kg_s * fuel_scale)
