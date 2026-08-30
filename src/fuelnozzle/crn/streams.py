"""Combustor air-side bookkeeping: how much air goes where, and at what state.

The liner air split is the shared-hardware compromise in a dual-fuel combustor. Jet-A
at LTO wants a low dome fraction if it runs rich-quench-lean; LNG at cruise wants a high
dome fraction if it runs lean premixed. Fixed hole areas cannot deliver both directly.

This module also carries the mechanism that partly reconciles them. Because the two fuel
circuits are separate hardware, air flows through *both* injector passages at all times,
including the one whose fuel is shut off. Sizing the two passages therefore shifts the
effective near-field air split between mission segments with completely fixed geometry.
See ``docs/CRN_PLAN.md`` Section 8.2.1.

That effect depends on how much of the idle passage's unfueled air reaches the near
field, which a reactor network cannot determine on its own. It is exposed here as an
explicit parameter, and results that depend on it must be reported with a sensitivity
across its range rather than at a single assumed value.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import cantera as ct
from pydantic import BaseModel, ConfigDict, Field, model_validator

from fuelnozzle.crn.chemistry import (
    DRY_AIR_MOLE_FRACTIONS,
    FuelKind,
    MechanismSpec,
    stoichiometric_air_fuel_ratio,
)
from fuelnozzle.models import ModelWarning, WarningSeverity

#: Fractions must sum to one within this tolerance before a split is accepted.
SPLIT_SUM_TOLERANCE = 1.0e-6


class CoolingAirDestination(StrEnum):
    """Where liner cooling film air rejoins the flow.

    The choice changes CO burnout and exit temperature uniformity, so it is stated
    explicitly rather than assumed.
    """

    EXIT = "exit"
    DILUTION = "dilution"
    PRIMARY = "primary"


class AirSplit(BaseModel):
    """Fractions of total combustor air delivered to each station.

    ``dome`` is further divided between the two injector air passages by
    ``jet_a_passage_share``, which is the design lever of Section 8.2.1.
    """

    model_config = ConfigDict(frozen=True)

    dome: float = Field(ge=0.0, le=1.0)
    primary: float = Field(ge=0.0, le=1.0)
    quench: float = Field(ge=0.0, le=1.0)
    dilution: float = Field(ge=0.0, le=1.0)
    cooling: float = Field(ge=0.0, le=1.0)
    cooling_destination: CoolingAirDestination = CoolingAirDestination.DILUTION

    jet_a_passage_share: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Share of dome air carried by the Jet-A injector passage; the LNG "
        "passage carries the remainder.",
    )
    idle_passage_mixing_fraction: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Fraction of the inactive circuit's unfueled passage air that reaches "
        "the near field rather than bypassing it. A reactor network cannot determine this; "
        "report results with a sensitivity across its range.",
    )

    @model_validator(mode="after")
    def validate_sum(self) -> AirSplit:
        total = self.dome + self.primary + self.quench + self.dilution + self.cooling
        if abs(total - 1.0) > SPLIT_SUM_TOLERANCE:
            raise ValueError(
                f"Air split fractions must sum to 1.0; they sum to {total:.6f}. "
                "The tool does not renormalize them silently."
            )
        return self

    @property
    def jet_a_passage(self) -> float:
        """Absolute fraction of total air through the Jet-A injector passage."""
        return self.dome * self.jet_a_passage_share

    @property
    def lng_passage(self) -> float:
        """Absolute fraction of total air through the LNG injector passage."""
        return self.dome * (1.0 - self.jet_a_passage_share)

    def active_passage(self, fuel: FuelKind) -> float:
        """Dome air through the passage of the fuel that is currently burning."""
        return self.jet_a_passage if fuel is FuelKind.JET_A else self.lng_passage

    def idle_passage(self, fuel: FuelKind) -> float:
        """Dome air through the passage whose fuel is shut off."""
        return self.lng_passage if fuel is FuelKind.JET_A else self.jet_a_passage

    def near_field_air_fraction(self, fuel: FuelKind) -> float:
        """Air fraction setting the near-field equivalence ratio for the active fuel.

        The active passage always contributes. The idle passage contributes only the
        portion that actually reaches the near field.
        """
        return self.active_passage(fuel) + (
            self.idle_passage(fuel) * self.idle_passage_mixing_fraction
        )


@dataclass(frozen=True)
class AirState:
    """Thermodynamic state of the combustor inlet air."""

    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    mean_molecular_weight_kg_kmol: float
    mole_fractions: dict[str, float]


@dataclass(frozen=True)
class AirStreams:
    """Resolved air mass flows for one operating point."""

    state: AirState
    total_mass_flow_kg_s: float
    dome_kg_s: float
    primary_kg_s: float
    quench_kg_s: float
    dilution_kg_s: float
    cooling_kg_s: float
    active_passage_kg_s: float
    idle_passage_kg_s: float
    near_field_air_kg_s: float
    overall_equivalence_ratio: float
    near_field_equivalence_ratio: float
    warnings: tuple[ModelWarning, ...]


def air_state(
    solution: ct.Solution,
    temperature_k: float,
    pressure_pa: float,
    mole_fractions: dict[str, float] | None = None,
) -> AirState:
    """Evaluate the inlet air state on the active mechanism."""
    composition = mole_fractions or DRY_AIR_MOLE_FRACTIONS
    missing = [name for name in composition if name not in solution.species_names]
    if missing:
        raise ValueError(f"Oxidizer species {', '.join(missing)} absent from the mechanism")
    solution.TPX = temperature_k, pressure_pa, composition
    return AirState(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        density_kg_m3=float(solution.density_mass),
        mean_molecular_weight_kg_kmol=float(solution.mean_molecular_weight),
        mole_fractions=dict(composition),
    )


def resolve_air_streams(
    solution: ct.Solution,
    spec: MechanismSpec,
    split: AirSplit,
    *,
    fuel: FuelKind,
    fuel_mass_flow_kg_s: float,
    temperature_k: float,
    pressure_pa: float,
    total_air_mass_flow_kg_s: float | None = None,
    overall_equivalence_ratio: float | None = None,
    oxidizer_mole_fractions: dict[str, float] | None = None,
) -> AirStreams:
    """Distribute combustor air across stations for one operating point.

    Supply exactly one of ``total_air_mass_flow_kg_s`` or
    ``overall_equivalence_ratio``; the other is derived from the stoichiometric
    air-fuel ratio of the active fuel.
    """
    if fuel_mass_flow_kg_s <= 0.0:
        raise ValueError("Air resolution requires a positive fuel mass flow")
    if (total_air_mass_flow_kg_s is None) == (overall_equivalence_ratio is None):
        raise ValueError(
            "Supply exactly one of total_air_mass_flow_kg_s or overall_equivalence_ratio"
        )

    state = air_state(solution, temperature_k, pressure_pa, oxidizer_mole_fractions)
    afr_stoich = stoichiometric_air_fuel_ratio(solution, spec, oxidizer_mole_fractions)

    if total_air_mass_flow_kg_s is None:
        assert overall_equivalence_ratio is not None
        if overall_equivalence_ratio <= 0.0:
            raise ValueError("Overall equivalence ratio must be positive")
        total_air = fuel_mass_flow_kg_s * afr_stoich / overall_equivalence_ratio
        phi_overall = overall_equivalence_ratio
    else:
        if total_air_mass_flow_kg_s <= 0.0:
            raise ValueError("Total air mass flow must be positive")
        total_air = total_air_mass_flow_kg_s
        phi_overall = fuel_mass_flow_kg_s * afr_stoich / total_air

    near_field_air = total_air * split.near_field_air_fraction(fuel)
    phi_near_field = (
        fuel_mass_flow_kg_s * afr_stoich / near_field_air if near_field_air > 0.0 else float("inf")
    )

    warnings: list[ModelWarning] = []
    if split.dome > 0.0 and split.idle_passage(fuel) > 0.0:
        warnings.append(
            ModelWarning(
                code="IDLE_PASSAGE_AIR_ASSUMPTION",
                severity=WarningSeverity.INFO,
                message=(
                    "Near-field equivalence ratio assumes "
                    f"{split.idle_passage_mixing_fraction:.0%} of the idle "
                    f"{'LNG' if fuel is FuelKind.JET_A else 'Jet-A'} passage air reaches the "
                    "near field. A reactor network cannot determine this value; report "
                    "conclusions with a sensitivity across its range."
                ),
            )
        )
    if near_field_air <= 0.0:
        warnings.append(
            ModelWarning(
                code="NO_NEAR_FIELD_AIR",
                severity=WarningSeverity.ERROR,
                message="No air reaches the near field; the equivalence ratio is undefined.",
            )
        )

    return AirStreams(
        state=state,
        total_mass_flow_kg_s=total_air,
        dome_kg_s=total_air * split.dome,
        primary_kg_s=total_air * split.primary,
        quench_kg_s=total_air * split.quench,
        dilution_kg_s=total_air * split.dilution,
        cooling_kg_s=total_air * split.cooling,
        active_passage_kg_s=total_air * split.active_passage(fuel),
        idle_passage_kg_s=total_air * split.idle_passage(fuel),
        near_field_air_kg_s=near_field_air,
        overall_equivalence_ratio=phi_overall,
        near_field_equivalence_ratio=phi_near_field,
        warnings=tuple(warnings),
    )
