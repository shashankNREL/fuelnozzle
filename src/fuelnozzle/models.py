"""Shared inputs and outputs for the fuel-nozzle models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import isfinite

from pydantic import BaseModel, ConfigDict, Field, field_validator

COOLPROP_TO_CANTERA_SPECIES = {
    "Methane": "CH4",
    "Ethane": "C2H6",
    "Propane": "C3H8",
    "Nitrogen": "N2",
}


class WarningSeverity(StrEnum):
    """Severity attached to an engineering model warning."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True)
class ModelWarning:
    """A model-limit, consistency, or design warning."""

    code: str
    severity: WarningSeverity
    message: str


class LNGComposition(BaseModel):
    """LNG composition expressed as normalized component mole fractions."""

    model_config = ConfigDict(frozen=True)

    mole_fractions: dict[str, float] = Field(min_length=1)

    @field_validator("mole_fractions")
    @classmethod
    def validate_and_normalize(cls, values: dict[str, float]) -> dict[str, float]:
        if any(not isfinite(value) or value < 0.0 for value in values.values()):
            raise ValueError("LNG mole fractions must be finite and non-negative")
        total = sum(values.values())
        if total <= 0.0:
            raise ValueError("At least one LNG component must have a positive fraction")
        return {name: value / total for name, value in values.items() if value > 0.0}

    @classmethod
    def pure_methane(cls) -> LNGComposition:
        return cls(mole_fractions={"Methane": 1.0})

    @property
    def components(self) -> tuple[str, ...]:
        return tuple(self.mole_fractions)

    @property
    def fractions(self) -> tuple[float, ...]:
        return tuple(self.mole_fractions.values())

    @property
    def is_pure(self) -> bool:
        return len(self.mole_fractions) == 1

    def cantera_mole_fractions(self) -> dict[str, float]:
        """Map thermodynamic-fluid names onto kinetic-mechanism species names."""
        unmapped = [
            component
            for component in self.mole_fractions
            if component not in COOLPROP_TO_CANTERA_SPECIES
        ]
        if unmapped:
            raise ValueError(
                "No Cantera species mapping is declared for LNG components "
                f"{', '.join(unmapped)}"
            )
        return {
            COOLPROP_TO_CANTERA_SPECIES[component]: fraction
            for component, fraction in self.mole_fractions.items()
        }


@dataclass(frozen=True)
class ThermodynamicState:
    """Mass-basis thermodynamic state returned by a property provider."""

    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float
    phase: str
    vapor_quality_mass: float | None
    cp_j_kg_k: float | None
    viscosity_pa_s: float | None
    conductivity_w_m_k: float | None
    surface_tension_n_m: float | None
    transport_model: str
    liquid_mole_fractions: dict[str, float] | None = None
    vapor_mole_fractions: dict[str, float] | None = None
