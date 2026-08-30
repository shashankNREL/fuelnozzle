"""Reactor specifications for the combustor network.

A reactor here is a zone of the combustor treated as internally uniform. Which idealized
type it is says what assumption is being made about mixing inside it:

- a **perfectly stirred reactor** assumes anything entering is instantly mixed
  throughout, which suits recirculation zones and flame zones where turbulence is
  violent;
- a **plug flow reactor** assumes the flow marches along without mixing back, which
  suits the post-flame and dilution length where the axial velocity dominates;
- an **evaporator** and a **mixer** are stirred reactors that additionally host droplets.

Plug flow is represented as a chain of stirred reactors rather than by Cantera's
``FlowReactor``, so that it participates in the same simultaneous network solve as
everything else. The number of links in the chain is a numerical parameter, and its
convergence is a verification task rather than an assumption.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, model_validator


class ReactorKind(StrEnum):
    """What mixing assumption a zone is being given."""

    PSR = "perfectly_stirred"
    PFR = "plug_flow"
    EVAPORATOR = "evaporator"
    MIXER = "mixer"

    @property
    def hosts_droplets(self) -> bool:
        """Whether liquid fuel may be present in this zone."""
        return self in (ReactorKind.EVAPORATOR, ReactorKind.MIXER)


class ReactorSpec(BaseModel):
    """One zone of the combustor.

    ``volume_m3`` sets the residence time through ``tau = rho * V / mdot``. Because a
    reactor stands for a coarse region with real internal variation, that residence time
    is an effective zone-scale quantity, not the transit time of any particular parcel.
    """

    model_config = ConfigDict(frozen=True)

    name: str = Field(min_length=1)
    kind: ReactorKind
    volume_m3: float = Field(gt=0.0)

    #: Number of stirred reactors used to represent a plug flow zone.
    plug_flow_segments: int = Field(default=8, ge=1, le=200)

    #: Heat removed from this zone, for liner loss. Positive means heat leaves the gas.
    heat_loss_w: float = Field(default=0.0, ge=0.0)
    heat_loss_basis: str | None = Field(default=None, min_length=1)

    #: Average distance a droplet travels inside this zone. Droplets and gas do not
    #: travel together, so this is not derivable from the reactor's residence time.
    spray_path_length_m: float | None = Field(default=None, gt=0.0)

    @model_validator(mode="after")
    def validate_kind_consistency(self) -> ReactorSpec:
        if self.spray_path_length_m is not None and not self.kind.hosts_droplets:
            raise ValueError(
                f"Reactor {self.name!r} is a {self.kind} but was given a spray path "
                "length. Only evaporator and mixer zones host droplets."
            )
        if self.kind is not ReactorKind.PFR and self.plug_flow_segments != 8:
            raise ValueError(
                f"Reactor {self.name!r} is a {self.kind}; plug_flow_segments applies "
                "only to plug flow zones."
            )
        if self.heat_loss_w > 0.0 and self.heat_loss_basis is None:
            raise ValueError(
                f"Reactor {self.name!r} has a prescribed heat loss without a calibration "
                "or physical-model identifier"
            )
        return self

    @property
    def segment_volume_m3(self) -> float:
        """Volume of one link in a plug flow chain."""
        if self.kind is not ReactorKind.PFR:
            return self.volume_m3
        return self.volume_m3 / self.plug_flow_segments


class InletSpec(BaseModel):
    """An external stream entering the network.

    Air stations and fuel vapor both arrive this way. Fuel vapor released by evaporating
    droplets enters at the *droplet* temperature, not the gas temperature, which is what
    keeps the latent heat from being counted twice.
    """

    model_config = ConfigDict(frozen=True)

    name: str = Field(min_length=1)
    target_reactor: str = Field(min_length=1)
    mass_flow_kg_s: float = Field(ge=0.0)
    temperature_k: float = Field(gt=0.0)
    mole_fractions: dict[str, float] = Field(min_length=1)
    at_reactor_exit: bool = False

    @model_validator(mode="after")
    def validate_composition(self) -> InletSpec:
        if any(value < 0.0 for value in self.mole_fractions.values()):
            raise ValueError(f"Inlet {self.name!r} has a negative mole fraction")
        if sum(self.mole_fractions.values()) <= 0.0:
            raise ValueError(f"Inlet {self.name!r} has an empty composition")
        return self


class OutletSpec(BaseModel):
    """Where exhaust leaves the network."""

    model_config = ConfigDict(frozen=True)

    source_reactor: str = Field(min_length=1)
    mass_flow_kg_s: float = Field(gt=0.0)
