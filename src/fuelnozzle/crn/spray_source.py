"""Bridge from the nozzle solvers to the droplet classes the network integrates.

This is where the existing package meets the new one. The Jet-A pressure-swirl solver
and the LNG flashing tiers already determine how the fuel leaves the injector; this
module turns that into an initial droplet population plus, for LNG, the fraction that
has already flashed to vapor.

Two decisions here differ from John et al. (2026) and are deliberate.

**Initial size comes from the nozzle model.** The paper starts every droplet at the
nozzle radius and lets aerodynamic breakup find the real size. This package already
solves the atomizer, so starting from its Sauter mean diameter uses information the
paper did not have. The paper's approach remains available for comparison.

**Breakup is skipped for flashing sprays.** The Taylor analogy describes a droplet torn
apart by aerodynamic forces. A flashing LNG droplet is burst from within by vapor
generation, which is a different mechanism that the analogy does not represent. In the
flashing regimes the size is taken from the Tier 3 result instead, and the decision is
recorded rather than left implicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import gamma, log, pi

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.jet_a import PressureSwirlResult
from fuelnozzle.models import COOLPROP_TO_CANTERA_SPECIES, ModelWarning, WarningSeverity
from fuelnozzle.spray import FlashSprayRegime, Tier3FlashSpray

#: Regimes in which breakup is driven by internal vaporization, not aerodynamics.
FLASHING_REGIMES = frozenset(
    {
        FlashSprayRegime.FULLY_FLASHING,
        FlashSprayRegime.TRANSITIONAL_FLASH,
        FlashSprayRegime.EXTERNAL_FLASH,
        FlashSprayRegime.UPSTREAM_TWO_PHASE,
    }
)

#: Below this spread parameter the Rosin-Rammler distribution has no finite Sauter
#: mean diameter, because Gamma(1 - 1/n) diverges as n approaches 1.
_ROSIN_RAMMLER_MIN_SPREAD = 1.2


class InitialSizePolicy(StrEnum):
    """Where the initial droplet radius comes from."""

    #: Sauter mean radius from the nozzle model. Requires a spray calibration.
    NOZZLE_SMD = "nozzle_smd"
    #: The paper's approach: start at the orifice radius and let breakup act.
    NOZZLE_RADIUS_TAB = "nozzle_radius_tab"
    #: Explicit user-supplied size.
    USER = "user"


@dataclass(frozen=True)
class DropletClass:
    """One droplet population entering the network."""

    radius_m: float
    temperature_k: float
    velocity_m_s: float
    mass_flow_kg_s: float
    number_rate_per_s: float
    origin: str

    @property
    def diameter_m(self) -> float:
        return 2.0 * self.radius_m


@dataclass(frozen=True)
class SprayBoundary:
    """Everything the network needs to know about the fuel leaving the injector."""

    fuel: FuelKind
    total_fuel_mass_flow_kg_s: float
    vapor_mass_flow_kg_s: float
    vapor_temperature_k: float
    droplet_classes: tuple[DropletClass, ...]
    injection_velocity_m_s: float
    cone_angle_deg: float | None
    apply_aerodynamic_breakup: bool
    size_policy: InitialSizePolicy
    calibration_id: str | None
    warnings: tuple[ModelWarning, ...]
    vapor_mole_fractions: dict[str, float] | None = None
    liquid_mole_fractions: dict[str, float] | None = None

    @property
    def vapor_mass_fraction(self) -> float:
        if self.total_fuel_mass_flow_kg_s <= 0.0:
            return 0.0
        return self.vapor_mass_flow_kg_s / self.total_fuel_mass_flow_kg_s

    @property
    def liquid_mass_flow_kg_s(self) -> float:
        return sum(cls.mass_flow_kg_s for cls in self.droplet_classes)


def number_rate(mass_flow_kg_s: float, radius_m: float, liquid_density_kg_m3: float) -> float:
    """Droplets injected per second, John et al. Eq. (6)."""
    if radius_m <= 0.0 or liquid_density_kg_m3 <= 0.0:
        raise ValueError("Droplet number rate requires a positive radius and density")
    droplet_mass = liquid_density_kg_m3 * (4.0 / 3.0) * pi * radius_m**3
    return mass_flow_kg_s / droplet_mass


def rosin_rammler_classes(
    sauter_mean_radius_m: float,
    mass_flow_kg_s: float,
    liquid_density_kg_m3: float,
    temperature_k: float,
    velocity_m_s: float,
    spread_parameter: float = 2.5,
    class_count: int = 5,
    origin: str = "rosin_rammler",
) -> tuple[DropletClass, ...]:
    """Split a spray into size classes following a Rosin-Rammler distribution.

    A single representative droplet understates the spread of a real spray: the largest
    droplets survive far longer than the mean and are what actually reach the flame.
    Classes carry equal mass but different sizes, so the surviving tail is represented.

    The paper lists a size distribution as future work; it is available here from the
    start, with a single class remaining the default for speed.
    """
    if class_count < 1:
        raise ValueError("At least one droplet class is required")
    if spread_parameter < _ROSIN_RAMMLER_MIN_SPREAD:
        raise ValueError(
            f"Rosin-Rammler spread parameter must be at least {_ROSIN_RAMMLER_MIN_SPREAD}; "
            "smaller values give a distribution with no finite Sauter mean diameter."
        )
    if class_count == 1:
        return (
            DropletClass(
                radius_m=sauter_mean_radius_m,
                temperature_k=temperature_k,
                velocity_m_s=velocity_m_s,
                mass_flow_kg_s=mass_flow_kg_s,
                number_rate_per_s=number_rate(
                    mass_flow_kg_s, sauter_mean_radius_m, liquid_density_kg_m3
                ),
                origin=origin,
            ),
        )

    # For a Rosin-Rammler distribution, D32 = X / Gamma(1 - 1/n), where X is the
    # characteristic diameter. Invert it so the classes reproduce the requested D32.
    characteristic_diameter = (
        2.0 * sauter_mean_radius_m * gamma(1.0 - 1.0 / spread_parameter)
    )

    # Equal-mass classes: sample the mass distribution at the midpoint of each interval.
    classes: list[DropletClass] = []
    class_mass_flow = mass_flow_kg_s / class_count
    for index in range(class_count):
        cumulative = (index + 0.5) / class_count
        # Inverting the cumulative mass distribution 1 - exp(-(D/X)^n).
        diameter = characteristic_diameter * (-log(1.0 - cumulative)) ** (
            1.0 / spread_parameter
        )
        radius = 0.5 * diameter
        classes.append(
            DropletClass(
                radius_m=radius,
                temperature_k=temperature_k,
                velocity_m_s=velocity_m_s,
                mass_flow_kg_s=class_mass_flow,
                number_rate_per_s=number_rate(
                    class_mass_flow, radius, liquid_density_kg_m3
                ),
                origin=f"{origin}_class_{index + 1}",
            )
        )
    return tuple(classes)


def jet_a_spray_boundary(
    result: PressureSwirlResult,
    fuel_mass_flow_kg_s: float,
    fuel_temperature_k: float,
    liquid_density_kg_m3: float,
    policy: InitialSizePolicy = InitialSizePolicy.NOZZLE_SMD,
    user_radius_m: float | None = None,
    class_count: int = 1,
    spread_parameter: float = 2.5,
) -> SprayBoundary:
    """Build the droplet population leaving a Jet-A pressure-swirl atomizer."""
    warnings = list(result.warnings)
    velocity = result.axial_exit_velocity_m_s
    radius, resolved_policy, extra = _resolve_initial_radius(
        policy,
        smd_m=result.smd_estimate_m,
        orifice_radius_m=0.5 * result.modeled_exit_diameter_m,
        user_radius_m=user_radius_m,
        fuel_label="Jet-A",
    )
    warnings.extend(extra)

    classes = rosin_rammler_classes(
        radius,
        fuel_mass_flow_kg_s,
        liquid_density_kg_m3,
        fuel_temperature_k,
        velocity,
        spread_parameter=spread_parameter,
        class_count=class_count,
        origin="jet_a_pressure_swirl",
    )
    return SprayBoundary(
        fuel=FuelKind.JET_A,
        total_fuel_mass_flow_kg_s=fuel_mass_flow_kg_s,
        vapor_mass_flow_kg_s=0.0,
        vapor_temperature_k=fuel_temperature_k,
        droplet_classes=classes,
        injection_velocity_m_s=velocity,
        cone_angle_deg=result.full_cone_angle_deg,
        apply_aerodynamic_breakup=resolved_policy is InitialSizePolicy.NOZZLE_RADIUS_TAB,
        size_policy=resolved_policy,
        calibration_id=None,
        warnings=tuple(warnings),
    )


def lng_spray_boundary(
    result: Tier3FlashSpray,
    fuel_mass_flow_kg_s: float,
    liquid_density_kg_m3: float,
    policy: InitialSizePolicy = InitialSizePolicy.NOZZLE_SMD,
    user_radius_m: float | None = None,
    class_count: int = 1,
    spread_parameter: float = 2.5,
) -> SprayBoundary:
    """Build the spray boundary for a flashing LNG injector.

    Part of the fuel has already vaporized inside the nozzle. Only the remaining liquid
    becomes droplets; treating all of it as liquid would double-count the atomization
    work the flash has already done.
    """
    warnings = list(result.warnings)
    boundary = result.cfd_boundary
    def cantera_composition(
        composition: dict[str, float] | None,
    ) -> dict[str, float] | None:
        if composition is None:
            return None
        return {
            COOLPROP_TO_CANTERA_SPECIES[name]: fraction
            for name, fraction in composition.items()
            if fraction > 0.0 and name in COOLPROP_TO_CANTERA_SPECIES
        } or None

    vapor_composition = cantera_composition(boundary.vapor_mole_fractions)
    liquid_composition = cantera_composition(boundary.liquid_mole_fractions)
    vapor_fraction = min(1.0, max(0.0, result.actual_exit_vapor_quality_mass))
    vapor_flow = fuel_mass_flow_kg_s * vapor_fraction
    liquid_flow = fuel_mass_flow_kg_s - vapor_flow
    flashing = result.regime in FLASHING_REGIMES

    if flashing:
        warnings.append(
            ModelWarning(
                code="AERODYNAMIC_BREAKUP_SKIPPED",
                severity=WarningSeverity.INFO,
                message=(
                    f"LNG regime is {result.regime}, where breakup is driven by internal "
                    "vaporization rather than aerodynamic distortion. The Taylor analogy "
                    "does not apply and is not used; size comes from the Tier 3 result."
                ),
            )
        )

    if liquid_flow <= 0.0:
        warnings.append(
            ModelWarning(
                code="LNG_FULLY_VAPORIZED_AT_INJECTION",
                severity=WarningSeverity.INFO,
                message=(
                    "LNG leaves the injector fully vaporized; no droplets enter the "
                    "network and evaporation modeling is unnecessary."
                ),
            )
        )
        return SprayBoundary(
            fuel=FuelKind.LNG,
            total_fuel_mass_flow_kg_s=fuel_mass_flow_kg_s,
            vapor_mass_flow_kg_s=fuel_mass_flow_kg_s,
            vapor_temperature_k=boundary.temperature_k,
            droplet_classes=(),
            injection_velocity_m_s=boundary.velocity_m_s,
            cone_angle_deg=result.full_cone_angle_estimate_deg,
            apply_aerodynamic_breakup=False,
            size_policy=policy,
            calibration_id=result.calibration_id,
            warnings=tuple(warnings),
            vapor_mole_fractions=vapor_composition,
            liquid_mole_fractions=liquid_composition,
        )

    radius, resolved_policy, extra = _resolve_initial_radius(
        policy,
        smd_m=result.smd_estimate_m,
        orifice_radius_m=0.5 * boundary.effective_orifice_diameter_m,
        user_radius_m=user_radius_m,
        fuel_label="LNG",
    )
    warnings.extend(extra)

    classes = rosin_rammler_classes(
        radius,
        liquid_flow,
        liquid_density_kg_m3,
        boundary.temperature_k,
        boundary.velocity_m_s,
        spread_parameter=spread_parameter,
        class_count=class_count,
        origin="lng_flash_spray",
    )
    return SprayBoundary(
        fuel=FuelKind.LNG,
        total_fuel_mass_flow_kg_s=fuel_mass_flow_kg_s,
        vapor_mass_flow_kg_s=vapor_flow,
        vapor_temperature_k=boundary.temperature_k,
        droplet_classes=classes,
        injection_velocity_m_s=boundary.velocity_m_s,
        cone_angle_deg=result.full_cone_angle_estimate_deg,
        apply_aerodynamic_breakup=(
            resolved_policy is InitialSizePolicy.NOZZLE_RADIUS_TAB and not flashing
        ),
        size_policy=resolved_policy,
        calibration_id=result.calibration_id,
        warnings=tuple(warnings),
        vapor_mole_fractions=vapor_composition,
        liquid_mole_fractions=liquid_composition,
    )


def _resolve_initial_radius(
    policy: InitialSizePolicy,
    *,
    smd_m: float | None,
    orifice_radius_m: float,
    user_radius_m: float | None,
    fuel_label: str,
) -> tuple[float, InitialSizePolicy, list[ModelWarning]]:
    """Choose the initial droplet radius, falling back only with an explicit warning.

    The existing package refuses to report a Sauter mean diameter without a hardware
    calibration. That discipline is inherited here: rather than inventing a size, the
    policy degrades to the paper's orifice-radius approach and says so.
    """
    warnings: list[ModelWarning] = []

    if policy is InitialSizePolicy.USER:
        if user_radius_m is None or user_radius_m <= 0.0:
            raise ValueError("InitialSizePolicy.USER requires a positive user_radius_m")
        return user_radius_m, policy, warnings

    if policy is InitialSizePolicy.NOZZLE_SMD:
        if smd_m is not None and smd_m > 0.0:
            return 0.5 * smd_m, policy, warnings
        warnings.append(
            ModelWarning(
                code="SPRAY_SIZE_POLICY_DOWNGRADED",
                severity=WarningSeverity.WARNING,
                message=(
                    f"No calibrated {fuel_label} Sauter mean diameter is available, so the "
                    "initial droplet size falls back to the orifice radius with "
                    "aerodynamic breakup. Atomization-quality results must not be "
                    "reported from this fallback; supply a spray calibration."
                ),
            )
        )
        return orifice_radius_m, InitialSizePolicy.NOZZLE_RADIUS_TAB, warnings

    return orifice_radius_m, policy, warnings
