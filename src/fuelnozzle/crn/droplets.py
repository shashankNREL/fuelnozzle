"""Droplet breakup, evaporation, and heating for the liquid-fueled reactor network.

Three processes act on a droplet between the injector exit and the flame: it is torn
into smaller droplets by aerodynamic forces, it is heated by the surrounding gas, and it
turns into vapor. Breakup happens roughly two orders of magnitude faster than
evaporation, so the two are treated in sequence rather than simultaneously, following
John et al. (2026).

Two departures from that paper are deliberate and are the reason this module exists in
its present form.

**A boiling branch.** The paper models Jet-A only, where evaporation is limited by how
fast vapor can diffuse away from the droplet surface. Flashing LNG is not diffusion
limited: the droplet is above its boiling point at chamber pressure, so vaporization is
limited instead by how fast heat arrives. The two branches meet exactly where the
surface vapor pressure reaches the ambient pressure, so the switch between them is
physically located, not tuned.

**Explicit energy accounting.** The gas loses the convective heat; the droplet gains it,
spends part on the phase change, and puts the rest into warming; and the vapor enters
the gas carrying the enthalpy it has *at the droplet temperature*, not at the gas
temperature. Written any other way, latent heat is silently counted twice. See
``docs/CRN_PLAN.md`` Section 3.6.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import log, pi, sqrt
from typing import Protocol

import numpy as np
from scipy.integrate import solve_ivp

# Taylor analogy breakup coefficients (O'Rourke and Amsden).
TAB_CK = 8.0
TAB_CF = 1.0 / 3.0
TAB_CD = 5.0
TAB_CB = 0.5
TAB_K = 10.0 / 3.0

#: Distortion at which a droplet breaks. y = 1 means the deformation has reached the
#: droplet radius, which is where the surface can no longer hold together.
TAB_BREAKUP_DISTORTION = 1.0

#: Guard against unbounded cascades if properties are pathological.
TAB_MAX_STAGES = 40


class EvaporationRegime(StrEnum):
    """Which physical process limits the vaporization rate."""

    DIFFUSION = "diffusion_limited"
    BOILING = "boiling_limited"
    INERT = "no_evaporation"


@dataclass(frozen=True)
class LiquidState:
    """Liquid fuel properties at one temperature and pressure."""

    density_kg_m3: float
    viscosity_pa_s: float
    surface_tension_n_m: float
    specific_heat_j_kg_k: float
    latent_heat_j_kg: float
    vapor_pressure_pa: float
    molecular_weight_kg_mol: float
    saturation_temperature_k: float


class LiquidPropertyProvider(Protocol):
    """Supplies liquid properties for whichever fuel is active.

    One protocol for both fuels keeps the droplet solver fuel-agnostic: Jet-A comes
    from a measured table, LNG from CoolProp.
    """

    def liquid_state(self, temperature_k: float, pressure_pa: float) -> LiquidState: ...

    @property
    def vapor_diffusivity_reference_m2_s(self) -> float: ...

    @property
    def vapor_diffusivity_exponent(self) -> float: ...


@dataclass(frozen=True)
class GasState:
    """Local gas conditions seen by a droplet."""

    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    viscosity_pa_s: float
    conductivity_w_m_k: float
    specific_heat_j_kg_k: float
    mean_molecular_weight_kg_mol: float
    fuel_vapor_mass_fraction: float = 0.0


@dataclass(frozen=True)
class BreakupResult:
    """Outcome of the aerodynamic breakup cascade."""

    initial_radius_m: float
    final_radius_m: float
    elapsed_time_s: float
    stages: int
    occurred: bool
    weber_number: float
    limited_by_stage_cap: bool = False


def reference_temperature_k(gas_temperature_k: float, droplet_temperature_k: float) -> float:
    """Film temperature at which gas transport properties are evaluated.

    Weighted two-thirds toward the droplet because most of the resistance sits in the
    thin layer hugging the surface. From John et al. Eq. (17).
    """
    return (gas_temperature_k + 2.0 * droplet_temperature_k) / 3.0


def gas_weber_number(gas: GasState, liquid: LiquidState, radius_m: float, velocity_m_s: float
                     ) -> float:
    """Ratio of disrupting aerodynamic force to restoring surface tension."""
    return gas.density_kg_m3 * velocity_m_s**2 * radius_m / liquid.surface_tension_n_m


def _tab_frequency_squared(liquid: LiquidState, radius_m: float) -> float:
    restoring = TAB_CK * liquid.surface_tension_n_m / (liquid.density_kg_m3 * radius_m**3)
    damping = TAB_CD * liquid.viscosity_pa_s / (2.0 * liquid.density_kg_m3 * radius_m**2)
    return restoring - damping**2


def taylor_analogy_breakup(
    gas: GasState,
    provider: LiquidPropertyProvider,
    radius_m: float,
    droplet_temperature_k: float,
    relative_velocity_m_s: float,
    max_stages: int = TAB_MAX_STAGES,
) -> BreakupResult:
    """Aerodynamic breakup by the Taylor analogy, cascaded until stable.

    The droplet is modeled as a forced, damped oscillator whose distortion ``y`` grows
    under aerodynamic forcing and is resisted by surface tension. It breaks when the
    distortion reaches the droplet radius. Each child becomes the new parent and the
    test repeats, because one pass rarely reaches a stable size.

    .. note::
       John et al. print the breakup criterion as ``A + We_g > 1``. That cannot be what
       they computed: for their own Fig. 3 case ``We_g`` is about 38, so the criterion
       would fire instantly and no timescale would be resolvable. The standard Taylor
       analogy criterion, used here, is ``We_c + A > 1`` where ``We_c = We_g/12`` is the
       forced equilibrium distortion. This reproduces a finite breakup time as their
       figure shows.
    """
    initial_radius = radius_m
    liquid = provider.liquid_state(droplet_temperature_k, gas.pressure_pa)
    weber = gas_weber_number(gas, liquid, radius_m, relative_velocity_m_s)

    if relative_velocity_m_s <= 0.0 or weber <= 0.0:
        return BreakupResult(initial_radius, radius_m, 0.0, 0, False, weber)

    elapsed = 0.0
    stages = 0
    current = radius_m

    while stages < max_stages:
        omega_squared = _tab_frequency_squared(liquid, current)
        if omega_squared <= 0.0:
            # Viscosity overdamps the oscillation; the droplet cannot break this way.
            break
        omega = sqrt(omega_squared)

        local_weber = gas_weber_number(gas, liquid, current, relative_velocity_m_s)
        forced_distortion = TAB_CF * local_weber / (TAB_CK * TAB_CB)

        # Starting each stage from an undistorted, motionless droplet gives amplitude
        # equal to the forced distortion, so the peak distortion is twice it.
        peak_distortion = 2.0 * forced_distortion
        if peak_distortion <= TAB_BREAKUP_DISTORTION:
            break

        # y(t) = We_c (1 - cos(omega t)) reaches 1 at this phase.
        cos_phase = 1.0 - TAB_BREAKUP_DISTORTION / forced_distortion
        phase = float(np.arccos(np.clip(cos_phase, -1.0, 1.0)))
        stage_time = phase / omega
        distortion_rate = forced_distortion * omega * float(np.sin(phase))

        child = _post_breakup_radius(current, liquid, distortion_rate)
        if not (0.0 < child < current):
            break

        elapsed += stage_time
        current = child
        stages += 1
    else:
        return BreakupResult(
            initial_radius, current, elapsed, stages, stages > 0, weber, True
        )

    return BreakupResult(initial_radius, current, elapsed, stages, stages > 0, weber)


def _post_breakup_radius(
    radius_m: float, liquid: LiquidState, distortion_rate_per_s: float
) -> float:
    """Child radius from the O'Rourke and Amsden energy balance, their Eq. (13)."""
    surface_term = (8.0 * TAB_K / 20.0) * TAB_BREAKUP_DISTORTION**2
    kinetic_term = (
        liquid.density_kg_m3
        * radius_m**3
        * distortion_rate_per_s**2
        / liquid.surface_tension_n_m
    ) * ((6.0 * TAB_K - 5.0) / 120.0)
    return radius_m / (1.0 + surface_term + kinetic_term)


def surface_vapor_mass_fraction(
    gas: GasState, liquid: LiquidState
) -> float:
    """Fuel vapor mass fraction at the droplet surface, John et al. Eq. (21).

    Warmer droplets have higher vapor pressure, so more fuel sits at the surface and
    evaporation accelerates. When the vapor pressure reaches the ambient pressure this
    returns 1, which is the boiling point and where the diffusion model stops applying.
    """
    if liquid.vapor_pressure_pa >= gas.pressure_pa:
        return 1.0
    pressure_ratio = gas.pressure_pa / liquid.vapor_pressure_pa
    fuel_mw = liquid.molecular_weight_kg_mol
    return fuel_mw / (fuel_mw + gas.mean_molecular_weight_kg_mol * (pressure_ratio - 1.0))


def spalding_mass_transfer_number(surface_fraction: float, far_field_fraction: float) -> float:
    """Driving force for evaporation, John et al. Eq. (15)."""
    if surface_fraction >= 1.0:
        return float("inf")
    return (surface_fraction - far_field_fraction) / (1.0 - surface_fraction)


def _spalding_correction(transfer_number: float) -> float:
    """The ln(1+B)/B factor that thins the boundary layer as blowing increases."""
    if transfer_number <= 0.0:
        return 1.0
    return log(1.0 + transfer_number) / transfer_number


def rho_d_product(
    provider: LiquidPropertyProvider, film_temperature_k: float
) -> float:
    """Density-diffusivity product, John et al. Eq. (18).

    Returned as a product because that is the grouping the evaporation rate needs; it
    avoids splitting and recombining two quantities that are only known together.
    """
    return (
        1.293
        * provider.vapor_diffusivity_reference_m2_s
        * (film_temperature_k / 273.0) ** (provider.vapor_diffusivity_exponent - 1.0)
    )


@dataclass(frozen=True)
class DropletRates:
    """Instantaneous rates for one droplet, with the gas-side sources it implies."""

    radius_rate_m_s: float
    temperature_rate_k_s: float
    mass_rate_kg_s: float
    convective_heat_w: float
    regime: EvaporationRegime
    transfer_number: float


def droplet_rates(
    gas: GasState,
    provider: LiquidPropertyProvider,
    radius_m: float,
    droplet_temperature_k: float,
    relative_velocity_m_s: float,
    heat_transfer_scaling: float = 1.0,
) -> DropletRates:
    """Evaporation and heating rates, choosing the branch the physics dictates.

    ``mass_rate_kg_s`` is negative while the droplet shrinks. ``convective_heat_w`` is
    the heat the *gas* gives up, which is the sign convention the reactor source terms
    need.
    """
    if radius_m <= 0.0:
        return DropletRates(0.0, 0.0, 0.0, 0.0, EvaporationRegime.INERT, 0.0)

    liquid = provider.liquid_state(droplet_temperature_k, gas.pressure_pa)
    film_temperature = reference_temperature_k(gas.temperature_k, droplet_temperature_k)
    diameter = 2.0 * radius_m
    area = 4.0 * pi * radius_m**2
    mass = (4.0 / 3.0) * pi * radius_m**3 * liquid.density_kg_m3

    reynolds = (
        gas.density_kg_m3 * relative_velocity_m_s * diameter / gas.viscosity_pa_s
        if gas.viscosity_pa_s > 0.0
        else 0.0
    )
    prandtl = gas.specific_heat_j_kg_k * gas.viscosity_pa_s / gas.conductivity_w_m_k
    convection = 0.6 * sqrt(max(reynolds, 0.0))

    boiling = droplet_temperature_k >= liquid.saturation_temperature_k or (
        liquid.vapor_pressure_pa >= gas.pressure_pa
    )

    if boiling:
        # Heat-transfer limited. The droplet cannot get hotter than its boiling point,
        # so every joule that arrives goes into the phase change.
        superheat = gas.temperature_k - liquid.saturation_temperature_k
        if superheat <= 0.0:
            return DropletRates(0.0, 0.0, 0.0, 0.0, EvaporationRegime.INERT, 0.0)
        transfer_number = gas.specific_heat_j_kg_k * superheat / liquid.latent_heat_j_kg
        nusselt = (2.0 + convection * prandtl ** (1.0 / 3.0)) * _spalding_correction(
            transfer_number
        )
        heat_w = (
            area
            * heat_transfer_scaling
            * nusselt
            * gas.conductivity_w_m_k
            * superheat
            / diameter
        )
        mass_rate = -heat_w / liquid.latent_heat_j_kg
        radius_rate = mass_rate / (area * liquid.density_kg_m3)
        return DropletRates(
            radius_rate, 0.0, mass_rate, heat_w, EvaporationRegime.BOILING, transfer_number
        )

    # Diffusion limited. Vapor must diffuse away before more can leave the surface.
    surface_fraction = surface_vapor_mass_fraction(gas, liquid)
    transfer_number = spalding_mass_transfer_number(
        surface_fraction, gas.fuel_vapor_mass_fraction
    )
    rho_d = rho_d_product(provider, film_temperature)
    schmidt = gas.viscosity_pa_s / rho_d if rho_d > 0.0 else 0.0

    if transfer_number <= 0.0:
        # Surrounding gas is already saturated; no net mass transfer, heating only.
        radius_rate = 0.0
        mass_rate = 0.0
        correction = 1.0
    else:
        correction = _spalding_correction(transfer_number)
        sherwood = (2.0 + convection * schmidt ** (1.0 / 3.0)) * correction
        radius_rate = -(
            rho_d * transfer_number * sherwood / (2.0 * liquid.density_kg_m3 * radius_m)
        )
        mass_rate = radius_rate * area * liquid.density_kg_m3

    nusselt = (2.0 + convection * prandtl ** (1.0 / 3.0)) * correction
    heat_w = (
        area
        * heat_transfer_scaling
        * nusselt
        * gas.conductivity_w_m_k
        * (gas.temperature_k - droplet_temperature_k)
        / diameter
    )

    # Energy split: arriving heat minus what the phase change consumes warms the liquid.
    # mass_rate is negative, so the latent term correctly subtracts.
    temperature_rate = (
        (heat_w + mass_rate * liquid.latent_heat_j_kg)
        / (liquid.specific_heat_j_kg_k * mass)
        if mass > 0.0
        else 0.0
    )
    regime = (
        EvaporationRegime.DIFFUSION if mass_rate < 0.0 else EvaporationRegime.INERT
    )
    return DropletRates(
        radius_rate, temperature_rate, mass_rate, heat_w, regime, transfer_number
    )


@dataclass(frozen=True)
class DropletHistory:
    """Trajectory of one droplet class through one reactor."""

    time_s: tuple[float, ...]
    radius_m: tuple[float, ...]
    temperature_k: tuple[float, ...]
    final_radius_m: float
    final_temperature_k: float
    evaporated_mass_fraction: float
    fully_evaporated: bool
    final_regime: EvaporationRegime


def integrate_droplet(
    gas: GasState,
    provider: LiquidPropertyProvider,
    initial_radius_m: float,
    initial_temperature_k: float,
    relative_velocity_m_s: float,
    residence_time_s: float,
    heat_transfer_scaling: float = 1.0,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-12,
) -> DropletHistory:
    """Advance one droplet through a reactor with the gas state held fixed.

    Holding the gas fixed is the inner half of the operator-split scheme: the outer
    iteration updates the gas from the evaporation these trajectories produce, then
    re-integrates. Splitting keeps the stiff droplet system away from the stiff
    chemistry system, which is what makes both solvable.
    """
    if initial_radius_m <= 0.0:
        raise ValueError("Initial droplet radius must be positive")
    if residence_time_s <= 0.0:
        raise ValueError("Residence time must be positive")

    def derivatives(_t: float, state: np.ndarray) -> list[float]:
        radius, temperature = float(state[0]), float(state[1])
        if radius <= 0.0:
            return [0.0, 0.0]
        rates = droplet_rates(
            gas, provider, radius, temperature, relative_velocity_m_s, heat_transfer_scaling
        )
        return [rates.radius_rate_m_s, rates.temperature_rate_k_s]

    def fully_evaporated(_t: float, state: np.ndarray) -> float:
        return float(state[0]) - 1.0e-9

    fully_evaporated.terminal = True
    fully_evaporated.direction = -1.0

    solution = solve_ivp(
        derivatives,
        (0.0, residence_time_s),
        [initial_radius_m, initial_temperature_k],
        method="LSODA",
        events=fully_evaporated,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not solution.success:  # pragma: no cover - solver failure path
        raise RuntimeError(f"Droplet integration failed: {solution.message}")

    radii = np.maximum(solution.y[0], 0.0)
    temperatures = solution.y[1]
    final_radius = float(radii[-1])
    volume_fraction_left = (final_radius / initial_radius_m) ** 3
    final_rates = droplet_rates(
        gas,
        provider,
        max(final_radius, 1.0e-12),
        float(temperatures[-1]),
        relative_velocity_m_s,
        heat_transfer_scaling,
    )

    return DropletHistory(
        time_s=tuple(float(value) for value in solution.t),
        radius_m=tuple(float(value) for value in radii),
        temperature_k=tuple(float(value) for value in temperatures),
        final_radius_m=final_radius,
        final_temperature_k=float(temperatures[-1]),
        evaporated_mass_fraction=float(min(1.0, max(0.0, 1.0 - volume_fraction_left))),
        fully_evaporated=bool(solution.status == 1 or final_radius <= 1.0e-9),
        final_regime=final_rates.regime,
    )
