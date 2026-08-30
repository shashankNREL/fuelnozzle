"""Operator-split coupling between the droplets and the reactor network.

The droplets and the gas depend on each other: evaporation is set by the gas
temperature, and the gas temperature is set partly by how much fuel has evaporated and
how much heat the droplets took to do it. Solving both at once would put a stiff droplet
system inside a stiff chemistry system, which is the quickest route to a solver that
fails for reasons nobody can diagnose.

Instead the two are alternated:

1. hold the gas fixed and march the droplets along the path they travel, accumulating
   how much fuel vapor each zone receives and how much heat it gives up;
2. hold those sources fixed and solve the reactor network to steady state;
3. repeat, moving only part of the way each time, until neither side is still changing.

Under-relaxation in step 3 matters. Evaporation and temperature reinforce each other --
hotter gas evaporates more fuel, which burns and makes the gas hotter still -- so a full
step invites oscillation between a flooded and a starved state.

The network is rebuilt on every iteration rather than mutated, because the evaporating
fuel changes the total gas mass flow and therefore the mass balance itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import cantera as ct

from fuelnozzle.crn.droplets import GasState, LiquidPropertyProvider, integrate_droplet
from fuelnozzle.crn.network import CombustorNetwork, NetworkSolution
from fuelnozzle.crn.reactors import InletSpec, OutletSpec, ReactorSpec
from fuelnozzle.crn.spray_source import SprayBoundary
from fuelnozzle.models import ModelWarning, WarningSeverity

#: Fraction of each update actually applied. Full steps oscillate.
DEFAULT_RELAXATION = 0.5

DEFAULT_TEMPERATURE_TOLERANCE_K = 1.0
DEFAULT_EVAPORATION_TOLERANCE = 1.0e-3


@dataclass(frozen=True)
class ZoneEvaporation:
    """What the droplets did inside one zone."""

    reactor_name: str
    evaporated_mass_flow_kg_s: float
    heat_drawn_from_gas_w: float
    mean_vapor_temperature_k: float
    exiting_liquid_mass_flow_kg_s: float


@dataclass(frozen=True)
class CoupledSolution:
    """Converged gas and droplet state."""

    network: NetworkSolution
    zones: tuple[ZoneEvaporation, ...]
    iterations: int
    converged: bool
    evaporated_fraction: float
    liquid_carryover_kg_s: float
    warnings: tuple[ModelWarning, ...]


def solve_coupled(
    reactors: tuple[ReactorSpec, ...] | list[ReactorSpec],
    air_inlets: tuple[InletSpec, ...] | list[InletSpec],
    outlet_reactor: str,
    internal_flows: dict[tuple[str, str], float],
    solution_factory,
    pressure_pa: float,
    spray: SprayBoundary,
    liquid_provider: LiquidPropertyProvider,
    spray_path: tuple[str, ...] | list[str],
    fuel_species: str,
    *,
    fixed_internal_flows: (
        frozenset[tuple[str, str]] | set[tuple[str, str]] | None
    ) = None,
    max_iterations: int = 20,
    relaxation: float = DEFAULT_RELAXATION,
    heat_transfer_scaling: float = 1.0,
) -> CoupledSolution:
    """Alternate droplet and network solves until both stop changing.

    ``spray_path`` is the ordered list of zones the droplets pass through. Droplets that
    survive it are reported as liquid carryover, which is a design failure rather than a
    rounding detail: fuel that never evaporates never burns.
    """
    path = tuple(spray_path)
    if not path:
        raise ValueError("A spray path with at least one reactor is required")
    known = {spec.name for spec in reactors}
    for name in path:
        if name not in known:
            raise ValueError(f"Spray path references unknown reactor {name!r}")
    if outlet_reactor not in known:
        raise ValueError(f"Unknown outlet reactor {outlet_reactor!r}")

    air_total = sum(inlet.mass_flow_kg_s for inlet in air_inlets)
    base_heat = {spec.name: spec.heat_loss_w for spec in reactors}

    # Start with every droplet already vaporized in the first zone. That gives the gas a
    # burning state to start from; the iteration then pulls it back toward reality.
    vapor: dict[str, tuple[float, float]] = {
        path[0]: (spray.liquid_mass_flow_kg_s, spray.vapor_temperature_k)
    }
    heat: dict[str, float] = {}

    warnings: list[ModelWarning] = []
    zones: tuple[ZoneEvaporation, ...] = ()
    previous_temperatures: dict[str, float] = {}
    previous_evaporated = -1.0
    evaporated_fraction = 0.0
    solution: NetworkSolution | None = None
    converged = False
    iteration = 0

    for iteration in range(1, max_iterations + 1):
        network = _build(
            reactors, air_inlets, outlet_reactor, internal_flows,
            vapor, heat, base_heat, spray, fuel_species, air_total, path[0],
            fixed_internal_flows,
        )
        solution = network.solve(solution_factory, pressure_pa)

        zones = _march_droplets(
            solution, network, spray, liquid_provider, path, pressure_pa,
            solution_factory, heat_transfer_scaling,
        )

        evaporated = sum(zone.evaporated_mass_flow_kg_s for zone in zones)
        evaporated_fraction = (
            evaporated / spray.liquid_mass_flow_kg_s
            if spray.liquid_mass_flow_kg_s > 0.0
            else 1.0
        )

        vapor = _relax_vapor(
            vapor,
            {
                zone.reactor_name: (
                    zone.evaporated_mass_flow_kg_s,
                    zone.mean_vapor_temperature_k,
                )
                for zone in zones
            },
            relaxation,
        )
        heat = _relax_heat(
            heat,
            {zone.reactor_name: zone.heat_drawn_from_gas_w for zone in zones},
            relaxation,
        )

        temperatures = {r.name: r.temperature_k for r in solution.reactors}
        moved = max(
            (
                abs(value - previous_temperatures.get(name, value - 1.0e9))
                for name, value in temperatures.items()
            ),
            default=float("inf"),
        )
        evaporation_moved = abs(evaporated_fraction - previous_evaporated)
        previous_temperatures = temperatures
        previous_evaporated = evaporated_fraction

        if (
            iteration > 1
            and moved < DEFAULT_TEMPERATURE_TOLERANCE_K
            and evaporation_moved < DEFAULT_EVAPORATION_TOLERANCE
        ):
            converged = True
            break

    assert solution is not None
    if not converged:
        warnings.append(
            ModelWarning(
                code="COUPLING_NOT_CONVERGED",
                severity=WarningSeverity.WARNING,
                message=(
                    f"Droplet and gas solutions were still moving after {iteration} "
                    "iterations. The reported state is the last iterate, not a converged "
                    "solution, and must not carry a design decision."
                ),
            )
        )

    carryover = zones[-1].exiting_liquid_mass_flow_kg_s if zones else 0.0
    total_fuel = max(spray.total_fuel_mass_flow_kg_s, 1.0e-12)
    if carryover > 1.0e-6 * total_fuel:
        warnings.append(
            ModelWarning(
                code="LIQUID_CARRYOVER",
                severity=WarningSeverity.ERROR,
                message=(
                    f"{carryover:.4g} kg/s of liquid fuel ({carryover / total_fuel:.1%} of "
                    "the total) leaves the spray path unevaporated. Fuel that never "
                    "evaporates never burns."
                ),
            )
        )

    warnings.extend(solution.warnings)
    warnings.extend(spray.warnings)
    return CoupledSolution(
        network=solution,
        zones=zones,
        iterations=iteration,
        converged=converged,
        evaporated_fraction=evaporated_fraction,
        liquid_carryover_kg_s=carryover,
        warnings=tuple(warnings),
    )


def _build(
    reactors,
    air_inlets,
    outlet_reactor: str,
    internal_flows: dict[tuple[str, str], float],
    vapor: dict[str, tuple[float, float]],
    heat: dict[str, float],
    base_heat: dict[str, float],
    spray: SprayBoundary,
    fuel_species: str,
    air_total: float,
    first_zone: str,
    fixed_internal_flows: (
        frozenset[tuple[str, str]] | set[tuple[str, str]] | None
    ),
) -> CombustorNetwork:
    """Assemble the network for one iteration.

    The convective heat the droplets absorbed is applied as a heat loss on each zone.
    Without it, injecting vapor at the droplet temperature would hand the gas back the
    latent heat it had just spent, and every zone would run hot.
    """
    inlets = list(air_inlets)
    if spray.vapor_mass_flow_kg_s > 0.0:
        inlets.append(
            InletSpec(
                name="fuel_vapor_flash",
                target_reactor=first_zone,
                mass_flow_kg_s=spray.vapor_mass_flow_kg_s,
                temperature_k=spray.vapor_temperature_k,
                mole_fractions={fuel_species: 1.0},
            )
        )
    gas_fuel = spray.vapor_mass_flow_kg_s
    for name, (flow, temperature) in vapor.items():
        if flow <= 0.0:
            continue
        gas_fuel += flow
        inlets.append(
            InletSpec(
                name=f"fuel_vapor_{name}",
                target_reactor=name,
                mass_flow_kg_s=flow,
                temperature_k=max(temperature, 1.0),
                mole_fractions={fuel_species: 1.0},
            )
        )

    specs = tuple(
        spec.model_copy(
            update={
                "heat_loss_w": base_heat.get(spec.name, 0.0)
                + max(heat.get(spec.name, 0.0), 0.0),
                "heat_loss_basis": (
                    spec.heat_loss_basis
                    or (
                        "coupled droplet sensible-plus-latent enthalpy"
                        if heat.get(spec.name, 0.0) > 0.0
                        else None
                    )
                ),
            }
        )
        for spec in reactors
    )
    return CombustorNetwork(
        specs,
        inlets,
        OutletSpec(source_reactor=outlet_reactor, mass_flow_kg_s=air_total + gas_fuel),
        internal_flows,
        fixed_internal_flows=fixed_internal_flows,
    )


def _relax_vapor(previous, proposed, relaxation):
    updated: dict[str, tuple[float, float]] = {}
    for name in set(previous) | set(proposed):
        old_flow, old_temperature = previous.get(name, (0.0, 0.0))
        new_flow, new_temperature = proposed.get(name, (0.0, old_temperature))
        updated[name] = (
            old_flow + relaxation * (new_flow - old_flow),
            new_temperature if new_temperature > 0.0 else old_temperature,
        )
    return updated


def _relax_heat(previous, proposed, relaxation):
    return {
        name: previous.get(name, 0.0)
        + relaxation * (proposed.get(name, 0.0) - previous.get(name, 0.0))
        for name in set(previous) | set(proposed)
    }


def _march_droplets(
    solution: NetworkSolution,
    network: CombustorNetwork,
    spray: SprayBoundary,
    liquid_provider: LiquidPropertyProvider,
    spray_path: tuple[str, ...],
    pressure_pa: float,
    solution_factory,
    heat_transfer_scaling: float,
) -> tuple[ZoneEvaporation, ...]:
    """Carry every droplet class along the spray path with the gas held fixed."""
    template = solution_factory()
    zones: list[ZoneEvaporation] = []

    surviving = [
        (cls.radius_m, cls.temperature_k, cls.mass_flow_kg_s, cls.velocity_m_s)
        for cls in spray.droplet_classes
    ]

    for reactor_name in spray_path:
        reactor = solution.by_name(reactor_name)
        spec = next(item for item in network.reactors if item.name == reactor_name)
        gas = _gas_state(template, reactor, pressure_pa)

        evaporated = 0.0
        heat_drawn = 0.0
        weighted_temperature = 0.0
        next_surviving: list[tuple[float, float, float, float]] = []

        for radius, temperature, mass_flow, velocity in surviving:
            if mass_flow <= 0.0 or radius <= 0.0:
                continue
            # Droplets and gas do not travel together, so droplet residence follows the
            # spray path length and droplet velocity, not the reactor residence time.
            path_length = spec.spray_path_length_m
            residence = (
                path_length / velocity
                if path_length is not None and velocity > 0.0
                else reactor.residence_time_s
            )
            residence = min(max(residence, 1.0e-9), 1.0)

            history = integrate_droplet(
                gas, liquid_provider, radius, temperature, velocity,
                residence_time_s=residence, heat_transfer_scaling=heat_transfer_scaling,
            )
            released = mass_flow * history.evaporated_mass_fraction
            remaining = mass_flow - released
            evaporated += released
            weighted_temperature += released * history.final_temperature_k

            # Heat the gas gave up, from the droplet energy balance integrated over the
            # zone: latent heat for what vaporized, plus sensible heat for what did not.
            liquid = liquid_provider.liquid_state(history.final_temperature_k, pressure_pa)
            heat_drawn += released * liquid.latent_heat_j_kg + remaining * (
                liquid.specific_heat_j_kg_k * (history.final_temperature_k - temperature)
            )

            if remaining > 0.0 and history.final_radius_m > 0.0:
                next_surviving.append(
                    (history.final_radius_m, history.final_temperature_k, remaining, velocity)
                )

        zones.append(
            ZoneEvaporation(
                reactor_name=reactor_name,
                evaporated_mass_flow_kg_s=evaporated,
                heat_drawn_from_gas_w=heat_drawn,
                mean_vapor_temperature_k=(
                    weighted_temperature / evaporated
                    if evaporated > 0.0
                    else gas.temperature_k
                ),
                exiting_liquid_mass_flow_kg_s=sum(item[2] for item in next_surviving),
            )
        )
        surviving = next_surviving

    return tuple(zones)


def _gas_state(template: ct.Solution, reactor, pressure_pa: float) -> GasState:
    """Gas conditions a droplet sees inside one reactor, with transport properties."""
    template.TPY = reactor.temperature_k, pressure_pa, reactor.mass_fractions
    return GasState(
        temperature_k=float(template.T),
        pressure_pa=pressure_pa,
        density_kg_m3=float(template.density_mass),
        viscosity_pa_s=float(template.viscosity),
        conductivity_w_m_k=float(template.thermal_conductivity),
        specific_heat_j_kg_k=float(template.cp_mass),
        mean_molecular_weight_kg_mol=float(template.mean_molecular_weight) / 1000.0,
        fuel_vapor_mass_fraction=0.0,
    )
