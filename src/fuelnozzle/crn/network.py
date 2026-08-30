"""Combustor reactor network: mass balance, assembly, and steady solution.

Two problems have to be solved before any chemistry runs.

**The flow splits will not balance.** Whether they come from a user's estimate or from
post-processed CFD, the reactor-to-reactor flows will not conserve mass exactly. John et
al. report imbalances up to 9.3% on one reactor. Rather than rescaling arbitrarily, the
smallest correction that restores conservation everywhere is found, so the user's intent
is disturbed as little as the physics allows.

**Recirculation makes the network a graph, not a chain.** Hot products flowing back
toward the injector are what stabilizes the flame, so a reactor's inlet can depend on its
own downstream. All reactors are therefore handed to a single Cantera ``ReactorNet`` and
solved simultaneously, which handles the recycle without any tearing or iteration on our
side.
"""

from __future__ import annotations

from dataclasses import dataclass

import cantera as ct
import numpy as np

from fuelnozzle.crn.reactors import InletSpec, OutletSpec, ReactorKind, ReactorSpec
from fuelnozzle.models import ModelWarning, WarningSeverity

#: A correction above this fraction of the flow means the supplied splits were badly
#: inconsistent, and the user should be told rather than quietly accommodated.
LARGE_CORRECTION_FRACTION = 0.05

#: Relative tolerance for declaring the global external balance closed.
GLOBAL_BALANCE_TOLERANCE = 1.0e-9

#: A reactor this close to its inlet temperature has not ignited.
EXTINCTION_TEMPERATURE_RISE_K = 50.0

#: How many residence times to march before giving up on reaching steady state.
STEADY_RESIDENCE_MULTIPLES = 200.0

#: Residence times per convergence check.
STEADY_CHUNK_MULTIPLES = 2.0

#: Temperature movement per chunk below which the network is called converged.
STEADY_TEMPERATURE_TOLERANCE_K = 0.05


class NetworkError(RuntimeError):
    """The network is malformed or cannot be solved."""


@dataclass(frozen=True)
class MassBalanceReport:
    """Outcome of the minimum-norm correction."""

    corrected_flows: dict[tuple[str, str], float]
    correction: dict[tuple[str, str], float]
    initial_residual_kg_s: float
    final_residual_kg_s: float
    correction_norm_kg_s: float
    largest_relative_correction: float
    warnings: tuple[ModelWarning, ...]


def minimum_norm_mass_correction(
    reactor_names: tuple[str, ...] | list[str],
    internal_flows: dict[tuple[str, str], float],
    external_inflows: dict[str, float],
    external_outflows: dict[str, float],
) -> MassBalanceReport:
    """Find the smallest change to the flows that conserves mass in every reactor.

    For each reactor, what enters must equal what leaves::

        external_in_i + sum_j m_ji  =  external_out_i + sum_j m_ij

    Writing the internal flows as a vector ``z``, that is a linear system ``A z = b``.
    The supplied flows generally miss it by a residual ``r = b - A z``, and we want the
    smallest ``dz`` with ``A dz = r``. Since there are usually more flows than reactors
    the system is underdetermined, and the minimum-norm solution is exactly what a
    least-squares solve returns.

    Solving it this way rather than by scaling means no single flow absorbs the whole
    error, and flows that were already consistent are barely touched.
    """
    names = list(reactor_names)
    if len(set(names)) != len(names):
        raise NetworkError("Reactor names must be unique")
    index = {name: position for position, name in enumerate(names)}

    for source, target in internal_flows:
        for label in (source, target):
            if label not in index:
                raise NetworkError(f"Flow references unknown reactor {label!r}")

    total_in = sum(external_inflows.values())
    total_out = sum(external_outflows.values())
    scale = max(total_in, total_out, 1.0e-30)
    if abs(total_in - total_out) > GLOBAL_BALANCE_TOLERANCE * scale:
        raise NetworkError(
            "External inflow and outflow do not balance globally "
            f"({total_in:.6g} vs {total_out:.6g} kg/s). No redistribution of internal "
            "flows can fix this; the boundary conditions themselves are inconsistent."
        )

    edges = list(internal_flows)
    matrix = np.zeros((len(names), len(edges)))
    for column, (source, target) in enumerate(edges):
        matrix[index[source], column] -= 1.0
        matrix[index[target], column] += 1.0

    target_vector = np.array(
        [
            external_outflows.get(name, 0.0) - external_inflows.get(name, 0.0)
            for name in names
        ]
    )
    flows = np.array([internal_flows[edge] for edge in edges])
    residual = target_vector - matrix @ flows

    warnings: list[ModelWarning] = []
    if edges:
        adjustment, *_ = np.linalg.lstsq(matrix, residual, rcond=None)
    else:
        adjustment = np.zeros(0)

    corrected = flows + adjustment
    final_residual = target_vector - matrix @ corrected

    negative = [edges[i] for i, value in enumerate(corrected) if value < 0.0]
    if negative:
        warnings.append(
            ModelWarning(
                code="MASS_CORRECTION_REVERSED_FLOW",
                severity=WarningSeverity.WARNING,
                message=(
                    "Mass correction drove "
                    f"{', '.join(f'{a}->{b}' for a, b in negative)} negative, meaning the "
                    "supplied splits imply flow in the opposite direction. Review the "
                    "topology rather than accepting the corrected values."
                ),
            )
        )

    relative = [
        abs(adjustment[i]) / abs(flows[i])
        for i in range(len(edges))
        if abs(flows[i]) > 0.0
    ]
    largest_relative = max(relative) if relative else 0.0
    if largest_relative > LARGE_CORRECTION_FRACTION:
        warnings.append(
            ModelWarning(
                code="LARGE_MASS_CORRECTION",
                severity=WarningSeverity.WARNING,
                message=(
                    f"Restoring mass conservation required changing one flow by "
                    f"{largest_relative:.1%}. The supplied splits were substantially "
                    "inconsistent; the corrected network may not represent the intent."
                ),
            )
        )

    return MassBalanceReport(
        corrected_flows={edge: float(corrected[i]) for i, edge in enumerate(edges)},
        correction={edge: float(adjustment[i]) for i, edge in enumerate(edges)},
        initial_residual_kg_s=float(np.linalg.norm(residual)),
        final_residual_kg_s=float(np.linalg.norm(final_residual)),
        correction_norm_kg_s=float(np.linalg.norm(adjustment)),
        largest_relative_correction=largest_relative,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True)
class ReactorSolution:
    """Converged state of one reactor."""

    name: str
    kind: ReactorKind
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    mass_flow_kg_s: float
    residence_time_s: float
    mole_fractions: dict[str, float]
    mass_fractions: dict[str, float]


@dataclass(frozen=True)
class NetworkSolution:
    """Converged state of the whole network."""

    reactors: tuple[ReactorSolution, ...]
    outlet: ReactorSolution
    converged: bool
    element_balance_error: float
    warnings: tuple[ModelWarning, ...]

    def by_name(self, name: str) -> ReactorSolution:
        for reactor in self.reactors:
            if reactor.name == name:
                return reactor
        raise KeyError(name)

    @property
    def peak_temperature_k(self) -> float:
        return max(reactor.temperature_k for reactor in self.reactors)


class CombustorNetwork:
    """A mass-balanced graph of reactors, ready to solve."""

    def __init__(
        self,
        reactors: tuple[ReactorSpec, ...] | list[ReactorSpec],
        inlets: tuple[InletSpec, ...] | list[InletSpec],
        outlet: OutletSpec,
        internal_flows: dict[tuple[str, str], float],
    ) -> None:
        self.reactors = tuple(reactors)
        self.inlets = tuple(inlets)
        self.outlet = outlet
        names = tuple(reactor.name for reactor in self.reactors)
        if len(set(names)) != len(names):
            raise NetworkError("Reactor names must be unique")
        known = set(names)
        for inlet in self.inlets:
            if inlet.target_reactor not in known:
                raise NetworkError(
                    f"Inlet {inlet.name!r} targets unknown reactor "
                    f"{inlet.target_reactor!r}"
                )
        if outlet.source_reactor not in known:
            raise NetworkError(f"Outlet leaves unknown reactor {outlet.source_reactor!r}")

        external_in: dict[str, float] = {}
        for inlet in self.inlets:
            external_in[inlet.target_reactor] = (
                external_in.get(inlet.target_reactor, 0.0) + inlet.mass_flow_kg_s
            )
        external_out = {outlet.source_reactor: outlet.mass_flow_kg_s}

        self.balance = minimum_norm_mass_correction(
            names, internal_flows, external_in, external_out
        )
        self.flows = self.balance.corrected_flows

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(reactor.name for reactor in self.reactors)

    def inflow_of(self, name: str) -> float:
        """Total mass entering a reactor, external plus internal."""
        external = sum(
            inlet.mass_flow_kg_s for inlet in self.inlets if inlet.target_reactor == name
        )
        internal = sum(
            flow for (_source, target), flow in self.flows.items() if target == name
        )
        return external + internal

    def characteristic_time_s(self, pressure_pa: float, temperature_k: float) -> float:
        """Longest zone residence time, estimated before the network is solved.

        Used to size the time-integration fallback. Estimating the density from a cold
        inlet temperature overstates it, which overstates the residence time, which is
        the safe direction for a bound.
        """
        density = pressure_pa * 0.029 / (8.314462618 * max(temperature_k, 1.0))
        times = [
            spec.volume_m3 * density / inflow
            for spec in self.reactors
            if (inflow := self.inflow_of(spec.name)) > 0.0
        ]
        return max(times) if times else 1.0e-3

    def solve(
        self,
        solution_factory,
        pressure_pa: float,
        *,
        max_time_s: float | None = None,
        use_preconditioner: bool = True,
        use_steady_solver: bool = False,
        ignition_temperature_k: float = 2000.0,
    ) -> NetworkSolution:
        """Assemble the Cantera network and drive it to steady state.

        ``solution_factory`` must return a *fresh* ``Solution`` on every call. Cantera
        reactors mutate the object they are given, so a shared one would alias every
        reactor's state together.
        """
        warnings: list[ModelWarning] = list(self.balance.warnings)
        residence_scale = self.characteristic_time_s(pressure_pa, ignition_temperature_k)
        if max_time_s is None:
            max_time_s = STEADY_RESIDENCE_MULTIPLES * residence_scale

        reservoirs: dict[str, ct.Reservoir] = {}
        for inlet in self.inlets:
            gas = solution_factory()
            gas.TPX = inlet.temperature_k, pressure_pa, inlet.mole_fractions
            reservoirs[inlet.name] = ct.Reservoir(gas, name=f"inlet_{inlet.name}", clone=False)

        # Start the reactors hot. A network initialized at inlet temperature converges
        # to the trivial unlit solution, which is a valid steady state of the equations
        # and a useless answer for a combustor.
        objects: dict[str, ct.Reactor] = {}
        for spec in self.reactors:
            gas = solution_factory()
            gas.TPX = ignition_temperature_k, pressure_pa, self._initial_composition()
            gas.equilibrate("HP")
            reactor = ct.IdealGasConstPressureMoleReactor(
                gas, name=spec.name, clone=False
            )
            reactor.volume = spec.volume_m3
            objects[spec.name] = reactor

        exhaust_gas = solution_factory()
        exhaust_gas.TPX = ignition_temperature_k, pressure_pa, self._initial_composition()
        exhaust = ct.Reservoir(exhaust_gas, name="exhaust", clone=False)

        controllers: list[ct.MassFlowController] = []
        for inlet in self.inlets:
            if inlet.mass_flow_kg_s <= 0.0:
                continue
            controllers.append(
                ct.MassFlowController(
                    reservoirs[inlet.name],
                    objects[inlet.target_reactor],
                    mdot=inlet.mass_flow_kg_s,
                )
            )
        for (source, target), flow in self.flows.items():
            if flow <= 0.0:
                continue
            controllers.append(
                ct.MassFlowController(objects[source], objects[target], mdot=flow)
            )
        controllers.append(
            ct.MassFlowController(
                objects[self.outlet.source_reactor],
                exhaust,
                mdot=self.outlet.mass_flow_kg_s,
            )
        )

        walls: list[ct.Wall] = []
        for spec in self.reactors:
            if spec.heat_loss_w <= 0.0:
                continue
            ambient_gas = solution_factory()
            ambient_gas.TPX = 300.0, pressure_pa, self._initial_composition()
            ambient = ct.Reservoir(ambient_gas, name=f"ambient_{spec.name}", clone=False)
            walls.append(
                ct.Wall(objects[spec.name], ambient, A=1.0, Q=spec.heat_loss_w)
            )

        network = ct.ReactorNet(list(objects.values()))
        # Zones with a long residence time need many internal steps to cross one chunk.
        # The default cap of 20000 is reached by volumes of order a cubic metre.
        network.max_steps = 200_000
        if use_preconditioner:
            try:
                network.preconditioner = ct.AdaptivePreconditioner()
            except Exception:  # pragma: no cover - depends on reactor type support
                warnings.append(
                    ModelWarning(
                        code="PRECONDITIONER_UNAVAILABLE",
                        severity=WarningSeverity.INFO,
                        message="Sparse preconditioning is unavailable; the solve is slower.",
                    )
                )

        # Time-march to steady state rather than calling advance_to_steady_state.
        # Cantera's steady solver proved unreliable on recirculating networks with large
        # mechanisms: it ground for over three minutes and then failed on a four-reactor,
        # 71-species case that time-marches to a converged answer in under a second.
        # Marching is bounded, observable, and its convergence is checked explicitly.
        # Note that advance() takes an absolute time, not a duration.
        converged = False
        previous = np.array([reactor.phase.T for reactor in objects.values()])
        chunk = max(residence_scale * STEADY_CHUNK_MULTIPLES, 1.0e-9)
        steps = max(int(max_time_s / chunk), 1)
        try:
            for _ in range(steps):
                network.advance(network.time + chunk)
                current = np.array([reactor.phase.T for reactor in objects.values()])
                if np.max(np.abs(current - previous)) < STEADY_TEMPERATURE_TOLERANCE_K:
                    converged = True
                    break
                previous = current
        except Exception as error:  # pragma: no cover - solver failure
            raise NetworkError(f"Network solution failed: {error}") from error

        if use_steady_solver:
            try:
                network.advance_to_steady_state()
                converged = True
            except Exception as error:
                warnings.append(
                    ModelWarning(
                        code="STEADY_SOLVER_FAILED",
                        severity=WarningSeverity.INFO,
                        message=(
                            f"The optional steady-state polish failed ({error}). The "
                            "time-marched solution is reported instead."
                        ),
                    )
                )

        if not converged:
            warnings.append(
                ModelWarning(
                    code="NETWORK_NOT_CONVERGED",
                    severity=WarningSeverity.WARNING,
                    message=(
                        f"Reactor temperatures were still moving after {max_time_s:.3g} s, "
                        f"about {STEADY_RESIDENCE_MULTIPLES:.0f} residence times. The "
                        "reported state is not a converged steady solution."
                    ),
                )
            )

        solutions = tuple(
            self._extract(spec, objects[spec.name]) for spec in self.reactors
        )
        outlet_solution = next(
            solution
            for solution in solutions
            if solution.name == self.outlet.source_reactor
        )
        warnings.extend(self._diagnose(solutions))
        return NetworkSolution(
            reactors=solutions,
            outlet=outlet_solution,
            converged=converged,
            element_balance_error=self._element_balance_error(objects, solution_factory),
            warnings=tuple(warnings),
        )

    def _initial_composition(self) -> dict[str, float]:
        """Composition used to seed reactors, taken from the aggregate inlet stream."""
        totals: dict[str, float] = {}
        for inlet in self.inlets:
            weight = max(inlet.mass_flow_kg_s, 0.0)
            total_fraction = sum(inlet.mole_fractions.values())
            for species, fraction in inlet.mole_fractions.items():
                totals[species] = totals.get(species, 0.0) + weight * fraction / (
                    total_fraction or 1.0
                )
        if not totals or sum(totals.values()) <= 0.0:
            return {"N2": 1.0}
        return totals

    def _extract(self, spec: ReactorSpec, reactor: ct.Reactor) -> ReactorSolution:
        gas = reactor.phase
        inflow = self.inflow_of(spec.name)
        density = float(gas.density_mass)
        residence = spec.volume_m3 * density / inflow if inflow > 0.0 else float("inf")
        return ReactorSolution(
            name=spec.name,
            kind=spec.kind,
            temperature_k=float(gas.T),
            pressure_pa=float(gas.P),
            density_kg_m3=density,
            mass_flow_kg_s=inflow,
            residence_time_s=residence,
            mole_fractions={
                name: float(value) for name, value in zip(gas.species_names, gas.X, strict=True)
            },
            mass_fractions={
                name: float(value) for name, value in zip(gas.species_names, gas.Y, strict=True)
            },
        )

    def _diagnose(self, solutions: tuple[ReactorSolution, ...]) -> list[ModelWarning]:
        warnings: list[ModelWarning] = []
        inlet_temperature = max(
            (inlet.temperature_k for inlet in self.inlets), default=0.0
        )
        peak = max(solution.temperature_k for solution in solutions)
        if peak < inlet_temperature + EXTINCTION_TEMPERATURE_RISE_K:
            warnings.append(
                ModelWarning(
                    code="NETWORK_EXTINGUISHED",
                    severity=WarningSeverity.ERROR,
                    message=(
                        f"Peak reactor temperature {peak:.1f} K is barely above the inlet "
                        f"temperature {inlet_temperature:.1f} K. The network converged to "
                        "an unlit solution, which is a valid steady state of the equations "
                        "but not a burning combustor. Emissions from it are meaningless."
                    ),
                )
            )
        return warnings

    def _element_balance_error(
        self, objects: dict[str, ct.Reactor], solution_factory
    ) -> float:
        """Largest relative element imbalance between inflow and outflow.

        Atoms cannot be created, so this is an unforgiving check on both the flow
        bookkeeping and the solver.
        """
        template = solution_factory()
        elements = template.element_names

        inflow = np.zeros(len(elements))
        for inlet in self.inlets:
            gas = solution_factory()
            gas.TPX = inlet.temperature_k, 101_325.0, inlet.mole_fractions
            for position, element in enumerate(elements):
                inflow[position] += inlet.mass_flow_kg_s * gas.elemental_mass_fraction(
                    element
                )

        outlet_gas = objects[self.outlet.source_reactor].phase
        outflow = np.array(
            [
                self.outlet.mass_flow_kg_s * outlet_gas.elemental_mass_fraction(element)
                for element in elements
            ]
        )
        scale = np.maximum(inflow, 1.0e-30)
        return float(np.max(np.abs(outflow - inflow) / scale))
