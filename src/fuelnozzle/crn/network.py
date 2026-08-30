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
from enum import StrEnum
from itertools import pairwise

import cantera as ct
import numpy as np
from scipy.linalg import qr
from scipy.optimize import Bounds, LinearConstraint, linprog, minimize

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

#: Maximum movement of any species mass fraction over one convergence chunk.
STEADY_SPECIES_TOLERANCE = 1.0e-8

#: Maximum relative movement of reactor mass or enthalpy inventory per chunk.
STEADY_INVENTORY_TOLERANCE = 1.0e-8

#: Maximum state-vector derivative accumulated over one characteristic residence time.
STEADY_DERIVATIVE_TOLERANCE = 1.0e-6

#: Maximum acceptable relative element and steady energy imbalance.
CONSERVATION_TOLERANCE = 1.0e-5


class NetworkError(RuntimeError):
    """The network is malformed or cannot be solved."""


class InitializationBranch(StrEnum):
    """Initial state used to find lit and unlit steady branches."""

    HOT = "hot"
    COLD = "cold"


@dataclass(frozen=True)
class MassBalanceReport:
    """Outcome of uncertainty-weighted, nonnegative mass closure."""

    corrected_flows: dict[tuple[str, str], float]
    correction: dict[tuple[str, str], float]
    initial_residual_kg_s: float
    final_residual_kg_s: float
    correction_norm_kg_s: float
    largest_relative_correction: float
    warnings: tuple[ModelWarning, ...]


@dataclass(frozen=True)
class ConvergenceReport:
    """Largest state movement observed in the final integration chunk."""

    temperature_change_k: float
    species_mass_fraction_change: float
    relative_mass_change: float
    relative_enthalpy_inventory_change: float
    relative_physical_volume_error: float
    scaled_state_derivative: float


def minimum_norm_mass_correction(
    reactor_names: tuple[str, ...] | list[str],
    internal_flows: dict[tuple[str, str], float],
    external_inflows: dict[str, float],
    external_outflows: dict[str, float],
    *,
    fixed_internal_flows: frozenset[tuple[str, str]] | set[tuple[str, str]] | None = None,
    flow_uncertainties_kg_s: dict[tuple[str, str], float] | None = None,
) -> MassBalanceReport:
    """Close every reactor balance without reversing a declared flow.

    For each reactor, what enters must equal what leaves::

        external_in_i + sum_j m_ji  =  external_out_i + sum_j m_ij

    Writing the internal flows as a vector ``z``, that is a linear system ``A z = b``.
    The correction minimizes changes normalized by each flow's declared uncertainty.
    Corrected flows are constrained nonnegative, and named fixed flows (normally measured
    or prescribed recirculation) cannot move. An infeasible directed topology is rejected
    rather than represented by a negative controller that the solver later omits.
    """
    names = list(reactor_names)
    if len(set(names)) != len(names):
        raise NetworkError("Reactor names must be unique")
    index = {name: position for position, name in enumerate(names)}

    for source, target in internal_flows:
        for label in (source, target):
            if label not in index:
                raise NetworkError(f"Flow references unknown reactor {label!r}")
    if any(flow < 0.0 for flow in internal_flows.values()):
        raise NetworkError("Supplied internal flows must be nonnegative")
    if any(flow < 0.0 for flow in external_inflows.values()) or any(
        flow < 0.0 for flow in external_outflows.values()
    ):
        raise NetworkError("External flows must be nonnegative")

    total_in = sum(external_inflows.values())
    total_out = sum(external_outflows.values())
    scale = max(total_in, total_out, 1.0e-30)
    if abs(total_in - total_out) > GLOBAL_BALANCE_TOLERANCE * scale:
        raise NetworkError(
            "External inflow and outflow do not balance globally "
            f"({total_in:.6g} vs {total_out:.6g} kg/s). No redistribution of internal "
            "flows can fix this; the boundary conditions themselves are inconsistent."
        )
    neighbors = {name: set() for name in names}
    for source, target in internal_flows:
        neighbors[source].add(target)
        neighbors[target].add(source)
    unseen = set(names)
    while unseen:
        pending = [unseen.pop()]
        component: set[str] = set()
        while pending:
            name = pending.pop()
            component.add(name)
            discovered = neighbors[name].intersection(unseen)
            unseen.difference_update(discovered)
            pending.extend(discovered)
        component_in = sum(external_inflows.get(name, 0.0) for name in component)
        component_out = sum(external_outflows.get(name, 0.0) for name in component)
        if abs(component_in - component_out) > GLOBAL_BALANCE_TOLERANCE * scale:
            raise NetworkError(
                "External flows do not balance within connected reactor component "
                f"{sorted(component)}"
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

    fixed = set(fixed_internal_flows or ())
    unknown_fixed = fixed.difference(edges)
    if unknown_fixed:
        raise NetworkError(f"Fixed-flow set contains unknown edges: {sorted(unknown_fixed)}")
    uncertainties = flow_uncertainties_kg_s or {}
    unknown_uncertainties = set(uncertainties).difference(edges)
    if unknown_uncertainties:
        raise NetworkError(
            f"Flow uncertainties contain unknown edges: {sorted(unknown_uncertainties)}"
        )
    if any(value <= 0.0 for value in uncertainties.values()):
        raise NetworkError("Flow uncertainties must be positive")

    adjustable_indices = [i for i, edge in enumerate(edges) if edge not in fixed]
    fixed_indices = [i for i, edge in enumerate(edges) if edge in fixed]
    corrected = flows.copy()
    if np.linalg.norm(residual) > GLOBAL_BALANCE_TOLERANCE * scale:
        if not adjustable_indices:
            raise NetworkError("Fixed internal flows do not close every reactor balance")
        adjustable_matrix = matrix[:, adjustable_indices]
        fixed_contribution = (
            matrix[:, fixed_indices] @ flows[fixed_indices]
            if fixed_indices
            else np.zeros(len(names))
        )
        adjustable_target = target_vector - fixed_contribution
        independent_rank = np.linalg.matrix_rank(adjustable_matrix)
        if independent_rank < len(adjustable_indices):
            missing_uncertainties = [
                edges[index]
                for index in adjustable_indices
                if edges[index] not in uncertainties
            ]
            if missing_uncertainties:
                raise NetworkError(
                    "Underdetermined flow closure requires an uncertainty for every "
                    f"adjustable edge; missing {missing_uncertainties}"
                )
        if independent_rank:
            _, _, pivots = qr(adjustable_matrix.T, pivoting=True, mode="economic")
            independent_rows = pivots[:independent_rank]
            equality_matrix = adjustable_matrix[independent_rows]
            equality_target = adjustable_target[independent_rows]
        else:
            equality_matrix = np.empty((0, len(adjustable_indices)))
            equality_target = np.empty(0)

        feasibility = linprog(
            np.zeros(len(adjustable_indices)),
            A_eq=equality_matrix,
            b_eq=equality_target,
            bounds=(0.0, None),
            method="highs",
        )
        if not feasibility.success:
            raise NetworkError(
                "No nonnegative internal-flow solution satisfies the directed topology"
            )

        nominal = flows[adjustable_indices]
        sigma = np.array(
            [
                uncertainties.get(
                    edges[index],
                    max(0.05 * abs(flows[index]), 1.0e-6 * scale),
                )
                for index in adjustable_indices
            ]
        )

        scaled_matrix = equality_matrix * sigma[np.newaxis, :]
        scaled_target = equality_target - equality_matrix @ nominal
        initial_scaled = (feasibility.x - nominal) / sigma
        lower_scaled = -nominal / sigma

        def objective(values: np.ndarray) -> float:
            return float(0.5 * np.sum(values**2))

        def gradient(values: np.ndarray) -> np.ndarray:
            return values

        constraints = (
            [LinearConstraint(scaled_matrix, scaled_target, scaled_target)]
            if independent_rank
            else []
        )
        optimization = minimize(
            objective,
            initial_scaled,
            jac=gradient,
            method="SLSQP",
            bounds=Bounds(lower_scaled, np.inf),
            constraints=constraints,
            options={"ftol": 1.0e-12, "maxiter": 1000},
        )
        if not optimization.success:
            raise NetworkError(
                f"Nonnegative uncertainty-weighted flow closure failed: {optimization.message}"
            )
        corrected[adjustable_indices] = nominal + sigma * optimization.x

    adjustment = corrected - flows
    final_residual = target_vector - matrix @ corrected
    if np.linalg.norm(final_residual) > GLOBAL_BALANCE_TOLERANCE * scale:
        raise NetworkError("Internal-flow closure left an unresolved reactor mass residual")

    warnings: list[ModelWarning] = []
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


# New code should use the physical name; the historical symbol remains API-compatible.
close_internal_flows = minimum_norm_mass_correction


@dataclass(frozen=True)
class ReactorSolution:
    """Converged state of one reactor."""

    name: str
    kind: ReactorKind
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    volume_m3: float
    mass_kg: float
    mass_flow_kg_s: float
    residence_time_s: float
    inventory_relative_error: float
    parent_zone: str
    mole_fractions: dict[str, float]
    mass_fractions: dict[str, float]


@dataclass(frozen=True)
class NetworkSolution:
    """Converged state of the whole network."""

    reactors: tuple[ReactorSolution, ...]
    outlet: ReactorSolution
    converged: bool
    element_balance_error: float
    energy_balance_error: float
    convergence: ConvergenceReport
    initialization_branch: InitializationBranch
    warnings: tuple[ModelWarning, ...]

    def by_name(self, name: str) -> ReactorSolution:
        for reactor in self.reactors:
            if reactor.name == name:
                return reactor
        raise KeyError(name)

    @property
    def peak_temperature_k(self) -> float:
        return max(reactor.temperature_k for reactor in self.reactors)

    def zone(self, name: str) -> tuple[ReactorSolution, ...]:
        """All numerical reactors representing one physical zone."""
        matches = tuple(reactor for reactor in self.reactors if reactor.parent_zone == name)
        if not matches:
            raise KeyError(name)
        return matches

    def zone_residence_time_s(self, name: str) -> float:
        return sum(reactor.residence_time_s for reactor in self.zone(name))


@dataclass(frozen=True)
class BranchSolutions:
    """Cold- and hot-initialized steady solutions at identical conditions."""

    cold: NetworkSolution
    hot: NetworkSolution

    @property
    def distinct(self) -> bool:
        return (
            abs(self.hot.outlet.temperature_k - self.cold.outlet.temperature_k)
            >= EXTINCTION_TEMPERATURE_RISE_K
        )


def _expand_plug_flow_reactors(
    reactors: tuple[ReactorSpec, ...],
    inlets: tuple[InletSpec, ...],
    outlet: OutletSpec,
    internal_flows: dict[tuple[str, str], float],
) -> tuple[
    tuple[ReactorSpec, ...],
    tuple[InletSpec, ...],
    OutletSpec,
    dict[tuple[str, str], float],
    dict[str, str],
]:
    """Replace each physical PFR zone by a volume-preserving PSR chain."""
    first: dict[str, str] = {}
    last: dict[str, str] = {}
    parent: dict[str, str] = {}
    expanded_reactors: list[ReactorSpec] = []

    for spec in reactors:
        if spec.kind is not ReactorKind.PFR or spec.plug_flow_segments == 1:
            first[spec.name] = last[spec.name] = spec.name
            parent[spec.name] = spec.name
            expanded_reactors.append(spec)
            continue
        segment_names = [
            (
                spec.name
                if index == spec.plug_flow_segments - 1
                else f"{spec.name}__pfr_{index + 1:03d}"
            )
            for index in range(spec.plug_flow_segments)
        ]
        first[spec.name], last[spec.name] = segment_names[0], segment_names[-1]
        for segment_name in segment_names:
            parent[segment_name] = spec.name
            expanded_reactors.append(
                spec.model_copy(
                    update={
                        "name": segment_name,
                        "volume_m3": spec.segment_volume_m3,
                        "plug_flow_segments": 1,
                        "heat_loss_w": spec.heat_loss_w / spec.plug_flow_segments,
                    }
                )
            )

    remapped_inlets: list[InletSpec] = []
    for inlet in inlets:
        target_spec = next(spec for spec in reactors if spec.name == inlet.target_reactor)
        if inlet.at_reactor_exit and target_spec.kind is not ReactorKind.PFR:
            raise NetworkError(
                f"Inlet {inlet.name!r} requests a reactor exit that is not a PFR"
            )
        target = (
            last[inlet.target_reactor]
            if inlet.at_reactor_exit
            else first[inlet.target_reactor]
        )
        remapped_inlets.append(inlet.model_copy(update={"target_reactor": target}))

    expanded_flows = {
        (last[source], first[target]): flow
        for (source, target), flow in internal_flows.items()
    }
    for spec in reactors:
        if spec.kind is not ReactorKind.PFR or spec.plug_flow_segments == 1:
            continue
        segment_names = [
            reactor.name
            for reactor in expanded_reactors
            if parent[reactor.name] == spec.name
        ]
        through_flow = sum(
            inlet.mass_flow_kg_s
            for inlet in inlets
            if inlet.target_reactor == spec.name and not inlet.at_reactor_exit
        ) + sum(
            flow for (_source, target), flow in internal_flows.items() if target == spec.name
        )
        for source, target in pairwise(segment_names):
            expanded_flows[(source, target)] = through_flow

    return (
        tuple(expanded_reactors),
        tuple(remapped_inlets),
        outlet.model_copy(update={"source_reactor": last[outlet.source_reactor]}),
        expanded_flows,
        parent,
    )


class CombustorNetwork:
    """A mass-balanced graph of reactors, ready to solve."""

    def __init__(
        self,
        reactors: tuple[ReactorSpec, ...] | list[ReactorSpec],
        inlets: tuple[InletSpec, ...] | list[InletSpec],
        outlet: OutletSpec,
        internal_flows: dict[tuple[str, str], float],
        *,
        fixed_internal_flows: (
            frozenset[tuple[str, str]] | set[tuple[str, str]] | None
        ) = None,
        flow_uncertainties_kg_s: dict[tuple[str, str], float] | None = None,
    ) -> None:
        coarse_reactors = tuple(reactors)
        coarse_inlets = tuple(inlets)
        names = tuple(reactor.name for reactor in coarse_reactors)
        if len(set(names)) != len(names):
            raise NetworkError("Reactor names must be unique")
        inlet_names = [inlet.name for inlet in coarse_inlets]
        if len(set(inlet_names)) != len(inlet_names):
            raise NetworkError("Inlet names must be unique")
        known = set(names)
        for inlet in coarse_inlets:
            if inlet.target_reactor not in known:
                raise NetworkError(
                    f"Inlet {inlet.name!r} targets unknown reactor "
                    f"{inlet.target_reactor!r}"
                )
        if outlet.source_reactor not in known:
            raise NetworkError(f"Outlet leaves unknown reactor {outlet.source_reactor!r}")

        external_in: dict[str, float] = {}
        for inlet in coarse_inlets:
            external_in[inlet.target_reactor] = (
                external_in.get(inlet.target_reactor, 0.0) + inlet.mass_flow_kg_s
            )
        external_out = {outlet.source_reactor: outlet.mass_flow_kg_s}

        self.balance = minimum_norm_mass_correction(
            names,
            internal_flows,
            external_in,
            external_out,
            fixed_internal_flows=fixed_internal_flows,
            flow_uncertainties_kg_s=flow_uncertainties_kg_s,
        )
        (
            self.reactors,
            self.inlets,
            self.outlet,
            self.flows,
            self._parent_by_name,
        ) = _expand_plug_flow_reactors(
            coarse_reactors,
            coarse_inlets,
            outlet,
            self.balance.corrected_flows,
        )
        expanded_names = tuple(reactor.name for reactor in self.reactors)
        expanded_external_in: dict[str, float] = {}
        for inlet in self.inlets:
            expanded_external_in[inlet.target_reactor] = (
                expanded_external_in.get(inlet.target_reactor, 0.0)
                + inlet.mass_flow_kg_s
            )
        minimum_norm_mass_correction(
            expanded_names,
            self.flows,
            expanded_external_in,
            {self.outlet.source_reactor: self.outlet.mass_flow_kg_s},
            fixed_internal_flows=set(self.flows),
        )
        stagnant = [name for name in expanded_names if self.inflow_of(name) <= 0.0]
        if stagnant:
            raise NetworkError(f"Reactors have no through-flow: {', '.join(stagnant)}")

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
        initialization_branch: InitializationBranch = InitializationBranch.HOT,
        initial_solution: NetworkSolution | None = None,
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
        initial_composition, initial_enthalpy = self._initial_mixture(
            solution_factory, pressure_pa
        )

        reservoirs: dict[str, ct.Reservoir] = {}
        for inlet in self.inlets:
            gas = solution_factory()
            gas.TPX = inlet.temperature_k, pressure_pa, inlet.mole_fractions
            reservoirs[inlet.name] = ct.Reservoir(gas, name=f"inlet_{inlet.name}", clone=False)

        objects: dict[str, ct.Reactor] = {}
        for spec in self.reactors:
            gas = solution_factory()
            if initial_solution is not None:
                try:
                    prior = initial_solution.by_name(spec.name)
                except KeyError as error:
                    raise NetworkError(
                        "Continuation requires identical expanded reactor names; "
                        f"{spec.name!r} is absent from the prior solution"
                    ) from error
                gas.TPY = prior.temperature_k, pressure_pa, prior.mass_fractions
            elif initialization_branch is InitializationBranch.HOT:
                gas.HPX = initial_enthalpy, pressure_pa, initial_composition
                gas.equilibrate("HP")
            else:
                gas.HPX = initial_enthalpy, pressure_pa, initial_composition
            reactor = ct.IdealGasConstPressureMoleReactor(
                gas, name=spec.name, clone=False
            )
            reactor.volume = spec.volume_m3
            objects[spec.name] = reactor

        exhaust_gas = solution_factory()
        if initialization_branch is InitializationBranch.HOT:
            exhaust_gas.HPX = initial_enthalpy, pressure_pa, initial_composition
            exhaust_gas.equilibrate("HP")
        else:
            exhaust_gas.HPX = initial_enthalpy, pressure_pa, initial_composition
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
            ambient_gas.TPX = 300.0, pressure_pa, initial_composition
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
        previous = self._state_snapshot(objects)
        convergence = ConvergenceReport(
            temperature_change_k=float("inf"),
            species_mass_fraction_change=float("inf"),
            relative_mass_change=float("inf"),
            relative_enthalpy_inventory_change=float("inf"),
            relative_physical_volume_error=float("inf"),
            scaled_state_derivative=float("inf"),
        )
        chunk = max(residence_scale * STEADY_CHUNK_MULTIPLES, 1.0e-9)
        steps = max(int(max_time_s / chunk), 1)
        try:
            for _ in range(steps):
                network.advance(network.time + chunk)
                volume_error = max(
                    abs(objects[spec.name].volume - spec.volume_m3) / spec.volume_m3
                    for spec in self.reactors
                )
                # Cantera's constant-pressure reactor otherwise keeps its initial mass and
                # lets volume drift. Re-imposing the physical control volume makes the
                # inventory at every steady iterate m = rho(T, P, Y) * V_physical.
                for spec in self.reactors:
                    objects[spec.name].volume = spec.volume_m3
                network.reinitialize()
                current = self._state_snapshot(objects)
                state = network.get_state()
                derivative = network.get_derivative(1)
                derivative_residual = (
                    residence_scale
                    * float(np.linalg.norm(derivative))
                    / max(float(np.linalg.norm(state)), 1.0e-30)
                )
                convergence = self._convergence_report(
                    previous, current, volume_error, derivative_residual
                )
                if (
                    convergence.temperature_change_k
                    < STEADY_TEMPERATURE_TOLERANCE_K
                    and convergence.species_mass_fraction_change
                    < STEADY_SPECIES_TOLERANCE
                    and convergence.relative_mass_change
                    < STEADY_INVENTORY_TOLERANCE
                    and convergence.relative_enthalpy_inventory_change
                    < STEADY_INVENTORY_TOLERANCE
                    and convergence.relative_physical_volume_error
                    < STEADY_INVENTORY_TOLERANCE
                    and convergence.scaled_state_derivative
                    < STEADY_DERIVATIVE_TOLERANCE
                ):
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
        element_error = self._element_balance_error(objects, solution_factory)
        energy_error = self._energy_balance_error(objects, solution_factory, pressure_pa)
        if element_error > CONSERVATION_TOLERANCE:
            converged = False
            warnings.append(
                ModelWarning(
                    code="ELEMENT_BALANCE_FAILED",
                    severity=WarningSeverity.ERROR,
                    message=f"Relative element imbalance is {element_error:.3g}.",
                )
            )
        if energy_error > CONSERVATION_TOLERANCE:
            converged = False
            warnings.append(
                ModelWarning(
                    code="ENERGY_BALANCE_FAILED",
                    severity=WarningSeverity.ERROR,
                    message=f"Relative steady enthalpy imbalance is {energy_error:.3g}.",
                )
            )
        return NetworkSolution(
            reactors=solutions,
            outlet=outlet_solution,
            converged=converged,
            element_balance_error=element_error,
            energy_balance_error=energy_error,
            convergence=convergence,
            initialization_branch=initialization_branch,
            warnings=tuple(warnings),
        )

    def solve_branches(
        self, solution_factory, pressure_pa: float, **solve_options
    ) -> BranchSolutions:
        """Solve cold and hot initial states without conflating their steady branches."""
        cold = self.solve(
            solution_factory,
            pressure_pa,
            initialization_branch=InitializationBranch.COLD,
            **solve_options,
        )
        hot = self.solve(
            solution_factory,
            pressure_pa,
            initialization_branch=InitializationBranch.HOT,
            **solve_options,
        )
        return BranchSolutions(cold=cold, hot=hot)

    def _initial_mixture(
        self, solution_factory, pressure_pa: float
    ) -> tuple[dict[str, float], float]:
        """Aggregate species molar flows and inlet enthalpy for a physical mixed seed."""
        molar_flows: dict[str, float] = {}
        enthalpy_flow_w = 0.0
        total_mass_flow = 0.0
        for inlet in self.inlets:
            if inlet.mass_flow_kg_s <= 0.0:
                continue
            gas = solution_factory()
            gas.TPX = inlet.temperature_k, pressure_pa, inlet.mole_fractions
            molar_rate_kmol_s = inlet.mass_flow_kg_s / gas.mean_molecular_weight
            for species, fraction in zip(gas.species_names, gas.X, strict=True):
                if fraction > 0.0:
                    molar_flows[species] = (
                        molar_flows.get(species, 0.0) + molar_rate_kmol_s * fraction
                    )
            enthalpy_flow_w += inlet.mass_flow_kg_s * gas.enthalpy_mass
            total_mass_flow += inlet.mass_flow_kg_s
        if not molar_flows or total_mass_flow <= 0.0:
            raise NetworkError("At least one positive external inlet is required")
        total_molar_flow = sum(molar_flows.values())
        return (
            {
                species: molar_flow / total_molar_flow
                for species, molar_flow in molar_flows.items()
            },
            enthalpy_flow_w / total_mass_flow,
        )

    @staticmethod
    def _state_snapshot(
        objects: dict[str, ct.Reactor],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        temperatures = np.array([reactor.phase.T for reactor in objects.values()])
        species = np.concatenate([reactor.phase.Y for reactor in objects.values()])
        masses = np.array([reactor.mass for reactor in objects.values()])
        enthalpy_inventories = np.array(
            [
                reactor.mass * reactor.phase.enthalpy_mass
                for reactor in objects.values()
            ]
        )
        return temperatures, species, masses, enthalpy_inventories

    @staticmethod
    def _convergence_report(
        previous: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        current: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        relative_physical_volume_error: float,
        scaled_state_derivative: float,
    ) -> ConvergenceReport:
        temperature_change = float(np.max(np.abs(current[0] - previous[0])))
        species_change = float(np.max(np.abs(current[1] - previous[1])))

        def relative_change(before: np.ndarray, after: np.ndarray) -> float:
            scale = np.maximum(np.maximum(np.abs(before), np.abs(after)), 1.0e-30)
            return float(np.max(np.abs(after - before) / scale))

        return ConvergenceReport(
            temperature_change_k=temperature_change,
            species_mass_fraction_change=species_change,
            relative_mass_change=relative_change(previous[2], current[2]),
            relative_enthalpy_inventory_change=relative_change(
                previous[3], current[3]
            ),
            relative_physical_volume_error=relative_physical_volume_error,
            scaled_state_derivative=scaled_state_derivative,
        )

    def _extract(self, spec: ReactorSpec, reactor: ct.Reactor) -> ReactorSolution:
        gas = reactor.phase
        inflow = self.inflow_of(spec.name)
        density = float(gas.density_mass)
        volume = float(reactor.volume)
        mass = float(reactor.mass)
        expected_mass = density * volume
        inventory_error = abs(mass - expected_mass) / max(abs(expected_mass), 1.0e-30)
        residence = mass / inflow if inflow > 0.0 else float("inf")
        return ReactorSolution(
            name=spec.name,
            kind=spec.kind,
            temperature_k=float(gas.T),
            pressure_pa=float(gas.P),
            density_kg_m3=density,
            volume_m3=volume,
            mass_kg=mass,
            mass_flow_kg_s=inflow,
            residence_time_s=residence,
            inventory_relative_error=inventory_error,
            parent_zone=self._parent_by_name[spec.name],
            mole_fractions={
                name: float(value)
                for name, value in zip(gas.species_names, gas.X, strict=True)
            },
            mass_fractions={
                name: float(value)
                for name, value in zip(gas.species_names, gas.Y, strict=True)
            },
        )

    def _energy_balance_error(
        self,
        objects: dict[str, ct.Reactor],
        solution_factory,
        pressure_pa: float,
    ) -> float:
        """Relative steady-flow enthalpy residual including prescribed wall heat."""
        enthalpy_in_w = 0.0
        for inlet in self.inlets:
            gas = solution_factory()
            gas.TPX = inlet.temperature_k, pressure_pa, inlet.mole_fractions
            enthalpy_in_w += inlet.mass_flow_kg_s * gas.enthalpy_mass
        outlet_gas = objects[self.outlet.source_reactor].phase
        enthalpy_out_w = self.outlet.mass_flow_kg_s * outlet_gas.enthalpy_mass
        heat_loss_w = sum(spec.heat_loss_w for spec in self.reactors)
        residual = enthalpy_in_w - enthalpy_out_w - heat_loss_w
        scale = max(abs(enthalpy_in_w), abs(enthalpy_out_w) + heat_loss_w, 1.0)
        return abs(residual) / scale

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


def solve_continuation(
    networks: tuple[CombustorNetwork, ...] | list[CombustorNetwork],
    solution_factory,
    pressure_pa: float,
    *,
    initialization_branch: InitializationBranch,
    **solve_options,
) -> tuple[NetworkSolution, ...]:
    """Follow one initialized steady branch across an ordered network sequence."""
    solved: list[NetworkSolution] = []
    prior: NetworkSolution | None = None
    for network in networks:
        prior = network.solve(
            solution_factory,
            pressure_pa,
            initialization_branch=initialization_branch,
            initial_solution=prior,
            **solve_options,
        )
        solved.append(prior)
    return tuple(solved)
