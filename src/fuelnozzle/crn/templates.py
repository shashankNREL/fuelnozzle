"""Combustor architecture templates.

Each builder returns a complete network description: the zones, where air enters, where
exhaust leaves, the flows between zones, and the path droplets take. The four
architectures differ in one decision -- how much air meets the fuel before it burns --
and that decision is what the dual-fuel study exists to make.

Rich-quench-lean burns the fuel rich, then dumps air in fast enough to skip past
stoichiometric before much nitric oxide forms. Lean premixed mixes fuel and air fully
before burning, so no part of the flame is hot enough to make much. Lean direct injection
sits between them, injecting into a lean dome without a premixing passage.

Because the two fuel circuits are separate hardware, air flows through *both* injector
passages at all times, including the one whose fuel is shut off. Every template routes
that idle-passage air explicitly, since sizing the two passages is what lets one liner
run rich for Jet-A and lean for LNG. See ``docs/CRN_PLAN.md`` Section 8.2.1.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from fuelnozzle.crn.chemistry import DRY_AIR_MOLE_FRACTIONS, FuelKind
from fuelnozzle.crn.reactors import InletSpec, ReactorKind, ReactorSpec
from fuelnozzle.crn.streams import AirSplit
from fuelnozzle.models import ModelWarning, WarningSeverity

#: Rich-zone equivalence ratios outside this band are unusual for RQL and are flagged.
RQL_RICH_BAND = (1.2, 1.8)

#: Quench residence times outside this band are flagged. Real quench jets penetrate and
#: mix in about half a millisecond to a millisecond -- that speed is the entire point of
#: the architecture. A measured sweep (implementation log, O-005) shows the NOx optimum
#: flattening and drifting once the quench is slower than a few milliseconds, and
#: disappearing altogether by twenty.
RQL_QUENCH_TIME_BAND_S = (2.0e-4, 2.0e-3)

#: Below this many stages the quench discretization is not converged. Measured: one stage
#: underpredicts NOx by more than half, five is about 5% low, and twelve is within about
#: 1% of twenty.
RQL_MINIMUM_CONVERGED_STAGES = 12


class QuenchSchedule(StrEnum):
    """How quench air is distributed along the mixing length.

    The schedule matters because nitric oxide is made while the mixture traverses
    stoichiometric on its way from rich to lean. Front-loading the air shortens that
    traverse but risks quenching the reaction; rear-loading does the opposite.
    """

    UNIFORM = "uniform"
    FRONT_LOADED = "front_loaded"
    REAR_LOADED = "rear_loaded"

    def weights(self, stages: int) -> tuple[float, ...]:
        if stages < 1:
            raise ValueError("At least one quench stage is required")
        if self is QuenchSchedule.UNIFORM:
            raw = [1.0] * stages
        elif self is QuenchSchedule.FRONT_LOADED:
            raw = [float(stages - index) for index in range(stages)]
        else:
            raw = [float(index + 1) for index in range(stages)]
        total = sum(raw)
        return tuple(value / total for value in raw)


@dataclass(frozen=True)
class Architecture:
    """A complete network description, ready for :func:`solve_coupled`."""

    name: str
    reactors: tuple[ReactorSpec, ...]
    air_inlets: tuple[InletSpec, ...]
    outlet_reactor: str
    internal_flows: dict[tuple[str, str], float]
    spray_path: tuple[str, ...]
    near_field_equivalence_ratio: float
    warnings: tuple[ModelWarning, ...]

    @property
    def reactor_names(self) -> tuple[str, ...]:
        return tuple(spec.name for spec in self.reactors)


@dataclass(frozen=True)
class ArchitectureInputs:
    """Everything a template needs that does not depend on which template it is."""

    fuel: FuelKind
    fuel_mass_flow_kg_s: float
    total_air_mass_flow_kg_s: float
    air_temperature_k: float
    air_split: AirSplit
    stoichiometric_air_fuel_ratio: float
    dome_volume_m3: float = 1.0e-3
    mixing_volume_m3: float = 1.0e-3
    flame_volume_m3: float = 1.5e-3
    post_volume_m3: float = 4.0e-3
    #: Sized for a quench residence time near one millisecond at combustor conditions.
    #: The previous default of 1.0e-3 gave about 4.6 ms, which flattened the NOx optimum
    #: and shifted it rich. See implementation log O-005.
    quench_volume_m3: float = 2.0e-4
    recirculation_ratio: float = 0.3
    spray_path_length_m: float = 0.03
    oxidizer_mole_fractions: dict[str, float] | None = None

    @property
    def oxidizer(self) -> dict[str, float]:
        return self.oxidizer_mole_fractions or DRY_AIR_MOLE_FRACTIONS

    def air(self, fraction: float) -> float:
        return self.total_air_mass_flow_kg_s * fraction

    def inlet(self, name: str, reactor: str, fraction: float) -> InletSpec:
        return InletSpec(
            name=name,
            target_reactor=reactor,
            mass_flow_kg_s=self.air(fraction),
            temperature_k=self.air_temperature_k,
            mole_fractions=self.oxidizer,
        )


def _near_field_inlets(
    inputs: ArchitectureInputs, dome_reactor: str, downstream_reactor: str
) -> tuple[list[InletSpec], float, list[ModelWarning]]:
    """Route dome air through both injector passages.

    The active circuit's passage air always reaches the near field. The idle circuit's
    passage air is split between the near field and downstream according to the mixing
    fraction, which is the assumption a reactor network cannot settle by itself.
    """
    split = inputs.air_split
    active = split.active_passage(inputs.fuel)
    idle = split.idle_passage(inputs.fuel)
    mixing = split.idle_passage_mixing_fraction

    inlets = [inputs.inlet("dome_active_passage", dome_reactor, active)]
    if idle > 0.0:
        if mixing > 0.0:
            inlets.append(
                inputs.inlet("dome_idle_passage", dome_reactor, idle * mixing)
            )
        if mixing < 1.0:
            inlets.append(
                inputs.inlet(
                    "idle_passage_bypass", downstream_reactor, idle * (1.0 - mixing)
                )
            )

    near_field_air = inputs.air(split.near_field_air_fraction(inputs.fuel))
    phi = (
        inputs.fuel_mass_flow_kg_s
        * inputs.stoichiometric_air_fuel_ratio
        / near_field_air
        if near_field_air > 0.0
        else float("inf")
    )

    warnings: list[ModelWarning] = []
    if idle > 0.0:
        warnings.append(
            ModelWarning(
                code="IDLE_PASSAGE_AIR_ROUTED",
                severity=WarningSeverity.INFO,
                message=(
                    f"{mixing:.0%} of the idle injector passage air is routed to the near "
                    "field. This fraction is an assumption a reactor network cannot "
                    "determine; report conclusions across its full range."
                ),
            )
        )
    return inlets, phi, warnings


def rql_architecture(
    inputs: ArchitectureInputs,
    quench_stages: int = RQL_MINIMUM_CONVERGED_STAGES,
    schedule: QuenchSchedule = QuenchSchedule.UNIFORM,
) -> Architecture:
    """Rich burn, rapid quench, lean burn.

    The quench is modeled as a chain of stages rather than one mixing point. That choice
    is not cosmetic, and it is now measured: a single perfectly mixed quench point jumps
    over the stoichiometric crossing where most of the nitric oxide is made, and
    **underpredicts NOx by more than half**. The result converges to within about 1% by
    twelve stages, which is the default.

    Quench *speed* matters more than either. Taking the quench from about one millisecond
    to eighteen roughly triples NOx and destroys the optimum in rich-zone equivalence
    ratio entirely. That speed is the whole point of the architecture, and it is why the
    quench volume default is sized for a short residence time.
    """
    if quench_stages < 1:
        raise ValueError("RQL requires at least one quench stage")

    dome, mixer, recirc = "rich_dome", "rich_mixer", "rich_recirc"
    lean, post = "lean_zone", "post_flame"
    stage_names = tuple(f"quench_{index + 1}" for index in range(quench_stages))

    inlets, phi_rich, warnings = _near_field_inlets(inputs, dome, stage_names[0])
    warnings = list(warnings)

    if not RQL_RICH_BAND[0] <= phi_rich <= RQL_RICH_BAND[1]:
        warnings.append(
            ModelWarning(
                code="RQL_RICH_ZONE_OUT_OF_BAND",
                severity=WarningSeverity.WARNING,
                message=(
                    f"Rich-zone equivalence ratio is {phi_rich:.2f}, outside the "
                    f"{RQL_RICH_BAND[0]}-{RQL_RICH_BAND[1]} band where rich-quench-lean "
                    "normally minimizes NOx. Check the dome air fraction."
                ),
            )
        )

    if quench_stages < RQL_MINIMUM_CONVERGED_STAGES:
        warnings.append(
            ModelWarning(
                code="RQL_QUENCH_UNDER_RESOLVED",
                severity=WarningSeverity.WARNING,
                message=(
                    f"{quench_stages} quench stages is below the {RQL_MINIMUM_CONVERGED_STAGES} "
                    "needed for a converged result. A single stage underpredicts NOx by "
                    "more than half; five is about 5% low."
                ),
            )
        )

    weights = schedule.weights(quench_stages)
    for index, (name, weight) in enumerate(zip(stage_names, weights, strict=True)):
        inlets.append(
            inputs.inlet(f"quench_air_{index + 1}", name, inputs.air_split.quench * weight)
        )
    inlets.append(inputs.inlet("dilution_air", post, inputs.air_split.dilution))
    if inputs.air_split.primary > 0.0:
        inlets.append(inputs.inlet("primary_air", mixer, inputs.air_split.primary))
    if inputs.air_split.cooling > 0.0:
        inlets.append(inputs.inlet("cooling_air", post, inputs.air_split.cooling))

    reactors = [
        ReactorSpec(
            name=dome, kind=ReactorKind.EVAPORATOR, volume_m3=inputs.dome_volume_m3,
            spray_path_length_m=inputs.spray_path_length_m,
        ),
        ReactorSpec(
            name=mixer, kind=ReactorKind.MIXER, volume_m3=inputs.mixing_volume_m3,
            spray_path_length_m=inputs.spray_path_length_m,
        ),
        ReactorSpec(name=recirc, kind=ReactorKind.PSR, volume_m3=inputs.flame_volume_m3),
        *[
            ReactorSpec(
                name=name, kind=ReactorKind.PSR,
                volume_m3=inputs.quench_volume_m3 / quench_stages,
            )
            for name in stage_names
        ],
        ReactorSpec(name=lean, kind=ReactorKind.PSR, volume_m3=inputs.flame_volume_m3),
        ReactorSpec(name=post, kind=ReactorKind.PFR, volume_m3=inputs.post_volume_m3),
    ]

    flows = _chain_flows(
        inlets,
        inputs,
        ordered=[dome, mixer, *stage_names, lean, post],
        recirculation=((mixer, recirc, dome), inputs.recirculation_ratio),
    )
    return Architecture(
        name="rql",
        reactors=tuple(reactors),
        air_inlets=tuple(inlets),
        outlet_reactor=post,
        internal_flows=flows,
        spray_path=(dome, mixer),
        near_field_equivalence_ratio=phi_rich,
        warnings=tuple(warnings),
    )


def ldi_architecture(inputs: ArchitectureInputs) -> Architecture:
    """Lean direct injection: fuel sprays straight into a lean dome.

    There is no premixing passage, so there is nowhere for the mixture to ignite before
    it is meant to. That immunity is the point of the architecture, and it is why LDI is
    the natural comparison for Jet-A at landing and take-off, where inlet air is hot
    enough to make a premixer risky.
    """
    dome, mixer, recirc, flame, post = (
        "lean_dome", "lean_mixer", "recirc", "flame", "post_flame"
    )
    inlets, phi_near, warnings = _near_field_inlets(inputs, dome, flame)
    warnings = list(warnings)

    if inputs.air_split.primary > 0.0:
        inlets.append(inputs.inlet("primary_air", mixer, inputs.air_split.primary))
    if inputs.air_split.quench > 0.0:
        inlets.append(inputs.inlet("quench_air", flame, inputs.air_split.quench))
    inlets.append(inputs.inlet("dilution_air", post, inputs.air_split.dilution))
    if inputs.air_split.cooling > 0.0:
        inlets.append(inputs.inlet("cooling_air", post, inputs.air_split.cooling))

    reactors = [
        ReactorSpec(
            name=dome, kind=ReactorKind.EVAPORATOR, volume_m3=inputs.dome_volume_m3,
            spray_path_length_m=inputs.spray_path_length_m,
        ),
        ReactorSpec(
            name=mixer, kind=ReactorKind.MIXER, volume_m3=inputs.mixing_volume_m3,
            spray_path_length_m=inputs.spray_path_length_m,
        ),
        ReactorSpec(name=recirc, kind=ReactorKind.PSR, volume_m3=inputs.flame_volume_m3),
        ReactorSpec(name=flame, kind=ReactorKind.PSR, volume_m3=inputs.flame_volume_m3),
        ReactorSpec(name=post, kind=ReactorKind.PFR, volume_m3=inputs.post_volume_m3),
    ]
    flows = _chain_flows(
        inlets,
        inputs,
        ordered=[dome, mixer, flame, post],
        recirculation=((flame, recirc, dome), inputs.recirculation_ratio),
    )
    return Architecture(
        name="ldi",
        reactors=tuple(reactors),
        air_inlets=tuple(inlets),
        outlet_reactor=post,
        internal_flows=flows,
        spray_path=(dome, mixer),
        near_field_equivalence_ratio=phi_near,
        warnings=tuple(warnings),
    )


def lpp_architecture(
    inputs: ArchitectureInputs, premixer_volume_m3: float = 0.4e-3
) -> Architecture:
    """Lean premixed and prevaporized: mix everything before burning any of it.

    .. note::
       Chemistry is **not** frozen in the premixer, although the idea of a premixer is
       that nothing reacts there. Freezing it would guarantee the assumption rather than
       test it. Left reacting, the premixer ignites in the model whenever residence time
       exceeds ignition delay -- which is the real failure mode this architecture has,
       and the one the autoignition screen exists to predict.
    """
    premix, recirc, flame, post = "premixer", "recirc", "flame", "post_flame"
    inlets, phi_near, warnings = _near_field_inlets(inputs, premix, flame)
    warnings = list(warnings)
    warnings.append(
        ModelWarning(
            code="LPP_PREMIXER_REACTS",
            severity=WarningSeverity.INFO,
            message=(
                "The premixer is modeled as a reacting zone, not a frozen one. If it "
                "ignites here, the premixed design is invalid at this condition rather "
                "than the model being wrong."
            ),
        )
    )

    if inputs.air_split.primary > 0.0:
        inlets.append(inputs.inlet("primary_air", premix, inputs.air_split.primary))
    if inputs.air_split.quench > 0.0:
        inlets.append(inputs.inlet("quench_air", flame, inputs.air_split.quench))
    inlets.append(inputs.inlet("dilution_air", post, inputs.air_split.dilution))
    if inputs.air_split.cooling > 0.0:
        inlets.append(inputs.inlet("cooling_air", post, inputs.air_split.cooling))

    reactors = [
        ReactorSpec(
            name=premix, kind=ReactorKind.EVAPORATOR, volume_m3=premixer_volume_m3,
            spray_path_length_m=inputs.spray_path_length_m,
        ),
        ReactorSpec(name=recirc, kind=ReactorKind.PSR, volume_m3=inputs.flame_volume_m3),
        ReactorSpec(name=flame, kind=ReactorKind.PSR, volume_m3=inputs.flame_volume_m3),
        ReactorSpec(name=post, kind=ReactorKind.PFR, volume_m3=inputs.post_volume_m3),
    ]
    flows = _chain_flows(
        inlets,
        inputs,
        ordered=[premix, flame, post],
        recirculation=((flame, recirc, flame), inputs.recirculation_ratio),
    )
    return Architecture(
        name="lpp",
        reactors=tuple(reactors),
        air_inlets=tuple(inlets),
        outlet_reactor=post,
        internal_flows=flows,
        spray_path=(premix,),
        near_field_equivalence_ratio=phi_near,
        warnings=tuple(warnings),
    )


def _chain_flows(
    inlets: list[InletSpec],
    inputs: ArchitectureInputs,
    ordered: list[str],
    recirculation: tuple[tuple[str, str, str], float] | None = None,
) -> dict[tuple[str, str], float]:
    """Accumulate flows along a chain, adding each zone's air where it enters.

    Building the flows by accumulation rather than by assignment means the chain is
    mass-consistent before the network's own correction ever runs, so any correction
    that does appear points at a real inconsistency rather than at bookkeeping.
    """
    entering: dict[str, float] = {}
    for inlet in inlets:
        entering[inlet.target_reactor] = (
            entering.get(inlet.target_reactor, 0.0) + inlet.mass_flow_kg_s
        )

    flows: dict[tuple[str, str], float] = {}
    carried = 0.0
    for index, name in enumerate(ordered):
        carried += entering.get(name, 0.0)
        if index + 1 < len(ordered):
            # Fuel is added by the coupling layer, so the gas carried here is air only.
            flows[(name, ordered[index + 1])] = carried

    if recirculation is not None:
        (source, loop, target), ratio = recirculation
        rate = ratio * inputs.total_air_mass_flow_kg_s
        if rate > 0.0:
            flows[(source, loop)] = rate
            flows[(loop, target)] = rate
            # Recirculated gas re-enters upstream and travels the chain again, so every
            # segment it passes through carries it a second time. Omitting this would
            # leave an imbalance for the network correction to absorb, which would hide
            # a bookkeeping error behind a physical-looking repair.
            if target in ordered and source in ordered:
                start, stop = ordered.index(target), ordered.index(source)
                for index in range(start, stop):
                    edge = (ordered[index], ordered[index + 1])
                    if edge in flows:
                        flows[edge] += rate
    return flows


def quench_residence_time_s(solution) -> float:
    """Total residence time across the quench stages of a solved RQL network.

    Reported rather than prescribed, because it depends on the density the solution finds.
    It is the single most influential quantity in an RQL NOx prediction.
    """
    return sum(
        reactor.residence_time_s
        for reactor in solution.reactors
        if reactor.name.startswith("quench_")
    )


def check_quench_residence_time(solution) -> tuple[ModelWarning, ...]:
    """Warn when the achieved quench is slower than real hardware achieves."""
    achieved = quench_residence_time_s(solution)
    low, high = RQL_QUENCH_TIME_BAND_S
    if achieved > high:
        return (
            ModelWarning(
                code="RQL_QUENCH_TOO_SLOW",
                severity=WarningSeverity.WARNING,
                message=(
                    f"Quench residence time is {1e3 * achieved:.2f} ms, above the "
                    f"{1e3 * high:.1f} ms that real quench jets achieve. A slow quench "
                    "dwells near stoichiometric, inflating NOx and flattening its "
                    "dependence on rich-zone equivalence ratio. Reduce the quench volume."
                ),
            ),
        )
    if achieved < low:
        return (
            ModelWarning(
                code="RQL_QUENCH_VERY_FAST",
                severity=WarningSeverity.INFO,
                message=(
                    f"Quench residence time is {1e3 * achieved:.2f} ms, faster than "
                    "typical hardware. Check that the quench volume is realistic."
                ),
            ),
        )
    return ()
