"""Numerical verification: limits, convergence, and the paper's trend check.

Verification asks whether the software solves the equations it claims to solve. It is
kept separate from validation, which asks whether those equations describe reality. These
tests need no experimental data.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import cantera as ct
import numpy as np
import pytest

from fuelnozzle.crn.chemistry import (
    FuelKind,
    MechanismRegistry,
    MechanismRole,
    MechanismSpec,
    equivalence_ratio,
)
from fuelnozzle.crn.coupling import solve_coupled
from fuelnozzle.crn.droplets import GasState, LiquidState, integrate_droplet
from fuelnozzle.crn.liquids import JetALiquidProvider
from fuelnozzle.crn.network import CombustorNetwork
from fuelnozzle.crn.reactors import InletSpec, OutletSpec, ReactorKind, ReactorSpec
from fuelnozzle.crn.spray_source import DropletClass, InitialSizePolicy, SprayBoundary
from fuelnozzle.jet_a import JetAPropertyTable

MECH_DIR = Path(__file__).resolve().parents[1] / "mech"
PRESSURE_PA = 20.0 * ct.one_atm


def gri30() -> ct.Solution:
    return ct.Solution("gri30.yaml")


@dataclass(frozen=True)
class ConstantLiquidProvider:
    """Fixed properties, so convergence is tested without interpolation noise."""

    state: LiquidState

    @property
    def vapor_diffusivity_reference_m2_s(self) -> float:
        return 4.16e-6

    @property
    def vapor_diffusivity_exponent(self) -> float:
        return 1.6

    def liquid_state(self, temperature_k: float, pressure_pa: float) -> LiquidState:
        return self.state


def air(temperature_k: float, pressure_pa: float = ct.one_atm) -> GasState:
    return GasState(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        density_kg_m3=pressure_pa / (287.05 * temperature_k),
        viscosity_pa_s=2.08e-5,
        conductivity_w_m_k=0.030,
        specific_heat_j_kg_k=1009.0,
        mean_molecular_weight_kg_mol=0.02897,
    )


def jet_a_provider() -> JetALiquidProvider:
    return JetALiquidProvider(
        JetAPropertyTable(
            temperature_k=(280.0, 300.0, 350.0, 400.0, 450.0, 470.0),
            density_kg_m3=(815.0, 800.0, 765.0, 730.0, 695.0, 680.0),
            viscosity_pa_s=(2.1e-3, 1.5e-3, 8.0e-4, 5.0e-4, 3.5e-4, 3.0e-4),
            surface_tension_n_m=(0.027, 0.0255, 0.0215, 0.0175, 0.0135, 0.012),
            liquid_cp_j_kg_k=(1950.0, 2000.0, 2150.0, 2300.0, 2450.0, 2500.0),
            latent_heat_j_kg=(3.6e5, 3.55e5, 3.4e5, 3.2e5, 3.0e5, 2.9e5),
            molecular_weight_kg_mol=0.15429,
            boiling_point_k=470.0,
            source="verification only",
        )
    )


# --- V1.5: the residence-time identity ------------------------------------------------


@pytest.mark.parametrize(
    ("flow_change", "expected_tau_change"),
    [(-0.264, 0.359), (0.258, -0.205)],
)
def test_residence_time_scaling_matches_the_paper(flow_change, expected_tau_change):
    """John et al. state that at fixed volume and density, tau scales as 1/mdot.

    They quote a -26.4% inflow change giving +35.9% in tau, and +25.8% giving -20.5%.
    Trivial arithmetic, but it pins our definition of residence time to theirs.
    """
    assert 1.0 / (1.0 + flow_change) - 1.0 == pytest.approx(expected_tau_change, abs=1e-3)


def test_residence_time_definition_follows_the_identity():
    """The same relation, exercised through the solver rather than by hand."""
    def solve(fuel_flow, air_flow):
        total = air_flow + fuel_flow
        network = CombustorNetwork(
            [ReactorSpec(name="psr", kind=ReactorKind.PSR, volume_m3=2.0e-3)],
            [
                InletSpec(
                    name="air", target_reactor="psr", mass_flow_kg_s=air_flow,
                    temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
                ),
                InletSpec(
                    name="fuel", target_reactor="psr", mass_flow_kg_s=fuel_flow,
                    temperature_k=300.0, mole_fractions={"CH4": 1.0},
                ),
            ],
            OutletSpec(source_reactor="psr", mass_flow_kg_s=total),
            {},
        )
        return network.solve(gri30, PRESSURE_PA).by_name("psr")

    base = solve(0.035, 1.0)
    doubled = solve(0.070, 2.0)
    ratio = base.residence_time_s / doubled.residence_time_s
    # Same composition and temperature, twice the flow, so half the residence time.
    assert ratio == pytest.approx(2.0, rel=0.02)


# --- Limiting behaviour ---------------------------------------------------------------


def test_long_residence_reactor_approaches_equilibrium():
    """A stirred reactor given unlimited time must reach the equilibrium state."""
    network = CombustorNetwork(
        [ReactorSpec(name="psr", kind=ReactorKind.PSR, volume_m3=0.5)],
        [
            InletSpec(
                name="air", target_reactor="psr", mass_flow_kg_s=1.0,
                temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
            ),
            InletSpec(
                name="fuel", target_reactor="psr", mass_flow_kg_s=0.035,
                temperature_k=750.0, mole_fractions={"CH4": 1.0},
            ),
        ],
        OutletSpec(source_reactor="psr", mass_flow_kg_s=1.035),
        {},
    )
    solved = network.solve(gri30, PRESSURE_PA).by_name("psr")

    reference = gri30()
    reference.TP = 750.0, PRESSURE_PA
    reference.set_equivalence_ratio(
        0.035 * reference.stoich_air_fuel_ratio("CH4:1", "O2:0.21,N2:0.79", basis="mole")
        / 1.0,
        "CH4:1",
        "O2:0.21,N2:0.79",
    )
    reference.equilibrate("HP")
    assert solved.temperature_k == pytest.approx(reference.T, rel=0.02)


def test_shrinking_the_reactor_extinguishes_it():
    """Below a critical residence time a stirred reactor cannot sustain combustion."""
    def peak(volume_m3):
        network = CombustorNetwork(
            [ReactorSpec(name="psr", kind=ReactorKind.PSR, volume_m3=volume_m3)],
            [
                InletSpec(
                    name="air", target_reactor="psr", mass_flow_kg_s=1.0,
                    temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
                ),
                InletSpec(
                    name="fuel", target_reactor="psr", mass_flow_kg_s=0.035,
                    temperature_k=300.0, mole_fractions={"CH4": 1.0},
                ),
            ],
            OutletSpec(source_reactor="psr", mass_flow_kg_s=1.035),
            {},
        )
        return network.solve(gri30, PRESSURE_PA).peak_temperature_k

    assert peak(2.0e-3) > 1800.0
    assert peak(1.0e-7) < 1200.0


def test_zero_fuel_leaves_the_air_unchanged():
    network = CombustorNetwork(
        [ReactorSpec(name="psr", kind=ReactorKind.PSR, volume_m3=2.0e-3)],
        [
            InletSpec(
                name="air", target_reactor="psr", mass_flow_kg_s=1.0,
                temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
            )
        ],
        OutletSpec(source_reactor="psr", mass_flow_kg_s=1.0),
        {},
    )
    solved = network.solve(gri30, PRESSURE_PA, ignition_temperature_k=750.0)
    assert solved.by_name("psr").temperature_k == pytest.approx(750.0, abs=1.0)


# --- Convergence ----------------------------------------------------------------------


def test_plug_flow_converges_with_segment_count():
    """A plug flow zone is a chain of stirred reactors; the answer must stop moving."""
    def exit_temperature(segments):
        network = CombustorNetwork(
            [
                ReactorSpec(
                    name="post",
                    kind=ReactorKind.PFR,
                    volume_m3=4.0e-3,
                    plug_flow_segments=segments,
                )
            ],
            [
                InletSpec(
                    name="air", target_reactor="post", mass_flow_kg_s=1.0,
                    temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
                ),
                InletSpec(
                    name="fuel", target_reactor="post", mass_flow_kg_s=0.035,
                    temperature_k=300.0, mole_fractions={"CH4": 1.0},
                ),
            ],
            OutletSpec(source_reactor="post", mass_flow_kg_s=1.035),
            {},
        )
        return network.solve(gri30, PRESSURE_PA).outlet.temperature_k

    coarse, medium, fine = (exit_temperature(n) for n in (2, 4, 8))
    assert abs(fine - medium) < abs(medium - coarse)
    assert abs(fine - medium) < 5.0


def test_droplet_integration_is_insensitive_to_solver_tolerance():
    """Tightening the ODE tolerance must not move the answer."""
    provider = ConstantLiquidProvider(
        LiquidState(
            density_kg_m3=800.0, viscosity_pa_s=1.5e-3, surface_tension_n_m=0.025,
            specific_heat_j_kg_k=2000.0, latent_heat_j_kg=3.5e5,
            vapor_pressure_pa=2.0e3, molecular_weight_kg_mol=0.15429,
            saturation_temperature_k=470.0,
        )
    )
    results = [
        integrate_droplet(
            air(900.0), provider, 25.0e-6, 350.0, 0.0, residence_time_s=2.0e-3,
            rtol=tolerance, atol=1.0e-14,
        ).evaporated_mass_fraction
        for tolerance in (1.0e-6, 1.0e-8, 1.0e-10)
    ]
    assert max(results) - min(results) < 1.0e-4


def test_ignition_table_converges_with_resolution():
    """A finer temperature grid must not change the interpolated delay materially."""
    from fuelnozzle.crn.autoignition import IgnitionDelayTable

    registry = MechanismRegistry(
        [
            MechanismSpec(
                path=str(MECH_DIR / "A2NTCfast_ske.yaml"), fuel=FuelKind.JET_A,
                role=MechanismRole.IGNITION_DELAY,
                fuel_mole_fractions={"POSF10325": 1.0}, provenance="verification",
            )
        ]
    )
    coarse = IgnitionDelayTable(
        registry, FuelKind.JET_A, (700.0, 800.0, 900.0), (PRESSURE_PA,), (0.5,)
    )
    fine = IgnitionDelayTable(
        registry, FuelKind.JET_A,
        tuple(np.arange(700.0, 901.0, 25.0)), (PRESSURE_PA,), (0.5,),
    )
    query = 837.0
    coarse_value, fine_value = coarse(query, PRESSURE_PA, 0.5), fine(query, PRESSURE_PA, 0.5)
    assert coarse_value is not None and fine_value is not None
    assert abs(coarse_value - fine_value) / fine_value < 0.15


# --- Tier V2: the paper's trend check --------------------------------------------------


def build_case(vaporized: bool):
    """One network, run either with droplets or as a fully premixed gaseous network.

    This is the paper's comparison. The gaseous case is what they call a GFRN: the fuel
    is assumed to arrive already gaseous *and perfectly mixed*, so every zone sees the
    same equivalence ratio. They report an 11-zone GFRN spanning only 0.652 to 0.653.

    Injecting the vapour at the dome instead would keep the air-staging heterogeneity and
    would not be their comparison; the two cases then differ by a few percent rather than
    by the order of magnitude they report.
    """
    fuel_flow, air_flow = 0.035, 1.0
    reactors = [
        ReactorSpec(
            name="evap", kind=ReactorKind.EVAPORATOR, volume_m3=0.6e-3,
            spray_path_length_m=0.03,
        ),
        ReactorSpec(
            name="mixer", kind=ReactorKind.MIXER, volume_m3=0.8e-3,
            spray_path_length_m=0.05,
        ),
        ReactorSpec(name="flame", kind=ReactorKind.PSR, volume_m3=1.5e-3),
        ReactorSpec(name="recirc", kind=ReactorKind.PSR, volume_m3=1.0e-3),
        ReactorSpec(name="post", kind=ReactorKind.PSR, volume_m3=4.0e-3),
    ]
    # Air is staged rather than all admitted at the dome. Without staging every zone
    # sees the same mixture and the network is spatially uniform, which leaves nowhere
    # for a locally rich pocket to exist and makes the comparison meaningless.
    dome_air, mixer_air, flame_air = 0.30, 0.20, 0.50
    inlets = [
        InletSpec(
            name="dome_air", target_reactor="evap", mass_flow_kg_s=air_flow * dome_air,
            temperature_k=800.0, mole_fractions={"O2": 0.21, "N2": 0.79},
        ),
        InletSpec(
            name="mixer_air", target_reactor="mixer",
            mass_flow_kg_s=air_flow * mixer_air, temperature_k=800.0,
            mole_fractions={"O2": 0.21, "N2": 0.79},
        ),
        InletSpec(
            name="flame_air", target_reactor="flame",
            mass_flow_kg_s=air_flow * flame_air, temperature_k=800.0,
            mole_fractions={"O2": 0.21, "N2": 0.79},
        ),
    ]
    # Mass-consistent by construction, including the recirculated stream that re-enters
    # at the mixer and traverses mixer->flame a second time.
    recirculated = 0.3
    flows = {
        ("evap", "mixer"): air_flow * dome_air,
        ("mixer", "flame"): air_flow * (dome_air + mixer_air) + recirculated,
        ("flame", "recirc"): recirculated,
        ("recirc", "mixer"): recirculated,
        ("flame", "post"): air_flow,
    }
    spray = SprayBoundary(
        fuel=FuelKind.JET_A,
        total_fuel_mass_flow_kg_s=fuel_flow,
        vapor_mass_flow_kg_s=fuel_flow if vaporized else 0.0,
        vapor_temperature_k=470.0 if vaporized else 300.0,
        droplet_classes=()
        if vaporized
        else (
            DropletClass(
                radius_m=12.0e-6, temperature_k=300.0, velocity_m_s=45.0,
                mass_flow_kg_s=fuel_flow, number_rate_per_s=1.0, origin="verification",
            ),
        ),
        injection_velocity_m_s=45.0, cone_angle_deg=90.0,
        apply_aerodynamic_breakup=False, size_policy=InitialSizePolicy.USER,
        calibration_id=None, warnings=(),
    )
    mechanism = lambda: ct.Solution(str(MECH_DIR / "A2NOx_skeletal.yaml"))  # noqa: E731

    if vaporized:
        # Distribute the fuel across every zone in proportion to the air it receives, so
        # the equivalence ratio is uniform. That is what "gaseous, perfectly mixed" means.
        air_by_zone = {
            "evap": air_flow * dome_air,
            "mixer": air_flow * mixer_air,
            "flame": air_flow * flame_air,
        }
        total_air = sum(air_by_zone.values())
        inlets = list(inlets) + [
            InletSpec(
                name=f"premixed_fuel_{zone}", target_reactor=zone,
                mass_flow_kg_s=fuel_flow * share / total_air, temperature_k=470.0,
                mole_fractions={"POSF10325": 1.0},
            )
            for zone, share in air_by_zone.items()
        ]
        network = CombustorNetwork(
            reactors, inlets,
            OutletSpec(source_reactor="post", mass_flow_kg_s=air_flow + fuel_flow),
            flows,
            fixed_internal_flows={("flame", "recirc"), ("recirc", "mixer")},
        )
        from fuelnozzle.crn.coupling import CoupledSolution

        return CoupledSolution(
            network=network.solve(mechanism, PRESSURE_PA), zones=(), iterations=1,
            converged=True, evaporated_fraction=1.0, liquid_carryover_kg_s=0.0,
            warnings=(),
        )

    return solve_coupled(
        reactors, inlets, "post", flows, mechanism, PRESSURE_PA, spray,
        jet_a_provider(), ("evap", "mixer"), "POSF10325",
        fixed_internal_flows={("flame", "recirc"), ("recirc", "mixer")},
        max_iterations=8,
    )


def equivalence_ratios(result) -> list[float]:
    solution = ct.Solution(str(MECH_DIR / "A2NOx_skeletal.yaml"))
    spec = MechanismSpec(
        path=str(MECH_DIR / "A2NOx_skeletal.yaml"), fuel=FuelKind.JET_A,
        role=MechanismRole.NETWORK, fuel_mole_fractions={"POSF10325": 1.0},
        provenance="verification",
    )
    values = []
    for reactor in result.network.reactors:
        solution.TPY = reactor.temperature_k, PRESSURE_PA, reactor.mass_fractions
        values.append(equivalence_ratio(solution, spec))
    return values


@pytest.fixture(scope="module")
def lfrn_and_gfrn():
    return build_case(vaporized=False), build_case(vaporized=True)


def test_premixed_case_has_a_uniform_equivalence_ratio(lfrn_and_gfrn):
    """The paper reports an 11-zone GFRN spanning only 0.652 to 0.653."""
    _, gfrn = lfrn_and_gfrn
    values = equivalence_ratios(gfrn)
    assert max(values) - min(values) < 0.01


def test_spray_case_spans_a_wide_range_of_equivalence_ratios(lfrn_and_gfrn):
    """Evaporating spray leaves the near-nozzle zone rich while the flame runs lean."""
    lfrn, _ = lfrn_and_gfrn
    values = equivalence_ratios(lfrn)
    assert max(values) - min(values) > 0.5
    assert max(values) > 1.0, "no rich pocket formed"


def test_both_cases_agree_closely_on_exit_temperature(lfrn_and_gfrn):
    """Half of the paper's result: the global energy balance barely notices the spray."""
    lfrn, gfrn = lfrn_and_gfrn
    assert lfrn.network.outlet.temperature_k == pytest.approx(
        gfrn.network.outlet.temperature_k, rel=0.01
    )


def test_nox_differs_by_orders_of_magnitude(lfrn_and_gfrn):
    """The other half, and the entire argument for modelling the spray.

    John et al. report their gaseous network underpredicting NO by 54-91%, and by more
    than an order of magnitude at the baseline. The same asymmetry appears here: exit
    temperature agrees to well under a percent while NO differs by a large factor.
    """
    lfrn, gfrn = lfrn_and_gfrn
    lfrn_no = lfrn.network.outlet.mole_fractions["NO"]
    gfrn_no = gfrn.network.outlet.mole_fractions["NO"]

    assert lfrn_no > 10.0 * gfrn_no, "the spray must change NOx by an order of magnitude"

    temperature_error = abs(
        lfrn.network.outlet.temperature_k - gfrn.network.outlet.temperature_k
    ) / gfrn.network.outlet.temperature_k
    nox_error = abs(lfrn_no - gfrn_no) / gfrn_no
    assert nox_error > 100.0 * temperature_error


def test_spray_produces_a_much_hotter_peak(lfrn_and_gfrn):
    """Local richness makes hot pockets, and NO responds exponentially to them.

    The paper reports peak temperatures above 2200 K with spray against about 1800 K
    without. The same separation appears here.
    """
    lfrn, gfrn = lfrn_and_gfrn
    assert lfrn.network.peak_temperature_k > gfrn.network.peak_temperature_k + 300.0


def test_premixed_network_is_spatially_uniform(lfrn_and_gfrn):
    """With one equivalence ratio everywhere, every zone reaches the same temperature."""
    _, gfrn = lfrn_and_gfrn
    assert gfrn.network.peak_temperature_k == pytest.approx(
        gfrn.network.outlet.temperature_k, rel=0.01
    )


def test_spray_case_fully_evaporates_so_the_comparison_is_not_confounded(lfrn_and_gfrn):
    """Unevaporated fuel would make the spray case cooler for a trivial reason."""
    lfrn, _ = lfrn_and_gfrn
    assert lfrn.evaporated_fraction > 0.99
    assert lfrn.liquid_carryover_kg_s < 1.0e-6


# --- O-005: quench sensitivity -------------------------------------------------------
#
# The rich-quench-lean NOx optimum initially fell at phi_rich = 1.22 rather than the
# 1.4-1.6 that classical practice indicates, and absolute levels were several times real
# hardware. The study below established that the cause was an unrealistically long quench
# residence time, not a defect in the reactor model.


def rql_case(dome: float, stages: int, quench_volume_m3: float):
    """Solve one RQL point with prevaporized fuel.

    Evaporation completes upstream of the quench, so prevaporizing isolates the quench
    and runs five times faster. Verified to reproduce the coupled solution to 2.7%.
    """
    from fuelnozzle.crn.chemistry import (
        emission_index_g_per_kg,
        stoichiometric_air_fuel_ratio,
    )
    from fuelnozzle.crn.streams import AirSplit
    from fuelnozzle.crn.templates import ArchitectureInputs, rql_architecture

    mechanism = lambda: ct.Solution(str(MECH_DIR / "A2NOx_skeletal.yaml"))  # noqa: E731
    spec = MechanismSpec(
        path=str(MECH_DIR / "A2NOx_skeletal.yaml"), fuel=FuelKind.JET_A,
        role=MechanismRole.NETWORK, fuel_mole_fractions={"POSF10325": 1.0},
        provenance="quench study",
    )
    afr = stoichiometric_air_fuel_ratio(mechanism(), spec)
    fuel, air = 0.035, 1.0

    architecture = rql_architecture(
        ArchitectureInputs(
            fuel=FuelKind.JET_A, fuel_mass_flow_kg_s=fuel, total_air_mass_flow_kg_s=air,
            air_temperature_k=800.0,
            air_split=AirSplit(
                dome=dome, primary=0.05, quench=0.30, dilution=1.0 - dome - 0.40,
                cooling=0.05, jet_a_passage_share=1.0,
                idle_passage_mixing_fraction=0.0,
            ),
            stoichiometric_air_fuel_ratio=afr, quench_volume_m3=quench_volume_m3,
        ),
        quench_stages=stages,
    )
    inlets = list(architecture.air_inlets) + [
        InletSpec(
            name="fuel_vapor", target_reactor=architecture.spray_path[0],
            mass_flow_kg_s=fuel, temperature_k=470.0,
            mole_fractions={"POSF10325": 1.0},
        )
    ]
    network = CombustorNetwork(
        architecture.reactors, inlets,
        OutletSpec(source_reactor=architecture.outlet_reactor, mass_flow_kg_s=air + fuel),
        architecture.internal_flows,
        fixed_internal_flows=architecture.fixed_internal_flows,
    )
    solved = network.solve(mechanism, PRESSURE_PA)
    gas = mechanism()
    gas.TPY = solved.outlet.temperature_k, PRESSURE_PA, solved.outlet.mass_fractions
    ei = emission_index_g_per_kg(
        float(gas.Y[gas.species_index("NO")]) * 46.0055 / 30.0061
        + float(gas.Y[gas.species_index("NO2")]),
        air + fuel, fuel,
    )
    return ei, solved, fuel * afr / (air * dome)


def test_single_stage_quench_badly_underpredicts_nox():
    """A single mixing point jumps over the stoichiometric crossing that makes the NOx.

    This is the measurement that justifies modelling the quench as a chain.
    """
    one_stage, _, _ = rql_case(0.38, 1, 1.0e-3)
    converged, _, _ = rql_case(0.38, 12, 1.0e-3)
    assert one_stage < 0.6 * converged


def test_quench_discretization_converges_with_stage_count():
    from fuelnozzle.crn.templates import RQL_MINIMUM_CONVERGED_STAGES

    coarse, _, _ = rql_case(0.38, 5, 1.0e-3)
    medium, _, _ = rql_case(0.38, RQL_MINIMUM_CONVERGED_STAGES, 1.0e-3)
    fine, _, _ = rql_case(0.38, 20, 1.0e-3)

    assert abs(fine - medium) / fine < 0.02
    assert abs(fine - medium) < abs(medium - coarse)


def test_slower_quench_produces_more_nox():
    """Quench speed is the dominant influence on an RQL NOx prediction."""
    fast, _, _ = rql_case(0.38, 12, 2.0e-4)
    slow, _, _ = rql_case(0.38, 12, 4.0e-3)
    assert slow > 2.0 * fast


def test_slow_quench_is_flagged():
    from fuelnozzle.crn.templates import check_quench_residence_time

    _, solved, _ = rql_case(0.38, 12, 4.0e-3)
    codes = {warning.code for warning in check_quench_residence_time(solved)}
    assert "RQL_QUENCH_TOO_SLOW" in codes


def test_under_resolved_quench_is_flagged():
    from fuelnozzle.crn.streams import AirSplit
    from fuelnozzle.crn.templates import ArchitectureInputs, rql_architecture

    architecture = rql_architecture(
        ArchitectureInputs(
            fuel=FuelKind.JET_A, fuel_mass_flow_kg_s=0.035, total_air_mass_flow_kg_s=1.0,
            air_temperature_k=800.0,
            air_split=AirSplit(
                dome=0.38, primary=0.05, quench=0.30, dilution=0.22, cooling=0.05,
                jet_a_passage_share=1.0, idle_passage_mixing_fraction=0.0,
            ),
            stoichiometric_air_fuel_ratio=14.69,
        ),
        quench_stages=3,
    )
    codes = {warning.code for warning in architecture.warnings}
    assert "RQL_QUENCH_UNDER_RESOLVED" in codes


def test_realistic_quench_puts_the_nox_optimum_in_the_expected_band():
    """O-005 resolved: with a realistic quench the optimum lands where practice says.

    At a quench of roughly one millisecond the minimum sits at phi_rich = 1.35 and the
    curve is sharply peaked. At the previous 4.6 ms default it drifted to 1.22 and the
    curve was almost flat, which is what raised the discrepancy in the first place.
    """
    results = [rql_case(dome, 12, 2.0e-4) for dome in (0.30, 0.34, 0.38, 0.42, 0.46)]
    emissions = [item[0] for item in results]
    ratios = [item[2] for item in results]

    best = ratios[emissions.index(min(emissions))]
    assert 1.25 <= best <= 1.65, f"optimum at phi_rich={best:.2f}, outside the expected band"

    # The minimum must be interior and pronounced, not a drift to the edge of the sweep.
    assert min(emissions) < 0.85 * max(emissions)
    assert emissions.index(min(emissions)) not in (0, len(emissions) - 1)
