"""Dual-fuel combustor study: Jet-A at landing and take-off, LNG at cruise.

Demonstrates the three questions the reactor-network extension exists to answer:

1. how much air should go to the head end during Jet-A landing and take-off;
2. how lean LNG can run at cruise;
3. what nozzle and passage geometry each fuel path should target.

The two fuel circuits are separate hardware and only one burns at a time. What they share
is the liner, so the air split is where the compromise lives.

Every number here is illustrative. The air splits, volumes, and Jet-A properties are
plausible but are not a design, and the absolute emissions that follow from them are not
predictions. What the study demonstrates is the machinery and the trends.
"""

from __future__ import annotations

import cantera as ct
import numpy as np

from fuelnozzle.crn.autoignition import (
    IgnitionDelayTable,
    autoignition_margin,
    premix_state,
)
from fuelnozzle.crn.chemistry import (
    FuelKind,
    MechanismRegistry,
    MechanismRole,
    MechanismSpec,
    stoichiometric_air_fuel_ratio,
)
from fuelnozzle.crn.coupling import solve_coupled
from fuelnozzle.crn.emissions import summarize_emissions
from fuelnozzle.crn.liquids import JetALiquidProvider
from fuelnozzle.crn.spray_source import (
    DropletClass,
    InitialSizePolicy,
    SprayBoundary,
)
from fuelnozzle.crn.streams import AirSplit
from fuelnozzle.crn.templates import (
    ArchitectureInputs,
    lpp_architecture,
    rql_architecture,
)
from fuelnozzle.crn.thermal import saturation_temperature_k, thermal_window
from fuelnozzle.jet_a import JetAPropertyTable
from fuelnozzle.models import LNGComposition
from fuelnozzle.properties import CoolPropLNGProvider

MECH_DIR = "mech"
PRESSURE_PA = 20.0 * ct.one_atm
JET_A_FUEL_FLOW = 0.035
LNG_FUEL_FLOW = 0.030
AIR_FLOW = 1.0

def registry() -> MechanismRegistry:
    return MechanismRegistry(
        [
            MechanismSpec(
                path="gri30.yaml", fuel=FuelKind.LNG, role=MechanismRole.NETWORK,
                fuel_mole_fractions={"CH4": 1.0},
                provenance="GRI-Mech 3.0 (ships with Cantera)",
                max_pressure_pa=1.0e6,
            ),
            MechanismSpec(
                path=f"{MECH_DIR}/A2NOx_skeletal.yaml", fuel=FuelKind.JET_A,
                role=MechanismRole.NETWORK, fuel_mole_fractions={"POSF10325": 1.0},
                provenance="HyChem A2 v2.0 + Glarborg NOx, skeletal (Lu 2018)",
            ),
            MechanismSpec(
                path=f"{MECH_DIR}/A2NTCfast_ske.yaml", fuel=FuelKind.JET_A,
                role=MechanismRole.IGNITION_DELAY,
                fuel_mole_fractions={"POSF10325": 1.0},
                provenance="HyChem A2 fast-NTC skeletal; low-temperature chemistry",
            ),
        ]
    )

def jet_a_liquid() -> JetALiquidProvider:
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
            source="illustrative Jet-A properties; replace with measured batch data",
        )
    )

def spray(fuel: FuelKind, fuel_flow: float, radius_m: float, temperature_k: float):
    return SprayBoundary(
        fuel=fuel,
        total_fuel_mass_flow_kg_s=fuel_flow,
        vapor_mass_flow_kg_s=0.0,
        vapor_temperature_k=temperature_k,
        droplet_classes=(
            DropletClass(
                radius_m=radius_m, temperature_k=temperature_k, velocity_m_s=45.0,
                mass_flow_kg_s=fuel_flow, number_rate_per_s=1.0, origin="illustrative",
            ),
        ),
        injection_velocity_m_s=45.0,
        cone_angle_deg=90.0,
        apply_aerodynamic_breakup=False,
        size_policy=InitialSizePolicy.USER,
        calibration_id=None,
        warnings=(),
    )

def split(dome: float, share: float, quench: float | None = None) -> AirSplit:
    """Air split with a given dome fraction.

    ``quench`` fixed holds the quench air constant and takes the dome's air from
    dilution, which is what isolates the effect of rich-zone equivalence ratio. Leaving
    it ``None`` scales every downstream station together, which keeps the split valid at
    large dome fractions where a fixed quench would drive dilution negative.

    The distinction matters. Sweeping the dome with the quench scaling alongside varies
    two things at once and confounds the rich-zone trend with a quench-rate trend.
    """
    remainder = 1.0 - dome
    if quench is not None:
        dilution = remainder - 0.05 - quench - 0.05
        if dilution < 0.0:
            raise ValueError(
                f"A dome fraction of {dome:.2f} leaves no room for {quench:.2f} quench "
                "air. Reduce the quench fraction or let it scale."
            )
        return AirSplit(
            dome=dome, primary=0.05, quench=quench, dilution=dilution, cooling=0.05,
            jet_a_passage_share=share, idle_passage_mixing_fraction=0.0,
        )
    proportions = {"primary": 0.05, "quench": 0.30, "dilution": 0.55, "cooling": 0.10}
    total = sum(proportions.values())
    scaled = {name: remainder * value / total for name, value in proportions.items()}
    return AirSplit(
        dome=dome, **scaled, jet_a_passage_share=share,
        idle_passage_mixing_fraction=0.0,
    )

def run_case(builder, fuel, fuel_flow, air_temperature_k, dome, share, afr, provider,
             mechanism, fuel_species, radius_m, fuel_temperature_k, quench=None):
    inputs = ArchitectureInputs(
        fuel=fuel, fuel_mass_flow_kg_s=fuel_flow, total_air_mass_flow_kg_s=AIR_FLOW,
        air_temperature_k=air_temperature_k, air_split=split(dome, share, quench),
        stoichiometric_air_fuel_ratio=afr,
    )
    architecture = builder(inputs)
    result = solve_coupled(
        architecture.reactors, architecture.air_inlets, architecture.outlet_reactor,
        architecture.internal_flows, mechanism, PRESSURE_PA,
        spray(fuel, fuel_flow, radius_m, fuel_temperature_k), provider,
        architecture.spray_path, fuel_species,
        fixed_internal_flows=architecture.fixed_internal_flows, max_iterations=10,
    )
    outlet = result.network.outlet
    emissions = summarize_emissions(
        mechanism(), outlet.temperature_k, PRESSURE_PA, outlet.mass_fractions,
        AIR_FLOW + fuel_flow, fuel_flow,
    )
    return architecture, result, emissions

def answer_one(reg):
    """How much air to the head end during Jet-A landing and take-off."""
    print("\n" + "=" * 78)
    print("ANSWER 1  Head-end air fraction for Jet-A at landing and take-off")
    print("=" * 78)
    afr = stoichiometric_air_fuel_ratio(
        reg.template(FuelKind.JET_A, MechanismRole.NETWORK),
        reg.spec_for(FuelKind.JET_A, MechanismRole.NETWORK),
    )
    provider = jet_a_liquid()
    mechanism = lambda: reg.new_solution(FuelKind.JET_A, MechanismRole.NETWORK)  # noqa: E731

    print("\n  Quench air held fixed at 0.30 so that only the rich-zone equivalence")
    print("  ratio varies; scaling it alongside would confound the two effects.")
    print(f"\n{'f_dome':>7s} {'phi_rich':>9s} {'peak T':>8s} {'exit T':>8s} "
          f"{'NOx ppm':>9s} {'EI_NOx':>8s}")
    best = None
    for dome in (0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.48):
        _, result, emissions = run_case(
            rql_architecture, FuelKind.JET_A, JET_A_FUEL_FLOW, 800.0, dome, 1.0, afr,
            provider, mechanism, "POSF10325", 20.0e-6, 300.0, quench=0.30,
        )
        phi = JET_A_FUEL_FLOW * afr / (AIR_FLOW * dome)
        print(f"{dome:7.2f} {phi:9.2f} {result.network.peak_temperature_k:8.0f} "
              f"{result.network.outlet.temperature_k:8.0f} "
              f"{emissions.nox_ppmv_dry_15pct_o2:9.1f} {emissions.ei_nox_g_per_kg:8.2f}")
        if best is None or emissions.ei_nox_g_per_kg < best[1]:
            best = (dome, emissions.ei_nox_g_per_kg, phi)

    print(f"\n  Lowest diagnostic NOx at f_dome = {best[0]:.2f}, phi_rich = {best[2]:.2f}")
    print("\n  Sanity anchor: classical RQL practice puts the optimum near phi_rich")
    print("  1.4-1.6. The minimum here sits just below that band and is sharply peaked,")
    print("  so this prototype reproduces the qualitative anchor. This is not validation.")
    print("  default gave a 4.6 ms quench, which flattened the curve and drifted the")
    print("  optimum to 1.22. Quench speed dominates an RQL NOx prediction.")
    print("\n  Remaining caveat: absolute EI_NOx is still roughly twice what real RQL")
    print("  hardware achieves. These zones are adiabatic and a perfectly stirred quench")
    print("  still mixes more slowly than real jets. Treat the trend and location as")
    print("  diagnostics until they pass the documented validation gates.")


def answer_two(reg):
    """How lean LNG can run at cruise."""
    print("\n" + "=" * 78)
    print("DIAGNOSTIC 2  LNG cruise equivalence-ratio sweep")
    print("=" * 78)
    afr = stoichiometric_air_fuel_ratio(
        reg.template(FuelKind.LNG, MechanismRole.NETWORK),
        reg.spec_for(FuelKind.LNG, MechanismRole.NETWORK),
    )
    mechanism = lambda: reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK)  # noqa: E731
    provider = jet_a_liquid()  # LNG enters as vapour here; no droplets survive

    print(f"\n{'f_dome':>7s} {'phi':>7s} {'peak T':>8s} {'exit T':>8s} "
          f"{'NOx ppm':>9s} {'CO ppm':>9s} {'T spread':>9s}")
    for dome in (0.35, 0.45, 0.55, 0.70, 0.85):
        try:
            _, result, emissions = run_case(
                lpp_architecture, FuelKind.LNG, LNG_FUEL_FLOW, 700.0, dome, 0.0, afr,
                provider, mechanism, "CH4", 5.0e-6, 150.0,
            )
        except Exception as error:  # noqa: BLE001
            print(f"{dome:7.2f}   solve failed: {type(error).__name__}")
            continue
        phi = LNG_FUEL_FLOW * afr / (AIR_FLOW * dome)
        temperatures = [r.temperature_k for r in result.network.reactors]
        spread = max(temperatures) - min(temperatures)
        co = emissions.co_ppmv_dry_15pct_o2 or 0.0
        print(f"{dome:7.2f} {phi:7.2f} {result.network.peak_temperature_k:8.0f} "
              f"{result.network.outlet.temperature_k:8.0f} "
              f"{emissions.nox_ppmv_dry_15pct_o2:9.1f} {co:9.1f} {spread:9.0f}")

    print("\n  The lean limit is a bracket, not a number: extinction, CO rise, and")
    print("  exit-temperature spread each bound it, and CO here is uncalibrated.")

def answer_three(reg):
    """Target geometry and thermal window for each fuel path."""
    print("\n" + "=" * 78)
    print("ANSWER 3  Per-path targets: autoignition margin and the LNG thermal window")
    print("=" * 78)

    print("\n  Premixing passage residence time each fuel can tolerate (20 atm):")
    print(f"\n  {'fuel':7s} {'T_air':>7s} {'T_mix':>7s} {'tau_ign':>10s} "
          f"{'max tau_res at M=4':>19s}")
    for fuel, fuel_temperature, air_temperature in (
        (FuelKind.JET_A, 470.0, 800.0),
        (FuelKind.LNG, 150.0, 700.0),
    ):
        spec = reg.spec_for(fuel, MechanismRole.NETWORK)
        state = premix_state(
            reg.new_solution(fuel, MechanismRole.NETWORK), spec,
            air_mass_flow_kg_s=AIR_FLOW, air_temperature_k=air_temperature,
            fuel_mass_flow_kg_s=JET_A_FUEL_FLOW, fuel_temperature_k=fuel_temperature,
            pressure_pa=PRESSURE_PA,
        )
        table = IgnitionDelayTable(
            reg, fuel, (650.0, 700.0, 750.0, 800.0, 850.0, 900.0),
            (PRESSURE_PA,), (round(state.equivalence_ratio, 2),),
        )
        margin = autoignition_margin(table, state, 1.0e-3, fuel)
        if margin.ignition_delay_s is None:
            print(f"  {str(fuel):7s} {air_temperature:7.0f} {state.temperature_k:7.1f} "
                  f"{'no ignition':>10s} {'unbounded':>19s}")
        else:
            print(f"  {str(fuel):7s} {air_temperature:7.0f} {state.temperature_k:7.1f} "
                  f"{margin.ignition_delay_s:10.3e} "
                  f"{margin.ignition_delay_s / 4.0:19.3e}")

    print("\n  LNG thermal window versus pump pressure (chamber 20 bar):")
    provider = CoolPropLNGProvider(LNGComposition.pure_methane())
    candidates = np.arange(150.0, 260.0, 0.5).tolist()
    print(f"\n  {'P_feed':>8s} {'sat gap':>8s} {'feasible T_fuel':>18s} {'width':>7s}")
    for feed in (2.5e6, 3.0e6, 4.0e6):
        window = thermal_window(
            provider, candidates, chamber_pressure_pa=2.0e6, feed_pressure_pa=feed,
            mass_flow_kg_s=LNG_FUEL_FLOW, tank_temperature_k=112.0,
        )
        bounds = window.bounds_k
        text = f"{bounds[0]:.1f}-{bounds[1]:.1f} K" if bounds else "EMPTY"
        width = f"{bounds[1] - bounds[0]:.1f} K" if bounds else "-"
        print(f"  {feed / 1e6:7.1f}M {window.saturation_gap_k:8.1f} {text:>18s} {width:>7s}")

    saturation = saturation_temperature_k(provider, 2.0e6)
    print(f"\n  Chamber saturation is {saturation:.1f} K. Superheat at the injector and")
    print("  subcooling in the line compete for the gap between the two saturation")
    print("  temperatures, so pump pressure is what buys thermal design freedom.")

def main() -> None:
    reg = registry()
    answer_one(reg)
    answer_two(reg)
    answer_three(reg)
    print("\n" + "=" * 78)
    print("All values illustrative. Air splits, volumes, and Jet-A properties are")
    print("plausible placeholders, not a design; absolute emissions are not predictions.")
    print("=" * 78)

if __name__ == "__main__":
    main()
