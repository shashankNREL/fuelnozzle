"""Verification of the mechanism registry, combustion bookkeeping, and air splits."""

from __future__ import annotations

import math
from pathlib import Path

import cantera as ct
import pytest

from fuelnozzle.crn.chemistry import (
    DRY_AIR_MOLE_FRACTIONS,
    FuelKind,
    MechanismError,
    MechanismRegistry,
    MechanismRole,
    MechanismSpec,
    bilger_mixture_fraction,
    corrected_ppmv,
    dry_mole_fraction,
    emission_index_g_per_kg,
    equivalence_ratio,
    nox_pathway_coverage,
    stoichiometric_air_fuel_ratio,
    validate_mechanism,
)
from fuelnozzle.crn.streams import AirSplit, CoolingAirDestination, resolve_air_streams
from fuelnozzle.models import LNGComposition

JET_A_MECH = Path(__file__).resolve().parents[1] / "mech" / "A2NOx_skeletal.yaml"
JET_A_NTC_MECH = Path(__file__).resolve().parents[1] / "mech" / "A2NTCfast_ske.yaml"


def lng_spec() -> MechanismSpec:
    return MechanismSpec(
        path="gri30.yaml",
        fuel=FuelKind.LNG,
        role=MechanismRole.NETWORK,
        fuel_mole_fractions={"CH4": 1.0},
        provenance="GRI-Mech 3.0, ships with Cantera",
        max_pressure_pa=1.0e6,
    )


def jet_a_spec() -> MechanismSpec:
    return MechanismSpec(
        path=str(JET_A_MECH),
        fuel=FuelKind.JET_A,
        role=MechanismRole.NETWORK,
        fuel_mole_fractions={"POSF10325": 1.0},
        provenance="HyChem A2 v2.0 + Glarborg NOx, skeletal (Lu 2018)",
    )


def registry() -> MechanismRegistry:
    return MechanismRegistry([lng_spec(), jet_a_spec()])


# --- Stoichiometry against hand calculation (plan Task 1.2 acceptance) ---------------


def test_methane_stoichiometric_afr_matches_hand_calculation():
    """CH4 + 2 O2, with air at 0.21/0.79 mole, gives AFR = 17.13 by hand."""
    spec = lng_spec()
    solution = registry().template(spec.fuel, spec.role)
    afr = stoichiometric_air_fuel_ratio(solution, spec)

    moles_air = 2.0 / 0.21
    mw_air = 0.21 * 31.998 + 0.79 * 28.014
    expected = moles_air * mw_air / 16.043
    assert afr == pytest.approx(expected, rel=2.0e-3)


def test_jet_a_surrogate_stoichiometric_afr_matches_hand_calculation():
    """POSF10325 is C11H22, so it needs 16.5 O2 per mole of fuel."""
    spec = jet_a_spec()
    solution = registry().template(spec.fuel, spec.role)
    afr = stoichiometric_air_fuel_ratio(solution, spec)

    moles_air = 16.5 / 0.21
    mw_air = 0.21 * 31.998 + 0.79 * 28.014
    mw_fuel = 11 * 12.011 + 22 * 1.008
    expected = moles_air * mw_air / mw_fuel
    assert afr == pytest.approx(expected, rel=2.0e-3)


@pytest.mark.parametrize("phi", [0.4, 0.6, 1.0, 1.5])
def test_equivalence_ratio_is_recovered_exactly(phi):
    """A mixture set to a known phi must report that phi back."""
    spec = lng_spec()
    solution = registry().template(spec.fuel, spec.role)
    solution.TP = 800.0, 2.0e6
    solution.set_equivalence_ratio(phi, spec.fuel_string, DRY_AIR_MOLE_FRACTIONS)
    assert equivalence_ratio(solution, spec) == pytest.approx(phi, rel=1.0e-9)


def test_bilger_mixture_fraction_is_consistent_with_equivalence_ratio():
    """Z and phi are two views of the same mixture and must agree."""
    spec = lng_spec()
    solution = registry().template(spec.fuel, spec.role)
    solution.TP = 800.0, 2.0e6
    solution.set_equivalence_ratio(0.6, spec.fuel_string, DRY_AIR_MOLE_FRACTIONS)

    z = bilger_mixture_fraction(solution, spec)
    z_stoich = 1.0 / (1.0 + stoichiometric_air_fuel_ratio(solution, spec))
    phi_from_z = (z / (1.0 - z)) * ((1.0 - z_stoich) / z_stoich)
    assert phi_from_z == pytest.approx(0.6, rel=1.0e-6)


def test_equivalence_ratio_survives_partial_burning():
    """The element-based ratio must stay correct after reaction consumes the fuel."""
    spec = lng_spec()
    solution = registry().new_solution(spec.fuel, spec.role)
    solution.TP = 1800.0, 2.0e6
    solution.set_equivalence_ratio(0.8, spec.fuel_string, DRY_AIR_MOLE_FRACTIONS)

    reactor = ct.IdealGasConstPressureReactor(solution, clone=False)
    ct.ReactorNet([reactor]).advance(0.05)

    assert solution["CH4"].X[0] < 1.0e-6, "fuel should be consumed"
    assert equivalence_ratio(solution, spec) == pytest.approx(0.8, rel=1.0e-6)


# --- Mechanism validation (plan Task 1.3 acceptance) --------------------------------


def test_mechanism_without_nitrogen_chemistry_raises_when_nox_requested():
    """A silent zero NO is worse than a failed run, so this must raise."""
    spec = MechanismSpec(
        path=str(JET_A_NTC_MECH),
        fuel=FuelKind.JET_A,
        role=MechanismRole.IGNITION_DELAY,
        fuel_mole_fractions={"POSF10325": 1.0},
        provenance="HyChem A2 fast-NTC skeletal; no N chemistry",
    )
    solution = ct.Solution(spec.path)

    with pytest.raises(MechanismError, match="lacks nitrogen chemistry"):
        validate_mechanism(spec, solution, require_nox=True)

    assert validate_mechanism(spec, solution, require_nox=False) == ()


def test_gri30_warns_above_its_validity_pressure():
    spec = lng_spec()
    solution = registry().template(spec.fuel, spec.role)
    warnings = validate_mechanism(spec, solution, require_nox=True, pressure_pa=3.0e6)
    codes = {warning.code for warning in warnings}
    assert "MECHANISM_PRESSURE_ABOVE_VALIDITY" in codes


def test_unrepresented_lng_component_raises_rather_than_being_dropped():
    spec = lng_spec()
    solution = registry().template(spec.fuel, spec.role)
    composition = LNGComposition(mole_fractions={"Methane": 0.9, "Ethane": 0.1})

    with pytest.raises(MechanismError, match="absent from mechanism"):
        validate_mechanism(spec, solution, require_nox=True, lng_composition=composition)


def test_missing_fuel_species_raises():
    spec = MechanismSpec(
        path="gri30.yaml",
        fuel=FuelKind.LNG,
        role=MechanismRole.NETWORK,
        fuel_mole_fractions={"POSF10325": 1.0},
        provenance="deliberately wrong fuel for this mechanism",
    )
    solution = ct.Solution("gri30.yaml")
    with pytest.raises(MechanismError, match="not present in mechanism"):
        validate_mechanism(spec, solution, require_nox=False)


def test_nox_pathway_coverage_differs_between_the_two_fuel_mechanisms():
    """Neither mechanism is a superset of the other; the code must report, not assume."""
    jet_a = nox_pathway_coverage(ct.Solution(str(JET_A_MECH)))
    lng = nox_pathway_coverage(ct.Solution("gri30.yaml"))

    assert jet_a.prompt_ncn and not jet_a.nnh_route
    assert lng.nnh_route and not lng.prompt_ncn
    assert jet_a.thermal and lng.thermal
    assert jet_a.n2o_route and lng.n2o_route


# --- Registry behaviour --------------------------------------------------------------


def test_ignition_role_falls_back_to_network_mechanism():
    reg = registry()
    assert not reg.has_dedicated_ignition_mechanism(FuelKind.JET_A)
    fallback = reg.spec_for(FuelKind.JET_A, MechanismRole.IGNITION_DELAY)
    assert fallback.role is MechanismRole.NETWORK


def test_new_solution_returns_independent_objects():
    """Reactors mutate their Solution, so the registry must not hand out one shared copy."""
    reg = registry()
    first = reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK)
    second = reg.new_solution(FuelKind.LNG, MechanismRole.NETWORK)
    first.TP = 500.0, ct.one_atm
    second.TP = 900.0, ct.one_atm
    assert first.T == pytest.approx(500.0)


def test_duplicate_mechanism_registration_raises():
    with pytest.raises(MechanismError, match="Duplicate mechanism"):
        MechanismRegistry([lng_spec(), lng_spec()])


# --- Emissions conversions -----------------------------------------------------------


def test_emission_index_is_grams_per_kilogram_of_fuel():
    assert emission_index_g_per_kg(1.0e-3, 20.0, 0.5) == pytest.approx(40.0)


def test_corrected_ppmv_is_identity_at_the_reference_oxygen_level():
    assert corrected_ppmv(1.0e-5, 0.15) == pytest.approx(10.0)


def test_corrected_ppmv_scales_up_a_more_dilute_sample():
    """Leaner exhaust has more O2, so the same raw ppm means more real emission."""
    lean = corrected_ppmv(1.0e-5, 0.18)
    assert lean > corrected_ppmv(1.0e-5, 0.15)
    assert lean == pytest.approx(1.0e6 * 1.0e-5 * (20.9 - 15.0) / (20.9 - 18.0))


def test_corrected_ppmv_rejects_oxygen_at_or_above_air():
    with pytest.raises(ValueError, match="undefined"):
        corrected_ppmv(1.0e-5, 0.209)


def test_dry_mole_fraction_removes_water():
    solution = ct.Solution("gri30.yaml")
    solution.TPX = 1500.0, ct.one_atm, {"H2O": 0.1, "NO": 0.9}
    assert dry_mole_fraction(solution, "NO") == pytest.approx(0.9 / 0.9)


# --- Air splits ----------------------------------------------------------------------


def base_split(**overrides) -> AirSplit:
    values = {
        "dome": 0.20,
        "primary": 0.15,
        "quench": 0.35,
        "dilution": 0.20,
        "cooling": 0.10,
        "cooling_destination": CoolingAirDestination.DILUTION,
    }
    values.update(overrides)
    return AirSplit(**values)


def test_air_split_must_sum_to_one():
    with pytest.raises(ValueError, match="must sum to 1.0"):
        AirSplit(dome=0.2, primary=0.2, quench=0.2, dilution=0.2, cooling=0.3)


def test_injector_passages_partition_the_dome_air():
    split = base_split(jet_a_passage_share=0.25)
    assert split.jet_a_passage == pytest.approx(0.05)
    assert split.lng_passage == pytest.approx(0.15)
    assert split.jet_a_passage + split.lng_passage == pytest.approx(split.dome)


def test_passage_sizing_shifts_near_field_air_between_fuels():
    """The Section 8.2.1 lever: fixed hardware, different effective near-field air.

    A small Jet-A passage and a large LNG passage give a rich Jet-A near field and a
    lean LNG near field without any moving parts.
    """
    split = base_split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0)
    jet_a_air = split.near_field_air_fraction(FuelKind.JET_A)
    lng_air = split.near_field_air_fraction(FuelKind.LNG)
    assert jet_a_air == pytest.approx(0.05)
    assert lng_air == pytest.approx(0.15)
    assert lng_air > jet_a_air


def test_idle_passage_mixing_fraction_bounds_the_near_field_air():
    """The lever's strength depends entirely on an assumption the CRN cannot settle."""
    segregated = base_split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0)
    mixed = base_split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=1.0)

    assert segregated.near_field_air_fraction(FuelKind.JET_A) == pytest.approx(0.05)
    assert mixed.near_field_air_fraction(FuelKind.JET_A) == pytest.approx(0.20)


def test_resolve_air_streams_from_overall_equivalence_ratio():
    spec = lng_spec()
    solution = registry().new_solution(spec.fuel, spec.role)
    split = base_split(jet_a_passage_share=0.25, idle_passage_mixing_fraction=0.0)

    streams = resolve_air_streams(
        solution,
        spec,
        split,
        fuel=FuelKind.LNG,
        fuel_mass_flow_kg_s=0.08,
        temperature_k=750.0,
        pressure_pa=2.0e6,
        overall_equivalence_ratio=0.35,
    )

    afr = stoichiometric_air_fuel_ratio(solution, spec)
    assert streams.total_mass_flow_kg_s == pytest.approx(0.08 * afr / 0.35)
    assert streams.overall_equivalence_ratio == pytest.approx(0.35)

    stations = (
        streams.dome_kg_s
        + streams.primary_kg_s
        + streams.quench_kg_s
        + streams.dilution_kg_s
        + streams.cooling_kg_s
    )
    assert stations == pytest.approx(streams.total_mass_flow_kg_s)
    assert streams.active_passage_kg_s == pytest.approx(streams.total_mass_flow_kg_s * 0.15)
    assert streams.near_field_equivalence_ratio > streams.overall_equivalence_ratio


def test_resolve_air_streams_round_trips_explicit_air_flow():
    spec = lng_spec()
    solution = registry().new_solution(spec.fuel, spec.role)
    streams = resolve_air_streams(
        solution,
        spec,
        base_split(),
        fuel=FuelKind.LNG,
        fuel_mass_flow_kg_s=0.08,
        temperature_k=750.0,
        pressure_pa=2.0e6,
        total_air_mass_flow_kg_s=4.0,
    )
    afr = stoichiometric_air_fuel_ratio(solution, spec)
    assert streams.overall_equivalence_ratio == pytest.approx(0.08 * afr / 4.0)


def test_resolve_air_streams_requires_exactly_one_air_specification():
    spec = lng_spec()
    solution = registry().new_solution(spec.fuel, spec.role)
    kwargs = dict(
        fuel=FuelKind.LNG,
        fuel_mass_flow_kg_s=0.08,
        temperature_k=750.0,
        pressure_pa=2.0e6,
    )
    with pytest.raises(ValueError, match="exactly one"):
        resolve_air_streams(solution, spec, base_split(), **kwargs)
    with pytest.raises(ValueError, match="exactly one"):
        resolve_air_streams(
            solution,
            spec,
            base_split(),
            total_air_mass_flow_kg_s=4.0,
            overall_equivalence_ratio=0.35,
            **kwargs,
        )


def test_air_state_density_follows_ideal_gas():
    spec = lng_spec()
    solution = registry().new_solution(spec.fuel, spec.role)
    streams = resolve_air_streams(
        solution,
        spec,
        base_split(),
        fuel=FuelKind.LNG,
        fuel_mass_flow_kg_s=0.08,
        temperature_k=750.0,
        pressure_pa=2.0e6,
        total_air_mass_flow_kg_s=4.0,
    )
    state = streams.state
    expected = state.pressure_pa * state.mean_molecular_weight_kg_kmol / (
        8314.462618 * state.temperature_k
    )
    assert state.density_kg_m3 == pytest.approx(expected, rel=1.0e-6)
    assert math.isfinite(state.density_kg_m3)
