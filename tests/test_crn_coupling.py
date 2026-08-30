"""Verification of the operator-split droplet/network coupling."""

from __future__ import annotations

from pathlib import Path

import cantera as ct
import pytest

from fuelnozzle.crn.chemistry import FuelKind
from fuelnozzle.crn.coupling import solve_coupled
from fuelnozzle.crn.liquids import JetALiquidProvider
from fuelnozzle.crn.reactors import InletSpec, ReactorKind, ReactorSpec
from fuelnozzle.crn.spray_source import DropletClass, InitialSizePolicy, SprayBoundary
from fuelnozzle.jet_a import JetAPropertyTable

MECH = str(Path(__file__).resolve().parents[1] / "mech" / "A2NOx_skeletal.yaml")
PRESSURE_PA = 10.0 * ct.one_atm
AIR_FLOW = 1.0
FUEL_FLOW = 0.035


def mechanism() -> ct.Solution:
    return ct.Solution(MECH)


def liquid_provider() -> JetALiquidProvider:
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
            source="illustrative Jet-A properties for verification only",
        )
    )


def spray(radius_m: float = 15.0e-6, velocity_m_s: float = 40.0) -> SprayBoundary:
    return SprayBoundary(
        fuel=FuelKind.JET_A,
        total_fuel_mass_flow_kg_s=FUEL_FLOW,
        vapor_mass_flow_kg_s=0.0,
        vapor_temperature_k=300.0,
        droplet_classes=(
            DropletClass(
                radius_m=radius_m,
                temperature_k=300.0,
                velocity_m_s=velocity_m_s,
                mass_flow_kg_s=FUEL_FLOW,
                number_rate_per_s=1.0,
                origin="test",
            ),
        ),
        injection_velocity_m_s=velocity_m_s,
        cone_angle_deg=90.0,
        apply_aerodynamic_breakup=False,
        size_policy=InitialSizePolicy.USER,
        calibration_id=None,
        warnings=(),
    )


def reactors(spray_path_m: float = 0.03) -> list[ReactorSpec]:
    return [
        ReactorSpec(
            name="evap", kind=ReactorKind.EVAPORATOR, volume_m3=0.6e-3,
            spray_path_length_m=spray_path_m,
        ),
        ReactorSpec(
            name="mixer", kind=ReactorKind.MIXER, volume_m3=0.8e-3,
            spray_path_length_m=spray_path_m,
        ),
        ReactorSpec(name="flame", kind=ReactorKind.PSR, volume_m3=1.5e-3),
        ReactorSpec(name="post", kind=ReactorKind.PSR, volume_m3=4.0e-3),
    ]


AIR = [
    InletSpec(
        name="air", target_reactor="evap", mass_flow_kg_s=AIR_FLOW,
        temperature_k=750.0, mole_fractions={"O2": 0.21, "N2": 0.79},
    )
]
FLOWS = {
    ("evap", "mixer"): AIR_FLOW + FUEL_FLOW + 0.3,
    ("mixer", "flame"): AIR_FLOW + FUEL_FLOW + 0.3,
    ("flame", "post"): AIR_FLOW + FUEL_FLOW,
    ("flame", "mixer"): 0.3,
    ("mixer", "evap"): 0.3,
}


def run(**overrides):
    settings = dict(
        spray=spray(),
        spray_path=("evap", "mixer", "flame"),
        reactor_list=reactors(),
    )
    settings.update(overrides)
    return solve_coupled(
        settings["reactor_list"],
        AIR,
        "post",
        FLOWS,
        mechanism,
        PRESSURE_PA,
        settings["spray"],
        liquid_provider(),
        settings["spray_path"],
        "POSF10325",
        max_iterations=12,
    )


def test_coupled_solution_converges():
    result = run()
    assert result.converged
    assert result.iterations <= 12


def test_small_droplets_fully_evaporate_in_a_hot_combustor():
    result = run()
    assert result.evaporated_fraction == pytest.approx(1.0, abs=1.0e-3)
    assert result.liquid_carryover_kg_s < 1.0e-9


def test_evaporation_happens_where_the_droplets_are():
    """Fuel must be released along the spray path, not everywhere at once."""
    result = run()
    released = {zone.reactor_name: zone.evaporated_mass_flow_kg_s for zone in result.zones}
    assert set(released) == {"evap", "mixer", "flame"}
    assert sum(released.values()) == pytest.approx(FUEL_FLOW, rel=1.0e-3)
    assert released["evap"] > 0.0


def test_evaporation_draws_heat_from_the_gas():
    """The latent heat has to come from somewhere; it comes from the gas."""
    result = run()
    total_heat = sum(zone.heat_drawn_from_gas_w for zone in result.zones)
    assert total_heat > 0.0

    liquid = liquid_provider().liquid_state(400.0, PRESSURE_PA)
    # At minimum, the latent heat of everything that vaporized.
    assert total_heat > 0.5 * FUEL_FLOW * liquid.latent_heat_j_kg


def test_combustion_proceeds_and_produces_nox():
    result = run()
    assert result.network.peak_temperature_k > 1500.0
    assert result.network.outlet.mole_fractions["NO"] > 0.0


def test_no_accumulates_downstream_in_the_coupled_solution():
    result = run()
    evap = result.network.by_name("evap").mole_fractions["NO"]
    post = result.network.by_name("post").mole_fractions["NO"]
    assert post > evap


def test_large_droplets_leave_unevaporated_fuel_and_are_flagged():
    """Liquid reaching the end of the spray path is a design failure, not a detail."""
    result = run(spray=spray(radius_m=400.0e-6, velocity_m_s=120.0))

    assert result.liquid_carryover_kg_s > 0.0
    assert result.evaporated_fraction < 1.0
    codes = {warning.code for warning in result.warnings}
    assert "LIQUID_CARRYOVER" in codes


def test_larger_droplets_evaporate_less_completely():
    small = run(spray=spray(radius_m=15.0e-6))
    large = run(spray=spray(radius_m=200.0e-6, velocity_m_s=120.0))
    assert large.evaporated_fraction < small.evaporated_fraction


def test_spray_path_must_reference_real_reactors():
    with pytest.raises(ValueError, match="unknown reactor"):
        run(spray_path=("evap", "ghost"))


def test_empty_spray_path_is_rejected():
    with pytest.raises(ValueError, match="at least one reactor"):
        run(spray_path=())


def test_prevaporized_fuel_bypasses_the_droplet_path():
    """LNG arriving already flashed should need no evaporation modeling."""
    prevaporized = SprayBoundary(
        fuel=FuelKind.JET_A,
        total_fuel_mass_flow_kg_s=FUEL_FLOW,
        vapor_mass_flow_kg_s=FUEL_FLOW,
        vapor_temperature_k=400.0,
        droplet_classes=(),
        injection_velocity_m_s=40.0,
        cone_angle_deg=None,
        apply_aerodynamic_breakup=False,
        size_policy=InitialSizePolicy.USER,
        calibration_id=None,
        warnings=(),
    )
    result = run(spray=prevaporized)

    assert result.liquid_carryover_kg_s == 0.0
    assert result.network.peak_temperature_k > 1500.0
