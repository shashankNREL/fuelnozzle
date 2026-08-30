"""Verification of droplet breakup, evaporation, and heating.

The quantitative anchors come from John et al. (2026) Fig. 3 and Section 3.1. That case
is underspecified in the paper -- neither the ambient pressure nor the Jet-A liquid
properties are stated -- so these tests assert the paper's *conclusion* and *trend*
rather than its absolute numbers. See the implementation log, deviation D-011.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from fuelnozzle.crn.droplets import (
    EvaporationRegime,
    GasState,
    LiquidState,
    droplet_rates,
    gas_weber_number,
    integrate_droplet,
    reference_temperature_k,
    spalding_mass_transfer_number,
    surface_vapor_mass_fraction,
    taylor_analogy_breakup,
)
from fuelnozzle.crn.liquids import JetALiquidProvider, LiquidPropertyError, LNGLiquidProvider
from fuelnozzle.jet_a import JetAPropertyTable
from fuelnozzle.models import LNGComposition
from fuelnozzle.properties import CoolPropLNGProvider

ATM = 101_325.0


def jet_a_table(**overrides) -> JetAPropertyTable:
    values = dict(
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
    values.update(overrides)
    return JetAPropertyTable(**values)


def air(temperature_k: float, pressure_pa: float = ATM, fuel_vapor: float = 0.0) -> GasState:
    return GasState(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        density_kg_m3=pressure_pa / (287.05 * temperature_k),
        viscosity_pa_s=2.08e-5,
        conductivity_w_m_k=0.030,
        specific_heat_j_kg_k=1009.0,
        mean_molecular_weight_kg_mol=0.02897,
        fuel_vapor_mass_fraction=fuel_vapor,
    )


@dataclass(frozen=True)
class ConstantLiquidProvider:
    """Fixed properties, so an analytic result can be compared without interpolation."""

    state: LiquidState

    @property
    def vapor_diffusivity_reference_m2_s(self) -> float:
        return 4.16e-6

    @property
    def vapor_diffusivity_exponent(self) -> float:
        return 1.6

    def liquid_state(self, temperature_k: float, pressure_pa: float) -> LiquidState:
        return self.state


def constant_liquid(**overrides) -> ConstantLiquidProvider:
    values = dict(
        density_kg_m3=800.0,
        viscosity_pa_s=1.5e-3,
        surface_tension_n_m=0.025,
        specific_heat_j_kg_k=2000.0,
        latent_heat_j_kg=3.5e5,
        vapor_pressure_pa=2.0e3,
        molecular_weight_kg_mol=0.15429,
        saturation_temperature_k=470.0,
    )
    values.update(overrides)
    return ConstantLiquidProvider(LiquidState(**values))


# --- Film temperature and dimensionless groups ---------------------------------------


def test_reference_temperature_is_weighted_toward_the_droplet():
    assert reference_temperature_k(900.0, 300.0) == pytest.approx(500.0)


def test_weber_number_matches_its_definition():
    gas = air(350.0)
    liquid = constant_liquid().state
    weber = gas_weber_number(gas, liquid, 75.0e-6, 112.0)
    expected = gas.density_kg_m3 * 112.0**2 * 75.0e-6 / liquid.surface_tension_n_m
    assert weber == pytest.approx(expected)
    assert weber == pytest.approx(37.95, rel=0.01)


def test_surface_vapor_fraction_matches_hand_calculation():
    """John et al. Eq. (21), evaluated by hand."""
    gas = air(350.0)
    liquid = constant_liquid(vapor_pressure_pa=2.0e3).state
    fraction = surface_vapor_mass_fraction(gas, liquid)
    ratio = ATM / 2.0e3
    expected = 0.15429 / (0.15429 + 0.02897 * (ratio - 1.0))
    assert fraction == pytest.approx(expected)


def test_surface_vapor_fraction_saturates_at_the_boiling_point():
    """When vapor pressure reaches ambient, the surface is pure fuel vapor."""
    gas = air(350.0)
    liquid = constant_liquid(vapor_pressure_pa=ATM).state
    assert surface_vapor_mass_fraction(gas, liquid) == 1.0
    assert spalding_mass_transfer_number(1.0, 0.0) == float("inf")


def test_transfer_number_vanishes_in_saturated_surroundings():
    gas = air(350.0, fuel_vapor=0.05)
    liquid = constant_liquid(vapor_pressure_pa=2.0e3).state
    surface = surface_vapor_mass_fraction(gas, liquid)
    assert spalding_mass_transfer_number(surface, surface) == pytest.approx(0.0)


# --- Breakup -------------------------------------------------------------------------


def test_breakup_reproduces_the_paper_order_of_magnitude():
    """Fig. 3 case: 150 um, 300 K droplet, 350 K gas, 112 m/s.

    The paper reports breakup at about 2e-5 s. Ambient pressure and liquid properties
    are not stated, so agreement is asserted only to an order of magnitude.
    """
    provider = JetALiquidProvider(jet_a_table())
    result = taylor_analogy_breakup(air(350.0), provider, 75.0e-6, 300.0, 112.0)

    assert result.occurred
    assert 1.0e-5 <= result.elapsed_time_s <= 1.0e-4
    assert result.final_radius_m < result.initial_radius_m / 5.0


def test_breakup_is_stronger_at_elevated_pressure():
    """Weber number scales with gas density, so pressure drives the cascade further."""
    provider = JetALiquidProvider(jet_a_table())
    low = taylor_analogy_breakup(air(350.0, ATM), provider, 75.0e-6, 300.0, 112.0)
    high = taylor_analogy_breakup(air(350.0, 20.0 * ATM), provider, 75.0e-6, 300.0, 112.0)

    assert high.weber_number == pytest.approx(20.0 * low.weber_number, rel=1.0e-6)
    assert high.final_radius_m < low.final_radius_m
    assert high.elapsed_time_s < low.elapsed_time_s


def test_no_breakup_without_relative_velocity():
    provider = JetALiquidProvider(jet_a_table())
    result = taylor_analogy_breakup(air(350.0), provider, 75.0e-6, 300.0, 0.0)
    assert not result.occurred
    assert result.final_radius_m == pytest.approx(75.0e-6)


def test_no_breakup_below_the_critical_weber_number():
    """A slow, small droplet is held together by surface tension."""
    provider = JetALiquidProvider(jet_a_table())
    result = taylor_analogy_breakup(air(350.0), provider, 5.0e-6, 300.0, 2.0)
    assert result.weber_number < 1.0
    assert not result.occurred


def test_breakup_cascade_terminates():
    provider = JetALiquidProvider(jet_a_table())
    result = taylor_analogy_breakup(air(350.0, 40.0 * ATM), provider, 200.0e-6, 300.0, 200.0)
    assert result.final_radius_m > 0.0
    assert result.stages >= 1


# --- Evaporation: the two branches ---------------------------------------------------


def test_d_squared_law_is_recovered_in_the_quiescent_limit():
    """With no relative velocity and fixed properties, d^2 must fall linearly in time.

    This is the classical result every droplet evaporation model must reproduce.
    """
    provider = constant_liquid()
    history = integrate_droplet(
        air(800.0), provider, 40.0e-6, 400.0, relative_velocity_m_s=0.0,
        residence_time_s=0.02,
    )
    times = np.array(history.time_s)
    temperatures = np.array(history.temperature_k)
    diameters_squared = (2.0 * np.array(history.radius_m)) ** 2

    # The D^2 law applies once the droplet has reached its wet-bulb temperature. Before
    # that, the film temperature is still moving and the vaporization rate with it, so
    # the early transient is legitimately non-linear and is excluded.
    settled = (np.abs(temperatures - temperatures[-1]) < 1.0) & (
        diameters_squared > 0.05 * diameters_squared[0]
    )
    assert settled.sum() > 10, "not enough settled points to test the D^2 law"

    slope, intercept = np.polyfit(times[settled], diameters_squared[settled], 1)
    predicted = slope * times[settled] + intercept
    residual = np.max(np.abs(diameters_squared[settled] - predicted))

    assert slope < 0.0
    assert residual < 0.01 * diameters_squared[0]


def test_boiling_branch_activates_above_the_saturation_temperature():
    provider = constant_liquid(saturation_temperature_k=380.0)
    rates = droplet_rates(air(900.0), provider, 20.0e-6, 400.0, 0.0)
    assert rates.regime is EvaporationRegime.BOILING


def test_boiling_branch_pins_the_droplet_temperature():
    """A boiling droplet cannot get hotter; all arriving heat drives the phase change."""
    provider = constant_liquid(saturation_temperature_k=380.0)
    rates = droplet_rates(air(900.0), provider, 20.0e-6, 400.0, 0.0)
    assert rates.temperature_rate_k_s == 0.0


def test_boiling_rate_equals_heat_divided_by_latent_heat():
    """The defining statement of the heat-transfer-limited branch, checked exactly."""
    provider = constant_liquid(saturation_temperature_k=380.0, latent_heat_j_kg=5.0e5)
    rates = droplet_rates(air(900.0), provider, 20.0e-6, 400.0, 0.0)
    assert rates.mass_rate_kg_s == pytest.approx(-rates.convective_heat_w / 5.0e5)
    assert rates.convective_heat_w > 0.0


def test_boiling_stops_when_the_gas_is_colder_than_saturation():
    provider = constant_liquid(saturation_temperature_k=380.0)
    rates = droplet_rates(air(370.0), provider, 20.0e-6, 400.0, 0.0)
    assert rates.regime is EvaporationRegime.INERT
    assert rates.mass_rate_kg_s == 0.0


def test_diffusion_branch_used_below_saturation():
    provider = constant_liquid(saturation_temperature_k=470.0, vapor_pressure_pa=2.0e3)
    rates = droplet_rates(air(600.0), provider, 20.0e-6, 350.0, 0.0)
    assert rates.regime is EvaporationRegime.DIFFUSION
    assert rates.mass_rate_kg_s < 0.0


def test_hotter_gas_evaporates_faster():
    provider = constant_liquid()
    cool = droplet_rates(air(500.0), provider, 20.0e-6, 350.0, 0.0)
    hot = droplet_rates(air(800.0), provider, 20.0e-6, 350.0, 0.0)
    assert abs(hot.mass_rate_kg_s) > abs(cool.mass_rate_kg_s)


def test_relative_velocity_enhances_evaporation():
    provider = constant_liquid()
    still = droplet_rates(air(700.0), provider, 20.0e-6, 350.0, 0.0)
    moving = droplet_rates(air(700.0), provider, 20.0e-6, 350.0, 50.0)
    assert abs(moving.mass_rate_kg_s) > abs(still.mass_rate_kg_s)


# --- Energy accounting (plan Section 3.6) --------------------------------------------


def test_droplet_energy_split_closes():
    """Arriving heat is split between the phase change and warming the liquid.

    Latent heat must be counted once. If the sensible term were computed from the
    convective heat alone, this balance would fail.
    """
    provider = constant_liquid()
    radius, temperature = 20.0e-6, 350.0
    rates = droplet_rates(air(700.0), provider, radius, temperature, 0.0)

    liquid = provider.liquid_state(temperature, ATM)
    mass = (4.0 / 3.0) * np.pi * radius**3 * liquid.density_kg_m3
    sensible = liquid.specific_heat_j_kg_k * mass * rates.temperature_rate_k_s
    latent = -rates.mass_rate_kg_s * liquid.latent_heat_j_kg

    assert sensible + latent == pytest.approx(rates.convective_heat_w, rel=1.0e-10)


def test_heat_transfer_scaling_is_linear_in_the_boiling_branch():
    provider = constant_liquid(saturation_temperature_k=380.0)
    base = droplet_rates(air(900.0), provider, 20.0e-6, 400.0, 0.0, heat_transfer_scaling=1.0)
    scaled = droplet_rates(air(900.0), provider, 20.0e-6, 400.0, 0.0, heat_transfer_scaling=2.0)
    assert scaled.convective_heat_w == pytest.approx(2.0 * base.convective_heat_w)


# --- Timescale separation: the paper's actual conclusion ------------------------------


def _time_to_one_tenth_diameter(gas_temperature_k: float) -> float:
    provider = JetALiquidProvider(jet_a_table())
    history = integrate_droplet(
        air(gas_temperature_k), provider, 75.0e-6, 300.0, 112.0, residence_time_s=5.0
    )
    radii = np.array(history.radius_m)
    times = np.array(history.time_s)
    reached = radii <= 7.5e-6
    assert reached.any(), "droplet never reached one-tenth of its initial size"
    return float(times[np.argmax(reached)])


def test_breakup_is_far_faster_than_evaporation_at_350_k():
    """The paper's conclusion: separation exceeding two orders of magnitude.

    The paper's own ratio is 133. We obtain a much larger separation because our
    evaporation is slower, but the conclusion that justifies treating breakup and
    evaporation sequentially is reproduced, which is what the implementation depends on.
    """
    provider = JetALiquidProvider(jet_a_table())
    breakup = taylor_analogy_breakup(air(350.0), provider, 75.0e-6, 300.0, 112.0)
    ratio = _time_to_one_tenth_diameter(350.0) / breakup.elapsed_time_s
    assert ratio > 100.0


def test_separation_survives_at_900_k_but_narrows():
    """Also the paper's conclusion: still above one order of magnitude at 900 K."""
    provider = JetALiquidProvider(jet_a_table())
    breakup = taylor_analogy_breakup(air(900.0), provider, 75.0e-6, 300.0, 112.0)
    ratio = _time_to_one_tenth_diameter(900.0) / breakup.elapsed_time_s
    assert ratio > 10.0


def test_timescale_separation_narrows_with_temperature():
    """Direction of the paper's 133 -> 21 trend: evaporation speeds up more than breakup."""
    provider = JetALiquidProvider(jet_a_table())
    cool = _time_to_one_tenth_diameter(350.0) / taylor_analogy_breakup(
        air(350.0), provider, 75.0e-6, 300.0, 112.0
    ).elapsed_time_s
    hot = _time_to_one_tenth_diameter(900.0) / taylor_analogy_breakup(
        air(900.0), provider, 75.0e-6, 300.0, 112.0
    ).elapsed_time_s
    assert hot < cool


# --- Integration behaviour -----------------------------------------------------------


def test_integration_terminates_on_full_evaporation():
    provider = constant_liquid()
    history = integrate_droplet(air(1200.0), provider, 10.0e-6, 400.0, 0.0, residence_time_s=1.0)
    assert history.fully_evaporated
    assert history.evaporated_mass_fraction == pytest.approx(1.0, abs=1.0e-6)
    assert history.time_s[-1] < 1.0


def test_partial_evaporation_reports_a_fraction_between_zero_and_one():
    provider = constant_liquid()
    history = integrate_droplet(
        air(600.0), provider, 50.0e-6, 350.0, 0.0, residence_time_s=1.0e-4
    )
    assert not history.fully_evaporated
    assert 0.0 < history.evaporated_mass_fraction < 1.0


def test_integration_rejects_nonpositive_inputs():
    provider = constant_liquid()
    with pytest.raises(ValueError, match="radius must be positive"):
        integrate_droplet(air(600.0), provider, 0.0, 350.0, 0.0, residence_time_s=1.0e-3)
    with pytest.raises(ValueError, match="Residence time must be positive"):
        integrate_droplet(air(600.0), provider, 1.0e-5, 350.0, 0.0, residence_time_s=0.0)


# --- Liquid property adapters ---------------------------------------------------------


def test_jet_a_provider_requires_the_droplet_columns():
    """The hydraulic table is valid without them; the droplet model must not guess."""
    table = JetAPropertyTable(
        temperature_k=(280.0, 300.0),
        density_kg_m3=(815.0, 800.0),
        viscosity_pa_s=(2.1e-3, 1.5e-3),
        surface_tension_n_m=(0.027, 0.0255),
        source="hydraulic-only table",
    )
    with pytest.raises(LiquidPropertyError, match="liquid_cp_j_kg_k"):
        JetALiquidProvider(table).liquid_state(300.0, ATM)


def test_jet_a_saturation_temperature_rises_with_pressure():
    provider = JetALiquidProvider(jet_a_table())
    low = provider.liquid_state(300.0, ATM).saturation_temperature_k
    high = provider.liquid_state(300.0, 20.0 * ATM).saturation_temperature_k
    assert low == pytest.approx(470.0, rel=1.0e-6)
    assert high > low


def test_lng_provider_supplies_a_usable_liquid_state():
    coolprop = CoolPropLNGProvider(LNGComposition.pure_methane())
    provider = LNGLiquidProvider(coolprop)
    state = provider.liquid_state(120.0, 2.0e5)

    assert state.density_kg_m3 > 300.0
    assert state.latent_heat_j_kg > 1.0e5
    assert state.molecular_weight_kg_mol == pytest.approx(0.016043, rel=1.0e-3)
    assert 100.0 < state.saturation_temperature_k < 160.0


def test_lng_droplet_above_saturation_uses_the_boiling_branch():
    """A flashing LNG droplet is not diffusion limited, which is why the branch exists."""
    coolprop = CoolPropLNGProvider(LNGComposition.pure_methane())
    provider = LNGLiquidProvider(coolprop)
    saturation = provider.liquid_state(120.0, 2.0e5).saturation_temperature_k

    gas = GasState(
        temperature_k=800.0,
        pressure_pa=2.0e5,
        density_kg_m3=2.0e5 / (287.05 * 800.0),
        viscosity_pa_s=3.5e-5,
        conductivity_w_m_k=0.055,
        specific_heat_j_kg_k=1100.0,
        mean_molecular_weight_kg_mol=0.02897,
    )
    rates = droplet_rates(gas, provider, 10.0e-6, saturation + 1.0, 0.0)
    assert rates.regime is EvaporationRegime.BOILING
    assert rates.temperature_rate_k_s == 0.0
    assert rates.mass_rate_kg_s < 0.0
