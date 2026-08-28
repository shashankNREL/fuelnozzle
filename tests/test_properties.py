import pytest

from fuelnozzle.models import LNGComposition
from fuelnozzle.properties import CoolPropLNGProvider


def test_pure_methane_direct_flashes_round_trip() -> None:
    provider = CoolPropLNGProvider(LNGComposition.pure_methane())
    inlet = provider.state_pt(2.0e6, 120.0)

    ph_state = provider.state_ph(8.0e5, inlet.enthalpy_j_kg)
    ps_state = provider.state_ps(8.0e5, inlet.entropy_j_kg_k)

    assert ph_state.enthalpy_j_kg == pytest.approx(inlet.enthalpy_j_kg, rel=1.0e-9)
    assert ps_state.entropy_j_kg_k == pytest.approx(inlet.entropy_j_kg_k, rel=1.0e-9)
    assert provider.coolprop_version.startswith("8.")


def test_mixture_pt_roots_round_trip_and_convert_quality() -> None:
    composition = LNGComposition(
        mole_fractions={"Methane": 0.90, "Ethane": 0.07, "Nitrogen": 0.03}
    )
    provider = CoolPropLNGProvider(composition)
    cold_liquid = provider.state_pt(2.3e6, 120.0)
    bubble = provider.bubble_state_at_pressure(5.0e5)
    dew = provider.dew_state_at_pressure(5.0e5)
    two_phase = provider.state_pt(5.0e5, 0.5 * (bubble.temperature_k + dew.temperature_k))

    ph_state = provider.state_ph(
        5.0e5,
        two_phase.enthalpy_j_kg,
        temperature_hint_k=two_phase.temperature_k,
    )
    ps_state = provider.state_ps(
        5.0e5,
        two_phase.entropy_j_kg_k,
        temperature_hint_k=two_phase.temperature_k,
    )

    assert bubble.temperature_k < dew.temperature_k
    assert two_phase.phase == "two_phase"
    assert two_phase.vapor_quality_mass is not None
    assert 0.0 < two_phase.vapor_quality_mass < 1.0
    assert cold_liquid.viscosity_pa_s is not None
    assert cold_liquid.conductivity_w_m_k is not None
    assert cold_liquid.transport_model == "component_mixing_fallback"
    assert ph_state.temperature_k == pytest.approx(two_phase.temperature_k, abs=1.0e-5)
    assert ps_state.temperature_k == pytest.approx(two_phase.temperature_k, abs=1.0e-5)


def test_composition_is_normalized() -> None:
    composition = LNGComposition(mole_fractions={"Methane": 90.0, "Ethane": 10.0})
    assert sum(composition.fractions) == pytest.approx(1.0)
