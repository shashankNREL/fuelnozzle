"""Liquid fuel property adapters behind one protocol.

The droplet solver must not know which fuel it is integrating. Jet-A properties come
from a measured table supplied by the user; LNG properties come from CoolProp through
the same provider the nozzle tiers already use, so a study cannot end up with the nozzle
and the combustor disagreeing about the fuel.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, log

from fuelnozzle.crn.droplets import LiquidState
from fuelnozzle.jet_a import JetAProperties, JetAPropertyTable
from fuelnozzle.properties import CoolPropLNGProvider, PropertyCalculationError

#: Vapor-air diffusivity constants for Jet-A, John et al. Eq. (18).
JET_A_DIFFUSIVITY_REFERENCE_M2_S = 4.16e-6
JET_A_DIFFUSIVITY_EXPONENT = 1.6

#: Methane-air equivalents. Methane is far more mobile than a kerosene vapor, so its
#: reference diffusivity is roughly five times larger.
METHANE_DIFFUSIVITY_REFERENCE_M2_S = 2.10e-5
METHANE_DIFFUSIVITY_EXPONENT = 1.75

#: Universal gas constant, for the Clausius-Clapeyron fallback.
_CLAUSIUS_CLAPEYRON_GAS_CONSTANT = 8.314462618

#: Reference pressure for the normal boiling point.
ATMOSPHERIC_PRESSURE_PA = 101_325.0


class LiquidPropertyError(RuntimeError):
    """A liquid property needed by the droplet models is unavailable."""


@dataclass(frozen=True)
class JetALiquidProvider:
    """Jet-A liquid properties from a measured table.

    The table must carry the optional droplet columns. They are optional on
    :class:`~fuelnozzle.jet_a.JetAPropertyTable` because the hydraulic nozzle
    calculation does not need them, but the droplet models cannot proceed without
    them and say so rather than substituting a guess.
    """

    table: JetAPropertyTable

    @property
    def vapor_diffusivity_reference_m2_s(self) -> float:
        return JET_A_DIFFUSIVITY_REFERENCE_M2_S

    @property
    def vapor_diffusivity_exponent(self) -> float:
        return JET_A_DIFFUSIVITY_EXPONENT

    def liquid_state(self, temperature_k: float, pressure_pa: float) -> LiquidState:
        properties, _ = self.table.at_temperature(temperature_k)
        self._require(properties)
        assert properties.liquid_cp_j_kg_k is not None
        assert properties.latent_heat_j_kg is not None
        assert properties.molecular_weight_kg_mol is not None
        assert properties.boiling_point_k is not None

        vapor_pressure = self._vapor_pressure(properties, temperature_k, pressure_pa)
        return LiquidState(
            density_kg_m3=properties.density_kg_m3,
            viscosity_pa_s=properties.viscosity_pa_s,
            surface_tension_n_m=properties.surface_tension_n_m,
            specific_heat_j_kg_k=properties.liquid_cp_j_kg_k,
            latent_heat_j_kg=properties.latent_heat_j_kg,
            vapor_pressure_pa=vapor_pressure,
            molecular_weight_kg_mol=properties.molecular_weight_kg_mol,
            saturation_temperature_k=self._saturation_temperature(properties, pressure_pa),
        )

    @staticmethod
    def _require(properties: JetAProperties) -> None:
        missing = [
            name
            for name, value in (
                ("liquid_cp_j_kg_k", properties.liquid_cp_j_kg_k),
                ("latent_heat_j_kg", properties.latent_heat_j_kg),
                ("molecular_weight_kg_mol", properties.molecular_weight_kg_mol),
                ("boiling_point_k", properties.boiling_point_k),
            )
            if value is None
        ]
        if missing:
            raise LiquidPropertyError(
                "Jet-A droplet models require "
                f"{', '.join(missing)} on the property table. These are optional for the "
                "hydraulic nozzle calculation but cannot be guessed for evaporation."
            )

    def _vapor_pressure(
        self, properties: JetAProperties, temperature_k: float, pressure_pa: float
    ) -> float:
        """Vapor pressure, from the table when present, otherwise Clausius-Clapeyron.

        The fallback is anchored on the measured normal boiling point, so it is a
        one-parameter extrapolation rather than a fitted correlation. It is adequate
        for screening and is documented as such.
        """
        if properties.vapor_pressure_pa is not None:
            return properties.vapor_pressure_pa
        assert properties.latent_heat_j_kg is not None
        assert properties.molecular_weight_kg_mol is not None
        assert properties.boiling_point_k is not None
        specific_gas_constant = (
            _CLAUSIUS_CLAPEYRON_GAS_CONSTANT / properties.molecular_weight_kg_mol
        )
        exponent = (properties.latent_heat_j_kg / specific_gas_constant) * (
            1.0 / properties.boiling_point_k - 1.0 / temperature_k
        )
        return ATMOSPHERIC_PRESSURE_PA * exp(max(-700.0, min(700.0, exponent)))

    def _saturation_temperature(
        self, properties: JetAProperties, pressure_pa: float
    ) -> float:
        """Boiling temperature at the local pressure, by inverting Clausius-Clapeyron."""
        assert properties.latent_heat_j_kg is not None
        assert properties.molecular_weight_kg_mol is not None
        assert properties.boiling_point_k is not None
        specific_gas_constant = (
            _CLAUSIUS_CLAPEYRON_GAS_CONSTANT / properties.molecular_weight_kg_mol
        )
        inverse = 1.0 / properties.boiling_point_k - specific_gas_constant * log(
            pressure_pa / ATMOSPHERIC_PRESSURE_PA
        ) / properties.latent_heat_j_kg
        return 1.0 / inverse if inverse > 0.0 else properties.boiling_point_k


@dataclass(frozen=True)
class LNGLiquidProvider:
    """LNG liquid properties from CoolProp, via the existing nozzle property provider.

    Reusing :class:`~fuelnozzle.properties.CoolPropLNGProvider` means the combustor and
    the nozzle tiers share one equation of state and one composition.
    """

    properties: CoolPropLNGProvider

    @property
    def vapor_diffusivity_reference_m2_s(self) -> float:
        return METHANE_DIFFUSIVITY_REFERENCE_M2_S

    @property
    def vapor_diffusivity_exponent(self) -> float:
        return METHANE_DIFFUSIVITY_EXPONENT

    def liquid_state(self, temperature_k: float, pressure_pa: float) -> LiquidState:
        try:
            saturated_liquid = self.properties.bubble_state_at_pressure(pressure_pa)
            saturated_vapor = self.properties.dew_state_at_pressure(pressure_pa)
        except PropertyCalculationError as error:
            raise LiquidPropertyError(
                f"CoolProp could not evaluate LNG saturation at {pressure_pa:.4g} Pa: {error}"
            ) from error

        latent_heat = saturated_vapor.enthalpy_j_kg - saturated_liquid.enthalpy_j_kg
        if latent_heat <= 0.0:
            raise LiquidPropertyError(
                "Non-positive LNG latent heat; the state is at or above the critical point "
                "and a droplet model does not apply."
            )

        # Below saturation the liquid is compressed; at or above it the droplet is
        # boiling and the saturated liquid is the correct state to report.
        compressed = None
        if temperature_k < saturated_liquid.temperature_k:
            try:
                compressed = self.properties.state_pt(pressure_pa, temperature_k)
            except PropertyCalculationError:
                compressed = None

        # Saturated liquid at the droplet's own temperature. Surface tension is defined
        # only along the saturation line, so a compressed-liquid state cannot supply it
        # and this is where it must come from.
        saturated_at_droplet = None
        try:
            saturated_at_droplet = self.properties.state_tq(temperature_k, 0.0)
        except PropertyCalculationError:
            saturated_at_droplet = None

        candidates = tuple(
            state
            for state in (compressed, saturated_at_droplet, saturated_liquid)
            if state is not None
        )

        try:
            vapor_pressure = self.properties.bubble_pressure_at_temperature(temperature_k)
        except PropertyCalculationError:
            vapor_pressure = pressure_pa

        density = candidates[0].density_kg_m3
        return LiquidState(
            density_kg_m3=density,
            viscosity_pa_s=self._first_available(candidates, "viscosity_pa_s", "viscosity"),
            surface_tension_n_m=self._first_available(
                candidates, "surface_tension_n_m", "surface tension"
            ),
            specific_heat_j_kg_k=self._first_available(
                candidates, "cp_j_kg_k", "specific heat"
            ),
            latent_heat_j_kg=latent_heat,
            vapor_pressure_pa=vapor_pressure,
            molecular_weight_kg_mol=self.properties.molar_mass_kg_mol,
            saturation_temperature_k=saturated_liquid.temperature_k,
        )

    @staticmethod
    def _first_available(states: tuple, attribute: str, label: str) -> float:
        """Take a property from the most representative state that reports it.

        CoolProp leaves some properties undefined for a compressed-liquid state --
        surface tension exists only along the saturation line -- so the search falls
        back to the saturated liquid at the droplet temperature and then at the local
        pressure, rather than failing on a property that is genuinely available nearby.
        """
        for state in states:
            value = getattr(state, attribute, None)
            if value is not None and value > 0.0:
                return float(value)
        raise LiquidPropertyError(
            f"CoolProp reported no {label} for any LNG liquid state near this condition; "
            "the droplet model cannot proceed without it."
        )
