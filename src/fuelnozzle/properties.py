"""CoolProp-backed thermodynamic properties for methane and compositional LNG."""

from __future__ import annotations

from collections.abc import Callable
from math import exp, isfinite, log

import CoolProp
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import brentq

from fuelnozzle.models import LNGComposition, ThermodynamicState


class PropertyCalculationError(RuntimeError):
    """Raised when a requested thermodynamic state cannot be established."""


class CoolPropLNGProvider:
    """Reusable CoolProp state for pure methane or a specified LNG mixture.

    Pure fluids use CoolProp's direct mass-enthalpy and mass-entropy flashes.
    CoolProp mixtures currently support PT, PQ, and TQ inputs, so PH and PS
    states are obtained with an outer, bracketed temperature root around PT
    flashes.
    """

    def __init__(self, composition: LNGComposition, backend: str = "HEOS") -> None:
        self.composition = composition
        self.backend = backend
        fluid_string = "&".join(composition.components)
        try:
            self._state = CP.AbstractState(backend, fluid_string)
            if not composition.is_pure:
                self._state.set_mole_fractions(list(composition.fractions))
        except (RuntimeError, ValueError) as exc:
            raise PropertyCalculationError(
                f"CoolProp could not initialize LNG components {composition.components}: {exc}"
            ) from exc

        self._component_molar_masses = tuple(
            CP.PropsSI("M", component) for component in composition.components
        )
        self._component_states = tuple(
            CP.AbstractState(backend, component) for component in composition.components
        )
        self.transport_fallback_used = False
        self.transport_saturation_fallback_components: set[str] = set()
        self.transport_omitted_components: set[str] = set()

    @property
    def coolprop_version(self) -> str:
        return CoolProp.__version__

    @property
    def fluid_label(self) -> str:
        return "&".join(self.composition.components)

    @property
    def molar_mass_kg_mol(self) -> float:
        """Mole-fraction-weighted molar mass of the LNG mixture.

        Needed by the reactor-network droplet models to convert between the mole basis
        CoolProp reports and the mass basis the evaporation equations use.
        """
        return float(
            sum(
                fraction * molar_mass
                for fraction, molar_mass in zip(
                    self.composition.fractions, self._component_molar_masses, strict=True
                )
            )
        )

    def state_pt(self, pressure_pa: float, temperature_k: float) -> ThermodynamicState:
        self._validate_pressure(pressure_pa)
        if not isfinite(temperature_k) or temperature_k <= 0.0:
            raise ValueError("Temperature must be finite and positive")
        try:
            self._state.update(CP.PT_INPUTS, pressure_pa, temperature_k)
        except (RuntimeError, ValueError) as exc:
            raise PropertyCalculationError(
                f"PT flash failed at P={pressure_pa:g} Pa, T={temperature_k:g} K: {exc}"
            ) from exc
        return self._current_state()

    def state_ph(
        self,
        pressure_pa: float,
        enthalpy_j_kg: float,
        *,
        temperature_hint_k: float | None = None,
    ) -> ThermodynamicState:
        self._validate_pressure(pressure_pa)
        if self.composition.is_pure:
            try:
                self._state.update(CP.HmassP_INPUTS, enthalpy_j_kg, pressure_pa)
            except (RuntimeError, ValueError) as exc:
                raise PropertyCalculationError(
                    f"PH flash failed at P={pressure_pa:g} Pa, h={enthalpy_j_kg:g} J/kg: {exc}"
                ) from exc
            return self._current_state()
        return self._mixture_target_state(
            pressure_pa,
            enthalpy_j_kg,
            self._state.hmass,
            "enthalpy",
            temperature_hint_k,
        )

    def state_ps(
        self,
        pressure_pa: float,
        entropy_j_kg_k: float,
        *,
        temperature_hint_k: float | None = None,
    ) -> ThermodynamicState:
        self._validate_pressure(pressure_pa)
        if self.composition.is_pure:
            try:
                self._state.update(CP.PSmass_INPUTS, pressure_pa, entropy_j_kg_k)
            except (RuntimeError, ValueError) as exc:
                raise PropertyCalculationError(
                    f"PS flash failed at P={pressure_pa:g} Pa, s={entropy_j_kg_k:g} J/kg/K: {exc}"
                ) from exc
            return self._current_state()
        return self._mixture_target_state(
            pressure_pa,
            entropy_j_kg_k,
            self._state.smass,
            "entropy",
            temperature_hint_k,
        )

    def state_pq(self, pressure_pa: float, molar_quality: float) -> ThermodynamicState:
        self._validate_pressure(pressure_pa)
        if not 0.0 <= molar_quality <= 1.0:
            raise ValueError("Molar vapor quality must be between zero and one")
        try:
            self._state.update(CP.PQ_INPUTS, pressure_pa, molar_quality)
        except (RuntimeError, ValueError) as exc:
            raise PropertyCalculationError(
                f"PQ flash failed at P={pressure_pa:g} Pa, Q={molar_quality:g}: {exc}"
            ) from exc
        return self._current_state()

    def state_tq(self, temperature_k: float, molar_quality: float) -> ThermodynamicState:
        if not 0.0 <= molar_quality <= 1.0:
            raise ValueError("Molar vapor quality must be between zero and one")
        try:
            self._state.update(CP.QT_INPUTS, molar_quality, temperature_k)
        except (RuntimeError, ValueError) as exc:
            raise PropertyCalculationError(
                f"TQ flash failed at T={temperature_k:g} K, Q={molar_quality:g}: {exc}"
            ) from exc
        return self._current_state()

    def bubble_state_at_pressure(self, pressure_pa: float) -> ThermodynamicState:
        return self.state_pq(pressure_pa, 0.0)

    def dew_state_at_pressure(self, pressure_pa: float) -> ThermodynamicState:
        return self.state_pq(pressure_pa, 1.0)

    def bubble_pressure_at_temperature(self, temperature_k: float) -> float:
        return self.state_tq(temperature_k, 0.0).pressure_pa

    def _mixture_target_state(
        self,
        pressure_pa: float,
        target: float,
        getter: Callable[[], float],
        quantity_name: str,
        temperature_hint_k: float | None,
    ) -> ThermodynamicState:
        def residual(temperature_k: float) -> float:
            self._state.update(CP.PT_INPUTS, pressure_pa, temperature_k)
            return getter() - target

        roots: list[tuple[float, float]] = []
        valid_points: list[tuple[float, float]] = []
        minimum_k = max(float(self._state.Tmin()) + 1.0e-4, 60.0)
        maximum_k = min(float(self._state.Tmax()) - 1.0e-4, 800.0)

        if temperature_hint_k is not None and minimum_k < temperature_hint_k < maximum_k:
            local_bracket = self._local_temperature_bracket(
                residual,
                temperature_hint_k,
                minimum_k,
                maximum_k,
            )
            if local_bracket is not None:
                roots.append(local_bracket)

        if not roots:
            candidates = np.linspace(minimum_k, maximum_k, 160)
            for temperature_k in candidates:
                try:
                    value = residual(float(temperature_k))
                except (RuntimeError, ValueError):
                    continue
                if isfinite(value):
                    valid_points.append((float(temperature_k), value))

            for left, right in zip(valid_points, valid_points[1:], strict=False):
                if left[1] == 0.0:
                    roots.append((left[0], left[0]))
                elif left[1] * right[1] < 0.0:
                    roots.append((left[0], right[0]))

        if not roots:
            raise PropertyCalculationError(
                f"Could not bracket mixture {quantity_name} target {target:g} at "
                f"P={pressure_pa:g} Pa"
            )

        if temperature_hint_k is None:
            bracket = roots[0]
        else:
            bracket = min(roots, key=lambda pair: abs(sum(pair) / 2.0 - temperature_hint_k))

        try:
            if bracket[0] == bracket[1]:
                temperature_k = bracket[0]
            else:
                temperature_k = brentq(residual, bracket[0], bracket[1], xtol=1.0e-8)
            self._state.update(CP.PT_INPUTS, pressure_pa, temperature_k)
        except (RuntimeError, ValueError) as exc:
            raise PropertyCalculationError(
                f"Mixture {quantity_name} root failed at P={pressure_pa:g} Pa: {exc}"
            ) from exc
        return self._current_state()

    @staticmethod
    def _local_temperature_bracket(
        residual: Callable[[float], float],
        temperature_hint_k: float,
        minimum_k: float,
        maximum_k: float,
    ) -> tuple[float, float] | None:
        points: dict[float, float] = {}
        for span_k in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0):
            for temperature_k in (
                max(minimum_k, temperature_hint_k - span_k),
                min(maximum_k, temperature_hint_k + span_k),
            ):
                if temperature_k in points:
                    continue
                try:
                    value = residual(temperature_k)
                except (RuntimeError, ValueError):
                    continue
                if isfinite(value):
                    points[temperature_k] = value

            ordered = sorted(points.items())
            for left, right in zip(ordered, ordered[1:], strict=False):
                if left[1] == 0.0:
                    return left[0], left[0]
                if left[1] * right[1] < 0.0:
                    return left[0], right[0]
        return None

    def _current_state(self) -> ThermodynamicState:
        phase_index = self._state.phase()
        phase = self._phase_name(phase_index)
        quality = self._mass_quality() if phase_index == CP.iphase_twophase else None
        liquid_composition: dict[str, float] | None = None
        vapor_composition: dict[str, float] | None = None
        if phase_index == CP.iphase_twophase:
            liquid_composition = dict(
                zip(
                    self.composition.components,
                    self._state.mole_fractions_liquid(),
                    strict=True,
                )
            )
            vapor_composition = dict(
                zip(
                    self.composition.components,
                    self._state.mole_fractions_vapor(),
                    strict=True,
                )
            )
        elif phase_index in (
            CP.iphase_gas,
            CP.iphase_supercritical_gas,
        ):
            vapor_composition = dict(self.composition.mole_fractions)
        elif phase_index in (
            CP.iphase_liquid,
            CP.iphase_supercritical_liquid,
        ):
            liquid_composition = dict(self.composition.mole_fractions)
        viscosity = self._optional_property(self._state.viscosity)
        conductivity = self._optional_property(self._state.conductivity)
        surface_tension = self._optional_property(self._state.surface_tension)
        transport_model = "coolprop_native"
        if not self.composition.is_pure and (viscosity is None or conductivity is None):
            fallback_viscosity, fallback_conductivity, fallback_surface_tension = (
                self._component_transport_fallback(phase_index)
            )
            viscosity = viscosity or fallback_viscosity
            conductivity = conductivity or fallback_conductivity
            surface_tension = surface_tension or fallback_surface_tension
            transport_model = "component_mixing_fallback"
            self.transport_fallback_used = True

        return ThermodynamicState(
            pressure_pa=float(self._state.p()),
            temperature_k=float(self._state.T()),
            density_kg_m3=float(self._state.rhomass()),
            enthalpy_j_kg=float(self._state.hmass()),
            entropy_j_kg_k=float(self._state.smass()),
            phase=phase,
            vapor_quality_mass=quality,
            cp_j_kg_k=self._optional_property(self._state.cpmass),
            viscosity_pa_s=viscosity,
            conductivity_w_m_k=conductivity,
            surface_tension_n_m=surface_tension,
            transport_model=transport_model,
            liquid_mole_fractions=liquid_composition,
            vapor_mole_fractions=vapor_composition,
        )

    def _component_transport_fallback(
        self,
        phase_index: int,
    ) -> tuple[float | None, float | None, float | None]:
        molar_quality = float(self._state.Q()) if phase_index == CP.iphase_twophase else -1.0
        use_vapor = phase_index in (CP.iphase_gas, CP.iphase_supercritical_gas) or (
            phase_index == CP.iphase_twophase and molar_quality >= 1.0 - 1.0e-10
        )
        if phase_index == CP.iphase_twophase:
            fractions = (
                self._state.mole_fractions_vapor()
                if use_vapor
                else self._state.mole_fractions_liquid()
            )
        else:
            fractions = list(self.composition.fractions)

        viscosities: list[tuple[float, float]] = []
        conductivities: list[tuple[float, float]] = []
        surface_tensions: list[tuple[float, float]] = []
        for component, fraction, component_state in zip(
            self.composition.components,
            fractions,
            self._component_states,
            strict=True,
        ):
            transport_state_available = False
            try:
                component_state.specify_phase(CP.iphase_gas if use_vapor else CP.iphase_liquid)
                component_state.update(CP.PT_INPUTS, self._state.p(), self._state.T())
                component_state.unspecify_phase()
                transport_state_available = True
            except (RuntimeError, ValueError):
                component_state.unspecify_phase()
                if not use_vapor:
                    try:
                        component_state.update(CP.QT_INPUTS, 0.0, self._state.T())
                        transport_state_available = True
                        self.transport_saturation_fallback_components.add(component)
                    except (RuntimeError, ValueError):
                        self.transport_omitted_components.add(component)

            if transport_state_available:
                component_viscosity = self._optional_property(component_state.viscosity)
                component_conductivity = self._optional_property(component_state.conductivity)
                if component_viscosity is not None:
                    viscosities.append((fraction, component_viscosity))
                if component_conductivity is not None:
                    conductivities.append((fraction, component_conductivity))

            if not use_vapor:
                try:
                    component_state.update(CP.QT_INPUTS, 0.0, self._state.T())
                    component_surface_tension = self._optional_property(
                        component_state.surface_tension
                    )
                    if component_surface_tension is not None:
                        surface_tensions.append((fraction, component_surface_tension))
                except (RuntimeError, ValueError):
                    pass

        viscosity_total = sum(fraction for fraction, _ in viscosities)
        conductivity_total = sum(fraction for fraction, _ in conductivities)
        surface_tension_total = sum(fraction for fraction, _ in surface_tensions)
        viscosity = (
            exp(
                sum(
                    fraction / viscosity_total * log(value)
                    for fraction, value in viscosities
                )
            )
            if viscosity_total > 0.0
            else None
        )
        conductivity = (
            sum(fraction / conductivity_total * value for fraction, value in conductivities)
            if conductivity_total > 0.0
            else None
        )
        surface_tension = (
            sum(
                fraction / surface_tension_total * value
                for fraction, value in surface_tensions
            )
            if surface_tension_total > 0.0
            else None
        )
        return viscosity, conductivity, surface_tension

    def _mass_quality(self) -> float:
        molar_quality = min(1.0, max(0.0, float(self._state.Q())))
        if self.composition.is_pure:
            return molar_quality

        liquid_molar_mass = sum(
            fraction * molar_mass
            for fraction, molar_mass in zip(
                self._state.mole_fractions_liquid(),
                self._component_molar_masses,
                strict=True,
            )
        )
        vapor_molar_mass = sum(
            fraction * molar_mass
            for fraction, molar_mass in zip(
                self._state.mole_fractions_vapor(),
                self._component_molar_masses,
                strict=True,
            )
        )
        vapor_mass = molar_quality * vapor_molar_mass
        liquid_mass = (1.0 - molar_quality) * liquid_molar_mass
        return vapor_mass / (vapor_mass + liquid_mass)

    @staticmethod
    def _optional_property(getter: Callable[[], float]) -> float | None:
        try:
            value = float(getter())
        except (RuntimeError, ValueError):
            return None
        return value if isfinite(value) else None

    @staticmethod
    def _validate_pressure(pressure_pa: float) -> None:
        if not isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise ValueError("Pressure must be finite and positive")

    @staticmethod
    def _phase_name(phase_index: int) -> str:
        phases = {
            CP.iphase_liquid: "liquid",
            CP.iphase_gas: "gas",
            CP.iphase_twophase: "two_phase",
            CP.iphase_supercritical_liquid: "supercritical_liquid",
            CP.iphase_supercritical_gas: "supercritical_gas",
            CP.iphase_supercritical: "supercritical",
        }
        return phases.get(phase_index, "unknown")
