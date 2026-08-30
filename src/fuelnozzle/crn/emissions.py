"""Emissions reporting: indices, corrected concentrations, and the ICAO LTO metric.

Raw exhaust concentrations cannot be compared between combustors, because a leaner one
dilutes its own output. Two conventions fix that. The emission index expresses output per
kilogram of fuel burned, and concentrations are quoted dry and corrected to a reference
oxygen level.

For certification the relevant number is not any single operating point but a weighted
total over the landing and take-off cycle, which is what this module assembles.

Only NOx carries an accuracy claim. CO and unburned hydrocarbons are computed and
reported because they are informative about the lean limit and about quench design, but
they were not part of the calibration set and are labelled accordingly wherever they
appear.
"""

from __future__ import annotations

from dataclasses import dataclass

import cantera as ct

from fuelnozzle.crn.chemistry import (
    corrected_ppmv,
    dry_mole_fraction,
    emission_index_g_per_kg,
)
from fuelnozzle.models import ModelWarning, WarningSeverity

#: Standard ICAO landing and take-off cycle: thrust setting and time in mode.
ICAO_LTO_MODES = (
    ("takeoff", 1.00, 0.7 * 60.0),
    ("climb_out", 0.85, 2.2 * 60.0),
    ("approach", 0.30, 4.0 * 60.0),
    ("idle", 0.07, 26.0 * 60.0),
)

#: Species reported without an accuracy claim.
UNCALIBRATED_SPECIES = ("CO", "OH", "CH4", "POSF10325")


@dataclass(frozen=True)
class EmissionsSummary:
    """Exhaust composition of one operating point, in reportable units."""

    temperature_k: float
    no_ppmv_dry_15pct_o2: float
    no2_ppmv_dry_15pct_o2: float
    nox_ppmv_dry_15pct_o2: float
    ei_nox_g_per_kg: float
    co_ppmv_dry_15pct_o2: float | None
    ei_co_g_per_kg: float | None
    dry_o2_mole_fraction: float
    warnings: tuple[ModelWarning, ...]


def summarize_emissions(
    solution: ct.Solution,
    temperature_k: float,
    pressure_pa: float,
    mass_fractions: dict[str, float],
    exhaust_mass_flow_kg_s: float,
    fuel_mass_flow_kg_s: float,
) -> EmissionsSummary:
    """Convert an exhaust state into emission indices and corrected concentrations.

    NOx is reported as the sum of NO and NO2, with NO2 counted as its NO2 mass, which is
    the certification convention and avoids the ambiguity of quoting them separately.
    """
    if fuel_mass_flow_kg_s <= 0.0:
        raise ValueError("Emissions reporting requires a positive fuel mass flow")

    solution.TPY = temperature_k, pressure_pa, mass_fractions
    warnings: list[ModelWarning] = []

    oxygen = dry_mole_fraction(solution, "O2")
    no = dry_mole_fraction(solution, "NO")
    no2 = dry_mole_fraction(solution, "NO2") if "NO2" in solution.species_names else 0.0

    no_ppm = corrected_ppmv(no, oxygen)
    no2_ppm = corrected_ppmv(no2, oxygen)

    def mass_fraction(name: str) -> float:
        if name not in solution.species_names:
            return 0.0
        return float(solution.Y[solution.species_index(name)])

    # NOx as equivalent NO2 mass, the certification convention.
    no2_molar_mass = 46.0055
    no_molar_mass = 30.0061
    ei_nox = emission_index_g_per_kg(
        mass_fraction("NO") * no2_molar_mass / no_molar_mass + mass_fraction("NO2"),
        exhaust_mass_flow_kg_s,
        fuel_mass_flow_kg_s,
    )

    co_ppm: float | None = None
    ei_co: float | None = None
    if "CO" in solution.species_names:
        co_ppm = corrected_ppmv(dry_mole_fraction(solution, "CO"), oxygen)
        ei_co = emission_index_g_per_kg(
            mass_fraction("CO"), exhaust_mass_flow_kg_s, fuel_mass_flow_kg_s
        )
        warnings.append(
            ModelWarning(
                code="CO_UNCALIBRATED",
                severity=WarningSeverity.INFO,
                message=(
                    "CO is reported as a diagnostic only. It was not part of the "
                    "calibration set and carries no accuracy claim; do not use it as a "
                    "quantitative lean-limit or efficiency prediction."
                ),
            )
        )

    return EmissionsSummary(
        temperature_k=temperature_k,
        no_ppmv_dry_15pct_o2=no_ppm,
        no2_ppmv_dry_15pct_o2=no2_ppm,
        nox_ppmv_dry_15pct_o2=no_ppm + no2_ppm,
        ei_nox_g_per_kg=ei_nox,
        co_ppmv_dry_15pct_o2=co_ppm,
        ei_co_g_per_kg=ei_co,
        dry_o2_mole_fraction=oxygen,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True)
class LTOMode:
    """One point of the landing and take-off cycle."""

    name: str
    thrust_fraction: float
    duration_s: float
    fuel_mass_flow_kg_s: float
    ei_nox_g_per_kg: float

    @property
    def nox_mass_g(self) -> float:
        return self.ei_nox_g_per_kg * self.fuel_mass_flow_kg_s * self.duration_s

    @property
    def fuel_mass_kg(self) -> float:
        return self.fuel_mass_flow_kg_s * self.duration_s


@dataclass(frozen=True)
class LTOResult:
    """Cycle total against the certification denominator."""

    modes: tuple[LTOMode, ...]
    rated_thrust_kn: float
    total_nox_g: float
    total_fuel_kg: float
    dp_foo_g_per_kn: float
    warnings: tuple[ModelWarning, ...]


def lto_dp_foo(
    modes: tuple[LTOMode, ...] | list[LTOMode], rated_thrust_kn: float
) -> LTOResult:
    """Cycle NOx mass per unit rated thrust, the ICAO reporting quantity.

    This is a faithful arithmetic assembly of the certification metric from the model's
    own emission indices. It is **not** a certification result: the underlying emission
    indices come from a reduced-order model that has not been validated against engine
    measurements, and certification requires measured data from the actual engine.
    """
    if rated_thrust_kn <= 0.0:
        raise ValueError("Rated thrust must be positive")
    if not modes:
        raise ValueError("At least one LTO mode is required")

    total_nox = sum(mode.nox_mass_g for mode in modes)
    total_fuel = sum(mode.fuel_mass_kg for mode in modes)

    warnings = [
        ModelWarning(
            code="LTO_NOT_A_CERTIFICATION_RESULT",
            severity=WarningSeverity.WARNING,
            message=(
                "Dp/Foo is assembled from modelled emission indices, not measured ones. "
                "It is useful for comparing designs against each other and must not be "
                "presented as a certification or compliance estimate."
            ),
        )
    ]
    covered = {mode.name for mode in modes}
    missing = [name for name, _, _ in ICAO_LTO_MODES if name not in covered]
    if missing:
        warnings.append(
            ModelWarning(
                code="LTO_CYCLE_INCOMPLETE",
                severity=WarningSeverity.WARNING,
                message=(
                    f"The cycle is missing {', '.join(missing)}. Dp/Foo computed from a "
                    "partial cycle is not comparable with a complete one; idle in "
                    "particular dominates the time in mode."
                ),
            )
        )

    return LTOResult(
        modes=tuple(modes),
        rated_thrust_kn=rated_thrust_kn,
        total_nox_g=total_nox,
        total_fuel_kg=total_fuel,
        dp_foo_g_per_kn=total_nox / rated_thrust_kn,
        warnings=tuple(warnings),
    )


def standard_lto_modes(
    fuel_flows_kg_s: dict[str, float], ei_nox_g_per_kg: dict[str, float]
) -> tuple[LTOMode, ...]:
    """Build the standard four-mode cycle from per-mode flows and emission indices."""
    modes: list[LTOMode] = []
    for name, thrust_fraction, duration in ICAO_LTO_MODES:
        if name not in fuel_flows_kg_s or name not in ei_nox_g_per_kg:
            continue
        modes.append(
            LTOMode(
                name=name,
                thrust_fraction=thrust_fraction,
                duration_s=duration,
                fuel_mass_flow_kg_s=fuel_flows_kg_s[name],
                ei_nox_g_per_kg=ei_nox_g_per_kg[name],
            )
        )
    return tuple(modes)
