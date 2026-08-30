"""Kinetic mechanism registry, fuel definitions, and combustion bookkeeping.

One fuel is active per operating point; Jet-A and LNG are never co-fired. Each fuel
carries its own mechanism, and Jet-A additionally carries a second mechanism used only
for ignition-delay calculations. That split is not stylistic: measured during Phase 0,
a high-temperature-only mechanism overpredicts Jet-A ignition delay by 510x at 700 K and
71x at 800 K (20 atm, phi=0.5), which would declare an unsafe premixer safe. See
``mech/README.md`` for the evidence and ``docs/CRN_PLAN.md`` Section 4.4 for the design.

Mixture fraction, equivalence ratio, and stoichiometric air-fuel ratio are delegated to
Cantera's own implementations rather than reimplemented here. Cantera's
``mixture_fraction`` uses the Bilger element-based definition, which stays valid in the
presence of partially burned gas and unevaporated fuel, which is exactly the condition
inside an evaporator reactor.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache

import cantera as ct
from pydantic import BaseModel, ConfigDict, Field, field_validator

from fuelnozzle.models import LNGComposition, ModelWarning, WarningSeverity

#: Volume fraction of O2 in dry air, used for the standard 15% O2 emission correction.
DRY_AIR_O2_PERCENT = 20.9

#: Reference O2 percentage that gas-turbine emissions are conventionally corrected to.
REFERENCE_O2_PERCENT = 15.0

#: Species that must exist before any NOx result may be reported.
REQUIRED_NOX_SPECIES = ("NO", "NO2", "N", "N2O")

#: Default oxidizer composition, dry air on a mole basis.
DRY_AIR_MOLE_FRACTIONS = {"O2": 0.21, "N2": 0.79}


class FuelKind(StrEnum):
    """Which fuel circuit an operating point is burning."""

    JET_A = "jet_a"
    LNG = "lng"


class MechanismRole(StrEnum):
    """What a mechanism is allowed to be used for.

    A mechanism valid for one role is not automatically valid for the other.
    """

    NETWORK = "network"
    IGNITION_DELAY = "ignition_delay"


class MechanismError(RuntimeError):
    """A mechanism is missing, unloadable, or unfit for the requested purpose."""


@dataclass(frozen=True)
class NOxPathwayCoverage:
    """Which NO formation routes a mechanism actually carries.

    Recorded on every result rather than enforced as a fixed species list, because
    neither supplied mechanism is a superset of the other: the Jet-A model has the
    modern NCN prompt route but no NNH, while GRI-Mech 3.0 has NNH but the older
    Fenimore HCN route.
    """

    thermal: bool
    n2o_route: bool
    prompt_ncn: bool
    prompt_hcn: bool
    nnh_route: bool

    @property
    def summary(self) -> str:
        present = [
            name
            for name, flag in (
                ("thermal", self.thermal),
                ("N2O", self.n2o_route),
                ("prompt-NCN", self.prompt_ncn),
                ("prompt-HCN", self.prompt_hcn),
                ("NNH", self.nnh_route),
            )
            if flag
        ]
        return ", ".join(present) if present else "none"


class MechanismSpec(BaseModel):
    """A kinetic mechanism, the fuel it represents, and where it is valid.

    ``fuel_mole_fractions`` describes the fuel stream as the mechanism represents it:
    a single surrogate species for Jet-A, or the real component mix for LNG.
    """

    model_config = ConfigDict(frozen=True)

    path: str = Field(min_length=1)
    fuel: FuelKind
    role: MechanismRole
    fuel_mole_fractions: dict[str, float] = Field(min_length=1)
    provenance: str = Field(min_length=1)
    phase_name: str | None = None
    min_pressure_pa: float | None = Field(default=None, gt=0.0)
    max_pressure_pa: float | None = Field(default=None, gt=0.0)
    min_temperature_k: float | None = Field(default=None, gt=0.0)
    max_temperature_k: float | None = Field(default=None, gt=0.0)

    @field_validator("fuel_mole_fractions")
    @classmethod
    def validate_fuel_fractions(cls, values: dict[str, float]) -> dict[str, float]:
        if any(value <= 0.0 for value in values.values()):
            raise ValueError("Fuel mole fractions must be positive")
        total = sum(values.values())
        return {name: value / total for name, value in values.items()}

    @property
    def fuel_string(self) -> str:
        """Fuel composition in the ``'A:1, B:2'`` form Cantera expects."""
        return ", ".join(f"{name}:{value!r}" for name, value in self.fuel_mole_fractions.items())


@lru_cache(maxsize=32)
def _load_template(path: str, phase_name: str | None) -> ct.Solution:
    """Load and cache one mechanism for read-only interrogation.

    The cached object is never handed to a reactor. Cantera reactors mutate the
    ``Solution`` they are given, so sharing one across reactors aliases their states.
    Use :meth:`MechanismRegistry.new_solution` to obtain an independent instance.
    """
    try:
        return ct.Solution(path, phase_name) if phase_name else ct.Solution(path)
    except Exception as error:  # pragma: no cover - depends on user file
        raise MechanismError(f"Could not load mechanism {path!r}: {error}") from error


class MechanismRegistry:
    """The mechanisms available to a study, indexed by fuel and role.

    A Jet-A entry with role ``IGNITION_DELAY`` is optional. When it is absent the
    network mechanism is used for ignition delay too, and a warning records that the
    result may be badly wrong at premixer temperatures.
    """

    def __init__(self, specs: tuple[MechanismSpec, ...] | list[MechanismSpec]) -> None:
        self._specs: dict[tuple[FuelKind, MechanismRole], MechanismSpec] = {}
        for spec in specs:
            key = (spec.fuel, spec.role)
            if key in self._specs:
                raise MechanismError(
                    f"Duplicate mechanism for fuel {spec.fuel} and role {spec.role}"
                )
            self._specs[key] = spec
        if not self._specs:
            raise MechanismError("At least one mechanism is required")

    @property
    def specs(self) -> tuple[MechanismSpec, ...]:
        return tuple(self._specs.values())

    def spec_for(self, fuel: FuelKind, role: MechanismRole) -> MechanismSpec:
        """Return the mechanism for a fuel and role, falling back for ignition delay."""
        key = (fuel, role)
        if key in self._specs:
            return self._specs[key]
        if role is MechanismRole.IGNITION_DELAY:
            fallback = self._specs.get((fuel, MechanismRole.NETWORK))
            if fallback is not None:
                return fallback
        raise MechanismError(f"No mechanism registered for fuel {fuel} and role {role}")

    def has_dedicated_ignition_mechanism(self, fuel: FuelKind) -> bool:
        return (fuel, MechanismRole.IGNITION_DELAY) in self._specs

    def template(self, fuel: FuelKind, role: MechanismRole) -> ct.Solution:
        """Cached, read-only solution. Do not attach this to a reactor."""
        spec = self.spec_for(fuel, role)
        return _load_template(spec.path, spec.phase_name)

    def new_solution(self, fuel: FuelKind, role: MechanismRole) -> ct.Solution:
        """A fresh, independent solution suitable for attaching to a reactor."""
        spec = self.spec_for(fuel, role)
        try:
            return ct.Solution(spec.path, spec.phase_name) if spec.phase_name else ct.Solution(
                spec.path
            )
        except Exception as error:  # pragma: no cover - depends on user file
            raise MechanismError(f"Could not load mechanism {spec.path!r}: {error}") from error


def nox_pathway_coverage(solution: ct.Solution) -> NOxPathwayCoverage:
    """Report which NO formation routes a mechanism carries."""
    names = set(solution.species_names)
    return NOxPathwayCoverage(
        thermal=all(name in names for name in ("NO", "N", "O")),
        n2o_route="N2O" in names,
        prompt_ncn="NCN" in names,
        prompt_hcn="HCN" in names,
        nnh_route="NNH" in names,
    )


def validate_mechanism(
    spec: MechanismSpec,
    solution: ct.Solution,
    *,
    require_nox: bool,
    pressure_pa: float | None = None,
    temperature_k: float | None = None,
    lng_composition: LNGComposition | None = None,
) -> tuple[ModelWarning, ...]:
    """Check a mechanism is fit for purpose, raising only where silence would mislead.

    A missing nitrogen submodel raises rather than warns: it produces zero NO with no
    other symptom, and a silent zero is worse than a failed run.
    """
    names = set(solution.species_names)

    for fuel_species in spec.fuel_mole_fractions:
        if fuel_species not in names:
            raise MechanismError(
                f"Fuel species {fuel_species!r} is not present in mechanism {spec.path!r}."
            )

    if require_nox:
        missing = [name for name in REQUIRED_NOX_SPECIES if name not in names]
        if missing:
            raise MechanismError(
                f"Mechanism {spec.path!r} lacks nitrogen chemistry (missing "
                f"{', '.join(missing)}), so NOx cannot be predicted. It would report zero "
                "NO with no other symptom. Supply a mechanism with an N submodel."
            )

    warnings: list[ModelWarning] = []

    if require_nox:
        coverage = nox_pathway_coverage(solution)
        if not (coverage.prompt_ncn or coverage.prompt_hcn):
            warnings.append(
                ModelWarning(
                    code="NOX_PROMPT_PATHWAY_ABSENT",
                    severity=WarningSeverity.WARNING,
                    message=(
                        f"Mechanism {spec.path!r} carries neither NCN nor HCN, so prompt NO "
                        "is not represented. Rich-zone NO will be underpredicted."
                    ),
                )
            )
        warnings.append(
            ModelWarning(
                code="NOX_PATHWAY_COVERAGE",
                severity=WarningSeverity.INFO,
                message=(
                    f"Mechanism {spec.path!r} NO pathways: {coverage.summary}. Mechanisms "
                    "for different fuels may cover different pathways, so cross-fuel "
                    "absolute NOx is trend-level only."
                ),
            )
        )

    if lng_composition is not None and spec.fuel is FuelKind.LNG:
        unrepresented = [name for name in lng_composition.components if name not in names]
        if unrepresented:
            raise MechanismError(
                f"LNG components {', '.join(unrepresented)} are absent from mechanism "
                f"{spec.path!r}. Dropping them would silently change the fuel; supply a "
                "mechanism that represents them or restate the composition."
            )

    warnings.extend(_range_warnings(spec, pressure_pa, temperature_k))
    return tuple(warnings)


def _range_warnings(
    spec: MechanismSpec,
    pressure_pa: float | None,
    temperature_k: float | None,
) -> list[ModelWarning]:
    warnings: list[ModelWarning] = []
    checks = (
        ("pressure", pressure_pa, spec.min_pressure_pa, spec.max_pressure_pa, "Pa"),
        ("temperature", temperature_k, spec.min_temperature_k, spec.max_temperature_k, "K"),
    )
    for label, value, low, high, unit in checks:
        if value is None:
            continue
        if low is not None and value < low:
            warnings.append(
                ModelWarning(
                    code=f"MECHANISM_{label.upper()}_BELOW_VALIDITY",
                    severity=WarningSeverity.WARNING,
                    message=(
                        f"{label.capitalize()} {value:.4g} {unit} is below the declared "
                        f"validity floor {low:.4g} {unit} of mechanism {spec.path!r}."
                    ),
                )
            )
        if high is not None and value > high:
            warnings.append(
                ModelWarning(
                    code=f"MECHANISM_{label.upper()}_ABOVE_VALIDITY",
                    severity=WarningSeverity.WARNING,
                    message=(
                        f"{label.capitalize()} {value:.4g} {unit} exceeds the declared "
                        f"validity ceiling {high:.4g} {unit} of mechanism {spec.path!r}. "
                        "Results are extrapolative."
                    ),
                )
            )
    return warnings


def stoichiometric_air_fuel_ratio(
    solution: ct.Solution,
    spec: MechanismSpec,
    oxidizer: dict[str, float] | None = None,
) -> float:
    """Stoichiometric air-fuel ratio on a mass basis."""
    return float(
        solution.stoich_air_fuel_ratio(
            spec.fuel_string,
            oxidizer or DRY_AIR_MOLE_FRACTIONS,
            basis="mole",
        )
    )


def equivalence_ratio(
    solution: ct.Solution,
    spec: MechanismSpec,
    oxidizer: dict[str, float] | None = None,
) -> float:
    """Element-based equivalence ratio of the solution's current state.

    Valid for partially burned mixtures, unlike a ratio formed from unburned
    reactant amounts.
    """
    return float(
        solution.equivalence_ratio(
            spec.fuel_string,
            oxidizer or DRY_AIR_MOLE_FRACTIONS,
            basis="mole",
        )
    )


def bilger_mixture_fraction(
    solution: ct.Solution,
    spec: MechanismSpec,
    oxidizer: dict[str, float] | None = None,
) -> float:
    """Bilger element-based mixture fraction of the solution's current state."""
    return float(
        solution.mixture_fraction(
            spec.fuel_string,
            oxidizer or DRY_AIR_MOLE_FRACTIONS,
            basis="mole",
            element="Bilger",
        )
    )


def emission_index_g_per_kg(
    species_mass_fraction: float,
    exhaust_mass_flow_kg_s: float,
    fuel_mass_flow_kg_s: float,
) -> float:
    """Emission index in grams of species per kilogram of fuel burned."""
    if fuel_mass_flow_kg_s <= 0.0:
        raise ValueError("Emission index requires a positive fuel mass flow")
    return 1000.0 * species_mass_fraction * exhaust_mass_flow_kg_s / fuel_mass_flow_kg_s


def dry_mole_fraction(solution: ct.Solution, species: str) -> float:
    """Mole fraction of a species with water removed, the convention for emissions."""
    water = solution.X[solution.species_index("H2O")] if "H2O" in solution.species_names else 0.0
    if water >= 1.0:
        raise ValueError("Mixture is entirely water; a dry basis is undefined")
    return float(solution.X[solution.species_index(species)] / (1.0 - water))


def corrected_ppmv(
    dry_species_mole_fraction: float,
    dry_o2_mole_fraction: float,
    reference_o2_percent: float = REFERENCE_O2_PERCENT,
) -> float:
    """Dry ppmv corrected to a reference O2 level, the standard emissions basis.

    Correcting to a reference oxygen level removes the effect of dilution, so that two
    combustors running at different overall leanness can be compared.
    """
    measured_o2_percent = 100.0 * dry_o2_mole_fraction
    denominator = DRY_AIR_O2_PERCENT - measured_o2_percent
    if denominator <= 0.0:
        raise ValueError(
            "Measured dry O2 is at or above the O2 content of air; the correction to a "
            "reference O2 level is undefined."
        )
    factor = (DRY_AIR_O2_PERCENT - reference_o2_percent) / denominator
    return 1.0e6 * dry_species_mole_fraction * factor


@dataclass(frozen=True)
class MergeConflict:
    """One reason two mechanisms cannot simply be concatenated."""

    kind: str
    severity: WarningSeverity
    detail: str


#: Species that are the same chemistry under different spellings. A merge that does not
#: reconcile these produces two independent pools of what is physically one species.
KNOWN_ALIAS_PAIRS = (
    ("CH2*", "CH2(S)"),
    ("CH2S", "CH2(S)"),
    ("C2H5OH", "ETOH"),
)

#: Relative disagreement in molar enthalpy above which shared thermo is inconsistent.
THERMO_TOLERANCE = 0.02


def check_mechanism_merge(
    first: ct.Solution, second: ct.Solution
) -> tuple[MergeConflict, ...]:
    """Report why two mechanisms would conflict if merged.

    This checks; it does not merge. Combining a fuel mechanism with a nitrogen submodel
    is a chemistry decision with real consequences for the answer, and it belongs to
    whoever is accountable for that answer.
    """
    conflicts: list[MergeConflict] = []
    shared = set(first.species_names) & set(second.species_names)
    # Comparing thermodynamic data requires setting a temperature, which mutates the
    # object. The registry hands out cached read-only templates, so the original state
    # is restored before returning.
    original = (first.TPX, second.TPX)

    for name in sorted(shared):
        left = first.species(name)
        right = second.species(name)
        if left.composition != right.composition:
            conflicts.append(
                MergeConflict(
                    kind="composition_collision",
                    severity=WarningSeverity.ERROR,
                    detail=(
                        f"Species {name!r} is {left.composition} in one mechanism and "
                        f"{right.composition} in the other. The same name means two "
                        "different molecules; merging would silently conflate them."
                    ),
                )
            )
            continue
        for temperature in (300.0, 1000.0, 2000.0):
            try:
                first.TP = temperature, ct.one_atm
                second.TP = temperature, ct.one_atm
                left_h = first.standard_enthalpies_RT[first.species_index(name)]
                right_h = second.standard_enthalpies_RT[second.species_index(name)]
            except Exception:  # pragma: no cover - depends on mechanism data
                break
            scale = max(abs(left_h), abs(right_h), 1.0)
            if abs(left_h - right_h) / scale > THERMO_TOLERANCE:
                conflicts.append(
                    MergeConflict(
                        kind="thermo_inconsistency",
                        severity=WarningSeverity.WARNING,
                        detail=(
                            f"Species {name!r} has different thermodynamic data in the two "
                            f"mechanisms at {temperature:.0f} K. A merge must choose one "
                            "set; mixing them makes the combined mechanism non-conservative."
                        ),
                    )
                )
                break

    left_reactions = {reaction.equation for reaction in first.reactions()}
    duplicates = sorted(
        {reaction.equation for reaction in second.reactions()} & left_reactions
    )
    if duplicates:
        conflicts.append(
            MergeConflict(
                kind="duplicate_reactions",
                severity=WarningSeverity.WARNING,
                detail=(
                    f"{len(duplicates)} reactions appear in both mechanisms, for example "
                    f"{duplicates[0]!r}. Keeping both doubles their rate unless they are "
                    "explicitly marked as duplicates."
                ),
            )
        )

    first.TPX, second.TPX = original

    for left_name, right_name in KNOWN_ALIAS_PAIRS:
        if (left_name in first.species_names and right_name in second.species_names) or (
            right_name in first.species_names and left_name in second.species_names
        ):
            conflicts.append(
                MergeConflict(
                    kind="species_alias",
                    severity=WarningSeverity.ERROR,
                    detail=(
                        f"One mechanism spells a species {left_name!r} and the other "
                        f"{right_name!r}. These are the same chemistry; unless they are "
                        "renamed to match, the merged mechanism carries two separate "
                        "pools of one species."
                    ),
                )
            )

    return tuple(conflicts)
