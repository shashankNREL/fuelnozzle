# Dual-Fuel Reactor-Network Model: Implementation Log

Companion to [`CRN_PLAN.md`](CRN_PLAN.md). The plan records **design intent**; this file
records **what was actually built, what was verified, and where reality diverged from the
plan**.

Conventions, following the existing `V_AND_V_ROADMAP.md` house style:

- Every task gets an entry when it is completed, with its acceptance evidence.
- Every divergence from the plan is recorded in Section 3 with a rationale and a
  consequence, never silently absorbed.
- Verification output is quoted, not paraphrased. If a test fails, that is recorded.
- Open items are carried forward in Section 4 until closed.

---

## 1. Environment record

| Item | Value | Recorded |
|---|---|---|
| Plan approved | 2026-08-28 | — |
| Implementation started | 2026-08-28, Phase 0 | — |
| Python | 3.11 | pixi.toml |
| Cantera | **3.2.0** | 2026-08-28 |
| CoolProp | >=8.0 (pypi) | pre-existing |
| matplotlib | 3.11.1 | 2026-08-28 |
| PyYAML | 6.0.3 | 2026-08-28 |
| pixi CLI in use | v0.67.0 | 2026-08-28 |
| Baseline test count before any CRN work | 21 passing | 2026-08-28 |

---

## 2. Task log

### Phase 0 — Environment and scaffolding

#### Task 0.1 — Add Cantera and supporting dependencies — **COMPLETE** (2026-08-28)

Added to `pixi.toml` `[dependencies]` via `pixi add`:

```
cantera    >=3.2.0,<4
matplotlib >=3.11.1,<4
pyyaml     >=6.0.3,<7
```

**Acceptance evidence.**

```
cantera 3.2.0
matplotlib 3.11.1
gri30: 53 species, 325 reactions
NOx species present: True     # NO, NO2, N, N2O, NNH all present
```

```
$ pixi run test
21 passed in 35.48s
```

Both acceptance criteria met: Cantera imports and runs on osx-arm64, and the 21
pre-existing tests still pass with no edits to existing source.

**Notes.**

- GRI-Mech 3.0 ships with Cantera as `gri30.yaml` and loads directly, confirming the
  Section 4.4 decision that no conversion or vendoring is needed for the LNG mechanism.
  Its version is pinned transitively by the Cantera version, which is recorded above.
- Cantera resolved to **3.2.0**, above the `>=3.1` floor the plan assumed. The mole-based
  reactors and adaptive preconditioning that Section 4.1 depends on are present.
- See deviation **D-002** regarding the pixi lock file.

#### Tasks 0.4 / 0.5 — Mechanism directory, download, conversion — **COMPLETE** (2026-08-28)

Created `mech/` with immutable `mech/source/` originals and converted YAML alongside.
Full provenance, SHA-256 hashes, and the evidence below are recorded in
[`../mech/README.md`](../mech/README.md).

| File | Bytes (CTI) | Species | Reactions |
|---|---|---|---|
| `A2NOx_skeletal` | 116,044 | 71 | 538 |
| `A2NTCfast_ske` | 68,515 | 47 | 247 |

Converted with `cti2yaml` from Cantera 3.2.0; both validated on conversion and confirmed
by loading through `ct.Solution`.

**Acceptance evidence.**

- Both YAML files load. Species counts are 71 and 47 as expected.
- `A2NOx_skeletal` produces NO: mole fraction **1.536e-02** after 50 ms at 1800 K, 20 atm,
  phi = 1. Non-zero, so the nitrogen chemistry is live.
- `A2NTCfast_ske` reproduces an NTC region — **but not where the plan said to look.**
  See deviation **D-004**.

**Finding 1 — the two-mechanism split is quantitatively necessary, not stylistic.**
Ignition delay at premixer-relevant conditions (phi = 0.5), high-temperature-only
mechanism versus the low-temperature one:

| P | 700 K | 800 K | 900 K | 1000 K |
|---|---|---|---|---|
| 20 atm | **510x** | **71x** | 9.3x | 2.0x |
| 40 atm | **490x** | **86x** | 13.5x | 2.9x |

LTO `T3` is roughly 700-900 K. Using `A2NOx_skeletal` for the autoignition screen there
would overpredict ignition delay by one to two orders of magnitude — **declaring a
premixer safe when it is not**. Plan risk R1b is confirmed as real and safety-relevant.

**Finding 2 — NTC is pressure-dependent, and its disappearance is not what it looks
like.** A clear non-monotonic NTC zone exists at 5 atm (800-950 K), weakens at 10 atm, and
is gone by 20-40 atm. **But the non-monotonic *shape* vanishing does not mean
low-temperature chemistry stops mattering** — Finding 1 is measured at 20 and 40 atm,
where the shape is already gone, and the effect is still 71x-510x. Two different
statements; conflating them would be a serious error. Recorded in `mech/README.md`.

**Finding 3 — NOx pathway coverage differs between the two fuels' mechanisms.**

| Marker | `A2NOx_skeletal` (Jet-A) | `gri30` (LNG) |
|---|---|---|
| Prompt via `NCN` (modern) | **Yes** | No |
| Prompt via `HCN` (Fenimore) | Yes | Yes |
| NNH route | **No** | **Yes** |

Neither is a superset of the other. Each carries the pathway most relevant to its own
regime (NCN prompt for rich Jet-A, NNH for lean premixed methane), which is fortunate but
must be stated, not relied upon. This sharpens plan risk **R1e**: cross-fuel *absolute*
NOx is not strictly comparable. The `dA`/`dL` metric is within-fuel and unaffected.

#### Task 0.2 — `fuelnozzle.crn` package scaffold — **COMPLETE** (2026-08-28)

Created `src/fuelnozzle/crn/__init__.py` documenting scope, the fixed assumptions, and
the planned module layout. Per-module stub files were deliberately not created (deviation
**D-006**); modules appear as they are implemented, so an importable name is always a
working one.

**Acceptance evidence.** `import fuelnozzle.crn` succeeds; `pixi run lint` clean.

#### Task 0.3 — LaTeX skeleton and build tasks — **COMPLETE** (2026-08-28)

Created `docs/crn_technical_reference.tex` reusing the nozzle reference's preamble,
`\code` / `\srcline` / `\warningbox` macros, and beginner-oriented voice. Sections are
present as an outline with each unwritten one marked *pending*, so the document's own
table of contents doubles as a progress indicator.

`pixi.toml` tasks restructured:

```
doc-nozzle = "tectonic docs/technical_reference.tex --keep-logs"
doc-crn    = "tectonic docs/crn_technical_reference.tex --keep-logs"
doc        = { depends-on = ["doc-nozzle", "doc-crn"] }
```

**Acceptance evidence.** `pixi run doc` builds both PDFs with no LaTeX errors
(`crn_technical_reference.pdf` 45,443 bytes; `technical_reference.pdf` 161,396 bytes).

### Phase 0 — **COMPLETE** (2026-08-28)

All five tasks done. Closing state: `21 passed`, `ruff` clean, both PDFs building,
three mechanisms registered and verified.

### Phase 1 — Chemistry and stream layer

#### Task 1.1 — Cantera API surface pinned — **COMPLETE** (2026-08-28)

Verified by introspecting the installed **Cantera 3.2.0** rather than reading the
documentation page, since the reactor interface changed between 2.6 and 3.x. Full table
recorded in the LaTeX appendix "Cantera Interface Used by This Work".

Two findings:

- **Mole-based reactors are required, not preferred.** Sparse Jacobian preconditioning
  (`AdaptivePreconditioner`) is supported *only* for `IdealGasMoleReactor` and
  `IdealGasConstPressureMoleReactor`. This is what makes a 70-species network solve in
  seconds. Confirms the plan's Section 4.1 default.
- **`PlugFlowReactor` does not exist** — see deviation **D-007**.

#### Task 1.2 — `chemistry.py` — **COMPLETE** (2026-08-28)

`MechanismSpec`, `MechanismRegistry`, `FuelKind`, `MechanismRole`, `NOxPathwayCoverage`,
cached template loading with independent `new_solution`, plus stoichiometry, equivalence
ratio, Bilger mixture fraction, emission index, dry-basis and 15% O2 corrections.

Mixture fraction, equivalence ratio, and stoichiometric AFR are delegated to Cantera
(deviation **D-008**) rather than reimplemented.

**Acceptance evidence** — hand-calculation checks required by the plan:

| Quantity | Hand calculation | Code | Agreement |
|---|---|---|---|
| CH4 stoichiometric AFR | 17.13 | **17.127** | within 0.02% |
| C11H22 (POSF10325) stoichiometric AFR | 14.69 | **14.691** | within 0.01% |
| phi recovered from a prescribed mixture | 0.600000 | **0.600000** | exact to 6 dp |
| phi from Bilger Z via the Z-phi identity | 0.6 | 0.6 | within 1e-6 |

Also verified that the element-based equivalence ratio survives combustion: after burning
a phi = 0.8 mixture to completion (CH4 mole fraction below 1e-6), the reported phi is
still 0.8 to 1e-6. A ratio formed from unburned reactant amounts would have failed this.

#### Task 1.3 — Mechanism validation — **COMPLETE** (2026-08-28)

**Acceptance evidence.** A mechanism without nitrogen chemistry (`A2NTCfast_ske`)
**raises** `MechanismError` when NOx is requested, and validates cleanly when it is not.
GRI-Mech 3.0 at 30 atm emits `MECHANISM_PRESSURE_ABOVE_VALIDITY`. An `LNGComposition`
containing Ethane against a mechanism lacking it raises rather than dropping the
component. A fuel species absent from the mechanism raises.

Pathway coverage is *recorded*, not required — see deviation **D-009**.

#### Task 1.4 — `streams.py` — **COMPLETE** (2026-08-28)

`AirSplit` (validated to sum to 1.0, never silently renormalized), `CoolingAirDestination`,
`AirState`, `AirStreams`, and `resolve_air_streams` accepting either an explicit air flow
or an overall equivalence ratio.

Implements the Section 8.2.1 passage lever: `a_J + a_L = f_dome`, with
`near_field_air_fraction = a_active + chi * a_idle`.

**Acceptance evidence.** Station flows sum to the total; air flow and overall phi round-trip
through the stoichiometric AFR; density matches the ideal gas law to 1e-6; a split not
summing to 1.0 raises. The lever is verified directly: with `a_J` share 0.25 and chi = 0,
near-field air is 0.05 for Jet-A and 0.15 for LNG from identical hardware. At chi = 1 the
Jet-A near-field air rises to 0.20 and the lever disappears — the R12 sensitivity, made
explicit in a test.

### Phase 1 — **COMPLETE** (2026-08-28)

Closing state: **50 passed** (21 pre-existing, 29 new), `ruff` clean, both PDFs building.

### Phase 2 — Droplet physics

#### Tasks 2.1-2.6 — **COMPLETE** (2026-08-28)

New modules `droplets.py` (breakup, both evaporation branches, heating, energy
accounting), `liquids.py` (Jet-A and LNG property adapters behind one protocol), and
`spray_source.py` (nozzle-result to droplet-class bridge, three size policies,
Rosin-Rammler discretization). `JetAPropertyTable` extended additively (Task 2.4).

**Acceptance evidence.**

- D2-law recovered in the quiescent limit: `d^2` linear in time to within 1% of the
  initial value, after excluding the heat-up transient (see D-012).
- Boiling branch verified exactly: `dm/dt = -Q/L_vap`, droplet temperature rate is
  identically zero, and the branch switches where surface vapor pressure reaches ambient.
- Droplet energy split closes to 1e-10 relative: convective heat equals the sum of the
  latent and sensible terms. This is the Section 3.6 accounting that prevents
  double-counting latent heat.
- Weber number, surface vapor mass fraction (Eq. 21), and droplet number rate (Eq. 6)
  each verified against hand calculation.
- Rosin-Rammler classes conserve mass and reproduce the requested Sauter mean diameter
  to within 5% at 400 classes.
- Existing 21 tests untouched and passing after the `JetAPropertyTable` extension.

**Finding — two real bugs were caught by these tests, not by inspection.**

1. **Rosin-Rammler inversion was backwards.** The code set the characteristic diameter
   to `D32 / Gamma(1 - 1/n)` when the correct relation is `D32 = X / Gamma(1 - 1/n)`,
   requiring multiplication. The recovered Sauter mean diameter was 18.2 um against a
   requested 40 um, off by exactly `Gamma(0.6)^2`. Fixed.
2. **CoolProp returns no surface tension for a compressed-liquid state**, because
   surface tension is defined only along the saturation line. The LNG provider now walks
   a fallback chain: compressed state, then saturated liquid at the droplet temperature,
   then saturated liquid at the local pressure.

#### Paper Fig. 3 comparison (Tier V1.1 / V1.2) — **PARTIALLY REPRODUCIBLE**

Measured, using illustrative Jet-A properties at 1 atm:

| Quantity | Paper | This implementation |
|---|---|---|
| Breakup time | ~2e-5 s | **3.4e-5 s** (same order) |
| Post-breakup diameter | ~1 um | 19.9 um at 1 atm; **1.2 um at 20 atm** |
| Evaporation to one-tenth size, 350 K | ~6e-3 s | 0.80 s |
| Timescale ratio at 350 K | 133 | ~2.4e4 |
| Timescale ratio at 900 K | 21 | ~8.1e2 |

**Breakup agrees to an order of magnitude. Evaporation does not, and the reason matters.**
The paper's stated 6 ms for near-complete evaporation of a 150 um Jet-A droplet in 350 K
air is difficult to reconcile with physically reasonable kerosene properties: Jet-A boils
between roughly 420 and 570 K, its vapor pressure at 350 K is of order 1 kPa, and the
resulting Spalding transfer number is about 0.05. Our 0.8 s is consistent with published
kerosene evaporation constants; 6 ms would require a far more volatile liquid.

The paper does not state the ambient pressure for this case, nor the Jet-A liquid
properties used, so the case cannot be closed either way. See deviation **D-011**.

**What is reproduced is the paper's conclusion**, which is what the implementation
actually depends on: breakup is faster than evaporation by more than two orders of
magnitude at 350 K, the separation narrows with gas temperature, and it still exceeds one
order of magnitude at 900 K. That justifies treating breakup and evaporation sequentially.
The tests assert the conclusion and the trend, not the paper's absolute numbers.

### Phase 2 — **COMPLETE** (2026-08-28)

Closing state: **97 passed** (21 pre-existing, 76 new), `ruff` clean.

### Phase 3 — Reactors and network

#### Tasks 3.1-3.4 — **COMPLETE** (2026-08-28)

New modules `reactors.py` (zone specifications), `network.py` (mass balance, assembly,
steady solve, extinction detection), and `coupling.py` (operator-split droplet/gas
iteration).

**Acceptance evidence.**

- **Minimum-norm correction verified against a hand solution.** Three reactors in a
  loop: the balances force `z1 = z2 = t` and `z3 = t - 1`, so minimizing the squared
  change from (1.10, 0.95, 0.30) gives `t = 3.35/3 = 1.116667`. The code returns exactly
  that. A separate test confirms every other feasible solution is a larger change.
- A deliberately inconsistent split is closed to a residual below 1e-12 and reported.
- Globally inconsistent boundary flows **raise**, since no redistribution can fix them.
- Element balance closes to **1.4e-7** relative on a 71-species four-reactor network
  (limited by the marching tolerance, not by conservation).
- Recirculation converges: a reactor whose inlet depends on its own downstream settles
  without any tearing iteration.
- Lean methane at phi = 0.6 preheated to 750 K burns at 1982 K, within the expected band
  for 10 atm, and NO accumulates downstream (28.7 -> 41.4 -> 47.9 ppm) as residence time
  builds.
- Coupled solve converges in 4 iterations with 100% evaporation, 10.2 kW drawn from the
  gas, and NO rising 20.5 -> 26.2 ppm along the path.
- Oversized droplets produce liquid carryover and are flagged `LIQUID_CARRYOVER`.

**Finding — two real bugs, both caught by tests rather than inspection.**

1. **`ReactorNet.advance()` takes an absolute time, not a duration.** The fallback path
   called `advance(5.0)` after `advance_to_steady_state` had already moved the clock to
   28,235 s, so it tried to integrate backwards and raised. This masked every
   extinguished-network case as a solver failure; fixing it made extinction detection
   work.
2. **The fallback integration horizon was catastrophically oversized.** Five seconds of
   physical time on a system with millisecond residence times is roughly a thousand
   times more work than reaching steady state requires.

#### Performance: the steady solver was replaced — **325 s to 0.86 s**

`ReactorNet.advance_to_steady_state` proved unreliable on recirculating networks with
large mechanisms. On a four-reactor, 71-species case it ground for **210 s and then
failed**, falling through to the (oversized) time integration for a total of 325 s. The
same case time-marches to a converged answer in **0.86 s** with an identical peak
temperature of 1881 K.

Time-marching is now the primary method (deviation **D-014**): integrate in chunks of two
residence times, stop when every reactor temperature moves less than 0.05 K per chunk,
give up after 200 residence times with a warning. The steady solver remains available as
an opt-in polish. A single PSR still solves fine under `advance_to_steady_state`, so the
difficulty is specific to recirculating multi-reactor networks.

### Phase 3 — **COMPLETE** (2026-08-28)

Closing state: **133 passed** (21 pre-existing, 112 new), `ruff` clean.

### Phase 4 — Architecture templates

#### Tasks 4.1-4.3 — **COMPLETE** (2026-08-28)

`templates.py` provides RQL (with a staged quench), LDI, and LPP builders, each returning
a complete network description ready for `solve_coupled`. `check_mechanism_merge` was
added to `chemistry.py`.

**Acceptance evidence.**

- **Every template is mass-consistent before the network's correction runs**
  (correction norm below 1e-12). This matters: a template that needed correcting would
  have a bookkeeping error hidden behind a physical-looking repair.
- All supplied air enters the network; every spray-path zone hosts droplets; every
  architecture terminates in a plug-flow zone.
- RQL quench air is distributed across stages and sums to the quench fraction. The
  front-loaded schedule puts more air in the first stage than the uniform one.
- A dome sized for `phi_rich` = 1.5 passes without complaint; the default split gives an
  absurdly rich dome and raises `RQL_RICH_ZONE_OUT_OF_BAND`.
- The Section 8.2.1 lever is verified at the architecture level: with identical total
  air, changing only the injector passage share changes the near-field equivalence
  ratio, and setting the idle mixing fraction to 1.0 removes the lever entirely.

**End-to-end architecture comparison** (Jet-A, 20 atm, T3 = 800 K, illustrative splits):

| Architecture | phi near-field | Peak T | Exit T | NO @ 15% O2 | Evaporated |
|---|---|---|---|---|---|
| RQL | 1.43 | 2602 K | 1924 K | 1703 ppm | 1.000 |
| LDI | 0.93 | 2406 K | 1926 K | 1319 ppm | 1.000 |
| LPP | 0.93 | 2456 K | 1922 K | **728 ppm** | 0.994 |

The ordering LPP < LDI < RQL is the expected one, and **exit temperature is essentially
identical across all three while NOx differs by more than a factor of two** -- the same
asymmetry John et al. report, and the one Tier V2 exists to check.

**The absolute values are not credible and should not be quoted.** These air splits are
illustrative, not a design: a near-field `phi` of 0.93 is nearly stoichiometric, giving
peak temperatures of 2400-2600 K and correspondingly enormous thermal NO. A real
lean-burn dome runs near `phi` = 0.5-0.6, which requires most of the air through the
swirler. The comparison demonstrates the machinery works and orders architectures
sensibly; it is not a result.

#### Mechanism merge checker

Run against the two mechanisms actually in use, it reports why they cannot be
concatenated: **6 species with inconsistent thermodynamic data** (including OH),
**223 duplicate reactions**, and the **CH2\* versus CH2(S) alias collision**. That last
one is an error rather than a warning, since a merged mechanism would carry two separate
pools of the same species.

This independently supports the two-mechanism design of Section 4.4: merging these files
naively would produce a broken mechanism, quietly.

### Phase 4 — **COMPLETE** (2026-08-28)

Closing state: **165 passed** (21 pre-existing, 144 new), `ruff` clean.

### Phase 5 — Design-question modules

#### Tasks 5.1-5.7 — **COMPLETE** (2026-08-28)

`autoignition.py`, `thermal.py`, `emissions.py`, the additive `OperatingPoint` fields,
the idle-circuit screen, and `examples/run_combustor_study.py`.

**Finding 1 — flash cooling is real, and mostly is not latent heat.** Premixed mixture
temperature with air at 800 K, 20 atm:

| Case | T_mix | Drop | Latent share |
|---|---|---|---|
| Jet-A vapour at 470 K | 774.1 K | 25.9 K | 0.0 K |
| Jet-A, 30% liquid | 771.5 K | 28.5 K | 2.6 K |
| LNG vapour at 150 K | 747.2 K | 52.8 K | 0.0 K |
| LNG, 30% liquid | 742.8 K | 57.2 K | 4.4 K |

LNG cools the mixture about twice as much as Jet-A. **The latent contribution is only
about 4 K of a 57 K drop** -- the benefit is overwhelmingly the sensible enthalpy of cold
vapour, not the phase change. Worth recording because the opposite is easy to assume.

**Finding 2 — the premixing margin, quantified.** At 800 K air, 20 atm:

| Fuel | T_mix | tau_ign | Margin at 1 ms | at 3 ms |
|---|---|---|---|---|
| Jet-A | 774 K | 4.04 ms | 4.0 (safe) | 1.3 (marginal) |
| LNG | 747 K | no ignition in 10 s | unbounded | unbounded |

Jet-A tolerates roughly 1 ms of premixing at LTO conditions; LNG is effectively inert on
premixer timescales. This is the central hypothesis of plan Section 8.4, now measured.

**Finding 3 — pump pressure is what buys thermal design freedom.** Superheat at the
injector and subcooling in the line compete for the temperature gap between saturation at
feed pressure and at chamber pressure. With a 20 bar chamber:

| Feed pressure | Saturation gap | Feasible fuel temperature | Width |
|---|---|---|---|
| 2.5 MPa | 6.1 K | **none** | - |
| 3.0 MPa | 11.4 K | 171.0-172.0 K | 1.0 K |
| 4.0 MPa | 20.2 K | 171.0-181.0 K | 10.0 K |

**Below about 3 MPa feed pressure there is no feasible fuel temperature at all.** No
amount of heating fixes it; the gap itself is too narrow. This connects directly to the
pressure budget the package already computes.

**Robustness gap found and fixed.** Above methane's critical pressure (4.6 MPa) there is
no saturation, and the window calculation raised a raw CoolProp message. It now raises
`SupercriticalFeedError` explaining that supercritical injection is a legitimate design
choice to which these constraints simply do not apply.

**Lean sweep (Answer 2) behaves cleanly.** As the LNG dome fraction rises and `phi` falls
from 1.47 to 0.60: NOx falls 1284 to 20 ppm, peak temperature 2397 to 1962 K, and the
exit temperature spread narrows 581 to 142 K. Textbook lean-burn behaviour.

#### Answer 1 sanity anchor — **NOT cleanly met; recorded as an open discrepancy**

Sweeping the head-end air fraction with **quench air held fixed** (see D-020) gives a
genuine NOx minimum, but at `phi_rich` = 1.22 rather than the 1.4-1.6 band that classical
RQL practice indicates. The curve is shallow between roughly 1.2 and 1.4.

Per the plan's own acceptance criterion for Task 8.7c, a minimum outside that band is to
be treated as a possible tool defect and investigated **before** being reported as a
result. It is therefore recorded here as open (**O-005**), not as a validation. Two
candidate causes, both plausible and both untested:

- a perfectly stirred quench zone mixes more slowly than real quench jets, so it dwells
  near stoichiometric and overproduces NO. Absolute EI_NOx is several times what real RQL
  hardware achieves, which is consistent with this.
- the optimum's location is sensitive to quench stage count and quench time, neither of
  which is calibrated.

### Phase 5 — **COMPLETE** (2026-08-28)

Closing state: **207 passed** (21 pre-existing, 186 new), `ruff` clean, example runs.

### Phase 7 — Verification and documentation

#### Tasks 7.1-7.9 — **COMPLETE** (2026-08-29)

`tests/test_crn_verification.py` adds limits, convergence studies, and the Tier V2 paper
comparison. The LaTeX reference is complete; `README.md`, `docs/modeling.md`, and
`docs/V_AND_V_ROADMAP.md` updated.

**Verification achieved.**

| Check | Result |
|---|---|
| Residence-time identity (V1.5) | Reproduced, both quoted cases |
| Long-residence reactor to equilibrium | Within 2% |
| Extinction below a critical volume | Reproduced |
| Zero-fuel limit | Air unchanged to 1 K |
| Plug-flow segment convergence | Under 5 K between 4 and 8 segments, monotone |
| Droplet solver tolerance | Under 1e-4 across rtol 1e-6 to 1e-10 |
| Ignition-table resolution | Under 15% under fourfold refinement |

#### Tier V2 — **REPRODUCED** (2026-08-29)

Identical topology and flows, run once with droplets and once with the fuel gaseous and
perfectly mixed:

| | phi range | Peak T | Exit T | NO |
|---|---|---|---|---|
| With spray (LFRN) | 0.514-1.714 | 2453 K | 1926.6 K | 1868 ppm |
| Gaseous premixed (GFRN) | 0.514-0.514 | 1928 K | 1928.0 K | 57 ppm |

**Exit temperature differs by 0.07%. NO differs by a factor of 33.** That is the paper's
asymmetry: both methods get the global energy balance right and only one gets the
emissions right. Their reported contrast — a gaseous network underpredicting NO by 54-91%,
peak temperatures near 1800 K against above 2200 K with spray — matches in direction and
character. **Per the plan's exit criterion, paper-comparison work is now closed.**

**Getting there took three corrections to my own test, each of which is instructive:**

1. **Hand-written flows ignored recirculated mass**, producing a 37% mass correction — the
   exact mistake `_chain_flows` was built to prevent, made by hand in a test. The
   comparison was running on heavily corrected flows.
2. **All air entered the first zone**, so every reactor sat at the same temperature and
   there was nowhere for a rich pocket to exist. The topology could not exercise the
   mechanism at all. Fixed by staging air across dome, mixer, and flame.
3. **The gaseous case was not the paper's gaseous case.** Injecting vapour only at the
   dome retains the air-staging heterogeneity; the two cases then differed by 3%, not 33x.
   The paper's GFRN is fuel premixed with *all* the air, giving one equivalence ratio
   everywhere (they report 0.652-0.653). Distributing fuel to every zone in proportion to
   its air is what produced the signature.

Also confounding an earlier attempt: 30 um droplets left 5.5% liquid carryover, so the
spray case ran cooler simply because less fuel burned. Reduced to 12 um for full
evaporation, and a test now asserts that the comparison is not confounded that way.

**Solver robustness fixed.** A zone with a long residence time exceeded Cantera's default
20000-step cap. `ReactorNet.max_steps` is now raised to 200000.

### Phase 7 — **COMPLETE** (2026-08-29)

Closing state: **223 passed** (21 pre-existing, 202 new), `ruff` clean, both PDFs build.

### O-005 resolved — quench sensitivity study (2026-08-29)

The rich-quench-lean NOx optimum fell at `phi_rich` = 1.22 rather than the expected
1.4-1.6, with absolute levels several times real hardware. The plan's own acceptance
criterion required investigating this before reporting any RQL optimum. **The cause was an
unrealistically long quench residence time, not a defect in the reactor model.**

Prevaporized fuel was used for the sweep after verifying it reproduces the fully coupled
solution to 2.7% and runs five times faster; evaporation completes upstream of the quench,
so this isolates it cleanly.

**Quench stage count**, at fixed quench volume and `phi_rich` = 1.35:

| Stages | tau_q (ms) | EI_NOx | Change |
|---|---|---|---|
| 1 | 4.07 | 58.99 | — |
| 2 | 4.31 | 112.29 | +90.3% |
| 5 | 4.49 | 126.15 | +4.0% |
| 12 | 4.58 | 131.49 | +1.5% |
| 20 | 4.60 | 133.08 | +1.2% |

**A single mixing point underpredicts NOx by more than half**, which is direct
confirmation that modelling the quench as a chain was necessary. The result is converged
to about 1% by twelve stages; the previous default of five was roughly 5% low.

**Quench residence time**, at twelve stages:

| Quench volume | tau_q (ms) | EI_NOx |
|---|---|---|
| 2.0e-4 | 0.92 | 68.5 |
| 1.0e-3 | 4.58 | 131.5 |
| 4.0e-3 | 18.34 | 198.6 |

**This is the dominant influence**, and it is physical rather than numerical. Real quench
jets penetrate and mix in roughly half a millisecond to a millisecond, which is the entire
point of the architecture. The previous default gave 4.6 ms.

**Effect on the optimum**, sweeping `phi_rich` at twelve stages:

| Quench | Optimum | Minimum EI | Curve |
|---|---|---|---|
| ~0.9 ms (realistic) | **phi = 1.35** | 68.5 | sharp: 141, 133, 108, 82, **69**, 77, 109, 149 |
| ~4.6 ms (old default) | phi = 1.22 | 130.6 | flat: 159, 165, 155, 141, **131**, 131, 142, 162 |
| ~18 ms (slow) | phi = 2.34 | 176.9 | no interior minimum at all |

With a realistic quench the optimum lands at `phi_rich` = 1.35 with a pronounced minimum,
which meets the sanity anchor. A slow quench flattens the curve and drifts the optimum
rich; a very slow one destroys it entirely.

**Changes made.** Default quench stages raised from 5 to 12; default quench volume reduced
from 1.0e-3 to 2.0e-4; `RQL_QUENCH_UNDER_RESOLVED` warns below twelve stages;
`check_quench_residence_time` warns when the achieved quench exceeds 2 ms. Six regression
tests added.

**What remains, and is now a stated limitation rather than an open item.** Absolute
EI_NOx is still roughly twice what real RQL hardware achieves. The zones are adiabatic and
a perfectly stirred quench still mixes more slowly than real jets even at one millisecond.
**Relative comparison and the optimum's location are supported; absolute level is not.**

### Phase 8 — Dual-fuel design optimization

#### Tasks 8.1-8.4, 8.7b, 8.9 — **COMPLETE** (2026-08-29)

`design.py` (Class A / A2 / B variables, bounds, mission points), `evaluate.py` (one
design across a mission set), `objectives.py` (four objectives, feasibility-first
constraints), `optimize.py` (sensitivity, Latin-hypercube sampling, Pareto front,
cost-of-shared-liner, focused sweeps). 26 tests.

**Stage 1 sensitivity — the most useful single output, as the plan predicted.** Elasticity
of each objective to each variable, from a one-at-a-time screen (19 designs, 209 s):

| Variable | Jet-A NOx | LNG NOx |
|---|---|---|
| `dome_air_fraction` | **-3.66** | **-12.09** |
| `quench_air_fraction` | -2.59 | -7.70 |
| `quench_volume_m3` | +1.20 | 0.00 |
| `flame_volume_m3` | +0.50 | +1.19 |
| `jet_a_passage_share` | -0.16 | -0.22 |
| fuel temperature | +0.055 | +0.046 |

Three things follow.

1. **The dome air fraction dominates both fuels**, by two to three times over the next
   variable and by twenty to seventy times over the injector variables. This is direct
   support for the plan's claim that the shared liner air split is the crux.
2. **Quench volume has exactly zero effect on LNG**, because the LNG path uses a lean
   premixed architecture with no quench zone. A correct structural result, and a useful
   check that the per-fuel architecture selection is actually taking effect.
3. **Fuel temperature barely moves NOx directly** (elasticity ~0.05). This *partly
   challenges* the emphasis the plan placed on thermal management as a NOx lever. Its real
   influence is through the autoignition *constraint*, which is not an objective, and
   through atomization. Worth stating because the opposite is easy to assume.

**The passage lever (8.7b), quantified across the full assumption range.** Near-field
equivalence ratio at fixed hardware, dome = 0.38:

| chi | Jet-A phi at share 0.15 -> 0.85 | LNG phi at share 0.15 -> 0.85 |
|---|---|---|
| 0.0 | 9.02 -> 1.59 | 1.59 -> 9.01 |
| 0.5 | 2.35 -> 1.46 | 1.46 -> 2.35 |
| 1.0 | 1.35 -> 1.35 | 1.35 -> 1.35 |

**At full mixing the lever disappears entirely** -- identical phi for both fuels at every
share. The Section 8.2.1 idea works only insofar as idle-passage air stays segregated,
exactly as risk R12 anticipated, and the tool now reports it across the whole range rather
than at one assumed value.

#### Cost of sharing a liner — machinery works, result **not yet reportable**

Two sweeps were run and **both were caught by robustness guards**, which is the intended
behaviour of Task 8.9 rather than a failure.

*First attempt*, 36-design Latin hypercube over 9 variables (330 s, 13 feasible): the best
"compromise" came out **identical to the LNG-optimal design**. With 36 samples in 9
dimensions the sample simply contained no genuine intermediate. Added
`SHARED_LINER_COST_DEGENERATE`, which refuses to present penalties as a cost of sharing
when the compromise is just one of the endpoints.

*Second attempt*, focused 10-point sweep over `dome_air_fraction` alone, using the
sensitivity ranking to reduce dimension (2970 s):

| dome | Jet-A EI | LNG EI |
|---|---|---|
| 0.22 | 139.79 | 71.30 |
| 0.38 | 62.53 | 9.04 |
| 0.50 | 21.01 | **2.03** (LNG optimum) |
| 0.54 | 13.42 | 2.30 (compromise) |
| 0.58 | **8.96** (Jet-A optimum) | 3.32 |

Not degenerate, giving dA = +49.7% and dL = +13.5%. **But the Jet-A optimum sits at the
last point of the sweep**: its curve is monotonic throughout, so the optimum was never
bracketed. The cause is that with `jet_a_passage_share` = 0.5 and chi = 0, only half the
dome air reaches the Jet-A near field, so effective `phi_rich` spans only 4.67 down to
1.77 and never reaches the ~1.35 optimum established when O-005 was resolved.

Added `SINGLE_FUEL_OPTIMUM_UNBRACKETED`, which the degeneracy check alone does not catch:
the compromise can be interior while a single-fuel optimum is pinned to the boundary.

**Therefore dA and dL are not reported as the cost of dual-fuel operation.** What is
supported is the *direction*, which is itself meaningful: the two fuels want different
dome air fractions, and the conflict the plan predicted is real and measurable.

**Practical constraint discovered.** Evaluation cost is not uniform. A design near the
baseline solves in about 11 s, but designs with a large dome fraction and almost no
dilution air become very stiff and take several minutes, so the focused sweep averaged
297 s per design. Any Phase 8 campaign must budget for this rather than assuming a flat
per-design cost.

### Phase 8 — **PARTIALLY COMPLETE** (2026-08-29)

Delivered: 8.1, 8.2, 8.3, 8.4, 8.7b, 8.9. Deferred: 8.5 (NSGA-II), 8.6 (architecture-pair
comparison), 8.7/8.7c-f (the reportable cost chart and per-path targets), 8.8, 8.10.

Closing state: **255 passed** (21 pre-existing, 234 new), `ruff` clean.

### Bracketed dome sweep — the bracketing failure had a physical cause (2026-08-29)

O-006 asked for a re-run at `jet_a_passage_share` = 1.0 to bracket the Jet-A optimum.
That exact setting is unusable: with chi = 0 it leaves the LNG passage with **zero**
near-field air, so LNG's equivalence ratio diverges and half the comparison is
meaningless. The sweep was instead run at **chi = 1.0**, where idle-passage air reaches
the near field and both fuels see the full dome air on identical hardware. That is also
the conservative case, since the passage lever contributes nothing.

Quench, primary, and cooling air were trimmed to leave dilution headroom. Evaluation cost
fell from about **297 s to 7-12 s per design**, which confirms the O-007 diagnosis:
networks with almost no dilution air are extremely stiff.

**The Jet-A optimum still would not bracket, and the reason is not the sweep range.**

First, two checks. The solver is exactly deterministic: the same design returned
`EI = 143.4792` on three consecutive runs. And a refined grid showed the curve is smooth,
not jagged -- the apparent jaggedness was misreading a coarse grid. What the refined grid
showed instead was a **maximum** at `phi_rich` = 1.43, where the O-005 study had found a
**minimum** at 1.35.

The only substantive difference between the two studies was quench air fraction. Testing
it directly:

| Quench air | Curve shape | Best point |
|---|---|---|
| 0.15 | monotone rising | rich edge, no interior optimum |
| 0.20 | **maximum** at `phi` = 1.43 | rich edge, no interior optimum |
| 0.30 | **interior minimum** | `phi` = 1.29, EI 67.6 |
| 0.38 | **interior minimum** | `phi` = 1.35, EI 43.2 |

**The rich-quench-lean optimum exists only when there is enough quench air.** Below
roughly 0.25 to 0.30 the curve inverts: instead of `phi_rich` near 1.4 being the best
choice, it becomes the worst, and nitric oxide rises monotonically as the dome is leaned
toward it.

This is physically coherent and it extends the O-005 finding rather than contradicting
it. Rich-quench-lean works by crossing stoichiometric quickly, and the crossing speed
depends on the quantity of quench air as well as on the time available. O-005 established
the time dependence; this establishes the quantity dependence. **The O-005 resolution
holds, but is conditional on quench air of about 0.30 or more**, which is what that study
used. A regression test now covers both branches.

**Consequence for the cost of sharing a liner.** The bracketing failure was never a
sampling problem, so widening the sweep would not have fixed it. At the quench-air
fraction used there is no interior Jet-A optimum to find. `dA` and `dL` therefore remain
unreportable, but for a now well-understood reason rather than an unexplained one.

---

## 3. Deviations from the plan

Recorded as they occur. "Consequence" states what a future reader must know.

| # | Plan said | What actually happened | Why | Consequence |
|---|---|---|---|---|
| **D-001** | `cantera >=3.1` | Resolved to **3.2.0** | Latest conda-forge build at install time | None adverse. Section 4.1 features (mole reactors, `AdaptivePreconditioner`, `ExtensibleReactor`, `advance_to_steady_state`) are all present in 3.2. Task 1.1 will pin the exact API surface against 3.2.0 rather than 3.1. |
| **D-002** | *(not anticipated)* | `pixi add` reported: *"Lock-file version 7 is newer than supported. Maximum supported version: 6 (pixi v0.67.0). The lock-file will be treated as missing and regenerated."* `pixi.lock` was **regenerated from v7 to v6 format**. | The committed `pixi.lock` had been written by a newer pixi than the v0.67.0 CLI available in this environment | **Action required by the user.** The regenerated lock is a downgrade in format. If other machines or CI use a newer pixi, either (a) run `pixi self-update` here and regenerate, or (b) accept the v6 lock repo-wide. Flagged as open item **O-001**. The environment itself resolves and works correctly either way. |
| **D-004** | Task 0.5 acceptance: *"`A2NTCfast_ske` reproduces an NTC region (non-monotonic ignition delay vs 1000/T) between roughly 650 and 900 K"* | The first check ran at **20 atm** and **failed** — delay was monotonic. Re-running across 5/10/20/40 atm found a clear NTC zone at **5 atm, 800-950 K**, weak at 10 atm, absent at 20-40 atm. | The plan's acceptance criterion did not specify a pressure, and NTC in jet fuels is strongly pressure-dependent — it washes out at combustor pressure | **Criterion amended:** the NTC check must be run at 5-10 atm. The original 20 atm failure was a defective test, not a defective mechanism. More importantly, this produced Finding 2 above, which corrects a conflation latent in the plan's own wording. |
| **D-005** | Plan Section 4.4 treated `A2NTCfast_ske` as *the* Jet-A low-temperature mechanism | Its header describes it as a **"hypothetical A2 fuel with fast NTC"** — a sensitivity variant, not a validated real-Jet-A model. Stanford publishes fast and slow NTC variants that bracket the real uncertainty. | Not noticed at plan time; only visible in the file header | Fast NTC gives the **shorter** delay and is therefore **conservative** for an autoignition safety screen, so current use is safe but possibly pessimistic. Ideally the margin is reported as a bracket across both variants. Tracked as **O-004**. |
| **D-006** | Task 0.2: *"Create `src/fuelnozzle/crn/` with `__init__.py` and module stubs"* | Package created with a documented `__init__.py`; **empty per-module stub files were not created**. Modules will be added as they are implemented. | Seventeen empty files that lint-pass but do nothing are navigational noise and invite half-written imports | None adverse. The planned layout is documented in the package docstring instead. |
| **D-007** | Plan Section 4.1 referred to Cantera's native steady **`PlugFlowReactor`** as a Phase 7 cross-check | No such class exists in Cantera 3.2.0. The class is **`FlowReactor`**; "Plug Flow Reactor" is only the documentation's prose name | Documentation prose was taken for an API name | Name corrected in the LaTeX appendix. No design impact: plug flow is represented by a series of stirred reactors, as planned, and `FlowReactor` remains the optional cross-check. |
| **D-008** | Task 1.2 listed Bilger mixture fraction and element-based equivalence ratio as things to implement | **Delegated to Cantera's own `mixture_fraction`, `equivalence_ratio`, and `stoich_air_fuel_ratio`** rather than reimplemented | Cantera 3.2 provides all three, with Bilger as an explicit option. Reimplementing well-tested library code adds bug surface for no benefit | Less code, less risk. Verified against hand calculation anyway (see Task 1.2). |
| **D-009** | Task 1.3: *"warn if `NNH` and `HCN` are absent, since prompt NO will then be missing"* | Changed to **recording** which pathways each mechanism carries, and warning only when *no* prompt route (neither NCN nor HCN) is present | Measured in Phase 0: the Jet-A mechanism has NCN but no NNH, and GRI-Mech has NNH but no NCN. A hard NNH requirement would wrongly flag the better Jet-A mechanism | Pathway coverage is attached to every result instead, which is more informative than a pass/fail and directly supports the R1e caveat. |
| **D-010** | Task 0.2 stated modules would appear as implemented | Confirmed working: `chemistry.py` and `streams.py` are the first two | — | None. |
| **D-011** | Plan Section 12.2 listed V1.1/V1.2 as *"exact quantitative unit tests from fully-specified paper data"*, asserting breakup at ~2e-5 s, evaporation at ~6e-3 s, and ratios of 133 and 21 | Breakup reproduces to an order of magnitude; **evaporation does not, and the paper's case is underspecified** (no ambient pressure, no liquid properties stated). The paper's 6 ms evaporation appears inconsistent with physically reasonable Jet-A volatility | The Fig. 3 case cannot be reconstructed from what the paper publishes | **Tests assert the paper's conclusion and trend rather than its numbers**: separation above two orders of magnitude at 350 K, above one at 900 K, and narrowing with temperature. This matches the 2026-08-28 direction to get trends and conclusions right. The absolute-number comparison is recorded here and closed. |
| **D-012** | Plan Section 12.1 listed the D2-law check as *"diffusion-limited branch recovers `d^2` linear decay in a quiescent, constant-property limit"* | Linearity holds only **after the heat-up transient**. The first test failed with a 10% residual because it fitted through the warm-up period, where the film temperature and hence the vaporization rate are still moving | The classical D2 law assumes the droplet has reached its wet-bulb temperature; the plan's wording omitted that condition | Test restricted to the settled region (droplet within 1 K of its final temperature), where the residual is under 1%. The transient is physically real, not numerical error. |
| **D-013** | Plan Section 3.3 quoted the paper's breakup criterion as `A + We_g > 1` | Implemented the standard Taylor analogy criterion `We_c + A > 1`, where `We_c = We_g/12` | The paper's printed form cannot be what they computed: for their own Fig. 3 case `We_g` is about 38, so the criterion would fire instantly and no breakup timescale would be resolvable. Their figure shows a finite ~2e-5 s | Documented in the `taylor_analogy_breakup` docstring. The standard criterion reproduces a finite breakup time of the right order, so it is almost certainly what the paper implemented. |
| **D-014** | Plan Section 4.2: *"Solution: `ReactorNet.advance_to_steady_state()`, with `AdaptivePreconditioner` ... and a fallback of long-time `advance()`"* | **Reversed**: bounded time-marching is primary, the steady solver is opt-in polish | Measured: `advance_to_steady_state` ground for 210 s and then failed on a four-reactor 71-species recirculating network that time-marches to the same answer in 0.86 s. Risk R3 materialized, in the form the plan anticipated | 375x faster and observable. Convergence is now checked explicitly against a temperature tolerance rather than delegated. Cost: element balance closes to ~1e-7 rather than ~1e-11, reflecting the marching tolerance. |
| **D-015** | Plan Section 4.1 listed `PressureController` for the network exit | Every connection, including the exit, uses `MassFlowController` | All flows are prescribed by the balanced split matrix, so there is nothing for a pressure controller to determine; prescribing the exit flow directly is simpler and deterministic | None observed. Mass and element balances close. Revisit if a future topology needs pressure-driven flow. |
| **D-016** | Plan Section 4.3 described the operator split abstractly | The network is **rebuilt** each outer iteration rather than mutated | Evaporating fuel changes the total gas mass flow, so the mass balance itself changes between iterations; mutating the inlet list would leave the balance stale | Slightly more work per iteration, but the balance is always consistent with the sources. |
| **D-017** | Plan Section 5.2: LPP *"Chemistry is frozen in the premixer by construction"* | The premixer is modeled as a **reacting** zone | Freezing chemistry would guarantee the assumption rather than test it. Left reacting, the premixer ignites in the model exactly when residence time exceeds ignition delay, which is the real failure mode of the architecture | More honest, and it makes the failure visible in the solution rather than only in a separate screen. Recorded on every LPP result as `LPP_PREMIXER_REACTS`. |
| **D-018** | *(not anticipated)* | `_chain_flows` initially omitted recirculated mass from the chain segments it re-traverses | Recirculated gas re-enters upstream and travels the chain a second time | Fixed. Without it, every template would have arrived at the network with an imbalance for the correction to absorb, which would have disguised a bookkeeping error as a physical repair. Templates now balance to machine precision. |
| **D-019** | *(not anticipated)* | `check_mechanism_merge` mutated the `Solution` objects passed to it, because comparing thermodynamic data requires setting a temperature | The registry hands out **cached** read-only templates, so this would have silently corrupted state shared by every other caller | Original state is captured and restored. Covered by a regression test. |
| **D-020** | *(not anticipated)* | The first head-end air sweep varied dome air and quench air **together**, because the example rescaled all downstream stations to keep the split summing to one. The result was monotonic with no minimum | Confounded experiment: the rich-zone trend and the quench-rate trend moved at once | Quench air is now held fixed for that sweep, taking the dome's air from dilution, which isolates the rich-zone effect and produces a genuine minimum. The scaling variant is retained for the lean sweep, where a fixed quench would drive dilution negative. A caution about the distinction is printed with the result. |
| **D-021** | Task 5.5 called for a `combustor_study.py` orchestration module | Orchestration is implemented in `examples/run_combustor_study.py` instead | A generic multi-point orchestrator is what Phase 8's `evaluate.py` must be, and building a second one now would duplicate machinery that Phase 8 would then replace | No capability is missing; the example demonstrates the full chain end to end. Phase 8 will provide the reusable form. |
| **D-022** | Task 5.1 planned a flashback screen | Implemented, but **suppressed unless the user supplies a laminar flame speed** | Computing it requires a flame calculation this tool does not perform, and inventing one would put an unmarked guess underneath a safety screen | Consistent with the package's refusal to invent an SMD. Reports nothing rather than something unfounded. |
| **D-023** | Plan Section 12.2 Tier V2 described building the gaseous comparison by introducing the fuel as vapour | The gaseous case distributes fuel to **every zone in proportion to its air**, giving a uniform equivalence ratio | Injecting vapour only at the dome keeps the air-staging heterogeneity and is not what the paper means by a gaseous network. It reproduces a 3% difference instead of the 33x one | The signature only appears with a genuinely premixed comparison. Documented in the test and in the LaTeX so the distinction is not lost. |
| **D-024** | *(not anticipated)* | `ReactorNet.max_steps` raised from the Cantera default of 20000 to 200000 | A zone of order a cubic metre has a residence time of seconds, and one marching chunk then exceeds the default cap | Long-residence zones now solve. No effect on ordinary combustor volumes. |
| **D-025** | Plan Task 8.5 specified NSGA-II via `pymoo` for the Pareto search | Implemented a **plain non-dominated sort over Latin-hypercube samples**, with no new dependency | At roughly 11 s to several minutes per design, an evolutionary run of the thousands of evaluations NSGA-II needs is hours to days. The affordable sample size is in the hundreds, where a quadratic non-dominated sort is exact and a population-based heuristic has no advantage | The Pareto machinery is complete and dependency-free. `pymoo` remains an option once a surrogate makes evaluation cheap enough to justify it. |
| **D-026** | *(not anticipated)* | `DesignVector.with_values` used `model_copy(update=...)`, which **does not re-run pydantic validators** | Discovered when a sweep produced a negative dilution fraction deep inside the network builder rather than a clean rejection at construction | Rebuilt through the constructor so validation always runs. This was a real correctness bug: any sweep could silently construct an impossible liner. Covered by a regression test. |
| **D-027** | Plan Task 8.9 described the robustness gate as re-checking conclusions against the sensitivity range and a mechanism substitution | The gate is additionally implemented as **two automatic guards inside the cost calculation**: degeneracy, and an unbracketed single-fuel optimum | Both failure modes appeared in the first two real sweeps. A gate applied only at reporting time would have let a plausible-looking dA/dL through | The guards fire before a number is presented, and both are error severity. |
| **D-003** | *(pre-implementation)* | `pypdf` was temporarily added to `pixi.toml` during planning to extract text from the reference paper, then removed and `pixi.toml`/`pixi.lock` restored via `git checkout` | No PDF text extractor was available in the environment | None. Recorded only so a future reader who sees it in shell history knows it was intentional and reverted. |

---

## 4. Open items

| # | Item | Raised | Status |
|---|---|---|---|
| **O-001** | `pixi.lock` was regenerated from format v7 to v6 (see D-002). Decide whether to `pixi self-update` and restore v7, or standardize on v6. | 2026-08-28 | **Open — needs user decision.** Does not block Phase 0. |
| **O-002** | No LNG-side spray calibration exists, so LNG SMD and all atomization-quality outputs that depend on it remain gated (plan risk R2). One flashing flow test would also pin the Tier 2 relaxation time `tau` and convert the `L/D` geometric target from a trend into a number (plan Section 8.11.2). | 2026-08-28 | Open — experimental dependency, not a code blocker. |
| **O-006** | Cost-of-shared-liner sweep would not bracket the Jet-A optimum. | 2026-08-29 | **DIAGNOSED 2026-08-29, still blocking.** Not a sampling problem: at quench air 0.20 there is **no interior optimum to bracket**, because the RQL curve inverts below about 0.25-0.30 quench air. Re-run at quench air >= 0.30 with the LNG side bounded by its stability limit rather than a NOx minimum. |
| **O-008** | The LNG lean-premixed optimum is a **stability boundary, not a stationary point** -- NOx falls monotonically toward lean until blowout. The `SINGLE_FUEL_OPTIMUM_UNBRACKETED` guard flags it, which is technically right but conceptually misleading. Teach the guard to distinguish "edge because the range was too narrow" from "edge because a constraint binds there". | 2026-08-29 | Open — affects reporting, not correctness. |
| **O-007** | Evaluation cost varies from ~11 s to ~300 s per design depending on how little dilution air remains. Investigate why low-dilution networks are so stiff before running a large campaign. | 2026-08-29 | Open — performance, not correctness. |
| **O-004** | `A2NTCfast_ske` is the *fast* NTC variant of a *hypothetical* A2 fuel (see D-005). Download the *slow* variant and report the autoignition margin as a bracket across both. Current use is conservative, so this is an accuracy improvement, not a safety gap. | 2026-08-28 | Open — low priority. |
| **O-005** | RQL NOx optimum location and absolute level. | 2026-08-28 | **RESOLVED 2026-08-29.** Cause was an unrealistically long quench residence time, not the reactor model. With a realistic quench the optimum is `phi_rich` = 1.35 with a sharp minimum. Defaults corrected, guards and six regression tests added. Absolute level remains about 2x real hardware and is now a stated limitation. |
| **O-003** | Decide whether CO is promoted from uncalibrated diagnostic to calibrated QoI. Affects whether the LNG lean limit (plan Section 8.10.2) is directional or quantitative. Not needed before Phase 8. | 2026-08-28 | Open — needs user decision. |

---

## 5. Verification record

Cumulative test counts, so regressions in the existing package are visible immediately.

| Date | Event | Existing tests | New CRN tests | Total |
|---|---|---|---|---|
| 2026-08-28 | Baseline before CRN work | 21 pass | 0 | 21 pass |
| 2026-08-28 | After Task 0.1 (Cantera added) | 21 pass | 0 | 21 pass |
| 2026-08-29 | Quench-air conditionality of the RQL optimum established; O-006 diagnosed | 21 pass | **235 pass** | **256 pass**; `ruff` clean |
| 2026-08-29 | **Phase 8 partial** (8.1-8.4, 8.7b, 8.9; sensitivity + guards) | 21 pass | **234 pass** | **255 pass**; `ruff` clean |
| 2026-08-29 | **O-005 resolved** (quench sensitivity study; defaults corrected) | 21 pass | **208 pass** | **229 pass**; `ruff` clean |
| 2026-08-29 | **Phase 7 complete** (7.1-7.9 verification, LaTeX, repo docs) | 21 pass | **202 pass** | **223 pass**; `ruff` clean; both PDFs build |
| 2026-08-28 | **Phase 5 complete** (5.1-5.7 autoignition, thermal, emissions, example) | 21 pass | **186 pass** | **207 pass**; `ruff` clean; example runs |
| 2026-08-28 | **Phase 4 complete** (4.1-4.3 templates, quench schedule, merge checker) | 21 pass | **144 pass** | **165 pass**; `ruff` clean |
| 2026-08-28 | **Phase 3 complete** (3.1-3.4 reactors, network, coupling) | 21 pass | **112 pass** | **133 pass**; `ruff` clean |
| 2026-08-28 | **Phase 2 complete** (2.1-2.6 droplets, liquids, spray bridge) | 21 pass | **76 pass** | **97 pass**; `ruff` clean |
| 2026-08-28 | **Phase 1 complete** (1.1 API, 1.2 chemistry, 1.3 validation, 1.4 streams) | 21 pass | **29 pass** | **50 pass**; `ruff` clean; both PDFs build |
| 2026-08-28 | **Phase 0 complete** (0.2 scaffold, 0.3 LaTeX, 0.4/0.5 mechanisms) | 21 pass | 0 | 21 pass; `ruff` clean; both PDFs build |

---

## 6. Independent physics and equation audit (2026-08-30)

### 6.1 Scope, intended use, and disposition

This audit reviewed commit `1f8ef3e7` with the requested use case held fixed:

- Jet-A at all four ICAO landing-and-take-off (LTO) modes;
- LNG at a representative cruise condition;
- separate fuel injectors, one active fuel at a time, and a shared liner; and
- a conceptual design result, not merely a demonstration that a CRN can run.

The review traced the production data path from `OperatingPoint` and the nozzle solvers
through spray, thermal, chemistry, network, emissions, objective, and optimization code. It
also compared implementation claims against the tests and CRN documentation. No source-code
fixes were made as part of this audit.

**Disposition: no-go for physics-based conceptual design in the current state.** The branch
is a useful research prototype for isolated equation studies and qualitative sensitivity
work. It is not yet suitable for selecting a dual-fuel architecture, sizing a combustor, or
reporting a Jet-A LTO optimum or LNG cruise lean limit. Several design-driving quantities are
either disconnected from the evaluator, represented by invalid proxies, or computed with a
different model than the documentation describes. The existing validation disclaimers in
the README and V&V roadmap therefore govern all current use.

Severity in this section has the following meaning:

- **Blocker:** can reverse an architecture/design decision, violate conservation, or return
  an unsafe state as acceptable.
- **Major:** important model-form omission or equation limitation that prevents quantitative
  use but does not always invalidate a qualitative trend.
- **Moderate:** robustness, traceability, or reporting issue that must be resolved before a
  defensible release.

### 6.2 What is physically sound or useful

The following foundations should be retained:

1. Keeping the fuels mutually exclusive at an operating point and providing separate fuel
   hardware is appropriate for the stated concept.
2. CoolProp-based LNG state calculations, explicit pressure budgets, and suppression of
   uncalibrated LNG and Jet-A SMD values are sound modeling practices.
3. Cantera is used correctly for stoichiometric AFR, element-based equivalence ratio, Bilger
   mixture fraction, and finite-rate homogeneous chemistry.
4. The NOx emission index correctly treats NO as NO2-equivalent mass. The standalone
   `lto_dp_foo()` arithmetic and the ICAO time-in-mode table are also correct.
5. The dedicated network and ignition-delay mechanisms acknowledge an important physical
   fact: a high-temperature Jet-A mechanism cannot safely screen low-temperature
   autoignition.
6. The hand-solved minimum-norm flow-correction fixture, analytical droplet limits, element
   accounting, and explicit warning objects are useful software-verification foundations.
7. The code correctly refuses to invent a flashback result when no laminar flame speed is
   supplied and correctly labels many spray correlations as requiring calibration.

These points do not offset the blockers below because most are not connected into one
conservative, mission-level design evaluation.

### 6.3 Blocking findings

| ID | Finding and physical consequence | Evidence |
|---|---|---|
| B-01 | **The design evaluator bypasses the advertised nozzle-to-combustor chain.** It injects both fuels as prescribed, fully gaseous streams. Injector geometry, LNG flash quality, droplet distribution, evaporation, TAB breakup, spray calibration, and liquid carryover therefore have no effect on the optimized design. A cold liquid fuel is represented as cold gas without paying the latent-heat requirement. | `crn/evaluate.py:181-198`; contrast `crn/__init__.py:3-7` and `crn/coupling.py` |
| B-02 | **A configured architecture pair is not closed as one physical shared liner.** The evaluator can assign a Jet-A topology and a separate LNG topology, but this is a logical pairing of templates. It does not demonstrate that their different zone connections, air-entry locations, effective areas, swirler/injector footprints, and recirculation structures coexist in one geometry. Topology is an evaluator setting rather than a design variable. | `crn/evaluate.py:114-137,156-171`; `crn/templates.py`; `CRN_PLAN.md:863-864,1345-1346` |
| B-03 | **A declared PFR is actually one fully back-mixed PSR.** `plug_flow_segments` and `segment_volume_m3` are never consumed by `CombustorNetwork.solve()`. The post-flame volume, where thermal NO can continue to form, consequently has the wrong residence-time distribution. The production PFR convergence test manually builds unrelated PSRs and cannot catch this defect. | `crn/reactors.py:54-83`; `crn/network.py:309-317`; `tests/test_crn_verification.py:207-237`; `crn_technical_reference.tex:535-538` |
| B-04 | **Nominal, simulated, and reported reactor residence times are inconsistent.** A constant-pressure Cantera reactor is seeded hot at the nominal volume; its initial mass is then retained while its volume changes with state. Extraction instead reports `rho_converged * volume_nominal / mdot`. The seed composition is also formed by mass-flow-weighting inlet mole fractions rather than combining inlet molar flows. Quench time, extinction, NO residence, and volume optimization therefore do not necessarily use the physical volume or inventory shown to the user. | `crn/network.py:305-317,442-468` |
| B-05 | **Adding fuel makes the air-only template internally inconsistent and invokes an unphysical flow repair.** The minimum-norm correction redistributes the missing fuel mass among forward and recirculation edges rather than constructing physically complete flows. It is unweighted, so edges with different uncertainty and physical roles receive equal mathematical treatment. | `crn/templates.py:398-440`; `crn/evaluate.py:181-198`; `crn/network.py:65-172` |
| B-06 | **A negative corrected flow is still accepted.** The solver omits its mass-flow controller, but `inflow_of()` and residence-time reporting retain the signed value. The graph that is solved is therefore different from the graph that is balanced and reported. Global balance also does not guarantee balance of disconnected components. | `crn/network.py:130-143,253-261,334-339` |
| B-07 | **The stated air-system and pressure-loss physics are absent.** `liner_pressure_loss_fraction` is unused; all reactors are assigned one pressure; air fractions are fixed inputs rather than results of effective area, discharge coefficient, density, and pressure ratio; and recirculation flows are prescribed rather than momentum- or pressure-driven. The same fractions are reused at LTO and cruise. This prevents a fixed liner from being evaluated consistently across the mission. | `operating.py:42-49`; `crn/network.py:300-345`; `crn/design.py:73-117`; `crn/templates.py` |
| B-08 | **Spray coupling omits the breakup it advertises.** `apply_aerodynamic_breakup` is created but never consumed, so the TAB fallback has no effect. The TAB step response itself omits the exponential/sine terms of the documented damped oscillator and rejects all overdamped breakup, although a sufficiently large constant forcing can produce monotonic distortion. | `crn/spray_source.py:213-224`; `crn/coupling.py:311-342`; `crn/droplets.py:137-219`; `crn_technical_reference.tex:374-419` |
| B-09 | **Droplet transport, mass, and integrated energy closure are incomplete.** Injection velocity is used as the gas-relative velocity through every zone; there is no drag, slip relaxation, trajectory, or gas velocity. Local fuel-vapor mass fraction is forced to zero, so vapor accumulation cannot inhibit evaporation. Coupling reconstructs heat from zone endpoints and has no explicit liquid outlet enthalpy ledger. Radius evolution uses instantaneous density without its temperature derivative, while evaporated fraction assumes constant density; the physical mass implied by `rho(T) * volume` can therefore change during heating while reported mass remains fixed, especially for LNG. | `crn/coupling.py:297-390`; `crn/droplets.py:279-300,332-357,450-472` |
| B-10 | **The LNG droplet treatment is not a flashing model.** Vaporization is powered only by gas convection, not by a metastable droplet's own superheat. A droplet initialized above saturation receives zero temperature derivative and remains superheated rather than being driven to saturation. Vapor and liquid phase compositions from a multicomponent flash are collapsed into one species. | `crn/droplets.py:310-335,399-472`; `crn/coupling.py:67-91,216-274` |
| B-11 | **Runtime chemistry validation never runs.** Mechanism species/pathway and range checks are called only by tests. The production evaluator can therefore use GRI-Mech 3.0 above its stated pressure basis or a mechanism with missing fuel/NOx coverage without the documented runtime warning. Species presence is also not evidence that the relevant NOx reactions or pathway rates are adequate. | `crn/chemistry.py:208-315`; no production callers |
| B-12 | **Real LNG composition cannot propagate into combustion.** CoolProp names such as `Methane` are compared literally with Cantera names such as `CH4`, so even represented C1-C3 components fail the guard. More importantly, `MechanismSpec.fuel_mole_fractions` is independent of `LNGComposition`; the main path always burns pure CH4 even if nozzle properties use a mixture. This makes phase behavior, stoichiometric AFR, phi, ignition, flame speed, and emissions compositionally inconsistent. | `models.py:36-55`; `crn/chemistry.py:266-273`; `examples/run_combustor_study.py:64-78` |
| B-13 | **Unknown autoignition states are returned as safe.** An out-of-table query and “did not ignite during the integration horizon” both become `None`, and `None` becomes `SAFE`. A query above the 1000 K table ceiling is thus labeled safe exactly where ignition delay is shortest. Pressure and phi are nearest-neighbor selections, with no stoichiometric point in the default grid. | `crn/autoignition.py:245-272,301-320`; `crn/evaluate.py:49-57,145-153` |
| B-14 | **Autoignition is disconnected from the physical passage.** The screen uses total dome air rather than the active-passage air plus the specified mixed share of idle-passage air. The design residence-time scalar does not change premixer volume, flow, mixing, evaporation, or network solution. The hypothetical fast-NTC A2 variant is conservative but is not a validated Jet-A uncertainty bound. | `crn/evaluate.py:262-293`; `crn/design.py:105-113`; `mech/README.md:33-44` |
| B-15 | **The Jet-A objective is not an ICAO LTO objective.** `weighted_ei_nox()` time-weights emission indices instead of weighting by fuel mass burned. The correct `sum(EI * mdot_f * time)` Dp/Foo helper exists but is disconnected from optimization. The example evaluates one takeoff-like point, not takeoff, climb-out, approach, and idle. | `crn/evaluate.py:98-108`; `crn/objectives.py:113-128`; `crn/emissions.py:130-231`; `examples/run_combustor_study.py:169-207` |
| B-16 | **The LNG cruise lean limit is not implemented.** There is no continuation to extinction, combustion-efficiency threshold, calibrated CO criterion, or dynamic/flameholding screen. A hot-initialized temperature threshold is not lean blowout. The example prints a five-point model sweep and calls the result a limit. | `crn/network.py:381-411`; `examples/run_combustor_study.py:222-242`; `CRN_PLAN.md:946-973` |
| B-17 | **Two optimization objectives are not the quantities their names imply.** “Mixing nonuniformity” is the range of mean phi among serial and functionally different zones; “exit temperature spread” is the range across every reactor. Neither represents unmixedness, an exit traverse, radial temperature distribution, or pattern factor. `rank_key()` then adds g/kg, kelvin, and dimensionless values, making the ranking depend on units and numerical scale. | `crn/evaluate.py:233-257`; `crn/objectives.py:113-144` |
| B-18 | **The example's LNG path is physically wrong.** It supplies `jet_a_liquid()` to the methane case, describes the fuel as vapor-fed while creating zero initial vapor and a liquid droplet population, and uses Jet-A fuel flow in the LNG premix calculation. It therefore cannot serve as a reference result. | `examples/run_combustor_study.py:97-115,210-263` |
| B-19 | **The model has no combustor geometry closure.** Absolute sector volumes and mass flows are accepted without an annulus, cup count, reference area, length, liner effective areas, jet momentum, dome swirler, or periodic-sector scaling check. It cannot produce the promised combustor dimensions or show that a selected set of volumes fits the shared hardware and pressure-loss budget. | `crn/design.py`; `crn/templates.py`; `operating.py:12-55` |
| B-20 | **No experimental validation supports design acceptance.** The 255 reported tests primarily verify code identities, monotonic trends, or internally generated regressions. The John et al. case does not reconstruct their unpublished topology or data, the RQL optimum was used to tune model defaults, and no experimental dataset is read by CRN tests. | `tests/test_crn_verification.py`; `docs/V_AND_V_ROADMAP.md:455-459` |
| B-21 | **The nozzle-layer “cruise” reference condition is not a representative combustor inlet.** `examples/run_study.py` and many LNG nozzle tests use `P3 = 0.1 MPa`, while the CRN example uses about `2 MPa`. This greatly changes LNG saturation margin, pressure ratio, choking, flash quality, and spray regime. No engine cycle deck or plausibility check prevents the 1-atm case from being read as a cruise-combustor design result. | `examples/run_study.py:45`; `examples/run_combustor_study.py:53-54`; `tests/test_lng_equilibrium.py`; `tests/test_lng_relaxation.py`; `tests/test_spray.py` |

### 6.4 Major model-form and equation limitations

#### 6.4.1 Network, mixing, heat transfer, and numerics

1. All reactor kinds map to the same constant-pressure stirred-reactor class. They are
   adiabatic by default, with only optional fixed-watt heat sources; “evaporator” and
   “mixer” are labels unless the separate coupling path is explicitly used.
2. Steady convergence checks temperature change only. Slow species, especially NO, can still
   change after temperature appears stationary. Species, energy, mass, and ODE residuals are
   not convergence criteria.
3. Hot equilibrium initialization finds a lit branch; it does not show that a flame can
   ignite, remain attached, or recover after a disturbance. Cold/hot multiplicity and
   hysteresis are not mapped.
4. Fixed heat-loss watts are not related to liner area, wall temperature, emissivity,
   coolant flow, or heat-transfer coefficients. Rich Jet-A zones can also radiate through
   soot, which is not modeled.
5. `CoolingAirDestination` is unused. Cooling air always enters the downstream post-flame
   zone, so primary-zone film cooling, dilution, and exhaust routing cannot be distinguished.
6. The quench is an arbitrary staged-PSR schedule. Refinement of that schedule demonstrates
   numerical consistency with the assumed mixer, not physical quench-jet mixing or scalar
   dissipation. The RQL optimum's strong dependence on quench fraction is therefore expected
   model-form sensitivity, not validation.
7. Recirculation ratio and zone volumes are fixed across operating conditions. Neither is
   linked to swirler aerodynamics, density, pressure loss, or flame state.
8. The solver's integration-horizon estimate is documented as using cold density but is
   normally evaluated from the hot ignition seed, which can understate the time required for
   a cold or extinguished state.
9. Mechanism extrapolation and some unexecuted safety checks do not prevent design
   acceptance. Network non-convergence itself is correctly propagated as an infeasible
   point.

#### 6.4.2 Jet-A atomization and evaporation

1. Jet-A is a single pseudo-component with one boiling point, molecular weight, diffusivity,
   and vapor-pressure relation. It cannot represent the roughly 100 K distillation range,
   preferential evaporation, or changing vapor composition of a real kerosene.
2. The pressure-swirl model is a screening construction, not a validated simplex-atomizer
   solution. Flow is sized with the full exit area and a prescribed discharge coefficient,
   after which an air core is inferred; swirl is capped by available energy. These steps need
   hardware data before air-core, film, cone, cavitation, or turndown decisions.
3. The default pressure-swirl discharge coefficient is `0.65`, while the repository's own
   cited Lacava case reports about `0.273`. Because predicted flow is linear in this
   coefficient, the default is not a transferable pressure-swirl value and can materially
   under-size area unless hardware data replace it.
4. The optional SMD relation lacks an external-air breakup model and is intentionally
   calibration-dependent. Ambient gas density is calculated but does not enter the SMD
   relation, so the calibrated SMD has no chamber-pressure dependence. A scalar calibration
   at one LTO mode cannot establish transfer to the others.
5. Rosin-Rammler midpoint classes only approximately reconstruct the requested D32 at
   practical class counts. The 400-class reconstruction test does not establish convergence
   for the default class count used in coupled runs.
6. Wall impingement, film formation, splash, multi-orifice interaction, spray cone
   dispersion, and circumferential maldistribution are absent. These omissions are
   especially important for idle CO/UHC and liner temperature.

#### 6.4.3 LNG nozzle, thermal management, and phase behavior

1. Tier 2 imposes a linear pressure path and equilibrium velocity, then relaxes quality in
   time. The reported bounded mass flux need not equal `actual_density * velocity`. Diameter
   is calculated afterward and does not participate in momentum, pressure, or residence-time
   closure. The current model can screen a relaxation length but cannot establish an `L/D`
   design rule.
2. For a choked case, the CFD boundary is reported at chamber pressure while the operating
   flux is taken at the upstream critical point. Boundary density, velocity, area, and mass
   flow therefore need not satisfy continuity.
3. The regime classifier can call a case fully flashing from pressure ratio or Jakob number
   even when the finite-rate exit quality remains small. This is a screening label, not a
   phase-distribution prediction.
4. `thermal_window()` checks feed-line subcooling at feed pressure but computes heating duty
   at chamber pressure. The heat-addition path must be evaluated at the actual pressurized
   feed state, followed by the nozzle expansion/flash.
5. For multicomponent LNG, a temperature just above the bubble point is reported with quality
   one; the bubble-to-dew enthalpy interval is labeled entirely latent even though it includes
   temperature glide. Total enthalpy and vapor quality should come directly from the
   equilibrium state rather than pure-fluid logic.
6. Target superheat is a temperature difference from chamber saturation, not an isenthalpic
   flash quality. The actual design variable must close feed heating, nozzle pressure drop,
   flashing, and downstream phase enthalpy in one path.
7. The advertised four-way thermal window omits autoignition from
   `ThermalWindowPoint`. Available aircraft/engine heat is a scalar ceiling rather than a
   heat-exchanger effectiveness, pinch, pressure-drop, mass, and mission analysis.
8. Idle Jet-A coking and idle-LNG vapor-lock are threshold screens and are not called by the
   design evaluator. They do not model thermal soak, wetted volume, purge, restart, deposits,
   or fuel switching.
9. The feed-line two-phase model is deliberately treated as an error state, but if retained
   for diagnosis it also needs acceleration pressure drop, minor-loss treatment, and
   two-phase heat-transfer/flow-regime validation.
10. When a feed line is present, its calculated outlet enthalpy is retained but its
    calculated outlet pressure is replaced by the independently requested nozzle-budget
    pressure. A warning does not close these incompatible boundaries; the feed and nozzle
    must be solved as one pressure/enthalpy path or rejected.

#### 6.4.4 Chemistry, emissions, ignition, and operability

1. GRI-Mech 3.0 is a reasonable methane screening mechanism, but its commonly stated
   optimization range extends only to roughly 10 atm. Cruise combustor calculations near
   20 atm are extrapolative; higher-pressure C1-C3 chemistry and NOx pathways require
   validation or a better-suited mechanism.
2. `A2NOx_skeletal` is explicitly high-temperature chemistry. It should not be relied on for
   reacting LPP premixer behavior, idle CO/UHC, or low-temperature ignition. Its declared
   mechanism metadata also lacks pressure, temperature, phi, reduction-target, and error
   ranges.
3. `A2NTCfast_ske` is a hypothetical fast-NTC A2 sensitivity variant. Using it alone is a
   conservative screen, not a quantified uncertainty band for Jet-A.
4. The 400 K ignition temperature-rise marker is reproducible but not universal. Two-stage
   ignition, cool flames, OH/CH2O markers, and maximum-temperature-rise-rate definitions can
   give different delays. Constant-pressure ignition should also be reconciled with the
   source experiments used to validate the mechanism.
5. The premix enthalpy calculation evaluates gaseous fuel thermochemistry below the NASA
   polynomial range: for example, LNG is evaluated at 150 K although GRI species
   polynomials generally begin at 200 K. Liquid/flash enthalpy should enter through the
   property model, not extrapolated ideal-gas methane.
6. Cross-fuel absolute NOx comparison is not defensible with mechanisms of different
   provenance, reduction error, and prompt/NNH pathway coverage. Relative conclusions must
   survive mechanism substitution and experimental comparison.
7. CO is explicitly uncalibrated, UHC is not actually aggregated, and no soot/nvPM model
   exists. The Jet-A mechanism lacks the larger PAH population needed for soot prediction;
   GRI-Mech has no aromatic/soot chemistry. These outputs cannot establish a quantitative
   lean limit or ICAO nvPM result.
8. The mean state of a 0-D reactor does not resolve the temperature/composition PDF that
   controls thermal NOx. Until CFD/rig calibration supplies residence, mixing, and
   unmixedness information, absolute NOx is a model output rather than a validated
   prediction.
9. Flashback is not connected to feasibility and uses an unvalidated generic turbulent
   multiplier when invoked. Lean blowout, combustion efficiency, ignition energy,
   light-around, altitude relight, fuel-transfer transients, thermoacoustics, and dynamic
   stability are absent.
10. A steady CRN is appropriate for stabilized point emissions but is not sufficient for
    safety-critical premixer design or fuel-switch operability.

#### 6.4.5 Mission aggregation, objectives, and optimization

1. `MissionPoint` duplicates rather than adapts `OperatingPoint`. Combustor air flow, overall
   phi, liner loss, thrust, wall temperature, pressure budget, sector multiplier, and nozzle
   result are not closed in one canonical state.
2. `thrust_fraction` is carried but does not affect physics or aggregation. Rated thrust is
   absent from `DesignEvaluator`, so Dp/Foo cannot be a native objective.
3. Autoignition margins between one and the requested design floor of four are correctly
   converted into optimizer constraint violations. Unknown ignition delay and uninvoked
   mechanism-range checks remain the unsafe gap.
4. No constraint enforces combustion efficiency, material temperature, pressure loss, liner
   cooling, pattern factor, nozzle turndown/cavitation, spray carryover, flashback, or
   idle-circuit thermal limits in the production objective. Packaging is only a user-supplied
   scalar budget, not a result of injector/liner geometry.
5. Shared-liner optimum bracketing assumes the tuple is a sorted one-dimensional sweep. It is
   invalid for unordered Latin-hypercube samples or a multidimensional design space.
6. Sensitivity and Pareto tools operate on deterministic point values. Model-form,
   mechanism, calibration, property, manufacturing, and operating-condition uncertainty are
   not propagated into feasibility or ranking.

### 6.5 Assessment of the stated assumptions

| Assumption | Assessment for the target design |
|---|---|
| One active fuel at a time; separate fuel hardware | **Appropriate**, but purge, trapped fuel, switching transients, and inactive-circuit heat soak remain required operability cases. |
| Shared liner | **Appropriate design question**, but the present model does not prove that its configurable Jet-A/LNG topologies and head-end routes coexist in one physical geometry. |
| Steady state | **Acceptable for stabilized emissions screening.** Not acceptable by itself for ignition, LBO, relight, switching, or thermoacoustic stability. |
| One representative sector | **Conditionally acceptable** when cup count, periodicity, mass-flow scaling, and sector volume are explicit. It cannot predict circumferential maldistribution or engine exit pattern factor. |
| Constant combustor pressure | **Acceptable inside a calibrated chemical submodel**, but not as a substitute for the diffuser/liner pressure-loss and effective-area calculation that determines air split and jet momentum. |
| Prescribed air split and recirculation | **Screening only.** It cannot identify a fixed shared-liner design across LTO and cruise without geometry/pressure coupling. |
| PSR flame zones | **Acceptable only after CFD or rig calibration** of volume, residence-time distribution, and mixing. |
| PFR post-flame zone | **Physically reasonable**, but not implemented in the production network. |
| Adiabatic zones | **Useful upper-temperature sensitivity**, not adequate for quantitative RQL NOx, liner cooling, wall temperature, or pattern-factor design. |
| Pure methane as LNG | **Not adequate** for a representative LNG design; C2/C3/N2 content and phase fractionation matter. |
| Single pseudo-component Jet-A liquid | **Not adequate** for evaporation, ignition, CO/UHC, and soot-sensitive design; useful only for first-order hydraulic screening. |
| Two different fuel mechanisms | **Scientifically necessary**, but cross-fuel absolute metrics need separate validation and uncertainty before comparison. |
| Hot equilibrium initialization | **Useful for locating a lit steady branch**, not evidence of ignition or stability. |
| Minimum ignition margin of four | **A tunable engineering screen, not a universal law.** It needs mechanism uncertainty, passage nonuniformity, wall effects, and rig evidence. |

### 6.6 Documentation and traceability corrections required with the physics fixes

1. Claims that PFR zones are PSR chains, TAB fallback is applied, every LNG
   cooling-to-mixing link is computed, and mechanism validity is enforced at runtime are not
   true for the production evaluator.
2. `docs/modeling.md` retains stale RQL status.
3. The technical reference has no complete bibliography or equation-to-source/validity-range
   map and still contains “Pending” material despite being marked complete elsewhere.
4. The John et al. exercise must be labeled a qualitative model-form comparison, not
   validation. Its unpublished topology, volumes, spray propagation, and calibration cannot
   be reconstructed.
5. The adjusted RQL default that recovers a classical optimum is a calibration/sanity
   choice. Tests using the same model are regression tests, not independent confirmation.
6. Every future result must record mechanism hashes, property backend/version, fuel
   composition, geometry revision, calibration identifiers, numerical settings, warning
   disposition, and uncertainty interval.

---

## 7. Proposed remediation plan — awaiting approval

No item below should be interpreted as implemented. Source changes should begin only after
the plan is approved.

### Stage 0 — protect users from unsupported conclusions

1. Mark the current design evaluator as prototype/screening-only. Preserve distinct
   **pass/fail/unknown** states and make unknown, non-converged, out-of-range, negative-flow,
   and unvalidated safety results ineligible for design acceptance, never safe.
2. Replace the example's design claims with reproducible diagnostics until the end-to-end
   gates in Section 8 pass.
3. Establish an equation/assumption register with source, validity range, uncertainty,
   implementation symbol, verification test, and validation dataset for every
   design-driving relation.

**Exit gate:** the software cannot emit an apparently feasible design when a required
physics or safety calculation was skipped.

### Stage 1 — establish one mission and hardware state model

1. Make `OperatingPoint` the canonical boundary condition and add an explicit adapter for
   the four ICAO LTO points and representative cruise points from an engine cycle deck.
2. Separate compressor discharge, dome, liner-zone, nozzle-inlet, critical, and chamber
   pressures. Close pump/feed/nozzle and diffuser/liner pressure budgets.
3. Define engine, annular, cup, and sector flow/volume scaling unambiguously.
4. Define one shared physical liner and separate Jet-A and LNG fuel passages. Permit
   fuel-specific head-end routing and schedules without silently changing shared geometry.
5. Derive air splits at every operating point from effective areas, discharge coefficients,
   local density, and pressure ratio; retain prescribed splits only as a documented
   calibration mode.
6. Carry representative LNG composition and phase compositions through CoolProp, Cantera,
   stoichiometry, thermal, spray, and emissions interfaces with an explicit species-name map.

**Exit gate:** every mass flow, pressure, temperature, composition, geometry, and sector
volume has one source and can be traced from mission input to CRN boundary.

### Stage 2 — correct reactor topology, conservation, and steady solution

1. Expand every `ReactorKind.PFR` into the requested PSR chain, preserving total physical
   volume, or use an independently cross-checked spatial PFR where network coupling permits.
2. Rework reactor formulation so the simulated gas inventory equals
   `rho_converged * physical_volume`; report actual Cantera mass and volume and reconcile
   both with the residence-time definition.
3. Construct initial compositions from species molar flows, not mass-flow-weighted mole
   fractions, and initialize each reactor from a physically defined local state.
4. Construct internal flows after all gas and droplet sources are known. Preserve
   physically specified recirculation and staging rather than using fuel addition as a flow
   “measurement error.”
5. Replace unconstrained minimum-norm repair with uncertainty-weighted, nonnegative,
   component-wise closure. Reject topology reversal and unresolved residuals.
6. Converge temperature, all relevant species, reactor derivatives, mass, elements, and
   energy. Cross-check time marching against Cantera steady solutions on tractable networks.
7. Add hot/cold continuation and branch tracking for PSR extinction; keep ignition,
   flameholding, and blowout conclusions distinct.
8. Route cooling to its declared destination and replace free heat-loss watts with a
   traceable wall/coolant model or an explicitly calibrated boundary.

**Exit gate:** the exact production topology passes the analytical and conservation tests in
Section 8 with no negative flow and grid-independent PFR/quench results.

### Stage 3 — make nozzle, phase, spray, and CRN coupling conservative

1. Feed the actual Jet-A and LNG nozzle results into every design evaluation; remove the
   all-vapor shortcut from design acceptance.
2. Replace the universal Jet-A discharge coefficient with a hardware/air-core closure and
   validate liquid-sheet plus gas-driven breakup across the full LTO density range.
3. Correct the damped TAB solution and invoke breakup when requested. Keep flash breakup and
   aerodynamic breakup separate and calibrated.
4. Add drag/slip relaxation and local gas velocity; propagate size, temperature, velocity,
   and trajectory by class and zone.
5. Integrate liquid mass directly, or include the density derivative consistently, so
   variable-density heating cannot create or destroy liquid.
6. Use local vapor composition in evaporation and close vapor inhibition, condensation
   policy, and high-pressure applicability.
7. Replace endpoint heat reconstruction with integrated gas-plus-liquid species, mass, and
   enthalpy ledgers, including carryover.
8. Add multicomponent Jet-A evaporation suitable for the required fidelity or explicitly
   gate the single-component approximation to screening.
9. Preserve LNG vapor and liquid phase compositions and use a flash/evaporation formulation
   driven by liquid enthalpy as well as gas convection.
10. Close choked nozzle and CRN boundary mass flux, density, velocity, pressure, quality, and
   area. Calibrate finite-rate flashing before making an `L/D` recommendation.
11. Add wall impingement/carryover bounds or declare geometries requiring those physics
   outside the model domain.

**Exit gate:** nozzle-to-exhaust mass, each element, and total enthalpy close within the
proposed tolerances for zero, partial, and complete evaporation cases.

### Stage 4 — correct LNG thermal management

1. Calculate tank-to-injector duty along the actual feed-pressure/heat-exchanger path, then
   perform the nozzle pressure/enthalpy flash separately.
2. Use equilibrium quality and phase compositions throughout the mixture bubble/dew glide;
   do not infer quality from a pure-fluid temperature test.
3. Couple available heat to effectiveness, pinch, pressure loss, mass, control schedule, and
   off-design heat source rather than a single scalar ceiling.
4. Include actual autoignition state in the thermal window and treat an infeasible inverse
   target as such rather than clamping it.
5. Integrate idle-path thermal soak, purge, coking/deposit, vapor-lock, and restart screens
   into mission feasibility.

**Exit gate:** thermal-window points reproduce independent enthalpy/flash calculations and
remain feasible at every required mission condition and uncertainty bound.

### Stage 5 — establish chemistry and operability credibility

1. Invoke mechanism validation before every solve and attach mechanism identity, range,
   pathway, and extrapolation status to the result.
2. Select and validate pressure-appropriate C1-C3 LNG chemistry and a representative LNG
   composition; use GRI-Mech 3.0 only where its evidence supports it.
3. Validate the Jet-A network mechanism against its source flame/speciation/NOx data and
   bracket low-temperature ignition with defensible Jet-A mechanisms/data rather than one
   hypothetical fast variant.
4. Use multi-dimensional, bounded ignition interpolation in log delay and preserve separate
   states for interpolated, censored-lower-bound, and unavailable results. Test multiple
   ignition markers.
5. Compute laminar flame speed where supported and calibrate flashback/turbulent-flame-speed
   treatment to the actual passage.
6. Implement a continuation-based LNG extinction/efficiency/CO bracket as a CRN lean-limit
   screen, then calibrate it to high-pressure rig data.
7. Add the pollutant fidelity needed for each claim. Do not claim quantitative CO, UHC,
   soot, or nvPM until suitable chemistry/modeling and data exist.
8. Keep transient ignition, relight, switching, and thermoacoustics as explicit external
   gates unless a validated transient model is added.

**Exit gate:** the lower confidence bound, not only the nominal mechanism, satisfies
autoignition/flashback/LBO criteria, and emissions errors meet held-out acceptance targets.

### Stage 6 — replace proxies and optimize only validated quantities

1. Optimize Jet-A LTO Dp/Foo using all four modes and fuel-flow/time weighting. Retain
   per-mode EI and constraint results so idle, approach, climb-out, and takeoff tradeoffs
   remain visible.
2. Define LNG cruise performance at named cruise points, not a generic “LNG” average.
3. Replace serial-reactor phi range with a calibrated unmixedness/PDF metric. Replace
   all-reactor temperature range with mass-weighted parallel exit-stream statistics or a
   CFD/rig exit pattern factor.
4. Make pressure loss, combustion efficiency, temperatures, cooling, carryover, nozzle
   hydraulics, autoignition, flashback, idle thermal state, packaging, and operability
   explicit constraints.
5. Use nondimensional normalization only for a documented scalar preference; preserve the
   Pareto front as the primary multi-objective result.
6. Make sweep bracketing aware of the varied coordinate and sort/group samples before
   endpoint tests. Use adaptive boundary expansion for single-fuel optima.
7. Propagate input, calibration, mechanism, manufacturing, numerical, and model-form
   uncertainty. Rank designs by robust feasibility and confidence intervals.
8. Compare the shared-liner design against separately optimized Jet-A and LNG liners only
   after all three designs use the same validation and uncertainty rules.

**Exit gate:** the recommended design remains feasible and non-dominated under uncertainty,
mechanism substitution, numerical refinement, and held-out validation cases.

### Stage 7 — release a conceptual-design result

1. Run the full LTO/cruise matrix, off-design corners, fuel-switch/idle-path external gates,
   and uncertainty ensemble.
2. Report dimensions, effective areas, pressure loss, flows, volumes, residence
   distributions, injector ranges, thermal schedule, constraint margins, Pareto position,
   and evidence grade.
3. State which quantities are verified, validated, calibrated, extrapolated, or unavailable.
4. Require independent technical review before lifting the no-go disposition.

---

## 8. Proposed verification and validation program

Acceptance values below are initial engineering gates, not substitutes for dataset-specific
measurement uncertainty. They should be finalized before calibration and never relaxed after
seeing holdout results.

### 8.1 Level V0 — equation, dimensional, and limit tests

| Area | Required independent test |
|---|---|
| Units and dimensions | Check every public equation for dimensional homogeneity; test unit conversions and scale invariance. Use the existing unit package rather than duplicating formulas in tests. |
| Pressure and flow | Verify SPI against a hand solution; recover the incompressible limit; demonstrate a back-pressure choking plateau; refine HEM pressure stations; require CFD-boundary `mdot = Cd * rho * u * A` consistency under the declared convention. |
| Feed thermal/hydraulic | Zero-heat, heat-only, friction-only, hydrostatic-only, and minor-loss limits; verify `delta h = Qdot / mdot`; refine axial segmentation; compare two-phase pressure loss with an independent implementation. |
| Thermodynamics | Compare methane and representative-LNG saturation, density, enthalpy, entropy, and PH/PS flash against NIST REFPROP/GERG-quality references, including points inside the bubble/dew glide and near critical conditions. |
| Reactor algebra | Zero-fuel pass-through, long-PSR equilibrium, residence-time/inventory identity, analytical first-order PSR/PFR cases, disconnected-graph rejection, and exact nonnegative mass closure. |
| PFR | Exercise `ReactorKind.PFR` itself; compare its segment chain with Cantera's independent PFR example and refine segment count. |
| Droplets | D2-law, zero-relative-velocity, zero-vapor-driving-force, saturation inhibition, condensation policy, damped-TAB under/critical/overdamped limits, drag relaxation, and class mass/number conservation. |
| Coupled energy | Nonreacting zero/partial/full evaporation with independent inlet/outlet gas-plus-liquid enthalpy; reacting adiabatic and prescribed-wall-heat cases; include carryover. |
| Autoignition | Exact grid points, all grid faces/corners, T/P/phi refinement, out-of-range state, censored nonignition, and multiple ignition markers. |
| Objectives | Hand-calculated four-mode Dp/Foo with unequal fuel flows; mass-weighted parallel exit statistics; invariance to unit rescaling; unordered/multidimensional sample bracketing. |

The Jet-A injector suite must also include a discharge-coefficient/air-core benchmark against
the cited pressure-swirl case and an LTO gas-density sweep. The latter must show the expected
external-breakup/SMD pressure response rather than a pressure-independent calibrated value.

Proposed numerical gates:

- algebraic identities: relative error below `1e-10`;
- full mass and elemental residual: below `1e-8` relative;
- coupled enthalpy residual: below `0.1%`;
- no negative or silently omitted flow;
- HEM, feed, Tier-2, PFR, quench, and coupling refinement: below `0.5%` change
  in hydraulic/quality outputs, `1%` in EI-NOx, and `2 K` in temperature; and
- ignition interpolation error below `5%` against directly integrated points.

Residuals must use declared scales. Use total positive inlet mass flow for total-mass
normalization, incoming elemental mass flow for each materially present element, and the sum
of absolute inlet/outlet/heat energy fluxes for enthalpy normalization. For a zero or trace
quantity, use a predeclared absolute tolerance based on solver precision instead of dividing
by a vanishing reference. Record both normalized and dimensional residuals.

### 8.2 Level V1 — component validation

1. **Properties:** methane reference EOS plus representative LNG mixtures. Initial targets:
   pure-fluid density/saturation/enthalpy within `0.5%`, mixture thermodynamics within `2%`,
   and transport within `10%`, unless source uncertainty is larger.
2. **Critical flow:** ingest the bundled Hammer spreadsheet and the De Lorenzo and
   Kim-O'Neal cases already identified in `V_AND_V_ROADMAP.md`. These validate model form;
   they are not substitutes for LNG data.
3. **Cryogenic spray:** use the identified liquid-nitrogen flash-spray datasets first, then
   obtain direct liquid-methane or representative-LNG onset, mass-flow, cone, and size data
   for a geometrically similar injector.
4. **Jet-A nozzle/spray:** use independent pressure-swirl flow number, air-core/film, cone,
   SMD/distribution, cavitation, and turndown data across LTO pressure/density conditions.
5. **Jet-A evaporation:** validate multicomponent vaporization and composition history
   against single-droplet or spray data at elevated pressure and temperature.
6. **Jet-A chemistry:** replay the source HyChem POSF10325 ignition, flame-speed,
   speciation, and NOx datasets, preserving train/holdout separation.
7. **LNG chemistry:** validate the selected mechanism and representative compositions
   against high-pressure ignition delay, laminar flame speed, extinction, JSR/PSR species,
   and NOx/CO data over the cruise domain.
8. **Mixing and cooling:** validate quench/dilution jet penetration and mixing with
   Holdeman-type confined-crossflow data and cold-flow sector measurements; validate wall
   heat/cooling against thermal data.

Initial holdout targets:

- nozzle mass flow within `5%`;
- critical mass flux within `5-10%`;
- cone angle within `5 deg` or measurement uncertainty;
- SMD within `20%` or measurement uncertainty, with uncertainty bands achieving their
  declared coverage;
- ignition delay within a factor of two and no false-safe classification; and
- PSR/JSR species and emissions within experiment-specific uncertainty targets.

### 8.3 Level V2 — subsystem CRN validation

1. Reconstruct published PSR/PFR/CRN cases without tuning their topology after inspecting the
   answer.
2. Calibrate zone volumes, flows, heat loss, and mixing to one CFD/rig subset; freeze them
   before evaluating held-out operating conditions, traverses, and emissions.
3. Use the John et al. liquid/gaseous comparison only as a qualitative model-form check
   unless their topology and calibration data become available.
4. Validate hot/cold branches, extinction residence time, and mechanism substitution.
5. Compare local temperature/phi PDFs, recirculation, residence-time distributions,
   pressure loss, exit traverse, and EI rather than calibrating only a single exit number.

### 8.4 Level V3 — mission and hardware validation

#### Jet-A LTO

1. Use an engine cycle deck to define P3, T3, air flow, fuel flow, fuel temperature, pressure
   loss, and scheduled geometry at takeoff, climb-out, approach, and idle.
2. Validate sector-rig RQL/LDI pressure loss, pattern factor, efficiency, NOx, CO/UHC, LBO,
   and wall temperatures at held-out conditions.
3. Compare complete four-mode EI and Dp/Foo trends with the ICAO Aircraft Engine Emissions
   Databank for an appropriate technology class. This is a model validation benchmark, not
   a certification claim.

#### LNG cruise

1. Define top-of-climb, nominal cruise, and end-of-cruise points, plus fuel-composition and
   heat-source envelopes.
2. Validate representative-LNG flame stability, combustion efficiency, NOx/CO, pressure
   loss, flash state, and exit traverse at pressure-matched conditions.
3. Use PRECCINSTA or similar methane data only as a topology/trend check; a pressure- and
   geometry-matched sector rig remains the acceptance evidence.

Proposed system gates are exit temperature within `2%`, EI-NOx within `20%` or two combined
standard uncertainties, and lean-limit phi within `0.03`, subject to tighter or looser
dataset-specific uncertainty. For held-out pairwise design comparisons, require the predicted
ordering of each named objective to match the measured ordering when the experimental 95%
intervals do not overlap; treat overlapping intervals as a tie. Require Pareto dominance to
use those same named objectives and uncertainty rule. An average error cannot hide a wrong
mode trend.

### 8.5 Level V4 — uncertainty, robustness, and credibility

1. Separate numerical, input, parameter, mechanism, measurement, and model-form uncertainty
   following ASME V&V 20 / NASA-STD-7009 principles.
2. Publish calibration priors/ranges and posterior coverage; do not calibrate and validate on
   the same operating points or geometry.
3. Run global sensitivity and uncertainty propagation over composition, P3/T3, flow,
   effective area, heat transfer, residence/mixing, mechanism, spray, and manufacturing
   variables.
4. Require safety margins at a predeclared one-sided 95% adverse confidence or credible
   bound, with the statistical construction recorded, and require the preferred design to
   remain non-dominated under plausible model substitution.
5. Add continuous-integration regression cases for every accepted validation case, with raw
   data checksum, source, preprocessing, software versions, and pass/fail metric.
6. Use independent review to approve any change in model-form choice or acceptance threshold.

### 8.6 Reference anchors for the remediation and V&V work

The implementation should trace each equation to the exact edition/page or archival data
record. The initial reference set is:

- Cantera `continuous_reactor.py`, `combustor.py`, and `pfr.py` examples for CSTR residence,
  PSR extinction continuation, and a reactor-chain PFR;
- Lefebvre and Ballal, *Gas Turbine Combustion*, for architecture, pressure loss, cooling,
  loading, and operability context;
- Holdeman, “Mixing of Multiple Jets With a Confined Subsonic Crossflow,” for
  quench/dilution mixing similarity;
- Abramzon and Sirignano, “Droplet Vaporization Model for Spray Combustion Calculations,”
  for a higher-fidelity evaporation reference;
- Downar-Zapolski et al., “The Non-Equilibrium Relaxation Model for One-Dimensional Flashing
  Liquid Flow,” for finite-rate flashing context;
- Wang et al. and Xu et al., HyChem Parts I and II, plus Saggese et al., “HyChem V - NOx
  Formation from a Typical Jet A,” for Jet-A chemistry provenance and validation;
- Glarborg et al., *Progress in Energy and Combustion Science* 67 (2018), for nitrogen
  chemistry assessment;
- GRI-Mech 3.0 primary documentation and pressure/composition validation basis;
- ICAO Annex 16 Volume II and the Aircraft Engine Emissions Databank for LTO definition and
  comparison data; and
- ASME V&V 20, AIAA G-077, Oberkampf and Roy, and NASA-STD-7009 for credibility and
  uncertainty practice.

---

## 9. Approval gate

The recommended implementation order is Stages 0-2 first, because unsafe status handling,
mission-state closure, physical volume/residence, PFR topology, and conservative flow
assembly can invalidate every later result. Stages 3-5 then establish the nozzle/spray,
thermal, chemistry, and operability evaluator. Optimization in Stage 6 should not be used to
select hardware until the relevant Section 8 component and subsystem gates pass.

**No source implementation should begin until this audit and staged plan are approved.**

---

## 10. Approved remediation implementation log (started 2026-08-30)

The user approved Stages 0--2 followed by the remaining stages. This section is the running
record of implementation details, verification, and deviations. Entries are appended as work
is completed; a stage is not complete merely because code exists.

### 10.1 Stage 0 -- fail-closed screening status (complete)

Implemented so far:

- `AutoignitionVerdict.UNKNOWN` distinguishes an unavailable ignition calculation from a
  safe one. `autoignition_margin()` now emits an error and makes the design ineligible for
  acceptance when the table cannot supply a delay.
- `evaluate_objectives()` creates an explicit constraint violation for unknown
  autoignition. The prototype evaluator and example now call their outputs screening
  diagnostics rather than validated design answers.
- `crn/status.py` provides general `PASS`, `FAIL`, and `UNKNOWN` gates. `PointResult` and
  `DesignResult` now report computational and acceptance status separately. The legacy
  `feasible` property is true only for acceptance `PASS`; the current prevaporized prototype
  remains `UNKNOWN` even after a successful numerical solve.
- `docs/CRN_EQUATION_REGISTER.md` records each design-driving relation, its implementation,
  applicability, verification, validation evidence, and current acceptance status.
- A regression test exercises an above-table ignition state and requires fail-closed
  behavior.

Deviations and constraints:

- Evidence grades and the external release gate remain in Stage 7; the general tri-state
  acceptance mechanism itself is now implemented.
- “Unknown” is not treated as physical failure. It is a distinct state that is ineligible
  for design acceptance until evidence resolves it.

Verification at closure: 61 focused tests passed, Ruff passed, and the changed files passed
secret scanning. Commit `4bdf83f` records the stage.

### 10.2 Stage 1 -- canonical mission state and shared hardware (complete)

Implemented:

- `OperatingPoint` remains the canonical user input. `resolve_pressure_stations()` now names
  compressor-discharge, dome, combustor-exit, and both fuel-nozzle inlet pressures, applies the
  declared liner pressure-loss fraction, and rejects an active-fuel pump that cannot supply the
  requested chamber pressure plus nozzle drop.
- `mission_point_from_operating()` is the only adapter from canonical operating inputs into the
  CRN `MissionPoint`. It retains the source point and resolved stations for traceability.
  `MissionProfile.from_icao_lto()` requires exactly one Jet-A takeoff, climb-out, approach, and
  idle point, including explicit durations; `from_cruise()` requires one or more LNG points.
- `SectorDefinition` explicitly converts engine-total flow to the modeled cup sector and back.
  The adapter applies that scale to both fuel and combustor air, preserving air/fuel ratio.
- `DualFuelHardware` installs separate Jet-A and LNG passage geometries in one immutable
  `SharedLinerGeometry`. The passage-air share is calculated from fixed effective areas instead
  of remaining an independent optimizer choice. Zone volumes and cooling destination likewise
  come from the shared liner when hardware is supplied.
- `AirAdmission.mass_flow_kg_s()` implements the compressible ideal-gas orifice relation with
  the subcritical and choked branches. `SharedLinerGeometry.area_derived_split()` normalizes
  flows through the five fixed effective areas. Prescribed fractions remain available only as a
  named calibration mode with a non-empty calibration identifier.
- `DesignEvaluator` accepts `MissionProfile`, canonical `OperatingPoint` objects, or legacy
  `MissionPoint` objects. Area-derived hardware replaces optimized station fractions; legacy
  direct mission input remains supported but receives an `UNKNOWN` canonical-hardware gate.
- `LNGComposition.cantera_mole_fractions()` maps CoolProp fluid names to mechanism species.
  `MechanismRegistry.with_lng_composition()` propagates that declared composition into both
  LNG mechanism roles. Unsupported components are rejected rather than silently discarded.
- New regression tests cover pressure closure and pump deficit, sector scaling, effective-area
  air splits, calibration traceability, passage shares, mission ordering/completeness, and LNG
  composition propagation/rejection.

Assumptions and limitations retained deliberately:

- The air-hole equation assumes steady, one-dimensional, isentropic, calorically perfect air
  with constant `gamma=1.4` and `R=287.05 J/(kg K)`. A discharge coefficient represents losses.
  All stations presently use the same prescribed downstream pressure. The model does not
  resolve liner axial pressure, swirl, crossflow jet penetration, or hole-to-hole interaction.
- Air-hole flows are normalized into fractions. They do not yet enforce that the sum of
  absolute hole flows equals the canonical combustor air flow, nor do they solve liner pressure
  loss from geometry. Consequently this is a more physical parameterization, not a validated
  pressure-loss model; acceptance remains `UNKNOWN`.
- Sector scaling assumes nominally identical cups and that canonical flows are engine totals.
  Distortion and circumferential maldistribution require external data.
- Real liner effective areas, volumes, discharge coefficients, cycle-deck states, and
  representative LNG compositions are unavailable. No values were invented.

Deviations from the approved plan:

- Stage 1 introduces the complete hardware and mission interfaces but cannot calibrate or
  validate them without the external geometry and cycle deck assigned to later evidence gates.
- Legacy `MissionPoint` and prescribed air fractions remain as compatibility and explicitly
  traceable calibration paths. They are not accepted as canonical design evidence.

Verification at closure: 67 focused Stage 1/chemistry/design tests passed and Ruff passed.

### 10.3 Stage 2 -- topology, conservation, inventory, and steady branches (complete)

Implemented:

- `ReactorKind.PFR` now expands during `CombustorNetwork` assembly into the requested
  equal-volume PSR chain. Incoming flows and ordinary inlets enter the first segment, the
  outlet and outgoing flows leave the last, heat loss is divided without changing its total,
  and cooling declared for the physical exit joins the last segment. The final segment keeps
  the physical-zone name, while `NetworkSolution.zone()` exposes every numerical segment.
- The historical `minimum_norm_mass_correction()` API now performs nonnegative,
  uncertainty-weighted closure. It checks global and connected-component boundary balances,
  preserves explicitly fixed recirculation flows, rejects negative inputs and infeasible
  directed graphs, and rejects underdetermined repairs unless every adjustable flow has a
  declared uncertainty. `close_internal_flows` is the preferred descriptive alias.
- Fuel vapor and staged droplet sources are included before closure. Architecture templates
  label their recirculation edges as fixed, so adding fuel changes the downstream through-flow
  rather than redistributing the fuel residual into a recirculation loop. Recirculation ratio
  is now based on total air-plus-fuel boundary flow rather than air alone.
- `_initial_mixture()` converts each inlet mass flow to species molar flows before combining
  mole fractions. It also computes the mass-flow-weighted inlet enthalpy used by cold and hot
  seeds. The earlier mass-flow-weighting of mole fractions was dimensionally incorrect.
- A constant-pressure Cantera reactor otherwise preserves its initial inventory and changes
  volume as temperature and composition change. The solver now re-imposes each declared
  physical control volume during time marching and reinitializes Cantera, so final reported
  `mass_kg`, `volume_m3`, and residence time obey `mass = density * volume` and
  `tau = mass / inflow`.
- Convergence now requires temperature, every species mass fraction, mass inventory, enthalpy
  inventory, physical volume, and a scaled state derivative to settle. Final element and
  steady enthalpy balances are computed and can invalidate convergence.
- `InitializationBranch`, `solve_branches()`, and `solve_continuation()` keep cold/unlit and
  hot/lit solutions distinct and allow a named branch to seed the next same-topology point.
  `NETWORK_EXTINGUISHED` remains a combustion-state diagnosis, not a lean-blowout or
  flameholding claim.
- `CoolingAirDestination` is now honored: primary cooling rejoins the primary/mixer zone,
  dilution cooling enters the start of the post-flame PFR, and exit cooling enters its last
  segment. `ReactorSpec` rejects prescribed heat-loss watts without a calibration or
  physical-model identifier; spray sensible/latent heat carries an explicit computed basis.
- Regression tests exercise exact PFR expansion and grid refinement, nonnegative/fixed and
  component-wise closure, fuel-source propagation, molar mixing, inventory/residence
  identities, cooling destinations, energy/element closure, heat-loss traceability, and
  hot/cold continuation.

Assumptions and limitations retained deliberately:

- A PFR remains a numerical chain of ideal stirred cells, not a resolved velocity profile.
  Segment refinement is required for each design-driving result.
- Uncertainty-weighted closure does not infer uncertainty. Underdetermined data-driven
  networks must supply it; exact template chains close uniquely after fixed recycle flows are
  declared.
- The control-volume reinitialization is a steady-state inventory constraint for Cantera's
  constant-pressure zero-dimensional reactor. It is not intended to represent physical
  transient volume forcing.
- Hot initialization uses the equilibrium state at the aggregate inlet enthalpy. Cold
  initialization uses the corresponding unreacted mixed state. These locate possible steady
  branches but do not prove real flameholding or lean blowout.
- Prescribed wall heat is now traceable but remains calibration-only. Stage 4 must replace it
  with a wall/coolant heat-transfer closure before thermal acceptance.

Implementation deviations:

- The public name `minimum_norm_mass_correction()` remains for compatibility, although its
  algorithm is no longer unconstrained minimum norm. New code should use
  `close_internal_flows`.
- Continuation requires identical expanded reactor names between adjacent points. A segment
  count or topology change is rejected clearly rather than interpolating states between
  unlike control volumes.

Verification at closure: 82 focused network, template, coupling, and design tests passed and
Ruff passed. The earlier long-running full physics suite was stopped under the user's time
limit; that deviation did not waive any Stage 7 validation gate.

### 10.4 Stage 3 -- conservative nozzle, phase, spray, and CRN coupling (implemented;
external validation open)

Implemented:

- `DesignEvaluator` now accepts an explicit `spray_model` callback. When present,
  `evaluate_point()` uses the returned `SprayBoundary` and liquid-property provider in
  `solve_coupled()` rather than inserting all fuel as vapor. Fuel identity and nozzle mass
  flow must match the mission point. Omitting the callback retains the quick gas-only
  diagnostic, but the model-fidelity gate remains `UNKNOWN`.
- `ThermodynamicState`, `CFDSprayBoundary`, and `SprayBoundary` preserve equilibrium liquid
  and vapor mole fractions. The LNG nozzle bridge maps each CoolProp component to the
  corresponding Cantera species before the phase streams enter the CRN. The finite-rate
  quality remains the mass split; equilibrium phase compositions are identified as such.
- `taylor_analogy_breakup()` now includes exponential viscous damping in the forced oscillator
  and finds the first threshold crossing with a bounded root solve. `_march_droplets()` invokes
  it only when `SprayBoundary.apply_aerodynamic_breakup` requests it. Internally driven LNG
  flash breakup remains separate.
- `integrate_droplet()` advances liquid mass directly, then obtains radius from current mass,
  density, and temperature. This removes the variable-density error caused by inferring mass
  from radius cubed. The same ODE advances droplet velocity with Schiller--Naumann drag and
  integrates gas-to-liquid heat directly.
- `_march_droplets()` estimates local gas velocity from declared spray-path length and gas
  residence time, propagates droplet velocity between zones, calculates local fuel-vapor mass
  fraction from the reactor composition, and therefore activates vapor-inhibition physics.
  Its gas heat sink is the integrated convective heat for every droplet class, not an endpoint
  latent/sensible reconstruction.
- `CoupledSolution` reports the nozzle-vapor, evaporated-fuel, and liquid-carryover mass ledger,
  its residual, and integrated droplet heat. `DesignEvaluator` adds explicit mass-closure and
  carryover gates and uses gas exhaust flow, excluding carryover, in emission-index arithmetic.

Verification:

- Existing coupling and spray tests passed (28 tests), and the updated droplet, thermal, and
  evidence-focused group passed (70 tests total). New checks verify direct liquid-mass closure,
  drag relaxation, integrated heat, phase glide, and coupled fuel closure.

Deviations and open evidence:

- Real installed nozzle geometry, Jet-A sheet/air-core calibration, LNG finite-rate flash
  calibration, high-pressure breakup data, chamber wall geometry, and impingement data were
  not available. No values were invented. Uncalibrated size, flash, single-surrogate Jet-A
  evaporation, dense-gas transfer, and possible wall impingement remain `UNKNOWN` for design
  acceptance.
- A complete multicomponent Jet-A distillation model was not added because no validated
  composition/property basis exists in the repository. The current single-surrogate model is
  explicitly screening-only.
- Trajectory is one-dimensional along a declared path; transverse position and wall films
  remain outside the model domain.

### 10.5 Stage 4 -- LNG thermal management (implemented; hardware thermal evidence open)

Implemented:

- `feed_path_heat_budget()` checks heat against enthalpy rise on the actual segmented
  pressure/enthalpy path produced by `solve_lng_feed_line()`. Every `FeedLinePoint` now retains
  equilibrium liquid and vapor compositions as well as quality.
- `heat_sink_budget()` obtains the target state from a pressure-temperature flash and divides
  its enthalpy rise at the bubble and dew enthalpies. A multicomponent target inside the
  temperature glide therefore has a quality between zero and one instead of being forced to
  all-liquid or all-vapor.
- `thermal_window()` calculates upstream heating duty at feed pressure. Chamber superheat is
  still evaluated at chamber pressure, while the nozzle flash is correctly understood as a
  later isenthalpic pressure drop.
- `HeatSourceModel` bounds recoverable heat by declared source mass flow, heat capacity,
  inlet temperature, exchanger effectiveness, minimum pinch, and fuel-side pressure loss.
  A missing evidence identifier makes the source status `UNKNOWN`; the old scalar ceiling is
  retained only as an unvalidated compatibility input.
- `fuel_temperature_for_target_superheat()` and `solve_temperature_for_duty()` now raise
  `InfeasibleThermalTargetError` when a target crosses feed saturation or lies outside the
  inverse bracket. They no longer silently return a bracket boundary.
- `ThermalWindowPoint` carries autoignition and heat-source statuses. A supplied
  `autoignition_evaluator` is evaluated at every candidate temperature.
- `idle_circuit_screen()` now distinguishes the steady wall-temperature screen from transient
  soak, purge, and restart evidence. `DesignEvaluator` attaches this screen when wall
  temperature is present; missing transient data remains an explicit `UNKNOWN` gate.

Deviations and open evidence:

- A specific aircraft/engine heat exchanger and off-design heat-rejection map were not
  available. The code accepts the required physical inputs but cannot create a control schedule
  or calibration from no data.
- Deposit accumulation, actual purge volume, soak-back transient, and restart performance need
  fuel-system geometry and rig evidence. The implemented inputs fail closed; they are not a
  fitted coking model.
- Equilibrium phase compositions are retained. Nonequilibrium phase fractionation inside the
  finite-rate nozzle still requires dedicated validation.

### 10.6 Stage 5 -- chemistry and operability credibility (implemented; validation gates open)

Implemented:

- `evaluate_point()` calls `validate_mechanism()` before a production solve and attaches
  mechanism path and provenance to `PointResult`. Separate applicability and held-out
  validation gates prevent a structurally loadable mechanism from being mistaken for a
  pressure-appropriate validated one.
- `IgnitionDelayTable.evaluate()` performs bounded interpolation in inverse temperature,
  pressure, and equivalence ratio on logarithmic delay. `IgnitionEvidenceState` separates a
  finite interpolated result, a no-ignition lower bound censored at the integration window, and
  an out-of-domain/unavailable result. A censored bound proves safety only when its lower-bound
  margin clears the requested minimum.
- `IgnitionMarker` supports both fixed temperature rise and maximum temperature-rate timing.
  Marker sensitivity remains evidence to report rather than a reason to choose the most
  favorable delay.
- `laminar_flame_speed()` provides a Cantera freely propagating flame calculation where the
  mechanism supports it. `flashback_screen()` now requires a passage-correlation calibration
  identifier before a numerically safe screen can pass an acceptance gate.
- `continuation_lean_limit_screen()` follows the hot CRN branch from rich to lean, reports the
  last lit and first extinguished points, and carries CO only as an uncalibrated diagnostic.
  Its warnings state explicitly that numerical extinction is not physical lean blowout.
- Lean blowout, transient ignition, relight, fuel switching, and thermoacoustics are materialized
  as `UNKNOWN` external gates at every mission point instead of existing only in prose.

Deviations and open evidence:

- No new chemical mechanism or external chemistry dataset was available. GRI-Mech remains an
  LNG baseline only; the Jet-A fast-NTC file remains one hypothetical sensitivity mechanism.
  The missing mechanism bracket and held-out flame, ignition, and emissions validations block
  acceptance.
- Flame speed can be computed, but turbulent flashback correlation constants still require the
  actual passage and rig. CRN extinction can be bracketed, but high-pressure sector calibration
  is required before calling the result LBO.
- Quantitative CO, UHC, soot, and nvPM claims remain unavailable. No transient or thermoacoustic
  model was fabricated.

### 10.7 Stage 6 -- mission objectives and robust feasibility (implemented; recommendation gate
open)

Implemented:

- `DesignResult.weighted_ei_nox()` now weights emission index by fuel mass
  (`fuel flow * duration`), not by time alone. `lto_emissions()` assembles every available
  Jet-A mode with rated thrust through the existing `lto_dp_foo()` arithmetic. Missing rated
  thrust or any of the four ICAO modes creates an objective constraint rather than a
  certification-looking number.
- `lng_cruise_ei_by_point()` and `ObjectiveVector.named_metrics` retain every named LNG cruise
  result. `MissionProfile.from_cruise()` rejects blank or duplicate point names.
- The serial-reactor equivalence-ratio range and all-reactor temperature range are no longer
  optimizer objectives. Their compatibility fields are `NaN`. A calibrated mixture-fraction
  PDF/unmixedness input and a mass-weighted parallel exit traverse remain unavailable, so the
  optimizer is not allowed to minimize misleading proxies.
- Every explicit failed point gate becomes a constraint. The primary objective set is Jet-A
  four-mode Dp/Foo and named-cruise LNG NOx. `rank_key()` uses the caller's declared
  lexicographic preference only as a tie breaker and does not sum incompatible units; the
  Pareto front remains primary.
- `cost_of_shared_liner()` identifies a one-coordinate sweep, sorts samples by that physical
  coordinate, and checks endpoints there. It no longer mistakes random insertion order for a
  swept range.
- `evaluate_uncertainty_ensemble()` accepts explicit input, calibration, mechanism,
  manufacturing, numerical, and model-form cases, returns empirical objective intervals, and
  reports robust `PASS` only when every category is covered and every case is feasible.

Deviations and open evidence:

- The code provides uncertainty propagation, not uncertainty values. Required distributions,
  correlations, mechanism alternatives, manufacturing tolerances, and model-form cases must
  come from data and engineering ownership.
- Adaptive range expansion was not made automatic because changing hardware bounds without an
  approved design envelope can be unsafe. An unbracketed optimum is reported and the caller
  must authorize expanded bounds.
- No calibrated unmixedness model, exit traverse, combustion-efficiency model, liner thermal
  limit, or validated pressure-loss closure was available. Their gates remain open rather than
  being replaced by CRN-internal spreads.

### 10.8 Stage 7 -- evidence-graded conceptual release (implemented; disposition NO-GO)

Implemented:

- `EvidenceGrade` and `EvidenceRecord` distinguish verified software, calibrated inputs,
  held-out validation, extrapolation, unavailable evidence, and missing calculations.
  Verification or calibration alone does not pass a predictive validation gate.
- `build_conceptual_design_report()` returns one schema for dimensions, effective areas,
  pressure loss, mission flows, volumes, residence times, injector ranges, thermal schedule,
  constraint gates, Pareto position, uncertainty status, evidence, and independent review.
  Missing geometry and ranges are represented by `None`, never a plausible placeholder.
- A release is `GO` only if mission acceptance, every required validation record, all
  uncertainty categories, and independent technical review pass. `UNKNOWN` is operationally
  `NO_GO`.

Current release disposition:

- **NO-GO.** The repository lacks the engine cycle deck, installed dimensions/effective-area
  calibration, nozzle and spray rig data, representative high-pressure LNG/Jet-A chemistry and
  emissions holdouts, exit traverse and liner thermal evidence, transient fuel-switch/relight
  evidence, uncertainty inputs, and independent technical review required by Stages 3--7.
- The implemented code makes these absences visible and machine-readable. It does not satisfy
  the external exit gates merely by having software paths for them.

Stage 3--7 focused verification at implementation time: Ruff passed and 70 focused droplet,
thermal, objective, uncertainty, and release tests passed. Additional existing coupling/spray
tests passed (28 tests). Full-suite validation and security checks are recorded below when
completed.
