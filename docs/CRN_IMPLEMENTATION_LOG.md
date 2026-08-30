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
