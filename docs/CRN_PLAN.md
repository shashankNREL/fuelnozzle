# Dual-Fuel Reactor-Network Combustor Model: Plan and Implementation Log

Plan date: 2026-08-28
Target package: `fuelnozzle` (extends the existing nozzle tool; new `fuelnozzle.crn` subpackage)
Status: **APPROVED 2026-08-28. Phase 0 in progress.**

---

## 0. Purpose

**Primary goal (2026-08-28 direction): produce an optimal fuel-nozzle and combustor
design for dual-fuel Jet-A + LNG operation.** That configuration is novel, nobody has a
design rule for it, and the purpose of this tool is to find one. Everything in this plan
is subordinate to that: the reactor network is the *evaluator*, Section 8 is the *design
method*, and the John et al. comparison is a directional check that the evaluator works,
nothing more.

The existing package answers *nozzle* questions: can the LNG flash, how big is the
orifice, what is the Jet-A cone angle. This extension answers the *combustor and design*
questions that follow:

1. **Architecture selection.** For a dual-fuel engine burning **LNG on cruise** and
   **Jet-A on landing/take-off (LTO)**, should the combustor use a **lean-burn**
   strategy (lean premixed/prevaporized, LPP; or lean direct injection, LDI) or a
   **rich-quench-lean** (RQL) strategy? One set of hardware must serve both fuels,
   so the answer is a compromise, and the trade must be quantified, not asserted.
2. **Nozzle design for flash atomization without autoignition.** LNG flashing is
   free atomization energy. Heating the LNG in the feed line increases the flash
   quality and improves the spray, but it also raises the mixture temperature and
   shortens the autoignition delay in any premixing passage. There is a design
   window. This tool must find it.
3. **Cryogenic thermal fuel management.** LNG is a large heat sink. How much heat
   *must* be put into the fuel to get useful flash atomization, how much heat is
   *available* to dump into it from engine/aircraft loads, and what does that heat
   do to combustion, are the same question asked three ways.
4. **The optimal dual-fuel design, and what it costs.** **Separate injector hardware per
   fuel** — they atomize too differently to share — feeding **one shared combustor**, one
   fuel at a time. Each injector is optimized for its own fuel with no compromise; the
   compromise lives entirely in the shared liner air split and zone volumes. Which
   architecture pair, which per-path injector geometry, which air split, and which LNG
   thermal schedule minimize NOx for both fuels while keeping atomization, mixing,
   autoignition margin, idle-circuit thermal limits, and cryogenic feasibility inside
   their bounds — and **how much worse is the shared-liner compromise than a liner
   dedicated to either fuel alone?** That last number is the headline deliverable
   (Section 8.8).

**The two concrete numbers the tool must produce** (2026-08-28 direction, detailed in
Section 8.10):

- **How much air should go to the head end during Jet-A LTO** — the dome air fraction and
  rich-zone equivalence ratio that minimize LTO NOx.
- **How lean can LNG run at cruise** — the minimum stable cruise equivalence ratio, the
  NOx it achieves, and the exit temperature uniformity it delivers.
- **What nozzle geometry to target for each fuel path** (Section 8.11) — orifice sizing,
  flow number, swirl groups for Jet-A; orifice sizing and `L/D` for LNG; premix passage
  dimensions for both. Reported as feasible regions with the driving constraint named.

The first two are questions about the head-end air fraction, asked from opposite ends,
and are coupled through the shared liner — so the real answer there is the triple
`(f_dome, a_J, a_L)`. The third largely decouples, because the paths are separate
hardware.

### 0.1 Confirmed scope decisions

| # | Question | Answer (2026-08-28) | Consequence for this plan |
|---|---|---|---|
| 1 | Are the fuels ever co-fired? | **No. Strictly one fuel per mission segment.** | Exactly one active fuel per operating point. Two separate kinetic mechanisms selected per point. No merged dual-fuel mechanism, no two concurrent spray paths. Data structures stay fuel-generic so co-firing is not *precluded*, but it is not built or tested. |
| 2 | Quantities of interest | **Atomization characteristics, fuel-air mixing, temperature, NOx.** | Calibration targets and acceptance criteria are built around SMD/evaporation, local and exit equivalence ratio, reactor and exit temperature, and NO/NO2. CO, UHC, and soot precursors are *reported as diagnostics* but are **not** calibration targets and carry no accuracy claim. |
| 3 | Will CFD define the reactor zones? | **Eventually yes.** | Phase 1 builds parametric architecture templates driven by user-specified flow splits. Phase 6 builds the k-means clustering + mass-balance + calibration ingestion path from the paper. A **CFD export contract** is written *now* (Section 8) so the CFD runs are set up to produce what the ingestion module will need. |

### 0.2 Additional assumptions (flag if wrong)

- **A1.** Steady state only. No transients, no thermoacoustics, no light-around, no
  altitude relight. Each operating point is an independent steady solve.
- **A2.** One representative fuel injector / one combustor sector is modeled. Multi-cup
  and circumferential non-uniformity are out of scope. `flow_multiplier` semantics
  from `OperatingPoint` carry over unchanged.
- **A3.** Constant pressure inside the combustor. Liner pressure loss is an input that
  sets jet velocities and mixing rates; it is not solved.
- **A4.** The user supplies the kinetic mechanisms as Cantera YAML. The tool validates
  them but does not author or reduce them.
- **A5.** Autoignition is treated as a **screening constraint** on the design, evaluated
  from homogeneous ignition-delay chemistry, not as a calibrated QoI.

---

## 1. Reference basis

**Primary inspiration:** John, P.; Chandrasekhar, H.; Saha, S.; Owoyele, O. (2026).
"A liquid-fueled reactor network model for enhanced NOx prediction in gas turbine
combustors." *Combustion and Flame* 292, 115199. DOI `10.1016/j.combustflame.2026.115199`.
Local copy: `docs/papers/1-s2.0-S0010218026004359-main.pdf`.

What we take from it:

- The core idea: a conventional gaseous CRN cannot predict NOx in a liquid-fueled
  combustor because it homogenizes the equivalence ratio field. Their GFRN
  underpredicted NO by 54–91%; their LFRN reached 0.29–26%. Spray-induced
  equivalence-ratio heterogeneity is the whole mechanism, through the exponential
  temperature sensitivity of the Zeldovich path.
- Two specialized reactor types: an **evaporator/breakup** reactor and a **mixer**
  reactor, with droplets that survive a reactor and convect downstream, carrying
  their own size/temperature history as separate **droplet classes**.
- Sub-models: TAB breakup (their Eqs. 7–13), Frossling evaporation (14–18),
  infinite-liquid-conductivity heating (19–20), vapor-pressure surface composition (21).
- The breakup/evaporation timescale separation that justifies sequential treatment
  (their Fig. 3: breakup ~2e-5 s, evaporation ~6e-3 s, ratio ~133 at 350 K gas and
  ~21 at 900 K gas). This is a directly reproducible verification target.
- k-means partitioning on `[T, phi]`, minimum-norm least-squares mass correction,
  and residence-time calibration (their Eqs. 22–24).

**What we deliberately change:**

| Paper | This tool | Why |
|---|---|---|
| Initial droplet radius = nozzle radius, then TAB breaks it up | Initial droplet size comes from the **existing nozzle models** (`PressureSwirlResult.smd_estimate_m`, `Tier3FlashSpray.smd_estimate_m` / regime), with the paper's r0=r_nozzle+TAB path retained as a selectable fallback | The package already solves the injector. Using r0=r_nozzle throws that away. TAB from a 250 um blob is a crude stand-in for an atomizer model we already have. |
| Jet-A only, diffusion-limited evaporation | Adds a **superheated/boiling branch** for LNG: when the droplet is above its saturation temperature at chamber pressure, vaporization is heat-transfer-limited and the droplet temperature is pinned near T_sat | Flashing LNG does not evaporate by Frossling diffusion. This is the single most important physics addition for the LNG side. |
| Fuel enters entirely as liquid | LNG enters partly as **vapor**: mass fraction `x` from Tier 2/3, only `(1-x)` as droplets | Tier 3 already computes the exit quality. Ignoring it would double-count the atomization work. |
| Topology from CFD k-means | Phase 1: parametric RQL / LPP / LDI / staged templates from user splits. Phase 6: k-means from CFD. | No CFD yet; and architecture *comparison* needs templates you can reason about, not clusters. |
| No autoignition treatment | Explicit premixer autoignition margin module | It is a stated design goal and it is the constraint that decides LPP vs LDI. |

Secondary: Turns, *An Introduction to Combustion*, for PSR/PFR formulation and the
extended Zeldovich mechanism; Lefebvre & Ballal, *Gas Turbine Combustion*, for RQL
and lean-burn architecture conventions and residence-time/loading correlations.

Cantera reactor documentation: <https://cantera.org/dev/reference/reactors/index.html>.
Task 1.1 includes reading this against the installed Cantera version and recording the
exact API surface used, because the reactor/flow-device API changed between 2.6 and 3.x.

---

## 2. How this attaches to what already exists

The new subpackage consumes the existing nozzle solvers rather than duplicating them.

```
OperatingPoint (P3, T3, mdot_fuel, pump pressure, dP_nozzle, fuel temperature)
        |
        +--> solve_lng_feed_line  ------> nozzle inlet enthalpy/temperature
        |          (feed.py)                        |
        |                                           v
        +--> screen_lng_flash / Tier1 / Tier2 / solve_lng_flash_spray
        |          (lng.py, spray.py)               |
        |                                  Tier3FlashSpray:
        |                                    - exit vapor quality x
        |                                    - exit velocity, density, T
        |                                    - regime, SMD (if calibrated)
        |                                    - CFDSprayBoundary
        |
        +--> solve_jet_a_pressure_swirl  -> PressureSwirlResult:
                   (jet_a.py)                  - cone angle, film thickness
                                               - exit velocities
                                               - SMD (if calibrated)
                                               - Re, We, Oh
                                                   |
=====================  NEW: fuelnozzle.crn  ========|===========================
                                                    v
                              SprayBoundary  (fuel-agnostic droplet + vapor source)
                                                    |
                                                    v
              CombustorNetwork  =  air split model + reactor graph + droplet transport
                                                    |
                    +-------------------------------+------------------------------+
                    v                               v                              v
              CombustorResult              AutoignitionScreen             HeatSinkBudget
         (T, phi, NO/NO2, EI, per-reactor    (tau_ign / tau_res margin)    (LNG cooling duty
          droplet/evaporation history)                                     vs heat required)
```

**Nothing in the existing modules is rewritten.** Three small, additive changes are
needed and are listed as explicit tasks (2.4, 2.5, 5.3):

- `JetAPropertyTable` gains optional liquid specific heat, latent heat, molecular
  weight, and boiling point columns (all `None`-able, existing constructors unaffected).
- `CFDSprayBoundary` is reused as-is; no change.
- `OperatingPoint` gains optional combustor-side fields (air mass flow or overall
  equivalence ratio, liner pressure loss, active fuel). All defaulted so existing
  studies and the 21 existing tests keep passing untouched.

---

## 3. Physical model

Every equation below goes into `docs/crn_technical_reference.tex` in full, derived
in plain English for a reader with no combustion background. Sketched here only so
the plan can be reviewed on physics grounds.

### 3.1 Air-side bookkeeping

Total combustor air `mdot_air` follows from either a user-supplied value or a
user-supplied overall equivalence ratio and the stoichiometric air-fuel ratio of the
active fuel. It is then split into user-specified fractions:

```
f_dome (swirler/nozzle air) + f_primary + f_quench + f_dilution + f_cooling = 1
```

Each fraction enters the network at a named reactor. Cooling-film air may be set to
rejoin the flow at a specified downstream station or to bypass to the exit; the choice
changes CO burnout and exit pattern factor and must be explicit, not implicit.

Local equivalence ratio in reactor *j* is computed from element balance (C, H, O)
rather than from an assumed reaction, so it is valid in the presence of unevaporated
fuel and partial products:

```
phi_j = (Z_j / (1 - Z_j)) * (1 - Z_st) / Z_st ,   Z = Bilger mixture fraction
```

### 3.2 Spray boundary (the bridge)

A `SprayBoundary` carries, per fuel:

- `vapor_mass_fraction_at_inlet` — 0 for Jet-A, `x` from Tier 2/3 for LNG;
- `droplet_classes` — initial radius, temperature, velocity, and number flow rate;
- `injection_velocity`, `cone_angle`, `regime`, and full provenance of where each
  number came from (nozzle model, calibration id, or user override).

Number flow rate follows the paper's Eq. 6, applied per class:

```
Ndot_d = mdot_fuel,liquid / ( rho_l * (4/3) * pi * r0^3 )
```

with `r0` chosen by a declared policy:

- `POLICY_NOZZLE_SMD` (default): `r0 = SMD/2` from the nozzle model. Requires a spray
  calibration; if absent, the tool **warns and refuses to report atomization QoIs**,
  matching the existing package's refusal to invent an SMD.
- `POLICY_NOZZLE_RADIUS_TAB` (paper-faithful): `r0 = r_nozzle`, then TAB breakup.
- `POLICY_USER`: explicit user-supplied distribution.

Rosin-Rammler discretization into N classes is supported from the start (the paper
lists this as future work); a single representative class is the default for speed.

### 3.3 Droplet breakup — TAB

Non-dimensional distortion `y = x/(C_b r0)`:

```
d2y/dt2 = (C_f/C_b) * rho_g u_d^2 / (rho_l r0^2)  -  C_k sigma/(rho_l r0^3) y  -  C_d mu_l/(rho_l r0^2) dy/dt
```

with `C_k=8, C_f=1/3, C_d=5, C_b=1/2`. Breakup when `A + We_g > 1`, where

```
We_g = rho_g u_d^2 r0 / sigma
A    = sqrt( (y - (C_f/(C_k C_b)) We_g)^2 + (ydot/omega)^2 )
omega^2 = C_k sigma/(rho_l r0^3) - (C_d mu_l/(2 rho_l r0^2))^2
```

Post-breakup radius (O'Rourke & Amsden), `K = 10/3`:

```
r = r0 / ( 1 + (8K/20) y^2 + (rho_l r0^3 / sigma) ydot^2 (6K-5)/120 )
```

**LNG note.** TAB is a mechanical-breakup model. In the `FULLY_FLASHING` and
`TRANSITIONAL_FLASH` regimes already classified by `spray.py`, breakup is driven by
internal bubble growth, not aerodynamic distortion, and TAB is not applicable. The
tool will **skip TAB for flashing regimes** and take the size directly from the Tier 3
flash-spray result, emitting an INFO warning that records the decision.

### 3.4 Evaporation — two branches

**Branch A, diffusion-limited (Frossling).** Used for Jet-A always, and for LNG when
the droplet is subcooled relative to chamber saturation.

```
dr/dt = - (rho_g D) / (2 rho_l r) * B_d * Sh_d
B_d   = (Y_s - Y_inf) / (1 - Y_s)
Sh_d  = (2 + 0.6 Re_d^0.5 Sc^(1/3)) * ln(1+B_d)/B_d
Tbar  = (T_gas + 2 T_d)/3
rho_g D = 1.293 D0 (Tbar/273)^(nD - 1)      [Jet-A: D0=4.16e-6, nD=1.6]
Y_s   = MW_F / ( MW_F + MW_mix (p_g/p_v(T_d) - 1) )
```

**Branch B, heat-transfer-limited boiling (new).** Used for LNG when
`T_d >= T_sat(P_chamber)`, i.e. a superheated droplet in a flashing spray. The droplet
temperature is pinned at `T_sat` and all convected heat goes to phase change:

```
dm_d/dt = - Q_conv / L_vap(T_sat)
Q_conv  = A_d * beta_spray * Nu_d * k_gas * (T_gas - T_sat) / (2r)
```

Switching between branches is continuous in `dm/dt` by construction and is logged per
droplet class per reactor. `L_vap` and `T_sat` for LNG come from `CoolPropLNGProvider`
(existing `bubble_state_at_pressure` / `dew_state_at_pressure`), so mixture LNG is
handled with the same equation of state as the nozzle tiers.

### 3.5 Droplet heating — infinite liquid conductivity

```
c_l m_d dT_d/dt = A_d beta_spray Nu_d k_air (T_gas - T_d) / (2r)  +  (dm_d/dt) L_vap
Nu_d = (2 + 0.6 Re_d^0.5 Pr_d^(1/3)) * ln(1+B_d)/B_d
```

`beta_spray` is the paper's user-defined heat-transfer scaling factor; it is one of the
few tunable parameters and is treated as a declared calibration constant, defaulting to 1.

### 3.6 Energy accounting between droplet and gas (implementation-critical)

This is the detail that most CRN implementations get wrong, so it is stated explicitly:

- The gas **loses** `Q_conv` (convective heat to the droplet).
- The droplet **gains** `Q_conv`, spends `|dm/dt| L_vap` on vaporization, and puts the
  remainder into sensible heating.
- The vapor **enters the gas carrying enthalpy `h_vapor(T_d)`**, i.e. at the droplet
  temperature, not at the gas temperature.

Implemented this way, total enthalpy is conserved and the latent heat is never
double-counted. Task 7.2 is an explicit energy-closure test asserting this to
round-off on a non-reacting case.

### 3.7 Droplet transport between reactors

Droplets leaving reactor *i* are distributed to downstream reactors by spray
propagation fractions `eta_ij` (paper Section 2.2.5), which are user inputs in Phase 1
and optimization variables in Phase 6. Residence of a droplet class in reactor *i* is
`s_i / u_d`, where `s_i` is the spray path length through that reactor — a declared
input, not the reactor's bulk residence time, because droplets and gas do not travel
together. Unevaporated droplets convect on. Droplets that reach the exit are reported
as **liquid carryover**, which is a hard design failure signal.

### 3.8 Reaction

Gas-phase chemistry is integrated by Cantera with the user-supplied mechanism.
NOx comes from the mechanism's own nitrogen chemistry (extended Zeldovich, N2O
pathway, prompt/Fenimore, NNH) — no separate post-processed NOx model, and no
Zeldovich bolt-on. See Section 4.4 for the mechanism validation that this requires.

---

## 4. Cantera implementation strategy

Reference: <https://cantera.org/dev/reference/reactors/index.html>, read 2026-08-28 and
re-verified against the installed version in Task 1.1.

Confirmed from that page: Cantera exposes eight homogeneous reactor variants
(control-volume / constant-pressure, each in ideal-gas and mole-based forms), an
`ExtensibleReactor` family for custom governing equations, a dedicated steady
**Plug Flow Reactor**, `Reservoir` boundary conditions, the flow devices
`MassFlowController` / `PressureController` / `Valve`, `Wall` for heat transfer and
volume change, `ReactorSurface` for heterogeneous chemistry, and `ReactorNet`
integrating with CVODES/IDAS. `ReactorNet.advance_to_steady_state` exists.
**Sparse Jacobian preconditioning is supported only for `IdealGasMoleReactor` and
`IdealGasConstPressureMoleReactor`** — which is the reason this plan selects the
mole-based constant-pressure reactor as the default PSR type.

### 4.1 Reactor type mapping

| Network element | Cantera object | Notes |
|---|---|---|
| PSR (flame, recirculation, quench sub-zone) | `IdealGasConstPressureMoleReactor` | Mole-based reactors are preferred for large mechanisms because they support the adaptive preconditioner. Constant pressure matches assumption A3. |
| PFR (post-flame, dilution) | Series of `IdealGasConstPressureMoleReactor` (default), or Cantera's native steady `PlugFlowReactor` for a true spatial solve | Series-of-PSRs converges to plug flow, reuses the same network machinery, and participates in the simultaneous steady solve; grid convergence in reactor count is a verification task. The native PFR is a Phase 7 cross-check, not the default, because it does not embed in a recirculating `ReactorNet`. |
| Evaporator / mixer (specialized) | `IdealGasConstPressureMoleReactor` plus an evaporation source | See 4.3. |
| Air inlets, fuel-vapor source | `Reservoir` | Fuel-vapor reservoir temperature is set to the droplet temperature, per 3.6. |
| Exhaust | `Reservoir` | |
| Reactor-to-reactor flow | `MassFlowController` | Mass flows from the split matrix. |
| Network exit | `PressureController` | Holds combustor pressure and closes the mass balance. |
| Liner heat loss | `Wall` with prescribed `U`/area or prescribed `Q` | Per-reactor heat-loss fraction is a user input; adiabatic is the default with a warning that exit temperature will be biased high. |

Solution: `ReactorNet.advance_to_steady_state()`, with
`AdaptivePreconditioner` enabled for mechanisms above ~50 species, and a fallback of
long-time `advance()` with a steady-state residual check if the steady solver stalls.

### 4.2 Network solve

All reactors and controllers live in **one** `ReactorNet`. Recirculation zones — which
are essential to flame stabilization and are the whole reason a CRN is not a chain —
are handled natively by the simultaneous solve, so no recycle-tearing iteration is
needed. A sequential-modular fallback with successive substitution and Aitken
acceleration is planned only if the simultaneous solve proves fragile (risk R3).

**Mass balance closure.** User-supplied splits will not close exactly. Before
constructing the network, the split matrix is corrected by the paper's minimum-norm
least-squares projection (their Eq. 23):

```
minimize ||dz||_2   subject to   A (z + dz) = b
```

solved as an equality-constrained least-squares problem. The correction magnitude per
reactor is reported; a correction above a threshold raises a WARNING because it means
the user's splits were substantially inconsistent.

### 4.3 Droplet-gas coupling: operator splitting (primary)

Embedding stiff droplet ODEs inside a stiff chemistry integration is the fastest way to
a fragile solver. The primary approach is **operator-split with outer iteration**:

1. Guess each reactor's gas state (initial guess from an unsteady ignition solve or from
   equilibrium).
2. Integrate all droplet classes through the network with the gas states frozen, using
   `scipy.integrate.solve_ivp` (LSODA/BDF) on the coupled `(r, T_d)` system.
   Result: per-reactor fuel-vapor source `mdot_evap,j`, source temperature, and `Q_conv,j`.
3. Rebuild the Cantera network with those sources and solve to steady state.
4. Repeat 2–3 until reactor temperatures and evaporation rates converge (relative
   tolerance ~1e-4, damped/under-relaxed updates).

This converges in a handful of outer iterations in practice, keeps each sub-problem
well-conditioned, and makes the droplet physics independently testable.

**Fully-coupled alternative (Phase 5, optional).** Cantera's `ExtensibleReactor`
family (confirmed present) allows Python-defined additional state variables and source
terms inside the reactor's own governing equations. If the
split iteration proves too slow or inaccurate for strongly coupled evaporator zones,
droplet state can be added as extra reactor state variables there. This is scoped as an
upgrade, not the first implementation.

### 4.4 Mechanism handling — and a warning worth reading

Fuel selection is per operating point (decision #1). A `MechanismSpec` carries the YAML
path, the fuel species and surrogate mole fractions, the phase name, a provenance
string, and a declared validity range in pressure, temperature, and equivalence ratio.

**Jet-A — resolved, with a caveat that changes the design.** Two Stanford HyChem
files were inspected on 2026-08-28:

| File | Species | N chemistry | Low-T / NTC | Use here |
|---|---|---|---|---|
| `A2NOx_skeletal.cti` | 71 | **Yes** — HyChem A2 v2.0 + Glarborg et al. NOx (*Prog. Energy Combust. Sci.* 67 (2018) 31–68); NO, NO2, N, N2O, HCN, HNO, CN, NH, NH2 all present | **No** — high-temperature chemistry only | **Reactor network / NOx** |
| `A2NTCfast_ske.cti` | 47 | **No** — only N2 as bath gas | **Yes** — 8 pyrolysis + 8 low-temperature oxidation steps, fast-NTC variant | **Autoignition delay tabulation** |

Both are skeletal reductions by T. Lu (UConn, 2018) of HyChem Jet-A POSF-10325
(`POSF10325`, C11H22), and both use the `CH2*` naming convention rather than `CH2(S)`.

**Neither file has both NOx and low-temperature chemistry, and this tool needs both —
for different sub-models.**

**Measured 2026-08-28 during Phase 0, and the result is stronger than anticipated.**
Ignition delay at premixer-relevant conditions (`phi` = 0.5), high-temperature-only
versus low-temperature mechanism:

| P | 700 K | 800 K | 900 K | 1000 K |
|---|---|---|---|---|
| 20 atm | **510x** | **71x** | 9.3x | 2.0x |
| 40 atm | **490x** | **86x** | 13.5x | 2.9x |

LTO `T3` is roughly 700–900 K. Using the high-temperature-only mechanism there would
overpredict ignition delay by **one to two orders of magnitude — declaring a premixing
passage safe when it is not.** That is a safety-relevant error in the worst possible
direction, and it settles the two-mechanism question definitively.

**A correction to an earlier conflation in this plan.** The non-monotonic NTC *shape*
turns out to be strongly pressure-dependent: clear at 5 atm (800–950 K), weak at 10 atm,
and absent by 20–40 atm. It would be easy to read that as "low-temperature chemistry does
not matter at combustor pressure." **It is not.** The 510x table above is measured at 20
and 40 atm, where the non-monotonic shape has already vanished, and low-temperature
chemistry still dominates the absolute delay. The shape disappearing and the chemistry
mattering are two different statements. Full evidence in `mech/README.md`.

The resolution is a **declared two-mechanism split**, each used strictly inside its own
validity:

- `A2NOx_skeletal` drives the reactor network, combustion, and all NOx output;
- `A2NTCfast_ske` drives the ignition-delay table consumed by `autoignition.py`.

`MechanismRegistry` therefore carries a `role` on each `MechanismSpec`
(`NETWORK` or `IGNITION_DELAY`), and every result records which mechanism produced which
number. A single mechanism serving both roles is permitted and is the simpler path if a
combined NOx+NTC A2 model becomes available. **Consistency between the two is not
assumed and must be checked** (Task 7.8): high-temperature ignition delays from both
files should agree where both are valid, and a disagreement above a declared tolerance
is reported rather than averaged away.

**Format.** Both files are **CTI, which Cantera 3.x no longer parses.** They must be
converted once with the `cti2yaml` script that ships with Cantera, and the converted
YAML, its SHA-256, the source URL, and the converter version recorded (Task 0.5).

**LNG / natural gas — confirmed.** **GRI-Mech 3.0**
(<http://combustion.berkeley.edu/gri-mech/version30/text30.html>), 53 species,
325 reactions, with nitrogen chemistry included (thermal, prompt/Fenimore, N2O, and
NNH pathways, plus reburn). It ships with Cantera as `gri30.yaml`, so no conversion is
needed and the version is pinned by the Cantera version.

Two validity caveats are recorded on the `MechanismSpec` and warned at runtime, not
buried:

1. **Pressure.** GRI-Mech 3.0 was optimized against data to roughly 10 atm. A combustor
   at 10-40 atm is an extrapolation, and NO predictions are the most pressure-sensitive
   output. The tool warns whenever `P3` exceeds the declared range. If dual-fuel NOx
   comparisons later hinge on absolute LNG values, FFCM-2 or an AramcoMech variant with
   an N submodel is the upgrade path; the registry accepts either without code changes.
2. **Composition.** GRI-Mech 3.0 covers methane, ethane, and propane chemistry. Higher
   alkanes or significant nitrogen diluent in the `LNGComposition` must be checked
   against the mechanism species list, and an unrepresented component **raises** rather
   than being silently dropped.

Because the Jet-A and LNG mechanisms differ in size, provenance, and validity, **absolute
NOx is not comparable across fuels without stating this.** Every dual-fuel comparison
reports both mechanisms and flags that the comparison is trend-level.

**General mechanism guard, still required.** For any supplied mechanism the tool will:

- verify `NO`, `NO2`, `N`, `N2O` are present when NOx is requested;
- **record which pathways each mechanism actually carries** rather than requiring a fixed
  species list. Measured 2026-08-28: `A2NOx_skeletal` has the modern `NCN` prompt route
  but **no `NNH`**; `gri30` has `NNH` but uses the older `HCN`/Fenimore prompt route.
  Neither is a superset of the other, so a hard `NNH` requirement would wrongly reject a
  good mechanism. Pathway coverage is attached to every result instead;
- **raise an error, not a warning,** if NOx is requested and no nitrogen chemistry is
  present. Silent zeros are worse than a failed run.

Merging mechanisms remains **the user's call**; the tool validates and reports, never
merges. Task 4.3 provides a merge *checker* (duplicate reactions, inconsistent thermo,
species-name collisions such as `CH2*` vs `CH2(S)`), not a merger.

**LNG / natural gas.** GRI-Mech 3.0 (53 species) includes NOx and is the obvious
default, but it is validated only to roughly 10 atm and to methane-dominant mixtures.
At combustor pressures of 10–40 atm the tool will emit a validity WARNING. FFCM-2 or an
AramcoMech variant with an N submodel are better choices at pressure; the mechanism
registry accepts any of them. Higher LNG alkanes (ethane, propane) must be present in
the mechanism if they are present in the `LNGComposition` — this is checked, and an
unrepresented component raises an error rather than being dropped.

### 4.5 Performance

Estimated per operating point: seconds to ~1 minute for a 7–12 reactor network with a
120-species mechanism, consistent with the paper's 12–36 s for 5–11 clusters. Mechanism
objects and `Solution` instances are cached and reused across operating points, the same
way `CoolPropLNGProvider` is reused across a study today.

---

## 5. Architecture templates (Phase 1 topologies)

Each template is a parametric builder returning a `CombustorNetwork`. All splits,
volumes, and path lengths are user inputs. Each template works for either fuel; the
comparison you actually want is **the same template evaluated with Jet-A at LTO
conditions and LNG at cruise conditions**, and then across templates.

### 5.1 RQL

```
[dome air f_dome] + [fuel spray]
   -> EVAPORATOR (rich, phi ~ 1.2-1.8)
   -> MIXER (rich)  <-> RECIRC-RICH (PSR, feeds back)
   -> QUENCH: N sub-PSRs with distributed jet air (f_quench spread over N stages)
   -> LEAN PSR
   -> PFR (post-flame + dilution f_dilution)
   -> exit
```

The quench is the design crux and is *not* modeled as a single mixing point. It is a
chain of `N_quench` sub-reactors with a user-specified air-addition schedule and a
quench time `tau_q`. NOx produced during quench — as the mixture traverses stoichiometric
on its way from rich to lean — is the dominant RQL NOx source and the whole reason RQL
succeeds or fails. A single perfectly-mixed quench point would hide it. Sensitivity to
`N_quench` and `tau_q` is a required output, not an option.

### 5.2 Lean premixed / prevaporized (LPP)

```
[premixer air f_dome] + [fuel spray]
   -> PREMIXER (evaporator + mixer, CHEMISTRY FROZEN, autoignition checked)
   -> FLAME PSR  <-> RECIRC PSR
   -> PFR (+ dilution)
   -> exit
```

Chemistry is frozen in the premixer by construction (that is what a premixer is *for*),
but the autoignition module then checks whether that assumption is defensible. If
`tau_ign < tau_res` in the premixer, the tool emits an **ERROR-severity** warning: the
LPP design is invalid at that operating point, and the reported emissions are
meaningless. This is expected to be the binding constraint for Jet-A at LTO `T3`.

### 5.3 Lean direct injection (LDI)

```
[dome air f_dome] + [fuel spray]
   -> EVAPORATOR (lean, no premixing passage)
   -> MIXER  <-> RECIRC PSR
   -> FLAME PSR
   -> PFR (+ dilution)
   -> exit
```

No premixing residence time, therefore no autoignition path — at the cost of a less
uniform equivalence ratio and higher NOx than ideal LPP. The LDI-vs-LPP comparison is
precisely the "how much premixing dare I use" question, and this tool answers it
quantitatively.

### 5.4 Staged pilot / main

A pilot circuit (richer, always lit, stability) and a main circuit (lean, high power),
with a user-specified fuel split. Both circuits burn the *same* fuel at any given
operating point, per decision #1. Needed to represent realistic lean-burn hardware,
which is never a single unstaged dome.

### 5.5 Custom

An explicit reactor list plus a connection matrix, for topologies not covered above and,
later, for CFD-derived networks.

---

## 6. Autoignition and flashback screening

Serves the stated goal: *exploit flash atomization, avoid autoignition*.

For each premixing or mixing passage:

1. **Mixture state at the passage.** Air at `T3`, `P3`, mixed with fuel vapor. For LNG,
   the vapor arrives cold — flash vaporization is an internal refrigerator, and the
   resulting mixture temperature `T_mix` can be tens of kelvin below `T3`. Computed by
   adiabatic mixing with the actual vapor enthalpy from CoolProp, including the latent
   heat drawn from the air. **This is a design lever, and quantifying it is a headline
   output of this tool.**
2. **Ignition delay** `tau_ign(T_mix, P3, phi_local)` from a Cantera constant-pressure
   homogeneous reactor, with the ignition criterion defined as max `dT/dt` (and OH
   peak recorded as a cross-check). Tabulated over a `(T, P, phi)` grid once per
   mechanism and interpolated, because it is called inside sweeps.
3. **Residence time** `tau_res` in the passage, from the passage volume and the flow,
   with the injector exit velocity taken from the existing nozzle result.
4. **Margin** `M = tau_ign / tau_res`. Design guidance and the literature convention put
   the acceptable floor around `M >= 3–5`; the threshold is a declared, editable input,
   not a hidden constant.
5. **Flashback** is screened separately by comparing the turbulent flame speed
   (correlation-based, declared) with the local passage velocity, and by a boundary-layer
   flashback criterion. Reported as a screen with explicit uncertainty; not calibrated.
6. **Evaporation completeness** at the premixer exit is reported alongside, because an
   LPP design that avoids autoignition by running cold but exits with 40% liquid is not
   an LPP design.

Expected and testable finding: LNG tolerates far more premixing than Jet-A at the same
`T3`, and the flash-cooling effect widens that margin further. That would be a concrete
argument for lean-premixed LNG cruise plus RQL or LDI Jet-A LTO on shared hardware —
exactly the question in the goal statement. The tool exists to test it, not to assume it.

---

## 7. Cryogenic thermal fuel management

Closes the loop that already half-exists in `feed.py`.

- **Heat sink capacity.** `Q_sink = mdot_LNG * (h(T_nozzle_inlet) - h(T_tank))` from
  CoolProp, decomposed into sensible-liquid, latent, and superheat contributions, and
  reported against the heat loads the user wants to dump (oil cooler, ECS, generator,
  Jet-A anti-coking margin).
- **Required heat for target flash quality.** Inverse problem: given a target exit
  vapor quality or a target flash-spray regime from `spray.py`, find the nozzle-inlet
  enthalpy that achieves it, then the feed-line heat load that delivers it. Solved by a
  bracketed root on the existing Tier 0–3 chain. This turns "how much should I heat the
  LNG" from a guess into a number.
- **Constraint set.** The answer must simultaneously satisfy: enough superheat for the
  desired atomization regime; enough subcooling margin *upstream* in the feed line to
  avoid vapor lock and two-phase pump feed (already computed by `screen_lng_flash` and
  `solve_lng_feed_line`); autoignition margin (Section 6); and the available heat loads.
- **Output.** A feasible-window plot and table in `(T_fuel_nozzle, dP_nozzle)` space,
  with the four constraint boundaries drawn on it, per operating point. This is the
  single most useful design artifact the tool will produce.

---

## 8. Dual-fuel design optimization — the primary deliverable

This section is the point of the tool. Sections 3–7 build an *evaluator*; this section
turns it into a *design* method.

**The three concrete outputs this tool must produce are stated in Sections 8.10 and 8.11.**
Read those first; everything else in Section 8 exists to support them.

### 8.1 The design problem, restated for separate fuel paths

**Per the 2026-08-28 direction: the Jet-A and LNG circuits are separate fuel paths with
separate injector hardware.** They atomize by entirely different mechanisms — a
pressure-swirl sheet versus flash-driven breakup — so forcing them through common
hardware would compromise both. Each injector is designed for its own fuel, with no
cross-fuel compromise.

This removes the injector-level conflict completely. **What remains shared is the
combustor**, and that is where the whole design tension now lives:

| Genuinely per-fuel — **no compromise** | Genuinely shared — **compromise required** |
|---|---|
| LNG orifice count, diameter, `L/D` | Liner air split: dome, primary, quench, dilution, cooling |
| Jet-A swirler ports, tangency radius, swirl chamber, exit orifice | Zone volumes and residence times (one liner length) |
| **Premixing passage length and area — per path** | Combustor pressure loss (set by total hole area) |
| Per-path swirler/air passage sizing | Dome packaging envelope (both injectors must fit) |
| Fuel temperature schedule | Quench port axial station and jet momentum |
| Nozzle `dP` and pump pressure | Thermal environment shared by both circuits |

**Consequence worth stating plainly: the premixing-length conflict I previously called
the crux is gone.** With separate paths, LNG gets a long premixing passage and Jet-A gets
a short one or none, and the autoignition constraint decouples per path. That was the
right thing to correct.

**The crux is now the air split.** A liner plumbed for RQL Jet-A (rich dome, `phi ~ 1.4-1.6`,
large quench jets) and one plumbed for lean-premixed LNG (`phi ~ 0.5` dome) differ in dome
air fraction by roughly a factor of three. Fixed hole areas cannot deliver both — unless
one of the mechanisms in Section 8.3 is used.

### 8.2 Two effects that separate paths introduce, and that must be modeled

These follow directly from having two circuits, and neither exists in a single-fuel
combustor. Both are new model requirements, not commentary.

#### 8.2.1 The inactive circuit still passes air — and this is a free design lever

Air flows through **both** injector air passages at all times, including the one whose
fuel is shut off. When LNG is burning, the Jet-A swirler is still flowing unfueled air
into the dome, and vice versa.

Let the Jet-A passage carry dome-air fraction `a_J` and the LNG passage `a_L`. Then:

- **Burning Jet-A:** near-field equivalence ratio is set by `a_J` alone; `a_L` enters as
  unfueled air that dilutes downstream.
- **Burning LNG:** near-field equivalence ratio is set by `a_L` alone; `a_J` dilutes.

So **sizing the two passages sets how the effective near-field air split shifts between
segments, with completely fixed hardware.** Making `a_J` small and `a_L` large yields a
rich Jet-A near field (RQL-like) and a lean LNG near field (LPP-like) *in the same
liner*. This is passive air staging obtained for free from the dual-path architecture.

**This is the single most promising idea in the plan and it is stated as a hypothesis for
the tool to test, not as a result.** The honest caveat: whether unfueled air from the
idle passage stays segregated from the near field long enough to preserve the intended
local `phi`, or short-circuits into it, is a mixing question. A CRN answers it only to the
extent that the user-supplied split and mixing rates are right, and it is exactly the
question the Phase 6 CFD should later settle. The tool will therefore report near-field
`phi` **with an explicit sensitivity to the idle-passage mixing fraction**, and will not
report a design conclusion that survives only at one assumed value.

#### 8.2.2 The idle circuit is a thermal and safety problem

A dual-fuel nozzle spends every mission segment with one circuit hot and stagnant:

- **Jet-A idle during LNG cruise:** stagnant fuel in a hot nozzle **cokes**. Wall
  temperature must stay below a declared coking limit, which is a real constraint on how
  the Jet-A circuit is packaged next to a burning dome.
- **LNG idle during Jet-A LTO:** stagnant cryogenic fuel warms, boils, and **vapor-locks**;
  restart then delivers vapour rather than liquid.
- **Cross-talk:** the cryogenic LNG line runs next to a hot Jet-A line inside one nozzle
  body. Heat leak into LNG is a loss to the Jet-A circuit and a gain to the LNG circuit,
  and it operates in both directions depending on which fuel is running.

This becomes an explicit screen (Task 5.7), reusing the conduction/insulation machinery
already in `feed.py`. It is a **safety-margin** output in the sense the direction asked
for, and it can veto an otherwise attractive packaging.

### 8.3 Reconciling the shared air split — four candidate mechanisms

Since the air split is now the binding compromise, the tool must be able to evaluate every
credible way around it. These are the design options, and comparing them is a headline
output:

1. **Matched-regime operation.** Run both fuels in the same dome-`phi` regime — Jet-A LDI
   plus LNG LPP, both lean. Air split is then compatible by construction. Costs Jet-A the
   RQL soot/stability advantages; buys simplicity.
2. **Passive air staging via passage sizing** (Section 8.2.1). No moving parts. The most
   attractive option if the mixing assumption holds.
3. **Fuel staging.** Shift effective near-field `phi` by re-allocating fuel between pilot
   and main circuits per segment, with fixed air. This is the practical lever real staged
   combustors use, and it is a Class B (schedulable) variable.
4. **Variable geometry / staged air valves.** The only mechanism that fully decouples the
   two, at real cost in complexity, weight, and certification. Evaluated so the benefit of
   the complexity can be quantified, not because it is recommended.

### 8.4 What cryogenic LNG buys, and how each effect is modeled

The direction asked to exploit *every* aspect of cryogenic LNG for mixing performance and
safety margin. Each item below is an explicitly modeled effect with a named home in the
code, so none of them stays a slogan:

| Effect | Mechanism | Where modeled |
|---|---|---|
| **Flash atomization** | Superheated discharge shatters the jet without needing high `dP` | Existing Tiers 0–3, `spray.py` |
| **Flash cooling of the premix** | Vaporization draws latent heat from the air; `T_mix` drops tens of K below `T3` | `autoignition.py` (Section 6) |
| **Longer ignition delay from that cooling** | Lower `T_mix` raises `tau_ign` exponentially | `autoignition.py` — **this is the safety margin that buys the long premixer** |
| **Therefore: a genuinely long LNG premixer** | More premixing length within the same margin | Per-path premix length, Class A2 variable |
| **Better mixing, lower NOx** | Uniform lean `phi` removes the hot streaks that drive Zeldovich NO | CRN unmixedness metric, objective 4 |
| **Denser cold fuel** | Smaller lines and orifices for the same mass flow | Existing `CoolPropLNGProvider` |
| **Heat sink capacity** | Cooling duty for oil, ECS, and the Jet-A circuit | `thermal.py` (Section 7) |
| **Methane's intrinsic chemistry** | Long ignition delay, high H/C, low soot | GRI-Mech 3.0 in the network |

The chain that matters runs left to right: **flash cooling → longer ignition delay →
longer premixer allowed → better mixing → lower NOx.** Quantifying that chain end to end
is the strongest single argument this tool can make for LNG, and every link in it is
computed rather than assumed.

### 8.5 Design variables

**Class A — shared combustor.** Chosen once; both fuels live with it.

- Liner air split: dome, primary, quench, dilution, cooling fractions
- Quench port axial station and jet momentum ratio
- Zone volumes: dome, primary-zone length, quench length, total length
- Combustor pressure loss
- Dome packaging envelope

**Class A2 — per-fuel injector hardware.** Independently optimized; **no cross-fuel
compromise.**

- LNG: orifice count, diameter, `L/D`, premixing passage length and area, air passage
  fraction `a_L`
- Jet-A: inlet port count and diameter, tangency radius, swirl chamber radius and length,
  exit orifice diameter, premix length (may be zero), air passage fraction `a_J`
- Constraint linking them: `a_J + a_L = ` dome air fraction (Class A), and both injectors
  must fit the dome envelope

**Class B — schedulable per segment.**

- Fuel temperature at the nozzle (LNG: set by feed-line heat pickup — the main lever)
- Nozzle pressure drop and pump pressure
- Pilot/main fuel split
- Fuel mass flow (set by the mission, not free)

### 8.6 Objectives and constraints

**Objectives** (multi-objective; a Pareto front, not a single scalar):

1. Minimize NOx at the Jet-A LTO points — ICAO `Dp/Foo` if thrust and time-in-mode are
   supplied, otherwise EI weighted by time-in-mode
2. Minimize NOx at the LNG cruise points
3. Maximize atomization quality for both fuels — SMD and evaporation completeness at the
   primary-zone entry
4. Maximize mixing uniformity — an equivalence-ratio unmixedness metric at the flame zone
5. **Maximize exit temperature uniformity** — minimize the temperature spread among the
   reactors feeding the exit (see the caveat in Section 8.10.3)

Objectives 3 and 4 are the confirmed QoIs from decision #2 and are first-class, not
constraints, because atomization and mixing quality **are** the nozzle design question.
Objective 5 was added on 2026-08-28 as a named deliverable.

**Hard constraints** (violation makes a design infeasible, not merely worse):

- Autoignition margin `M = tau_ign/tau_res >= M_min`, evaluated **per path** — Jet-A in
  its own premixer (or trivially satisfied if it has none) and LNG in its own
- LNG thermal window: superheat sufficient for the target flash regime at the nozzle
  **and** subcooling preserved upstream (no vapour lock, no two-phase pump feed) — the
  four-constraint window of Section 7
- **Idle-circuit limits (Section 8.2.2): Jet-A wall temperature below the coking limit
  while LNG burns; LNG circuit restartable after Jet-A operation**
- Zero liquid carryover at the combustor exit, both fuels
- Pressure budget closes at every operating point (existing `fuel_pressure_budget`)
- No Jet-A cavitation, no LNG upstream two-phase
- Turndown: the same hardware flows takeoff and idle within available pump pressure
- Lean blowout margin at LNG cruise — a screen, with stated uncertainty
- Dome packaging: both injectors and their air passages fit

Constraints use **feasibility-first ranking** (feasible dominates infeasible; among
infeasible, rank by total violation) rather than weighted penalties, so a design is never
rewarded for trading a safety constraint against an emissions objective.

### 8.7 Optimization strategy — staged, deliberately not "press the optimize button"

One evaluation = build geometry, run the existing nozzle solvers for both circuits, build
spray boundaries, run the CRN at every mission point, aggregate. At seconds to a minute
per point and 4–8 points, that is **minutes per design**.

**Stage 1 — Sweeps and sensitivity.** Sobol/Morris screening over Classes A, A2, and B.
Produces the ranking that says which five variables actually matter. Expected to be the
most valuable single output, and it reduces dimension for everything after.

**Stage 2 — Multi-objective Pareto search.** NSGA-II via `pymoo` over the reduced set,
feasibility-first constraints, evaluations parallel across cores, checkpointed.

**Stage 3 — Surrogate-assisted refinement** near the Pareto knee. Optional.

**Stage 4 — Architecture-pair comparison.** Architecture is now a **pair** — one per fuel
path, sharing a liner. Credible pairs: (Jet-A RQL, LNG LPP), (Jet-A LDI, LNG LPP),
(Jet-A RQL, LNG LDI), (Jet-A LDI, LNG LDI). Each pair is evaluated against each air-split
mechanism from Section 8.3.

Staged on purpose: a single global optimization over a partly-validated reduced-order
model yields a confident number nobody should believe. Sweep first, optimize second.

### 8.8 Headline output: the cost of sharing a combustor

With separate injectors, the compromise is no longer in the nozzle — it is **entirely in
the shared liner**. So the headline metric sharpens accordingly.

Run three optimizations: **Jet-A-optimal liner**, **LNG-optimal liner**, and the
**shared-liner compromise**, each with its own fully-optimized per-fuel injectors. Then:

```
                         NOx (Jet-A LTO)     NOx (LNG cruise)
Jet-A-optimal liner            best              (poor)
LNG-optimal liner             (poor)              best
Shared-liner compromise      best + dA          best + dL
```

`dA` and `dL` are **the quantified price of sharing one combustor between two fuels** —
the number nobody has for this configuration. Reported per architecture pair and per
air-split mechanism, so the deliverable reads: *"a shared liner costs dA on Jet-A NOx and
dL on LNG NOx; passive passage sizing recovers most of it; variable geometry recovers the
rest but is not worth the complexity"* — or whatever the model actually says.

Because `dA` and `dL` are each computed **within one fuel and one mechanism**, this metric
is immune to the cross-fuel mechanism-comparability problem of R1e.

Supporting artifacts: the Section 7 four-constraint thermal window per segment; an
autoignition margin map in `(premix length, T_fuel)` drawn **separately per path**, since
they no longer share a passage; the idle-circuit thermal screen; and the Stage 1
sensitivity ranking.

### 8.9 Honesty constraints on the output

Non-negotiable, and stated in the report the tool emits:

- The Pareto front is a front **of this model**, not of the hardware. Its value is
  *relative* ranking and *trend* identification.
- Cross-fuel absolute NOx is trend-level only (R1e); `dA` and `dL` are within-fuel by
  construction and are the defensible metric.
- The passive-air-staging result (Section 8.2.1) depends on an assumed idle-passage mixing
  fraction and **must be reported with that sensitivity**, never at a single value.
- Any conclusion that does not survive the Stage 1 sensitivity range, or that flips under
  a mechanism substitution, is **not reported as design guidance**.
- Optimal designs are candidates for CFD and rig testing, never final designs.

### 8.10 Primary design answers 1 and 2

The 2026-08-28 direction named the two concrete numbers this tool exists to produce.
A third — target nozzle geometry per fuel path — follows in Section 8.11.
Both are questions about the **head-end (dome) air fraction**, asked from opposite ends —
which is the clearest possible confirmation that the shared liner air split is the crux
(Section 8.1).

#### 8.10.1 Answer 1 — How much air to the head end during Jet-A LTO?

**Output:** the dome/head-end air fraction `f_dome` (and the resulting rich-zone
equivalence ratio `phi_rich`) that minimizes Jet-A NOx over the LTO cycle.

**Method.** Sweep `f_dome` across the RQL range, and at each value run the full LTO
mission set (takeoff, climb-out, approach, idle) with the Jet-A path active and the LNG
path passing unfueled air (Section 8.2.1). Report NOx as ICAO `Dp/Foo` when thrust and
time-in-mode are supplied, and as EI otherwise.

**Expected shape of the answer.** A minimum in NOx at intermediate `phi_rich`. Too rich
and the quench traverse through stoichiometric is long and hot; too lean and the dome
itself approaches stoichiometric. Classical RQL practice puts the optimum near
`phi_rich ~ 1.4-1.6`, and **if the tool does not reproduce a minimum in that neighbourhood,
that is a signal to distrust the tool, not the literature.** This is a built-in sanity
anchor, already listed in Section 12.3.

**Constraints that bound the answer:** quench-zone residence time (a slow quench destroys
the whole benefit — hence the staged quench of Section 5.1); combustor exit temperature;
Jet-A autoignition margin in its own path; idle-circuit limits; and turndown to idle.

**Reported with:** sensitivity to quench stage count `N_quench` and quench time `tau_q`,
because those set how much NOx is made in transit through stoichiometric and are the
largest modelling uncertainty in the RQL answer.

#### 8.10.2 Answer 2 — How lean can LNG run at cruise?

**Output:** the minimum stable cruise equivalence ratio `phi_lean` (equivalently the
maximum dome air fraction for the LNG path), the NOx achieved there, and the exit
temperature uniformity that comes with it.

**Method.** Sweep `phi_lean` downward at cruise conditions with the LNG path active and
the Jet-A path passing unfueled air, tracking NOx, exit temperature spread, combustion
efficiency, and approach to extinction.

**How the lean limit is determined — and its honest uncertainty.** The limit is bracketed
by three independent criteria rather than asserted from one:

1. **PSR extinction by continuation.** Reduce the flame-zone residence time (or lower
   `phi`) until the Cantera PSR extinguishes. This is a genuine physics-based limit
   computed from the mechanism, not a correlation, and it is the most trustworthy of the
   three. **Caveat:** it is the extinction limit of an *idealized perfectly stirred
   reactor*, not a real combustor blowout, which additionally involves flow dynamics,
   recirculation strength, and unsteadiness that a steady CRN cannot represent.
2. **CO and combustion-efficiency rise.** The practical lean limit is usually reached
   when CO climbs before extinction does. **But CO is an uncalibrated diagnostic in this
   tool** (decision #2), so this criterion is directional only.
3. **Exit temperature uniformity degradation.** As the flame weakens, the temperature
   spread among exit-feeding reactors grows.

**The lean limit is therefore reported as a bracket with the governing criterion named,
never as a single number**, and it is flagged as the output most in need of rig
validation. Promising a precise LNG lean limit from this tool would be overclaiming.

**Why LNG can go leaner than Jet-A, and why that is the point.** The chain of Section 8.4
runs: flash cooling lowers `T_mix`, which raises `tau_ign`, which permits a genuinely long
premixing passage, which delivers a uniform lean mixture, which suppresses the hot streaks
that drive Zeldovich NO. Every link is computed. **Quantifying how much leaner LNG can run
than Jet-A, and how much NOx that buys, is the single most valuable output of this tool.**

#### 8.10.3 The two answers are coupled — and that coupling is the deliverable

Answer 1 wants a **low** dome air fraction for Jet-A (rich head end). Answer 2 wants a
**high** dome air fraction for LNG (lean head end). One liner, fixed holes.

The reconciliation mechanisms of Section 8.3 are what make both achievable, and the
passive lever of Section 8.2.1 is the most promising: by sizing the two injector air
passages `a_J` and `a_L`, the *effective* head-end air differs between segments even
though total hole area is fixed. **So the real output is not two independent numbers but a
triple:** `(f_dome, a_J, a_L)` that simultaneously delivers rich Jet-A LTO and lean LNG
cruise — plus the residual `dA`/`dL` penalty of Section 8.8 that says how much the sharing
still costs.

**Caveat on "uniform temperature profile" — read before relying on it.** A chemical
reactor network **cannot predict pattern factor**. Pattern factor is a two-dimensional
exit-plane temperature map governed by dilution-jet penetration and mixing that a
zero-dimensional network does not resolve. What this tool reports is the **temperature
spread among the reactors feeding the exit**, which is a *proxy* for exit uniformity and
is useful for ranking designs against each other. It is not a pattern-factor number, must
not be quoted as one, and cannot substitute for CFD or a rig traverse. Ranking designs by
it is defensible; predicting turbine inlet distortion from it is not.

### 8.11 Answer 3 — target nozzle geometry for each fuel path

**Good news first: this is the easiest of the three answers.** Because the paths are
separate hardware (Section 8.1), each injector is a *single-fuel* design problem with no
cross-fuel compromise. The only couplings are the dome packaging envelope and the air
passage fractions `a_J` / `a_L`, which tie back to the liner. Everything else is
independent.

Geometry is reported as a **feasible region with the driving constraint named for each
feature**, not as a single point — the same philosophy as the thermal window of Section 7.
A bare number like "0.42 mm" hides which assumption it rests on.

#### 8.11.1 Jet-A pressure-swirl path

Driven by `solve_jet_a_pressure_swirl`, which already computes most of this. Task 8.7f
adds the standard atomizer design groups that the existing result does not yet report.

| Feature | Symbol | Primary driver | What bounds it |
|---|---|---|---|
| Exit orifice diameter | `d_o` | Max-power `mdot`, `Cd`, `dP` | Turndown at idle; manufacturability; clogging |
| **Flow number** | `FN = mdot/sqrt(rho dP)` | Mission max flow vs available `dP` | Pump pressure budget (existing `fuel_pressure_budget`) |
| **Atomizer constant** | `K = A_p/(D_s d_o)` | Target cone angle and `Cd` | Film thickness; air-core stability |
| Swirl chamber ratio | `D_s/d_o` | Cone angle, discharge coefficient | Dome packaging |
| Exit orifice length ratio | `l_o/d_o` | Film development, `Cd` | Too long collapses the cone |
| Inlet port count | `n_p` | Flow uniformity, manufacturability | Minimum drillable port diameter |
| Inlet port diameter | `d_p` | Inlet velocity, hence swirl momentum | Port cavitation (existing screen) |
| *Resulting:* cone angle | `2 theta` | `u_t/u_a` at exit | Dome impingement, wall wetting |
| *Resulting:* film thickness | `t` | Air-core solution | Sets SMD — **gated on calibration** |

`FN` and `K` are the two groups that atomizer design charts are actually written in, and
neither is currently an output. Adding them is the difference between a result an
injector designer can use and one they cannot.

**Turndown is the geometric driver most likely to force the architecture.** A simplex
atomizer's flow scales as `sqrt(dP)`, so a wide LTO thrust range may not be coverable by
one orifice within the available pump pressure — which is precisely what forces duplex or
pilot/main staging in real hardware. The tool evaluates atomization quality (`Re`, `We`,
`Oh`, film thickness) at **idle as well as takeoff**, and reports if a single simplex
cannot cover the range. That is a genuine architecture finding, not a detail.

#### 8.11.2 LNG flashing path

Driven by Tiers 1–3. The critical difference from Jet-A: **`L/D` is the primary design
knob, and the existing Tier 2 machinery is uniquely able to inform it.**

| Feature | Driver | What bounds it |
|---|---|---|
| Orifice diameter | `mdot/(Cd G)` with `G` from the Tier 1 SPI/HEM bracket and Tier 2 | Choking; manufacturability |
| Orifice count | Total area split; spray distribution and mixing | Packaging; jet-to-jet interaction |
| **`L/D`** | **Flash onset location** — Tier 2 residence time versus relaxation time `tau` | Too long: flash migrates upstream, risking vapour in the feed. Too short: flash happens externally, wasting the atomization |
| Premix passage length | Autoignition margin (Section 6) | `tau_ign/tau_res >= M_min` |
| Premix passage area | Velocity, hence `tau_res` and flashback margin | Flashback screen |
| Air passage fraction `a_L` | Near-field `phi` at cruise (Section 8.2.1) | `a_J + a_L =` dome fraction |
| *Resulting:* exit quality `x` | Tier 2 | Target flash regime |
| *Resulting:* regime | Tier 3 classifier | Fully-flashing target |

**The `L/D` answer is the most valuable geometric output for LNG and also the most
fragile.** Tier 2 predicts where flash onset falls along the hole as a function of `L/D`,
which is exactly the design question. But it does so through the relaxation time `tau`,
which the existing documentation is explicit about: **`tau` is a calibration parameter,
not a universal LNG constant.** So the tool delivers the *trend* and the *sensitivity of
onset location to `L/D`* reliably, while the absolute `L/D` target moves with `tau`.
**One flashing flow test on representative hardware would pin `tau` and convert this from
a trend into a number.** That is the highest-value single experiment this plan can
recommend, and it is cheap relative to a rig campaign.

#### 8.11.3 What this tool cannot tell you about geometry

Stated plainly so these are sourced from CFD or rig rather than over-read from here:

- **Inlet edge radius, chamfer, and surface roughness.** These strongly affect discharge
  coefficient and, for LNG, nucleation onset. Not modeled at all.
- **Manufacturing tolerance stack-up** and its effect on cup-to-cup flow scatter.
- **Spray-swirler air interaction** — the aerodynamic field that actually disperses the
  spray is not resolved; the CRN receives mixing rates as inputs.
- **Film breakup length and ligament structure** for the Jet-A sheet.
- **Jet-to-jet interaction** between multiple LNG orifices.
- **Any SMD at all without a hardware calibration.** The existing package refuses to
  invent one, and every geometric target that depends on droplet size inherits that
  gating. Geometry driven by *hydraulics* (diameter, `FN`, `L/D`, cone angle, film
  thickness) is available without calibration; geometry driven by *atomization quality*
  is not.
- The pressure-swirl internal closure is simplified relative to a full
  Giffen-Muraszewski solution — `V_AND_V_ROADMAP.md` already flags this and schedules the
  Lacava verification. Geometric targets from it are screening-level until that closes.

## 9. CFD ingestion contract (write now, implement Phase 6)

Because CFD will eventually be available (decision #3), the export requirements are
fixed **now**, so the CFD is run in a way that is usable later. Re-running CFD because
the export lacked face fluxes would be an expensive, avoidable mistake.

Required per-cell export: cell id, centroid `(x,y,z)`, volume, density, velocity vector,
temperature, pressure, mixture fraction and/or equivalence ratio, major species mass
fractions, and NO/NO2 mass fractions.

Required per-face export (**the part that is easy to forget**): face id, owner and
neighbour cell ids, face area, outward normal, and face mass flux
`mdot_f = rho_f (u_f . n_f) A_f`. Reconstructing fluxes from cell-centred velocities
after the fact is lossy and will not close the mass balance.

Required Lagrangian spray export, if the CFD is liquid-fueled: parcel positions,
diameters, temperatures, velocities, and evaporation rates, binned per zone.

Ingestion pipeline (paper Section 2.2.5): interpolate to a structured auxiliary mesh →
k-means on the `[T, phi]` feature vector → map clusters to space → build the
reactor-to-reactor mass-flow matrix `M` by summing positive face fluxes across cluster
boundaries → minimum-norm least-squares mass correction → assign reactor types
(evaporator at the fuel inlet, PFR at the exit/largest-axial-extent zone, PSR elsewhere)
→ emit a `CombustorNetwork`.

Calibration (their Eq. 24): weighted relative-error objective over
`[T_out, Y_NO,out, Y_NO2,out, phi_mixer]`, with mass flows, spray-path reactor volumes,
spray path lengths, and spray propagation fractions as variables, bounds of ±30% on
flows and volumes and ±50% on path lengths, mass re-balanced after every optimizer
proposal, and infeasible evaluations penalized. Solved with CMA-ES and a genetic
algorithm, best-of-two retained. Cluster-count sensitivity (5–11) is reproduced as a
convergence study.

---

## 10. Public API sketch

```python
from fuelnozzle.crn import (
    ActiveFuel, AirSplit, CombustorGeometry, MechanismRegistry, MechanismSpec,
    RQLTemplate, LPPTemplate, LDITemplate, StagedTemplate,
    SprayBoundaryPolicy, DropletSettings, AutoignitionSettings,
    run_combustor_study,
)

result = run_combustor_study(
    operating_points=points,                 # existing OperatingPoint + combustor fields
    architecture=RQLTemplate(
        rich_zone_equivalence_ratio=1.5,
        quench_stages=5,
        quench_time_s=1.5e-3,
        air_split=AirSplit(dome=0.20, primary=0.15, quench=0.35,
                           dilution=0.20, cooling=0.10),
        reactor_volumes_m3={...},
        recirculation_fractions={...},
    ),
    mechanisms=MechanismRegistry(
        jet_a=MechanismSpec(path="mech/hychem_posf10325_nox.yaml", ...),
        lng=MechanismSpec(path="mech/gri30.yaml", ...),
    ),
    lng_geometry=..., jet_a_geometry=..., jet_a_properties=...,   # existing objects
    composition=..., feed_line=...,                               # existing objects
    droplet_settings=DropletSettings(spray_policy=SprayBoundaryPolicy.NOZZLE_SMD),
    autoignition=AutoignitionSettings(minimum_margin=4.0),
)
```

Returned per operating point: per-reactor temperature, pressure, equivalence ratio,
composition, residence time, and heat loss; droplet class histories with evaporation
fraction per reactor; exit temperature, NO/NO2 in ppmv dry corrected to 15% O2 and as
EI (g/kg fuel); CO/UHC as **uncalibrated diagnostics**; liquid carryover; autoignition
margin per passage; heat-sink budget; and the full warning chain, preserving the
existing package's provenance discipline. Every result nests the nozzle result that
produced it, exactly as `Tier3FlashSpray` nests Tiers 2/1/0.

---

## 11. Task list

Tasks are numbered for selective approval. Each has an acceptance criterion.

### Phase 0 — Environment and scaffolding

- **0.1** Add `cantera >=3.1` (conda-forge) to `pixi.toml`. Add `matplotlib` and
  `pyyaml`. Defer `scikit-learn` to Phase 6 and `pymoo` to Phase 8.
  *Accept:* `pixi run python -c "import cantera; print(cantera.__version__)"` works on
  osx-arm64; `pixi run test` still passes all 21 existing tests.
- **0.2** Create `src/fuelnozzle/crn/` with `__init__.py` and module stubs. Extend
  `pyproject.toml` dependencies. Keep `ruff` clean at line-length 100.
- **0.3** Create `docs/crn_technical_reference.tex` skeleton reusing the existing
  preamble, `\code`/`\srcline`/`\warningbox` macros, and beginner-oriented voice.
  Add a `doc-crn` pixi task and make `doc` build both documents with tectonic.
  *Accept:* `pixi run doc` produces both PDFs with no LaTeX errors.
- **0.4** Add a `mech/` directory with a README stating the mechanism contract, and
  `.gitignore` entries for large mechanism files if the user prefers them out of git.
- **0.5** Convert `A2NOx_skeletal.cti` and `A2NTCfast_ske.cti` to YAML with `cti2yaml`.
  Record source URL, download date, SHA-256, Cantera converter version, and species/
  reaction counts in `mech/README.md`.
  *Accept:* both YAML files load via `ct.Solution`; species counts are 71 and 47;
  `A2NOx_skeletal` reports a non-zero NO production rate for a stoichiometric
  POSF10325/air mixture at 1800 K; and `A2NTCfast_ske` reproduces an NTC region
  **at 5-10 atm** (amended — see deviation D-004; NTC washes out at combustor pressure,
  so testing at 20 atm gives a false negative).

### Phase 1 — Chemistry and stream layer

- **1.1** Read the Cantera reactor documentation against the installed version; record
  the exact classes, constructor signatures, and methods used in a
  `docs/crn_technical_reference.tex` appendix, with the version pinned.
  *Accept:* every Cantera call used later appears in that table.
- **1.2** `chemistry.py`: `MechanismSpec`, `MechanismRegistry`, `ActiveFuel`, cached
  `Solution` loading, surrogate composition handling, stoichiometric AFR, Bilger
  mixture fraction, element-based equivalence ratio, ppmv-dry-corrected-to-15%-O2 and
  EI conversions.
  *Accept:* unit tests for stoichiometric AFR of CH4 and a Jet-A surrogate against hand
  calculations; `phi` recovered exactly for a known premixed state.
- **1.3** Mechanism validation: required species present, nitrogen chemistry present
  when NOx is requested (**error** if absent), LNG components represented, declared
  validity range vs the operating point (warning if outside).
  *Accept:* a mechanism with no `NO` raises; GRI-3.0 at 30 atm warns.
- **1.4** `streams.py`: `AirSplit`, air property state at `T3`/`P3`, split validation
  (fractions sum to 1 within tolerance), and the mass-flow bookkeeping.

### Phase 2 — Droplet physics

- **2.1** `droplets.py`: `DropletClass` state, TAB breakup integration, Frossling
  diffusion-limited evaporation, boiling-limited evaporation, infinite-conductivity
  heating, with the energy accounting of Section 3.6.
- **2.2** Liquid property adapters: Jet-A from the extended property table; LNG from
  `CoolPropLNGProvider`. Both behind one protocol so the droplet solver is fuel-agnostic.
- **2.3** Gas-side transport property evaluation at the reference temperature
  `Tbar = (T_gas + 2 T_d)/3` from the Cantera solution.
- **2.4** Extend `JetAPropertyTable` with optional `liquid_cp_j_kg_k`,
  `latent_heat_j_kg`, `molecular_weight_kg_mol`, `boiling_point_k`. All optional;
  existing constructors unchanged.
  *Accept:* all 21 existing tests pass with no edits.
- **2.5** `spray_source.py`: build `SprayBoundary` from `PressureSwirlResult` or
  `Tier3FlashSpray`, implement the three `r0` policies, honour the existing
  calibration-gating discipline (no invented SMD), skip TAB for flashing regimes.
- **2.6** Rosin-Rammler discretization into N droplet classes.

### Phase 3 — Reactors and network

- **3.1** `reactors.py`: thin wrappers over the Cantera reactor types in Section 4.1,
  carrying volume, residence time, heat-loss specification, resident droplet classes,
  and spray path length.
- **3.2** `network.py`: `CombustorNetwork` graph, split-matrix assembly, minimum-norm
  mass-balance correction, Cantera `ReactorNet` construction with `MassFlowController`s
  and a `PressureController`, steady-state solve with preconditioning and a documented
  fallback.
  *Accept:* mass closes to solver tolerance on every reactor; a deliberately
  inconsistent split is corrected minimally and reported.
- **3.3** Operator-split outer iteration (Section 4.3) with under-relaxation, an
  iteration cap, and non-convergence reported as a warning rather than a silent result.
- **3.4** Initialization strategy: equilibrium or unsteady-ignition first guess, so the
  network does not converge to the trivial extinguished solution. Extinguished
  solutions are detected and reported explicitly.

### Phase 4 — Architecture templates

- **4.1** `templates.py`: RQL with a staged quench, LPP, LDI, staged pilot/main, custom.
- **4.2** Quench sub-model with an air-addition schedule; `N_quench` and `tau_q`
  sensitivity as a first-class output.
- **4.3** Mechanism merge *checker* (duplicate reactions, thermo inconsistency, species
  collisions). Reports; does not merge.

### Phase 5 — Design-question modules

- **5.1** `autoignition.py`: ignition-delay tabulation and interpolation, mixture-state
  computation including LNG flash cooling, residence time, margin, flashback screen,
  evaporation completeness.
  *Accept:* CH4 and a Jet-A surrogate ignition delays reproduce published Arrhenius
  trends; the LNG flash-cooling temperature drop is verified by an independent enthalpy
  balance.
- **5.2** `thermal.py`: heat-sink budget, inverse solve for the fuel temperature that
  achieves a target flash regime, four-constraint feasible-window map.
- **5.3** Add optional combustor fields to `OperatingPoint` (air flow or overall `phi`,
  liner pressure loss, active fuel), all defaulted.
  *Accept:* existing tests and `examples/run_study.py` unchanged and passing.
- **5.7** Idle-circuit thermal screen (Section 8.2.2): Jet-A coking-limit wall
  temperature while LNG burns, LNG vapour-lock and restartability while Jet-A burns, and
  cross-talk heat leak between the two circuits in one nozzle body. Reuses the
  conduction/insulation machinery in `feed.py`.
  *Accept:* a packaging that violates the Jet-A coking limit is flagged as infeasible.
- **5.4** `emissions.py`: EI, ppmv corrections, ICAO LTO `Dp/Foo` when thrust and
  time-in-mode are supplied. CO/UHC clearly labeled uncalibrated.
- **5.5** `combustor_study.py`: orchestration across operating points, provider and
  mechanism reuse, warning aggregation, architecture-comparison report.
- **5.6** `examples/run_combustor_study.py`: LNG cruise on LPP and LDI, Jet-A LTO on
  RQL and LDI, printing the architecture trade table and the feasible thermal window.

### Phase 6 — CFD ingestion and calibration (deferred until CFD exists)

- **6.1** Write and circulate the CFD export contract of Section 9 **before** the CFD
  campaign runs. This task is due early even though the rest of Phase 6 is deferred.
- **6.2** `cfd_ingest.py`: structured interpolation, k-means on `[T, phi]`, cluster-to-
  space mapping, face-flux-based `M` assembly, mass correction, reactor-type assignment.
- **6.3** `calibrate.py`: the Eq. 24 objective, sensitivity-based mass-flow variable
  selection, bounds, re-balancing after each proposal, CMA-ES and GA, best-of-two.
- **6.4** Cluster-count convergence study (5–11), reproducing the paper's saturation
  behaviour beyond 9 clusters.

### Phase 7 — Verification, validation, documentation

- **7.1** Unit verification (Section 12).
- **7.2** Energy- and mass-closure tests, including the Section 3.6 enthalpy accounting.
- **7.3** Implement **Tier V1** paper tests (V1.1-V1.6, Section 12.2): droplet timescale
  separation at 350 K and 900 K, D2-law shape, the Table 1 mass-correction fixture, and
  the residence-time scaling identity.
  *Accept:* V1.1 and V1.2 pass within a declared tolerance; this is the gate for Phase 2.
- **7.8** Implement **Tier V2** trend check: the seven-reactor LFRN-vs-GFRN signature
  test of Section 12.2, with assumed volumes explicitly declared in the fixture. This is
  the exit criterion for all paper-comparison work.
  Also cross-check the two Jet-A mechanisms (Section 4.4) against each other in their
  overlapping high-temperature validity range.
  *Accept:* the exit-temperature-agrees / NOx-diverges asymmetry is reproduced.
- **7.9** Record in `mech/README.md` and the LaTeX document that our chemistry differs
  from the paper's, that the paper does not state its NOx submodel, and that absolute
  NOx agreement with the paper is therefore not claimed.
- **7.4** Grid/parameter convergence: PFR sub-reactor count, droplet ODE tolerance,
  outer-iteration tolerance, ignition-delay table resolution.
- **7.5** Complete `docs/crn_technical_reference.tex` — every equation, every algorithm,
  every source-line map, plain English throughout, following the existing document's
  structure and its habit of separating physical law from numerical approximation from
  empirical calibration.
- **7.6** Update `docs/modeling.md`, `docs/V_AND_V_ROADMAP.md`, and `README.md`.
- **7.7** Maintain the implementation log in Section 14 of this file as work proceeds.

---

### Phase 8 — Dual-fuel design optimization (**the payoff phase**)

Depends on Phases 1-5. Independent of Phase 6 — this does **not** wait for CFD.

- **8.1** `design.py`: `DesignVariable` with the Class A (shared combustor) / Class A2
  (per-fuel injector) / Class B (schedulable) distinction of Section 8.5, bounds, the
  `a_J + a_L` air-passage coupling, the dome packaging constraint, and a `DesignVector`
  mapping to both nozzle geometries, the liner air split, the architecture *pair*, and
  the per-segment schedule.
  *Accept:* a design vector round-trips to a fully specified study and back.
- **8.2** `evaluate.py`: one design vector plus a mission set to a full multi-point
  result. Caches property providers, `Solution` objects, and ignition-delay tables across
  evaluations. Parallel across designs with `multiprocessing`/`joblib`.
  *Accept:* repeated evaluation of one design is bitwise identical; N designs on N cores
  scale near-linearly.
- **8.3** `objectives.py`: the four objectives and the hard-constraint set of Section 8.3,
  with **feasibility-first ranking** rather than weighted penalties.
  *Accept:* a design violating autoignition margin never outranks a feasible design,
  whatever its NOx.
- **8.4** Stage 1 sweeps: one-at-a-time and Sobol/Morris screening, producing a
  sensitivity ranking over Class A and Class B variables.
  *Accept:* the ranking is stable under sample-size doubling.
- **8.5** Stage 2 `optimize.py`: NSGA-II via `pymoo` over the reduced variable set,
  parallel evaluation, feasibility-first constraints, checkpoint/restart (runs will take
  hours).
  *Accept:* the front is non-dominated and stable under re-run with a different seed.
- **8.6** Stage 4 architecture-**pair** comparison: (Jet-A RQL, LNG LPP), (Jet-A LDI, LNG
  LPP), (Jet-A RQL, LNG LDI), (Jet-A LDI, LNG LDI), each crossed with the four air-split
  mechanisms of Section 8.3. Compare the resulting Pareto fronts.
- **8.7** **Cost-of-shared-liner chart** (Section 8.8): Jet-A-optimal, LNG-optimal, and
  shared-liner optima each with fully optimized per-fuel injectors, reporting `dA` and
  `dL` per architecture pair and air-split mechanism. Plus **per-path** autoignition
  margin maps in `(premix length, T_fuel)`.
  *Accept:* `dA` and `dL` are each computed within one fuel and one mechanism, and the
  report states this.
- **8.7c** **Answer 1 (Section 8.10.1):** head-end air sweep for Jet-A LTO, reporting
  `f_dome` and `phi_rich` at minimum LTO NOx, with sensitivity to `N_quench` and `tau_q`.
  *Accept:* a NOx minimum appears near `phi_rich ~ 1.4-1.6`; if not, the result is treated
  as a tool defect and investigated before being reported.
- **8.7d** **Answer 2 (Section 8.10.2):** LNG cruise lean sweep, reporting the lean limit
  as a **bracket** from PSR extinction continuation, CO rise, and exit-uniformity
  degradation, with the governing criterion named.
  *Accept:* the limit is never reported as a single number; PSR extinction is computed by
  continuation, not from a correlation.
- **8.7e** **Coupled answer (Section 8.10.3):** solve for the triple `(f_dome, a_J, a_L)`
  that simultaneously delivers rich Jet-A LTO and lean LNG cruise, with the residual
  `dA`/`dL` penalty.
- **8.7f** **Answer 3 (Section 8.11):** per-path target geometry reporting. Add flow
  number `FN` and atomizer constant `K` to `PressureSwirlResult`; add an `L/D` sweep
  reporting flash-onset location with explicit sensitivity to the Tier 2 relaxation time
  `tau`; evaluate Jet-A atomization quality at idle as well as takeoff and report whether
  a single simplex covers the turndown range; emit feasible regions with the driving
  constraint named per feature.
  *Accept:* every geometric target carries its driver and bounding constraint; no SMD-
  dependent target is emitted without a calibration; the `L/D` result is always
  accompanied by its `tau` sensitivity.
- **8.7b** Passive-air-staging study (Section 8.2.1): sweep `a_J`/`a_L` and the
  idle-passage mixing fraction; report near-field `phi` per fuel with explicit sensitivity.
  *Accept:* no conclusion is reported that holds only at one assumed mixing fraction.
- **8.8** Stage 3 surrogate-assisted refinement near the Pareto knee. *Optional.*
- **8.9** Robustness gate: re-check every reported conclusion against the Stage 1
  sensitivity range and against a mechanism substitution. Conclusions that flip are
  removed from the design guidance and recorded as unresolved.
- **8.10** `examples/run_dual_fuel_design.py`: end-to-end sweep, Pareto search, and
  cost-of-dual-fuel chart for one mission set.

## 12. Verification and validation plan

Follows the discipline already established in `docs/V_AND_V_ROADMAP.md`: verification
that the software solves the equations correctly, kept strictly separate from validation
that the equations describe reality.

### 12.1 Verification (no experiments needed)

| Check | Criterion |
|---|---|
| Single PSR, long residence time | Approaches Cantera equilibrium to < 1 K and < 1% on major species |
| PSR extinction | Reproduces the classical S-curve; blowout residence time is grid-independent |
| PFR as N sub-PSRs | Exit state converges with N; report the N needed for 0.5% |
| Mass balance | Closes on every reactor to solver tolerance; element balance C/H/O/N closes |
| Energy balance | Total enthalpy in = out + wall loss, to round-off, on a non-reacting evaporating case |
| D2-law | Diffusion-limited branch recovers `d^2` linear decay in a quiescent, constant-property limit |
| Boiling branch | Analytical `Q/L` rate recovered for a fixed-property sphere |
| TAB | Paper Fig. 3: 150 um, `T_d`=300 K, `T_gas`=350 K, `u_d`=112 m/s gives breakup near 2e-5 s and evaporation near 6e-3 s; `tau_evap/tau_breakup` ~133 at 350 K and ~21 at 900 K |
| Mass correction | Returns the true minimum-norm solution on a constructed problem with a known answer |
| Zero-fuel limit | Network reduces to air passing through unchanged |
| Zero-evaporation limit | LFRN reduces to a gaseous CRN, reproducing the GFRN behaviour |
| Reproducibility | Bitwise-stable results for fixed settings and mechanism version |

### 12.2 Validation against John et al. (the LFRN paper)

**Purpose of this comparison, per the 2026-08-28 direction: get the trends and the
conclusions right, to prove the implementation is going in the right direction.** It is
a directional check on our implementation of their method, not a validation campaign and
not a figure-reproduction exercise. Reproducing their absolute NOx is neither required
nor possible — they publish reactor *mass flows* (Tables 1 and 2) but not reactor
volumes, not the connectivity matrix, not spray path lengths `s_i`, and not propagation
fractions `eta_ij`, so their calibrated network is not reconstructable.

Two tiers, both cheap, and then we move on to the dual-fuel work that is the actual
point of this tool.

#### Tier V1 — Sub-model unit tests from fully-specified paper data

These are not figure reproductions; they are regression tests with published numbers on
both sides, validating the droplet sub-models and the mass-correction algorithm — the
majority of the new physics, and the part most likely to be silently wrong.

| Test | Paper source | Given | Assert |
|---|---|---|---|
| **V1.1 Breakup/evaporation timescales** | Fig. 3 and Sec. 3.1 text | `d0`=150 um, `T_d`=300 K, `T_gas`=350 K, `u_d`=112 m/s | Breakup collapses 150 um to ~1 um at `t ~ 2e-5` s; evaporation near-complete at `t ~ 6e-3` s; `tau_evap/tau_breakup ~ 133` measured at one-tenth of initial size |
| **V1.2 Timescale ratio at high temperature** | Sec. 3.1 text | Same, `T_gas`=900 K | Ratio falls to `~21`, still exceeding one order of magnitude |
| **V1.3 D2-law shape** | Sec. 3.1 ("following the classical D2-law") | Evaporation-only case | `d^2` decays linearly in time to within a declared tolerance |
| **V1.4 Mass-correction fixture** | Table 1 | Seven reactors with published `mdot_in`, `mdot_out`, and imbalances `eps` from +0.23 to -4.58 g/s (up to -9.31% on PSR-6) | Minimum-norm correction closes every reactor balance; the correction is the true min-norm solution; corrected flows sit within the published raw-to-balanced spread |
| **V1.5 Residence-time scaling** | Sec. 2.2.5 text under Table 2 | `tau = rho V / mdot` | A -26.4% inflow change gives +35.9% in `tau`; +25.8% gives -20.5%. Trivial arithmetic, but it pins our `tau` definition to theirs |
| **V1.6 Calibration bounds and weights** | Sec. 2.2.5 | `w = [5, 4, 2, 2]` on `[T_out, Y_NO, Y_NO2, phi_m]`; ±30% on flows/volumes, ±50% on path lengths; failed evaluations penalized at `1e6` | Phase 6 optimizer reproduces these settings exactly as defaults |

V1.1 and V1.2 are the strongest tests available anywhere in this plan: fully specified
inputs, published outputs, and they exercise TAB, Frossling, and the heating model
*together*. They are the acceptance gate for Phase 2.

#### Tier V2 — Trend and conclusion check (the one that matters)

The paper's central claim is a *mechanism*, and the mechanism is testable from numbers
quoted in the running text rather than read off figures. Build a seven-reactor network
matching their Fig. 2 / Fig. 5 description — evaporator, mixer, flame-1, flame-2,
recirculation-1, recirculation-2, PFR — using the Table 2 mass flows, with volumes
assumed and declared. Then run the identical network twice: once with the spray models
active (LFRN) and once with the fuel introduced as vapour at the inlet (GFRN).

Assert the signature, not the absolute values:

| Quantity | Paper (text-quoted) | Test |
|---|---|---|
| GFRN equivalence-ratio spread | 0.652 to 0.653 across an 11-cluster network | GFRN `phi` variance is near zero |
| LFRN equivalence-ratio spread | Locally rich near the nozzle, diluting toward the flame | LFRN `phi` spans a wide range, rich in evaporator/mixer |
| GFRN peak temperature | `~1800` K | GFRN peak well below LFRN peak |
| LFRN peak temperature | `>2200` K (CFD peak 2180 K) | LFRN peak exceeds GFRN by several hundred K |
| Outlet NO | GFRN "more than an order of magnitude lower" than CFD | GFRN/LFRN outlet NO ratio below ~0.1 |
| Exit temperature | Both within 3% of CFD | LFRN and GFRN exit temperatures agree with each other within a few percent |
| Cluster-count trend | Error decreases 5 to 9 clusters, saturates beyond 9 | Monotone-then-flat trend reproduced qualitatively |

That last row of the temperature comparison is the sharpest part of the paper: **both
methods get exit temperature right and only one gets NOx right.** If our implementation
reproduces that asymmetry, the spray-heterogeneity mechanism is working. If it does not,
something is wrong regardless of what the absolute numbers say.

Their reported error bands, for reference in the test docstrings: LFRN NO 0.54% at
baseline growing to 26.3% at 500 K, and 0.29–0.67% across the fuel-flow sweep; GFRN NO
54–91% low throughout. NO2 is weaker for both (LFRN within 46% / 29–41%, GFRN 52–87%),
which is itself worth encoding as an expectation that **NO2 is the less trustworthy
output**.

#### Tier V3 — Digitized figure comparison: **dropped**

Figs. 6, 7, 8, and 11 exist only as plots and would require digitization.
**Explicitly descoped on 2026-08-28.** If it is ever wanted, it is done under the
transcription discipline already written in `docs/V_AND_V_ROADMAP.md` Section 7.

#### What is explicitly **not** claimed

- We will **not** claim to reproduce their absolute NOx, because volumes, connectivity,
  spray lengths, and propagation fractions are unpublished.
- We will **not** match their chemistry exactly. They used Cantera 2.6 with HyChem +
  USC Mech II — and USC Mech II has no nitrogen chemistry, while **the paper never states
  which NOx submodel they appended.** We use HyChem A2 + Glarborg (Section 4.4), which is
  a defensible and probably better choice, but it is *not* their mechanism.
- Runtime comparison (their Table 3) is hardware-dependent and is recorded as an
  order-of-magnitude sanity note, never as a test.

**Exit criterion for paper work.** Once V1 passes and V2 reproduces the
exit-temperature-agrees / NOx-diverges asymmetry, the paper comparison is **done**. No
further effort is spent on it. Everything after that goes into Section 8.

### 12.3 Validation against other evidence

Ordered by availability:

1. **Against your own CFD** once available (decision #3): the calibration and holdout
   protocol of Section 9, with train/holdout splits declared before fitting.
2. **Against rig or engine data** if obtainable: exit temperature, NOx EI, and
   pattern factor.
3. **Against the cryogenic spray datasets** already registered in
   `V_AND_V_ROADMAP.md` Section 6.3 (DaRUS-2527, DaRUS-2076), for the flash-boiling
   droplet branch, which the John paper cannot validate because it is Jet-A only.
4. **Anchoring checks** that cost nothing and catch gross errors: `EI_NOx` scaling
   roughly as `P^0.5` and exponentially in flame temperature; RQL NOx minimizing near
   `phi_rich` ~1.4–1.6; lean-burn NOx falling steeply with `phi` until CO rises.

**Claim boundary, to be stated in the README and the LaTeX document:**

> The reactor-network extension is a reduced-order design-space screening tool for
> comparing combustor architectures and fuel-nozzle thermal strategies. Absolute NOx
> predictions are not validated for this hardware and must not be used as certification
> or compliance estimates. Its intended use is *relative* comparison — RQL vs lean-burn,
> LNG vs Jet-A, one fuel temperature vs another — under a single set of stated
> assumptions. CO, UHC, and soot outputs are uncalibrated diagnostics.

---

## 13. Risks

| # | Risk | Mitigation |
|---|---|---|
| R1 | ~~Supplied Jet-A mechanism lacks nitrogen chemistry~~ **RESOLVED 2026-08-28**: `A2NOx_skeletal.cti` (71 sp, HyChem A2 v2.0 + Glarborg NOx) confirmed to carry full N chemistry | Guard retained anyway (hard error in Task 1.3) because it now protects the *LNG* mechanism, which is not yet supplied. |
| R1b | **CONFIRMED AND QUANTIFIED 2026-08-28.** The NOx-bearing A2 mechanism is high-temperature only. Measured error in Jet-A ignition delay at premixer conditions: **510x at 700 K and 71x at 800 K (20 atm, phi=0.5)** — it would declare an unsafe premixer safe | Declared two-mechanism split by `role` (Section 4.4): `A2NOx_skeletal` for the network, `A2NTCfast_ske` for ignition delay. Cross-checked in the overlapping high-temperature range (Task 7.8); disagreement reported, never averaged. |
| R1f | **New.** `A2NTCfast_ske` is the *fast* NTC variant of a **hypothetical** A2 fuel, not a validated real-Jet-A low-temperature model; Stanford publishes fast and slow variants that bracket the uncertainty | Fast NTC gives the *shorter* delay and is therefore **conservative** for a safety screen, so current use is safe but may be pessimistic. Report the margin as a bracket across both variants once the slow variant is obtained (open item O-004). |
| R1c | **New.** Both supplied files are CTI, which Cantera 3.x cannot parse | One-time `cti2yaml` conversion with URL, date, SHA-256, and converter version recorded (Task 0.5). |
| R1d | ~~No LNG mechanism supplied~~ **RESOLVED 2026-08-28**: GRI-Mech 3.0, ships with Cantera as `gri30.yaml` | Pressure-validity warning retained (R6). |
| R1e | **New.** Jet-A and LNG use mechanisms of different size, provenance, and validity, so cross-fuel *absolute* NOx is not strictly comparable | Every dual-fuel comparison reports both mechanisms and is labeled trend-level. Cross-fuel conclusions must survive a mechanism-substitution check before being reported as design guidance. |
| R2 | No spray calibration exists, so SMD is unavailable and atomization QoIs cannot be reported | Inherit the existing package's refusal to invent an SMD. Offer the paper-faithful `r0 = r_nozzle` + TAB path as an explicitly-labeled fallback so the network can still run, with the atomization QoI suppressed. |
| R3 | Simultaneous `ReactorNet` steady solve fails to converge with recirculation and stiff chemistry | Preconditioning; unsteady-ignition initialization; documented sequential-modular fallback; explicit non-convergence reporting rather than a silent bad answer. |
| R4 | Operator-split iteration converges slowly in strongly coupled evaporator zones | Under-relaxation; fall back to the extensible-reactor coupled formulation (Section 4.3) for those reactors only. |
| R5 | Template flow splits are guesses, so architecture comparison is guess-vs-guess | Every comparison ships with a split-sensitivity study; conclusions are only reported where they survive the sensitivity range. This is why Phase 6 CFD ingestion matters. |
| R6 | GRI-Mech 3.0 used far outside its validated pressure range | Validity-range warnings; recommend FFCM-2 or Aramco+N at pressure; record the mechanism in every result. |
| R7 | The flash-boiling droplet branch has no cryogenic validation data | Label it clearly as model-form extension; connect it to the LN2 datasets already registered in `V_AND_V_ROADMAP.md` Section 6.3 (DaRUS-2527, DaRUS-2076). |
| R8 | Scope creep into a general CFD-replacement CRN framework | The four QoIs of decision #2 bound the work. Anything not serving atomization, mixing, temperature, or NOx is deferred. Paper comparison has an explicit exit criterion (Section 12.2). |
| R9 | **Highest design-level risk.** An optimizer run on a partly-validated reduced-order model produces a confident "optimal design" that nobody should believe | Staged strategy (Section 8.4): sensitivity first, optimization second. Feasibility-first constraints. Task 8.9 robustness gate deletes any conclusion that does not survive the sensitivity range or a mechanism substitution. Outputs are labeled candidates for CFD and rig test, never final designs. |
| R10 | Optimization cost: minutes per design times thousands of designs | Stage 1 dimension reduction before Stage 2; parallel evaluation; checkpoint/restart; surrogate refinement only near the knee. If still infeasible, reduce the mission set to the two binding points (Jet-A takeoff, LNG cruise). |
| R11 | The dual-fuel compromise turns out to be dominated by one architecture pair at every point, making the study anticlimactic | That is a *result*, not a failure, and worth knowing early. Stage 1 sweeps reveal it before the expensive Stage 2 runs. |
| R12 | **The passive-air-staging idea (Section 8.2.1) depends on unfueled air from the idle passage staying segregated from the near field.** A CRN cannot settle this; it only propagates the user's assumed mixing fraction | Always reported with a sensitivity across that fraction, never at one value. Flagged as the highest-priority question for the Phase 6 CFD. If it fails, options 1, 3, and 4 of Section 8.3 remain. |
| R13 | Idle-circuit coking or vapour lock kills an otherwise optimal packaging late in the study | Screened early as a hard constraint (Task 5.7), not discovered at the end. |

---

## 14. Implementation log

**The detailed running log lives in [`CRN_IMPLEMENTATION_LOG.md`](CRN_IMPLEMENTATION_LOG.md).**
It records per-task completion with acceptance evidence, a numbered deviations table
(plan said / what happened / why / consequence), open items, and a cumulative
verification record. This file stays the statement of design intent; that file is the
record of execution.

The summary table below tracks plan-level decisions only.

*(Appended as work proceeds — date, task number, what was built, what was verified,
what was deviated from and why. Empty until the plan is approved.)*

| Date | Task | Change | Verification | Notes |
|---|---|---|---|---|
| 2026-08-28 | — | Plan created. | — | **Approved 2026-08-28; implementation started at Phase 0.** |
| 2026-08-28 | — | Confirmed scope answers 1-3; assumptions A1-A3 confirmed by user. | — | Steady + constant pressure approved. |
| 2026-08-28 | 4.4 | Inspected `A2NOx_skeletal.cti` (71 sp, HyChem A2 v2.0 + Glarborg NOx) and `A2NTCfast_ske.cti` (47 sp, NTC, no N). Adopted two-mechanism split by role. | Species lists read from source files. | R1 resolved; R1b/R1c/R1d opened. |
| 2026-08-28 | 4.4 | GRI-Mech 3.0 confirmed for LNG (ships with Cantera as `gri30.yaml`). Pressure-validity caveat recorded. | — | R1d resolved; R1e opened on cross-fuel comparability. |
| 2026-08-28 | 8.11 | Added Answer 3: per-path target nozzle geometry. Jet-A gains `FN` and `K` design groups plus a turndown check; LNG gains an `L/D` flash-onset sweep. Recorded that `L/D` accuracy hinges on the Tier 2 relaxation time and that one flashing test would pin it. | — | Easiest of the three answers: separate hardware means no cross-fuel compromise. |
| 2026-08-28 | 8.10 | Named the two primary outputs: head-end air fraction for Jet-A LTO, and the LNG cruise lean limit. Added exit temperature uniformity as objective 5, the three-criterion lean-limit bracket, and the pattern-factor limitation caveat. | — | Confirms air split is the crux; the two answers are coupled through `(f_dome, a_J, a_L)`. |
| 2026-08-28 | 8 | **Separate fuel paths per fuel, shared combustor.** Rewrote Section 8: injector-level compromise removed, air split becomes the crux; added inactive-circuit air lever (8.2.1), idle-circuit thermal screen (8.2.2), four air-split reconciliation mechanisms (8.3), and the LNG-exploitation chain (8.4). Headline metric became cost-of-shared-liner. | — | Supersedes the earlier "premixing length is the crux" framing, which assumed shared hardware. R12/R13 opened. |
| 2026-08-28 | 8 | **Reframed: primary deliverable is the optimal dual-fuel design, not architecture evaluation.** Added Section 8 (design variables, objectives, staged optimization, cost-of-dual-fuel chart) and Phase 8 tasks. | — | Paper validation demoted to a directional check with an explicit exit criterion; Tier V3 dropped. |
| 2026-08-28 | 12.2 | Added tiered John et al. validation plan (V1 exact / V2 behavioural / V3 digitized) after establishing that reactor volumes, connectivity, spray lengths, and propagation fractions are unpublished. | — | Absolute NOx reproduction explicitly not claimed. |
