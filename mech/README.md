# Kinetic mechanisms

Mechanisms used by `fuelnozzle.crn`. See `docs/CRN_PLAN.md` Section 4.4 for the role
each one plays and `docs/CRN_IMPLEMENTATION_LOG.md` for the verification record.

`source/` holds the files exactly as downloaded, immutable. The `.yaml` files at this
level are the converted, Cantera-loadable versions and are what the code reads.

## Registry

| File | Role | Species | Reactions | N chemistry | Low-T / NTC |
|---|---|---|---|---|---|
| `A2NOx_skeletal.yaml` | Jet-A **network + NOx** | 71 | 538 | Yes (Glarborg) | **No** |
| `A2NTCfast_ske.yaml` | Jet-A **ignition delay** | 47 | 247 | No (N2 only) | Yes (fast variant) |
| `gri30.yaml` | LNG **network + NOx** | 53 | 325 | Yes | n/a |

`gri30.yaml` ships with Cantera and is not vendored here; its version is pinned by the
Cantera version recorded in the implementation log.

## Provenance

### A2NOx_skeletal

- Source: <https://web.stanford.edu/group/haiwanglab/HyChem/download/cantera/A2NOx_skeletal.cti>
- Downloaded: 2026-08-28, 116,044 bytes
- SHA-256 (source CTI): `5618ae19b95d7f82a8814ddde5f7f00fcceb8f2f980372f993246efbb76aeeaa`
- Skeletal reduction by Tianfeng Lu (University of Connecticut, 2018)
- Combines the HyChem model for Jet-A POSF10325 (Version 2.0) with the NOx chemistry of
  Glarborg et al., *Prog. Energy Combust. Sci.* **67** (2018) 31-68
- Fuel species: `POSF10325` (C11H22)
- Converted with `cti2yaml` from Cantera 3.2.0

### A2NTCfast_ske

- Source: <https://web.stanford.edu/group/haiwanglab/HyChem/download/cantera/A2NTCfast_ske.cti>
- Downloaded: 2026-08-28, 68,515 bytes
- SHA-256 (source CTI): `cd2a40126fcac116032518c55ed28f49a684c181ba673dca2ddc009122c47427`
- Skeletal reduction by Tianfeng Lu (2018). Header describes it as *"a 47-species skeletal
  model for a **hypothetical** A2 fuel with **fast** NTC"* — 8 pyrolysis steps and 8
  low-temperature oxidation steps on a USC Mech IIa base.
- **This is a sensitivity variant, not a validated real-Jet-A low-temperature model.**
  Stanford publishes fast and slow NTC variants that together bracket the real
  uncertainty in Jet-A low-temperature chemistry. See "Open question" below.
- Converted with `cti2yaml` from Cantera 3.2.0

## Why two Jet-A mechanisms

Neither file has both NOx and low-temperature chemistry, and the tool needs both for
different sub-models. The split is **not** a stylistic choice; it is quantitatively
necessary, and this is the evidence.

### Ignition delay at premixer-relevant conditions (phi = 0.5)

Constant-pressure homogeneous ignition, criterion `T > T0 + 400 K`:

| P | T | `A2NOx_skeletal` (high-T only) | `A2NTCfast_ske` (with low-T) | Ratio |
|---|---|---|---|---|
| 20 atm | 700 K | 6.57e+00 s | 1.29e-02 s | **510x** |
| 20 atm | 800 K | 2.15e-01 s | 3.01e-03 s | **71x** |
| 20 atm | 900 K | 1.52e-02 s | 1.64e-03 s | 9.3x |
| 20 atm | 1000 K | 2.04e-03 s | 1.02e-03 s | 2.0x |
| 40 atm | 700 K | 4.17e+00 s | 8.51e-03 s | **490x** |
| 40 atm | 800 K | 1.33e-01 s | 1.55e-03 s | **86x** |
| 40 atm | 900 K | 9.14e-03 s | 6.76e-04 s | 13.5x |
| 40 atm | 1000 K | 1.18e-03 s | 4.08e-04 s | 2.9x |

**LTO `T3` sits at roughly 700-900 K.** Using the high-temperature-only mechanism there
would overpredict ignition delay by one to two orders of magnitude — declaring a
premixing passage safe when it is not. That is a safety-relevant error in the worst
possible direction, and it is why `autoignition.py` must use `A2NTCfast_ske`.

### NTC is pressure-dependent, and "washed out" does not mean "irrelevant"

Ignition delay for `A2NTCfast_ske`, phi = 1, showing where the non-monotonic NTC region
actually exists:

| P | 600 K | 700 K | 800 K | 850 K | 900 K | 950 K | 1000 K |
|---|---|---|---|---|---|---|---|
| 5 atm | 1.05e-1 | 1.62e-2 | **7.37e-3** | 7.76e-3 | 9.97e-3 | 1.02e-2 | 5.99e-3 |
| 10 atm | 6.30e-2 | 8.42e-3 | 2.74e-3 | **2.31e-3** | 2.33e-3 | 2.69e-3 | 2.20e-3 |
| 20 atm | 4.12e-2 | 4.58e-3 | 1.23e-3 | 8.58e-4 | 7.32e-4 | 7.06e-4 | 6.96e-4 |
| 40 atm | 2.95e-2 | 2.64e-3 | 6.12e-4 | 3.79e-4 | 2.79e-4 | 2.36e-4 | 2.15e-4 |

A clear NTC zone (delay *rising* with temperature) exists at 5 atm between 800 and 950 K,
weakens at 10 atm, and disappears by 20-40 atm.

**Do not read the disappearance as "low-temperature chemistry does not matter at
combustor pressure."** The non-monotonic *shape* vanishes, but low-temperature chemistry
still dominates the *absolute* ignition delay at 20-40 atm, as the 71x-510x table above
shows. These are two different statements and conflating them would be a serious error.

## NOx pathway coverage differs between fuels

| Pathway marker | `A2NOx_skeletal` (Jet-A) | `gri30` (LNG) |
|---|---|---|
| Thermal (Zeldovich): `NO`, `NO2`, `N` | Yes | Yes |
| N2O route: `N2O` | Yes | Yes |
| Prompt via `NCN` (modern Glarborg route) | **Yes** | No |
| Prompt via `HCN` (older Fenimore route) | Yes | Yes |
| NNH route: `NNH` | **No** | **Yes** |

**Neither mechanism is a superset of the other.** The Jet-A model carries the modern NCN
prompt-NO route but lacks NNH; GRI-Mech 3.0 carries NNH but uses the older Fenimore
prompt route.

Each happens to carry the pathway most relevant to its own regime — NCN prompt matters in
rich Jet-A zones, NNH matters in lean premixed methane — which is fortunate, but it must
be stated rather than relied on silently. **Consequence: cross-fuel absolute NOx is not
strictly comparable** (plan risk R1e). The `dA`/`dL` metric of plan Section 8.8 is
computed within one fuel and one mechanism and is unaffected.

## Validity caveats enforced at runtime

- **GRI-Mech 3.0** was optimized against data to roughly 10 atm. At 10-40 atm combustor
  pressure this is an extrapolation and NO is the most pressure-sensitive output. The
  mechanism registry warns whenever `P3` exceeds the declared range. Upgrade path if
  absolute LNG NOx becomes decision-critical: FFCM-2 or an AramcoMech variant with an N
  submodel.
- **GRI-Mech 3.0** covers C1-C3. An `LNGComposition` component absent from the mechanism
  **raises** rather than being silently dropped.
- Both Jet-A files use `CH2*` rather than `CH2(S)` naming. Any future merge must
  reconcile this.

## Open question

`A2NTCfast_ske` is the **fast** NTC variant. Stanford also publishes a **slow** NTC
variant. Because the two bracket the real uncertainty in Jet-A low-temperature chemistry,
the autoignition margin should ideally be reported as a **bracket** across both rather
than a single curve.

Fast NTC gives the **shorter** ignition delay and is therefore the **conservative** choice
for an autoignition safety screen, so using it alone is safe but may be pessimistic.
Tracked as open item **O-004** in the implementation log.
