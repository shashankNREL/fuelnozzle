# Fuel-Nozzle Tool: Implementation and V&V Roadmap

Status date: 2026-08-26  
Package version: 0.1.0  
Current regression state: 21 tests passing; Ruff clean before this documentation update.

> Source line numbers in this document refer to the repository state on the status date.
> They are navigation aids, not permanent identifiers. Function names are the durable references.

## 1. Purpose and Claim Boundary

The package is a steady, reduced-order engineering tool for:

- Jet-A simplex pressure-swirl sizing and screening;
- liquid-LNG feed-line thermal/hydraulic calculations;
- LNG phase-risk screening;
- equilibrium and finite-rate flashing calculations;
- preliminary flash-spray classification;
- repeated evaluation over user-specified flight-envelope operating points.

The current defensible claim is:

> The software implementation is regression tested and ready for numerical verification.
> Model-form validation will use CO2/refrigerant critical-flow data and cryogenic LN2 spray
> data. LNG predictions remain extrapolative until direct liquid-methane or representative
> LNG nozzle data are available.

The package is not a certification model, a reacting combustor model, or a substitute for
hardware testing.

## 2. User-Provided Inputs

Each `OperatingPoint` supplies:

- flight-stage name and optional duration;
- `p3_pa` and `t3_k`;
- LNG and Jet-A mass flow rates;
- LNG pump outlet pressure and temperature;
- Jet-A pump outlet pressure and nozzle inlet temperature;
- pressure drop across each fuel nozzle;
- optional multiplier for identical engines/nozzle circuits.

`p3_pa` is currently treated as both combustor inlet pressure and the receiving pressure at
the fuel-nozzle exit. A future interface should split compressor discharge pressure from local
combustor/nozzle back pressure when liner pressure loss matters.

The pressure budget is

```text
P_nozzle,in = P3 + DeltaP_nozzle
DeltaP_feed,available = P_pump,out - P_nozzle,in
```

Implemented by `fuel_pressure_budget` in `src/fuelnozzle/operating.py:55`, with the equations
at lines 62-63. Inputs are never silently altered to close an inconsistent pressure budget.

## 3. What Has Been Implemented

### 3.1 Property layer

`CoolPropLNGProvider` (`src/fuelnozzle/properties.py:20`) provides:

- one reusable CoolProp `AbstractState` per pure fluid or composition;
- direct PT, PH, PS, PQ, and TQ flashes for pure methane;
- outer temperature roots for mixture PH and PS states;
- bubble/dew states and bubble pressure;
- conversion from molar vapor fraction to mass quality;
- a documented component-mixing fallback when CoolProp has no native mixture transport.

Mixture interaction parameters are not invented. Unsupported compositions fail explicitly.

### 3.2 LNG feed line

`solve_lng_feed_line` (`src/fuelnozzle/feed.py:63`) implements a uniform, segmented line:

- single-phase Darcy-Weisbach pressure loss using ChEDL `friction_factor`;
- two-phase friction using ChEDL `Kim_Mudawar` or `Chisholm`;
- static head and distributed minor loss;
- imposed W/m heat leak or an `ht` cylindrical wall/insulation calculation;
- enthalpy marching and detection of the first equilibrium two-phase station;
- comparison of calculated outlet pressure with the available pressure budget.

### 3.3 LNG Tier 0

`screen_lng_flash` (`src/fuelnozzle/lng.py:151`) returns:

- pump and nozzle-inlet states;
- temperature and pressure subcooling margins;
- equilibrium flash fraction at P3;
- equilibrium flash-onset pressure;
- predicted onset category: none, upstream, internal, or exit;
- warnings for pressure and phase risks.

### 3.4 LNG Tier 1

`solve_lng_equilibrium_flow` (`src/fuelnozzle/lng.py:256`) returns:

- the single-phase incompressible (SPI) mass-flux bound;
- an isentropic homogeneous-equilibrium (HEM) pressure path;
- maximum mass flux and corresponding critical pressure;
- an equilibrium choking decision;
- required total area and per-orifice diameter;
- predicted mass flow when physical geometry is supplied.

### 3.5 LNG Tier 2

`solve_lng_relaxation_flow` (`src/fuelnozzle/lng.py:363`) returns:

- delayed nucleation pressure;
- finite-rate movement of vapor quality toward equilibrium quality;
- residence time along the short hole;
- actual flash onset and exit quality;
- mass flux bounded between local HEM and SPI estimates;
- required area/diameter and choking classification.

The delay and relaxation time are calibration parameters, not universal LNG constants.

### 3.6 LNG Tier 3

`solve_lng_flash_spray` (`src/fuelnozzle/spray.py:102`) returns:

- saturation-to-chamber pressure ratio;
- superheat and Jakob number when available;
- mechanical, external-flash, transitional, fully-flashing, or upstream-two-phase class;
- optional calibrated SMD and cone-angle ranges;
- a CFD boundary record.

Without a `FlashSprayCalibration`, SMD and cone angle are deliberately `None`.

### 3.7 Jet-A pressure swirl

`solve_jet_a_pressure_swirl` (`src/fuelnozzle/jet_a.py:143`) returns:

- required exit area and diameter;
- predicted mass flow and effective discharge coefficient;
- inlet-port, axial, and tangential velocities;
- air-core radius and film thickness;
- nominal full cone angle;
- a centrifugal-pressure cavitation screen;
- Reynolds, Weber, and Ohnesorge numbers;
- optional calibrated SMD range.

Jet-A properties come from a measured table or declared source. Jet-A is not represented as a
native CoolProp fluid.

### 3.8 Flight-envelope study

`run_nozzle_study` (`src/fuelnozzle/study.py:75`) evaluates both fuel circuits at every operating
point. It reuses one LNG property provider across the study, integrates fuel mass when all phase
durations are present, compares those masses with externally supplied mission totals, and
preserves warning provenance.

## 4. Existing Verification Tests

The current 21 regression tests cover:

| Area | Existing checks |
|---|---|
| CoolProp | Pure methane PH/PS round trips; mixture PH/PS roots; quality conversion; composition normalization; transport fallback |
| Feed line | Imposed heat-leak energy increase; positive ChEDL loss; insulation-derived heat ingress |
| Tier 0 | No-flash case; internal equilibrium crossing; pressure-budget deficit |
| Tier 1 | Positive SPI/HEM fluxes; geometry sizing; path length |
| Tier 2 | Fast internal flashing; delayed external flashing; no-flash limit; HEM/SPI bounds |
| Tier 3 | Calibration gating; calibrated ranges; cold mechanical-breakup case |
| Jet-A | Sized-flow closure; hollow-cone result; property interpolation/extrapolation; SMD gating |
| Study | Both circuits; stage-mass integration; mission mismatch warning; compositional LNG end to end |

These are useful regression and behavior tests. They do not yet constitute quantitative
experimental validation.

## 5. Deviations from the Original Plan

| Planned direction | Current implementation | Consequence / next action |
|---|---|---|
| Full Abramovich/Giffen-Muraszewski pressure-swirl internal solution | Simplified pressure-energy and angular-momentum closure | Verify algebra with Lacava; later add a named literature model and compare air-core/film predictions |
| Literature pressure-swirl SMD tiers | Calibration coefficient multiplying a capillary sheet scale | Keep outputs labeled calibrated; add literature correlation adapters with validity metadata |
| HEM critical nozzle model | Isentropic thermodynamic path sampled over prescribed pressure stations | Good equilibrium bound; does not solve area/friction-coupled nozzle momentum |
| Homogeneous relaxation model | Quality relaxation imposed on the HEM velocity/path | Screening model only; not a fully coupled HRM conservation solver |
| Internal pressure profile solved from nozzle equations | Pressure stations are linear between inlet and P3 and mapped linearly to axial position | Flash-onset location is approximate and must not be over-interpreted |
| Two-phase choking from a sonic/eigenvalue condition | Choking inferred from an upstream maximum in sampled mass flux | Add convergence tests and compare critical pressure/mass flux against Hammer and De Lorenzo data |
| Detailed LNG cavitation and nucleation physics | Bubble-pressure crossing plus user-calibrated pressure delay | Cannot independently predict nucleation delay without data |
| Universal flash-spray prediction | Threshold classifier plus optional calibration scaling | Correctly prevents unsupported deterministic SMD claims |
| Native LNG mixture transport | CoolProp native where available, otherwise component mixing | Validate viscosity/conductivity; uncertainty is not yet propagated |
| Mission repository integration | Generic stage inputs and optional total-mass check | No direct adapter to `lng_aviation`; P3/T3 remain user inputs as intended |
| Correlation registry with machine-enforced validity | Warnings and local selection logic only | Add structured correlation provenance/validity objects |
| Automated uncertainty propagation | Fixed warning bands and calibration uncertainty | Add Monte Carlo or polynomial-chaos propagation after data schemas stabilize |
| CFD export files | In-memory `CFDSprayBoundary` record | Add OpenFOAM/CSV writers after boundary convention is selected |
| CLI/config study runner | Python public API and example script | Add validated YAML/CSV schemas if batch use is required |

## 6. Planned Paper and Dataset Set

### 6.1 Thermodynamic verification

1. **Setzmann, U.; Wagner, W. (1991).** “A New Equation of State and Tables of
   Thermodynamic Properties for Methane Covering the Range from the Melting Line to
   625 K at Pressures up to 1000 MPa.” *Journal of Physical and Chemical Reference
   Data*, 20(6), 1061-1155. DOI: `10.1063/1.555898`.

   Use: methane saturation and caloric-property reference values. This is the basis of the
   CoolProp methane EOS; fixed regression points should still be recorded.

2. **Bell, I. H. et al. (2014).** “Pure and Pseudo-pure Fluid Thermophysical Property
   Evaluation and the Open-Source Thermophysical Property Library CoolProp.”
   *Industrial & Engineering Chemistry Research*, 53. DOI: `10.1021/ie4033999`.

   Use: property-library provenance and independent API/version checks.

### 6.2 Mass flow, choking, and non-equilibrium model form

3. **Hammer, M.; Deng, H.; Austegard, A.; Log, A. M.; Munkejord, S. T. (2022).**
   “Experiments and modelling of choked flow of CO2 in orifices and nozzles.”
   *International Journal of Multiphase Flow*, 156, 104201.
   DOI: `10.1016/j.ijmultiphaseflow.2022.104201`.

   Use: primary Tier 1/Tier 2 model-form validation. The paper reports that HEM
   underpredicts mass flux and compares orifice/nozzle behavior.

   Transcribe: restriction geometry, upstream state, back pressure, measured mass flux,
   uncertainty, and model curves for one complete orifice series and one nozzle series.

4. **De Lorenzo, M.; Lafon, Ph.; Seynhaeve, J.-M.; Bartosiewicz, Y. (2017).**
   “Benchmark of Delayed Equilibrium Model (DEM) and classic two-phase critical flow
   models against experimental data.” *International Journal of Multiphase Flow*, 92,
   112-130. DOI: `10.1016/j.ijmultiphaseflow.2017.03.004`.

   Use: HEM/frozen/delayed-equilibrium ordering, critical mass flux, and relaxation trends.

   Transcribe: selected Super Moby Dick case geometry/conditions, experimental mass flux,
   model predictions, fitted DEM parameters, and uncertainties.

5. **Kim, Y.; O'Neal, D. L. (1995).** “A comparison of critical flow models for
   estimating two-phase flow of HCFC22 and HFC134a through short-tube orifices.”
   *International Journal of Refrigeration*. DOI: `10.1016/0140-7007(95)93785-I`.

   Use: short-residence metastability and model comparison for orifice-like restrictions.

   Transcribe: individual HFC-134a geometry and inlet state, measured critical flow, and
   the model-comparison figures.

6. **Payne, W. V.; O'Neal, D. L.** “A Mass Flowrate Correlation for Refrigerants
   and Refrigerant Mixtures Flowing Through Short Tubes.” NIST-hosted report.

   Use: supporting geometry/range information and aggregate uncertainty. It includes
   sharp-edged diameters about 1.09-1.94 mm, lengths 9.5-25.4 mm, single- and
   two-phase inlet conditions, and over 1200 points. It is model-form evidence, not LNG
   validation.

### 6.3 Cryogenic flash regime and spray validation

7. **Rees, A.; Oschwald, M. (2022).** “Evolution of flash boiling liquid nitrogen
   sprays by means of high-speed shadowgraphy.” DaRUS V1.
   DOI: `10.18419/DARUS-2527`.

   Use: Tier 3 regime and spray-angle validation over 22 LN2 conditions. The repository
   includes an injection-condition CSV and shadowgraph images, so manual transcription is
   not initially required.

8. **Rees, A.; Salzmann, H.; Sender, J.; Oschwald, M. (2020).** “About the
   Morphology of Flash Boiling Liquid Nitrogen Sprays.” *Atomization and Sprays*,
   30(10), 713-740. DOI: `10.1615/AtomizSpr.2020035265`.

   Use: breakup-regime terminology and quantitative spray-angle interpretation associated
   with DaRUS-2527.

9. **Rees, A.; Araneo, L.; Salzmann, H.; Lamanna, G.; Sender, J.; Oschwald, M.
   (2020).** “Droplet velocity and diameter distributions in flash boiling liquid
   nitrogen jets by means of phase Doppler diagnostics.” *Experiments in Fluids*,
   61, article 182. DOI: `10.1007/s00348-020-03020-7`.

   Associated raw data: Rees, A.; Oschwald, M., DaRUS IN-1,
   DOI: `10.18419/DARUS-2076`.

   Use: local D10 and velocity-field validation and evidence of two droplet populations.
   IN-1 provides raw CSV/XLSX files. Published conditions include approximately 89.7 K,
   440 kPa injection pressure, 7.3 kPa chamber pressure, 1.00 mm diameter, and L/D=2.9.

10. **Gaertner, J. W.; Kronenburg, A.; Rees, A.; Sender, J.; Oschwald, M.;
    Lamanna, G. (2020).** “Numerical and experimental analysis of flashing cryogenic
    nitrogen.” *International Journal of Multiphase Flow*, 130, 103360.
    DOI: `10.1016/j.ijmultiphaseflow.2020.103360`.

    Use: qualitative/quantitative comparison of a numerical flashing model with the DLR
    cryogenic experiment.

    Important correction: DOI `10.1016/j.ijmultiphaseflow.2020.103275` is a different
    automotive in-nozzle simulation and must not be cited as the DLR LN2 paper.

### 6.4 Pressure-swirl validation

11. **Lacava, P. T.; Bastos-Netto, D.; Pimenta, A. P. (2004).** “Design Procedure
    and Experimental Evaluation of Pressure-Swirl Atomizers.” 24th ICAS.

    Use: reproducible hydraulic and spray benchmark for one four-port water atomizer.
    It is an algebra/model-form benchmark, not Jet-A property validation.

    Transcribe:

    - Tables 1, 2, and 4;
    - Figure 4(a), mass flow versus pressure drop;
    - Figure 4(b), discharge coefficient versus pressure drop;
    - Figure 6, spray semi-angle versus pressure drop;
    - Figure 7, experimental/theoretical SMD versus pressure drop;
    - optionally Figure 8, droplet-volume fractions.

    At the 4 atm design point the paper reports 6.13 g/s, Cd=0.2728,
    semi-angle=34.5 degrees, and SMD=80.83 micrometers.

12. **Dodge, L. G.; Biaglow, J. A. (1986).** “Effect of elevated temperature and
    pressure on sprays from simplex swirl atomizers.” *Journal of Engineering for
    Gas Turbines and Power*, 108; ASME 85-GT-58. NASA NTRS ID `19860037997`.

    Use: direct Jet-A SMD, Rosin-Rammler width, axial evolution, ambient density,
    and air-temperature effects. The NTRS record is metadata-only.

    Transcribe from a full copy: atomizer geometry/flow number, Jet-A temperature and
    properties, fuel pressure drop and flow, ambient pressure/temperature/density,
    measurement coordinates, SMD, distribution parameter, and uncertainty.

### 6.5 Sources deliberately excluded from the core quantitative set

- NASA “Liquid Methane/Oxygen Injector Study for Mars Ascent Engines,” NTRS
  `20000013281`: available record is a hot-fire abstract without reconstructable
  injector mass-flow data.
- NASA discharge-coefficient uncertainty report `20190030455`: useful measurement
  methodology, but the identified tests are water cold-flow and not direct LCH4 flashing.
- Automotive flash-boiling studies: useful later for breakup model comparisons, but not
  first-line cryogenic validation.

## 7. Transcription Contracts

Preserve published units and mark gauge versus absolute pressure. Never convert a half-angle
into a full cone angle without retaining the original field.

### 7.1 Restriction/choking row

```text
source, case_id, fluid, geometry_id, diameter, length, L_over_D, inlet_shape,
P0_abs, T0, h0_or_x0, Pback_abs, mdot_or_mass_flux, measurement_uncertainty,
critical_pressure, choked, model_name_if_curve, figure_or_table, notes
```

### 7.2 Pressure-swirl row

```text
source, case_id, fluid, deltaP, deltaP_basis, mdot, Cd,
half_angle_deg, full_angle_deg, D10, D32_or_SMD,
measurement_axial_position, measurement_radial_position,
uncertainty, figure_or_table, notes
```

### 7.3 Provenance rules

- One CSV file per publication/dataset.
- Include a README naming the person and date of transcription.
- Store the page, table, figure, panel, curve color/marker, and digitization tool.
- Keep raw transcribed values immutable; perform unit conversion in a derived file or loader.
- Record uncertainty as published. If absent, write `not_reported`, not zero.
- Mark digitized values separately from table values.

## 8. Verification and Validation Implementation Plan

### Phase A: software and numerical verification

1. Add closed-form SPI mass-flux tests over a parameter grid.
2. Add pressure-grid convergence for Tier 1 using 45, 90, 180, 360, and 720 stations.
3. Require critical mass flux and pressure to converge within selected tolerances.
4. Add Tier 2 axial-grid/time-step convergence.
5. Test limiting behavior:
   - zero superheat gives no flash;
   - relaxation time approaching zero approaches equilibrium quality;
   - very large relaxation time approaches metastable/frozen behavior;
   - zero pressure delay starts relaxation at bubble crossing;
   - increasing back pressure removes the choked plateau;
   - mass flux remains within the implemented HEM/SPI bracket.
6. Add conservation diagnostics for feed-line pressure and enthalpy marching.
7. Add independent CoolProp reference tables and version-pinned regressions.

### Phase B: critical-flow model-form validation

1. Implement a generic validation-case schema decoupled from `OperatingPoint` and methane.
2. Generalize the property-provider protocol so CO2, R134a, and water can exercise the same
   Tier 1 numerical kernel.
3. Load one Hammer orifice series and one nozzle series.
4. Compare HEM mass flux and critical pressure using bias, MAE, MAPE, RMSE, and uncertainty-
   normalized residuals.
5. Calibrate Tier 2 only on a declared training subset.
6. Hold out at least one geometry and one inlet-state series.
7. Repeat with selected De Lorenzo/Super Moby Dick and Kim-O'Neal cases.

Initial acceptance targets, subject to published uncertainty:

- numerical repeatability below 0.1% for fixed settings;
- grid-converged Tier 1 mass flux within 0.5%;
- HEM results reproduce the source paper's HEM curve within digitization/property tolerance;
- Tier 2 improves held-out mass-flux error over HEM without violating physical bounds;
- choking/non-choking classification is correct for at least 95% of clearly classified cases.

### Phase C: cryogenic regime and spray validation

1. Download DaRUS-2527 directly and compute all dimensionless quantities using CoolProp N2.
2. Define image-derived spray angle with a reproducible threshold/edge method.
3. Compare Tier 3 classes and predicted angle trends over all 22 cases.
4. Use DaRUS-2076 raw D10/U/V fields as a separate validation set.
5. Do not fit and validate to the same spatial points.
6. Compare the CFD boundary velocity and scale with IN-1 injection conditions.

Initial targets:

- correct broad regime family for at least 18 of 22 morphology cases;
- spray full angle within experimental/digitization uncertainty plus 10 degrees;
- velocity-scale error reported with the published +/-3 m/s injection uncertainty;
- no claim that calibrated LN2 D10 is LNG D32.

### Phase D: pressure-swirl verification and calibration

1. Add a Lacava fixture containing the full 2-6 atm curves.
2. Verify mass-flow/Cd/cone-angle trends and reproduce the 4 atm table.
3. Split pressure points into calibration and holdout sets for SMD.
4. Add the Dodge/Biaglow Jet-A environmental data as a separate validation set.
5. Replace or supplement the present capillary calibration with named literature adapters.

Initial targets:

- mass flow within 5% and Cd within 0.02 for the Lacava holdout points;
- full cone angle within 5 degrees after consistently converting the published semi-angle;
- SMD uncertainty band covers held-out Lacava and Jet-A points at the declared coverage rate;
- ambient-density direction and axial SMD trend agree with Dodge/Biaglow.

### Phase E: LNG-specific validation

1. Seek direct LCH4 or representative LNG data with composition, inlet state, geometry,
   back-pressure sweep, mass flow, and uncertainty.
2. Freeze Tier 2 parameters before evaluating direct LNG holdout cases.
3. Validate pure methane first, then composition sensitivity.
4. Validate the mixture transport fallback independently or replace it with measured data.
5. Keep all surrogate-fluid conclusions separate from LNG conclusions.

## 9. Highest-Priority Data Request

The first useful transcription batch is:

1. Hammer et al.: full geometry table plus one complete CO2 orifice back-pressure/state series,
   including measured uncertainty and all model curves.
2. Lacava et al.: Figures 4, 6, and 7 at every plotted pressure.
3. Dodge and Biaglow: one complete Jet-A ambient-pressure series once the full paper is available.

This sequence enables quantitative choking validation, immediate pressure-swirl regression tests,
and a direct Jet-A spray check without mixing calibration and validation evidence.

## 10. Documentation and Reproducibility Deliverables

Planned repository additions during V&V implementation:

```text
validation/
  README.md
  schemas/
  raw/<source-key>/
  derived/<source-key>/
  cases/
  reports/
tests/verification/
tests/validation/
```

Each validation report should record package version, CoolProp version, source data checksum,
model settings, calibrated parameters, train/holdout split, metrics, plots, and warning output.
