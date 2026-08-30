# CRN Equation and Assumption Register

Status: implementation register; validation evidence is not implied by a populated row.

Every design-driving relation must identify its implementation, applicability, uncertainty,
verification oracle, and validation evidence before it can support design acceptance.
`missing` means the relation is deliberately unavailable rather than silently assumed.

| Relation or assumption | Implementation | Source / basis | Applicability and uncertainty | Verification | Validation evidence | Acceptance status |
|---|---|---|---|---|---|---|
| Fuel pressure budget | `operating.fuel_pressure_budget` | Bernoulli pressure accounting | Steady, prescribed losses | Hand closure and deficit limits | Feed/nozzle rig | screening |
| Combustor pressure stations | `operating.resolve_pressure_stations` | Prescribed fractional liner pressure loss | Treats P3 as both compressor discharge and dome pressure; active pump must supply chamber plus nozzle drop | Pressure-path closure and deficit rejection | Engine cycle deck and pressure survey | screening |
| Sector scaling | `crn.hardware.SectorDefinition` | Geometric cup-count ratio | Requires nominally identical cups and engine-total source flows | Forward/reverse scaling identity | Annular/sector flow survey | arithmetic only |
| Liner air admission | `crn.hardware.AirAdmission.mass_flow_kg_s` | Quasi-steady, one-dimensional, isentropic ideal-gas orifice relation with choked limit | Constant gamma/R; one effective discharge coefficient and common downstream pressure; does not predict jet penetration or solve pressure loss | Subcritical/choked analytical limits and area-ratio tests | Cold-flow pressure/flow survey | unvalidated |
| Fixed-hardware air split | `crn.hardware.SharedLinerGeometry.area_derived_split` | Normalized effective-area flows | Uses prescribed P3/P4; absolute areas do not yet close total combustor air flow or determine pressure loss | Fraction sum and effective-area ratio | Sector effective-area calibration | unvalidated |
| LNG thermodynamic state and flash | `properties.CoolPropLNGProvider` | CoolProp EOS | Backend and mixture dependent | PT/PH/PS round trips | REFPROP/GERG/NIST holdout | screening |
| LNG equilibrium critical flow | `lng.solve_lng_equilibrium_flow` | Isentropic HEM | Homogeneous equilibrium bound | SPI limit and grid refinement | Critical-flow datasets | screening |
| LNG finite-rate flashing | `lng.solve_lng_relaxation_flow` | Empirical relaxation | Requires geometry-matched calibration | Fast/slow limits and grid refinement | Direct LCH4/LNG flow test | unvalidated |
| Jet-A pressure-swirl hydraulics | `jet_a.solve_jet_a_pressure_swirl` | Reduced-order Bernoulli/angular momentum | Hardware-specific discharge and swirl closures | Mass/energy identities | Flow, air-core, cone and cavitation rig | unvalidated |
| Spray-size distributions | `crn.spray_source` | Rosin-Rammler | Requires calibrated D32 and distribution exponent | Number/mass/D32 reconstruction | Phase-Doppler or imaging data | unvalidated |
| TAB aerodynamic breakup | `crn.droplets.taylor_analogy_breakup` | Taylor analogy breakup | Aerodynamic breakup only | Damped-oscillator limits | Elevated-pressure breakup data | unvalidated |
| Droplet evaporation/heating | `crn.droplets.integrate_droplet` | Convective mass/heat transfer | Single pseudo-component; no wall interaction | D2 and full enthalpy limits | Jet-A/LNG droplet data | unvalidated |
| Reactor mass balance | `crn.network.close_internal_flows` | Node-wise steady mass conservation | Directed connected graph; fixed measured flows and uncertainties must be declared | Independent graph fixtures, infeasible topology, and nonnegative bounds | not applicable | verified arithmetic |
| PSR residence time | `crn.network.ReactorSolution` | `tau = mass / mass flow`, with `mass = rho * physical volume` | Constant-pressure steady control volume; physical volume is re-imposed during time marching | Cantera mass/volume/inventory identity | Tracer residence distribution | verified arithmetic |
| PFR approximation | `crn.network._expand_plug_flow_reactors` | Series of equal-volume stirred reactors | Segment-converged per case; numerical axial mixing remains | Exact volume/flow topology and 2/4/8 segment refinement | Tracer/species axial profiles | verified implementation, unvalidated physics |
| Quench/dilution mixing | `crn.templates` | Staged PSR approximation | Calibration-specific schedule and volume | Stage refinement | Confined-jet/sector data | unvalidated |
| Wall heat transfer and cooling | `crn.reactors`, `crn.templates._cooling_inlet` | Declared cooling re-entry plus prescribed calibrated heat | Free heat-loss watts are rejected; physical wall/coolant closure remains missing | Destination routing and steady enthalpy balance | Liner/coolant thermal data | calibrated screening only |
| Jet-A chemistry and NOx | `crn.chemistry` | HyChem A2 + nitrogen chemistry | Mechanism-range and reduction dependent | Source-case replay | Held-out flame/rig data | unvalidated |
| LNG chemistry and NOx | `crn.chemistry.MechanismRegistry.with_lng_composition` | GRI-Mech 3.0 baseline with explicit property-name-to-species mapping | Pressure/composition extrapolation limited; unsupported components are rejected rather than dropped | Composition propagation and missing-species tests | High-pressure representative-LNG data | unvalidated |
| Autoignition delay | `crn.autoignition` | Homogeneous kinetics | Mechanism, marker, P/T/phi dependent | Direct integration and interpolation | Shock-tube/RCM holdout | screening |
| Flashback | `crn.autoignition.flashback_screen` | Velocity/flame-speed comparison | Requires passage flow and validated flame speed | Limit tests | Premixer rig | unavailable by default |
| Lean blowout | missing | Damkohler/extinction plus flameholding | Steady CRN alone is insufficient | Continuation/hysteresis | Pressure-matched sector rig | missing |
| NOx emission index | `crn.emissions` | Exhaust mass per fuel mass | Mechanism/network validation required | Hand species accounting | ICAO/sector-rig holdout | unvalidated |
| ICAO LTO aggregation | `crn.emissions.lto_dp_foo` | ICAO time in mode and fuel mass | All four modes required | Hand cycle | ICAO EDB comparison | arithmetic only |
| Exit temperature quality | missing | Mass-weighted exit traverse/pattern factor | Parallel exit representation required | Synthetic traverse | Sector/engine traverse | missing |
| Optimization ranking | `crn.objectives`, `crn.optimize` | Feasibility-first Pareto search | Only validated quantities may drive selection | Synthetic dominance and scale tests | End-to-end holdouts | screening |

The authoritative remediation and V&V requirements are in
`docs/CRN_IMPLEMENTATION_LOG.md`, Sections 6–9. Update this register in the same change as
any equation, correlation, validity range, or evidence status.
