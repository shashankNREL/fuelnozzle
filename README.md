# fuelnozzle

Reduced-order design models for a Jet-A pressure-swirl injector and a
liquid-LNG injector. LNG flashing is an output of the calculation, not an input
mode selected by the user.

The package evaluates every user-provided flight-envelope operating point with:

- user-specified $P_3$, $T_3$, fuel mass flows, pump pressures, fuel
  temperatures, and nozzle pressure drops;
- a Jet-A simplex pressure-swirl model;
- an optional heated LNG feed-line model;
- four LNG tiers from phase screening through calibrated spray estimates;
- mission fuel-mass consistency checks when stage durations are provided.

All numerical inputs and outputs use base SI units. Field names include units.

## Environment

The project uses `pixi`; do not install its dependencies with the system Python.

```shell
cd /Users/syellapa/Documents/Research/2026/SAF/LNG/fuelnozzle
pixi install
pixi run test
pixi run lint
pixi run example
pixi run doc
```

CoolProp is installed from PyPI at version 8 or newer because the compositional
LNG implementation depends on CoolProp 8 mixture stability and PT phase-split
behavior. The remaining dependencies are resolved through conda-forge.

## LNG Tiers

1. **Tier 0, phase screen:** subcooling margins, equilibrium flash fraction,
	and predicted saturation-crossing location.
2. **Tier 1, SPI and HEM:** single-phase incompressible and homogeneous
	equilibrium mass-flux bounds, critical pressure, choking, and orifice sizing.
3. **Tier 2, finite-rate relaxation:** metastable pressure delay, phase-change
	relaxation time, actual flash-onset location, and a mass flux bounded by SPI
	and HEM.
4. **Tier 3, spray screen:** mechanical, external-flash, transitional, or fully
	flashing regime plus CFD boundary data. SMD and cone angle are returned only
	when a traceable hardware calibration is supplied.

Each Tier 3 result retains its nested Tier 2, Tier 1, and Tier 0 results, so no
intermediate assumptions or warnings are lost.

## Pressure Contract

For each fuel, the nozzle inlet pressure required by the user inputs is

$$
P_{\mathrm{nozzle,in}} = P_3 + \Delta P_{\mathrm{nozzle}}.
$$

The remaining pump pressure is the allowance for the feed system:

$$
\Delta P_{\mathrm{feed,available}}
=P_{\mathrm{pump,out}}-P_{\mathrm{nozzle,in}}.
$$

The tool never silently changes $P_3$, pump pressure, or nozzle pressure drop to
close this budget. It emits a warning or error when the inputs are inconsistent.
In the current model, `p3_pa` is also the receiving pressure at the nozzle exit.
If combustor liner loss must be represented separately, supply the downstream
pressure appropriate to the nozzle as `p3_pa`.

## Mission Coupling

The separate `lng_aviation` repository can provide mission-level fuel masses
and phase durations. This package does not infer $P_3$, $T_3$, or injector
conditions from that model. Use `MissionFuelMasses` for a consistency check and
provide the station schedule explicitly with `OperatingPoint` objects.

`flow_multiplier` states how many identical engine/nozzle circuits a point's
mass flows represent. It affects mission mass integration, not individual
nozzle sizing.

See `examples/run_study.py` for an executable study.

## Combustor Reactor Network

The `fuelnozzle.crn` subpackage extends the nozzle models past the injector exit into the
combustor, for dual-fuel design studies: **Jet-A at landing and take-off, LNG at cruise**,
one fuel at a time, with **separate injector hardware per fuel** feeding a **shared
liner**. It consumes the nozzle results rather than duplicating them.

It answers three questions:

1. how much air should go to the head end during Jet-A landing and take-off;
2. how lean LNG can run at cruise;
3. what nozzle and passage geometry each fuel path should target.

Run `pixi run example-combustor` for an end-to-end demonstration.

### Scope and claim boundary

Steady state, constant combustor pressure, one sector, one fuel per operating point.
Quantities of interest are atomization, fuel-air mixing, temperature, and NOx. **CO, UHC,
and soot are uncalibrated diagnostics and carry no accuracy claim.**

> The extension is a reduced-order design-space screening tool for comparing combustor
> architectures and nozzle thermal strategies. Absolute NOx is not validated for this
> hardware and must not be used as a certification or compliance estimate. Its intended
> use is *relative* comparison under one stated set of assumptions.

Jet-A and LNG use mechanisms of different provenance and validity, so **cross-fuel
absolute NOx is trend-level only**.

### Mechanisms

See `mech/README.md`. Jet-A needs **two**: `A2NOx_skeletal` for the network and NOx, and
`A2NTCfast_ske` for ignition delay. That split is not stylistic — a high-temperature-only
mechanism overpredicts Jet-A ignition delay by 510x at 700 K, which would declare an
unsafe premixer safe. LNG uses GRI-Mech 3.0, which ships with Cantera.

## Documentation

- `docs/modeling.md`: concise model equations, CoolProp integration details,
	assumptions, and validation limits.
- `docs/V_AND_V_ROADMAP.md`: implementation inventory, deviations from the
	original plan, paper and dataset register, and staged V&V plan.
- `docs/technical_reference.tex`: beginner-oriented derivation of the equations,
	algorithms, source-function line map, limitations, and planned evidence.
- `docs/technical_reference.pdf`: compiled technical reference generated by
	Tectonic with `pixi run doc`.
- `docs/CRN_PLAN.md`: design intent, task list, and risk register for the reactor-network
	extension.
- `docs/CRN_IMPLEMENTATION_LOG.md`: build record, numbered deviations from the plan, open
	items, and cumulative verification results.
- `docs/crn_technical_reference.tex`: beginner-oriented derivation of the reactor-network
	equations, with the same separation of physical law, numerical approximation, and
	empirical calibration.
- `mech/README.md`: kinetic mechanism registry, provenance, and validity caveats.

## Status

This is a reduced-order research and screening tool. It is not a certification,
combustor operability, or detailed injector-CFD method. Validate final hardware
against cryogenic flow and spray tests with comparable fluid composition,
orifice geometry, pressure ratio, and thermal boundary conditions.
