# Modeling Basis

## Scope

`fuelnozzle` is a steady, reduced-order screening tool. It answers hydraulic,
thermal, phase-regime, and preliminary spray questions across a supplied flight
envelope. It does not resolve the internal three-dimensional vortex, bubble
nucleation field, primary breakup, combustor reaction, or thermoacoustic
response.

## CoolProp Integration

`CoolPropLNGProvider` owns one reusable CoolProp `AbstractState` for each LNG
composition. Reusing the state avoids repeated backend construction and keeps
all LNG tiers on one equation of state and reference state.

Pure methane uses direct CoolProp flashes:

- PT for user and feed-line states;
- PH for adiabatic pressure changes and heated-line marching;
- PS for the isentropic HEM path;
- PQ/TQ for saturation properties.

CoolProp mixtures support PT, PQ, and TQ inputs but not general PH/PS input
pairs. The provider therefore solves an outer, bracketed temperature root around
repeated PT flashes:

$$
f_h(T)=h(P,T,\mathbf z)-h^*=0,
$$

or

$$
f_s(T)=s(P,T,\mathbf z)-s^*=0.
$$

The previous pressure station's temperature supplies a continuation bracket;
a broad scan is used only if continuation fails. CoolProp reports mixture vapor
fraction on a molar basis, so the provider converts it to mass quality before
passing it to hydraulic models.

When CoolProp has no native transport model for an LNG mixture, the provider
uses reusable pure-component states at the local pressure and temperature. It
uses logarithmic mole-fraction mixing for viscosity, linear mole-fraction
mixing for thermal conductivity, and a liquid-component surface-tension mix.
Every affected state is labeled `component_mixing_fallback`, and feed-line
results emit a warning. If a dissolved component has no pure forced-liquid PT
solution, its same-temperature saturated-liquid transport is used. If that is
also unavailable, available component weights are renormalized and the omitted
component is named in the warning. These estimates require validation against
mixture data.

CoolProp mixture interaction parameters are not automatically estimated. A
composition that CoolProp cannot initialize fails explicitly. Mixture runs are
slower than pure methane because every PH/PS state requires an outer root.

## Pressure And Thermal Boundary

The user supplies $P_3$, $T_3$, pump outlet pressure and fuel temperature,
nozzle pressure drop, and fuel mass flow at each operating point. The required
nozzle inlet pressure is $P_3+\Delta P_n$. Pump pressure above that value is
available for feed-line loss.

When `LNGFeedLine` is supplied, pressure and enthalpy are marched from the pump
to the injector. Single-phase pressure loss uses Darcy-Weisbach with ChEDL's
friction factor. Two-phase friction uses ChEDL `Kim_Mudawar` within its published
diameter range and `Chisholm` for larger lines. Static head is added separately.

Heat ingress can be supplied as measured W/m. Otherwise, the model uses ChEDL
`ht` internal convection and cylindrical wall/insulation conduction. The line
outlet enthalpy, rather than pump temperature, becomes the nozzle inlet thermal
condition.

## LNG Tier 0

The principal phase margins are

$$
\Delta T_{sub}=T_{bubble}(P,\mathbf z)-T
$$

and

$$
\Delta P_{sub}=P-P_{bubble}(T,\mathbf z).
$$

The equilibrium flash fraction at $P_3$ comes from an isenthalpic PH state.
Assuming a linear internal pressure profile, the bubble-pressure crossing gives
an equilibrium onset location. This is a risk screen, not proof that bubbles
nucleate there.

## LNG Tier 1

The single-phase incompressible bound is

$$
G_{SPI}=\sqrt{2\rho_0\Delta P_n}.
$$

The HEM path is isentropic. At each pressure,

$$
u(P)=\sqrt{2[h_0-h(P,s_0)]}, \qquad G(P)=\rho(P,s_0)u(P).
$$

The critical mass flux and pressure are the maximum of $G(P)$ between nozzle
inlet and $P_3$. If that maximum occurs upstream of the exit, the equilibrium
model reports choking. Required geometric area is

$$
A=\frac{\dot m}{C_dG}.
$$

HEM assumes instantaneous thermal and mechanical equilibrium. For a short
cryogenic orifice it normally predicts phase change too early and is treated as
one side of a model bracket, not truth.

## LNG Tier 2

Tier 2 uses a homogeneous relaxation relation:

$$
\frac{Dx}{Dt}=\frac{x_{eq}-x}{\tau}.
$$

Relaxation begins only below

$$
P_{nuc}=P_{bubble}-\Delta P_{delay}.
$$

The pressure delay and relaxation time are explicit calibration parameters.
The path is integrated with local residence time. The resulting mass flux is
bounded at each station by the HEM and SPI limits. Tier 2 determines whether
actual flashing is upstream, internal, at the exit, external, or absent.

## LNG Tier 3

Tier 3 reports saturation-pressure ratio, equilibrium and actual vapor
fractions, superheat, Jakob number when available, regime, and a CFD boundary
record. The regime thresholds in `FlashSpraySettings` are visible assumptions.

There is no default universal LNG SMD equation. Without a
`FlashSprayCalibration`, SMD and cone-angle fields are `None` and a warning is
returned. A calibration supplies reference values, scaling exponents, and
uncertainty. It should be based on data with comparable composition,
temperature, pressure ratio, $L/D$, inlet edge, roughness, and ambient gas.

## Jet-A Pressure Swirl

Jet-A is represented by measured property tables or declared values, never as a
native CoolProp fluid. Exit area is sized with

$$
\dot m=C_d A_o\sqrt{2\rho\Delta P_n}.
$$

Tangential port momentum is propagated to the exit with an explicit angular
momentum efficiency and constrained by pressure energy and continuity. The
liquid annulus gives air-core radius and film thickness. Axial and tangential
velocities give the nominal cone angle. Centrifugal pressure depression is used
for a cavitation screen when Jet-A vapor pressure is supplied.

The solver reports Reynolds, Weber, and Ohnesorge numbers. As for LNG, the
Jet-A SMD field remains empty unless a hardware-specific coefficient is supplied.

## Validation Priorities

1. CoolProp states against NIST/REFPROP values for the selected LNG composition.
2. Isothermal line pressure drop and imposed heat-load energy balance.
3. Nozzle mass flow versus inlet state and back pressure, including the choked
   plateau.
4. Flash-onset location from transparent-nozzle or imaging data.
5. Jet-A flow number, air-core size, cone angle, and film thickness.
6. LNG and Jet-A droplet distributions only after hardware calibration.

Every correlation or fitted parameter used for a design decision should retain
its source, validity range, and uncertainty outside this numerical core.