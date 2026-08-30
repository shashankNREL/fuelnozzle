"""Chemical reactor network models for dual-fuel Jet-A and LNG combustors.

This subpackage extends the nozzle models in :mod:`fuelnozzle` past the injector
exit and into the combustor. It consumes, rather than duplicates, the existing
solvers: :func:`fuelnozzle.solve_lng_flash_spray` and
:func:`fuelnozzle.solve_jet_a_pressure_swirl` supply the spray boundary that the
reactor network starts from.

The design intent is recorded in ``docs/CRN_PLAN.md`` and the build record in
``docs/CRN_IMPLEMENTATION_LOG.md``. Equations and derivations are in
``docs/crn_technical_reference.tex``.

Scope, per the approved plan:

- One fuel is active per operating point. Jet-A and LNG are never co-fired.
- Jet-A and LNG have separate injector hardware; only the combustor is shared.
- Steady state, constant combustor pressure, one sector.
- Quantities of interest: atomization, fuel-air mixing, temperature, and NOx.
  CO, UHC, and soot are uncalibrated diagnostics and carry no accuracy claim.

Planned module layout. Modules appear here as they are implemented rather than
as empty stubs, so an importable name is always a working one::

    chemistry.py        mechanism registry, fuel definitions, equivalence ratio
    streams.py          air split bookkeeping and inlet states
    droplets.py         TAB breakup, evaporation, droplet heating
    spray_source.py     bridge from the nozzle results to droplet classes
    reactors.py         Cantera reactor wrappers
    network.py          reactor graph, mass balance, steady solve
    templates.py        RQL, LPP, LDI, and staged architectures
    autoignition.py     ignition delay, residence time, margin, flashback
    thermal.py          LNG heat-sink budget and feasible thermal window
    emissions.py        emission indices and corrected concentrations
    combustor_study.py  orchestration across operating points
    design.py           design variables and bounds
    evaluate.py         one design vector to a full multi-point result
    objectives.py       objectives and feasibility-first constraints
    optimize.py         sweeps, sensitivity, and Pareto search
    cfd_ingest.py       clustering and mass correction from CFD (deferred)
    calibrate.py        network calibration against CFD (deferred)
"""

from __future__ import annotations

__all__: list[str] = []
