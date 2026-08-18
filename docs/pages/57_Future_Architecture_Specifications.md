@page 57_Future_Architecture_Specifications Future Architecture Specifications

@anchor _Future_Architecture_Specifications

This page records phased architectural direction and is retained for explicit
review and approval.
It prevents design decisions from being lost between independently reviewed
branches. A specification on this page is not a statement that its YAML or file
format is currently accepted by PICurv.

@tableofcontents

@section p57_status_sec 1. Status and Dependency Order

| Specification | Status | Dependency |
| --- | --- | --- |
| @ref 56_Field_Identity_and_Layout_Catalog | implemented on `main` | none |
| @ref 58_Turbulence_Statistics_Pipeline_Specification | Phase 1 implemented, awaiting review and merge; Phase 2 specified but not implemented | field catalog |
| @ref 59_Function_Identity_and_Observability_Specification | deferred, independently benchmarked | none |

The field-catalog phase established typed Eulerian and particle identities,
layout metadata, and existing-vector views without changing vector ownership.
That is the required correctness foundation for statistics on shifted,
face-centered, and component-staggered fields.

@section p57_rules_sec 2. Rules Shared by Future Work

Future architectural work must:

- use a separate branch for each reviewable phase;
- inspect relevant implementation history before changing layout, periodic, I/O,
  restart, or monitoring behavior;
- reuse existing setup, geometry, mask, field, I/O, logger, and postprocessing
  surfaces before adding a new implementation;
- retain the high-level setup and run-loop shape unless a proven requirement
  makes a local orchestration change necessary;
- add regression tests that fail without the intended change;
- update public templates, configuration validation, ingestion documentation,
  runtime documentation, and developer documentation in the same phase; and
- run the full serial, MPI, periodic, restart, postprocessing, ingress, and
  documentation gates before merge when those surfaces are affected.

Strings remain correct at true ingress and presentation boundaries. Runtime
systems resolve them once into the typed identity appropriate to that system;
PICurv does not use one universal ID namespace for fields, functions, boundary
handlers, and postprocessing operations.

@section p57_stats_summary_sec 3. Statistics Direction

The statistics work will add reusable completed-state scheduling and
field-layout-aware reducers shared where useful by scientific statistics and
rolling physical-solution monitoring. It will preserve existing logger formats
and PETSc monitors, accumulate numerically stable centered primitive moments
online, and derive Reynolds stresses, RMS, TKE, turbulent fluxes, and normalized
outputs in postprocessing.

The current Phase 1 foundation keeps `control` as the single generated C-ingress
artifact, moves convergence policy to
`monitor.yml -> solution_monitoring.convergence`, and removes rather than wraps
the legacy `-averaging`/`su*` system. It deliberately exposes no scientific
window YAML before the accumulator, checkpoint, and postprocessor path is
usable. The committed-checkpoint contract and its coordinator, validator, and
restart discovery are implemented, which completes Phase 1.

The authoritative design, YAML placement, checkpoint naming, restart rules, and
implementation phases are in
**@subpage 58_Turbulence_Statistics_Pipeline_Specification**.

@section p57_functions_summary_sec 4. Function Identity Direction

Function logging/profiling identity is deliberately separate. It may replace
repeated function-name scans with compile-time IDs or cached handles only after
benchmarks show that the lookup overhead is material. It does not block the
statistics pipeline. See
**@subpage 59_Function_Identity_and_Observability_Specification**.

@section p57_branch_sec 5. Branch Policy

Specification-only work uses a documentation branch. Implementation begins only
after explicit plan approval, from current `main`, and each completed phase is
merged before the next phase branch is created. A later phase must not rely on
unreviewed changes in another long-lived branch.
