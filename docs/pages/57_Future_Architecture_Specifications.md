@page 57_Future_Architecture_Specifications Future Architecture Specifications

@anchor _Future_Architecture_Specifications

This page indexes architectural work that is proposed rather than built, and
records the dependency order between the pieces. It prevents design decisions
from being lost between independently reviewed branches.

Being listed here is not a statement that a YAML or file format is accepted. The
status column is authoritative: a specification is proposed until its own page
says otherwise.

@tableofcontents

@section p57_status_sec 1. Status and Dependency Order

| Specification | Status | Dependency |
| --- | --- | --- |
| @ref 56_Field_Identity_and_Layout_Catalog | implemented | none |
| @ref 58_Field_Statistics | implemented | field catalog |
| @ref 60_Field_Statistics_Planned_Extensions | proposed | field statistics |
| @ref 59_Function_Identity_and_Observability_Specification | deferred, benchmark-gated | none |

The typed field catalog established Eulerian and particle identities, layout
metadata, and non-owning views over existing vectors without changing vector
ownership. That is the correctness foundation statistics rely on for shifted,
face-centered, and component-staggered fields, which is why it had to land first.

Field statistics build directly on it and are implemented; what remains proposed
for that subsystem is indexed at
**@subpage 60_Field_Statistics_Planned_Extensions**, which carries its own
dependency order.

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

Field statistics are implemented. Windows accumulate numerically stable centered
moments online, ride in the committed checkpoint bundle, resume on continuation,
and are derived into Reynolds stresses, RMS, turbulent kinetic energy, and fluxes
in post-processing. The contract is at **@subpage 58_Field_Statistics**.

The direction the remaining work follows is a shared, field-layout-aware reducer:
the same traversal serving spatial targets, reduced statistics, and eventually the
rolling physical-solution monitor, instead of one implementation per consumer.
That replacement is permitted only where tests demonstrate identical results.
Existing logger formats and PETSc monitors are preserved throughout.

The proposed extensions and their dependency order are at
**@subpage 60_Field_Statistics_Planned_Extensions**.

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
