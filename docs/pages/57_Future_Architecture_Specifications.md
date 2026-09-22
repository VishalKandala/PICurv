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
| Immersed boundaries, moving bodies, moving frames (@ref p57_ibm_sec) | planned; switches refused | none |
| Multi-block coupling (@ref p57_multiblock_sec) | planned; `blocks` other than 1 refused | none |
| RANS closures (@ref p57_rans_sec) | planned; `rans` block and `-rans` refused | none |
| Newton-Krylov Jacobian types and modes (@ref p57_nk_jacobian_sec) | planned; refused at validation | Newton-Krylov solver |
| @ref 17_Workflow_Extensibility | proposed extension directions | none |

The typed field catalog established Eulerian and particle identities, layout
metadata, and non-owning views over existing vectors without changing vector
ownership. That is the correctness foundation statistics rely on for shifted,
face-centered, and component-staggered fields, which is why it had to land first.

Field statistics build directly on it and are implemented; what remains proposed
for that subsystem is indexed at
**@subpage 60_Field_Statistics_Planned_Extensions**, which carries its own
dependency order.

Proposed extensions to the conductor workflow - further grid geometries, grid-quality
gates, sweep retry, convergence-based run completion, and data-driven particle
closures - are collected at **@subpage 17_Workflow_Extensibility**.

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

The direction the remaining work follows is that averaging commutes with linear
operations and not with anything else. Spatial reduction over accumulated
pointwise state is exact after the fact, so profiles, regions, and bins are
post-processing operations; only what the reduction cannot reach — the moment
order, the choice of products, statistics of interpolated or nonlinear
quantities, and conditional sampling — has to be resolved while the solver runs.
The same reduction traversal could eventually serve the rolling
physical-solution monitor, permitted only where tests demonstrate identical
results. Existing logger formats and PETSc monitors are preserved throughout.

The proposed extensions and their dependency order are at
**@subpage 60_Field_Statistics_Planned_Extensions**.

@section p57_functions_summary_sec 4. Function Identity Direction

Function logging/profiling identity is deliberately separate. It may replace
repeated function-name scans with compile-time IDs or cached handles only after
benchmarks show that the lookup overhead is material. It does not block the
statistics pipeline. See
**@subpage 59_Function_Identity_and_Observability_Specification**.

@section p57_ibm_sec 5. Immersed Boundaries, Moving Bodies, and Moving Frames

**Status: planned, not implemented.** The legacy CURVIB code this solver descends from
imposed immersed bodies on the curvilinear grid, moved them, and coupled their motion to
the flow. None of that was ported. What reached this tree is the switches and some
dormant plumbing:

- `models.physics.fsi.immersed` and `moving_fsi` (`-imm`, `-fsi`), plus the
  passthrough-only `-rfsi`, `-mframe`, `-rframe`, `-mhv` and `-lv`;
- a Poisson pre-solve branch for immersed cases (solid-aware restriction of `Nvert`
  through the multigrid hierarchy and a check for fully blocked regions), and `Nvert`
  solid markers that the momentum and Poisson kernels honour;
- commented-out calls to `ibm_interpolation_advanced`, a function defined nowhere, in
  both momentum solvers and the flow-solver orchestration; and a commented-out
  moving-frame convection branch in `ComputeRHS`.

Each switch used to be accepted, and each ran a different problem without saying so:
`immersed` reconfigured the Poisson solve around a body that was never loaded,
`moving_fsi` did nothing, and `-mframe`/`-rframe` skipped the convective term
altogether. They are now refused - the YAML switches by `picurv validate`, and the flags
by `CreateSimulationContext` with `PETSC_ERR_SUP` - and
`tests/c/test_setup_lifecycle.c` holds that refusal in place.

**What must exist before any switch is accepted again.** Loading a body surface and
classifying cells against it into `Nvert`; the interpolation that sets velocity at the
immersed-boundary nodes each stage, in every momentum solver that claims support;
an end-to-end check that the dormant Poisson branch conserves mass around a real body;
and, for moving bodies or frames, the motion model and the frame-relative convection
term. The Newton-Krylov solver refuses all of these independently and would need its
own treatment.

**Design owner:** the repository owner. Nothing here is scheduled.

@section p57_multiblock_sec 6. Multi-Block Coupling

**Status: planned, not implemented.** The configuration and setup layers are
multi-block aware: `models.domain.blocks` reaches `-nblk`, PICGRID files carry several
blocks, `boundary_conditions` accepts one face list per block and stages a `bcs` file
for each, and the solvers loop over blocks. What is missing is the coupling. Every
inter-block exchange (`Block_Interface_U`) is commented out - in the explicit and
dual-time momentum solvers and in the initial condition - and nothing couples the
pressure solve across blocks, so a multi-block case was a set of isolated single-block
solves that reported themselves as one domain. No test or shipped example ever ran one.

Validation now refuses `blocks` other than 1, and setup refuses `-nblk` other than 1.

**What must exist before more than one block is accepted.** Interface exchange of the
velocity fields at shared faces in every momentum solver and in the initial
condition; pressure coupling across blocks; particle hand-off between blocks; and a
runtime harness showing that a domain split into blocks reproduces the single-block
solution. The Newton-Krylov solver refuses multiple blocks independently. Field
statistics payloads are already block scoped and follow the same natural ordering
the Eulerian payloads do; the harness that proves multi-block equivalence should
cover them at the same time.

**Design owner:** the repository owner. Nothing here is scheduled.

@section p57_rans_sec 6a. RANS Closures

**Status: planned, not implemented.** PICurv models turbulence with the LES closures at
**@subpage 72_LES_Turbulence_Closure**; there is no Reynolds-averaged path. A `k_omega`
selector existed until 2026-09-22 and was removed: nothing behind it was ever built. Setup
allocated an eddy viscosity but never the `K_Omega` fields, the transport update in
`FlowSolver` was commented out, and the function it called was defined nowhere, so a case
that enabled it copied a null vector and aborted at the end of the first timestep. It was
recorded known-defective on 2026-09-18 and, since no implementation had ever existed,
returned to planned when the dead hooks came out. `src/guide.md` lists exactly what was
removed.

`models.physics.turbulence.rans` is refused at validation, and `-rans` is refused at setup,
so nothing is silently ignored.

**What must exist before a RANS selector returns.** Storage and ghost exchange for the
turbulence variables; a transport equation for each, discretized on the same curvilinear
metrics as momentum, with their production, dissipation and cross-diffusion terms; wall
treatment matched to the closure, including which wall function is admissible with it;
checkpoint and restart of the turbulence state; and a validation case with a reference
profile - a channel at a published `Re_tau` - showing the mean profile and the eddy
viscosity the closure is supposed to produce. A wall-modelled LES path is the nearer
alternative for the same engineering questions, and the wall functions already exist.

**Design owner:** the repository owner. Nothing here is scheduled.

@section p57_nk_jacobian_sec 7. Newton-Krylov Jacobian Types and Modes

**Status: planned, not implemented.** The Newton-Krylov solver's `jacobian` block is a
discriminated configuration so that further constructions can be added beside the one
that exists, `type: finite_difference` with `mode: matrix_free`. Two are designed; both
are refused at validation until they exist.

A second finite-difference mode, `colored_sparse`, would assemble a sparse numerical
Jacobian using coloring:

```yaml
    jacobian:
      type: finite_difference
      finite_difference:
        mode: colored_sparse
```

Frozen-momentum approximations are a different Jacobian type, not storage modes of
finite difference:

```yaml
    jacobian:
      type: frozen_momentum_approximation
      frozen_momentum_approximation:
        structure: diagonal      # or: full_sparse
```

Where each would be added is recorded with the solver's extension points at
@ref p55_precond_sec.

**Design owner:** the repository owner. Nothing here is scheduled.

@section p57_branch_sec 8. Branch Policy

Specification-only work uses a documentation branch. Implementation begins only
after explicit plan approval, from current `main`, and each completed phase is
merged before the next phase branch is created. A later phase must not rely on
unreviewed changes in another long-lived branch.
