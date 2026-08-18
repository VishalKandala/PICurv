@page 18_Changelog Changelog

@anchor _Changelog

# Changelog

## Unreleased

- Restored future-architecture documentation and added the authoritative
  proposed field-statistics specification, including centered moment state,
  layout/mask rules, future spatial targets, restart behavior, postprocessing
  responsibilities, legacy removal, and phased implementation plan. Also
  restored the separate benchmark-gated function-identity specification.

- Field-statistics foundation (Phase 1 implementation):
  - moved physical-solution convergence configuration to
    `monitor.yml -> solution_monitoring.convergence`, preserving its existing
    flags, every-completed-step behavior, and log format while adding an
    effective enable switch;
  - kept the generated master `control` as the single C-ingress artifact and
    removed the premature observation sidecar, schema/version envelope,
    exposed cadence, and unsupported-window startup gate;
  - kept unfinished scientific field-statistics keys out of active user YAML
    until accumulation, checkpointing, and postprocessing work end to end;
  - removed `SimCtx::averaging`, legacy sum vectors, `su0`/`su1`/`su2`/`sp`
    read/write hooks, the dead averaging call, old templates, and the reserved
    averaged-field postprocessor passthrough; and
  - replaced flat restart files with immutable committed checkpoint bundles
    containing a text manifest, manifest SHA-256 commit marker, geometry/layout
    identity, catalog-selected PETSc-vector payloads, optional particle/model
    state, and the previous contravariant velocity needed for exact BDF2 restart;
  - routed solve output, continuation, restart staging, live post discovery, and
    final/graceful-shutdown output through the bundle contract while retaining
    the existing generic vector and swarm readers/writers;
  - kept solver/convergence/profiling and other persistent logs independent of
    checkpoint state; and
  - updated Python/C regression coverage and ingress/runtime/reference
    documentation. Scientific statistics accumulation remains Phase 2 work.

- Routed the remaining out-of-factory allocations through the designated
  infrastructure. The coarsest-level immersed-boundary `KSKE` workspace was
  allocated and freed on every Poisson solve from inside the solver, with a
  defensive second free in teardown; it is now allocated once by
  `CreateAndInitializeAllVectors` and freed once in teardown, which is safe
  because `FullyBlocked` guards every read with its own flag array and so carries
  no state between solves. The corner-staging workspace was created lazily inside
  the interpolation routine and rebuilt whenever the block size changed; it is now
  two explicitly typed pairs, `CellScalarAtCorner` and `CellVectorAtCorner`,
  created by the same factory.
- Gave the corner-staging workspace typed field identities. Because the vectors
  are now fixed named members with global and local pairs, they satisfy the
  catalog's `offsetof` requirement, so the hand-rolled global-to-local scatter in
  the interpolation path is replaced by `UpdateLocalGhosts`. Corner anatomy
  logging takes the field identity from its caller rather than inferring the
  degree of freedom from a cached vector's block size.

- Added pointwise spatial target resolution for the field-statistics pipeline:
  `include/statistics_target.h` and `src/statistics_target.c` resolve, for one
  field on one block, an iteration domain that excludes PETSc halo storage and
  solver-layout boundary, dummy, and duplicate-periodic indices, which are
  distinct categories. Layout classification comes from the typed field catalog:
  cell-like dimensions span the shifted interior, node-like dimensions carry one
  more entry, and each face family is node-like only in its own direction.
  Periodic directions drop the wrapped duplicate plane, which leaves cell-like
  spans unchanged because their duplicates already sat outside. Component-
  staggered fields are rejected, since their components live on different face
  families and cannot share one pointwise domain. `SpatialTargetPlan` exists with
  only the pointwise identity mapping so later spatial bins extend rather than
  retrofit it. Covered by a new `unit-statistics-target` suite across cell, node,
  and I/J/K-face layouts in nonperiodic, mixed, and fully periodic domains, plus
  a multi-rank case asserting the resolved domain is decomposition independent.

- Added the weighted centered-moment kernels backing the field-statistics
  accumulator contract: `include/statistics_moments.h` and
  `src/statistics_moments.c` implement the stable weighted Welford update, the
  compatible co-moment update, and the weighted parallel-merge formula, plus
  weighted variance, covariance, and Kish effective sample size. The module is
  pure: no configuration, no PETSc vectors, and no window or scheduling
  knowledge. Non-positive sample weights are rejected rather than silently
  corrupting an accumulator. Covered by a new `unit-statistics` suite asserting
  bitwise-zero variance for constant signals, known scalar and two-sample
  results under equal and unequal weights, all six symmetric velocity
  components of a three-sample self-product, exact agreement between a
  self-paired co-moment and the scalar second moment, high-mean/low-fluctuation
  precision where a raw sum-of-squares would cancel, and merge-equals-sequential
  including empty-partition no-ops.

- Added the field-statistics Phase 2 implementation specification, settling the
  contracts page 58 deferred: right-rectangle endpoint quadrature for
  physical-time weighting (chosen because the difference from trapezoidal is an
  endpoint artifact decaying as 1/N, far below the statistical sampling floor,
  while trapezoidal would force a checkpointed field snapshot per accumulated
  field); final-interval clipping at a bounded window's end; and an interval
  convention under which a zero-length interval is not a sample, making
  initial-state handling identical under both weightings and removing the
  `include_initial` control. Also fixes the user contract, cross-field
  covariance spelling, window identity hash inputs, the indexed control-option
  scheme with its ingress-audit extension, the `statistics/` checkpoint
  namespace, the postprocessing contract, and the seven implementation stages.
  Amended page 58 accordingly.

- Field-statistics Phase 1 closeout:
  - added the same-step checkpoint acceptance test, covering both the idempotent
    no-op rewrite of an already-committed step and the rejection of a same-step
    write whose physical time disagrees with the committed bundle;
  - added a multi-rank smoke scenario that restarts one committed bundle at two
    different rank counts, asserting that load-and-store payloads survive the
    decomposition change byte identically and that the resumed run keeps
    marching;
  - documented that particle payloads are deliberately not block scoped in
    version 1 of the checkpoint format, since PICurv carries a single swarm; and
  - corrected stale status prose that still described the checkpoint contract as
    outstanding Phase 1 work; and
  - characterized restart fidelity and recorded it in the maintenance backlog:
    re-applying boundary conditions to loaded fields means a restart resumes from
    a boundary-consistent state that differs from the stored one at the outlet.
    The offset is outlet-local, decays, sits below solver stopping tolerances, and
    is invariant to them, so it cannot be tightened away. Softened the checkpoint
    contract's BDF2 wording accordingly: preserving `Ucont_rm1` preserves the order
    of the restart, not bitwise identity.

- Field identity infrastructure (phase 1):
  - added a 38-entry typed `FieldId`/`FieldDescriptor` catalog with canonical
    names, active aliases, DM family, degree of freedom, shifted/staggered
    layout, synchronization class, availability, capabilities, and existing
    `UserCtx` vector bindings;
  - added non-owning runtime `FieldView` resolution without changing PETSc
    vector allocation or teardown;
  - replaced string dispatch in ghost updates, periodic field synchronization,
    initial interior-field routing, and field-layout diagnostics with catalog
    metadata while preserving their numerical loops;
  - added a separate ten-entry `ParticleFieldId` catalog for DMSwarm component
    count, PETSc data type, registration ownership, initialization, model-update,
    and Eulerian-scatter metadata;
  - converted particle registration, initialization, scalar physics, analytical
    assignment, and Eulerian/particle interpolation/scatter orchestration to
    typed identities;
  - centralized fixed DMSwarm access names and made generic particle gather and
    restart/output casting use PETSc's registered field datatype;
  - retained runtime configuration, function logging/profiling, generic swarm
    I/O, postprocessing recipe dispatch, and all file formats as textual
    boundaries; and
  - expanded setup, periodic, runtime, solver, and logging regression coverage.

- Solver continuation now rejects `run_control.start_step: 0` instead of treating
  a fresh start as an in-place continuation and appending misleading step-zero
  separators to existing logs.

- Picard-Jameson spectral-radius pseudo-time step (semantic change):
  - changed the pseudo-time step from `dtau = pseudo_cfl × dt` to
    `dtau = pseudo_cfl / lambda_max`, where `lambda_max` is the global maximum
    convective spectral radius computed once per physical timestep. This makes
    `pseudo_cfl.*` true dimensionless Courant numbers, independent of `dt`, grid
    size, and flow speed.
  - changed shipped defaults to `pseudo_cfl.initial: 0.5` and
    `pseudo_cfl.maximum: 2.0` (stable 4-stage Jameson range ~0–2.83); existing
    configs with `0.1`/`1.0` still run but now produce physically smaller, safer
    `dtau` values.
  - added a BDF2 lower bound `lambda_max >= COEF_TIME_ACCURACY/dt` so zero-flow
    startup yields a finite `dtau`.
  - added the diagnostic field `mom_last_lambda_max` to the simulation context.
  - changed the per-trial momentum convergence log fields from `Pseudo-cfl`/
    `cfl_after` to `dtau`/`cfl_eff` and `dtau_after`/`cfl_eff_after`; `picurv
    summarize` parsing, JSON payload, console output, and plot series were
    updated to match.

- Initial-condition contract and workflow:
  - replaced user-facing legacy IC selectors with canonical `generated` and `file` modes.
  - added built-in zero, Cartesian constant, streamwise constant, and Poiseuille generators.
  - added single-block file-backed `Ucat` and `Ucont` startup through the existing field reader.
  - added `ic_gen` run/precompute orchestration, defaulting to `generators/ic.gen` with an optional compatible script override.
  - added the repository `generators/ic.gen` expression engine for grid-aware `Ucat` and staggered `Ucont` PETSc vectors.
  - made initial conditions subordinate to analytical, load, and restart Eulerian-state authority.
  - aligned grid, inlet-profile, and IC generators around repository defaults plus optional case-relative or absolute overrides.

- Picard-Jameson momentum-controller hardening:
  - changed pseudo-time trial acceptance and rollback to one global
    transactional decision across blocks and MPI ranks.
  - bounded runtime by counting rejected and accepted trials against
    `max_iterations`, and guaranteed nonfatal exits retain the last accepted
    finite state.
  - unified pseudo-CFL adaptation and carry-over, removed the extra hard-coded
    reduction multiplier and independent low-CFL rebound, and added controller
    bounds validation.
  - added optional `residual_absolute_tol` and `residual_relative_tol`;
    configurations without enabled residual tolerances retain legacy
    update-only convergence semantics.
  - deprecated unused `step_tol` while retaining compatibility ingestion.

- Search robustness observability:
  - added always-on aggregate search instrumentation for particle-enabled runs.
  - added `logs/search_metrics.csv` with timestep-level search, traversal, tie-break, boundary-clamp, bbox-guess, pass-depth, per-step loss, run-local cumulative loss, V2 population/outcome counters, and derived `search_failure_fraction`, `search_work_index`, and `re_search_fraction` signals.
  - added `LOG_SEARCH_METRICS` for compact DEBUG-gated console summaries when explicitly allow-listed.
  - added the `examples/search_robustness/` example family with Brownian Cartesian/curvilinear baselines plus deterministic Cartesian/curvilinear `UNIFORM_FLOW` migration-stress variants, a study starter, and a dedicated metrics-reference docs page.

- Profiling/logging cleanup:
  - removed the `LOG_PROFILE` log level from the C logging enum and all code paths.
  - per-step profiling summaries no longer use a special log level.
  - profiling is now configured explicitly via `profiling.timestep_output` (`off`, `selected`, `all`) and writes timestep rows to a dedicated profiling log file.
  - `profiling.final_summary.enabled` now controls whether the end-of-run `ProfilingSummary_*.log` file is written.
  - `LOG_LEVEL=PROFILE` is no longer a supported runtime setting.
  - removed the temporary `profiling.critical_functions` compatibility shorthand; monitors must now use `profiling.timestep_output`.

- Executable naming refresh:
  - renamed the conductor from `pic.flow` to `picurv`.
  - renamed the main solver executable from `picsolver` to `simulator`.
  - `init` now creates config-only case directories; binaries are resolved from `bin/` via PATH.
  - `init --pin-binaries` copies `simulator`/`postprocessor` into the case for version-pinning (protects running jobs from concurrent rebuilds).
  - `bin/picurv` is now a launcher for the `picurv_cli/picurv` source-tree entrypoint.
  - `sync-binaries` pins specific binary versions into a case directory (optional, equivalent to `--pin-binaries` after init).
  - removed the temporary `--copy-binaries` init flag.

- Grid contract correction:
  - `programmatic_c.im/jm/km` are now treated as cell counts as documented.
  - `picurv` now converts those values to node counts before emitting `-im/-jm/-km`.
  - `grid_gen` remains unchanged: `grid.gen` still accepts cell counts and writes node counts into `.picgrid`.
  - compatibility note: a historical `programmatic_c` case with `im=32` previously yielded 31 physical cells; it now yields the documented 32 physical cells.

- Interface contract cleanup:
  - removed deprecated numeric and alternate-string selector forms for field initialization, particle initialization, and analytical type selection; docs now use canonical YAML values only.
  - renamed the dual-time smoother from `Dual Time Picard RK4` to the more precise `Dual Time Picard Jameson RK`; deprecated RK4 selector, YAML-block, noise-control, C API, and C CLI spellings remain accepted as compatibility aliases.
  - removed the exposed placeholder Newton-Krylov momentum solver modes from parser and C runtime enums.
  - canonical `PICGRID` headers are required for file-based grids in C runtime ingestion.
  - added `grid.gen legacy1d` converter and optional `grid.legacy_conversion` wrapper in `picurv` for headerless 1D-axis legacy payload migration.

- Workflow launch policy update:
  - `run --num-procs` now applies to solver stage sizing.
  - post stage is forced to single-rank execution in local mode.
  - generated `post.sbatch` now uses forced single-task resources (`nodes=1`, `ntasks_per_node=1`).
  - manifests/dry-run plans now expose stage-specific MPI counts.

- Cluster orchestration and sweeps:
  - added `cluster.yml` Slurm contract support to `picurv run` (`--cluster`, `--scheduler`, `--no-submit`).
  - added `picurv submit` as the delayed-submit counterpart to `--no-submit` for existing run/study artifacts.
  - added `picurv cancel` so Slurm jobs can be stopped by `--run-dir` instead of manual job-id lookup; `picurv cancel --stage solve --graceful` requests solver final-output shutdown with `SIGUSR1`.
  - added scheduler artifact generation and submission metadata (`solver.sbatch`, `post.sbatch`, `submission.json`, `manifest.json`).
  - added `picurv sweep` for parameter studies using Slurm job arrays with post-stage dependency chaining.
  - added `picurv sweep --continue --study-dir <path>` for resuming partially-completed studies: detects per-case completion status, prepares checkpoint restarts via `resolve_restart_source`, and submits sparse solver arrays for incomplete cases only.
  - added `picurv sweep --reaggregate --study-dir <path>` for manual metrics re-aggregation on existing study outputs.
  - added automatic post-completion metrics aggregation via a chained Slurm job (`metrics_aggregate.sbatch`, `afterany` dependency on post array).
  - `detect_last_checkpoint_step` now falls back to particle checkpoint files (`position*.dat`) for analytical-mode cases with no eulerian output.
  - added study aggregation outputs (`metrics_table.csv`, `results/plots`, `summary.json`).
  - added templates: `master_cluster.yml`, `master_study.yml`.
  - added docs pages: cluster run guide and sweep/study guide.
  - documented signal-triggered final snapshot writes for impending walltime/termination handling (`SIGUSR1`, `SIGTERM`, `SIGINT`) with launcher-specific Slurm signal guidance.
  - added automatic runtime walltime guarding for generated Slurm solver jobs, with tunable `execution.walltime_guard` policy and exported batch metadata (`PICURV_JOB_START_EPOCH`, `PICURV_WALLTIME_LIMIT_SECONDS`).

- Documentation enforcement and CI:
  - enforced function-level Doxygen-compatible coverage across C product code, C tests, Python product scripts, and Python tests.
  - added `tests/tooling/audit_function_docs.py` as the repository-wide audit gate.
  - wired the audit into repository consistency tests and the GitHub Actions quality workflow as an explicit pre-pytest step.

- Run inspection and local log routing:
  - added `picurv summarize` for read-only per-step health summaries derived from existing run artifacts.
  - expanded `picurv summarize` with additive `--overview`, `--case`, `--solver`, and `--monitor` configuration views that work before timestep artifacts exist.
  - added `summarize --list-plot-series` and `--plot` time-history workflows backed by standalone `generators/plot.gen`, with append-order windows, automatic residual/norm log scaling, explicit saves, and headless fallback.
  - report missing `matplotlib` as an actionable `DEPENDENCY_MISSING` plotting error instead of a configuration-value error.
  - install `matplotlib` by default in the bootstrap-managed Python environment so plotting works after a standard install.
  - preserve the existing pip version during normal bootstrap runs; add explicit `--upgrade-pip` for environments that need it.
  - verify Python dependencies under the launcher's isolated environment and report the specific missing plotting dependency.
  - local wrapper stream logs for solver/post stages now write under `runs/<run_id>/scheduler/` instead of `runs/<run_id>/logs/`.
  - summary lookup can use continuity, particle metrics, momentum, Poisson, profiling, and sampled particle snapshot artifacts when present.
  - sampled particle snapshot summaries now include compact speed/bounds/rank/weight diagnostics plus sampled deltas against the previous snapshot when possible.

- Initial Doxygen docs scaffold, architecture, and developer guide pages.
- Main docs refresh:
  - updated tutorials/references for current YAML -> `picurv` -> C contract.
  - added method overview pages (CurvIB, fractional-step, dual-time Jameson RK, pressure Poisson/multigrid, walking search, interpolation/projection, IEM/averaging).
  - added non-dimensionalization page and linked it into nav/reference flow.
  - updated landing page visuals (`docs/assets/curv.gif`, `docs/assets/paraview_flat_channel.png`).
  - added maintenance backlog page for low-priority warning/refactor tracking.
- Repository organization updates:
  - moved build configs to `config/build/`.
  - moved shared `grid.gen` profile to `config/grids/coarse_square_tube_curved.cfg`.
  - converted `sandbox/` into explicit developer-sandbox documentation.
  - standardized internal directory guides to `guide.md` (non-root).
  - moved documentation media into `docs/assets/`.
  - removed non-root `README.md` files (single top-level `README.md` retained).
  - routed Doxygen warning output to `logs/doxygen.warnings`.
- Hardened YAML -> C config contract across case/solver/monitor/post:
  - scalar-only `da_processors_*` contract with validation.
  - structured mappings for analytical type and monitor subdirectories.
  - post input extension keys and statistics pipeline wiring.
- Added config contract documentation set:
  - `14_Config_Contract.md`
  - `15_Config_Ingestion_Map.md`
  - `16_Config_Extension_Playbook.md`
- Added ingress drift guard:
  - `tests/tooling/audit_ingress.py`
  - `tests/tooling/audit_ingress_manifest.json`
- Added documentation note on data-driven particle-closure expansion:
  - offline workflows supported now via solver/post artifacts.
  - tightly coupled inference requires runtime-selectable closure interface extension.
- Removed obsolete `stubs/` archive from repository after extracting useful documentation content.
