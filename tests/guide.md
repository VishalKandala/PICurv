# Tests Guide

PICurv testing is intentionally layered: Python control-plane validation, PETSc installation checks, focused C unit suites, executable smoke runs, MPI variants, and coverage gates. Choose the narrowest layer that answers your current question.

Function-level docs in `tests/c/` are part of the repository contract. For new
or touched C tests/helpers, keep the Doxygen blocks concise but descriptive:
the `@brief` should state what the routine verifies or sets up, not just that it
is a test-local routine.

## Canonical Targets

- repository/doc contract:
  - `python3 tests/tooling/audit_function_docs.py`
- Python and coverage:
  - `make test-python` (`make test` alias)
  - `make coverage-python`
  - `make coverage-c`
  - `make coverage`
- installation/toolchain:
  - `make doctor` (`make install-check` alias)
- C unit suites:
  - `make unit`
  - `make unit-geometry`
  - `make unit-setup`
  - `make unit-solver`
  - `make unit-particles`
  - `make unit-io`
  - `make unit-logging`
  - `make unit-post`
  - `make unit-post-compute-mpi`
  - `make unit-post-eulerian-vtk-mpi`
  - `make unit-post-particle-vtk-mpi`
  - `make unit-grid`
  - `make unit-metric`
  - `make unit-boundaries`
  - `make unit-poisson-rhs`
  - `make unit-runtime`
  - `make unit-simulation`
  - `make unit-mpi`
  - `make unit-periodic`
  - `make unit-periodic-dev`
- momentum solver suites:
  - `make unit-newton-krylov`
  - `make unit-momentum-candidates`
  - `make unit-momentum-newton-boundary-fixedpoint`
- field-statistics suites:
  - `make unit-statistics`
  - `make unit-statistics-target`
  - `make unit-statistics-window`
  - `make unit-statistics-accumulator`
  - `make unit-statistics-config`
- smoke/integration:
  - `make smoke`
  - `make smoke-mpi`
  - `make smoke-mpi-matrix`
  - `make smoke-stress`
  - `make smoke-periodic`
  - `make smoke-periodic-dev`
  - `make smoke-driven-periodic`
- documentation gates:
  - `make audit-agent-setup`
  - `make audit-ingress`
  - `make audit-docs-expansion`
  - `make certify-docs-fast`
  - `make certify-docs`
- aggregate gates:
  - `make check`
  - `make check-mpi`
  - `make check-mpi-matrix`
  - `make check-full`
  - `make check-stress`

## Python Test Files (`tests/test_*.py`)

- `test_cli_smoke.py`
  - CLI help and argument contract checks
  - dry-run plan schema checks (text/json)
  - staged Slurm workflow coverage for `submit`, `cancel`, and `sweep`
  - run summary coverage for `summarize` JSON/text output, time-history discovery/plot orchestration, standalone `plot.gen`, and failure paths
  - restart path resolution checks
  - cluster no-submit manifest/script checks
  - grid-gen/PICGRID header and node-count translation checks
  - case-local binary preference behavior for copied or path-resolved `picurv`
- `test_case_maintenance.py`
  - `init` origin metadata behavior
  - real CLI wrapper coverage for `build`, `sync-config`, `status-source`, and `pull-source`
  - source-root resolution for build/sync/pull commands
  - template sync behavior (`overwrite`, `prune`)
  - source/case drift reporting (`status-source`)
- `test_config_regressions.py`
  - ingress-manifest drift checks
  - post recipe alias compatibility
  - post validation guards and statistics artifact pathing
- `test_repo_consistency.py`
  - validates example bundles and study bundles via `picurv validate`
  - scans docs/examples/tests for stale/forbidden contract literals
  - wraps the repository-wide function documentation audit script
- `test_review_packet.py`
  - page mode renders every published page, including pages owning no invariant contract
  - page mode reports only the freshness surfaces that route review to that page
  - page-mode output does not depend on the order contracts are iterated
  - contract mode renders every registered contract and rejects unknown ids
- `test_storage.py`
  - rclone backend round trip, protect/offload/restore, and remote-only recovery
  - incomplete-checkpoint archival refusal and partial-restore markers
  - cold-run guards so `post` and `submit` refuse before missing data is misdiagnosed
- `test_run_directory_containment.py`
  - escaping log/output directories are errors; restart directories are exempt
  - override semantics, reserved-directory collisions, and run-root rejection
- `test_c_guard_behaviour.py`
  - C write-authorization guard: contained paths accepted, run-root, ancestor, and
    relative-symlink escapes refused, and the escapes authorization never waives
- `test_repo_file_enumeration.py`
  - enumeration used by every scanner: deletions, staged deletions, untracked files,
    honoured skip directories, and the non-Git fallback
  - each scanner runs against the current dirty tree
- `test_post_spectra_config.py`
  - spectra recipe normalization, task identity, and stage selection
  - periodic-box preconditions, and refusal of wall-bounded or multi-block domains
  - recipe signature, predicted artifacts, and staged follow commands
- `test_spectra_shell_spectrum.py`
  - determinism, prescribed-envelope recovery, and Parseval closure for both symbols
  - DMDA interior extraction, and refusal of stretched or curvilinear grids
- `test_spectral_random_velocity.py`
  - solenoidal projection, component-energy balance, and directional isotropy
  - PETSc axis mapping, Nyquist reality, and subordination to a restarted state
- `test_profile_field_slice.py`
  - cell-centered DMDA slice length, dummy layers never sampled, and each face
    reading its adjacent interior cell
- `test_picgrid.py`
  - PICGRID writers share round-trip-safe precision
- `test_grid_generator.py`
  - shared generator stages: axis construction, `origin` placement, and the ordered
    placement/similarity transform list
  - every transform kind is reachable, and every malformed token names its reason
  - handedness is reported separately from Jacobian sign consistency, which a
    uniformly inverted grid satisfies
  - `box` piecewise wall-height fields: segment continuity, a step moving only its
    own wall, a hill closing on itself, and the spanwise envelope giving an
    obstacle finite width
  - a sub-cell corner is refused, a coarse one is reported, and a sharper corner
    costs orthogonality
  - `sweep` cross-sections and centrelines: a straight rectangular sweep is a box,
    the disc map fills the circle with no hole, and swept grids are right-handed
  - parallel transport keeps an out-of-plane path untwisted, and an arc about its
    own tangent is refused
- `test_capability_tooling.py`
  - capability parity, value metadata, and coverage-entry shape
  - contract registry and artifact-topology snapshot validation
  - documentation scaffolding output, and repository-wide fragment-link resolution
- `test_field_catalog.py`
  - the published field inventory matches the compiled Eulerian and particle catalogs
  - missing fields, wrong layout groups, undocumented layout values, and renamed
    inventory groups are each rejected
  - the invariant contract registers this audit as its blocking checker
- `test_subsystem_lifecycle.py`
  - what each status owes, illegal transitions, and reasoned "not applicable"
  - dangling obligation page and anchor references are rejected
- `test_freshness.py`
  - manifest structure, hard-blocking versus soft-advisory drift, and promotion rules
  - attestation digests are path-sensitive, order-independent, and never a commit sha
- `test_inline_choices.py`
  - unnamed inline choice sets and argparse `choices` lists are detected
  - every waiver is typed, reasoned, and claims no pending family
- `test_path_notation.py`
  - bare run-owned prefixes rejected, logical identities accepted
  - runnable command blocks keep concrete paths; repository directories survive
- `test_family_census_classifications.py`
  - every census entry typed and reasoned, with parameter and alias entries naming
    a real owner
  - discovery covers every authoritative constant
- `test_documentation_counts.py`
  - counts stated in prose match the registries they describe: contracts, logical
    identities, containment points, and census status lines
  - capability, subsystem, and freshness-surface modes render every declared identifier
  - optional xrefs are accepted only when schema and dirty-byte stamps are current
  - changed-path classification distinguishes routed production paths from declared gaps

## C Unit Files (`tests/c/test_*.c`)

- `test_install_check.c`: PETSc environment and basic object viability (`doctor`: `environment-visible`, `basic-petsc-objects`)
- `test_geometry.c`: interpolation and geometric signed-distance helpers
- `test_setup_lifecycle.c`: setup/cleanup lifecycle, RNG, and initialized particle-settlement contracts
- `test_solver_kernels.c`: analytical geometry/particle dispatch, LES filter/eddy-viscosity, FlowSolver guardrails, and analytical source helpers
- `test_particle_kernels.c`: walking-search helper kernels
- `test_io.c`: I/O path checks, parser helpers, scaling-ingestion contracts, and startup-banner summary contracts
- `test_logging.c`: log-level, allow-list, string-conversion, continuity/min-max/interpolation diagnostics, profiling, and snapshot-cadence contracts
- `test_postprocessing.c`: post-processing kernel contracts (specific-KE, displacement, nodal average, normalization, dimensionalization, Q-criterion)
- `test_vtk_io.c`: VTK writer and data-preparation contracts (coordinates, field gather/subsampling, particle prep)
- `test_postprocessor.c`: postprocessing orchestration contracts (swarm setup, pipeline dispatch, eulerian/particle output, statistics dispatch)
- `test_post_compute_mpi.c`: analytic serial/MPI checks for Eulerian, particle, resizing, and derived field-statistics compute
- `test_post_eulerian_vtk_mpi.c`: exhaustive serial/MPI Eulerian VTK value and byte-equivalence checks
- `test_post_particle_vtk_mpi.c`: exhaustive serial/MPI particle VTK value and byte-equivalence checks
- `test_statistics.c`: particle statistics kernel contracts (MSD CSV output and empty-swarm behavior)
- `test_statistics_moments.c`: weighted centered moment, co-moment, and parallel-merge kernels
- `test_statistics_window.c`: window scheduling, weighting, clipping, and identity hashing
- `test_statistics_accumulator.c`: per-window storage, accumulation, masks, and derived quantities
- `test_statistics_target.c`: spatial target plan across cell, node, face, and periodic layouts
- `test_statistics_config.c`: field-statistics control resolution from the generated control file
- `test_grid.c`: local/global bounding-box helpers
- `test_metric.c`: metric inversion, contravariant velocity, face geometry helpers
- `test_boundaries.c`: boundary factory, direct handler-behavior checks, and wall-model dispatch across the log law, Werner-Wengle, and Cabot
- `test_les.c`: LES closure kernels - symmetric-tensor algebra, strain rate, filter width, the Germano model tensor and Leonard stress, coefficient clipping and the viscosity floor, and homogeneous-direction resolution
- `test_periodic_dev.c`: gating geometric-periodic boundary and synchronization harness
- `test_poisson_rhs.c`: pressure update, RHS, projection, body-force and diffusivity helpers
- `test_runtime_kernels.c`: setup/runloop/particle/interpolation/scatter/wall/walltime-guard/LES helper contracts
- `test_mpi_kernels.c`: multi-rank particle distribution, bbox collectives, and restart migration behavior
- `test_momentum_newton_krylov.c`: matrix-free Newton-Krylov momentum callbacks and helper contracts
- `test_momentum_convective_candidates.c`: finite-difference Jacobian study of the production convection-only residual
- `test_momentum_newton_boundary_fixedpoint.c`: opt-in regression for the Newton-Krylov Cartesian-seed correction on a production-sized straight duct
- shared fixture layer:
  - `test_support.c`
  - `test_support.h`
  - exposes both a fast minimal fixture and a richer tiny-runtime fixture built through the real setup path
  - mirrors the production `da/fda/swarm` contract (`da = IM+1/JM+1/KM+1`, coordinate-DM `fda`, production swarm fields)

## Smoke Harness (`tests/smoke/run_smoke.sh`)

Single-rank smoke (`make smoke`) verifies:

- binary `-help` launch viability (`simulator`, `postprocessor`)
  - banner presence is authoritative; local PETSc may exit with code `62` (`PETSC_ERR_ARG_WRONG`) and still satisfy the smoke contract
- `picurv init` self-contained case creation and metadata
- template matrix `init + validate + dry-run` checks (`flat_channel`, `bent_channel`, `brownian_motion`)
- dry-run plan schema and restart-source resolution
- a tiny PETSc diagnostics solve with `malloc_debug`, `malloc_dump`, `malloc_view`, and `log_view`
- tiny real solve+post for flat and bent channels
- tiny particle solve+post and restart branches (`load`, `init`)
- restart-equivalence continuity check
- tiny analytical Brownian run with VTP + MSD CSV checks

Opt-in stress smoke (`make smoke-stress`) additionally verifies:

- longer particle-cycle runtime sequences
- chained restart workflows over more than one split point
- BC-focused runtime variants for the stable path (parabolic inlet) plus periodic constant-flux validate/dry-run coverage
- an extra-rank MPI particle runtime stress case

Geometric-periodic smoke (`make smoke-periodic`) is part of the standard `check` gate:

- real runtime solve+post coverage for the current periodic constant-flux development path
- verifies the supported Eulerian geometric-periodic runtime path

Multi-rank smoke (`make smoke-mpi`, `make smoke-mpi-matrix`) additionally verifies:

- rank-dependent runtime launch behavior
- flat/bent multi-rank tiny solves
- particle restart branches under multi-rank execution

Useful env knobs:

- `TEST_MPI_NPROCS` for `unit-mpi`
- `SMOKE_MPI_NPROCS` for `smoke-mpi`
- `SMOKE_MPI_MATRIX_NPROCS` for `smoke-mpi-matrix`
- `KEEP_SMOKE_TMP=1` to preserve smoke temp workspace for debugging

## Suggested Command Cadence

- editing `picurv_cli/` or YAML contracts:
  - `make test-python`
- editing C/Python functions or helper/test docstrings:
  - `python3 tests/tooling/audit_function_docs.py`
- editing shared agent instructions or canonical skills:
  - `make sync-agent-skills`
  - `make audit-agent-setup`
- editing one C subsystem:
  - targeted `make unit-<area>`
- editing setup/teardown lifecycle code:
  - `make unit-setup`
- editing PETSc object ownership, monitors, or memory diagnostics:
  - `make smoke`
- editing core simulation numerics or particle orchestration:
  - `make unit-simulation`
- editing periodic BC code under development:
  - `make unit-periodic`
  - `make smoke-periodic`
- editing runtime orchestration, restart, or output contracts:
  - `make smoke` plus MPI variant if rank behavior is involved
- running the opt-in medium-budget extension tier:
  - `make smoke-stress`
- pre-merge:
  - `make check` (or `make check-mpi`)
- pre-release:
  - `make check-full`
  - `make coverage`

## Current targeted backlog

- walking search:
  - add direct unit coverage for `LocateParticleOrFindMigrationTarget` boundary-clamp, ghost-handoff, tie-breaker, `LOST`, and `MIGRATING_OUT` branches
  - add explicit direction-complete and failure-path checks for the `GuessParticleOwnerWithBBox` heuristic
- particle migration:
  - add non-restart MPI migration tests for multi-pass handoff, newcomer flagging, and particle-count conservation
- momentum:
  - add direct positive-path tests for `MomentumSolver_Explicit_RungeKutta4`
  - add one small direct invariant harness for `MomentumSolver_DualTime_Picard_JamesonRK`
- pressure/Poisson:
  - add deeper `PoissonSolver_MG` and periodic/IBM stencil checks beyond the current `Projection`/`PoissonLHSNew` helper surface
- grid/metrics/setup:
  - broaden the richer runtime fixture to more geometry/topology variants; current setup coverage is production-faithful but still mostly tiny Cartesian cases
- coverage follow-up:
  - `2026-03-20` audit snapshot: `src/AnalyticalSolutions.c` (`23.93%`) and `src/BodyForces.c` (`11.69%`) still need direct unit harnesses
  - `src/les.c` was rewritten after the `2026-03-20` snapshot, so that figure no longer describes the file; `test_les.c` covers the kernels directly, and the assembled model path has not been re-measured
  - `2026-03-20` audit snapshot: `src/poisson.c` (`45.66%`) and `src/rhs.c` (`61.95%`) still need deeper direct branch coverage
  - `2026-03-20` audit snapshot: `src/runloop.c` (`64.76%`) still needs more runtime-orchestration branch coverage
  - `2026-03-20` audit snapshot: `src/Boundaries.c` (`68.10%`) still needs additional non-periodic boundary edge-case coverage
- periodic BC:
  - add nontrivial curvilinear seam and deeper periodic Poisson coverage

## Notes

- Python tests do not require PETSc.
- GitHub Actions quality CI runs `python tests/tooling/audit_function_docs.py`, then `pytest -q`, then markdown link checks.
- `doctor`, `unit-*`, `smoke*`, `check*`, and `coverage-c` require PETSc/MPI tooling.
- `check-full` is the single-command comprehensive gate (`check` + `unit-mpi` + `smoke-mpi` + `smoke-mpi-matrix`).
- `check-stress` extends `check-full` with the opt-in `smoke-stress` layer.
- geometric-periodic harnesses are included in the standard gates.
- compatibility aliases (`install-check`, `ctest-*`) still exist, but canonical names are preferred in docs and CI.

## Authoritative Docs

- https://vishalkandala.me/docs/picurv/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/docs/picurv/51_C_Test_Suite_Developer_Guide.html
