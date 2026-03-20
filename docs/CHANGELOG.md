@page 18_Changelog Changelog

@anchor _Changelog

# Changelog

## Unreleased

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
  - `bin/picurv` is now a symlink to `scripts/picurv` (single source of truth).
  - `sync-binaries` pins specific binary versions into a case directory (optional, equivalent to `--pin-binaries` after init).
  - removed the temporary `--copy-binaries` init flag.

- Grid contract correction:
  - `programmatic_c.im/jm/km` are now treated as cell counts as documented.
  - `picurv` now converts those values to node counts before emitting `-im/-jm/-km`.
  - `grid_gen` remains unchanged: `grid.gen` still accepts cell counts and writes node counts into `.picgrid`.
  - compatibility note: a historical `programmatic_c` case with `im=32` previously yielded 31 physical cells; it now yields the documented 32 physical cells.

- Interface contract cleanup:
  - removed deprecated numeric and alternate-string selector forms for field initialization, particle initialization, and analytical type selection; docs now use canonical YAML values only.
  - removed deprecated direct/internal momentum solver selector forms; `strategy.momentum_solver` now accepts only `Explicit RK4` and `Dual Time Picard RK4`.
  - removed the exposed placeholder Newton-Krylov momentum solver modes from parser and C runtime enums.
  - canonical `PICGRID` headers are required for file-based grids in C runtime ingestion.
  - added `grid.gen legacy1d` converter and optional `grid.legacy_conversion` wrapper in `picurv` for headerless 1D-axis legacy payload migration.

- Workflow launch policy update:
  - `run --num-procs` now applies to solver stage sizing.
  - post stage defaults to single-rank execution in local mode.
  - generated `post.sbatch` now defaults to single-task resources (`nodes=1`, `ntasks_per_node=1`).
  - manifests/dry-run plans now expose stage-specific MPI counts.

- Cluster orchestration and sweeps:
  - added `cluster.yml` Slurm contract support to `picurv run` (`--cluster`, `--scheduler`, `--no-submit`).
  - added `picurv submit` as the delayed-submit counterpart to `--no-submit` for existing run/study artifacts.
  - added `picurv cancel` so Slurm jobs can be stopped by `--run-dir` instead of manual job-id lookup.
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
  - added `scripts/audit_function_docs.py` as the repository-wide audit gate.
  - wired the audit into repository consistency tests and the GitHub Actions quality workflow as an explicit pre-pytest step.

- Run inspection and local log routing:
  - added `picurv summarize` for read-only per-step health summaries derived from existing run artifacts.
  - local wrapper stream logs for solver/post stages now write under `runs/<run_id>/scheduler/` instead of `runs/<run_id>/logs/`.
  - summary lookup can use continuity, particle metrics, momentum, Poisson, profiling, and sampled particle snapshot artifacts when present.
  - sampled particle snapshot summaries now include compact speed/bounds/rank/weight diagnostics plus sampled deltas against the previous snapshot when possible.

- Initial Doxygen docs scaffold, architecture, and developer guide pages.
- Main docs refresh:
  - updated tutorials/references for current YAML -> `picurv` -> C contract.
  - added method overview pages (CurvIB, fractional-step, dual-time RK4, pressure Poisson/multigrid, walking search, interpolation/projection, IEM/averaging).
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
  - `scripts/audit_ingress.py`
  - `scripts/audit_ingress_manifest.json`
- Added documentation note on data-driven particle-closure expansion:
  - offline workflows supported now via solver/post artifacts.
  - tightly coupled inference requires runtime-selectable closure interface extension.
- Removed obsolete `stubs/` archive from repository after extracting useful documentation content.
