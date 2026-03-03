@page 18_Changelog Changelog

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
  - `init` now produces self-contained case directories built around `picurv`, `simulator`, and `postprocessor`.
  - removed the temporary `--copy-binaries` init flag; copying built executables is now unconditional.

- Grid contract correction:
  - `programmatic_c.im/jm/km` are now treated as cell counts as documented.
  - `picurv` now converts those values to node counts before emitting `-im/-jm/-km`.
  - `grid_gen` remains unchanged: `grid.gen` still accepts cell counts and writes node counts into `.picgrid`.
  - compatibility note: a historical `programmatic_c` case with `im=32` previously yielded 31 physical cells; it now yields the documented 32 physical cells.

- Interface contract cleanup:
  - removed deprecated numeric and alternate-string selector forms for field initialization, particle initialization, and analytical type selection; docs now use canonical YAML values only.
  - removed deprecated direct/internal momentum solver selector forms; `strategy.momentum_solver` now accepts only `Explicit RK4` and `Dual Time Picard RK4`.
  - removed the exposed placeholder Newton-Krylov momentum solver modes from parser and C runtime enums.
  - canonical `PICGRID` headers are now required for file-based grids in both `picurv` and the C runtime; convert older grids before use.

- Workflow launch policy update:
  - `run --num-procs` now applies to solver stage sizing.
  - post stage defaults to single-rank execution in local mode.
  - generated `post.sbatch` now defaults to single-task resources (`nodes=1`, `ntasks_per_node=1`).
  - manifests/dry-run plans now expose stage-specific MPI counts.

- Cluster orchestration and sweeps:
  - added `cluster.yml` Slurm contract support to `picurv run` (`--cluster`, `--scheduler`, `--no-submit`).
  - added scheduler artifact generation and submission metadata (`solver.sbatch`, `post.sbatch`, `submission.json`, `manifest.json`).
  - added `picurv sweep` for parameter studies using Slurm job arrays with post-stage dependency chaining.
  - added study aggregation outputs (`metrics_table.csv`, `results/plots`, `summary.json`).
  - added templates: `master_cluster.yml`, `master_study.yml`.
  - added docs pages: cluster run guide and sweep/study guide.

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
