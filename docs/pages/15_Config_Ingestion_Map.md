@page 15_Config_Ingestion_Map Developer Ingestion Map

@anchor _Config_Ingestion_Map

This page maps configuration flow from YAML schema to generated artifacts and C ingestion/use sites.

@tableofcontents

@section p15_pipeline_sec 1. End-to-End Flow

1. `picurv` validates YAML inputs.
2. `picurv` generates `*.control`, `bcs*.run`, and `post.run`, plus optional
   `whitelist.run` / `profile.run` sidecars when those features are enabled.
3. `CreateSimulationContext()` in `src/setup.c` loads `-control_file` and ingests options.
4. Secondary file parsers ingest BC/post/logging profile inputs.
5. Runtime modules consume parsed values during setup, solve, and postprocess.

@section p15_map_sec 2. Mapping Matrix

| YAML / Source Key | Generated Artifact Key | C Ingestion Site | Primary Runtime Consumers |
| :--- | :--- | :--- | :--- |
| `case.run_control.*` | `-start_step`, `-totalsteps`, `-dt` | `src/setup.c` (`CreateSimulationContext`) | `src/runloop.c`, setup/timestep logic |
| `case.grid.programmatic_settings.im/jm/km/...` | `-im/-jm/-km`, bounds, stretch (`im/jm/km` translated from YAML cell counts to C node counts) | `src/io.c` (`ReadGridGenerationInputs`, `PopulateFinestUserGridResolutionFromOptions`) | `src/grid.c`, analytical geometry setup |
| `case.grid.da_processors_*` | `-da_processors_x/y/z` | `src/setup.c` | `src/grid.c` DMDA creation |
| `case.models.domain.blocks` | `-nblk` | `src/setup.c` | `src/grid.c` |
| paired `case.boundary_conditions` periodic faces | BC rows; periodic axes are derived in C | `src/Boundaries.c`, `src/io.c` (`DeterminePeriodicity`) | `src/grid.c` DMDA creation, periodic field synchronization |
| `case.models.physics.particles.*` | `-numParticles`, `-pinit`, `-particle_restart_mode`, `-psrc_*` | `src/setup.c` | `src/ParticleSwarm.c`, `src/ParticleMotion.c`, statistics kernels |
| `case.models.physics.turbulence.*` | `-les`, `-les_constant_cs`, `-les_dynamic_frequency`, `-les_filter_width`, `-les_test_filter_kernel`, `-les_test_filter_width_ratio`, `-les_averaging_mode`, `-les_averaging_directions`, `-les_clip_mode`, `-les_clip_max_cs`, `-les_min_viscosity_ratio`, `-les_yoshizawa_ci`, `-les_gradient_model`, `-les_diagnostics`, `-les_diagnostics_cadence`, `-rans`, `-wallfunction`, `-wall_roughness` | `src/setup.c` | `src/solvers.c`, `src/les.c`, `src/Filter.c`, `src/rhs.c`, `src/Boundaries.c`, `src/wallfunction.c` |
| `solver.interpolation.method` | `-interpolation_method` | `src/setup.c` | `src/interpolation.c` (dispatch in @ref InterpolateEulerFieldToSwarm) |
| `case.boundary_conditions` | `bcs*.run` rows | BC parser path (`src/Boundaries.c` + helpers) | BC handler factory and boundary application |
| `solver.operation_mode.*` | `-euler_field_source`, `-analytical_type`, `-analytical_uniform_u/-v/-w` | `src/setup.c` | `src/initialcondition.c`, `src/grid.c`, `src/AnalyticalSolutions.c` |
| `solver.verification.sources.diffusivity.*` | `-verification_diffusivity_mode/-profile/-gamma0/-slope_x` | `src/setup.c` | `src/verification_sources.c`, `src/rhs.c` |
| `solver.verification.sources.scalar.*` | `-verification_scalar_mode/-profile/-value/-phi0/-slope_x/-amplitude/-kx/-ky/-kz` | `src/setup.c` | `src/verification_sources.c`, `src/AnalyticalSolutions.c`, `src/ParticlePhysics.c`, `src/logging.c` |
| `solver.strategy/tolerances/momentum_solver.*` | solver flags (`-mom_*`, pseudo-CFL, etc.) | `src/setup.c` | `src/momentumsolvers.c` |
| `solver.momentum_solver.newton_krylov.*` | application selectors `-mom_nk_jacobian_*`, `-mom_nk_preconditioner_*`; prefixed `-mom_nk_snes_*`, `-mom_nk_ksp_*` | application parsing plus PETSc options db via `SNESSetFromOptions()` | `src/momentum_newton_krylov.c` |
| `solver.scalar_transport.*` | `-schmidt_number`, `-turb_schmidt_number` | `src/setup.c` | `src/rhs.c`, `src/particle_statistics.c`, scalar/particle transport |
| `solver.poisson_solver.*` / legacy `solver.pressure_solver.*` | `-poisson_tol`, `-ps_ksp_*`, `-ps_pc_type`, `-mg_*`, `-ps_mg_levels_*` | `src/setup.c` + PETSc options db | `src/poisson.c` |
| `solver.petsc_passthrough_options` | raw flags in control | PETSc options db | PETSc KSP/PC stack, mostly in `src/poisson.c` |
| `monitor.io.data_output_frequency` | `-tio` | `src/setup.c` | `src/io.c`, `src/setup.c`, `src/runloop.c` |
| `monitor.io.particle_console_output_frequency` | `-particle_console_output_freq` | `src/setup.c` | `src/io.c`, `src/setup.c`, particle console logging |
| `monitor.io.particle_log_interval` | `-logfreq` | `src/setup.c` | particle console row subsampling |
| `monitor.io.statistics_console_output_frequency` | `-statistics_console_output_freq` | `src/statistics_config.c` | `src/logging.c` statistics console snapshot, `src/runloop.c` |
| `monitor.io.directories.*` | `-output_dir`, `-restart_dir`, `-log_dir` | `src/setup.c` | `src/io.c`, `src/setup.c`, `src/runloop.c` |
| `monitor.logging.*` | `whitelist.run`, `LOG_LEVEL` env | `src/setup.c` + `src/logging.c` | logging macros/system |
| `monitor.profiling.*` | `profile.run` (selected-mode only) + explicit profiling flags in `*.control` | `src/setup.c` + profiling init | profiler summaries |
| `monitor.diagnostics.petsc.*` | solver/postprocessor executable arguments (`-malloc_*`, `-log_view`, `-objects_dump`, etc.) | PETSc initialization | PETSc memory/log/object diagnostics |
| `monitor.diagnostics.runtime_memory_log.*` | `-runtime_memory_log_enabled`, `-runtime_memory_log_file` | `src/setup.c` | `src/logging.c`, `src/runloop.c`, `src/postprocessor.c`, `src/simulator.c` |
| `monitor.field_statistics.*` | `-field_statistics_enabled`, `-field_statistics_window_count`, `-field_statistics_window_<i>_*` | `src/statistics_config.c` (`ParseFieldStatisticsConfig`) | `src/statistics_window.c`, `src/statistics_accumulator.c`, `src/runloop.c`, `src/io.c` checkpoint state |
| `monitor.solution_monitoring.convergence.*` | `-solution_convergence_enabled/-solution_convergence_mode/-solution_convergence_period_steps/-solution_convergence_window_steps` in `*.control` | `src/setup.c` | existing `src/logging.c` and `src/runloop.c` convergence path; every completed timestep |
| `cluster.execution.walltime_guard.*` | `-walltime_guard_*` in solver `*.control` | `src/setup.c` | `src/runloop.c` runtime walltime estimator and graceful final-write cutoff |
| `monitor.solver_monitoring.poisson.*` | prefixed Poisson monitor flags in control (`-ps_ksp_*`) | PETSc options db | `src/poisson.c` KSP monitor setup |
| `monitor.solver_monitoring.momentum.*` | bare Newton history/SNES/KSP monitor flags (`-mom_nk_*`) | `src/setup.c` + PETSc options db | `src/momentum_newton_krylov.c` |
| `monitor.solver_monitoring.petsc_passthrough_options` | raw flags in control | PETSc options db | PETSc monitors/debug options |
| `post.run_control.*` | `startTime/endTime/timeStep` in `post.run` | `src/io.c` (`ParsePostProcessingSettings`) | `src/postprocessor.c` main loop |
| `post.io.input_extensions.*` | `eulerianExt`, `particleExt` in `post.run` | `src/io.c` | `ReadSimulationFields`, `ReadAllSwarmFields`, swarm precheck |
| `post.statistics_pipeline.*` | `statistics_pipeline`, `statistics_output_prefix` | `src/io.c` | `GlobalStatisticsPipeline` dispatch |
| `post.field_statistics.*` | `field_statistics_windows`, `field_statistics_outputs`, `field_statistics_formats`, `field_statistics_source_step` | `src/io.c` (`ParsePostProcessingSettings`) | `FieldStatisticsPipeline` in `src/postprocessor.c`, derivation in `src/statistics_accumulator.c` |
| `post.spectra.*` (digest only) | `spectra_signature` | `src/io.c` (`ParsePostProcessingSettings`) | accepted and ignored; exists so a changed spectra recipe reaches the `--continue` recipe fingerprint |

@section p15_python_only_sec 3. Python-Only Orchestration Mapping (No C Ingestion)

These keys are consumed by `picurv` orchestration only:

| YAML / Source Key | Generated Artifact Key | Ingestion Site | Runtime Consumer |
| :--- | :--- | :--- | :--- |
| `cluster.scheduler.type` | scheduler selector | `picurv_cli/core.py` (`validate_cluster_config`) | `run`/`sweep` scheduler dispatch |
| `cluster.resources.*` | `#SBATCH` directives | `picurv_cli/core.py` (Slurm renderers) | Slurm scheduler |
| `cluster.notifications.*` | `#SBATCH --mail-*` | `picurv_cli/core.py` | Slurm scheduler |
| `cluster.execution.module_setup` | pre-launch shell lines in `*.sbatch` | `picurv_cli/core.py` | batch script runtime env |
| `cluster.execution.launcher*` | launch command (`srun`/`mpirun`) | `picurv_cli/core.py` | solver and field postprocessor launches |
| `study.base_configs.*` | per-case config materialization | `picurv_cli/core.py` (`sweep_workflow`) | case generation pipeline |
| `study.parameters` | case matrix expansion | `picurv_cli/core.py` (`expand_parameter_matrix`) | study case synthesis |
| `study.execution.max_concurrent_array_tasks` | Slurm array `%N` throttle | `picurv_cli/core.py` | Slurm scheduler |
| `study.metrics` | `metrics_table.csv` extraction contract | `picurv_cli/core.py` metric extractors | study aggregation/reporting |
| `study.plotting` | `results/plots/*` output controls | `picurv_cli/core.py` plotting pipeline | study reporting |
| `post.spectra.tasks[*]` | `spectra.gen shell-spectrum` arguments | `picurv_cli/core.py` (`normalize_post_spectra_config`) | `run_post_spectra_stage`, `generators/spectra.gen` |
| `post.spectra.output_prefix` | spectra CSV basenames under `<monitor output>/spectra/` | `picurv_cli/core.py` (`post_spectra_task_basename`) | `run_post_spectra_stage` |

@section p15_exceptions_sec 4. Important Exceptions

- PETSc runtime option consumption is not only explicit `PetscOptionsGet*`.
- `KSPSetFromOptions` in `src/poisson.c` ingests prefixed PETSc options dynamically.
- `monitor.diagnostics.petsc.*` is intentionally emitted as executable
  arguments rather than control-file rows so PETSc initialization-time options
  are available before runtime setup reads `*.control`.
- `LOG_LEVEL` is environment-driven (`src/logging.c`) and intentionally outside control-file parsing.

@section p15_mode_dependent_sec 5. Mode-Dependent Normalization in `picurv`

Some launcher behaviors depend on other config selections before values ever reach C:

- `case.properties.initial_conditions.mode: generated` resolves `generator` and `params` into
  built-in C flags or a staged `ic_gen` PETSc vector.
- `case.properties.initial_conditions.mode: file` validates and stages one `Ucat` or `Ucont`
  PETSc vector in the existing @ref ReadFieldData naming layout.
- file-backed ICs are rejected for multi-block cases in the first implementation.
- `solver.operation_mode.eulerian_field_source` and `case.run_control.start_step` supersede IC
  materialization when load, analytical, or restart state has authority.
- `solver.operation_mode.eulerian_field_source: analytical`
  routes through the analytical grid-ingestion split:
  `TGV3D` still requires `case.grid.mode: programmatic_c`, while `ZERO_FLOW` and `UNIFORM_FLOW`
  accept `case.grid.mode: programmatic_c` or `case.grid.mode: file`.
- `case.properties.initial_conditions.mode` is an explicit launcher requirement.

This means the YAML contract is intentionally stricter than "whatever C would do with missing
options" in several places.

@section p15_safety_sec 6. Safety-Critical Ingress

`-allow_unsafe_log_dir` is generated **only** when
`monitor.io.directories.allow_unsafe_paths` is a real YAML boolean `true`. It carries an
explicit authorization across the Python/C boundary: the runtime performs its own
containment check before recursively deleting the log directory, and waives the
external-location restriction only when this option is present. It never waives a
run-root target, a reserved run directory, a log/output overlap, or a malformed value.

It is a reserved flag: raw PETSc passthrough cannot set it, so authorization cannot be
forged. See @ref p71_isolation_sub.

@section p15_maintenance_sec 7. Drift Prevention

- Use `tests/tooling/audit_ingress.py` to compare PETSc option ingress with the maintained
  manifest. The manifest's own `sources` list is authoritative for which files are scanned;
  it currently covers `src/setup.c`, `src/io.c`, and `src/statistics_config.c`.
- Keep this map and the manifest updated whenever new options are introduced.
- Run manually with:
  - `python3 tests/tooling/audit_ingress.py`
  - or `make audit-ingress`

For roadmap-oriented workflow extensions built on this contract, see **@subpage 17_Workflow_Extensibility**.
For selector-specific contributor hook points, see **@subpage 50_Modular_Selector_Extension_Guide**.
