@page 15_Config_Ingestion_Map Developer Ingestion Map

@anchor _Config_Ingestion_Map

This page maps configuration flow from YAML schema to generated artifacts and C ingestion/use sites.

@tableofcontents

@section p15_pipeline_sec 1. End-to-End Flow

1. `picurv` validates YAML inputs.
2. `picurv` generates `*.control`, `bcs*.run`, and `post.run`, plus optional `whitelist.run` / `profile.run` sidecars when those features are enabled.
3. `CreateSimulationContext()` in `src/setup.c` loads `-control_file` and ingests options.
4. Secondary file parsers ingest BC/post/logging profile inputs.
5. Runtime modules consume parsed values during setup, solve, and postprocess.

@section p15_map_sec 2. Mapping Matrix

| YAML / Source Key | Generated Artifact Key | C Ingestion Site | Primary Runtime Consumers |
| :--- | :--- | :--- | :--- |
| `case.run_control.*` | `-start_step`, `-totalsteps`, `-dt` | `src/setup.c` (`CreateSimulationContext`) | `src/runloop.c`, setup/timestep logic |
| `case.grid.programmatic_settings.im/jm/km/...` | `-im/-jm/-km`, bounds, stretch (`im/jm/km` translated from YAML cell counts to C node counts) | `src/io.c` (`ReadGridGenerationInputs`, `PopulateFinestUserGridResolutionFromOptions`) | `src/grid.c`, analytical geometry setup |
| `case.grid.da_processors_*` | `-da_processors_x/y/z` | `src/setup.c` | `src/grid.c` DMDA creation |
| `case.models.domain.*` | `-nblk`, periodic flags | `src/setup.c` | `src/grid.c`, BC setup |
| `case.models.physics.particles.*` | `-numParticles`, `-pinit`, `-particle_restart_mode`, `-psrc_*` | `src/setup.c` | `src/ParticleSwarm.c`, `src/ParticleMotion.c`, statistics kernels |
| `solver.interpolation.method` | `-interpolation_method` | `src/setup.c` | `src/interpolation.c` (dispatch in @ref InterpolateEulerFieldToSwarm) |
| `case.boundary_conditions` | `bcs*.run` rows | BC parser path (`src/Boundaries.c` + helpers) | BC handler factory and boundary application |
| `solver.operation_mode.*` | `-euler_field_source`, `-analytical_type`, `-analytical_uniform_u/-v/-w` | `src/setup.c` | `src/initialcondition.c`, `src/grid.c`, `src/AnalyticalSolutions.c` |
| `solver.verification.sources.diffusivity.*` | `-verification_diffusivity_mode/-profile/-gamma0/-slope_x` | `src/setup.c` | `src/verification_sources.c`, `src/rhs.c` |
| `solver.strategy/tolerances/momentum_solver.*` | solver flags (`-mom_*`, pseudo-CFL, etc.) | `src/setup.c` | `src/momentumsolvers.c` |
| `solver.pressure_solver.*` | `-poisson_tol`, `-mg_*`, prefixed PETSc flags | `src/setup.c` + PETSc options db | `src/poisson.c` |
| `solver.petsc_passthrough_options` | raw flags in control | PETSc options db | PETSc KSP/PC stack, mostly in `src/poisson.c` |
| `monitor.io.data_output_frequency` | `-tio` | `src/setup.c` | `src/io.c`, `src/setup.c`, `src/runloop.c` |
| `monitor.io.particle_console_output_frequency` | `-particle_console_output_freq` | `src/setup.c` | `src/io.c`, `src/setup.c`, particle console logging |
| `monitor.io.particle_log_interval` | `-logfreq` | `src/setup.c` | particle console row subsampling |
| `monitor.io.directories.*` | `-output_dir`, `-restart_dir`, `-log_dir`, `-euler_subdir`, `-particle_subdir` | `src/setup.c` | `src/io.c`, `src/setup.c`, `src/runloop.c` |
| `monitor.logging.*` | `whitelist.run`, `LOG_LEVEL` env | `src/setup.c` + `src/logging.c` | logging macros/system |
| `monitor.profiling.*` | `profile.run` (selected-mode only) + explicit profiling flags in `*.control` | `src/setup.c` + profiling init | profiler summaries |
| `cluster.execution.walltime_guard.*` | `-walltime_guard_*` in solver `*.control` | `src/setup.c` | `src/runloop.c` runtime walltime estimator and graceful final-write cutoff |
| `monitor.solver_monitoring` | raw flags in control | PETSc options db | PETSc monitors/convergence output |
| `post.run_control.*` | `startTime/endTime/timeStep` in `post.run` | `src/io.c` (`ParsePostProcessingSettings`) | `src/postprocessor.c` main loop |
| `post.io.input_extensions.*` | `eulerianExt`, `particleExt` in `post.run` | `src/io.c` | `ReadSimulationFields`, `ReadAllSwarmFields`, swarm precheck |
| `post.io.eulerian_fields_averaged` | `output_fields_averaged` in `post.run` | `src/io.c` | reserved/no-op in current writer path |
| `post.statistics_pipeline.*` | `statistics_pipeline`, `statistics_output_prefix` | `src/io.c` | `GlobalStatisticsPipeline` dispatch |

@section p15_python_only_sec 3. Python-Only Orchestration Mapping (No C Ingestion)

These keys are consumed by `picurv` orchestration only:

| YAML / Source Key | Generated Artifact Key | Ingestion Site | Runtime Consumer |
| :--- | :--- | :--- | :--- |
| `cluster.scheduler.type` | scheduler selector | `scripts/picurv` (`validate_cluster_config`) | `run`/`sweep` scheduler dispatch |
| `cluster.resources.*` | `#SBATCH` directives | `scripts/picurv` (Slurm renderers) | Slurm scheduler |
| `cluster.notifications.*` | `#SBATCH --mail-*` | `scripts/picurv` | Slurm scheduler |
| `cluster.execution.module_setup` | pre-launch shell lines in `*.sbatch` | `scripts/picurv` | batch script runtime env |
| `cluster.execution.launcher*` | launch command (`srun`/`mpirun`) | `scripts/picurv` | solver launch plus forced single-rank post launch |
| `study.base_configs.*` | per-case config materialization | `scripts/picurv` (`sweep_workflow`) | case generation pipeline |
| `study.parameters` | case matrix expansion | `scripts/picurv` (`expand_parameter_matrix`) | study case synthesis |
| `study.execution.max_concurrent_array_tasks` | Slurm array `%N` throttle | `scripts/picurv` | Slurm scheduler |
| `study.metrics` | `metrics_table.csv` extraction contract | `scripts/picurv` metric extractors | study aggregation/reporting |
| `study.plotting` | `results/plots/*` output controls | `scripts/picurv` plotting pipeline | study reporting |

@section p15_exceptions_sec 4. Important Exceptions

- PETSc runtime option consumption is not only explicit `PetscOptionsGet*`.
- `KSPSetFromOptions` in `src/poisson.c` ingests prefixed PETSc options dynamically.
- `LOG_LEVEL` is environment-driven (`src/logging.c`) and intentionally outside control-file parsing.

@section p15_mode_dependent_sec 5. Mode-Dependent Normalization in `picurv`

Some launcher behaviors depend on other config selections before values ever reach C:

- `case.properties.initial_conditions.mode: Zero`
  allows `u_physical`, `v_physical`, and `w_physical` to be omitted; `picurv` writes zeros.
- `case.properties.initial_conditions.mode: Poiseuille`
  accepts `peak_velocity_physical` as the preferred YAML input and maps it to the inlet-aligned
  `-ucont_*` component that C later interprets as `Vmax`.
- `solver.operation_mode.eulerian_field_source: analytical`
  routes through the analytical grid-ingestion split:
  `TGV3D` still requires `case.grid.mode: programmatic_c`, while `ZERO_FLOW` and `UNIFORM_FLOW`
  accept `case.grid.mode: programmatic_c` or `case.grid.mode: file`.
- `case.properties.initial_conditions.mode`
  is now an explicit launcher requirement even though raw C has its own internal default.

This means the YAML contract is intentionally stricter than "whatever C would do with missing
options" in several places.

@section p15_maintenance_sec 6. Drift Prevention

- Use `scripts/audit_ingress.py` to compare PETSc option ingress in `setup.c/io.c` with the maintained manifest.
- Keep this map and the manifest updated whenever new options are introduced.
- Run manually with:
  - `python3 scripts/audit_ingress.py`
  - or `make audit-ingress`

For roadmap-oriented workflow extensions built on this contract, see **@subpage 17_Workflow_Extensibility**.
For selector-specific contributor hook points, see **@subpage 50_Modular_Selector_Extension_Guide**.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Developer Ingestion Map** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
