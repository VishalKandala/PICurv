@page 15_Config_Ingestion_Map Developer Ingestion Map

This page maps configuration flow from YAML schema to generated artifacts and C ingestion/use sites.

@tableofcontents

@section pipeline_sec 1. End-to-End Flow

1. `pic.flow` validates YAML inputs.
2. `pic.flow` generates `*.control`, `bcs*.run`, `whitelist.run`, `profile.run`, `post.run`.
3. `CreateSimulationContext()` in `src/setup.c` loads `-control_file` and ingests options.
4. Secondary file parsers ingest BC/post/logging profile inputs.
5. Runtime modules consume parsed values during setup, solve, and postprocess.

@section map_sec 2. Mapping Matrix

| YAML / Source Key | Generated Artifact Key | C Ingestion Site | Primary Runtime Consumers |
| :--- | :--- | :--- | :--- |
| `case.run_control.*` | `-start_step`, `-totalsteps`, `-dt` | `src/setup.c` (`CreateSimulationContext`) | `src/simulation.c`, setup/timestep logic |
| `case.grid.programmatic_settings.im/jm/km/...` | `-im/-jm/-km`, bounds, stretch (`im/jm/km` translated from YAML cell counts to C node counts) | `src/io.c` (`ReadGridGenerationInputs`, `PopulateFinestUserGridResolutionFromOptions`) | `src/grid.c`, analytical geometry setup |
| `case.grid.programmatic_settings.da_processors_*` | `-da_processors_x/y/z` | `src/setup.c` | `src/grid.c` DMDA creation |
| `case.models.domain.*` | `-nblk`, periodic flags | `src/setup.c` | `src/grid.c`, BC setup |
| `case.models.physics.particles.*` | `-numParticles`, `-pinit`, `-particle_restart_mode`, `-psrc_*` | `src/setup.c` | `src/ParticleSwarm.c`, `src/ParticleMotion.c`, statistics kernels |
| `case.boundary_conditions` | `bcs*.run` rows | BC parser path (`src/Boundaries.c` + helpers) | BC handler factory and boundary application |
| `solver.operation_mode.*` | `-euler_field_source`, `-analytical_type` | `src/setup.c` | `src/initialcondition.c`, `src/grid.c`, `src/AnalyticalSolutions.c` |
| `solver.strategy/tolerances/momentum_solver.*` | solver flags (`-mom_*`, pseudo-CFL, etc.) | `src/setup.c` | `src/momentumsolvers.c` |
| `solver.pressure_solver.*` | `-poisson_tol`, `-mg_*`, prefixed PETSc flags | `src/setup.c` + PETSc options db | `src/poisson.c` |
| `solver.petsc_passthrough_options` | raw flags in control | PETSc options db | PETSc KSP/PC stack, mostly in `src/poisson.c` |
| `monitor.io.*` | `-tio`, `-logfreq`, `-output_dir`, `-restart_dir`, `-log_dir`, `-euler_subdir`, `-particle_subdir` | `src/setup.c` | `src/io.c`, `src/setup.c` env setup |
| `monitor.logging.*` | `whitelist.run`, `LOG_LEVEL` env | `src/setup.c` + `src/logging.c` | logging macros/system |
| `monitor.profiling.*` | `profile.run` | `src/setup.c` + profiling init | profiler summaries |
| `monitor.solver_monitoring` | raw flags in control | PETSc options db | PETSc monitors/convergence output |
| `post.run_control.*` | `startTime/endTime/timeStep` in `post.run` | `src/io.c` (`ParsePostProcessingSettings`) | `src/postprocessor.c` main loop |
| `post.io.input_extensions.*` | `eulerianExt`, `particleExt` in `post.run` | `src/io.c` | `ReadSimulationFields`, `ReadAllSwarmFields`, swarm precheck |
| `post.io.eulerian_fields_averaged` | `output_fields_averaged` in `post.run` | `src/io.c` | reserved/no-op in current writer path |
| `post.statistics_pipeline.*` | `statistics_pipeline`, `statistics_output_prefix` | `src/io.c` | `GlobalStatisticsPipeline` dispatch |

@section python_only_sec 3. Python-Only Orchestration Mapping (No C Ingestion)

These keys are consumed by `pic.flow` orchestration only:

| YAML / Source Key | Generated Artifact Key | Ingestion Site | Runtime Consumer |
| :--- | :--- | :--- | :--- |
| `cluster.scheduler.type` | scheduler selector | `scripts/pic.flow` (`validate_cluster_config`) | `run`/`sweep` scheduler dispatch |
| `cluster.resources.*` | `#SBATCH` directives | `scripts/pic.flow` (Slurm renderers) | Slurm scheduler |
| `cluster.notifications.*` | `#SBATCH --mail-*` | `scripts/pic.flow` | Slurm scheduler |
| `cluster.execution.module_setup` | pre-launch shell lines in `*.sbatch` | `scripts/pic.flow` | batch script runtime env |
| `cluster.execution.launcher*` | launch command (`srun`/`mpirun`) | `scripts/pic.flow` | solver/post executable launch |
| `study.base_configs.*` | per-case config materialization | `scripts/pic.flow` (`sweep_workflow`) | case generation pipeline |
| `study.parameters` | case matrix expansion | `scripts/pic.flow` (`expand_parameter_matrix`) | study case synthesis |
| `study.execution.max_concurrent_array_tasks` | Slurm array `%N` throttle | `scripts/pic.flow` | Slurm scheduler |
| `study.metrics` | `metrics_table.csv` extraction contract | `scripts/pic.flow` metric extractors | study aggregation/reporting |
| `study.plotting` | `results/plots/*` output controls | `scripts/pic.flow` plotting pipeline | study reporting |

@section exceptions_sec 4. Important Exceptions

- PETSc runtime option consumption is not only explicit `PetscOptionsGet*`.
- `KSPSetFromOptions` in `src/poisson.c` ingests prefixed PETSc options dynamically.
- `LOG_LEVEL` is environment-driven (`src/logging.c`) and intentionally outside control-file parsing.

@section mode_dependent_sec 5. Mode-Dependent Normalization in `pic.flow`

Some launcher behaviors depend on other config selections before values ever reach C:

- `case.properties.initial_conditions.mode: Zero`
  allows `u_physical`, `v_physical`, and `w_physical` to be omitted; `pic.flow` writes zeros.
- `case.properties.initial_conditions.mode: Poiseuille`
  accepts `peak_velocity_physical` as the preferred YAML input and maps it to the inlet-aligned
  `-ucont_*` component that C later interprets as `Vmax`.
- `solver.operation_mode.eulerian_field_source: analytical`
  currently requires `case.grid.mode: programmatic_c` at launcher validation time because the
  active C ingestion path for analytical runs does not consume file/grid-gen geometry in the same
  way as standard solve mode.
- `case.properties.initial_conditions.mode`
  is now an explicit launcher requirement even though raw C has its own internal default.

This means the YAML contract is intentionally stricter than "whatever C would do with missing
options" in several places.

@section maintenance_sec 6. Drift Prevention

- Use `scripts/audit_ingress.py` to compare PETSc option ingress in `setup.c/io.c` with the maintained manifest.
- Keep this map and the manifest updated whenever new options are introduced.
- Run manually with:
  - `python3 scripts/audit_ingress.py`
  - or `make audit-ingress`

For roadmap-oriented workflow extensions built on this contract, see **@subpage 17_Workflow_Extensibility**.
