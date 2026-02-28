@page 14_Config_Contract Configuration Contract (YAML -> Generated Artifacts -> Runtime)

This page is the user-facing source of truth for the configuration contract implemented by `pic.flow`.
It describes the launcher-level contract, which may be stricter or more explicit than the raw C defaults because `pic.flow` validates and normalizes inputs before runtime.

@tableofcontents

@section inputs_sec 1. Required Input Roles

`pic.flow` composes a standard single-run workflow from five logical inputs, with two additional files for cluster/sweep modes:

1. `case.yml`: physics, grid, BC definitions, run control.
2. `solver.yml`: numerical strategy and solver parameters.
3. `monitor.yml`: I/O and logging/profiling controls.
4. `post.yml`: post-processing recipe.
5. MPI launch settings (`-n`, executable stage selection).
   - `-n/--num-procs` sizes the solver stage launch.
   - post stage defaults to single-rank execution by workflow policy.
6. (Cluster mode) `cluster.yml`: scheduler/resource/launcher contract.
7. (Sweep mode) `study.yml`: parameter matrix + metrics/plot contract.

You can name files however you want. File names are not hardcoded on the C side; `pic.flow` resolves paths and emits generated artifacts.

These roles are intentionally modular:

- `case.yml` describes physical setup and geometry contract.
- `solver.yml` describes numerical strategy.
- `monitor.yml` describes logging and I/O behavior.
- `post.yml` describes post-processing outputs.

In normal use, you reuse and mix these files instead of cloning one monolithic config for every run.

@section artifacts_sec 2. Generated Artifacts

For each run, `pic.flow` generates:

- `<run_id>.control`: master PETSc/control flags for solver/post setup.
- `bcs.run` or `bcs_block*.run`: boundary condition definitions.
- `whitelist.run`: logging function allow-list.
- `profile.run`: profiling critical-function list.
- `post.run`: key=value post-processing recipe consumed by C post parser.

@section case_sec 3. Case Contract Highlights

- `grid.mode` supports: `file`, `programmatic_c`, `grid_gen`.
- For `programmatic_c`, per-block arrays are supported for geometry (`im/jm/km`, bounds, stretching).
- `da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block MPI decomposition is not currently supported.
- `boundary_conditions` supports single-block list or multi-block list-of-lists.
- `solver_parameters` is an advanced passthrough map for raw flags not yet modeled in schema.
- `properties.initial_conditions.mode` is required explicitly by the launcher.
- `properties.initial_conditions.mode: Zero` may omit velocity components.
- `properties.initial_conditions.mode: Poiseuille` supports:
  - `peak_velocity_physical` (preferred scalar centerline speed), or
  - `u_physical/v_physical/w_physical` (legacy explicit component form).

@section solver_sec 4. Solver Contract Highlights

- `operation_mode.eulerian_field_source` -> `-euler_field_source`
- `operation_mode.analytical_type` -> `-analytical_type`
- `strategy.momentum_solver` -> `-mom_solver_type` via normalized names.
- Solver-specific block support currently includes `momentum_solver.dual_time_picard_rk4`.
- `petsc_passthrough_options` remains the escape hatch for advanced PETSc/C flags.

Analytical-mode compatibility rule:

- when `operation_mode.eulerian_field_source: analytical` is selected, the current launcher contract requires `case.yml -> grid.mode: programmatic_c`.
- this reflects the current C analytical ingestion path, which does not consume `file`/`grid_gen` geometry in the standard way.

@section monitor_sec 5. Monitor Contract Highlights

- `io.data_output_frequency` -> `-tio`
- `io.particle_log_interval` -> `-logfreq`
- `io.directories.output/restart/log` -> `-output_dir/-restart_dir/-log_dir`
- `io.directories.eulerian_subdir/particle_subdir` -> `-euler_subdir/-particle_subdir`
- `solver_monitoring` maps raw flags directly into control output.

@section post_sec 6. Post Contract Highlights

- Pipelines are serialized into semicolon-delimited C pipeline strings.
- `io.eulerian_fields` -> `output_fields_instantaneous`
- `io.eulerian_fields_averaged` -> `output_fields_averaged` (reserved/no-op in current writer path)
- `io.particle_fields` -> `particle_fields_instantaneous`
- `io.input_extensions.eulerian/particle` -> `eulerianExt/particleExt` for post input readers
- `source_data.directory` -> `source_directory`

@section cluster_sec 7. Cluster Contract Highlights (cluster.yml)

- `scheduler.type` currently supports `slurm` only.
- `resources.account/nodes/ntasks_per_node/mem/time` are required.
- `resources.partition` is optional.
- `notifications.mail_user/mail_type` are optional; email is validated when provided.
- `execution.module_setup` injects shell lines before launch.
- `execution.launcher` controls launch style (`srun`, `mpirun`, custom).
- `execution.launcher_args` provides site-specific launch flags.
- `execution.extra_sbatch` supports scheduler-specific pass-through flags.

`pic.flow run --cluster ...` generates:
- `runs/<run_id>/scheduler/solver.sbatch`
- `runs/<run_id>/scheduler/post.sbatch`
- `runs/<run_id>/scheduler/submission.json`

@section study_sec 8. Study Contract Highlights (study.yml)

- `base_configs` provides case/solver/monitor/post template paths.
- `study_type` is one of:
  - `grid_independence`
  - `timestep_independence`
  - `sensitivity`
- `parameters` defines cartesian-product sweeps using keys of form:
  - `case.<yaml.path>`
  - `solver.<yaml.path>`
  - `monitor.<yaml.path>`
  - `post.<yaml.path>`
- `metrics` defines CSV/log extractors for aggregate tables.
- `plotting` controls whether plots are generated and output format.
- `execution.max_concurrent_array_tasks` maps to Slurm array throttling `%N`.

`pic.flow sweep --study ... --cluster ...` generates:
- `studies/<study_id>/scheduler/case_index.tsv`
- `studies/<study_id>/scheduler/solver_array.sbatch`
- `studies/<study_id>/scheduler/post_array.sbatch`
- `studies/<study_id>/results/metrics_table.csv`
- `studies/<study_id>/results/plots/*`

@section passthrough_sec 9. Escape Hatches and Defaults

- Escape hatches stay supported:
  - `case.solver_parameters`
  - `solver.petsc_passthrough_options`
  - `monitor.solver_monitoring`
- Default post input/output extension is `dat` unless overridden.

Launcher defaults vs C defaults:

- some omitted keys intentionally preserve C defaults,
- some omitted keys are now rejected because the launcher requires explicit selection,
- some omissions produce warnings because the launcher knows the downstream C fallback is important.

Examples:

- omitting `properties.initial_conditions.mode` is rejected at the launcher level,
- omitting velocity components for `mode: Zero` is accepted and defaults to zero,
- omitting `models.physics.particles.restart_mode` on a particle restart emits a warning that C will default to `load`.

For workflow growth patterns (grid generation orchestration, multi-run studies, and ML coupling paths), see **@subpage 17_Workflow_Extensibility**.
For worked examples and profile-composition patterns, see **@subpage 49_Workflow_Recipes_and_Config_Cookbook**.
