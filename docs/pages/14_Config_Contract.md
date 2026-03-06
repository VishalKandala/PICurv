@page 14_Config_Contract Configuration Contract (YAML -> Generated Artifacts -> Runtime)

@anchor _Config_Contract

This page is the user-facing source of truth for the configuration contract implemented by `picurv`.
It describes the launcher-level contract, which may be stricter or more explicit than the raw C defaults because `picurv` validates and normalizes inputs before runtime.

@tableofcontents

@section p14_inputs_sec 1. Required Input Roles

`picurv` composes a standard single-run workflow from five logical inputs, with two additional files for cluster/sweep modes:

1. `case.yml`: physics, grid, BC definitions, run control.
2. `solver.yml`: numerical strategy and solver parameters.
3. `monitor.yml`: I/O and logging/profiling controls.
4. `post.yml`: post-processing recipe.
5. MPI launch settings (`-n`, executable stage selection).
   - `-n/--num-procs` sizes the solver stage launch.
   - post stage defaults to single-rank execution by workflow policy.
6. (Cluster mode) `cluster.yml`: scheduler/resource/launcher contract.
7. (Sweep mode) `study.yml`: parameter matrix + metrics/plot contract.

You can name files however you want. File names are not hardcoded on the C side; `picurv` resolves paths and emits generated artifacts.

These roles are intentionally modular:

- `case.yml` describes physical setup and geometry contract.
- `solver.yml` describes numerical strategy.
- `monitor.yml` describes logging and I/O behavior.
- `post.yml` describes post-processing outputs.

In normal use, you reuse and mix these files instead of cloning one monolithic config for every run.

@section p14_artifacts_sec 2. Generated Artifacts

For each run, `picurv` generates:

- `<run_id>.control`: master PETSc/control flags for solver/post setup.
- `bcs.run` or `bcs_block*.run`: boundary condition definitions.
- `whitelist.run`: logging function allow-list.
- `profile.run`: selected per-step profiling function list (only when `profiling.timestep_output.mode: selected`).
- `post.run`: key=value post-processing recipe consumed by C post parser.

@section p14_case_sec 3. Case Contract Highlights

- `grid.mode` supports: `file`, `programmatic_c`, `grid_gen`.
- For `programmatic_c`, per-block arrays are supported for geometry (`im/jm/km`, bounds, stretching).
- `programmatic_c.im/jm/km` are cell counts in YAML; `picurv` converts them to node counts before writing `-im/-jm/-km`.
- `da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block MPI decomposition is not currently supported.
- For `grid_gen`, `grid.generator.config_file` is required today. `grid.gen` consumes cell counts and writes node counts into `.picgrid`.
- For `file`, optional `grid.legacy_conversion` can call `grid.gen legacy1d` to convert headerless 1D-axis legacy payloads before standard validation/non-dimensionalization.
- `boundary_conditions` supports single-block list or multi-block list-of-lists.
- `solver_parameters` is an advanced passthrough map for raw flags not yet modeled in schema.
- `properties.initial_conditions.mode` is required explicitly by the launcher.
- `properties.initial_conditions.mode: Zero` may omit velocity components.
- `properties.initial_conditions.mode: Poiseuille` supports:
  - `peak_velocity_physical` (scalar centerline speed), or
  - `u_physical/v_physical/w_physical` (explicit component override).

@section p14_solver_sec 4. Solver Contract Highlights

- `operation_mode.eulerian_field_source` -> `-euler_field_source`
- `operation_mode.analytical_type` -> `-analytical_type`
- `strategy.momentum_solver` -> `-mom_solver_type` via normalized names.
- Solver-specific block support currently includes `momentum_solver.dual_time_picard_rk4`.
- `petsc_passthrough_options` remains the escape hatch for advanced PETSc/C flags.

Analytical-mode compatibility rule:

- when `operation_mode.eulerian_field_source: analytical` is selected, the current launcher contract requires `case.yml -> grid.mode: programmatic_c`.
- this reflects the current C analytical ingestion path, which does not consume `file`/`grid_gen` geometry in the standard way.

@section p14_monitor_sec 5. Monitor Contract Highlights

- `io.data_output_frequency` -> `-tio`
- `io.particle_console_output_frequency` -> `-particle_console_output_freq` (defaults to `data_output_frequency` when omitted)
- `io.particle_log_interval` -> `-logfreq`
- `io.directories.output/restart/log` -> `-output_dir/-restart_dir/-log_dir`
- `io.directories.eulerian_subdir/particle_subdir` -> `-euler_subdir/-particle_subdir`
- `solver_monitoring` maps raw flags directly into control output.

@section p14_post_sec 6. Post Contract Highlights

- Pipelines are serialized into semicolon-delimited C pipeline strings.
- `io.eulerian_fields` -> `output_fields_instantaneous`
- `io.eulerian_fields_averaged` -> `output_fields_averaged` (reserved/no-op in current writer path)
- `io.particle_fields` -> `particle_fields_instantaneous`
- `io.input_extensions.eulerian/particle` -> `eulerianExt/particleExt` for post input readers
- `source_data.directory` -> `source_directory`

@section p14_cluster_sec 7. Cluster Contract Highlights (cluster.yml)

- `scheduler.type` currently supports `slurm` only.
- `resources.account/nodes/ntasks_per_node/mem/time` are required.
- `resources.partition` is optional.
- `notifications.mail_user/mail_type` are optional; email is validated when provided.
- `execution.module_setup` injects shell lines before launch.
- `execution.launcher` controls launch style (`srun`, `mpirun`, custom). A multi-word launcher string is accepted for site compatibility, but keeping the executable here and extra flags in `execution.launcher_args` is the preferred portable form.
- `execution.launcher_args` provides site-specific launch flags and is appended after any inline tokens parsed from `execution.launcher`.
- `execution.extra_sbatch` supports scheduler-specific pass-through flags.

`picurv run --cluster ...` generates:
- `runs/<run_id>/scheduler/solver.sbatch`
- `runs/<run_id>/scheduler/post.sbatch`
- `runs/<run_id>/scheduler/submission.json`

@section p14_study_sec 8. Study Contract Highlights (study.yml)

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

`picurv sweep --study ... --cluster ...` generates:
- `studies/<study_id>/scheduler/case_index.tsv`
- `studies/<study_id>/scheduler/solver_array.sbatch`
- `studies/<study_id>/scheduler/post_array.sbatch`
- `studies/<study_id>/results/metrics_table.csv`
- `studies/<study_id>/results/plots/*`

@section p14_passthrough_sec 9. Escape Hatches and Defaults

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
For selector-specific contributor hook points, see **@subpage 50_Modular_Selector_Extension_Guide**.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Configuration Contract (YAML -> Generated Artifacts -> Runtime)** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
