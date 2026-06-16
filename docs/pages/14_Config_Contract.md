@page 14_Config_Contract Configuration Contract (YAML -> Generated Artifacts -> Runtime)

@anchor _Config_Contract

This page is the user-facing source of truth for the configuration contract implemented by `picurv`.
It describes the launcher-level contract, which may be stricter or more explicit than the raw C defaults because `picurv` validates and normalizes inputs before runtime.

@tableofcontents

@section p14_inputs_sec 1. Required Input Roles

`picurv` composes a standard single-run workflow from five logical inputs, with two additional files for cluster/sweep modes:

1. `case.yml`: physics, grid, BC definitions, run control.
2. `solver.yml`: numerical strategy and solver parameters.
3. `monitor.yml`: I/O, logging/profiling, and diagnostics controls.
4. `post.yml`: post-processing recipe.
5. MPI launch settings (`-n`, executable stage selection).
   - `-n/--num-procs` sizes the solver stage launch.
   - post stage is forced to single-rank execution by workflow policy, even when the solver stage uses more ranks.
   - optional site execution defaults can be supplied via nearest `.picurv-execution.yml`.
6. (Cluster mode) `cluster.yml`: scheduler/resource/launcher contract.
7. (Sweep mode) `study.yml`: parameter matrix + metrics/plot contract.

You can name files however you want. File names are not hardcoded on the C side; `picurv` resolves paths and emits generated artifacts.

These roles are intentionally modular:

- `case.yml` describes physical setup and geometry contract.
- `solver.yml` describes numerical strategy.
- `monitor.yml` describes logging, diagnostics, profiling, and I/O behavior.
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
- `grid.da_processors_x/y/z` optionally set the global DMDA layout for any grid mode.
- Legacy `grid.programmatic_settings.da_processors_*` is still accepted for compatibility.
- `da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block MPI decomposition is not currently supported.
- For `grid_gen`, `grid.generator.config_file` is required today. `grid.gen` consumes cell counts and writes node counts into `.picgrid`.
- For `file`, optional `grid.legacy_conversion` can call `grid.gen legacy1d` to convert headerless 1D-axis legacy payloads before standard validation/non-dimensionalization.
- `boundary_conditions` supports single-block list or multi-block list-of-lists.
- `INLET` + `prescribed_flow` supports `source.type: file`, `source.type: generated`,
  and `source.type: field_slice`. Generated square-duct Poiseuille profiles and
  old-field slices are produced by Python, written under `config/`, summarized
  in `profile.info`, then converted to the existing C-side `source_file` contract.
- `solver_parameters` is an advanced passthrough map for raw flags not yet modeled in schema.
- `properties.initial_conditions.mode` is `generated` or `file`.
- generated built-ins are `zero`, `constant`, `streamwise_constant`, and `poiseuille`; their
  inputs live under `params`.
- `generator: ic_gen` requires `params.field` and `params.config_file`, defaults to
  `generators/ic.gen`, accepts optional `params.script`, and stages its PETSc vector output exactly like file mode.
  The repository generator evaluates `[expression]` configs on a staged nondimensional single-block PICGRID.
- generated and field-sliced `prescribed_flow` sources default to `generators/profile.gen`
  and accept optional `source.script`.
- `mode: file` requires `field: Ucat|Ucont` and `source_file`.
- file-backed ICs currently support single-block cases only.
- `solver.operation_mode.eulerian_field_source` and `run_control.start_step` determine whether
  `initial_conditions` has authority.

@section p14_solver_sec 4. Solver Contract Highlights

- `operation_mode.eulerian_field_source` -> `-euler_field_source`
- `operation_mode.analytical_type` -> `-analytical_type`
- `operation_mode.uniform_flow.{u,v,w}` -> `-analytical_uniform_u/-analytical_uniform_v/-analytical_uniform_w` for `UNIFORM_FLOW`
- `verification.sources.diffusivity.*` -> `-verification_diffusivity_*`
- `verification.sources.scalar.*` -> `-verification_scalar_*`
- `strategy.momentum_solver` -> `-mom_solver_type` via normalized names.
- Solver-specific block support currently includes `momentum_solver.dual_time_picard_jameson_rk`.
- Deprecated `dual_time_picard_rk4` and `rk4_residual_noise_allowance_factor`
  spellings are accepted only as compatibility aliases and normalize to Jameson controls.
- `momentum_solver.dual_time_picard_jameson_rk` controls and their C-side flags:
  - `max_pseudo_steps` -> `-mom_max_pseudo_steps` (default 50): maximum accepted pseudo-time
    iterations per physical step; rejected iterations are counted separately.
  - `absolute_tol` -> `-mom_atol` (default 1e-7): stop when `||ΔU||_∞ < tol`.
  - `relative_tol` -> `-mom_rtol` (default 1e-4): stop when `||ΔU||_∞ / ||ΔU₀||_∞ < tol`.
  - `jameson_residual_noise_allowance_factor` -> `-mom_dt_jameson_residual_norm_noise_allowance_factor`
    (default 1.1, must be ≥ 1): EMA-smoothed residual ratio threshold above which a
    pseudo-time trial is rejected and the pseudo-CFL is reduced. Raise toward 1.2–1.5 for
    convection-dominated or non-monotone convergence histories; lower toward 1.05 for strict
    monotone enforcement.
  - `ratio_ema_alpha` -> `-mom_ratio_ema_alpha` (default 0.3, range [0, 1]): exponential
    moving-average coefficient applied to the step-to-step residual ratio before the rejection
    decision. `smoothed = α × raw + (1−α) × prev`. Setting `α = 1.0` recovers the original
    raw-ratio behavior; `α = 0.3` requires ~3–4 consecutive bad trials to trigger rejection.
  - `pseudo_cfl.*` values are dimensionless Courant numbers (Phase 3+). The solver computes the
    pseudo-time step as `dtau = pseudo_cfl / lambda_max`, where `lambda_max` is the global maximum
    convective spectral radius of the current velocity field. This makes `pseudo_cfl` independent
    of the physical timestep `dt`, grid size, and flow speed. Stable range for the 4-stage Jameson
    RK smoother is ~0–2.83.
  - `pseudo_cfl.initial` -> `-pseudo_cfl` (default 0.5): starting pseudo-CFL (Courant number).
  - `pseudo_cfl.minimum` -> `-min_pseudo_cfl` (default 0.001): floor below which no further
    CFL reduction is attempted; the solver breaks the retry loop at this point.
  - `pseudo_cfl.maximum` -> `-max_pseudo_cfl` (default 2.0): ceiling on CFL growth (stability limit ~2.83).
  - `pseudo_cfl.growth_factor` -> `-pseudo_cfl_growth_factor` (default 1.1, must be ≥ 1):
    factor applied to the pseudo-CFL after a well-converging accepted trial. Set to 1.0 to
    disable CFL growth entirely.
  - `pseudo_cfl.reduction_factor` -> `-pseudo_cfl_reduction_factor` (default 0.75, must be
    in (0,1)): factor applied to the pseudo-CFL after a rejected trial.
- `solution_convergence.*` -> `-solution_convergence_*` for physical solution drift logging.
- `interpolation.method` -> `-interpolation_method`. Defaults to `Trilinear` (direct cell-center, second-order). Set to `CornerAveraged` for the legacy two-stage path.
- `petsc_passthrough_options` remains the escape hatch for advanced PETSc/C flags.
- `scalar_transport.schmidt_number` and `scalar_transport.turbulent_schmidt_number`
  are the structured scalar/Brownian transport controls; do not use passthrough
  for ordinary Schmidt-number tuning.
- `case.yml -> models.physics.turbulence` is the structured turbulence control surface.
  LES uses `les.enabled/model` plus `constant_cs`, `max_cs`, `dynamic_frequency`,
  and `test_filter`; wall functions use `wall_function.enabled/model/roughness_height`.
  Legacy `les: true` remains constant Smagorinsky (`-les 1`), while `les: 2`
  selects dynamic Smagorinsky.

Analytical-mode compatibility rule:

- when `operation_mode.eulerian_field_source: analytical` is selected, `TGV3D` still requires `case.yml -> grid.mode: programmatic_c`.
- `ZERO_FLOW` and `UNIFORM_FLOW` support `case.yml -> grid.mode: programmatic_c` and `case.yml -> grid.mode: file`.
- `grid_gen` remains outside the current documented analytical contract.

Verification-pathway rule:

- `solver.yml -> verification.sources.diffusivity` and `solver.yml -> verification.sources.scalar` are reserved for verification-only source overrides when no cleaner end-to-end path exists.
- they are only valid for analytical solver runs.
- `verification.sources.scalar` prescribes particle `Psi` and drives the runtime diagnostic `logs/scatter_metrics.csv` while leaving ordinary production runs unchanged when absent.
- new verification source overrides belong in `verification_sources.*`, with production call sites kept as thin delegation points.

Solution-convergence rule:

- `solver.yml -> solution_convergence` is an optional override block for the always-on runtime solution-convergence logger.
- the runtime always writes `logs/solution_convergence.log` once per physical timestep, defaulting to `steady_deterministic` mode when the block is omitted.
- it complements, rather than replaces, the inner solver-health logs:
  - `logs/Momentum_Solver_Convergence_History_Block_*.log`
  - `logs/Poisson_Solver_Convergence_History_Block_*.log`
  - `logs/Continuity_Metrics.log`
- supported modes are:
  - `steady_deterministic`
  - `periodic_deterministic`
  - `statistical_steady`
  - `transient`
- `periodic_deterministic.period_steps` is required only when `mode: periodic_deterministic`.
- `statistical_steady.window_steps` is required only when `mode: statistical_steady`.
- nested mode blocks are rejected when they do not match the selected `mode`.

@section p14_monitor_sec 5. Monitor Contract Highlights

- `io.data_output_frequency` -> `-tio`
- `io.particle_console_output_frequency` -> `-particle_console_output_freq` (defaults to `data_output_frequency` when omitted)
- `io.particle_log_interval` -> `-logfreq`
- `io.directories.output/restart/log` -> `-output_dir/-restart_dir/-log_dir`
- `io.directories.eulerian_subdir/particle_subdir` -> `-euler_subdir/-particle_subdir`
- `profiling.timestep_output` -> `profile.run` when `mode: selected`, plus profiling control flags
- `diagnostics.petsc` -> PETSc startup arguments on solver/postprocessor commands
- `diagnostics.runtime_memory_log` -> `-runtime_memory_log_enabled/-runtime_memory_log_file`
- `solver_monitoring.poisson.*` maps readable monitor keys into prefixed Poisson KSP flags.
- `solver_monitoring.petsc_passthrough_options` maps raw PETSc flags directly into control output.

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
- when `execution.launcher` / `execution.launcher_args` are omitted, `picurv` falls back to nearest `.picurv-execution.yml` (`cluster_execution`, then `default_execution`) before using the built-in default `srun`.
- `execution.walltime_guard` optionally tunes the automatic runtime walltime estimator for generated Slurm solver jobs. When omitted, generated solver jobs still use the built-in default policy (`enabled: true`, `warmup_steps: 10`, `multiplier: 2.0`, `min_seconds: 60`, `estimator_alpha: 0.35`).
- `execution.extra_sbatch` supports scheduler-specific pass-through flags.
- `cluster.yml` does not currently define run naming. `picurv` derives `run_id` from `<case_basename>_<timestamp>` and uses that same run ID to name generated scheduler jobs.

Optional shared runtime execution file:

- `picurv init` writes `.picurv-execution.yml` into each new case with inert defaults.
- nearest `.picurv-execution.yml` may define:
  - `default_execution`
  - `local_execution`
  - `cluster_execution`
- local precedence is: `PICURV_MPI_LAUNCHER` -> `MPI_LAUNCHER` -> `.picurv-execution.yml` -> legacy `.picurv-local.yml` -> default `mpiexec`
- cluster precedence is: `cluster.yml.execution` -> `.picurv-execution.yml cluster_execution` -> `.picurv-execution.yml default_execution` -> default `srun`

`picurv run --cluster ...` generates:
- `runs/<run_id>/scheduler/solver.sbatch`
- `runs/<run_id>/scheduler/post.sbatch`
- `runs/<run_id>/scheduler/solver_<jobid>.out/.err` and `post_<jobid>.out/.err` after submission
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
- `studies/<study_id>/scheduler/solver_<array_jobid>_<taskid>.out/.err` and `post_<array_jobid>_<taskid>.out/.err` after submission
- `studies/<study_id>/results/metrics_table.csv`
- `studies/<study_id>/results/plots/*`

@section p14_passthrough_sec 9. Escape Hatches and Defaults

- Escape hatches stay supported:
  - `case.solver_parameters`
  - `solver.petsc_passthrough_options`
  - `monitor.solver_monitoring.petsc_passthrough_options`
- Prefer structured keys first: use `solver.poisson_solver` for pressure-solver
  KSP/MG setup and `monitor.diagnostics.petsc` for PETSc startup diagnostics.
- Default post input/output extension is `dat` unless overridden.

Launcher defaults vs C defaults:

- some omitted keys intentionally preserve C defaults,
- some omitted keys are now rejected because the launcher requires explicit selection,
- some omissions produce warnings because the launcher knows the downstream C fallback is important.

Examples:

- omitting `properties.initial_conditions.mode` is rejected at the launcher level,
- `generator: zero` requires no `params`,
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
