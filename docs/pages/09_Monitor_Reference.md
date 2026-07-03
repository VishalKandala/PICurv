@page 09_Monitor_Reference Configuration Reference: Monitor YAML

@anchor _Monitor_Reference

For the full commented template, see:

@verbinclude master_template/master_monitor.yml

`monitor.yml` controls observability, diagnostics, and run I/O behavior.

@tableofcontents

@section p09_io_sec 1. io

```yaml
io:
  data_output_frequency: 100
  particle_console_output_frequency: 100
  particle_log_interval: 10
  directories:
    output: "output"
    restart: "restart"
    log: "logs"
    eulerian_subdir: "eulerian"
    particle_subdir: "particles"
```

Mappings:
- `data_output_frequency` -> `-tio`
- `particle_console_output_frequency` -> `-particle_console_output_freq`
- `particle_log_interval` -> `-logfreq`
- `directories.output` -> `-output_dir`
- `directories.restart` -> `-restart_dir`
- `directories.log` -> `-log_dir`
- `directories.eulerian_subdir` -> `-euler_subdir`
- `directories.particle_subdir` -> `-particle_subdir`

Semantics:
- `data_output_frequency` controls file/restart output cadence.
- `particle_console_output_frequency` controls how often particle snapshots are printed to the main log.
- `particle_log_interval` controls row subsampling within each particle snapshot.
- If `particle_console_output_frequency` is omitted, `picurv` defaults it to `data_output_frequency`.
- If `particle_console_output_frequency` is `0`, periodic particle console snapshots are disabled.

@section p09_logging_sec 2. logging

```yaml
logging:
  verbosity: "WARNING"
  enabled_functions: []
```

- `verbosity` maps to environment variable `LOG_LEVEL` via `picurv` launcher.
- `enabled_functions` is serialized into `whitelist.run` only when non-empty.
- If `enabled_functions` is empty, `picurv` omits `whitelist.run` and the C runtime falls back to its default allow-list.
- An explicitly provided `whitelist.run` must contain at least one function name; an empty whitelist file is invalid.
- `config/monitors/Standard_Output.yml` uses `WARNING` with an empty allow-list for quiet production runs; the startup banner still reports the walltime-guard status.
- Some runtime artifacts are independent of console verbosity. For particle-enabled runs, `logs/search_metrics.csv` is written automatically and includes both raw search counters and derived signals such as `search_failure_fraction`, `search_work_index`, and `re_search_fraction`; allow-listing `LOG_SEARCH_METRICS` only affects the optional compact console summary.
- Use **@subpage 53_Search_Robustness_Metrics_Reference** for the exact metric definitions and formulas.

Supported verbosity strings:
- `ERROR`
- `WARNING`
- `INFO`
- `DEBUG`
- `TRACE`
- `VERBOSE`

@section p09_profiling_sec 3. profiling

```yaml
profiling:
  timestep_output:
    mode: "selected"
    functions:
      - Flow_Solver
      - AdvanceSimulation
    file: "Profiling_Timestep_Summary.csv"
  final_summary:
    enabled: true
```

Rules:
- `timestep_output.mode`:
  - `off`: disable per-step profiling output
  - `selected`: write only the listed functions each timestep
  - `all`: write all instrumented functions seen in a timestep
- `timestep_output.functions` is required only when `mode: selected`
- `timestep_output.file` sets the filename written under the run `logs/` directory
- `final_summary.enabled` controls the end-of-run `ProfilingSummary_*.log` file

@section p09_diagnostics_sec 4. diagnostics

Structured diagnostics for PETSc memory/object/function debugging plus a compact
PICurv runtime memory log:

```yaml
diagnostics:
  petsc:
    malloc_debug: false
    malloc_test: false
    malloc_dump: false
    malloc_view: false
    malloc_view_threshold: null
    memory_view: false
    log_view: false
    log_view_memory: false
    log_all: false
    log_trace: false
    objects_dump: false
    options_left: null
  runtime_memory_log:
    enabled: true
    file: "Runtime_Memory.log"
```

Rules:
- PETSc initialization-time diagnostics such as `malloc_debug` and `malloc_test`
  are passed on the executable command line, not only through the generated
  `.control` file.
- For example, this YAML:

```yaml
diagnostics:
  petsc:
    malloc_debug: true
    log_view: true
    objects_dump: "all"
```

adds PETSc startup arguments like:

```text
-malloc_debug -log_view :runs/<run_id>/logs/PETSc_LogView_Solver.log -objects_dump all
```

for the solver stage, with analogous `PostProcessor` log names for post runs.
- `malloc_view`, `log_view`, and `log_trace` accept `false`, `true`, or a
  non-empty PETSc viewer/path string. When set to `true`, PICurv writes
  run-local defaults such as `logs/PETSc_MallocView_Solver.log`,
  `logs/PETSc_LogView_Solver.log`, and matching `PostProcessor` files.
- `objects_dump` accepts `false`, `true`, or `all`.
- `options_left` accepts `true`, `false`, or `null`; use `null` to omit the
  PETSc option entirely.
- PETSc diagnostics that support output files use run-local defaults under
  `logs/`, with solver/postprocessor-specific filenames. Boolean-only PETSc
  diagnostics remain in the captured solver/post stream logs.
- `runtime_memory_log` writes a rank-reduced, terminal-readable log with max
  process/PETSc allocation signals per step.
- `picurv summarize` reports the latest runtime memory signals when the log is
  present.
- `picurv summarize --list-plot-series` exposes available plottable scalar monitor histories,
  and `--plot <qualified-series>` renders full or last-N append-order histories
  through standalone `generators/plot.gen`.

@section p09_solver_monitoring_sec 5. solver_monitoring

Human-readable solver monitor controls. PICurv maps these keys to the raw
C/PETSc flags written into the generated `.control` file:

```yaml
solver_monitoring:
  momentum:
    newton_krylov_history: false
    snes_monitor: true
    snes_converged_reason: true
    ksp_monitor: true
    ksp_converged_reason: true
  poisson:
    pic_true_residual: true
    true_residual: false
    converged_reason: true
    view: false
  petsc_passthrough_options:
    -ps_pc_svd_monitor: false
```

Mappings:
- `momentum.newton_krylov_history` -> `-mom_nk_pic_monitor`
- `momentum.snes_monitor` -> `-mom_nk_snes_monitor`
- `momentum.snes_converged_reason` -> `-mom_nk_snes_converged_reason`
- `momentum.ksp_monitor` -> `-mom_nk_ksp_monitor`
- `momentum.ksp_converged_reason` -> `-mom_nk_ksp_converged_reason`
- `poisson.pic_true_residual` -> `-ps_ksp_pic_monitor_true_residual`
- `poisson.true_residual` -> `-ps_ksp_monitor_true_residual`
- `poisson.converged_reason` -> `-ps_ksp_converged_reason`
- `poisson.view` -> `-ps_ksp_view`

Rules:
- The structured `momentum.*` and `poisson.*` keys are booleans. `true` emits
  a bare switch and `false` emits nothing; no Boolean monitor is serialized as
  a numeric viewer argument.
- `petsc_passthrough_options` remains available for raw PETSc flags not yet
  exposed as structured YAML. In passthrough, `true` emits a switch-only flag,
  `false` omits the flag, and non-boolean values emit `flag value`.
- `solver_monitoring` is written into the generated `.control` file and
  consumed during solver setup, while `diagnostics.petsc` is placed on the
  executable command line for PETSc initialization-time diagnostics.
- Legacy direct flag entries under `solver_monitoring` are still accepted for
  compatibility, but new profiles should prefer the structured form above.

@section p09_next_steps_sec 6. Next Steps

Proceed to **@subpage 10_Post_Processing_Reference**.

Also see:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 50_Modular_Selector_Extension_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Configuration Reference: Monitor YAML** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
