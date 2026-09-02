@page 09_Monitor_Reference Configuration Reference: Monitor YAML

@anchor _Monitor_Reference

For the full commented template, see:

@verbinclude master_template/master_monitor.yml

`monitor.yml` controls runtime monitoring, diagnostics, profiling, run I/O
behavior, and scientific Eulerian field statistics.

@tableofcontents

@section p09_solution_monitoring_sec 1. solution_monitoring

Physical-solution convergence uses the existing runtime monitor through this
contract:

```yaml
solution_monitoring:
  convergence:
    enabled: true
    mode: statistical_steady
    statistical_steady:
      window_steps: 100
```

PICurv writes the normalized settings directly into the generated master
`.control` file as the existing `-solution_convergence_*` options. It retains
`<run.runtime_logs>/solution_convergence.log`, its columns, and its current observable loops.
`enabled: false` disables that writer and its private history allocation. The
default remains enabled `steady_deterministic` monitoring. Convergence is
observed after every completed timestep; there is no user-facing cadence key.

Mode-specific rules:

- `periodic_deterministic` requires
  `periodic_deterministic.period_steps > 0`;
- `statistical_steady` requires `statistical_steady.window_steps > 0`;
- `steady_deterministic` and `transient` take no nested mode block; and
- a nested block for a different mode is rejected.

There is no `observations.run` sidecar and no observation-plan schema/version
at this boundary. `control` remains the single generated C-ingress artifact for
these settings.

Scientific field statistics are configured separately, under `field_statistics`;
see @ref p09_field_statistics_sec. Solution monitoring answers whether the run has
converged, field statistics answer what the converged flow is, and they share no
state.

The removed `case.yml -> models.statistics.time_averaging` and `-averaging`
surface is not compatible input and is not translated into a replacement
workflow.

@section p09_io_sec 2. io

`io` controls output cadence and console sampling. Directory selection is not a
monitor option: every run uses the fixed homes declared by the workspace contract.
Solver state goes under `<run.solver_output>`, checkpoints under `<run.checkpoints>`,
runtime logs under `<run.runtime_logs>`, and restart materializations under
`<run.inputs>/restart`. A `directories` block is rejected as unsupported input.

```yaml
io:
  data_output_frequency: 100
  particle_console_output_frequency: 100
  statistics_console_output_frequency: 100
  particle_log_interval: 10
```

Mappings:
- `data_output_frequency` -> `-tio`
- `particle_console_output_frequency` -> `-particle_console_output_freq`
- `statistics_console_output_frequency` -> `-statistics_console_output_freq`
- `particle_log_interval` -> `-logfreq`

Semantics:
- `data_output_frequency` controls committed-checkpoint cadence. Initial and
  final completed states are also checkpointed, including an off-cadence final
  state.
- Checkpoint-internal `eulerian/`, `particles/`, and `statistics/` names are
  fixed and are not user configuration. See @ref p58_checkpoint_sec.
- `statistics_console_output_frequency` controls how often field-statistics window
  snapshots are printed to the main log. It defaults to `data_output_frequency`, and
  `0` disables the snapshot without disabling accumulation. It is reporting only:
  changing it cannot change any accumulated result.
- `particle_console_output_frequency` controls how often particle snapshots are printed to the main log.
- `particle_log_interval` controls row subsampling within each particle snapshot.
- If `particle_console_output_frequency` is omitted, `picurv` defaults it to `data_output_frequency`.
- If `particle_console_output_frequency` is `0`, periodic particle console snapshots are disabled.

@section p09_field_statistics_sec 3. field_statistics

Accumulates weighted centered moments of Eulerian fields while the solver runs.
Omit the block entirely to accumulate nothing: no storage is allocated and no
statistics output is produced.

```yaml
field_statistics:
  enabled: true
  windows:
    - name: production
      start_time: 50.0
      end_time: 250.0          # optional; open-ended when omitted
      weighting: physical_time # or: sample
      step_cadence: 10         # exactly one of step_cadence / time_cadence
      fields:
        - {field: Ucat, moments: [first, second]}
        - {field: P, moments: [first]}
      covariances:
        - [Ucat, P]
```

Mappings:
- `enabled` -> `-field_statistics_enabled`
- window count -> `-field_statistics_window_count`
- per window `<i>` -> `-field_statistics_window_<i>_name`, `_start_time`,
  `_end_time`, `_weighting`, `_step_cadence` or `_time_cadence`, `_field_count`,
  `_field_<j>_name`, `_field_<j>_moments`, `_covariance_count`, `_covariance_<k>`

`end_time` is omitted from the control file rather than given a sentinel, because
its absence is what makes a window open-ended.

@subsection p09_fs_windows_sub 3.1 Windows

Each window is an independent average. Windows may overlap; each owns its state.
The name identifies a window's saved state across a restart, so it must be unique
within the file and stable across runs.

- `start_time` is when the window begins. The first state at or after it anchors
  the window without being sampled; see @ref p58_bounds_sub.
- `end_time` is optional. A bounded window completes when it passes it and stops
  changing; an open window accumulates until the run ends.
- `weighting` is `sample` (every accepted state counts equally) or
  `physical_time` (a state is weighted by the interval it represents). Use
  `physical_time` whenever the timestep can vary. On a constant timestep the two
  agree exactly.
- Exactly one of `step_cadence` (positive integer) and `time_cadence` (positive
  number) must be set.

@subsection p09_fs_fields_sub 3.2 Fields and covariances

`fields` lists what to accumulate. The first moment is always kept, because every
centered product is measured against it; adding `second` keeps the centered second
moment that Reynolds stresses, RMS, and TKE are derived from.

| Field | Components | Requires |
| --- | --- | --- |
| `Ucat` | 3 | — |
| `P` | 1 | — |
| `Nvert` | 1 | — |
| `Phi` | 1 | — |
| `Psi` | 1 | particles |
| `ParticleCount` | 1 | particles |
| `Nu_t` | 1 | a turbulence model |
| `CS` | 1 | LES |

`covariances` lists cross-field co-moments, such as `[Ucat, P]` for a
pressure-velocity flux. A vector's own Reynolds stresses are **not** a covariance;
they come from `moments: [second]` on that field.

@subsection p09_fs_validation_sub 3.3 Validation

`picurv validate` rejects, naming the offending window:

- a duplicate window name;
- `start_time` greater than or equal to `end_time`;
- neither or both of `step_cadence` and `time_cadence`, or a non-positive value;
- `weighting` outside `{sample, physical_time}`;
- an empty `fields` list, a repeated field, or an empty or unknown `moments` entry;
- a field name that is not accumulable, or whose subsystem is inactive for the run;
- a covariance whose members are not both listed in `fields`, that names one field
  twice, that pairs two vector fields, or that is requested twice.

Combining `--restart-from` with enabled statistics is also rejected: a branch may
follow a different trajectory than the samples already collected, so the windows
would silently restart from zero. Use `--continue` to resume in place, or disable
statistics to branch without them.

@subsection p09_fs_reading_sub 3.4 What a run produces

Accumulated state is written into each committed checkpoint under
`statistics/window_<n>/block_<n>/`, and is restored on `--continue`. Nothing is
derived at solve time; @ref 10_Post_Processing_Reference turns saved state into
Reynolds stresses, RMS, TKE, and fluxes.

While the run is live, `io.statistics_console_output_frequency` prints one line
per window with its state, sample count, accumulated weight, represented time,
progress, and mask coverage. @ref 58_Field_Statistics explains the semantics
behind each of those numbers, and
@ref 60_Field_Statistics_Planned_Extensions records what is not built yet.

@section p09_logging_sec 4. logging

```yaml
logging:
  verbosity: "WARNING"
  enabled_functions: []
```

- `verbosity` maps to environment variable `LOG_LEVEL` via `picurv` launcher.
- `enabled_functions` is serialized into `whitelist.run` only when non-empty.
- If `enabled_functions` is empty, `picurv` omits `whitelist.run` and the C runtime falls back to its default allow-list, which is `main` and `CreateSimulationContext` only.
- Console output is two-tiered by design. **Operator reporting** — the startup banner, progress, graceful-shutdown notices, and the periodic particle and statistics console snapshots — is gated by `verbosity` alone, so a run always reports what it is doing. **Code-level diagnostics** are dual-gated on `verbosity` *and* this allow list, so instrumenting one function for debugging does not drag every other function's output along with it.
- Raising `verbosity` therefore does not surface a function's diagnostics on its own; name the function here as well. See **@subpage 11_User_How_To_Guides** section 3.4 for the workflow.
- An explicitly provided `whitelist.run` must contain at least one function name; an empty whitelist file is invalid.
- `config/monitors/Standard_Output.yml` uses `WARNING` with an empty allow-list for quiet production runs; the startup banner still reports the walltime-guard status.
- Some runtime artifacts are independent of console verbosity. For particle-enabled runs, `<run.runtime_logs>/search_metrics.csv` is written automatically and includes both raw search counters and derived signals such as `search_failure_fraction`, `search_work_index`, and `re_search_fraction`; allow-listing `LOG_SEARCH_METRICS` only affects the optional compact console summary.
- LES runs with `case.yml -> models.physics.turbulence.les.diagnostics.enabled` write `<run.runtime_logs>/les_coefficient.csv`, one row per step or per configured cadence. It carries the effective model coefficient reported as `Cs`, its spatial spread, eddy-viscosity levels, the modelled subgrid kinetic energy, and the volume fractions that were backscattering or limited before clipping. Those values are instantaneous volume statistics, not time averages; window-averaged statistics of the same model's fields come from `field_statistics` instead. Column definitions are at @ref p72_diagnostics_sec, and the history is plottable without leaving the CLI: `picurv summarize --run-dir <run> --plot les.cs_effective`.
- Use **@subpage 53_Search_Robustness_Metrics_Reference** for the exact metric definitions and formulas.

Supported verbosity strings:
- `ERROR`
- `WARNING`
- `INFO`
- `DEBUG`
- `TRACE`
- `VERBOSE`

@section p09_startup_banner_sec Startup banner and solver-specific output

The rank-zero `CASE SUMMARY` is an effective-configuration report, rather than
a dump of every possible option. It is printed independently of `logging.verbosity`
so a quiet production profile still records the run identity, I/O cadence,
walltime-guard status, and active state source.

- Particle cadence, row sampling, restart mode, initialization mode, and
  interpolation are shown only when particles are configured.
- Eulerian source details distinguish a fresh initial condition, a restart, and
  `load` authority.
- `Initial Pseudo-CFL (Courant)` appears only for `Dual Time Picard Jameson RK`.
  It is intentionally absent for `Newton Krylov` and explicit RK, where that
  controller is not active.
- Newton--Krylov iteration detail is controlled by `solver_monitoring.momentum`
  below and its structured solver controls are documented in
  **@subpage 08_Solver_Reference** and **@subpage 55_Newton_Krylov_Momentum_Solver**.
- The banner also records the effective console log level, profiling mode and
  final-summary setting, runtime-memory-log state, and the selected solution
  convergence mode. Periodic and statistical convergence modes additionally
  show their active period or window; those fields are omitted otherwise.
- `Statistics Console Cadence` is always present. It reports the configured
  cadence, or `DISABLED` distinguishing a silenced console over accumulating
  windows from a run that configured none, so a log is never ambiguous about
  whether statistics were being collected.

The C `unit-io` contract test verifies these conditional banner fields, while
`unit-logging` verifies the logging behavior described on this page. Both are
included in the runtime layer of `make certify-docs`; the GitHub documentation
workflow runs the dependency-free structural documentation checks on every main
push and every third day.

The structural certification also runs a user-facing reporting audit. It rejects
raw C console emissions outside the approved banner/logging/bootstrap/progress
surfaces and verifies that the C banner plus CLI and `picurv summarize` retain
their solver-specific reporting branches. The checked-in
`tests/tooling/user_facing_reporting_contract.json` also declares every CLI
command's parser, handler, and context marker, including submission and plot
paths.

@section p09_profiling_sec 5. profiling

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
- `timestep_output.file` sets the filename written under `<run.runtime_logs>/`
- `final_summary.enabled` controls the end-of-run `ProfilingSummary_*.log` file

@section p09_diagnostics_sec 6. diagnostics

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
  run-local defaults such as `<run.runtime_logs>/PETSc_MallocView_Solver.log`,
  `<run.runtime_logs>/PETSc_LogView_Solver.log`, and matching `PostProcessor` files.
- `objects_dump` accepts `false`, `true`, or `all`.
- `options_left` accepts `true`, `false`, or `null`; use `null` to omit the
  PETSc option entirely.
- PETSc diagnostics that support output files use run-local defaults under
  `<run.runtime_logs>/`, with solver/postprocessor-specific filenames. Boolean-only PETSc
  diagnostics remain in the captured solver/post stream logs.
- `runtime_memory_log` writes a rank-reduced, terminal-readable log with max
  process/PETSc allocation signals per step.
- `picurv summarize` reports the latest runtime memory signals when the log is
  present.
- `picurv summarize --list-plot-series` exposes available plottable scalar monitor histories,
  and `--plot <qualified-series>` renders full or last-N append-order histories
  through standalone `generators/plot.gen`.

@section p09_solver_monitoring_sec 7. solver_monitoring

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

@section p09_next_steps_sec 8. Next Steps

Proceed to **@subpage 10_Post_Processing_Reference**.

Also see:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 50_Modular_Selector_Extension_Guide**
- **@subpage 58_Field_Statistics**
