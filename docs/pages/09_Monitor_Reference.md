@page 09_Monitor_Reference Configuration Reference: Monitor YAML

For the full commented template, see:

@verbinclude master_template/master_monitor.yml

`monitor.yml` controls observability and run I/O behavior.

@tableofcontents

@section io_sec 1. io

```yaml
io:
  data_output_frequency: 100
  particle_console_output_frequency: 100
  particle_log_interval: 10
  directories:
    output: "results"
    restart: "results"
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

@section logging_sec 2. logging

```yaml
logging:
  verbosity: "INFO"
  enabled_functions:
    - AdvanceSimulation
```

- `verbosity` maps to environment variable `LOG_LEVEL` via `picurv` launcher.
- `enabled_functions` is serialized into `whitelist.run` only when non-empty.
- If `enabled_functions` is empty, `picurv` omits `whitelist.run` and the C runtime falls back to its default allow-list.
- An explicitly provided `whitelist.run` must contain at least one function name; an empty whitelist file is invalid.

Supported verbosity strings:
- `ERROR`
- `WARNING`
- `INFO`
- `DEBUG`
- `TRACE`
- `VERBOSE`

@section profiling_sec 3. profiling

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

@section solver_monitoring_sec 4. solver_monitoring

Raw flag passthrough for PETSc monitors/debug options:

```yaml
solver_monitoring:
  -ps_ksp_pic_monitor_true_residual: true
  -ps_ksp_converged_reason: true
```

Rules:
- Keys must be full flags (with leading `-`).
- `true` emits switch-only flag.
- `false` omits flag.
- Non-boolean values emit `flag value`.

@section next_steps_sec 5. Next Steps

Proceed to **@subpage 10_Post_Processing_Reference**.

Also see:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 50_Modular_Selector_Extension_Guide**
