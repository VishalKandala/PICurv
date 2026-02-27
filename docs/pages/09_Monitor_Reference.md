@page 09_Monitor_Reference Configuration Reference: monitor.yml

For the full commented template, see:

@verbinclude master_template/master_monitor.yml

`monitor.yml` controls observability and run I/O behavior.

@tableofcontents

@section io_sec 1. io

```yaml
io:
  data_output_frequency: 100
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
- `particle_log_interval` -> `-logfreq`
- `directories.output` -> `-output_dir`
- `directories.restart` -> `-restart_dir`
- `directories.log` -> `-log_dir`
- `directories.eulerian_subdir` -> `-euler_subdir`
- `directories.particle_subdir` -> `-particle_subdir`

@section logging_sec 2. logging

```yaml
logging:
  verbosity: "INFO"
  enabled_functions:
    - AdvanceSimulation
```

- `verbosity` maps to environment variable `LOG_LEVEL` via `pic.flow` launcher.
- `enabled_functions` is serialized into `whitelist.run` and loaded by C logging filters.

Supported verbosity strings:
- `ERROR`
- `WARNING`
- `PROFILE`
- `INFO`
- `DEBUG`
- `TRACE`
- `VERBOSE`

@section profiling_sec 3. profiling

```yaml
profiling:
  critical_functions:
    - Flow_Solver
    - AdvanceSimulation
```

This list is serialized into `profile.run` for critical-path profiling summaries.

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
