# PICurv Monitor Profiles Library

This directory stores reusable `monitor.yml` profiles for solver observability. Monitor settings do not change physics, but they strongly affect your ability to diagnose convergence, instability, and runtime regressions.

## 1. What Monitor Profiles Control

- logging verbosity and optional function whitelist filters,
- timestep progress reporting and end-of-run summaries,
- file output cadence and output directory behavior,
- particle console reporting cadence,
- structured diagnostics under `diagnostics`:
  - PETSc initialization diagnostics such as malloc/log/object reporting,
  - PICurv's compact `Runtime_Memory.log`,
- PETSc monitor passthrough options.

## 2. Included Profiles

- `Standard_Output.yml`: recommended default for most runs.
- `examples/search_robustness/Search_Robustness_Output.yml`: search-focused monitor profile that keeps runtime quiet by default while preserving `search_metrics.csv`.

For exhaustive schema coverage, see:
- `examples/master_template/master_monitor.yml`

## 3. Typical Usage

```bash
./bin/picurv run --solve -n 8 \
  --case my_study/case.yml \
  --solver my_study/solver.yml \
  --monitor config/monitors/Standard_Output.yml
```

## 4. Verbosity Levels

`ERROR`, `WARNING`, `INFO`, `DEBUG`, `TRACE`, `VERBOSE`

`Standard_Output.yml` now uses `WARNING` plus an empty function allow-list as the
production default. That keeps solver logs quiet while still surfacing runtime
warnings, graceful-shutdown notices, and the startup banner summary.

Use `INFO` when you want more live progress detail, and reserve high verbosity
(`DEBUG+`) for short diagnostic runs only.

For particle-enabled runs, note that some runtime artifacts are written
independently of console verbosity. In particular, `logs/search_metrics.csv` is
always emitted by the runtime when particle tracking is active; allow-listing
`LOG_SEARCH_METRICS` only controls the optional compact console summary.

## 5. Diagnostics and PETSc Monitoring

Use structured `diagnostics.petsc` for PETSc options that must be visible at
`PetscInitialize()` time, such as `malloc_debug`, `malloc_view`, `log_view`,
`log_trace`, `objects_dump`, and `options_left`. Use
`diagnostics.runtime_memory_log` for PICurv's rank-reduced runtime memory log.

Use `solver_monitoring.momentum.newton_krylov_history` to enable PICurv's
optional one-row-per-Newton-iteration nonlinear norm history. The compact
one-row-per-timestep Newton summary is always written, independent of console
verbosity and function allow-lists. The same `momentum` block exposes Boolean
`snes_monitor`, `snes_converged_reason`, `ksp_monitor`, and
`ksp_converged_reason` switches; `true` emits a bare PETSc flag and `false`
emits nothing. Use `solver_monitoring.poisson` for readable
solver-time Poisson monitor controls. PICurv maps keys such as
`pic_true_residual`, `converged_reason`, and `view` to the prefixed PETSc/C flags
consumed from the generated control file. Keep raw one-off PETSc flags under
`solver_monitoring.petsc_passthrough_options`.

## 6. CFD-Oriented Monitoring Advice

- Track continuity and Poisson convergence together; divergence without KSP degradation often points to BC/profile mismatches.
- Use small output intervals during solver tuning, then coarsen intervals for production throughput.
- Keep monitor profiles versioned with studies when publishing benchmark results.
