# PICurv Monitor Profiles Library

This directory stores reusable `monitor.yml` profiles for solver observability. Monitor settings do not change physics, but they strongly affect your ability to diagnose convergence, instability, and runtime regressions.

## 1. What Monitor Profiles Control

- logging verbosity and optional function whitelist filters,
- timestep progress reporting and end-of-run summaries,
- file output cadence and output directory behavior,
- particle console reporting cadence,
- PETSc monitor passthrough options.

## 2. Included Profiles

- `Standard_Output.yml`: recommended default for most runs.

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

Use high verbosity (`DEBUG+`) for short diagnostic runs only; large production runs are easier to manage with `INFO` plus targeted PETSc monitor flags.

## 5. CFD-Oriented Monitoring Advice

- Track continuity and Poisson convergence together; divergence without KSP degradation often points to BC/profile mismatches.
- Use small output intervals during solver tuning, then coarsen intervals for production throughput.
- Keep monitor profiles versioned with studies when publishing benchmark results.
