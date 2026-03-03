# PICurv Monitor Profiles Library

## 1. Purpose

This directory stores reusable `monitor.yml` profiles for solver observability.

Monitor profiles control:
- logging verbosity and function whitelist,
- profiling timestep-reporting and final-summary behavior,
- file output cadence, particle console cadence, and directories,
- raw PETSc monitor passthrough.

## 2. Included Profiles

- `Standard_Output.yml`: recommended default for most runs.

For full schema coverage, use:
- `examples/master_template/master_monitor.yml`

## 3. How to Use

```bash
./bin/picurv run --solve -n 8 \
  --case my_study/case.yml \
  --solver my_study/solver.yml \
  --monitor config/monitors/Standard_Output.yml
```

## 4. Current Verbosity Levels

`ERROR`, `WARNING`, `INFO`, `DEBUG`, `TRACE`, `VERBOSE`
