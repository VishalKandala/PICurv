# Scheduler Config Guide

This directory stores reusable scheduler profiles for `picurv` cluster runs and parameter sweeps. Scheduler profiles are the bridge between logical PICurv workflow stages and actual batch system execution policy (resources, queueing, job dependencies, and runtime commands).

## What This Config Should Capture

- partition/account settings,
- walltime and resource requests,
- launcher behavior (`srun`, `mpirun`, etc. as expected by site policy),
- output/logging file conventions,
- array/dependency settings used by `picurv sweep`.

## Baseline Profile

- `slurm_default.yml`: conservative default Slurm profile suitable as a starting template.

## How To Use

- Single run:
  - `./bin/picurv run --cluster config/schedulers/slurm_default.yml ...`
- Sweep run:
  - `./bin/picurv sweep --cluster config/schedulers/slurm_default.yml ...`

## Operational Guidance

1. Validate scheduler settings on a tiny dry-run or short queue job before production.
2. Keep cluster profile changes separate from physics/numerics changes when debugging failures.
3. Preserve one site-approved baseline profile and derive specialized profiles from it.

## Reference Schema

- `examples/master_template/master_cluster.yml`
