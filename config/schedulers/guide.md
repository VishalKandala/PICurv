# Scheduler Config Guide

This directory stores reusable scheduler profiles for `picurv` cluster runs and parameter sweeps. Scheduler profiles are the bridge between logical PICurv workflow stages and actual batch system execution policy (resources, queueing, job dependencies, and runtime commands).

## What This Config Should Capture

- partition/account settings,
- walltime and resource requests,
- batch launcher behavior (`srun`, `mpirun`, etc. as expected by site policy),
- output/logging file conventions,
- array/dependency settings used by `picurv sweep`.

Keep launcher executable and site flags separated when possible:

- prefer `execution.launcher: "srun"` or `"mpirun"`
- put site flags in `execution.launcher_args`
- multi-word `execution.launcher` strings are accepted for site compatibility, but they are not the preferred portable form

If your site wants the same launcher tokens for login-node and batch runs, put those shared defaults in `.picurv-execution.yml` and only override them here when batch jobs need something different.

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
4. `cluster.yml` does not name the run directory; `picurv` generates `runs/<case_basename>_<timestamp>/` automatically and derives Slurm job names from that run ID.
5. For interactive multi-rank runs on login nodes, prefer `.picurv-execution.yml` (or `PICURV_MPI_LAUNCHER` for one-off runs) instead of overloading scheduler profiles.

## Reference Schema

- `examples/master_template/master_cluster.yml`
