# Runtime Config Guide

This directory documents optional runtime-side knobs that are not part of the core case/solver/monitor/post contracts.

## Preferred Shared Execution Config

`picurv run` defaults to `mpiexec` for local multi-rank solver launches and `srun` for generated Slurm jobs.

If your site needs different MPI launcher tokens, create an optional `.picurv-execution.yml` in either:

- the repository root, or
- a case directory (or any parent directory above the case file).

Nearest config wins.

Start from [execution.example.yml](execution.example.yml):

```yaml
default_execution:
  launcher: "mpirun"
  launcher_args:
    - -mca
    - pml
    - ucx

local_execution: {}
cluster_execution: {}
```

Behavior by environment:

1. Local workstation: no file needed; local multi-rank solve uses `mpiexec`.
2. Cluster login node: add `.picurv-execution.yml`; local solve uses `local_execution`, otherwise `default_execution`.
3. Cluster batch job: generated job scripts use `cluster.yml.execution` first, otherwise `cluster_execution`, otherwise `default_execution`.

Local precedence:

1. `PICURV_MPI_LAUNCHER`
2. `MPI_LAUNCHER`
3. nearest `.picurv-execution.yml`
4. nearest legacy `.picurv-local.yml`
5. built-in default `mpiexec`

Cluster precedence:

1. `cluster.yml -> execution.launcher` / `execution.launcher_args`
2. nearest `.picurv-execution.yml -> cluster_execution`
3. nearest `.picurv-execution.yml -> default_execution`
4. built-in default `srun`

## Legacy Local-Only Compatibility

Existing `.picurv-local.yml` files still work for local/login-node runs. They are ignored for cluster job generation.
