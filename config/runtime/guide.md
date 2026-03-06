# Runtime Config Guide

This directory documents optional runtime-side knobs that are not part of the core case/solver/monitor/post contracts.

## Local MPI Launcher Override

`picurv run` in local mode defaults to `mpiexec` for multi-rank solver launches.

If your site needs a different interactive/login-node launcher, create an optional `.picurv-local.yml` in either:

- the repository root, or
- a case directory (or any parent directory above the case file).

Nearest config wins.

Example:

```yaml
local_execution:
  launcher: "mpirun"
  launcher_args:
    - -mca
    - pml
    - ucx
    - -mca
    - btl
    - ^uct,ofi
    - -mca
    - mtl
    - ^ofi
```

Priority order:

1. `PICURV_MPI_LAUNCHER`
2. `MPI_LAUNCHER`
3. nearest `.picurv-local.yml`
4. built-in default `mpiexec`

Use this only for local multi-rank runs. Cluster batch launches should continue to use `cluster.yml`.
