
- Supported (no longer experimental): the `point_variance` field-statistics output, verified against the exact
  within-cell variance of a known particle ensemble; Newton-Krylov momentum in 2D mode, matching Picard through a
  mixing-layer roll-up on the cluster with every step converged.
- Supported (no longer experimental), verified on the cluster: the MPI decomposition from `--ni/--nj/--nk`,
  `cluster.yml` `execution.decomposition` or the automatic choice, through dry run, `--no-submit` and `submit`; the
  scalar-transport keys in `case.yml` and the LES diagnostics in `monitor.yml`, with their deprecated locations; and
  Newton-Krylov momentum monitors passed only to a Newton-Krylov solve.
