- `monitor.yml` function lists (`logging.enabled_functions`,
  `profiling.timestep_output.functions`) that name a Poisson function the rewrite renamed,
  such as `PoissonSolver_MG` or `Projection`, now translate it to the current name with a
  warning instead of selecting nothing.
