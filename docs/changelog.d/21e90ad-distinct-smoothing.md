- `poisson_solver.multigrid.pre_sweeps` and `post_sweeps` are now applied separately.
  They were previously collapsed to the larger of the two, with a warning that wrongly
  blamed PETSc. Configurations with equal counts, including every shipped profile, are
  unchanged.
