- Experimental: `picurv run` takes the MPI decomposition from optional `--ni/--nj/--nk`, from `cluster.yml`
  `execution.decomposition`, or chooses it from the grid, its multigrid levels and periodicity when only `-n` is
  given; the choice is printed before launch and in `--dry-run`, the postprocessor gets its own, and `case.yml`
  `grid.da_processors_x/y/z` is deprecated but still read.
