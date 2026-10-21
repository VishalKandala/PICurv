- `picurv validate` now refuses a `grid.da_processors_x/y/z` layout that leaves a rank fewer
  points than the stencil width at the coarsest multigrid level, as `picurv run --dry-run` does; with
  `--cluster` it also checks the rank product against the profile's task count.
