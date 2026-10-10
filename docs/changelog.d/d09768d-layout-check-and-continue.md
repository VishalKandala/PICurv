- `picurv run` (including `--dry-run`) now refuses a `grid.da_processors_x/y/z` layout that
  would leave the coarsest multigrid level fewer points per rank than the stencil width,
  naming the axis and the most ranks it can take, instead of letting the job abort in grid
  setup after the queue wait. `--dry-run` also catches a layout whose product differs from
  the rank count. `--continue` now accepts a changed rank layout or rank count; only a
  change to the physical case is refused.
