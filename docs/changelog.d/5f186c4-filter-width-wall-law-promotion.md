- `filter_width: geometric_mean`, `wall_function.model: log_law` and
  `test_filter.kernel: simpson_ik` are now supported, on the HOM02 benchmark and the
  wall-modelled Re_tau ~ 1000 channel. `wall_function.model: cabot` stays experimental: on
  that channel its pressure-gradient term put the friction velocity 22% low, so use `werner`
  or `log_law` for results.
