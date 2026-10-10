- The dynamic Smagorinsky procedure now weights the `simpson_ik` test filter by the two
  directions it filters: `alpha = width_ratio^(4/3)` instead of `width_ratio^2`. The squared
  ratio had halved the dynamic coefficient (0.091 against the box filter's 0.194 on HOM02).
  Runs with `test_filter.kernel: simpson_ik` get a larger coefficient than before; the
  `volume_weighted_box` path is unchanged.
