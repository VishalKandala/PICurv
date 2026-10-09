- The LES subgrid models and the Werner–Wengle wall model are now `supported`.
  `werner` reproduces Lee & Moser's Re_tau = 1000 channel within criteria fixed before the
  run (friction velocity 2.9% high, mean velocity 1.1% RMS). All four LES models
  (`constant_smagorinsky`, `dynamic_smagorinsky`, `vreman`, `wale`), the `local`,
  `homogeneous` and `global` averaging modes, the `clamp`, `clip_negative` and `none`
  clipping modes, the `volume_weighted_box` test filter, and the `cube_root_volume`,
  `max_edge` and `scotti` filter widths were characterized on the AGARD HOM02
  decaying-isotropic benchmark: every model removes the grid-cutoff energy pile-up, and the
  dynamic coefficient settles at 0.18-0.19. Their accuracy is characterized, not validated:
  no model keeps the resolved energy decay within 10% of the DNS during start-up, and
  Vreman's spectrum misses its criterion narrowly. A wall-modelled LES of the turbulent
  Re_D = 40,000 duct bend ran to completion with them. `simpson_ik`, `geometric_mean`,
  `log_law` and `cabot` remain experimental.
