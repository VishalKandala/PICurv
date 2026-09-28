- `models.physics.particles.fields` gives particle-carried fields an initial value - today
  `Psi` - as a number, an expression in physical `x y z`, normalized `xn yn zn`, `pid`, `t`
  and per-particle `uniform()`/`normal()` draws, or regions (`half_space`, `slab`, `box`,
  `ball`, `cylinder`, with optional smoothed edges) painted over a background. Values are
  applied before the first scatter and summarized in `particle_initial_fields.csv`. The
  expression language is shared with `ic_gen`, which gains `and`, `or`, `not`, and chained
  comparisons.
