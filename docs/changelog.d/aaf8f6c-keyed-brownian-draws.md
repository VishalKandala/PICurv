- Experimental: Brownian increments are drawn from `random_seed`, the particle ID and the step, so a restart with
  diffusion continues the uninterrupted run exactly and the result no longer depends on the rank count. Earlier
  stochastic runs do not reproduce under this build.
