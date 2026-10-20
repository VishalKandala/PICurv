- Experimental: a particle that leaves the domain through a removing boundary face now ends
  its search as `EXITED` and is reported separately from a search failure.
  - `search_metrics.csv` and `Particle_Metrics.log` add `exited`, `exited_cumulative` and
    `exited_psi_sum`. `lost` and `lost_psi_sum` now count only search failures.
  - `picurv summarize` reports the exited counts.
  - Settlement no longer searches a lost or exited particle again in later passes. Before,
    each loss was counted once per pass, and on curved boundaries a departed particle could
    migrate between ranks until the run aborted with "Particle migration failed to converge
    after 50 passes".
  - A walk from a particle's previous cell that ends more than three cells beyond the
    boundary is treated as blocked and searched again from scratch in the same pass.
