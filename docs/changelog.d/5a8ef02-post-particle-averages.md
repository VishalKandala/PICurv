- Post-processing now loads the saved particle averages (`Psi`, `ParticleCount`) for a case
  with particle `restart_mode: init`; it previously skipped them, so `Psi_nodal` was written
  as zero. The skip now applies only to a solver restart, which reseeds the particles.
