- Experimental: `Volume` particle seeding is now uniform per unit physical volume on stretched
  and curvilinear grids. Before, it was uniform per cell, which over-filled small cells. Cells
  are chosen in proportion to their volume, and each rank's share follows its owned volume. On
  uniform grids the distribution is unchanged, but the particle positions are a different random
  realisation. `MetricJacobian` now returns the physical Jacobian; it returned J/8 and had no
  callers.
