- Newton--Krylov momentum and its unpreconditioned matrix-free path are now `supported`
  within the documented single-block scope. A retained 96 x 96 x 256-cell,
  144-rank channel campaign records 10,000 converged and committed steps, the
  restart statistics window, and an exploratory Lee--Moser DNS comparison.
  The run's inferred `Re_tau` is 202.876 against 182.088 and its skin friction
  is 29.3% higher, so the campaign establishes production execution rather than
  quantitative DNS agreement. The frozen-Jacobian preconditioner, cluster
  scheduling, and workspace asset lifecycle retain their experimental status.
- The physical-units particle smoke fixture uses a grid large enough for the
  three-rank multigrid decomposition; the former 8-cells-per-axis grid left a
  coarse-grid partition narrower than the PETSc stencil.
