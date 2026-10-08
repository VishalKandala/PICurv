
- The frozen-momentum point-block preconditioner for the Newton–Krylov momentum solver
  (`preconditioner.model: frozen_momentum_jacobian`) is now `supported`. Every one of 13,500
  steps converged and committed across a three-grid laminar study of the Humphrey, Taylor &
  Whitelaw (1977) square-duct bend at Re = 790, with divergence below 5e-10. The claim is
  correct execution, not speed: it was not compared against running without a
  preconditioner, and it is untested on wall-resolved, high-Reynolds-number grids. The
  campaign's comparison with the 1977 measurements and its grid-convergence analysis are
  recorded in the Newton–Krylov page.
