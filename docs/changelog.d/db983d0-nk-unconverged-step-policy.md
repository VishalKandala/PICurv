- The Newton-Krylov momentum solver no longer stops the run whenever SNES fails. A failed
  step whose residual is finite and no larger than at entry is committed as not converged,
  with a warning and `state: committed_unconverged` in the summary log, as the dual-time
  solver already does. A non-finite or grown residual still rolls the step back and stops
  the run. Both families' policy is described on the momentum-solvers page.
