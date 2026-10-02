
- `picurv init` now organizes a workspace created inside a directory named `runs` or
  `studies`, such as a cluster's `.../runs/` scratch directory. Such a workspace used to keep
  its template files at its root with no `<workspace>/config/case.yml`, although `init` reported success;
  `sync-config` and `status-source`, which render templates the same way, were affected too.
- IEM micromixing of the particle scalar `Psi` is now `supported`. In a zero-flow test from a
  random 0/1 start, the within-cell variance decays at the model rate `exp(-2 Omega t)` to
  within 0.71% over its first e-fold for `C_IEM` of 2, 20 and 200, the whole decay curve
  matches an independent implementation of the same step, and a run continued from a
  checkpoint carries on the same decay. Later, the decay slows where neighbouring cells hold
  different means, because particles crossing faces carry those differences in; this is the
  model's variance production, and it falls as particles per cell increase. The IEM model
  itself, including `C_IEM = 2`, is not yet validated against a turbulent flow. Page 28 now
  covers verification, setting up a mixing case, restart, and troubleshooting.
- The `raw-output` storage retention component is now `supported`, so every part of
  `picurv storage` is supported. It is promoted on the owner's decision: the storage campaign
  runs held no raw output, so it was never exercised with data.
