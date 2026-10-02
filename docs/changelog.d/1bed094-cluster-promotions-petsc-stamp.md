
- Cluster scheduling, parameter studies, and the workspace asset store are now `supported`,
  as are all three study types (`grid_independence`, `timestep_independence`, `sensitivity`)
  and all four input import modes (`copy`, `hardlink`, `reference`, `reflink`). On a Slurm
  cluster, `picurv cancel` stopped a running solver, `cancel --graceful` made it write an
  off-cadence checkpoint that `run --continue` resumed from, post jobs chose their steps when
  they started, array studies of every type honoured `max_concurrent_array_tasks` and
  aggregated a row per member, and each import mode behaved as documented on Lustre. A
  successful `reflink` has not been observed; the filesystems tested refuse it.
- `sweep --continue` works for studies inside a workspace. It looked for the study's base
  configurations in the workspace rather than in the study, so it failed for every such study.
- `submit --stage post-process` no longer waits on a solve job that has already finished.
  Slurm forgets a finished job and then rejects any dependency on it; a completed solve is
  now not waited for, and a failed one is refused unless `--force` is given.
- `simulator --version` and `postprocessor --version` now also name the PETSc they were built
  against (version, debug or optimized, arch, and directory). `picurv version` shows it, even
  for an executable that cannot start, together with the `libpetsc` the current shell would
  load, and marks a mismatch such as a debug build run against an optimized PETSc
  (`undefined symbol: petscstack`). Each run's software lock records it, and a job now stops
  if its executable was rebuilt against a different PETSc after staging.
- Study metric tables and plots are documented at their real location,
  `studies/<study_id>/output/analysis/`. New troubleshooting entries cover an `sbatch` that
  rejects a job while exiting 0, and the `petscstack` error.
