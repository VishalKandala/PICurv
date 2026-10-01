
- `picurv run --post-process` is incremental: it processes only the requested steps whose
  output is missing or stale (no output, a changed recipe, or a changed checkpoint) and
  launches nothing when every step is up to date. Each produced step is recorded in the
  recipe's `state.json`; output written before these records existed is kept when it is
  newer than its checkpoint. Steps need not be contiguous, so a single missing file is
  reprocessed alone.
  - Post jobs decide their steps when they start, under the post lock, so a Slurm post job
    staged before the solver ran, or resubmitted later, and every study post array task
    process only what is missing then. A job that fails part-way keeps the steps it finished.
  - `--recompute` regenerates every requested step, for every selected stage. Rebuilding the
    postprocessor does not invalidate output; the run reports how many steps another build made.
  - `--continue` now belongs to `--solve`; on a post-only command it is ignored with a warning.
  - Spectra follow the same rule and keep the rows of steps they do not re-measure.
  - Branch runs post-process only the steps from their fork onward.
- The postprocessor no longer crashes on two or more ranks for periodic cases (Q-criterion and
  field-statistics staging); outputs are byte-identical across rank counts.
- A field-statistics window that opens late no longer makes earlier steps count as incomplete.
- An in-place `--solve --continue` from a `start_step` before, or past, the run's last committed
  checkpoint is refused; branch with `--restart-from` to redo a stretch.
- Identical asset inputs reuse one published object, generated initial conditions are
  byte-reproducible, and a restart no longer builds an initial condition it does not use.
- A run or study moved or copied to another directory is pointed at its new location the next
  time PICurv opens it.
- `picurv submit` accepts a post-only staging, submits a study's latest staged set (including
  `sweep --continue --no-submit`), and chains the study's metrics aggregation job.
- Conductor messages use consistent `[WARNING]`/`[FATAL]` prefixes, several stale or misleading
  messages are corrected, and post runs no longer list the solver's PETSc options as unused.
