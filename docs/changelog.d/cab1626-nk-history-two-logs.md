- `solver_monitoring.momentum.newton_krylov_history` now writes two rank-zero files, and
  said so in only one of the three places a reader looks. The master monitor template
  still described it as the "nonlinear history", and the monitor reference listed the
  key-to-flag mapping without saying what the switch produces; only the Newton-Krylov
  guide named both files. Both now state that it writes one row per Newton iteration and
  one row per Krylov iteration, what the per-Krylov row is for - separating an
  Eisenstat-Walker tolerance change from the linear solve degrading - and that it is
  flushed every iteration, so it is a switch for diagnosing a solve rather than one to
  leave on for a production run.
