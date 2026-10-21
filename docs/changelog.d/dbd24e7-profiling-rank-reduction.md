- Experimental: `ProfilingSummary_*.log` now reports each function's total time on the slowest
  and fastest rank (with the rank) and the mean over ranks, sorted by the slowest rank; it was rank
  0's timers alone. The per-step profiling file still reports rank 0.
