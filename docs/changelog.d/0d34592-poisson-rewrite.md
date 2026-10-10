- The pressure-Poisson solver is built once per run and reused every step instead of
  being rebuilt, including its coarse-level LU factor, which profiling had put at 60-65%
  of each step on a wall-resolved LES duct. Results are bitwise identical to the previous
  solver. The continuity log's `Sum(RHS)` column is now `Poisson Source Imbalance`. The
  previous `src/poisson.c`, with its dormant immersed-boundary routines, is available at
  `53ba654` and earlier.
