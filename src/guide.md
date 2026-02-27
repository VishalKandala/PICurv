# Source Guide

This directory contains the C implementation for the solver, postprocessor, and shared runtime subsystems.
It is organized by functional module rather than by executable.

## Runtime Entry Points

- `picsolver.c`: main entry for the flow solver executable.
- `postprocessor.c`: main entry for postprocessing executable.

## Major Module Files

- setup and lifecycle: `setup.c`, `simulation.c`
- flow stepping: `solvers.c`, `momentumsolvers.c`, `rhs.c`, `poisson.c`
- grid and metrics: `grid.c`, `Metric.c`
- BC handling: `Boundaries.c`, `BC_Handlers.c`
- particle system: `ParticleSwarm.c`, `ParticleMotion.c`, `interpolation.c`, `ParticlePhysics.c`
- postprocessing/statistics: `postprocessing_kernels.c`, `particle_statistics.c`, `vtk_io.c`

## Development Notes

1. Shared state is centralized in `SimCtx`/`UserCtx` from `variables.h`.
2. Prefer extending existing module boundaries before adding new files.
3. Keep high-level orchestration in `solvers.c`/`simulation.c`; keep kernels in module files.
4. Update matching headers in `include/` for every public symbol change.

## Related Docs

- `include/guide.md`
- `docs/pages/13_Code_Architecture.md`
- `docs/pages/21_Methods_Overview.md`
