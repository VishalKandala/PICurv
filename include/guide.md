# Include Guide

This directory contains public header files for PICurv solver and postprocessor modules.
Headers in this folder define interfaces consumed by translation units in `src/`.

## How To Read This Directory

- `variables.h`: core shared types (`SimCtx`, `UserCtx`, enums, BC structs).
- `setup.h`, `simulation.h`, `solvers.h`: top-level runtime orchestration interfaces.
- `grid.h`, `Metric.h`, `Boundaries.h`, `poisson.h`, `rhs.h`: Eulerian solver subsystems.
- `ParticleSwarm.h`, `ParticleMotion.h`, `interpolation.h`, `ParticlePhysics.h`: Lagrangian and coupling subsystems.
- `postprocessor.h`, `postprocessing_kernels.h`, `particle_statistics.h`: postprocessing/statistics APIs.

## Maintenance Rules

1. Keep function signatures in headers synchronized with definitions in `src/*.c`.
2. Keep Doxygen `@param` names exact; mismatches produce warning noise.
3. Put cross-module shared enums/structs in `variables.h` unless module-local.
4. When adding a new public API, document call order assumptions and ownership semantics.

## Related Docs

- `src/guide.md`
- https://vishalkandala.me/picurv-docs/13_Code_Architecture.html
- https://vishalkandala.me/picurv-docs/35_API_Documentation_Status.html
