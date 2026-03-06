# Include Guide

This directory contains public header files for PICurv solver and postprocessor modules. Headers in this folder define cross-module interfaces and are the contract between translation units in `src/`.

For maintainers, this directory is a stability boundary: changes here can impact many modules, tests, and documentation pages simultaneously.

## How To Read This Directory

- `variables.h`: shared core types (`SimCtx`, `UserCtx`, enums, BC structs).
- `setup.h`, `runloop.h`, `solvers.h`: top-level runtime orchestration APIs.
- `grid.h`, `Metric.h`, `Boundaries.h`, `poisson.h`, `rhs.h`: Eulerian solver subsystem APIs.
- `ParticleSwarm.h`, `ParticleMotion.h`, `interpolation.h`, `ParticlePhysics.h`: Lagrangian/coupling subsystem APIs.
- `postprocessor.h`, `postprocessing_kernels.h`, `particle_statistics.h`: post/statistics interfaces.

## Maintenance Rules

1. Keep function signatures synchronized with definitions in `src/*.c`.
2. Keep Doxygen `@param` names exact to avoid warning noise and stale docs.
3. Place cross-module shared enums/structs in `variables.h` unless module-local.
4. Document call-order assumptions and ownership semantics for every new public API.

## API Change Guidance

- Prefer additive changes over breaking signature rewrites when possible.
- If a signature must change, update all call sites and related docs in one commit.
- Include migration notes in changelog/docs when behavior changes are user-visible.

## Related Docs

- `src/guide.md`
- https://vishalkandala.me/picurv-docs/13_Code_Architecture.html
- https://vishalkandala.me/picurv-docs/35_API_Documentation_Status.html
