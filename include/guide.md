# Include Guide

This directory contains public header files for PICurv solver and postprocessor modules. Headers in this folder define cross-module interfaces and are the contract between translation units in `src/`.

For maintainers, this directory is a stability boundary: changes here can impact many modules, tests, and documentation pages simultaneously.

## How To Read This Directory

- `variables.h`: shared core types (`SimCtx`, `UserCtx`, enums, BC structs).
- `field_catalog.h`: typed identities, immutable layout/storage metadata, and
  non-owning runtime views for persistent Eulerian fields.
- `particle_field_catalog.h`: separate DMSwarm field identities, storage types,
  registration ownership, initialization, and Eulerian-coupling metadata.
- `setup.h`, `runloop.h`, `solvers.h`: top-level runtime orchestration APIs.
- `grid.h`, `Metric.h`, `Boundaries.h`, `poisson.h`, `rhs.h`: Eulerian solver subsystem APIs.
- `ParticleSwarm.h`, `ParticleMotion.h`, `interpolation.h`, `ParticlePhysics.h`: Lagrangian/coupling subsystem APIs.
- `postprocessor.h`, `postprocessing_kernels.h`, `particle_statistics.h`: post/statistics interfaces.
- `statistics_moments.h`, `statistics_window.h`, `statistics_accumulator.h`,
  `statistics_target.h`, `statistics_config.h`: Eulerian field-statistics kernels,
  window lifecycle, per-window storage, spatial targeting, and control ingress.

## Maintenance Rules

1. Keep function signatures synchronized with definitions in `src/*.c`.
2. Keep Doxygen `@param` names exact to avoid warning noise and stale docs.
3. Place broad simulation state in `variables.h`; extend `field_catalog.h` for
   persistent Eulerian metadata and `particle_field_catalog.h` for persistent
   solver-particle metadata.
4. Document call-order assumptions and ownership semantics for every new public API.
5. Keep `tests/tooling/audit_function_docs.py` passing when changing public headers or their matching C implementations.

## API Change Guidance

- Prefer additive changes over breaking signature rewrites when possible.
- If a signature must change, update all call sites and related docs in one commit.
- Include migration notes in changelog/docs when behavior changes are user-visible.

## Related Docs

- `src/guide.md`
- https://vishalkandala.me/docs/picurv/13_Code_Architecture.html
- https://vishalkandala.me/docs/picurv/35_API_Documentation_Status.html
