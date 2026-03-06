@page 51_C_Test_Suite_Developer_Guide C Test Suite Developer Guide

This page documents how the PETSc-backed C testing layer is structured and how to extend it safely.

@tableofcontents

@section layout_sec 1. Test Layout

Key directories:

- `tests/c/`: C test sources and shared fixture helpers
- `tests/smoke/`: executable smoke runner for help/init/dry-run/restart workflow checks
- `tests/fixtures/`: Python/control-plane YAML fixtures

Current C files:

- `test_support.h`
- `test_support.c`
- `test_install_check.c`
- `test_geometry.c`
- `test_solver_kernels.c`
- `test_particle_kernels.c`
- `test_io.c`
- `test_postprocessing.c`
- `test_grid.c`
- `test_metric.c`
- `test_boundaries.c`
- `test_poisson_rhs.c`

@section targets_sec 2. Canonical Targets

Canonical user-facing targets:

- `make doctor`
- `make unit`
- `make unit-geometry`
- `make unit-solver`
- `make unit-particles`
- `make unit-io`
- `make unit-post`
- `make unit-grid`
- `make unit-metric`
- `make unit-boundaries`
- `make unit-poisson-rhs`
- `make smoke`
- `make check`

Compatibility aliases (`install-check`, `ctest*`) should be preserved for continuity, but documentation and new contributor guidance should always prefer the canonical names.

@section choosing_sec 3. Where A Test Belongs

Use this rule set:

- `doctor`:
  - only for installation/toolchain/PETSc viability checks
- `unit-*`:
  - for isolated functions, kernels, and small component seams
- `smoke`:
  - for real executable launch or tiny end-to-end entrypoint checks
- `test-python`:
  - for `scripts/picurv`, schema, and repository metadata behavior

If a test can run against a small local fixture without invoking the full executables, it almost always belongs in `unit-*`.

@section fixtures_sec 4. Fixture Helpers

`tests/c/test_support.*` provides the shared PETSc-backed fixture layer.

Current responsibilities:

- create minimal `SimCtx` and `UserCtx` fixtures
- create deterministic `DMDA` objects
- create deterministic `DMSwarm` objects
- seed identity metric fields
- create unique temporary directories under `/tmp`
- provide assertion helpers for reals, ints, vectors, and file existence

Future C tests should reuse these helpers instead of rebuilding fixtures inline unless a suite truly needs a specialized setup.

@section style_sec 5. Conventions

Naming:

- keep file names aligned with the target name (`test_io.c` -> `make unit-io`)
- use function names that state the behavior under test

Assertions:

- fail fast
- print the failing context plus expected and actual values

Tolerances:

- `1e-12` for algebraic invariants
- `1e-10` to `1e-9` for interpolation and I/O round-trips
- exact equality for integer/state checks

@section fs_sec 6. Temporary Files

Rules:

- write temporary artifacts under `/tmp`
- never write repo-tracked test outputs
- keep paths easy to inspect (`/tmp/picurv-test-*`) so a failed test can be debugged manually

@section extending_sec 7. Adding A New C Unit Suite

Recommended workflow:

1. add or reuse a `tests/c/test_<area>.c` file
2. build on `test_support.*`
3. add a dedicated executable in the `Makefile`
4. expose a canonical `unit-<area>` target
5. add a compatibility alias only if it matches an existing naming family
6. document the new target in the README and testing guide if it is user-facing

@section refs_sec 8. Related Pages

- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 01_Installation**
