# Branch Sync Divergence Matrix

This matrix freezes `main` vs `hotfix/geometric-periodic-modern-helpers` divergence during the periodic-isolation sync.

## Classification Rules

- `PERIODIC_CORE`: geometric periodic implementation behavior and direct orchestration logic.
- `PERIODIC_SUPPORT`: periodic-only validation/test-plan content.
- `SHARED_NON_PERIODIC`: tooling, tests, docs, build UX, maintenance flows.

## Commit-Level Classification

| Commit / Area | Classification | Sync Action |
|---|---|---|
| `7ad6bc8` Hotfix geometric periodic orchestration and cluster validation plan | `PERIODIC_CORE` + `PERIODIC_SUPPORT` | Keep hotfix-only |
| `7469459` test: add unified validation suite | `SHARED_NON_PERIODIC` | Sync both branches |
| `ea8144a` test: add temporary local PETSc 3.19 wrapper | `SHARED_NON_PERIODIC` | Sync both branches |
| Installation automation/docs updates (bootstrap + install guide alignment) | `SHARED_NON_PERIODIC` | Sync both branches |

## File-Level Approval Table

| File/Area | Classification | Approved Action |
|---|---|---|
| `include/Boundaries.h`, `src/Boundaries.c`, `src/poisson.c`, periodic-related `scripts/picurv` sections | `PERIODIC_CORE` | Keep isolated to hotfix |
| `docs/GeometricPeriodic_Hotfix_Cluster_TestPlan.md` | `PERIODIC_SUPPORT` | Keep isolated to hotfix |
| `Makefile` doctor/unit/smoke/check targets and test wiring | `SHARED_NON_PERIODIC` | Sync both branches |
| `Makefile.local-petsc19` | `SHARED_NON_PERIODIC` | Sync both branches |
| `tests/c/*`, `tests/smoke/*`, `tests/guide.md` | `SHARED_NON_PERIODIC` | Sync both branches |
| `docs/pages/51_C_Test_Suite_Developer_Guide.md` | `SHARED_NON_PERIODIC` | Sync both branches |
| `README.md`, `docs/pages/01_Installation.md`, `docs/pages/05_The_Conductor_Script.md`, `docs/pages/40_Testing_and_Quality_Guide.md` | `SHARED_NON_PERIODIC` | Sync both branches |
| `scripts/bootstrap_install.sh` | `SHARED_NON_PERIODIC` | Sync both branches |
| `.gitignore` generated-artifact policy (`docs_build/`, `bin/`) | `SHARED_NON_PERIODIC` | Sync both branches |
| `tmp/programmatic_grid_translation_check/*`, `tmp/tiny_grid.picgrid` | `SHARED_NON_PERIODIC` (temporary assets) | Keep isolated to hotfix unless later promoted |

## Outcome Expectation

After sync, `main..hotfix` should differ only by intentionally isolated periodic items and temporary hotfix-only assets listed above.
