# Tests Guide

PICurv now exposes a layered local testing model:

- `make test-python`: Python-only CLI/config regression suite
- `make test`: backward-compatible alias to `test-python`
- `make coverage-python`: dependency-free Python line-coverage gate for core runtime scripts
- `make coverage-c`: gcov-backed C line-coverage gate (`unit + smoke` with coverage flags)
- `make coverage`: runs Python and C coverage gates
- `make doctor`: installation and PETSc provisioning validation
- `make unit`: all isolated C unit/component suites
- `make unit-geometry`
- `make unit-solver`
- `make unit-particles`
- `make unit-io`
- `make unit-post`
- `make unit-grid`
- `make unit-metric`
- `make unit-boundaries`
- `make unit-poisson-rhs`
- `make unit-runtime`
- `make unit-mpi`: dedicated multi-rank MPI consistency tests (default 2 ranks)
- `make smoke`: executable end-to-end smoke checks (template matrix + tiny runtime sequences, including flat/bent/brownian)
- `make smoke-mpi`: multi-rank runtime smoke checks for tiny flat+bent solve/post plus flat particle+restart workflows
- `make smoke-mpi-matrix`: multi-rank runtime smoke checks across a rank matrix (`SMOKE_MPI_MATRIX_NPROCS`)
- `make check`: full local validation sweep
- `make check-mpi`: `make check` plus multi-rank MPI tests
- `make check-mpi-matrix`: `make check` plus rank-matrix MPI smoke and `unit-mpi`

## Layout

- `test_cli_smoke.py`: CLI help/validate/dry-run smoke coverage for `scripts/picurv`
- `test_case_maintenance.py`: case-origin and sync/build maintenance regressions
- `test_config_regressions.py`: ingress/schema drift guards
- `test_repo_consistency.py`: example validation and repository-wide consistency checks
- `tests/c/`: PETSc-backed C unit binaries used by `make doctor` and `make unit-*`
- `tests/smoke/`: executable smoke runner for template matrix init/validate/dry-run plus tiny end-to-end solve/post/restart/analytical workflows
- `fixtures/valid/`: canonical valid YAML input sets
- `fixtures/invalid/`: intentionally broken YAML input sets

## Canonical Commands

Python-only:

```bash
make test
make coverage-python
```

Installation validation:

```bash
make doctor
```

Subsystem-only C tests:

```bash
make unit-io
make unit-particles
make unit-runtime
make unit-mpi
```

Everything:

```bash
make check
```

Coverage + exhaustive-gate checks:

```bash
make coverage
```

## Notes

- The Python suite does not require PETSc.
- `doctor`, `unit-*`, `unit`, `smoke`, `smoke-mpi`, `smoke-mpi-matrix`, `check`, `check-mpi`, `check-mpi-matrix`, and `coverage-c` assume a working PETSc/MPI toolchain.
- Compatibility aliases such as `make install-check` and `make ctest-*` remain available, but the canonical user-facing names are `doctor` and `unit-*`.
- `make check` is the gate command for pre-merge confidence; `make unit-<area>` commands are the development-loop commands for isolated subsystem work.
- `make check-mpi` is the extended gate command when MPI multi-rank behavior is in scope.
- `make coverage` is the line-coverage gate for the core Python runtime scripts and C solver sources.

## Authoritative Docs

- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/picurv-docs/51_C_Test_Suite_Developer_Guide.html
