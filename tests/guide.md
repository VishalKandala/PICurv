# Tests Guide

PICurv now exposes a layered local testing model:

- `make test-python`: Python-only CLI/config regression suite
- `make test`: backward-compatible alias to `test-python`
- `make doctor`: installation and PETSc provisioning validation
- `make unit`: all isolated C unit/component suites
- `make unit-geometry`
- `make unit-solver`
- `make unit-particles`
- `make unit-io`
- `make unit-post`
- `make smoke`: executable entrypoint smoke checks for `bin/simulator` and `bin/postprocessor`
- `make check`: full local validation sweep

## Layout

- `test_cli_smoke.py`: CLI help/validate/dry-run smoke coverage for `scripts/picurv`
- `test_case_maintenance.py`: case-origin and sync/build maintenance regressions
- `test_config_regressions.py`: ingress/schema drift guards
- `test_repo_consistency.py`: example validation and repository-wide consistency checks
- `tests/c/`: PETSc-backed C unit binaries used by `make doctor` and `make unit-*`
- `tests/smoke/`: lightweight executable smoke runner and reserved fixture directories
- `fixtures/valid/`: canonical valid YAML input sets
- `fixtures/invalid/`: intentionally broken YAML input sets

## Canonical Commands

Python-only:

```bash
make test
```

Installation validation:

```bash
make doctor
```

Subsystem-only C tests:

```bash
make unit-io
make unit-particles
```

Everything:

```bash
make check
```

## Notes

- The Python suite does not require PETSc.
- `doctor`, `unit-*`, `unit`, `smoke`, and `check` assume a working PETSc/MPI toolchain.
- Compatibility aliases such as `make install-check` and `make ctest-*` remain available, but the canonical user-facing names are `doctor` and `unit-*`.

## Authoritative Docs

- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/picurv-docs/51_C_Test_Suite_Developer_Guide.html
