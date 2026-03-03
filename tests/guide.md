# Tests Guide

This directory contains automated smoke and regression checks for `picurv`.

## Layout

- `test_cli_smoke.py`: CLI help/validate/dry-run smoke coverage
- `test_case_maintenance.py`: regression checks for case-origin metadata, build-from-case, sync/pull maintenance commands, status reporting, and safe config pruning
- `test_config_regressions.py`: ingress/schema regression checks for config translation and drift guards
- `test_repo_consistency.py`: repository-wide example validation and stale-doc contract scans
- `fixtures/valid/`: canonical valid YAML input sets
- `fixtures/invalid/`: intentionally broken YAML sets with expected failure patterns

## Local Run

```bash
pytest -q
```

If local Python dependency installation is restricted, run tests in CI (`.github/workflows/quality.yml`).

## Authoritative Docs

- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/picurv-docs/39_Common_Fatal_Errors.html
