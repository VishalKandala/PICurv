# Tests Guide

This directory contains automated smoke and regression checks for `pic.flow`.

## Layout

- `test_cli_smoke.py`: CLI help/validate/dry-run smoke coverage
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
