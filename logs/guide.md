# Logs Guide

This directory stores generated logs and warning outputs used during development and CI checks.

## What Lives Here

- `logs/doxygen.warnings`: Doxygen warning output configured by `docs/Doxyfile`.
- ad-hoc local logs created by scripts/tools during diagnostics.

Most solver/postprocessor runtime logs are generated under run/study directories rather than here.

## Recommended Workflow

1. Regenerate docs or run tooling.
2. Inspect warning/log files in this directory.
3. Fix root causes in source/docs rather than suppressing warnings.

## Related Docs

- `docs/pages/35_API_Documentation_Status.md`
- `docs/pages/40_Testing_and_Quality_Guide.md`
