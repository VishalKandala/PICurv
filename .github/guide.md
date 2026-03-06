# GitHub Workflows Guide

This directory stores CI workflow definitions used to protect documentation quality and command-level runtime behavior before merge.

CI here is intentionally scoped to high-signal checks that fail fast when contract or docs hygiene drifts.

## Workflow Files

- `workflows/docs.yml`: builds Doxygen docs and mirrors published HTML artifacts.
- `workflows/quality.yml`: runs CLI smoke tests and markdown link checks.

## How To Evolve Workflows Safely

1. Keep workflow changes minimal and explicit.
2. Prefer adding one new check at a time with clear failure messages.
3. Keep local command parity (`make` targets or scripts) so failures are reproducible outside CI.
4. Update related documentation whenever gate behavior changes.

## Related Docs

- `docs/pages/40_Testing_and_Quality_Guide.md`
- `docs/pages/35_API_Documentation_Status.md`
