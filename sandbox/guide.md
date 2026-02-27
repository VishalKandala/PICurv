# Sandbox Guide

## Purpose

`sandbox/` is an optional developer scratch area for experiments that are not yet production-ready.

Use it for:
- prototype scripts
- temporary notes
- exploratory config snippets
- one-off data transforms before promoting to tracked workflows

Do not use it for:
- canonical run configurations
- release artifacts
- required build/runtime dependencies

## Directory Hygiene

Recommended structure inside `sandbox/`:
- `notes/` for temporary design notes
- `prototypes/` for throwaway code
- `scratch-configs/` for trial YAML or `.cfg`
- `tmp/` for generated transient files

Keep these lightweight and delete stale content regularly.

## Promotion Rules

When an experiment stabilizes, promote it out of `sandbox/`:
- reusable config -> `config/`
- runnable sample case -> `examples/`
- automation logic -> `scripts/`
- user/developer explanation -> `docs/pages/`
- production C implementation -> `src/` + `include/`

## Review Checklist Before Promotion

1. Is the feature/config reproducible without sandbox context?
2. Is there a stable interface or schema contract?
3. Are docs updated in `docs/pages/` and relevant `guide.md` files?
4. Is ownership clear (module/file location, maintenance path)?

If all are yes, move it out of sandbox.

## Related Docs

- https://vishalkandala.me/picurv-docs/30_Repository_Navigation.html
- https://vishalkandala.me/picurv-docs/29_Maintenance_Backlog.html
