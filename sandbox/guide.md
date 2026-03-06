# Sandbox Guide

## Purpose

`sandbox/` is a developer scratch area for experiments that are not yet production-ready. Use this directory to iterate quickly without polluting canonical workflow assets in `config/`, `examples/`, `scripts/`, or `docs/pages/`.

## Good Uses

- prototype scripts,
- temporary implementation notes,
- exploratory config snippets,
- one-off data transforms before promotion.

## Bad Uses

- canonical run configurations,
- release artifacts,
- required build/runtime dependencies,
- anything another developer must rely on without additional context.

## Suggested Internal Structure

- `notes/` for temporary design notes.
- `prototypes/` for throwaway code.
- `scratch-configs/` for trial YAML or `.cfg`.
- `tmp/` for generated transient files.

Keep these lightweight and delete stale content regularly.

## Promotion Rules

When an experiment stabilizes, promote it out of `sandbox/`:

- reusable config -> `config/`
- runnable sample case -> `examples/`
- automation logic -> `scripts/`
- user/developer explanation -> `docs/pages/`
- production C implementation -> `src/` + `include/`

## Promotion Checklist

1. Is the feature/config reproducible without hidden sandbox context?
2. Is there a stable interface or schema contract?
3. Are docs updated in `docs/pages/` and relevant `guide.md` files?
4. Is ownership clear (module location and maintenance path)?

If all are yes, move it out of sandbox.

## Related Docs

- https://vishalkandala.me/picurv-docs/30_Repository_Navigation.html
- https://vishalkandala.me/picurv-docs/29_Maintenance_Backlog.html
