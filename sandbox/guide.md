# Sandbox Guide

## Purpose

`sandbox/` is a developer scratch area for experiments that are not yet production-ready. Use this directory to iterate quickly without polluting canonical workflow assets in `<repo>/config/`, `examples/`, `picurv_cli/`, `generators/`, or `docs/pages/`.

## Good Uses

- prototype scripts,
- temporary implementation notes,
- exploratory config snippets,
- one-off data transforms before promotion.
- throwaway cases initialized with `picurv init --dest sandbox/<name>`, whose
  runs and studies stay inside the case directory and are deleted with it.
- cluster handoff material: the analysis scripts and note that accompany a
  staged long run, and the derived results brought back from it.

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

- reusable config -> `<repo>/config/`
- runnable sample case -> `examples/`
- conductor logic -> `picurv_cli/`; standalone generation -> `generators/`
- user/developer explanation -> `docs/pages/`
- production C implementation -> `src/` + `include/`

## Promotion Checklist

1. Is the feature/config reproducible without hidden sandbox context?
2. Is there a stable interface or schema contract?
3. Are docs updated in `docs/pages/` and relevant `guide.md` files?
4. Is ownership clear (module location and maintenance path)?

If all are yes, move it out of sandbox.

## Related Docs

- https://vishalkandala.me/docs/picurv/30_Repository_Navigation.html
- https://vishalkandala.me/docs/picurv/29_Maintenance_Backlog.html
