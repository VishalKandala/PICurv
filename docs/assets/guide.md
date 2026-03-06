# Assets Guide

This folder stores media used by authored docs pages and the Doxygen main page. Assets in this directory are part of the user-facing product, not internal scratch files, so they should be treated with the same stability expectations as public docs pages.

## Asset Usage Principles

- Prefer relative references such as `@image html assets/<file>` in docs pages.
- Use descriptive filenames that encode scenario or method context.
- Keep filenames stable once published to reduce stale cached links and broken references.
- Prefer deterministic exports (consistent dimensions/aspect ratio) when replacing existing media.

## Current Primary Assets

- `curv.gif`: landing-page simulation preview used to convey flow/particle coupling behavior.
- `paraview_flat_channel.png`: baseline postprocessing visualization example.

## Maintenance Workflow

1. Add or replace media in this directory.
2. Update references in `docs/mainpage.md` and any page-specific usage.
3. Regenerate docs and verify rendering size/clarity.
4. Confirm alt-text and captions remain accurate.

## Review Checklist

- Is the asset still representative of current solver output?
- Is resolution high enough for desktop reading but not unnecessarily large?
- Are all references updated (including old tutorial pages)?
- Does the new image preserve interpretability for CFD newcomers?
