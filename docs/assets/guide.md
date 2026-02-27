# Assets Guide

This folder stores media used by authored docs pages and the Doxygen main page.

## Usage

- Prefer relative references such as `@image html assets/<file>` in docs pages.
- Keep filenames stable once published to avoid stale cached links on the website.
- Use descriptive names tied to case/method context.

## Current Primary Assets

- `curv.gif`: landing-page simulation preview.
- `paraview_flat_channel.png`: sample postprocessing visualization.

## Maintenance

If you replace an asset with a new revision, prefer re-exporting to the same dimensions and filename where possible, or update all references in `docs/mainpage.md` and any pages using it.
