# Docs Guide

This directory contains authored documentation source, Doxygen configuration, and documentation assets.

## Core Files

- `Doxyfile`: documentation build configuration.
- `DoxygenLayout.xml`: top-level nav and tab layout.
- `mainpage.md`: landing page content (`@mainpage`).
- `CHANGELOG.md`: release/change notes (`@page 18_Changelog`).
- `custom.css`: site-level visual overrides.

## Subdirectories

- `pages/`: authored long-form docs pages.
- `assets/`: media used by landing page and tutorials.

## Maintenance Checklist

1. Keep cross-page IDs stable (`@page`, `@subpage`, `@ref`).
2. Run markdown link checker before pushing:
   - `python3 scripts/check_markdown_links.py`
3. Keep Doxygen nav structure in sync with new index pages and docs organization.
4. If equations are added, verify they render with configured MathJax settings.

## Entry Pages

- `docs/pages/41_Getting_Started_Index.md`
- `docs/pages/42_User_Guide_Index.md`
- `docs/pages/43_Developer_Portal_Index.md`
