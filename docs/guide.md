# Docs Guide

This directory is the authored source of the PICurv documentation site. It contains long-form user and developer pages, the Doxygen build configuration, and visual assets used by landing/tutorial pages.

The goal of this guide is to explain not only where files live, but how each documentation artifact participates in the release workflow. If you are updating solver behavior, configuration contracts, or runtime diagnostics, this directory is where you make those changes visible and actionable for users.

## Core Files and Their Roles

- `Doxyfile`: controls documentation generation options, warning behavior, math rendering, and source indexing policy.
- `DoxygenLayout.xml`: defines navigation hierarchy and tab ordering in the generated site.
- `mainpage.md`: root landing page (`@mainpage`) that frames the project and links users into key entry points.
- `CHANGELOG.md`: release-oriented behavior notes (`@page 18_Changelog`) tied to contract/runtime evolution.
- `custom.css`: site-level presentation overrides for generated HTML.
- `Documentation_Code_Alignment_Audit_2026-03-06.md`: code-vs-doc gap audit artifact covering `src/`, `scripts/`, and test suites.

## Subdirectories and Intent

- `pages/`: canonical long-form documentation pages. These are the primary location for method details, user procedures, and developer extension guidance.
- `assets/`: images/media referenced by pages, tutorial walkthroughs, and the main landing page.

## Documentation Maintenance Workflow

1. Identify the behavior change and map it to impacted docs pages.
2. Update user-facing pages first (what changed and how to use it).
3. Update developer/reference pages second (runtime mapping, implementation details).
4. Regenerate docs locally and resolve warnings.
5. Run link checks and verify navigation placement.

## Quality Checklist Before Commit

- Keep cross-page IDs stable (`@page`, `@subpage`, `@ref`) so existing links do not break.
- Run markdown link checking:
  - `python3 scripts/check_markdown_links.py`
- Keep `docs/DoxygenLayout.xml` and `docs/mainpage.md` synchronized with new/renamed pages.
- For equations or symbols, verify MathJax rendering in the generated HTML.
- Avoid introducing orphan pages that are discoverable only by URL.

## Entry Pages

- https://vishalkandala.me/picurv-docs/41_Getting_Started_Index.html
- https://vishalkandala.me/picurv-docs/42_User_Guide_Index.html
- https://vishalkandala.me/picurv-docs/43_Developer_Portal_Index.html
