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
2. Update user-facing pages first (what changed and how to use it), including `README.md` when entry workflows or top-level commands changed.
3. Update developer/reference pages second (runtime mapping, implementation details).
4. Run repository doc-contract checks (`python3 scripts/audit_function_docs.py` when function comments/docstrings changed).
5. Regenerate docs locally and resolve warnings.
6. Run link checks and verify navigation placement.

## Quality Checklist Before Commit

- Keep cross-page IDs stable (`@page`, `@subpage`, `@ref`) so existing links do not break.
- Run markdown link checking:
  - `python3 scripts/check_markdown_links.py`
- Run the function-doc audit when C/Python executable APIs or test helpers changed:
  - `python3 scripts/audit_function_docs.py`
- Keep `docs/DoxygenLayout.xml` and `docs/mainpage.md` synchronized with new/renamed pages.
- When changing run/cluster lifecycle behavior, keep `README.md`, `05_The_Conductor_Script.md`, `36_Cluster_Run_Guide.md`, `52_Run_Lifecycle_Guide.md`, and the relevant `config/*/guide.md` pages aligned.
- For equations or symbols, verify MathJax rendering in the generated HTML.
- Avoid introducing orphan pages that are discoverable only by URL.

## Entry Pages

- https://vishalkandala.me/picurv-docs/41_Getting_Started_Index.html
- https://vishalkandala.me/picurv-docs/42_User_Guide_Index.html
- https://vishalkandala.me/picurv-docs/43_Developer_Portal_Index.html
