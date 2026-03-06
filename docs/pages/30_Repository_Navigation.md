@page 30_Repository_Navigation Repository Navigation and Directory Guides

@anchor _Repository_Navigation

This page is the documentation map of the repository layout.
Use it as the single directory-guide index so internal `guide.md` files do not need to appear as standalone documentation pages.

@tableofcontents

@section p30_top_nav_sec 1. Top-Level Repository Areas

- Configuration library: `config/`
- Documentation source: `docs/`
- Runnable templates: `examples/`
- Public headers: `include/`
- C source implementation: `src/`
- Automation scripts: `scripts/`
- Tests and fixtures: `tests/`
- Developer scratch area: `sandbox/`
- Logs and warnings: `logs/`
- GitHub workflow definitions: `.github/`

@section p30_top_links_sec 2. Top-Level Guide Links

- [Configuration Guide](https://github.com/VishalKandala/PICurv/blob/main/config/guide.md)
- [Documentation Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/guide.md)
- [Examples Guide](https://github.com/VishalKandala/PICurv/blob/main/examples/guide.md)
- [Include Guide](https://github.com/VishalKandala/PICurv/blob/main/include/guide.md)
- [Source Guide](https://github.com/VishalKandala/PICurv/blob/main/src/guide.md)
- [Scripts Guide](https://github.com/VishalKandala/PICurv/blob/main/scripts/guide.md)
- [Tests Guide](https://github.com/VishalKandala/PICurv/blob/main/tests/guide.md)
- [Sandbox Guide](https://github.com/VishalKandala/PICurv/blob/main/sandbox/guide.md)
- [Logs Guide](https://github.com/VishalKandala/PICurv/blob/main/logs/guide.md)
- [GitHub Workflows Guide](https://github.com/VishalKandala/PICurv/blob/main/.github/guide.md)

@section p30_config_nav_sec 3. Configuration Subdirectory Guides

- [Build Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/build/guide.md)
- [Grid Config Library Guide](https://github.com/VishalKandala/PICurv/blob/main/config/grids/guide.md)
- [Monitor Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/monitors/guide.md)
- [Solver Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/solvers/guide.md)
- [Postprocessor Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/postprocessors/guide.md)
- [Scheduler Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/schedulers/guide.md)
- [Study Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/studies/guide.md)

@section p30_docs_nav_sec 4. Documentation-Local Guides

- [Pages Authoring Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/pages/guide.md)
- [Assets Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/assets/guide.md)

@section p30_usage_sec 5. Documentation Organization Rules

1. Repository-level conceptual docs should live in `docs/pages/` as `@page` entries.
2. Directory orientation notes should live in local `guide.md` files.
3. Directory guides are indexed from this page, not promoted to top-level docs navigation.
4. When adding a major new directory, add a `guide.md` and add a link here.

@section p30_related_sec 6. Related Documentation

- **@subpage 41_Getting_Started_Index**
- **@subpage 42_User_Guide_Index**
- **@subpage 43_Developer_Portal_Index**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Repository Navigation and Directory Guides** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.

