@page 30_Repository_Navigation Repository Navigation and Directory Guides

This page is the documentation map of the repository layout.
Use it as the single directory-guide index so internal `guide.md` files do not need to appear as standalone documentation pages.

@tableofcontents

@section top_nav_sec 1. Top-Level Repository Areas

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

@section top_links_sec 2. Top-Level Guide Links

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

@section config_nav_sec 3. Configuration Subdirectory Guides

- [Build Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/build/guide.md)
- [Grid Config Library Guide](https://github.com/VishalKandala/PICurv/blob/main/config/grids/guide.md)
- [Monitor Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/monitors/guide.md)
- [Solver Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/solvers/guide.md)
- [Postprocessor Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/postprocessors/guide.md)
- [Scheduler Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/schedulers/guide.md)
- [Study Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/studies/guide.md)

@section docs_nav_sec 4. Documentation-Local Guides

- [Pages Authoring Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/pages/guide.md)
- [Assets Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/assets/guide.md)

@section usage_sec 5. Documentation Organization Rules

1. Repository-level conceptual docs should live in `docs/pages/` as `@page` entries.
2. Directory orientation notes should live in local `guide.md` files.
3. Directory guides are indexed from this page, not promoted to top-level docs navigation.
4. When adding a major new directory, add a `guide.md` and add a link here.

@section related_sec 6. Related Documentation

- **@subpage 41_Getting_Started_Index**
- **@subpage 42_User_Guide_Index**
- **@subpage 43_Developer_Portal_Index**
