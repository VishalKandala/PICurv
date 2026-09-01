@page 30_Repository_Navigation Repository Navigation and Directory Guides

@anchor _Repository_Navigation

This page is the documentation map of the repository layout.
Use it as the single directory-guide index so internal `guide.md` files do not need to appear as standalone documentation pages.

@tableofcontents

@section p30_top_nav_sec 1. Top-Level Repository Areas

- Configuration library: `<repo>/config/`
- Documentation source: `docs/`
- Runnable templates: `examples/`
- Public headers: `include/`
- C source implementation: `src/`
- Python conductor package and source-tree entrypoint: `picurv_cli/`
- Tests and fixtures: `tests/`
- Developer scratch area: `sandbox/`
- Logs and warnings: `<repo>/logs/`
- GitHub workflow definitions: `.github/`

@section p30_top_links_sec 2. Top-Level Guide Links

- [Configuration Guide](https://github.com/VishalKandala/PICurv/blob/main/config/guide.md)
- [Documentation Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/guide.md)
- [Examples Guide](https://github.com/VishalKandala/PICurv/blob/main/examples/guide.md)
- [Include Guide](https://github.com/VishalKandala/PICurv/blob/main/include/guide.md)
- [Source Guide](https://github.com/VishalKandala/PICurv/blob/main/src/guide.md)
- [Generators Guide](https://github.com/VishalKandala/PICurv/blob/main/generators/guide.md)
- [PICurv CLI Package Guide](https://github.com/VishalKandala/PICurv/blob/main/picurv_cli/guide.md)
- [Tests Guide](https://github.com/VishalKandala/PICurv/blob/main/tests/guide.md)
- [Smoke Test Guide](https://github.com/VishalKandala/PICurv/blob/main/tests/smoke/guide.md)
- [Test Tooling Guide](https://github.com/VishalKandala/PICurv/blob/main/tests/tooling/guide.md)
- [Sandbox Guide](https://github.com/VishalKandala/PICurv/blob/main/sandbox/guide.md)
- [Logs Guide](https://github.com/VishalKandala/PICurv/blob/main/logs/guide.md)
- [GitHub Workflows Guide](https://github.com/VishalKandala/PICurv/blob/main/.github/guide.md)

@section p30_config_nav_sec 3. Configuration Subdirectory Guides

- [Build Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/build/guide.md)
- [Grid Config Library Guide](https://github.com/VishalKandala/PICurv/blob/main/config/grids/guide.md)
- [Initial Condition Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/initial_conditions/guide.md)
- [Monitor Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/monitors/guide.md)
- [Solver Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/solvers/guide.md)
- [Postprocessor Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/postprocessors/guide.md)
- [Composable Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/profiles/guide.md)
- [Runtime Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/runtime/guide.md)
- [Scheduler Profiles Guide](https://github.com/VishalKandala/PICurv/blob/main/config/schedulers/guide.md)
- [Study Config Guide](https://github.com/VishalKandala/PICurv/blob/main/config/studies/guide.md)

@section p30_examples_nav_sec 4. Example Case Guides

- [Bent Channel Guide](https://github.com/VishalKandala/PICurv/blob/main/examples/bent_channel/guide.md)
- [Flat Channel Guide](https://github.com/VishalKandala/PICurv/blob/main/examples/flat_channel/guide.md)
- [Master Template Guide](https://github.com/VishalKandala/PICurv/blob/main/examples/master_template/guide.md)

@section p30_docs_nav_sec 5. Documentation-Local Guides

- [Pages Authoring Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/pages/guide.md)
- [Assets Guide](https://github.com/VishalKandala/PICurv/blob/main/docs/assets/guide.md)

@section p30_usage_sec 6. Documentation Organization Rules

1. Repository-level conceptual docs should live in `docs/pages/` as `@page` entries.
2. Directory orientation notes should live in local `guide.md` files.
3. Directory guides are indexed from this page, not promoted to top-level docs navigation.
4. When adding a major new directory, add a `guide.md` and add a link here.

@section p30_related_sec 7. Related Documentation

- **@subpage 41_Getting_Started_Index**
- **@subpage 42_User_Guide_Index**
- **@subpage 43_Developer_Portal_Index**
