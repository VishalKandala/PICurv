@page 40_Testing_and_Quality_Guide Testing, Smoke, and Quality Gates Guide

This page explains how to verify the Iteration A CLI features end-to-end, how the smoke tests are structured, and how CI enforces quality gates.

@tableofcontents

@section scope_sec 1. What This Guide Covers

Iteration A introduced:
- `picurv validate`
- `picurv run --dry-run` and `--format json`
- standardized error envelope:
  - `ERROR <CODE> | key=<...> | file=<...> | message=<...> | hint=<...>`
- CI workflow for CLI smoke tests and markdown link checks

This guide documents:
1. manual smoke commands for users and maintainers,
2. automated smoke tests in `tests/test_cli_smoke.py`,
3. workflow behavior in `.github/workflows/quality.yml`,
4. extension patterns when adding future CLI features.

@section quick_smoke_sec 2. Quick Smoke Matrix (Manual)

Run from repository root.

@subsection help_smoke_ssec 2.1 Command Discovery Smoke

```bash
./bin/picurv --help
./bin/picurv run --help
./bin/picurv validate --help
```

Expected:
- top-level help lists `validate`,
- run help lists `--dry-run` and `--format {text,json}`,
- validate help lists file-role flags and `--strict`.

@subsection validate_smoke_ssec 2.2 Config-Only Validation Smoke

Valid fixture set:

```bash
./bin/picurv validate \
  --case tests/fixtures/valid/case.yml \
  --solver tests/fixtures/valid/solver.yml \
  --monitor tests/fixtures/valid/monitor.yml \
  --post tests/fixtures/valid/post.yml \
  --cluster tests/fixtures/valid/cluster.yml \
  --study tests/fixtures/valid/study.yml
```

Expected:
- exit code `0`,
- `[SUCCESS] Validation completed for ... file(s)`.

Invalid fixture example:

```bash
./bin/picurv validate \
  --case tests/fixtures/invalid/case_missing_properties.yml \
  --solver tests/fixtures/valid/solver.yml \
  --monitor tests/fixtures/valid/monitor.yml
```

Expected:
- exit code `1`,
- at least one structured error line with `ERROR CFG_MISSING_SECTION`.

@subsection dryrun_smoke_ssec 2.3 Dry-Run Smoke (No File Writes)

Human-readable plan:

```bash
./bin/picurv run --solve --post-process \
  --case tests/fixtures/valid/case.yml \
  --solver tests/fixtures/valid/solver.yml \
  --monitor tests/fixtures/valid/monitor.yml \
  --post tests/fixtures/valid/post.yml \
  --dry-run
```

Machine-readable plan:

```bash
./bin/picurv run --solve --post-process \
  --case tests/fixtures/valid/case.yml \
  --solver tests/fixtures/valid/solver.yml \
  --monitor tests/fixtures/valid/monitor.yml \
  --post tests/fixtures/valid/post.yml \
  --dry-run --format json
```

Expected:
- no new run directories/files are created,
- JSON output parses cleanly and includes `mode`, `launch_mode`, `stages`, and `artifacts`.

@section code_map_sec 3. Code Map (Where Behavior Lives)

Primary implementation in `scripts/picurv`:
- standardized error emitter: `emit_structured_error(...)`
- CLI usage failure exit code (`2`): `fail_cli_usage(...)`
- dry-run plan builder: `build_run_dry_plan(args)`
- dry-run renderer (text/json): `render_run_dry_plan(plan, output_format=...)`
- validate command implementation: `validate_workflow(args)`
- parser composition: `build_main_parser()`
- dispatch and argument combination checks: `dispatch_command(args)`

Related docs:
- **@subpage 05_The_Conductor_Script**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 14_Config_Contract**

@section automated_sec 4. Automated Smoke Tests

Automated smoke tests live in:
- `tests/test_cli_smoke.py`

Fixture packs:
- valid configs: `tests/fixtures/valid/`
- invalid configs: `tests/fixtures/invalid/`

Current test categories:
1. help output smoke,
2. validate pass/fail behavior,
3. structured error presence for invalid inputs,
4. dry-run no-write guarantee,
5. dry-run JSON schema sanity,
6. markdown link checker pass.

Run locally (when `pytest` is available):

```bash
pytest -q
```

@section ci_sec 5. CI Quality Gate Behavior

Workflow file:
- `.github/workflows/quality.yml`

Pipeline steps:
1. checkout,
2. setup Python,
3. install `pytest`, `pyyaml`, `numpy`,
4. run `pytest -q`,
5. run `python scripts/check_markdown_links.py`.

If either step fails, the workflow fails and blocks merge (for repositories using required checks).

Separate docs build/mirror workflow:
- `.github/workflows/docs.yml`
- handles Doxygen generation and docs artifact mirroring.

@section linkcheck_sec 6. Markdown Link Checking

Link checker script:
- `scripts/check_markdown_links.py`

Coverage:
- `README.md`
- `docs/**/*.md` (excluding generated `docs_build/`)

It checks local relative links only and intentionally skips:
- `http://...`
- `https://...`
- `mailto:...`
- in-page anchors (`#...`)

@section extend_sec 7. How to Extend Smoke Coverage

When adding a new CLI option or subcommand:
1. add a help-output assertion in `tests/test_cli_smoke.py`,
2. add one success-path and one failure-path test,
3. add/update fixtures under `tests/fixtures/`,
4. document user-facing behavior in **@subpage 05_The_Conductor_Script**,
5. update this page’s smoke matrix if the command is user-critical.

When adding new docs pages:
1. add nav entries in `docs/DoxygenLayout.xml`,
2. run `python scripts/check_markdown_links.py`,
3. ensure links from at least one entry page (`mainpage`/conductor/how-to).

@section troubleshooting_sec 8. Local Troubleshooting

Common local issues:
- `pytest: command not found`:
  - install pytest in your local Python environment,
  - or run tests in CI where dependencies are provisioned.
- `No module named pytest`:
  - use the Python interpreter/environment intended for development.
- network-restricted environments:
  - local package install may fail; rely on CI run for authoritative smoke status.

@section refs_sec 9. External References

- Pytest docs: https://docs.pytest.org/
- GitHub Actions workflow syntax: https://docs.github.com/actions/using-workflows/workflow-syntax-for-github-actions
- CommonMark link syntax reference: https://spec.commonmark.org/
