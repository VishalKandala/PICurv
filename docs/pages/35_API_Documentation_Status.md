@page 35_API_Documentation_Status API Documentation Status

@anchor _API_Documentation_Status

This page describes current API-doc quality, warning sources, and the expected standard for new changes.

@tableofcontents

@section p35_status_sec 1. Current Status

Function-level documentation coverage is now enforced across:

- public C headers in `include/`,
- C implementations in `src/`,
- C tests in `tests/c/`,
- Python product code in `picurv_cli/` and `generators/`,
- Python tests and documentation tooling in `tests/`.

The current repository contract requires every executable function to have an
attached descriptive comment/docstring. Public C headers own the rendered
Doxygen API contract, including `@brief`, exact parameter coverage, and
documented return values where applicable; implementation comments preserve
local intent without duplicating public parameter documentation.

@section p35_warning_sec 2. Warning Log and Build Path

Doxygen warnings are written to:

- `logs/doxygen.warnings`

Configured in `docs/Doxyfile` via:

- `WARN_LOGFILE = ../logs/doxygen.warnings`

Repository consistency checks now also enforce function-level documentation via:

- `tests/tooling/audit_function_docs.py`

GitHub Actions quality CI also runs the audit explicitly before `pytest -q` on pull requests and pushes to `main`.

The audit requires a meaningful explanation, not merely a tag or symbol label:
public C declarations must have a specific Doxygen brief, private C helpers
must have a specific attached implementation comment, and Python functions must
have a specific brief docstring. It scans production code, generators, tests,
and documentation tooling.

`make build-docs` remains useful for rendered output verification, but the audit
script is the primary completeness gate for function comments. For a
commit-scoped claim that the documentation is current, use `make certify-docs`.
It requires a clean worktree, runs link/API/ingress/configuration/Doxygen gates
plus `make check-full`, and writes an ignored certificate named
`logs/documentation-certificate-<full-sha>.md`. The certificate is valid only
through that exact Git commit.

Every rendered Doxygen page also carries a bottom banner linking to the commit
whose source tree produced it. The banner identifies the documentation revision;
the full certificate records whether the PETSc/MPI runtime gate was included.

@section p35_expected_sec 3. Expected Standard For New APIs

For newly added or modified functions and test helpers:

1. each parameter must have exactly one matching `@param`,
2. public-header summaries must describe the result, state change, or numerical
   role; labels such as “helper function,” “public interface,” “implementation,”
   and “internal helper” are rejected by the API audit,
3. return/CHKERRQ semantics should be documented,
4. cross-module dependencies should be explicit when non-obvious,
5. Python functions must use Doxygen-compatible docstrings, not plain one-line docstrings.

For C test files, concise briefs are acceptable, but they must still say what
the test/helper verifies or sets up. Placeholder summaries like `Test-local
routine.` are rejected.

Minimum acceptable quality is interface correctness and discoverability, even when deep theoretical derivation is documented elsewhere.

For nontrivial numerical APIs, the header must also identify the coordinate or
field convention, ownership/ghost-state preconditions, and output meaning when
those facts determine correct use. A parameter name by itself is not a contract.

Inside a function, comment only the information that code cannot communicate
reliably by itself: numerical rationale, invariant or unit convention,
distributed-memory/ownership boundary, intentionally surprising early return,
or error-handling constraint. Do not add line-by-line narration of obvious
assignments or assertions. These in-body explanations are reviewed with the
code change; the automated gate enforces the function-level baseline and
rejects known stub language.

@section p35_workflow_sec 4. Practical Cleanup Workflow

1. run `python3 tests/tooling/audit_function_docs.py`,
2. patch one module at a time (`include/*.h` + matching `src/*.c`, plus affected Python helpers/tests),
3. run the repo consistency tests that wrap the audit,
4. run docs build (`make build-docs`) to verify rendered output where needed.

Batching by subsystem still keeps review scope manageable, but the audit gate now
prevents silent regressions between cleanup passes.

@section p35_refs_sec 5. Related Pages

- **@subpage 13_Code_Architecture**
- **@subpage 29_Maintenance_Backlog**
- **@subpage 40_Testing_and_Quality_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **API Documentation Status** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
