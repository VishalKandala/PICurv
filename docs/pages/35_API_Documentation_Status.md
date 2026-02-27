@page 35_API_Documentation_Status API Documentation Status

This page describes current API-doc quality, warning sources, and the expected standard for new changes.

@tableofcontents

@section status_sec 1. Current Status

API documentation quality is improving but not yet uniformly clean across all modules.
Known issues in some headers/source pairs include:

- stale parameter names in Doxygen blocks,
- duplicate/missing `@param` entries,
- partial function-level docs on legacy-heavy paths.

Treat these as active maintenance items rather than hidden debt.

@section warning_sec 2. Warning Log and Build Path

Doxygen warnings are written to:

- `logs/doxygen.warnings`

Configured in `docs/Doxyfile` via:

- `WARN_LOGFILE = ../logs/doxygen.warnings`

CI currently enforces Markdown-link integrity; Doxygen warning cleanup remains a repository maintenance objective.

@section expected_sec 3. Expected Standard For New APIs

For newly added or modified public functions:

1. each parameter must have exactly one matching `@param`,
2. function summary should describe behavior and side effects,
3. return/CHKERRQ semantics should be documented,
4. cross-module dependencies should be explicit when non-obvious.

Minimum acceptable quality is interface correctness and discoverability, even when deep theoretical derivation is documented elsewhere.

@section workflow_sec 4. Practical Cleanup Workflow

1. run docs build (`make build-docs`) where available,
2. inspect `logs/doxygen.warnings`,
3. patch one module at a time (`include/*.h` + matching `src/*.c`),
4. re-run warning pass and track net reduction.

Batching by module avoids regressions and keeps review scope manageable.

@section refs_sec 5. Related Pages

- **@subpage 13_Code_Architecture**
- **@subpage 29_Maintenance_Backlog**
- **@subpage 40_Testing_and_Quality_Guide**
