@page 35_API_Documentation_Status API Documentation Status

This page defines how API-doc quality is checked and where warning output is tracked.

@tableofcontents

@section status_sec 1. Current Status

Function/parameter documentation is partially complete across the codebase.
There are known mismatches in some headers (duplicate `@param`, stale parameter names, and missing docs for some signatures).
This should be treated as a known quality gap until a clean Doxygen warning pass is achieved.

@section warning_log_sec 2. Warning Log Location

Doxygen warning output is configured to:
- `logs/doxygen.warnings`

This is controlled by `docs/Doxyfile` (`WARN_LOGFILE = ../logs/doxygen.warnings`).

@section workflow_sec 3. Verification Workflow

1. Run `make build-docs`.
2. Inspect `logs/doxygen.warnings`.
3. Fix warnings in small batches by module (`include/*.h` + matching `src/*.c`).
4. Re-run and track deltas.

@section quality_gate_sec 4. Practical Quality Gate

For new/modified APIs:
- all arguments should have exactly one matching `@param` entry
- no stale argument names in doc blocks
- no undocumented public function declarations in `include/`

This does not enforce full theoretical derivations; it enforces interface correctness and discoverability.
