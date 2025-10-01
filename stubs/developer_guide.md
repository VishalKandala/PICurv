/** @page devguide Developer Guide {#devguide}
@tableofcontents

@section build Build & Toolchain
- Compiler: GCC/Clang, C17
- Dependencies: (PETSc/MPI/etc. if applicable)
- Targets: `picsolver`, `inttest`, `postprocess` (example)

@section style Coding Style
- Headers in `include/` define the public API; keep comments Doxygen-ready.
- Use `///` for brief one-line docs and `/** ... */` for blocks.

@section testing Testing
- `make test` runs unit/integration tests (describe how to add a test).

@section perf Performance Notes
- Log events around walking-search, interpolation, solver phases.
- Prefer AoS/SoA decisions documented here with rationale.

@section contrib Contributing
- Branches: `main`, feature branches
- Commit style: imperative, scoped
- PR checklist: build, tests, docs, clang-format
*/
