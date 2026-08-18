@page 59_Function_Identity_and_Observability_Specification Function Identity and Observability Specification

@anchor _Function_Identity_and_Observability_Specification

This deferred specification records the possible optimization of function-name
filtering and profiling. It is not part of the field-statistics implementation
and currently has no implementation branch.

@tableofcontents

@section p59_problem_sec 1. Problem to Measure

Logging allowlists and profiling currently identify instrumented functions by
name. Emitted logging is normally dominated by formatting and I/O, but filtered
calls in hot loops and repeated profiler start/end name searches may have
measurable cost. That cost must be benchmarked before changing the system.

@section p59_direction_sec 2. Recorded Direction

If benchmarks justify the work:

- give instrumented application functions a compile-time `FunctionId` or a
  call-site-cached registry handle;
- resolve configured names once into a bitset/indexed allow table;
- make profiler start/end update directly indexed state;
- check log level before function filtering in every hot-path macro;
- parse profiling mode once at ingress; and
- retain canonical function names for YAML, diagnostics, errors, and output.

The implementation should preserve existing logging macros, log formats,
profiling files, and default behavior as far as possible. A cached handle may be
less intrusive than a manually maintained enum of every C function; the
benchmark/prototype phase chooses between them.

@section p59_boundaries_sec 3. Explicit Boundaries

This work does not:

- reuse `FieldId` or create one universal identity namespace;
- assign IDs to filenames, arbitrary configuration values, PETSc object names,
  DMSwarm field names, or user-defined postprocessing labels;
- replace boundary-condition enums/function pointers already resolved at
  ingress;
- redesign runtime physical-solution monitoring or scientific statistics; or
- block @ref 58_Turbulence_Statistics_Pipeline_Specification.

One-time string comparisons at external ingress and presentation boundaries are
correct and remain.

@section p59_plan_sec 4. Future Plan

1. Benchmark disabled, filtered, selected, and all-function logging/profiling on
   representative serial and MPI runs.
2. Inventory call sites and macros; identify only repeated runtime lookups.
3. Prototype cached handles and compile-time IDs without changing output.
4. Select the simpler design that demonstrates a material improvement.
5. Add behavior, concurrency, and overhead regressions, then update monitoring
   and contributor documentation.

This work receives its own branch series only after explicit approval.

