@page 59_Function_Identity_and_Observability_Specification Function Identity and Observability Specification

@anchor _Function_Identity_and_Observability_Specification

How a C function is named and timed, and the one part of the original deferred work,
the logging allowlist, that is still open.

@tableofcontents

@section p59_identity_sec 1. Function Identity

A function's name everywhere is C99's predefined `__func__`: profiling events, logging
filters (`LOG_ALLOW`, `is_function_allowed(__func__)`) and PETSc's error traces (PETSc
3.20's `PETSC_FUNCTION_NAME` is `__func__`). No file defines or reads `__FUNCT__`; a
hand-written name could drift from the function it labels, which is how
`ComputeKFaceMetrics` was once reported as `ComputeJFaceMetrics`. Names written in
`monitor.yml` (`logging.enabled_functions`, `profiling.timestep_output.functions`,
`profiling.helpers`) are these function names.

@section p59_profiling_sec 2. Profiling Contract

Every function defined in `src/*.c` that is not `static inline` times itself:

- it begins one timer, `PROFILE_FUNCTION_BEGIN` or `PROFILE_HELPER_BEGIN`, as its first
  statement after `PetscFunctionBeginUser` (an early return placed before the timer, such
  as a once-only guard, needs no END);
- every explicit `return` and `PetscFunctionReturn` after it is immediately preceded by
  `PROFILE_FUNCTION_END` (or its alias `PROFILE_HELPER_END`), and a body that can fall off
  its end ends with one. Error returns inside `SETERRQ`, `PetscCall` and `CHKERRQ` end the
  run and need none.

A missed END leaves the function's PETSc event open; its later calls are then not recorded.
`make audit-profiling` (also run by `make check`) enforces both rules and the absence of
`__FUNCT__`, against the declared exclusions in `tests/tooling/audit_profiling_coverage.py`:
`main`, the asynchronous signal handler, the logging filter itself (`get_log_level`,
`is_function_allowed`, `set_allowed_functions`), the profiler's own functions, and two
dependency-free leaf modules, `src/checksum.c` and `src/statistics_moments.c`, which include
only `petscsys.h`.

Each timer is a PETSc log event named by `__func__`, registered on the function's first call
and cached at the call site. Registration waits for `PetscInitialize` and stops after
`PetscFinalize`, so a function that runs outside them (argument handling in `main`) stays
untimed instead of failing. There are two kinds:

- `PROFILE_FUNCTION_BEGIN`, class `PICurv`: functions called a handful of times per step,
  solver pass or post-processing step. Always timed.
- `PROFILE_HELPER_BEGIN`, class `PICurvHelper`: functions called per particle, per cell or
  in inner loops. Off unless `profiling.helpers` names them or is `all`
  (@ref p09_profiling_sec).

The split is set by measured calls per rank per step with every helper on, over a particle
run, an implicit-momentum flow run, a Newton–Krylov run and a post-processing step: at most 20
is an aggregate. A function no measured case reached starts as a helper and is promoted once a
case that exercises it shows it is called rarely. A new function starts as a helper unless it is
plainly once per step.

@section p59_cost_sec 3. Measured Cost

Per call, PETSc 3.20.3 built optimized (`--with-debugging=0`), one rank on Grace:

| Timer | Cost |
|---|---|
| empty loop iteration (reference) | 2.5 ns |
| helper, off | 3.5 ns |
| PETSc event, active | 98 ns |
| former PICurv registry timer (removed) | 54 ns |

A particle step makes about 19 helper calls per particle, so the former always-on timers cost
1.4 s per step at 1.64M particles per rank (20–33% of step time); off, helpers cost a few
milliseconds. `helpers: all` costs about 100 ns per call and roughly doubles a particle
step, so it is for short diagnostic runs. A debug PETSc build pays about 26 ns per switched-off
helper, which is visible in debug runs.

@section p59_logging_sec 4. Still Deferred: Logging Filter

Logging still filters by name: `LOG_ALLOW` checks the level first, then `is_function_allowed`
compares `__func__` against the allowlist. Emitted logging is dominated by formatting and I/O,
and the level check keeps filtered calls cheap at production levels, so this is unchanged. If a
benchmark shows the name comparison is material at debug levels, the recorded direction is to
resolve the configured names once into an indexed allow table and cache each call site's entry,
as the profiler does, without changing log formats or the YAML names.

This work does not reuse `FieldId` or create one universal identity namespace, assign IDs to
filenames, configuration values, PETSc object names, DMSwarm field names or post-processing
labels, or replace boundary-condition dispatch resolved at ingress. One-time string comparisons
at external ingress and presentation boundaries are correct and remain.
