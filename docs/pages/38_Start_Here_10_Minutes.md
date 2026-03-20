@page 38_Start_Here_10_Minutes Start Here in 10 Minutes

@anchor _Start_Here_10_Minutes

This page is the shortest reliable path from repository clone to a successful solver and post-processing run.
It intentionally mirrors the beginning of **@subpage 02_Tutorial_Programmatic_Grid**.
Use this page when you want speed; use the tutorial page when you want deeper explanation.

@tableofcontents

@section p38_intent_sec 1. What This Quickstart Covers

In one pass you will:

1. build binaries,
2. create a runnable case directory,
3. validate all YAML files before launch,
4. run solver and postprocessor,
5. confirm expected output structure.

@section p38_prereq_sec 2. Minimal Prerequisites

You need:

- PETSc configured (`PETSC_DIR`, `PETSC_ARCH`),
- MPI runtime,
- Python 3.

If any of these are missing, follow **@subpage 01_Installation** first.

@section p38_build_sec 3. Build Binaries

From repository root:

```bash
./scripts/picurv build
source etc/picurv.sh
```

The first command builds `bin/simulator`, `bin/postprocessor`, and creates
`bin/picurv` (a symlink to `scripts/picurv`). The second adds `bin/` to your PATH
so `picurv` works from any directory. Add the source line to `~/.bashrc` to make it permanent.

If build fails, go to **@subpage 01_Installation** and verify PETSc/MPI toolchain setup.

@section p38_init_sec 4. Initialize a Starter Case

```bash
./bin/picurv init flat_channel --dest my_case
```

You should get:

```text
my_case/
|- flat_channel.yml
|- Imp-MG-Standard.yml
|- Standard_Output.yml
`- standard_analysis.yml
```

`init` creates the case directory with config files. Runtime binaries (`simulator`, `postprocessor`) are resolved from the project `bin/` directory via PATH — source `etc/picurv.sh` to set it up.

These are the four runtime roles:

- case physics and grid,
- solver numerics,
- monitor/logging controls,
- post-processing pipeline.

@section p38_validate_sec 5. Validate Before Running

```bash
./bin/picurv validate \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

Why this matters:

- catches schema and contract errors early,
- avoids wasting cluster/local runtime on bad config,
- produces structured error messages with key and file context.

Optional launch preview with no execution:

```bash
./bin/picurv run --solve --post-process \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml \
  --dry-run
```

@section p38_run_sec 6. Run Solver and Post

```bash
./bin/picurv run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

This command does:

1. validate + normalize config,
2. create `runs/<run_id>/config/` artifacts,
3. launch `simulator`,
4. launch `postprocessor`.

@section p38_check_sec 7. Confirm Success

Check these locations:

- `runs/<run_id>/config/` (generated control artifacts)
- `runs/<run_id>/logs/` (run logs)
- `output/` (solver outputs, if configured)
- `viz/` or configured post output directory (VTK files)

Quick sanity check in ParaView:

- open `Field_*.vts`,
- add `Slice`,
- color by `Ucat_nodal`.

@section p38_next_sec 8. Continue With Full Guidance

- Full walkthrough with deeper explanation: **@subpage 02_Tutorial_Programmatic_Grid**
- File-based-grid tutorial: **@subpage 03_Tutorial_File-Based_Grid**
- Command/reference guide for orchestration: **@subpage 05_The_Conductor_Script**
- Operational run reuse/restart flow: **@subpage 52_Run_Lifecycle_Guide**
- Troubleshooting map: **@subpage 39_Common_Fatal_Errors**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Start Here in 10 Minutes** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
