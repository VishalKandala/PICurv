@page 38_Start_Here_10_Minutes Start Here in 10 Minutes

This page is the shortest reliable path from repository clone to a successful solver and post-processing run.
It intentionally mirrors the beginning of **@subpage 02_Tutorial_Programmatic_Grid**.
Use this page when you want speed; use the tutorial page when you want deeper explanation.

@tableofcontents

@section intent_sec 1. What This Quickstart Covers

In one pass you will:

1. build binaries,
2. create a runnable case directory,
3. validate all YAML files before launch,
4. run solver and postprocessor,
5. confirm expected output structure.

@section prereq_sec 2. Minimal Prerequisites

You need:

- PETSc configured (`PETSC_DIR`, `PETSC_ARCH`),
- MPI runtime,
- Python 3.

If any of these are missing, follow **@subpage 01_Installation** first.

@section build_sec 3. Build Binaries

From repository root:

```bash
./scripts/pic.flow build
```

Expected binaries:

- `bin/picsolver`
- `bin/postprocessor`

If build fails, go to **@subpage 01_Installation** and verify PETSc/MPI toolchain setup.

@section init_sec 4. Initialize a Starter Case

```bash
./scripts/pic.flow init flat_channel --dest my_case
```

You should get:

```text
my_case/
|- flat_channel.yml
|- Imp-MG-Standard.yml
|- Standard_Output.yml
`- standard_analysis.yml
```

These are the four runtime roles:

- case physics and grid,
- solver numerics,
- monitor/logging controls,
- post-processing pipeline.

@section validate_sec 5. Validate Before Running

```bash
./scripts/pic.flow validate \
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
./scripts/pic.flow run --solve --post-process \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml \
  --dry-run
```

@section run_sec 6. Run Solver and Post

```bash
./scripts/pic.flow run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

This command does:

1. validate + normalize config,
2. create `runs/<run_id>/config/` artifacts,
3. launch `picsolver`,
4. launch `postprocessor`.

@section check_sec 7. Confirm Success

Check these locations:

- `runs/<run_id>/config/` (generated control artifacts)
- `runs/<run_id>/logs/` (run logs)
- `results/` (solver outputs, if configured)
- `viz/` or configured post output directory (VTK files)

Quick sanity check in ParaView:

- open `Field_*.vts`,
- add `Slice`,
- color by `Ucat_nodal`.

@section next_sec 8. Continue With Full Guidance

- Full walkthrough with deeper explanation: **@subpage 02_Tutorial_Programmatic_Grid**
- File-based-grid tutorial: **@subpage 03_Tutorial_File-Based_Grid**
- Command/reference guide for orchestration: **@subpage 05_The_Conductor_Script**
- Troubleshooting map: **@subpage 39_Common_Fatal_Errors**
