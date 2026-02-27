@page 38_Start_Here_10_Minutes Start Here in 10 Minutes

This is the fastest path from clone to first solver+post output.

@tableofcontents

@section step1_sec 1. Install Prerequisites

Follow **@subpage 01_Installation** for full setup details. For quick local runs you need:
- PETSc environment (`PETSC_DIR`, `PETSC_ARCH`)
- MPI, Make, Python 3

@section step2_sec 2. Build Binaries

From repository root:

```bash
./scripts/pic.flow build
```

Expected artifacts:
- `bin/picsolver`
- `bin/postprocessor`

@section step3_sec 3. Initialize a Template Case

```bash
./scripts/pic.flow init flat_channel --dest my_case
```

@section step4_sec 4. Validate Configs (No Run Yet)

```bash
./scripts/pic.flow validate \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

Optional launch preview (no files created):

```bash
./scripts/pic.flow run --solve --post-process \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml \
  --dry-run
```

@section step5_sec 5. Run Solver + Post

```bash
./scripts/pic.flow run --solve --post-process -n 4 \
  --case my_case/flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/Standard_Output.yml \
  --post my_case/standard_analysis.yml
```

@section step6_sec 6. Inspect Outputs

Check generated directories:
- `runs/<run_id>/config/`
- `runs/<run_id>/logs/`
- solver output path from `monitor.yml` (typically `results/`)
- post output path from `post.yml` (`io.output_directory`)

@section step7_sec 7. First Troubleshooting Actions

If validation or run fails:
1. Re-run `pic.flow validate ...` for the same files.
2. Check structured error lines (`ERROR <CODE> | key=... | file=...`).
3. Compare against `examples/master_template/*.yml`.
4. Use **@subpage 39_Common_Fatal_Errors** for mapped fixes.

Next:
- **@subpage 11_User_How_To_Guides**
- **@subpage 05_The_Conductor_Script**
