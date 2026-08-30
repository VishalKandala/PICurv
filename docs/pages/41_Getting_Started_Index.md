@page 41_Getting_Started_Index Getting Started

@anchor _Getting_Started_Index

@pagemeta{Tutorial, New users, Verified end to end}

This is the fastest reliable path from a fresh clone to a finished run. It is
ordered to minimize first-run failure: build first, run a known-good case, then
inspect what it produced. PICurv's YAML files are modular by design, so the goal
is to learn to recombine profiles rather than to write one large configuration
per run.

Work through it in one pass. When you want deeper explanation of any step, the
full walkthrough is **@subpage 02_Tutorial_Programmatic_Grid**.

@note **Verification status.** Every command on this page has been executed against
current `main`, including the full solve-and-post path, which completed in about ten
seconds on a debug build and produced the visualization file named in section 6. It is
not yet gated by CI, so it is verified-by-execution rather than continuously enforced.

@tableofcontents

@section p41_prereq_sec 1. Prerequisites

You need:

- PETSc configured, with `PETSC_DIR` set. `PETSC_ARCH` is required only for old-style
  in-tree PETSc builds; prefix installs do not use it.
- an MPI runtime,
- **Python 3.10 or newer** if you use the managed CLI bootstrap
  (`bootstrap_install.sh` prefers 3.12/3.11/3.10). The repository tooling itself runs
  on 3.8, but the bootstrap path expects a newer interpreter.

If any of these are missing, work through **@subpage 01_Installation** first and
return here.

@section p41_build_sec 2. Build

From the repository root:

```bash
./picurv_cli/picurv build
source etc/picurv.sh
```

The first command builds `bin/simulator` and `bin/postprocessor`, and creates
`bin/picurv` as a launcher. The second puts `bin/` on your PATH, exports the
managed CLI Python when bootstrap created one, and exposes `picurv_cli/` as a
fallback so `picurv` works from any directory even before the launcher is
rebuilt. Add the `source` line to your shell profile to make it permanent.

The build stream is written to `<repo>/logs/build.log`. If the build fails, return to
**@subpage 01_Installation** and verify the PETSc/MPI toolchain.

@section p41_init_sec 3. Create a case

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

Those four files are the four runtime roles: case physics and grid, solver
numerics, monitor and logging controls, and the post-processing pipeline. Runtime
binaries are resolved from `bin/` through your PATH, which is what sourcing
`etc/picurv.sh` set up.

@section p41_validate_sec 4. Validate before running

```bash
./bin/picurv validate \
  --case my_case/quickstart_flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/quickstart_Standard_Output.yml \
  --post my_case/quickstart_standard_analysis.yml
```

Validation catches schema and contract errors before any compute is spent, and
reports them with the offending key and file. To preview the launch without
executing anything, add `--dry-run` to the `run` command below; it prints the
resolved commands and the artifacts the run would produce.

@section p41_run_sec 5. Run

```bash
./bin/picurv run --solve --post-process -n 2 \
  --case my_case/quickstart_flat_channel.yml \
  --solver my_case/Imp-MG-Standard.yml \
  --monitor my_case/quickstart_Standard_Output.yml \
  --post my_case/quickstart_standard_analysis.yml
```

The `quickstart_*` files are a deliberately small variant: a 9x9x17 grid over 20
timesteps, which completes **solve and post-processing in about ten seconds** even in
a PETSc debug build. The unprefixed `flat_channel.yml` files are the full case - a
25x25x97 grid over 1000 steps - and take substantially longer in a debug build.

This validates and normalizes the configuration, writes the generated control
artifacts, launches the solver, then launches the postprocessor.

@section p41_check_sec 6. Confirm it worked

Look for:

- `<run.config>/` - generated runtime control artifacts, including `<run_id>.control`,
- `<run.runtime_logs>/` - solver runtime logs,
- `<run.scheduler>/` - the solver and postprocessor **stream** logs,
- `<run.solver_output>/checkpoints/` - committed checkpoint bundles,
- `<run.visualization>/standard_analysis/` - the VTK output.

The quickstart writes exactly one visualization file:

```text
<run>/visualization/standard_analysis/eulerian_data_00020.vts
```

A quick visual check in ParaView: open that file, add a `Slice`, and colour by
`Ucat_nodal`. **@subpage 04_Visualization_Tutorial** covers this properly.

@note The post recipe targets a specific step. `quickstart_standard_analysis.yml`
targets step 20 to match the quickstart run; the full `standard_analysis.yml` targets
step 1000. A recipe asking for a step the run never reached is skipped with
"first requested source step ... is incomplete", and no visualization is written.

If something failed, **@subpage 39_Common_Fatal_Errors** maps common first-run
failures to corrective actions.

@section p41_path_sec 7. Recommended read order

Once the run above succeeds:

1. **@subpage 02_Tutorial_Programmatic_Grid** - the same path with full
   explanation and artifact inspection.
2. **@subpage 03_Tutorial_File-Based_Grid** - using an externally supplied grid.
3. **@subpage 04_Visualization_Tutorial** - getting meaningful pictures out.
4. **@subpage 05_The_Conductor_Script** - the command and option model.
5. **@subpage 49_Workflow_Recipes_and_Config_Cookbook** - recombining profiles.

@section p41_outputs_sec 8. What you should be able to do next

- Build `bin/simulator` and `bin/postprocessor`.
- Generate valid runtime control artifacts from YAML.
- Execute `run --solve --post-process` locally.
- Inspect VTK output in ParaView.
- Recombine `case.yml`, `solver.yml`, `monitor.yml`, and `post.yml` intentionally.
- Map a first-run failure to a corrective action.

@section p41_next_sec 9. Where to go next

- Authoring your own runs: **@subpage 42_User_Guide_Index**
- Grid generation and `grid.mode: grid_gen`: **@subpage 48_Grid_Generator_Guide**
- Running on a cluster, restarting, and reusing runs: **@subpage 52_Run_Lifecycle_Guide**
- CI and smoke-test contract: **@subpage 40_Testing_and_Quality_Guide**
- Full page index: **@subpage Documentation_Map**
