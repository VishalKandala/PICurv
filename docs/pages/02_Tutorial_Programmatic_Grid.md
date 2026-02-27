@page 02_Tutorial_Programmatic_Grid Tutorial: Your First Simulation (Flat Channel)

This is the complete first simulation walkthrough for PICurv using the `flat_channel` template.
Compared to **@subpage 38_Start_Here_10_Minutes**, this page explains what each step does and how to verify outcomes.

@tableofcontents

@section goal_sec 1. Tutorial Goal

By the end of this tutorial you will have:

- generated a run directory from YAML,
- executed a solver step pipeline,
- produced postprocessed VTK output,
- understood where each configuration file is used.

@section init_sec 2. Initialize the Study Directory

From repository root:

```bash
./scripts/pic.flow init flat_channel --dest my_first_run
```

Expected files:

```text
my_first_run/
|- flat_channel.yml
|- Imp-MG-Standard.yml
|- Standard_Output.yml
`- standard_analysis.yml
```

@section roles_sec 3. Understand File Roles

- `flat_channel.yml`
  case physics, scaling, grid, BC handlers, model flags.
- `Imp-MG-Standard.yml`
  momentum strategy, tolerances, pressure multigrid controls.
- `Standard_Output.yml`
  logging verbosity, profiling, output directories, monitor flags.
- `standard_analysis.yml`
  postprocessing pipeline tasks and VTK/statistics output settings.

@section validate_sec 4. Validate Inputs Before Launch

```bash
./scripts/pic.flow validate \
  --case my_first_run/flat_channel.yml \
  --solver my_first_run/Imp-MG-Standard.yml \
  --monitor my_first_run/Standard_Output.yml \
  --post my_first_run/standard_analysis.yml
```

Validation confirms the YAML-to-runtime contract and catches mis-typed keys before execution.

@section run_sec 5. Run Solver and Postprocessor

```bash
./scripts/pic.flow run \
  --case my_first_run/flat_channel.yml \
  --solver my_first_run/Imp-MG-Standard.yml \
  --monitor my_first_run/Standard_Output.yml \
  --post my_first_run/standard_analysis.yml \
  -n 4 --solve --post-process
```

What happens internally:

1. all input files are loaded and validated,
2. run ID directory is created under `runs/`,
3. C-consumed artifacts are generated (`*.control`, `bcs*.run`, `whitelist.run`, `profile.run`, `post.run`),
4. `bin/picsolver` is launched,
5. `bin/postprocessor` is launched.

@section artifacts_sec 6. Inspect Generated Artifacts

Typical output structure:

```text
runs/
`- flat_channel_YYYYMMDD-HHMMSS/
   |- config/
   |- logs/
   |- results/
   `- viz/
```

Interpretation:

- `config/`: exact runtime artifact snapshot for reproducibility,
- `logs/`: function-level and solver monitor traces,
- `results/`: solver field outputs,
- `viz/`: postprocessed VTK files for ParaView.

@section checks_sec 7. First Validation Checks

After run completion, verify:

- no fatal errors in `runs/<run_id>/logs/`,
- non-empty VTK outputs in visualization directory,
- expected time sequence naming (`Field_*.vts`, optionally `Particle_*.vtp`).

If results are missing, check `post.yml` output toggles and directory settings first.

@section viz_sec 8. Visualize in ParaView

1. Open `Field_*.vts` time series.
2. Click `Apply`.
3. Add `Slice` filter.
4. Color by `Ucat_nodal` magnitude or streamwise component.
5. Optionally add `Stream Tracer` from inlet region.

Expected behavior: smooth channel flow with profile development consistent with setup.

@section troubleshooting_sec 9. Common First-Run Issues

- Missing binary: rerun `./scripts/pic.flow build` and check PETSc vars.
- Validation failure: compare key paths against examples and **@subpage 14_Config_Contract**.
- Empty post output: verify `--post-process` flag and `post.yml` output settings.

Mapped error guide:

- **@subpage 39_Common_Fatal_Errors**

@section next_steps_sec 10. Next Steps

- File-grid workflow: **@subpage 03_Tutorial_File-Based_Grid**
- Visualization details: **@subpage 04_Visualization_Tutorial**
- Operational recipes: **@subpage 11_User_How_To_Guides**
