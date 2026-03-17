@page 02_Tutorial_Programmatic_Grid Tutorial: Your First Simulation (Flat Channel)

@anchor _Tutorial_Programmatic_Grid

This is the complete first simulation walkthrough for PICurv using the `flat_channel` template.
This is the canonical first-simulation walkthrough with both commands and reasoning.

@tableofcontents

@section p02_goal_sec 1. Tutorial Goal

By the end of this tutorial you will have:

- generated a run directory from YAML,
- executed a solver step pipeline,
- produced postprocessed VTK output,
- understood where each configuration file is used.

@section p02_init_sec 2. Initialize the Study Directory

From repository root:

```bash
./bin/picurv init flat_channel --dest my_first_run
```

Expected files:

```text
my_first_run/
|- flat_channel.yml
|- Imp-MG-Standard.yml
|- Standard_Output.yml
`- standard_analysis.yml
```

@section p02_roles_sec 3. Understand File Roles

- `flat_channel.yml`
  case physics, scaling, grid, BC handlers, model flags.
- `Imp-MG-Standard.yml`
  momentum strategy, tolerances, pressure multigrid controls.
- `Standard_Output.yml`
  logging verbosity, profiling, output directories, monitor flags.
- `standard_analysis.yml`
  postprocessing pipeline tasks and VTK/statistics output settings.

@section p02_validate_sec 4. Validate Inputs Before Launch

```bash
./bin/picurv validate \
  --case my_first_run/flat_channel.yml \
  --solver my_first_run/Imp-MG-Standard.yml \
  --monitor my_first_run/Standard_Output.yml \
  --post my_first_run/standard_analysis.yml
```

Validation confirms the YAML-to-runtime contract and catches mis-typed keys before execution.

@section p02_run_sec 5. Run Solver and Postprocessor

```bash
./bin/picurv run \
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
4. `bin/simulator` is launched,
5. `bin/postprocessor` is launched.

@section p02_artifacts_sec 6. Inspect Generated Artifacts

Typical output structure:

```text
runs/
`- flat_channel_YYYYMMDD-HHMMSS/
   |- config/
   |- logs/
   |- output/
   `- viz/
```

Interpretation:

- `config/`: exact runtime artifact snapshot for reproducibility,
- `logs/`: function-level and solver monitor traces,
- `output/`: solver field outputs,
- `viz/`: postprocessed VTK files for ParaView.

@section p02_checks_sec 7. First Validation Checks

After run completion, verify:

- no fatal errors in `runs/<run_id>/logs/`,
- non-empty VTK outputs in visualization directory,
- expected time sequence naming (`Field_*.vts`, optionally `Particle_*.vtp`).

If results are missing, check `post.yml` output toggles and directory settings first.

@section p02_viz_sec 8. Visualize in ParaView

1. Open `Field_*.vts` time series.
2. Click `Apply`.
3. Add `Slice` filter.
4. Color by `Ucat_nodal` magnitude or streamwise component.
5. Optionally add `Stream Tracer` from inlet region.

Expected behavior: smooth channel flow with profile development consistent with setup.

@section p02_troubleshooting_sec 9. Common First-Run Issues

- Missing binary: rerun `./scripts/picurv build` and check PETSc vars.
- Validation failure: compare key paths against examples and **@subpage 14_Config_Contract**.
- Empty post output: verify `--post-process` flag and `post.yml` output settings.

Mapped error guide:

- **@subpage 39_Common_Fatal_Errors**

@section p02_next_steps_sec 10. Next Steps

- File-grid workflow: **@subpage 03_Tutorial_File-Based_Grid**
- Visualization details: **@subpage 04_Visualization_Tutorial**
- Operational recipes: **@subpage 11_User_How_To_Guides**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Tutorial: Your First Simulation (Flat Channel)** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

