@page 02_Tutorial_Programmatic_Grid Tutorial: Your First Simulation (Flat Channel)

This tutorial runs the `flat_channel` example end-to-end with the current YAML pipeline.

@tableofcontents

@section init_sec 1. Initialize a Study Directory

From project root:

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

@section config_sec 2. What Each File Does

- `flat_channel.yml`: case physics, grid, BCs, runtime control
- `Imp-MG-Standard.yml`: solver strategy and numerics
- `Standard_Output.yml`: I/O, logging, profiling, monitoring
- `standard_analysis.yml`: post-processing pipeline recipe

@section run_sec 3. Run Solver + Postprocessor

```bash
./scripts/pic.flow run \
    --case my_first_run/flat_channel.yml \
    --solver my_first_run/Imp-MG-Standard.yml \
    --monitor my_first_run/Standard_Output.yml \
    --post my_first_run/standard_analysis.yml \
    -n 4 --solve --post-process
```

What `pic.flow` does:

1. Validates YAML files.
2. Creates `runs/<run_id>/`.
3. Generates C-consumed config artifacts:
   - `<run_id>.control`
   - `bcs.run` (or `bcs_block*.run`)
   - `whitelist.run`
   - `profile.run`
   - `post.run`
4. Launches `bin/picsolver`.
5. Launches `bin/postprocessor` (because `--post-process` was set).

@section output_sec 4. Inspect Output Layout

```text
runs/
`- flat_channel_YYYYMMDD-HHMMSS/
   |- config/
   |- logs/
   |- results/
   `- viz/
```

- `config/`: generated runtime artifacts and copied input YAMLs
- `results/`: solver binary field dumps
- `viz/`: postprocessed `.vts` / `.vtp`

@section viz_sec 5. Visualize in ParaView

1. Open `runs/<run_id>/viz/Field_*.vts` as a time series.
2. Click `Apply`.
3. Add a `Slice` filter.
4. Color by `Ucat_nodal` and select the flow-direction component.

You should observe the expected Poiseuille profile development.

@section next_steps_sec 6. Next Steps

Proceed to **@subpage 03_Tutorial_File-Based_Grid** for file-based grid input workflow.
