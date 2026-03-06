@page 03_Tutorial_File-Based_Grid Tutorial: Using a File-Based Grid (Bent Channel)

@anchor _Tutorial_File-Based_Grid

This tutorial demonstrates `grid.mode: file` using the `bent_channel` template.

@tableofcontents

@section p03_init_sec 1. Initialize a Study

```bash
./bin/picurv init bent_channel --dest my_bent_channel_run
```

Expected files:

```text
my_bent_channel_run/
|- bent_channel.yml
|- Imp-MG-Standard.yml
|- Standard_Output.yml
|- standard_analysis.yml
|- bent_channel_coarse.picgrid
|- bent_channel_coarse.vts
`- bent_channel_coarse.info
```

@section p03_config_sec 2. File-Grid Configuration

In `bent_channel.yml`:

```yaml
grid:
  mode: file
  source_file: bent_channel_coarse.picgrid
```

Behavior:
- `picurv` validates the source grid file exists.
- Coordinates are non-dimensionalized before C execution using `properties.scaling.length_ref`.
- Generated normalized grid is staged into run config artifacts.

@section p03_run_sec 3. Run the Case

```bash
./bin/picurv run \
    --case my_bent_channel_run/bent_channel.yml \
    --solver my_bent_channel_run/Imp-MG-Standard.yml \
    --monitor my_bent_channel_run/Standard_Output.yml \
    --post my_bent_channel_run/standard_analysis.yml \
    -n 4 --solve --post-process
```

@section p03_viz_sec 4. Visualize Results

1. Open `runs/<run_id>/viz/Field_*.vts`.
2. Use `Stream Tracer` seeded near inlet.
3. Color by `Ucat_nodal` magnitude.

This provides a fast sanity check of bend-flow behavior.

@section p03_next_steps_sec 5. Next Steps

Proceed to **@subpage 04_Visualization_Tutorial**.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Tutorial: Using a File-Based Grid (Bent Channel)** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

