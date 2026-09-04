@page 03_Tutorial_File-Based_Grid Tutorial: Using a File-Based Grid (Search Robustness)

@anchor _Tutorial_File-Based_Grid

This tutorial demonstrates `grid.mode: file` using the `search_robustness` template, which
carries the only `.picgrid` the examples still ship.

Use `file` mode when the geometry already exists as a grid file and should stay fixed and
explicit - a mesh from another tool, or one whose exact node positions are the thing under
test. When the geometry is something the bundled generator can express, `grid_gen` says it
in a few lines instead of a few megabytes; `bent_channel` moved that way once `sweep`
could express a square duct through a quarter turn. See **@subpage 48_Grid_Generator_Guide**.

@tableofcontents

@section p03_init_sec 1. Initialize a Study

```bash
./bin/picurv init search_robustness --dest my_file_grid_run
```

Among the files it places:

```text
my_file_grid_run/
|- search_robustness_curvilinear.yml
|- Imp-MG-Standard.yml
|- Search_Robustness_Output.yml
|- search_robustness_analysis.yml
`- bent_channel_coarse.picgrid
```

The grid is a unit square duct through a 90-degree bend, 21 x 21 x 145 nodes.

@section p03_config_sec 2. File-Grid Configuration

In `search_robustness_curvilinear.yml`:

```yaml
grid:
  mode: file
  source_file: bent_channel_coarse.picgrid
```

Behavior:
- `picurv` validates the source grid file exists.
- Coordinates are non-dimensionalized before C execution using `properties.scaling.length_ref`.
- Generated normalized grid is staged into run config artifacts.

For geometric-periodic file grids, the runtime additionally verifies that each
paired surface matches pointwise under one nonzero constant translation and
that the periodic axis has at least four physical nodes.

@section p03_run_sec 3. Run the Case

```bash
./bin/picurv run \
    --case my_file_grid_run/config/case.yml \
    --solver my_file_grid_run/config/solver.yml \
    --monitor my_file_grid_run/config/monitor.yml \
    --post my_file_grid_run/config/post.yml \
    -n 4 --solve --post-process
```

@section p03_viz_sec 4. Visualize Results

1. Open `<run.visualization>/<recipe_id>/Field_*.vts`.
2. Use `Stream Tracer` seeded near inlet.
3. Color by `Ucat_nodal` magnitude.

This provides a fast sanity check of bend-flow behavior.

@section p03_next_steps_sec 5. Next Steps

Proceed to **@subpage 04_Visualization_Tutorial**.
