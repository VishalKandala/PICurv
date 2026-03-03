@page 03_Tutorial_File-Based_Grid Tutorial: Using a File-Based Grid (Bent Channel)

This tutorial demonstrates `grid.mode: file` using the `bent_channel` template.

@tableofcontents

@section init_sec 1. Initialize a Study

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

@section config_sec 2. File-Grid Configuration

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

@section run_sec 3. Run the Case

```bash
./bin/picurv run \
    --case my_bent_channel_run/bent_channel.yml \
    --solver my_bent_channel_run/Imp-MG-Standard.yml \
    --monitor my_bent_channel_run/Standard_Output.yml \
    --post my_bent_channel_run/standard_analysis.yml \
    -n 4 --solve --post-process
```

@section viz_sec 4. Visualize Results

1. Open `runs/<run_id>/viz/Field_*.vts`.
2. Use `Stream Tracer` seeded near inlet.
3. Color by `Ucat_nodal` magnitude.

This provides a fast sanity check of bend-flow behavior.

@section next_steps_sec 5. Next Steps

Proceed to **@subpage 04_Visualization_Tutorial**.
