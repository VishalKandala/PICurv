# Profile Generation Guide

This directory documents reusable inlet-profile patterns for the data sources of a
`prescribed_velocity` inlet. The actual profile choice lives in `case.yml` because the
BC face and grid dimensions determine the required `PICSLICE` shape.

## Supported Patterns

File-backed profile:

```yaml
boundary_conditions:
  - face: "-Zeta"
    type: INLET
    handler: prescribed_velocity
    source:
      type: file
      path: profiles/inlet.picslice
```

Field-sliced profile from an old PICurv `ufield*.dat`:

```yaml
boundary_conditions:
  - face: "-Zeta"
    type: INLET
    handler: prescribed_velocity
    source:
      type: field_slice
      field_file: ../old_run/output/checkpoints/step_000000010000/eulerian/block_0000/Ucat.dat
      grid_file: ../old_run/inputs/grid/grid.run
      source_case: ../old_run/config/case.yml
      slice:
        face: "+Zeta"
        orientation: opposite
```

`field_slice` reads Cartesian `Ucat`, computes geometry normals from the old
and target `PICGRID` files, projects to positive normal speed magnitudes, and
writes the same canonical `PICSLICE` format as the other source types. V1
requires exact slice dimension matches; it does not interpolate or resample.

`picurv run --solve` automatically materializes sliced profiles under
`<run.config>/`, writes `profile.info`, stages a nondimensional solver
copy, and passes the staged `source_file` to the C runtime.

Use `picurv precompute --case case.yml --only inlet-profiles` when you want to
generate and inspect profile artifacts before launching the solver; they are
published in the workspace asset store and reused by the run that follows.

## Standalone Engine

`generators/profile.gen` is the low-level generator:

```bash
python3 generators/profile.gen field-slice \
  --output inlet_from_old_run.picslice \
  --field-file old_run/output/checkpoints/step_000000010000/eulerian/block_0000/Ucat.dat \
  --source-grid old_run/inputs/grid/grid.run \
  --target-grid new_run/inputs/grid/grid.run \
  --slice-face=+Zeta \
  --target-face=-Zeta \
  --orientation opposite \
  --velocity-scale 1.0
```

Prefer the `picurv` conductor for normal cases, because it derives dimensions
from `case.yml`, `grid.mode`, block index, and BC face.

## Design Intent

Profile sources intentionally converge on one artifact contract: a dimensional
canonical `PICSLICE` that `picurv` can validate, nondimensionalize, stage, and
wire into `bcs.run` as `source_file`. This keeps the C handler source-agnostic:
hand-authored files and old-field slices both reach the
solver through the same path.

`profile.gen` is the standalone engine for profile-producing operations.
`picurv precompute` is the orchestration command for materializing configured
grid/profile artifacts before a solve, and `picurv run --solve` performs the
same materialization automatically when a sliced source is selected.
Future derived-profile sources should extend
`profile.gen` first, then be exposed as a `prescribed_velocity` data source
using the same artifact pipeline. A profile with a closed form needs neither: write it as
an expression source.
