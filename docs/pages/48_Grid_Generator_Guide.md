@page 48_Grid_Generator_Guide Grid Generator Guide: scripts/grid.gen

@anchor _Grid_Generator_Guide

`scripts/grid.gen` is the standalone structured-grid utility used by the `grid.mode: grid_gen` workflow and legacy grid conversion paths.
This page documents what it can generate or convert directly, how its config files work, and how `picurv` wraps it.

@tableofcontents

@section p48_role_sec 1. What `grid.gen` Is For

Use `grid.gen` when you want:

- a reproducible Python-side grid-generation step before solver launch,
- generated `.picgrid` geometry instead of C-side `programmatic_c` geometry,
- reusable generator configs that can be versioned and shared across cases,
- optional quick mesh diagnostics (`.info`) and visualization output (`.vts`).
- conversion of legacy headerless 1D-axis grid payloads into canonical PICGRID files.

This is distinct from:

- `grid.mode: programmatic_c`: geometry is generated inside the C runtime,
- `grid.mode: file`: an existing `.picgrid` file is staged and used directly.

@section p48_invocation_sec 2. Direct CLI Usage

General form:

```bash
python3 scripts/grid.gen [--config <cfgfile>] <grid_type> [grid-specific options...]
```

Examples from the current script interface:

```bash
python3 scripts/grid.gen --config config/grids/coarse_square_tube_curved.cfg pipe --ncells-path 200
python3 scripts/grid.gen cpipe --orientation xz --vts cpipe_xz.vts
python3 scripts/grid.gen cpipe --ncells-i 32 --ncells-j 32 --stretch-i 2.5 --stretch-j 2.5
python3 scripts/grid.gen cpipe --ncells-i 32 --ncells-j 32 --first-cell-i-start 0.01 --first-cell-j-start 0.01
python3 scripts/grid.gen legacy1d --input legacy_flat.grid --output legacy_flat.picgrid --no-write-vtk
```

`grid.gen` accepts both:

- a config file (`--config`) for reusable defaults,
- direct CLI flags for overrides.

CLI values override config-file values.

@subsection p48_option_matrix_ssec 2.1 Option Families

Common options (all grid types):

- `--config <file>`
- `--output <picgrid_path>`
- `--vts <vtk_path>`
- `--stats-file <info_path>`
- `--origin X Y Z`
- `--show-stats` / `--no-show-stats`
- `--write-vtk` / `--no-write-vtk`

Global stretching controls:

- `--stretch-i`, `--stretch-j`, `--stretch-k`
- first-cell sizing controls:
  - `--first-cell-i-start`, `--first-cell-i-end`
  - `--first-cell-j-start`, `--first-cell-j-end`
  - `--first-cell-k-start`, `--first-cell-k-end`

Grid-type-specific options:

- `cpipe`: `--ncells-i`, `--ncells-j`, `--ncells-k`, `--side-lengths`, `--rc-factor`, `--straight-factor`, `--bend-angle`, `--orientation`
- `pipe`: `--ncells-phi`, `--ncells-r`, `--ncells-path`, `--diameter`, `--pinhole-factor`, `--rc-factor`, `--straight-factor`, `--bend-angle`, `--orientation`
- `warp`: `--ncells-i`, `--ncells-j`, `--ncells-k`, `--bounds-x`, `--bounds-y`, `--bounds-z`, `--amp-A`, `--amp-B`, `--amp-C`
- `legacy1d`: `--input`, `--axis-columns`, `--strict-trailing` / `--allow-trailing`

Notes for CFD users:

- `ncells-*` are cell counts at input.
- PICGRID header dimensions are node counts (`cells + 1`) after conversion.
- this conversion behavior is regression-tested in `tests/test_cli_smoke.py`.

@section p48_types_sec 3. Supported Grid Types

Current subcommands:

- `cpipe`: bent pipe with a rectangular cross-section
- `pipe`: bent pipe with a circular cross-section (O-grid style)
- `warp`: generic warped Cartesian block
- `legacy1d`: legacy 1D-axis payload converter (headerless block + dims + x/y/z axis lists -> canonical PICGRID)

These map to different geometric parameter sets inside `grid.gen`.

@section p48_config_sec 4. Config File Model

Reusable config files live naturally in `config/grids/`, but can also live beside a study for reproducibility snapshots.

Current shared example:

- `config/grids/coarse_square_tube_curved.cfg`

The config file is INI-style:

- the section name typically matches the grid type (`[cpipe]`, `[pipe]`, `[warp]`)
- values in that section provide defaults
- CLI arguments can override selected keys at run time

Typical contents include:

- cell counts,
- geometric size/bounds,
- bend/orientation settings,
- stretching controls,
- output toggles and optional filenames.

@section p48_outputs_sec 5. Outputs

Depending on settings, `grid.gen` can emit:

- `.picgrid`: solver-ingestible structured grid file
- `.info`: grid quality/statistics summary
- `.vts`: visualization helper for ParaView

`grid.gen` also prints mesh-quality information to the console when enabled, including normalized warpage checks.

@section p48_case_sec 6. Using `grid.gen` Through `case.yml`

`picurv` wraps the generator through:

```yaml
grid:
  mode: grid_gen
  generator:
    config_file: config/grids/coarse_square_tube_curved.cfg
    grid_type: cpipe
    cli_args: ["--ncells-i", "96", "--ncells-j", "96"]
```

Current contract notes:

- `grid.generator.config_file` is required today.
- `picurv` does not generate a temporary `grid.cfg` for `grid.gen`.
- `grid.gen` input dimensions (`ncells_*`, `--ncells-*`) are cell counts.
- `grid.gen` converts those to node counts before writing the `.picgrid` header.

Behavior in `picurv`:

1. validate generator config path and wrapper settings,
2. call `scripts/grid.gen`,
3. validate and non-dimensionalize the generated grid,
4. stage the runtime `grid.run` artifact for C-side ingestion.

This means `grid_gen` remains a Python-side preprocessing workflow even though the final solver still receives a normalized file grid.

@subsection p48_file_legacy_ssec 6.1 Legacy File Conversion Through `case.yml`

For `grid.mode: file`, `picurv` can call `grid.gen legacy1d` before normal grid staging:

```yaml
grid:
  mode: file
  source_file: legacy_flat.grid
  legacy_conversion:
    enabled: true
    format: legacy1d
    axis_columns: [0, 1, 2]
    strict_trailing: true
```

Execution sequence:

1. run `grid.gen legacy1d --input <source_file> --output <run>/config/grid.converted.picgrid`
2. validate converted PICGRID payload
3. non-dimensionalize to `grid.run`
4. pass `-grid_file <run>/config/grid.run` to the C runtime

@section p48_choose_sec 7. When To Use Which Grid Path

Use `programmatic_c` when:

- you want the simplest structured box-style runtime setup,
- you do not need an external generated mesh artifact,
- you want direct C-side control of the structured grid box inputs.

Use `file` when:

- you already have a stable `.picgrid`,
- geometry should be fixed and explicit,
- you want no generator step during launch.

Use `grid_gen` when:

- you want generated geometry but still want file-based staging,
- you want a reusable mesh recipe with controlled overrides,
- you want generated `.info`/`.vts` outputs alongside the mesh.

@section p48_related_sec 8. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 17_Workflow_Extensibility**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Grid Generator Guide: scripts/grid.gen** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
