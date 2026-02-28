@page 48_Grid_Generator_Guide Grid Generator Guide: scripts/grid.gen

`scripts/grid.gen` is the standalone structured-grid generation utility used by the `grid.mode: grid_gen` workflow.
This page documents what it can generate directly, how its config files work, and how `pic.flow` wraps it.

@tableofcontents

@section role_sec 1. What `grid.gen` Is For

Use `grid.gen` when you want:

- a reproducible Python-side grid-generation step before solver launch,
- generated `.picgrid` geometry instead of C-side `programmatic_c` geometry,
- reusable generator configs that can be versioned and shared across cases,
- optional quick mesh diagnostics (`.info`) and visualization output (`.vts`).

This is distinct from:

- `grid.mode: programmatic_c`: geometry is generated inside the C runtime,
- `grid.mode: file`: an existing `.picgrid` file is staged and used directly.

@section invocation_sec 2. Direct CLI Usage

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
```

`grid.gen` accepts both:

- a config file (`--config`) for reusable defaults,
- direct CLI flags for overrides.

CLI values override config-file values.

@section types_sec 3. Supported Grid Types

Current subcommands:

- `cpipe`: bent pipe with a rectangular cross-section
- `pipe`: bent pipe with a circular cross-section (O-grid style)
- `warp`: generic warped Cartesian block

These map to different geometric parameter sets inside `grid.gen`.

@section config_sec 4. Config File Model

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

@section outputs_sec 5. Outputs

Depending on settings, `grid.gen` can emit:

- `.picgrid`: solver-ingestible structured grid file
- `.info`: grid quality/statistics summary
- `.vts`: visualization helper for ParaView

`grid.gen` also prints mesh-quality information to the console when enabled, including normalized warpage checks.

@section case_sec 6. Using `grid.gen` Through `case.yml`

`pic.flow` wraps the generator through:

```yaml
grid:
  mode: grid_gen
  generator:
    config_file: config/grids/coarse_square_tube_curved.cfg
    grid_type: cpipe
    cli_args: ["--ncells-i", "96", "--ncells-j", "96"]
```

Behavior in `pic.flow`:

1. validate generator config path and wrapper settings,
2. call `scripts/grid.gen`,
3. validate and non-dimensionalize the generated grid,
4. stage the runtime `grid.run` artifact for C-side ingestion.

This means `grid_gen` remains a Python-side preprocessing workflow even though the final solver still receives a normalized file grid.

@section choose_sec 7. When To Use Which Grid Path

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

@section related_sec 8. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 17_Workflow_Extensibility**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
