@page 48_Grid_Generator_Guide Grid Generator Guide: generators/grid.gen

@anchor _Grid_Generator_Guide

`generators/grid.gen` is the standalone structured-grid utility used by the `grid.mode: grid_gen` workflow and legacy grid conversion paths.
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
python3 generators/grid.gen [--config <cfgfile>] <grid_type> [grid-specific options...]
```

Examples from the current script interface:

```bash
python3 generators/grid.gen --config config/grids/coarse_square_tube_curved.cfg pipe --ncells-path 200
python3 generators/grid.gen cpipe --orientation xz --vts cpipe_xz.vts
python3 generators/grid.gen cpipe --ncells-i 32 --ncells-j 32 --stretch-i 2.5 --stretch-j 2.5
python3 generators/grid.gen cpipe --ncells-i 32 --ncells-j 32 --first-cell-i-start 0.01 --first-cell-j-start 0.01
python3 generators/grid.gen legacy1d --input legacy_flat.grid --output legacy_flat.picgrid --no-write-vtk
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
2. call `generators/grid.gen`,
3. validate and non-dimensionalize the generated grid,
4. stage the runtime `grid.run` artifact for C-side ingestion.

This means `grid_gen` remains a Python-side preprocessing workflow even though the final solver still receives a normalized file grid.
Fresh solve runs and `--restart-from` runs stage a grid in the target run directory.
For `--continue` solve runs, an existing `<run.config>/grid.run` is reused when present
so the generator is not called again. Post-only runs read the existing control
file and do not regenerate grids.

@subsection p48_promote_generated_ssec 6.1 Promoting a Generated Grid to File Mode

Use `grid_gen` when you want PICurv to run `grid.gen` as part of staging a case.
Use `file` mode for later runs when you want to reuse a generated `.picgrid`
without invoking the generator again.

A practical workflow is:

1. stage the original case with `grid.mode: grid_gen` and `--no-submit`,
2. inspect the generated PICGRID artifact under the run config directory,
3. copy or reference that `.picgrid` from a later case,
4. switch the later case to `grid.mode: file`.

By default the generated artifact is:

```text
runs/<run_id>/config/grid.generated.picgrid
```

The solver does not ingest that dimensional file directly. PICurv validates and
non-dimensionalizes it into:

```text
runs/<run_id>/config/grid.run
```

and the generated solver control file points the C runtime at `<run.config>/grid.run`.
When reusing the generated mesh in another case, point `grid.source_file` at the
`.picgrid` artifact, not at `grid.run`.

Example follow-up case:

```yaml
grid:
  mode: file
  source_file: ../runs/<run_id>/config/grid.generated.picgrid
```

@subsection p48_file_legacy_ssec 6.2 Legacy File Conversion Through `case.yml`

For `grid.mode: file`, `picurv` can call `grid.gen legacy1d` before normal grid staging:

```yaml
grid:
  mode: file
  source_file: legacy_flat.grid
  legacy_conversion:
    enabled: true
    format: legacy1d            # or column_text; aliases: legacy_1d, les_flat_1d, les-flat-1d
    output_file: null           # optional: override the generated .picgrid output path
    script: null                # optional: override the conversion script (default: generators/grid.gen)
    axis_columns: [0, 1, 2]    # preferred source columns for X/Y/Z axis rows
    strict_trailing: true
    cli_args: []                # additional raw tokens forwarded to the conversion script
```

Execution sequence:

1. run `grid.gen legacy1d --input <source_file> --output <run.config>/grid.converted.picgrid`
2. validate converted PICGRID payload
3. non-dimensionalize to `grid.run`
4. pass `-grid_file <run.config>/grid.run` to the C runtime

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

For a generated grid used with geometric-periodic BCs, each opposite periodic
surface must be a pointwise translated copy of its partner. The translation
must be nonzero and constant over the entire surface; rotating, shearing, and
nonmatching surface pairs are rejected at runtime.

@section p48_cap_geom_sec 7.1 Grid Generator Geometry Entries

@htmlinclude generated/capability_inventory_grid_generator_type.html

@subsection p48_cap_geom_warp_sub warp

@anchor p48_cap_geom_warp

**Identity.** `grid.gen warp`, or `grid.generator.grid_type: warp`.

**What it does.** Generates a Cartesian block optionally deformed by three independent sinusoidal amplitudes. With all amplitudes at zero it is a plain stretched box, which is how the shipped periodic cases use it.

**When to choose it.** The default choice for channels, ducts, and any box-like domain generated rather than declared - it is the geometry every shipped `grid_gen` case selects. Turn the amplitudes up only to exercise curvilinear machinery on a deformation you know analytically.

**Parameters it owns.** `--ncells-i`, `--ncells-j`, `--ncells-k`, `--bounds-x`, `--bounds-y`, `--bounds-z`, and the deformation amplitudes `--amp-A`, `--amp-B`, `--amp-C`.

**Interactions.** Reads defaults from the `[warp]` config section - see `config/grids/` for the plane-channel and square-duct configurations. With zero amplitudes the result is equivalent to a `programmatic_c` box, and the choice between them is whether the geometry lives in a generator config or in the case file.

**Diagnostics.** `grid.gen` prints mesh-quality information including normalized warpage checks; `--vts` writes a viewable mesh for inspection before a solve.

**Evidence.** Production exercised - `examples/periodic_test` driven-channel and driven-duct cases all generate their meshes with this type.

**Limitations.** The deformation is a fixed three-amplitude sinusoid; arbitrary warping is not expressible, and no case establishes which amplitudes stay well-conditioned.

@subsection p48_cap_geom_cpipe_sub cpipe

@anchor p48_cap_geom_cpipe

**Identity.** `grid.gen cpipe`, or `grid.generator.grid_type: cpipe`.

**What it does.** Generates a bent pipe with a **rectangular** cross-section as a single curvilinear block. The bend is parameterised by a radius factor and an angle rather than drawn.

**When to choose it.** For bent-duct geometries with a square or rectangular cross-section. Choose `pipe` instead when the cross-section must be circular, and `warp` for anything unbent.

**Parameters it owns.** `--ncells-i`, `--ncells-j`, `--ncells-k`, `--side-lengths`, `--rc-factor`, `--straight-factor`, `--bend-angle`, `--orientation`, plus the stretching and first-cell controls in **@ref p48_invocation_sec**.

**Interactions.** Reads defaults from the `[cpipe]` section - `config/grids/coarse_square_tube_curved.cfg` ships one. `--orientation` fixes which pair of axes the bend lies in, which must match the periodic-axis choice if the case is periodic.

**Diagnostics.** Mesh-quality output reports normalized warpage, which is where a too-tight bend radius shows up first.

**Evidence.** Implemented only. A generator config ships, but no shipped case selects this type, so nothing here is exercised end to end.

**Limitations.** One block only, and the bend is a single constant-radius arc; compound or variable-radius bends are not expressible. Experimental until a case exercises it.

@subsection p48_cap_geom_pipe_sub pipe

@anchor p48_cap_geom_pipe

**Identity.** `grid.gen pipe`, or `grid.generator.grid_type: pipe`.

**What it does.** Generates a bent pipe with a **circular** cross-section, meshed O-grid style so the centre is not a singular point.

**When to choose it.** For round-pipe geometries, where a rectangular cross-section would change the physics rather than merely the mesh. The O-grid topology is what makes a circular section tractable on a structured block.

**Parameters it owns.** `--ncells-phi`, `--ncells-r`, `--ncells-path`, `--diameter`, `--pinhole-factor`, `--rc-factor`, `--straight-factor`, `--bend-angle`, `--orientation`.

**Interactions.** Reads defaults from the `[pipe]` config section. `--pinhole-factor` controls the O-grid core and interacts with near-axis cell quality; it is the first control to revisit when metric diagnostics complain near the centreline.

**Diagnostics.** Same mesh-quality reporting as the other geometries. Near-axis cells are where warpage warnings appear first.

**Evidence.** Implemented only. No shipped case or study selects it.

**Limitations.** O-grid cell quality degrades toward the axis, and nothing establishes what `--pinhole-factor` should be for a given Reynolds number. Experimental until a case exercises it.

@section p48_related_sec 8. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 17_Workflow_Extensibility**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@ref p54_geometric_periodic "Periodic Boundaries and Driven Flows"**
