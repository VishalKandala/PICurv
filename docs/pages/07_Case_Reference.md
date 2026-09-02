@page 07_Case_Reference Configuration Reference: Case YAML

@anchor _Case_Reference

For the full commented template, see:

@verbinclude master_template/master_case.yml

`case.yml` defines physical setup, grid source, domain topology, and boundary conditions.
It is intentionally modular: the same `case.yml` can be paired with different `solver.yml`, `monitor.yml`, and `post.yml` profiles when the combination remains contract-compatible.

@tableofcontents

@section p07_properties_sec 1. properties

```yaml
properties:
  scaling:
    length_ref: 0.1
    velocity_ref: 1.5
  fluid:
    density: 1000.0
    viscosity: 0.001
  initial_conditions:
    mode: generated
    generator: constant
    params:
      u_physical: 1.5
      v_physical: 0.0
      w_physical: 0.0
```

Alternative IC forms:

```yaml
# Single streamwise speed
initial_conditions:
  mode: generated
  generator: streamwise_constant
  params:
    velocity_physical: 1.5
    flow_direction: "+Zeta"

# File-backed Cartesian velocity
initial_conditions:
  mode: file
  field: Ucat
  source_file: initial_conditions/velocity.dat
```

Key mappings:
- `scaling.length_ref` -> `-scaling_L_ref`
- `scaling.velocity_ref` -> `-scaling_U_ref`
- `fluid.density` and `fluid.viscosity` are used by `picurv` to compute Reynolds number -> `-ren`
- `generator: zero|constant|poiseuille|streamwise_constant` -> the corresponding built-in `-finit` mode
- `mode: file` and `generator: ic_gen` -> `-finit 4`, `-ic_field`, and staged `-ic_dir`
- `params.u_physical/v_physical/w_physical` -> `-ucont_x/-ucont_y/-ucont_z`
- `params.velocity_physical` and `params.peak_velocity_physical` -> `-ic_velocity_physical`
- `flow_direction` -> `-flow_direction <int>` (`+Xi=0,-Xi=1,+Eta=2,-Eta=3,+Zeta=4,-Zeta=5`)

For the scaling model and conversion logic, see **@subpage 19_Nondimensionalization**.
For detailed startup behavior of field initialization modes, see **@subpage 33_Initial_Conditions**.

Practical contract notes:

- `initial_conditions.mode` is `generated` or `file`.
- generated built-ins are `zero`, `constant`, `streamwise_constant`, and `poiseuille`.
- `generator: ic_gen` defaults to `generators/ic.gen`; optional `params.script` selects a compatible override.
- file-backed ICs accept one PETSc binary `Ucat` or `Ucont` vector and currently require a single-block case.
- `flow_direction` is required for curvilinear Constant and Poiseuille when no INLET face exists.
- `eulerian_field_source` and restart selection supersede `initial_conditions`.

@section p07_run_control_sec 2. run_control

```yaml
run_control:
  dt_physical: 0.0001
  start_step: 0
  total_steps: 2000
```

Mappings:
- `dt_physical` -> `-dt` (non-dimensionalized)
- `start_step` -> `-start_step`
- `total_steps` -> `-totalsteps`

Restart is handled via CLI flags rather than `case.yml` keys:
- `--restart-from <previous_run_dir>` -> `picurv` resolves the prior run's actual restart source directory and emits `-restart_dir <absolute_previous_source>`
- `--continue` -> shorthand for resuming from the most recent run of the same case

@section p07_grid_sec 3. grid

Supported modes:
- `programmatic_c`
- `file`
- `grid_gen`

Mode compatibility note:

- for normal `solve` and `load` workflows, all three grid modes are supported.
- for `solver.yml -> operation_mode.eulerian_field_source: analytical`, `TGV3D` requires `grid.mode: programmatic_c`.
- `ZERO_FLOW` and `UNIFORM_FLOW` support `grid.mode: programmatic_c` and `grid.mode: file`.

Optional global DMDA layout hints apply to all grid modes:

```yaml
grid:
  da_processors_x: 4
  da_processors_y: 2
  da_processors_z: 2
```

These are scalar global values, not per-block vectors. Legacy placement under
`grid.programmatic_settings.da_processors_*` is still accepted for compatibility,
but the shared top-level `grid.da_processors_*` form is preferred.

@subsection p07_grid_prog_ssec 3.1 mode: programmatic_c

`programmatic_settings` supports per-block lists for geometry arrays:
- `im/jm/km`
- `xMins/xMaxs`, `yMins/yMaxs`, `zMins/zMaxs`
- `rxs/rys/rzs`
- `cgrids` — per-block integer grid type: `0` = Cartesian (default), `1` = O-type curvilinear → `-cgrids`

Dimension contract:
- `im/jm/km` in YAML are cell counts.
- `picurv` converts them to node counts before emitting `-im/-jm/-km` for the C runtime.

Important constraint:
- `grid.da_processors_x/y/z` are scalar integers only (global DMDA layout). Per-block processor decomposition is not implemented.

@subsection p07_grid_file_ssec 3.2 mode: file

```yaml
grid:
  mode: file
  source_file: my_grid.picgrid
```

`picurv` validates existence and prepares normalized grid data for C-side ingestion.

Optional legacy conversion path (headerless 1D-axis payloads):

```yaml
grid:
  mode: file
  source_file: legacy_flat.grid
  legacy_conversion:
    enabled: true
    format: legacy1d            # aliases: legacy_1d, les_flat_1d, les-flat-1d; or column_text
    output_file: null           # optional: override generated .picgrid output path
    script: null                # optional: override conversion script path (default generators/grid.gen)
    axis_columns: [0, 1, 2]    # preferred source columns for X/Y/Z axis rows
    strict_trailing: true
    cli_args: []                # additional raw tokens forwarded to the conversion script
```

When enabled, `picurv` first calls `generators/grid.gen legacy1d` to create a canonical PICGRID file in the run config, then runs the normal validation/non-dimensionalization staging path.

@subsection p07_grid_gen_ssec 3.3 mode: grid_gen

```yaml
grid:
  mode: grid_gen
  generator:
    config_file: config/grids/coarse_square_tube_curved.cfg
    grid_type: cpipe
    cli_args: ["--ncells-i", "96", "--ncells-j", "96"]
```

This runs `generators/grid.gen` before solver launch and stages generated grid artifacts into the run config.
`grid.generator.config_file` is required today; `picurv` does not synthesize a temporary grid config.
`grid.gen` accepts cell-count inputs (`ncells_*` / `--ncells-*`) and writes node counts into the generated `.picgrid` header.

For direct `grid.gen` usage, generator types, and config-file structure, see **@subpage 48_Grid_Generator_Guide**.

@section p07_cap_gridmode_sec 3.4 Grid Ingestion Mode Entries

@htmlinclude generated/capability_inventory_grid_mode.html

@subsection p07_cap_gridmode_file_sub file

@anchor p07_cap_gridmode_file

**Identity.** `grid.mode: file` with `grid.source_file` naming a `.picgrid` file.

**What it does.** Reads an externally produced curvilinear mesh. `picurv` validates that the file exists and stages it into `<run.config>/`, so the run carries the mesh it was computed with rather than a path that may later change.

**When to choose it.** When the geometry came from `grid.gen`, from a mesher, or from a previous study, and you want the exact same mesh across runs. Prefer it over `grid_gen` once a mesh is settled: regenerating on every launch is a way to silently change geometry.

**Parameters it owns.** `grid.source_file` -> the staged `.picgrid` path. Global DMDA hints `grid.da_processors_x/y/z` apply here as to every mode.

**Interactions.** Mutually exclusive with `programmatic_c` and `grid_gen`. Compatible with the `ZERO_FLOW` and `UNIFORM_FLOW` analytical sources; `TGV3D` is not, because it needs the programmatic Cartesian box.

**Diagnostics.** The startup banner reports the resolved grid source. A missing or unreadable file is a fatal `picurv validate` error naming the path, not a runtime failure.

**Evidence.** Production exercised - `examples/bent_channel` runs from a staged `.picgrid`.

**Limitations.** `picurv` does not inspect the mesh beyond existence and header validity, so a geometrically wrong but well-formed file is accepted here and only shows up in metric or Jacobian diagnostics.

@subsection p07_cap_gridmode_programmatic_c_sub programmatic_c

@anchor p07_cap_gridmode_programmatic_c

**Identity.** `grid.mode: programmatic_c` with a `grid.programmatic_settings` block.

**What it does.** Builds a structured block directly in the C runtime from cell counts, extents, and stretching ratios. No mesh file exists or is read.

**When to choose it.** For Cartesian boxes and simple O-type blocks - channels, ducts, isotropic boxes - and for anything where the geometry is a handful of numbers you want in the case file rather than in a binary alongside it. It is also the only mode `TGV3D` accepts.

**Parameters it owns.** `grid.programmatic_settings.im/jm/km` (cell counts, converted to node counts before `-im/-jm/-km`), `xMins/xMaxs`, `yMins/yMaxs`, `zMins/zMaxs`, `rxs/rys/rzs`, and `cgrids` (`0` Cartesian, `1` O-type curvilinear) -> `-cgrids`.

**Interactions.** Mutually exclusive with `file` and `grid_gen`. Required by the `TGV3D` analytical source. Per-block lists are supported for the geometry arrays, but `grid.da_processors_*` stay scalar: per-block processor decomposition is not implemented.

**Diagnostics.** The startup banner reports the resolved block extents and counts. A count/extent mismatch across per-block lists is a validation error naming the offending key.

**Evidence.** Production exercised - `examples/flat_channel`; integration verified - `make unit-grid`.

**Limitations.** Restricted to what the built-in generator can express: rectangular blocks and the O-type variant. Anything bent, branched, or externally meshed needs `file` or `grid_gen`.

@subsection p07_cap_gridmode_grid_gen_sub grid_gen

@anchor p07_cap_gridmode_grid_gen

**Identity.** `grid.mode: grid_gen` with a `grid.generator` block naming a config file and a `grid_type`.

**What it does.** Runs `generators/grid.gen` before the solver launches and stages the generated `.picgrid` into `<run.config>/`. The generated mesh, not the generator invocation, is what the solver reads.

**When to choose it.** While a geometry is still being iterated - sweeping cell counts or bend angles across a study, where regenerating per case is the point. Once the mesh is settled, switch to `file` so the geometry stops depending on the generator's behaviour.

**Parameters it owns.** `grid.generator.config_file` (required; `picurv` does not synthesize one), `grid.generator.grid_type` -> **@ref p48_cap_geom_sec**, and `grid.generator.cli_args` for per-run overrides such as `--ncells-i`.

**Interactions.** Mutually exclusive with `file` and `programmatic_c`. The generator accepts cell counts (`ncells_*`) and writes node counts into the `.picgrid` header, so the count convention matches `programmatic_c` at the YAML boundary.

**Diagnostics.** Generator stdout is captured into the run's scheduler log, and the staged `.picgrid` appears under `<run.config>/`. A generator failure aborts before the solver launches.

**Evidence.** Production exercised - `examples/periodic_test` cases stage grids this way.

**Limitations.** Adds a build step to every launch, and the run's geometry depends on the generator's current behaviour rather than on a fixed artefact. `config_file` is mandatory today.

@section p07_models_sec 4. models

```yaml
models:
  domain:
    blocks: 1
  physics:
    dimensionality: "3D"
    turbulence:
      les:
        enabled: false
      rans:
        enabled: false
      wall_function:
        enabled: false
    particles:
      count: 0
      init_mode: "Surface"
      restart_mode: "init"
```

Common mappings:
- `domain.blocks` -> `-nblk`
- periodic axes are derived from paired `PERIODIC` boundary conditions before
  DMDA creation; `models.domain` does not accept periodic flags
- `physics.dimensionality: "2D"` -> `-TwoD 1`
- `physics.turbulence.les.enabled/model` -> `-les` (`0` none, `1` constant Smagorinsky, `2` dynamic Smagorinsky)
- `physics.turbulence.les.constant_cs` -> `-les_constant_cs`
- `physics.turbulence.les.dynamic_frequency` -> `-les_dynamic_frequency`
- `physics.turbulence.les.filter_width` -> `-les_filter_width`
- `physics.turbulence.les.test_filter.kernel/width_ratio` -> `-les_test_filter_kernel`, `-les_test_filter_width_ratio`
- `physics.turbulence.les.averaging.mode/directions` -> `-les_averaging_mode`, `-les_averaging_directions`
- `physics.turbulence.les.clipping.mode/max_cs/min_viscosity_ratio` -> `-les_clip_mode`, `-les_clip_max_cs`, `-les_min_viscosity_ratio`
- `physics.turbulence.les.gradient_model.enabled` -> `-les_gradient_model`
- `physics.turbulence.les.diagnostics.enabled/cadence/yoshizawa_ci` -> `-les_diagnostics`, `-les_diagnostics_cadence`, `-les_yoshizawa_ci`
- `physics.turbulence.rans.enabled/model` -> `-rans` (`k_omega` accepted; runtime update currently incomplete)
- `physics.turbulence.wall_function.enabled/model` -> `-wallfunction` (`1` log law, `2` Werner-Wengle, `3` Cabot)
- `physics.turbulence.wall_function.roughness_height` -> `-wall_roughness` (`log_law` only)
- `physics.particles.count` -> `-numParticles`
- `physics.particles.init_mode` -> `-pinit` (`Surface`, `Volume`, `PointSource`, `SurfaceEdges`)
- `physics.particles.restart_mode` -> `-particle_restart_mode`
- point source coordinates -> `-psrc_x/-psrc_y/-psrc_z`

Legacy turbulence shorthand remains valid:

- `les: false` -> `-les 0`
- `les: true` or `les: 1` -> constant Smagorinsky (`-les 1`)
- `les: 2` -> dynamic Smagorinsky (`-les 2`)

LES and RANS are mutually exclusive in one case. `wall_function` is a sibling option because wall functions can be used with wall-modeled LES or RANS boundary treatments.

Restart note:

- if `run_control.start_step > 0`, particles are enabled, and `restart_mode` is omitted, `picurv` warns that C will default to `load`.

For mode-specific particle behavior and restart flow, see **@subpage 45_Particle_Initialization_and_Restart**.

@section p07_les_sec 5. LES Subgrid Models

The subgrid model is chosen by `models.physics.turbulence.les.model`.

@htmlinclude generated/capability_inventory_turbulence_les_model.html

@subsection p07_cap_les_none_sub none

@anchor p07_cap_les_none

**Identity.** `les.enabled: false`, or `model: none` -> `-les 0` ->
`NO_LES_MODEL`. Accepted spellings: `off`, `disabled`, `no_les`. The legacy
shorthand `les: false` is equivalent.

**What it does.** Disables subgrid modelling entirely. The momentum equations are
solved with molecular viscosity alone.

**When to choose it.** Laminar cases, and any DNS where the grid resolves the
dissipative scales. Also the correct choice while verifying something else: a
subgrid model is one more thing to be wrong when you are chasing a discrepancy.

**Parameters it owns.** None. Setting `constant_cs` or dynamic-model controls
alongside `none` has no effect.

**Interactions.** LES and RANS are mutually exclusive within one case. Wall
functions are configured separately and are not implied by disabling LES.

**Diagnostics.** The startup banner reports the resolved turbulence model. With no
model active, no eddy-viscosity field is computed or logged.

**Evidence.** Production exercised - `examples/flat_channel` and
`examples/bent_channel` both run with no subgrid model.

**Limitations.** None; this is the absence of a model rather than a model.

@subsection p07_cap_les_constant_smagorinsky_sub constant_smagorinsky

@anchor p07_cap_les_constant_smagorinsky

**Identity.** `model: constant_smagorinsky` -> `-les 1` -> `CONSTANT_SMAGORINSKY`.
Accepted spellings: `constant`, `smagorinsky`. The legacy shorthand `les: true` and
`les: 1` both select this model.

**What it does.** Applies a Smagorinsky eddy viscosity built from a **fixed**
coefficient and the local strain-rate magnitude, added to the molecular viscosity:
`nu_t = Cs^2 Delta^2 |S|`.

**When to choose it.** When you want subgrid dissipation, know a coefficient
appropriate to your flow, and would rather not pay for the dynamic procedure. It is
also the right choice on a geometry with no homogeneous direction, where a local
dynamic coefficient is at its noisiest.

**Parameters it owns.** `constant_cs` -> `-les_constant_cs`, the fixed coefficient, and
`filter_width` -> `-les_filter_width`, which sets how `Delta` is derived per cell.

@warning `constant_cs` defaults to **0.03**, which is a wall-bounded value. Isotropic
turbulence conventionally uses Lilly's 0.16-0.17, but only where the grid cutoff sits
in an inertial range; at low `Re_lambda` a smaller value near 0.1 is more appropriate.
The default is not a universal choice - set it deliberately for your flow.

**Interactions.** Mutually exclusive with RANS. The test-filter, averaging, and
clipping controls belong to the dynamic procedure and are rejected here rather than
silently ignored.

**Storage.** This model allocates no coefficient field. The coefficient is a number
from the configuration, so `nu_t` is built from it directly; nothing is synchronized,
checkpointed, or written to disk on its behalf.

**Diagnostics.** The startup banner reports the resolved turbulence model. With
`diagnostics.enabled` set, `<run.runtime_logs>/les_coefficient.csv` records the eddy-viscosity levels
and the modelled subgrid energy each step. The eddy-viscosity field is available for
output.

**Evidence.** Implemented, with unit coverage: `tests/c/test_les.c` case
`constant-model-needs-no-coefficient-field` runs the model with no coefficient field
allocated and checks the resulting `nu_t` against the closed form.

**Limitations.** The coefficient is constant in space and time, so it cannot adapt near
walls and does not vanish in laminar regions the way a dynamic procedure does.

@subsection p07_cap_les_dynamic_smagorinsky_sub dynamic_smagorinsky

@anchor p07_cap_les_dynamic_smagorinsky

**Identity.** `model: dynamic_smagorinsky` -> `-les 2` -> `DYNAMIC_SMAGORINSKY`.
Accepted spelling: `dynamic`. The legacy shorthand `les: 2` selects it.

**What it does.** Measures the model coefficient from the resolved field each update
rather than prescribing it, using the Germano identity and Lilly's least-squares
contraction. The resolved velocity is filtered a second time at a wider test scale; the
stress carried between the two scales is computable, and matching the model against it
determines the coefficient. The full derivation is in
**@subpage 72_LES_Turbulence_Closure**.

**When to choose it.** When the coefficient should respond to the flow: transitional
regions, near walls, or decaying turbulence, where a single constant is wrong somewhere
in the domain.

**Parameters it owns.** `dynamic_frequency` -> `-les_dynamic_frequency`, plus the
`filter_width`, `test_filter`, `averaging`, and `clipping` blocks documented below.

**Interactions.** Mutually exclusive with RANS. `constant_cs` belongs to the constant
model and is rejected here.

**Storage.** The coefficient field is stored in `CS` and checkpointed. It holds `C`, the
factor multiplying `Delta^2 |S|`, which is `Cs^2` in the classical notation - not `Cs`.
Under `clipping.mode: none` it is signed, because a negative coefficient is backscatter.

**Diagnostics.** With `diagnostics.enabled`, `<run.runtime_logs>/les_coefficient.csv` records the
effective coefficient converted to `Cs`, its spatial spread, eddy-viscosity levels, and
the fractions of the domain that were backscattering or limited before clipping. Those
last two describe a pre-clipping state that no stored field preserves.
`diagnostics.cadence` sets the interval, and `diagnostics.yoshizawa_ci` sets the
constant in the reported subgrid kinetic energy; neither affects the solution.
`picurv summarize --plot les.cs_effective` renders the coefficient history directly.

**Evidence.** Implemented, with unit coverage in `tests/c/test_les.c`: the model tensor
is checked against its closed form on constant strain, the filtered product is pinned
apart from the product of filtered factors, the procedure returns exactly zero on
uniform flow, and global averaging is checked to give one coefficient per block. The
averaging reduction is checked to be decomposition-independent in
`tests/c/test_mpi_kernels.c`. `examples/decaying_isotropic_turbulence` exercises the
model end to end; on a shortened 32^3 run of it the coefficient field comes out uniform
to roundoff under `averaging.mode: homogeneous`, as averaging over a homogeneous set
must, and `Cs` climbs monotonically from a random initial field. That is a trend check,
not a validation - the run was not carried far enough for `Cs` to settle.

@note **Status: experimental.** The formulation is correct and unit-tested, but the
coefficient has not yet been validated against a reference flow. The check that would
settle it is decaying isotropic turbulence with homogeneous averaging, where `Cs(t)`
should settle near 0.16-0.17; until that run is recorded, treat the magnitude as
uncharacterized.

**Limitations.** The procedure needs a developed field to sample, so the coefficient is
held at zero for the first two steps of a run started from rest. With
`averaging.mode: local` the coefficient is noisy and the least-squares closure is
formally inconsistent, since it assumes the coefficient is constant over the averaging
set; prefer `homogeneous` wherever the flow has a homogeneous direction.

@section p07_les_width_sec 5.1 Grid Filter Width Entries

`filter_width` sets how each cell's filter width `Delta` is derived from its metrics.
It applies to both LES models, since both scale the eddy viscosity by `Delta^2`.

@htmlinclude generated/capability_inventory_turbulence_les_filter_width.html

@subsection p07_cap_width_cube_root_volume_sub cube_root_volume

@anchor p07_cap_width_cube_root_volume

**Identity.** `les.filter_width: cube_root_volume` -> `-les_filter_width 0`. The default.

**What it does.** Sets `Delta` to the cube root of the cell volume.

**When to choose it.** Near-isotropic cells, where it is exact and cheapest.

**Parameters it owns.** None.

**Interactions.** None beyond scaling `nu_t`.

**Diagnostics.** No dedicated output; its effect appears in the eddy-viscosity level.

**Evidence.** Implemented, covered by `tests/c/test_les.c` case
`filter-width-models-separate-on-stretched-cell`.

**Limitations.** Underestimates the width on a stretched cell, because a geometric mean
is dominated by the short directions.

@subsection p07_cap_width_geometric_mean_sub geometric_mean

@anchor p07_cap_width_geometric_mean

**Identity.** `les.filter_width: geometric_mean` -> `-les_filter_width 1`.

**What it does.** Sets `Delta` to the geometric mean of the Cartesian cell extents
resolved from the cell metrics rather than from the volume alone.

**When to choose it.** Curvilinear grids where the cell volume and the metric-derived
extents disagree, so the width should follow the actual cell shape.

**Parameters it owns.** None.

**Interactions.** None beyond scaling `nu_t`.

**Diagnostics.** As above.

**Evidence.** Implemented, covered by the same unit case.

**Limitations.** Still a geometric mean, so it shares the cube-root model's optimism on
strongly stretched cells.

@subsection p07_cap_width_max_edge_sub max_edge

@anchor p07_cap_width_max_edge

**Identity.** `les.filter_width: max_edge` -> `-les_filter_width 2`.

**What it does.** Sets `Delta` to the longest Cartesian cell extent.

**When to choose it.** Strongly stretched grids, such as a wall-normal channel mesh,
where the largest unresolved scale is set by the long direction.

**Parameters it owns.** None.

**Interactions.** None beyond scaling `nu_t`.

**Diagnostics.** As above.

**Evidence.** Implemented, covered by the same unit case.

**Limitations.** The most dissipative of the three; on a near-isotropic grid it
overestimates the width and adds subgrid dissipation the flow does not need.

@section p07_les_avg_sec 5.2 Coefficient Averaging Entries

`averaging.mode` selects the set the two Germano contractions are averaged over before
they are divided. Lilly's closure assumes the coefficient is constant across that set,
so the set should span directions in which the flow really is statistically
homogeneous. Dynamic model only.

@htmlinclude generated/capability_inventory_turbulence_les_averaging_mode.html

@subsection p07_cap_avg_local_sub local

@anchor p07_cap_avg_local

**Identity.** `les.averaging.mode: local` -> `-les_averaging_mode 0`. The default.

**What it does.** Divides the two contractions cell by cell, giving a coefficient field
that varies pointwise.

**When to choose it.** Geometry with no homogeneous direction, where no larger averaging
set can be justified.

**Parameters it owns.** None.

**Interactions.** `averaging.directions` does not apply and is rejected.

**Diagnostics.** The coefficient spread columns in `<run.runtime_logs>/les_coefficient.csv` show how
noisy the field is.

**Evidence.** Implemented, covered by `tests/c/test_les.c` case
`average-ratio-local-is-pointwise`.

**Limitations.** The denominator collapses wherever the resolved strain is briefly
small, so the coefficient is noisy, and the least-squares derivation is formally
inconsistent because the coefficient it assumes constant is not.

@subsection p07_cap_avg_homogeneous_sub homogeneous

@anchor p07_cap_avg_homogeneous

**Identity.** `les.averaging.mode: homogeneous` -> `-les_averaging_mode 1`.

**What it does.** Averages both contractions over the homogeneous directions, then
divides. The directions come from `averaging.directions` when given, and otherwise from
the case's resolved periodicity, which the loader derives from the boundary pairs and
rejects the case unless opposite faces, and all blocks, agree on. A triply periodic box therefore yields one coefficient for the
whole domain, and a channel periodic in `i` and `k` yields a wall-normal profile, in
both cases with no extra configuration.

**When to choose it.** Whenever the flow has a homogeneous direction. This is the mode
the dynamic procedure was derived for.

**Parameters it owns.** `averaging.directions` -> `-les_averaging_directions`, a subset
of `[i, j, k]` overriding the periodic default.

**Interactions.** Naming a direction the case does not declare `PERIODIC` is accepted
with a warning, since homogeneity is the user's claim to make. Requesting this mode with
neither periodic pairs nor an explicit list is rejected.

**Diagnostics.** `cs_effective` in `<run.runtime_logs>/les_coefficient.csv` is the whole-domain value
regardless of mode, so it stays comparable across modes.

**Evidence.** Implemented, covered by `tests/c/test_les.c` cases
`homogeneous-averaging-derives-periodic-axes` and
`average-ratio-retains-unaveraged-direction`.

**Limitations.** Averaging over a direction the flow is not homogeneous in smears real
spatial variation into a single number.

@subsection p07_cap_avg_global_sub global

@anchor p07_cap_avg_global

**Identity.** `les.averaging.mode: global` -> `-les_averaging_mode 2`.

**What it does.** Averages both contractions over the entire block, giving one
coefficient per block per update.

**When to choose it.** A domain you know is homogeneous but have not declared periodic.
Where the domain is periodic in all three directions, `homogeneous` reaches the same
answer without asserting anything the boundary conditions do not already say.

**Parameters it owns.** None.

**Interactions.** `averaging.directions` does not apply and is rejected. In a
multi-block case the average is per block.

**Diagnostics.** As above.

**Evidence.** Implemented, covered by `tests/c/test_les.c` cases
`average-ratio-divides-summed-fields` and
`dynamic-procedure-global-average-is-uniform`, and by the decomposition-independence
case in `tests/c/test_mpi_kernels.c`.

**Limitations.** Discards all spatial variation in the coefficient, which is wrong
wherever the flow is inhomogeneous.

@section p07_les_clip_sec 5.3 Coefficient Limiting Entries

`clipping.mode` bounds the coefficient the contraction produces. Dynamic model only.

@htmlinclude generated/capability_inventory_turbulence_les_clip_mode.html

@subsection p07_cap_clip_clamp_sub clamp

@anchor p07_cap_clip_clamp

**Identity.** `les.clipping.mode: clamp` -> `-les_clip_mode 0`. The default.

**What it does.** Restricts the coefficient to `[0, max_cs^2]`.

**When to choose it.** The conservative default, and the safe choice when a run is at
risk of a diverging coefficient.

**Parameters it owns.** `clipping.max_cs` -> `-les_clip_max_cs`, a ceiling on `Cs` and
not on the stored coefficient. Defaults to 0.3, roughly twice the physical value, so it
catches divergence without shaping the ordinary distribution.

**Interactions.** Removes backscatter, so the total-viscosity floor never engages from
below.

**Diagnostics.** `limited_fraction` in `<run.runtime_logs>/les_coefficient.csv` reports the volume
fraction the clip modified. A ceiling that is doing nothing reads near zero.

**Evidence.** Implemented, covered by `tests/c/test_les.c` case
`clip-model-coefficient-modes`.

**Limitations.** Discarding the negative tail is not neutral: the positives that would
have cancelled it survive, so the mean subgrid dissipation is biased upward.

@subsection p07_cap_clip_clip_negative_sub clip_negative

@anchor p07_cap_clip_clip_negative

**Identity.** `les.clipping.mode: clip_negative` -> `-les_clip_mode 1`.

**What it does.** Discards negative coefficients but imposes no ceiling.

**When to choose it.** When backscatter is not wanted but a ceiling would distort a
coefficient you have reason to trust at large values.

**Parameters it owns.** None; `clipping.max_cs` imposes no bound here, and naming it under this mode is rejected rather than ignored.

**Interactions.** As for `clamp`.

**Diagnostics.** As for `clamp`.

**Evidence.** Implemented, covered by the same unit case.

**Limitations.** Carries the same upward bias in mean dissipation as `clamp`, without
the divergence guard.

@subsection p07_cap_clip_none_sub none

@anchor p07_cap_clip_none

**Identity.** `les.clipping.mode: none` -> `-les_clip_mode 2`.

**What it does.** Keeps the signed coefficient, so cells where the resolved scales are
receiving energy from the unresolved ones - backscatter - are represented rather than
zeroed.

**When to choose it.** When the physics of backscatter matters, or when the upward bias
the clipping modes introduce in mean dissipation is itself the thing under study.

**Parameters it owns.** None directly, but it is the mode that makes
`clipping.min_viscosity_ratio` -> `-les_min_viscosity_ratio` load-bearing: the total
viscosity is floored at `ratio * nu`, which bounds the quantity that actually has to
stay positive for the momentum operator to remain well posed.

**Interactions.** The stored coefficient field becomes signed. Anything reading `CS`
must expect negative values.

**Diagnostics.** `backscatter_fraction` in `<run.runtime_logs>/les_coefficient.csv` reports the volume
fraction with a negative coefficient, which is exactly what the other two modes discard.

**Evidence.** Implemented, covered by `tests/c/test_les.c` cases
`clip-model-coefficient-modes` and `eddy-viscosity-floor-bounds-total-viscosity`.

**Limitations.** A locally averaged coefficient is negative at a large fraction of
points in developed turbulence, so this mode leans hard on the viscosity floor. Pair it
with `averaging.mode: homogeneous`, where the averaged coefficient is negative only when
the whole homogeneous set is backscattering on balance.

@section p07_rans_filter_sec 6. RANS and Test-Filter Entries

@htmlinclude generated/capability_inventory_turbulence_rans_model.html
@htmlinclude generated/capability_inventory_turbulence_les_test_filter.html

@subsection p07_cap_rans_none_sub none

@anchor p07_cap_rans_none

**Identity.** `turbulence.rans.enabled: false`, or `model: none` -> `-rans 0`. Accepted
spellings: `off`, `disabled`.

**What it does.** Disables RANS modelling. Momentum is closed by molecular viscosity, or
by LES if that is enabled instead.

**When to choose it.** Whenever you are not running RANS - which, given the status of the
alternative below, is currently always.

**Parameters it owns.** None.

**Interactions.** RANS and LES are mutually exclusive within one case.

**Diagnostics.** The startup banner reports the resolved turbulence model.

**Evidence.** Production exercised - `examples/flat_channel` and `examples/bent_channel`
both run with RANS off.

**Limitations.** None; this is the absence of a model.

@subsection p07_cap_rans_k_omega_sub k_omega

@anchor p07_cap_rans_k_omega

**Identity.** `turbulence.rans.model: k_omega` -> `-rans 1`. Accepted spelling: `komega`.

**What it does.** Intended to close the momentum equations with a two-equation k-omega
model.

**When to choose it.** Not currently - see the status below.

**Parameters it owns.** The RANS block in `case.yml`.

**Interactions.** Mutually exclusive with LES. Wall functions are configured separately
and are not implied by enabling RANS.

**Diagnostics.** The startup banner reports the model as resolved even though the update
is incomplete, so the banner alone is not evidence that it is working.

**Evidence.** Implemented only. No facet is claimed.

@warning **Status: experimental - the runtime update is incomplete.** The configuration
layer accepts `k_omega` and the flag reaches the runtime, but the transport equations are
not fully updated each step. Do not treat RANS results from this tree as meaningful.

**Limitations.** Beyond the incomplete update, only `k_omega` is exposed; no other
closure is selectable.

@subsection p07_cap_filter_volume_weighted_box_sub volume_weighted_box

@anchor p07_cap_filter_volume_weighted_box

**Identity.** `les.test_filter.kernel: volume_weighted_box` -> `-les_test_filter_kernel 0`.
Accepted spelling: `box`. The default.

**What it does.** Averages the 3x3x3 stencil with cell-volume weights, so neighbouring
cells of different size contribute in proportion to the volume they represent. Makes no
assumption about which directions are homogeneous, and excludes solid cells.

**When to choose it.** The general-purpose choice, and the only defensible one on a
geometry with no homogeneous direction.

**Parameters it owns.** None. `test_filter.width_ratio` belongs to the block, not to the
kernel, and applies to either kernel.

**Interactions.** Consulted only by @ref p07_cap_les_dynamic_smagorinsky "dynamic_smagorinsky".

**Diagnostics.** No dedicated output; its effect reaches
`<run.runtime_logs>/les_coefficient.csv` through the coefficient.

**Evidence.** Implemented, covered by `tests/c/test_solver_kernels.c`, which checks that
it preserves a constant field and returns zero when the stencil is entirely solid.

**Limitations.** Less accurate than a Simpson stencil where one is admissible.

@subsection p07_cap_filter_simpson_ik_sub simpson_ik

@anchor p07_cap_filter_simpson_ik

**Identity.** `les.test_filter.kernel: simpson_ik` -> `-les_test_filter_kernel 1`.

**What it does.** Applies a two-dimensional Simpson stencil across the central eta-plane
of the 3x3x3 neighbourhood, weighting the nine samples 1, 4, and 16. Volume weights are
not used: the stencil assumes uniform spacing in the plane it averages over.

**When to choose it.** A plane channel or any flow homogeneous in xi and zeta, where the
higher-order stencil is both admissible and more accurate.

**Parameters it owns.** None.

**Interactions.** The homogeneity assumption is now checked rather than trusted:
selecting this kernel without declaring both xi and zeta `PERIODIC` is rejected during
validation.

**Diagnostics.** As above.

**Evidence.** Implemented, covered by `tests/c/test_solver_kernels.c` alongside the box
filter.

**Limitations.** Averages over only the central eta-plane, so it ignores variation in
the wall-normal direction entirely; that is the point on a channel and wrong anywhere
else.

@htmlinclude generated/capability_inventory_turbulence_wall_function.html

A wall model supplies the stress of a boundary layer the mesh does not resolve, so it is
only meaningful when the unresolved motions are modelled somewhere. Three pairings are
therefore rejected before a run starts.

- **No turbulence model.** A wall model with LES and RANS both off is refused. This
  solver has no implicit-LES scheme to stand in for the missing closure: its convection
  is QUICK, whose numerical dissipation is linear and upwind-biased rather than a
  limiter-based truncation error standing in for a subgrid stress. Enable LES, or
  resolve the wall and switch the wall function off.
- **`cabot` with RANS.** Cabot closes the wall layer with its own mixing-length eddy
  viscosity. Under a RANS model that layer would carry two closures with no matching
  between them.
- **`werner` with RANS.** Werner-Wengle applies its power law to an instantaneous
  filtered velocity, which is a large-eddy quantity. A RANS field is already averaged and
  wants a law derived for the mean profile.

A wall model on a laminar case is refused for the same reason from the other direction:
these laws describe a turbulent boundary layer, and below transition there is no inertial
region for them to represent.

Where the model's stress actually enters is worth stating, because it is not obvious. The
correction sets the near-wall velocity, and the momentum equation reaches the wall through
a viscosity times a velocity gradient - so the modelled stress arrives only if that
viscosity is the one which reproduces it, `nu_eff = tau_w y / u`. Molecular viscosity
alone delivers a fraction `u+/y+` of the stress, which at `y+ = 267` is about a
fourteenth of it. The wall model therefore installs its own effective eddy viscosity at
its wall face, in place of the subgrid value, which is zero there in a wall-resolved run.
`nu_wall_over_nu_mean` in `<run.runtime_logs>/wall_model.csv` reports it.

@subsection p07_cap_wall_log_law_sub log_law

@anchor p07_cap_wall_log_law

**Identity.** `turbulence.wall_function.model: log_law` -> `-wallfunction`. Accepted
spelling: `loglaw`. It is also the default when `wall_function.enabled` is true and no
model is named.

**What it does.** Applies a logarithmic law-of-the-wall boundary treatment instead of
resolving the viscous sublayer.

**When to choose it.** When the near-wall grid is too coarse to resolve the sublayer and
you accept a modelled wall stress rather than a computed one.

**Parameters it owns.** `wall_function.roughness_height` -> `-wall_roughness`.

**Interactions.** Configured independently of LES and RANS - enabling a turbulence model
does not imply a wall function, and vice versa. The wall handler itself remains
@ref p44_cap_noslip "noslip"; the wall function changes how the stress is imposed, not the
BC selection.

The correction also reaches the LES closure. `Contra2Cart` rebuilds the Cartesian field
from the contravariant one at the top of each timestep and so discards the previous
solve's near-wall correction; the wall model is re-applied before the subgrid strain
rates are formed, so the closure and the momentum equation agree about the near-wall
velocity rather than disagreeing within one timestep.

**Diagnostics.** The startup banner reports whether wall functions are enabled and the
configured roughness height.

**Evidence.** Implemented, with unit coverage of the velocity laws and their
friction-velocity root-finds in `tests/c/test_runtime_kernels.c` - each is checked to
invert its own law - and of the integration in `tests/c/test_boundaries.c`, which checks
that the corrected first interior velocity is the modelled profile at that cell's wall
distance, for the friction velocity the integration itself stored. No shipped example
enables it and no reference-flow comparison has been run.

**Limitations.** Validity depends on the first cell falling in the logarithmic region,
which nothing checks for you.

@subsection p07_cap_wall_werner_sub werner

@anchor p07_cap_wall_werner

**Identity.** `turbulence.wall_function.model: werner` -> `-wallfunction 2`. Accepted
spelling: `werner_wengle`.

**What it does.** Applies the Werner-Wengle power-law wall model. The two-layer profile
- linear below `y+ = 11.81`, `u+ = 8.3 (y+)^(1/7)` above - is assumed to hold across the
first cell and integrated, which inverts in closed form: the wall stress comes straight
from the resolved speed with no root-find at all, in either region. The corrected
velocity is produced by the exact inverse of that same relation, so a cell nearer the
wall is always assigned the slower speed.

**When to choose it.** When the first cell may fall inside or near the viscous sublayer,
where a pure log law has no valid branch. Its explicit form also avoids the inner
iteration the log law needs, which matters for a matrix-free Jacobian: an iterative inner
solve makes the residual only piecewise smooth.

**Parameters it owns.** None. It has no roughness formulation, so `roughness_height` is
rejected rather than accepted and ignored - a rough wall is not representable through
this model.

**Interactions.** As for @ref p07_cap_wall_log_law "log_law".

**Diagnostics.** As for @ref p07_cap_wall_log_law "log_law"; the banner reports the
selected model by name.

**Evidence.** Implemented, with unit coverage of `find_utau_Werner` and `u_Werner` in
`tests/c/test_runtime_kernels.c`, each checked to invert the other. The integration that
dispatches to it is covered in `tests/c/test_boundaries.c`. No reference-flow comparison
has been run.

**Limitations.** No roughness. The power-law constants are fixed, and the model carries
no pressure-gradient term, so it describes an equilibrium layer only - `cabot` is the
selection when the near-wall layer is not in equilibrium.


@subsection p07_cap_wall_cabot_sub cabot

@anchor p07_cap_wall_cabot

**Identity.** `turbulence.wall_function.model: cabot` -> `-wallfunction 3`.

**What it does.** Applies Cabot's wall model, which integrates a mixing-length eddy
viscosity across the wall layer rather than assuming a profile shape.

**When to choose it.** When the near-wall flow departs from equilibrium enough that a
fixed profile is a poor fit, and you are prepared to pay for the integration.

**Parameters it owns.** None currently.

**What it needs.** Cabot's departure from an equilibrium profile is driven by the
wall-parallel pressure gradient, which is taken from the resolved pressure at the
near-wall cell through @ref ComputeScalarFieldDerivatives. The pressure is the previous
projection's, the same lag every other explicit use of it carries.

**Interactions.** As for @ref p07_cap_wall_log_law "log_law".

**Diagnostics.** As for @ref p07_cap_wall_log_law "log_law".

**Evidence.** Implemented, with unit coverage of `find_utau_Cabot` and `u_Cabot` in
`tests/c/test_runtime_kernels.c`. The integration that dispatches to it is covered in
`tests/c/test_boundaries.c`. No reference-flow comparison has been run.

**Limitations.** No roughness formulation, so `roughness_height` is rejected rather than
accepted and ignored. The mixing-length constant is fixed.

@section p07_bc_sec 7. boundary_conditions

Single-block syntax: list of 6 face entries.
Multi-block syntax: list-of-lists, one 6-face list per block.

Supported face names:
- `-Xi`, `+Xi`, `-Eta`, `+Eta`, `-Zeta`, `+Zeta`

Supported type/handler combinations:
- `INLET` + `constant_velocity` (`vx/vy/vz`)
- `INLET` + `parabolic` (`v_max`)
- `INLET` + `prescribed_flow`:
  - file-backed: `params.source.type: file`, `params.source.path`
  - generated: `params.source.type: generated`, `params.source.generator: square_duct_poiseuille`
  - field-sliced: `params.source.type: field_slice`, `params.source.field_file`, `params.source.grid_file`, `params.source.slice`
  - generated and field-sliced sources default to `generators/profile.gen`; optional
    `params.source.script` selects a compatible override
- `OUTLET` + `conservation`
- `WALL` + `noslip`
- `PERIODIC` + `geometric`
- `PERIODIC` + `constant_flux` (`target_flux`, optional `enforce_seam_flux`; `apply_trim` is a deprecated alias)
- `PERIODIC` + `initial_flux` (no params; optional `enforce_seam_flux`) - target measured from the initial condition

All six faces must be explicitly provided for each block.
For detailed handler semantics, validation constraints, and C dispatch path, see **@subpage 44_Boundary_Conditions_Guide**.

Generated profile example:

```yaml
- face: "-Zeta"
  type: INLET
  handler: prescribed_flow
  params:
    source:
      type: generated
      generator: square_duct_poiseuille
      params:
        bulk_velocity: 1.0
        n_terms: 101
```

`picurv run --solve` generates the dimensional `.picslice`, writes
`profile.info`, stages the solver-scale `.picslice`, and passes the existing
`source_file` key to the C runtime. Use `picurv precompute --case ...` to create
the same deterministic artifacts without launching the solver.

Field-sliced profile example:

```yaml
- face: "-Zeta"
  type: INLET
  handler: prescribed_flow
  params:
    source:
      type: field_slice
      field_file: ../old_run/output/checkpoints/step_000000010000/eulerian/block_0000/Ucat.dat
      grid_file: ../old_run/config/grid.run
      source_case: ../old_run/config/case.yml
      slice:
        face: "+Zeta"
        orientation: opposite
```

`field_slice` uses Python preprocessing to write a normal dimensional PICSLICE;
the C runtime still sees only the staged `source_file`.

@section p07_passthrough_sec 8. solver_parameters (Advanced)

Optional escape hatch for flags not yet exposed in structured schema:

```yaml
solver_parameters:
  -read_fields: true
  -some_new_flag: 123
```

Use sparingly and prefer structured keys when available.

@section p07_modular_sec 9. Mixing With Other Profiles

`case.yml` is designed to be combined with reusable profiles for the other config roles.

Common patterns:

- one `case.yml` + multiple `monitor.yml` files (debug vs production output),
- one `case.yml` + multiple `post.yml` recipes (quick scalar check vs heavy VTK/statistics),
- one `solver.yml` reused across many `case.yml` files,
- one `cluster.yml` reused across many runs and sweeps.

For worked combinations, see **@subpage 49_Workflow_Recipes_and_Config_Cookbook**.

@section p07_next_steps_sec 8. Next Steps

Proceed to **@subpage 08_Solver_Reference**.

Cross-file contract/mapping:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 48_Grid_Generator_Guide**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 32_Analytical_Solutions**
- **@subpage 33_Initial_Conditions**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**
