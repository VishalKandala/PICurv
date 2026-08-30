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
- `physics.turbulence.les.constant_cs` -> `-const_cs`
- `physics.turbulence.les.max_cs` -> `-max_cs`
- `physics.turbulence.les.dynamic_frequency` -> `-dynamic_freq`
- `physics.turbulence.les.test_filter` -> `-testfilter_ik` (`volume_weighted_box` = 0, `homogeneous_ik` = 1)
- `physics.turbulence.rans.enabled/model` -> `-rans` (`k_omega` accepted; runtime update currently incomplete)
- `physics.turbulence.wall_function.enabled` -> `-wallfunction`
- `physics.turbulence.wall_function.roughness_height` -> `-wall_roughness`
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
coefficient and the local strain-rate magnitude, added to the molecular viscosity.

**When to choose it.** **Not currently** - see the status below. On a fresh run the
coefficient is never computed, so this model applies no subgrid dissipation at all.
It is the model the shipped LES cases are configured for, which is why those cases
inherit the defect rather than escaping it. If the defect is fixed, this becomes the
model to choose whenever you want subgrid dissipation and are prepared to pick a
coefficient appropriate to your flow.

**Parameters it owns.** `constant_cs` -> `-const_cs`, the fixed Smagorinsky
coefficient.

@warning `constant_cs` defaults to **0.03**, which is a wall-bounded value. Isotropic
turbulence conventionally uses Lilly's 0.16-0.17, but only where the grid cutoff sits
in an inertial range; at low `Re_lambda` a smaller value near 0.1 is more appropriate.
The default is not a universal choice - set it deliberately for your flow.

**Interactions.** Mutually exclusive with RANS. The test-filter and dynamic-frequency
controls belong to the dynamic model, not this one, and have no effect here.

**Diagnostics.** The startup banner reports the resolved turbulence **model**; it does
not print the coefficient. `ComputeSmagorinskyConstant` logs the coefficient at
`LOG_INFO` when it actually runs - and its absence from the log is exactly the
symptom of the defect below. The eddy-viscosity field is available for output.

**Evidence.** Implemented only. `examples/periodic_test/driven_channel/les_retau180`
selects this model, but that case inherits the defect below, so its selection is not
evidence that the model works. `decaying_isotropic_turbulence` selects the dynamic
model instead.

@warning **Status: known-defective on fresh runs. Do not use for published work.**
On a fresh run (`StartStep == 0`), `FlowSolver` calls `ComputeSmagorinskyConstant`
only at `step == StartStep + 1`, i.e. step 1. That call hits the early-return branch
that sets `CS = 0` for `step < 2`, and for the constant model the function is never
called again. The result is a run that reports LES as enabled while the Smagorinsky
coefficient stays zero - no subgrid dissipation is applied. The existing unit test
invokes the helper directly at step 2 and therefore does not exercise this lifecycle.
This is a runtime defect tracked outside the documentation work; it has not been
fixed here.

**Limitations.** Beyond the defect above: the coefficient is constant in space and
time, so it cannot adapt near walls and does not vanish in laminar regions the way a
dynamic procedure would.

@subsection p07_cap_les_dynamic_smagorinsky_sub dynamic_smagorinsky

@anchor p07_cap_les_dynamic_smagorinsky

**Identity.** `model: dynamic_smagorinsky` -> `-les 2` -> `DYNAMIC_SMAGORINSKY`.
Accepted spelling: `dynamic`. The legacy shorthand `les: 2` selects it.

**What it does.** Intended to compute the Smagorinsky coefficient dynamically from a
test-filtered field via the Germano identity, rather than from a fixed constant.

**When to choose it.** Not currently. See the status below.

**Parameters it owns.** `dynamic_frequency` -> `-dynamic_freq`, `test_filter` ->
`-testfilter_ik`, and `max_cs` -> `-max_cs` (the coefficient clip). The
`-i_homo_filter`/`-j_homo_filter`/`-k_homo_filter` options are parsed but inert while
the averaging they gate remains unimplemented.

**Interactions.** Mutually exclusive with RANS. `constant_cs` belongs to the constant
model and has no effect here.

**Diagnostics.** `ComputeSmagorinskyConstant` logs the maximum computed coefficient
and the clip value when it runs.

**Evidence.** Implemented only. No facet is claimed.

@warning **Status: known-defective. Do not use for scientific work.** The model tensor
is built from separately filtered quantities rather than the test filter of the
product, so both terms share an identical factor and `M_ij` collapses to a fixed
scalar multiple of the filtered strain-rate tensor. The Germano scale-similarity
relation is therefore never evaluated - the procedure has no dynamic content. The
homogeneous averaging it assumes is also an unimplemented placeholder.

**Limitations.** Beyond the defect: the detailed formulation, averaging strategy, and
validation evidence for this model are **deliberately deferred** until the
implementation is corrected. The scoped record of what that documentation must cover
is in `tests/tooling/capability_scope_records.json`.

@warning **Neither LES model is currently safe for scientific use.**
`dynamic_smagorinsky` is selectable but known-defective and deliberately undocumented
(record: `tests/tooling/capability_scope_records.json`).
`constant_smagorinsky` is inert on fresh runs, as described above. The shipped
`decaying_isotropic_turbulence` example selects the dynamic model and inherits its
defect. Treat LES results from this tree as uncharacterized until both issues are
fixed and tested end to end.

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

**Identity.** `les.test_filter: volume_weighted_box` -> `-testfilter_ik 0`. Accepted
spellings: `box`, `general_box`.

**What it does.** Applies a volume-weighted box test filter over the local stencil,
making no assumption about which directions are homogeneous.

**When to choose it.** The general-purpose choice, and the only defensible one on a
geometry with no homogeneous direction.

**Parameters it owns.** None; it selects the filter shape used by the dynamic procedure.

**Interactions.** Only consulted by @ref p07_cap_les_dynamic_smagorinsky "dynamic_smagorinsky",
which is known-defective. With the constant model or with LES off, this setting has no
effect.

**Diagnostics.** No dedicated output; its effect is visible only through the dynamic
coefficient, which the defect renders meaningless.

**Evidence.** Implemented only.

**Limitations.** Inherits the dynamic model's defect: selecting a filter cannot rescue a
procedure that never evaluates the Germano identity.

@subsection p07_cap_filter_homogeneous_ik_sub homogeneous_ik

@anchor p07_cap_filter_homogeneous_ik

**Identity.** `les.test_filter: homogeneous_ik` -> `-testfilter_ik 1`. Accepted
spellings: `ik_homogeneous`, `simpson_ik`.

**What it does.** Applies Simpson-weighted filtering in the i and k directions, assuming
those are statistically homogeneous.

**When to choose it.** A plane channel, where the streamwise and spanwise directions are
homogeneous and filtering along them is physically justified.

**Parameters it owns.** None.

**Interactions.** The homogeneity assumption is the user's to make; nothing validates that
i and k actually are homogeneous for your geometry. Same dynamic-model dependency as
above.

**Diagnostics.** As above.

**Evidence.** Implemented only.

**Limitations.** Wrong on any geometry where i and k are not homogeneous, and inherits the
dynamic model's defect regardless.

@htmlinclude generated/capability_inventory_turbulence_wall_function.html

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

**Diagnostics.** The startup banner reports whether wall functions are enabled and the
configured roughness.

**Evidence.** Implemented only. No shipped example enables it.

**Limitations.** The only exposed model. Its validity depends on the first cell falling in
the logarithmic region, which nothing checks for you.

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
