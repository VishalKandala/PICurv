@page 33_Initial_Conditions Initial Condition Modes

@anchor _Initial_Conditions

This page documents how PICurv initializes Eulerian velocity fields and particles at startup.
It covers both user-facing YAML inputs and the C implementation path that applies them.

@tableofcontents

@section p33_overview_sec 1. Where Initialization Happens

Startup sequence:

1. `picurv_cli/core.py` resolves built-in, file, or external-generator YAML.
2. File-backed ICs are staged in the filename layout expected by @ref ReadFieldData.
3. @ref InitializeEulerianState chooses fresh solve, restart load, or analytical initialization.
4. On a fresh solve, @ref PopulateInitialUcont populates `Ucont`; existing finalization then applies
   boundary values and derives the remaining velocity state.
5. If particles are enabled, @ref InitializeParticleSwarm runs independently after Eulerian setup.

@section p33_euler_sec 2. Eulerian Field Initialization (`properties.initial_conditions`)

The canonical contract has two modes: `generated` and `file`. Generated ICs may use a built-in
C generator or the repository `generators/ic.gen` utility. Both produce the same solver-facing result: one
initial velocity field.

Built-in zero:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: zero
```

Built-in Cartesian constant:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: constant
    params:
      u_physical: 1.5
      v_physical: 0.0
      w_physical: 0.0
```

Built-in streamwise constant:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: streamwise_constant
    params:
      velocity_physical: 1.5
      flow_direction: "+Zeta"
```

Built-in Poiseuille:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: poiseuille
    params:
      peak_velocity_physical: 1.5
      flow_direction: "+Zeta"
```

File-backed `Ucat` or `Ucont`:

```yaml
properties:
  initial_conditions:
    mode: file
    field: Ucat
    source_file: initial_conditions/velocity.dat
```

The input must be one PETSc binary `.dat` vector readable by @ref ReadFieldData. `Ucat` inputs
are converted by @ref Cart2Contra; `Ucont` inputs are used directly. The first implementation
supports one block only because it intentionally reuses the existing single-field
@ref ReadFieldData path.

Repository generator:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: ic_gen
    params:
      field: Ucat
      script: tools/custom_ic.py  # optional; defaults to generators/ic.gen
      config_file: config/initial_conditions/expression.cfg
```

`picurv` chooses where the payload goes: it is always written to
`<run.inputs>/initial_condition/initial_condition.generated.dat`. The retired
destination keys `output_file`, `summary_json`, and `spectrum_csv` are refused
with an error rather than ignored.

The launcher invokes `generators/ic.gen` by default, or the optional case-relative/absolute
`params.script` override, as:

```text
python <ic-generator> -c <config_file> --field Ucat|Ucont --output <run.inputs>/initial_condition/initial_condition.generated.dat --grid <run.inputs>/grid/grid.run [cli_args...]
```

`picurv run --solve` materializes the result after grid preparation. `picurv precompute --case ...`
materializes and stages the same artifact without running the solver.

The repository `ic.gen` accepts an INI `[expression]` section. `Ucat` configs
define `u`, `v`, and `w`, evaluated at actual cell centers with extrapolated
dummy-cell centers. `Ucont` configs define `u_xi`, `u_eta`, and `u_zeta`,
evaluated at their corresponding geometric face centers. Expressions may use
`x/y/z`, normalized logical `xi/eta/zeta`, storage `i/j/k`, `pi`, and the
documented numerical functions. The first implementation supports one block.
The repository generator requires a staged PICGRID. `grid.mode: file` and
`grid.mode: grid_gen` provide that grid directly; for single-block
`grid.mode: programmatic_c`, the launcher materializes a nondimensional
`<run.inputs>/grid/grid.run` bridge from scalar `programmatic_settings` before invoking
`ic.gen`.

Repository spectral provider, for a triply periodic box:

```yaml
properties:
  initial_conditions:
    mode: generated
    generator: spectral_random_velocity
    params:
      field: Ucat
      seed: 12345
      spectrum: {type: k4_exponential, k0: 4.0, k_cut: 20.0}
      projection: {type: solenoidal, operator: picurv_discrete}
      normalization: {type: component_rms, target: 1.0}
      remove_mean: true
```

It draws a seeded random field, shapes it to the spectrum envelope
(`k4_exponential` is `E(k) ~ k^4 exp(-2 (k/k0)^2)`, nothing above `k_cut`), projects it
to be divergence-free under the chosen operator, and scales it to the requested
component RMS. It requires every face geometric-periodic and a fresh 3D run, and writes
`<run.analysis>/metrics/initial_condition_summary.json` and the staged spectrum beside
the field. `examples/decaying_isotropic_turbulence` configures it; its README covers the
keys. With `operator: picurv_discrete` the solver's own step-0 divergence is at round-off;
`initial-conditions-2026-09-21` measured 3.1e-14 against a flux scale of 4.3, component RMS 1.000000, the shell
spectrum within 4.7% (energy-weighted) of the envelope, and the same seed reproducing
the same field.

@section p33_entries_sec 3. Field Initialization Mode Entries

@htmlinclude generated/capability_inventory_initial_field_mode.html

Every entry below, and the `streamwise_constant`, `ic_gen` and `file` sources above, was
checked against the field it claims to produce by `initial-conditions-2026-09-21`: on a
programmatic duct the staged interior face fluxes equal their definitions to 3.5e-18 or
exactly, and a `file` source reproduces the generated run it was taken from bitwise.
Step-0 `Ucat` next to walls and inlets reflects the boundary faces, because it is
reconstructed from the fluxes after boundary conditions apply.

@subsection p33_cap_zero_sub Zero

@anchor p33_cap_zero

**Identity.** Field mode `Zero` -> `-finit 0` -> `IC_MODE_ZERO`.

**What it does.** Initializes the velocity field to zero everywhere.

**When to choose it.** When the flow is driven entirely by its boundaries and you want no
assumption about the interior - an inlet-driven channel starting from rest. Choose
`Poiseuille` instead when you want a developed profile immediately, and `Constant` when
you want uniform motion.

**Parameters it owns.** None.

**Interactions.** A zero field with periodic boundaries and no driving force stays zero;
that combination needs a driven handler to do anything.

**Diagnostics.** The startup banner reports the resolved initial-condition mode; step 0
output shows a zero field.

**Evidence.** Analytically verified - `initial-conditions-2026-09-21`: every staged face flux is exactly zero.

**Limitations.** Transition from rest can be slow, and for turbulence it will not occur
at all without a finite-amplitude perturbation.

@subsection p33_cap_constant_sub Constant

@anchor p33_cap_constant

**Identity.** Field mode `Constant` -> `-finit 1` -> `IC_MODE_CONSTANT_CARTESIAN`.

**What it does.** Initializes the whole field to one uniform Cartesian velocity.

**When to choose it.** Plug flow, or seeding a bulk motion the boundaries will then shape.

**Parameters it owns.** The constant velocity components in the initial-condition block.

**Interactions.** A uniform field will not satisfy no-slip at walls; the first steps
resolve that discontinuity, which can be abrupt.

**Diagnostics.** As above.

**Evidence.** Regression verified - `make smoke-periodic` asserts the streamwise-constant
initial-condition banner. Production exercised - `examples/flat_channel` starts from a
constant field. Analytically verified - `initial-conditions-2026-09-21`: `(0.3, -0.2, 0.7)` and a streamwise constant
of 1.3 along `+Zeta` stage the face fluxes they define to round-off.

**Limitations.** Not divergence-free in general on a curvilinear grid, so the first
pressure solve does real work.

@subsection p33_cap_poiseuille_sub Poiseuille

@anchor p33_cap_poiseuille

**Identity.** Field mode `Poiseuille` -> `-finit 2` -> `IC_MODE_POISEUILLE`.

**What it does.** Sets the streamwise face flux to `peak * f1 * f2` times the face area,
where each cross-stream factor is `1 - s^2` and `s` runs from -1 to 1 across the axis in
logical (index) space: cell `c` of an axis with `n` nodes sits at `s = (c - 1/2 - (n-1)/2) /
((n-1)/2)`, so the profile vanishes on the walls. A periodic cross-stream axis has no walls,
and its factor is 1.

**When to choose it.** A wall-bounded channel or duct where you want a developed-looking
start rather than a transient from rest.

**Parameters it owns.** The profile scaling in the initial-condition block.

**Interactions.** The coordinate is logical, so the profile is a parabola in physical space
only on a grid uniform across the section. On a channel periodic in the spanwise axis it is
the exact laminar profile; on a rectangular duct the product of two parabolas is a start,
not the duct solution.

**Diagnostics.** As above; inspect the step-0 field before committing to a long run.

**Evidence.** Unit verified -
`set-initial-interior-field-poiseuille-profile` and
`set-initial-interior-field-poiseuille-uniform-along-periodic-axis` in `make unit-runtime`
check the wall-vanishing product profile and the periodic-axis rule;
`initial-conditions-2026-09-21` found the staged interior face fluxes equal to this formula
to 3.5e-18 in a solver run.

**Limitations.** Logical-space only: on a stretched or curved cross-section the profile is
not a physical parabola. Before 2026-09-18 the parabola vanished at the first cell centres
rather than on the walls, and it was applied across periodic axes too.

@htmlinclude generated/capability_inventory_initial_target_field.html

@subsection p33_cap_field_ucat_sub Ucat

@anchor p33_cap_field_ucat

**Identity.** `initial_conditions.field: Ucat` -> the Cartesian velocity field is the
initialization target (internal code `0`).

**What it does.** Applies the chosen field mode to the Cartesian velocity, from which the
contravariant fluxes are then derived.

**When to choose it.** The normal choice, and the one that matches how the field modes are
expressed - a uniform or parabolic profile is naturally stated in Cartesian components.

**Parameters it owns.** None; it selects which field the mode writes to.

**Interactions.** `Ucat` is the derived representation, so initializing it requires a
conversion to the evolved `Ucont` fluxes before the first step. See
**@subpage 20_Grid_Cell_Architecture_Guide**.

**Diagnostics.** Step-0 output shows the initialized field.

**Evidence.** Production exercised - `examples/flat_channel`. Analytically verified -
`initial-conditions-2026-09-21`: `ic_gen` `Ucat` expressions and every built-in mode stage the fluxes their Cartesian
values define.

**Limitations.** On a strongly curvilinear grid the conversion from a Cartesian profile to
fluxes is not exact in the sense of preserving the intended profile shape.

@subsection p33_cap_field_ucont_sub Ucont

@anchor p33_cap_field_ucont

**Identity.** `initial_conditions.field: Ucont` -> the contravariant flux field is the
initialization target (internal code `1`).

**What it does.** Applies the chosen field mode directly to the evolved flux variable,
skipping the Cartesian-to-flux conversion.

**When to choose it.** When you want exact control of the initial fluxes - for instance
setting a precise volumetric flux - rather than of the Cartesian velocity.

**Parameters it owns.** None.

**Interactions.** `Ucont` is what the solver actually evolves, so this path avoids one
conversion. The field modes are still expressed in the same components, which is what
makes this the less intuitive of the two.

**Diagnostics.** Step-0 output; check the derived `Ucat` looks as intended.

**Evidence.** Analytically verified - `initial-conditions-2026-09-21`: `ic_gen` `Ucont` expressions stage exactly,
and a `file` source of `Ucont` reproduces the generated run bitwise. No shipped example
selects it.

**Limitations.** Easier to get wrong than `Ucat`, because the values are fluxes rather
than velocities: a face flux carries the face area, so a velocity typed where a flux is
expected is off by that area. Only uniform Cartesian grids were checked.

@section p33_euler_modes_sec 4. C Runtime Modes and Entry Points

The launcher maps the contract to one `InitialConditionMode` enum value:

| Initial-condition selection | C mode |
|-----------------------------|--------|
| `generator: zero` | `IC_MODE_ZERO` |
| `generator: constant` | `IC_MODE_CONSTANT_CARTESIAN` |
| `generator: streamwise_constant` | `IC_MODE_CONSTANT_STREAMWISE` |
| `generator: poiseuille` | `IC_MODE_POISEUILLE` |
| `mode: file` or `generator: ic_gen` | `IC_MODE_FILE` |

@ref PopulateInitialUcont is the fresh-solve dispatcher. Built-in modes reuse
@ref SetInitialInteriorField and @ref UniformCart2Contra. File mode reuses @ref ReadFieldData;
when its field selector is `Ucat`, @ref Cart2Contra converts the loaded vector field to `Ucont`.
After that point, the existing finalization path treats every IC source identically.

@section p33_restart_modes_sec 5. Authority and Restart Branches

In @ref InitializeEulerianState "InitializeEulerianState":

- `eulerian_field_source=solve` and `StartStep == 0` consumes `initial_conditions`.
- `eulerian_field_source=solve` and `StartStep > 0` consumes the restart path.
- `eulerian_field_source=load` consumes the existing load path.
- `eulerian_field_source=analytical` consumes the analytical initializer.

Thus `eulerian_field_source` supersedes `initial_conditions`; the launcher does not materialize
a configured file or `ic_gen` artifact when another source has authority.

Operational note:

- `StartStep` identifies the saved restart state being loaded, not the first new step to compute.
- If a run completed through step `N`, restart with `start_step: N`; the first newly advanced step will be `N+1`.

@section p33_particle_link_sec 6. Particle Initialization Relation

Particle initialization is configured in `case.yml -> models.physics.particles`, but executed by a separate subsystem.

For full particle mode and restart details, use:

- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**

@section p33_checks_sec 7. Practical Checks

After startup, confirm:

- banner line shows expected field initialization mode,
- no warnings about missing inlet face (for inlet-driven setups),
- first output step has non-empty `Ucat`/`Ucont` fields.

Common pitfalls:

- using `Poiseuille` in strongly non-rectangular topology and expecting a textbook cylindrical profile,
- supplying a bulk/mean velocity to Poiseuille mode when the current implementation expects `Vmax`,
- forgetting that initialization sets the interior only; boundary handlers then overwrite face values,
- providing a PETSc vector whose size does not match the target DM,
- omitting `flow_direction` when the domain is fully periodic (no INLET face),
- comparing `u_physical` directly to `Ucont` without accounting for metric-face scaling.

@section p33_refs_sec 7. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 32_Analytical_Solutions**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 50_Modular_Selector_Extension_Guide**

@section p33_wall_spectral_sec Wall-bounded spectral seeds

@htmlinclude generated/capability_inventory_initial_file_generator.html

`channel_spectral_velocity` uses one wall axis; `duct_spectral_velocity` uses two.
The remaining axes are periodic, with driving allowed only on `streamwise_axis`.
Both use `Ucat` on one Cartesian block. Uniform spacing is required along periodic
axes; wall-normal axes may stretch. `grid_gen` provides an inspectable staged grid.

```yaml
mode: generated
generator: channel_spectral_velocity
params:
  field: Ucat
  seed: 12345
  wall_axes: [j]
  streamwise_axis: k
  bulk_velocity: 1.0
  perturbation_rms: 0.1
  wall_modes: 3
  spectrum: {type: k4_exponential, k0: 4.0, k_cut: 8.0}
  initial_spectra:
    - {task: plane_spectrum, axes: [i, k], fixed_indices: {j: 32}, subtract_mean: sample}
```

Parameters use solver nondimensional units. `perturbation_rms` is
`sqrt(volume_mean(u'^2+v'^2+w'^2)/3)`, not the RMS of each component independently.
The spectral envelope shapes a seeded vector potential along periodic directions;
wall basis functions are `sin(n*pi*t)*sin(pi*t)`, evaluated at physical cell
centers: the least envelope that keeps the potential and its wall-normal derivative
zero at the wall, so the seed has the physical near-wall behaviour, wall-parallel
`u' ~ y` and wall-normal `v' ~ y^2`. A curl using the separable face-average divergence operators produces a
perturbation in the discrete divergence nullspace, including zero wall-face flux.
Periodic modes at or above one third of the cell count are removed. Wall mode count
plus one (the highest sine index added by the envelope) must fit below half the wall-axis cell count; stretched
meshes still need a physical resolution check. The envelope is not a prescribed final
isotropic velocity spectrum: the curl, wall functions, and stretching change it.

The mean is a parabola on each wall axis, multiplied for a duct and normalized with
cell-volume weights to the requested bulk velocity. The duct product is a startup
profile, not the exact laminar rectangular-duct solution. Perturbation net flux is
zero. Match the driven boundary's target flux to bulk velocity times cross-section.
Wall dummy values are odd reflections, periodic dummies wrap. Runtime retains
ownership of boundary enforcement and Cartesian/contravariant reconstruction.

`initial_spectra` is a required nonempty list using the plane/line task syntax from
@ref p10_spectra_sec. Initial measurements cannot use a statistics-window mean.
Staging measures the generated Ucat, before runtime reconstruction. Compare against
step-zero checkpoint spectra to distinguish initialization reconstruction from
subsequent physical evolution. The summary and spectra use the canonical metrics and
spectra directories, including through precomputed asset publication and reuse.
There is no per-timestep injection, and transition/sustained turbulence is not guaranteed.

@subsection p33_cap_gen_ic_gen_sub ic_gen

@anchor p33_cap_gen_ic_gen

**Identity.** `properties.initial_conditions.generator: ic_gen` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the existing file IC runtime path.

**What it does.** Expression-based Ucat or Ucont generation on the staged grid.

**When to choose it.** Choose explicit analytic expressions; use a spectral provider for randomized velocity seeds.

**Parameters it owns.** Requires config_file; field selects Ucat/Ucont.

**Interactions.** The existing PICGRID coordinates, expression evaluator and PETSc writer own this route.

**Diagnostics.** Malformed expressions and missing inputs fail before solver launch.

**Evidence.** Unit verified - `make test-python` covers the generator and the CLI file
route. Analytically verified - `initial-conditions-2026-09-21`: staged `Ucat` and `Ucont`
expressions equal their definitions exactly, and a `file` source re-reading the staged
vectors reproduces the run bitwise.

**Limitations.** No divergence or turbulence property is implied by arbitrary
expressions: whatever the expression says is what is staged. Measured on uniform
Cartesian grids; expressions on curved grids are not covered.

@subsection p33_cap_gen_spectral_random_velocity_sub spectral_random_velocity

@anchor p33_cap_gen_spectral_random_velocity

**Identity.** `properties.initial_conditions.generator: spectral_random_velocity` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the existing file IC runtime path.

**What it does.** The existing seeded triply periodic random spectral velocity.

**When to choose it.** Choose this for DIT; channel and duct providers supply wall-compatible seeds.

**Parameters it owns.** seed, random, spectrum, projection, normalization and remove_mean retain their existing contracts.

**Interactions.** Requires six geometric-periodic faces and a uniform Cartesian grid; restart remains authoritative.

**Diagnostics.** The IC summary reports selected-operator divergence and realized energy; shell spectra measure the generated field.

**Evidence.** Unit verified - `make test-python`. Production exercised in
`examples/decaying_isotropic_turbulence`. Analytically verified -
`initial-conditions-2026-09-21`: fluctuation RMS 1.000000 and kinetic energy 1.500000 as
requested, the solver's own step-0 divergence 3.1e-14 against a flux scale of 4.3 under
`operator: picurv_discrete`, the shell spectrum within 4.7% (energy-weighted) of the
`k4_exponential` envelope with nothing above `k_cut`, and identical seeds reproducing
identical fields.

**Limitations.** A continuum solenoidal field need not be discretely solenoidal - choose
`operator: picurv_discrete` when the solver's own divergence must vanish - and runtime
reconstruction changes the energy from the staged value, which the summary predicts
separately. A startup field, not developed turbulence.

@subsection p33_cap_gen_channel_spectral_velocity_sub channel_spectral_velocity

@anchor p33_cap_gen_channel_spectral_velocity

**Identity.** `properties.initial_conditions.generator: channel_spectral_velocity` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the existing file IC runtime path.

**What it does.** A flux-normalized channel mean plus a discrete-curl spectral perturbation.

**When to choose it.** Choose for two no-slip walls; use duct_spectral_velocity for four walls.

**Parameters it owns.** The wall-bounded parameter block above applies with one wall_axes entry; initial_spectra is required.

**Interactions.** Single-block 3D Cartesian grid, two periodic axes, one no-slip pair, Ucat file staging.

**Diagnostics.** Summary reports bulk velocity, perturbation RMS and discrete divergence; selected plane/line spectra measure the seed.

**Evidence.** Unit verified - `make test-python`. Production exercised in
`examples/turbulent_channel`. Analytically verified - `wall-spectral-ic-2026-09-24`: on a
stretched 32 x 32 x 64 channel the volume-averaged streamwise velocity equals the
requested bulk to 2e-16, the mean profile equals the flux-normalized `4t(1-t)` parabola to
2.0e-15, the perturbation RMS is exactly the requested 0.1 with zero component means, the
dummy layers are exact odd reflections so the wall faces carry zero velocity, the
runtime's step-0 flux field has a cell divergence of 1.4e-17 against a flux scale of
4.1e-2, and the wall-parallel perturbation grows linearly from the wall (first-to-second
cell RMS ratio 2.42, the discrete curl's value for a potential ~ y^2).

**Limitations.** A startup construction: amplitude and resolution do not establish
sustained turbulence, and the spectrum shapes the vector potential rather than the
velocity. Measured on a reduced grid on one rank, not at the shipped production size.

@subsection p33_cap_gen_duct_spectral_velocity_sub duct_spectral_velocity

@anchor p33_cap_gen_duct_spectral_velocity

**Identity.** `properties.initial_conditions.generator: duct_spectral_velocity` dispatches through `GENERATED_IC_PROVIDERS` into `generators/ic.gen`, then the existing file IC runtime path.

**What it does.** A flux-normalized product-parabola mean plus a discrete-curl spectral perturbation.

**When to choose it.** Choose for rectangular/square ducts with four walls; channel_spectral_velocity has only one wall axis.

**Parameters it owns.** The wall-bounded parameter block above applies with two wall_axes entries; initial_spectra selects streamwise lines.

**Interactions.** Single-block 3D Cartesian grid, one periodic axis, two no-slip pairs, Ucat file staging.

**Diagnostics.** Summary reports bulk velocity and discrete divergence; actual selected line spectra report sample energy and Parseval residual.

**Evidence.** Unit verified - `make test-python`. Production exercised in
`examples/periodic_test/driven_duct`. Analytically verified - `wall-spectral-ic-2026-09-24`:
on a doubly stretched 32 x 32 x 64 duct the volume-averaged streamwise velocity equals the
requested bulk exactly, the mean equals the flux-normalized product of `4t(1-t)` parabolas
to 3.8e-15, the perturbation RMS is exactly the requested 0.1, all four walls carry zero
velocity through exact odd reflection, the runtime's step-0 flux field has a cell
divergence of 1.0e-17 against a flux scale of 2.7e-2, and near-wall growth is linear on
both wall axes (ratios 2.42 and 2.44).

**Limitations.** The product-parabola startup mean is not the exact laminar duct solution,
and a startup construction establishes no sustained turbulence. Circular pipes and
immersed boundaries are outside this provider. Measured on a reduced grid on one rank.
