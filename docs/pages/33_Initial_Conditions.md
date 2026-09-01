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
      output_file: config/initial_condition.generated.dat
```

The launcher invokes `generators/ic.gen` by default, or the optional case-relative/absolute
`params.script` override, as:

```text
python <ic-generator> -c <config_file> --field Ucat|Ucont --output <output_file> --grid <grid.run> [cli_args...]
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
`<run.config>/grid.run` bridge from scalar `programmatic_settings` before invoking
`ic.gen`.

@section p33_entries_sec 3. Field Initialization Mode Entries

@htmlinclude generated/capability_inventory_initial_field_mode.html

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

**Evidence.** Implemented only.

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

**Evidence.** Implemented only.

**Limitations.** Not divergence-free in general on a curvilinear grid, so the first
pressure solve does real work.

@subsection p33_cap_poiseuille_sub Poiseuille

@anchor p33_cap_poiseuille

**Identity.** Field mode `Poiseuille` -> `-finit 2` -> `IC_MODE_POISEUILLE`.

**What it does.** Initializes with an analytic parabolic profile.

**When to choose it.** A wall-bounded channel or duct where you want a developed-looking
start rather than a transient from rest.

**Parameters it owns.** The profile scaling in the initial-condition block.

**Interactions.** Assumes a wall-bounded geometry; in a domain without the expected walls
the profile will not vanish where you expect.

**Diagnostics.** As above; inspect the step-0 field before committing to a long run.

**Evidence.** Implemented only.

**Limitations.** An analytic approximation, not a solution of the case being run. The C
enum also carries `IC_MODE_CONSTANT_STREAMWISE` and `IC_MODE_FILE`, which this selector
does not expose.

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

**Evidence.** Production exercised - `examples/flat_channel`.

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

**Evidence.** Implemented only. No shipped example selects it.

**Limitations.** Unverified by any in-tree case, and easier to get wrong because the
values are fluxes rather than velocities.

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
- attempting a file-backed IC for a multi-block case in the first implementation,
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
