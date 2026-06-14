@page 33_Initial_Conditions Initial Condition Modes

@anchor _Initial_Conditions

This page documents how PICurv initializes Eulerian velocity fields and particles at startup.
It covers both user-facing YAML inputs and the C implementation path that applies them.

@tableofcontents

@section p33_overview_sec 1. Where Initialization Happens

Startup sequence:

1. `scripts/picurv` resolves built-in, file, or external-generator YAML.
2. File-backed ICs are staged in the filename layout expected by @ref ReadFieldData.
3. @ref InitializeEulerianState chooses fresh solve, restart load, or analytical initialization.
4. On a fresh solve, @ref PopulateInitialUcont populates `Ucont`; existing finalization then applies
   boundary values and derives the remaining velocity state.
5. If particles are enabled, @ref InitializeParticleSwarm runs independently after Eulerian setup.

@section p33_euler_sec 2. Eulerian Field Initialization (`properties.initial_conditions`)

The canonical contract has two modes: `generated` and `file`. Generated ICs may use a built-in
C generator or the repository `scripts/ic.gen` utility. Both produce the same solver-facing result: one
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
      script: tools/custom_ic.py  # optional; defaults to scripts/ic.gen
      config_file: config/initial_conditions/expression.cfg
      output_file: config/initial_condition.generated.dat
```

The launcher invokes `scripts/ic.gen` by default, or the optional case-relative/absolute
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
The repository generator requires a staged PICGRID, so normal use selects
`grid.mode: file` or `grid.mode: grid_gen`. A custom `params.script` may define
its own behavior for `programmatic_c`.

@section p33_euler_modes_sec 3. C Runtime Modes and Entry Points

The launcher maps the contract to one @ref InitialConditionMode:

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

@section p33_restart_modes_sec 4. Authority and Restart Branches

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

@section p33_particle_link_sec 5. Particle Initialization Relation

Particle initialization is configured in `case.yml -> models.physics.particles`, but executed by a separate subsystem.

For full particle mode and restart details, use:

- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**

@section p33_checks_sec 6. Practical Checks

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

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Initial Condition Modes** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
