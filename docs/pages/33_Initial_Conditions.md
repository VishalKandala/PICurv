@page 33_Initial_Conditions Initial Condition Modes

@anchor _Initial_Conditions

This page documents how PICurv initializes Eulerian velocity fields and particles at startup.
It covers both user-facing YAML inputs and the C implementation path that applies them.

@tableofcontents

@section p33_overview_sec 1. Where Initialization Happens

Startup sequence:

1. `scripts/picurv` converts YAML to control flags (`-finit`, `-ucont_*`, `-pinit`, etc.).
2. @ref InitializeEulerianState chooses fresh-start, load, or analytical branch.
3. For fresh starts, @ref SetInitialInteriorField initializes interior `Ucont`.
4. Boundary handlers then enforce face values and ghost-cell consistency.
5. If particles are enabled, @ref InitializeParticleSwarm runs independently after Eulerian setup.

@section p33_euler_sec 2. Eulerian Field Initialization (`properties.initial_conditions`)

YAML:

```yaml
properties:
  initial_conditions:
    mode: "Constant"   # Zero | Constant | Poiseuille
    u_physical: 1.5
    v_physical: 0.0
    w_physical: 0.0
```

Preferred Poiseuille form:

```yaml
properties:
  initial_conditions:
    mode: "Poiseuille"
    peak_velocity_physical: 1.5
```

`picurv` mapping:

- `mode` -> `-finit` via `normalize_field_init_mode`:
  - `Zero` (or `0`) -> `0`
  - `Constant` (or `1`) -> `1`
  - `Poiseuille` (or `2`) -> `2`
- `u_physical/v_physical/w_physical` -> `-ucont_x/-ucont_y/-ucont_z` after non-dimensionalization by `U_ref`.
- `peak_velocity_physical` (Poiseuille only) -> `picurv` infers the unique inlet axis from `boundary_conditions` and maps the scalar peak speed onto the matching `-ucont_*` component.

Launcher-side contract:

- `mode` should be set explicitly.
- for `Zero`, velocity components may be omitted and default to zero.
- for `Constant`, explicit component values are required.
- for `Poiseuille`, use either:
  - `peak_velocity_physical`, or
  - `u_physical/v_physical/w_physical`
  but not both in the same block.

C-side entry points:

- @ref InitializeEulerianState
- @ref SetInitialInteriorField
- @ref FieldInitializationToString

@section p33_euler_modes_sec 3. Eulerian Mode Details

`-finit = 0` (Zero):

- sets interior contravariant velocity to zero.

`-finit = 1` (Constant Normal Velocity):

- uses inlet-face orientation to select the normal component,
- uses `InitialConstantContra` as target normal speed,
- applies sign according to inlet on min or max face.

`-finit = 2` (Poiseuille-like Normal Velocity):

- computes a separable parabolic profile using index-space normalized cross-stream coordinates,
- applies profile in the two cross-stream directions,
- clamps small negative roundoff to zero.
- treats the inlet-aligned input component as the profile peak (`Vmax` / centerline speed), not bulk-average velocity.

This means `peak_velocity_physical` is the clearest user-facing representation for Poiseuille mode.

@section p33_euler_formula_sec 4. Contravariant Initialization Note

Initialization is applied to contravariant velocity components, scaled by face-area metrics:

\f[
U_n = v_n\,A_n,
\f]

where face areas are derived from metric vectors (`Csi`, `Eta`, `Zet`) and sign is aligned with inlet face orientation.

This is why equivalent physical inflow speed can map to different `Ucont` magnitudes on stretched/curved meshes.

@section p33_restart_modes_sec 5. Eulerian Restart Branches

In @ref InitializeEulerianState "InitializeEulerianState":

- `StartStep == 0`:
  - `eulerian_field_source=solve`: fresh initialization path.
  - `eulerian_field_source=load`: reads initial fields from restart files.
  - `eulerian_field_source=analytical`: uses analytical initializer.
- `StartStep > 0`:
  - `load` source reloads restart fields for the requested step.
  - `analytical` source regenerates analytical field at current `(t, step)`.

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

- using `Poiseuille` in strongly non-rectangular topology and expecting textbook cylindrical profile,
- supplying a bulk/mean velocity to Poiseuille mode when the current implementation expects `Vmax`,
- forgetting that initialization sets interior only; boundary handlers then overwrite face behavior,
- comparing `u_physical` directly to `Ucont` without accounting for metric-face scaling.

@section p33_refs_sec 8. Related Pages

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
