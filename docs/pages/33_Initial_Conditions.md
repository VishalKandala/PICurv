@page 33_Initial_Conditions Initial Condition Modes

This page documents how PICurv initializes Eulerian velocity fields and particles at startup.
It covers both user-facing YAML inputs and the C implementation path that applies them.

@tableofcontents

@section overview_sec 1. Where Initialization Happens

Startup sequence:

1. `scripts/pic.flow` converts YAML to control flags (`-finit`, `-ucont_*`, `-pinit`, etc.).
2. @ref InitializeEulerianState chooses fresh-start, load, or analytical branch.
3. For fresh starts, @ref SetInitialInteriorField initializes interior `Ucont`.
4. Boundary handlers then enforce face values and ghost-cell consistency.
5. If particles are enabled, @ref InitializeParticleSwarm runs independently after Eulerian setup.

@section euler_sec 2. Eulerian Field Initialization (`properties.initial_conditions`)

YAML:

```yaml
properties:
  initial_conditions:
    mode: "Constant"   # Zero | Constant | Poiseuille
    u_physical: 1.5
    v_physical: 0.0
    w_physical: 0.0
```

`pic.flow` mapping:

- `mode` -> `-finit` via `normalize_field_init_mode`:
  - `Zero` (or `0`) -> `0`
  - `Constant` (or `1`) -> `1`
  - `Poiseuille` (or `2`) -> `2`
- `u_physical/v_physical/w_physical` -> `-ucont_x/-ucont_y/-ucont_z` after non-dimensionalization by `U_ref`.

C-side entry points:

- @ref InitializeEulerianState
- @ref SetInitialInteriorField
- @ref FieldInitializationToString

@section euler_modes_sec 3. Eulerian Mode Details

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

@section euler_formula_sec 4. Contravariant Initialization Note

Initialization is applied to contravariant velocity components, scaled by face-area metrics:

\f[
U_n = v_n\,A_n,
\f]

where face areas are derived from metric vectors (`Csi`, `Eta`, `Zet`) and sign is aligned with inlet face orientation.

This is why equivalent physical inflow speed can map to different `Ucont` magnitudes on stretched/curved meshes.

@section restart_modes_sec 5. Eulerian Restart Branches

In @ref InitializeEulerianState:

- `StartStep == 0`:
  - `eulerian_field_source=solve`: fresh initialization path.
  - `eulerian_field_source=load`: reads initial fields from restart files.
  - `eulerian_field_source=analytical`: uses analytical initializer.
- `StartStep > 0`:
  - `load` source reloads restart fields for the requested step.
  - `analytical` source regenerates analytical field at current `(t, step)`.

@section particle_link_sec 6. Particle Initialization Relation

Particle initialization is configured in `case.yml -> models.physics.particles`, but executed by a separate subsystem.

For full particle mode and restart details, use:

- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**

@section checks_sec 7. Practical Checks

After startup, confirm:

- banner line shows expected field initialization mode,
- no warnings about missing inlet face (for inlet-driven setups),
- first output step has non-empty `Ucat`/`Ucont` fields.

Common pitfalls:

- using `Poiseuille` in strongly non-rectangular topology and expecting textbook cylindrical profile,
- forgetting that initialization sets interior only; boundary handlers then overwrite face behavior,
- comparing `u_physical` directly to `Ucont` without accounting for metric-face scaling.

@section refs_sec 8. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 32_Analytical_Solutions**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 39_Common_Fatal_Errors**
