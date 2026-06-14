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

There are three IC modes (Zero, Constant, Poiseuille) and two Constant sub-modes (cartesian
and curvilinear). The sub-mode is inferred from which velocity keys are present — no explicit
`coordinate_system` key is needed.

**Zero** — set everything to zero:

```yaml
properties:
  initial_conditions:
    mode: "Zero"
```

**Constant, cartesian** — three Cartesian components; works for any domain:

```yaml
properties:
  initial_conditions:
    mode: "Constant"
    u_physical: 1.5   # x-component (physical units)
    v_physical: 0.0
    w_physical: 0.0
```

**Constant, curvilinear** — single streamwise speed; requires `flow_direction` or an INLET face:

```yaml
properties:
  initial_conditions:
    mode: "Constant"
    velocity_physical: 1.5        # streamwise speed (physical units)
    flow_direction: "+Zeta"       # required when no INLET face is present
```

Valid `flow_direction` tokens: `+Xi`, `-Xi`, `+Eta`, `-Eta`, `+Zeta`, `-Zeta`.

**Poiseuille** — separable parabolic profile in the two cross-stream directions:

```yaml
properties:
  initial_conditions:
    mode: "Poiseuille"
    peak_velocity_physical: 1.5   # centerline speed (physical units)
    flow_direction: "+Zeta"       # required when no INLET face is present
```

`picurv` mapping:

| YAML key | CLI flag | Notes |
|----------|----------|-------|
| `mode: "Zero"` | `-finit 0` | velocity components ignored |
| `mode: "Constant"` + `u/v/w_physical` | `-finit 1 -ucont_x … -ucont_y … -ucont_z …` | cartesian sub-mode; `flow_direction` forbidden |
| `mode: "Constant"` + `velocity_physical` | `-finit 1 -ic_coordinate_system 1 -ic_velocity_physical …` | curvilinear sub-mode |
| `mode: "Poiseuille"` | `-finit 2 -ic_velocity_physical …` | `peak_velocity_physical` is the centerline speed |
| `flow_direction: "+Zeta"` | `-flow_direction 4` | token-to-int map: `+Xi=0,-Xi=1,+Eta=2,-Eta=3,+Zeta=4,-Zeta=5` |

`flow_direction` resolution order (C side):
1. Derived from identified INLET face (if present).
2. Explicit `-flow_direction` CLI flag (set by `flow_direction` key).
3. Error — one of the above is required for curvilinear Constant and Poiseuille.

Launcher-side contract:

- `mode` must be set explicitly.
- For `Zero`, velocity keys may be omitted.
- For `Constant`, provide either `u/v/w_physical` (cartesian) **or** `velocity_physical` (curvilinear), not both.
- For cartesian `Constant`, `flow_direction` is forbidden.
- For curvilinear `Constant` and `Poiseuille`, `flow_direction` is required when no INLET face exists.
- When both `flow_direction` and an INLET face are present (Poiseuille only), their axes must agree.

C-side entry points:

- @ref InitializeEulerianState
- @ref SetInitialInteriorField
- @ref Cart2Contra
- @ref FieldInitializationToString

@section p33_euler_modes_sec 3. Eulerian Mode Details

`-finit = 0` (Zero):

- sets all interior contravariant velocity components to zero.

`-finit = 1` (Constant Velocity) — two sub-modes selected by `-ic_coordinate_system`:

- **Cartesian** (`-ic_coordinate_system 0`, default): `Cart2Contra` dots the Cartesian vector
  `(-ucont_x, -ucont_y, -ucont_z)` with face-area metric vectors (`Csi`, `Eta`, `Zet`) at every
  interior node, filling all three contravariant components correctly.  The old
  `identifiedInletBCFace`-based path is no longer used for this mode.

- **Curvilinear / streamwise** (`-ic_coordinate_system 1`): sets only the normal contravariant
  component along the streamwise axis (derived from `flow_direction`), scaled by the face-area
  magnitude.  The two cross-stream components are left at zero.

`-finit = 2` (Poiseuille-like Normal Velocity):

- computes a separable parabolic profile in the two cross-stream index directions,
- applies zero-clamping to suppress roundoff negatives near corners,
- `peak_velocity_physical` is the centerline speed (`Vmax`), not the bulk/mean average.

@section p33_euler_formula_sec 4. Contravariant Initialization Note

Contravariant velocity (flux) at a face equals the physical velocity dotted with the face-area
metric vector:

\f[
U_\xi = \mathbf{v} \cdot \mathbf{A}_\xi, \quad
U_\eta = \mathbf{v} \cdot \mathbf{A}_\eta, \quad
U_\zeta = \mathbf{v} \cdot \mathbf{A}_\zeta.
\f]

`Cart2Contra` performs this dot product over all owned interior nodes for the cartesian Constant
mode and is also used by the analytical solver.  For curvilinear / Poiseuille modes the
streamwise normal component is computed directly from the face-area magnitude and the physical
speed.

Because face areas scale with mesh spacing and curvature, the same physical speed produces
different `Ucont` magnitudes on stretched or curved meshes — this is expected and correct.

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

- using `Poiseuille` in strongly non-rectangular topology and expecting a textbook cylindrical profile,
- supplying a bulk/mean velocity to Poiseuille mode when the current implementation expects `Vmax`,
- forgetting that initialization sets the interior only; boundary handlers then overwrite face values,
- mixing `u/v/w_physical` and `velocity_physical` in the same `initial_conditions` block,
- omitting `flow_direction` when the domain is fully periodic (no INLET face),
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
