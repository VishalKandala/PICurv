@page 44_Boundary_Conditions_Guide Boundary Conditions Guide

@anchor _Boundary_Conditions_Guide

This page is the detailed reference for boundary-condition authoring in `case.yml`.
It documents what is currently supported by the Python conductor and what is currently wired in C.

@tableofcontents

@section p44_grammar_sec 1. Boundary-Condition Grammar

`case.yml -> boundary_conditions` supports:

- single-block form: one list with 6 faces,
- multi-block form: list-of-lists, one 6-face list per block.

Face names:

- `-Xi`, `+Xi`, `-Eta`, `+Eta`, `-Zeta`, `+Zeta`

Every block must provide all six faces exactly once.

@section p44_supported_sec 2. Supported User-Facing Combinations (`picurv`)

The validator accepts these type/handler pairs:

| Type | Handler | Required params | Optional params | Units expected in YAML |
|---|---|---|---|---|
| `WALL` | `noslip` | none | none | n/a |
| `INLET` | `constant_velocity` | `vx`, `vy`, `vz` | none | m/s |
| `INLET` | `parabolic` | `v_max` | none | m/s |
| `OUTLET` | `conservation` | none | none | n/a |
| `PERIODIC` | `geometric` | none | none | n/a |
| `PERIODIC` | `constant_flux` | `target_flux` | `apply_trim` | m^3/s |

Notes:

- legacy `params.vector` and `params.velocity` are rejected; use scalar `vx/vy/vz`,
- unknown params are rejected per handler,
- `constant_flux` must be configured on both faces of the periodic pair.

@section p44_nondim_sec 3. Non-Dimensionalization Before C Input

`picurv` converts physical BC values to solver-scale values before writing `bcs.run`:

\f[
v^* = \frac{v}{U_{ref}}, \qquad
Q^* = \frac{Q}{U_{ref} L_{ref}^{2}}.
\f]

Specifically:

- `vx/vy/vz/v_max` are divided by `properties.scaling.velocity_ref`,
- `target_flux` is divided by `velocity_ref * length_ref^2`.

@section p44_periodic_rules_sec 4. Periodicity Consistency Rules

Validator checks:

1. each axis pair must match periodicity:
   - `(-Xi,+Xi)`, `(-Eta,+Eta)`, `(-Zeta,+Zeta)` are both periodic or both non-periodic.
2. driven periodic handler (`constant_flux`) must match on both faces of a periodic pair.
3. incomplete face sets fail validation immediately.

The C side also performs periodic consistency checks while deriving global periodicity in @ref DeterminePeriodicity.

@section p44_c_pipeline_sec 5. C-Side Parsing and Dispatch

Generated `bcs.run` lines use:

```text
<face> <type> <handler> key=value key=value ...
```

C ingestion path:

1. @ref ParseAllBoundaryConditions parses and broadcasts face configs.
2. `StringToBCFace`, `StringToBCType`, `StringToBCHandlerType` convert tokens.
3. `ValidateBCHandlerForBCType` enforces type-handler compatibility.
4. @ref BoundaryCondition_Create maps handler enums to concrete implementations.

`boundary_faces` is the canonical in-memory representation. Do not restore or
maintain duplicate persistent legacy BC arrays in `UserCtx`. If a future
low-level integration genuinely needs a translated view, add a narrow
on-demand adapter only when a concrete consumer exists.

Current factory-wired handlers include:

- wall no-slip,
- inlet constant velocity,
- inlet parabolic,
- outlet conservation,
- periodic geometric,
- periodic driven constant flux.

@section p44_c_gap_sec 6. Exposed vs Latent Options

Important contributor note:

- C enums include additional BC categories and handlers (for example symmetry-family and `initial_flux` enum entries),
- but the current end-to-end user contract intentionally exposes the validated subset listed above.

If you add a new BC mode, update all three layers in one change:

1. Python validator/writer (`scripts/picurv`),
2. C parser/factory/handler implementation,
3. docs and tests (fixtures + smoke checks).

@section p44_examples_sec 7. Authoring Examples

Constant-velocity inlet + outlet conservation:

```yaml
boundary_conditions:
  - face: "-Xi"
    type: INLET
    handler: constant_velocity
    params: { vx: 1.5, vy: 0.0, vz: 0.0 }
  - face: "+Xi"
    type: OUTLET
    handler: conservation
  - face: "-Eta"
    type: WALL
    handler: noslip
  - face: "+Eta"
    type: WALL
    handler: noslip
  - face: "-Zeta"
    type: WALL
    handler: noslip
  - face: "+Zeta"
    type: WALL
    handler: noslip
```

Driven periodic flux on Xi:

```yaml
boundary_conditions:
  - face: "-Xi"
    type: PERIODIC
    handler: constant_flux
    params: { target_flux: 0.02, apply_trim: true }
  - face: "+Xi"
    type: PERIODIC
    handler: constant_flux
    params: { target_flux: 0.02, apply_trim: true }
  - face: "-Eta"
    type: WALL
    handler: noslip
  - face: "+Eta"
    type: WALL
    handler: noslip
  - face: "-Zeta"
    type: WALL
    handler: noslip
  - face: "+Zeta"
    type: WALL
    handler: noslip
```

@section p44_troubleshoot_sec 8. Common Failure Modes

- `missing required key 'face'/'type'/'handler'`: malformed BC entry object.
- `Duplicate face`: same face declared twice in one block.
- `incomplete ... Missing faces`: not all 6 faces provided.
- `Inconsistent periodicity`: one side periodic and opposite side non-periodic.
- `Unknown params`: extra keys not allowed for selected handler.

@section p44_refs_sec 9. Related Pages

- **@subpage 07_Case_Reference**
- **@subpage 33_Initial_Conditions**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 13_Code_Architecture**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 50_Modular_Selector_Extension_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Boundary Conditions Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

