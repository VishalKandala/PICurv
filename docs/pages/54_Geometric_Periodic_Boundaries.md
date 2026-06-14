/**
@page p54_geometric_periodic Geometric Periodic Boundaries

@section p54_scope_sec 1. Scope

PICurv supports geometric periodicity for Eulerian fields. Periodic axes are
derived exclusively from paired `PERIODIC` boundary conditions before the
DMDA is created. The resulting internal periodic flags select PETSc periodic
index wrapping; they are not user configuration options.

Particle positions are not currently wrapped across periodic boundaries.

@section p54_grid_sec 2. Grid Contract

For every active periodic axis, the two opposite nodal surfaces must match
pointwise under one nonzero constant Cartesian translation:

\f[
\mathbf{x}_{+}(p) - \mathbf{x}_{-}(p) = \mathbf{T}
\f]

for every corresponding surface point \f$p\f$. The translation may have
components in all Cartesian directions, but it may not vary over the surface.
Each periodic axis must contain at least four physical nodes. These checks are
performed independently, so mixed periodic/non-periodic domains are supported.

This admits translated curvilinear periods, including a repeating wavy
channel, while rejecting rotating, shearing, or otherwise nonmatching paired
surfaces.

@section p54_fields_sec 3. Field Synchronization

PETSc `DM_BOUNDARY_PERIODIC` wraps indices but does not translate coordinates
or resolve the ownership conventions of every staggered field. PICurv therefore
handles periodic Eulerian data by layout:

- cell-centered fields,
- persistent fields belonging to one face family,
- component-staggered fields such as `Ucont`,
- local-only staggered work fields.

Face-center and cell-center coordinate ghosts are wrapped and shifted by the
validated translation. Higher-depth QUICK stencil preparation is handled
separately from ordinary local ghost refresh.

@section p54_config_sec 4. Configuration

Configure both faces of an axis as periodic:

```yaml
boundary_conditions:
  blocks:
    - id: 0
      faces:
        - {face: "-Xi", type: PERIODIC, handler: geometric}
        - {face: "+Xi", type: PERIODIC, handler: geometric}
```

Do not add `i_periodic`, `j_periodic`, or `k_periodic` under `models.domain`;
those legacy keys are rejected.

@section p54_diag_sec 5. Diagnostics and Tests

At startup, the banner reports BC-derived periodic axes and each validated
translation. Use `make unit-periodic` for focused synchronization and geometry
coverage and `make smoke-periodic` for the real runtime periodic workflow.
Both are included in the standard `unit` and `check` gates.

@section p54_refs_sec 6. Related Pages

- @subpage 44_Boundary_Conditions_Guide
- @subpage 07_Case_Reference
- @subpage 40_Testing_and_Quality_Guide
*/
