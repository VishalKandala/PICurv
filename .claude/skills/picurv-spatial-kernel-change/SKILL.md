---
name: picurv-spatial-kernel-change
description: Design or refactor PICurv stencils, interpolation/scatter, and distributed field loops across layouts, metrics, ghost/periodic synchronization, and MPI ownership. Not configuration-only work or initial diagnosis.
---

# PICurv spatial kernel change

Establish the spatial data contract before editing any kernel that reads neighbors,
writes distributed fields, interpolates between layouts, or assembles local values.
Use `picurv-solver-debugging` first when the cause of an observed anomaly is still
unknown; switch here once evidence localizes it to indexing, layout, stencil, metric,
or synchronization behavior. Combine this skill with `picurv-capability-change` when
the spatial work also adds or changes a user-selectable capability.

## Resolve current layout facts

Run `make review-packet CONTRACT=field.identity_and_layout`. Use
`docs/pages/56_Field_Identity_and_Layout_Catalog.md` sections 3 through 5 to locate
the field's descriptor, layout, and synchronization concepts, then confirm them in
the declared catalog sources, the runtime DM/setup path, and the affected producer
and consumers. The contract is tracked and report-only: neither the page nor the
catalog alone proves the live allocation, ownership, or call order.

Apply the documentation-as-index rule in `AGENTS.md`: use the packet and guides to
form the smallest plausible source-and-test inspection set, not as a reason to read
the whole subsystem. An empty route, missing symbol, suspect or never-attested surface,
or report-only claim requires targeted code tracing before design; report the routing
gap rather than filling it with inference.

When a freshness surface covers the localized implementation, also run
`make review-packet SURFACE=<surface-id>`; for a solver-owned spatial path, add the
applicable `SUBSYSTEM=<subsystem-id>` packet. Use a current optional xref section to
inspect the bounded direct and two-hop callers around the declared field or
synchronization symbols. Doxygen edges are supporting source evidence, not semantic
dispatch proof; trace registry tables, callbacks, PETSc operations, and runtime selection
manually when the packet warns that they are indirect. Do not build the optional index
merely to make an otherwise complete registry route look more authoritative.

Use `docs/pages/20_Grid_Cell_Architecture_Guide.md`,
`docs/pages/54_Geometric_Periodic_Boundaries.md`, `src/guide.md`, and
`include/guide.md` only to narrow exploration. Resolve current field inventories,
shifted indices, stencil widths, metric formulas, and periodic repairs from code
rather than copying them into this skill.

## Reuse spatial infrastructure first

Before adding loop-bound, layout, periodicity, scatter, assembly, or field-access code,
apply the reuse gate in `AGENTS.md`. Inspect current callers of shared operations such as
`GetOwnedCellRange`, `FieldGetView`, `UpdateLocalGhosts`, and the layout-specific periodic
synchronizers; reuse them only when their ownership and vector-authority contracts match.
Also locate the established boundary/periodicity state rather than creating a
subsystem-local interpretation of an axis.

If a shared operation is too narrow, prefer generalizing it at its owning module and
migrating every applicable caller. Keep a specialized implementation only when stencil
semantics or measured hot-path cost makes the shared path incorrect; record that reason.
Do not hide communication or allocation inside a supposedly general numerical kernel.

## Write the transient spatial contract

For every affected field and operator, record these six decisions in the working
notes or implementation plan:

1. the exact physical entries this rank may write;
2. the entries it may read, including the representation that owns the current value;
3. the per-axis stencil reach and every physical-boundary variant;
4. the scatter, boundary repair, periodic repair, and assembly order that makes reads
   fresh and writes visible;
5. the result that must remain invariant under a different MPI decomposition; and
6. the applicable block and multigrid-level scope.

Include local-to-global semantics in the write decision: use insertion only for a
single-writer contract and accumulation only when contributions are intentionally
combined.

> Never infer a physical iteration domain directly from PETSc ownership, and never
> infer ghost freshness from the existence of a local vector.

Keep geometric indices, solver-layout indices, PETSc-owned indices, and halo indices
distinct. A physical boundary or dummy slot is not an MPI halo. Outputs normally
write owned physical entries; neighbor reads may use ghosted local entries only after
the required refresh and repairs.

If the field identity/layout, authoritative vector, or stencil reach is still unknown,
name those gaps and retain explicit branches in the plan. Do not select a synchronizer
by inference.

## Trace the execution order

Trace the live path end to end:

`producer -> authoritative global/local value -> scatter -> boundary/periodic repair -> kernel -> assembly/sync -> consumer`

At each edge, verify the actual Vec and DM, ownership range, ghost width, and call
ordering. Check that the DM halo covers the stated stencil reach for the smallest
legal local subdomain, including boundary variants. If the operator uses curvilinear
geometry, identify the correct center- or face-located metric/Jacobian at the point of
application and verify the resulting units.

Before choosing a communication helper, classify each required value as an interior
decomposition halo, a physical periodic endpoint/duplicate, a staggered normal-face
ghost, or a wider-stencil ghost. These are different obligations. A standard refresh
publishes current canonical global values across interior seams; layout-specific
periodic synchronizers also repair endpoint or dummy planes. If the producer writes a
local-only or noncanonical trial vector, publish that representation first—a catalogued
global-to-local refresh may otherwise overwrite new values with stale canonical state.

For particle-grid interpolation or scatter, also establish support ownership,
normalization, conservation, migration timing, and whether contributions cross rank
boundaries. For matrix or `MatStencil` assembly, verify that logical indices map to
the intended layout and that off-rank reach is representable. Do not assume
multiblock or multigrid-level equivalence when the current runtime and tests do not
demonstrate it.

## Implement and verify

Derive loop bounds from the stated physical domain and layout. Keep owned writes
separate from ghosted reads, make synchronization explicit at the owning orchestration
layer, and preserve existing boundary variants unless the requested change alters
their contract. Avoid a helper whose name hides a scatter, repair, or collective that
changes call-order requirements.

Read `tests/guide.md` and choose the narrowest applicable checks. A spatial change
normally needs a subset of this evidence matrix, selected by the contract above:

- constant, linear, or manufactured fields with an analytic result;
- index-unique sentinel values that expose transposition, offset, or ownership errors;
- interior points and each affected physical-boundary variant;
- nonperiodic, mixed-periodic, or fully periodic topology as applicable;
- one rank and a decomposition whose seam crosses the stencil, plus an alternate
  decomposition when the result should be identical;
- stretched or curvilinear geometry when metrics participate;
- masks or immersed boundaries when they alter support; and
- conservation, constant preservation, or normalization for particle-grid coupling.

For a stale-ghost regression, poison the local halos, invoke the production
orchestration path, and assert the downstream kernel result as well as the repaired
values. A helper-only test does not establish correct call ordering.

Compare global results rather than rank-local ordering when decomposition invariance
is the claim. State which matrix cells were applicable, which ran, and which remain
unverified.

When permitted and useful, delegate independent read-only traces for field/layout,
synchronization/order, and test/oracle selection. Give every lane the same spatial
contract; keep one lead responsible for edits and integration, and serialize shared
PETSc/MPI runs.

## Documentation and handoff

If the behavior or a declared contract changes, use the applicable review packet and
the documentation rules in `AGENTS.md`; never hand-edit generated documentation. Keep
exact layout tables and formulas in their existing code or documentation owner rather
than creating an agent-only copy.

Report the six-part spatial contract, reuse decisions, any specialized path retained,
changed execution order, tests and rank/topology matrix run, applicable cases skipped,
and residual uncertainty. Distinguish a value verified by execution from a layout fact
merely found in a report-only record.

Before handoff, run `make review-packet CHANGED=working-tree`. Any unrouted production
path requires the packet's nearest-guide fallback and targeted code search; the changed
route is advisory and never replaces the spatial execution-order trace.
