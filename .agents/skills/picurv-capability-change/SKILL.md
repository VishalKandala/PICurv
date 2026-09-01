---
name: picurv-capability-change
description: Route new or changed PICurv solvers, preconditioners, forces, boundaries, postprocessors, CLI features, and other capabilities. Use for extension work; not narrow fixes or diagnosis.
---

# PICurv capability change

Implement the smallest complete vertical slice while keeping runtime behavior, enforced
contracts, tests, and documentation aligned. Do not create a feature-specific checklist
when an existing registry, generated artifact, guide, or audit already owns the facts.

## Classify before searching broadly

Choose the smallest applicable change class:

- a value in an existing public selector family;
- a new public selector key or family;
- a new subsystem or a subsystem lifecycle change;
- a non-selector YAML/configuration field;
- a CLI-only command, flag, or workflow; or
- an internal implementation with no new public ingress.

Use `docs/pages/64_Documentation_Extension_Framework.md` to distinguish value, family,
and subsystem work. Use `docs/pages/16_Config_Extension_Playbook.md` for non-selector
configuration and the current CLI parser/generated CLI reference for CLI-only work.
Internal work must not invent a public selector or capability entry.

Use the documentation-as-index rule in `AGENTS.md` before searching broadly. Let the
applicable registry, generated inventory, guide, and review packet produce the smallest
plausible set of source symbols, files, contracts, and tests, then verify that set in
code. If the route is empty, suspect, never attested, report-only, or names missing
symbols, switch immediately to targeted code search and report the indexing gap.

For a value in an existing selector family, begin with
`make review-packet CAPABILITY=<family-id>`. For a subsystem or lifecycle change, also
run `make review-packet SUBSYSTEM=<subsystem-id>`. These joins are complete only over
declared registry data; they do not prove code behavior. If `docs_build/xref.json` is
current, use its bounded direct and two-hop edges to choose caller sites to inspect, while
respecting the printed Doxygen call-graph caveat. A missing optional index is normal;
continue with the base route or run `make docs-xref` when caller discovery will repay the
build cost.

## Trace a live sibling

Find the nearest implemented sibling and trace the path that exists at the current
`HEAD`. Use the trust hierarchy in `AGENTS.md`; guides route exploration but do not prove
current behavior.

- For a public selector, start with its record in
  `tests/tooling/capability_families.json`. Follow `public_surface`, `parity_sources`, and
  `family_page`, then confirm acceptance, generated mapping, runtime ingestion,
  dispatch/consumer behavior, and diagnostics in code. If the registry cannot express a
  link in that chain, verify it manually and report the audit coverage gap; extend parity
  enforcement only when the tooling can represent it honestly.
- For configuration ingress, route through `config/guide.md`, the role-specific guide it
  names, `picurv_cli/guide.md`, and `docs/pages/15_Config_Ingestion_Map.md`. Update the
  ingress manifest when the enforced PETSc option surface changes.
- For CLI work, inspect the assembled parser and dispatch rather than regexing one source
  file. Treat `docs/generated/cli_reference.*` as generated output.
- For runtime work, use `src/guide.md` and `include/guide.md` to locate ownership,
  setup/teardown, dispatch, and consumers.
- For reusable configuration assets, follow `config/guide.md` and its starter-content
  contract.
- For a postprocessing product, trace the compute task and the independent output-name
  surface through producer, selectable field name, and writer mapping.

Resolve identifiers and paths from the repository instead of relying on values remembered
from an earlier session.

## Pass the reuse gate

Apply the repository-wide reuse rule in `AGENTS.md` before designing new code. Build a
small reuse map from each required operation to the existing APIs and sibling call sites
inspected. For every proposed helper, wrapper, state field, or duplicate branch, state why
direct reuse cannot satisfy the semantics or performance and whether the current owner can
be generalized instead. A second implementation of an existing responsibility is not a
valid vertical slice.

Prefer extending the authoritative owner so all applicable callers gain the improvement.
Keep feature-specific mathematics at its subsystem boundary, but route reusable mechanics
such as parsing, field access, lifecycle, I/O, synchronization, and diagnostics through
their established shared infrastructure. Do not introduce a generic abstraction unless
its contract is stable and it serves the current change without extra hot-path work. If
the reuse map finds existing bespoke callers that the generalized operation could replace,
present the consolidation scope and tradeoffs and ask the user before changing those
callers, as required by `AGENTS.md`.

## Build an impact map

Activate only the concern modules that apply in
`docs/pages/64_Documentation_Extension_Framework.md`. Explicitly decide whether the change
touches:

- user ingress, validation, normalization, or generated controls;
- runtime state, allocation, ownership, setup, teardown, or dispatch;
- units, numerical behavior, observability, or scientific evidence;
- restart/persistent state, artifacts, file formats, or postprocessing exposure;
- MPI decomposition, collectives, deterministic behavior, or rank equivalence;
- compatibility, lifecycle status, examples, or user documentation.

For every new field, first distinguish process-lifetime storage from durable/restart state.
Then decide allocation/ownership, field identity/cataloguing, checkpoint read/write,
postprocessing/output exposure, teardown, and MPI layout independently. A derived or
temporary kernel field does not automatically require all of them. For derivative kernels,
also verify units and pipeline ordering against the in-memory field state at the point of
computation.

Concern modules are an analysis checklist. Record them only where the current capability or
subsystem contract requires; do not manufacture empty sections. Treat lifecycle status as a
scientific/product claim supported by current evidence, not as a default inferred from how
new the code is.

## Execute and coordinate

Implement the narrowest end-to-end path that makes the capability real. Do not leave a
selectable placeholder or a runtime branch that users cannot reach.

When permitted and two or more lanes are genuinely independent, delegate bounded read-only
exploration such as runtime/state tracing, configuration/contract tracing, and test/document
routing. Keep one lead responsible for the impact map, edits, and integration. If write
delegation is useful, assign disjoint files. Serialize PETSc/MPI runs and any command
sharing generated outputs.

## Verify and document

Read `tests/guide.md` and select the narrowest test that answers each affected layer. Add a
real executable smoke check when unit tests cannot exercise ingress or lifecycle behavior,
and an applicable MPI check when rank behavior can differ. Run broader gates in proportion
to the change; report applicable checks that were not run.

Route documentation work through its existing owner instead of copying its procedure here:

- Selector value/family/subsystem: follow
  `docs/pages/64_Documentation_Extension_Framework.md`, then the hook map in
  `docs/pages/50_Modular_Selector_Extension_Guide.md` and the owning registry. A new
  selector key/family, rather than a value in an existing family, also invokes the family
  census.
- Non-selector configuration: follow `docs/pages/16_Config_Extension_Playbook.md`.
- CLI contract: regenerate with `make docs-cli-reference` and run
  `make audit-cli-reference`; change explanatory prose only when semantics or workflow
  changed.
- Subsystem status or obligations: update `tests/tooling/subsystem_records.json` and run
  `make audit-subsystems`.
- Before editing a documentation page or contract, run the applicable
  `make review-packet PAGE=<n>` or `CONTRACT=<id>`.

Use the scaffold/generator/audit commands named by those owners and `AGENTS.md`; never
hand-edit `docs/generated/`. If an audit, generator, or review-packet command fails, report
the tooling failure and its partial scope; do not silently replace machine evidence with a
prose claim. If visual preview is warranted, distinguish a successful build/gate from an
interactive serving process.

## Handoff

Run `make review-packet CHANGED=working-tree` before the final broad verification pass.
Treat its production-path coverage as advisory: inspect every unrouted production path
through the nearest guide and a targeted responsibility/symbol search, then report the
registry gap rather than claiming the change is fully routed.

State the change class, integration surfaces changed, tests and audits run, applicable
checks skipped, reuse decisions, any net-new abstraction and its justification, and
remaining risks. Separate behavior actually exercised from evidence merely declared in a
registry. Call out any code/document disagreement instead of choosing the more convenient
source.
