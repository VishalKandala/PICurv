@page 64_Documentation_Extension_Framework Documentation Extension Framework

@anchor _Documentation_Extension_Framework

@pagemeta{Reference, Contributors adding features, Derived from the boundary-condition pilot}

Adding a feature to PICurv has a defined documentation shape. This page is that
shape. It was not designed in the abstract: the capability-entry contract below was
extracted from real entries written for the boundary-handler family and then
validated against the momentum-solver family, and the sections it asks for are the
ones those entries actually needed.

The framework covers three scales of change and applies the same lifecycle and
concern modules at each.

@tableofcontents

@section p64_scales_sec 1. Three Scales of Change

| Scale | Example | What you write |
|---|---|---|
| New **value** in an existing family | a boundary handler, a momentum solver, a post kernel | one capability entry |
| New **family** in an existing subsystem | a family of preconditioners, particle-wall models | a family page, then one entry per value |
| New **subsystem** | immersed boundaries, a grid-generation framework, ML integration | a subsystem charter, then families beneath it |

Pick the smallest scale that fits. Most changes are the first.

@section p64_entry_sec 2. The Capability Entry Contract

A **canonical** value gets an entry with these eight parts, in this order. Each part
is a bold run-in label, and `make audit-capability` checks that every one is present.

Worked examples: @ref p44_cap_prescribed_flow (complex, nested input),
@ref p44_cap_noslip (minimal), @ref p08_cap_newton_krylov (deep, promoted).

1. **Identity** - the full chain from YAML selector through generated token and C
   enum to the implementing function. This part carries the implementation pointer:
   a reader must be able to get from the value they typed to the code that runs.
2. **What it does** - the behavior, in a short paragraph.
3. **When to choose it** - explicitly against its siblings. This is the part
   readers most need and the part most often missing; an entry that does not
   compare is a description, not guidance.
4. **Parameters it owns** - the options exclusive to this value, with required and
   optional separated. Options shared across the family belong on the family page.
5. **Interactions** - what it requires, conflicts with, or implies, including
   anything enforced by validation.
6. **Diagnostics** - what it logs, and how a reader tells it is working. Where
   possible, name the failure signature rather than only the success one.
7. **Evidence** - which facets are established, each naming its test, example, or
   report. Use the facet vocabulary in **@subpage 62_Capability_Status_Vocabulary**.
8. **Limitations** - what it cannot do. Where the value has been promoted to its own
   page, this part defers there rather than duplicating it.

@note An earlier draft of this contract listed nine parts, separating implementation
pointers from Identity. Writing the pilot entries showed the two always belonged
together - the chain in part 1 *is* the implementation pointer - so the contract was
revised to match the reference implementation rather than leaving the two disagreeing.

@subsection p64_alias_sub 2.1 Deprecated Alias Stubs

A **deprecated alias** is not a capability and does not get the full contract. It
gets three parts, and the audit enforces exactly those:

1. **Identity** - the alias, and the canonical value it normalizes to.
2. **Status** - that it is deprecated, and why it still parses.
3. **Migration** - the exact edits to move to the canonical value.

@ref p08_cap_dual_time_picard_rk4 is the worked example. Declare the alias in
`tests/tooling/capability_families.json` with `canonical: false` and `alias_of`
naming its target; the audit then requires the stub shape instead of the full one.

Anchor each entry as `@anchor <prefix><value_slug>` using the family's declared
prefix, on its own line and outside any code block. `make audit-capability` uses
that anchor to verify coverage, so an entry without one does not count as written -
and an anchor with the family prefix that no current value claims is reported as a
stale entry.

@subsection p64_promotion_sub 2.1 The Promotion Rule

When an entry outgrows its family page, promote it to its own page and leave a stub
plus a link behind. The entry keeps all eight parts; parts 5 and 8 defer to the
promoted page rather than duplicating it.

@ref p44_cap_constant_flux is the worked example: the driven-flux controller has a
control law, a cadence, a restart contract, and worked validation cases, so those
live on **@subpage p54_geometric_periodic** while the entry stays short and points
there. Promotion is the reason a value gets its own page - not the accident of who
was working on it recently.

@section p64_family_sec 3. The Family Contract

A family page owns what the values have in common:

- what the family controls, and the selector that chooses within it,
- the public selector grammar,
- a **generated** inventory of values (never hand-maintained - see @ref p64_generated_sec),
- rules that apply across all values,
- the capability entries themselves,
- compatibility and precedence between values,
- shared diagnostics,
- the extension procedure and testing obligation.

**@subpage 44_Boundary_Conditions_Guide** is the reference implementation.

@section p64_subsystem_sec 4. The Subsystem Contract

A new subsystem needs a charter, but not all of it at once. Requirements grow with
the status the subsystem claims, so experimental work is never blocked by
documentation it cannot yet honestly write.

| Status | Minimum documentation obligation |
|---|---|
| `planned` | Purpose, intended scope, design owner, explicit not-implemented status |
| `internal` | Scope boundary, architecture boundary, dependencies, developer entry points |
| `experimental` | Configuration, selection guidance, observability, limitations, safe-use boundaries |
| `supported` | Evidence, lifecycle and restart behavior, operations, troubleshooting, examples, complete reference |
| `deprecated` | Migration path, replacement, compatibility period, removal policy |
| `known-defective` | Defect disclosure at every selector surface, a scope record, safe-use boundaries, limitations |
| `removed` | Changelog and history record, plus the rejection behavior |

The gate applies when a feature **claims the next status**, not before.

The first four statuses form a ladder and are **cumulative**: claiming `supported`
owes everything below it as well, so a subsystem cannot reach the top by skipping a
rung. `known-defective` and `deprecated` additionally owe the ladder up to the
`peak_status` they reached; `removed` documents an absence and owes only its own row.
Visibility is a separate axis: a subsystem may be `experimental` and public, but
`internal` and public is a contradiction the audit rejects.

@subsection p64_lifecycle_enforcement_sub 4.1 How the Gate Is Enforced

Each subsystem has a record in `tests/tooling/subsystem_records.json` naming its
status, its previous status, the capability families it owns, and where each
obligation is satisfied. `make audit-subsystems` fails when an obligation is
undeclared, when a cited page or anchor does not resolve, when a publicly visible
subsystem cites an unpublished page, when a concern falls outside the vocabulary in
section 5, when a `planned` subsystem already owns a `supported` capability value, or
when the declared transition is not one the lifecycle allows - a removed subsystem
cannot be revived in place, and nothing that was ever built reaches `removed` without
first being deprecated or declared defective.

`planned -> removed` is the single deliberate exception: a cancelled design was never
built, so it owes a history record and the rejection behaviour, not a migration path
for users who never had it.

An obligation is satisfied by a resolving documentation reference, by a literal value
where the obligation is a fact rather than prose, or by a **reasoned** not-applicable
declaration. A bare `N/A` is rejected. `tests/test_subsystem_lifecycle.py` holds the
negative tests: each one writes a record a careless author could plausibly produce and
asserts the audit refuses it.

@subsection p64_status_promotion_sub 4.2 A Promotion Is a Human Decision

The gate proves that the documentation a status owes exists. It cannot prove that a
subsystem is suitable for production or that its results are scientifically correct,
and those are what `supported` actually claims.

A record may therefore declare a `proposed_status` above the one it claims, with a
`promotion_rationale` saying what the reviewer is being asked to confirm. Obligations
are then checked against the **proposal**, so the writing is held to the higher bar
while the published status stays conservative until a human agrees. `make
audit-subsystems` lists every pending promotion after its summary.

@warning A promotion recorded by tooling is not a promotion. If a subsystem's status
changed because a gate passed, the gate was measuring the wrong thing.

@section p64_concerns_sec 5. Concern Modules

A feature declares which concerns apply to it. Each activated concern adds a
checklist. This is what lets the framework cover changes nobody has anticipated:
there is no "ML template", but an ML integration activates external service,
credentials, determinism, and artifact-topology concerns and inherits their
questions.

| Concern | Activates when the feature... |
|---|---|
| Numerical method | changes discretization, a solver, or a model |
| User selector | adds or changes something selectable from YAML or CLI |
| Units and nondimensionalization | introduces a physical quantity crossing the YAML/C boundary |
| Artifact topology and storage lifecycle | reads, writes, moves, archives, restores, or deletes artifacts |
| Persistent/restart state | carries state across a checkpoint |
| Determinism and reproducibility | can change results run to run, or break restart equivalence |
| MPI/distributed execution | has rank-dependent behavior or collective calls |
| External service | talks to anything outside the process |
| Security/credentials | handles secrets or authenticated access |
| Generated artifact | produces something a downstream stage consumes |
| File format | defines or changes a durable on-disk format |
| Concurrency, permissions, destructive scope | can remove or overwrite data, or shares a location |
| Backward compatibility | changes or retires an existing contract |
| Scientific verification and validation | makes a claim about physical correctness |

@subsection p64_units_sub 5.1 Why Units Is Its Own Concern

Unit and scaling errors are the classic silent failure in CFD: wrong by a constant
factor, plausible-looking output, no crash. Any feature introducing a physical
quantity must state the units expected in YAML and where non-dimensionalization
happens. @ref p44_nondim_sec is the worked example.

@subsection p64_na_sub 5.2 Declaring a Concern Not Applicable

An activated concern may be declared not applicable **with a stated reason**. That
is a real answer and is accepted as one.

@warning What is not accepted is an empty heading or filler prose. Sections written
to satisfy a template rather than a reader are exactly what produced the generic
"CFD Reader Guidance" tails removed from 51 pages during this overhaul, and
`make audit-docs-expansion` now rejects that specific failure by signature.

@section p64_generated_sec 6. Generated Surfaces Are Not Hand-Maintained

Any enumeration that can be extracted is extracted. The capability inventory is
generated from the Python validation layer - the code that actually decides what a
user may select - and never from C enums, which carry values that are not exposed
end to end. The boundary system has 16 handler constants in its enum and **8**
public handlers; documenting the enum would advertise eight capabilities that do
not exist.

```bash
make docs-inventory     # regenerate the inventory and per-family fragments
make audit-capability   # verify parity, coverage, and that the inventory is current
```

Include a family's table with the HTML include command and the generated **`.html`**
fragment:

```text
@htmlinclude generated/capability_inventory_<family>.html
```

The `.md` snapshot in the same directory is a diff-friendly record, not the include
target: Doxygen's plain include command inserts Markdown verbatim as a code block,
so including it renders the table as unformatted text.

@section p64_adding_sec 7. Adding a Capability Value: Procedure

1. Implement the value end to end - validator, generated mapping, C parser, factory
   dispatch, runtime behavior.
2. Add it to the family's public surface in `picurv_cli/core.py`.
3. Run `make docs-inventory`. The value appears in the generated table.
4. Write its capability entry on the family page, with the required anchor.
5. Declare which concern modules apply; answer them or record a reasoned N/A.
6. Add tests, and name them in the entry's Evidence part.
7. Run `make audit-capability`. Parity and coverage must pass.
8. Run `make preview-docs` and read the entry as a reader would.

Skipping step 4 fails step 7 - coverage is enforced for every registered family. That
is the intent: the documentation obligation lands at the same moment the capability
becomes selectable. Adding a value to a family that does not exist yet fails earlier
still, at `make audit-family-census`, which refuses a public selector surface no
family covers.

@section p64_related_sec 8. Related Documentation

- **@subpage 62_Capability_Status_Vocabulary** - status words and evidence facets
- **@subpage 63_Page_Type_Contract** - which page type you are writing
- **@subpage 50_Modular_Selector_Extension_Guide** - the code-side hook points
