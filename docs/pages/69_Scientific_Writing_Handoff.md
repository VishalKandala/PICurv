@page 69_Scientific_Writing_Handoff Scientific Writing Handoff

@anchor _Scientific_Writing_Handoff

@pagemeta{Reference, Owner and collaborators, Open work package}

Six foundational pages describe the numerical core of PICurv and are substantially
thinner than the pages describing its newest features. That imbalance is deliberate
to leave in place for now: the physics belongs to the owner, and agent-authored
derivations would be worse than honest gaps.

This page records what each of those pages owes, so the eventual writing session is a
writing session rather than a scoping one.

@tableofcontents

@section p69_state_sec 1. Current State

| Page | Lines | Owes |
|---|---:|---|
| **@subpage 19_Nondimensionalization** | ~89 | The scaling set, which quantities are scaled where, and the YAML/C boundary |
| **@subpage 22_CURVIB_Method** | ~95 | Discrete curvilinear formulation; the immersed-boundary scope actually implemented |
| **@subpage 23_Fractional_Step_Method** | ~84 | Discrete form of each stage; accuracy and stability properties |
| **@subpage 26_Walking_Search_Method** | ~104 | The algorithm, its failure modes, and its termination guarantees |
| **@subpage 28_IEM_and_Statistical_Averaging** | ~114 | The mixing model and the averaging operators actually applied |
| **@subpage 34_Particle_Model_Overview** | ~130 | Physical model behind each displacement term |

For comparison, the three most recently written feature pages total roughly 1,400
lines. The imbalance is in depth, not page count.

@section p69_contract_sec 2. What Each Page Owes

Per the content contract agreed for foundational pages, each should carry:

1. Scope and assumptions — what the method covers and what it does not.
2. Governing equation in continuous form.
3. The **discrete** form actually implemented, not the textbook form.
4. Field and coordinate conventions, tied to
   **@subpage 56_Field_Identity_and_Layout_Catalog**.
5. Timestep or call sequence, tied to the diagrams already in place.
6. Boundary treatment.
7. Configuration controls that affect it.
8. Diagnostics that reveal whether it is behaving.
9. Validation evidence, using the facets in **@subpage 62_Capability_Status_Vocabulary**.
10. Known limitations.
11. Literature references and implementation anchors.

@section p69_ready_sec 3. Scaffolding Already In Place

The structural work has been done, so none of it needs doing during the writing pass:

- **Diagrams.** Fractional-step sequence, field layout, multigrid hierarchy, and
  particle chronology are drawn, embedded, and marked for owner validation of their
  scientific interpretation.
- **Page metadata and status vocabulary** exist, so a partially-written page can
  declare its status honestly rather than implying completeness.
- **The evidence matrix** (**@subpage 66_Evidence_Matrix**) now covers every public
  capability family. What it records is *declared sources*, not verified results, and
  no capability yet claims an analytical or reference-validation facet. It shows where
  the scientific evidence is thin; it does not supply any.
- **Cross-reference targets** are stable: the field catalog, capability entries, and
  glossary can be linked to without inventing anchors.
- **Corrected factual baseline.** Page 34's displacement equation and step ordering
  were verified against `src/runloop.c` and `src/ParticleMotion.c` and corrected, so
  the writing pass starts from accurate mechanics.

@section p69_notdo_sec 4. What Was Deliberately Not Done

No derivations, formulations, or physical explanations were written into these pages
by an agent. No placeholder TODO sections were published. Where a page is thin, it is
visibly thin rather than padded — that is the honest state, and the generic filler
removed at the start of this overhaul is the cautionary example.

@warning Do not let an agent fill these pages unsupervised. The failure mode is not
obvious nonsense; it is plausible, well-formatted prose that is subtly wrong about a
discretization, which is far harder to catch in review than an empty section.

@section p69_open_sec 5. Adjacent Open Items

These are not documentation tasks but block or shape the documentation:

- **Dynamic Smagorinsky.** Implementation correction pending; documentation scoped and
  deliberately deferred. Record: `tests/tooling/capability_scope_records.json`.
- **Periodic wall-bounded convergence.** Requires re-characterization at current
  `HEAD`; see @ref p54_driven_limits_sub.
- **Explicit RK4 evidence gap.** No positive-path harness; recorded in `src/guide.md`
  and visible in the evidence matrix.
- **Example cost measurement.** The catalog carries the column; values are unmeasured.
- **Analytical evidence mapping.** The particle verification examples do compare
  against exact solutions, but that evidence has not been mapped onto capability rows.

@section p69_related_sec 6. Related Documentation

- **@subpage 21_Methods_Overview** — the map across all method pages
- **@subpage 64_Documentation_Extension_Framework** — the contracts these pages follow
