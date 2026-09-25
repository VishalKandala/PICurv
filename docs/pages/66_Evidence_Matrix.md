@page 66_Evidence_Matrix Capability Evidence Matrix

@anchor _Evidence_Matrix

@pagemeta{Reference, Scientists assessing credibility, Declared sources for every public capability}

What confidence this project claims for each capability in the families covered so far.

The table is generated from the capability registry and now covers every public
capability family the census recognises - 39 families, 129 canonical values.

@warning **Coverage is not credibility.** A complete table means every capability has
been *asked* what evidence stands behind it, not that the answers are strong. Many
rows carry a single facet, and some carry none. Read the gaps in section 3 before
citing anything here.

A method being *documented* says nothing about whether it has been *verified*, which
is what this table is for.

@tableofcontents

@section p66_reading_sec 1. How to Read It

@warning **A tick is a declared evidence source, not a verified scientific result.**
Tooling checks that each declared source exists and that the capability entry cites
the same source the registry does. It cannot check that the test or example actually
establishes the claimed property — that is a human review judgement. Read a tick as
"someone recorded this source against this claim", not as "this has been proven".

Each column is an independent evidence facet, not a rung on a ladder. Production use
does not imply analytical verification, and a benchmark comparison does not imply
restart equivalence. A capability can be heavily exercised and still never checked
against an exact solution.

A row with **no ticks is not an error**. It means *implemented only*: the code path
exists and runs, and nobody has yet established more than that. That is a legitimate
state to record and a dishonest one to hide.

Facet definitions are in **@subpage 62_Capability_Status_Vocabulary**.

@section p66_matrix_sec 2. The Matrix

@htmlinclude generated/evidence_matrix.html

Accepted spellings and latent values are omitted: they carry no evidence of their own
beyond the canonical value they resolve to.

@section p66_gaps_sec 3. Reading the Gaps

@subsection p66_pvd_evidence Physical-Time ParaView Collections (post.pipeline)

PVD indexing is a presentation option within `post.pipeline`, so it has no separate
selector-value row in the generated table. Its direct regression evidence is
`file:tests/test_cli_smoke.py`; run
`python3 -m pytest -q tests/test_cli_smoke.py -k 'paraview or lineage_enabled_post_plan'`.

| Checked behavior | Direct test |
|---|---|
| Toggling indexing preserves computational recipe identity | `test_paraview_series_is_presentation_only_for_recipe_identity` |
| Checkpoint times, including a nonzero initial time, survive repeated collection growth | `test_paraview_series_uses_checkpoint_physical_time_and_refreshes_live_output` |
| Nested A/B/C ancestry clips parents and selects the child at a fork | `test_paraview_lineage_flattens_nested_branches_and_child_wins_fork` |
| Reinitialized particles start a new collection | `test_paraview_particle_series_starts_at_branch_when_particles_reinitialize` |
| A child's post plan starts at its first owned cadence step | `test_lineage_enabled_post_plan_starts_at_child_owned_cadence` |

These tests construct checkpoint metadata and visualization fixtures; they do not
run a live concurrent solver or open ParaView. The nested test supplies nonuniform
times and frame spacing under one recipe ID. It does not establish stitching across
different post strides. `make smoke` exercises executable solver/post workflows and
the supplied templates, but has no dedicated PVD-content assertions. PVD-specific
Slurm execution, statistics carry/reset, interrupted-write recovery, pruned-checkpoint
fallback, and a real ParaView reader session remain **not verified by this evidence**.
The historical kernel measurements in the generated table do not establish these
indexing properties. See @ref p10_io_sec for the implementation's current contract.

@subsection p66_scientific_gaps Scientific Evidence Limits

Four gaps in the current table are worth naming, because they are the ones most
likely to matter:

- **Analytical facets rest on single recorded measurements, and none is `reference`.**
  87 values cite `measurement:` records - 25 distinct ones, taken between 2026-09-18 and
  2026-09-24 - covering the Picard and Explicit RK4 solvers' orders, duct, channel and
  pipe Poiseuille flow, the Poisson options, every initial-condition mode, every
  particle seeding and restart mode, both interpolation methods, the post-processing
  kernels, field statistics, shell, plane and line spectra, metric closure on every
  grid-generator feature, the workspace input import modes, and the four generated
  initial-condition providers.
  Each record states what it does not establish, and none is gated in CI. Six records
  are cited by no value: the two `not-met` ones (the Q-criterion output placement and the
  generated inlet flux on `programmatic_c` grids, both since fixed and re-measured), the
  superseded `inconclusive` drift measurement, the 2026-09-22 wall-seed measurement
  taken on the former `sin(pi t)^4` envelope (superseded by `wall-spectral-ic-2026-09-24`), and the paired drift and scalar-scatter
  measurements, which verify subsystems with no selector value to cite them. No
  capability has been compared against external reference data.
- **The grid generator's closure is metric consistency, not accuracy on every shape.**
  Every geometry, section, wall segment, path segment and transform closes a uniform
  flow to round-off in the solver's metrics, but solves have been run only on flat
  boxes, the swept circle, and a mirrored hill channel. Metric quality at a resolved
  step corner is reported by the generator and not validated against a solution.
- **The turbulence closures carry no measured facet.** Every LES and wall-function
  value is experimental; at most they carry unit coverage and, for the dynamic model, a
  production example. The LES models are detailed next.

- **All four LES models are `experimental`.** All are implemented and carry unit
  coverage in `tests/c/test_les.c`, including an analytic check of the Germano model
  tensor, a decomposition-independence check of the coefficient averaging, and checks of
  the Vreman and WALE kernels against independent evaluations. None has a validated
  coefficient magnitude: no reference-flow comparison has been run and
  gated. The check that would close the gap is decaying isotropic turbulence with
  homogeneous averaging, where `Cs(t)` should settle near 0.16-0.17. The dynamic model
  declares `examples/decaying_isotropic_turbulence` as production evidence, which records
  that the model runs there and reproduces the trends a correct implementation must show
  - not that the coefficient has been checked against a reference. The constant model
  declares none, because no shipped example selects it.

@section p66_updating_sec 4. Keeping It Current

The matrix is generated from `value_metadata[*].evidence` in
`tests/tooling/capability_families.json` and regenerated by `make docs-inventory`.

Each facet maps to source identifiers (`make:<target>`, `example:<dir>`,
`file:<path>`, or `measurement:<id>`). `make audit-capability` enforces two things:
every declared source must exist, and the capability entry's **Evidence** part must
cite the same source. It does not and cannot verify that the source establishes the
claim, so adding a facet remains a deliberate human assertion.

@subsection p66_measurements_sub 4.1 Citing a Measurement That Already Happened

The first three source kinds name something anyone can re-run on demand. That works
for unit and integration coverage, and for a shipped example. It does not work for
the evidence the `analytical`, `benchmark`, and `reference` facets ask for: a
refinement study, a decaying-turbulence coefficient check, a comparison against a
published result. Those run once, on a cluster, at a resolution nobody will repeat
casually, and the answer is a number rather than a command.

`measurement:<id>` cites such a run, recorded in
`tests/tooling/measurement_records.json`. A record states the question asked, the
date and revision it was taken at, the machine and rank layout that produced it, the
configuration that determines the result, what was measured against what threshold,
and — required, with `none` rejected — what the measurement does not establish. A
verdict of `not-met` or `inconclusive` is a legitimate record: a campaign that failed
its threshold is evidence, and deleting it invites the same run again.

`make review-packet CAPABILITY=<family>` prints the records cited by that family
before its verification commands, so a reader planning an expensive run sees what has
already been answered first. That is the point of recording them: the registry exists
to stop a question being paid for twice, and to give a status promotion something
specific to rest on.

@section p66_related_sec 5. Related Documentation

- **@subpage 62_Capability_Status_Vocabulary** — the facet vocabulary
- **@subpage 65_Example_Catalog** — which example establishes what
- **@subpage 40_Testing_and_Quality_Guide** — the test suites behind the facets
