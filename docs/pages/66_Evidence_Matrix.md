@page 66_Evidence_Matrix Capability Evidence Matrix

@anchor _Evidence_Matrix

@pagemeta{Reference, Scientists assessing credibility, Declared sources for every public capability}

What confidence this project claims for each capability in the families covered so far.

The table is generated from the capability registry and now covers every public
capability family the census recognises - 39 families, 130 canonical values.

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

- **Most facets rest on single recorded measurements, and three records are `reference`.**
  125 values cite `measurement:` records - 38 distinct ones, taken between 2026-09-18 and
  2026-10-09 - covering the Picard and Explicit RK4 solvers' orders, duct, channel and
  pipe Poiseuille flow, the Poisson options, every initial-condition mode, every
  particle seeding and restart mode, both interpolation methods, the post-processing
  kernels, field statistics, shell, plane and line spectra, metric closure on every
  grid-generator feature, the workspace input import modes, and the four generated
  initial-condition providers, plus Newton--Krylov production execution on a
  144-rank turbulent channel and, with the point-block preconditioner, on a three-grid
  laminar curved duct and a wall-modelled turbulent one, and the storage compression
  levels, offload policies, and retention components against a Google Drive remote, the
  three study types run as Slurm arrays, the import of an external periodic velocity
  field, and the LES models and options against the HOM02 benchmark.
  Each record states what it does not establish, and none is gated in CI. Twelve records
  are cited by no value: six `not-met` ones (the Q-criterion output placement and the
  generated inlet flux on `programmatic_c` grids, both since fixed and re-measured,
  asset reuse for identical inputs, fixed and since confirmed on the cluster by
  `asset-lifecycle-grace-2026-10-01`, `hom02-les-accuracy-2026-10-09`, cited in the
  LES entries' limitations, and the `cabot` channel run and the second `werner` window,
  cited in those wall functions' entries), the
  superseded `inconclusive` drift measurement, the 2026-09-22 wall-seed measurement
  taken on the former `sin(pi t)^4` envelope (superseded by
  `wall-spectral-ic-2026-09-24`), and the paired drift, scalar-scatter, pending-job cancellation, and IEM
  variance-decay measurements, which verify subsystems with no selector value to cite them.
  The three `reference` records are all the wall-modelled Re_tau = 1000 channel against
  Lee & Moser (2015), judged against criteria fixed before the runs, on one grid:
  `wmles-channel-retau1000-werner-2026-10-09` backs `werner` and the dynamic-Smagorinsky
  path, `wmles-channel-retau1000-loglaw-2026-10-09` backs `log_law`, and
  `wmles-channel-retau1000-simpson-2026-10-09` backs `simpson_ik`. The earlier Lee--Moser comparison in `turbulent-channel-nk-2026-09-29`
  retains an 11.4% friction-Reynolds-number discrepancy and 29.3% excess skin friction;
  its `met` verdict concerns converged production execution, so it contributes a
  `production` facet rather than a `reference` facet. See @ref p55_channel_evidence_sub.
- **The grid generator's closure is metric consistency, not accuracy on every shape.**
  Every geometry, section, wall segment, path segment and transform closes a uniform
  flow to round-off in the solver's metrics, but solves have been run only on flat
  boxes, the swept circle, a mirrored hill channel, and the swept-square 90-degree bend
  of `humphrey-laminar-bend-nk-pointblock-2026-10-07` and
  `humphrey-turbulent-bend-wmles-2026-10-09`, whose comparisons with the Humphrey, Taylor
  & Whitelaw and Taylor, Whitelaw & Yianneskis measurements are exploratory. Metric
  quality at a resolved step corner is reported by the generator and not validated
  against a solution.
- **The LES `benchmark` facets characterize; they do not pass an accuracy test.** The
  four LES models and the averaging, clipping and width options rest on
  `hom02-les-execution-2026-10-09`: every model ran HOM02 decaying isotropic turbulence
  from Wray's DNS field to the end and removed the grid-cutoff pile-up, and the dynamic
  coefficient settled at 0.18-0.19. The accuracy criteria fixed before those runs were
  not met (`hom02-les-accuracy-2026-10-09`): every model ran 18-29% above the DNS energy
  during start-up from the DNS field, and Vreman missed the low-band spectrum criterion.
  `geometric_mean` is supported on `hom02-geometric-mean-equivalence-2026-10-09`, which
  shows it equal to `cube_root_volume` on orthogonal cells, the only cells run.
  `simpson_ik` rests on `hom02-simpson-ik-2026-10-09` and the channel record below; both
  showed its former weighting under-dissipative, and its validation is of the corrected one.
- **The wall functions rest on one channel.** `werner` and `log_law` are validated on the
  wall-modelled Re_tau ~ 1000 channel (`wmles-channel-retau1000-werner-2026-10-09`,
  `wmles-channel-retau1000-loglaw-2026-10-09`), and `werner` also runs a production bend.
  A second window of `werner` (`wmles-channel-retau1000-werner-window2-2026-10-09`, not met)
  put its friction velocity 3.2% high against the 3% bound, so its bias straddles the
  criterion. `cabot` failed the same channel (`wmles-channel-retau1000-cabot-2026-10-09`) and
  stays experimental. Unit coverage alone said nothing about any of this: until 2026-10-08
  `werner` passed it while applying a cell-averaged relation to a point velocity and
  over-predicting the channel's friction velocity by 14.5%.

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
