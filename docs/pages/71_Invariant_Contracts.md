@page 71_Invariant_Contracts Invariant Contracts

@anchor _Invariant_Contracts

@pagemeta{Reference, Maintainers and contributors, 13 enforced of 20 registered}

PICurv documentation covers two kinds of truth, and only one of them is a choice.

**Capabilities** are things a user selects — a boundary handler, a momentum solver.
They are governed by the capability system in
**@subpage 64_Documentation_Extension_Framework**.

**Invariants** are guarantees the system maintains regardless of what anyone selects
— where a run writes its files, what the startup banner prints, which PETSc options
the runtime ingests. Nobody chooses them, so no selector changes when they change,
and every capability gate stays green while the documentation silently goes stale.

This page is the register of the second kind.

@tableofcontents

@section p71_why_sec 1. Why This Exists

A stress test asked what would happen if run output moved to a centralized location.
The answer was uncomfortable: no enum changes, no selector changes, no capability
entry changes — and roughly 120 path literals across 23 files quietly become wrong,
with every gate still passing.

That is the failure mode this register closes. The lesson generalized: **documentation
contracts must cover what the system guarantees, not only what the user chooses.**

@section p71_status_sec 2. Contract Status

| Status | Meaning |
|---|---|
| `enforced` | An extractor and a blocking checker exist; drift fails the build |
| `tracked` | Sources, owner, and dependent pages are recorded, so a change is *locatable*; no automated check |
| `planned` | Identified as important, but no reliable extraction is possible yet |
| `retired` | Kept only for migration and changelog purposes |

@warning `tracked` is not `verified`. It means a change is findable, not that anyone
would be told. Most invariants in this project are currently tracked.

@section p71_register_sec 3. The Register

Run `make audit-contracts` for the live list. Currently **13 enforced, 6 tracked,
1 planned**.

**Enforced.** PETSc option ingress; user-facing reporting grammar; run artifact
**path** topology; logical run-path locators in prose; the public capability surface;
the capability family census; published-site integrity; the generated CLI reference;
page-type coverage; named public choice sets; subsystem lifecycle obligations;
documentation freshness attestation; the field identity and layout inventory.

**Tracked.** Run artifact **lifecycle semantics**; checkpoint bundle schema; units and
non-dimensionalization; run/study/submission/storage manifest schemas; CLI structured
output schemas; durable file formats such as PICGRID and PICSLICE.

@warning `field.identity_and_layout` is enforced over the catalog's **identity and
layout** only: `audit_field_catalog.py` compares the field names and their layout
assignment in `56_Field_Identity_and_Layout_Catalog` section 6, and the layout
vocabulary in section 4, against `src/field_catalog.c`,
`src/particle_field_catalog.c`, and `include/field_catalog.h`. Degrees of freedom, DM
family, synchronization class, availability conditions, capability flags, and whether
the runtime actually allocates a field are **not** checked and remain human review.

@warning The artifact contract is deliberately **split**. `run.artifact_topology` is
enforced and covers only the planned path set: the extractor compares what the CLI
plans against a recorded snapshot across default, solve-only, and custom-directory
scenarios, and every planned path must map to a declared logical identity.
`run.artifact_lifecycle` is **tracked, not verified** — ownership, mutability,
atomicity, retention, restart role, containment, and destructive scope are declared by
hand and reviewed by humans. Nothing checks them against the sources, and the
certificate must not describe them as verified.

**Planned.** Stable runtime error identifiers, which need runtime work before the
documentation can point at anything.

@subsection p71_freshness_sub 3.1 Freshness Is a Separate Axis

Being enforced says a contract is checked against its sources. It says nothing about
whether the *prose* describing it was read since those sources last changed. That is
what `make audit-freshness` answers, from `tests/tooling/freshness_manifest.json`.

| Tier | What is fingerprinted | On change |
|---|---|---|
| Hard | A normalized semantic artifact - generated argparse structure, the capability inventory, the topology snapshot, a registry | Blocking. The owning and dependent pages are named and must be reviewed |
| Soft | Raw sources where no semantic extractor exists, such as a solver implementation file | Advisory. Reported with the pages to review; does not fail the build |
| Integrity | A malformed manifest, or a watched path that no longer exists | Always blocking - coverage is broken, so the absence of suspicion proves nothing |

A soft surface may be **promoted** to blocking with a stated reason. One is:
`run.containment_guard`, because the code it watches decides whether a directory is
recursively deleted.

Hard surfaces fingerprint normalized artifacts precisely so that a comment or
docstring change does not fire them. Soft surfaces have no such filter, which is the
price of watching something no extractor understands.

Re-attestation is `make attest-freshness ARGS="<surface-id>"`. It records that a
person or agent compared the owning pages against the sources - **not** that the
comparison found nothing, and not that the page is correct. Nothing stores a commit
SHA inside a page; git already records who re-attested and when.

@warning A surface may also be **never attested**. That is reported and does not
block for soft surfaces: it states honestly that no page has yet been compared against
that source. Eight implementation surfaces are in that state now, most of them the
physics pages whose scientific review is still outstanding.

@note The highest-value promotion from `tracked` to `enforced` is
**units and non-dimensionalization**. Unit errors are silent — wrong by a constant
factor, plausible output, no crash — which makes them exactly the class where an
automated check earns the most.

@section p71_topology_sec 4. Artifact Topology

The one invariant given a full contract so far is where a run puts its files.
`tests/tooling/artifact_topology.json` names 26
logical identities — `run.root`, `run.config`, `run.config.history`, `run.config.active`, `run.post_recipes`, `run.post_recipe.files`, `run.control`, `run.config.inputs`, `run.inputs`, `run.input.roles`, `run.asset_lock`, `run.input.files`, `run.runtime_logs`, `run.runtime_logs.files`, `run.scheduler`, `run.scheduler.files`, `run.solver_output`, `run.checkpoints`, `run.analysis`, `run.analysis.metrics`, `run.analysis.statistics`, `run.analysis.spectra`, `run.analysis.plots`, `run.visualization`, `run.visualization.recipe`, `run.manifest` — covering the run root, configuration revisions,
materialized input roles, runtime logs, scheduler state, solver output, checkpoints,
analysis families, recipe-specific visualization, and the manifest. Each identity
records its resolution rule, owner, writers and readers, mutability, retention,
restart role, and destructive scope.

Narrative pages should refer to these **logical identities**, and use `<run.config>`
notation rather than repeating physical paths. The concrete layout is extracted from
the CLI's own dry-run planner and normalized, so the recorded shape cannot drift from
what the tool actually plans.

@subsection p71_notation_sub 4.1 How Prose Names These Paths

Narrative documentation uses logical names because they express ownership and storage
policy, not because users may redirect the directories. `make audit-path-literals`
enforces that notation, and these counter-examples sit in a code block because the
audit would otherwise reject the page that defines the rule:

```text
write                                     instead of
<run.runtime_logs>/Runtime_Memory.log     logs/Runtime_Memory.log
<run.solver_output>/checkpoints/          output/checkpoints/
<run.config>/case.yml                     runs/<run_id>/config/case.yml
<run.runtime_logs>/...                    <run_dir>/logs/...
<repo>/logs/build.log                     logs/build.log
```

The logical identity lets storage and lifecycle documentation say whether an artifact
is essential, restart-bearing, raw, or derived without duplicating a physical path.
The `<repo>/` prefix distinguishes repository build/test logs from run-owned logs.

Runnable command blocks may keep concrete paths - a reader has to type something - and
the changelog is exempt, because rewriting history would falsify it.

@subsection p71_isolation_sub 4.2 Fixed Run Trees and Shared Immutable Assets

`isolated_run_tree` is the enforced mutable-output layout. Each initialized run owns
its complete `config`, `inputs`, `output`, `logs`, and `scheduler` tree. YAML cannot
redirect those homes to peer `diagnostics` or `results` directories. Isolation between
runs comes from the run id namespace and fixed run-local paths.

| Model | Meaning |
|---|---|
| `isolated_run_tree` | Mutable run state is namespaced and self-contained. **Enforced.** |
| `shared_content_addressed` | Immutable workspace assets are shared by digest and materialized read-only into a run. |

@note **The fixed names are enforced, not merely defaulted.** A configurable
`io.directories` surface once let a monitor point the log directory anywhere, and on a
fresh solve the C runtime calls `PetscRMTree` on that directory — so an escaping value
meant the solver recursively deleted whatever was there. The surface was removed rather
than guarded further, but the guards it produced were kept: a hand-edited control file,
or one staged by an older version, still has to get past them.

Containment is checked at six points:

1. configuration validation (`validate`, `run`, `submit`), which rejects a monitor
   still carrying the removed `io.directories` surface;
2. the reserved-flag guard that stops raw PETSc passthrough in any of case, solver, or
   monitor from setting `-log_dir`, `-output_dir`, `-restart_dir`, `-analysis_dir`, or
   `-allow_unsafe_log_dir`;
3. the control-file emitter, which writes the canonical values last and
   unconditionally — PETSc takes the final occurrence, so a validated directory written
   there wins over anything earlier, and writing them always removes the omitted-default
   gap where an environment variable could choose the directory instead;
4. the dry-run planner, so a plan reports what a real run would refuse;
5. a submission preflight that re-reads the staged `.control` file rather than trusting
   that staging validated it, and refuses any run-owned flag whose value is not the
   canonical one; and
6. independently of all of those, the C runtime's own guard immediately before the
   delete.

The last one matters because the first five live in `picurv`, and the solver can be
launched directly on a hand-written control file. `ClassifyLogDirectory`
(`src/setup.c`) re-derives the decision from the options actually in force, resolves
the path physically — including a directory that does not exist yet, through its
longest existing prefix — and returns a typed verdict.

@warning **Residual risk.** The restart input is not containment-checked as a
destination: it is a read path materialized from a source that was itself validated.
A simulator invoked by hand with a non-canonical `-log_dir` is protected by the C
guard alone, since no preflight ran. And the canonical values are asserted in one
place, `CANONICAL_RUN_PATHS`; a future layout change that updated the constant without
re-reading its consumers would move the guards with it silently.

External source files use a different contract: `picurv inputs import <kind> <source>
--mode reference` stores a small `.reference.yml` descriptor in the workspace.
Validation fails loudly when its target is unavailable. Mutable runtime output never
uses that indirection.

@section p71_change_sec 5. When a Contract Changes

1. Change the authoritative source.
2. Regenerate: `make docs-topology`, `make docs-inventory`.
3. The relevant checker fails, naming the contract rather than the file.
4. `make review-packet CONTRACT=<id>` prints the contract's scope, sources, canonical
   and dependent pages, logical identities, uncommitted changes, and how to verify it.
   `make review-packet PAGE=<n>` does the same from a page's point of view, listing
   source *symbols* rather than whole files.
5. `make audit-path-literals` catches prose that still hardcodes a run-relative path
   anywhere in the repository — not only in the pages the registry happens to list.
5. Update the canonical documentation, then the dependent pages.
6. Re-run `make audit-contracts`.

The tooling makes a change **locatable and bounded**. It does not write the migration
note, decide whether old runs stay discoverable, or judge whether a layout change is
safe. Those remain human decisions.

@section p71_related_sec 6. Related Documentation

- **@subpage 64_Documentation_Extension_Framework** — the capability side
- **@subpage 52_Run_Artifact_Lifecycle_Contract** — what the run topology means operationally
- **@subpage 61_Storage_Management_Guide** — archiving and restoring runs
