@page 71_Invariant_Contracts Invariant Contracts

@anchor _Invariant_Contracts

@pagemeta{Reference, Maintainers and contributors, 12 enforced of 20 registered}

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

Run `make audit-contracts` for the live list. Currently **12 enforced, 7 tracked,
1 planned**.

**Enforced.** PETSc option ingress; user-facing reporting grammar; run artifact
**path** topology; logical run-path locators in prose; the public capability surface;
the capability family census; published-site integrity; the generated CLI reference;
page-type coverage; named public choice sets; subsystem lifecycle obligations;
documentation freshness attestation.

**Tracked.** Run artifact **lifecycle semantics**; checkpoint bundle schema; field
identity and layout; units and non-dimensionalization; run/study/submission/storage
manifest schemas; CLI structured output schemas; durable file formats such as PICGRID
and PICSLICE.

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
`tests/tooling/artifact_topology.json` names 13 logical identities — `run.root`, `run.config`, `run.control`, `run.runtime_logs`, `run.scheduler`, `run.solver_output`, `run.checkpoints`, `run.manifest`, `run.config.inputs`, `run.post_control`, `run.runtime_logs.files`, `run.scheduler.files`, `run.visualization` — each recording its
resolution rule, owner, writers and readers, mutability, retention, restart role, and
destructive scope.

Narrative pages should refer to these **logical identities**, and use `<run.config>`
notation rather than repeating physical paths. The concrete layout is extracted from
the CLI's own dry-run planner and normalized, so the recorded shape cannot drift from
what the tool actually plans.

@subsection p71_notation_sub 4.1 How Prose Must Name These Paths

Two of the run-owned directories - the log and output directories - are
**configurable**. A page that names one concretely is therefore wrong the moment
`io.directories` is customized, which is exactly the storage-layout case this contract
exists for. `make audit-path-literals` enforces the notation, and the counter-examples
below sit in a code block because the audit would otherwise reject the page that
defines the rule:

```text
write                                     instead of
<run.runtime_logs>/Runtime_Memory.log     logs/Runtime_Memory.log
<run.solver_output>/checkpoints/          output/checkpoints/
<run.config>/case.yml                     runs/<run_id>/config/case.yml
<run.runtime_logs>/...                    <run_dir>/logs/...
<repo>/logs/build.log                     logs/build.log
```

The first three are the ordinary case: a logical identity instead of a path string.
The fourth matters because a placeholder run root still hardcodes which subdirectory
follows it. The fifth is the one that bites: the repository's own build and test log
directory is a **different** directory that no run configuration moves, and in prose it
used to be indistinguishable from a run's runtime logs. The `<repo>/` prefix is not
decoration.

Runnable command blocks may keep concrete paths - a reader has to type something - and
the changelog is exempt, because rewriting history would falsify it.

@subsection p71_isolation_sub 4.2 Default Layout, Not Enforced Isolation

`isolated_run_tree` is the **default layout** a run gets when monitor
`io.directories` values are left alone. The run-owned directories (`log`, `output`)
are containment-validated, so they cannot be pointed outside the run tree without an
explicit override. Isolation between *separate runs* still rests on distinct run ids
rather than on an enforced storage model, which is why this is described as a layout
rather than a guarantee.

| Model | Meaning |
|---|---|
| `isolated_run_tree` | Each run owns a complete tree; nothing shared. **Default, not enforced.** |
| `shared_root_namespaced` | Runs share a parent root, namespaced by run id |
| `shared_pooled` | Shared storage with no per-run namespacing |

@note **Run-owned directories are containment-validated.** `io.directories.log` and
`io.directories.output` must resolve inside the run directory. Absolute paths, values
starting with `~`, and parent traversal are rejected by `picurv validate` and therefore by `run` and
`submit`. This closes a real data-loss path: on a fresh solve the C runtime calls
`PetscRMTree` on its configured log directory, so an escaping value meant the solver
recursively deleted whatever was there.

An **absolute** path outside the run tree can still be used deliberately by setting
`io.directories.allow_unsafe_paths: true`, which downgrades that one rejection to a
prominent warning naming the consequence. It is never the right default, and it waives
nothing else: a relative escape, the run root itself, a reserved directory name, a
`log`/`output` overlap, and a value carrying whitespace, quotes or `#` stay fatal with
the override set. A relative escape is refused outright because it lands among sibling
runs and study members; if an external location is genuinely intended, name it
absolutely. A value starting with `~` is refused outright rather than
expanded: nothing expands it anywhere in the pipeline, so it names a literal `~`
directory inside the run tree. An earlier version of the physical containment check
did expand it and classified the result as an authorizable external location, while
every layer that consumed the value treated it literally - the two disagreed about
which directory was about to be deleted. Refusing the value ends that disagreement;
expanding it in one layer would only have moved it.

Containment is checked at six points:

1. configuration validation (`validate`, `run`, `submit`);
2. the reserved-flag guard that stops raw PETSc passthrough from setting
   `-log_dir`/`-output_dir`/`-restart_dir`;
3. the control-file emitter itself;
4. the dry-run planner, so a plan reports what a real run would refuse;
5. a submission preflight that re-reads the staged `.control` file rather than
   trusting that staging validated it; and
6. independently of all of those, the C runtime's own guard immediately before the
   delete.

The last one matters because the first five live in `picurv`, and the solver can be
launched directly on a hand-written control file. `ClassifyLogDirectory`
(`src/setup.c`) re-derives the decision from the options actually in force, resolves
the path physically — including a directory that does not exist yet, through its
longest existing prefix — and returns a typed verdict. `LogDirectoryIsSafeToWipe`
then applies the waiver to exactly one of those verdicts. The classifier never sees
the authorization, so no ordering mistake can let it short-circuit a physical check.
`-allow_unsafe_log_dir` is reserved, so raw PETSc passthrough cannot forge it.

@warning **Residual risk.** `io.directories.restart` is deliberately not
containment-checked, because it is a read path that may legitimately point at another
run's checkpoints. And an authorization is exactly what it says: with
`allow_unsafe_paths: true` and an absolute external `log`, the runtime *will*
recursively delete that directory on every fresh solve. The guard establishes that the
deletion was asked for, not that it is wise.

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
- **@subpage 52_Run_Lifecycle_Guide** — what the run topology means operationally
- **@subpage 61_Storage_Management_Guide** — archiving and restoring runs
