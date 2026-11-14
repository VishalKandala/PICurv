# PICurv — working agreement for AI coding sessions

Applies to Claude Code, Codex, and any other agent working in this repo.
`AGENTS.md` is the canonical shared instruction file. Tool-specific instruction
files may import it, but must not duplicate it.

## Never create new top-level directories

The repository root is a fixed, curated set. Do not add directories to it —
not for scratch work, not for results, not for notes, not "temporarily".
If you believe a new top-level directory is genuinely required, stop and ask
first; do not create it and mention it afterwards.

## Where work goes

| Kind of file | Location |
|---|---|
| Transient files with no lasting value (intermediate output, one-off scripts, parsing buffers) | the session scratchpad **outside** the repo |
| In-repo scratch: prototypes, exploratory notes, trial configs | `sandbox/` (see `sandbox/guide.md`) |
| Generated transient files that must sit in-repo | `sandbox/tmp/` |
| Solver run output | `runs/` — already gitignored |
| Sweep/study output | `studies/<study_id>/` — created by `picurv sweep`/`submit`, already gitignored |
| Solver and tooling logs | `<repo>/logs/` — already gitignored |

Prefer the session scratchpad. Only put scratch in `sandbox/` when it needs to
outlive the session or be shared with a human.

`runs/`, `studies/`, and `<repo>/logs/` are the product's own output paths, owned by
`picurv_cli`. Let the tooling create them. Do not hand-create them, write
unrelated files into them, or invent parallel alternatives.

Those rows describe the project's own output. Your own trial runs, pilots, and
sweeps are scratch and do not belong there. Initialize a throwaway case inside
`sandbox/` and work inside it:

```bash
./bin/picurv init flat_channel --dest sandbox/<name>
```

`picurv init` places the case relative to the current working directory, and every
artifact a run produces is written under the case it belongs to, so that case's
own `runs/` and `studies/` sit inside `sandbox/<name>/`. Deleting that one
directory removes the experiment whole, which is the point: delete it once it has
answered its question.

Do not instead redirect a repository-root run's directories into `sandbox/`. The
run-directory containment guard rejects paths that escape their run root, and it
is right to.

Never delete, move, or reorganize the user's own runs; `picurv storage` owns their
lifecycle.

## Sandbox hygiene

`sandbox/` is scratch, not a staging area for deliverables. Its contents are
gitignored except `guide.md`, so nothing there reaches a commit by accident —
which also means nothing there is backed up. Anything that matters gets
promoted out per the promotion rules in `sandbox/guide.md`, or it is lost.

Clean up after yourself: delete scratch you created once it has served its
purpose, rather than leaving it for the next session to inherit.

## Long jobs go to an executor, not to your session

Estimate wall-clock cost before starting any solve, sweep, or post-processing
run. The cheapest estimate is a pilot: run a handful of timesteps at the real
grid and rank layout, take the per-step time, and extrapolate to the requested
step count. This is the same warmup-and-extrapolate estimator the generated
scheduler jobs already use for their walltime guard. Do not guess from grid size
alone, and do not start a long run to see how it goes.

**Local limit: 8 ranks and 1 hour.** A local agent runs locally only what fits 8
ranks and 1 hour by a measured estimate: builds, unit and smoke tests,
development checks, and evidence runs for the measurement registry. Run it in the
background and report. The local PETSc is a debug build, so a local run is
evidence of behaviour, never of cost or speed.

This limit is narya's, on the development machine. An executor follows its own
cluster's handbook for what may run on its login node; everything else there is a
job.

**Beyond the limit:** do not start it, detach it, or chop it into pieces that fit.
The work goes to a cluster executor under the agent protocol on Drive
(`Research/PICurv/agents/protocol.md`): the dispatcher (narya) sends a `RUN` naming
the specification or its record, the build, the stop conditions and the reports
wanted. A co-developer (Codex) hands such work to narya and never instructs an
executor. Within the limit, narya may still send work to an executor when it needs
an optimized build, a cost measurement, data that lives on a cluster, or narya's
time for something else.

The executor stages and submits with picurv (`picurv run|sweep ... --cluster
<cluster.yml> --no-submit`, then `picurv submit --run-dir|--study-dir`), using its
cluster's handbook (`agents/clusters/<name>.md`) for accounts, partitions and the
environment. Do not invent site policy.

The `RUN`, or the record it points to, carries what the tooling does not
generate: the analysis scripts that turn raw solver output into the answer, and
exactly which files to bring back. Ask for derived data, not raw state: metrics,
summaries, spectra, profiles, solver and scheduler text output, post-processing
results. Request checkpoint or field binaries only when the set is small and
nothing derived can answer the question. Returned data lands in `sandbox/`, and
the analysis continues there under the same hygiene and cleanup rules as any
other scratch.

## Re-run what the change touched, not everything

While iterating on a failure, the question is which tests the fix could have
changed the answer for. Re-running only the test that failed is not enough: a fix
can break a test that was passing. Re-running the whole suite after every fix is
the other extreme, and on a change confined to one module most of it is answering
a question nobody asked.

Re-run the tests that cover the modules the fix touched. `tests/guide.md` maps
subsystems to their narrowest targets, and `make review-packet CHANGED=working-tree`
names the test targets declared for the paths you actually changed. When the packet
reports a changed path as unrouted, fall back to the nearest directory `guide.md`
and choose by responsibility; do not treat an unrouted path as "nothing to re-run".

The pre-push hook runs the full certification on every push to `main`, so the full suite is
not re-run by hand before committing. Run it yourself only when its result is the question -
another index width, another MPI layout, the cluster - or before a push that is not to
`main`; branch pushes are not gated. `picurv-test-design` covers which checks answer which
question and when a check becomes a permanent test.

Do not pre-run the merge gates. `make certify-docs` and `certify-docs-fast` are
commit-scoped by construction: they refuse a dirty tree, because a certificate is
claimed for one commit. The pre-push hook and CI run them on the way to `main`
anyway. Running them mid-iteration returns an error, not an answer.

A test that was passing before your change is not evidence that it still passes.
State which targets you re-ran and which you did not, so the difference is visible
rather than assumed.

## Plan before acting

Before any code change, case, run or campaign step, load the workflow skill that owns it -
`picurv-capability-change`, `picurv-spatial-kernel-change`, `picurv-solver-debugging`,
`picurv-case-design`, `picurv-test-design` - as the first step, not after reading the code.
Skill loads made for an earlier change do not carry over to a new one. Read what the
repository and records already answer, write the plan, and when there is a real choice put
the options to the owner with a recommendation. An approval of compute or of a goal is not
approval of a design.

## Technical communication

For substantive PICurv explanations, plans, findings, observations, recommendations,
and user-facing or agent-to-agent technical handoffs, apply
`picurv-technical-communication` alongside every applicable workflow skill. This includes
case design and execution, output or plot interpretation, scientific comparison,
mathematics, code behavior, debugging, reviews, and feature design. Trivial file lookups
or commands do not require it.

Workflow skills control investigation, actions, and verification. The communication skill
controls explanation depth, organization, terminology, and claim labeling. Use the
repository documentation as an index and verify behavioral claims according to the trust
hierarchy below. Do not let presentation guidance weaken a workflow skill's evidence or
handoff requirements.

## Code-change workflow

The standard path for a code change (owner, 2026-10-19):

1. Work on a separate branch. A feature branch (any branch but `main`) may be pushed without
   asking, by narya or Codex, so that an executor can build it.
2. Test locally as far as the local limit allows.
3. Merge into `main` with the commit stating that the change is still `EXPERIMENTAL`.
   A push to `main` needs the owner's go or a grant recorded in the thread (agent protocol §7).
4. Executor verification: narya names the commit in a `BUILD` (agent protocol §6) and an
   executor tests it on its cluster. Iterate (branch, test, merge, verify) until done.
5. Update documentation, guides and registries.
6. Delete the branch.

Follow the existing commenting and naming conventions, and the reuse rules below.

## Committing

Do not commit scratch, temp files, run output, or generated artifacts. Before
committing, check `git status` and confirm every path is one you intend to
track. Commit only when explicitly asked. A feature branch may be pushed without
asking; a push to `main` needs the owner's go or a grant recorded in the thread
(agent protocol §7).

After a push, read the workflow results for that commit (`gh run list --commit
<sha>`) before reporting the push as done. The local Doxygen is older than CI's,
so a green local gate does not predict CI.

## Reuse before writing code

The default is no new production code. Before adding a function, helper, loop,
wrapper, state field, or special-case path, search for the operation that already
owns the behavior and inspect its callers. Search by responsibility as well as by
the name you expect. In the plan or handoff, account for each new abstraction with:
the existing candidates inspected, why direct reuse is or is not correct, and why
generalizing the current owner is or is not preferable.

A new capability is not by itself a reason to create a new production code file.
Extend the existing owning module by default. A new source, header, or production
module requires a distinct and durable responsibility, or a real build, lifecycle,
dependency, or measured performance boundary. Before creating one, identify the
plausible existing homes and explain why none is coherent. Unless the user explicitly
authorized a new module, ask before creating the file. Do not create one-feature files
or generic `utils`/`helpers` dumping grounds.

In particular, do not reimplement physical/owned range selection, field lookup,
I/O, allocation/lifecycle, PETSc scatter or assembly, periodic-axis detection,
periodic or boundary repair, metric access, or logging when the repository already
has the operation. This list is representative, not exhaustive. Resolve current API
names and contracts from code; similar-looking operations are not interchangeable
when layout, ownership, vector authority, or call ordering differs.

When an existing operation is almost sufficient, prefer a coherent generalization at
its owning module and migrate applicable callers instead of adding a parallel helper
or compatibility wrapper. A new shared helper is justified only when no existing
contract can satisfy the required semantics and performance. Put it at the narrowest
stable owner, give it explicit inputs and side effects, and make its contract reusable
by other current or foreseeable callers when that does not weaken semantics or add
cost. Separate reusable mechanics from feature-specific mathematics; do not bake a
selector, solver, or physics term into a general operation that does not require it.
Do not create speculative abstractions for hypothetical reuse.

Identify every current bespoke call site that the new or generalized operation could
replace. Unless the user already authorized consolidation, show them the candidates,
semantic and performance differences, migration risk, and verification needed, then
ask whether those callers should be migrated. Do not modify them without approval. If
consolidation remains out of scope, record the remaining bespoke paths and do not
pretend the duplication was eliminated.

`UserCtx` has a back-pointer to `SimCtx`; `SimCtx` owns the multigrid hierarchy of
per-level, per-block `UserCtx` objects. Reuse that existing state instead of copying
it, but pass the narrowest context that identifies the required ownership. Do not
pass `SimCtx` merely because it can reach everything when the operation belongs to
one block or level.

Reuse must not hide or add costly work. Do not introduce extra allocation, copying,
collectives, scatters, synchronization, or hot-loop indirection for architectural
tidiness. Preserve a specialized fast path when measurement or the numerical contract
requires it, and document why the shared path is unsuitable.

## Performance changes

Most PICurv cases are coupled: the flow is solved and particles are carried. Judge a
performance change on a coupled case whenever it touches a shared path, and do not tune
to particle-only runs or to one multigrid level count or rank layout. A particle-only
setup is a development harness; a fix judged only on it is not done. Keep a change whose
results are proven identical even when the gain is small, unless it adds complexity.
Identical means Eulerian fields bit-identical (flag a change that legitimately alters
summation order) and particles identical one by one: the storage order of particles in
the swarm need not be preserved, so compare checkpoints sorted by particle id, not byte
for byte.

## Runnables and modules

Everything a person runs is extensionless (`bin/picurv`, `generators/*.gen`).
Importable Python modules live only in `picurv_cli/`; a generator loads one by path and
lists it in `GENERATOR_SHARED_MODULES`, so the run lock hashes it. Tests stay `.py`.

## Documentation: what to trust, and what you owe

### Use documentation as an index

Start with the applicable registry or generated inventory, directory `guide.md`, and
review packet to reduce the repository to the smallest plausible set of source symbols,
files, contracts, and tests. Then inspect that code and those tests before relying on a
behavioral claim or editing the implementation. Documentation routes; runtime and code
remain authoritative.

Treat an empty review packet, a suspect or never-attested freshness surface, a
report-only contract, or a missing named symbol as an incomplete index. Fall back
immediately to targeted code search and live-path tracing; do not compensate by reading
more prose. Report the routing gap and update its existing documentation owner when the
task's documentation obligations authorize or require that work.

### Trust hierarchy

When sources disagree, believe them in this order. Do not "fix" a lower tier to
match a higher one without checking which is actually wrong.

1. **Runtime and code behaviour** — what the program does.
2. **Enforced contracts** — schemas, manifests, and the registries under
   `tests/tooling/`. These are machine-checked, so they are rarely wrong for long.
3. **Generated documentation** — anything under `docs/generated/`. Never edit these
   by hand; regenerate with `make docs-inventory` or `make docs-topology`.
4. **Directory `guide.md` files** — local routing. Useful, not authoritative.
5. **Explanatory prose** — long-form pages. Most likely to have drifted.
6. **Proposals and history** — `57_Future_Architecture_Specifications`, changelog
   entries. These describe intent or the past, never current behaviour.

A page being well-written is not evidence that it is current. Verify any enumeration,
default, or ordering against tier 1 or 2 before repeating it.

### How pages are written

- A page states the current behaviour. Do not frame it against the previous behaviour
  ("now", "used to", "no longer"); that history belongs in changelog fragments and commit
  messages. An alias may be listed as an accepted spelling without calling it old.
- Planned work is described only on page 57 (`57_Future_Architecture_Specifications`) and
  the pages it indexes (59, 60, 17). Everywhere else (reference and how-to pages, examples,
  `master_template`, `<repo>/config/`, READMEs) it is absent, including "planned, see page 57"
  pointers and keys that exist only for a planned feature. Stating that a capability is
  absent is allowed; source stubs and conductor error messages are fine.
- The foundational physics pages (19, 22, 23, 26, 28, 34) are updated together with the
  owner, never authored unsupervised.
- Before writing that the reason for a design choice is unknown, or asking for it, search
  `git log` (subject, body, `-S`); commit messages often hold rationale that pages and
  comments lack. Put this check in the brief of every subagent that writes about the code.

### Before changing documentation

Run `make review-packet PAGE=<n>` or `make review-packet CONTRACT=<id>`. It assembles
the declared sources, owning and dependent contracts, capability families, evidence
sources, uncommitted changes, and how to verify — so you read the right few hundred
lines instead of grepping.

### If you add or change a subsystem's status

Statuses are `planned`, `internal`, `experimental`, `supported`, `known-defective`,
`deprecated`, `removed`. Each owes more documentation than the one below it, and the
ladder is cumulative. Update the record in `tests/tooling/subsystem_records.json` and
run `make audit-subsystems`; it will name what the new status owes and refuse a
transition the lifecycle does not allow. A reasoned "not applicable" is a real answer.
Empty prose written to pass the gate is not.

### If you add or change a user-selectable option

The documentation obligation is enforced, not advisory. See
`docs/pages/50_Modular_Selector_Extension_Guide.md` §1.1. In short: declare it in
`tests/tooling/capability_families.json`, scaffold its entry with
`tests/tooling/scaffold_documentation.py`, then `make docs-inventory` and
`make audit-capability`.

### Verification

```bash
make audit-capability      # capability parity, coverage, evidence, value lifecycle
make audit-subsystems      # what each subsystem's claimed status owes
make audit-page-types      # every published page is typed
make audit-freshness       # hard (blocking) and soft (advisory) staleness suspicion
make audit-contracts       # invariant contracts (run after build-docs)
make audit-path-literals   # no hardcoded run paths in prose
make preview-docs          # build, gate, and serve the site
make certify-docs          # everything, for a clean tree
```

`make audit-freshness` reports which pages need review because something they describe
changed. After you have actually compared a page with its sources, record it:

```bash
make attest-freshness ARGS="<surface-id>"
```

Attesting records that the comparison happened. It is not a claim that the page is
correct, and it is never a substitute for making the page correct.

### Documentation trail

Work on the code or in a campaign is done when its trail is complete, so that both the
code reasoning and the scientific reasoning can be found later:

- a changelog fragment for every push to `main`, in a shape CI's Doxygen accepts: one
  top-level bullet with everything nested, plain prose, or several top-level bullets after a
  leading blank line; never de-dent back to column 0 after nesting
  (`docs/changelog.d/README.md`);
- repository documentation, directory `guide.md` files and code comments that describe the
  behaviour as it now is;
- the research log on Overleaf (local clone `../picurv-research-log`): Part I describes the
  method as it now is, with its reasoning and limitations; Part II records the dated study
  or decision, citing record, product, run and archive ids. Pull before editing and push
  after, following `picurv-technical-communication`, `references/research-log.md`;
- dated records and products on the research drive (`research-data` skill);
- a measurement record for evidence that backs a status.

### Two standing rules

**Never claim more than was checked.** "Built from commit X" is not "verified". A
declared evidence source is not a verified result. If automation cannot establish
something, say so rather than implying it did.

For reviews, a statement that no issue was found applies only to the files, symbols,
contracts, tests, ranks, configurations, and runtime paths actually inspected. Mark a
claim **not verified** when its check was skipped or unavailable, it lies outside the
declared scope, or the necessary evidence is absent. Passing a schema, registry, or
documentation audit does not establish runtime or numerical correctness.

Apply the relevant PICurv skill during native agent review. Review-only requests are
read-only unless the user separately authorizes changes. Report confirmed findings with
file and line locations, and keep plausible risks explicitly separate from demonstrated
defects. Do not turn an absence of evidence into an all-clear.

**Never write filler.** Guidance that would read the same on any page is not guidance.
`make audit-docs-expansion` rejects the specific filler that was removed from 51 pages;
the rule is broader than the check. A section you cannot fill with something specific
should not exist.
