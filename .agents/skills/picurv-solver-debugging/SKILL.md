---
name: picurv-solver-debugging
description: Diagnose PICurv convergence, PETSc or MPI, boundary, restart, and numerical anomalies. Use for evidence-first investigation or requested fixes; not new capability work.
---

# PICurv solver debugging

Establish a reproducible symptom and discriminate among causes before changing solver
settings or code. If the user asked only for diagnosis, stop at an evidence-backed cause
and proposed next step; implement a fix only when requested.

## Capture the actual run

Record the current commit, dirty paths, exact run directory and segment, rank count, binary
hash/provenance, PETSc/MPI versions, launcher and thread settings, rank decomposition, and
the smallest inputs needed to reproduce the symptom.
Inspect the copied run configuration and generated control files, not only the source YAML
that may have produced them.

Start with `docs/pages/67_Troubleshooting.md`, the current `picurv summarize --help`, and
`docs/generated/cli_reference.json`. Prefer read-only JSON summaries for initial triage;
confirm output paths before using plotting commands that may write files. If the user has
not supplied a run or artifacts, state what cannot yet be diagnosed and return a precise
collection plus bounded-experiment plan instead of inventing a cause.

## Route evidence, do not copy lore

Apply the documentation-as-index rule in `AGENTS.md`. Use the relevant registry,
generated inventory, guide, freshness record, and review packet to identify the
smallest plausible code, test, and runtime-evidence set, then verify the live path.
If routing is empty, suspect, never attested, report-only, or names missing symbols,
fall back immediately to targeted code search and report the indexing gap; do not
substitute additional prose for implementation evidence.

- Use `src/guide.md` to trace the live run path from orchestration into momentum, pressure,
  boundary, restart, and monitoring code.
- Use `tests/tooling/capability_families.json` for selector reachability and declared
  evidence, `tests/tooling/subsystem_records.json` for subsystem maturity and its canonical
  observability/troubleshooting pages, and `tests/tooling/freshness_manifest.json` for
  staleness signals.
- Use `docs/pages/09_Monitor_Reference.md` for current monitor keys and log grammar,
  `docs/pages/66_Evidence_Matrix.md` for what has actually been declared as evidence, and
  the solver-specific page routed by the subsystem record for algorithmic invariants.
- Use `tests/tooling/contract_registry.json` to determine whether an applicable invariant is
  enforced, tracked, or planned. Do not call a report-only record enforced, and remember
  that a schema check does not prove numerical restart equivalence.

Do not copy dated iteration counts, timing fractions, convergence floors, or known-defect
measurements into this skill or a new side record. Locate the current canonical prose and
remeasure the claim at the checked-out `HEAD`. Keep selector status, subsystem maturity,
declared evidence, freshness, and an observed run result as distinct facts.

Start the implementation trace with `make review-packet SUBSYSTEM=<subsystem-id>` when a
registered subsystem owns the symptom, and add `SURFACE=<surface-id>` for the watched
implementation or `CONTRACT=<contract-id>` for the suspected invariant. A current optional
xref section can bound caller and reference inspection, including one intermediate hop.
Treat those Doxygen edges only as supporting source references: function pointers,
callbacks, PETSc dispatch, macros, and runtime-selected routes still require a live code
trace. Missing xrefs do not block diagnosis; run `make docs-xref` only when the extra
caller index is worth its build cost.

## Measure and discriminate

Classify the failure quantitatively using the observability that applies:

- continuity or maximum divergence;
- momentum updates/residuals and accepted/rejected trials, or SNES/KSP reasons;
- Poisson tracked and true residuals, iteration counts, and multigrid behavior;
- boundary handler resolution and state changes;
- restart source/step, loaded fields, first-step behavior, and later trajectory; or
- single-rank versus multi-rank behavior with otherwise identical inputs.

For restart anomalies, compare three aligned states: the source checkpoint, the restarted
branch-start state before advancing, and the uninterrupted state at the same step. Inspect
both checkpointed fields and stateful controllers or histories that influence the next step
but may not be serialized. Passing tests for MPI, periodicity, and restart separately do not
establish their combined behavior.

Create the smallest reproduction that preserves the defining condition and a negative
control that removes one suspected condition. Maintain a short hypothesis table with the
predicted signature, supporting/refuting evidence, and cheapest discriminating test.
Change one variable per experiment.

Do not respond to uncertainty with an unbounded parameter sweep, a larger case, or more
ranks. Bound experiments in advance and stop when the evidence separates the leading
hypotheses or when the remaining uncertainty can be stated precisely.

## Coordinate independent evidence lanes

When permitted and useful, delegate independent read-only lanes: one agent traces
code/config resolution, one extracts quantitative log/run evidence, and one challenges the
leading hypothesis or maps the narrowest tests. Give each the same commit and reproduction
identity. The lead owns causal synthesis, edits, and the final claim. Serialize expensive
solver runs and all runs sharing output locations.

## Verify a requested fix

Read `tests/guide.md` and select the narrowest relevant unit target, then exercise the exact
reproducer and its negative control. Add a real runtime smoke test when the defect depends
on ingress, lifecycle, boundary ordering, pressure coupling, or restart. Repeat at another
rank count when decomposition could matter.

Compare before/after metrics under identical inputs. A changed outcome is not yet a causal
explanation; show that the predicted signature changed and that nearby behavior did not
regress.

Before adding a fix-specific helper or call-site workaround, apply the reuse gate in
`AGENTS.md`. Repair or generalize the owning shared path when the defect is common to its
contract; retain a specialized path only when differing semantics or performance are
demonstrated.

If the fix changes documented behavior or invalidates an empirical characterization, use
the applicable review packet before editing the canonical page, follow the freshness and
documentation audits in `AGENTS.md`, and update that page rather than creating a parallel
diagnostic ledger. Report any review/audit tooling failure explicitly.

## Handoff

For an implemented fix, run `make review-packet CHANGED=working-tree` before handoff and
investigate any unrouted production path using its nearest guide and targeted search. The
changed-set result is advisory and does not turn declared evidence into an observed fix.

Provide the reproduction identity, measured symptom, leading cause and causal evidence,
excluded alternatives, before/after metrics when a fix was made, tests/rank matrix run,
reuse or new-helper decisions, and residual risks. Label inference as inference. Never turn
a passing schema, declared evidence source, or one numerical case into a broader
verification claim.
