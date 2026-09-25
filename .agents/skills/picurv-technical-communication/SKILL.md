---
name: picurv-technical-communication
description: Explain and report substantive PICurv work consistently for technically trained users and developers. Apply to case design or execution, output and plot interpretation, scientific comparison, mathematics, code, debugging findings, reviews, and feature plans; compose with the applicable PICurv workflow skill rather than replacing it.
---

# PICurv technical communication

Make PICurv discussions accessible to readers who may be expert in fluid mechanics,
numerical analysis, scientific computing, or software development without assuming that
they are expert in all four. This skill controls communication and handoff quality. It
does not replace the workflow that constructs or runs a case, analyzes results, diagnoses
a defect, reviews code, or designs a feature.

## Compose with the current task

Apply every applicable PICurv workflow skill first. In particular, retain the investigation,
implementation, and verification requirements of `picurv-solver-debugging`,
`picurv-capability-change`, and `picurv-spatial-kernel-change`. Any future case, analysis,
or validation skill likewise owns its operational workflow.

Use `AGENTS.md` and the selected workflow skill to locate the smallest relevant registry,
generated inventory, guide, review packet, source path, test, and runtime evidence. Treat
documentation as an index and follow the repository trust hierarchy. Current runtime and
code behavior remain authoritative.

Do not copy changing PICurv facts into this skill. Resolve selectors, defaults, equations,
field layouts, paths, and execution ordering from their current owners. The workflow skill
controls what to inspect or do; this skill controls how to explain the resulting knowledge.

## Calibrate without talking down

Infer the reader's perspective from the request and conversation. Assume technical
literacy and undergraduate mathematics, but do not assume familiarity with PICurv,
PETSc, MPI, curvilinear coordinates, a particular discretization, turbulence model, or
solver family.

- For a fluid-mechanics reader, connect modeled physics to discretization and observable
  behavior.
- For a numerical-methods reader, emphasize discrete operators, approximation assumptions,
  stability, conservation, conditioning, and convergence.
- For a scientific-computing developer, emphasize control flow, data layout, ownership,
  synchronization, lifecycle, and cost.
- For a mixed or unknown audience, begin with the shared physical or mathematical idea and
  introduce specialist details only where they affect the conclusion.

Define unfamiliar terms locally. Preserve real equations, mechanisms, and caveats. An
analogy may establish intuition, but state where it stops matching the implementation or
mathematics.

## Lead with meaning and connect only useful layers

Start with the result, practical consequence, or decision. Then provide enough mechanism
and evidence for the reader to assess it independently. When relevant, connect:

`physics -> mathematical model -> discrete method -> PICurv implementation -> configuration and execution -> logs, fields, statistics, or plots -> scientific interpretation`

Use only the layers needed for the question. Do not force every response through the full
chain, repeat the conclusion in several sections, or supply background the reader already
demonstrates. Prefer one representative live path or example over an exhaustive survey.

## Keep claim types distinct

Use these concepts consistently, with labels when ambiguity would matter:

- **Observed:** directly present in inspected output or an executed experiment.
- **Confirmed:** established from current code or an enforced contract.
- **Documented:** stated in current documentation but not independently established.
- **Inferred:** supported by evidence but not directly demonstrated.
- **Proposed:** a suggested interpretation, experiment, design, or next step.
- **Not verified:** applicable or plausible, but the required check was unavailable or skipped.

Separate what was found, what it means, and what should happen next. State the comparison
basis, assumptions, and important alternatives when they can change the conclusion. Never
turn declared evidence, a schema check, or one numerical case into a broader correctness
claim.

## Route conditional detail

Read only the references needed for the current discussion:

- Read [references/science-and-mathematics.md](references/science-and-mathematics.md) for
  equations, models, discretization, units, numerical properties, or physical meaning.
- Read [references/code-and-runtime.md](references/code-and-runtime.md) for source paths,
  PETSc/MPI behavior, ownership, synchronization, lifecycle, or performance.
- Read [references/evidence-and-results.md](references/evidence-and-results.md) for runs,
  plots, statistics, comparisons, validation, literature, plans, or technical handoffs.

Read more than one only when the question genuinely crosses those concerns.

## Scale the response

For a focused question, answer directly in a few paragraphs. For a substantive finding,
plan, or handoff, preserve the following semantic fields without necessarily turning each
one into a heading:

- bottom line or objective;
- claim status and evidence location;
- mechanism or reasoning;
- assumptions, constraints, and scope checked;
- recommended next action; and
- unresolved uncertainty or decision.

For agent-to-agent handoffs, transmit findings, evidence locations, constraints, open
questions, and the next useful action. Omit the chronological narrative of the search
unless it changes how the evidence should be interpreted.
