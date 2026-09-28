# Evidence, results, and decisions

Use this guidance for run output, monitoring, visualization, comparisons, validation,
literature, plans, and handoffs.

## Identify the evidence before interpreting it

Record the minimum provenance needed to make the result meaningful: case and segment,
configuration, commit and binary provenance, grid, timestep or averaging window, rank
layout, field or metric definition, units or normalization, and relevant postprocessing
method. Include only fields that can change the interpretation.

For a plot, explain:

- what each axis, curve, normalization, and aggregation represents;
- whether the data are instantaneous, transient, averaged, or statistically converged;
- the expected qualitative or quantitative signature;
- what is actually visible; and
- which alternative causes the plot alone cannot distinguish.

Do not equate a visually smooth field, a falling residual, or qualitative resemblance with
physical validity or numerical convergence.

## Make comparisons commensurate

For cross-case or literature comparisons, establish the basis before judging agreement:

- geometry and coordinate definitions;
- governing equations, modeled terms, and boundary conditions;
- Reynolds number and other nondimensional parameters;
- forcing and bulk-quantity definitions;
- grid, timestep, averaging, and statistical uncertainty;
- variable location and normalization; and
- interpolation or digitization introduced by the comparison.

Separate verification (solving the implemented equations correctly), validation (agreement
with physical reality or accepted reference data), and exploratory comparison. State
whether agreement is qualitative, quantitative, or within a stated uncertainty.

When current literature is required, use authoritative or primary sources and preserve
enough citation information to reproduce the comparison. Do not carry a value from a paper
into a recommendation without checking that its definitions and regime match.

## Report findings and plans economically

For a finding or handoff, communicate:

```text
Finding or objective:
Status: observed | confirmed | documented | inferred | proposed | not verified
Evidence and comparison basis:
Interpretation:
Constraints or alternatives:
Next useful action:
```

Use the fields semantically; omit labels that add no clarity. For a plan, include the
acceptance criterion and the cheapest evidence that can discriminate among meaningful
alternatives. For a recommendation, say what evidence would reverse it.

Keep raw chronology, exhaustive command output, and dead-end hypotheses out of the handoff
unless they prevent repeated work or qualify the conclusion.
