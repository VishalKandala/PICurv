---
name: picurv-test-design
description: Choose, size and write PICurv tests. Use whenever a test is written, changed, run to verify a fix or feature, or proposed as a regression guard - unit, smoke, MPI, restart, convergence or validation - and whenever deciding which suites to run.
---

# PICurv test design

Run the checks that answer the question in front of you, and write tests that cover the
class of situations a defect lives in, not the one that exposed it. Compose with the
workflow skill that owns the change.

## Run what answers the question

- While iterating, run the narrow targets the change touches: `tests/guide.md` and
  `make review-packet CHANGED=working-tree` name them.
- The pre-push hook runs the full certification on every push to `main`; leave the full
  suite to it unless the run itself is the question (another index width, another MPI
  layout, the cluster). Branch pushes are not gated.
- A run outside the gate - a sandbox case, a scratch script, a cluster check - is evidence
  too. Record what ran and what it showed in the commit message or the dated record.

## Decide whether a check becomes a permanent test

A one-off feature check stays bespoke unless it guards a contract no gated case exercises.
Promote it when it does, at the cheapest tier that covers it (unit, then smoke, then MPI).
When a class of error reached a run past the gate, extend the gate for that class, not only
for the instance.

## Cover the class, not the instance

When a defect is found, write down its class boundary: which inputs, layouts, faces,
topologies or regimes trigger the same mechanism, and which do not. Then test:

- at least two distinct representatives inside the boundary, and one control outside it;
- that each new test fails on the code before the fix, passes after it, and fails again
  if the fix alone is reverted.

For each test, name the production change that would make it fail. Mutate the code
mentally - a wrong constant, a swapped branch, a missing boundary or zero case - and make
sure some test fails for each.

## Use relations when there is no exact answer

Most solver results have no closed-form value. Assert relations instead:

- **Decomposition.** One rank and a layout whose seam crosses the stencil must agree; two
  layouts must agree with each other.
- **Restart.** A run split at its midpoint and resumed must match the uninterrupted run, to
  the documented restart floor; compare particles by ID, not by order.
- **Orientation.** A boundary condition placed on each of the six faces gives the mirrored
  or rotated result.
- **Symmetry and shifts.** Periodic shifts, reflections, free-slip half-domains.
- **Order of accuracy.** For a discretisation change, a manufactured or analytic solution on
  a skewed or curved grid at three or more resolutions, asserting the observed slope;
  round-off and iteration error at least 100 times below the discretisation error.

## New boundary condition or wall treatment

Cover: a kernel unit test against an analytic value; each face it can sit on; an edge or
corner cell shared with another boundary; a rank split that cuts the face; a midpoint
restart; interaction with a periodic axis; and, where the scheme allows, an order test with
the condition on one face.

## Tests are part of the trail

Name what each test verifies in its docstring, and cite it from the capability or
subsystem evidence it backs (`tests/tooling/capability_families.json`). A test weakened or
deleted to make a change pass is a defect in the change.
