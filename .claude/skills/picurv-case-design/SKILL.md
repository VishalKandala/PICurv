---
name: picurv-case-design
description: Design PICurv solver runs before staging them. Use whenever a case, pilot, run, extension, sweep, study, diagnostic matrix or campaign is being set up, changed or submitted, locally or on a cluster - before `picurv init`, `run`, `sweep`, `submit` or any job script.
---

# PICurv case design

Decide what each run is for, what it shares with other runs, and what must happen before
it, then stage it once. Every rule here exists because a run was staged without it.
Compose with `research-data` for filing, `picurv-test-design` when a run is a check of the
code, and `picurv-technical-communication` for the write-up.

## 1. Look, then lay out approaches

Read before asking: the campaign README and dated records, the measurement registry
(`tests/tooling/measurement_records.json`), existing restarts and archives, and the case's
own `grid.info`. Never ask what the repository or records already answer.

Then put two or three approaches to the owner in a message of their own - what each runs,
what it can and cannot show, its cost - with a recommendation. Stage nothing until one is
chosen. A blanket compute approval covers the budget, not the design.

## 2. State the question and its anchor

For each run, write one sentence: the question it answers and the decision its result will
drive. Name the reference it is judged against (experiment preferred, then DNS/LES in the
literature, then an analytic result) and the criteria, fixed now, before any data. Name what
the run does not test, so a later reader cannot over-read it. A new campaign starts with the
literature round (`research-data`, "Starting a campaign"); the owner picks the conditions.

## 3. Find shared states and order the phases

- **One state per physics.** Runs that differ only downstream of a field - scalar
  orientation, particles, post-processing, a statistics window - restart from one shared
  run of that field. Two scalar experiments on the same flow are one flow and two restarts.
- **Phase order.** Flow spin-up, then a stationarity test fixed in advance, then particles,
  scalars or statistics from a steady restart (`restart_mode: init`,
  `--statistics-state reset`). Particles injected during a transient give no usable data.
- **Preconditions.** A variant may be staged only when its precursor's stationarity
  record exists. Extend in blocks (for flow: two flow-throughs, then one at a time) with
  checkpoints sparse until steady.

## 4. Fix the case specification before any run

Write down, per run: initial condition, boundary conditions, grid and rank layout, time
step and CFL, the discard window (transient) and how its end is detected, the averaging
window, the stopping rule, outputs (derived data first), and checkpoint cadence. Every
reported statistic carries an error bar from correlated-sample statistics; a grid
difference smaller than the error bars supports no convergence claim.

## 5. Check the grid and the resolution

- Read `grid.info`: expansion ratio (aim at or below about 1.2), non-orthogonality, aspect
  ratio, Jacobian, nominal wall units.
- Compare y+, Δx+, Δz+ and cells across the key features with the reference's grid and with
  the wall treatment's valid range (page 07).
- For LES, estimate the subgrid energy fraction (Pope's M) and plan at least three
  resolutions before calling a result grid-independent.
- Check the multigrid rank-layout rule (page 25) and that every rank owns at least three
  nodes normal to a modelled wall.

## 6. Estimate cost from a pilot that reaches the physics

Pilot at the real grid and layout; take the per-step cost and extrapolate (AGENTS.md: over
45 minutes goes to the cluster). A pilot must reach the time or state where the physics of
interest - and its failure modes - appear, or say plainly what it cannot show. Scale between
Reynolds numbers only with a stated law. Pack related runs into one allocation: stage each
with `picurv run --no-submit`, then run their `solver.sbatch` in sequence, with walltime for
the sum. Large campaigns need the owner's approval; approval does not carry across campaigns.

## 7. Design diagnostic matrices

Change one setting per run, from the same start. Fix the success measure before the runs.
Name every confound a run cannot separate (a coarser cross-section also moves y+). Choose
the smallest set that discriminates the candidate causes, and include the control.

## 8. Record the decision, review, then stage

Write the plan and its decisions as a dated record (`research-data`). For runs costing more
than a few node-hours, have it reviewed by a reader given the plan, not the conclusion.
Stage the whole package once; do not revise staged runs piecemeal.

## 9. After the run

Compare against the criteria fixed in step 2 before interpreting. Then complete the
documentation trail in `AGENTS.md`, "Documentation trail": products and records
(`research-data`), the research-log study entry, and any measurement record.
