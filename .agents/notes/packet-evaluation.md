# Packet-first vs baseline: routing evaluation

Question: does `make review-packet` reduce the exploration an agent does before it
starts working, and by how much?

## Why this cannot be run from inside a working session

The measurement is contaminated by knowledge of the answer. A session that has
already read the registries knows where everything is, so its "baseline" arm is not
a baseline. Each arm has to be a session that starts cold, on the same commit, with
no memory of the other arm. Run them separately, ideally on different days or with
memory disabled.

## Arms

Same task, same commit, two independent sessions.

- **A — baseline.** Prompt states only the task.
- **B — packet-first.** Prompt states the task and adds one line:
  "Start by running the applicable `make review-packet` mode before opening files."

Do not hint at the answer in either prompt. Do not run B first: if you must order
them, run A first, since B benefits less from prior knowledge than A loses.

## Tasks

Start with two. Add the rest only if the first two disagree.

1. **Add a value to an existing selector family.** Cheapest, highest frequency,
   exercises `CAPABILITY=`. Suggested: a new `turbulence.wall_function` value.
2. **Diagnose a rank-dependent restart discrepancy.** Exercises `SUBSYSTEM=` and
   `SURFACE=`, and the failure-signature routing in `picurv-solver-debugging`.

Later, if wanted: add a body force; refactor a spatial kernel; add a CLI feature;
generalise a shared operation. The last is the one no skill currently owns.

## Metrics

Primary, and the only one that needs no judgement:

    python3 <script> <session.jsonl> --repo .

The prototype that computed this (`sandbox/prototypes/session_files.py`) was scratch
and has been deleted; rewrite it before running this evaluation. It reported total
tool calls and every repository file the session touched, ranked by how often.
Transcripts are at `~/.claude/projects/-root-projects-PICurv/<session-id>.jsonl`.
The script reports total tool calls and every repository file the session touched,
ranked by how often. Compare:

- **repository files touched** - the headline number;
- **irrelevant files touched** - files not in the correct answer's change set, judged
  once by you after both arms finish, against the same list for both;
- **relevant files missed** - files the change should have touched and neither arm
  opened;
- **tool calls** - a proxy for wall-clock and token cost.

Secondary, recorded by hand:

- did the arm choose the narrowest correct test targets, or run the full suite blindly;
- did it find the existing owner, or write a parallel helper;
- did it surface the applicable contract, freshness state, and documentation obligation.

`/cost` gives per-session token totals if you want them, but files touched is the
metric the packet was built to move and it needs no normalisation.

## Reading the result

The packet is worth its maintenance if arm B touches materially fewer files without
missing relevant ones. If B touches fewer files *and* misses relevant ones, the packet
is over-narrowing and the routing is wrong, not the idea.

A null result is a real result: it would say the deterministic modes are not paying
for themselves yet, most likely because registry coverage is too thin - 24 of 38
`src/*.c` files are currently declared - and the next move is declaring paths rather
than building more tooling.

## Housekeeping

This note now lives in `.agents/notes/` rather than sandbox scratch. Keep it here
while the question is open; if the evaluation runs, promote only the resulting
numbers into `tests/tooling/measurement_records.json` if they justify a status claim
— this methodology writeup does not itself belong in the certified doc site.
