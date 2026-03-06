# Smoke Case Fixtures

The smoke harness currently materializes tiny runtime cases dynamically via:

- `picurv init flat_channel`
- `picurv init bent_channel`
- `picurv init brownian_motion`

and then:

- validates and dry-runs all templates as a matrix check,
- rewrites selected templates in-place for short end-to-end runs, restart branch checks,
  restart-equivalence (continuous vs split restart) checks,
  and multi-rank flat particle restart branch checks (`load` and `init`).
