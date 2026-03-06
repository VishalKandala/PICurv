# Smoke Case Fixtures

The smoke harness currently materializes tiny runtime cases dynamically via:

- `picurv init flat_channel`
- `picurv init brownian_motion`

and then rewrites them in-place for short end-to-end runs and restart branch checks.
