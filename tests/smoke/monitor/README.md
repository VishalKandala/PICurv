# Smoke Monitor Fixtures

Monitor profiles are rewritten by the smoke runner so tiny runtime sequences produce deterministic, assertion-friendly output.

## Rewrite Intent

- force output cadence to every step so restart artifacts are guaranteed,
- align particle console cadence with output cadence in particle scenarios,
- retain enough logging detail to diagnose failures without overwhelming CI logs.

## Why Deterministic Monitor Cadence Matters

Smoke assertions rely on specific artifact availability (for example restart checkpoints and postprocessable outputs). If monitor cadence drifts, tests become flaky and failures no longer indicate real runtime regressions.

## Maintenance Guidance

- Keep smoke monitor rewrites minimal and purpose-driven.
- If monitor schema keys change, update rewrites and assertions together.
- Avoid introducing non-deterministic output intervals in smoke-only profiles.

## Contributor Reminder

Keep smoke fixture docs synchronized with run script intent so scenario coverage stays explicit.
