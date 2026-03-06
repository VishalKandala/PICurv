# Smoke Post Fixtures

Post profiles are rewritten in the smoke runner to guarantee deterministic output checks on tiny runtime paths.

## Assertions Enabled By Post Rewrites

- Eulerian VTS output on tiny solve runs,
- particle VTP output on particle-enabled runs,
- MSD CSV statistics output on analytical Brownian smoke sequences.

## Why This Layer Exists

Smoke testing must verify that postprocessing contracts still hold when solver output structure, naming, or run orchestration changes. Deterministic post profiles make these regressions visible quickly.

## Maintenance Notes

- If post schema keys or output naming conventions change, update rewrites and assertions in the same commit.
- Keep smoke post tasks focused on contract coverage, not expensive analysis.
- Use larger analytical studies outside smoke for in-depth physics validation.

## Contributor Reminder

Keep this fixture documentation aligned with smoke assertions whenever post schema or naming changes.
