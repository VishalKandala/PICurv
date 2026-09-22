- Post-processing accepts a single checkpoint with `start_step = end_step` and
  `step_interval: 1`, regardless of the monitor output cadence. Zero and negative
  intervals are rejected in Python and C instead of allowing a non-advancing loop.
