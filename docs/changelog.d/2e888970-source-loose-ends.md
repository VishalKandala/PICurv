
- Initial conditions: a source's `pressure` expression [Pa] sets the initial `P` in a solved case.
- Fixed: `ber` and `bei` silently lost accuracy at large arguments (1e-9 relative at 60, 4e-4 at 100, meaningless by 150). They use an asymptotic expansion from `|x| = 30`.
- Fixed: changing a frame of a `field_slice` recording did not invalidate a precomputed profile.
- Fixed: `GridSpace.z` averaged one corner of the lower k-face twice. Nothing reads the field, so no result changes.
- Examples: `lid_driven_cavity`, a square cavity at `Re = 100` driven by `handler: moving`. Its 64 x 64 grid has not been run.
- Tests: a test failing on some ranks aborts the job rather than hanging it, and `unit-mpi` passes at 4 ranks with two two-rank tests skipped.
