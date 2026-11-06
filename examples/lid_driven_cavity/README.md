# Lid-driven cavity

A square cavity whose top wall moves at unit speed, the example of `handler: moving`. At
`Re = 100` the flow settles to one primary vortex with two small corner eddies; Ghia, Ghia and
Shin (1982, J. Comput. Phys. 48) tabulate its centreline velocities, the usual reference.

`grid.gen` reads `grid.cfg`: 64 x 64 cells on a unit square, with 8 periodic cells along `k`, so
the flow is two-dimensional. Multigrid coarsens only `i` and `j`, because a periodic axis keeps at
least the stencil's width on every level. The lid is `+Eta`, moving along `+x`; the other walls
are at rest and the velocity starts from zero.

The velocity is discontinuous at the lid's two top corners, where the moving wall meets a wall at
rest. No grid resolves that, so the error near them falls more slowly than second order, and the
pressure there grows as the grid refines. Centreline comparisons are the better measure.

## Run and compare

```bash
picurv init lid_driven_cavity --dest my_cavity
cd my_cavity
picurv run --solve --post-process -n 4 --case config/case.yml --solver config/solver.yml \
    --monitor config/monitor.yml --post config/post.yml
```

The shipped 3000 steps reach `t = 30`, about 3 viscous times of the cavity. Check steadiness in
`solution_convergence.csv`, then cut `Ucat_nodal` along `x = 0.5` and compare `u` with Table I
of Ghia et al.

## What was checked

On a reduced 32 x 32 grid, run once with `dt = 0.02` to `t = 24`, `u` along `x = 0.5` stayed within
0.003 of Ghia's values, typed from memory of the table, with a minimum of -0.2088 against
-0.2109. This grid has not been compared with the paper, so treat it as an example of the moving
wall rather than a validated benchmark.
