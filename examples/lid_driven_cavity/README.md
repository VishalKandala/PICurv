# Lid-driven cavity

A square cavity whose top wall moves at unit speed, the example of `handler: moving`. At
`Re = 100` the flow settles to one primary vortex with two small corner eddies. Botella and
Peyret (1998, Computers and Fluids 27) give its centreline velocity extrema to seven digits; the
older Ghia, Ghia and Shin (1982) table is 1.5 to 3.3% from them at this Reynolds number.

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
`solution_convergence.csv`, then take the extrema of `u` along `x = 0.5` and of `v` along
`y = 0.5` and compare them with Table 2 of Botella and Peyret.

## What was checked

Measurement `lid-driven-cavity-re100-2026-11-06`: on 64 x 64 and 128 x 128 cells the extrema
converge to Botella and Peyret's at second order, within 1.5e-3 and 3.2e-4.
