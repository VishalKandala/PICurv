# Driven Periodic Plane Channel

Streamwise-periodic channel flow held at a prescribed bulk flux by a driven
boundary handler. Five cases, in increasing cost:

| Directory | What it is | Grid | Cells |
|---|---|---|---|
| `laminar/` | exact closed-form verification | `plane_channel_laminar.cfg` | 16 x 32 x 16 |
| `les_retau180/` | constant-Smagorinsky LES | `plane_channel_les_retau180.cfg` | 32 x 64 x 32 |
| `les_wallmodel_retau1000/` | wall-modelled dynamic-Smagorinsky LES, `Re_tau ~ 1000` | `plane_channel_wm_retau1000.cfg` | 32 x 20 x 48 |
| `turbulent_retau180/` | DNS, `Re_tau = 180` | `plane_channel_retau180.cfg` | 128^3 (2.1M) |
| `turbulent_retau395/` | DNS, `Re_tau = 395` | `plane_channel_retau395.cfg` | 128 x 256 x 256 (8.4M) |

Every directory ships `case_constant_flux.yml`, `solver.yml`, `monitor.yml` and
`post.yml`; `laminar/` and `turbulent_retau180/` also ship
`case_initial_flux.yml` and `solver_newton_krylov.yml`.

## 1. Geometry and axis convention

All five cases use the same axis assignment:

- `i` (Xi) — spanwise, periodic (`geometric`)
- `j` (Eta) — wall-normal, no-slip walls at `y = 0` and `y = 2`
  (`les_wallmodel_retau1000/` applies a wall function on both)
- `k` (Zeta) — streamwise, periodic and **driven**

The half-height is `h = 1`, so the centreline is at `y = 1`. `length_ref` is `h`
and `velocity_ref` is the target bulk velocity `U_b = 1`, which makes the case
Reynolds number `U_b h / nu` and the nondimensional grid identical to the
physical one.

## 2. Grids are generated, not shipped

No `.picgrid` files live here. Each case sets `grid.mode: grid_gen` pointing at a
checked-in `.cfg` under `config/grids/`, so the grid is reproducible and
adjustable rather than an opaque binary.

Wall spacing is the knob you will actually turn. `first_cell_j_start` and
`first_cell_j_end` are expressed as a **fraction of the axis length**, not in
physical units, so for a wall-normal axis of length `Ly = 2`:

    first_cell_j = (y_plus_target / Re_tau) / 2

Retarget `y+` for a different `Re_tau` by editing those two numbers alone.
`stretch_j` is only the initial guess for the tanh fit; the generator solves for
the factor that actually realizes the requested first cell and reports it.

You do not have to trust that arithmetic by eye. Give the generator the case's
reference scales plus a target `--re-tau` and its `.info` report states the `y+`
the spacing actually realizes - `First_Cell_j_Start_Plus` alongside the solver-
and wall-unit sections. `picurv run`/`precompute` already pass `length_ref` and
`nu` from the case automatically; `re_tau` is a design target, so it stays a
generator argument. See **@subpage 48_Grid_Generator_Guide**, section 5.2.

Node counts (cells + 1) must stay **odd at every multigrid level**, because each
level coarsens as `IM -> (IM+1)/2`. The shipped node counts sit on the ladder
`5 -> 9 -> 17 -> 33 -> 65 -> 129 -> 257`, so the cell counts are 16, 32, 64, 128
and 256. Each config declares `mg_levels`, so the generator refuses a count that
breaks the ladder.

## 3. Choosing a handler

`case_constant_flux.yml` prescribes the flux; `case_initial_flux.yml` measures
the flux of the initial condition and holds that. They are otherwise identical.
For the shipped cases both routes target the same number, so running the pair is
a direct cross-check of the two handlers. Semantics, the control law, and the
restart contract: `docs/pages/54_Geometric_Periodic_Boundaries.md`, section 5.

`target_flux` is a **volumetric** flux — `U_b` times the cross-sectional area —
not a velocity:

| Case | Cross-section | `target_flux` |
|---|---|---|
| `laminar/` | `1 x 2` | `2.0` |
| `les_retau180/`, `turbulent_retau180/` | `2*pi x 2` | `4*pi = 12.566370614359172` |
| `turbulent_retau395/` | `pi x 2` | `2*pi = 6.283185307179586` |

## 4. Running

```bash
./bin/picurv run --solve -n 8 \
  --case    examples/periodic_test/driven_channel/laminar/case_constant_flux.yml \
  --solver  examples/periodic_test/driven_channel/laminar/solver.yml \
  --monitor examples/periodic_test/driven_channel/laminar/monitor.yml
```

Swap `solver.yml` for `solver_newton_krylov.yml` to run the same case under the
matrix-free Newton--Krylov momentum solver.

Rank count is bounded by multigrid depth, not just by available cores: every
level must leave each rank at least 3 nodes per axis (the stencil width when any
axis is periodic). The shipped level counts assume a moderate decomposition; see
`docs/pages/25_Pressure_Poisson_GMRES_Multigrid.md` for the formula.

## 5. Acceptance criteria

### 5.1 Laminar (`laminar/`) — exact verification

For half-height `h`, viscosity `nu` and the converged body force `f`:

    u(y)  = (f / 2nu) (h^2 - (y - h)^2)
    U_b   = f h^2 / 3nu      ->   f = 3 nu U_b / h^2
    u_tau = sqrt(f h)                    (exact, from the mean force balance)

With `h = 1`, `nu = 0.01`, `U_b = 1`:

| Quantity | Expected |
|---|---|
| `f` | `0.03` |
| `u_tau` | `sqrt(0.03) = 1.7320508e-01` |

Accept when the converged `f` matches `3 nu U_b / h^2` to solver tolerance,
`u_tau` from the computed wall shear matches `sqrt(f h)`, and the profile
matches the parabola at the expected order under refinement. The wall-normal
mesh is deliberately uniform here so the order check is clean.

This verifies the source term and the flux controller with no turbulence model
involved, which is why it comes first.

### 5.2 Turbulent DNS (`turbulent_retau180/`, `turbulent_retau395/`)

DNS resolution, no SGS model, so a flux-controller error is not conflated with
an SGS-model error. Reference: Moser, Kim & Mansour (1999), Phys. Fluids 11, 943.

Accept on the mean profile in wall units (`U+` vs `y+`) against the log law and
the DNS, the RMS fluctuation profiles, and the realized `u_tau = sqrt(f h)`:

| Case | Target `u_tau = Re_tau * nu` | Implied `f = u_tau^2 / h` |
|---|---|---|
| `Re_tau = 180` | `6.4286e-02` | `4.1327e-03` |
| `Re_tau = 395` | `5.7455e-02` | `3.3011e-03` |

### 5.3 LES (`les_retau180/`)

The coarse repeat of the `Re_tau = 180` case, run **after** the DNS so the SGS
contribution is assessed against an in-tree DNS rather than against literature
alone.

> **Experimental LES.** All four LES models (`constant_smagorinsky`,
> `dynamic_smagorinsky`, `vreman`, `wale`) are implemented and unit-tested, but none
> has a validated coefficient magnitude. This case selects `constant_smagorinsky`,
> which now applies its configured coefficient from the first step. The channel is
> periodic in xi and zeta, so `dynamic_smagorinsky` with `averaging.mode: homogeneous`
> is also available here and derives those two directions from the boundary pairs,
> giving a wall-normal coefficient profile. See `docs/pages/07_Case_Reference.md` and
> `docs/pages/72_LES_Turbulence_Closure.md`.

### 5.4 Wall-modelled LES (`les_wallmodel_retau1000/`)

`Re = U_b h / nu = 20,000`, which Dean's correlation puts at `Re_tau ~ 1016`. The
grid is uniform, so the wall layer is modelled rather than resolved: the first cell
centre sits at `y+ ~ 51` once the flow is developed, with `dx+ ~ 100` and
`dz+ ~ 133`. The closure is `dynamic_smagorinsky` averaged over `i` and `k`,
with a `werner` wall function and LES coefficient diagnostics every 10 steps.

`werner` rather than `log_law` is deliberate. The startup wall stress is near
laminar, so the first cell reads `y+ ~ 26` for the first steps (a 20-step pilot on
2026-09-24 measured `u_tau = 0.031`, mean `y+ = 31`), and the runtime stops a
`log_law` run after 10 samples below `y+ = 30`. Werner's two-layer law is valid
into the viscous sublayer. The same pilot took 4.4 s per step on 2 ranks of a debug
build, so the 50,000-step run is a cluster job.

Acceptance: realized `u_tau` against Dean (`0.0508 U_b`), `U+` against the log law
above the first cell, and resolved RMS profiles against Lee & Moser (2015) away from
the wall. None of these has been run.

## 6. Before you launch a campaign

One thing still stands between the turbulent case files and a finished validation;
see `docs/pages/54_Geometric_Periodic_Boundaries.md`, section 5.7.

**Momentum convergence** is no longer one of them. The pseudo-time stall observed on
2026-08-24 was re-characterized on 2026-09-18 and does not reproduce: the laminar case
converges every step, and at `Re = 10` reproduces the exact parabola at second order
under both driven handlers.

**Initial-condition seeding.** The turbulent and LES cases seed with
`channel_spectral_velocity`: a discretely divergence-free perturbation, zero on
the walls, of RMS `0.1 U_b` on a parabolic mean normalized to `U_b`. Precompute
(`picurv precompute --case ...`) reports the realized bulk velocity, RMS, and
divergence, and records plane spectra at three wall-normal stations. The DNS cases
use `k0 = 4`, `k_cut = 8`; the LES grid's coarse streamwise spacing bounds `k_cut`
at 5.3, so it uses `k0 = 2`, `k_cut = 5`.

Transition from this seed is not guaranteed. The perturbation grows linearly from
the wall, like a developed near-wall field: at `Re_tau = 180` its wall-parallel RMS
is about `0.013 U_b` at `y+ ~ 5` and `0.031 U_b` at `y+ ~ 13`, against `0.13 U_b` at
the centreline, whereas developed turbulence peaks near `0.17 U_b` at `y+ ~ 13`. Monitor fluctuation energy and
Reynolds shear stress, and if the flow relaminarizes, develop a field at coarse
resolution or from a precursor, carry it in with `mode: file`, and use
`initial_flux` to hold its flux.

## 7. Extracting profiles for DNS comparison

`monitor.yml` accumulates first and second moments of `Ucat` over a stationary
window. Turning those into `U+` vs `y+` needs an average over the two
homogeneous directions, and the postprocessor has no homogeneous-averaging task,
so that reduction is done outside it — see `tools/`.
