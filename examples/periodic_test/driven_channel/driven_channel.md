# Driven Periodic Plane Channel

Streamwise-periodic channel flow held at a prescribed bulk flux by a driven
boundary handler. Four cases, in increasing cost:

| Directory | What it is | Grid | Cells |
|---|---|---|---|
| `laminar/` | exact closed-form verification | `plane_channel_laminar.cfg` | 17 x 33 x 17 |
| `les_retau180/` | constant-Smagorinsky LES | `plane_channel_les_retau180.cfg` | 33 x 65 x 33 |
| `turbulent_retau180/` | DNS, `Re_tau = 180` | `plane_channel_retau180.cfg` | 129^3 (2.15M) |
| `turbulent_retau395/` | DNS, `Re_tau = 395` | `plane_channel_retau395.cfg` | 129 x 257 x 257 (8.5M) |

Every directory ships `case_constant_flux.yml`, `solver.yml`, `monitor.yml` and
`post.yml`; `laminar/` and `turbulent_retau180/` also ship
`case_initial_flux.yml` and `solver_newton_krylov.yml`.

## 1. Geometry and axis convention

All four cases use the same axis assignment:

- `i` (Xi) — spanwise, periodic (`geometric`)
- `j` (Eta) — wall-normal, no-slip walls at `y = 0` and `y = 2`
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

Cell counts must stay **odd at every multigrid level**, because each level
coarsens as `IM -> (IM+1)/2`. The shipped counts sit on the ladder
`5 -> 9 -> 17 -> 33 -> 65 -> 129 -> 257`.

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
alone. Use **constant** Smagorinsky only — the dynamic model in this tree has a
known defect (incorrect `M_ij`, no homogeneous averaging).

## 6. Before you launch a campaign

Two things currently stand between these case files and a finished validation.
Both are described in full in `docs/pages/54_Geometric_Periodic_Boundaries.md`,
section 5.7.

**Momentum convergence.** Under the Dual Time Picard Jameson RK solver these
periodic wall-bounded cases do not reach the pseudo-time tolerance: `dtau`
collapses to its floor and every step logs "reached N total attempts without
convergence", continuing from the last accepted state. This is *not* caused by
the driven handlers — an otherwise identical case with plain `geometric`
handlers behaves the same, and the Poisson side stays healthy throughout
(maximum divergence around 1e-14). It must be resolved before the turbulent
runs are meaningful.

**Initial-condition seeding.** The shipped cases seed with
`streamwise_constant`, which is a laminar profile. Transition from it is slow
and is not guaranteed without a finite-amplitude perturbation. There is no
in-tree generator that produces a perturbed divergence-free field respecting
no-slip: `spectral_random_velocity` explicitly requires `PERIODIC`/`geometric`
on all six faces and is rejected here, correctly, because it cannot respect a
wall. Until such a generator exists, develop a field at coarse resolution or
from a precursor, carry it in, and use `initial_flux` to hold its flux.

## 7. Extracting profiles for DNS comparison

`monitor.yml` accumulates first and second moments of `Ucat` over a stationary
window. Turning those into `U+` vs `y+` needs an average over the two
homogeneous directions, and the postprocessor has no homogeneous-averaging task,
so that reduction is done outside it — see `tools/`.
