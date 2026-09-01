# Driven Periodic Square Duct

Streamwise-periodic turbulent flow through a square duct at `Re_b = 4410`, held
at a prescribed bulk flux by a driven boundary handler.

A square duct is the case worth running after the plane channel, because it is
not a plane channel with corners. The anisotropy of the Reynolds stresses in the
cross-plane drives a **secondary flow of the second kind**: a set of eight
counter-rotating cells that carry fluid along the corner bisectors toward the
corners and back along the wall bisectors. A plane channel has no such motion.
It is small — 1–3% of `U_b` — which is exactly why it is a demanding test: it is
a weak residual of a long average, and a flux controller with a systematic bias
will drown it.

## 1. Geometry and axis convention

- `i` (Xi) — cross-stream, no-slip walls at `x = 0` and `x = 2`
- `j` (Eta) — cross-stream, no-slip walls at `y = 0` and `y = 2`
- `k` (Zeta) — streamwise, periodic and **driven**, length `4*pi`

The duct half-width is `a = 1`, so the side length and hydraulic diameter are
`D = 2a = 2`. Grid: `config/grids/square_duct_reb4410.cfg`, 129 x 129 x 257
cells (4.3M), with two-sided tanh clustering on **both** walled directions.
Clustering the two identically keeps the corner region symmetric, which matters
because the corner is where the secondary flow lives.

## 2. Reynolds number bookkeeping

Read this carefully; the factor of two is easy to lose.

| Quantity | Value |
|---|---|
| Case `Re` (what `case.yml` produces) | `U_b a / nu = 2205` |
| Literature `Re_b` | `U_b D / nu = 4410` |
| `viscosity` in `case.yml` | `4.5351474e-04` |
| Nominal `Re_tau` | `~300` |

`length_ref` is the half-width `a`, so the case Reynolds number is built on `a`,
while Gavrilakis quotes it on the hydraulic diameter `D = 2a`.

## 3. References

- Gavrilakis (1992), JFM **244**, 101 — `Re_b = 4410`, `Re_tau ≈ 300`
- Huser & Biringen (1993), JFM **257**, 65 — `Re_tau ≈ 600`

## 4. Running

```bash
./bin/picurv run --solve -n 8 \
  --case    examples/periodic_test/driven_duct/case_constant_flux.yml \
  --solver  examples/periodic_test/driven_duct/solver.yml \
  --monitor examples/periodic_test/driven_duct/monitor.yml
```

`case_initial_flux.yml` runs the same case under the `initial_flux` handler, and
`solver_newton_krylov.yml` under the matrix-free Newton--Krylov momentum solver.
`target_flux` is a volumetric flux: `U_b * Lx * Ly = 1 * 2 * 2 = 4.0`.

## 5. Acceptance criteria

1. **Mean streamwise profile** along the wall bisector and along the corner
   bisector, in wall units, against Gavrilakis. The two bisectors differ, and
   that difference is the point.
2. **Secondary flow** present in the cross-plane with the correct eight-cell
   topology and a magnitude of **1–3% of `U_b`**. Absence of the secondary flow,
   or a magnitude an order of magnitude off, is a failure even if the mean
   streamwise profile looks right.
3. **Wall shear stress distribution** along a wall, which should peak near the
   wall-bisector and fall toward the corners.

The statistics window in `monitor.yml` is deliberately long. The secondary flow
is a small residual of a long average; a short window will show noise where the
cells should be.

## 6. Before you launch

The two blockers described in
`examples/periodic_test/driven_channel/driven_channel.md` section 6 —
pseudo-time momentum convergence on periodic wall-bounded flow, and the absence
of a wall-respecting perturbed initial-condition generator — apply here too, and
bite harder: a duct will not develop turbulence from a laminar seed at this
Reynolds number without a finite-amplitude perturbation.

The momentum-convergence blocker carries a status caveat: it was measured on
2026-08-24 and **requires re-characterization at current `HEAD`** after the
2026-08-25 convergence-criterion change. See
`docs/pages/54_Geometric_Periodic_Boundaries.md`, section 5.7, for the canonical
dated statement. The initial-condition gap is unaffected and still stands.

Validate the plane channel first. It is cheaper, it has an exact laminar answer,
and it isolates the flux controller from everything the duct adds.
