@page 28_IEM_and_Statistical_Averaging IEM Mixing and Statistical Averaging

@anchor _IEM_and_Statistical_Averaging

PICurv couples a particle scalar micromixing model (IEM-style) with particle
statistics postprocessing. Eulerian field statistics are a separate system with
its own configuration and state; see @ref 58_Field_Statistics.

@tableofcontents

@section p28_iem_sec 1. IEM Mixing Update In Current Code

For `Psi`, @ref UpdateParticleField uses:

\f[
\frac{d\Psi}{dt} = -\Omega(\Psi-\tilde\Psi),
\qquad
\Omega = C_{IEM}\,\frac{\Gamma_{eff}}{\Delta^2},
\qquad
\Delta^2\approx V^{2/3}.
\f]

Closed-form update implemented in code:

\f[
\Psi^{n+1}=\tilde\Psi + (\Psi^n-\tilde\Psi)e^{-\Omega\Delta t},
\qquad
\tilde\Psi_p = m_c + \zeta_c\,\nabla\langle\Psi\rangle\cdot(\mathbf{x}_p-\mathbf{x}_c),
\qquad
m_c = \frac{\sum_{q\in c} w_q\Psi_q}{\sum_{q\in c} w_q},\quad
\mathbf{x}_c = \frac{\sum_{q\in c} w_q\mathbf{x}_q}{\sum_{q\in c} w_q},\quad
w_q = 1-e^{-\Omega_q\Delta t}.
\f]

The target \f$\tilde\Psi_p\f$ is the host cell's limited linear reconstruction (Jenny et al.
2001, Eq. 48; Viswanathan, Wang & Pope 2011). Its parts:

- \f$m_c\f$ and \f$\mathbf{x}_c\f$ are the cell's mean `Psi` and particle centroid, each weighted by
  the particle's relaxation weight \f$w_q\f$, the fraction of the way it moves toward its target
  in one step. With that weighting the update keeps every cell's `Psi` total exactly for any
  slope, even when particles in one cell relax at different rates (\f$\Gamma_{eff}\f$ varies
  between particles under LES).
- \f$\nabla\langle\Psi\rangle\f$ is the gradient of the plain cell means from the previous
  scatter: the central difference of neighbouring cell means, one-sided into the domain at a
  non-periodic boundary, mapped with the cell's metrics, and zero along an axis whose stencil
  reaches a cell with no particles.
- \f$\zeta_c\in[0,1]\f$ is the largest factor that keeps every target of the cell's particles
  within the minimum and maximum `Psi` of all particles at the start of the step (a global
  reduction). The unlimited slope carries a front's gradient into neighbouring pure cells and
  pushes their targets outside the scalar's range; in the duct-bend R2 pilot 9-18% of particles
  left [0, 1], by up to 0.07. With the limiter no particle leaves the range.

Relaxing toward the stepped mean alone would pull every
particle toward its own cell's value, which transports the resolved mean as extra diffusion:
by a factor \f$x\coth x\f$, \f$x = \tfrac{1}{2}\sqrt{C_{IEM}}\,h/\Delta\f$ for cell width
\f$h\f$ across the gradient (+16% at the default constant on cubic cells), which refinement does
not remove.

`C_IEM` is `case.yml -> models.physics.particles.iem_constant` (`-iem_constant`), default 2.0,
the conventional value; it is exposed for sensitivity studies, not because another value
has been validated. `C_IEM = 0` switches micromixing off: @ref UpdateAllParticleFields
skips the update rather than running it at a zero rate, because
\f$\langle\Psi\rangle + (\Psi^n-\langle\Psi\rangle)\f$ need not round to
\f$\Psi^n\f$ exactly, and a passive label must stay bit-identical. With a 0/1 label the
scattered cell mean is then the local fraction of particles that carried it, and the
`lost_psi_sum` column of `search_metrics.csv` counts the labelled particles leaving the
domain each step.

Code touchpoints:

- particle kernel: @ref UpdateParticleField
- per-field particle loop: @ref UpdateFieldForAllParticles
- stage wrapper: @ref UpdateAllParticleFields

The scalar to mix comes from `models.physics.particles.fields` (@ref p45_particle_values_sec),
which sets each particle's `Psi` before the first scatter; without it, every particle
starts at the catalog default 0.0 and the update has nothing to relax. IEM relaxes each
particle toward its own cell's reconstruction, so it only mixes values that share a cell: a
double-delta start inside every cell is `"where(uniform() < 0.5, 0, 1)"`, whereas a sharp
`half_space` region mixes only in the cells its boundary cuts, and spreads beyond them only
as fast as particles move.

@subsection p28_verification_ssec 1.1 Verification

**Verified** (measurement `iem-variance-decay-2026-10-01` in
`tests/tooling/measurement_records.json`): in zero flow from a random 0/1 start on 8^3 uniform
cells, the within-cell variance of `Psi` follows \f$e^{-2\Omega t}\f$ to 0.71% over its first
e-fold for `C_IEM` = 2, 20 and 200. The solver's whole decay curve, over up to 2.9 decades,
matches an independent NumPy emulation of the same step order to within 2.1 seed standard
deviations, and a run continued from a checkpoint with `restart_mode: load` stays within
0.52% of the uninterrupted run. The record's runs relaxed toward the stepped cell mean. On its
configuration, with `C_IEM` = 2, the reconstruction's variance ratio stays within 0.1% of the
stepped mean's through the first e-fold (local runs, not a registered measurement): with no
mean gradient, the reconstruction's slope is the gradient of sampling noise.

Below the first e-fold, the solver decays more slowly than \f$e^{-2\Omega t}\f$ (+4% at a
tenth of the initial variance with `C_IEM` = 20), and the emulation does the same. This is
the model working as intended. \f$e^{-2\Omega t}\f$ is the decay of a field whose cell means
are all equal; IEM never changes a cell mean, so any difference between neighbouring cell
means survives the mixing, and particles diffusing across a face carry that difference into
the next cell as new within-cell variance - the discrete form of scalar-variance production
by a mean gradient. In this test the differences are sampling noise of a random 0/1 start,
whose cell-mean variance is \f$1/(4N_{pc})\f$ for \f$N_{pc}\f$ particles per cell, so the
excess falls as \f$1/N_{pc}\f$: about 1% at 256 per cell and 0.2% at 1024 in the emulation.
Freezing the particles removes it exactly, and scattering the mean after the move instead
of using the previous step's leaves it unchanged, so neither the kernel nor the step order
causes it.

**Measured** (measurement `iem-cell-reconstruction-2026-11-08`, threshold not met): in a two-stream plug flow at 64
ranks, relaxing toward the stepped mean added transport of the resolved mean, D_eff/D 1.044, 1.089 and
1.171 at `C_IEM` = 1, 2 and 4 against a `C_IEM` = 0 reference, the x coth x factor within 0.9 of its
excess. With the reconstruction the excess over the reference is +0.0011 +- 0.0016 (`C_IEM` = 2) and
+0.0025 +- 0.0016 (`C_IEM` = 4), and the variance per squared mean gradient stays within 1.5% of the
stepped mean's. Against a tighter two-seed reference a residual is resolved: a plateau excess of
+0.0025 +- 0.0009 (`C_IEM` = 2) and +0.0038 +- 0.0009 (`C_IEM` = 4), decaying along the layer as a 1D
model of the reconstruction predicts (+0.012 at z = 0.5 to +0.005 at 3.5 for `C_IEM` = 2), though the
last station lies 2.9 standard errors below it. The reconstruction removes about 97% of the stepped
mean's excess, not all of it.

The unit tests in `tests/c/test_setup_lifecycle.c` pin the kernel: `TestConfiguredIEMUpdatesSwarm`
checks relaxation toward a prescribed mean for the default and an overridden constant;
`TestIEMRelaxesTowardOwnCellMean` checks that each particle relaxes toward its own cell's
scattered mean, conserving the scalar total; `TestIEMRelaxesTowardCellReconstruction` checks
that a linear field stays far closer to its line than relaxing to the stepped mean leaves it,
with each cell's sum kept, and `iem-reconstruction-across-seams` in `make unit-mpi` checks the same
at 2 and 4 ranks; and `TestConfiguredParticleInitialValue`
checks that a configured value reaches every particle and the t=0 cell mean.

**Verified** (measurement `iem-aprime-update-2026-11-15`): the limited, weighted target above.
In the zero-flow case of `iem-variance-decay-2026-10-01` (8^3 cells, 64 particles per cell,
2 ranks) the within-cell variance follows \f$e^{-2\Omega t}\f$ within 0.52% (`C_IEM` = 20) and
0.37% (`C_IEM` = 200) over its first e-fold, and matches a 64-seed NumPy emulation of this step to
within 2.66 and 0.91 seed standard deviations at every output step, over 2.1 and 2.5 decades.
`Psi` stays within [0, 1] to about 1e-18 and the `Psi` budget closes to about 1e-13 per block in
every duct-bend run. The unit tests `iem-keeps-cell-totals-with-varying-rates` and
`iem-stays-within-particle-extrema` in `make unit-mpi` fail without the weighting and the
limiter respectively.

**Cost of the limiter:** in the two-stream plug flow above, effective diffusivity rises by
+0.17% (`C_IEM` = 2) and +0.25% (`C_IEM` = 4) against the unlimited target, the same for every
seed, and falls only slowly with particles per cell (to +0.12% and +0.21% at 60 per cell). The
reading is that near a resolved front the scalar is curved within a cell, so the linear target
overshoots at the front's edges and the limiter flattens the slope there; it should fall with
front resolution, which is not tested. The fraction of cells with \f$\zeta_c<1\f$ is not
reported by any output. On a coupled 128^3 case with 105 M particles the step costs 6.8% more,
nearly all in the mixing update.

**Not validated:** the model. Two reference flows were run with this update
(records/2026-11-14_duct-bend_iem-model-checks-and-r1-c-sensitivity, ADDENDUM_A_results.md):
for a line source in grid turbulence (Warhaft 1984), whose scalar starts finer than a cell, the
centreline intensity is 3.4 times the measured value at `C_IEM` = 2 and still twice it at 3, so the
model under-mixes; for the temporal mixing layer of Colucci et al. (1998), where resolved scales
carry the variance, `C_IEM` = 2 leaves about half the DNS subgrid variance, so it over-mixes,
while 0.75 comes within about 10%. No single `C_IEM` fits both: the frequency
\f$C_{IEM}\Gamma_{eff}/V^{2/3}\f$ is tied to the cell size, not to the scalar's own sub-filter
scale. `C_IEM = 2` remains the conventional default, not a calibrated value. In both flows the
limited target changes the result by under 1% against the unlimited one.

@subsection p28_running_ssec 1.2 Running a Mixing Case

IEM needs only particles and a non-uniform `Psi`. This is the configuration of the
verification runs, on any case with particles:

```yaml
# case.yml
models:
  physics:
    particles:
      particles_per_cell: 64
      init_mode: Volume
      fields:
        Psi: {expression: "where(uniform() < 0.5, 0, 1)"}   # double delta inside every cell
      iem_constant: 2.0                       # the default; 0 makes Psi a passive label
```

Checkpoints store each particle's `Psi` (`particles/Psi.dat`). To see it, list `Psi` in
`post.yml` `io.particle_fields` for the particle `.vtp` files, and add a `nodal_average`
task from `Psi` to `Psi_nodal` for the scattered cell mean in the `.vts` files
(@ref p10_cap_eul_nodal_average_sub describes the task). `Psi` is also a statistics-window field
(@ref p58_scope_sec).

@subsection p28_restart_ssec 1.3 Restart

With `restart_mode: load`, each particle resumes with its saved `Psi`, the loaded state is
scattered before the first step so step one relaxes toward a mean that has seen it, and
`fields` is ignored with a warning. With `restart_mode: init`, the population is reseeded
and `fields` applies again, so any mixing already done is discarded.
@ref p45_restart_matrix_sec has the full matrix.

@subsection p28_troubleshoot_ssec 1.4 Troubleshooting

- **`Psi` never changes.** Either every particle in a cell holds the same value (the
  default 0.0 when `fields` is absent, or a region whose boundary cuts no cell), or
  `iem_constant` is 0, or `solver.yml` `verification.sources.scalar` is set, which
  prescribes `Psi` and bypasses the update.
- **Variance decays more slowly than \f$e^{-2\Omega t}\f$.** Cell means differ, and
  particles crossing faces turn those differences into within-cell variance
  (@ref p28_verification_ssec). In a resolved mixing problem that is real variance
  production; when the differences are only sampling noise, more particles per cell
  reduce it in proportion.
- **The global mean of `Psi` drifts.** IEM conserves each cell's total, so a drift comes from
  particles carrying `Psi` out of the domain; `lost_psi_sum` in `search_metrics.csv` counts
  it each step.

@section p28_dataflow_sec 2. Required Dataflow For IEM

IEM update requires:

- per-particle diffusivity from swarm fields,
- host-cell IDs for indexing,
- Eulerian mean field (`user->lPsi`) and Jacobian (`user->lAj`) for volume scaling.

This means scatter/interpolation order matters: stale Eulerian means produce stale IEM forcing.

@section p28_stats_sec 3. Statistics Pipeline (Postprocessor)

Current primary reduction kernel:

- @ref ComputeParticleMSD

Implemented MSD physics includes:

\f[
D = \frac{1}{Re\,Sc},
\qquad
r_{theory} = \sqrt{6Dt},
\f]

with global MPI reductions and CSV output per statistics call.

@section p28_field_observations_sec 4. Relationship To Eulerian Field Statistics

Eulerian field statistics are configured at `monitor.yml -> field_statistics`,
accumulate weighted centered moments while the solver runs, and are derived into
Reynolds stresses, RMS, turbulent kinetic energy, and fluxes by
`post.yml -> field_statistics`. @ref 58_Field_Statistics is the full contract.

Particle MSD stays under the postprocessor `statistics_pipeline` described above.
It is not an Eulerian field window, shares no accumulator state, and is
configured separately. The two are kept apart deliberately so a particle
reduction and a resolved-turbulence average are never conflated.

`Psi` is the field where the two subsystems meet: it is a particle-carried scalar
projected onto the grid, and the projected field can itself be accumulated as a
statistics window when particles are active.

@section p28_terminology_sec 5. Averaging Terminology In PICurv

"Averaging" appears in multiple contexts:

- particle-to-grid count-normalized scatter,
- Eulerian field-statistics windows,
- postprocessing global statistical reductions.

Treat these as distinct workflows with different configuration points.

@section p28_refs_sec 6. Related Pages

- **@subpage 27_Trilinear_Interpolation_and_Projection**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 58_Field_Statistics**
