@page 28_IEM_and_Statistical_Averaging IEM Mixing and Statistical Averaging

@anchor _IEM_and_Statistical_Averaging

PICurv couples a particle scalar micromixing model (IEM-style) with particle
statistics postprocessing. Eulerian field statistics are a separate system with
its own configuration and state; see @ref 58_Field_Statistics.

@tableofcontents

@section p28_iem_sec 1. IEM Mixing Update In Current Code

For `Psi`, @ref UpdateParticleField uses:

\f[
\frac{d\Psi}{dt} = -\Omega(\Psi-\langle\Psi\rangle),
\qquad
\Omega = C_{IEM}\,\frac{\Gamma_{eff}}{\Delta^2},
\qquad
\Delta^2\approx V^{2/3}.
\f]

Closed-form update implemented in code:

\f[
\Psi^{n+1}=\langle\Psi\rangle + (\Psi^n-\langle\Psi\rangle)e^{-\Omega\Delta t}.
\f]

`C_IEM` is `solver.yml -> scalar_transport.iem_constant` (`-iem_constant`), default 2.0,
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
particle toward its own cell's mean, so it only mixes values that share a cell: a
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
0.52% of the uninterrupted run.

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

The unit tests in `tests/c/test_setup_lifecycle.c` pin the kernel: `TestConfiguredIEMUpdatesSwarm`
checks relaxation toward a prescribed mean for the default and an overridden constant;
`TestIEMRelaxesTowardOwnCellMean` checks that each particle relaxes toward its own cell's
scattered mean, conserving the scalar total; and `TestConfiguredParticleInitialValue`
checks that a configured value reaches every particle and the t=0 cell mean.

**Not validated:** the model. No turbulent flow has been run against a measured
scalar-variance decay, so `C_IEM = 2` and the mixing time scale
\f$V^{2/3}/(C_{IEM}\Gamma_{eff})\f$ are conventional choices, not calibrated ones.

@subsection p28_running_ssec 1.2 Running a Mixing Case

IEM needs only particles and a non-uniform `Psi`. This is the configuration of the
verification runs, on any case with particles:

```yaml
# case.yml
models:
  physics:
    particles:
      count: 32768
      init_mode: Volume
      fields:
        Psi: "where(uniform() < 0.5, 0, 1)"   # double delta inside every cell
# solver.yml
scalar_transport:
  iem_constant: 2.0                            # the default; 0 makes Psi a passive label
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
