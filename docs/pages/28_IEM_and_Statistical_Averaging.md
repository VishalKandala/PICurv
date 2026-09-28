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
starts at the catalog default 0.0 and the update has nothing to relax. A double-delta
start, for example, is a sharp `half_space` region with `value: 1` over `background: 0`.

@warning **Coupled scalar-variance decay is not yet verified.** The verification scalar
source (@ref p08_verification_sec) prescribes `Psi` exactly and bypasses this update, and
the scatter of `Psi` to the grid is verified through it. `TestConfiguredIEMUpdatesSwarm`
in `tests/c/test_setup_lifecycle.c` checks relaxation toward a prescribed mean for the
default and an overridden constant, and `TestConfiguredParticleInitialValue` checks that a
configured value reaches every particle and the t=0 cell mean. Neither establishes the
decay of scalar variance in a production flow.

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
