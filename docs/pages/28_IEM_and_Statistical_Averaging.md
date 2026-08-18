@page 28_IEM_and_Statistical_Averaging IEM Mixing and Statistical Averaging

@anchor _IEM_and_Statistical_Averaging

PICurv couples a particle scalar micromixing model (IEM-style) with particle
statistics postprocessing. The planned Eulerian field-statistics pipeline is a
separate system.

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

Code touchpoints:

- particle kernel: @ref UpdateParticleField
- per-field particle loop: @ref UpdateFieldForAllParticles
- stage wrapper: @ref UpdateAllParticleFields

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

@section p28_field_observations_sec 4. Eulerian Field Statistics Status

Eulerian scientific statistics are not currently exposed in active YAML. The
old `models.statistics.time_averaging`, raw `su0`/`su1`/`su2`/`sp` state, and
reserved averaged-field post output were removed rather than retained as a
second workflow. Until the replacement is complete, users can compute offline
statistics only from the instantaneous states they elected to save.

The replacement design uses field-catalog-aware pointwise weighted centered
first and second moments with independent state per named window. Reynolds
stress, RMS, and TKE remain postprocessor-derived quantities. See the
authoritative future specification for the full contract and staged status.

Particle MSD remains under the existing postprocessor `statistics_pipeline`.
It is not an Eulerian field window and does not share accumulator state.

@section p28_terminology_sec 5. Averaging Terminology In PICurv

"Averaging" appears in multiple contexts:

- particle->grid count-normalized scatter,
- future Eulerian field-statistics windows,
- postprocessing global statistical reductions.

Treat these as distinct workflows with different configuration points.

@section p28_refs_sec 6. Related Pages

- **@subpage 27_Trilinear_Interpolation_and_Projection**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 58_Turbulence_Statistics_Pipeline_Specification**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **IEM Mixing and Statistical Averaging** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
