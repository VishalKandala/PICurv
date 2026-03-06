@page 27_Trilinear_Interpolation_and_Projection Trilinear Interpolation and Particle-Grid Projection

@anchor _Trilinear_Interpolation_and_Projection

Eulerian-Lagrangian coupling in PICurv is built from interpolation (grid -> particle) and scatter/normalization (particle -> grid).

@tableofcontents

@section p27_g2p_sec 1. Grid -> Particle Interpolation

For a particle inside a host cell with local coordinates \f$(a_1,a_2,a_3)\in[0,1]^3\f$, trilinear interpolation uses 8 corner weights:

\f[
\phi_p = \sum_{m=1}^{8} w_m(a_1,a_2,a_3)\,\phi_m,
\qquad
\sum_m w_m = 1.
\f]

Code touchpoints:

- scalar kernel: @ref TrilinearInterpolation_Scalar
- vector kernel: @ref TrilinearInterpolation_Vector
- field-by-field wrapper: @ref InterpolateEulerFieldToSwarm
- stage wrapper: @ref InterpolateAllFieldsToSwarm

Implementation note: interpolation path includes center-to-corner staging before final particle evaluation.

@section p27_p2g_sec 2. Particle -> Grid Scatter and Normalization

Scatter computes per-cell sums, then normalizes by particle count:

\f[
\bar{\phi}_{i,j,k} =
\begin{cases}
\dfrac{\sum_{p\in cell(i,j,k)} \phi_p}{N_{i,j,k}}, & N_{i,j,k}>0, \\
0, & N_{i,j,k}=0.
\end{cases}
\f]

Code touchpoints:

- particle count accumulation: @ref CalculateParticleCountPerCell
- scatter orchestrator: @ref ScatterAllParticleFieldsToEulerFields

Current standard scatter path actively maps particle `Psi` to Eulerian `Psi`; additional fields are scaffolded and can be enabled with matching vector/DM contracts.

@section p27_coupling_sec 3. Accuracy and Stability Considerations

Coupling quality depends on:

- accurate host-cell IDs and interpolation weights,
- ghost synchronization before interpolation/scatter,
- consistent field DOF/DM association,
- avoiding stale particle-count vectors between scatter calls.

Inconsistency in any of these usually appears as noisy particle statistics or nonphysical reconstructed Eulerian fields.

@section p27_refs_sec 4. Related Pages

- **@subpage 26_Walking_Search_Method**
- **@subpage 28_IEM_and_Statistical_Averaging**
- **@subpage 34_Particle_Model_Overview**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Trilinear Interpolation and Particle-Grid Projection** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

