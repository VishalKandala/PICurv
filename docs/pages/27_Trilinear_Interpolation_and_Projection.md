@page 27_Trilinear_Interpolation_and_Projection Trilinear Interpolation and Particle-Grid Projection

Eulerian-Lagrangian coupling in PICurv is built from interpolation (grid -> particle) and scatter/normalization (particle -> grid).

@tableofcontents

@section g2p_sec 1. Grid -> Particle Interpolation

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

@section p2g_sec 2. Particle -> Grid Scatter and Normalization

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

@section coupling_sec 3. Accuracy and Stability Considerations

Coupling quality depends on:

- accurate host-cell IDs and interpolation weights,
- ghost synchronization before interpolation/scatter,
- consistent field DOF/DM association,
- avoiding stale particle-count vectors between scatter calls.

Inconsistency in any of these usually appears as noisy particle statistics or nonphysical reconstructed Eulerian fields.

@section refs_sec 4. Related Pages

- **@subpage 26_Walking_Search_Method**
- **@subpage 28_IEM_and_Statistical_Averaging**
- **@subpage 34_Particle_Model_Overview**
