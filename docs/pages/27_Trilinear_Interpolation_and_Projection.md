@page 27_Trilinear_Interpolation_and_Projection Trilinear Interpolation and Particle-Grid Projection

@anchor _Trilinear_Interpolation_and_Projection

Eulerian-Lagrangian coupling in PICurv is built from interpolation (grid -> particle) and scatter/normalization (particle -> grid).

@tableofcontents

@section p27_g2p_sec 1. Grid -> Particle Interpolation

For a particle inside a host cell with local coordinates \f$(a_1,a_2,a_3)\in[0,1]^3\f$, trilinear interpolation uses 8 source-point weights:

\f[
\phi_p = \sum_{m=1}^{8} w_m(a_1,a_2,a_3)\,\phi_m,
\qquad
\sum_m w_m = 1.
\f]

PICurv supports two grid-to-particle interpolation methods, selectable at runtime via the `interpolation_method` YAML key (or `-interpolation_method` flag).

@subsection p27_trilinear_ssec 1.1 Trilinear (Direct Cell-Center) — Default

The recommended method. For each particle, the 8 nearest cell centers are identified via octant detection (based on existing host-cell weights), forming a "dual cell." Trilinear interpolation is performed directly from these cell-centered field values to the particle position.

Key properties:

- **Second-order accurate** on both uniform and curvilinear grids.
- **No intermediate staging**: operates directly on the ghosted cell-centered field (`lUcat`, `lDiffusivity`, etc.) and ghosted cell center coordinates (`lCent`). No extra ghost exchange is needed beyond what the solver already provides.
- **Boundary handling**: at non-periodic boundaries, the dual-cell octant is clamped to use only physical cell centers, and trilinear weights are left unclamped (may exceed [0,1]) to provide second-order linear extrapolation. Weights still sum to 1.0.
- **Periodic boundaries**: PETSc supplies wrapped field indices, while
  `ApplyPeriodicCorrectionsToCellCentersAndSpacing()` supplies translated
  coordinate images. This interpolation path is Eulerian-only; particles are
  not currently wrapped across periodic boundaries.

Code path: @ref InterpolateEulerFieldToSwarm dispatches to `InterpolateEulerFieldFromCenterToSwarm`.

@subsection p27_corner_averaged_ssec 1.2 Corner-Averaged (Legacy)

The original two-stage method:

1. **Center -> Corner**: average surrounding cell-center values to each grid node via @ref InterpolateFieldFromCenterToCorner (unweighted `sum/count`).
2. **Corner -> Particle**: standard trilinear interpolation from the 8 corner values of the host cell.

This is second-order on uniform Cartesian grids but degrades to first-order on curvilinear grids because the unweighted averaging does not account for the asymmetric placement of cell centers around nodes.

Code path: @ref InterpolateEulerFieldToSwarm dispatches to `InterpolateEulerFieldFromCornerToSwarm`.

@subsection p27_shared_ssec 1.3 Shared Kernels

Both methods share the same underlying trilinear kernels:

- scalar kernel: @ref TrilinearInterpolation_Scalar
- vector kernel: @ref TrilinearInterpolation_Vector
- weight computation: @ref ComputeTrilinearWeights
- field-by-field dispatcher: @ref InterpolateEulerFieldToSwarm
- stage wrapper: @ref InterpolateAllFieldsToSwarm

@subsection p27_config_ssec 1.4 Configuration

YAML (`solver.yml`):

```yaml
interpolation:
  method: "Trilinear"     # default; or "CornerAveraged"
```

C flag: `-interpolation_method 0` (Trilinear) or `-interpolation_method 1` (CornerAveraged).

Enum: `InterpolationMethod` in `include/variables.h`. Stored in `SimCtx.interpolationMethod`.

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
