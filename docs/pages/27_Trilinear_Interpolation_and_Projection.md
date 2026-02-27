@page 27_Trilinear_Interpolation_and_Projection Trilinear Interpolation and Particle-Grid Projection

PICurv performs bidirectional coupling between Eulerian and Lagrangian representations.

@tableofcontents

@section g2p_sec 1. Grid -> Particle Interpolation

- Corner/ghost-aware interpolation paths are used to evaluate fields at particle locations.
- Trilinear interpolation with eight-corner weights is the primary path for velocity/field transfer.

Code touchpoints:
- @ref TrilinearInterpolation_Scalar
- @ref TrilinearInterpolation_Vector
- @ref InterpolateEulerFieldToSwarm
- @ref InterpolateAllFieldsToSwarm

@section p2g_sec 2. Particle -> Grid Projection (Scatter/Averaging)

- Particle fields are accumulated to Eulerian storage.
- Per-cell normalization/averaging yields physically usable grid fields.
- This supports reconstructed Eulerian particle statistics and diagnostics.

Code touchpoints:
- Particle-to-grid scatter entry: @ref ScatterAllParticleFieldsToEulerFields
- Particle count accumulation: @ref CalculateParticleCountPerCell

@section refs_sec 3. References

- Trilinear interpolation background: https://en.wikipedia.org/wiki/Trilinear_interpolation
