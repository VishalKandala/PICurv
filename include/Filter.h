#ifndef FILTER_H
#define FILTER_H

#include  "variables.h"
#include  "logging.h"

/**
 * @brief Applies a numerical "test filter" to a 3x3x3 stencil of data points.
 *
 * EXPLANATION of "Test Filter":
 * In the dynamic Smagorinsky LES model the model coefficient is not prescribed. It
 * is calculated dynamically by sampling the turbulent energy cascade at two
 * different scales. The first scale is the grid filter itself, implicit in the
 * discretization. The second, larger scale is this explicit "test filter". By
 * comparing how the resolved flow behaves at the two scales, the model determines
 * how much dissipation the unresolved motions require.
 *
 * The kernel is selected by the caller rather than read from the simulation
 * context, so a filter can be exercised on 27 numbers without a live run.
 *
 * @param[in] kernel  Discrete stencil to apply.
 * @param[in] values  The 3x3x3 array of scalar values at the stencil points,
 *                    indexed `[k][j][i]` in the DMDA convention.
 * @param[in] weights The 3x3x3 array of cell volume weights. Ignored by
 *                    ::LES_TEST_FILTER_SIMPSON_IK, which assumes uniform spacing
 *                    in its homogeneous plane.
 * @return The filtered value, or zero when every weight vanishes because the
 *         stencil lies entirely inside a solid body.
 */
double ApplyLESTestFilter(LESTestFilterKernel kernel, double values[3][3][3], double weights[3][3][3]);

/**
 * @brief Applies the test filter to all six components of a symmetric tensor.
 *
 * Filtering is linear, so applying it componentwise is exact; this wrapper exists
 * so callers state the intent once instead of repeating six near-identical calls
 * whose only difference is easy to mistype.
 *
 * @param[in]  kernel   Discrete stencil to apply.
 * @param[in]  values   Per-stencil-point tensors, indexed `[k][j][i]`.
 * @param[in]  weights  Cell volume weights over the same stencil.
 * @param[out] filtered Componentwise filtered tensor.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyLESTestFilterSymTensor(LESTestFilterKernel kernel, SymTensor values[3][3][3],
                                           double weights[3][3][3], SymTensor *filtered);

#endif // FILTER_H
