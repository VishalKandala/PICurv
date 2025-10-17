#ifndef FILTER_H
#define FILTER_H

#include  "variables.h"
#include  "logging.h"
/**
 * @brief Applies a numerical "test filter" to a 3x3x3 stencil of data points.
 *
 * EXPLANATION of "Test Filter":
 * In the dynamic Smagorinsky LES model, the model coefficient (Cs) is not constant.
 * It is calculated "dynamically" by sampling the turbulent energy cascade at two
 * different scales. The first scale is the grid filter itself (implicit in the discretization).
 * The second, larger scale is this explicit "test filter". By comparing the behavior
 * of the resolved flow at these two scales, the model can determine the appropriate
 * amount of dissipation needed, and thus calculate the Smagorinsky coefficient.
 *
 * This function is the primary interface for this process. It selects the appropriate
 * filtering algorithm based on simulation settings.
 *
 * @param simCtx The global simulation context, containing flags like `testfilter_ik`.
 * @param values The 3x3x3 array of scalar values at the stencil points.
 * @param weights The 3x3x3 array of weights (typically 1/cell volume) for averaging.
 * @return The resulting scalar value after applying the test filter.
 */
double ApplyLESTestFilter(const SimCtx *simCtx, double values[3][3][3], double weights[3][3][3]);
#endif // FILTER_H