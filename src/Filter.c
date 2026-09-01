/**
 * @file Filter.c
 * @brief Implements numerical filtering schemes for Large Eddy Simulation (LES).
 *
 * This file contains the functions necessary to apply a "test filter" to the resolved
 * velocity field, which is a core component of the dynamic Smagorinsky turbulence model.
 * The choice of filter (e.g., a general 3D box filter or a specialized 2D
 * homogeneous filter) is determined by the simulation's configuration.
 */

#include "Filter.h"
#include <math.h> // Required for fabs() in the safe division check.

/*================================================================================*
 *                       INTERNAL (STATIC) FUNCTIONS                              *
 *================================================================================*/

/**
 * @brief Apply the two-dimensional Simpson stencil on the central homogeneous-plane samples.
 */
static double ApplySimpsonRuleHomogeneousFilter(double values[3][3][3])
{
	// The stencil only uses the central j-plane (j-index = 1, corresponding to the y-direction).
	// The formula is a weighted sum of the 9 points in that plane.
	const double corners = values[0][1][0] + values[2][1][0] + values[0][1][2] + values[2][1][2]; // 4 corner points
	const double edges   = values[0][1][1] + values[1][1][0] + values[2][1][1] + values[1][1][2]; // 4 edge-center points
	const double center  = values[1][1][1];                                                        // 1 center point

	// The weights (1, 4, 16) and normalization factor (36) come from the 2D Simpson's rule.
	return (corners + 4.0 * edges + 16.0 * center) / 36.0;
}


/**
 * @brief Average a 3x3x3 stencil with cell-volume weights while excluding solid cells.
 */
static double ApplyVolumeWeightedBoxFilter(double values[3][3][3], double weights[3][3][3])
{
    // v1...v8 store the sum of (value * weight) for each of the 8 sub-cubes.
    // w1...w8 store the sum of (weight) for each of the 8 sub-cubes.
	double v1, v2, v3, v4, v5, v6, v7, v8;
	double w1, w2, w3, w4, w5, w6, w7, w8;

	// --- Calculations for the 4 sub-cubes on the bottom layer (k-indices 0 and 1) ---

	// Bottom-Back-Left sub-cube (i-indices: 0,1; j-indices: 0,1; k-indices: 0,1)
	v1 = ( values[0][0][0]*weights[0][0][0] + values[1][0][0]*weights[1][0][0] + values[0][1][0]*weights[0][1][0] + values[1][1][0]*weights[1][1][0] +
           values[0][0][1]*weights[0][0][1] + values[1][0][1]*weights[1][0][1] + values[0][1][1]*weights[0][1][1] + values[1][1][1]*weights[1][1][1] );
	w1 = ( weights[0][0][0] + weights[1][0][0] + weights[0][1][0] + weights[1][1][0] +
           weights[0][0][1] + weights[1][0][1] + weights[0][1][1] + weights[1][1][1] );

	// Bottom-Back-Right sub-cube (i-indices: 1,2; j-indices: 0,1; k-indices: 0,1)
	v2 = ( values[1][0][0]*weights[1][0][0] + values[2][0][0]*weights[2][0][0] + values[1][1][0]*weights[1][1][0] + values[2][1][0]*weights[2][1][0] +
           values[1][0][1]*weights[1][0][1] + values[2][0][1]*weights[2][0][1] + values[1][1][1]*weights[1][1][1] + values[2][1][1]*weights[2][1][1] );
	w2 = ( weights[1][0][0] + weights[2][0][0] + weights[1][1][0] + weights[2][1][0] +
           weights[1][0][1] + weights[2][0][1] + weights[1][1][1] + weights[2][1][1] );

	// Bottom-Front-Left sub-cube (i-indices: 0,1; j-indices: 1,2; k-indices: 0,1)
	v3 = ( values[0][1][0]*weights[0][1][0] + values[1][1][0]*weights[1][1][0] + values[0][2][0]*weights[0][2][0] + values[1][2][0]*weights[1][2][0] +
           values[0][1][1]*weights[0][1][1] + values[1][1][1]*weights[1][1][1] + values[0][2][1]*weights[0][2][1] + values[1][2][1]*weights[1][2][1] );
	w3 = ( weights[0][1][0] + weights[1][1][0] + weights[0][2][0] + weights[1][2][0] +
           weights[0][1][1] + weights[1][1][1] + weights[0][2][1] + weights[1][2][1] );

	// Bottom-Front-Right sub-cube (i-indices: 1,2; j-indices: 1,2; k-indices: 0,1)
	v4 = ( values[1][1][0]*weights[1][1][0] + values[2][1][0]*weights[2][1][0] + values[1][2][0]*weights[1][2][0] + values[2][2][0]*weights[2][2][0] +
           values[1][1][1]*weights[1][1][1] + values[2][1][1]*weights[2][1][1] + values[1][2][1]*weights[1][2][1] + values[2][2][1]*weights[2][2][1] );
	w4 = ( weights[1][1][0] + weights[2][1][0] + weights[1][2][0] + weights[2][2][0] +
           weights[1][1][1] + weights[2][1][1] + weights[1][2][1] + weights[2][2][1] );


	// --- Calculations for the 4 sub-cubes on the top layer (k-indices 1 and 2) ---

	// Top-Back-Left sub-cube (i-indices: 0,1; j-indices: 0,1; k-indices: 1,2)
	v5 = ( values[0][0][1]*weights[0][0][1] + values[1][0][1]*weights[1][0][1] + values[0][1][1]*weights[0][1][1] + values[1][1][1]*weights[1][1][1] +
           values[0][0][2]*weights[0][0][2] + values[1][0][2]*weights[1][0][2] + values[0][1][2]*weights[0][1][2] + values[1][1][2]*weights[1][1][2] );
	w5 = ( weights[0][0][1] + weights[1][0][1] + weights[0][1][1] + weights[1][1][1] +
           weights[0][0][2] + weights[1][0][2] + weights[0][1][2] + weights[1][1][2] );

	// Top-Back-Right sub-cube (i-indices: 1,2; j-indices: 0,1; k-indices: 1,2)
	v6 = ( values[1][0][1]*weights[1][0][1] + values[2][0][1]*weights[2][0][1] + values[1][1][1]*weights[1][1][1] + values[2][1][1]*weights[2][1][1] +
           values[1][0][2]*weights[1][0][2] + values[2][0][2]*weights[2][0][2] + values[1][1][2]*weights[1][1][2] + values[2][1][2]*weights[2][1][2] );
	w6 = ( weights[1][0][1] + weights[2][0][1] + weights[1][1][1] + weights[2][1][1] +
           weights[1][0][2] + weights[2][0][2] + weights[1][1][2] + weights[2][1][2] );

	// Top-Front-Left sub-cube (i-indices: 0,1; j-indices: 1,2; k-indices: 1,2)
	v7 = ( values[0][1][1]*weights[0][1][1] + values[1][1][1]*weights[1][1][1] + values[0][2][1]*weights[0][2][1] + values[1][2][1]*weights[1][2][1] +
           values[0][1][2]*weights[0][1][2] + values[1][1][2]*weights[1][1][2] + values[0][2][2]*weights[0][2][2] + values[1][2][2]*weights[1][2][2] );
	w7 = ( weights[0][1][1] + weights[1][1][1] + weights[0][2][1] + weights[1][2][1] +
           weights[0][1][2] + weights[1][1][2] + weights[0][2][2] + weights[1][2][2] );

	// Top-Front-Right sub-cube (i-indices: 1,2; j-indices: 1,2; k-indices: 1,2)
	v8 = ( values[1][1][1]*weights[1][1][1] + values[2][1][1]*weights[2][1][1] + values[1][2][1]*weights[1][2][1] + values[2][2][1]*weights[2][2][1] +
           values[1][1][2]*weights[1][1][2] + values[2][1][2]*weights[2][1][2] + values[1][2][2]*weights[1][2][2] + values[2][2][2]*weights[2][2][2] );
	w8 = ( weights[1][1][1] + weights[2][1][1] + weights[1][2][1] + weights[2][2][1] +
           weights[1][1][2] + weights[2][1][2] + weights[1][2][2] + weights[2][2][2] );

    // Sum the contributions from all 8 octants.
	double total_weighted_value = v1+v2+v3+v4+v5+v6+v7+v8;
	double total_weight         = w1+w2+w3+w4+w5+w6+w7+w8;

    // Production safety check: avoid division by zero if all weights are zero
    // (e.g., if the stencil is entirely inside a solid body).
    if (fabs(total_weight) < 1.0e-12) {
        return 0.0;
    }

	return total_weighted_value / total_weight;
}


/*================================================================================*
 *                          PUBLIC FUNCTIONS                                      *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "ApplyLESTestFilter"
/**
 * @brief Implementation of \ref ApplyLESTestFilter().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Filter.h`.
 * @see ApplyLESTestFilter()
 */
double ApplyLESTestFilter(LESTestFilterKernel kernel, double values[3][3][3], double weights[3][3][3])
{
	// This function acts as a dispatcher over the configured stencil.
	switch (kernel) {
	case LES_TEST_FILTER_SIMPSON_IK:
		// The caller has declared the i and k directions homogeneous, so the more
		// accurate Simpson stencil applies. Volume weights are irrelevant there.
		return ApplySimpsonRuleHomogeneousFilter(values);
	case LES_TEST_FILTER_VOLUME_WEIGHTED_BOX:
	default:
		// Default for general, non-uniform, curvilinear grids. The volume-weighted
		// box filter keeps the average correct when neighbouring cells differ in size.
		return ApplyVolumeWeightedBoxFilter(values, weights);
	}
}

#undef __FUNCT__
#define __FUNCT__ "ApplyLESTestFilterSymTensor"
/**
 * @brief Implementation of \ref ApplyLESTestFilterSymTensor().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Filter.h`.
 * @see ApplyLESTestFilterSymTensor()
 */
PetscErrorCode ApplyLESTestFilterSymTensor(LESTestFilterKernel kernel, SymTensor values[3][3][3],
                                           double weights[3][3][3], SymTensor *filtered)
{
	// One scratch plane per component keeps the scalar filter unaware of tensors.
	double component[6][3][3][3];

	PetscFunctionBeginUser;
	PetscCheck(filtered != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
	           "Filtered tensor destination cannot be NULL.");

	for (PetscInt r = 0; r < 3; r++)
	for (PetscInt q = 0; q < 3; q++)
	for (PetscInt p = 0; p < 3; p++) {
		component[0][r][q][p] = values[r][q][p].xx;
		component[1][r][q][p] = values[r][q][p].xy;
		component[2][r][q][p] = values[r][q][p].xz;
		component[3][r][q][p] = values[r][q][p].yy;
		component[4][r][q][p] = values[r][q][p].yz;
		component[5][r][q][p] = values[r][q][p].zz;
	}

	filtered->xx = ApplyLESTestFilter(kernel, component[0], weights);
	filtered->xy = ApplyLESTestFilter(kernel, component[1], weights);
	filtered->xz = ApplyLESTestFilter(kernel, component[2], weights);
	filtered->yy = ApplyLESTestFilter(kernel, component[3], weights);
	filtered->yz = ApplyLESTestFilter(kernel, component[4], weights);
	filtered->zz = ApplyLESTestFilter(kernel, component[5], weights);

	PetscFunctionReturn(0);
}
