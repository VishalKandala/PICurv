#ifndef LES_H
#define LES_H

#include "variables.h"
#include "logging.h"
#include "Metric.h"
#include "setup.h"
#include "Filter.h"

/**
 * @brief Computes the dynamic Smagorinsky constant (Cs) for the LES model.
 *
 * This function implements the dynamic procedure for calculating the Smagorinsky
 * constant at each grid point. It involves test-filtering the velocity field and
 * strain rate tensor to determine the constant locally. The calculation can be
 * computationally expensive and is typically run less frequently than every
 * time step, controlled by `simCtx->dynamic_freq`.
 *
 *
 * @param user The user context for a single computational block.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeSmagorinskyConstant(UserCtx *user);

/**
 * @brief Computes the turbulent eddy viscosity (Nu_t) for the LES model.
 *
 * Using the Smagorinsky constant (Cs) computed by `ComputeSmagorinskyConstant`,
 * this function calculates the eddy viscosity at each grid point based on the
 * local strain rate magnitude and filter width.
 *
 * nu_t = (Cs * filter_width)^2 * |S|*
 *
 * @param user The user context for a single computational block.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeEddyViscosityLES(UserCtx *user);
#endif // LES_H