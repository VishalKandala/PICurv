#ifndef RHS_H
#define RHS_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h"
#include "Metric.h" 

PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc);

PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv);
/**
 * @brief Applies a body force term to drive flow in channel simulations.
 *
 * This function calculates and adds a source term to the z-component of the
 * momentum equation's right-hand side (RHS). This term acts as a constant
 * pressure gradient to drive the flow, which is necessary for periodic channel
 * flow setups (controlled by the `channelz` flag in SimCtx).
 *
 * The force is only applied in fluid cells (where `nvert` is low). The function
 * modifies the input vector `Rct` in-place and updates the time-smoothed state
 * variable.
 *
 * @param user The UserCtx containing block-specific data and state variables.
 * @param Rct  The PETSc Vec for the contravariant RHS, which will be modified.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeBodyForce(UserCtx *user, Vec Rct);
/**
 * @brief Computes the Right-Hand Side (RHS) of the momentum equations.
 *
 * This function calculates the contribution of the convective and diffusive terms.
 * It is called by the momentum solvers (e.g., RungeKutta).
 *
 * @param user The UserCtx for a single block.
 * @param Rhs  The PETSc Vec where the RHS result will be stored.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode ComputeRHS(UserCtx *user, Vec Rhs);


#endif // RHS_H
