#ifndef RHS_H
#define RHS_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h"
#include "Metric.h"
#include  "BodyForces.h"

PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc);

PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv);

/**
 * @brief General dispatcher for applying all active body forces (momentum sources).
 *
 * This function serves as a central hub for adding momentum source terms to the
 * contravariant right-hand-side (Rct) of the momentum equations. It is called once per RHS
 * evaluation (e.g., once per Runge-Kutta stage).
 *
 * It introspects the simulation configuration to determine which, if any, body
 * forces are active and calls their specific implementation functions.
 *
 * @param user The UserCtx containing the simulation state for a single block.
 * @param Rhs  The PETSc Vec for the contravariant RHS, which will be modified in-place.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeBodyForces(UserCtx *user, Vec Rct);

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
