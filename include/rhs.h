#ifndef RHS_H
#define RHS_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h"
#include "Metric.h"
#include  "BodyForces.h"

/**
 * @brief Computes the viscous contribution to the contravariant momentum RHS.
 *
 * This routine evaluates diffusive fluxes on the curvilinear grid and writes the
 * resulting term into `Visc`. The caller is responsible for providing compatible
 * vectors and for assembling any additional source terms afterwards.
 *
 * @param[in]  user  Block-level solver context containing metrics and model parameters.
 * @param[in]  Ucont Contravariant velocity field used by the discretization.
 * @param[in]  Ucat  Cartesian velocity field used for derivative evaluation.
 * @param[out] Visc  Output vector receiving the viscous RHS contribution.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc);

/**
 * @brief Computes the convective contribution to the contravariant momentum RHS.
 *
 * This routine evaluates the advection operator on the current velocity state and
 * stores the contribution in `Conv` for subsequent combination with viscous and
 * body-force terms.
 *
 * @param[in]  user  Block-level solver context containing metrics and numerics settings.
 * @param[in]  Ucont Contravariant velocity field used in face-normal flux construction.
 * @param[in]  Ucat  Cartesian velocity field used by the convective stencil.
 * @param[out] Conv  Output vector receiving the convection RHS contribution.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeConvection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv);

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
 * @param Rct  The PETSc Vec for the contravariant RHS, modified in place.
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
