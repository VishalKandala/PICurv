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
/**
 * @brief Computes the effective diffusivity scalar field (Gamma_eff) on the Eulerian grid.
 *
 * This function calculates the total diffusivity used to drive the stochastic 
 * motion of particles (Scalar FDF). It combines molecular diffusion and 
 * turbulent diffusion.
 *
 * Formula:
 *   Gamma_eff = nu/Sc + nu_t/Sc_t
 *
 * Where:
 * - nu = 1/Re (kinematic viscosity)
 * - nu_t (eddy viscosity from LES/RANS model)
 * - Sc (molecular Schmidt number)
 * - Sc_t (turbulent Schmidt number)
 *
 * @note If turbulence models are disabled, nu_t is assumed to be 0.
 * @note This function updates the local ghost values of lDiffusivity at the end 
 *       to ensure gradients can be computed correctly at subdomain boundaries.
 *
 * @param[in,out] user  Pointer to the user context containing grid data and simulation parameters.
 * 
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeEulerianDiffusivity(UserCtx *user);

/**
 * @brief Computes the Eulerian gradient of the effective diffusivity field.
 *
 * Reads the scalar diffusivity field and writes a vector gradient field used by
 * particle stochastic transport updates.
 *
 * @param[in,out] user Pointer to user context containing diffusivity vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeEulerianDiffusivityGradient(UserCtx *user);

#endif // RHS_H
