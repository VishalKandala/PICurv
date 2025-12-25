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
/**
 * @brief Computes the effective diffusivity scalar field ($\Gamma_{eff}$) on the Eulerian grid.
 *
 * This function calculates the total diffusivity used to drive the stochastic 
 * motion of particles (Scalar FDF). It combines molecular diffusion and 
 * turbulent diffusion.
 *
 * **Formula:**
 * \f[
 *    \Gamma_{eff} = \underbrace{\frac{\nu}{Sc}}_{\text{Molecular}} + \underbrace{\frac{\nu_t}{Sc_t}}_{\text{Turbulent}}
 * \f]
 *
 * Where:
 * - \f$ \nu = 1/Re \f$ (Kinematic Viscosity)
 * - \f$ \nu_t \f$ (Eddy Viscosity from LES/RANS model)
 * - \f$ Sc \f$ (Molecular Schmidt Number, user-defined)
 * - \f$ Sc_t \f$ (Turbulent Schmidt Number, user-defined)
 *
 * @note If turbulence models are disabled, \f$ \nu_t \f$ is assumed to be 0.
 * @note This function updates the local ghost values of lDiffusivity at the end 
 *       to ensure gradients can be computed correctly at subdomain boundaries.
 *
 * @param[in,out] user  Pointer to the user context containing grid data and simulation parameters.
 * 
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeEulerianDiffusivity(UserCtx *user);

#endif // RHS_H
