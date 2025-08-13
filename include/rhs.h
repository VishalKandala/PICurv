#ifndef RHS_H
#define RHS_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h" 
/*================================================================================*
 *                       CORE NUMERICAL KERNELS                                   *
 *================================================================================*/

PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc);

PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv);

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
extern PetscErrorCode FormFunction1(UserCtx *user, Vec Rhs);


#endif // RHS_H
