#ifndef POISSON_H
#define POISSON_H

#include "variables.h" // Provides definitions for UserCtx, UserMG, Vec, Mat, etc.
#include "rhs.h"       // Provides some primitives necessary.

/*================================================================================*
 *                     HIGH-LEVEL POISSON SOLVER                                  *
 *================================================================================*/

/**
 * @brief Solves the pressure-Poisson equation using a geometric multigrid method.
 *
 * This function orchestrates the entire multigrid V-cycle for the pressure
 * correction equation. It assembles the Laplacian matrix on all grid levels,
 * sets up the KSP solvers, smoothers, restriction/interpolation operators,
 * and executes the solve.
 *
 * @param usermg The UserMG context containing the entire multigrid hierarchy.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode PoissonSolver_MG(UserMG *usermg);


/*================================================================================*
 *                  CORE COMPONENTS OF THE POISSON SOLVER                         *
 *================================================================================*/

/**
 * @brief Assembles the Left-Hand-Side (LHS) matrix (Laplacian operator) for the
 *        Poisson equation on a single grid level.
 *
 * @param user The UserCtx for the grid level on which to assemble the matrix.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode PoissonLHSNew(UserCtx *user);

/**
 * @brief Computes the Right-Hand-Side (RHS) of the Poisson equation, which is
 *        the divergence of the intermediate velocity field.
 *
 * @param user The UserCtx for the grid level.
 * @param B    The PETSc Vec where the RHS result will be stored.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode PoissonRHS(UserCtx *user, Vec B);

/**
 * @brief Updates the pressure field `P` with the pressure correction `Phi`
 *        computed by the Poisson solver. (P = P + Phi)
 *
 * @param user The UserCtx containing the P and Phi vectors.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode UpdatePressure(UserCtx *user);

/**
 * @brief Corrects the contravariant velocity field `Ucont` to be divergence-free
 *        using the gradient of the pressure correction field `Phi`.
 *
 * @param user The UserCtx containing the Ucont and Phi vectors.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode Projection(UserCtx *user);


/*================================================================================*
 *                    MULTIGRID & NULL SPACE HELPERS                              *
 *================================================================================*/

/**
 * @brief The callback function for PETSc's MatNullSpace object.
 *
 * This function removes the null space from the Poisson solution vector by
 * ensuring the average pressure is zero, which is necessary for problems with
 * pure Neumann boundary conditions.
 *
 * @param nullsp The MatNullSpace context.
 * @param X      The vector to be corrected.
 * @param ctx    A void pointer to the UserCtx.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode PoissonNullSpaceFunction(MatNullSpace nullsp, Vec X, void *ctx);

/**
 * @brief The callback function for the multigrid restriction operator (MatShell).
 *
 * Defines the fine-to-coarse grid transfer for the Poisson residual.
 *
 * @param A The shell matrix context.
 * @param X The fine-grid source vector.
 * @param F The coarse-grid destination vector.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode MyRestriction(Mat A, Vec X, Vec F);

/**
 * @brief The callback function for the multigrid interpolation operator (MatShell).
 *
 * Defines the coarse-to-fine grid transfer for the pressure correction.
 *
 * @param A The shell matrix context.
 * @param X The coarse-grid source vector.
 * @param F The fine-grid destination vector.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode MyInterpolation(Mat A, Vec X, Vec F);


/*================================================================================*
 *                 IMMERSED BOUNDARY RELATED HELPERS (Optional)                   *
 *================================================================================*/

// These functions are called by the Poisson solver but are specific to the
// immersed boundary method. It is good practice to declare them here as they
// are part of the Poisson module's dependencies.

/**
 * @brief Calculates the net flux across the immersed boundary surface.
 * @param user      The UserCtx for the grid level.
 * @param ibm_Flux  (Output) The calculated net flux.
 * @param ibm_Area  (Output) The total surface area of the IB.
 * @param flg       A flag controlling the correction behavior.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg);

/**
 * @brief A specialized version of VolumeFlux, likely for reversed normals.
 * @param user      The UserCtx for the grid level.
 * @param ibm_Flux  (Output) The calculated net flux.
 * @param ibm_Area  (Output) The total surface area of the IB.
 * @param flg       A flag controlling the correction behavior.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode VolumeFlux_rev(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg);


#endif // POISSON_H
