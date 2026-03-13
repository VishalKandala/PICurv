#ifndef POISSON_H
#define POISSON_H

#include "variables.h" // Provides definitions for UserCtx, UserMG, Vec, Mat, etc.
#include "Metric.h"       // Provides some primitives necessary.
#include "Boundaries.h"
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
 *
 * @note Testing status:
 *       This routine is exercised in runtime smoke, but still needs deeper
 *       direct bespoke coverage for debugging and branch isolation.
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
 *
 * @note Testing status:
 *       Direct unit coverage exists for core operator assembly, but periodic and
 *       immersed-boundary stencil branches remain thinner than the Cartesian baseline.
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
 * @brief Enforces a constant volumetric flux profile along the entire length of a driven periodic channel.
 *
 * This function is a "hard" corrector, called at the end of the projection step.
 * The projection ensures the velocity field is divergence-free (3D continuity), but this
 * function enforces a stricter 1D continuity condition (`Flux(plane) = constant`)
 * required for physically realistic, fully-developed periodic channel/pipe flow.
 *
 * The process is as follows:
 * 1.  Introspects the boundary condition handlers to detect if a `DRIVEN_` flow is active
 *     and in which direction ('X', 'Y', or 'Z'). If none is found, it exits.
 * 2.  Measures the current volumetric flux through *every single cross-sectional plane*
 *     in the driven direction.
 * 3.  For each plane, it calculates the velocity correction required to make its flux
 *     match the global `targetVolumetricFlux` (which was set by the controller).
 * 4.  It applies this spatially-uniform (but plane-dependent) velocity correction directly
 *     to the `ucont` field, ensuring `Flux(plane) = TargetFlux` for all planes.
 *
 * @param user The UserCtx containing the simulation state for a single block.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CorrectChannelFluxProfile(UserCtx *user);

/**
 * @brief Corrects the contravariant velocity field `Ucont` to be divergence-free
 *        using the gradient of the pressure correction field `Phi`.
 *
 * @param user The UserCtx containing the Ucont and Phi vectors.
 * @return PetscErrorCode 0 on success.
 *
 * @note Testing status:
 *       Direct unit coverage exists for basic projection invariants; periodic
 *       and immersed-boundary correction branches remain part of the next-gap backlog.
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
