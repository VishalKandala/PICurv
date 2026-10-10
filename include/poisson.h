#ifndef POISSON_H
#define POISSON_H

#include "variables.h"
#include "Boundaries.h"

/**
 * @file poisson.h
 * @brief Pressure-Poisson projection: the multigrid pressure-correction solve, the
 *        pressure update, and the velocity correction that makes `Ucont` divergence-free.
 *
 * The discrete operator, the right-hand side, and the projection gradient are built from
 * one face-gradient definition, so the Laplacian the solver inverts is exactly the
 * divergence of the gradient the projection applies.
 *
 * Every masked operation reads the solid field from the `UserCtx` of the level it acts
 * on (`lNvert`). Solid cells are excluded from the operator, the right-hand side, the
 * projection, and the multigrid transfers; cells next to a solid face use one-sided
 * transverse differences.
 */

/**
 * @brief Solves the pressure-correction equation for every block with geometric multigrid.
 *
 * On the first call for a block, assembles the operator on every multigrid level and
 * builds the outer Krylov solver with its `PCMG` preconditioner, grid transfers, level
 * smoothers, coarse solve, and null space. The solver is kept in the finest level's
 * `UserCtx::ksp` and reused on every later call; it is destroyed with the context.
 *
 * Each call refreshes the ghosted contravariant flux, forms the right-hand side, solves
 * for `Phi` on the finest level, and appends the iteration history to
 * `Poisson_Solver_Convergence_History_Block_<bi>.log` in the run's log directory.
 *
 * Solver controls are read from PETSc options under the `ps_` prefix
 * (`-ps_ksp_*`, `-ps_mg_levels_N_*`, `-ps_mg_coarse_*`) when the solver is built.
 *
 * @param[in,out] usermg Multigrid hierarchy; the finest level's `Phi` receives the solution.
 * @return PETSc error code. A non-finite residual or a preconditioner that could not be
 *         built stops the run with PETSC_ERR_NOT_CONVERGED.
 */
extern PetscErrorCode PoissonSolver_Multigrid(UserMG *usermg);

/**
 * @brief Assembles the pressure-correction operator on one multigrid level.
 *
 * Allocates `user->A` on the first call and reassembles its entries on later calls into
 * the same nonzero structure. Fluid rows carry the 19-point curvilinear Laplacian;
 * dummy rows on the domain boundary and solid rows are identities. Non-periodic faces
 * are homogeneous Neumann; periodic faces wrap to the opposite interior layer.
 *
 * @param[in,out] user Level context supplying metrics, `lNvert`, and boundary types.
 * @return PETSc error code.
 */
extern PetscErrorCode AssemblePoissonOperator(UserCtx *user);

/**
 * @brief Forms the right-hand side of the pressure-correction equation.
 *
 * Writes the scaled divergence of the ghosted contravariant flux `lUcont` into @p B,
 * zero on dummy and solid cells, and stores its domain integral in
 * `SimCtx::poissonSourceImbalance`.
 *
 * @param[in]  user Finest-level context.
 * @param[out] B    Right-hand-side vector on `user->da`.
 * @return PETSc error code.
 */
extern PetscErrorCode ComputePoissonRHS(UserCtx *user, Vec B);

/**
 * @brief Adds the pressure correction to the pressure, `P += Phi`, and refreshes both
 *        fields' periodic images and ghosts.
 * @param[in,out] user Block context holding `P` and `Phi`.
 * @return PETSc error code.
 */
extern PetscErrorCode UpdatePressure(UserCtx *user);

/**
 * @brief Corrects the contravariant flux with the gradient of `Phi`.
 *
 * Subtracts `dt / COEF_TIME_ACCURACY` times the face pressure-gradient flux from every
 * fluid face of `Ucont`, then refreshes the periodic images, reconstructs the Cartesian
 * velocity, and finalizes the cell fields that depend on it.
 *
 * @param[in,out] user Block context; reads the ghosted `lPhi`.
 * @return PETSc error code.
 */
extern PetscErrorCode ProjectVelocity(UserCtx *user);

#endif /* POISSON_H */
