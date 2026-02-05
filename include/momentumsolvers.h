#ifndef MOMENTUNSOLVERS_H
#define MOMENTUMSOLVERS_H

#include "variables.h" // Provides definitions for UserCtx, SimCtx, IBMNodes, etc.
#include "logging.h"
#include "rhs.h"
#include "Boundaries.h"

/*================================================================================*
 *                        MOMENTUM EQUATION SOLVERS                               *
 *================================================================================*/

/**
 * @brief Advances the momentum equations using an explicit 4th-order Runge-Kutta scheme.
 * @param user Array of UserCtx structs for all blocks.
 * @param ibm  (Optional) Pointer to IBM data. Pass NULL if disabled.
 * @param fsi  (Optional) Pointer to FSI data. Pass NULL if disabled.
 * @return PetscErrorCode 0 on success.
 */
extern PetscErrorCode MomentumSolver_Explicit_RungeKutta4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

/**
 * @brief Solves the Momentum Equations using Dual-Time Stepping with a Fixed-Point RK4 Smoother.
 *
 * =================================================================================================
 * GLOSSARY & THEORETICAL BASIS
 * =================================================================================================
 * 
 * 1. METHODOLOGY: Dual-Time Stepping (Pseudo-Time Integration)
 *    We aim to solve the implicit BDF equation:  R_spatial(U) + dU/dt_physical = 0.
 *    We do this by introducing a fictitious "Pseudo-Time" (tau) and iterating to steady state:
 *       dU/d(tau) = - [ R_spatial(U) + BDF_Terms(U) ]
 *    When dU/d(tau) -> 0, the physical time step is satisfied.
 *
 * 2. ALGORITHM: Fixed-Point Iteration with Explicit Runge-Kutta
 *    This is technically a Fixed-Point iteration on the operator:
 *       U_new = U_old + pseudo_dt_scaling * dt_pseudo * Total_Residual(U_old)
 *    We use a 4-Stage Explicit RK scheme (Jameson-Schmidt-Turkel coeffs) to smooth errors.
 *
 * 3. STABILITY: Backtracking Line Search
 *    If a pseudo-time step causes the Residual or Solution Error to GROW (Divergence),
 *    the solver "Backtracks": it restores the previous solution, cuts the pseudo-time step
 *    scaling factor (pseudo_dt_scaling) in half, and retries the iteration.
 *
 * =================================================================================================
 * VARIABLE MAPPING
 * =================================================================================================
 * -- Physics Variables (Legacy Names Kept) --
 * ti   : Physical Time Step Index.
 * dt   : Physical Time Step size (Delta t).
 * st   : Pseudo-Time Step size (Delta tau).
 * alfa : Runge-Kutta stage coefficients {1/4, 1/3, 1/2, 1}.
 *
 * -- Convergence & Solver Control (Renamed) --
 * pseudo_iter       : Counter for the inner dual-time loop.
 * pseudo_dt_scaling : Adaptive scalar for the pseudo-time step (formerly lambda).
 * delta_sol_norm    : The L_inf norm of the change in solution (dU).
 * resid_norm        : The L_inf norm of the Total Residual (RHS).
 * =================================================================================================
 */
PetscErrorCode MomentumSolver_DualTime_Picard_RK4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

#endif // MOMENTUMSOLVERS_H
