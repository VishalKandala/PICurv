#ifndef MOMENTUMSOLVERS_H
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
 *
 * @note Testing status:
 *       The explicit RK path remains on the near-term backlog for direct
 *       positive-path bespoke coverage; today it is weaker than the dual-time
 *       path in the test surface.
 */
extern PetscErrorCode MomentumSolver_Explicit_RungeKutta4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

/**
 * @brief Solves one physical momentum step with matrix-free Newton--Krylov.
 *
 * Version one uses a finite-difference matrix-free Jacobian, GMRES, and no
 * preconditioner. All PETSc solver objects are local to this call. Rows removed
 * by legacy boundary residual enforcement are made explicit: conditioned normal
 * rows use X-Uconditioned, untouched dummy/tangential rows use X, and periodic
 * duplicates use Xdup-Xrepresentative. Unsupported masked, interface, and
 * component-disabled rows are rejected before setup.
 *
 * @param user Single-block momentum context.
 * @param ibm  Must be NULL; immersed boundaries are not supported in version one.
 * @param fsi  Must be NULL; moving-body coupling is not supported in version one.
 * @return 0 on convergence, PETSC_ERR_CONV_FAILED after rollback on nonconvergence.
 */
PetscErrorCode MomentumSolver_NewtonKrylov(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

/**
 * @brief Solves the momentum equations using dual-time Picard iteration with Jameson RK smoothing.
 *
 * =================================================================================================
 * GLOSSARY & THEORETICAL BASIS
 * =================================================================================================
 * 1. METHODOLOGY: Dual-Time Stepping (Pseudo-Time Integration)
 * We aim to solve the implicit BDF equation:  R_spatial(U) + dU/dt_physical = 0.
 * We do this by introducing a fictitious "Pseudo-Time" (tau) and iterating to steady state:
 * dU/d(tau) = - [ R_spatial(U) + BDF_Terms(U) ]
 * When dU/d(tau) -> 0, the physical time step is satisfied.
 * 2. ALGORITHM: Fixed-Point Iteration with Explicit Runge-Kutta
 * This is technically a Fixed-Point iteration on the operator:
 * U_new = U_old + pseudo_dtau * alfa_stage * Total_Residual(U_old)
 * where pseudo_dtau = pseudo_cfl / lambda_max is the spectral-radius-based pseudo-time step
 * (lambda_max = global max convective spectral radius). This makes pseudo_cfl a true
 * dimensionless Courant number, independent of the physical time step dt.
 * We use a 4-Stage Explicit RK scheme (Jameson-Schmidt-Turkel coeffs) to smooth errors.
 * 3. STABILITY: Adaptive Pseudo-CFL Trial Acceptance and Rollback
 * If a pseudo-time trial causes excessive residual growth, the solver restores the
 * previous accepted state, reduces the global pseudo-CFL, and retries.
 * =================================================================================================
 * VARIABLE MAPPING
 * =================================================================================================
 * -- Physics Variables (Legacy Names Kept) --
 * ti   : Physical Time Step Index.
 * dt   : Physical Time Step size (Delta t).
 * st   : Pseudo-Time Step size (Delta tau).
 * alfa : Runge-Kutta stage coefficients {1/4, 1/3, 1/2, 1}.
 * -- Convergence & Solver Control (Renamed) --
 * pseudo_iter       : Counter for the inner dual-time loop.
 * pseudo_dtau       : Adaptive pseudo-time step [physical time], = pseudo_cfl / lambda_max.
 * lambda_max        : Global max convective spectral radius [1/s] from the current field.
 * delta_sol_norm    : The L_inf norm of the change in solution (dU).
 * resid_norm        : The L_inf norm of the Total Residual (RHS).
 * =================================================================================================
 *
 * @param user Primary `UserCtx` input for the operation.
 * @param ibm Parameter `ibm` passed to `MomentumSolver_DualTime_Picard_JamesonRK()`.
 * @param fsi Parameter `fsi` passed to `MomentumSolver_DualTime_Picard_JamesonRK()`.
 * @return PetscErrorCode 0 on success.
 *
 * @note Testing status:
 *       This solver is covered primarily through runtime smoke and orchestration
 *       tests. A smaller direct invariant-style positive-path harness remains
 *       part of the next-gap backlog.
 */
PetscErrorCode MomentumSolver_DualTime_Picard_JamesonRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi);

/** @deprecated Use MomentumSolver_DualTime_Picard_JamesonRK(). */
#define MomentumSolver_DualTime_Picard_RK4 MomentumSolver_DualTime_Picard_JamesonRK

/*================================================================================*
 *               SHARED PHYSICAL-TIME (BDF) COEFFICIENT PLUMBING                  *
 *================================================================================*/

/**
 * @brief Returns whether the current physical step uses the BDF2 discretization.
 *
 * Single source of truth for the BDF1/BDF2 selection. The predicate is identical
 * to the one historically inlined in ComputeTotalResidual():
 *   BDF2 when COEF_TIME_ACCURACY > 1.1 AND step != StartStep AND step != 1,
 *   otherwise BDF1 (startup step and the first step after a restart).
 *
 * @param simCtx Master simulation context (reads step, StartStep).
 * @return PETSC_TRUE for BDF2, PETSC_FALSE for BDF1.
 */
PetscBool MomentumUsesBDF2(SimCtx *simCtx);

/**
 * @brief Returns the BDF physical-time coefficient a0 for the current step.
 *
 * a0 = 1.5 (== COEF_TIME_ACCURACY) for BDF2, a0 = 1.0 for BDF1. Used both as the
 * leading coefficient of the physical-time term in the residual and as the
 * additive temporal contribution lambda_t = a0/dt in the momentum stability
 * estimate, keeping the two numerically consistent.
 *
 * @param simCtx Master simulation context.
 * @return a0 in {1.0, 1.5}.
 */
PetscReal MomentumBDFCoefficient(SimCtx *simCtx);

/**
 * @brief Computes the shared spatial-plus-BDF momentum residual in user->Rhs.
 * @param user Block context with an allocated Rhs vector.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeTotalResidual(UserCtx *user);

/*================================================================================*
 *                 MOMENTUM PSEUDO-TIME STABILITY ESTIMATE (SHADOW)               *
 *================================================================================*/

/**
 * @brief Convective-estimate candidate selector. See ComputeMomentumStabilityEstimate().
 *
 * B: six-face transport scale (f_c * Aj * sum|U_f| / 2).
 * C: B + frozen-advector discrete-divergence diagonal term.
 * D: C + nonlinear velocity-gradient row-norm term (lambda_grad_u).
 */
typedef enum {
    MOM_STAB_CAND_B = 0,
    MOM_STAB_CAND_C = 1,
    MOM_STAB_CAND_D = 2
} MomStabCandidate;

/**
 * @brief Dominant stiffness contributor at the controlling cell.
 */
typedef enum {
    MOM_STAB_LIMITER_TIME       = 0,
    MOM_STAB_LIMITER_CONVECTION = 1,
    MOM_STAB_LIMITER_VISCOSITY  = 2
} MomStabLimiter;

/**
 * @brief Diagnostic report produced by ComputeMomentumStabilityEstimate().
 *
 * This is a PRACTICAL CONSERVATIVE STABILITY ESTIMATE (operator-scaled pseudo-time
 * estimate), not a proven spectral radius. All lambda_* are global maxima in [1/s].
 */
typedef struct {
    PetscReal lambda;        /* selected-candidate global max estimate [1/s] */
    PetscReal lambda_t;      /* temporal term a0/dt (uniform across cells)   */
    PetscReal lambda_c;      /* convective part at the controlling cell      */
    PetscReal lambda_v;      /* viscous part at the controlling cell         */
    PetscReal lambda_B;      /* global max of (lambda_t + lambda_c^B + lambda_v) */
    PetscReal lambda_C;      /* global max with candidate C convective term  */
    PetscReal lambda_D;      /* global max with candidate D convective term  */
    PetscInt  ci, cj, ck;    /* controlling-cell global index (selected cand) */
    PetscInt  cblock;        /* controlling-cell block                        */
    PetscInt  cclass;        /* 0=interior, 1=physical-boundary, 2=IB-adjacent */
    PetscInt  one_sided;     /* controlling cell used the one-sided viscous x2 */
    PetscInt  active_cells;  /* global count of active (non-masked) cells       */
    PetscBool estimate_incomplete; /* true if Clark/RANS/vel-dependent force is active (uncovered) */
    MomStabLimiter limiter;  /* dominant contributor at the controlling cell  */
} MomStabilityReport;

/**
 * @brief Compute the momentum pseudo-time stability estimate (shadow/diagnostic).
 *
 * Conservative, operator-scaled estimate: lambda = max_cell (a0/dt + lambda_c + lambda_nu),
 * over active, non-solid cells, blocks, and MPI ranks, where lambda_c already includes the
 * per-direction QUICK scheme factors. Read-only; performs no halo exchange, but does perform
 * global scalar collectives (see implementation). This is a PRACTICAL CONSERVATIVE estimate,
 * not a proven spectral radius. See the Workstream-A design.
 *
 * Call-site preconditions (NOT enforced internally): lUcont, lUcat, lNu_t, lNvert
 * must be fresh; lAj, face Jacobians and face metrics are static after grid init.
 *
 * @param[in]  user         Array of UserCtx (one per block).
 * @param[in]  block_number Number of blocks.
 * @param[in]  dt           Physical time step.
 * @param[in]  candidate    Convective candidate driving rep->lambda (B, C or D).
 * @param[out] rep          Filled diagnostic report (global maxima + breakdown).
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeMomentumStabilityEstimate(UserCtx *user, PetscInt block_number,
                                                PetscReal dt, MomStabCandidate candidate,
                                                MomStabilityReport *rep);

/**
 * @brief Active staggered-momentum row mask for a cell (exposed for unit testing).
 *
 * Returns a 3-bit mask (xi=1, eta=2, zeta=4). A solid cell yields 0 (all inactive); a
 * positive solid neighbour or a positive non-periodic physical face clears the corresponding
 * normal row; TwoD (1/2/3) clears the homogeneous direction's row.
 *
 * @param nvert Local nvert array (ghosted).
 * @param k Cell k index.
 * @param j Cell j index.
 * @param i Cell i index.
 * @param mx Global x dimension.
 * @param my Global y dimension.
 * @param mz Global z dimension.
 * @param np_x1 True if the positive-x face is non-periodic.
 * @param np_y1 True if the positive-y face is non-periodic.
 * @param np_z1 True if the positive-z face is non-periodic.
 * @param twoD TwoD homogeneous-direction selector (0 none, 1 xi, 2 eta, 3 zeta).
 * @return 3-bit active-row mask (0 when the location carries no active unknown).
 */
PetscInt MomCellActiveRows(PetscReal ***nvert, PetscInt k, PetscInt j, PetscInt i,
                           PetscInt mx, PetscInt my, PetscInt mz,
                           PetscBool np_x1, PetscBool np_y1, PetscBool np_z1, PetscInt twoD);

#endif // MOMENTUMSOLVERS_H
