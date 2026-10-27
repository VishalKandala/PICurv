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

/* ===========================================================================================
 *  Face-gradient operator, shared with the momentum residual's pressure gradient.
 *  poisson.c documents the face gradient; the per-face primitives are inline so the
 *  momentum residual's hot loop does not pay a call per face.
 * =========================================================================================== */

/** Cells whose `nvert` exceeds this value are solid. */
#define POISSON_SOLID_THRESHOLD 0.1

/**
 * @brief Transverse difference used at one face.
 *
 * The difference is `weight * (sum of Phi on row hi - sum of Phi on row lo)`, where a row
 * is the pair of cells on either side of the face, displaced by the given offset along the
 * transverse axis.
 */
typedef struct {
    PetscInt  lo;      /**< Transverse offset of the subtracted row pair. */
    PetscInt  hi;      /**< Transverse offset of the added row pair. */
    PetscReal weight;  /**< 0.25 central, 0.5 one-sided, 0 when no fluid side remains. */
} PoissonTransverseDifference;

/** @brief Everything needed to evaluate the gradient flux through one face. */
typedef struct {
    PetscReal                   dot[3];   /**< F_b . F_n on the face, b = 0..2. */
    PetscReal                   aj;       /**< Inverse Jacobian on the face. */
    PoissonTransverseDifference diff[3];  /**< Transverse differences; diff[n] is unused. */
} PoissonFaceGradient;

/** @brief Read-only face metric arrays: `metric[n][b]` is F_b on the n-faces. */
typedef struct {
    const Cmpnts    ***metric[3][3];
    const PetscReal ***aj[3];
} PoissonFaceMetrics;

/**
 * @brief Value of a cell-centred array at cell @p c displaced by @p d.
 *
 * @param[in] field Ghosted cell-centred array, indexed [k][j][i].
 * @param[in] c     Cell (i, j, k).
 * @param[in] d     Displacement (di, dj, dk).
 * @return The value at `c + d`.
 */
static inline PetscReal PoissonOperator_At(const PetscReal ***field, const PetscInt c[3],
                                          const PetscInt d[3])
{
    return field[c[2] + d[2]][c[1] + d[1]][c[0] + d[0]];
}

/**
 * @brief Chooses the transverse difference along axis @p t at the face between cell
 *        @p c and `c + e_n`.
 *
 * The central difference averages the two cells beside the face over rows -1 and +1.
 * When the +1 row is a non-periodic boundary layer or touches a solid cell, the difference
 * falls back to rows -1 and 0, and symmetrically to rows 0 and +1; with neither side
 * available the transverse term vanishes.
 *
 * @param[in] nvert    Ghosted solid field.
 * @param[in] c        Cell on the low side of the face.
 * @param[in] n        Face orientation (0, 1, 2 for i, j, k).
 * @param[in] t        Transverse axis, not @p n.
 * @param[in] m        Grid sizes (mx, my, mz).
 * @param[in] periodic Whether each axis is periodic.
 * @return The rows and weight of the difference.
 */
static inline PoissonTransverseDifference PoissonOperator_TransverseDifference(
    const PetscReal ***nvert, const PetscInt c[3], PetscInt n, PetscInt t,
    const PetscInt m[3], const PetscBool periodic[3])
{
    PoissonTransverseDifference diff = {0, 0, 0.0};
    PetscInt plus[3] = {0, 0, 0}, plus_n[3] = {0, 0, 0};
    PetscInt minus[3] = {0, 0, 0}, minus_n[3] = {0, 0, 0};
    const PetscInt s = c[t];

    plus[t] = 1;  plus_n[t] = 1;  plus_n[n] = 1;
    minus[t] = -1; minus_n[t] = -1; minus_n[n] = 1;
    const PetscReal solid_plus  = PoissonOperator_At(nvert, c, plus)  + PoissonOperator_At(nvert, c, plus_n);
    const PetscReal solid_minus = PoissonOperator_At(nvert, c, minus) + PoissonOperator_At(nvert, c, minus_n);

    if ((s == m[t] - 2 && !periodic[t]) || solid_plus > POISSON_SOLID_THRESHOLD) {
        if (solid_minus < POISSON_SOLID_THRESHOLD && (s != 1 || periodic[t])) {
            diff.lo = -1; diff.hi = 0; diff.weight = 0.5;
        }
    } else if ((s == 1 && !periodic[t]) || solid_minus > POISSON_SOLID_THRESHOLD) {
        if (solid_plus < POISSON_SOLID_THRESHOLD) {
            diff.lo = 0; diff.hi = 1; diff.weight = 0.5;
        }
    } else {
        diff.lo = -1; diff.hi = 1; diff.weight = 0.25;
    }
    return diff;
}

/**
 * @brief Gathers the metric coefficients and transverse differences of the gradient flux
 *        through the face between cell @p c and `c + e_n`.
 *
 * @param[in] metrics  Face metric arrays from PoissonOperator_GetFaceMetrics().
 * @param[in] nvert    Ghosted solid field.
 * @param[in] c        Cell on the low side of the face.
 * @param[in] n        Face orientation (0, 1, 2 for i, j, k).
 * @param[in] m        Grid sizes (mx, my, mz).
 * @param[in] periodic Whether each axis is periodic.
 * @return The face's metric products, inverse Jacobian and transverse differences.
 */
static inline PoissonFaceGradient PoissonOperator_FaceGradientStencil(
    const PoissonFaceMetrics *metrics, const PetscReal ***nvert, const PetscInt c[3],
    PetscInt n, const PetscInt m[3], const PetscBool periodic[3])
{
    PoissonFaceGradient face;
    const Cmpnts normal = metrics->metric[n][n][c[2]][c[1]][c[0]];

    face.aj = metrics->aj[n][c[2]][c[1]][c[0]];
    for (PetscInt b = 0; b < 3; b++) {
        const Cmpnts base = metrics->metric[n][b][c[2]][c[1]][c[0]];
        face.dot[b] = base.x * normal.x + base.y * normal.y + base.z * normal.z;
        if (b == n) {
            face.diff[b].lo = 0; face.diff[b].hi = 0; face.diff[b].weight = 0.0;
        } else {
            face.diff[b] = PoissonOperator_TransverseDifference(nvert, c, n, b, m, periodic);
        }
    }
    return face;
}

/**
 * @brief Gradient flux of a cell-centred @p field through the face between cell @p c and
 *        `c + e_n`: `sum_b g_nb D_b(field)` of the face gradient defined in poisson.c.
 *
 * The projection applies it to `Phi`; the momentum residual applies it to `P`, so both use
 * the one discrete gradient.
 *
 * @param[in] metrics  Face metric arrays from PoissonOperator_GetFaceMetrics().
 * @param[in] nvert    Ghosted solid field.
 * @param[in] field    Ghosted cell-centred field whose gradient is taken.
 * @param[in] c        Cell on the low side of the face.
 * @param[in] n        Face orientation (0, 1, 2 for i, j, k).
 * @param[in] m        Grid sizes (mx, my, mz).
 * @param[in] periodic Whether each axis is periodic.
 * @return The gradient flux through the face.
 */
static inline PetscReal PoissonOperator_FaceGradientFlux(
    const PoissonFaceMetrics *metrics, const PetscReal ***nvert, const PetscReal ***field,
    const PetscInt c[3], PetscInt n, const PetscInt m[3], const PetscBool periodic[3])
{
    const PoissonFaceGradient face = PoissonOperator_FaceGradientStencil(metrics, nvert, c, n, m, periodic);
    PetscInt  across[3] = {0, 0, 0};
    PetscReal difference[3];

    across[n] = 1;
    for (PetscInt b = 0; b < 3; b++) {
        if (b == n) {
            difference[b] = PoissonOperator_At(field, c, across) - field[c[2]][c[1]][c[0]];
            continue;
        }
        const PoissonTransverseDifference diff = face.diff[b];
        PetscInt hi_lower[3] = {0, 0, 0}, hi_upper[3] = {0, 0, 0};
        PetscInt lo_lower[3] = {0, 0, 0}, lo_upper[3] = {0, 0, 0};
        hi_lower[b] = diff.hi; hi_upper[b] = diff.hi; hi_upper[n] = 1;
        lo_lower[b] = diff.lo; lo_upper[b] = diff.lo; lo_upper[n] = 1;
        difference[b] = (diff.weight == 0.0) ? 0.0 :
            (PoissonOperator_At(field, c, hi_lower) + PoissonOperator_At(field, c, hi_upper) -
             PoissonOperator_At(field, c, lo_lower) - PoissonOperator_At(field, c, lo_upper)) * diff.weight;
    }
    return difference[0] * face.dot[0] * face.aj +
           difference[1] * face.dot[1] * face.aj +
           difference[2] * face.dot[2] * face.aj;
}

/**
 * @brief Records which axes are periodic, from the negative face of each axis.
 *
 * @param[in]  user     Block context.
 * @param[out] periodic Whether each axis is periodic.
 */
void PoissonOperator_PeriodicAxes(const UserCtx *user, PetscBool periodic[3]);

/**
 * @brief Borrows read access to the face metric arrays of @p user.
 *
 * @param[in]  user    Block context whose ghosted face metrics are read.
 * @param[out] metrics Face metric arrays; return them with PoissonOperator_RestoreFaceMetrics().
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode PoissonOperator_GetFaceMetrics(UserCtx *user, PoissonFaceMetrics *metrics);

/**
 * @brief Returns the arrays borrowed by PoissonOperator_GetFaceMetrics().
 *
 * @param[in]     user    Block context the arrays were borrowed from.
 * @param[in,out] metrics Face metric arrays to return.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode PoissonOperator_RestoreFaceMetrics(UserCtx *user, PoissonFaceMetrics *metrics);

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
