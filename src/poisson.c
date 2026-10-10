/**
 * @file poisson.c
 * @brief Pressure-Poisson projection: operator, right-hand side, multigrid solve,
 *        pressure update, and velocity correction.
 *
 * Grid layout. Cell-centred quantities live at DMDA indices 1..m-2 on each axis; indices 0
 * and m-1 are dummy layers that carry no unknown. On a periodic axis the interior wraps
 * from m-2 to 1. The contravariant flux `Ucont[k][j][i].x` sits on the face between cells
 * i and i+1, and likewise for the other two components.
 *
 * Face gradient. The flux of grad(Phi) through the face between cell c and c + e_n is
 *
 *     sum_b g_nb D_b(Phi),   g_nb = (F_b . F_n) * aj_n,
 *
 * where F_b are the contravariant base vectors stored on the n-faces, aj_n is the inverse
 * Jacobian there, D_n is the difference across the face, and D_b (b != n) is the transverse
 * difference chosen by PoissonOperator_TransverseDifference(). The operator assembles the
 * divergence of this flux, and the projection subtracts it from `Ucont`, so the two are
 * consistent by construction.
 */
#include "poisson.h"
#include "logging.h"
#include "setup.h"

/** Cells whose `nvert` exceeds this value are solid. */
#define POISSON_SOLID_THRESHOLD 0.1

/** Offsets of the 19 stencil points, in the column order used to insert each row. */
static const PetscInt POISSON_STENCIL_OFFSETS[19][3] = {
    { 0,  0,  0},                                              /* centre */
    { 1,  0,  0}, {-1,  0,  0}, { 0,  1,  0}, { 0, -1,  0},    /* faces  */
    { 0,  0,  1}, { 0,  0, -1},
    { 1,  1,  0}, { 1, -1,  0}, {-1,  1,  0}, {-1, -1,  0},    /* edges  */
    { 0,  1,  1}, { 0,  1, -1}, { 0, -1,  1}, { 0, -1, -1},
    { 1,  0,  1}, { 1,  0, -1}, {-1,  0,  1}, {-1,  0, -1},
};

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

/** @brief Field IDs of the face metric vectors, indexed as PoissonFaceMetrics. */
static const FieldId POISSON_FACE_METRIC_FIELDS[3][3] = {
    {FIELD_ID_ICSI, FIELD_ID_IETA, FIELD_ID_IZET},
    {FIELD_ID_JCSI, FIELD_ID_JETA, FIELD_ID_JZET},
    {FIELD_ID_KCSI, FIELD_ID_KETA, FIELD_ID_KZET},
};
static const FieldId POISSON_FACE_AJ_FIELDS[3] = {FIELD_ID_IAJ, FIELD_ID_JAJ, FIELD_ID_KAJ};

/**
 * @brief Maps the 19 stencil offsets to their slots; corners of the 3x3x3 block are -1.
 */
static PetscInt PoissonOperator_StencilSlot(const PetscInt d[3])
{
    for (PetscInt s = 0; s < 19; s++) {
        if (POISSON_STENCIL_OFFSETS[s][0] == d[0] && POISSON_STENCIL_OFFSETS[s][1] == d[1] &&
            POISSON_STENCIL_OFFSETS[s][2] == d[2]) return s;
    }
    return -1;
}

/** @brief Value of a cell-centred array at cell @p c displaced by @p d. */
static inline PetscReal PoissonOperator_At(const PetscReal ***field, const PetscInt c[3],
                                          const PetscInt d[3])
{
    return field[c[2] + d[2]][c[1] + d[1]][c[0] + d[0]];
}

/** @brief Records which axes are periodic, from the negative face of each axis. */
static void PoissonOperator_PeriodicAxes(const UserCtx *user, PetscBool periodic[3])
{
    periodic[0] = (PetscBool)(user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC);
    periodic[1] = (PetscBool)(user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC);
    periodic[2] = (PetscBool)(user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC);
}

/**
 * @brief Chooses the transverse difference along axis @p t at the face between cell
 *        @p c and `c + e_n`.
 *
 * The central difference averages the two cells beside the face over rows -1 and +1.
 * When the +1 row is a non-periodic boundary layer or touches a solid cell, the difference
 * falls back to rows -1 and 0, and symmetrically to rows 0 and +1; with neither side
 * available the transverse term vanishes.
 */
static PoissonTransverseDifference PoissonOperator_TransverseDifference(
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
 */
static PoissonFaceGradient PoissonOperator_FaceGradientStencil(
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

/** @brief Borrows read access to the face metric arrays of @p user. */
static PetscErrorCode PoissonOperator_GetFaceMetrics(UserCtx *user, PoissonFaceMetrics *metrics)
{
    FieldView view;

    PetscFunctionBeginUser;
    for (PetscInt n = 0; n < 3; n++) {
        for (PetscInt b = 0; b < 3; b++) {
            PetscCall(FieldGetView(user, POISSON_FACE_METRIC_FIELDS[n][b], &view));
            PetscCall(DMDAVecGetArrayRead(view.dm, view.local_vec, (void *)&metrics->metric[n][b]));
        }
        PetscCall(FieldGetView(user, POISSON_FACE_AJ_FIELDS[n], &view));
        PetscCall(DMDAVecGetArrayRead(view.dm, view.local_vec, (void *)&metrics->aj[n]));
    }
    PetscFunctionReturn(0);
}

/** @brief Returns the arrays borrowed by PoissonOperator_GetFaceMetrics(). */
static PetscErrorCode PoissonOperator_RestoreFaceMetrics(UserCtx *user, PoissonFaceMetrics *metrics)
{
    FieldView view;

    PetscFunctionBeginUser;
    for (PetscInt n = 0; n < 3; n++) {
        for (PetscInt b = 0; b < 3; b++) {
            PetscCall(FieldGetView(user, POISSON_FACE_METRIC_FIELDS[n][b], &view));
            PetscCall(DMDAVecRestoreArrayRead(view.dm, view.local_vec, (void *)&metrics->metric[n][b]));
        }
        PetscCall(FieldGetView(user, POISSON_FACE_AJ_FIELDS[n], &view));
        PetscCall(DMDAVecRestoreArrayRead(view.dm, view.local_vec, (void *)&metrics->aj[n]));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Index of the neighbour at offset @p d (-1, 0, +1) from @p v on an axis of @p m
 *        points, wrapping between the interior layers 1 and m-2 when periodic.
 */
static inline PetscInt PoissonOperator_NeighborIndex(PetscInt v, PetscInt d, PetscInt m, PetscBool periodic)
{
    if (periodic && d == 1 && v == m - 2) return 1;
    if (periodic && d == -1 && v == 1) return m - 2;
    return v + d;
}

/**
 * @brief Adds the signed gradient flux through one face to a row's stencil coefficients.
 *
 * @param[in]     face  Face gradient from PoissonOperator_FaceGradientStencil().
 * @param[in]     own   Offset of the face's lower cell from the row cell.
 * @param[in]     n     Normal axis of the face.
 * @param[in]     sign  +1 for the row's upper face on the axis, -1 for its lower face.
 * @param[in,out] coefficients The row's 19 coefficients.
 */
static void PoissonOperator_AddFaceFlux(const PoissonFaceGradient *face, const PetscInt own[3],
                                        PetscInt n, PetscReal sign, PetscScalar coefficients[19])
{
    for (PetscInt b = 0; b < 3; b++) {
        const PetscReal g = face->dot[b] * face->aj;
        PetscInt lower[3] = {own[0], own[1], own[2]};
        PetscInt upper[3] = {own[0], own[1], own[2]};

        upper[n] += 1;
        if (b == n) {
            coefficients[PoissonOperator_StencilSlot(lower)] += sign * (-g);
            coefficients[PoissonOperator_StencilSlot(upper)] += sign * g;
            continue;
        }

        const PoissonTransverseDifference diff = face->diff[b];
        if (diff.weight == 0.0) continue;
        const PetscReal term = g * diff.weight;
        PetscInt hi_lower[3] = {lower[0], lower[1], lower[2]}, hi_upper[3] = {upper[0], upper[1], upper[2]};
        PetscInt lo_lower[3] = {lower[0], lower[1], lower[2]}, lo_upper[3] = {upper[0], upper[1], upper[2]};
        hi_lower[b] += diff.hi; hi_upper[b] += diff.hi;
        lo_lower[b] += diff.lo; lo_upper[b] += diff.lo;
        coefficients[PoissonOperator_StencilSlot(hi_lower)] += sign * term;
        coefficients[PoissonOperator_StencilSlot(hi_upper)] += sign * term;
        coefficients[PoissonOperator_StencilSlot(lo_lower)] += sign * (-term);
        coefficients[PoissonOperator_StencilSlot(lo_upper)] += sign * (-term);
    }
}

#undef __FUNCT__
#define __FUNCT__ "AssemblePoissonOperator"
/**
 * @brief Implementation of \ref AssemblePoissonOperator().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/poisson.h`.
 * @see AssemblePoissonOperator()
 */
PetscErrorCode AssemblePoissonOperator(UserCtx *user)
{
    const DMDALocalInfo info = user->info;
    const PetscInt      m[3] = {info.mx, info.my, info.mz};
    PetscBool           periodic[3];
    PoissonFaceMetrics  metrics;
    const PetscReal  ***nvert, ***aj;
    AO                  ao;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!user->A) {
        PetscInt local_rows;
        PetscCall(VecGetLocalSize(user->Phi, &local_rows));
        PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, local_rows, local_rows, m[0] * m[1] * m[2], m[0] * m[1] * m[2],
                               19, NULL, 19, NULL, &user->A));
    }
    PetscCall(MatZeroEntries(user->A));

    PoissonOperator_PeriodicAxes(user, periodic);
    PetscCall(DMDAGetAO(user->da, &ao));
    PetscCall(PoissonOperator_GetFaceMetrics(user, &metrics));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));

    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                const PetscInt c[3] = {i, j, k};
                PetscInt       row = i + j * m[0] + k * m[0] * m[1];
                PetscInt       columns[19];
                PetscScalar    coefficients[19] = {0.0};

                PetscCall(AOApplicationToPetsc(ao, 1, &row));
                if (i == 0 || i == m[0] - 1 || j == 0 || j == m[1] - 1 || k == 0 || k == m[2] - 1) {
                    const PetscScalar one = 1.0;
                    PetscCall(MatSetValues(user->A, 1, &row, 1, &row, &one, INSERT_VALUES));
                    continue;
                }

                for (PetscInt s = 0; s < 19; s++) {
                    const PetscInt *d = POISSON_STENCIL_OFFSETS[s];
                    columns[s] = PoissonOperator_NeighborIndex(i, d[0], m[0], periodic[0]) +
                                 PoissonOperator_NeighborIndex(j, d[1], m[1], periodic[1]) * m[0] +
                                 PoissonOperator_NeighborIndex(k, d[2], m[2], periodic[2]) * m[0] * m[1];
                }
                PetscCall(AOApplicationToPetsc(ao, 19, columns));

                if (nvert[k][j][i] > POISSON_SOLID_THRESHOLD) {
                    /* Solid rows keep the fluid row's structure with zero couplings, so a
                       solid field that changes during a run reassembles in place. */
                    coefficients[0] = 1.0;
                    PetscCall(MatSetValues(user->A, 1, &row, 19, columns, coefficients, INSERT_VALUES));
                    continue;
                }

                /* Faces in the order east, west, north, south, top, bottom. Non-periodic
                   boundary faces carry no flux (homogeneous Neumann). */
                for (PetscInt n = 0; n < 3; n++) {
                    const PetscInt first = periodic[n] ? 0 : 1;
                    const PetscInt last  = periodic[n] ? m[n] - 1 : m[n] - 2;
                    for (PetscInt side = 0; side < 2; side++) {
                        PetscInt across[3] = {0, 0, 0}, own[3] = {0, 0, 0}, face_cell[3] = {i, j, k};
                        const PetscBool upper = (PetscBool)(side == 0);

                        across[n] = upper ? 1 : -1;
                        if (PoissonOperator_At(nvert, c, across) >= POISSON_SOLID_THRESHOLD) continue;
                        if (c[n] == (upper ? last : first)) continue;
                        if (!upper) { own[n] = -1; face_cell[n] -= 1; }

                        const PoissonFaceGradient face =
                            PoissonOperator_FaceGradientStencil(&metrics, nvert, face_cell, n, m, periodic);
                        PoissonOperator_AddFaceFlux(&face, own, n, upper ? 1.0 : -1.0, coefficients);
                    }
                }

                for (PetscInt s = 0; s < 19; s++) coefficients[s] *= -aj[k][j][i];
                PetscCall(MatSetValues(user->A, 1, &row, 19, columns, coefficients, INSERT_VALUES));
            }
        }
    }

    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(PoissonOperator_RestoreFaceMetrics(user, &metrics));
    PetscCall(MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY));

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Poisson operator assembled on level %d.\n", user->thislevel);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputePoissonRHS"
/**
 * @brief Implementation of \ref ComputePoissonRHS().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/poisson.h`.
 * @see ComputePoissonRHS()
 */
PetscErrorCode ComputePoissonRHS(UserCtx *user, Vec B)
{
    SimCtx            *simCtx = user->simCtx;
    const DMDALocalInfo info = user->info;
    const PetscInt     mx = info.mx, my = info.my, mz = info.mz;
    const PetscReal    dt = simCtx->dt;
    const Cmpnts    ***ucont;
    const PetscReal ***nvert, ***aj;
    PetscReal       ***rhs;
    PetscReal          local_sum = 0.0, global_sum = 0.0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(DMDAVecGetArray(user->da, B, &rhs));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, (void *)&ucont));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));

    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                if (i == 0 || i == mx - 1 || j == 0 || j == my - 1 || k == 0 || k == mz - 1 ||
                    nvert[k][j][i] > POISSON_SOLID_THRESHOLD) {
                    rhs[k][j][i] = 0.0;
                } else {
                    rhs[k][j][i] = -(ucont[k][j][i].x - ucont[k][j][i-1].x +
                                     ucont[k][j][i].y - ucont[k][j-1][i].y +
                                     ucont[k][j][i].z - ucont[k-1][j][i].z) / dt * aj[k][j][i] * COEF_TIME_ACCURACY;
                }
            }
        }
    }

    /* The integral of the right-hand side is the net volume flux into the domain carried
       by the uncorrected velocity. With Neumann pressure boundaries it must vanish for
       the equation to have a solution. */
    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                local_sum += rhs[k][j][i] / aj[k][j][i] * dt / COEF_TIME_ACCURACY;
            }
        }
    }
    PetscCallMPI(MPI_Allreduce(&local_sum, &global_sum, 1, MPIU_REAL, MPI_SUM, PetscObjectComm((PetscObject)B)));
    simCtx->poissonSourceImbalance = global_sum;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Poisson source imbalance: %le\n", (double)global_sum);

    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, (void *)&ucont));
    PetscCall(DMDAVecRestoreArray(user->da, B, &rhs));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdatePressure"
/**
 * @brief Implementation of \ref UpdatePressure().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/poisson.h`.
 * @see UpdatePressure()
 */
PetscErrorCode UpdatePressure(UserCtx *user)
{
    const FieldId periodic_fields[] = {FIELD_ID_P, FIELD_ID_PHI};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(VecAXPY(user->P, 1.0, user->Phi));
    PetscCall(SynchronizePeriodicCellFields(user, 2, periodic_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_P));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_PHI));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ProjectVelocity"
/**
 * @brief Implementation of \ref ProjectVelocity().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/poisson.h`.
 * @see ProjectVelocity()
 */
PetscErrorCode ProjectVelocity(UserCtx *user)
{
    SimCtx             *simCtx = user->simCtx;
    const DMDALocalInfo info = user->info;
    const PetscInt      m[3] = {info.mx, info.my, info.mz};
    const PetscInt      start[3] = {info.xs, info.ys, info.zs};
    const PetscInt      end[3] = {info.xs + info.xm, info.ys + info.ym, info.zs + info.zm};
    const PetscReal     scale = simCtx->dt / COEF_TIME_ACCURACY;
    const FieldId       staggered_fields[] = {FIELD_ID_UCONT};
    PetscBool           periodic[3];
    PetscInt            interior_start[3], interior_end[3];
    PoissonFaceMetrics  metrics;
    const PetscReal  ***nvert, ***phi;
    Cmpnts           ***ucont;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PoissonOperator_PeriodicAxes(user, periodic);
    for (PetscInt a = 0; a < 3; a++) {
        interior_start[a] = (start[a] == 0) ? 1 : start[a];
        interior_end[a]   = (end[a] == m[a]) ? m[a] - 1 : end[a];
    }

    PetscCall(PoissonOperator_GetFaceMetrics(user, &metrics));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lPhi, (void *)&phi));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));

    /* One pass per face orientation. Faces between two interior cells are corrected on
       every axis; a periodic axis also corrects its seam face at index 0. */
    for (PetscInt n = 0; n < 3; n++) {
        PetscInt lo[3], hi[3];
        for (PetscInt a = 0; a < 3; a++) { lo[a] = interior_start[a]; hi[a] = interior_end[a]; }
        if (periodic[n] && start[n] == 0) lo[n] = 0;
        hi[n] = PetscMin(hi[n], periodic[n] ? m[n] - 1 : m[n] - 2);

        for (PetscInt k = lo[2]; k < hi[2]; k++) {
            for (PetscInt j = lo[1]; j < hi[1]; j++) {
                for (PetscInt i = lo[0]; i < hi[0]; i++) {
                    const PetscInt c[3] = {i, j, k};
                    PetscInt       across[3] = {0, 0, 0};
                    PetscReal      difference[3];

                    across[n] = 1;
                    if (nvert[k][j][i] > POISSON_SOLID_THRESHOLD ||
                        PoissonOperator_At(nvert, c, across) > POISSON_SOLID_THRESHOLD) continue;

                    const PoissonFaceGradient face =
                        PoissonOperator_FaceGradientStencil(&metrics, nvert, c, n, m, periodic);
                    for (PetscInt b = 0; b < 3; b++) {
                        if (b == n) {
                            difference[b] = PoissonOperator_At(phi, c, across) - phi[k][j][i];
                            continue;
                        }
                        const PoissonTransverseDifference diff = face.diff[b];
                        PetscInt hi_lower[3] = {0, 0, 0}, hi_upper[3] = {0, 0, 0};
                        PetscInt lo_lower[3] = {0, 0, 0}, lo_upper[3] = {0, 0, 0};
                        hi_lower[b] = diff.hi; hi_upper[b] = diff.hi; hi_upper[n] = 1;
                        lo_lower[b] = diff.lo; lo_upper[b] = diff.lo; lo_upper[n] = 1;
                        difference[b] = (diff.weight == 0.0) ? 0.0 :
                            (PoissonOperator_At(phi, c, hi_lower) + PoissonOperator_At(phi, c, hi_upper) -
                             PoissonOperator_At(phi, c, lo_lower) - PoissonOperator_At(phi, c, lo_upper)) * diff.weight;
                    }

                    const PetscReal flux = difference[0] * face.dot[0] * face.aj +
                                           difference[1] * face.dot[1] * face.aj +
                                           difference[2] * face.dot[2] * face.aj;
                    const PetscReal correction = flux * scale;
                    if (n == 0)      ucont[k][j][i].x -= correction;
                    else if (n == 1) ucont[k][j][i].y -= correction;
                    else             ucont[k][j][i].z -= correction;
                }
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lPhi, (void *)&phi));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(PoissonOperator_RestoreFaceMetrics(user, &metrics));

    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered_fields));
    PetscCall(Contra2Cart(user));
    PetscCall(FinalizePostProjectionCellFields(user));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Removes the null space of the Neumann pressure problem from a level vector.
 *
 * PETSc first removes the global constant. This callback then subtracts the mean over
 * the interior fluid cells and zeroes the dummy layers and solid cells, which carry no
 * unknown.
 */
static PetscErrorCode PoissonMultigrid_RemoveNullSpace(MatNullSpace nullsp, Vec X, void *ctx)
{
    UserCtx            *user = (UserCtx *)ctx;
    const DMDALocalInfo info = user->info;
    const PetscInt      mx = info.mx, my = info.my, mz = info.mz;
    const PetscInt      xs = info.xs, xe = info.xs + info.xm;
    const PetscInt      ys = info.ys, ye = info.ys + info.ym;
    const PetscInt      zs = info.zs, ze = info.zs + info.zm;
    const PetscInt      lxs = (xs == 0) ? 1 : xs, lxe = (xe == mx) ? mx - 1 : xe;
    const PetscInt      lys = (ys == 0) ? 1 : ys, lye = (ye == my) ? my - 1 : ye;
    const PetscInt      lzs = (zs == 0) ? 1 : zs, lze = (ze == mz) ? mz - 1 : ze;
    const PetscReal  ***nvert;
    PetscReal        ***x;
    PetscReal           local[2] = {0.0, 0.0}, global[2];
    MPI_Comm            comm = PetscObjectComm((PetscObject)X);

    PetscFunctionBeginUser;
    (void)nullsp;
    PetscCall(DMDAVecGetArray(user->da, X, &x));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));

    for (PetscInt k = lzs; k < lze; k++) {
        for (PetscInt j = lys; j < lye; j++) {
            for (PetscInt i = lxs; i < lxe; i++) {
                if (nvert[k][j][i] < POISSON_SOLID_THRESHOLD) {
                    local[0] += x[k][j][i];
                    local[1] += 1.0;
                }
            }
        }
    }
    PetscCallMPI(MPI_Allreduce(&local[0], &global[0], 1, MPIU_REAL, MPI_SUM, comm));
    PetscCallMPI(MPI_Allreduce(&local[1], &global[1], 1, MPIU_REAL, MPI_SUM, comm));
    const PetscReal shift = global[0] / (-1.0 * global[1]);
    for (PetscInt k = lzs; k < lze; k++) {
        for (PetscInt j = lys; j < lye; j++) {
            for (PetscInt i = lxs; i < lxe; i++) {
                if (nvert[k][j][i] < POISSON_SOLID_THRESHOLD) x[k][j][i] += shift;
            }
        }
    }

    for (PetscInt k = zs; k < ze; k++) {
        for (PetscInt j = ys; j < ye; j++) {
            for (PetscInt i = xs; i < xe; i++) {
                if (i == 0 || i == mx - 1 || j == 0 || j == my - 1 || k == 0 || k == mz - 1 ||
                    nvert[k][j][i] > POISSON_SOLID_THRESHOLD) x[k][j][i] = 0.0;
            }
        }
    }

    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArray(user->da, X, &x));
    PetscFunctionReturn(0);
}

/**
 * @brief Coarse cell and interpolation direction along one axis for fine index @p f.
 *
 * Each fine cell takes 3/4 of its parent coarse cell and 1/4 of the coarse neighbour on
 * the side it lies towards. The first and last interior cells and semi-coarsened axes use
 * the parent only, as does a neighbour that is solid on the coarse grid.
 */
static void PoissonMultigrid_InterpolationParent(PetscInt f, PetscInt m, PetscInt semi,
                                                 PetscInt *coarse, PetscInt *direction)
{
    if (semi) {
        *coarse = f;
        *direction = 0;
        return;
    }
    *coarse = (f + 1) / 2;
    *direction = (f - 2 * (*coarse)) == 0 ? 1 : -1;
    if (f == 1 || f == m - 2) *direction = 0;
}

/**
 * @brief Prolongs a coarse-level correction to the next finer level (MatShell multiply).
 * @details The shell context is the fine-level UserCtx.
 */
static PetscErrorCode PoissonMultigrid_Interpolate(Mat P, Vec X, Vec F)
{
    UserCtx            *user, *coarse;
    DMDALocalInfo       info;
    Vec                 lX;
    const PetscReal  ***x, ***nvert, ***nvert_c;
    PetscReal        ***f;

    PetscFunctionBeginUser;
    PetscCall(MatShellGetContext(P, &user));
    coarse = user->user_c;
    info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
    const PetscInt lxs = (xs == 0) ? 1 : xs, lxe = (xe == mx) ? mx - 1 : xe;
    const PetscInt lys = (ys == 0) ? 1 : ys, lye = (ye == my) ? my - 1 : ye;
    const PetscInt lzs = (zs == 0) ? 1 : zs, lze = (ze == mz) ? mz - 1 : ze;

    PetscCall(DMGetLocalVector(coarse->da, &lX));
    PetscCall(DMGlobalToLocalBegin(coarse->da, X, INSERT_VALUES, lX));
    PetscCall(DMGlobalToLocalEnd(coarse->da, X, INSERT_VALUES, lX));
    PetscCall(DMDAVecGetArrayRead(coarse->da, lX, (void *)&x));
    PetscCall(DMDAVecGetArrayRead(coarse->da, coarse->lNvert, (void *)&nvert_c));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArray(user->da, F, &f));

    for (PetscInt k = lzs; k < lze; k++) {
        for (PetscInt j = lys; j < lye; j++) {
            for (PetscInt i = lxs; i < lxe; i++) {
                PetscInt ic, jc, kc, ia, ja, ka;

                PoissonMultigrid_InterpolationParent(i, mx, user->isc, &ic, &ia);
                PoissonMultigrid_InterpolationParent(j, my, user->jsc, &jc, &ja);
                PoissonMultigrid_InterpolationParent(k, mz, user->ksc, &kc, &ka);
                if (ka == -1 && nvert_c[kc-1][jc][ic] > POISSON_SOLID_THRESHOLD) ka = 0;
                else if (ka == 1 && nvert_c[kc+1][jc][ic] > POISSON_SOLID_THRESHOLD) ka = 0;
                if (ja == -1 && nvert_c[kc][jc-1][ic] > POISSON_SOLID_THRESHOLD) ja = 0;
                else if (ja == 1 && nvert_c[kc][jc+1][ic] > POISSON_SOLID_THRESHOLD) ja = 0;
                if (ia == -1 && nvert_c[kc][jc][ic-1] > POISSON_SOLID_THRESHOLD) ia = 0;
                else if (ia == 1 && nvert_c[kc][jc][ic+1] > POISSON_SOLID_THRESHOLD) ia = 0;

                f[k][j][i] = (x[kc   ][jc   ][ic   ] * 9 +
                              x[kc   ][jc+ja][ic   ] * 3 +
                              x[kc   ][jc   ][ic+ia] * 3 +
                              x[kc   ][jc+ja][ic+ia]) * 3./64. +
                             (x[kc+ka][jc   ][ic   ] * 9 +
                              x[kc+ka][jc+ja][ic   ] * 3 +
                              x[kc+ka][jc   ][ic+ia] * 3 +
                              x[kc+ka][jc+ja][ic+ia]) / 64.;
            }
        }
    }
    for (PetscInt k = zs; k < ze; k++) {
        for (PetscInt j = ys; j < ye; j++) {
            for (PetscInt i = xs; i < xe; i++) {
                if (i == 0 || i == mx - 1 || j == 0 || j == my - 1 || k == 0 || k == mz - 1 ||
                    nvert[k][j][i] > POISSON_SOLID_THRESHOLD) f[k][j][i] = 0.0;
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->da, F, &f));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(coarse->da, coarse->lNvert, (void *)&nvert_c));
    PetscCall(DMDAVecRestoreArrayRead(coarse->da, lX, (void *)&x));
    PetscCall(DMRestoreLocalVector(coarse->da, &lX));
    PetscFunctionReturn(0);
}

/**
 * @brief Restricts a fine-level residual to the next coarser level (MatShell multiply).
 *
 * Each coarse cell averages the eight fine cells it covers, weighting each by its fluid
 * fraction; coarse dummy and solid cells receive zero. The shell context is the
 * coarse-level UserCtx.
 */
static PetscErrorCode PoissonMultigrid_Restrict(Mat R, Vec X, Vec F)
{
    UserCtx            *user, *fine;
    DMDALocalInfo       info;
    Vec                 lX;
    const PetscReal  ***x, ***nvert, ***nvert_f;
    PetscReal        ***f;

    PetscFunctionBeginUser;
    PetscCall(MatShellGetContext(R, &user));
    fine = user->user_f;
    info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    const PetscInt ia = user->isc ? 0 : 1, ja = user->jsc ? 0 : 1, ka = user->ksc ? 0 : 1;

    PetscCall(DMGetLocalVector(fine->da, &lX));
    PetscCall(DMGlobalToLocalBegin(fine->da, X, INSERT_VALUES, lX));
    PetscCall(DMGlobalToLocalEnd(fine->da, X, INSERT_VALUES, lX));
    PetscCall(DMDAVecGetArrayRead(fine->da, lX, (void *)&x));
    PetscCall(DMDAVecGetArrayRead(fine->da, fine->lNvert, (void *)&nvert_f));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArray(user->da, F, &f));

    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                if (i == 0 || i == mx - 1 || j == 0 || j == my - 1 || k == 0 || k == mz - 1 ||
                    nvert[k][j][i] > POISSON_SOLID_THRESHOLD) {
                    f[k][j][i] = 0.0;
                    continue;
                }
                const PetscInt ih = user->isc ? i : 2 * i;
                const PetscInt jh = user->jsc ? j : 2 * j;
                const PetscInt kh = user->ksc ? k : 2 * k;
                f[k][j][i] = 0.125 *
                    (x[kh   ][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih   ]) +
                     x[kh   ][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih-ia]) +
                     x[kh   ][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih   ]) +
                     x[kh-ka][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih   ]) +
                     x[kh   ][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih-ia]) +
                     x[kh-ka][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih   ]) +
                     x[kh-ka][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih-ia]) +
                     x[kh-ka][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih-ia]));
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->da, F, &f));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(fine->da, fine->lNvert, (void *)&nvert_f));
    PetscCall(DMDAVecRestoreArrayRead(fine->da, lX, (void *)&x));
    PetscCall(DMRestoreLocalVector(fine->da, &lX));
    PetscFunctionReturn(0);
}

/**
 * @brief Gives each block factor of a block-Jacobi level solver a small diagonal shift, so
 *        the factorization of a nearly singular Neumann block does not fail on a zero
 *        pivot. Does nothing for any other preconditioner.
 */
static PetscErrorCode PoissonMultigrid_ShiftBlockFactors(KSP level_ksp)
{
    PC        level_pc;
    PCType    level_pc_type;
    PetscBool is_bjacobi = PETSC_FALSE;
    KSP      *block_ksp;
    PetscInt  nblocks;

    PetscFunctionBeginUser;
    PetscCall(KSPGetPC(level_ksp, &level_pc));
    PetscCall(PCGetType(level_pc, &level_pc_type));
    if (level_pc_type) PetscCall(PetscStrcmp(level_pc_type, PCBJACOBI, &is_bjacobi));
    if (!is_bjacobi) PetscFunctionReturn(0);

    PetscCall(KSPSetUp(level_ksp));
    PetscCall(PCBJacobiGetSubKSP(level_pc, &nblocks, NULL, &block_ksp));
    for (PetscInt b = 0; b < nblocks; b++) {
        PC block_pc;
        PetscCall(KSPGetPC(block_ksp[b], &block_pc));
        PetscCall(PCFactorSetShiftAmount(block_pc, 1.e-10));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Builds the multigrid solver for block @p bi and stores it in the finest level.
 *
 * Assembles the operator on every level, then configures the outer `ps_` Krylov solver
 * with a multiplicative V-cycle `PCMG`: shell restriction and interpolation between
 * levels, block-Jacobi smoothers by default, a coarse solve limited to 40 iterations at
 * relative tolerance 1e-8, and the Neumann null space on every level. PETSc options
 * override these defaults.
 *
 * Each smoother runs `pre_sweeps` iterations before the coarse correction and
 * `post_sweeps` after it. When the two differ, the post-smoother is a separate solver
 * that starts as a copy of the configured pre-smoother and reads further options under
 * `ps_mg_levels_N_up_`. Everything built here depends only on the grid metrics, the
 * solid field, and the boundary types; calling it again after one of those changes
 * rebuilds the solver.
 */
static PetscErrorCode PoissonMultigrid_Build(UserMG *usermg, PetscInt bi)
{
    MGCtx          *mgctx = usermg->mgctx;
    const PetscInt  levels = usermg->mglevels;
    UserCtx        *finest = &mgctx[levels - 1].user[bi];
    SimCtx         *simCtx = finest->simCtx;
    DualMonitorCtx *monitor;
    KSP             ksp;
    PC              pc;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Block %d: building the multigrid Poisson solver on %d levels.\n", bi, levels);

    for (PetscInt l = levels - 1; l >= 0; l--) PetscCall(AssemblePoissonOperator(&mgctx[l].user[bi]));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPAppendOptionsPrefix(ksp, "ps_"));

    /* The convergence log is opened and closed around every solve; see
       PoissonMultigrid_OpenConvergenceLog(). The monitor owns its context. */
    PetscCall(PetscNew(&monitor));
    monitor->block_id = bi;
    monitor->file_handle = NULL;
    PetscCall(KSPMonitorSet(ksp, DualKSPMonitor, monitor, DualMonitorDestroy));

    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCMG));
    PetscCall(PCMGSetLevels(pc, levels, NULL));
    PetscCall(PCMGSetCycleType(pc, PC_MG_CYCLE_V));
    PetscCall(PCMGSetType(pc, PC_MG_MULTIPLICATIVE));
    PetscCall(PCMGSetNumberSmooth(pc, simCtx->mg_preItr));

    for (PetscInt l = levels - 1; l > 0; l--) {
        UserCtx *fine = &mgctx[l].user[bi];
        UserCtx *coarse = &mgctx[l - 1].user[bi];
        const PetscInt m_c = coarse->info.xm * coarse->info.ym * coarse->info.zm;
        const PetscInt m_f = fine->info.xm * fine->info.ym * fine->info.zm;
        const PetscInt M_c = coarse->info.mx * coarse->info.my * coarse->info.mz;
        const PetscInt M_f = fine->info.mx * fine->info.my * fine->info.mz;

        PetscCall(MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, coarse, &fine->MR));
        PetscCall(MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, fine, &fine->MP));
        PetscCall(MatShellSetOperation(fine->MR, MATOP_MULT, (void (*)(void))PoissonMultigrid_Restrict));
        PetscCall(MatShellSetOperation(fine->MP, MATOP_MULT, (void (*)(void))PoissonMultigrid_Interpolate));
        PetscCall(PCMGSetRestriction(pc, l, fine->MR));
        PetscCall(PCMGSetInterpolation(pc, l, fine->MP));
    }

    for (PetscInt l = levels - 1; l >= 0; l--) {
        UserCtx *level = &mgctx[l].user[bi];
        KSP      level_ksp;
        PC       level_pc;

        if (l > 0) {
            PetscCall(PCMGGetSmoother(pc, l, &level_ksp));
        } else {
            PetscCall(PCMGGetCoarseSolve(pc, &level_ksp));
            PetscCall(KSPSetTolerances(level_ksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 40));
        }
        PetscCall(KSPSetOperators(level_ksp, level->A, level->A));
        PetscCall(KSPGetPC(level_ksp, &level_pc));
        PetscCall(PCSetType(level_pc, PCBJACOBI));
        PetscCall(KSPSetFromOptions(level_ksp));
        PetscCall(PoissonMultigrid_ShiftBlockFactors(level_ksp));

        PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &level->nullsp));
        PetscCall(MatNullSpaceSetFunction(level->nullsp, PoissonMultigrid_RemoveNullSpace, level));
        PetscCall(MatSetNullSpace(level->A, level->nullsp));
        PetscCall(PCMGSetResidual(pc, l, PCMGResidualDefault, level->A));
        PetscCall(KSPSetUp(level_ksp));

        if (l > 0 && simCtx->mg_preItr != simCtx->mg_poItr) {
            KSP post_smoother;

            /* PETSc creates the post-smoother as a copy of the pre-smoother's type,
               preconditioner and tolerances as configured so far. */
            PetscCall(PCMGGetSmootherUp(pc, l, &post_smoother));
            PetscCall(KSPAppendOptionsPrefix(post_smoother, "up_"));
            PetscCall(KSPSetOperators(post_smoother, level->A, level->A));
            PetscCall(KSPSetTolerances(post_smoother, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, simCtx->mg_poItr));
            PetscCall(KSPSetFromOptions(post_smoother));
            PetscCall(PoissonMultigrid_ShiftBlockFactors(post_smoother));
            PetscCall(KSPSetUp(post_smoother));
        }

        if (l < levels - 1) {
            PetscCall(MatCreateVecs(level->A, &level->R, NULL));
            PetscCall(PCMGSetRhs(pc, l, level->R));
        }
    }

    PetscCall(KSPSetOperators(ksp, finest->A, finest->A));
    PetscCall(MatSetNullSpace(finest->A, finest->nullsp));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSetUp(ksp));
    PetscCall(VecDuplicate(finest->P, &finest->B));
    finest->ksp = ksp;

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Prepares the convergence monitor for one solve and opens its log file on rank 0.
 *
 * The first step of a fresh run truncates the log; every other step appends, and the
 * first step of a continued run records where it resumed.
 */
static PetscErrorCode PoissonMultigrid_OpenConvergenceLog(KSP ksp, SimCtx *simCtx, PetscInt bi)
{
    DualMonitorCtx *monitor = NULL;
    const PetscBool first_step = (PetscBool)(simCtx->step == simCtx->StartStep + 1);

    PetscFunctionBeginUser;
    PetscCall(KSPGetMonitorContext(ksp, &monitor));
    monitor->step = simCtx->step;
    monitor->log_to_console = simCtx->ps_ksp_pic_monitor_true_residual;
    monitor->file_handle = NULL;
    if (simCtx->rank == 0) {
        char filename[PETSC_MAX_PATH_LEN + 128];

        PetscCall(PetscSNPrintf(filename, sizeof(filename),
                                "%s/Poisson_Solver_Convergence_History_Block_%d.log", simCtx->log_dir, bi));
        monitor->file_handle = fopen(filename, (first_step && !simCtx->continueMode) ? "w" : "a");
        PetscCheck(monitor->file_handle, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "Could not open KSP monitor log file: %s", filename);
        if (simCtx->continueMode && first_step) {
            PetscCall(PetscFPrintf(PETSC_COMM_SELF, monitor->file_handle,
                                   "# Continuation from step %" PetscInt_FMT "\n", simCtx->StartStep));
        }
        PetscCall(PetscFPrintf(PETSC_COMM_SELF, monitor->file_handle,
                               "--- Convergence for Timestep %d, Block %d ---\n", (int)simCtx->step, bi));
    }
    PetscFunctionReturn(0);
}

/** @brief Closes the log file opened by PoissonMultigrid_OpenConvergenceLog(). */
static PetscErrorCode PoissonMultigrid_CloseConvergenceLog(KSP ksp)
{
    DualMonitorCtx *monitor = NULL;

    PetscFunctionBeginUser;
    PetscCall(KSPGetMonitorContext(ksp, &monitor));
    if (monitor->file_handle) {
        fclose(monitor->file_handle);
        monitor->file_handle = NULL;
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PoissonSolver_Multigrid"
/**
 * @brief Implementation of \ref PoissonSolver_Multigrid().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/poisson.h`.
 * @see PoissonSolver_Multigrid()
 */
PetscErrorCode PoissonSolver_Multigrid(UserMG *usermg)
{
    SimCtx       *simCtx = usermg->mgctx[0].user[0].simCtx;
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting Multigrid Poisson Solve...\n");

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        UserCtx           *user = &usermg->mgctx[usermg->mglevels - 1].user[bi];
        KSPConvergedReason reason;

        if (!user->ksp) PetscCall(PoissonMultigrid_Build(usermg, bi));

        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered_fields));
        PetscCall(ComputePoissonRHS(user, user->B));

        PetscCall(PoissonMultigrid_OpenConvergenceLog(user->ksp, simCtx, bi));
        PetscCall(KSPSolve(user->ksp, user->B, user->Phi));
        PetscCall(PoissonMultigrid_CloseConvergenceLog(user->ksp));

        /* A non-finite residual or a preconditioner that could not be built leaves Phi
           meaningless, and the projection would carry it into the next momentum step,
           where it surfaces as a failure of the wrong solver. Stopping at max_it is this
           solve's normal mode and stays silent; any other divergence is reported, because
           the projection proceeds on that Phi. */
        PetscCall(KSPGetConvergedReason(user->ksp, &reason));
        PetscCheck(reason != KSP_DIVERGED_NANORINF && reason != KSP_DIVERGED_PC_FAILED,
                   PETSC_COMM_WORLD, PETSC_ERR_NOT_CONVERGED,
                   "Pressure Poisson solve on block %" PetscInt_FMT " failed at step %" PetscInt_FMT
                   " (KSP reason %s). Known causes: a multigrid hierarchy coarsened too far "
                   "(reduce poisson_solver.multigrid.levels or refine the grid), or a momentum "
                   "field that has already diverged, such as an explicit time step beyond its "
                   "stability limit.",
                   bi, simCtx->step, KSPConvergedReasons[reason]);
        if (reason < 0 && reason != KSP_DIVERGED_ITS) {
            LOG(GLOBAL, LOG_WARNING, "Pressure Poisson solve on block %" PetscInt_FMT
                " diverged at step %" PetscInt_FMT " (KSP reason %s); the projection uses the last iterate.\n",
                bi, simCtx->step, KSPConvergedReasons[reason]);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Multigrid Poisson Solve complete.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
