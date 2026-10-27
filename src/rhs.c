#include "rhs.h"
#include "poisson.h"   /* the face gradient shared with the projection */
#include "verification_sources.h"
#include "les.h"   /* SymTensor kernels for the gradient (Clark) closure term */

/** @brief Owned-cell bounds of the residual: interior cells, excluding the dummy layers. */
static inline void RHSInteriorBounds(const DMDALocalInfo *info, PetscInt lo[3], PetscInt hi[3])
{
    const PetscInt m[3] = {info->mx, info->my, info->mz};
    const PetscInt start[3] = {info->xs, info->ys, info->zs};
    const PetscInt end[3] = {info->xs + info->xm, info->ys + info->ym, info->zs + info->zm};

    for (PetscInt a = 0; a < 3; a++) {
        lo[a] = (start[a] == 0) ? 1 : start[a];
        hi[a] = (end[a] == m[a]) ? m[a] - 1 : end[a];
    }
}

/** Cells whose `nvert` exceeds this are blanked inner-body cells, which convect as fluid. */
#define RHS_INNER_BLANK 7.0

/** @brief Whether a cell takes part in convection: fluid, or a blanked inner-body cell. */
static inline PetscBool RHSConvectsAsFluid(PetscReal nvert)
{
    return (PetscBool)(nvert < 0.1 || nvert > RHS_INNER_BLANK);
}

/** @brief Cartesian velocity at cell @p c displaced by @p d along axis @p n. */
static inline Cmpnts RHSVelocityAt(const Cmpnts ***ucat, const PetscInt c[3], PetscInt n, PetscInt d)
{
    PetscInt cc[3] = {c[0], c[1], c[2]};

    cc[n] += d;
    return ucat[cc[2]][cc[1]][cc[0]];
}

/**
 * @brief QUICK flux through the face between cell @p c and `c + e_n`, upwinded by the sign of
 *        @p ucon. @p m1 and @p p2 are the offsets used for the far upwind cells, -1 and +2 in
 *        the interior; a one-sided face substitutes 0 or +1 for the side it cannot reach.
 */
static inline Cmpnts RHSQuickFlux(const Cmpnts ***ucat, const PetscInt c[3], PetscInt n,
                                  PetscReal ucon, PetscInt m1, PetscInt p2)
{
    const PetscReal coef = 0.125;
    const PetscReal up = ucon + fabs(ucon), um = ucon - fabs(ucon);
    const Cmpnts    u0 = RHSVelocityAt(ucat, c, n, 0), u1 = RHSVelocityAt(ucat, c, n, 1);
    const Cmpnts    um1 = RHSVelocityAt(ucat, c, n, m1), up2 = RHSVelocityAt(ucat, c, n, p2);
    Cmpnts          f;

    f.x = um * (coef * (-up2.x -2. * u1.x +3. * u0.x) +u1.x) +
          up * (coef * (-um1.x -2. * u0.x +3. * u1.x) +u0.x);
    f.y = um * (coef * (-up2.y -2. * u1.y +3. * u0.y) +u1.y) +
          up * (coef * (-um1.y -2. * u0.y +3. * u1.y) +u0.y);
    f.z = um * (coef * (-up2.z -2. * u1.z +3. * u0.z) +u1.z) +
          up * (coef * (-um1.z -2. * u0.z +3. * u1.z) +u0.z);
    return f;
}

/** @brief Central flux through the face between cell @p c and `c + e_n`. */
static inline Cmpnts RHSCentralFlux(const Cmpnts ***ucat, const PetscInt c[3], PetscInt n, PetscReal ucon)
{
    const Cmpnts u0 = RHSVelocityAt(ucat, c, n, 0), u1 = RHSVelocityAt(ucat, c, n, 1);
    Cmpnts       f;

    f.x = ucon * ( u0.x + u1.x );
    f.y = ucon * ( u0.y + u1.y );
    f.z = ucon * ( u0.z + u1.z );
    return f;
}

/**
 * @brief Convective flux of Cartesian momentum through the face between cell @p c and
 *        `c + e_n`, carried by half the contravariant flux on that face.
 *
 * Interior faces use the central flux when @p central (central differencing or LES) and QUICK
 * otherwise; with central differencing the two end faces of an axis are central as well. A face
 * at the low end of an axis or behind a solid cell is QUICK one-sided toward +n, and one at the
 * high end or ahead of a solid cell one-sided toward -n, except across a periodic seam, which
 * keeps the full stencil. A face meeting none of these carries no flux.
 */
static inline Cmpnts RHSConvectiveFaceFlux(const Cmpnts ***ucat, const PetscReal ***nvert,
                                           const PetscInt c[3], PetscInt n, PetscReal ucon,
                                           PetscInt m, PetscBool periodic, PetscBool central)
{
    const PetscInt  s = c[n];
    const PetscReal nv_0 = nvert[c[2]][c[1]][c[0]];
    const PetscReal nv_p = nvert[c[2] + (n == 2)][c[1] + (n == 1)][c[0] + (n == 0)];
    const Cmpnts    zero = {0.0, 0.0, 0.0};
    /* The cell below exists in the local array everywhere but at s = 0 on a bounded axis. */
    const PetscBool has_m = (PetscBool)(s > 0 || periodic);
    const PetscReal nv_m = has_m ? nvert[c[2] - (n == 2)][c[1] - (n == 1)][c[0] - (n == 0)] : 0.0;
    const PetscBool full_seam = (PetscBool)(periodic && nv_m < 0.1 && nv_p < 0.1);

    if (s > 0 && s < m - 2 && RHSConvectsAsFluid(nv_p) && RHSConvectsAsFluid(nv_m)) {
        return central ? RHSCentralFlux(ucat, c, n, ucon) : RHSQuickFlux(ucat, c, n, ucon, -1, 2);
    }
    if (central && (s == 0 || s == m - 2) && RHSConvectsAsFluid(nv_p) && RHSConvectsAsFluid(nv_0)) {
        return RHSCentralFlux(ucat, c, n, ucon);
    }
    if (s < m - 2 && (s == 0 || nv_m > 0.1)) {
        return RHSQuickFlux(ucat, c, n, ucon, (full_seam && s == 0) ? -1 : 0, 2);
    }
    if (s > 0 && (s == m - 2 || nv_p > 0.1)) {
        return RHSQuickFlux(ucat, c, n, ucon, -1, (full_seam && s == m - 2) ? 2 : 1);
    }
    return zero;
}

/**
 * @brief Implementation of \ref ComputeConvection().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/rhs.h`.
 * @see ComputeConvection()
 */
PetscErrorCode ComputeConvection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv)
{
    const SimCtx        *simCtx = user->simCtx;
    const DMDALocalInfo *info = &user->info;
    const PetscInt       m[3] = {info->mx, info->my, info->mz};
    const PetscBool      central = (PetscBool)(simCtx->les || simCtx->central);
    PetscBool            periodic[3];
    PetscInt             lo[3], hi[3];
    Vec                  Fp[3];
    Cmpnts            ***fp[3], ***conv;
    const Cmpnts      ***ucont, ***ucat;
    const PetscReal   ***nvert;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(info, lo, hi);
    PoissonOperator_PeriodicAxes(user, periodic);
    PetscCall(PreparePeriodicQuickStencilFields(user, Ucat, user->lNvert));
    PetscCall(VecSet(Conv, 0.0));
    /* Face fluxes per axis, from fda's pool and zeroed as fresh vectors would be: faces that
       carry no flux must read as zero. */
    for (PetscInt n = 0; n < 3; n++) {
        PetscCall(DMGetLocalVector(user->fda, &Fp[n]));
        PetscCall(VecZeroEntries(Fp[n]));
        PetscCall(DMDAVecGetArray(user->fda, Fp[n], &fp[n]));
    }
    PetscCall(DMDAVecGetArrayRead(user->fda, Ucont, (void *)&ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, Ucat, (void *)&ucat));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArray(user->fda, Conv, &conv));

    /* Fluxes on every face of the owned cells: along n, from the face below the first cell. */
    for (PetscInt n = 0; n < 3; n++) {
        PetscInt flo[3] = {lo[0], lo[1], lo[2]};

        flo[n] -= 1;
        for (PetscInt k = flo[2]; k < hi[2]; k++) {
            for (PetscInt j = flo[1]; j < hi[1]; j++) {
                for (PetscInt i = flo[0]; i < hi[0]; i++) {
                    const PetscInt  c[3] = {i, j, k};
                    const PetscReal ucon = (&ucont[k][j][i].x)[n] * 0.5;

                    fp[n][k][j][i] = RHSConvectiveFaceFlux(ucat, (const PetscReal ***)nvert, c, n, ucon,
                                                           m[n], periodic[n], central);
                }
            }
        }
    }

    /* Convective term in Cartesian components: the net flux out of each cell. */
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                conv[k][j][i].x =
                    fp[0][k][j][i].x - fp[0][k][j][i-1].x +
                    fp[1][k][j][i].x - fp[1][k][j-1][i].x +
                    fp[2][k][j][i].x - fp[2][k-1][j][i].x;
                conv[k][j][i].y =
                    fp[0][k][j][i].y - fp[0][k][j][i-1].y +
                    fp[1][k][j][i].y - fp[1][k][j-1][i].y +
                    fp[2][k][j][i].y - fp[2][k-1][j][i].y;
                conv[k][j][i].z =
                    fp[0][k][j][i].z - fp[0][k][j][i-1].z +
                    fp[1][k][j][i].z - fp[1][k][j-1][i].z +
                    fp[2][k][j][i].z - fp[2][k-1][j][i].z;
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, Conv, &conv));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Ucat, (void *)&ucat));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Ucont, (void *)&ucont));
    for (PetscInt n = 0; n < 3; n++) {
        PetscCall(DMDAVecRestoreArray(user->fda, Fp[n], &fp[n]));
        PetscCall(DMRestoreLocalVector(user->fda, &Fp[n]));
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Convective term calculated .\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/** @brief Whether a cell is a blanked solid cell, which one-sided viscous differences avoid. */
static PICURV_ALWAYS_INLINE PetscBool RHSIsBlankedSolid(PetscReal nvert)
{
    return (PetscBool)(nvert > 0.5 && nvert < RHS_INNER_BLANK);
}

/** @brief Velocity at `c + p e_b + q e_n`, for the transverse stencils below. */
static PICURV_ALWAYS_INLINE Cmpnts RHSVelocityAtBN(const Cmpnts ***ucat, const PetscInt c[3],
                                                   PetscInt b, PetscInt p, PetscInt n, PetscInt q)
{
    return ucat[c[2] + (b == 2) * p + (n == 2) * q][c[1] + (b == 1) * p + (n == 1) * q][c[0] + (b == 0) * p + (n == 0) * q];
}

/** @brief Solid field at `c + p e_b + q e_n`. */
static PICURV_ALWAYS_INLINE PetscReal RHSNvertAtBN(const PetscReal ***nvert, const PetscInt c[3],
                                                   PetscInt b, PetscInt p, PetscInt n, PetscInt q)
{
    return nvert[c[2] + (b == 2) * p + (n == 2) * q][c[1] + (b == 1) * p + (n == 1) * q][c[0] + (b == 0) * p + (n == 0) * q];
}

/** @brief `(u(c+hi e_b+e_n) + u(c+hi e_b) - u(c+lo e_b+e_n) - u(c+lo e_b)) * weight`. */
static PICURV_ALWAYS_INLINE Cmpnts RHSTransverseVelocityDifference(const Cmpnts ***ucat, const PetscInt c[3],
                                                                   PetscInt n, PetscInt b, PetscInt lo,
                                                                   PetscInt hi, PetscReal weight)
{
    const Cmpnts a1 = RHSVelocityAtBN(ucat, c, b, hi, n, 1), a0 = RHSVelocityAtBN(ucat, c, b, hi, n, 0);
    const Cmpnts b1 = RHSVelocityAtBN(ucat, c, b, lo, n, 1), b0 = RHSVelocityAtBN(ucat, c, b, lo, n, 0);
    Cmpnts       d;

    d.x = (a1.x + a0.x - b1.x - b0.x) * weight;
    d.y = (a1.y + a0.y - b1.y - b0.y) * weight;
    d.z = (a1.z + a0.z - b1.z - b0.z) * weight;
    return d;
}

/**
 * @brief Velocity difference along transverse axis @p b on the face between cell @p c and
 *        `c + e_n`: the average over the two cells beside the face, one-sided away from a
 *        blanked solid cell. Each branch has constant offsets once @p n and @p b are constants.
 */
static PICURV_ALWAYS_INLINE Cmpnts RHSTransverseFaceDifference(const Cmpnts ***ucat, const PetscReal ***nvert,
                                                               const PetscInt c[3], PetscInt n, PetscInt b)
{
    if (RHSIsBlankedSolid(RHSNvertAtBN(nvert, c, b, 1, n, 0)) || RHSIsBlankedSolid(RHSNvertAtBN(nvert, c, b, 1, n, 1)))
        return RHSTransverseVelocityDifference(ucat, c, n, b, -1, 0, 0.5);
    if (RHSIsBlankedSolid(RHSNvertAtBN(nvert, c, b, -1, n, 0)) || RHSIsBlankedSolid(RHSNvertAtBN(nvert, c, b, -1, n, 1)))
        return RHSTransverseVelocityDifference(ucat, c, n, b, 0, 1, 0.5);
    return RHSTransverseVelocityDifference(ucat, c, n, b, -1, 1, 0.25);
}

/**
 * @brief Velocity differences on the face between cell @p c and `c + e_n`, one per
 *        computational direction: across the face along n, and along each transverse axis.
 */
static PICURV_ALWAYS_INLINE void RHSFaceVelocityDifferences(const Cmpnts ***ucat, const PetscReal ***nvert,
                                                            const PetscInt c[3], PetscInt n, Cmpnts d[3])
{
    const Cmpnts u0 = RHSVelocityAtBN(ucat, c, n, 0, n, 0), u1 = RHSVelocityAtBN(ucat, c, n, 0, n, 1);

    d[n].x = u1.x - u0.x; d[n].y = u1.y - u0.y; d[n].z = u1.z - u0.z;
    d[(n + 1) % 3] = RHSTransverseFaceDifference(ucat, nvert, c, n, (n + 1) % 3);
    d[(n + 2) % 3] = RHSTransverseFaceDifference(ucat, nvert, c, n, (n + 2) % 3);
}

/**
 * @brief Viscous stress through the face between cell @p c and `c + e_n`, in Cartesian
 *        components, from the face differences @p d and the face metrics.
 *
 * With F_b the contravariant base vectors on the face and N = F_n its normal, the flux of each
 * velocity component U is `(sum_b (F_b . N) dU/dxi_b + sum_c (sum_b dU_c/dxi_b F_b) . e_U N_c)
 * aj`: the full stress, gradient plus transpose. The eddy viscosity @p nu_t carries it when
 * LES is on, and the molecular viscosity @p nu is added; with the gradient (Clark) model the
 * tensor `sum_b d_b d_b^T` contracted with N, over 12, is subtracted.
 */
static PICURV_ALWAYS_INLINE Cmpnts RHSViscousFaceFlux(const PoissonFaceMetrics *metrics, const Cmpnts d[3],
                                        const PetscInt c[3], PetscInt n, PetscReal nu,
                                        PetscReal nu_t, PetscBool les, PetscBool clark)
{
    const Cmpnts    F0 = metrics->metric[n][0][c[2]][c[1]][c[0]];
    const Cmpnts    F1 = metrics->metric[n][1][c[2]][c[1]][c[0]];
    const Cmpnts    F2 = metrics->metric[n][2][c[2]][c[1]][c[0]];
    const Cmpnts    N = (n == 0) ? F0 : (n == 1) ? F1 : F2;
    const PetscReal ajc = metrics->aj[n][c[2]][c[1]][c[0]];
    const PetscReal g0 = F0.x * N.x + F0.y * N.y + F0.z * N.z;
    const PetscReal g1 = F1.x * N.x + F1.y * N.y + F1.z * N.z;
    const PetscReal g2 = F2.x * N.x + F2.y * N.y + F2.z * N.z;
    /* r_<velocity component><Cartesian axis> = sum_b dU/dxi_b F_b[axis] */
    const PetscReal r11 = d[0].x * F0.x + d[1].x * F1.x + d[2].x * F2.x;
    const PetscReal r21 = d[0].y * F0.x + d[1].y * F1.x + d[2].y * F2.x;
    const PetscReal r31 = d[0].z * F0.x + d[1].z * F1.x + d[2].z * F2.x;
    const PetscReal r12 = d[0].x * F0.y + d[1].x * F1.y + d[2].x * F2.y;
    const PetscReal r22 = d[0].y * F0.y + d[1].y * F1.y + d[2].y * F2.y;
    const PetscReal r32 = d[0].z * F0.y + d[1].z * F1.y + d[2].z * F2.y;
    const PetscReal r13 = d[0].x * F0.z + d[1].x * F1.z + d[2].x * F2.z;
    const PetscReal r23 = d[0].y * F0.z + d[1].y * F1.z + d[2].y * F2.z;
    const PetscReal r33 = d[0].z * F0.z + d[1].z * F1.z + d[2].z * F2.z;
    const PetscReal ex = g0 * d[0].x + g1 * d[1].x + g2 * d[2].x + r11 * N.x + r21 * N.y + r31 * N.z;
    const PetscReal ey = g0 * d[0].y + g1 * d[1].y + g2 * d[2].y + r12 * N.x + r22 * N.y + r32 * N.z;
    const PetscReal ez = g0 * d[0].z + g1 * d[1].z + g2 * d[2].z + r13 * N.x + r23 * N.y + r33 * N.z;
    Cmpnts          f = {0.0, 0.0, 0.0};

    if (les) {
        f.x = ex * ajc * (nu_t);
        f.y = ey * ajc * (nu_t);
        f.z = ez * ajc * (nu_t);
    }
    f.x += ex * ajc * (nu);
    f.y += ey * ajc * (nu);
    f.z += ez * ajc * (nu);
    if (clark) {
        /* d_b is the velocity difference between neighbouring cells along direction b, so it
           already carries that direction's spacing: sum_b d_b d_b^T is Delta_k^2 du_i/dx_k
           du_j/dx_k as it stands, and no further cell extent belongs in it. */
        const SymTensor gradient_tensor =
            SymTensorCombine(1.0, SymTensorSelfOuter(d[0]), 1.0,
                             SymTensorCombine(1.0, SymTensorSelfOuter(d[1]),
                                              1.0, SymTensorSelfOuter(d[2])));
        const Cmpnts gradient_flux = SymTensorTimesVector(gradient_tensor, N);

        f.x -= gradient_flux.x / 12.;
        f.y -= gradient_flux.y / 12.;
        f.z -= gradient_flux.z / 12.;
    }
    return f;
}

/**
 * @brief Viscous stresses on every face of orientation @p n of the owned cells, from the face
 *        below the first cell; @p n is a constant at each call, so each expands to its own loop.
 */
static PICURV_ALWAYS_INLINE void RHSViscousFluxesAlong(const UserCtx *user, const PoissonFaceMetrics *metrics,
                                                       const Cmpnts ***ucat, const PetscReal ***nvert,
                                                       const PetscReal ***lnu_t, const PetscReal ***lnu_wall,
                                                       const PetscInt lo[3], const PetscInt hi[3],
                                                       const PetscInt m[3], PetscInt n, PetscBool les,
                                                       PetscBool clark, Cmpnts ***fp)
{
    const PetscBool wall_lo = BoundaryFaceHasWallShear(&user->boundary_faces[2 * n]);
    const PetscBool wall_hi = BoundaryFaceHasWallShear(&user->boundary_faces[2 * n + 1]);
    const PetscReal ren = user->simCtx->ren;
    PetscInt        flo[3] = {lo[0], lo[1], lo[2]};
    PetscInt        across[3] = {0, 0, 0};

    flo[n] -= 1;
    across[n] = 1;
    for (PetscInt k = flo[2]; k < hi[2]; k++) {
        for (PetscInt j = flo[1]; j < hi[1]; j++) {
            for (PetscInt i = flo[0]; i < hi[0]; i++) {
                const PetscInt c[3] = {i, j, k};
                double         nu = 1./ren, nu_t = 0;
                Cmpnts         d[3];

                if (les) {
                    nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k + across[2]][j + across[1]][i + across[0]]);
                    /* Zero suits a wall-resolved run, where the subgrid stress vanishes at the
                       wall; with a wall model this face is where the modelled stress enters. */
                    if ((wall_lo && c[n] == 0) || (wall_hi && c[n] == m[n] - 2)) {
                        PetscInt w[3] = {i, j, k};
                        w[n] = (c[n] == 0) ? 1 : m[n] - 2;
                        nu_t = lnu_wall ? lnu_wall[w[2]][w[1]][w[0]] : 0.0;
                    }
                }
                RHSFaceVelocityDifferences(ucat, nvert, c, n, d);
                fp[k][j][i] = RHSViscousFaceFlux(metrics, d, c, n, nu, nu_t, les, clark);
            }
        }
    }
}

/**
 * @brief Implementation of \ref ComputeViscous().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/rhs.h`.
 * @see ComputeViscous()
 */
PetscErrorCode ComputeViscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc)
{
    const SimCtx        *simCtx = user->simCtx;
    const DMDALocalInfo *info = &user->info;
    const PetscInt       m[3] = {info->mx, info->my, info->mz};
    const PetscBool      les = (PetscBool)(simCtx->les != 0);
    const PetscBool      clark = (PetscBool)(simCtx->les_gradient_model != 0);
    PetscInt             lo[3], hi[3];
    PoissonFaceMetrics   metrics;
    Vec                  Fp[3];
    Cmpnts            ***fp[3], ***visc;
    const Cmpnts      ***ucat;
    const PetscReal   ***nvert, ***lnu_t = NULL, ***lnu_wall = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    (void)Ucont;
    RHSInteriorBounds(info, lo, hi);
    PetscCall(VecSet(Visc, 0.0));
    for (PetscInt n = 0; n < 3; n++) {
        PetscCall(DMGetLocalVector(user->fda, &Fp[n]));
        PetscCall(VecZeroEntries(Fp[n]));
        PetscCall(DMDAVecGetArray(user->fda, Fp[n], &fp[n]));
    }
    PetscCall(PoissonOperator_GetFaceMetrics(user, &metrics));
    PetscCall(DMDAVecGetArrayRead(user->fda, Ucat, (void *)&ucat));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    if (les) PetscCall(DMDAVecGetArrayRead(user->da, user->lNu_t, (void *)&lnu_t));
    /* A wall model replaces the subgrid viscosity at its own face rather than adding to it, so
       the flux there carries the stress the model computed instead of the molecular fraction. */
    if (simCtx->wallfunction) PetscCall(DMDAVecGetArrayRead(user->da, user->lNu_Wall, (void *)&lnu_wall));
    PetscCall(DMDAVecGetArray(user->fda, Visc, &visc));

    /* Stresses on every face of the owned cells, one constant direction per call. */
    RHSViscousFluxesAlong(user, &metrics, ucat, nvert, lnu_t, lnu_wall, lo, hi, m, 0, les, clark, fp[0]);
    RHSViscousFluxesAlong(user, &metrics, ucat, nvert, lnu_t, lnu_wall, lo, hi, m, 1, les, clark, fp[1]);
    RHSViscousFluxesAlong(user, &metrics, ucat, nvert, lnu_t, lnu_wall, lo, hi, m, 2, les, clark, fp[2]);

    /* Viscous term in Cartesian components: the net stress out of each cell. */
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                visc[k][j][i].x =
                    (fp[0][k][j][i].x - fp[0][k][j][i-1].x +
                     fp[1][k][j][i].x - fp[1][k][j-1][i].x +
                     fp[2][k][j][i].x - fp[2][k-1][j][i].x);
                visc[k][j][i].y =
                    (fp[0][k][j][i].y - fp[0][k][j][i-1].y +
                     fp[1][k][j][i].y - fp[1][k][j-1][i].y +
                     fp[2][k][j][i].y - fp[2][k-1][j][i].y);
                visc[k][j][i].z =
                    (fp[0][k][j][i].z - fp[0][k][j][i-1].z +
                     fp[1][k][j][i].z - fp[1][k][j-1][i].z +
                     fp[2][k][j][i].z - fp[2][k-1][j][i].z);
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, Visc, &visc));
    if (lnu_wall) PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNu_Wall, (void *)&lnu_wall));
    if (lnu_t) PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNu_t, (void *)&lnu_t));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Ucat, (void *)&ucat));
    PetscCall(PoissonOperator_RestoreFaceMetrics(user, &metrics));
    for (PetscInt n = 0; n < 3; n++) {
        PetscCall(DMDAVecRestoreArray(user->fda, Fp[n], &fp[n]));
        PetscCall(DMRestoreLocalVector(user->fda, &Fp[n]));
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Viscous terms calculated .\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeBodyForces()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeBodyForces(UserCtx *user, Vec Rct)
{
	PetscErrorCode ierr;
    PetscFunctionBeginUser;
	PROFILE_HELPER_BEGIN;

    // --- 1. Apply momentum source for driven channel/pipe flows ---
    // This function will internally check if a driven flow BC is active.
    ierr = ComputeDrivenChannelFlowSource(user, Rct); CHKERRQ(ierr);

    // --- 2. (Future Extension) Apply gravitational force ---
    // if (user->simCtx->gravityEnabled) {
    //     ierr = ApplyGravitationalForce(user, Rhs); CHKERRQ(ierr);
    // }
    //
    // Adding a body force here: see the contract in include/BodyForces.h.
    // In particular, this function runs once per RESIDUAL EVALUATION, not once
    // per timestep, so any force carrying state across calls (a filter, ramp,
    // moving average, or integral term) must gate its update on simCtx->step.

    PROFILE_HELPER_END;
    PetscFunctionReturn(0); 
}

/**
 * @brief Contravariant cell-centre components of the spatial terms, Rct = aj F_n . (Visc - Conv),
 *        with each base vector F_n averaged from its two faces to the cell centre.
 */
static PetscErrorCode RHSContravariantSpatialTerms(UserCtx *user, Vec Conv, Vec Visc, Vec Rct)
{
    PetscInt         lo[3], hi[3];
    const Cmpnts  ***csi, ***eta, ***zet, ***conv, ***visc;
    const PetscReal ***aj;
    Cmpnts        ***rct;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(&user->info, lo, hi);
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (void *)&csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecGetArrayRead(user->fda, Conv, (void *)&conv));
    PetscCall(DMDAVecGetArrayRead(user->fda, Visc, (void *)&visc));
    PetscCall(DMDAVecGetArray(user->fda, Rct, &rct));
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                const Cmpnts rc = {visc[k][j][i].x - conv[k][j][i].x,
                                   visc[k][j][i].y - conv[k][j][i].y,
                                   visc[k][j][i].z - conv[k][j][i].z};
                rct[k][j][i].x = aj[k][j][i] *
                    (0.5 * (csi[k][j][i].x + csi[k][j][i-1].x) * rc.x +
                     0.5 * (csi[k][j][i].y + csi[k][j][i-1].y) * rc.y +
                     0.5 * (csi[k][j][i].z + csi[k][j][i-1].z) * rc.z);
                rct[k][j][i].y = aj[k][j][i] *
                    (0.5 * (eta[k][j][i].x + eta[k][j-1][i].x) * rc.x +
                     0.5 * (eta[k][j][i].y + eta[k][j-1][i].y) * rc.y +
                     0.5 * (eta[k][j][i].z + eta[k][j-1][i].z) * rc.z);
                rct[k][j][i].z = aj[k][j][i] *
                    (0.5 * (zet[k][j][i].x + zet[k-1][j][i].x) * rc.x +
                     0.5 * (zet[k][j][i].y + zet[k-1][j][i].y) * rc.y +
                     0.5 * (zet[k][j][i].z + zet[k-1][j][i].z) * rc.z);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, Rct, &rct));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Visc, (void *)&visc));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Conv, (void *)&conv));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (void *)&csi));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/** @brief Face residual of component @p n at cell @p c: the two-cell average of @p rct less
 *         the pressure-gradient flux through the face. */
static PICURV_ALWAYS_INLINE PetscReal RHSFaceTerm(const PoissonFaceMetrics *metrics, const PetscReal ***nvert,
                                    const PetscReal ***p, const Cmpnts ***rct, const PetscInt c[3],
                                    PetscInt n, const PetscInt m[3], const PetscBool periodic[3])
{
    const PetscReal *own = &rct[c[2]][c[1]][c[0]].x;
    const PetscReal *next = &rct[c[2] + (n == 2)][c[1] + (n == 1)][c[0] + (n == 0)].x;

    return 0.5 * (own[n] + next[n]) - PoissonOperator_FaceGradientFlux(metrics, nvert, p, c, n, m, periodic);
}

/**
 * @brief Face residual from the ghosted cell-centre terms: the average of @p Rct over the two
 *        cells beside each face, less the pressure-gradient flux through it.
 *
 * The flux is the projection's face gradient (PoissonOperator_FaceGradientFlux()) applied to
 * `P`. Every interior face is set, and on a periodic axis also the seam face at index 0.
 */
static PetscErrorCode RHSFaceTermsWithPressureGradient(UserCtx *user, Vec Rct, Vec Rhs)
{
    const DMDALocalInfo *info = &user->info;
    const PetscInt       m[3] = {info->mx, info->my, info->mz};
    const PetscInt       start[3] = {info->xs, info->ys, info->zs};
    PetscBool            periodic[3];
    PetscInt             interior_lo[3], interior_hi[3];
    PoissonFaceMetrics   metrics;
    const PetscReal   ***nvert, ***p;
    const Cmpnts      ***rct;
    Cmpnts            ***rhs;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(info, interior_lo, interior_hi);
    PoissonOperator_PeriodicAxes(user, periodic);
    PetscCall(PoissonOperator_GetFaceMetrics(user, &metrics));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lP, (void *)&p));
    PetscCall(DMDAVecGetArrayRead(user->fda, Rct, (void *)&rct));
    PetscCall(DMDAVecGetArray(user->fda, Rhs, &rhs));
    /* One pass over the owned cells sets all three face components; with each direction a
       constant, every call below compiles to its own straight-line stencil. */
    for (PetscInt k = interior_lo[2]; k < interior_hi[2]; k++) {
        for (PetscInt j = interior_lo[1]; j < interior_hi[1]; j++) {
            for (PetscInt i = interior_lo[0]; i < interior_hi[0]; i++) {
                const PetscInt c[3] = {i, j, k};

                rhs[k][j][i].x = RHSFaceTerm(&metrics, nvert, p, rct, c, 0, m, periodic);
                rhs[k][j][i].y = RHSFaceTerm(&metrics, nvert, p, rct, c, 1, m, periodic);
                rhs[k][j][i].z = RHSFaceTerm(&metrics, nvert, p, rct, c, 2, m, periodic);
            }
        }
    }
    /* On a periodic axis the seam face at index 0 also carries an unknown. */
    if (periodic[0] && start[0] == 0) {
        for (PetscInt k = interior_lo[2]; k < interior_hi[2]; k++)
            for (PetscInt j = interior_lo[1]; j < interior_hi[1]; j++) {
                const PetscInt c[3] = {0, j, k};
                rhs[k][j][0].x = RHSFaceTerm(&metrics, nvert, p, rct, c, 0, m, periodic);
            }
    }
    if (periodic[1] && start[1] == 0) {
        for (PetscInt k = interior_lo[2]; k < interior_hi[2]; k++)
            for (PetscInt i = interior_lo[0]; i < interior_hi[0]; i++) {
                const PetscInt c[3] = {i, 0, k};
                rhs[k][0][i].y = RHSFaceTerm(&metrics, nvert, p, rct, c, 1, m, periodic);
            }
    }
    if (periodic[2] && start[2] == 0) {
        for (PetscInt j = interior_lo[1]; j < interior_hi[1]; j++)
            for (PetscInt i = interior_lo[0]; i < interior_hi[0]; i++) {
                const PetscInt c[3] = {i, j, 0};
                rhs[0][j][i].z = RHSFaceTerm(&metrics, nvert, p, rct, c, 2, m, periodic);
            }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, Rhs, &rhs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Rct, (void *)&rct));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lP, (void *)&p));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(PoissonOperator_RestoreFaceMetrics(user, &metrics));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Zeroes the residual rows that carry no unknown: the component a 2D run disables, every
 *        component of a solid cell, and each face component whose far cell is solid.
 */
static PetscErrorCode RHSZeroInactiveRows(UserCtx *user, Vec Rhs)
{
    const PetscInt     two_d = user->simCtx->TwoD;
    PetscInt           lo[3], hi[3];
    const PetscReal ***nvert;
    Cmpnts          ***rhs;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(&user->info, lo, hi);
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArray(user->fda, Rhs, &rhs));
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                if (two_d == 1)      rhs[k][j][i].x = 0.0;
                else if (two_d == 2) rhs[k][j][i].y = 0.0;
                else if (two_d == 3) rhs[k][j][i].z = 0.0;
                if (nvert[k][j][i] > 0.1) {
                    rhs[k][j][i].x = 0.0; rhs[k][j][i].y = 0.0; rhs[k][j][i].z = 0.0;
                }
                if (nvert[k][j][i+1] > 0.1) rhs[k][j][i].x = 0.0;
                if (nvert[k][j+1][i] > 0.1) rhs[k][j][i].y = 0.0;
                if (nvert[k+1][j][i] > 0.1) rhs[k][j][i].z = 0.0;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, Rhs, &rhs));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ComputeRHS().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/rhs.h`.
 * @see ComputeRHS()
 */
PetscErrorCode ComputeRHS(UserCtx *user, Vec Rhs)
{
    const FieldId cell_fields[] = {FIELD_ID_UCAT};
    Vec           Conv, Visc, Rct;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %" PetscInt_FMT ": Computing RHS (FormFunction1)...\n",
              user->simCtx->rank, user->_this);
    /* Work vectors from fda's pool, zeroed as fresh vectors would be: the residual is evaluated
       many times per step, and allocating them each time costs more than the zeroing. */
    PetscCall(DMGetLocalVector(user->fda, &Conv)); PetscCall(VecZeroEntries(Conv));
    PetscCall(DMGetLocalVector(user->fda, &Visc)); PetscCall(VecZeroEntries(Visc));
    PetscCall(DMGetLocalVector(user->fda, &Rct));  PetscCall(VecZeroEntries(Rct));

    /* Cartesian velocity, with its periodic copies and ghosts, from the contravariant flux. */
    PetscCall(Contra2Cart(user));
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));

    /* Spatial terms in Cartesian components at cell centres, then contravariant. */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Calculating convective terms...\n");
    PetscCall(ComputeConvection(user, user->lUcont, user->lUcat, Conv));
    if (!user->simCtx->invicid) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "  Calculating viscous terms...\n");
        PetscCall(ComputeViscous(user, user->lUcont, user->lUcat, Visc));
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Converting Cartesian RHS to Contravariant RHS...\n");
    PetscCall(RHSContravariantSpatialTerms(user, Conv, Visc, Rct));
    PetscCall(ComputeBodyForces(user, Rct));
    PetscCall(SynchronizePeriodicLocalStaggeredField(user, Rct));

    /* Face residual: averaged spatial terms less the pressure gradient, then inactive rows. */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Adding pressure gradient term to RHS...\n");
    PetscCall(RHSFaceTermsWithPressureGradient(user, Rct, Rhs));
    PetscCall(RHSZeroInactiveRows(user, Rhs));
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Final cleanup and edge-cases complete .\n");

    PetscCall(DMRestoreLocalVector(user->fda, &Rct));
    PetscCall(DMRestoreLocalVector(user->fda, &Visc));
    PetscCall(DMRestoreLocalVector(user->fda, &Conv));
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %" PetscInt_FMT ": RHS computation complete.\n",
              user->simCtx->rank, user->_this);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
