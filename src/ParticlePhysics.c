#include "ParticlePhysics.h"
#include "verification_sources.h"
#include "ParticleMotion.h"
#include "interpolation.h"

#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256
#endif

static PetscErrorCode ComputeCellScalarGradient(UserCtx *user, Vec lscalar, Vec lcount, Vec gradient);

/**
 * @brief Squared IEM mixing length of a cell: its volume to the power 2/3, from the reciprocal
 *        volume Aj, the volume floored at 1e-14 so a degenerate or cut cell cannot divide by zero.
 */
static inline PetscReal IEMMixingLength2(PetscReal aj)
{
    PetscReal vol = (aj > 1.0e-14) ? (1.0 / aj) : 1.0e-14;
    if (vol < 1.0e-14) vol = 1.0e-14;
    return PetscPowReal(vol, 0.6666667);
}

/**
 * @brief IEM decay factor exp(-omega dt) over one step, omega = C D / delta2.
 */
static inline PetscReal IEMDecayFactor(PetscReal C_model, PetscReal diffusivity, PetscReal delta2, PetscReal dt)
{
    PetscReal omega = C_model * diffusivity / delta2;
    return PetscExpReal(-omega * dt);
}

/**
 * @brief Implementation of \ref UpdateParticleField().
 * @see UpdateParticleField()
 */
PetscErrorCode UpdateParticleField(ParticleFieldId field_id,
                                   PetscReal dt,
                                   PetscReal *psi_io,
                                   PetscReal diffusivity,
                                   PetscReal mean_val,
                                   PetscReal delta2,
                                   PetscReal C_model)
{
    PetscFunctionBeginUser;

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        // The IEM model relaxes Psi exponentially toward its target, mean_val, over dt.
        PetscReal decay = IEMDecayFactor(C_model, diffusivity, delta2, dt);

        PetscReal psi_old = *psi_io;
        *psi_io = mean_val + (psi_old - mean_val) * decay;
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Per-cell quantities of the bounded, conservative IEM target.
 * @details The target of particle p in cell c is t_p = m_c + zeta_c g_c . (x_p - x_c), with
 *          - m_c and x_c the cell's Psi mean and particle centroid, each weighted by the particle's
 *            relaxation weight c_p = 1 - exp(-omega_p dt). With omega_p the particle's own rate, the
 *            cell's Psi total changes by sum_p c_p (t_p - psi_p), which these weights make zero for
 *            any zeta_c (Jenny et al. 2001, Eq. 48; Viswanathan, Wang & Pope 2011, Eqs. 20 and 29);
 *          - g_c the gradient of the plain cell means (`lslope`); any slope keeps the total;
 *          - zeta_c in [0, 1] the largest factor that keeps every target of the cell's particles
 *            within the extrema of all particles' Psi at the start of the step (Viswanathan et al.
 *            2011, Eqs. 17 and 28), so the update, a convex combination of psi_p and t_p, cannot
 *            leave them.
 *          A cell's particles may sit on a neighbouring rank, so the sums and the offset extrema are
 *          accumulated into ghosted copies and reduced into the owner (ADD, then MAX through the
 *          DMDA scatter), and the owner's results are read back through ghosted copies. A cell whose
 *          weights sum to zero (every particle with zero rate) keeps its plain mean; its particles do
 *          not change. Collective.
 * @param[in]  user         Context; `lPsi`, `lAj`, `ParticleCount` current.
 * @param[in]  dt           Step.
 * @param[in]  C_model      IEM constant.
 * @param[in]  lslope       Ghosted gradient of the plain cell means, its ghosts current.
 * @param[out] lmoments_out Ghosted `fda` vector: (m_c, zeta_c, 0) per cell; restore with DMRestoreLocalVector.
 * @param[out] lcentroid_out Ghosted `fda` vector: x_c per cell; restore with DMRestoreLocalVector.
 */
static PetscErrorCode ComputeIEMTargetMoments(UserCtx *user, PetscReal dt, PetscReal C_model, Vec lslope,
                                              Vec *lmoments_out, Vec *lcentroid_out)
{
    DM               swarm = user->swarm, fda = user->fda;
    DMDALocalInfo    info = user->info;
    PetscInt         n_local = 0, gxs, gys, gzs, gxm, gym, gzm;
    const PetscReal *psi = NULL, *diff = NULL, *pos = NULL;
    const PetscInt  *cell = NULL;
    PetscReal        bounds[2] = {PETSC_MAX_REAL, PETSC_MAX_REAL}; /* min psi, -max psi */
    Vec              lw = NULL, lwx = NULL, gw = NULL, gwx = NULL, gmom = NULL, gcen = NULL, lext = NULL, gext = NULL;
    Vec              lmom = NULL, lcen = NULL;
    VecScatter       gtol = NULL;
    const PetscReal ***aj = NULL, ***mean = NULL, ***count = NULL;
    Cmpnts        ***w = NULL, ***wx = NULL, ***mom = NULL, ***cen = NULL, ***ext = NULL;
    const Cmpnts  ***slope = NULL, ***lcen_arr = NULL, ***cent = NULL;
    /* Owned physical cells: indices 1 to M-2; 0 and M-1 hold boundary values. */
    const PetscInt   lo[3] = {PetscMax(info.xs, 1), PetscMax(info.ys, 1), PetscMax(info.zs, 1)};
    const PetscInt   hi[3] = {PetscMin(info.xs + info.xm, info.mx - 1), PetscMin(info.ys + info.ym, info.my - 1),
                              PetscMin(info.zs + info.zm, info.mz - 1)};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(DMDAGetGhostCorners(fda, &gxs, &gys, &gzs, &gxm, &gym, &gzm));
#define IN_PATCH(i, j, k) ((i) >= gxs && (i) < gxs + gxm && (j) >= gys && (j) < gys + gym && (k) >= gzs && (k) < gzs + gzm)
    PetscCall(DMSwarmGetLocalSize(swarm, &n_local));
    PetscCall(DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PSI), NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void **)&diff));
    PetscCall(DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void **)&cell));

    /* Pass 1: the Psi extrema, and per cell the weight sums (sum c, sum c psi) and sum c x. */
    PetscCall(DMGetLocalVector(fda, &lw));
    PetscCall(DMGetLocalVector(fda, &lwx));
    PetscCall(VecSet(lw, 0.0));
    PetscCall(VecSet(lwx, 0.0));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecGetArray(fda, lw, &w));
    PetscCall(DMDAVecGetArray(fda, lwx, &wx));
    for (PetscInt p = 0; p < n_local; ++p) {
        const PetscInt i = cell[3*p] + 1, j = cell[3*p+1] + 1, k = cell[3*p+2] + 1;
        bounds[0] = PetscMin(bounds[0], psi[p]);
        bounds[1] = PetscMin(bounds[1], -psi[p]);
        if (!IN_PATCH(i, j, k)) continue;
        const PetscReal c = 1.0 - IEMDecayFactor(C_model, diff[p], IEMMixingLength2(aj[k][j][i]), dt);
        w[k][j][i].x  += c;
        w[k][j][i].y  += c * psi[p];
        wx[k][j][i].x += c * pos[3*p];
        wx[k][j][i].y += c * pos[3*p+1];
        wx[k][j][i].z += c * pos[3*p+2];
    }
    PetscCall(DMDAVecRestoreArray(fda, lwx, &wx));
    PetscCall(DMDAVecRestoreArray(fda, lw, &w));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, bounds, 2, MPIU_REAL, MPIU_MIN, PETSC_COMM_WORLD));
    const PetscReal psi_min = bounds[0], psi_max = -bounds[1];

    PetscCall(DMGetGlobalVector(fda, &gw));
    PetscCall(DMGetGlobalVector(fda, &gwx));
    PetscCall(VecSet(gw, 0.0));
    PetscCall(VecSet(gwx, 0.0));
    PetscCall(DMLocalToGlobal(fda, lw, ADD_VALUES, gw));
    PetscCall(DMLocalToGlobal(fda, lwx, ADD_VALUES, gwx));
    PetscCall(DMRestoreLocalVector(fda, &lwx));
    PetscCall(DMRestoreLocalVector(fda, &lw));

    /* The owner forms the weighted mean and centroid. */
    PetscCall(DMGetGlobalVector(fda, &gmom));
    PetscCall(DMGetGlobalVector(fda, &gcen));
    PetscCall(VecSet(gmom, 0.0));
    PetscCall(VecSet(gcen, 0.0));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lPsi, (void *)&mean));
    PetscCall(DMDAVecGetArrayRead(user->da, user->ParticleCount, (void *)&count));
    PetscCall(DMDAVecGetArrayRead(fda, user->lCent, (void *)&cent));
    PetscCall(DMDAVecGetArrayRead(fda, gw, (void *)&w));
    PetscCall(DMDAVecGetArrayRead(fda, gwx, (void *)&wx));
    PetscCall(DMDAVecGetArray(fda, gmom, &mom));
    PetscCall(DMDAVecGetArray(fda, gcen, &cen));
    for (PetscInt k = lo[2]; k < hi[2]; ++k) {
        for (PetscInt j = lo[1]; j < hi[1]; ++j) {
            for (PetscInt i = lo[0]; i < hi[0]; ++i) {
                const PetscReal sw = w[k][j][i].x;
                if (count[k][j][i] > 0.5 && sw > 0.0) {
                    mom[k][j][i].x = w[k][j][i].y / sw;
                    cen[k][j][i].x = wx[k][j][i].x / sw;
                    cen[k][j][i].y = wx[k][j][i].y / sw;
                    cen[k][j][i].z = wx[k][j][i].z / sw;
                } else {
                    mom[k][j][i].x = mean[k][j][i];
                    cen[k][j][i]   = cent[k][j][i];
                }
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(fda, gcen, &cen));
    PetscCall(DMDAVecRestoreArray(fda, gmom, &mom));
    PetscCall(DMDAVecRestoreArrayRead(fda, gwx, (void *)&wx));
    PetscCall(DMDAVecRestoreArrayRead(fda, gw, (void *)&w));
    PetscCall(DMDAVecRestoreArrayRead(fda, user->lCent, (void *)&cent));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->ParticleCount, (void *)&count));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lPsi, (void *)&mean));
    PetscCall(DMRestoreGlobalVector(fda, &gwx));
    PetscCall(DMRestoreGlobalVector(fda, &gw));
    PetscCall(DMGetLocalVector(fda, &lcen));
    PetscCall(DMGlobalToLocal(fda, gcen, INSERT_VALUES, lcen));
    PetscCall(DMGetLocalVector(fda, &lmom));
    PetscCall(DMGlobalToLocal(fda, gmom, INSERT_VALUES, lmom));

    /* Pass 2: per cell the largest offset g . (x - x_c) above and below zero. Both are >= 0 over a
       cell's particles, since x_c is their weighted centroid, so an untouched 0 never wins the MAX. */
    PetscCall(DMGetLocalVector(fda, &lext));
    PetscCall(VecSet(lext, 0.0));
    PetscCall(DMDAVecGetArrayRead(fda, lslope, (void *)&slope));
    PetscCall(DMDAVecGetArrayRead(fda, lcen, (void *)&lcen_arr));
    PetscCall(DMDAVecGetArray(fda, lext, &ext));
    for (PetscInt p = 0; p < n_local; ++p) {
        const PetscInt i = cell[3*p] + 1, j = cell[3*p+1] + 1, k = cell[3*p+2] + 1;
        if (!IN_PATCH(i, j, k)) continue;
        const Cmpnts    g = slope[k][j][i], xc = lcen_arr[k][j][i];
        const PetscReal d = g.x * (pos[3*p] - xc.x) + g.y * (pos[3*p+1] - xc.y) + g.z * (pos[3*p+2] - xc.z);
        ext[k][j][i].x = PetscMax(ext[k][j][i].x, d);
        ext[k][j][i].y = PetscMax(ext[k][j][i].y, -d);
    }
    PetscCall(DMDAVecRestoreArray(fda, lext, &ext));
    PetscCall(DMDAVecRestoreArrayRead(fda, lcen, (void *)&lcen_arr));
    PetscCall(DMDAVecRestoreArrayRead(fda, lslope, (void *)&slope));
    /* DMLocalToGlobal takes only ADD and INSERT on a DMDA; its own scatter, reversed, takes MAX. */
    PetscCall(DMGetGlobalVector(fda, &gext));
    PetscCall(VecSet(gext, 0.0));
    PetscCall(DMDAGetScatter(fda, &gtol, NULL));
    PetscCall(VecScatterBegin(gtol, lext, gext, MAX_VALUES, SCATTER_REVERSE));
    PetscCall(VecScatterEnd(gtol, lext, gext, MAX_VALUES, SCATTER_REVERSE));
    PetscCall(DMRestoreLocalVector(fda, &lext));

    /* The owner sets zeta: the largest factor in [0, 1] keeping m_c + zeta d within [psi_min, psi_max]. */
    PetscCall(DMDAVecGetArrayRead(fda, gext, (void *)&ext));
    PetscCall(DMDAVecGetArray(fda, gmom, &mom));
    for (PetscInt k = lo[2]; k < hi[2]; ++k) {
        for (PetscInt j = lo[1]; j < hi[1]; ++j) {
            for (PetscInt i = lo[0]; i < hi[0]; ++i) {
                const PetscReal m = mom[k][j][i].x;
                PetscReal       zeta = 1.0;
                if (ext[k][j][i].x > 0.0) zeta = PetscMin(zeta, (psi_max - m) / ext[k][j][i].x);
                if (ext[k][j][i].y > 0.0) zeta = PetscMin(zeta, (m - psi_min) / ext[k][j][i].y);
                mom[k][j][i].y = PetscMax(zeta, 0.0);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(fda, gmom, &mom));
    PetscCall(DMDAVecRestoreArrayRead(fda, gext, (void *)&ext));
    PetscCall(DMRestoreGlobalVector(fda, &gext));
    PetscCall(DMGlobalToLocal(fda, gmom, INSERT_VALUES, lmom));
    PetscCall(DMRestoreGlobalVector(fda, &gcen));
    PetscCall(DMRestoreGlobalVector(fda, &gmom));
#undef IN_PATCH

    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void **)&diff));
    PetscCall(DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PSI), NULL, NULL, (void **)&psi));
    *lmoments_out  = lmom;
    *lcentroid_out = lcen;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref UpdateFieldForAllParticles().
 * @see UpdateFieldForAllParticles()
 */
PetscErrorCode UpdateFieldForAllParticles(UserCtx *user, ParticleFieldId field_id)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    DM             da = user->da;
    PetscInt       n_local;
    PetscReal      dt = user->simCtx->dt;
    PetscReal      C_IEM = user->simCtx->iem_constant;

    PetscReal       *psi_arr = NULL;
    PetscReal       *diff_arr = NULL;
    PetscReal       *pos_arr = NULL;
    PetscInt        *cell_arr = NULL;
    PetscInt64      *pid_arr = NULL;

    PetscReal       ***grid_aj = NULL;
    Cmpnts          ***grid_slope = NULL, ***grid_centroid = NULL, ***grid_moments = NULL;
    Vec             lslope = NULL, lcentroid = NULL, lmoments = NULL;

    PetscBool       accessed_eulerian = PETSC_FALSE;
    const ParticleFieldDescriptor *descriptor = NULL;
    const char      *fieldName = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = ParticleFieldGetDescriptor(field_id, &descriptor); CHKERRQ(ierr);
    PetscCheck((descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_MODEL_UPDATE) != 0,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Particle field '%s' has no registered model-update kernel.",
               descriptor->canonical_name);
    fieldName = descriptor->canonical_name;

    // Relax toward the mean of the particles where they are now. The Eulerian Psi otherwise
    // holds the previous step's scatter, taken before the particles moved and migrated; the
    // scatter recounts the particles per cell itself. Collective, so it comes before the early
    // return of ranks with no particles.
    if (field_id == PARTICLE_FIELD_ID_PSI) {
        Vec slope = NULL;

        // Psi relaxation requires ghosted Eulerian mean and Jacobian fields.
        PetscCheck(user->lPsi && user->lAj, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx lPsi or lAj not initialized.");
        ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
        // The target is the cell's limited linear reconstruction m + zeta g . (x - xc), with m and
        // xc weighted by each particle's relaxation weight so each cell's total is kept, and zeta
        // keeping every target within the particles' Psi extrema (ComputeIEMTargetMoments).
        // Relaxing to the stepped mean alone adds transport that does not vanish with refinement.
        // A particle's cell is stored one index past its origin node, so the cell fields are read
        // through ghosted copies.
        // The scatter's count holds only this rank's particles in its ghost cells; the empty-cell
        // test reads neighbours across subdomain seams, so their totals are fetched first.
        ierr = UpdateLocalGhosts(user, FIELD_ID_PARTICLE_COUNT); CHKERRQ(ierr);
        ierr = DMGetGlobalVector(user->fda, &slope); CHKERRQ(ierr);
        ierr = ComputeCellScalarGradient(user, user->lPsi, user->lParticleCount, slope); CHKERRQ(ierr);
        ierr = DMGetLocalVector(user->fda, &lslope); CHKERRQ(ierr);
        ierr = DMGlobalToLocal(user->fda, slope, INSERT_VALUES, lslope); CHKERRQ(ierr);
        ierr = DMRestoreGlobalVector(user->fda, &slope); CHKERRQ(ierr);
        ierr = ComputeIEMTargetMoments(user, dt, C_IEM, lslope, &lmoments, &lcentroid); CHKERRQ(ierr);
    }

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        if (lslope) { ierr = DMRestoreLocalVector(user->fda, &lslope); CHKERRQ(ierr); }
        if (lcentroid) { ierr = DMRestoreLocalVector(user->fda, &lcentroid); CHKERRQ(ierr); }
        if (lmoments) { ierr = DMRestoreLocalVector(user->fda, &lmoments); CHKERRQ(ierr); }
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void**)&pid_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);

        ierr = DMDAVecGetArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, lslope, &grid_slope); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, lcentroid, &grid_centroid); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, lmoments, &grid_moments); CHKERRQ(ierr);
        accessed_eulerian = PETSC_TRUE;
    }

    // The filter width depends only on the cell, and the particles are stored in cell order,
    // so it is recomputed only when the cell changes.
    PetscInt  width_cell[3] = {-1, -1, -1};
    PetscReal p_delta2      = 1.0;
    for (PetscInt p = 0; p < n_local; ++p) {
        PetscReal p_diff = 0.0;
        PetscReal p_mean = 0.0;

        if (field_id == PARTICLE_FIELD_ID_PSI) {
            // Cell IDs are 0-based geometric indices; cell-centred storage starts at index 1
            // (index 0 is the boundary slot), the shift the scatter applies when it writes lPsi.
            PetscInt i = cell_arr[3*p + 0] + 1;
            PetscInt j = cell_arr[3*p + 1] + 1;
            PetscInt k = cell_arr[3*p + 2] + 1;

            const Cmpnts g = grid_slope[k][j][i], c = grid_centroid[k][j][i], mz = grid_moments[k][j][i];
            p_diff = diff_arr[p];
            p_mean = mz.x + mz.y * (g.x * (pos_arr[3*p + 0] - c.x) + g.y * (pos_arr[3*p + 1] - c.y)
                                    + g.z * (pos_arr[3*p + 2] - c.z));

            if (i != width_cell[0] || j != width_cell[1] || k != width_cell[2]) {
                // Aj is the reciprocal cell volume in the curvilinear-grid representation.
                p_delta2 = IEMMixingLength2(grid_aj[k][j][i]);
                width_cell[0] = i; width_cell[1] = j; width_cell[2] = k;
            }
        }

        const PetscReal psi_old = psi_arr[p];
        ierr = UpdateParticleField(field_id, dt, &psi_arr[p], p_diff, p_mean, p_delta2, C_IEM);
        CHKERRQ(ierr);
        LOG_PARTICLE(user->simCtx, pid_arr[p], "mix",
                     "%s %.10e -> %.10e | target %.10e | diffusivity %.6e mixing length^2 %.6e C %.4f dt %.6e\n",
                     fieldName, psi_old, psi_arr[p], p_mean, p_diff, p_delta2, C_IEM, dt);
    }

    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void**)&pid_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);
    }

    if (accessed_eulerian) {
        ierr = DMDAVecRestoreArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, lslope, &grid_slope); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, lcentroid, &grid_centroid); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, lmoments, &grid_moments); CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(user->fda, &lslope); CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(user->fda, &lcentroid); CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(user->fda, &lmoments); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updated particle physics for field '%s'.\n", fieldName);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref UpdateAllParticleFields().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticlePhysics.h`.
 * @see UpdateAllParticleFields()
 */
PetscErrorCode UpdateAllParticleFields(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updating all particle physical properties...\n");

    if (VerificationScalarOverrideActive(user->simCtx)) {
        // Verification profiles define Psi exactly, so bypass the model-driven update.
        ierr = ApplyVerificationScalarOverrideToParticles(user); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO, "Verification scalar override active; skipped model-driven Psi update.\n");
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    /* C_IEM = 0 switches micromixing off. The update is skipped rather than run with a
     * zero rate: mean + (psi - mean) * 1 need not round-trip to psi exactly, and a passive
     * label must stay bit-identical. */
    if (user->simCtx->iem_constant == 0.0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "IEM constant is zero; Psi is carried unchanged.\n");
    } else {
        ierr = UpdateFieldForAllParticles(user, PARTICLE_FIELD_ID_PSI); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "All particle physical properties updated.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ComputeEulerianDiffusivity().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticlePhysics.h`.
 * @see ComputeEulerianDiffusivity()
 */
PetscErrorCode ComputeEulerianDiffusivity(UserCtx *user)
{
    const SimCtx       *simCtx = user->simCtx;
    const DMDALocalInfo info = user->info;
    // Grid uses indices 1 to M-2 for physical cells; index 0 and M-1 are boundary holders.
    const PetscInt      lo[3] = {PetscMax(info.xs, 1), PetscMax(info.ys, 1), PetscMax(info.zs, 1)};
    const PetscInt      hi[3] = {PetscMin(info.xs + info.xm, info.mx - 1), PetscMin(info.ys + info.ym, info.my - 1),
                                 PetscMin(info.zs + info.zm, info.mz - 1)};
    PetscReal        ***diff_arr = NULL;
    const PetscReal  ***nut_arr = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (VerificationDiffusivityOverrideActive(user->simCtx)) {
        PetscCall(ApplyVerificationDiffusivityOverride(user));
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // Gamma = nu / Sc + nu_t / Sc_t, with nu = 1/Re (0 when Re is unset or not positive) and the
    // Schmidt numbers 1 and 0.7 when unset. The eddy viscosity enters only with LES; a dynamic
    // model can give a slightly negative one, which is clipped to 0 so the Brownian amplitude
    // sqrt(2 Gamma dt) stays real.
    const PetscReal nu_molecular    = (simCtx->ren > 1.0e-12) ? 1.0 / simCtx->ren : 0.0;
    const PetscReal Sc              = (simCtx->schmidt_number > 1.0e-6) ? simCtx->schmidt_number : 1.0;
    const PetscReal Sct             = (simCtx->Turbulent_schmidt_number > 1.0e-6) ? simCtx->Turbulent_schmidt_number : 0.7;
    const PetscReal gamma_molecular = nu_molecular / Sc;
    const PetscBool use_turbulence_model = simCtx->les ? PETSC_TRUE : PETSC_FALSE;

    PetscCall(DMDAVecGetArray(user->da, user->Diffusivity, &diff_arr));
    if (use_turbulence_model) PetscCall(DMDAVecGetArrayRead(user->da, user->Nu_t, (void *)&nut_arr));
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                PetscReal gamma_turbulent = 0.0;
                if (use_turbulence_model) gamma_turbulent = PetscMax(nut_arr[k][j][i], 0.0) / Sct;
                diff_arr[k][j][i] = gamma_molecular + gamma_turbulent;
            }
        }
    }
    if (use_turbulence_model) PetscCall(DMDAVecRestoreArrayRead(user->da, user->Nu_t, (void *)&nut_arr));
    PetscCall(DMDAVecRestoreArray(user->da, user->Diffusivity, &diff_arr));

    // The diffusivity gradient and the particle interpolation read the periodic dummy cells and
    // the ghosts.
    const FieldId periodic_fields[] = {FIELD_ID_DIFFUSIVITY};
    PetscCall(SynchronizePeriodicCellFields(user, 1, periodic_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Cartesian gradient of a cell-centred scalar at every owned physical cell.
 * @details Each computational derivative is the central difference, or at a non-periodic
 *          boundary cell the second-order one-sided difference into the domain, mapped to
 *          Cartesian components with the cell's metrics. With `lcount`, a derivative whose
 *          stencil reaches a cell holding no particles is zero, so an empty cell's zero mean
 *          does not enter. Writes owned entries only; local to this translation unit.
 * @param[in]  user     Grid context.
 * @param[in]  lscalar  Ghosted scalar, its ghosts current.
 * @param[in]  lcount   Ghosted particle count, its ghosts current; or NULL.
 * @param[out] gradient Global vector on `fda`.
 */
static PetscErrorCode ComputeCellScalarGradient(UserCtx *user, Vec lscalar, Vec lcount, Vec gradient)
{
    DMDALocalInfo    info = user->info;
    const PetscInt   n[3] = {info.mx, info.my, info.mz};
    PetscBool        periodic[3];
    GetPeriodicAxes(user->simCtx, periodic);
    // Grid uses indices 1 to M-2 for physical cells; index 0 and M-1 are boundary holders.
    const PetscInt   lo[3] = {PetscMax(info.xs, 1), PetscMax(info.ys, 1), PetscMax(info.zs, 1)};
    const PetscInt   hi[3] = {PetscMin(info.xs + info.xm, n[0] - 1), PetscMin(info.ys + info.ym, n[1] - 1),
                              PetscMin(info.zs + info.zm, n[2] - 1)};
    const PetscReal  ***f = NULL, ***count = NULL, ***aj = NULL;
    const Cmpnts     ***csi = NULL, ***eta = NULL, ***zet = NULL;
    Cmpnts           ***grad = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(DMDAVecGetArrayRead(user->da, lscalar, (void *)&f));
    if (lcount) PetscCall(DMDAVecGetArrayRead(user->da, lcount, (void *)&count));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (void *)&csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecGetArray(user->fda, gradient, &grad));
#define CELL(a, d, o) ((d) == 0 ? (a)[k][j][i + (o)] : (d) == 1 ? (a)[k][j + (o)][i] : (a)[k + (o)][j][i])
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                const PetscInt idx[3] = {i, j, k};
                PetscReal      dG[3];
                for (PetscInt d = 0; d < 3; ++d) {
                    PetscInt  o[3] = {-1, 1, 0};
                    PetscReal c[3] = {-1.0, 1.0, 0.0};   /* central: [-1, 0, 1] / 2 */
                    if (!periodic[d] && idx[d] == 1) {
                        o[0] = 0; o[1] = 1; o[2] = 2; c[0] = -3.0; c[1] = 4.0; c[2] = -1.0;   /* forward: [-3, 4, -1] / 2 */
                    } else if (!periodic[d] && idx[d] == n[d] - 2) {
                        o[0] = 0; o[1] = -1; o[2] = -2; c[0] = 3.0; c[1] = -4.0; c[2] = 1.0;  /* backward: [1, -4, 3] / 2 */
                    }
                    dG[d] = (c[0] * CELL(f, d, o[0]) + c[1] * CELL(f, d, o[1]) + c[2] * CELL(f, d, o[2])) * 0.5;
                    if (count) {
                        for (PetscInt q = 0; q < 3; ++q) if (c[q] != 0.0 && CELL(count, d, o[q]) < 0.5) dG[d] = 0.0;
                    }
                }
                TransformScalarDerivativesToPhysical(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i],
                                                     dG[0], dG[1], dG[2], &grad[k][j][i]);
            }
        }
    }
#undef CELL
    PetscCall(DMDAVecRestoreArray(user->fda, gradient, &grad));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (void *)&csi));
    if (lcount) PetscCall(DMDAVecRestoreArrayRead(user->da, lcount, (void *)&count));
    PetscCall(DMDAVecRestoreArrayRead(user->da, lscalar, (void *)&f));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ComputeEulerianDiffusivityGradient().
 * @see ComputeEulerianDiffusivityGradient()
 */
PetscErrorCode ComputeEulerianDiffusivityGradient(UserCtx *user)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // Central differences read the neighbours across subdomain and periodic seams.
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY));
    PetscCall(ComputeCellScalarGradient(user, user->lDiffusivity, NULL, user->DiffusivityGradient));
    // A particle within half a cell of a periodic seam interpolates across it, reading the
    // periodic dummy cells, which the gradient's owned-cell loop does not write.
    const FieldId periodic_fields[] = {FIELD_ID_DIFFUSIVITY_GRADIENT};
    PetscCall(SynchronizePeriodicCellFields(user, 1, periodic_fields));
    // Particles near subdomain boundaries interpolate the gradient, so its ghosts are filled.
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY_GRADIENT));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
