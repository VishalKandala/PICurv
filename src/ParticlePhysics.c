#include "ParticlePhysics.h"
#include "verification_sources.h"
#include "ParticleMotion.h"
#include "interpolation.h"

#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256
#endif

static PetscErrorCode ComputeCellScalarGradient(UserCtx *user, Vec lscalar, Vec lcount, Vec gradient);

/**
 * @brief Internal helper implementation: `UpdateParticleField()`.
 * @details Local to this translation unit.
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
        PetscReal omega = C_model * diffusivity / delta2;
        PetscReal decay = PetscExpReal(-omega * dt);

        PetscReal psi_old = *psi_io;
        *psi_io = mean_val + (psi_old - mean_val) * decay;
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `UpdateFieldForAllParticles()`.
 * @details Local to this translation unit.
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

    PetscReal       ***grid_mean = NULL;
    PetscReal       ***grid_aj = NULL;
    Cmpnts          ***grid_slope = NULL, ***grid_centroid = NULL;
    Vec             lslope = NULL, lcentroid = NULL;

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
        Vec slope = NULL, centroid = NULL;

        // Psi relaxation requires ghosted Eulerian mean and Jacobian fields.
        PetscCheck(user->lPsi && user->lAj, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx lPsi or lAj not initialized.");
        ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
        // The target is the cell's linear reconstruction m + g . (x - xc), with xc the cell's
        // particle centroid so each cell mean is kept. Relaxing to the stepped mean m alone
        // adds transport that does not vanish with refinement. A particle's cell is stored one
        // index past its origin node, so both fields are read through ghosted copies.
        // The scatter's count holds only this rank's particles in its ghost cells; the empty-cell
        // test reads neighbours across subdomain seams, so their totals are fetched first.
        ierr = UpdateLocalGhosts(user, FIELD_ID_PARTICLE_COUNT); CHKERRQ(ierr);
        ierr = DMGetGlobalVector(user->fda, &slope); CHKERRQ(ierr);
        ierr = ComputeCellScalarGradient(user, user->lPsi, user->lParticleCount, slope); CHKERRQ(ierr);
        ierr = DMGetLocalVector(user->fda, &lslope); CHKERRQ(ierr);
        ierr = DMGlobalToLocal(user->fda, slope, INSERT_VALUES, lslope); CHKERRQ(ierr);
        // The scatter above has just counted the particles, so the centroid is not recounted.
        ierr = DMGetGlobalVector(user->fda, &centroid); CHKERRQ(ierr);
        ierr = ScatterParticleFieldCellMean(user, PARTICLE_FIELD_ID_POSITION, user->fda, user->ParticleCount, PETSC_FALSE,
                                            centroid); CHKERRQ(ierr);
        ierr = DMGetLocalVector(user->fda, &lcentroid); CHKERRQ(ierr);
        ierr = DMGlobalToLocal(user->fda, centroid, INSERT_VALUES, lcentroid); CHKERRQ(ierr);
        ierr = DMRestoreGlobalVector(user->fda, &centroid); CHKERRQ(ierr);
        ierr = DMRestoreGlobalVector(user->fda, &slope); CHKERRQ(ierr);
    }

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        if (lslope) { ierr = DMRestoreLocalVector(user->fda, &lslope); CHKERRQ(ierr); }
        if (lcentroid) { ierr = DMRestoreLocalVector(user->fda, &lcentroid); CHKERRQ(ierr); }
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void**)&pid_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);

        ierr = DMDAVecGetArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, lslope, &grid_slope); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, lcentroid, &grid_centroid); CHKERRQ(ierr);
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

            const Cmpnts g = grid_slope[k][j][i], c = grid_centroid[k][j][i];
            p_diff = diff_arr[p];
            p_mean = grid_mean[k][j][i] + g.x * (pos_arr[3*p + 0] - c.x) + g.y * (pos_arr[3*p + 1] - c.y)
                                         + g.z * (pos_arr[3*p + 2] - c.z);

            if (i != width_cell[0] || j != width_cell[1] || k != width_cell[2]) {
                // Aj is the reciprocal cell volume in the curvilinear-grid representation.
                PetscReal jac   = grid_aj[k][j][i];
                PetscReal p_vol = (jac > 1.0e-14) ? (1.0 / jac) : 1.0e-14;
                // Guard the LES mixing time scale against degenerate or cut-cell volumes.
                if (p_vol < 1.0e-14) p_vol = 1.0e-14;
                p_delta2 = PetscPowReal(p_vol, 0.6666667);
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
        ierr = DMDAVecRestoreArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, lslope, &grid_slope); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, lcentroid, &grid_centroid); CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(user->fda, &lslope); CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(user->fda, &lcentroid); CHKERRQ(ierr);
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
    PetscErrorCode ierr;
    DM             da = user->da;
    PetscInt       i, j, k, xs, ys, zs, xm, ym, zm, xe, ye, ze;
	PetscInt       lxs, lys, lzs, lxe, lye, lze;

    // Pointers for 3D grid access
    PetscReal      ***diff_arr; // Output: Diffusivity field
    PetscReal      ***nut_arr;  // Input:  Eddy Viscosity field (optional)

    // Physics parameters
    PetscReal      nu_molecular, gamma_molecular;
    PetscReal      nu_turbulent, gamma_turbulent;
    PetscReal      Sc, Sct;
    PetscBool      use_turbulence_model;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (VerificationDiffusivityOverrideActive(user->simCtx)) {
        ierr = ApplyVerificationDiffusivityOverride(user); CHKERRQ(ierr);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // ------------------------------------------------------------------------
    // 1. Parameter Setup & Safety Checks
    // ------------------------------------------------------------------------
    
    // Determine Molecular Viscosity (nu = 1/Re)
    // Guard against division by zero if Re is not set or infinite (inviscid)
    if (user->simCtx->ren > 1.0e-12) {
        nu_molecular = 1.0 / user->simCtx->ren;
    } else {
        nu_molecular = 0.0;
    }

    // Set Schmidt Numbers (Default to 1.0 if not provided to prevent NaN)
    Sc  = (user->simCtx->schmidt_number > 1.0e-6) ? user->simCtx->schmidt_number : 1.0;
    Sct = (user->simCtx->Turbulent_schmidt_number > 1.0e-6) ? user->simCtx->Turbulent_schmidt_number : 0.7;

    // Pre-calculate molecular component
    gamma_molecular = nu_molecular / Sc;

    // Check if a turbulence model is active
    use_turbulence_model = user->simCtx->les ? PETSC_TRUE : PETSC_FALSE;

    // ------------------------------------------------------------------------
    // 2. Data Access
    // ------------------------------------------------------------------------
    
    // Get local grid boundaries
	DMDALocalInfo info;
    info = user->info;

	xs = info.xs; ys = info.ys; zs = info.zs;
	xm = info.xm; ym = info.ym; zm = info.zm;
	xe = xs + xm; ye = ys + ym; ze = zs + zm;

	lxs = (xs == 0)? xs + 1 : xs;
	lys = (ys == 0)? ys + 1 : ys;
	lzs = (zs == 0)? zs + 1 : zs;
	lxe = (xe == info.mx)? xe - 1 : xe;
	lye = (ye == info.my)? ye - 1 : ye;
	lze = (ze == info.mz)? ze - 1 : ze;

    // Get write access to the output Diffusivity array
    ierr = DMDAVecGetArray(da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);

    // Get read access to Eddy Viscosity only if turbulence is active
    if (use_turbulence_model) {
        ierr = DMDAVecGetArrayRead(da, user->Nu_t, &nut_arr); CHKERRQ(ierr);
    }

    // ------------------------------------------------------------------------
    // 3. Calculation Loop
    // ------------------------------------------------------------------------

    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                
                gamma_turbulent = 0.0;

                if (use_turbulence_model) {
                    // Fetch local eddy viscosity
                    nu_turbulent = nut_arr[k][j][i];

                    // NUMERICAL SAFETY: 
                    // Some turbulence models (dynamic SGS) can locally produce 
                    // slightly negative viscosity. We clamp this to 0 to prevent
                    // negative diffusivity, which crashes the Langevin sqrt().
                    if (nu_turbulent < 0.0) {
                        nu_turbulent = 0.0;
                    }

                    gamma_turbulent = nu_turbulent / Sct;
                }

                // Sum components
                diff_arr[k][j][i] = gamma_molecular + gamma_turbulent;
            }
        }
    }

    // ------------------------------------------------------------------------
    // 4. Cleanup & Synchronization
    // ------------------------------------------------------------------------

    // Restore arrays
    if (use_turbulence_model) {
        ierr = DMDAVecRestoreArrayRead(da, user->Nu_t, &nut_arr); CHKERRQ(ierr);
    }
    ierr = DMDAVecRestoreArray(da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);

    // Update Ghost Points
    // This is required because downstream operations (Drift Gradient Calculation 
    // and Particle Interpolation) will need access to the halo regions of this field.
    const FieldId periodic_fields[] = {FIELD_ID_DIFFUSIVITY};
    ierr = SynchronizePeriodicCellFields(user, 1, periodic_fields); CHKERRQ(ierr);
	ierr = UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY); CHKERRQ(ierr);
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
    const PetscBool  periodic[3] = {(PetscBool)(user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC),
                                    (PetscBool)(user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC),
                                    (PetscBool)(user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC)};
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
 * @brief Internal helper implementation: `ComputeEulerianDiffusivityGradient()`.
 * @details Local to this translation unit.
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
