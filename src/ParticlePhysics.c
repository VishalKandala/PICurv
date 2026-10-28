#include "ParticlePhysics.h"
#include "verification_sources.h"
#include "ParticleMotion.h"
#include "interpolation.h"

#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256
#endif

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
        // The IEM model relaxes Psi exponentially toward the cell mean over dt.
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
    PetscInt        *cell_arr = NULL;

    PetscReal       ***grid_mean = NULL;
    PetscReal       ***grid_aj = NULL;

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
        ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
    }

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);

        // Psi relaxation requires ghosted Eulerian mean and Jacobian fields.
        if (!user->lPsi || !user->lAj) {
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx lPsi or lAj not initialized.");
        }

        ierr = DMDAVecGetArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
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

            p_diff = diff_arr[p];
            p_mean = grid_mean[k][j][i];

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

        ierr = UpdateParticleField(field_id, dt, &psi_arr[p], p_diff, p_mean, p_delta2, C_IEM);
        CHKERRQ(ierr);
    }

    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);
    }

    if (accessed_eulerian) {
        ierr = DMDAVecRestoreArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updated particle physics for field '%s'.\n", fieldName);

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
 * @brief Internal helper implementation: `ComputeEulerianDiffusivityGradient()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeEulerianDiffusivityGradient(UserCtx *user)
{
    PetscErrorCode ierr;
    DM             da = user->da, fda = user->fda;
    DMDALocalInfo  info = user->info;

    // 1. Determine Global Dimensions
    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    // 2. Determine Local Loop Bounds (Skip Ghosts/Unused Indices)
    //    Grid uses indices 1 to M-2 for physical cells.
    //    Index 0 and M-1 are ghost/boundary holders.
    
    // Start: If we own the global start (0), skip it and start at 1.
    PetscInt lxs = (info.xs == 0) ? 1 : info.xs;
    // End:   If we own the global end (mx), stop before it (mx-1), so loop covers mx-2.
    PetscInt lxe = (info.xs + info.xm == mx) ? mx - 1 : info.xs + info.xm;

    PetscInt lys = (info.ys == 0) ? 1 : info.ys;
    PetscInt lye = (info.ys + info.ym == my) ? my - 1 : info.ys + info.ym;

    PetscInt lzs = (info.zs == 0) ? 1 : info.zs;
    PetscInt lze = (info.zs + info.zm == mz) ? mz - 1 : info.zs + info.zm;

    PetscInt i, j, k;

    // Pointers
    PetscReal ***diff;      // Input (Scalar)
    Cmpnts    ***grad_diff; // Output (Vector)
    Cmpnts    ***csi, ***eta, ***zet;
    PetscReal ***aj;

    // Boundary Flags (Check if Periodic)
    PetscBool p_x = (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC);
    PetscBool p_y = (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC);
    PetscBool p_z = (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC);

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // 3. Update Ghosts for Diffusivity 
    //    Required so that Central Differences at i=2 can read i=1,
    //    and Central Differences at Periodic Boundaries work correctly.
    ierr = UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY); CHKERRQ(ierr);

    // 4. Get Arrays (Read Only)
    ierr = DMDAVecGetArrayRead(da,  user->lDiffusivity, &diff); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lCsi, &csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lEta, &eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lZet, &zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da,  user->lAj,  &aj);  CHKERRQ(ierr);
    
	LOG_ALLOW(GLOBAL,LOG_DEBUG,"Diffusivity and Metrics arrays accessed successfully! .\n");

    // 5. Get Output Array (Read/Write)
    ierr = DMDAVecGetArray(fda, user->DiffusivityGradient, &grad_diff); CHKERRQ(ierr);

	LOG_ALLOW(GLOBAL,LOG_DEBUG,"Diffusivity Gradient array accessed successfully! .\n");

    // 6. Loop over Physical Domain
    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                
                PetscReal dGdCsi, dGdEta, dGdZet;

                // ---------------------------------------------------------
                // I-Direction (Csi)
                // ---------------------------------------------------------
                if (!p_x && i == 1) {
                    // Physical Start: 2nd Order Forward Difference
                    // Stencil: [-3, 4, -1] / 2
                    dGdCsi = (-3.0 * diff[k][j][i] + 4.0 * diff[k][j][i+1] - diff[k][j][i+2]) * 0.5;
                } 
                else if (!p_x && i == mx - 2) {
                    // Physical End: 2nd Order Backward Difference
                    // Stencil: [1, -4, 3] / 2
                    dGdCsi = (3.0 * diff[k][j][i] - 4.0 * diff[k][j][i-1] + diff[k][j][i-2]) * 0.5;
                } 
                else {
                    // Interior / Periodic: Central Difference
                    // Stencil: [-1, 0, 1] / 2
                    dGdCsi = (diff[k][j][i+1] - diff[k][j][i-1]) * 0.5;
                }

                // ---------------------------------------------------------
                // J-Direction (Eta)
                // ---------------------------------------------------------
                if (!p_y && j == 1) {
                    // Physical Start: Forward
                    dGdEta = (-3.0 * diff[k][j][i] + 4.0 * diff[k][j+1][i] - diff[k][j+2][i]) * 0.5;
                } 
                else if (!p_y && j == my - 2) {
                    // Physical End: Backward
                    dGdEta = (3.0 * diff[k][j][i] - 4.0 * diff[k][j-1][i] + diff[k][j-2][i]) * 0.5;
                } 
                else {
                    // Interior: Central
                    dGdEta = (diff[k][j+1][i] - diff[k][j-1][i]) * 0.5;
                }

                // ---------------------------------------------------------
                // K-Direction (Zet)
                // ---------------------------------------------------------
                if (!p_z && k == 1) {
                    // Physical Start: Forward
                    dGdZet = (-3.0 * diff[k][j][i] + 4.0 * diff[k+1][j][i] - diff[k+2][j][i]) * 0.5;
                } 
                else if (!p_z && k == mz - 2) {
                    // Physical End: Backward
                    dGdZet = (3.0 * diff[k][j][i] - 4.0 * diff[k-1][j][i] + diff[k-2][j][i]) * 0.5;
                } 
                else {
                    // Interior: Central
                    dGdZet = (diff[k+1][j][i] - diff[k-1][j][i]) * 0.5;
                }

                // ---------------------------------------------------------
                // Transform to Physical Space (Cartesian Gradient)
                // ---------------------------------------------------------
                TransformScalarDerivativesToPhysical(aj[k][j][i],
                                                     csi[k][j][i], eta[k][j][i], zet[k][j][i],
                                                     dGdCsi, dGdEta, dGdZet,
                                                     &grad_diff[k][j][i]);
            }
        }
    }

    // 7. Restore Arrays
    ierr = DMDAVecRestoreArrayRead(da,  user->lDiffusivity, &diff); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fda, user->lCsi, &csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fda, user->lEta, &eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fda, user->lZet, &zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da,  user->lAj,  &aj);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, user->DiffusivityGradient, &grad_diff); CHKERRQ(ierr);

    // 8. Update Ghosts for the Result
    //    Important: Particles near subdomain boundaries will need to interpolate
    //    this gradient vector, so the ghosts of the vector field must be filled.
    ierr = UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY_GRADIENT); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
