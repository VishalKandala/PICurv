#include "BodyForces.h"


//////////////////////////////////////////////
// DRIVEN CHANNEL FLOW FORCE(EQUIVALENT) TERM 
/////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "ComputeDrivenChannelFlowSource"
/**
 * @brief Applies a momentum source term to drive flow in a periodic channel or pipe.
 *
 * This function is the "engine" of the driven flow control system. It operates by:
 * 1.  Introspecting the boundary condition handlers to see if a `DRIVEN_` flow
 *     handler is active on any face. This determines if a driven flow is enabled
 *     and in which direction ('X', 'Y', or 'Z').
 * 2.  If a driven flow is active, it reads the `bulkVelocityCorrection` value that
 *     was computed by the handler's `PreStep` method and stored in the `SimCtx`.
 * 3.  It translates this velocity correction into a momentum source term.
 * 4.  It adds this source term to the appropriate component of the contravariant
 *     RHS vector (`Rct`) for all fluid cells in the domain.
 *
 * If no driven flow handler is found, this function does nothing.
 *
 * @param user The UserCtx containing the simulation state for a single block.
 * @param Rct  The PETSc Vec for the contravariant RHS, which will be modified in-place.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeDrivenChannelFlowSource(UserCtx *user, Vec Rct)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    PetscFunctionBeginUser;

    // --- Step 1: Discover if and where a driven flow is active ---
    char drivenDirection = ' '; // Use space as a null/not-found indicator
    for (int i = 0; i < 6; i++) {
        BCHandlerType handler_type = user->boundary_faces[i].handler_type;
        if (handler_type == BC_HANDLER_DRIVEN_CONSTANT_FLUX ||
            handler_type == BC_HANDLER_DRIVEN_INITIAL_FLUX)
        {
            switch (user->boundary_faces[i].face_id) {
                case BC_FACE_NEG_X: case BC_FACE_POS_X: drivenDirection = 'X'; break;
                case BC_FACE_NEG_Y: case BC_FACE_POS_Y: drivenDirection = 'Y'; break;
                case BC_FACE_NEG_Z: case BC_FACE_POS_Z: drivenDirection = 'Z'; break;
            }
            break; // Found it, no need to check other faces
        }
    }

    // --- Step 2: Early exit if no driven flow is configured ---
    if (drivenDirection == ' ') {
        PetscFunctionReturn(0);
    }

    // --- Step 3: Get the control signal and exit if no correction is needed ---
    PetscReal bulkVelocityCorrection = simCtx->bulkVelocityCorrection;
    if (PetscAbsReal(bulkVelocityCorrection) < 1.0e-12) {
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %d: Applying driven flow momentum source in '%c' direction.\n",
              simCtx->rank, user->_this, drivenDirection);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  - Received Bulk Velocity Correction: %le\n", bulkVelocityCorrection);

    // --- Step 4: Setup for calculation ---
    DMDALocalInfo info = user->info;
    PetscInt i, j, k;
    PetscInt lxs = (info.xs == 0) ? 1 : info.xs;
    PetscInt lys = (info.ys == 0) ? 1 : info.ys;
    PetscInt lzs = (info.zs == 0) ? 1 : info.zs;
    PetscInt lxe = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    PetscInt lye = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    PetscInt lze = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    Cmpnts ***rct, ***csi, ***eta, ***zet;
    PetscReal ***nvert;
    ierr = DMDAVecGetArray(user->fda, Rct, &rct); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // Calculate the driving force magnitude for the current timestep, smoothed
    // with the value from the previous step for stability.
    const PetscReal forceScalingFactor  = simCtx->forceScalingFactor;
    PetscReal drivingForceMagnitude = (bulkVelocityCorrection / simCtx->dt / simCtx->st * COEF_TIME_ACCURACY);
    drivingForceMagnitude = (simCtx->drivingForceMagnitude * 0.5) + (drivingForceMagnitude * 0.5); 
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "  - Previous driving force:            %le\n", simCtx->drivingForceMagnitude);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "  - New smoothed driving force:        %le\n", drivingForceMagnitude);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "  - Force scaling factor:              %f\n",  simCtx->forceScalingFactor);

    PetscBool hasLoggedApplication = PETSC_FALSE; // Flag to log details only once per rank.
    // --- Step 5: Apply the momentum source to the correct RHS component ---
    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                if (nvert[k][j][i] < 0.1) { // Apply only to fluid cells
                    PetscReal faceArea = 0.0;
                    PetscReal momentumSource = 0.0;

                    switch (drivenDirection) {
                        case 'X':
                            faceArea = sqrt(csi[k][j][i].x * csi[k][j][i].x + csi[k][j][i].y * csi[k][j][i].y + csi[k][j][i].z * csi[k][j][i].z);
                            momentumSource = drivingForceMagnitude * forceScalingFactor * faceArea;
                            rct[k][j][i].x += momentumSource;

                            // Log details for the very first point where force is applied on this rank.
                            if (!hasLoggedApplication) {
                                LOG_ALLOW(LOCAL, LOG_DEBUG,"Body Force %le added at (%d,%d,%d)\n",force_z, k, j, i);
                                hasLoggedApplication = PETSC_TRUE;
                            }
                                break;
                        case 'Y':
                            faceArea = sqrt(eta[k][j][i].x * eta[k][j][i].x + eta[k][j][i].y * eta[k][j][i].y + eta[k][j][i].z * eta[k][j][i].z);
                            momentumSource = drivingForceMagnitude * forceScalingFactor * faceArea;
                            rct[k][j][i].y += momentumSource;
                            
                            // Log details for the very first point where force is applied on this rank.
                            if (!hasLoggedApplication) {
                                LOG_ALLOW(LOCAL, LOG_DEBUG,"Body Force %le added at (%d,%d,%d)\n",force_z, k, j, i);
                                hasLoggedApplication = PETSC_TRUE;
                            }
                            break;
                        case 'Z':
                            faceArea = sqrt(zet[k][j][i].x * zet[k][j][i].x + zet[k][j][i].y * zet[k][j][i].y + zet[k][j][i].z * zet[k][j][i].z);
                            momentumSource = drivingForceMagnitude * forceScalingFactor * faceArea;
                            rct[k][j][i].z += momentumSource;

                            // Log details for the very first point where force is applied on this rank.
                            if (!hasLoggedApplication) {
                                LOG_ALLOW(LOCAL, LOG_DEBUG,"Body Force %le added at (%d,%d,%d)\n",force_z, k, j, i);
                                hasLoggedApplication = PETSC_TRUE;
                            }
                            break;
                    }
                }
            }
        }
    }

    // Update simCtx with latest Driving Force
	simCtx->drivingForceMagnitude = drivingForceMagnitude;

    // --- Step 6: Restore arrays ---
    ierr = DMDAVecRestoreArray(user->fda, Rct, &rct); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}