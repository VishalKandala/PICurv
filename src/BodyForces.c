#include "BodyForces.h"
#include "io.h"


//////////////////////////////////////////////
// DRIVEN CHANNEL FLOW FORCE(EQUIVALENT) TERM 
/////////////////////////////////////////////

/**
 * @brief The axis a driven periodic handler drives, or ' ' when none is configured.
 * @param[in] user Block whose boundary faces are inspected.
 * @return 'X', 'Y', 'Z', or ' '.
 */
static char DrivenFlowDirection(const UserCtx *user)
{
    for (int i = 0; i < 6; i++) {
        const BCHandlerType handler_type = user->boundary_faces[i].handler_type;
        if (handler_type == BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX ||
            handler_type == BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX) {
            switch (user->boundary_faces[i].face_id) {
                case BC_FACE_NEG_X: case BC_FACE_POS_X: return 'X';
                case BC_FACE_NEG_Y: case BC_FACE_POS_Y: return 'Y';
                case BC_FACE_NEG_Z: case BC_FACE_POS_Z: return 'Z';
            }
            return ' ';
        }
    }
    return ' ';
}

#undef __FUNCT__
#define __FUNCT__ "ComputeDrivenChannelFlowSource"
/**
 * @brief Internal helper implementation: `ComputeDrivenChannelFlowSource()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeDrivenChannelFlowSource(UserCtx *user, Vec Rct)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    PetscFunctionBeginUser;

    // --- Step 1: Discover if and where a driven flow is active ---
    const char drivenDirection = DrivenFlowDirection(user); // ' ' when none is configured

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
    //
    // ONCE PER PHYSICAL STEP, NOT ONCE PER CALL. This function runs from
    // ComputeRHS, which executes once per Jameson RK stage under the Picard
    // solver and once per residual evaluation under Newton-Krylov. The smoothing
    // below carries state in simCtx across calls, so advancing it every call
    // would walk the applied force toward its target within a single timestep
    // (0.5, then 0.75, then 0.875 ... of the way there). The force would then
    // depend on how many residual evaluations preceded it - history dependence
    // that MomentumNewtonKrylov_FormResidual() explicitly forbids, and that
    // breaks the constant-forcing assumption behind the Picard shadow-Jacobian
    // estimate. Resolve it once for the step and reuse it thereafter.
    const PetscReal forceScalingFactor  = simCtx->forceScalingFactor;
    PetscReal drivingForceMagnitude;

    if (simCtx->drivingForceStep != simCtx->step) {
        const PetscReal targetForce = (bulkVelocityCorrection / simCtx->dt / 1.0 * COEF_TIME_ACCURACY); // replaced simCtx->st with 1.0.
        drivingForceMagnitude = (simCtx->drivingForceMagnitude * 0.5) + (targetForce * 0.5);
        simCtx->drivingForceMagnitude = drivingForceMagnitude;
        simCtx->drivingForceStep = simCtx->step;
    } else {
        drivingForceMagnitude = simCtx->drivingForceMagnitude;
    }
    
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
                                LOG_ALLOW(LOCAL, LOG_DEBUG,"Body Force %le added at (%d,%d,%d)\n",momentumSource, k, j, i);
                                hasLoggedApplication = PETSC_TRUE;
                            }
                                break;
                        case 'Y':
                            faceArea = sqrt(eta[k][j][i].x * eta[k][j][i].x + eta[k][j][i].y * eta[k][j][i].y + eta[k][j][i].z * eta[k][j][i].z);
                            momentumSource = drivingForceMagnitude * forceScalingFactor * faceArea;
                            rct[k][j][i].y += momentumSource;
                            
                            // Log details for the very first point where force is applied on this rank.
                            if (!hasLoggedApplication) {
                                LOG_ALLOW(LOCAL, LOG_DEBUG,"Body Force %le added at (%d,%d,%d)\n",momentumSource, k, j, i);
                                hasLoggedApplication = PETSC_TRUE;
                            }
                            break;
                        case 'Z':
                            faceArea = sqrt(zet[k][j][i].x * zet[k][j][i].x + zet[k][j][i].y * zet[k][j][i].y + zet[k][j][i].z * zet[k][j][i].z);
                            momentumSource = drivingForceMagnitude * forceScalingFactor * faceArea;
                            rct[k][j][i].z += momentumSource;

                            // Log details for the very first point where force is applied on this rank.
                            if (!hasLoggedApplication) {
                                LOG_ALLOW(LOCAL, LOG_DEBUG,"Body Force %le added at (%d,%d,%d)\n",momentumSource, k, j, i);
                                hasLoggedApplication = PETSC_TRUE;
                            }
                            break;
                    }
                }
            }
        }
    }

    // --- Step 6: Restore arrays ---
    ierr = DMDAVecRestoreArray(user->fda, Rct, &rct); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LogDrivenFlowDiagnostics"
/**
 * @brief Implementation of \ref LogDrivenFlowDiagnostics().
 * @details Full API contract is documented with the header declaration in
 *          `include/BodyForces.h`.
 */
PetscErrorCode LogDrivenFlowDiagnostics(UserCtx *user)
{
    SimCtx    *simCtx = user->simCtx;
    const char direction = DrivenFlowDirection(user);
    FILE      *file = NULL;

    PetscFunctionBeginUser;

    /* The controller state is global: every rank holds the same values, set from the
       controller's collective reductions, so rank 0 reports it once per step. */
    if (direction == ' ' || user->_this != 0 || simCtx->rank != 0) PetscFunctionReturn(0);

    /* What the momentum equation actually received this step. The source is resolved
       once per step and skipped entirely when the correction is negligible, in which
       case the smoothed magnitude is stale and nothing was applied. */
    const PetscBool applied = (PetscBool)(simCtx->drivingForceStep == simCtx->step &&
                                          PetscAbsReal(simCtx->bulkVelocityCorrection) >= 1.0e-12);
    const PetscReal acceleration = applied ? simCtx->drivingForceMagnitude * simCtx->forceScalingFactor : 0.0;
    const PetscReal area = simCtx->drivenFluxArea;
    const PetscReal bulk_velocity = (area > 0.0) ? simCtx->drivenFluxMeasured / area : 0.0;
    PetscReal physical_time = 0.0;

    PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
    PetscCall(PicurvOpenDiagnosticsCsv(simCtx, "driven_flow.csv",
                                       "step,time,direction,target_flux,measured_flux,cross_section_area,"
                                       "bulk_velocity,bulk_velocity_correction,driving_acceleration,physical_time",
                                       &file));
    fprintf(file, "%d,%.6e,%c,%.10e,%.10e,%.10e,%.10e,%.6e,%.10e,%.6e\n",
            (int)simCtx->step, (double)simCtx->ti, direction,
            (double)simCtx->targetVolumetricFlux, (double)simCtx->drivenFluxMeasured, (double)area,
            (double)bulk_velocity, (double)simCtx->bulkVelocityCorrection, (double)acceleration,
            (double)physical_time);
    PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
               "Unable to close the driven-flow diagnostics file.");
    PetscFunctionReturn(0);
}
