/**
 * @file test_solver_kernels.c
 * @brief C unit tests for solver-side analytical and LES helper kernels.
 */

#include "test_support.h"

#include "AnalyticalSolutions.h"
#include "BodyForces.h"
#include "Filter.h"
#include "les.h"
#include "solvers.h"
/**
 * @brief Tests LES test-filter helper paths for representative cases.
 */

static PetscErrorCode TestLESTestFilterPaths(void)
{
    SimCtx simCtx;
    double values[3][3][3];
    double weights[3][3][3];

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    for (PetscInt k = 0; k < 3; ++k) {
        for (PetscInt j = 0; j < 3; ++j) {
            for (PetscInt i = 0; i < 3; ++i) {
                values[k][j][i] = 2.0;
                weights[k][j][i] = 1.0;
            }
        }
    }

    simCtx.testfilter_ik = 1;
    PetscCall(PicurvAssertRealNear(2.0, ApplyLESTestFilter(&simCtx, values, weights), 1.0e-12,
                                   "Simpson-rule filter should preserve a constant field"));

    simCtx.testfilter_ik = 0;
    PetscCall(PicurvAssertRealNear(2.0, ApplyLESTestFilter(&simCtx, values, weights), 1.0e-12,
                                   "box filter should preserve a constant field"));

    for (PetscInt k = 0; k < 3; ++k) {
        for (PetscInt j = 0; j < 3; ++j) {
            for (PetscInt i = 0; i < 3; ++i) {
                weights[k][j][i] = 0.0;
            }
        }
    }
    PetscCall(PicurvAssertRealNear(0.0, ApplyLESTestFilter(&simCtx, values, weights), 1.0e-12,
                                   "box filter should return zero when all weights are zero"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests analytical geometry selection for supported analytical solutions.
 */

static PetscErrorCode TestAnalyticalGeometrySelection(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscErrorCode ierr_grid = 0;
    PetscErrorCode ierr_non_square = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 8, 8));
    PetscCall(PicurvAssertBool(AnalyticalTypeRequiresCustomGeometry("TGV3D"),
                               "TGV3D should require custom geometry"));
    PetscCall(PicurvAssertBool((PetscBool)!AnalyticalTypeRequiresCustomGeometry("ZERO_FLOW"),
                               "ZERO_FLOW should not require custom geometry"));

    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "ZERO_FLOW", sizeof(simCtx->AnalyticalSolutionType)));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_grid = SetAnalyticalGridInfo(user);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_grid != 0),
                               "SetAnalyticalGridInfo should reject analytical types without custom geometry"));

    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "TGV3D", sizeof(simCtx->AnalyticalSolutionType)));
    simCtx->block_number = 2;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_non_square = SetAnalyticalGridInfo(user);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_non_square != 0),
                               "TGV3D multi-block setup should reject non-square block counts"));

    simCtx->block_number = 1;
    PetscCall(SetAnalyticalGridInfo(user));
    PetscCall(PicurvAssertRealNear(0.0, user->Min_X, 1.0e-12, "TGV3D single-block xmin"));
    PetscCall(PicurvAssertRealNear(2.0 * PETSC_PI, user->Max_X, 1.0e-12, "TGV3D single-block xmax"));

    simCtx->block_number = 4;
    user->_this = 3;
    PetscCall(SetAnalyticalGridInfo(user));
    PetscCall(PicurvAssertRealNear(PETSC_PI, user->Min_X, 1.0e-12, "TGV3D multi-block xmin should reflect the block column"));
    PetscCall(PicurvAssertRealNear(2.0 * PETSC_PI, user->Max_X, 1.0e-12, "TGV3D multi-block xmax should reflect the block column"));
    PetscCall(PicurvAssertRealNear(PETSC_PI, user->Min_Y, 1.0e-12, "TGV3D multi-block ymin should reflect the block row"));
    PetscCall(PicurvAssertRealNear(2.0 * PETSC_PI, user->Max_Y, 1.0e-12, "TGV3D multi-block ymax should reflect the block row"));
    PetscCall(PicurvAssertRealNear(2.0 * PETSC_PI, user->Max_Z, 1.0e-12, "TGV3D multi-block zmax should span the full domain"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests analytical scalar verification helper routines.
 */

static PetscErrorCode TestAnalyticalScalarVerificationHelpers(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec target = NULL;
    PetscReal value = 0.0;
    PetscReal ***target_arr = NULL;
    PetscReal *positions = NULL;
    PetscReal *psi = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));
    simCtx->verificationScalar.enabled = PETSC_TRUE;
    PetscCall(PetscStrncpy(simCtx->verificationScalar.mode,
                           "analytical",
                           sizeof(simCtx->verificationScalar.mode)));

    PetscCall(PetscStrncpy(simCtx->verificationScalar.profile,
                           "CONSTANT",
                           sizeof(simCtx->verificationScalar.profile)));
    simCtx->verificationScalar.value = 2.5;
    PetscCall(EvaluateAnalyticalScalarProfile(simCtx, 0.2, 0.3, 0.4, 0.0, &value));
    PetscCall(PicurvAssertRealNear(2.5, value, 1.0e-12,
                                   "constant scalar profile should evaluate to the configured value"));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    positions[0] = 0.1; positions[1] = 0.2; positions[2] = 0.3;
    positions[3] = 0.8; positions[4] = 0.6; positions[5] = 0.4;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(SetAnalyticalScalarFieldOnParticles(user, "Psi"));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(PicurvAssertRealNear(2.5, psi[0], 1.0e-12,
                                   "SetAnalyticalScalarFieldOnParticles should overwrite the first particle scalar"));
    PetscCall(PicurvAssertRealNear(2.5, psi[1], 1.0e-12,
                                   "SetAnalyticalScalarFieldOnParticles should overwrite the second particle scalar"));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));

    PetscCall(PetscStrncpy(simCtx->verificationScalar.profile,
                           "LINEAR_X",
                           sizeof(simCtx->verificationScalar.profile)));
    simCtx->verificationScalar.phi0 = 1.0;
    simCtx->verificationScalar.slope_x = 2.0;
    PetscCall(EvaluateAnalyticalScalarProfile(simCtx, 0.25, 0.0, 0.0, 0.0, &value));
    PetscCall(PicurvAssertRealNear(1.5, value, 1.0e-12,
                                   "linear-x scalar profile should evaluate phi0 + slope_x * x"));

    PetscCall(VecDuplicate(user->Psi, &target));
    PetscCall(SetAnalyticalScalarFieldAtCellCenters(user, target));
    PetscCall(DMDAVecGetArrayRead(user->da, target, &target_arr));
    PetscCall(PicurvAssertRealNear(1.25, target_arr[1][1][1], 1.0e-12,
                                   "cell-center scalar fill should use the physical x center coordinate"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, target, &target_arr));
    PetscCall(VecDestroy(&target));

    PetscCall(PetscStrncpy(simCtx->verificationScalar.profile,
                           "SIN_PRODUCT",
                           sizeof(simCtx->verificationScalar.profile)));
    simCtx->verificationScalar.amplitude = 3.0;
    simCtx->verificationScalar.kx = PETSC_PI;
    simCtx->verificationScalar.ky = PETSC_PI;
    simCtx->verificationScalar.kz = PETSC_PI;
    PetscCall(EvaluateAnalyticalScalarProfile(simCtx, 0.5, 0.5, 0.5, 0.0, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12,
                                   "sin-product scalar profile should peak at pi/2 in each coordinate"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests analytical solution engine ZERO_FLOW, UNIFORM_FLOW, and unknown-type dispatch.
 */

static PetscErrorCode TestAnalyticalSolutionEngineDispatch(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscErrorCode ierr_unknown = 0;
    Cmpnts ***ucat = NULL;
    Cmpnts ***ubcs = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->Ucat, 3.0));
    PetscCall(VecSet(user->P, 5.0));
    PetscCall(VecSet(user->Bcs.Ubcs, 7.0));

    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "ZERO_FLOW", sizeof(simCtx->AnalyticalSolutionType)));
    PetscCall(AnalyticalSolutionEngine(simCtx));
    PetscCall(PicurvAssertVecConstant(user->Ucat, 0.0, 1.0e-12, "ZERO_FLOW should zero the Eulerian velocity field"));
    PetscCall(PicurvAssertVecConstant(user->P, 0.0, 1.0e-12, "ZERO_FLOW should zero the pressure field"));
    PetscCall(PicurvAssertVecConstant(user->Bcs.Ubcs, 0.0, 1.0e-12, "ZERO_FLOW should zero boundary-condition velocity data"));

    /* UNIFORM_FLOW now works in curvilinear form (sets Ucont via metric dot products,
       derives Ucat via Contra2Cart).  The test grid needs identity metrics so the
       transformation is invertible and Ucat recovers the physical velocity exactly. */
    {
        Cmpnts ***l_csi, ***l_eta, ***l_zet;
        DMDALocalInfo linfo;
        PetscCall(DMDAGetLocalInfo(user->fda, &linfo));
        PetscCall(DMDAVecGetArray(user->fda, user->lCsi, &l_csi));
        PetscCall(DMDAVecGetArray(user->fda, user->lEta, &l_eta));
        PetscCall(DMDAVecGetArray(user->fda, user->lZet, &l_zet));
        for (PetscInt k = linfo.zs; k < linfo.zs + linfo.zm; k++)
            for (PetscInt j = linfo.ys; j < linfo.ys + linfo.ym; j++)
                for (PetscInt i = linfo.xs; i < linfo.xs + linfo.xm; i++) {
                    l_csi[k][j][i] = (Cmpnts){1.0, 0.0, 0.0};
                    l_eta[k][j][i] = (Cmpnts){0.0, 1.0, 0.0};
                    l_zet[k][j][i] = (Cmpnts){0.0, 0.0, 1.0};
                }
        PetscCall(DMDAVecRestoreArray(user->fda, user->lZet, &l_zet));
        PetscCall(DMDAVecRestoreArray(user->fda, user->lEta, &l_eta));
        PetscCall(DMDAVecRestoreArray(user->fda, user->lCsi, &l_csi));
    }

    simCtx->AnalyticalUniformVelocity.x = 1.25;
    simCtx->AnalyticalUniformVelocity.y = -0.5;
    simCtx->AnalyticalUniformVelocity.z = 0.75;
    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "UNIFORM_FLOW", sizeof(simCtx->AnalyticalSolutionType)));
    PetscCall(AnalyticalSolutionEngine(simCtx));
    PetscCall(PicurvAssertVecConstant(user->P, 0.0, 1.0e-12, "UNIFORM_FLOW should keep the pressure field zero"));

    Cmpnts ***ucont = NULL;
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(PicurvAssertRealNear(1.25, ucat[1][1][1].x, 1.0e-12, "UNIFORM_FLOW should impose the configured x velocity"));
    PetscCall(PicurvAssertRealNear(-0.5, ucat[1][1][1].y, 1.0e-12, "UNIFORM_FLOW should impose the configured y velocity"));
    PetscCall(PicurvAssertRealNear(0.75, ucat[1][1][1].z, 1.0e-12, "UNIFORM_FLOW should impose the configured z velocity"));
    PetscCall(PicurvAssertRealNear(1.25, ucont[1][1][1].x, 1.0e-12, "UNIFORM_FLOW should set contravariant x flux (identity metric)"));
    PetscCall(PicurvAssertRealNear(-0.5, ucont[1][1][1].y, 1.0e-12, "UNIFORM_FLOW should set contravariant y flux (identity metric)"));
    PetscCall(PicurvAssertRealNear(0.75, ucont[1][1][1].z, 1.0e-12, "UNIFORM_FLOW should set contravariant z flux (identity metric)"));
    PetscCall(PicurvAssertRealNear(1.25, ubcs[0][1][1].x, 1.0e-12, "UNIFORM_FLOW should populate boundary x velocity"));
    PetscCall(PicurvAssertRealNear(-0.5, ubcs[0][1][1].y, 1.0e-12, "UNIFORM_FLOW should populate boundary y velocity"));
    PetscCall(PicurvAssertRealNear(0.75, ubcs[0][1][1].z, 1.0e-12, "UNIFORM_FLOW should populate boundary z velocity"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));

    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "NOT_A_REAL_ANALYTICAL_TYPE", sizeof(simCtx->AnalyticalSolutionType)));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_unknown = AnalyticalSolutionEngine(simCtx);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_unknown != 0),
                               "AnalyticalSolutionEngine should reject unknown analytical type strings"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests exact Taylor-Green samples on selected Eulerian interior and boundary points.
 */
static PetscErrorCode TestAnalyticalSolutionEngineTaylorGreenSamples(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***cent = NULL;
    Cmpnts ***cent_x = NULL;
    Cmpnts ***cent_y = NULL;
    Cmpnts ***cent_z = NULL;
    Cmpnts ***ucat = NULL;
    Cmpnts ***ubcs = NULL;
    PetscReal ***p = NULL;
    const PetscReal vel_decay = PetscExpReal(-0.5);

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->ren = 2.0;
    simCtx->ti = 0.5;
    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "TGV3D", sizeof(simCtx->AnalyticalSolutionType)));

    PetscCall(DMDAVecGetArray(user->fda, user->Cent, &cent));
    PetscCall(DMDAVecGetArray(user->fda, user->Centx, &cent_x));
    PetscCall(DMDAVecGetArray(user->fda, user->Centy, &cent_y));
    PetscCall(DMDAVecGetArray(user->fda, user->Centz, &cent_z));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                cent[k][j][i].x = 0.0;
                cent[k][j][i].y = 0.0;
                cent[k][j][i].z = 0.0;
                cent_x[k][j][i] = cent[k][j][i];
                cent_y[k][j][i] = cent[k][j][i];
                cent_z[k][j][i] = cent[k][j][i];
            }
        }
    }
    cent[1][1][1].x = 0.5 * PETSC_PI;
    cent[1][1][1].y = 0.0;
    cent[1][1][1].z = 0.0;
    cent[1][2][1].x = 0.0;
    cent[1][2][1].y = 0.5 * PETSC_PI;
    cent[1][2][1].z = 0.0;
    cent_z[0][1][1].x = 0.5 * PETSC_PI;
    cent_z[0][1][1].y = 0.0;
    cent_z[0][1][1].z = 0.0;
    PetscCall(DMDAVecRestoreArray(user->fda, user->Centz, &cent_z));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Centy, &cent_y));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Centx, &cent_x));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Cent, &cent));

    PetscCall(AnalyticalSolutionEngine(simCtx));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecGetArrayRead(user->da, user->P, &p));
    PetscCall(PicurvAssertRealNear(vel_decay, ucat[1][1][1].x, 1.0e-12, "TGV sample should set the expected interior x velocity"));
    PetscCall(PicurvAssertRealNear(0.0, ucat[1][1][1].y, 1.0e-12, "TGV sample should keep the paired interior y velocity at zero"));
    PetscCall(PicurvAssertRealNear(-vel_decay, ucat[1][2][1].y, 1.0e-12, "TGV sample should set the expected interior y velocity"));
    PetscCall(PicurvAssertRealNear(0.0, p[1][1][1], 1.0e-12, "Chosen TGV sample should produce zero pressure"));
    PetscCall(PicurvAssertRealNear(vel_decay, ubcs[0][1][1].x, 1.0e-12, "TGV sample should set the expected boundary x velocity"));
    PetscCall(PicurvAssertRealNear(0.0, ubcs[0][1][1].y, 1.0e-12, "TGV boundary sample should keep the paired y velocity at zero"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle analytical-solution dispatch for TGV3D, UNIFORM_FLOW, and non-analytical no-op paths.
 */

static PetscErrorCode TestAnalyticalSolutionForParticlesDispatch(void)
{
    SimCtx simCtx;
    Vec tempVec = NULL;
    PetscReal *data = NULL;
    const PetscReal vel_decay = PetscExpReal(-0.5);

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    simCtx.ren = 2.0;
    simCtx.ti = 0.5;
    PetscCall(PetscStrncpy(simCtx.AnalyticalSolutionType, "TGV3D", sizeof(simCtx.AnalyticalSolutionType)));

    PetscCall(VecCreateSeq(PETSC_COMM_SELF, 6, &tempVec));
    PetscCall(VecGetArray(tempVec, &data));
    data[0] = 0.5 * PETSC_PI; data[1] = 0.0;          data[2] = 0.0;
    data[3] = 0.0;           data[4] = 0.5 * PETSC_PI; data[5] = 0.0;
    PetscCall(VecRestoreArray(tempVec, &data));

    PetscCall(SetAnalyticalSolutionForParticles(tempVec, &simCtx));
    PetscCall(VecGetArray(tempVec, &data));
    PetscCall(PicurvAssertRealNear(vel_decay, data[0], 1.0e-12,
                                   "TGV3D particle dispatch should populate the x velocity at x=pi/2"));
    PetscCall(PicurvAssertRealNear(0.0, data[1], 1.0e-12,
                                   "TGV3D particle dispatch should leave the first particle y velocity at zero"));
    PetscCall(PicurvAssertRealNear(0.0, data[2], 1.0e-12,
                                   "TGV3D particle dispatch should leave the first particle z velocity at zero"));
    PetscCall(PicurvAssertRealNear(0.0, data[3], 1.0e-12,
                                   "TGV3D particle dispatch should leave the second particle x velocity at zero"));
    PetscCall(PicurvAssertRealNear(-vel_decay, data[4], 1.0e-12,
                                   "TGV3D particle dispatch should populate the y velocity at y=pi/2"));
    PetscCall(PicurvAssertRealNear(0.0, data[5], 1.0e-12,
                                   "TGV3D particle dispatch should leave the second particle z velocity at zero"));
    PetscCall(VecRestoreArray(tempVec, &data));

    PetscCall(PetscStrncpy(simCtx.AnalyticalSolutionType, "ZERO_FLOW", sizeof(simCtx.AnalyticalSolutionType)));
    PetscCall(VecGetArray(tempVec, &data));
    data[0] = 3.0;
    data[1] = 4.0;
    data[2] = 5.0;
    PetscCall(VecRestoreArray(tempVec, &data));
    PetscCall(SetAnalyticalSolutionForParticles(tempVec, &simCtx));
    PetscCall(VecGetArray(tempVec, &data));
    PetscCall(PicurvAssertRealNear(3.0, data[0], 1.0e-12,
                                   "Non-TGV particle dispatch should leave the vector untouched"));
    PetscCall(PicurvAssertRealNear(4.0, data[1], 1.0e-12,
                                   "Non-TGV particle dispatch should preserve the y component"));
    PetscCall(PicurvAssertRealNear(5.0, data[2], 1.0e-12,
                                   "Non-TGV particle dispatch should preserve the z component"));
    PetscCall(VecRestoreArray(tempVec, &data));

    simCtx.AnalyticalUniformVelocity.x = 0.125;
    simCtx.AnalyticalUniformVelocity.y = -0.25;
    simCtx.AnalyticalUniformVelocity.z = 0.375;
    PetscCall(PetscStrncpy(simCtx.AnalyticalSolutionType, "UNIFORM_FLOW", sizeof(simCtx.AnalyticalSolutionType)));
    PetscCall(SetAnalyticalSolutionForParticles(tempVec, &simCtx));
    PetscCall(VecGetArray(tempVec, &data));
    PetscCall(PicurvAssertRealNear(0.125, data[0], 1.0e-12,
                                   "UNIFORM_FLOW particle dispatch should populate the x velocity"));
    PetscCall(PicurvAssertRealNear(-0.25, data[1], 1.0e-12,
                                   "UNIFORM_FLOW particle dispatch should populate the y velocity"));
    PetscCall(PicurvAssertRealNear(0.375, data[2], 1.0e-12,
                                   "UNIFORM_FLOW particle dispatch should populate the z velocity"));
    PetscCall(PicurvAssertRealNear(0.125, data[3], 1.0e-12,
                                   "UNIFORM_FLOW particle dispatch should use the same x velocity for each particle"));
    PetscCall(PicurvAssertRealNear(-0.25, data[4], 1.0e-12,
                                   "UNIFORM_FLOW particle dispatch should use the same y velocity for each particle"));
    PetscCall(PicurvAssertRealNear(0.375, data[5], 1.0e-12,
                                   "UNIFORM_FLOW particle dispatch should use the same z velocity for each particle"));
    PetscCall(VecRestoreArray(tempVec, &data));

    PetscCall(VecDestroy(&tempVec));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests deterministic LES eddy-viscosity computation on a linear velocity field.
 */

static PetscErrorCode TestComputeEddyViscosityLESDeterministicField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucat = NULL;
    PetscReal ***nu_t = NULL;
    const PetscReal expected_nu_t = 0.25 * PetscSqrtReal(2.0);

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 5, 5, 5));
    PetscCall(DMCreateGlobalVector(user->da, &user->Nu_t));
    PetscCall(DMCreateLocalVector(user->da, &user->lNu_t));
    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    PetscCall(DMCreateLocalVector(user->da, &user->lCs));

    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(VecSet(user->Nu_t, 0.0));
    PetscCall(VecSet(user->CS, 0.5));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                ucat[k][j][i].x = (PetscReal)i;
                ucat[k][j][i].y = 0.0;
                ucat[k][j][i].z = 0.0;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalBegin(user->da, user->CS, INSERT_VALUES, user->lCs));
    PetscCall(DMGlobalToLocalEnd(user->da, user->CS, INSERT_VALUES, user->lCs));

    PetscCall(ComputeEddyViscosityLES(user));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Nu_t, &nu_t));
    PetscCall(PicurvAssertRealNear(expected_nu_t, nu_t[2][2][2], 1.0e-6,
                                   "linear velocity field should yield deterministic LES eddy viscosity"));
    PetscCall(PicurvAssertRealNear(0.0, nu_t[0][2][2], 1.0e-12,
                                   "boundary cells should remain untouched by the interior LES loop"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Nu_t, &nu_t));

    PetscCall(VecDestroy(&user->CS));
    PetscCall(VecDestroy(&user->lCs));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests FlowSolver guardrails for unsupported momentum solver selections.
 */

static PetscErrorCode TestFlowSolverRejectsUnsupportedMomentumSolverType(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscErrorCode ierr_flow = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->mom_solver_type = (MomentumSolverType)999;

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_flow = FlowSolver(simCtx);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_flow != 0),
                               "FlowSolver should reject unsupported momentum solver selectors"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests driven-channel flow source-term evaluation.
 */

static PetscErrorCode TestDrivenChannelFlowSource(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec rct = NULL;
    Cmpnts ***rct_arr = NULL;
    Cmpnts ***l_csi = NULL;
    Cmpnts ***l_eta = NULL;
    Cmpnts ***l_zet = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecDuplicate(user->Ucont, &rct));
    PetscCall(VecZeroEntries(rct));

    user->boundary_faces[BC_FACE_NEG_X].face_id = BC_FACE_NEG_X;
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    user->boundary_faces[BC_FACE_POS_X].face_id = BC_FACE_POS_X;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    simCtx->bulkVelocityCorrection = 2.0;
    simCtx->dt = 1.0;
    simCtx->forceScalingFactor = 1.0;
    simCtx->drivingForceMagnitude = 0.0;
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMDAVecGetArray(user->fda, user->lCsi, &l_csi));
    PetscCall(DMDAVecGetArray(user->fda, user->lEta, &l_eta));
    PetscCall(DMDAVecGetArray(user->fda, user->lZet, &l_zet));
    l_csi[1][1][1].x = 1.0; l_csi[1][1][1].y = 0.0; l_csi[1][1][1].z = 0.0;
    l_eta[1][1][1].x = 0.0; l_eta[1][1][1].y = 1.0; l_eta[1][1][1].z = 0.0;
    l_zet[1][1][1].x = 0.0; l_zet[1][1][1].y = 0.0; l_zet[1][1][1].z = 1.0;
    PetscCall(DMDAVecRestoreArray(user->fda, user->lZet, &l_zet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lEta, &l_eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lCsi, &l_csi));

    PetscCall(ComputeDrivenChannelFlowSource(user, rct));
    PetscCall(DMDAVecGetArrayRead(user->fda, rct, &rct_arr));
    PetscCall(PicurvAssertRealNear(1.5, simCtx->drivingForceMagnitude, 1.0e-12,
                                   "driven flow source should update the controller magnitude"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].y, 1.0e-12,
                                   "driven flow source should leave the y component unchanged"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].z, 1.0e-12,
                                   "driven flow source should leave the z component unchanged"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, rct, &rct_arr));

    PetscCall(VecDestroy(&rct));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-solver PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"les-filter-paths", TestLESTestFilterPaths},
        {"analytical-geometry-selection", TestAnalyticalGeometrySelection},
        {"analytical-scalar-verification-helpers", TestAnalyticalScalarVerificationHelpers},
        {"analytical-solution-engine-dispatch", TestAnalyticalSolutionEngineDispatch},
        {"analytical-solution-engine-taylor-green-samples", TestAnalyticalSolutionEngineTaylorGreenSamples},
        {"analytical-solution-for-particles-dispatch", TestAnalyticalSolutionForParticlesDispatch},
        {"compute-eddy-viscosity-les-deterministic-field", TestComputeEddyViscosityLESDeterministicField},
        {"flow-solver-rejects-unsupported-momentum-solver-type", TestFlowSolverRejectsUnsupportedMomentumSolverType},
        {"driven-channel-flow-source", TestDrivenChannelFlowSource},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv solver utility tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-solver", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
