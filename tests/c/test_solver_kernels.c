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
#include "momentumsolvers.h"
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

    simCtx->block_number = 1;
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
    PetscCall(DMDAVecGetArray(user->fda, user->lCentx, &cent_x));
    PetscCall(DMDAVecGetArray(user->fda, user->lCenty, &cent_y));
    PetscCall(DMDAVecGetArray(user->fda, user->lCentz, &cent_z));
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
    PetscCall(DMDAVecRestoreArray(user->fda, user->lCentz, &cent_z));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lCenty, &cent_y));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lCentx, &cent_x));
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
/* ===================================================================== *
 *   Stage A3: momentum pseudo-time stability estimate (shadow) tests     *
 *   Unit Cartesian grid via identity metrics: Csi=(1,0,0), aj=1, so the  *
 *   analytic targets are exact (lambda_t=a0/dt, viscous=4*6*nu=24nu,     *
 *   centered conv reduces to sum|U_f|, QUICK interior=4/3, boundary=2.5).*
 * ===================================================================== */

/**
 * @brief Sets the mathematical BC type on all six faces of a test UserCtx.
 * @param user Test user context.
 * @param t    BCType to assign to every face (e.g. PERIODIC, INLET).
 */
static void MomSetAllBC(UserCtx *user, BCType t)
{
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = t;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = t;
    user->boundary_faces[BC_FACE_NEG_Y].mathematical_type = t;
    user->boundary_faces[BC_FACE_POS_Y].mathematical_type = t;
    user->boundary_faces[BC_FACE_NEG_Z].mathematical_type = t;
    user->boundary_faces[BC_FACE_POS_Z].mathematical_type = t;
}

/**
 * @brief Fills a local Cmpnts vector (ghosts included) with a uniform vector value.
 * @param fda  Vector DM owning lvec.
 * @param lvec Local Cmpnts vector to fill.
 * @param x,y,z Uniform component values.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomFillLocalCmpnts(DM fda, Vec lvec, PetscReal x, PetscReal y, PetscReal z)
{
    Cmpnts ***a;
    PetscInt gxs, gys, gzs, gxm, gym, gzm;
    PetscFunctionBeginUser;
    PetscCall(DMDAGetGhostCorners(fda, &gxs, &gys, &gzs, &gxm, &gym, &gzm));
    PetscCall(DMDAVecGetArray(fda, lvec, &a));
    for (PetscInt k = gzs; k < gzs+gzm; ++k)
        for (PetscInt j = gys; j < gys+gym; ++j)
            for (PetscInt i = gxs; i < gxs+gxm; ++i) { a[k][j][i].x = x; a[k][j][i].y = y; a[k][j][i].z = z; }
    PetscCall(DMDAVecRestoreArray(fda, lvec, &a));
    PetscFunctionReturn(0);
}

/**
 * @brief Fills a local scalar vector (ghosts included) with a uniform value.
 * @param da   Scalar DM owning lvec.
 * @param lvec Local scalar vector to fill.
 * @param v    Uniform value.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomFillLocalScalar(DM da, Vec lvec, PetscReal v)
{
    PetscReal ***a;
    PetscInt gxs, gys, gzs, gxm, gym, gzm;
    PetscFunctionBeginUser;
    PetscCall(DMDAGetGhostCorners(da, &gxs, &gys, &gzs, &gxm, &gym, &gzm));
    PetscCall(DMDAVecGetArray(da, lvec, &a));
    for (PetscInt k = gzs; k < gzs+gzm; ++k)
        for (PetscInt j = gys; j < gys+gym; ++j)
            for (PetscInt i = gxs; i < gxs+gxm; ++i) a[k][j][i] = v;
    PetscCall(DMDAVecRestoreArray(da, lvec, &a));
    PetscFunctionReturn(0);
}

/**
 * @brief Sets a single local-scalar cell value (e.g. to mark one nvert solid).
 * @param da   Scalar DM owning lvec.
 * @param lvec Local scalar vector.
 * @param ci,cj,ck Cell indices to set.
 * @param v    Value to assign.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomSetLocalScalarCell(DM da, Vec lvec, PetscInt ci, PetscInt cj, PetscInt ck, PetscReal v)
{
    PetscReal ***a;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(da, lvec, &a));
    a[ck][cj][ci] = v;
    PetscCall(DMDAVecRestoreArray(da, lvec, &a));
    PetscFunctionReturn(0);
}

/* 1. Shared BDF coefficient: BDF1 on cold start (ti==1) and restart (ti==tistart), BDF2 otherwise. */
static PetscErrorCode TestMomentumBDFCoefficient(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL;
    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    simCtx->StartStep = 0; simCtx->step = 1;
    PetscCall(PicurvAssertRealNear(1.0, MomentumBDFCoefficient(simCtx), 1e-12, "a0 BDF1 on cold-start step 1"));
    PetscCall(PicurvAssertBool((PetscBool)(!MomentumUsesBDF2(simCtx)), "step 1 is BDF1"));
    simCtx->step = 5;
    PetscCall(PicurvAssertRealNear(1.5, MomentumBDFCoefficient(simCtx), 1e-12, "a0 BDF2 on interior step"));
    PetscCall(PicurvAssertBool(MomentumUsesBDF2(simCtx), "step 5 (StartStep 0) is BDF2"));
    simCtx->StartStep = 5; simCtx->step = 5;
    PetscCall(PicurvAssertRealNear(1.0, MomentumBDFCoefficient(simCtx), 1e-12, "a0 BDF1 on restart step (ti==tistart)"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* Common setup: unit-Cartesian fixture with identity metrics; caller tweaks flags/fields. */
static PetscErrorCode MomMakeUnitGrid(SimCtx **simCtx, UserCtx **user, PetscInt n, BCType bc)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(simCtx, user, n, n, n));
    MomSetAllBC(*user, bc);
    (*simCtx)->dt = 0.1; (*simCtx)->step = 1; (*simCtx)->StartStep = 0;   /* a0=1 -> lambda_t=10 */
    (*simCtx)->ren = 1.0; (*simCtx)->les = 0; (*simCtx)->rans = 0;
    (*simCtx)->central = 0; (*simCtx)->invicid = 0; (*simCtx)->block_number = 1;
    /* The minimal fixture does not allocate lNu_t (LES off by default); create a
       zeroed one so the estimator can read it when a test enables LES/RANS. */
    if (!(*user)->lNu_t) PetscCall(DMCreateLocalVector((*user)->da, &(*user)->lNu_t));
    PetscCall(MomFillLocalCmpnts((*user)->fda, (*user)->lUcont, 0.0, 0.0, 0.0));
    PetscCall(MomFillLocalCmpnts((*user)->fda, (*user)->lUcat,  0.0, 0.0, 0.0));
    PetscCall(MomFillLocalScalar((*user)->da, (*user)->lNvert, 0.0));
    PetscCall(MomFillLocalScalar((*user)->da, (*user)->lNu_t,  0.0));
    PetscFunctionReturn(0);
}

/* 2. Centered convection: f_c=1, and Aj*U^xi = u/dx scaling (identity grid -> lambda_c = sum|u|). */
static PetscErrorCode TestMomentumStabilityCenteredConvection(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->central = 1; simCtx->invicid = 1;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 2.0, 0.0, 0.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(10.0, rep.lambda_t, 1e-9, "lambda_t = a0/dt"));
    PetscCall(PicurvAssertRealNear(2.0,  rep.lambda_c, 1e-9, "centered conv f_c=1 (Aj*U^xi scaling)"));
    PetscCall(PicurvAssertRealNear(0.0,  rep.lambda_v, 1e-9, "inviscid: no viscous term"));
    PetscCall(PicurvAssertRealNear(12.0, rep.lambda,   1e-9, "total lambda"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)rep.cclass, "interior cell class"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 3. Interior QUICK convective factor = 4/3. */
static PetscErrorCode TestMomentumStabilityQuickInterior(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->central = 0; simCtx->invicid = 1;       /* QUICK branch */
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 3.0, 0.0, 0.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(4.0, rep.lambda_c, 1e-9, "QUICK interior f_c=4/3 (4/3*3=4)"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)rep.cclass, "interior cell class"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 4. Boundary/IB-modified QUICK convective factor = 2.5 (conservative branch class). */
static PetscErrorCode TestMomentumStabilityQuickBoundary(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, INLET));   /* non-periodic -> boundary band active */
    simCtx->central = 0; simCtx->invicid = 1;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 8.0, 0.0, 0.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(20.0, rep.lambda_c, 1e-9, "boundary QUICK f_c=2.5 (2.5*8=20)"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)rep.cclass, "physical-boundary cell class"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 5. Viscous longitudinal full-stress factor: 4*6*nu = 24nu (= 8nu*(1/dx^2+1/dy^2+1/dz^2)). */
static PetscErrorCode TestMomentumStabilityViscousLongitudinal(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 1.0;          /* nu = 1, no convection */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(24.0, rep.lambda_v, 1e-9, "viscous = 4*6*nu = 24nu (factor-8 longitudinal)"));
    PetscCall(PicurvAssertRealNear(34.0, rep.lambda,   1e-9, "total lambda_t+lambda_v"));
    PetscCall(PicurvAssertIntEqual((PetscInt)MOM_STAB_LIMITER_VISCOSITY, (PetscInt)rep.limiter, "viscosity-limited"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 6. Viscous estimate scales linearly with molecular viscosity (1/Re). */
static PetscErrorCode TestMomentumStabilityViscousScalesWithNu(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 4.0;          /* nu = 0.25 -> 24*0.25 = 6 */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(6.0, rep.lambda_v, 1e-9, "viscous scales with nu=1/Re"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 7. LES eddy viscosity adds into nu_eff,f (face average): nu_eff = 1/Re + nu_t. */
static PetscErrorCode TestMomentumStabilityLESEddyViscosity(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 1.0; simCtx->les = 1;   /* nu_eff = 1 + nu_t */
    PetscCall(MomFillLocalScalar(user->da, user->lNu_t, 1.0)); /* nu_t = 1 -> nu_eff = 2 -> 24*2 = 48 */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(48.0, rep.lambda_v, 1e-9, "LES nu_t adds: 4*6*(1+1) = 48"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 8. One-sided viscous branch near a solid: conservative x2 applied once per cell. */
static PetscErrorCode TestMomentumStabilityOneSidedViscous(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 1.0;
    /* Mark one interior cell solid; its fluid cross-neighbours get the one-sided x2 (24 -> 48). */
    PetscCall(MomSetLocalScalarCell(user->da, user->lNvert, 4, 5, 4, 1.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(48.0, rep.lambda_v, 1e-9, "one-sided viscous x2 (24 -> 48)"));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)rep.one_sided, "one-sided multiplier flagged"));
    PetscCall(PicurvAssertIntEqual(2, (PetscInt)rep.cclass, "IB-adjacent cell class"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* 9. Read-only: the estimator must not mutate any solver field. */
static PetscErrorCode TestMomentumStabilityReadOnly(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscReal n_ucont_0, n_ucat_0, n_nvert_0, n_ucont_1, n_ucat_1, n_nvert_1;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 2.0;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 1.5, -0.7, 0.3));
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcat,  0.9,  0.4, -0.2));
    PetscCall(VecNorm(user->lUcont, NORM_2, &n_ucont_0));
    PetscCall(VecNorm(user->lUcat,  NORM_2, &n_ucat_0));
    PetscCall(VecNorm(user->lNvert, NORM_2, &n_nvert_0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_D, &rep));
    PetscCall(VecNorm(user->lUcont, NORM_2, &n_ucont_1));
    PetscCall(VecNorm(user->lUcat,  NORM_2, &n_ucat_1));
    PetscCall(VecNorm(user->lNvert, NORM_2, &n_nvert_1));
    PetscCall(PicurvAssertRealNear(n_ucont_0, n_ucont_1, 1e-14, "lUcont unchanged by estimator"));
    PetscCall(PicurvAssertRealNear(n_ucat_0,  n_ucat_1,  1e-14, "lUcat unchanged by estimator"));
    PetscCall(PicurvAssertRealNear(n_nvert_0, n_nvert_1, 1e-14, "lNvert unchanged by estimator"));
    PetscCall(PicurvAssertBool((PetscBool)(rep.lambda > 0.0), "estimate finite and positive"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Sets uniform Cartesian metrics for spacing (dx,dy,dz) on all cell and face arrays.
 * @param user Test user context.
 * @param dx,dy,dz Cartesian cell spacings.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomSetCartesianMetrics(UserCtx *user, PetscReal dx, PetscReal dy, PetscReal dz)
{
    const PetscReal sx = dy*dz, sy = dx*dz, sz = dx*dy, aj = 1.0/(dx*dy*dz);
    Vec csiv[] = {user->lCsi, user->lICsi, user->lJCsi, user->lKCsi};
    Vec etav[] = {user->lEta, user->lIEta, user->lJEta, user->lKEta};
    Vec zetv[] = {user->lZet, user->lIZet, user->lJZet, user->lKZet};
    Vec ajv[]  = {user->lAj,  user->lIAj,  user->lJAj,  user->lKAj};
    PetscFunctionBeginUser;
    for (int t = 0; t < 4; ++t) {
        PetscCall(MomFillLocalCmpnts(user->fda, csiv[t], sx, 0, 0));
        PetscCall(MomFillLocalCmpnts(user->fda, etav[t], 0, sy, 0));
        PetscCall(MomFillLocalCmpnts(user->fda, zetv[t], 0, 0, sz));
        PetscCall(MomFillLocalScalar(user->da, ajv[t], aj));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Fills lUcont.x as a linear ramp slope*i (nonzero discrete contravariant divergence).
 * @param fda Vector DM.
 * @param lvec Local Cmpnts vector.
 * @param slope Ramp slope.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomFillUcontXRamp(DM fda, Vec lvec, PetscReal slope)
{
    Cmpnts ***a;
    PetscInt gxs, gys, gzs, gxm, gym, gzm;
    PetscFunctionBeginUser;
    PetscCall(DMDAGetGhostCorners(fda, &gxs, &gys, &gzs, &gxm, &gym, &gzm));
    PetscCall(DMDAVecGetArray(fda, lvec, &a));
    for (PetscInt k = gzs; k < gzs+gzm; ++k)
        for (PetscInt j = gys; j < gys+gym; ++j)
            for (PetscInt i = gxs; i < gxs+gxm; ++i) { a[k][j][i].x = slope*i; a[k][j][i].y = 0; a[k][j][i].z = 0; }
    PetscCall(DMDAVecRestoreArray(fda, lvec, &a));
    PetscFunctionReturn(0);
}

/**
 * @brief Fills lUcat as a linear shear u=(gamma*j,0,0) so |grad u|_inf = gamma.
 * @param fda Vector DM.
 * @param lvec Local Cmpnts vector.
 * @param gamma Shear rate.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomFillUcatShearY(DM fda, Vec lvec, PetscReal gamma)
{
    Cmpnts ***a;
    PetscInt gxs, gys, gzs, gxm, gym, gzm;
    PetscFunctionBeginUser;
    PetscCall(DMDAGetGhostCorners(fda, &gxs, &gys, &gzs, &gxm, &gym, &gzm));
    PetscCall(DMDAVecGetArray(fda, lvec, &a));
    for (PetscInt k = gzs; k < gzs+gzm; ++k)
        for (PetscInt j = gys; j < gys+gym; ++j)
            for (PetscInt i = gxs; i < gxs+gxm; ++i) { a[k][j][i].x = gamma*j; a[k][j][i].y = 0; a[k][j][i].z = 0; }
    PetscCall(DMDAVecRestoreArray(fda, lvec, &a));
    PetscFunctionReturn(0);
}

/* A3.5-1. Anisotropic Cartesian viscous: detects metric-normal argument-order mistakes
 *         that the identity-metric test cannot. Target = 8nu*(1/dx^2+1/dy^2+1/dz^2). */
static PetscErrorCode TestMomentumStabilityAnisotropicViscous(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 1.0;
    PetscCall(MomSetCartesianMetrics(user, 1.0, 2.0, 4.0));   /* 8*(1 + 1/4 + 1/16) = 10.5 */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(10.5, rep.lambda_v, 1e-9, "anisotropic viscous = 8nu*sum(1/d^2)"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-2. Directional QUICK: only the boundary-modified direction gets 2.5; others stay 4/3. */
static PetscErrorCode TestMomentumStabilityDirectionalQuick(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->central = 0; simCtx->invicid = 1;
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = INLET;   /* x non-periodic only */
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = INLET;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 3.0, 3.0, 3.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    /* boundary cell: 0.5*(2.5*6 + (4/3)*6 + (4/3)*6) = 15.5 ; NOT 22.5 (all-2.5). */
    PetscCall(PicurvAssertRealNear(15.5, rep.lambda_c, 1e-9, "only x-direction modified to 2.5"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-3. One-sided viscous trigger by x-, y-, z-neighbor and applied only once for multiples. */
static PetscErrorCode TestMomentumStabilityOneSidedDirections(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    const PetscInt offs[3][3] = {{1,0,0},{0,1,0},{0,0,1}};   /* x,y,z neighbors of (4,4,4) */
    PetscFunctionBeginUser;
    for (int d = 0; d < 3; ++d) {
        PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
        simCtx->invicid = 0; simCtx->ren = 1.0;
        PetscCall(MomSetLocalScalarCell(user->da, user->lNvert,
                                        4+offs[d][0], 4+offs[d][1], 4+offs[d][2], 1.0));
        PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
        PetscCall(PicurvAssertRealNear(48.0, rep.lambda_v, 1e-9, "single-direction one-sided x2"));
        PetscCall(PicurvAssertIntEqual(1, (PetscInt)rep.one_sided, "one-sided flagged"));
        PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    }
    /* multiple solid neighbors still apply x2 only once (48, not 96). */
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 1.0;
    PetscCall(MomSetLocalScalarCell(user->da, user->lNvert, 5, 4, 4, 1.0));
    PetscCall(MomSetLocalScalarCell(user->da, user->lNvert, 4, 5, 4, 1.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertRealNear(48.0, rep.lambda_v, 1e-9, "multiple triggers -> x2 once (not 96)"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-4. Active-row mask: a cell with an inactive normal row but active tangential rows is
 *         retained with the full estimate; a fully-inactive location is excluded. */
static PetscErrorCode TestMomentumStabilityActiveRowMask(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    /* Positive non-periodic x face (i=mx-2) disables the xi row only; cell stays active. */
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, INLET));
    simCtx->central = 1; simCtx->invicid = 1;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 2.0, 2.0, 2.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertBool((PetscBool)(rep.active_cells > 0), "boundary cells remain active"));
    PetscCall(PicurvAssertBool((PetscBool)(rep.lambda > rep.lambda_t), "tangential rows retained"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    /* All cells solid -> zero active -> estimate falls back to lambda_t only. */
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 1;
    PetscCall(MomFillLocalScalar(user->da, user->lNvert, 1.0));   /* every cell solid */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    PetscCall(PicurvAssertIntEqual(0, rep.active_cells, "all solid -> zero active cells"));
    PetscCall(PicurvAssertRealNear(rep.lambda_t, rep.lambda, 1e-12, "empty active set -> lambda_t only"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-5. Explicit input validation: invalid block_number and a non-finite metric both error. */
static PetscErrorCode TestMomentumStabilityValidation(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscErrorCode e1 = 0, e2 = 0;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 1.0;

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    e1 = ComputeMomentumStabilityEstimate(user, 0, simCtx->dt, MOM_STAB_CAND_C, &rep);  /* block_number=0 */
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(e1 != 0), "block_number<=0 rejected"));

    PetscCall(MomSetLocalScalarCell(user->da, user->lAj, 4, 4, 4, PETSC_INFINITY));     /* bad metric */
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    e2 = ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(e2 != 0), "non-finite inverse-Jacobian rejected"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-6a. Candidate C exceeds B by the discrete-divergence term when div(Ucont) != 0. */
static PetscErrorCode TestMomentumStabilityCandidateBCDivergence(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->central = 1; simCtx->invicid = 1;                      /* f_c=1 isolates the div term */
    PetscCall(MomFillUcontXRamp(user->fda, user->lUcont, 2.0));    /* div = slope = 2 */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    /* C - B = 0.5*Aj*|div| = 0.5*2 = 1.0 at the (shared) controlling cell. */
    PetscCall(PicurvAssertRealNear(1.0, rep.lambda_C - rep.lambda_B, 1e-9, "C-B = divergence term"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-6b. Candidate D adds |grad u|_inf on a linear shear, and equals C for uniform velocity. */
static PetscErrorCode TestMomentumStabilityCandidateDShear(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->central = 1; simCtx->invicid = 1;
    PetscCall(MomFillUcatShearY(user->fda, user->lUcat, 2.0));     /* |grad u|_inf = gamma = 2 */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_D, &rep));
    PetscCall(PicurvAssertRealNear(2.0, rep.lambda_D - rep.lambda_C, 1e-9, "D-C = |grad u| on shear"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));

    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->central = 1; simCtx->invicid = 1;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcat, 0.7, -0.3, 0.5));  /* uniform velocity */
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_D, &rep));
    PetscCall(PicurvAssertRealNear(rep.lambda_C, rep.lambda_D, 1e-12, "D == C for uniform velocity"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-6c. Wall faces suppress eddy viscosity exactly as in Viscous(). */
static PetscErrorCode TestMomentumStabilityWallSuppression(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    PetscFunctionBeginUser;
    /* n=3 -> DMDA mx=4 -> interior cells {1,2}; every interior cell touches 3 wall faces
       (no wall-free interior cell exists, so the global max reflects suppression). */
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 3, WALL));
    simCtx->invicid = 0; simCtx->ren = 1.0; simCtx->les = 1;
    PetscCall(MomFillLocalScalar(user->da, user->lNu_t, 1.0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    /* 3 wall faces (nu_t suppressed, nu_eff=1) + 3 interior faces (nu_eff=2): 4*(3*1+3*2)=36. */
    PetscCall(PicurvAssertRealNear(36.0, rep.lambda_v, 1e-9, "wall-face eddy viscosity suppressed"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/* A3.5-7. Strengthened read-only: exact per-vector equality (not just norms). */
static PetscErrorCode TestMomentumStabilityReadOnlyExact(void)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; MomStabilityReport rep;
    Vec ucont0, ucat0, nvert0, nut0, aj0;
    PetscBool eq;
    PetscFunctionBeginUser;
    PetscCall(MomMakeUnitGrid(&simCtx, &user, 8, PERIODIC));
    simCtx->invicid = 0; simCtx->ren = 2.0; simCtx->les = 1;
    PetscCall(MomFillLocalCmpnts(user->fda, user->lUcont, 1.5, -0.7, 0.3));
    PetscCall(MomFillUcatShearY(user->fda, user->lUcat, 1.1));
    PetscCall(MomFillLocalScalar(user->da, user->lNu_t, 0.4));
    PetscCall(VecDuplicate(user->lUcont, &ucont0)); PetscCall(VecCopy(user->lUcont, ucont0));
    PetscCall(VecDuplicate(user->lUcat,  &ucat0));  PetscCall(VecCopy(user->lUcat,  ucat0));
    PetscCall(VecDuplicate(user->lNvert, &nvert0)); PetscCall(VecCopy(user->lNvert, nvert0));
    PetscCall(VecDuplicate(user->lNu_t,  &nut0));   PetscCall(VecCopy(user->lNu_t,  nut0));
    PetscCall(VecDuplicate(user->lAj,    &aj0));    PetscCall(VecCopy(user->lAj,    aj0));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_D, &rep));
    PetscCall(VecEqual(user->lUcont, ucont0, &eq)); PetscCall(PicurvAssertBool(eq, "lUcont bit-identical"));
    PetscCall(VecEqual(user->lUcat,  ucat0,  &eq)); PetscCall(PicurvAssertBool(eq, "lUcat bit-identical"));
    PetscCall(VecEqual(user->lNvert, nvert0, &eq)); PetscCall(PicurvAssertBool(eq, "lNvert bit-identical"));
    PetscCall(VecEqual(user->lNu_t,  nut0,   &eq)); PetscCall(PicurvAssertBool(eq, "lNu_t bit-identical"));
    PetscCall(VecEqual(user->lAj,    aj0,    &eq)); PetscCall(PicurvAssertBool(eq, "lAj bit-identical"));
    PetscCall(VecDestroy(&ucont0)); PetscCall(VecDestroy(&ucat0)); PetscCall(VecDestroy(&nvert0));
    PetscCall(VecDestroy(&nut0));   PetscCall(VecDestroy(&aj0));
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
        {"momentum-bdf-coefficient", TestMomentumBDFCoefficient},
        {"momentum-stability-centered-convection", TestMomentumStabilityCenteredConvection},
        {"momentum-stability-quick-interior", TestMomentumStabilityQuickInterior},
        {"momentum-stability-quick-boundary", TestMomentumStabilityQuickBoundary},
        {"momentum-stability-viscous-cartesian-factor8-estimate", TestMomentumStabilityViscousLongitudinal},
        {"momentum-stability-viscous-scales-with-nu", TestMomentumStabilityViscousScalesWithNu},
        {"momentum-stability-les-eddy-viscosity", TestMomentumStabilityLESEddyViscosity},
        {"momentum-stability-one-sided-viscous", TestMomentumStabilityOneSidedViscous},
        {"momentum-stability-read-only", TestMomentumStabilityReadOnly},
        {"momentum-stability-anisotropic-viscous", TestMomentumStabilityAnisotropicViscous},
        {"momentum-stability-directional-quick", TestMomentumStabilityDirectionalQuick},
        {"momentum-stability-one-sided-directions", TestMomentumStabilityOneSidedDirections},
        {"momentum-stability-active-row-mask", TestMomentumStabilityActiveRowMask},
        {"momentum-stability-validation", TestMomentumStabilityValidation},
        {"momentum-stability-candidate-bc-divergence", TestMomentumStabilityCandidateBCDivergence},
        {"momentum-stability-candidate-d-shear", TestMomentumStabilityCandidateDShear},
        {"momentum-stability-wall-suppression", TestMomentumStabilityWallSuppression},
        {"momentum-stability-read-only-exact", TestMomentumStabilityReadOnlyExact},
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
