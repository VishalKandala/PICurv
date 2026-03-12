/**
 * @file test_solver_kernels.c
 * @brief C unit tests for solver-side analytical and LES helper kernels.
 */

#include "test_support.h"

#include "AnalyticalSolutions.h"
#include "BodyForces.h"
#include "Filter.h"
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

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 8, 8));
    PetscCall(PicurvAssertBool(AnalyticalTypeRequiresCustomGeometry("TGV3D"),
                               "TGV3D should require custom geometry"));
    PetscCall(PicurvAssertBool((PetscBool)!AnalyticalTypeRequiresCustomGeometry("ZERO_FLOW"),
                               "ZERO_FLOW should not require custom geometry"));

    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "TGV3D", sizeof(simCtx->AnalyticalSolutionType)));
    simCtx->block_number = 1;
    PetscCall(SetAnalyticalGridInfo(user));
    PetscCall(PicurvAssertRealNear(0.0, user->Min_X, 1.0e-12, "TGV3D single-block xmin"));
    PetscCall(PicurvAssertRealNear(2.0 * PETSC_PI, user->Max_X, 1.0e-12, "TGV3D single-block xmax"));

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

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecDuplicate(user->Ucont, &rct));
    PetscCall(VecZeroEntries(rct));

    user->boundary_faces[0].face_id = BC_FACE_NEG_X;
    user->boundary_faces[0].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    simCtx->bulkVelocityCorrection = 2.0;
    simCtx->dt = 1.0;
    simCtx->forceScalingFactor = 1.0;
    simCtx->drivingForceMagnitude = 0.0;

    PetscCall(ComputeDrivenChannelFlowSource(user, rct));
    PetscCall(DMDAVecGetArrayRead(user->fda, rct, &rct_arr));
    PetscCall(PicurvAssertRealNear(1.5, rct_arr[1][1][1].x, 1.0e-12,
                                   "driven flow source should update the x component"));
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
