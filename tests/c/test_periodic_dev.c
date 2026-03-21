/**
 * @file test_periodic_dev.c
 * @brief Non-gating periodic-boundary development harnesses.
 */

#include "test_support.h"

#include "Boundaries.h"

/**
 * @brief Appends one key/value pair to a linked list of boundary-condition parameters.
 */
static PetscErrorCode AppendBCParam(BC_Param **head, const char *key, const char *value)
{
    BC_Param *node = NULL;
    BC_Param **cursor = head;

    PetscFunctionBeginUser;
    PetscCheck(head != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BC param list head cannot be NULL.");
    PetscCall(PetscCalloc1(1, &node));
    PetscCall(PetscStrallocpy(key, &node->key));
    PetscCall(PetscStrallocpy(value, &node->value));
    while (*cursor) {
        cursor = &(*cursor)->next;
    }
    *cursor = node;
    PetscFunctionReturn(0);
}
/**
 * @brief Destroys one boundary-condition handler allocated by a periodic test.
 */
static PetscErrorCode DestroyBoundaryHandler(BoundaryCondition **bc_ptr)
{
    BoundaryCondition *bc = NULL;

    PetscFunctionBeginUser;
    if (!bc_ptr || !*bc_ptr) PetscFunctionReturn(0);

    bc = *bc_ptr;
    if (bc->Destroy) {
        PetscCall(bc->Destroy(bc));
    }
    PetscCall(PetscFree(bc));
    *bc_ptr = NULL;
    PetscFunctionReturn(0);
}
/**
 * @brief Marks the x faces as periodic for periodic-transfer harnesses.
 */
static void MarkXPeriodic(UserCtx *user)
{
    user->boundary_faces[BC_FACE_NEG_X].face_id = BC_FACE_NEG_X;
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;
    user->boundary_faces[BC_FACE_POS_X].face_id = BC_FACE_POS_X;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_X].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;
}
/**
 * @brief Tests periodic geometric factory construction in the non-gating periodic harness.
 */
static PetscErrorCode TestPeriodicGeometricFactoryAssignment(void)
{
    BoundaryCondition *bc = NULL;

    PetscFunctionBeginUser;
    PetscCall(BoundaryCondition_Create(BC_HANDLER_PERIODIC_GEOMETRIC, &bc));
    PetscCall(PicurvAssertBool((PetscBool)(bc != NULL), "periodic geometric factory should allocate a handler"));
    PetscCall(PicurvAssertIntEqual(BC_PRIORITY_WALL, bc->priority, "periodic geometric handler priority"));
    PetscCall(PicurvAssertBool((PetscBool)(bc->Apply == NULL), "periodic geometric handler should not expose Apply"));
    PetscCall(PicurvAssertBool((PetscBool)(bc->Initialize == NULL), "periodic geometric handler should not expose Initialize"));
    PetscCall(DestroyBoundaryHandler(&bc));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests direct periodic face transfer on the staggered velocity field.
 */
static PetscErrorCode TestTransferPeriodicFaceFieldCopiesXFaces(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucont = NULL;
    PetscReal expected_neg_face = 0.0;
    PetscReal expected_pos_face = 0.0;
    PetscReal expected_neg_ghost = 0.0;
    PetscReal expected_pos_ghost = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);

    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                ucont[k][j][i].x = (PetscReal)i;
                ucont[k][j][i].y = 10.0 + (PetscReal)i;
                ucont[k][j][i].z = 20.0 + (PetscReal)i;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont));
    expected_neg_face = ucont[2][2][user->info.mx - 2].x;
    expected_pos_face = ucont[2][2][1].x;
    expected_neg_ghost = ucont[2][2][user->info.mx - 3].x;
    expected_pos_ghost = ucont[2][2][2].x;
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont));

    PetscCall(TransferPeriodicFaceField(user, "Ucont"));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont));
    PetscCall(PicurvAssertRealNear(expected_neg_face, ucont[2][2][0].x, 1.0e-12, "NEG_X periodic face should copy from the opposite interior face"));
    PetscCall(PicurvAssertRealNear(expected_pos_face, ucont[2][2][user->info.mx - 1].x, 1.0e-12, "POS_X periodic face should copy from the leading interior face"));
    PetscCall(PicurvAssertRealNear(expected_neg_ghost, ucont[2][2][-1].x, 1.0e-12, "NEG_X ghost face should copy the second-to-last interior face"));
    PetscCall(PicurvAssertRealNear(expected_pos_ghost, ucont[2][2][user->info.mx].x, 1.0e-12, "POS_X ghost face should copy the second interior face"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests periodic metric transfer through the aggregate periodic-metrics helper.
 */
static PetscErrorCode TestApplyMetricsPeriodicBCsCopiesAjFaces(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***aj = NULL;
    PetscReal expected_neg_face = 0.0;
    PetscReal expected_pos_face = 0.0;
    PetscReal expected_neg_ghost = 0.0;
    PetscReal expected_pos_ghost = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);

    PetscCall(DMDAVecGetArray(user->da, user->Aj, &aj));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                aj[k][j][i] = 100.0 + (PetscReal)i;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->Aj, &aj));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
    expected_neg_face = aj[2][2][user->info.mx - 2];
    expected_pos_face = aj[2][2][1];
    expected_neg_ghost = aj[2][2][user->info.mx - 3];
    expected_pos_ghost = aj[2][2][2];
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));

    PetscCall(ApplyMetricsPeriodicBCs(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
    PetscCall(PicurvAssertRealNear(expected_neg_face, aj[2][2][0], 1.0e-12, "NEG_X periodic metric face should copy from the opposite interior face"));
    PetscCall(PicurvAssertRealNear(expected_pos_face, aj[2][2][user->info.mx - 1], 1.0e-12, "POS_X periodic metric face should copy from the leading interior face"));
    PetscCall(PicurvAssertRealNear(expected_neg_ghost, aj[2][2][-1], 1.0e-12, "NEG_X periodic metric ghost should copy the second-to-last interior face"));
    PetscCall(PicurvAssertRealNear(expected_pos_ghost, aj[2][2][user->info.mx], 1.0e-12, "POS_X periodic metric ghost should copy the second interior face"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests periodic driven-flow controller initialization, sensing, and trim application.
 */
static PetscErrorCode TestPeriodicDrivenConstantHandlerBehavior(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundaryCondition *bc = NULL;
    BCContext ctx;
    PetscReal dummy_inflow = 0.0;
    PetscReal dummy_outflow = 0.0;
    Cmpnts ***ucont = NULL;
    Cmpnts ***uch = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 6, 6, PETSC_FALSE, PETSC_FALSE, PETSC_TRUE));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Ucont, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));

    user->boundary_faces[BC_FACE_NEG_Z].face_id = BC_FACE_NEG_Z;
    user->boundary_faces[BC_FACE_NEG_Z].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_Z].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    user->boundary_faces[BC_FACE_POS_Z].face_id = BC_FACE_POS_Z;
    user->boundary_faces[BC_FACE_POS_Z].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_Z].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    PetscCall(AppendBCParam(&user->boundary_faces[BC_FACE_NEG_Z].params, "target_flux", "30.0"));
    PetscCall(AppendBCParam(&user->boundary_faces[BC_FACE_NEG_Z].params, "apply_trim", "true"));

    ctx.user = user;
    ctx.face_id = BC_FACE_NEG_Z;
    PetscCall(BoundaryCondition_Create(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX, &bc));
    PetscCall(bc->Initialize(bc, &ctx));
    PetscCall(bc->PreStep(bc, &ctx, &dummy_inflow, &dummy_outflow));
    PetscCall(PicurvAssertRealNear(30.0, simCtx->targetVolumetricFlux, 1.0e-12, "periodic driven initialization should store the target flux"));
    PetscCall(PicurvAssertRealNear(0.2, simCtx->bulkVelocityCorrection, 1.0e-12, "periodic driven controller should compute the bulk correction from the measured flux"));

    PetscCall(bc->Apply(bc, &ctx));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Uch, &uch));
    PetscCall(PicurvAssertRealNear(1.2, ucont[0][3][3].z, 1.0e-12, "periodic driven trim should update the boundary face flux"));
    PetscCall(PicurvAssertRealNear(0.2, uch[0][3][3].z, 1.0e-12, "periodic driven trim should be recorded in Uch"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Uch, &uch));

    PetscCall(DestroyBoundaryHandler(&bc));
    FreeBC_ParamList(user->boundary_faces[BC_FACE_NEG_Z].params);
    user->boundary_faces[BC_FACE_NEG_Z].params = NULL;
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the periodic driven handler rejects non-periodic faces during initialization.
 */
static PetscErrorCode TestPeriodicDrivenConstantRejectsNonPeriodicFace(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundaryCondition *bc = NULL;
    BCContext ctx;
    PetscErrorCode ierr_init = 0;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    user->boundary_faces[BC_FACE_NEG_Z].face_id = BC_FACE_NEG_Z;
    user->boundary_faces[BC_FACE_NEG_Z].mathematical_type = WALL;
    user->boundary_faces[BC_FACE_NEG_Z].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    PetscCall(AppendBCParam(&user->boundary_faces[BC_FACE_NEG_Z].params, "target_flux", "5.0"));
    ctx.user = user;
    ctx.face_id = BC_FACE_NEG_Z;

    PetscCall(BoundaryCondition_Create(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX, &bc));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_init = bc->Initialize(bc, &ctx);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_init != 0), "periodic driven initialization should reject non-periodic faces"));

    PetscCall(DestroyBoundaryHandler(&bc));
    FreeBC_ParamList(user->boundary_faces[BC_FACE_NEG_Z].params);
    user->boundary_faces[BC_FACE_NEG_Z].params = NULL;
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the non-gating periodic development PETSc test binary.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"periodic-geometric-factory-assignment", TestPeriodicGeometricFactoryAssignment},
        {"transfer-periodic-face-field-copies-x-faces", TestTransferPeriodicFaceFieldCopiesXFaces},
        {"apply-metrics-periodic-bcs-copies-aj-faces", TestApplyMetricsPeriodicBCsCopiesAjFaces},
        {"periodic-driven-constant-handler-behavior", TestPeriodicDrivenConstantHandlerBehavior},
        {"periodic-driven-constant-rejects-non-periodic-face", TestPeriodicDrivenConstantRejectsNonPeriodicFace},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv periodic development tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-periodic-dev", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
