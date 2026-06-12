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
 * @brief Marks the y faces as periodic for periodic-transfer harnesses.
 */
static void MarkYPeriodic(UserCtx *user)
{
    user->boundary_faces[BC_FACE_NEG_Y].face_id = BC_FACE_NEG_Y;
    user->boundary_faces[BC_FACE_NEG_Y].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_Y].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;
    user->boundary_faces[BC_FACE_POS_Y].face_id = BC_FACE_POS_Y;
    user->boundary_faces[BC_FACE_POS_Y].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_Y].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;
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
 * @brief Tests ordered cell-field synchronization across two periodic directions.
 */
static PetscErrorCode TestSynchronizePeriodicCellFieldsCopiesMixedAxes(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***p = NULL;
    PetscReal ***cs = NULL;
    PetscReal ***diffusivity = NULL;
    Cmpnts ***ucat = NULL;
    PetscInt mx, my;
    const char *fields[] = {"Ucat", "P", "Phi", "Nvert", "Nu_t", "CS", "Diffusivity"};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_TRUE, PETSC_FALSE));
    MarkXPeriodic(user);
    MarkYPeriodic(user);

    PetscCall(DMCreateGlobalVector(user->da, &user->Nu_t));
    PetscCall(DMCreateLocalVector(user->da, &user->lNu_t));
    PetscCall(VecSet(user->Nu_t, 5.0));
    PetscCall(VecSet(user->lNu_t, 0.0));
    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    PetscCall(DMCreateLocalVector(user->da, &user->lCs));
    PetscCall(VecSet(user->Phi, 3.0));
    PetscCall(VecSet(user->Nvert, 4.0));

    PetscCall(DMDAVecGetArray(user->da, user->P, &p));
    PetscCall(DMDAVecGetArray(user->da, user->CS, &cs));
    PetscCall(DMDAVecGetArray(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                PetscReal value = (PetscReal)(i + 10 * j + 100 * k);
                p[k][j][i] = value;
                cs[k][j][i] = value + 4000.0;
                diffusivity[k][j][i] = value + 5000.0;
                ucat[k][j][i].x = value + 1000.0;
                ucat[k][j][i].y = value + 2000.0;
                ucat[k][j][i].z = value + 3000.0;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArray(user->da, user->CS, &cs));
    PetscCall(DMDAVecRestoreArray(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));

    PetscCall(SynchronizePeriodicCellFields(user, 7, fields));

    mx = user->info.mx;
    my = user->info.my;
    PetscCall(DMDAVecGetArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecGetArrayRead(user->da, user->CS, &cs));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 20 + 200), p[2][2][0], 1.0e-12,
                                    "x-periodic negative endpoint should copy the opposite interior value"));
    PetscCall(PicurvAssertRealNear((PetscReal)(1 + 20 + 200), p[2][2][mx - 1], 1.0e-12,
                                    "x-periodic positive endpoint should copy the leading interior value"));
    PetscCall(PicurvAssertRealNear((PetscReal)(2 + 10 * (my - 2) + 200), p[2][0][2], 1.0e-12,
                                    "y-periodic negative endpoint should copy the opposite interior value"));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 10 * (my - 2) + 200), p[2][0][0], 1.0e-12,
                                    "ordered synchronization should propagate the opposite periodic corner"));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 10 * (my - 2) + 1200), ucat[2][0][0].x, 1.0e-12,
                                    "vector fields should use the same ordered periodic-corner protocol"));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 10 * (my - 2) + 4200), cs[2][0][0], 1.0e-12,
                                    "CS should use the ordered periodic-corner protocol"));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 10 * (my - 2) + 5200), diffusivity[2][0][0], 1.0e-12,
                                    "Diffusivity should use the ordered periodic-corner protocol"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->CS, &cs));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(PicurvAssertVecConstant(user->Nu_t, 5.0, 1.0e-12, "Nu_t should be accepted by periodic cell synchronization"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests mixed-boundary post-projection cell finalization and Ucont preservation.
 */
static PetscErrorCode TestFinalizePostProjectionCellFieldsMixedBoundaries(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec ucont_before = NULL;
    Vec ucont_difference = NULL;
    Cmpnts ***ucat = NULL;
    Cmpnts ***ucont = NULL;
    PetscReal ***p = NULL;
    PetscReal difference_norm = 0.0;
    PetscInt mx;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(VecSet(user->Bcs.Ubcs, 0.0));

    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecGetArray(user->da, user->P, &p));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                PetscReal value = (PetscReal)(i + 10 * j + 100 * k);
                ucat[k][j][i].x = value + 1000.0;
                ucat[k][j][i].y = value + 2000.0;
                ucat[k][j][i].z = value + 3000.0;
                ucont[k][j][i].x = value + 4000.0;
                ucont[k][j][i].y = value + 5000.0;
                ucont[k][j][i].z = value + 6000.0;
                p[k][j][i] = value;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));

    PetscCall(VecDuplicate(user->Ucont, &ucont_before));
    PetscCall(VecDuplicate(user->Ucont, &ucont_difference));
    PetscCall(VecCopy(user->Ucont, ucont_before));

    PetscCall(FinalizePostProjectionCellFields(user));

    PetscCall(VecWAXPY(ucont_difference, -1.0, user->Ucont, ucont_before));
    PetscCall(VecNorm(ucont_difference, NORM_INFINITY, &difference_norm));
    PetscCall(PicurvAssertRealNear(0.0, difference_norm, 1.0e-12,
                                    "post-projection cell finalization must preserve mixed-boundary Ucont"));

    mx = user->info.mx;
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArrayRead(user->da, user->P, &p));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 20 + 200 + 1000), ucat[2][2][0].x, 1.0e-12,
                                    "mixed finalization should restore the negative x-periodic Ucat endpoint"));
    PetscCall(PicurvAssertRealNear((PetscReal)(-(2 + 10 + 200 + 1000)), ucat[2][0][2].x, 1.0e-12,
                                    "mixed finalization should extrapolate the non-periodic y dummy face"));
    PetscCall(PicurvAssertRealNear((PetscReal)(-((mx - 2) + 10 + 200 + 1000)), ucat[2][0][0].x, 1.0e-12,
                                    "periodic synchronization should win at a mixed x-periodic/y-physical edge"));
    PetscCall(PicurvAssertRealNear((PetscReal)((mx - 2) + 20 + 200), p[2][2][0], 1.0e-12,
                                    "mixed finalization should restore the negative x-periodic pressure endpoint"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, &p));

    PetscCall(VecDestroy(&ucont_difference));
    PetscCall(VecDestroy(&ucont_before));
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
        {"synchronize-periodic-cell-fields-copies-mixed-axes", TestSynchronizePeriodicCellFieldsCopiesMixedAxes},
        {"finalize-post-projection-cell-fields-mixed-boundaries", TestFinalizePostProjectionCellFieldsMixedBoundaries},
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
