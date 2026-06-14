/**
 * @file test_periodic_dev.c
 * @brief Focused geometric-periodic boundary tests.
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
 * @brief Tests periodic configuration rejects an unpaired geometric periodic face.
 */
static PetscErrorCode TestPeriodicConfigurationRequiresPairedFaces(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscErrorCode ierr_validate = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_validate = BoundarySystem_Validate(user);
    PetscCall(PetscPopErrorHandler());

    PetscCall(PicurvAssertBool((PetscBool)(ierr_validate != 0), "unpaired periodic faces should be rejected"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests constant translational periodic geometry is accepted and stored.
 */
static PetscErrorCode TestPeriodicGeometryStoresTranslation(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec gcoor = NULL;
    const Cmpnts ***coor = NULL;
    PetscReal expected_x = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);

    PetscCall(DMGetCoordinates(user->da, &gcoor));
    PetscCall(DMDAVecGetArrayRead(user->fda, gcoor, &coor));
    expected_x = coor[0][0][user->info.mx - 2].x - coor[0][0][0].x;
    PetscCall(DMDAVecRestoreArrayRead(user->fda, gcoor, &coor));

    PetscCall(ValidatePeriodicGeometry(user));
    PetscCall(PicurvAssertBool(user->periodic_translation_valid[0], "X-periodic translation should be marked valid"));
    PetscCall(PicurvAssertBool((PetscBool)!user->periodic_translation_valid[1], "Y translation should remain invalid"));
    PetscCall(PicurvAssertBool((PetscBool)!user->periodic_translation_valid[2], "Z translation should remain invalid"));
    PetscCall(PicurvAssertRealNear(expected_x, user->periodic_translation[0].x, 1.0e-12, "stored X translation"));
    PetscCall(PicurvAssertRealNear(0.0, user->periodic_translation[0].y, 1.0e-12, "stored X translation y component"));
    PetscCall(PicurvAssertRealNear(0.0, user->periodic_translation[0].z, 1.0e-12, "stored X translation z component"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests mixed periodic directions are validated and stored independently.
 */
static PetscErrorCode TestPeriodicGeometryStoresMixedTranslations(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 6, 4, PETSC_TRUE, PETSC_TRUE, PETSC_FALSE));
    MarkXPeriodic(user);
    MarkYPeriodic(user);

    PetscCall(ValidatePeriodicGeometry(user));
    PetscCall(PicurvAssertBool(user->periodic_translation_valid[0], "mixed-periodic X translation should be valid"));
    PetscCall(PicurvAssertBool(user->periodic_translation_valid[1], "mixed-periodic Y translation should be valid"));
    PetscCall(PicurvAssertBool((PetscBool)!user->periodic_translation_valid[2], "mixed-periodic Z translation should remain invalid"));
    PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(user->periodic_translation[0].x) > 0.0),
                               "mixed-periodic X translation should have nonzero x component"));
    PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(user->periodic_translation[1].y) > 0.0),
                               "mixed-periodic Y translation should have nonzero y component"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests non-translational seam geometry fails before metric construction.
 */
static PetscErrorCode TestPeriodicGeometryRejectsVaryingTranslation(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec gcoor = NULL;
    Cmpnts ***coor = NULL;
    PetscErrorCode ierr_validate = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);

    PetscCall(DMGetCoordinates(user->da, &gcoor));
    PetscCall(DMDAVecGetArray(user->fda, gcoor, &coor));
    coor[1][1][user->info.mx - 2].y += 0.25;
    PetscCall(DMDAVecRestoreArray(user->fda, gcoor, &coor));
    PetscCall(UpdateLocalGhosts(user, "Coordinates"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_validate = ValidatePeriodicGeometry(user);
    PetscCall(PetscPopErrorHandler());

    PetscCall(PicurvAssertBool((PetscBool)(ierr_validate != 0), "varying seam translation should be rejected"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests face-center coordinates remain geometrically continuous across a periodic seam.
 */
static PetscErrorCode TestPeriodicFaceCenterCoordinateSynchronization(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***centx = NULL;
    const Cmpnts ***lcentx = NULL;
    const char *fields[] = {"Centx"};
    PetscReal translation, spacing;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(ValidatePeriodicGeometry(user));
    translation = user->periodic_translation[0].x;
    spacing = translation / (PetscReal)(user->info.mx - 2);

    PetscCall(DMDAVecGetArray(user->fda, user->Centx, &centx));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) {
                centx[k][j][i] = (Cmpnts){spacing * i, 2.0 + j, 3.0 + k};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Centx, &centx));

    PetscCall(SynchronizePeriodicFaceFields(user, 'i', 1, fields));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCentx, &lcentx));
    PetscCall(PicurvAssertRealNear(-spacing, lcentx[2][2][-1].x, 1.0e-12,
                                   "translated Centx negative adjacent ghost"));
    PetscCall(PicurvAssertRealNear(0.0, lcentx[2][2][0].x, 1.0e-12,
                                   "translated Centx negative endpoint"));
    PetscCall(PicurvAssertRealNear(translation + spacing, lcentx[2][2][user->info.mx - 1].x, 1.0e-12,
                                   "translated Centx positive endpoint"));
    PetscCall(PicurvAssertRealNear(translation + 2.0 * spacing, lcentx[2][2][user->info.mx].x, 1.0e-12,
                                   "translated Centx positive adjacent ghost"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCentx, &lcentx));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests the dedicated QUICK outer-ghost repair for cell-centered inputs.
 */
static PetscErrorCode TestPeriodicQuickStencilPreparation(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucat = NULL;
    PetscReal ***nvert = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 6, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);

    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) {
                ucat[k][j][i] = (Cmpnts){100.0 + i, 200.0 + i, 300.0 + i};
                nvert[k][j][i] = 400.0 + i;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(UpdateLocalGhosts(user, "Ucat"));
    PetscCall(UpdateLocalGhosts(user, "Nvert"));
    PetscCall(PreparePeriodicQuickStencilFields(user, user->lUcat, user->lNvert));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcat, &ucat));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, &nvert));
    PetscCall(PicurvAssertRealNear(100.0 + user->info.mx - 3, ucat[2][2][-1].x, 1.0e-12,
                                   "QUICK Ucat negative outer ghost"));
    PetscCall(PicurvAssertRealNear(102.0, ucat[2][2][user->info.mx].x, 1.0e-12,
                                   "QUICK Ucat positive outer ghost"));
    PetscCall(PicurvAssertRealNear(400.0 + user->info.mx - 3, nvert[2][2][-1], 1.0e-12,
                                   "QUICK Nvert negative outer ghost"));
    PetscCall(PicurvAssertRealNear(402.0, nvert[2][2][user->info.mx], 1.0e-12,
                                   "QUICK Nvert positive outer ghost"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, &nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcat, &ucat));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests periodic geometric factory construction.
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
static PetscErrorCode TestApplyMetricsPeriodicBCsSynchronizesAj(void)
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
    PetscCall(ValidatePeriodicGeometry(user));

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
    expected_neg_ghost = expected_pos_face;
    expected_pos_ghost = expected_neg_face;
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));

    PetscCall(ApplyMetricsPeriodicBCs(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
    PetscCall(PicurvAssertRealNear(expected_neg_face, aj[2][2][0], 1.0e-12, "NEG_X periodic metric face should copy from the opposite interior face"));
    PetscCall(PicurvAssertRealNear(expected_pos_face, aj[2][2][user->info.mx - 1], 1.0e-12, "POS_X periodic metric face should copy from the leading interior face"));
    PetscCall(PicurvAssertRealNear(expected_neg_ghost, aj[2][2][-1], 1.0e-12, "cell-centered Aj negative ghost should retain PETSc wraparound"));
    PetscCall(PicurvAssertRealNear(expected_pos_ghost, aj[2][2][user->info.mx], 1.0e-12, "cell-centered Aj positive ghost should retain PETSc wraparound"));
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
 * @brief Tests ordered persistent synchronization for an I-face scalar metric.
 */
static PetscErrorCode TestSynchronizePeriodicFaceFieldsCopiesMixedAxes(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***iaj = NULL;
    Cmpnts ***csi = NULL;
    const char *fields[] = {"Csi", "IAj"};
    PetscInt mx, my;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_TRUE, PETSC_FALSE));
    MarkXPeriodic(user);
    MarkYPeriodic(user);
    mx = user->info.mx;
    my = user->info.my;

    PetscCall(DMDAVecGetArray(user->da, user->IAj, &iaj));
    PetscCall(DMDAVecGetArray(user->fda, user->Csi, &csi));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) {
                PetscReal value = (PetscReal)(i + 10 * j + 100 * k);
                iaj[k][j][i] = value;
                csi[k][j][i] = (Cmpnts){value + 1000.0, value + 2000.0, value + 3000.0};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->IAj, &iaj));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Csi, &csi));

    PetscCall(SynchronizePeriodicFaceFields(user, 'i', 2, fields));

    PetscCall(DMDAVecGetArrayRead(user->da, user->IAj, &iaj));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Csi, &csi));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 2 + 20 + 200), iaj[2][2][0], 1.0e-12,
                                   "I-face negative seam should copy the opposite physical seam face"));
    PetscCall(PicurvAssertRealNear((PetscReal)(1 + 20 + 200), iaj[2][2][mx - 1], 1.0e-12,
                                   "I-face positive dummy should copy the leading physical face"));
    PetscCall(PicurvAssertRealNear((PetscReal)(2 + 10 * (my - 2) + 200), iaj[2][0][2], 1.0e-12,
                                   "I-face field should use cell-style synchronization tangentially"));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 2 + 10 * (my - 2) + 200), iaj[2][0][0], 1.0e-12,
                                   "ordered face synchronization should propagate mixed-axis corners"));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 2 + 10 * (my - 2) + 1200), csi[2][0][0].x, 1.0e-12,
                                   "vector I-face fields should use the same ordered synchronization"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->IAj, &iaj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Csi, &csi));

    PetscCall(DMDAVecGetArrayRead(user->da, user->lIAj, &iaj));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 3 + 20 + 200), iaj[2][2][-1], 1.0e-12,
                                   "I-face negative adjacent ghost should be shifted to the prior face"));
    PetscCall(PicurvAssertRealNear((PetscReal)(2 + 20 + 200), iaj[2][2][mx], 1.0e-12,
                                   "I-face positive adjacent ghost should be shifted to the next face"));
    PetscCall(PicurvAssertRealNear((PetscReal)(2 + 10 + 200), iaj[2][-1][2], 1.0e-12,
                                   "I-face tangential ghosts should retain cell-style wraparound"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lIAj, &iaj));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests ordered persistent synchronization for component-staggered Ucont.
 */
static PetscErrorCode TestSynchronizePeriodicStaggeredFieldsCopiesMixedAxes(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucont = NULL;
    const char *fields[] = {"Ucont"};
    PetscInt mx, my;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 4, 4, 4, PETSC_TRUE, PETSC_TRUE, PETSC_FALSE));
    MarkXPeriodic(user);
    MarkYPeriodic(user);
    mx = user->info.mx;
    my = user->info.my;

    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) {
                PetscReal value = (PetscReal)(i + 10 * j + 100 * k);
                ucont[k][j][i] = (Cmpnts){value + 1000.0, value + 2000.0, value + 3000.0};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));

    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 2 + 20 + 1200), ucont[2][2][0].x, 1.0e-12,
                                   "Ucont.x negative X seam should copy the opposite physical seam"));
    PetscCall(PicurvAssertRealNear((PetscReal)(1 + 20 + 2200), ucont[2][2][mx - 1].y, 1.0e-12,
                                   "Ucont.y positive X dummy should copy the leading physical value"));
    PetscCall(PicurvAssertRealNear((PetscReal)(2 + 10 * (my - 2) + 3200), ucont[2][0][2].z, 1.0e-12,
                                   "Ucont.z should synchronize tangentially in Y"));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 2 + 10 * (my - 2) + 2200), ucont[2][0][0].y, 1.0e-12,
                                   "ordered staggered synchronization should propagate mixed-axis corners"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 3 + 20 + 1200), ucont[2][2][-1].x, 1.0e-12,
                                   "Ucont.x should use face-normal X ghost repair"));
    PetscCall(PicurvAssertRealNear((PetscReal)(1 + 20 + 2200), ucont[2][2][-1].y, 1.0e-12,
                                   "Ucont.y should retain cell-style wraparound tangentially in X"));
    PetscCall(PicurvAssertRealNear((PetscReal)(2 + 20 + 1200), ucont[2][2][mx].x, 1.0e-12,
                                   "Ucont.x positive adjacent ghost should use the next face"));
    PetscCall(PicurvAssertRealNear((PetscReal)(mx - 2 + 20 + 2200), ucont[2][2][mx].y, 1.0e-12,
                                   "Ucont.y positive X ghost should retain cell-style wraparound"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont));

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
 * @brief Runs the focused geometric-periodic PETSc test binary.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"periodic-geometric-factory-assignment", TestPeriodicGeometricFactoryAssignment},
        {"periodic-configuration-requires-paired-faces", TestPeriodicConfigurationRequiresPairedFaces},
        {"periodic-geometry-stores-translation", TestPeriodicGeometryStoresTranslation},
        {"periodic-geometry-stores-mixed-translations", TestPeriodicGeometryStoresMixedTranslations},
        {"periodic-geometry-rejects-varying-translation", TestPeriodicGeometryRejectsVaryingTranslation},
        {"periodic-face-center-coordinate-synchronization", TestPeriodicFaceCenterCoordinateSynchronization},
        {"periodic-quick-stencil-preparation", TestPeriodicQuickStencilPreparation},
        {"transfer-periodic-face-field-copies-x-faces", TestTransferPeriodicFaceFieldCopiesXFaces},
        {"apply-metrics-periodic-bcs-synchronizes-aj", TestApplyMetricsPeriodicBCsSynchronizesAj},
        {"synchronize-periodic-cell-fields-copies-mixed-axes", TestSynchronizePeriodicCellFieldsCopiesMixedAxes},
        {"synchronize-periodic-face-fields-copies-mixed-axes", TestSynchronizePeriodicFaceFieldsCopiesMixedAxes},
        {"synchronize-periodic-staggered-fields-copies-mixed-axes", TestSynchronizePeriodicStaggeredFieldsCopiesMixedAxes},
        {"finalize-post-projection-cell-fields-mixed-boundaries", TestFinalizePostProjectionCellFieldsMixedBoundaries},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv geometric-periodic tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-periodic", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
