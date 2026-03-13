/**
 * @file test_boundaries.c
 * @brief C unit tests for boundary factories, inlet ownership, and face helpers.
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
 * @brief Destroys one boundary-condition handler allocated by a boundary test.
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
 * @brief Chooses a stable interior node index for face-sample assertions on tiny test grids.
 */
static PetscInt InteriorSampleIndex(PetscInt nodes)
{
    PetscInt index = 3;

    if (index >= nodes - 1) {
        index = nodes - 2;
    }
    if (index < 1) {
        index = 1;
    }
    return index;
}

/**
 * @brief Extracts the velocity component aligned with the supplied boundary-face normal.
 */
static PetscReal GetFaceNormalComponent(Cmpnts value, BCFace face)
{
    switch (face) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X:
            return value.x;
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y:
            return value.y;
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z:
            return value.z;
    }
    return 0.0;
}

/**
 * @brief Returns the sign convention used for face-normal flux expectations.
 */
static PetscReal GetFaceOrientationSign(BCFace face)
{
    switch (face) {
        case BC_FACE_NEG_X:
        case BC_FACE_NEG_Y:
        case BC_FACE_NEG_Z:
            return 1.0;
        case BC_FACE_POS_X:
        case BC_FACE_POS_Y:
        case BC_FACE_POS_Z:
            return -1.0;
    }
    return 0.0;
}

/**
 * @brief Maps an inlet face to the matching configuration key used by the handler parser.
 */
static const char *GetInletParamKey(BCFace face)
{
    switch (face) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X:
            return "vx";
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y:
            return "vy";
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z:
            return "vz";
    }
    return "";
}

/**
 * @brief Computes the number of face-interior sample points for a given boundary face.
 */
static PetscReal GetFaceInteriorPointCount(const UserCtx *user, BCFace face)
{
    switch (face) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X:
            return (PetscReal)(user->info.my - 2) * (PetscReal)(user->info.mz - 2);
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y:
            return (PetscReal)(user->info.mx - 2) * (PetscReal)(user->info.mz - 2);
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z:
            return (PetscReal)(user->info.mx - 2) * (PetscReal)(user->info.my - 2);
    }
    return 0.0;
}

/**
 * @brief Selects representative `Ucont` and `Ubcs` slots for face-matrix assertions.
 */
static PetscErrorCode GetRepresentativeFaceSlots(const UserCtx *user,
                                                 BCFace face,
                                                 PetscInt *ucont_k,
                                                 PetscInt *ucont_j,
                                                 PetscInt *ucont_i,
                                                 PetscInt *ubcs_k,
                                                 PetscInt *ubcs_j,
                                                 PetscInt *ubcs_i)
{
    const PetscInt sample_i = InteriorSampleIndex(user->info.mx);
    const PetscInt sample_j = InteriorSampleIndex(user->info.my);
    const PetscInt sample_k = InteriorSampleIndex(user->info.mz);

    PetscFunctionBeginUser;
    switch (face) {
        case BC_FACE_NEG_X:
            *ucont_k = sample_k; *ucont_j = sample_j; *ucont_i = 0;
            *ubcs_k = sample_k; *ubcs_j = sample_j; *ubcs_i = 0;
            break;
        case BC_FACE_POS_X:
            *ucont_k = sample_k; *ucont_j = sample_j; *ucont_i = user->info.mx - 2;
            *ubcs_k = sample_k; *ubcs_j = sample_j; *ubcs_i = user->info.mx - 1;
            break;
        case BC_FACE_NEG_Y:
            *ucont_k = sample_k; *ucont_j = 0; *ucont_i = sample_i;
            *ubcs_k = sample_k; *ubcs_j = 0; *ubcs_i = sample_i;
            break;
        case BC_FACE_POS_Y:
            *ucont_k = sample_k; *ucont_j = user->info.my - 2; *ucont_i = sample_i;
            *ubcs_k = sample_k; *ubcs_j = user->info.my - 1; *ubcs_i = sample_i;
            break;
        case BC_FACE_NEG_Z:
            *ucont_k = 0; *ucont_j = sample_j; *ucont_i = sample_i;
            *ubcs_k = 0; *ubcs_j = sample_j; *ubcs_i = sample_i;
            break;
        case BC_FACE_POS_Z:
            *ucont_k = user->info.mz - 2; *ucont_j = sample_j; *ucont_i = sample_i;
            *ubcs_k = user->info.mz - 1; *ubcs_j = sample_j; *ubcs_i = sample_i;
            break;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Selects center, off-center, and wall-adjacent slots for parabolic-face checks.
 */
static PetscErrorCode GetParabolicSampleSlots(const UserCtx *user,
                                              BCFace face,
                                              PetscInt *center_k,
                                              PetscInt *center_j,
                                              PetscInt *center_i,
                                              PetscInt *off_k,
                                              PetscInt *off_j,
                                              PetscInt *off_i,
                                              PetscInt *wall_k,
                                              PetscInt *wall_j,
                                              PetscInt *wall_i,
                                              PetscInt *ubcs_center_k,
                                              PetscInt *ubcs_center_j,
                                              PetscInt *ubcs_center_i)
{
    PetscFunctionBeginUser;
    switch (face) {
        case BC_FACE_NEG_X:
            *center_k = 3; *center_j = 3; *center_i = 0;
            *off_k = 2; *off_j = 2; *off_i = 0;
            *wall_k = 1; *wall_j = 1; *wall_i = 0;
            *ubcs_center_k = 3; *ubcs_center_j = 3; *ubcs_center_i = 0;
            break;
        case BC_FACE_POS_X:
            *center_k = 3; *center_j = 3; *center_i = user->info.mx - 2;
            *off_k = 2; *off_j = 2; *off_i = user->info.mx - 2;
            *wall_k = 1; *wall_j = 1; *wall_i = user->info.mx - 2;
            *ubcs_center_k = 3; *ubcs_center_j = 3; *ubcs_center_i = user->info.mx - 1;
            break;
        case BC_FACE_NEG_Y:
            *center_k = 3; *center_j = 0; *center_i = 3;
            *off_k = 2; *off_j = 0; *off_i = 2;
            *wall_k = 1; *wall_j = 0; *wall_i = 1;
            *ubcs_center_k = 3; *ubcs_center_j = 0; *ubcs_center_i = 3;
            break;
        case BC_FACE_POS_Y:
            *center_k = 3; *center_j = user->info.my - 2; *center_i = 3;
            *off_k = 2; *off_j = user->info.my - 2; *off_i = 2;
            *wall_k = 1; *wall_j = user->info.my - 2; *wall_i = 1;
            *ubcs_center_k = 3; *ubcs_center_j = user->info.my - 1; *ubcs_center_i = 3;
            break;
        case BC_FACE_NEG_Z:
            *center_k = 0; *center_j = 3; *center_i = 3;
            *off_k = 0; *off_j = 2; *off_i = 2;
            *wall_k = 0; *wall_j = 1; *wall_i = 1;
            *ubcs_center_k = 0; *ubcs_center_j = 3; *ubcs_center_i = 3;
            break;
        case BC_FACE_POS_Z:
            *center_k = user->info.mz - 2; *center_j = 3; *center_i = 3;
            *off_k = user->info.mz - 2; *off_j = 2; *off_i = 2;
            *wall_k = user->info.mz - 2; *wall_j = 1; *wall_i = 1;
            *ubcs_center_k = user->info.mz - 1; *ubcs_center_j = 3; *ubcs_center_i = 3;
            break;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Resets one boundary-face configuration entry to a neutral test-local baseline.
 */
static void ResetBoundaryFaceConfig(BoundaryFaceConfig *cfg)
{
    cfg->mathematical_type = INTERFACE;
    cfg->handler_type = BC_HANDLER_UNDEFINED;
    cfg->face_id = BC_FACE_NEG_X;
    if (cfg->params) {
        FreeBC_ParamList(cfg->params);
        cfg->params = NULL;
    }
}
/**
 * @brief Tests that face-service detection matches a defined inlet face.
 */

static PetscErrorCode TestCanRankServiceFaceMatchesInletWhenDefined(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscBool generic_face = PETSC_FALSE;
    PetscBool inlet_face = PETSC_FALSE;
    PetscBool any_face_serviceable = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    for (BCFace face = BC_FACE_NEG_X; face <= BC_FACE_POS_Z; face++) {
        PetscCall(CanRankServiceFace(&user->info, user->IM, user->JM, user->KM, face, &generic_face));
        if (generic_face) {
            any_face_serviceable = PETSC_TRUE;
        }

        user->inletFaceDefined = PETSC_TRUE;
        user->identifiedInletBCFace = face;
        PetscCall(CanRankServiceInletFace(user, &user->info, user->IM, user->JM, user->KM, &inlet_face));
        PetscCall(PicurvAssertBool((PetscBool)(inlet_face == generic_face), "inlet face service check should match generic face check"));
    }
    PetscCall(PicurvAssertBool(any_face_serviceable, "at least one global face should be serviceable in a non-degenerate domain"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that inlet-face service requires an inlet face to be defined.
 */

static PetscErrorCode TestCanRankServiceInletFaceRequiresDefinition(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscBool can_service = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    user->inletFaceDefined = PETSC_FALSE;
    for (BCFace face = BC_FACE_NEG_X; face <= BC_FACE_POS_Z; face++) {
        user->identifiedInletBCFace = face;
        PetscCall(CanRankServiceInletFace(user, &user->info, user->IM, user->JM, user->KM, &can_service));
        PetscCall(PicurvAssertBool((PetscBool)!can_service, "undefined inlet face should not be serviceable"));
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the boundary-condition factory assignments for representative handlers.
 */

static PetscErrorCode TestBoundaryConditionFactoryAssignments(void)
{
    BoundaryCondition *wall = NULL;
    BoundaryCondition *inlet = NULL;

    PetscFunctionBeginUser;
    PetscCall(BoundaryCondition_Create(BC_HANDLER_WALL_NOSLIP, &wall));
    PetscCall(PicurvAssertBool((PetscBool)(wall != NULL), "factory should return wall handler"));
    PetscCall(PicurvAssertIntEqual(BC_PRIORITY_WALL, wall->priority, "wall handler priority"));
    PetscCall(PicurvAssertBool((PetscBool)(wall->Apply != NULL), "wall handler should expose Apply"));
    PetscCall(PicurvAssertBool((PetscBool)(wall->Destroy == NULL), "wall handler should not require destroy hook"));
    PetscCall(DestroyBoundaryHandler(&wall));

    PetscCall(BoundaryCondition_Create(BC_HANDLER_INLET_CONSTANT_VELOCITY, &inlet));
    PetscCall(PicurvAssertBool((PetscBool)(inlet != NULL), "factory should return inlet handler"));
    PetscCall(PicurvAssertIntEqual(BC_PRIORITY_INLET, inlet->priority, "inlet handler priority"));
    PetscCall(PicurvAssertBool((PetscBool)(inlet->Initialize != NULL), "inlet handler should expose Initialize"));
    PetscCall(PicurvAssertBool((PetscBool)(inlet->Apply != NULL), "inlet handler should expose Apply"));
    PetscCall(PicurvAssertBool((PetscBool)(inlet->Destroy != NULL), "inlet handler should expose its destroy hook"));
    PetscCall(DestroyBoundaryHandler(&inlet));

    PetscFunctionReturn(0);
}
/**
 * @brief Tests the implemented-handler matrix exposed by the factory.
 */

static PetscErrorCode TestBoundaryConditionFactoryImplementedHandlerMatrix(void)
{
    struct HandlerExpectation {
        BCHandlerType handler;
        BCPriorityType priority;
        PetscBool expect_apply;
        PetscBool expect_initialize;
    };
    const struct HandlerExpectation expectations[] = {
        {BC_HANDLER_WALL_NOSLIP, BC_PRIORITY_WALL, PETSC_TRUE, PETSC_FALSE},
        {BC_HANDLER_INLET_CONSTANT_VELOCITY, BC_PRIORITY_INLET, PETSC_TRUE, PETSC_TRUE},
        {BC_HANDLER_INLET_PARABOLIC, BC_PRIORITY_INLET, PETSC_TRUE, PETSC_TRUE},
        {BC_HANDLER_OUTLET_CONSERVATION, BC_PRIORITY_OUTLET, PETSC_TRUE, PETSC_FALSE},
    };

    PetscFunctionBeginUser;
    for (size_t i = 0; i < sizeof(expectations) / sizeof(expectations[0]); ++i) {
        BoundaryCondition *bc = NULL;
        PetscCall(BoundaryCondition_Create(expectations[i].handler, &bc));
        PetscCall(PicurvAssertBool((PetscBool)(bc != NULL), "factory should allocate a handler object"));
        PetscCall(PicurvAssertIntEqual(expectations[i].priority, bc->priority, "handler priority should match expectation"));
        PetscCall(PicurvAssertBool((PetscBool)((bc->Apply != NULL) == expectations[i].expect_apply), "Apply hook expectation mismatch"));
        PetscCall(PicurvAssertBool((PetscBool)((bc->Initialize != NULL) == expectations[i].expect_initialize), "Initialize hook expectation mismatch"));
        PetscCall(DestroyBoundaryHandler(&bc));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that unsupported handler types are rejected by the factory.
 */

static PetscErrorCode TestBoundaryConditionFactoryRejectsUnsupportedHandler(void)
{
    BoundaryCondition *bc = NULL;
    PetscErrorCode ierr_create;

    PetscFunctionBeginUser;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_create = BoundaryCondition_Create(BC_HANDLER_OUTLET_PRESSURE, &bc);
    PetscCall(PetscPopErrorHandler());

    PetscCall(PicurvAssertBool((PetscBool)(ierr_create != 0), "unsupported handler should return a non-zero error code"));
    if (bc) {
        PetscCall(PetscFree(bc));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Tests deterministic inlet-face grid-location helpers across all faces.
 */

static PetscErrorCode TestGetDeterministicFaceGridLocationFaceMatrix(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 8, 8));
    simCtx->np = 128;

    for (BCFace face = BC_FACE_NEG_X; face <= BC_FACE_POS_Z; ++face) {
        PetscInt ci = -1, cj = -1, ck = -1;
        PetscReal xi = -1.0, eta = -1.0, zta = -1.0;
        PetscBool placed = PETSC_FALSE;

        user->identifiedInletBCFace = face;
        PetscCall(GetDeterministicFaceGridLocation(
            user, &user->info,
            user->info.xs, user->info.ys, user->info.zs,
            user->info.mx, user->info.my, user->info.mz,
            0,
            &ci, &cj, &ck, &xi, &eta, &zta, &placed));

        PetscCall(PicurvAssertBool(placed, "single-rank deterministic face placement should succeed"));
        PetscCall(PicurvAssertBool((PetscBool)(ci >= user->info.xs && ci < user->info.xs + user->info.xm), "ci must map to owned node window"));
        PetscCall(PicurvAssertBool((PetscBool)(cj >= user->info.ys && cj < user->info.ys + user->info.ym), "cj must map to owned node window"));
        PetscCall(PicurvAssertBool((PetscBool)(ck >= user->info.zs && ck < user->info.zs + user->info.zm), "ck must map to owned node window"));
        PetscCall(PicurvAssertBool((PetscBool)(xi >= 0.0 && xi < 1.0), "xi should be in [0,1)"));
        PetscCall(PicurvAssertBool((PetscBool)(eta >= 0.0 && eta < 1.0), "eta should be in [0,1)"));
        PetscCall(PicurvAssertBool((PetscBool)(zta >= 0.0 && zta < 1.0), "zta should be in [0,1)"));

        if (face == BC_FACE_NEG_X || face == BC_FACE_POS_X) {
            PetscCall(PicurvAssertRealNear(0.5, xi, 1.0e-10, "deterministic x-face placement should sit halfway into boundary-adjacent cell"));
        } else if (face == BC_FACE_NEG_Y || face == BC_FACE_POS_Y) {
            PetscCall(PicurvAssertRealNear(0.5, eta, 1.0e-10, "deterministic y-face placement should sit halfway into boundary-adjacent cell"));
        } else {
            PetscCall(PicurvAssertRealNear(0.5, zta, 1.0e-10, "deterministic z-face placement should sit halfway into boundary-adjacent cell"));
        }
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests random inlet-face cell selection across all supported faces.
 */

static PetscErrorCode TestGetRandomCellAndLogicalCoordsOnInletFaceMatrix(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscRandom rand_i = NULL, rand_j = NULL, rand_k = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 8, 8));

    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_i));
    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_j));
    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rand_k));
    PetscCall(PetscRandomSetInterval(rand_i, 0.0, 1.0));
    PetscCall(PetscRandomSetInterval(rand_j, 0.0, 1.0));
    PetscCall(PetscRandomSetInterval(rand_k, 0.0, 1.0));
    PetscCall(PetscRandomSetFromOptions(rand_i));
    PetscCall(PetscRandomSetFromOptions(rand_j));
    PetscCall(PetscRandomSetFromOptions(rand_k));

    for (BCFace face = BC_FACE_NEG_X; face <= BC_FACE_POS_Z; ++face) {
        PetscInt ci = -1, cj = -1, ck = -1;
        PetscReal xi = -1.0, eta = -1.0, zta = -1.0;
        const PetscReal boundary_eps = 5.0e-4;

        user->identifiedInletBCFace = face;
        PetscCall(GetRandomCellAndLogicalCoordsOnInletFace(
            user, &user->info,
            user->info.xs, user->info.ys, user->info.zs,
            user->info.mx, user->info.my, user->info.mz,
            &rand_i, &rand_j, &rand_k,
            &ci, &cj, &ck,
            &xi, &eta, &zta));

        PetscCall(PicurvAssertBool((PetscBool)(ci >= user->info.xs && ci < user->info.xs + user->info.xm), "ci must map to owned node window"));
        PetscCall(PicurvAssertBool((PetscBool)(cj >= user->info.ys && cj < user->info.ys + user->info.ym), "cj must map to owned node window"));
        PetscCall(PicurvAssertBool((PetscBool)(ck >= user->info.zs && ck < user->info.zs + user->info.zm), "ck must map to owned node window"));
        PetscCall(PicurvAssertBool((PetscBool)(xi >= 0.0 && xi <= 1.0), "xi should be in [0,1]"));
        PetscCall(PicurvAssertBool((PetscBool)(eta >= 0.0 && eta <= 1.0), "eta should be in [0,1]"));
        PetscCall(PicurvAssertBool((PetscBool)(zta >= 0.0 && zta <= 1.0), "zta should be in [0,1]"));

        switch (face) {
            case BC_FACE_NEG_X:
                PetscCall(PicurvAssertBool((PetscBool)(xi <= boundary_eps), "NEG_X should pin xi near 0"));
                break;
            case BC_FACE_POS_X:
                PetscCall(PicurvAssertBool((PetscBool)(xi >= 1.0 - boundary_eps), "POS_X should pin xi near 1"));
                break;
            case BC_FACE_NEG_Y:
                PetscCall(PicurvAssertBool((PetscBool)(eta <= boundary_eps), "NEG_Y should pin eta near 0"));
                break;
            case BC_FACE_POS_Y:
                PetscCall(PicurvAssertBool((PetscBool)(eta >= 1.0 - boundary_eps), "POS_Y should pin eta near 1"));
                break;
            case BC_FACE_NEG_Z:
                PetscCall(PicurvAssertBool((PetscBool)(zta <= boundary_eps), "NEG_Z should pin zta near 0"));
                break;
            case BC_FACE_POS_Z:
                PetscCall(PicurvAssertBool((PetscBool)(zta >= 1.0 - boundary_eps), "POS_Z should pin zta near 1"));
                break;
        }
    }

    PetscCall(PetscRandomDestroy(&rand_i));
    PetscCall(PetscRandomDestroy(&rand_j));
    PetscCall(PetscRandomDestroy(&rand_k));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests constant inlet handler initialization and face writes on a Z inlet.
 */

static PetscErrorCode TestInletConstantVelocityHandlerBehavior(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundaryCondition *bc = NULL;
    BCContext ctx;
    Cmpnts ***ucont = NULL;
    Cmpnts ***ubcs = NULL;
    PetscReal local_inflow = 0.0;
    PetscReal local_outflow = 0.0;
    PetscReal expected_flux = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Bcs.Ubcs, 0.0));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    user->boundary_faces[BC_FACE_NEG_Z].face_id = BC_FACE_NEG_Z;
    user->boundary_faces[BC_FACE_NEG_Z].mathematical_type = INLET;
    user->boundary_faces[BC_FACE_NEG_Z].handler_type = BC_HANDLER_INLET_CONSTANT_VELOCITY;
    PetscCall(AppendBCParam(&user->boundary_faces[BC_FACE_NEG_Z].params, "vz", "2.5"));

    ctx.user = user;
    ctx.face_id = BC_FACE_NEG_Z;
    PetscCall(BoundaryCondition_Create(BC_HANDLER_INLET_CONSTANT_VELOCITY, &bc));
    PetscCall(bc->Initialize(bc, &ctx));
    PetscCall(bc->PostStep(bc, &ctx, &local_inflow, &local_outflow));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(PicurvAssertRealNear(2.5, ucont[0][3][3].z, 1.0e-12, "constant inlet should set Ucont normal component"));
    PetscCall(PicurvAssertRealNear(2.5, ubcs[0][3][3].z, 1.0e-12, "constant inlet should set Ubcs normal component"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));

    expected_flux = (PetscReal)(user->info.mx - 2) * (PetscReal)(user->info.my - 2) * 2.5;
    PetscCall(PicurvAssertRealNear(expected_flux, local_inflow, 1.0e-12, "constant inlet PostStep should sum the face flux"));
    PetscCall(PicurvAssertRealNear(0.0, local_outflow, 1.0e-12, "constant inlet should not contribute to outflow accumulation"));

    PetscCall(DestroyBoundaryHandler(&bc));
    FreeBC_ParamList(user->boundary_faces[BC_FACE_NEG_Z].params);
    user->boundary_faces[BC_FACE_NEG_Z].params = NULL;
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests wall no-slip application across the full face matrix.
 */
static PetscErrorCode TestWallNoSlipHandlerFaceMatrix(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundaryCondition *bc = NULL;
    BCContext ctx;
    Cmpnts ***ucont = NULL;
    Cmpnts ***ubcs = NULL;
    const BCFace faces[] = {
        BC_FACE_NEG_X, BC_FACE_POS_X,
        BC_FACE_NEG_Y, BC_FACE_POS_Y,
        BC_FACE_NEG_Z, BC_FACE_POS_Z
    };

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(BoundaryCondition_Create(BC_HANDLER_WALL_NOSLIP, &bc));

    for (size_t n = 0; n < sizeof(faces) / sizeof(faces[0]); ++n) {
        const BCFace face = faces[n];
        PetscInt ucont_k, ucont_j, ucont_i, ubcs_k, ubcs_j, ubcs_i;

        PetscCall(VecSet(user->Ucont, 5.0));
        PetscCall(VecSet(user->Bcs.Ubcs, 7.0));
        ctx.user = user;
        ctx.face_id = face;
        PetscCall(bc->Apply(bc, &ctx));

        PetscCall(GetRepresentativeFaceSlots(user, face, &ucont_k, &ucont_j, &ucont_i, &ubcs_k, &ubcs_j, &ubcs_i));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertRealNear(0.0, GetFaceNormalComponent(ucont[ucont_k][ucont_j][ucont_i], face), 1.0e-12, "wall no-slip should zero the normal flux component"));
        PetscCall(PicurvAssertRealNear(0.0, ubcs[ubcs_k][ubcs_j][ubcs_i].x, 1.0e-12, "wall no-slip should zero Ubcs.x"));
        PetscCall(PicurvAssertRealNear(0.0, ubcs[ubcs_k][ubcs_j][ubcs_i].y, 1.0e-12, "wall no-slip should zero Ubcs.y"));
        PetscCall(PicurvAssertRealNear(0.0, ubcs[ubcs_k][ubcs_j][ubcs_i].z, 1.0e-12, "wall no-slip should zero Ubcs.z"));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
    }

    PetscCall(DestroyBoundaryHandler(&bc));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests constant-inlet initialization, flux accounting, and face writes across all faces.
 */
static PetscErrorCode TestInletConstantVelocityHandlerFaceMatrix(void)
{
    struct ConstantInletCase {
        BCFace face;
        PetscReal velocity;
        const char *value_text;
    };
    const struct ConstantInletCase cases[] = {
        {BC_FACE_NEG_X, 1.5, "1.5"},
        {BC_FACE_POS_X, 2.0, "2.0"},
        {BC_FACE_NEG_Y, 2.5, "2.5"},
        {BC_FACE_POS_Y, 3.0, "3.0"},
        {BC_FACE_NEG_Z, 3.5, "3.5"},
        {BC_FACE_POS_Z, 4.0, "4.0"},
    };
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BCContext ctx;
    Cmpnts ***ucont = NULL;
    Cmpnts ***ubcs = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    for (size_t n = 0; n < sizeof(cases) / sizeof(cases[0]); ++n) {
        const struct ConstantInletCase *test_case = &cases[n];
        BoundaryCondition *bc = NULL;
        PetscReal local_inflow = 0.0;
        PetscReal local_outflow = 0.0;
        PetscInt ucont_k, ucont_j, ucont_i, ubcs_k, ubcs_j, ubcs_i;
        const PetscReal expected_value = GetFaceOrientationSign(test_case->face) * test_case->velocity;
        const PetscReal expected_flux = GetFaceInteriorPointCount(user, test_case->face) * expected_value;

        PetscCall(VecSet(user->Ucont, 0.0));
        PetscCall(VecSet(user->Bcs.Ubcs, 0.0));
        ResetBoundaryFaceConfig(&user->boundary_faces[test_case->face]);
        user->boundary_faces[test_case->face].face_id = test_case->face;
        user->boundary_faces[test_case->face].mathematical_type = INLET;
        user->boundary_faces[test_case->face].handler_type = BC_HANDLER_INLET_CONSTANT_VELOCITY;
        PetscCall(AppendBCParam(&user->boundary_faces[test_case->face].params,
                                GetInletParamKey(test_case->face),
                                test_case->value_text));

        ctx.user = user;
        ctx.face_id = test_case->face;
        PetscCall(BoundaryCondition_Create(BC_HANDLER_INLET_CONSTANT_VELOCITY, &bc));
        PetscCall(bc->PreStep(bc, &ctx, &local_inflow, &local_outflow));
        PetscCall(PicurvAssertRealNear(0.0, local_inflow, 1.0e-12, "constant inlet PreStep should leave inflow unchanged"));
        PetscCall(PicurvAssertRealNear(0.0, local_outflow, 1.0e-12, "constant inlet PreStep should leave outflow unchanged"));
        PetscCall(bc->Initialize(bc, &ctx));
        PetscCall(bc->PostStep(bc, &ctx, &local_inflow, &local_outflow));

        PetscCall(GetRepresentativeFaceSlots(user, test_case->face, &ucont_k, &ucont_j, &ucont_i, &ubcs_k, &ubcs_j, &ubcs_i));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertRealNear(expected_value, GetFaceNormalComponent(ucont[ucont_k][ucont_j][ucont_i], test_case->face), 1.0e-12, "constant inlet should set the face-normal Ucont component"));
        PetscCall(PicurvAssertRealNear(expected_value, GetFaceNormalComponent(ubcs[ubcs_k][ubcs_j][ubcs_i], test_case->face), 1.0e-12, "constant inlet should set the face-normal Ubcs component"));
        PetscCall(PicurvAssertRealNear((test_case->face == BC_FACE_NEG_X || test_case->face == BC_FACE_POS_X) ? expected_value : 0.0, ubcs[ubcs_k][ubcs_j][ubcs_i].x, 1.0e-12, "constant inlet should only set the expected Ubcs axis"));
        PetscCall(PicurvAssertRealNear((test_case->face == BC_FACE_NEG_Y || test_case->face == BC_FACE_POS_Y) ? expected_value : 0.0, ubcs[ubcs_k][ubcs_j][ubcs_i].y, 1.0e-12, "constant inlet should only set the expected Ubcs axis"));
        PetscCall(PicurvAssertRealNear((test_case->face == BC_FACE_NEG_Z || test_case->face == BC_FACE_POS_Z) ? expected_value : 0.0, ubcs[ubcs_k][ubcs_j][ubcs_i].z, 1.0e-12, "constant inlet should only set the expected Ubcs axis"));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertRealNear(expected_flux, local_inflow, 1.0e-12, "constant inlet PostStep should integrate the face flux with orientation"));
        PetscCall(PicurvAssertRealNear(0.0, local_outflow, 1.0e-12, "constant inlet should not add to outflow"));

        PetscCall(DestroyBoundaryHandler(&bc));
        ResetBoundaryFaceConfig(&user->boundary_faces[test_case->face]);
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests parabolic inlet handler shape on a tiny Z-face.
 */

static PetscErrorCode TestInletParabolicProfileHandlerBehavior(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundaryCondition *bc = NULL;
    BCContext ctx;
    Cmpnts ***ucont = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    user->boundary_faces[BC_FACE_NEG_Z].face_id = BC_FACE_NEG_Z;
    user->boundary_faces[BC_FACE_NEG_Z].mathematical_type = INLET;
    user->boundary_faces[BC_FACE_NEG_Z].handler_type = BC_HANDLER_INLET_PARABOLIC;
    PetscCall(AppendBCParam(&user->boundary_faces[BC_FACE_NEG_Z].params, "v_max", "4.0"));

    ctx.user = user;
    ctx.face_id = BC_FACE_NEG_Z;
    PetscCall(BoundaryCondition_Create(BC_HANDLER_INLET_PARABOLIC, &bc));
    PetscCall(bc->Initialize(bc, &ctx));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear(4.0, ucont[0][3][3].z, 1.0e-12, "parabolic inlet should peak at the face centerline"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[0][1][1].z, 1.0e-12, "parabolic inlet should vanish at wall-adjacent locations"));
    PetscCall(PicurvAssertBool((PetscBool)(ucont[0][3][3].z > ucont[0][2][2].z), "parabolic inlet centerline should exceed off-center velocity"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(DestroyBoundaryHandler(&bc));
    FreeBC_ParamList(user->boundary_faces[BC_FACE_NEG_Z].params);
    user->boundary_faces[BC_FACE_NEG_Z].params = NULL;
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests parabolic-inlet centerline, wall, and flux behavior across all faces.
 */
static PetscErrorCode TestInletParabolicProfileHandlerFaceMatrix(void)
{
    const BCFace faces[] = {
        BC_FACE_NEG_X, BC_FACE_POS_X,
        BC_FACE_NEG_Y, BC_FACE_POS_Y,
        BC_FACE_NEG_Z, BC_FACE_POS_Z
    };
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BCContext ctx;
    Cmpnts ***ucont = NULL;
    Cmpnts ***ubcs = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    for (size_t n = 0; n < sizeof(faces) / sizeof(faces[0]); ++n) {
        const BCFace face = faces[n];
        BoundaryCondition *bc = NULL;
        PetscReal local_inflow = 0.0;
        PetscReal local_outflow = 0.0;
        PetscInt center_k, center_j, center_i;
        PetscInt off_k, off_j, off_i;
        PetscInt wall_k, wall_j, wall_i;
        PetscInt ubcs_center_k, ubcs_center_j, ubcs_center_i;
        const PetscReal expected_center = GetFaceOrientationSign(face) * 4.0;

        PetscCall(VecSet(user->Ucont, 0.0));
        PetscCall(VecSet(user->Bcs.Ubcs, 0.0));
        ResetBoundaryFaceConfig(&user->boundary_faces[face]);
        user->boundary_faces[face].face_id = face;
        user->boundary_faces[face].mathematical_type = INLET;
        user->boundary_faces[face].handler_type = BC_HANDLER_INLET_PARABOLIC;
        PetscCall(AppendBCParam(&user->boundary_faces[face].params, "v_max", "4.0"));

        ctx.user = user;
        ctx.face_id = face;
        PetscCall(BoundaryCondition_Create(BC_HANDLER_INLET_PARABOLIC, &bc));
        PetscCall(bc->PreStep(bc, &ctx, &local_inflow, &local_outflow));
        PetscCall(PicurvAssertRealNear(0.0, local_inflow, 1.0e-12, "parabolic inlet PreStep should leave inflow unchanged"));
        PetscCall(PicurvAssertRealNear(0.0, local_outflow, 1.0e-12, "parabolic inlet PreStep should leave outflow unchanged"));
        PetscCall(bc->Initialize(bc, &ctx));
        PetscCall(bc->PostStep(bc, &ctx, &local_inflow, &local_outflow));

        PetscCall(GetParabolicSampleSlots(user, face,
                                          &center_k, &center_j, &center_i,
                                          &off_k, &off_j, &off_i,
                                          &wall_k, &wall_j, &wall_i,
                                          &ubcs_center_k, &ubcs_center_j, &ubcs_center_i));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertRealNear(expected_center, GetFaceNormalComponent(ucont[center_k][center_j][center_i], face), 1.0e-12, "parabolic inlet should peak at the face centerline"));
        PetscCall(PicurvAssertRealNear(expected_center, GetFaceNormalComponent(ubcs[ubcs_center_k][ubcs_center_j][ubcs_center_i], face), 1.0e-12, "parabolic inlet should write the centerline boundary velocity"));
        PetscCall(PicurvAssertRealNear(0.0, GetFaceNormalComponent(ucont[wall_k][wall_j][wall_i], face), 1.0e-12, "parabolic inlet should vanish at the wall"));
        PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(GetFaceNormalComponent(ucont[center_k][center_j][center_i], face)) >
                                               PetscAbsReal(GetFaceNormalComponent(ucont[off_k][off_j][off_i], face))),
                                   "parabolic inlet centerline magnitude should exceed the off-center magnitude"));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(local_inflow) > 0.0), "parabolic inlet PostStep should accumulate a non-zero flux"));
        PetscCall(PicurvAssertBool((PetscBool)((local_inflow > 0.0) == (GetFaceOrientationSign(face) > 0.0)),
                                   "parabolic inlet PostStep should preserve the face orientation"));
        PetscCall(PicurvAssertRealNear(0.0, local_outflow, 1.0e-12, "parabolic inlet should not add to outflow"));

        PetscCall(DestroyBoundaryHandler(&bc));
        ResetBoundaryFaceConfig(&user->boundary_faces[face]);
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests outlet conservation handler correction and post-step flux accounting.
 */

static PetscErrorCode TestOutletConservationHandlerBehavior(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundaryCondition *bc = NULL;
    BCContext ctx;
    PetscReal global_inflow = 30.0;
    PetscReal global_farfield_in = 0.0;
    PetscReal global_farfield_out = 0.0;
    PetscReal global_outflow = 25.0;
    PetscReal local_inflow = 0.0;
    PetscReal local_outflow = 0.0;
    Cmpnts ***ucont = NULL;
    Cmpnts ***ucat = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Ucat, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    simCtx->AreaOutSum = (PetscReal)(user->info.mx - 2) * (PetscReal)(user->info.my - 2);

    user->boundary_faces[BC_FACE_POS_Z].face_id = BC_FACE_POS_Z;
    user->boundary_faces[BC_FACE_POS_Z].mathematical_type = OUTLET;
    user->boundary_faces[BC_FACE_POS_Z].handler_type = BC_HANDLER_OUTLET_CONSERVATION;

    ctx.user = user;
    ctx.face_id = BC_FACE_POS_Z;
    ctx.global_inflow_sum = &global_inflow;
    ctx.global_farfield_inflow_sum = &global_farfield_in;
    ctx.global_farfield_outflow_sum = &global_farfield_out;
    ctx.global_outflow_sum = &global_outflow;

    PetscCall(BoundaryCondition_Create(BC_HANDLER_OUTLET_CONSERVATION, &bc));
    PetscCall(bc->Apply(bc, &ctx));
    PetscCall(bc->PostStep(bc, &ctx, &local_inflow, &local_outflow));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcat, &ucat));
    PetscCall(PicurvAssertRealNear(1.2, ucont[5][3][3].z, 1.0e-12, "outlet correction should add the expected flux trim"));
    PetscCall(PicurvAssertRealNear(1.0, ucat[5][3][3].z, 1.0e-12, "outlet handler should preserve the interior Ucat reference"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcat, &ucat));

    PetscCall(PicurvAssertRealNear(30.0, local_outflow, 1.0e-12, "outlet PostStep should report corrected flux"));
    PetscCall(PicurvAssertRealNear(0.0, local_inflow, 1.0e-12, "outlet handler should not accumulate inflow"));

    PetscCall(DestroyBoundaryHandler(&bc));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests outlet-conservation correction and flux accounting across all outlet faces.
 */
static PetscErrorCode TestOutletConservationHandlerFaceMatrix(void)
{
    const BCFace faces[] = {
        BC_FACE_NEG_X, BC_FACE_POS_X,
        BC_FACE_NEG_Y, BC_FACE_POS_Y,
        BC_FACE_NEG_Z, BC_FACE_POS_Z
    };
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BCContext ctx;
    Cmpnts ***ucont = NULL;
    Cmpnts ***ubcs = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Ucat, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));

    for (size_t n = 0; n < sizeof(faces) / sizeof(faces[0]); ++n) {
        const BCFace face = faces[n];
        BoundaryCondition *bc = NULL;
        PetscReal global_inflow = 0.0;
        PetscReal global_farfield_in = 0.0;
        PetscReal global_farfield_out = 0.0;
        PetscReal global_outflow = 0.0;
        PetscReal local_inflow = 0.0;
        PetscReal local_outflow = 0.0;
        PetscReal area = GetFaceInteriorPointCount(user, face);
        PetscInt ucont_k, ucont_j, ucont_i, ubcs_k, ubcs_j, ubcs_i;

        PetscCall(VecSet(user->Ucont, 0.0));
        PetscCall(VecSet(user->Bcs.Ubcs, 0.0));
        ResetBoundaryFaceConfig(&user->boundary_faces[face]);
        user->boundary_faces[face].face_id = face;
        user->boundary_faces[face].mathematical_type = OUTLET;
        user->boundary_faces[face].handler_type = BC_HANDLER_OUTLET_CONSERVATION;
        user->simCtx->AreaOutSum = area;

        ctx.user = user;
        ctx.face_id = face;
        ctx.global_inflow_sum = &global_inflow;
        ctx.global_farfield_inflow_sum = &global_farfield_in;
        ctx.global_farfield_outflow_sum = &global_farfield_out;
        ctx.global_outflow_sum = &global_outflow;

        PetscCall(BoundaryCondition_Create(BC_HANDLER_OUTLET_CONSERVATION, &bc));
        PetscCall(bc->PreStep(bc, &ctx, &local_inflow, &local_outflow));
        PetscCall(PicurvAssertRealNear(area, local_outflow, 1.0e-12, "outlet PreStep should measure the uncorrected face flux"));

        global_inflow = area + 0.2 * area;
        global_outflow = local_outflow;
        local_outflow = 0.0;
        PetscCall(bc->Apply(bc, &ctx));
        PetscCall(bc->PostStep(bc, &ctx, &local_inflow, &local_outflow));

        PetscCall(GetRepresentativeFaceSlots(user, face, &ucont_k, &ucont_j, &ucont_i, &ubcs_k, &ubcs_j, &ubcs_i));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertRealNear(1.2, GetFaceNormalComponent(ucont[ucont_k][ucont_j][ucont_i], face), 1.0e-12, "outlet correction should add the configured flux trim"));
        PetscCall(PicurvAssertRealNear(1.0, ubcs[ubcs_k][ubcs_j][ubcs_i].x, 1.0e-12, "outlet handler should copy Ucat into Ubcs.x"));
        PetscCall(PicurvAssertRealNear(1.0, ubcs[ubcs_k][ubcs_j][ubcs_i].y, 1.0e-12, "outlet handler should copy Ucat into Ubcs.y"));
        PetscCall(PicurvAssertRealNear(1.0, ubcs[ubcs_k][ubcs_j][ubcs_i].z, 1.0e-12, "outlet handler should copy Ucat into Ubcs.z"));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(PicurvAssertRealNear(0.0, local_inflow, 1.0e-12, "outlet handler should not add to inflow"));
        PetscCall(PicurvAssertRealNear(1.2 * area, local_outflow, 1.0e-12, "outlet PostStep should report the corrected flux"));

        PetscCall(DestroyBoundaryHandler(&bc));
        ResetBoundaryFaceConfig(&user->boundary_faces[face]);
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-boundaries PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"can-rank-service-face-matches-inlet-when-defined", TestCanRankServiceFaceMatchesInletWhenDefined},
        {"can-rank-service-inlet-face-requires-definition", TestCanRankServiceInletFaceRequiresDefinition},
        {"boundary-condition-factory-assignments", TestBoundaryConditionFactoryAssignments},
        {"boundary-condition-factory-implemented-handler-matrix", TestBoundaryConditionFactoryImplementedHandlerMatrix},
        {"boundary-condition-factory-rejects-unsupported-handler", TestBoundaryConditionFactoryRejectsUnsupportedHandler},
        {"deterministic-face-grid-location-matrix", TestGetDeterministicFaceGridLocationFaceMatrix},
        {"random-inlet-face-location-matrix", TestGetRandomCellAndLogicalCoordsOnInletFaceMatrix},
        {"wall-no-slip-handler-face-matrix", TestWallNoSlipHandlerFaceMatrix},
        {"inlet-constant-velocity-handler-behavior", TestInletConstantVelocityHandlerBehavior},
        {"inlet-constant-velocity-handler-face-matrix", TestInletConstantVelocityHandlerFaceMatrix},
        {"inlet-parabolic-profile-handler-behavior", TestInletParabolicProfileHandlerBehavior},
        {"inlet-parabolic-profile-handler-face-matrix", TestInletParabolicProfileHandlerFaceMatrix},
        {"outlet-conservation-handler-behavior", TestOutletConservationHandlerBehavior},
        {"outlet-conservation-handler-face-matrix", TestOutletConservationHandlerFaceMatrix},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv boundary tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-boundaries", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
