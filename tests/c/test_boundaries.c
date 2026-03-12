/**
 * @file test_boundaries.c
 * @brief C unit tests for boundary factories, inlet ownership, and face helpers.
 */

#include "test_support.h"

#include "Boundaries.h"
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
    BoundaryCondition *driven = NULL;

    PetscFunctionBeginUser;
    PetscCall(BoundaryCondition_Create(BC_HANDLER_WALL_NOSLIP, &wall));
    PetscCall(PicurvAssertBool((PetscBool)(wall != NULL), "factory should return wall handler"));
    PetscCall(PicurvAssertIntEqual(BC_PRIORITY_WALL, wall->priority, "wall handler priority"));
    PetscCall(PicurvAssertBool((PetscBool)(wall->Apply != NULL), "wall handler should expose Apply"));
    PetscCall(PicurvAssertBool((PetscBool)(wall->Destroy == NULL), "wall handler should not require destroy hook"));
    PetscCall(DestroyBoundaryHandler(&wall));

    PetscCall(BoundaryCondition_Create(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX, &driven));
    PetscCall(PicurvAssertBool((PetscBool)(driven != NULL), "factory should return driven periodic handler"));
    PetscCall(PicurvAssertIntEqual(BC_PRIORITY_INLET, driven->priority, "driven handler priority"));
    PetscCall(PicurvAssertBool((PetscBool)(driven->Initialize != NULL), "driven handler should expose Initialize"));
    PetscCall(PicurvAssertBool((PetscBool)(driven->PreStep != NULL), "driven handler should expose PreStep"));
    PetscCall(PicurvAssertBool((PetscBool)(driven->Apply != NULL), "driven handler should expose Apply"));
    PetscCall(PicurvAssertBool((PetscBool)(driven->Destroy != NULL), "driven handler should expose Destroy"));
    PetscCall(DestroyBoundaryHandler(&driven));

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
        {BC_HANDLER_PERIODIC_GEOMETRIC, BC_PRIORITY_WALL, PETSC_FALSE, PETSC_FALSE},
        {BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX, BC_PRIORITY_INLET, PETSC_TRUE, PETSC_TRUE},
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
