#include "test_support.h"

#include "Boundaries.h"

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

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"can-rank-service-face-matches-inlet-when-defined", TestCanRankServiceFaceMatchesInletWhenDefined},
        {"can-rank-service-inlet-face-requires-definition", TestCanRankServiceInletFaceRequiresDefinition},
        {"boundary-condition-factory-assignments", TestBoundaryConditionFactoryAssignments},
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
