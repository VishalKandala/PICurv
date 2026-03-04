#include "test_support.h"

#include <stdlib.h>

#ifdef PETSC_ERR_ENVIRONMENT
#define PICURV_TEST_ERR_ENVIRONMENT PETSC_ERR_ENVIRONMENT
#else
#define PICURV_TEST_ERR_ENVIRONMENT PETSC_ERR_USER_INPUT
#endif

static PetscErrorCode TestEnvironmentVisible(void)
{
    const char *petsc_dir = getenv("PETSC_DIR");

    PetscFunctionBeginUser;
    PetscCheck(petsc_dir && petsc_dir[0] != '\0', PETSC_COMM_SELF, PICURV_TEST_ERR_ENVIRONMENT,
               "PETSC_DIR is not set in the environment.");
    PetscFunctionReturn(0);
}

static PetscErrorCode TestBasicPetscObjects(void)
{
    DM da = NULL, swarm = NULL;
    Vec vec = NULL;

    PetscFunctionBeginUser;
    PetscCall(DMDACreate3d(PETSC_COMM_WORLD,
                           DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DMDA_STENCIL_BOX,
                           3, 3, 3,
                           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                           1, 1,
                           NULL, NULL, NULL,
                           &da));
    PetscCall(DMSetUp(da));
    PetscCall(DMCreateGlobalVector(da, &vec));

    PetscCall(DMCreate(PETSC_COMM_WORLD, &swarm));
    PetscCall(DMSetType(swarm, DMSWARM));
    PetscCall(DMSetDimension(swarm, 3));
    PetscCall(DMSwarmSetType(swarm, DMSWARM_BASIC));
    PetscCall(DMSwarmSetCellDM(swarm, da));
    PetscCall(DMSwarmRegisterPetscDatatypeField(swarm, "position", 3, PETSC_REAL));
    PetscCall(DMSwarmFinalizeFieldRegister(swarm));
    PetscCall(DMSwarmSetLocalSizes(swarm, 1, 0));

    PetscCall(VecDestroy(&vec));
    PetscCall(DMDestroy(&swarm));
    PetscCall(DMDestroy(&da));
    PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"environment-visible", TestEnvironmentVisible},
        {"basic-petsc-objects", TestBasicPetscObjects},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv installation validation");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("doctor", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
