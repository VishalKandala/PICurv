/**
 * @file test_particle_kernels.c
 * @brief C test module for PICurv.
 */

#include "test_support.h"

#include "walkingsearch.h"
/**
 * @brief Test-local routine.
 */

static PetscErrorCode TestCheckCellWithinLocalGrid(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscBool within = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(CheckCellWithinLocalGrid(user, 1, 1, 1, &within));
    PetscCall(PicurvAssertBool(within, "cell (1,1,1) should be within the serial local grid"));

    PetscCall(CheckCellWithinLocalGrid(user, 3, 1, 1, &within));
    PetscCall(PicurvAssertBool((PetscBool)!within, "cell (3,1,1) should fall outside the valid local cell range"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Test-local routine.
 */

static PetscErrorCode TestInitializeTraversalParameters(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Particle particle;
    PetscInt idx = -1, idy = -1, idz = -1, steps = -1;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMemzero(&particle, sizeof(particle)));

    particle.PID = 42;
    particle.cell[0] = 1;
    particle.cell[1] = 1;
    particle.cell[2] = 1;

    PetscCall(InitializeTraversalParameters(user, &particle, &idx, &idy, &idz, &steps));
    PetscCall(PicurvAssertIntEqual(1, idx, "InitializeTraversalParameters should preserve prior i"));
    PetscCall(PicurvAssertIntEqual(1, idy, "InitializeTraversalParameters should preserve prior j"));
    PetscCall(PicurvAssertIntEqual(1, idz, "InitializeTraversalParameters should preserve prior k"));
    PetscCall(PicurvAssertIntEqual(0, steps, "InitializeTraversalParameters should reset traversal steps"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Test-local routine.
 */

static PetscErrorCode TestRetrieveCurrentCell(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cell cell;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMemzero(&cell, sizeof(cell)));

    PetscCall(RetrieveCurrentCell(user, 1, 1, 1, &cell));
    PetscCall(PicurvAssertRealNear(1.0, cell.vertices[0].x, 1.0e-12, "vertex 0 x coordinate"));
    PetscCall(PicurvAssertRealNear(1.0, cell.vertices[0].y, 1.0e-12, "vertex 0 y coordinate"));
    PetscCall(PicurvAssertRealNear(1.0, cell.vertices[0].z, 1.0e-12, "vertex 0 z coordinate"));
    PetscCall(PicurvAssertRealNear(2.0, cell.vertices[5].x, 1.0e-12, "vertex 5 x coordinate"));
    PetscCall(PicurvAssertRealNear(2.0, cell.vertices[5].y, 1.0e-12, "vertex 5 y coordinate"));
    PetscCall(PicurvAssertRealNear(2.0, cell.vertices[5].z, 1.0e-12, "vertex 5 z coordinate"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Entry point for this unit-test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"check-cell-within-local-grid", TestCheckCellWithinLocalGrid},
        {"initialize-traversal-parameters", TestInitializeTraversalParameters},
        {"retrieve-current-cell", TestRetrieveCurrentCell},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv particle kernel tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-particles", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
