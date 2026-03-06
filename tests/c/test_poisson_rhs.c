#include "test_support.h"

#include "poisson.h"
#include "rhs.h"

static PetscErrorCode EnsurePoissonAndRhsVectors(UserCtx *user)
{
    PetscFunctionBeginUser;
    if (!user->Phi) PetscCall(DMCreateGlobalVector(user->da, &user->Phi));
    if (!user->lPhi) PetscCall(DMCreateLocalVector(user->da, &user->lPhi));
    if (!user->Aj) PetscCall(DMCreateGlobalVector(user->da, &user->Aj));
    if (!user->lAj) PetscCall(DMCreateLocalVector(user->da, &user->lAj));
    if (!user->Diffusivity) PetscCall(DMCreateGlobalVector(user->da, &user->Diffusivity));
    if (!user->lDiffusivity) PetscCall(DMCreateLocalVector(user->da, &user->lDiffusivity));
    PetscFunctionReturn(0);
}

static PetscErrorCode TestUpdatePressureAddsPhi(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(EnsurePoissonAndRhsVectors(user));
    PetscCall(VecSet(user->P, 1.0));
    PetscCall(VecSet(user->Phi, 0.25));

    PetscCall(UpdatePressure(user));
    PetscCall(PicurvAssertVecConstant(user->P, 1.25, 1.0e-12, "UpdatePressure should add Phi into P"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

static PetscErrorCode TestPoissonRHSZeroDivergence(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec B = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(EnsurePoissonAndRhsVectors(user));

    simCtx->dt = 0.5;
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));

    PetscCall(VecDuplicate(user->P, &B));
    PetscCall(PoissonRHS(user, B));
    PetscCall(PicurvAssertVecConstant(B, 0.0, 1.0e-12, "zero velocity divergence should produce zero Poisson RHS"));
    PetscCall(PicurvAssertRealNear(0.0, simCtx->summationRHS, 1.0e-12, "global Poisson RHS sum should be zero"));

    PetscCall(VecDestroy(&B));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

static PetscErrorCode TestComputeBodyForcesDispatcher(void)
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

    PetscCall(ComputeBodyForces(user, rct));
    PetscCall(DMDAVecGetArrayRead(user->fda, rct, &rct_arr));
    PetscCall(PicurvAssertRealNear(1.5, rct_arr[1][1][1].x, 1.0e-12, "ComputeBodyForces should apply driven-channel source"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].y, 1.0e-12, "ComputeBodyForces should keep y unchanged"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].z, 1.0e-12, "ComputeBodyForces should keep z unchanged"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, rct, &rct_arr));

    PetscCall(VecDestroy(&rct));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

static PetscErrorCode TestComputeEulerianDiffusivityMolecularOnly(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***diff_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(EnsurePoissonAndRhsVectors(user));

    simCtx->ren = 2.0;
    simCtx->schmidt_number = 4.0;
    simCtx->les = PETSC_FALSE;
    simCtx->rans = PETSC_FALSE;
    PetscCall(VecZeroEntries(user->Diffusivity));

    PetscCall(ComputeEulerianDiffusivity(user));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diff_arr));
    PetscCall(PicurvAssertRealNear(0.125, diff_arr[1][1][1], 1.0e-12, "molecular diffusivity interior value"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diff_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"update-pressure-adds-phi", TestUpdatePressureAddsPhi},
        {"poisson-rhs-zero-divergence", TestPoissonRHSZeroDivergence},
        {"compute-body-forces-dispatcher", TestComputeBodyForcesDispatcher},
        {"compute-eulerian-diffusivity-molecular-only", TestComputeEulerianDiffusivityMolecularOnly},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv Poisson/RHS tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-poisson-rhs", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
