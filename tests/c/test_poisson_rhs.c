/**
 * @file test_poisson_rhs.c
 * @brief C unit tests for Poisson, RHS, body-force, and diffusivity helpers.
 */

#include "test_support.h"

#include "poisson.h"
#include "rhs.h"
#include "verification_sources.h"
/**
 * @brief Allocates Poisson/RHS support vectors required by the tests.
 */

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
/**
 * @brief Tests that pressure updates add the correction potential.
 */

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
/**
 * @brief Tests that the Poisson RHS is zero for zero-divergence velocity.
 */

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
/**
 * @brief Tests the body-force dispatcher across supported source terms.
 */

static PetscErrorCode TestComputeBodyForcesDispatcher(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec rct = NULL;
    Cmpnts ***rct_arr = NULL;
    Cmpnts ***l_csi = NULL;
    Cmpnts ***l_eta = NULL;
    Cmpnts ***l_zet = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecDuplicate(user->Ucont, &rct));
    PetscCall(VecZeroEntries(rct));

    user->boundary_faces[BC_FACE_NEG_X].face_id = BC_FACE_NEG_X;
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    user->boundary_faces[BC_FACE_POS_X].face_id = BC_FACE_POS_X;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    simCtx->bulkVelocityCorrection = 2.0;
    simCtx->dt = 1.0;
    simCtx->forceScalingFactor = 1.0;
    simCtx->drivingForceMagnitude = 0.0;
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMDAVecGetArray(user->fda, user->lCsi, &l_csi));
    PetscCall(DMDAVecGetArray(user->fda, user->lEta, &l_eta));
    PetscCall(DMDAVecGetArray(user->fda, user->lZet, &l_zet));
    l_csi[1][1][1].x = 1.0; l_csi[1][1][1].y = 0.0; l_csi[1][1][1].z = 0.0;
    l_eta[1][1][1].x = 0.0; l_eta[1][1][1].y = 1.0; l_eta[1][1][1].z = 0.0;
    l_zet[1][1][1].x = 0.0; l_zet[1][1][1].y = 0.0; l_zet[1][1][1].z = 1.0;
    PetscCall(DMDAVecRestoreArray(user->fda, user->lZet, &l_zet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lEta, &l_eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->lCsi, &l_csi));

    PetscCall(ComputeBodyForces(user, rct));
    PetscCall(DMDAVecGetArrayRead(user->fda, rct, &rct_arr));
    PetscCall(PicurvAssertRealNear(1.5, simCtx->drivingForceMagnitude, 1.0e-12,
                                   "ComputeBodyForces should update the driven-flow controller magnitude"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].y, 1.0e-12, "ComputeBodyForces should keep y unchanged"));
    PetscCall(PicurvAssertRealNear(0.0, rct_arr[1][1][1].z, 1.0e-12, "ComputeBodyForces should keep z unchanged"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, rct, &rct_arr));

    PetscCall(VecDestroy(&rct));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests Eulerian diffusivity for the molecular-only configuration.
 */

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
/**
 * @brief Tests that convection vanishes for a quiescent field.
 */

static PetscErrorCode TestConvectionZeroField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec conv = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecDuplicate(user->lUcont, &conv));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));

    PetscCall(Convection(user, user->lUcont, user->lUcat, conv));
    PetscCall(PicurvAssertVecConstant(conv, 0.0, 1.0e-12, "Convection should vanish for a quiescent field"));

    PetscCall(VecDestroy(&conv));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that viscous terms vanish for a uniform field.
 */

static PetscErrorCode TestViscousUniformField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec visc = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecDuplicate(user->lUcont, &visc));
    PetscCall(VecSet(user->Ucont, 1.0));
    PetscCall(VecSet(user->Ucat, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));

    PetscCall(Viscous(user, user->lUcont, user->lUcat, visc));
    PetscCall(PicurvAssertVecConstant(visc, 0.0, 1.0e-12, "Viscous should vanish for a uniform field"));

    PetscCall(VecDestroy(&visc));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the full RHS remains zero without forcing on a quiescent field.
 */

static PetscErrorCode TestComputeRHSZeroFieldNoForcing(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(ComputeRHS(user, user->Rhs));
    PetscCall(PicurvAssertVecConstant(user->Rhs, 0.0, 1.0e-12, "ComputeRHS should remain zero without forcing on a quiescent field"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests diffusivity-gradient computation on a constant field.
 */

static PetscErrorCode TestComputeEulerianDiffusivityGradientConstantField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Diffusivity, 0.75));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Diffusivity, INSERT_VALUES, user->lDiffusivity));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Diffusivity, INSERT_VALUES, user->lDiffusivity));

    PetscCall(ComputeEulerianDiffusivityGradient(user));
    PetscCall(PicurvAssertVecConstant(user->DiffusivityGradient, 0.0, 1.0e-12, "constant diffusivity should yield zero diffusivity gradient"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests verification-driven linear diffusivity override and its gradient.
 */

static PetscErrorCode TestComputeEulerianDiffusivityVerificationLinearX(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***diff_arr = NULL;
    Cmpnts ***grad_arr = NULL;
    Cmpnts ***cent = NULL;
    const PetscReal gamma0 = 0.5;
    const PetscReal slope_x = 0.25;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "analytical", sizeof(simCtx->eulerianSource)));
    simCtx->verificationDiffusivity.enabled = PETSC_TRUE;
    PetscCall(PetscStrncpy(simCtx->verificationDiffusivity.mode, "analytical", sizeof(simCtx->verificationDiffusivity.mode)));
    PetscCall(PetscStrncpy(simCtx->verificationDiffusivity.profile, "LINEAR_X", sizeof(simCtx->verificationDiffusivity.profile)));
    simCtx->verificationDiffusivity.gamma0 = gamma0;
    simCtx->verificationDiffusivity.slope_x = slope_x;

    PetscCall(PicurvAssertBool(VerificationDiffusivityOverrideActive(simCtx),
                               "verification diffusivity override should report as active"));
    PetscCall(ComputeEulerianDiffusivity(user));
    PetscCall(ComputeEulerianDiffusivityGradient(user));

    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diff_arr));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Cent, &cent));
    PetscCall(PicurvAssertRealNear(gamma0 + slope_x * cent[2][2][2].x, diff_arr[2][2][2], 1.0e-12,
                                   "verification override should populate the linear diffusivity field"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diff_arr));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->DiffusivityGradient, &grad_arr));
    PetscCall(PicurvAssertRealNear(slope_x / user->IM, grad_arr[2][2][2].x, 1.0e-12,
                                   "linear verification diffusivity should yield the current finite-difference x gradient"));
    PetscCall(PicurvAssertRealNear(0.0, grad_arr[2][2][2].y, 1.0e-12,
                                   "linear verification diffusivity should keep y gradient zero"));
    PetscCall(PicurvAssertRealNear(0.0, grad_arr[2][2][2].z, 1.0e-12,
                                   "linear verification diffusivity should keep z gradient zero"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->DiffusivityGradient, &grad_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the Poisson null-space operator removes the mean from a constant field.
 */

static PetscErrorCode TestPoissonNullSpaceFunctionRemovesMean(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec x = NULL;
    PetscReal sum = 0.0;
    PetscReal ***x_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecDuplicate(user->P, &x));
    PetscCall(VecSet(x, 3.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(PoissonNullSpaceFunction(NULL, x, user));
    PetscCall(VecSum(x, &sum));
    PetscCall(PicurvAssertRealNear(0.0, sum, 1.0e-10, "PoissonNullSpaceFunction should remove the global mean"));
    PetscCall(DMDAVecGetArrayRead(user->da, x, &x_arr));
    PetscCall(PicurvAssertRealNear(0.0, x_arr[1][1][1], 1.0e-10, "PoissonNullSpaceFunction should zero a uniform interior field"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, x, &x_arr));

    PetscCall(VecDestroy(&x));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that Poisson matrix assembly produces a populated operator on a tiny Cartesian grid.
 */
static PetscErrorCode TestPoissonLHSNewAssemblesOperator(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt rows = 0, cols = 0;
    PetscInt ncols = 0;
    const PetscInt *col_idx = NULL;
    const PetscScalar *values = NULL;
    PetscInt interior_row = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(PoissonLHSNew(user));
    PetscCall(PicurvAssertBool((PetscBool)(user->A != NULL), "PoissonLHSNew should allocate the Poisson operator"));

    PetscCall(MatGetSize(user->A, &rows, &cols));
    PetscCall(PicurvAssertIntEqual(user->info.mx * user->info.my * user->info.mz, rows, "Poisson operator row count should match the DA node count"));
    PetscCall(PicurvAssertIntEqual(rows, cols, "Poisson operator should be square"));

    interior_row = (2 * user->info.my + 2) * user->info.mx + 2;
    PetscCall(MatGetRow(user->A, interior_row, &ncols, &col_idx, &values));
    PetscCall(PicurvAssertBool((PetscBool)(ncols > 1), "Interior Poisson row should couple to neighboring nodes"));
    PetscCall(PicurvAssertBool((PetscBool)(values != NULL), "Interior Poisson row should expose non-null coefficients"));
    PetscCall(MatRestoreRow(user->A, interior_row, &ncols, &col_idx, &values));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that projection leaves a zero pressure-correction field unchanged.
 */

static PetscErrorCode TestProjectionZeroPhiLeavesVelocityUnchanged(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Phi, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(Projection(user));
    PetscCall(PicurvAssertVecConstant(user->Ucont, 0.0, 1.0e-12, "Projection should leave a zero-velocity field unchanged when Phi is zero"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that projection applies the expected x-direction correction for a linear pressure field.
 */
static PetscErrorCode TestProjectionLinearPhiCorrectsVelocity(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***phi = NULL;
    Cmpnts ***ucont = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    simCtx->dt = COEF_TIME_ACCURACY;
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMDAVecGetArray(user->da, user->Phi, &phi));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                phi[k][j][i] = (PetscReal)i;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->Phi, &phi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Phi, INSERT_VALUES, user->lPhi));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(Projection(user));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertRealNear(-1.0, ucont[2][2][2].x, 1.0e-10, "Projection should subtract the x pressure gradient under identity metrics"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[2][2][2].y, 1.0e-10, "Projection should leave the y component unchanged for an x-only gradient"));
    PetscCall(PicurvAssertRealNear(0.0, ucont[2][2][2].z, 1.0e-10, "Projection should leave the z component unchanged for an x-only gradient"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-poisson-rhs PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"update-pressure-adds-phi", TestUpdatePressureAddsPhi},
        {"poisson-rhs-zero-divergence", TestPoissonRHSZeroDivergence},
        {"compute-eulerian-diffusivity-molecular-only", TestComputeEulerianDiffusivityMolecularOnly},
        {"convection-zero-field", TestConvectionZeroField},
        {"viscous-uniform-field", TestViscousUniformField},
        {"compute-rhs-zero-field-no-forcing", TestComputeRHSZeroFieldNoForcing},
        {"compute-eulerian-diffusivity-gradient-constant-field", TestComputeEulerianDiffusivityGradientConstantField},
        {"compute-eulerian-diffusivity-verification-linear-x", TestComputeEulerianDiffusivityVerificationLinearX},
        {"poisson-null-space-function-removes-mean", TestPoissonNullSpaceFunctionRemovesMean},
        {"poisson-lhs-new-assembles-operator", TestPoissonLHSNewAssemblesOperator},
        {"projection-zero-phi-leaves-velocity-unchanged", TestProjectionZeroPhiLeavesVelocityUnchanged},
        {"projection-linear-phi-corrects-velocity", TestProjectionLinearPhiCorrectsVelocity},
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
