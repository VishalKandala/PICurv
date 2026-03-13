/**
 * @file test_postprocessing.c
 * @brief C unit tests for post-processing kernel functions.
 */

#include "test_support.h"

#include "postprocessing_kernels.h"
/**
 * @brief Tests specific turbulent kinetic energy computation on particle data.
 */

static PetscErrorCode TestComputeSpecificKE(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal (*vel_arr)[3] = NULL;
    PetscScalar *ske_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    vel_arr[0][0] = 1.0;
    vel_arr[0][1] = 2.0;
    vel_arr[0][2] = 2.0;
    vel_arr[1][0] = 0.0;
    vel_arr[1][1] = 3.0;
    vel_arr[1][2] = 4.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(ComputeSpecificKE(user, "velocity", "ske"));
    PetscCall(DMSwarmGetField(user->post_swarm, "ske", NULL, NULL, (void *)&ske_arr));
    PetscCall(PicurvAssertRealNear(4.5, ske_arr[0], 1.0e-12, "ComputeSpecificKE first particle"));
    PetscCall(PicurvAssertRealNear(12.5, ske_arr[1], 1.0e-12, "ComputeSpecificKE second particle"));
    PetscCall(DMSwarmRestoreField(user->post_swarm, "ske", NULL, NULL, (void *)&ske_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle displacement computation against reference positions.
 */

static PetscErrorCode TestComputeDisplacement(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal (*pos_arr)[3] = NULL;
    PetscScalar *disp_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "disp"));

    simCtx->psrc_x = 1.0;
    simCtx->psrc_y = 2.0;
    simCtx->psrc_z = 3.0;

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 1.0; pos_arr[0][1] = 2.0; pos_arr[0][2] = 3.0;
    pos_arr[1][0] = 4.0; pos_arr[1][1] = 6.0; pos_arr[1][2] = 3.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(ComputeDisplacement(user, "disp"));
    PetscCall(DMSwarmGetField(user->post_swarm, "disp", NULL, NULL, (void *)&disp_arr));
    PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(disp_arr[0]), 1.0e-12, "ComputeDisplacement first particle"));
    PetscCall(PicurvAssertRealNear(5.0, PetscRealPart(disp_arr[1]), 1.0e-12, "ComputeDisplacement second particle"));
    PetscCall(DMSwarmRestoreField(user->post_swarm, "disp", NULL, NULL, (void *)&disp_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests nodal averaging of scalar cell-centered data.
 */

static PetscErrorCode TestComputeNodalAverageScalar(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    const PetscScalar ***p_nodal_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->P, 7.0));
    PetscCall(VecSet(user->P_nodal, -1.0));

    PetscCall(ComputeNodalAverage(user, "P", "P_nodal"));

    PetscCall(DMDAVecGetArrayRead(user->da, user->P_nodal, (void *)&p_nodal_arr));
    PetscCall(PicurvAssertRealNear(7.0, PetscRealPart(p_nodal_arr[0][0][0]), 1.0e-12, "ComputeNodalAverage interior node"));
    PetscCall(PicurvAssertRealNear(-1.0, PetscRealPart(p_nodal_arr[user->KM][user->JM][user->IM]), 1.0e-12,
                                   "ComputeNodalAverage untouched non-physical boundary node"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P_nodal, (void *)&p_nodal_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests normalization of one field relative to a reference field.
 */

static PetscErrorCode TestNormalizeRelativeField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt ref_idx = 0;
    PetscScalar ref_value = 0.0;
    PetscReal vmin = 0.0, vmax = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscCalloc1(1, &simCtx->pps));
    simCtx->pps->reference[0] = 1;
    simCtx->pps->reference[1] = 1;
    simCtx->pps->reference[2] = 1;

    PetscCall(VecSet(user->P, 10.0));
    ref_idx = simCtx->pps->reference[2] * (user->IM * user->JM) +
              simCtx->pps->reference[1] * user->IM +
              simCtx->pps->reference[0];
    PetscCall(VecSetValue(user->P, ref_idx, 4.0, INSERT_VALUES));
    PetscCall(VecAssemblyBegin(user->P));
    PetscCall(VecAssemblyEnd(user->P));

    PetscCall(NormalizeRelativeField(user, "P"));

    PetscCall(VecGetValues(user->P, 1, &ref_idx, &ref_value));
    PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(ref_value), 1.0e-12, "NormalizeRelativeField reference value"));
    PetscCall(VecMin(user->P, NULL, &vmin));
    PetscCall(VecMax(user->P, NULL, &vmax));
    PetscCall(PicurvAssertRealNear(0.0, vmin, 1.0e-12, "NormalizeRelativeField minimum"));
    PetscCall(PicurvAssertRealNear(6.0, vmax, 1.0e-12, "NormalizeRelativeField maximum"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests dimensionalization of pressure data from nondimensional values.
 */

static PetscErrorCode TestDimensionalizePressureField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->scaling.P_ref = 3.0;

    PetscCall(VecSet(user->P, 2.0));
    PetscCall(DimensionalizeField(user, "P"));
    PetscCall(PicurvAssertVecConstant(user->P, 6.0, 1.0e-12, "DimensionalizeField should scale pressure by P_ref"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests Q-criterion computation for a quiescent velocity field.
 */

static PetscErrorCode TestComputeQCriterionZeroFlow(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(VecSet(user->Aj, 1.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));

    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));

    PetscCall(ComputeQCriterion(user));
    PetscCall(PicurvAssertVecConstant(user->Qcrit, 0.0, 1.0e-12, "ComputeQCriterion should be zero for uniform flow"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-post PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"compute-specific-ke", TestComputeSpecificKE},
        {"compute-displacement", TestComputeDisplacement},
        {"compute-nodal-average-scalar", TestComputeNodalAverageScalar},
        {"normalize-relative-field", TestNormalizeRelativeField},
        {"dimensionalize-pressure-field", TestDimensionalizePressureField},
        {"compute-qcriterion-zero-flow", TestComputeQCriterionZeroFlow},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv post-processing tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-post", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
