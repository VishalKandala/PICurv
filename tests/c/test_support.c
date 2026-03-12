/**
 * @file test_support.c
 * @brief Shared C test fixtures, assertions, and PETSc helper utilities.
 */

#include "test_support.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
/**
 * @brief Destroys a PETSc vector only when the handle is non-null.
 */

static PetscErrorCode DestroyVecIfSet(Vec *vec)
{
    PetscFunctionBeginUser;
    if (vec && *vec) {
        PetscCall(VecDestroy(vec));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Destroys a PETSc DM only when the handle is non-null.
 */

static PetscErrorCode DestroyDMIfSet(DM *dm)
{
    PetscFunctionBeginUser;
    if (dm && *dm) {
        PetscCall(DMDestroy(dm));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Registers one DMSwarm field used by the C test fixtures.
 */

static PetscErrorCode RegisterSwarmFieldForTests(DM swarm, const char *field_name, PetscInt field_dim, PetscDataType dtype)
{
    PetscFunctionBeginUser;
    PetscCall(DMSwarmRegisterPetscDatatypeField(swarm, field_name, field_dim, dtype));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs a named C test suite and prints pass/fail progress markers.
 */

PetscErrorCode PicurvRunTests(const char *suite_name, const PicurvTestCase *cases, size_t case_count)
{
    PetscFunctionBeginUser;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==> Running %s (%zu tests)\n", suite_name, case_count));
    for (size_t i = 0; i < case_count; ++i) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  -> %s\n", cases[i].name));
        PetscCall(cases[i].fn());
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "     [PASS] %s\n", cases[i].name));
    }

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==> %s complete\n", suite_name));
    PetscFunctionReturn(0);
}
/**
 * @brief Ensures a directory exists for test output.
 */

PetscErrorCode PicurvEnsureDir(const char *path)
{
    PetscFunctionBeginUser;
    if (mkdir(path, 0777) != 0 && errno != EEXIST) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create directory '%s': %s", path, strerror(errno));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Creates a unique temporary directory for one test case.
 */

PetscErrorCode PicurvMakeTempDir(char *path, size_t path_len)
{
    PetscFunctionBeginUser;
    if (!path || path_len < 24) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Temp directory buffer is missing or too small.");
    }

    PetscCall(PetscSNPrintf(path, path_len, "/tmp/picurv-test-XXXXXX"));
    if (!mkdtemp(path)) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "mkdtemp failed for '%s': %s", path, strerror(errno));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Builds minimal SimCtx and UserCtx fixtures for C unit tests.
 */

PetscErrorCode PicurvCreateMinimalContexts(SimCtx **simCtx_out, UserCtx **user_out, PetscInt mx, PetscInt my, PetscInt mz)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    if (!simCtx_out || !user_out) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointers cannot be NULL.");
    }

    PetscCall(PetscCalloc1(1, &simCtx));
    PetscCall(PetscCalloc1(1, &user));

    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &simCtx->rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &simCtx->size));
    simCtx->block_number = 1;
    simCtx->dt = 0.1;
    simCtx->tiout = 2;
    simCtx->forceScalingFactor = 1.0;
    simCtx->LoggingFrequency = 1;
    simCtx->exec_mode = EXEC_MODE_SOLVER;
    simCtx->mom_solver_type = MOMENTUM_SOLVER_EXPLICIT_RK;
    simCtx->poisson = 0;
    simCtx->ren = 1.0;
    simCtx->schmidt_number = 1.0;
    simCtx->StartStep = 0;
    simCtx->StepsToRun = 1;
    simCtx->step = 1;
    simCtx->np = 0;
    PetscCall(PetscStrncpy(simCtx->euler_subdir, "euler", sizeof(simCtx->euler_subdir)));
    PetscCall(PetscStrncpy(simCtx->particle_subdir, "particles", sizeof(simCtx->particle_subdir)));
    PetscCall(PetscStrncpy(simCtx->output_dir, "/tmp", sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, "/tmp", sizeof(simCtx->restart_dir)));
    simCtx->mglevels = 1;
    simCtx->usermg.mglevels = 1;
    PetscCall(PetscCalloc1(1, &simCtx->usermg.mgctx));
    simCtx->usermg.mgctx[0].thislevel = 0;
    PetscCall(PetscCalloc1(simCtx->size, &simCtx->bboxlist));

    user->simCtx = simCtx;
    user->_this = 0;
    user->thislevel = 0;
    user->mglevels = 1;
    user->IM = mx;
    user->JM = my;
    user->KM = mz;
    simCtx->usermg.mgctx[0].user = user;
    for (PetscMPIInt rank_idx = 0; rank_idx < simCtx->size; ++rank_idx) {
        simCtx->bboxlist[rank_idx].min_coords.x = 0.0;
        simCtx->bboxlist[rank_idx].min_coords.y = 0.0;
        simCtx->bboxlist[rank_idx].min_coords.z = 0.0;
        simCtx->bboxlist[rank_idx].max_coords.x = (PetscReal)(mx - 1);
        simCtx->bboxlist[rank_idx].max_coords.y = (PetscReal)(my - 1);
        simCtx->bboxlist[rank_idx].max_coords.z = (PetscReal)(mz - 1);
    }

    PetscCall(DMDACreate3d(PETSC_COMM_WORLD,
                           DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DMDA_STENCIL_BOX,
                           mx, my, mz,
                           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                           1, 1,
                           NULL, NULL, NULL,
                           &user->da));
    PetscCall(DMSetUp(user->da));
    PetscCall(DMDASetUniformCoordinates(user->da, 0.0, (PetscReal)(mx - 1), 0.0, (PetscReal)(my - 1), 0.0, (PetscReal)(mz - 1)));

    PetscCall(DMDACreate3d(PETSC_COMM_WORLD,
                           DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           DMDA_STENCIL_BOX,
                           mx, my, mz,
                           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                           3, 1,
                           NULL, NULL, NULL,
                           &user->fda));
    PetscCall(DMSetUp(user->fda));
    PetscCall(DMDASetUniformCoordinates(user->fda, 0.0, (PetscReal)(mx - 1), 0.0, (PetscReal)(my - 1), 0.0, (PetscReal)(mz - 1)));
    PetscCall(DMDAGetLocalInfo(user->fda, &user->info));

    PetscCall(DMCreateGlobalVector(user->da, &user->P));
    PetscCall(DMCreateLocalVector(user->da, &user->lP));
    PetscCall(DMCreateGlobalVector(user->da, &user->Nvert));
    PetscCall(DMCreateLocalVector(user->da, &user->lNvert));
    PetscCall(DMCreateGlobalVector(user->da, &user->ParticleCount));
    PetscCall(DMCreateLocalVector(user->da, &user->lParticleCount));
    PetscCall(DMCreateGlobalVector(user->da, &user->Psi));
    PetscCall(DMCreateLocalVector(user->da, &user->lPsi));
    PetscCall(DMCreateGlobalVector(user->da, &user->Qcrit));
    PetscCall(DMCreateGlobalVector(user->da, &user->P_nodal));
    PetscCall(DMCreateGlobalVector(user->da, &user->Psi_nodal));

    PetscCall(DMCreateGlobalVector(user->fda, &user->Ucat));
    PetscCall(DMCreateLocalVector(user->fda, &user->lUcat));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Ucont));
    PetscCall(DMCreateLocalVector(user->fda, &user->lUcont));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Csi));
    PetscCall(DMCreateLocalVector(user->fda, &user->lCsi));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Eta));
    PetscCall(DMCreateLocalVector(user->fda, &user->lEta));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Zet));
    PetscCall(DMCreateLocalVector(user->fda, &user->lZet));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Cent));
    PetscCall(DMCreateLocalVector(user->fda, &user->lCent));
    PetscCall(DMCreateGlobalVector(user->fda, &user->Ucat_nodal));

    PetscCall(VecZeroEntries(user->P));
    PetscCall(VecZeroEntries(user->Nvert));
    PetscCall(VecZeroEntries(user->ParticleCount));
    PetscCall(VecZeroEntries(user->Psi));
    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(VecZeroEntries(user->Ucont));
    PetscCall(VecZeroEntries(user->Cent));

    PetscCall(PicurvPopulateIdentityMetrics(user));

    *simCtx_out = simCtx;
    *user_out = user;
    PetscFunctionReturn(0);
}
/**
 * @brief Populates identity metric vectors on the minimal grid fixture.
 */

PetscErrorCode PicurvPopulateIdentityMetrics(UserCtx *user)
{
    Cmpnts ***csi = NULL;
    Cmpnts ***eta = NULL;
    Cmpnts ***zet = NULL;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecGetArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecGetArray(user->fda, user->Zet, &zet));

    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                csi[k][j][i].x = 1.0; csi[k][j][i].y = 0.0; csi[k][j][i].z = 0.0;
                eta[k][j][i].x = 0.0; eta[k][j][i].y = 1.0; eta[k][j][i].z = 0.0;
                zet[k][j][i].x = 0.0; zet[k][j][i].y = 0.0; zet[k][j][i].z = 1.0;
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Zet, &zet));

    PetscCall(DMGlobalToLocalBegin(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Psi, INSERT_VALUES, user->lPsi));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Psi, INSERT_VALUES, user->lPsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Cent, INSERT_VALUES, user->lCent));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Cent, INSERT_VALUES, user->lCent));
    PetscFunctionReturn(0);
}
/**
 * @brief Creates matched solver and post-processing swarms for tests.
 */

PetscErrorCode PicurvCreateSwarmPair(UserCtx *user, PetscInt nlocal, const char *post_field_name)
{
    PetscFunctionBeginUser;
    if (!user) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    }

    PetscCall(DMCreate(PETSC_COMM_WORLD, &user->swarm));
    PetscCall(DMSetType(user->swarm, DMSWARM));
    PetscCall(DMSetDimension(user->swarm, 3));
    PetscCall(DMSwarmSetType(user->swarm, DMSWARM_BASIC));
    PetscCall(DMSwarmSetCellDM(user->swarm, user->da));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "position", 3, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "velocity", 3, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "weight", 3, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "Psi", 1, PETSC_REAL));
    PetscCall(DMSwarmFinalizeFieldRegister(user->swarm));
    PetscCall(DMSwarmSetLocalSizes(user->swarm, nlocal, 0));

    PetscCall(DMCreate(PETSC_COMM_WORLD, &user->post_swarm));
    PetscCall(DMSetType(user->post_swarm, DMSWARM));
    PetscCall(DMSetDimension(user->post_swarm, 3));
    PetscCall(DMSwarmSetType(user->post_swarm, DMSWARM_BASIC));
    PetscCall(DMSwarmSetCellDM(user->post_swarm, user->da));
    PetscCall(RegisterSwarmFieldForTests(user->post_swarm, post_field_name, 1, PETSC_REAL));
    PetscCall(DMSwarmFinalizeFieldRegister(user->post_swarm));
    PetscCall(DMSwarmSetLocalSizes(user->post_swarm, nlocal, 0));
    PetscFunctionReturn(0);
}
/**
 * @brief Destroys minimal SimCtx/UserCtx fixtures and all owned PETSc objects.
 */

PetscErrorCode PicurvDestroyMinimalContexts(SimCtx **simCtx_ptr, UserCtx **user_ptr)
{
    UserCtx *user = NULL;
    SimCtx *simCtx = NULL;

    PetscFunctionBeginUser;
    if (user_ptr) {
        user = *user_ptr;
    }
    if (simCtx_ptr) {
        simCtx = *simCtx_ptr;
    }

    if (user) {
        PetscCall(DestroyDMIfSet(&user->swarm));
        PetscCall(DestroyDMIfSet(&user->post_swarm));

        PetscCall(DestroyVecIfSet(&user->P));
        PetscCall(DestroyVecIfSet(&user->lP));
        PetscCall(DestroyVecIfSet(&user->Phi));
        PetscCall(DestroyVecIfSet(&user->lPhi));
        PetscCall(DestroyVecIfSet(&user->Nvert));
        PetscCall(DestroyVecIfSet(&user->lNvert));
        PetscCall(DestroyVecIfSet(&user->ParticleCount));
        PetscCall(DestroyVecIfSet(&user->lParticleCount));
        PetscCall(DestroyVecIfSet(&user->Psi));
        PetscCall(DestroyVecIfSet(&user->lPsi));
        PetscCall(DestroyVecIfSet(&user->Qcrit));
        PetscCall(DestroyVecIfSet(&user->P_nodal));
        PetscCall(DestroyVecIfSet(&user->Psi_nodal));
        PetscCall(DestroyVecIfSet(&user->Ucat));
        PetscCall(DestroyVecIfSet(&user->lUcat));
        PetscCall(DestroyVecIfSet(&user->Ucont));
        PetscCall(DestroyVecIfSet(&user->lUcont));
        PetscCall(DestroyVecIfSet(&user->Csi));
        PetscCall(DestroyVecIfSet(&user->lCsi));
        PetscCall(DestroyVecIfSet(&user->Eta));
        PetscCall(DestroyVecIfSet(&user->lEta));
        PetscCall(DestroyVecIfSet(&user->Zet));
        PetscCall(DestroyVecIfSet(&user->lZet));
        PetscCall(DestroyVecIfSet(&user->Aj));
        PetscCall(DestroyVecIfSet(&user->lAj));
        PetscCall(DestroyVecIfSet(&user->IAj));
        PetscCall(DestroyVecIfSet(&user->lIAj));
        PetscCall(DestroyVecIfSet(&user->JAj));
        PetscCall(DestroyVecIfSet(&user->lJAj));
        PetscCall(DestroyVecIfSet(&user->KAj));
        PetscCall(DestroyVecIfSet(&user->lKAj));
        PetscCall(DestroyVecIfSet(&user->Diffusivity));
        PetscCall(DestroyVecIfSet(&user->lDiffusivity));
        PetscCall(DestroyVecIfSet(&user->DiffusivityGradient));
        PetscCall(DestroyVecIfSet(&user->lDiffusivityGradient));
        PetscCall(DestroyVecIfSet(&user->Nu_t));
        PetscCall(DestroyVecIfSet(&user->lNu_t));
        PetscCall(DestroyVecIfSet(&user->Cent));
        PetscCall(DestroyVecIfSet(&user->lCent));
        PetscCall(DestroyVecIfSet(&user->Ucat_nodal));

        PetscCall(DestroyDMIfSet(&user->da));
        PetscCall(DestroyDMIfSet(&user->fda));
        PetscCall(PetscFree(user));
        if (user_ptr) {
            *user_ptr = NULL;
        }
    }

    if (simCtx) {
        PetscCall(PetscFree(simCtx->bboxlist));
        PetscCall(PetscFree(simCtx->usermg.mgctx));
        PetscCall(PetscFree(simCtx->pps));
        PetscCall(PetscFree(simCtx));
        if (simCtx_ptr) {
            *simCtx_ptr = NULL;
        }
    }

    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that two real values agree within tolerance.
 */

PetscErrorCode PicurvAssertRealNear(PetscReal expected, PetscReal actual, PetscReal tol, const char *context)
{
    PetscFunctionBeginUser;
    if (PetscAbsReal(expected - actual) > tol) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "[FAIL] %s | expected=%0.12e actual=%0.12e tol=%0.12e\n",
                              context, (double)expected, (double)actual, (double)tol));
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Assertion failed.");
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that two integer values are equal.
 */

PetscErrorCode PicurvAssertIntEqual(PetscInt expected, PetscInt actual, const char *context)
{
    PetscFunctionBeginUser;
    if (expected != actual) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "[FAIL] %s | expected=%" PetscInt_FMT " actual=%" PetscInt_FMT "\n",
                              context, expected, actual));
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Assertion failed.");
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that one boolean condition is true.
 */

PetscErrorCode PicurvAssertBool(PetscBool value, const char *context)
{
    PetscFunctionBeginUser;
    if (!value) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "[FAIL] %s\n", context));
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Assertion failed.");
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that a filesystem path exists as a readable file.
 */

PetscErrorCode PicurvAssertFileExists(const char *path, const char *context)
{
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscTestFile(path, 'r', &exists));
    if (!exists) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "[FAIL] %s | missing file: %s\n", context, path));
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Expected file is missing.");
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that a PETSc vector is spatially constant within tolerance.
 */

PetscErrorCode PicurvAssertVecConstant(Vec vec, PetscScalar expected, PetscReal tol, const char *context)
{
    PetscReal vmin = 0.0;
    PetscReal vmax = 0.0;

    PetscFunctionBeginUser;
    PetscCall(VecMin(vec, NULL, &vmin));
    PetscCall(VecMax(vec, NULL, &vmax));
    PetscCall(PicurvAssertRealNear((PetscReal)expected, vmin, tol, context));
    PetscCall(PicurvAssertRealNear((PetscReal)expected, vmax, tol, context));
    PetscFunctionReturn(0);
}
