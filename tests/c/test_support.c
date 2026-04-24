/**
 * @file test_support.c
 * @brief Shared C test fixtures, assertions, and PETSc helper utilities.
 */

#include "test_support.h"

#include "grid.h"
#include "io.h"
#include "setup.h"

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
 * @brief Destroys a PETSc matrix only when the handle is non-null.
 */

static PetscErrorCode DestroyMatIfSet(Mat *mat)
{
    PetscFunctionBeginUser;
    if (mat && *mat) {
        PetscCall(MatDestroy(mat));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Destroys a PETSc KSP only when the handle is non-null.
 */

static PetscErrorCode DestroyKSPIfSet(KSP *ksp)
{
    PetscFunctionBeginUser;
    if (ksp && *ksp) {
        PetscCall(KSPDestroy(ksp));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Destroys a PETSc nullspace only when the handle is non-null.
 */

static PetscErrorCode DestroyNullSpaceIfSet(MatNullSpace *nullsp)
{
    PetscFunctionBeginUser;
    if (nullsp && *nullsp) {
        PetscCall(MatNullSpaceDestroy(nullsp));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Destroys a PETSc random generator only when the handle is non-null.
 */

static PetscErrorCode DestroyRandomIfSet(PetscRandom *rand_ctx)
{
    PetscFunctionBeginUser;
    if (rand_ctx && *rand_ctx) {
        PetscCall(PetscRandomDestroy(rand_ctx));
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
 * @brief Allocates and zeroes a global vector from the provided DM.
 */

static PetscErrorCode CreateZeroedGlobalVector(DM dm, Vec *vec)
{
    PetscFunctionBeginUser;
    PetscCall(DMCreateGlobalVector(dm, vec));
    PetscCall(VecSet(*vec, 0.0));
    PetscFunctionReturn(0);
}
/**
 * @brief Allocates and zeroes a local vector from the provided DM.
 */

static PetscErrorCode CreateZeroedLocalVector(DM dm, Vec *vec)
{
    PetscFunctionBeginUser;
    PetscCall(DMCreateLocalVector(dm, vec));
    PetscCall(VecSet(*vec, 0.0));
    PetscFunctionReturn(0);
}
/**
 * @brief Duplicates and zeroes a vector.
 */

static PetscErrorCode CreateZeroedDuplicate(Vec src, Vec *vec)
{
    PetscFunctionBeginUser;
    PetscCall(VecDuplicate(src, vec));
    PetscCall(VecSet(*vec, 0.0));
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
 * @brief Recursively removes a temporary directory created by PicurvMakeTempDir.
 */

PetscErrorCode PicurvRemoveTempDir(const char *path)
{
    char cmd[512];

    PetscFunctionBeginUser;
    if (!path || path[0] == '\0') PetscFunctionReturn(0);
    /* Safety: only remove paths under /tmp/picurv-test- */
    if (strncmp(path, "/tmp/picurv-test-", 17) != 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Refusing to remove path outside /tmp/picurv-test-*: '%s'", path);
    }
    PetscCall(PetscSNPrintf(cmd, sizeof(cmd), "rm -rf '%s'", path));
    if (system(cmd) != 0) {
        PetscCall(PetscPrintf(PETSC_COMM_SELF, "Warning: failed to remove temp dir '%s'\n", path));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Writes one small temporary text file used by the richer runtime fixtures.
 */

static PetscErrorCode WriteTextFileForTests(const char *path, const char *contents)
{
    FILE *file = NULL;

    PetscFunctionBeginUser;
    file = fopen(path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s' for writing.", path);
    fputs(contents, file);
    fclose(file);
    PetscFunctionReturn(0);
}
/**
 * @brief Creates a tiny control-file bundle used by richer runtime fixtures built through the setup path.
 */

static PetscErrorCode PrepareTinyRuntimeConfig(const char *bcs_contents,
                                               PetscBool enable_particles,
                                               char *tmpdir,
                                               size_t tmpdir_len,
                                               char *control_path,
                                               size_t control_path_len)
{
    char bcs_path[PETSC_MAX_PATH_LEN];
    char post_path[PETSC_MAX_PATH_LEN];
    char output_dir[PETSC_MAX_PATH_LEN];
    char log_dir[PETSC_MAX_PATH_LEN];
    char control_buffer[8192];
    const char *particle_block = NULL;
    const char *default_bcs =
        "-Xi WALL noslip\n"
        "+Xi WALL noslip\n"
        "-Eta WALL noslip\n"
        "+Eta WALL noslip\n"
        "-Zeta INLET constant_velocity vx=0.0 vy=0.0 vz=1.5\n"
        "+Zeta OUTLET conservation\n";

    PetscFunctionBeginUser;
    PetscCall(PicurvMakeTempDir(tmpdir, tmpdir_len));
    PetscCall(PetscSNPrintf(bcs_path, sizeof(bcs_path), "%s/bcs.run", tmpdir));
    PetscCall(PetscSNPrintf(post_path, sizeof(post_path), "%s/post.run", tmpdir));
    PetscCall(PetscSNPrintf(output_dir, sizeof(output_dir), "%s/results", tmpdir));
    PetscCall(PetscSNPrintf(log_dir, sizeof(log_dir), "%s/logs", tmpdir));
    PetscCall(PetscSNPrintf(control_path, control_path_len, "%s/test.control", tmpdir));

    PetscCall(WriteTextFileForTests(bcs_path, bcs_contents ? bcs_contents : default_bcs));
    PetscCall(WriteTextFileForTests(
        post_path,
        "startTime = 0\n"
        "endTime = 1\n"
        "timeStep = 1\n"
        "output_particles = false\n"));

    if (enable_particles) {
        particle_block =
            "-numParticles 8\n"
            "-pinit 2\n"
            "-psrc_x 0.5\n"
            "-psrc_y 0.5\n"
            "-psrc_z 0.5\n"
            "-particle_restart_mode init\n";
    } else {
        particle_block =
            "-numParticles 0\n"
            "-pinit 2\n";
    }

    PetscCall(PetscSNPrintf(
        control_buffer,
        sizeof(control_buffer),
        "-start_step 0\n"
        "-totalsteps 2\n"
        "-ren 100.0\n"
        "-dt 0.001\n"
        "-finit 1\n"
        "-ucont_x 0.0\n"
        "-ucont_y 0.0\n"
        "-ucont_z 1.5\n"
        "-bcs_files %s\n"
        "-profiling_timestep_mode off\n"
        "-profiling_final_summary true\n"
        "-postprocessing_config_file %s\n"
        "-grid\n"
        "-im 6\n"
        "-jm 6\n"
        "-km 6\n"
        "-xMins 0.0\n"
        "-xMaxs 1.0\n"
        "-yMins 0.0\n"
        "-yMaxs 1.0\n"
        "-zMins 0.0\n"
        "-zMaxs 1.0\n"
        "-rxs 1.0\n"
        "-rys 1.0\n"
        "-rzs 1.0\n"
        "-cgrids 0\n"
        "-nblk 1\n"
        "-euler_field_source solve\n"
        "-mom_solver_type EXPLICIT_RK\n"
        "-mg_level 1\n"
        "-poisson 0\n"
        "-tio 0\n"
        "-particle_console_output_freq 0\n"
        "-logfreq 1\n"
        "-output_dir %s\n"
        "-restart_dir %s\n"
        "-log_dir %s\n"
        "%s",
        bcs_path,
        post_path,
        output_dir,
        output_dir,
        log_dir,
        particle_block));
    PetscCall(WriteTextFileForTests(control_path, control_buffer));
    PetscFunctionReturn(0);
}
/**
 * @brief Builds minimal SimCtx and UserCtx fixtures for C unit tests with configurable periodicity.
 */
PetscErrorCode PicurvCreateMinimalContextsWithPeriodicity(SimCtx **simCtx_out,
                                                          UserCtx **user_out,
                                                          PetscInt mx,
                                                          PetscInt my,
                                                          PetscInt mz,
                                                          PetscBool x_periodic,
                                                          PetscBool y_periodic,
                                                          PetscBool z_periodic)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundingBox *boxes = NULL;
    PetscInt da_mx = mx + 1;
    PetscInt da_my = my + 1;
    PetscInt da_mz = mz + 1;
    DMBoundaryType x_boundary = x_periodic ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType y_boundary = y_periodic ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType z_boundary = z_periodic ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    PetscInt stencil_width = (x_periodic || y_periodic || z_periodic) ? 3 : 1;

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
    simCtx->i_periodic = x_periodic ? 1 : 0;
    simCtx->j_periodic = y_periodic ? 1 : 0;
    simCtx->k_periodic = z_periodic ? 1 : 0;
    PetscCall(PetscStrncpy(simCtx->euler_subdir, "euler", sizeof(simCtx->euler_subdir)));
    PetscCall(PetscStrncpy(simCtx->particle_subdir, "particles", sizeof(simCtx->particle_subdir)));
    PetscCall(PetscStrncpy(simCtx->output_dir, "/tmp", sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, "/tmp", sizeof(simCtx->restart_dir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, "/tmp", sizeof(simCtx->log_dir)));
    simCtx->mglevels = 1;
    simCtx->usermg.mglevels = 1;
    PetscCall(PetscCalloc1(1, &simCtx->usermg.mgctx));
    simCtx->usermg.mgctx[0].thislevel = 0;
    PetscCall(PetscCalloc1(simCtx->size, &simCtx->bboxlist));
    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &simCtx->BrownianMotionRNG));
    PetscCall(PetscRandomSetType(simCtx->BrownianMotionRNG, PETSCRAND48));
    PetscCall(PetscRandomSetInterval(simCtx->BrownianMotionRNG, 0.0, 1.0));
    PetscCall(PetscRandomSetSeed(simCtx->BrownianMotionRNG, 12345));
    PetscCall(PetscRandomSeed(simCtx->BrownianMotionRNG));

    user->simCtx = simCtx;
    user->_this = 0;
    user->thislevel = 0;
    user->mglevels = 1;
    user->IM = mx;
    user->JM = my;
    user->KM = mz;
    user->Min_X = 0.0;
    user->Min_Y = 0.0;
    user->Min_Z = 0.0;
    user->Max_X = 1.0;
    user->Max_Y = 1.0;
    user->Max_Z = 1.0;
    user->rx = 1.0;
    user->ry = 1.0;
    user->rz = 1.0;
    user->bbox.min_coords.x = 0.0;
    user->bbox.min_coords.y = 0.0;
    user->bbox.min_coords.z = 0.0;
    user->bbox.max_coords.x = 1.0;
    user->bbox.max_coords.y = 1.0;
    user->bbox.max_coords.z = 1.0;
    simCtx->usermg.mgctx[0].user = user;
    for (PetscMPIInt rank_idx = 0; rank_idx < simCtx->size; ++rank_idx) {
        simCtx->bboxlist[rank_idx].min_coords.x = 0.0;
        simCtx->bboxlist[rank_idx].min_coords.y = 0.0;
        simCtx->bboxlist[rank_idx].min_coords.z = 0.0;
        simCtx->bboxlist[rank_idx].max_coords.x = 1.0;
        simCtx->bboxlist[rank_idx].max_coords.y = 1.0;
        simCtx->bboxlist[rank_idx].max_coords.z = 1.0;
    }

    PetscCall(DMDACreate3d(PETSC_COMM_WORLD,
                           x_boundary, y_boundary, z_boundary,
                           DMDA_STENCIL_BOX,
                           da_mx, da_my, da_mz,
                           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                           1, stencil_width,
                           NULL, NULL, NULL,
                           &user->da));
    PetscCall(DMSetUp(user->da));
    PetscCall(DMGetCoordinateDM(user->da, &user->fda));
    PetscCall(PetscObjectReference((PetscObject)user->fda));
    PetscCall(DMDASetUniformCoordinates(user->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0));
    PetscCall(DMDAGetLocalInfo(user->da, &user->info));
    PetscCall(ComputeLocalBoundingBox(user, &user->bbox));
    PetscCall(GatherAllBoundingBoxes(user, &boxes));
    PetscCall(BroadcastAllBoundingBoxes(user, &boxes));
    for (PetscMPIInt rank_idx = 0; rank_idx < simCtx->size; ++rank_idx) {
        simCtx->bboxlist[rank_idx] = boxes[rank_idx];
    }
    free(boxes);

    PetscCall(CreateZeroedGlobalVector(user->da, &user->P));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lP));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Phi));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lPhi));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Nvert));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lNvert));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->ParticleCount));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lParticleCount));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Psi));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lPsi));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Qcrit));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->P_nodal));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Psi_nodal));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Aj));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lAj));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->Diffusivity));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lDiffusivity));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->B));
    PetscCall(CreateZeroedGlobalVector(user->da, &user->R));

    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Ucat));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->lUcat));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Ucont));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->lUcont));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Csi));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->lCsi));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->Eta));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lEta));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->Zet));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lZet));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Cent));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->lCent));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->GridSpace));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->lGridSpace));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->Centx));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->Centy));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->Centz));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Ucat_nodal));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->DiffusivityGradient));
    PetscCall(CreateZeroedLocalVector(user->fda, &user->lDiffusivityGradient));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Rhs));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->dUcont));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->pUcont));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Bcs.Ubcs));
    PetscCall(CreateZeroedGlobalVector(user->fda, &user->Bcs.Uch));

    PetscCall(CreateZeroedDuplicate(user->Ucont, &user->Ucont_o));
    PetscCall(CreateZeroedDuplicate(user->lUcont, &user->lUcont_o));
    PetscCall(CreateZeroedDuplicate(user->Ucont, &user->Ucont_rm1));
    PetscCall(CreateZeroedDuplicate(user->lUcont, &user->lUcont_rm1));
    PetscCall(CreateZeroedDuplicate(user->Ucat, &user->Ucat_o));
    PetscCall(CreateZeroedDuplicate(user->P, &user->P_o));
    PetscCall(CreateZeroedDuplicate(user->Nvert, &user->Nvert_o));
    PetscCall(CreateZeroedLocalVector(user->da, &user->lNvert_o));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->ICsi));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lICsi));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->IEta));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lIEta));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->IZet));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lIZet));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->JCsi));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lJCsi));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->JEta));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lJEta));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->JZet));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lJZet));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->KCsi));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lKCsi));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->KEta));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lKEta));
    PetscCall(CreateZeroedDuplicate(user->Csi, &user->KZet));
    PetscCall(CreateZeroedDuplicate(user->lCsi, &user->lKZet));
    PetscCall(CreateZeroedDuplicate(user->Aj, &user->IAj));
    PetscCall(CreateZeroedDuplicate(user->lAj, &user->lIAj));
    PetscCall(CreateZeroedDuplicate(user->Aj, &user->JAj));
    PetscCall(CreateZeroedDuplicate(user->lAj, &user->lJAj));
    PetscCall(CreateZeroedDuplicate(user->Aj, &user->KAj));
    PetscCall(CreateZeroedDuplicate(user->lAj, &user->lKAj));

    PetscCall(PicurvPopulateUniformCellCenters(user));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    *simCtx_out = simCtx;
    *user_out = user;
    PetscFunctionReturn(0);
}

/**
 * @brief Builds minimal SimCtx and UserCtx fixtures for C unit tests.
 */
PetscErrorCode PicurvCreateMinimalContexts(SimCtx **simCtx_out, UserCtx **user_out, PetscInt mx, PetscInt my, PetscInt mz)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(simCtx_out, user_out, mx, my, mz, PETSC_FALSE, PETSC_FALSE, PETSC_FALSE));
    PetscFunctionReturn(0);
}

/**
 * @brief Populates cell center coordinates for a uniform grid on [0,1]^3.
 *
 * Uses the shifted-index convention: cell (i,j,k) is stored at array index
 * (i+1, j+1, k+1). For mx cells on [0,1], cell i has center at (i+0.5)/mx.
 */
PetscErrorCode PicurvPopulateUniformCellCenters(UserCtx *user)
{
    Cmpnts ***cent = NULL;
    PetscInt mx = user->IM;
    PetscInt my = user->JM;
    PetscInt mz = user->KM;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Cent, &cent));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; i++) {
                cent[k][j][i].x = (i - 0.5) / (PetscReal)mx;
                cent[k][j][i].y = (j - 0.5) / (PetscReal)my;
                cent[k][j][i].z = (k - 0.5) / (PetscReal)mz;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Cent, &cent));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Cent, INSERT_VALUES, user->lCent));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Cent, INSERT_VALUES, user->lCent));
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
    Cmpnts ***icsi = NULL;
    Cmpnts ***ieta = NULL;
    Cmpnts ***izet = NULL;
    Cmpnts ***jcsi = NULL;
    Cmpnts ***jeta = NULL;
    Cmpnts ***jzet = NULL;
    Cmpnts ***kcsi = NULL;
    Cmpnts ***keta = NULL;
    Cmpnts ***kzet = NULL;
    PetscReal ***aj = NULL;
    PetscReal ***iaj = NULL;
    PetscReal ***jaj = NULL;
    PetscReal ***kaj = NULL;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecGetArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecGetArray(user->fda, user->Zet, &zet));
    PetscCall(DMDAVecGetArray(user->fda, user->ICsi, &icsi));
    PetscCall(DMDAVecGetArray(user->fda, user->IEta, &ieta));
    PetscCall(DMDAVecGetArray(user->fda, user->IZet, &izet));
    PetscCall(DMDAVecGetArray(user->fda, user->JCsi, &jcsi));
    PetscCall(DMDAVecGetArray(user->fda, user->JEta, &jeta));
    PetscCall(DMDAVecGetArray(user->fda, user->JZet, &jzet));
    PetscCall(DMDAVecGetArray(user->fda, user->KCsi, &kcsi));
    PetscCall(DMDAVecGetArray(user->fda, user->KEta, &keta));
    PetscCall(DMDAVecGetArray(user->fda, user->KZet, &kzet));
    PetscCall(DMDAVecGetArray(user->da, user->Aj, &aj));
    PetscCall(DMDAVecGetArray(user->da, user->IAj, &iaj));
    PetscCall(DMDAVecGetArray(user->da, user->JAj, &jaj));
    PetscCall(DMDAVecGetArray(user->da, user->KAj, &kaj));

    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                csi[k][j][i].x = 1.0; csi[k][j][i].y = 0.0; csi[k][j][i].z = 0.0;
                eta[k][j][i].x = 0.0; eta[k][j][i].y = 1.0; eta[k][j][i].z = 0.0;
                zet[k][j][i].x = 0.0; zet[k][j][i].y = 0.0; zet[k][j][i].z = 1.0;
                icsi[k][j][i] = csi[k][j][i];
                ieta[k][j][i] = eta[k][j][i];
                izet[k][j][i] = zet[k][j][i];
                jcsi[k][j][i] = csi[k][j][i];
                jeta[k][j][i] = eta[k][j][i];
                jzet[k][j][i] = zet[k][j][i];
                kcsi[k][j][i] = csi[k][j][i];
                keta[k][j][i] = eta[k][j][i];
                kzet[k][j][i] = zet[k][j][i];
                aj[k][j][i] = 1.0;
                iaj[k][j][i] = 1.0;
                jaj[k][j][i] = 1.0;
                kaj[k][j][i] = 1.0;
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Zet, &zet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->ICsi, &icsi));
    PetscCall(DMDAVecRestoreArray(user->fda, user->IEta, &ieta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->IZet, &izet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->JCsi, &jcsi));
    PetscCall(DMDAVecRestoreArray(user->fda, user->JEta, &jeta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->JZet, &jzet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->KCsi, &kcsi));
    PetscCall(DMDAVecRestoreArray(user->fda, user->KEta, &keta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->KZet, &kzet));
    PetscCall(DMDAVecRestoreArray(user->da, user->Aj, &aj));
    PetscCall(DMDAVecRestoreArray(user->da, user->IAj, &iaj));
    PetscCall(DMDAVecRestoreArray(user->da, user->JAj, &jaj));
    PetscCall(DMDAVecRestoreArray(user->da, user->KAj, &kaj));

    PetscCall(DMGlobalToLocalBegin(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->ICsi, INSERT_VALUES, user->lICsi));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->ICsi, INSERT_VALUES, user->lICsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->IEta, INSERT_VALUES, user->lIEta));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->IEta, INSERT_VALUES, user->lIEta));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->IZet, INSERT_VALUES, user->lIZet));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->IZet, INSERT_VALUES, user->lIZet));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->JCsi, INSERT_VALUES, user->lJCsi));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->JCsi, INSERT_VALUES, user->lJCsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->JEta, INSERT_VALUES, user->lJEta));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->JEta, INSERT_VALUES, user->lJEta));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->JZet, INSERT_VALUES, user->lJZet));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->JZet, INSERT_VALUES, user->lJZet));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->KCsi, INSERT_VALUES, user->lKCsi));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->KCsi, INSERT_VALUES, user->lKCsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->KEta, INSERT_VALUES, user->lKEta));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->KEta, INSERT_VALUES, user->lKEta));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->KZet, INSERT_VALUES, user->lKZet));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->KZet, INSERT_VALUES, user->lKZet));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Nvert, INSERT_VALUES, user->lNvert));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalBegin(user->da, user->IAj, INSERT_VALUES, user->lIAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->IAj, INSERT_VALUES, user->lIAj));
    PetscCall(DMGlobalToLocalBegin(user->da, user->JAj, INSERT_VALUES, user->lJAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->JAj, INSERT_VALUES, user->lJAj));
    PetscCall(DMGlobalToLocalBegin(user->da, user->KAj, INSERT_VALUES, user->lKAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->KAj, INSERT_VALUES, user->lKAj));
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
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "DMSwarm_CellID", 3, PETSC_INT));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "weight", 3, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "Diffusivity", 1, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "DiffusivityGradient", 3, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "Psi", 1, PETSC_REAL));
    PetscCall(RegisterSwarmFieldForTests(user->swarm, "DMSwarm_location_status", 1, PETSC_INT));
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
 * @brief Builds a tiny runtime context through the real setup path for behavior-level tests.
 */

PetscErrorCode PicurvBuildTinyRuntimeContext(const char *bcs_contents,
                                             PetscBool enable_particles,
                                             SimCtx **simCtx_out,
                                             UserCtx **user_out,
                                             char *tmpdir,
                                             size_t tmpdir_len)
{
    char control_path[PETSC_MAX_PATH_LEN];
    SimCtx *simCtx = NULL;

    PetscFunctionBeginUser;
    PetscCheck(simCtx_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx output cannot be NULL.");

    PetscCall(PetscOptionsClear(NULL));
    PetscCall(PrepareTinyRuntimeConfig(bcs_contents, enable_particles, tmpdir, tmpdir_len, control_path, sizeof(control_path)));
    PetscCall(PetscOptionsSetValue(NULL, "-control_file", control_path));
    PetscCall(CreateSimulationContext(0, NULL, &simCtx));
    simCtx->exec_mode = EXEC_MODE_SOLVER;
    PetscCall(SetupSimulationEnvironment(simCtx));
    PetscCall(SetupGridAndSolvers(simCtx));
    PetscCall(SetupBoundaryConditions(simCtx));
    PetscCall(SetupDomainRankInfo(simCtx));

    *simCtx_out = simCtx;
    if (user_out) {
        *user_out = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    }
    PetscFunctionReturn(0);
}
/**
 * @brief Finalizes and frees a runtime context built by `PicurvBuildTinyRuntimeContext`.
 */

PetscErrorCode PicurvDestroyRuntimeContext(SimCtx **simCtx_ptr)
{
    PetscFunctionBeginUser;
    if (simCtx_ptr && *simCtx_ptr) {
        PetscCall(FinalizeSimulation(*simCtx_ptr));
        PetscCall(PetscFree(*simCtx_ptr));
        *simCtx_ptr = NULL;
    }
    PetscCall(PetscOptionsClear(NULL));
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

    if (simCtx) {
        PetscCall(DestroySolutionConvergenceState(simCtx));
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
        PetscCall(DestroyVecIfSet(&user->CS));
        PetscCall(DestroyVecIfSet(&user->lCs));
        PetscCall(DestroyVecIfSet(&user->Nu_t));
        PetscCall(DestroyVecIfSet(&user->lNu_t));
        PetscCall(DestroyVecIfSet(&user->Ucont_o));
        PetscCall(DestroyVecIfSet(&user->lUcont_o));
        PetscCall(DestroyVecIfSet(&user->Ucat_o));
        PetscCall(DestroyVecIfSet(&user->P_o));
        PetscCall(DestroyVecIfSet(&user->Nvert_o));
        PetscCall(DestroyVecIfSet(&user->lNvert_o));
        PetscCall(DestroyVecIfSet(&user->Ucont_rm1));
        PetscCall(DestroyVecIfSet(&user->lUcont_rm1));
        PetscCall(DestroyVecIfSet(&user->Rhs));
        PetscCall(DestroyVecIfSet(&user->dUcont));
        PetscCall(DestroyVecIfSet(&user->pUcont));
        PetscCall(DestroyVecIfSet(&user->CellFieldAtCorner));
        PetscCall(DestroyVecIfSet(&user->lCellFieldAtCorner));
        PetscCall(DestroyVecIfSet(&user->B));
        PetscCall(DestroyVecIfSet(&user->R));
        PetscCall(DestroyVecIfSet(&user->Cent));
        PetscCall(DestroyVecIfSet(&user->lCent));
        PetscCall(DestroyVecIfSet(&user->GridSpace));
        PetscCall(DestroyVecIfSet(&user->lGridSpace));
        PetscCall(DestroyVecIfSet(&user->Centx));
        PetscCall(DestroyVecIfSet(&user->Centy));
        PetscCall(DestroyVecIfSet(&user->Centz));
        PetscCall(DestroyVecIfSet(&user->Ucat_nodal));
        PetscCall(DestroyVecIfSet(&user->Bcs.Ubcs));
        PetscCall(DestroyVecIfSet(&user->Bcs.Uch));
        PetscCall(DestroyMatIfSet(&user->A));
        PetscCall(DestroyMatIfSet(&user->C));
        PetscCall(DestroyKSPIfSet(&user->ksp));
        PetscCall(DestroyNullSpaceIfSet(&user->nullsp));
        PetscCall(PetscFree(user->RankCellInfoMap));
        user->RankCellInfoMap = NULL;
        for (PetscInt face = BC_FACE_NEG_X; face <= BC_FACE_POS_Z; ++face) {
            if (user->boundary_faces[face].params) {
                FreeBC_ParamList(user->boundary_faces[face].params);
                user->boundary_faces[face].params = NULL;
            }
        }

        PetscCall(DestroyDMIfSet(&user->fda));
        PetscCall(DestroyDMIfSet(&user->da));
        PetscCall(PetscFree(user));
        if (user_ptr) {
            *user_ptr = NULL;
        }
    }

    if (simCtx) {
        PetscCall(DestroyRandomIfSet(&simCtx->BrownianMotionRNG));
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
