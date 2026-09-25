/**
 * @file test_statistics.c
 * @brief C unit tests for particle statistics kernels and CSV output.
 */

#include "test_support.h"

#include "particle_statistics.h"

#include <stdio.h>
#include <string.h>
/**
 * @brief Tests MSD CSV generation for a populated particle swarm.
 */

static PetscErrorCode TestComputeParticleMSDWritesCSV(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char csv_prefix[PETSC_MAX_PATH_LEN];
    char csv_path[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;
    char header[512];
    char row[512];
    PetscReal (*pos_arr)[3] = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));
    simCtx->dt = 1.0;
    simCtx->ren = 1.0;
    simCtx->schmidt_number = 1.0;
    simCtx->psrc_x = 0.0;
    simCtx->psrc_y = 0.0;
    simCtx->psrc_z = 0.0;

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(csv_prefix, sizeof(csv_prefix), "%s/stats", tmpdir));
    PetscCall(PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", csv_prefix));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 1.0;
    pos_arr[0][1] = 0.0;
    pos_arr[0][2] = 0.0;
    pos_arr[1][0] = -1.0;
    pos_arr[1][1] = 0.0;
    pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(ComputeParticleMSD(user, csv_prefix, 1));
    PetscCall(PicurvAssertFileExists(csv_path, "ComputeParticleMSD should emit a CSV summary"));

    file = fopen(csv_path, "r");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open generated stats CSV '%s'.", csv_path);
    PetscCheck(fgets(header, sizeof(header), file) != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV header is missing.");
    PetscCheck(fgets(row, sizeof(row), file) != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV data row is missing.");
    fclose(file);

    PetscCall(PicurvAssertBool((PetscBool)(strstr(header, "MSD_total") != NULL),
                               "MSD CSV header should contain MSD_total"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(row, "2") != NULL),
                               "MSD CSV row should include the particle count"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that a dimensionalizing post-processor reports MSD columns physically.
 * @details Two particles at x = +/-1 (solver units) about the origin, step 1, dt = 1. At
 *          L = 2 and U = 4 the time column is scaled by L/U = 0.5, MSD by L^2 = 4, and the
 *          RMS radius and centre of mass by L; the relative error is unchanged.
 */
static PetscErrorCode TestComputeParticleMSDDimensionalizedColumns(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char csv_prefix[PETSC_MAX_PATH_LEN];
    char csv_path[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;
    char header[512];
    int step = 0;
    double values[15];
    PetscReal (*pos_arr)[3] = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));
    PetscCall(PetscCalloc1(1, &simCtx->pps));
    simCtx->pps->dimensionalize = PETSC_TRUE;
    simCtx->scaling.L_ref = 2.0;
    simCtx->scaling.U_ref = 4.0;
    simCtx->scaling.rho_ref = 1.0;
    simCtx->dt = 1.0;
    simCtx->ren = 1.0;
    simCtx->schmidt_number = 1.0;

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(csv_prefix, sizeof(csv_prefix), "%s/stats", tmpdir));
    PetscCall(PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", csv_prefix));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 1.0;  pos_arr[0][1] = 0.0; pos_arr[0][2] = 0.0;
    pos_arr[1][0] = -1.0; pos_arr[1][1] = 0.0; pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(ComputeParticleMSD(user, csv_prefix, 1));
    file = fopen(csv_path, "r");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s'.", csv_path);
    PetscCheck(fgets(header, sizeof(header), file) != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV header is missing.");
    PetscCheck(fscanf(file, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &step,
                      &values[0], &values[1], &values[2], &values[3], &values[4], &values[5], &values[6],
                      &values[7], &values[8], &values[9], &values[10], &values[11], &values[12],
                      &values[13], &values[14]) == 16,
               PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "MSD row has an unexpected layout.");
    fclose(file);

    /* Solver units: t = 1, MSD_x = 1, r = 1, r_theory = sqrt(6), com = 0. */
    PetscCall(PicurvAssertRealNear(0.5, values[0], 1.0e-6, "time is scaled by L/U"));
    PetscCall(PicurvAssertRealNear(4.0, values[2], 1.0e-6, "MSD_x is scaled by L^2"));
    PetscCall(PicurvAssertRealNear(4.0, values[5], 1.0e-6, "MSD_total is scaled by L^2"));
    PetscCall(PicurvAssertRealNear(2.0, values[6], 1.0e-6, "the RMS radius is scaled by L"));
    PetscCall(PicurvAssertRealNear(2.0 * PetscSqrtReal(6.0), values[7], 1.0e-5, "the theory radius is scaled by L"));
    PetscCall(PicurvAssertRealNear(100.0 * (PetscSqrtReal(6.0) - 1.0) / PetscSqrtReal(6.0), values[8], 1.0e-3,
                                   "the relative error is the same in either unit system"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that empty swarms do not emit MSD CSV output.
 */

static PetscErrorCode TestComputeParticleMSDEmptySwarmNoOutput(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char csv_prefix[PETSC_MAX_PATH_LEN];
    char csv_path[PETSC_MAX_PATH_LEN];
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 0, "ske"));

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(csv_prefix, sizeof(csv_prefix), "%s/stats", tmpdir));
    PetscCall(PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", csv_prefix));

    PetscCall(ComputeParticleMSD(user, csv_prefix, 5));
    PetscCall(VerifyPathExistence(csv_path, PETSC_FALSE, PETSC_TRUE, "MSD CSV for empty swarm", &exists));
    PetscCall(PicurvAssertBool((PetscBool)!exists,
                               "ComputeParticleMSD should not write output when no particles are present"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-post-statistics PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"compute-particle-msd-writes-csv", TestComputeParticleMSDWritesCSV},
        {"compute-particle-msd-dimensionalized-columns", TestComputeParticleMSDDimensionalizedColumns},
        {"compute-particle-msd-empty-swarm-no-output", TestComputeParticleMSDEmptySwarmNoOutput},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv statistics tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-post-statistics", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
