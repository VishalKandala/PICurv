#include "test_support.h"

#include "particle_statistics.h"
#include "postprocessing_kernels.h"

#include <stdio.h>
#include <string.h>

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

static PetscErrorCode TestComputeParticleMSD(void)
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
    pos_arr[0][0] = 1.0; pos_arr[0][1] = 0.0; pos_arr[0][2] = 0.0;
    pos_arr[1][0] = -1.0; pos_arr[1][1] = 0.0; pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(ComputeParticleMSD(user, csv_prefix, 1));
    PetscCall(PicurvAssertFileExists(csv_path, "ComputeParticleMSD should emit a CSV summary"));

    file = fopen(csv_path, "r");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open generated stats CSV '%s'.", csv_path);
    PetscCheck(fgets(header, sizeof(header), file) != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV header is missing.");
    PetscCheck(fgets(row, sizeof(row), file) != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV data row is missing.");
    fclose(file);

    PetscCall(PicurvAssertBool((PetscBool)(strstr(header, "MSD_total") != NULL), "MSD CSV header should contain MSD_total"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(row, "2") != NULL), "MSD CSV row should include the particle count"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"compute-specific-ke", TestComputeSpecificKE},
        {"compute-particle-msd", TestComputeParticleMSD},
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
