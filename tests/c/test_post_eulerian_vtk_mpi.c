/**
 * @file test_post_eulerian_vtk_mpi.c
 * @brief Focused multirank regression for instantaneous Eulerian VTK output.
 */

#include "test_support.h"

#include "postprocessor.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** @brief Compare one scalar from the raw VTK payload with its deterministic truth. */
static PetscBool ScalarMatches(PetscScalar actual, PetscReal expected)
{
    return (PetscBool)(PetscAbsReal(PetscRealPart(actual) - expected) <= 1.0e-12);
}

/** @brief Read and exhaustively validate every appended value in a generated VTS file. */
static PetscBool ValidateEulerianVTKFile(const char *path, PetscInt mx, PetscInt my, PetscInt mz)
{
    FILE *file = NULL;
    char line[1024];
    char footer[128];
    char expected_extent[128];
    PetscBool saw_pressure = PETSC_FALSE;
    PetscBool saw_velocity = PETSC_FALSE;
    PetscBool saw_qcrit = PETSC_FALSE;
    PetscBool saw_extent = PETSC_FALSE;
    uint32_t block_size = 0;
    PetscInt npoints = mx * my * mz;
    PetscScalar *coords = NULL;
    PetscScalar *pressure = NULL;
    PetscScalar *velocity = NULL;
    PetscScalar *qcrit = NULL;
    size_t footer_size = 0;
    PetscBool valid = PETSC_FALSE;

    file = fopen(path, "rb");
    if (!file) return PETSC_FALSE;
    snprintf(expected_extent, sizeof(expected_extent),
             "WholeExtent=\"0 %d 0 %d 0 %d\"", (int)(mx - 1), (int)(my - 1), (int)(mz - 1));

    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Name=\"P_nodal\"")) saw_pressure = PETSC_TRUE;
        if (strstr(line, "Name=\"Ucat_nodal\"")) saw_velocity = PETSC_TRUE;
        if (strstr(line, "Name=\"Qcrit\"")) saw_qcrit = PETSC_TRUE;
        if (strstr(line, expected_extent)) saw_extent = PETSC_TRUE;
        if (strstr(line, "<AppendedData encoding=\"raw\">")) break;
    }
    if (!saw_pressure || !saw_velocity || !saw_qcrit || !saw_extent || fgetc(file) != '_') goto cleanup;

    /* Coordinates. */
    if (fread(&block_size, sizeof(block_size), 1, file) != 1 ||
        block_size != (uint32_t)(3 * npoints * (PetscInt)sizeof(PetscScalar))) goto cleanup;
    coords = malloc((size_t)(3 * npoints) * sizeof(*coords));
    if (!coords ||
        fread(coords, sizeof(PetscScalar), (size_t)(3 * npoints), file) != (size_t)(3 * npoints)) goto cleanup;

    /* P_nodal. */
    if (fread(&block_size, sizeof(block_size), 1, file) != 1 ||
        block_size != (uint32_t)(npoints * (PetscInt)sizeof(PetscScalar))) goto cleanup;
    pressure = malloc((size_t)npoints * sizeof(*pressure));
    if (!pressure ||
        fread(pressure, sizeof(PetscScalar), (size_t)npoints, file) != (size_t)npoints) goto cleanup;

    /* Ucat_nodal. */
    if (fread(&block_size, sizeof(block_size), 1, file) != 1 ||
        block_size != (uint32_t)(3 * npoints * (PetscInt)sizeof(PetscScalar))) goto cleanup;
    velocity = malloc((size_t)(3 * npoints) * sizeof(*velocity));
    if (!velocity ||
        fread(velocity, sizeof(PetscScalar), (size_t)(3 * npoints), file) != (size_t)(3 * npoints)) goto cleanup;

    /* Qcrit. */
    if (fread(&block_size, sizeof(block_size), 1, file) != 1 ||
        block_size != (uint32_t)(npoints * (PetscInt)sizeof(PetscScalar))) goto cleanup;
    qcrit = malloc((size_t)npoints * sizeof(*qcrit));
    if (!qcrit ||
        fread(qcrit, sizeof(PetscScalar), (size_t)npoints, file) != (size_t)npoints) goto cleanup;

    for (PetscInt k = 0; k < mz; ++k) {
        for (PetscInt j = 0; j < my; ++j) {
            for (PetscInt i = 0; i < mx; ++i) {
                const PetscInt point = (k * my + j) * mx + i;
                const PetscReal base = (PetscReal)(10000 * k + 100 * j + i);

                if (!ScalarMatches(coords[3 * point + 0], (PetscReal)i / (PetscReal)mx) ||
                    !ScalarMatches(coords[3 * point + 1], (PetscReal)j / (PetscReal)my) ||
                    !ScalarMatches(coords[3 * point + 2], (PetscReal)k / (PetscReal)mz) ||
                    !ScalarMatches(pressure[point], base + 0.25) ||
                    !ScalarMatches(velocity[3 * point + 0], base + 1.0) ||
                    !ScalarMatches(velocity[3 * point + 1], base + 2.0) ||
                    !ScalarMatches(velocity[3 * point + 2], base + 3.0) ||
                    !ScalarMatches(qcrit[point], -base - 0.5)) goto cleanup;
            }
        }
    }

    footer_size = fread(footer, 1, sizeof(footer) - 1, file);
    footer[footer_size] = '\0';
    if (!strstr(footer, "</AppendedData>") || !strstr(footer, "</VTKFile>")) goto cleanup;

    valid = PETSC_TRUE;

cleanup:
    fclose(file);
    free(coords);
    free(pressure);
    free(velocity);
    free(qcrit);
    return valid;
}

/**
 * @brief Ensures every Eulerian field and coordinate is correct after collective output.
 */
static PetscErrorCode TestWriteEulerianFileCollectiveMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char vtk_path[PETSC_MAX_PATH_LEN];
    char vtk_tmp_path[PETSC_MAX_PATH_LEN];
    PetscReal ***pressure = NULL;
    PetscReal ***qcrit = NULL;
    Cmpnts ***velocity = NULL;
    PetscMPIInt rank = 0, size = 1;
    PetscMPIInt local_owner = -1, probe_owner = -1;
    PetscBool parsed = PETSC_FALSE;
    PetscBool tmp_exists = PETSC_FALSE;
    PetscBool prefix_set = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-vtk_test_output_prefix",
                                    pps.output_prefix, sizeof(pps.output_prefix), &prefix_set));
    PetscCall(PicurvAssertBool(prefix_set,
                               "the focused Eulerian VTK test requires -vtk_test_output_prefix"));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s_00003.vts", pps.output_prefix));
    PetscCall(PetscSNPrintf(vtk_tmp_path, sizeof(vtk_tmp_path), "%s.tmp", vtk_path));

    PetscCall(DMDAVecGetArray(user->da, user->P_nodal, &pressure));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat_nodal, &velocity));
    PetscCall(DMDAVecGetArray(user->da, user->Qcrit, &qcrit));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                const PetscReal base = (PetscReal)(10000 * k + 100 * j + i);

                pressure[k][j][i] = base + 0.25;
                velocity[k][j][i] = (Cmpnts){base + 1.0, base + 2.0, base + 3.0};
                qcrit[k][j][i] = -base - 0.5;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->P_nodal, &pressure));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat_nodal, &velocity));
    PetscCall(DMDAVecRestoreArray(user->da, user->Qcrit, &qcrit));

    if (user->IM - 1 >= user->info.xs && user->IM - 1 < user->info.xs + user->info.xm &&
        user->JM - 1 >= user->info.ys && user->JM - 1 < user->info.ys + user->info.ym &&
        user->KM - 1 >= user->info.zs && user->KM - 1 < user->info.zs + user->info.zm) {
        local_owner = rank;
    }
    PetscCallMPI(MPI_Allreduce(&local_owner, &probe_owner, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool((PetscBool)(probe_owner == ((size > 1) ? size - 1 : 0)),
                               "the final physical VTK point should have the expected owner"));

    PetscCall(PetscStrncpy(pps.output_fields_instantaneous,
                           "P_nodal,Ucat_nodal,Qcrit",
                           sizeof(pps.output_fields_instantaneous)));

    PetscCall(WriteEulerianFile(user, &pps, 3));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(PicurvAssertFileExists(vtk_path,
                                     "multirank WriteEulerianFile should emit one combined .vts file"));
    PetscCall(PetscTestFile(vtk_tmp_path, 'r', &tmp_exists));
    PetscCall(PicurvAssertBool((PetscBool)!tmp_exists,
                               "successful multirank Eulerian VTK output should not leave a temporary file"));

    if (rank == 0) {
        parsed = ValidateEulerianVTKFile(vtk_path, user->IM, user->JM, user->KM);
    }
    PetscCallMPI(MPI_Bcast(&parsed, 1, MPIU_BOOL, 0, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool(parsed,
                               "every coordinate and Eulerian field value in the generated VTK should match"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Runs the focused serial/MPI Eulerian VTK regression. */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"write-eulerian-file-collective-multi-rank", TestWriteEulerianFileCollectiveMultiRank},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv multirank Eulerian VTK tests");
    if (ierr) return (int)ierr;

    ierr = PicurvRunTests("unit-post-eulerian-vtk-mpi", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
