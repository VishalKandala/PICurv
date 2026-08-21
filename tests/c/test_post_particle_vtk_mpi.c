/**
 * @file test_post_particle_vtk_mpi.c
 * @brief Focused multirank regression for particle VTK output.
 */

#include "test_support.h"

#include "postprocessor.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum {
    TEST_PARTICLES = 7,
    TEST_OUTPUT_PARTICLES = 4
};

/** @brief Compare one scalar from the raw VTK payload with its deterministic truth. */
static PetscBool ScalarMatches(PetscScalar actual, PetscReal expected)
{
    return (PetscBool)(PetscAbsReal(PetscRealPart(actual) - expected) <= 1.0e-12);
}

/** @brief Read one length-prefixed raw appended-data block. */
static PetscBool ReadBlock(FILE *file, void *data, size_t count, size_t element_size)
{
    uint32_t block_size = 0;
    const size_t expected_size = count * element_size;

    if (expected_size > UINT32_MAX) return PETSC_FALSE;
    if (fread(&block_size, sizeof(block_size), 1, file) != 1 || block_size != (uint32_t)expected_size) {
        return PETSC_FALSE;
    }
    return (PetscBool)(fread(data, element_size, count, file) == count);
}

/** @brief Read and exhaustively validate every appended value in a generated VTP file. */
static PetscBool ValidateParticleVTKFile(const char *path)
{
    FILE *file = NULL;
    char line[1024];
    char footer[128];
    PetscScalar coords[3 * TEST_OUTPUT_PARTICLES];
    PetscScalar velocity[3 * TEST_OUTPUT_PARTICLES];
    PetscScalar cell_id[3 * TEST_OUTPUT_PARTICLES];
    PetscScalar migration_status[TEST_OUTPUT_PARTICLES];
    PetscScalar pid[TEST_OUTPUT_PARTICLES];
    PetscScalar ske[TEST_OUTPUT_PARTICLES];
    PetscInt connectivity[TEST_OUTPUT_PARTICLES];
    PetscInt offsets[TEST_OUTPUT_PARTICLES];
    PetscBool saw_piece = PETSC_FALSE;
    PetscBool saw_velocity = PETSC_FALSE;
    PetscBool saw_cell_id = PETSC_FALSE;
    PetscBool saw_migration_status = PETSC_FALSE;
    PetscBool saw_pid = PETSC_FALSE;
    PetscBool saw_ske = PETSC_FALSE;
    PetscBool valid = PETSC_FALSE;
    size_t footer_size = 0;

    file = fopen(path, "rb");
    if (!file) return PETSC_FALSE;

    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "NumberOfPoints=\"4\"") && strstr(line, "NumberOfVerts=\"4\"")) saw_piece = PETSC_TRUE;
        if (strstr(line, "Name=\"velocity\"") && strstr(line, "NumberOfComponents=\"3\"")) saw_velocity = PETSC_TRUE;
        if (strstr(line, "Name=\"CellID\"") && strstr(line, "NumberOfComponents=\"3\"")) saw_cell_id = PETSC_TRUE;
        if (strstr(line, "Name=\"Migration Status\"") && strstr(line, "NumberOfComponents=\"1\"")) saw_migration_status = PETSC_TRUE;
        if (strstr(line, "Name=\"pid\"") && strstr(line, "NumberOfComponents=\"1\"")) saw_pid = PETSC_TRUE;
        if (strstr(line, "Name=\"ske\"") && strstr(line, "NumberOfComponents=\"1\"")) saw_ske = PETSC_TRUE;
        if (strstr(line, "<AppendedData encoding=\"raw\">")) break;
    }
    if (!saw_piece || !saw_velocity || !saw_cell_id || !saw_migration_status || !saw_pid || !saw_ske ||
        fgetc(file) != '_') goto cleanup;

    if (!ReadBlock(file, coords, 3 * TEST_OUTPUT_PARTICLES, sizeof(PetscScalar)) ||
        !ReadBlock(file, velocity, 3 * TEST_OUTPUT_PARTICLES, sizeof(PetscScalar)) ||
        !ReadBlock(file, cell_id, 3 * TEST_OUTPUT_PARTICLES, sizeof(PetscScalar)) ||
        !ReadBlock(file, migration_status, TEST_OUTPUT_PARTICLES, sizeof(PetscScalar)) ||
        !ReadBlock(file, pid, TEST_OUTPUT_PARTICLES, sizeof(PetscScalar)) ||
        !ReadBlock(file, ske, TEST_OUTPUT_PARTICLES, sizeof(PetscScalar)) ||
        !ReadBlock(file, connectivity, TEST_OUTPUT_PARTICLES, sizeof(PetscInt)) ||
        !ReadBlock(file, offsets, TEST_OUTPUT_PARTICLES, sizeof(PetscInt))) goto cleanup;

    for (PetscInt p = 0; p < TEST_OUTPUT_PARTICLES; ++p) {
        const PetscInt source = 2 * p;

        if (!ScalarMatches(coords[3 * p + 0], (PetscReal)source + 0.1) ||
            !ScalarMatches(coords[3 * p + 1], (PetscReal)source + 0.2) ||
            !ScalarMatches(coords[3 * p + 2], (PetscReal)source + 0.3) ||
            !ScalarMatches(velocity[3 * p + 0], (PetscReal)source + 10.0) ||
            !ScalarMatches(velocity[3 * p + 1], (PetscReal)source + 20.0) ||
            !ScalarMatches(velocity[3 * p + 2], (PetscReal)source + 30.0) ||
            !ScalarMatches(cell_id[3 * p + 0], (PetscReal)source + 100.0) ||
            !ScalarMatches(cell_id[3 * p + 1], (PetscReal)source + 200.0) ||
            !ScalarMatches(cell_id[3 * p + 2], (PetscReal)source + 300.0) ||
            !ScalarMatches(migration_status[p], (PetscReal)source + 400.0) ||
            !ScalarMatches(pid[p], (PetscReal)source + 100000.0) ||
            !ScalarMatches(ske[p], (PetscReal)source + 0.5) ||
            connectivity[p] != p || offsets[p] != p + 1) goto cleanup;
    }

    footer_size = fread(footer, 1, sizeof(footer) - 1, file);
    footer[footer_size] = '\0';
    if (!strstr(footer, "</AppendedData>") || !strstr(footer, "</VTKFile>")) goto cleanup;

    valid = PETSC_TRUE;

cleanup:
    fclose(file);
    return valid;
}

/** @brief Fill each local swarm entry from its rank-ordered global particle index. */
static PetscErrorCode FillParticleFields(UserCtx *user, PetscInt nlocal, PetscInt global_offset)
{
    PetscReal (*real3)[3] = NULL;
    PetscInt (*int3)[3] = NULL;
    PetscInt *int1 = NULL;
    PetscInt64 *int64 = NULL;
    PetscReal *real1 = NULL;

    PetscFunctionBeginUser;
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&real3));
    for (PetscInt p = 0; p < nlocal; ++p) {
        const PetscReal source = (PetscReal)(global_offset + p);
        real3[p][0] = source + 0.1;
        real3[p][1] = source + 0.2;
        real3[p][2] = source + 0.3;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&real3));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&real3));
    for (PetscInt p = 0; p < nlocal; ++p) {
        const PetscReal source = (PetscReal)(global_offset + p);
        real3[p][0] = source + 10.0;
        real3[p][1] = source + 20.0;
        real3[p][2] = source + 30.0;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&real3));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void *)&int3));
    for (PetscInt p = 0; p < nlocal; ++p) {
        const PetscInt source = global_offset + p;
        int3[p][0] = source + 100;
        int3[p][1] = source + 200;
        int3[p][2] = source + 300;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void *)&int3));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void *)&int1));
    for (PetscInt p = 0; p < nlocal; ++p) int1[p] = global_offset + p + 400;
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void *)&int1));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void *)&int64));
    for (PetscInt p = 0; p < nlocal; ++p) int64[p] = (PetscInt64)(global_offset + p + 100000);
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void *)&int64));

    PetscCall(DMSwarmGetField(user->post_swarm, "ske", NULL, NULL, (void *)&real1));
    for (PetscInt p = 0; p < nlocal; ++p) real1[p] = (PetscReal)(global_offset + p) + 0.5;
    PetscCall(DMSwarmRestoreField(user->post_swarm, "ske", NULL, NULL, (void *)&real1));
    PetscFunctionReturn(0);
}

/** @brief Ensures particle fields are collectively gathered and written without corruption. */
static PetscErrorCode TestWriteParticleFileCollectiveMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char vtk_path[PETSC_MAX_PATH_LEN];
    char vtk_tmp_path[PETSC_MAX_PATH_LEN];
    PetscMPIInt rank = 0, size = 1;
    PetscInt nlocal = 0, global_offset = 0;
    PetscMPIInt local_owner = -1, final_owner = -1;
    PetscBool parsed = PETSC_FALSE;
    PetscBool tmp_exists = PETSC_FALSE;
    PetscBool prefix_set = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    nlocal = TEST_PARTICLES / size + (rank < TEST_PARTICLES % size ? 1 : 0);
    PetscCallMPI(MPI_Exscan(&nlocal, &global_offset, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    if (rank == 0) global_offset = 0;

    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, nlocal, "ske"));
    PetscCall(FillParticleFields(user, nlocal, global_offset));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-vtk_test_output_prefix",
                                    pps.particle_output_prefix, sizeof(pps.particle_output_prefix), &prefix_set));
    PetscCall(PicurvAssertBool(prefix_set,
                               "the focused particle VTK test requires -vtk_test_output_prefix"));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s_00004.vtp", pps.particle_output_prefix));
    PetscCall(PetscSNPrintf(vtk_tmp_path, sizeof(vtk_tmp_path), "%s.tmp", vtk_path));

    if (global_offset <= TEST_PARTICLES - 1 && TEST_PARTICLES - 1 < global_offset + nlocal) local_owner = rank;
    PetscCallMPI(MPI_Allreduce(&local_owner, &final_owner, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool((PetscBool)(final_owner == ((size < TEST_PARTICLES) ? size - 1 : TEST_PARTICLES - 1)),
                               "the final exported particle should have the expected owner"));

    pps.outputParticles = PETSC_TRUE;
    pps.particle_output_freq = 2;
    PetscCall(PetscStrncpy(pps.particle_fields,
                           "position,velocity,CellID,Migration Status,pid,ske",
                           sizeof(pps.particle_fields)));

    PetscCall(WriteParticleFile(user, &pps, 4));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(PicurvAssertFileExists(vtk_path,
                                     "multirank WriteParticleFile should emit one combined .vtp file"));
    PetscCall(PetscTestFile(vtk_tmp_path, 'r', &tmp_exists));
    PetscCall(PicurvAssertBool((PetscBool)!tmp_exists,
                               "successful multirank particle VTK output should not leave a temporary file"));

    if (rank == 0) parsed = ValidateParticleVTKFile(vtk_path);
    PetscCallMPI(MPI_Bcast(&parsed, 1, MPIU_BOOL, 0, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool(parsed,
                               "every coordinate, particle field, and vertex index in the generated VTK should match"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"write-particle-file-collective-multi-rank", TestWriteParticleFileCollectiveMultiRank},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv multirank particle VTK tests");
    if (ierr) return (int)ierr;

    ierr = PicurvRunTests("unit-post-particle-vtk-mpi", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
