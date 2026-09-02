/**
 * @file io.c
 * @brief Implementation of data input/output routines, focusing on grid configuration.
 *
 * This module provides functions to parse grid geometry information, either from
 * command-line options for programmatically generated grids or by reading the
 * header of a grid definition file.
 */

#include "io.h"
#include "checksum.h"
#include "field_catalog.h"
#include "particle_field_catalog.h"
#include "statistics_accumulator.h"

#include <errno.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#define PICURV_CHECKPOINT_FORMAT "picurv-checkpoint"
#define PICURV_CHECKPOINT_VERSION 1
#define PICURV_CHECKPOINTS_DIRECTORY "checkpoints"
#define PICURV_EULERIAN_DIRECTORY "eulerian"
#define PICURV_PARTICLE_DIRECTORY "particles"
#define PICURV_STATISTICS_DIRECTORY "statistics"
#define PICURV_CHECKPOINT_STEP_WIDTH 12

// =============================================================================
//          STATIC (PRIVATE) VARIABLES FOR ONE-TIME FILE READ
// =============================================================================

/** @brief Stores the number of blocks read from the grid file. */
static PetscInt  g_nblk_from_file = 0;
/** @brief Caches the IM dimensions for all blocks read from the grid file. */
static PetscInt* g_IMs_from_file = NULL;
/** @brief Caches the JM dimensions for all blocks read from the grid file. */
static PetscInt* g_JMs_from_file = NULL;
/** @brief Caches the KM dimensions for all blocks read from the grid file. */
static PetscInt* g_KMs_from_file = NULL;
/** @brief A flag to ensure the grid file is read only once. */
static PetscBool g_file_has_been_read = PETSC_FALSE;

/**
 * @brief Copies the owned entries of a ghosted scalar DMDA vector to its global vector.
 *
 * DMLocalToGlobal with INSERT_VALUES is unsupported for multidirection-periodic
 * DMDAs in PETSc.  Explicitly copying only the uniquely owned region preserves
 * INSERT semantics without reducing duplicate periodic ghost entries.
 */
static PetscErrorCode CopyOwnedLocalScalarToGlobal(DM dm, Vec local_vec, Vec global_vec)
{
    DMDALocalInfo info;
    const PetscScalar ***local_array = NULL;
    PetscScalar ***global_array = NULL;

    PetscFunctionBeginUser;

    PetscCall(DMDAGetLocalInfo(dm, &info));
    PetscCall(DMDAVecGetArrayRead(dm, local_vec, &local_array));
    PetscCall(DMDAVecGetArray(dm, global_vec, &global_array));
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i)
                global_array[k][j][i] = local_array[k][j][i];
    PetscCall(DMDAVecRestoreArray(dm, global_vec, &global_array));
    PetscCall(DMDAVecRestoreArrayRead(dm, local_vec, &local_array));
    PetscFunctionReturn(0);
}

/** @brief Return whether a catalogued field belongs in the current checkpoint. */
static PetscBool CheckpointFieldIsEnabled(const SimCtx *simCtx, const FieldDescriptor *descriptor)
{
    const unsigned int availability = descriptor ? descriptor->availability : FIELD_AVAILABILITY_ALWAYS;

    if (!simCtx || !descriptor || !(descriptor->capabilities & FIELD_CAPABILITY_CHECKPOINT)) return PETSC_FALSE;
    if ((availability & FIELD_AVAILABILITY_TURBULENCE) && !(simCtx->les || simCtx->rans)) return PETSC_FALSE;
    if ((availability & FIELD_AVAILABILITY_LES_DYNAMIC) && simCtx->les != DYNAMIC_SMAGORINSKY) return PETSC_FALSE;
    if ((availability & FIELD_AVAILABILITY_RANS) && !simCtx->rans) return PETSC_FALSE;
    if ((availability & FIELD_AVAILABILITY_PARTICLES) && simCtx->np <= 0) return PETSC_FALSE;
    return PETSC_TRUE;
}

/**
 * @brief Format any level of the statistics subtree, from the root down to one payload.
 *
 * This is the only place the statistics layout is written down. The directory
 * creator, the payload writer, the manifest inventory, and the restart reader all
 * ask for the level they need, so a path cannot be built one way and looked for
 * another.
 *
 * Statistics payloads are block scoped exactly as Eulerian payloads are, with a
 * window level above the block level, so a multiblock run keeps one accumulator tree
 * per block under each window.
 *
 * Pass @p root NULL for a bundle-relative path, and stop early by passing a negative
 * @p window or @p block, or a NULL @p payload_name.
 */
static PetscErrorCode FormatStatisticsPath(const char *root, PetscInt window, PetscInt block,
                                           const char *payload_name, char *path, size_t path_size)
{
    size_t used = 0;

    PetscFunctionBeginUser;
    PetscCheck(path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output path is required.");
    if (root) PetscCall(PetscSNPrintf(path, path_size, "%s/%s", root, PICURV_STATISTICS_DIRECTORY));
    else PetscCall(PetscSNPrintf(path, path_size, "%s", PICURV_STATISTICS_DIRECTORY));
    if (window < 0) PetscFunctionReturn(0);

    PetscCall(PetscStrlen(path, &used));
    PetscCall(PetscSNPrintf(path + used, path_size - used, "/window_%04" PetscInt_FMT, window));
    if (block < 0) PetscFunctionReturn(0);

    PetscCall(PetscStrlen(path, &used));
    PetscCall(PetscSNPrintf(path + used, path_size - used, "/block_%04" PetscInt_FMT, block));
    if (!payload_name) PetscFunctionReturn(0);

    PetscCall(PetscStrlen(path, &used));
    PetscCall(PetscSNPrintf(path + used, path_size - used, "/%s.dat", payload_name));
    PetscFunctionReturn(0);
}

/** @brief Gather a vector in decomposition-independent natural ordering onto rank zero. */
static PetscErrorCode GatherVectorToRankZero(Vec field_vec, Vec *sequential_vec)
{
    DM dm = NULL;
    const char *dm_type = NULL;
    Vec natural_vec = NULL;
    VecScatter scatter = NULL;

    PetscFunctionBeginUser;
    PetscCheck(field_vec != NULL && sequential_vec != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "A source vector and output vector pointer are required.");
    *sequential_vec = NULL;

    PetscCall(VecGetDM(field_vec, &dm));
    if (dm) PetscCall(DMGetType(dm, &dm_type));
    if (dm_type && !strcmp(dm_type, DMDA)) {
        PetscCall(DMDACreateNaturalVector(dm, &natural_vec));
        PetscCall(DMDAGlobalToNaturalBegin(dm, field_vec, INSERT_VALUES, natural_vec));
        PetscCall(DMDAGlobalToNaturalEnd(dm, field_vec, INSERT_VALUES, natural_vec));
    }

    PetscCall(VecScatterCreateToZero(natural_vec ? natural_vec : field_vec, &scatter, sequential_vec));
    PetscCall(VecScatterBegin(scatter, natural_vec ? natural_vec : field_vec, *sequential_vec,
                              INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, natural_vec ? natural_vec : field_vec, *sequential_vec,
                            INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterDestroy(&scatter));
    PetscCall(VecDestroy(&natural_vec));
    PetscFunctionReturn(0);
}

/** @brief Compute and cache a rank-count-independent hash of the active grid geometry. */
static PetscErrorCode ComputeCheckpointGeometrySHA256(SimCtx *simCtx, UserCtx *user, char digest_hex[65])
{
    PicurvSHA256Context hash_context;
    char metadata[256];

    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL && user != NULL && digest_hex != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Geometry hash inputs cannot be NULL.");
    if (simCtx->checkpointGeometryHashReady) {
        PetscCall(PetscStrncpy(digest_hex, simCtx->checkpointGeometrySHA256, 65));
        PetscFunctionReturn(0);
    }

    if (simCtx->rank == 0) PicurvSHA256Init(&hash_context);
    for (PetscInt block = 0; block < simCtx->block_number; ++block) {
        FieldView coordinates;
        Vec sequential_vec = NULL;

        PetscCall(FieldGetView(&user[block], FIELD_ID_COORDINATES, &coordinates));
        PetscCall(GatherVectorToRankZero(coordinates.global_vec, &sequential_vec));
        if (simCtx->rank == 0) {
            const PetscScalar *values = NULL;
            PetscInt value_count = 0;
            size_t metadata_length = 0;

            PetscCall(PetscSNPrintf(metadata, sizeof(metadata),
                                    "block=%" PetscInt_FMT ";im=%" PetscInt_FMT ";jm=%" PetscInt_FMT
                                    ";km=%" PetscInt_FMT ";periodic=%d,%d,%d;",
                                    block, user[block].IM, user[block].JM, user[block].KM,
                                    (int)simCtx->i_periodic, (int)simCtx->j_periodic, (int)simCtx->k_periodic));
            PetscCall(PetscStrlen(metadata, &metadata_length));
            PicurvSHA256Update(&hash_context, metadata, metadata_length);
            PetscCall(VecGetSize(sequential_vec, &value_count));
            PicurvSHA256Update(&hash_context, &value_count, sizeof(value_count));
            PetscCall(VecGetArrayRead(sequential_vec, &values));
            PicurvSHA256Update(&hash_context, values, (size_t)value_count * sizeof(*values));
            PetscCall(VecRestoreArrayRead(sequential_vec, &values));
        }
        PetscCall(VecDestroy(&sequential_vec));
    }

    if (simCtx->rank == 0) PicurvSHA256FinalHex(&hash_context, digest_hex);
    PetscCallMPI(MPI_Bcast(digest_hex, 65, MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx->checkpointGeometrySHA256, digest_hex, 65));
    simCtx->checkpointGeometryHashReady = PETSC_TRUE;
    PetscFunctionReturn(0);
}

/** @brief Format the canonical directory name for one completed step. */
static PetscErrorCode FormatCheckpointStepDirectory(const char *root, PetscInt step, char *path, size_t path_size)
{
    PetscFunctionBeginUser;
    PetscCheck(root != NULL && path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Checkpoint root and output path are required.");
    PetscCall(PetscSNPrintf(path, path_size, "%s/%s/step_%0*" PetscInt_FMT,
                            root, PICURV_CHECKPOINTS_DIRECTORY, PICURV_CHECKPOINT_STEP_WIDTH, step));
    PetscFunctionReturn(0);
}

/** @brief Resolve either an exact bundle or a run/output root to one step bundle. */
static PetscErrorCode ResolveCheckpointStepDirectory(const char *source_root, PetscInt step,
                                                     char *path, size_t path_size)
{
    char metadata_path[PETSC_MAX_PATH_LEN];
    PetscBool exact_bundle = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCheck(source_root != NULL && path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Checkpoint source and output path are required.");
    PetscCall(PetscSNPrintf(metadata_path, sizeof(metadata_path), "%s/checkpoint.meta", source_root));
    PetscCall(PetscTestFile(metadata_path, 'r', &exact_bundle));
    if (exact_bundle) PetscCall(PetscStrncpy(path, source_root, path_size));
    else PetscCall(FormatCheckpointStepDirectory(source_root, step, path, path_size));
    PetscFunctionReturn(0);
}

/** @brief Create one directory on rank zero and report failures collectively. */
static PetscErrorCode CreateCheckpointDirectoryCollective(const SimCtx *simCtx, const char *path)
{
    PetscMPIInt create_failed = 0;

    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL && path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Simulation context and directory path are required.");
    if (simCtx->rank == 0 && mkdir(path, 0777) != 0) {
        struct stat directory_stat;
        if (errno != EEXIST || stat(path, &directory_stat) != 0 || !S_ISDIR(directory_stat.st_mode)) {
            create_failed = 1;
        }
    }
    PetscCallMPI(MPI_Bcast(&create_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD));
    PetscCheck(!create_failed, PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
               "Unable to create checkpoint directory '%s'.", path);
    PetscFunctionReturn(0);
}

/** @brief Validate a committed bundle and return selected authoritative metadata. */
static PetscErrorCode ValidateCheckpointBundle(SimCtx *simCtx, UserCtx *user,
                                               const char *checkpoint_directory,
                                               PetscInt expected_step,
                                               PetscReal *physical_time,
                                               PetscInt *particle_count,
                                               PetscBool *particles_saved,
                                               PetscBool *les_saved,
                                               PetscBool *rans_saved)
{
    char metadata_path[PETSC_MAX_PATH_LEN];
    char commit_path[PETSC_MAX_PATH_LEN];
    char expected_digest[65] = "";
    char actual_digest[65] = "";
    char format[64] = "";
    char saved_geometry_digest[65] = "";
    char current_geometry_digest[65] = "";
    PetscOptions options = NULL;
    PetscInt version = 0;
    PetscInt saved_step = -1;
    PetscInt payload_count = 0;
    PetscInt saved_particle_count = 0;
    PetscReal saved_time = 0.0;
    PetscBool saved_particles = PETSC_FALSE;
    PetscBool saved_les = PETSC_FALSE;
    PetscBool saved_rans = PETSC_FALSE;
    PetscBool found = PETSC_FALSE;
    FILE *commit_file = NULL;

    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL && user != NULL && checkpoint_directory != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Checkpoint validation inputs cannot be NULL.");
    PetscCall(PetscSNPrintf(metadata_path, sizeof(metadata_path), "%s/checkpoint.meta", checkpoint_directory));
    PetscCall(PetscSNPrintf(commit_path, sizeof(commit_path), "%s/COMMITTED", checkpoint_directory));

    PetscCall(PetscTestFile(metadata_path, 'r', &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
               "Checkpoint metadata is missing: %s", metadata_path);
    PetscCall(PetscTestFile(commit_path, 'r', &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
               "Checkpoint is not committed: %s", checkpoint_directory);

    commit_file = fopen(commit_path, "r");
    PetscCheck(commit_file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
               "Unable to read checkpoint commit marker '%s'.", commit_path);
    PetscCheck(fgets(expected_digest, sizeof(expected_digest), commit_file) != NULL,
               PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
               "Checkpoint commit marker '%s' is empty.", commit_path);
    PetscCheck(fclose(commit_file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
               "Unable to close checkpoint commit marker '%s'.", commit_path);
    TrimWhitespace(expected_digest);
    PetscCheck(strlen(expected_digest) == 64, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint commit marker '%s' does not contain a SHA-256 digest.", commit_path);
    PetscCall(PicurvSHA256File(metadata_path, actual_digest));
    PetscCheck(!strcmp(expected_digest, actual_digest), PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint metadata hash mismatch in '%s'.", checkpoint_directory);

    PetscCall(PetscOptionsCreate(&options));
    PetscCall(PetscOptionsInsertFile(PETSC_COMM_WORLD, options, metadata_path, PETSC_TRUE));
    PetscCall(PetscOptionsGetString(options, NULL, "-checkpoint_format", format, sizeof(format), &found));
    PetscCheck(found && !strcmp(format, PICURV_CHECKPOINT_FORMAT), PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Unsupported checkpoint format '%s' in '%s'.", found ? format : "<missing>", metadata_path);
    PetscCall(PetscOptionsGetInt(options, NULL, "-checkpoint_version", &version, &found));
    PetscCheck(found && version == PICURV_CHECKPOINT_VERSION, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Unsupported checkpoint version %" PetscInt_FMT " in '%s'; expected %d.",
               version, metadata_path, PICURV_CHECKPOINT_VERSION);
    PetscCall(PetscOptionsGetInt(options, NULL, "-checkpoint_step", &saved_step, &found));
    PetscCheck(found && saved_step == expected_step, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint step is %" PetscInt_FMT ", expected %" PetscInt_FMT ".",
               saved_step, expected_step);
    PetscCall(PetscOptionsGetReal(options, NULL, "-checkpoint_time", &saved_time, &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' does not record physical time.", metadata_path);
    PetscCall(PetscOptionsGetString(options, NULL, "-checkpoint_geometry_sha256",
                                    saved_geometry_digest, sizeof(saved_geometry_digest), &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' does not record geometry identity.", metadata_path);
    PetscCall(ComputeCheckpointGeometrySHA256(simCtx, user, current_geometry_digest));
    PetscCheck(!strcmp(saved_geometry_digest, current_geometry_digest), PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
               "Checkpoint geometry/layout does not match the active grid.");

    PetscCall(PetscOptionsGetInt(options, NULL, "-checkpoint_payload_count", &payload_count, &found));
    PetscCheck(found && payload_count > 0, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' has no payload inventory.", metadata_path);
    PetscCall(PetscOptionsGetInt(options, NULL, "-checkpoint_particle_count",
                                 &saved_particle_count, &found));
    PetscCheck(found && saved_particle_count >= 0, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' has no valid particle count.", metadata_path);
    PetscCall(PetscOptionsGetBool(options, NULL, "-checkpoint_particles", &saved_particles, &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' does not record particle-state availability.", metadata_path);
    PetscCall(PetscOptionsGetBool(options, NULL, "-checkpoint_les", &saved_les, &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' does not record LES-state availability.", metadata_path);
    PetscCall(PetscOptionsGetBool(options, NULL, "-checkpoint_rans", &saved_rans, &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint '%s' does not record RANS-state availability.", metadata_path);
    for (PetscInt payload = 0; payload < payload_count; ++payload) {
        char option_name[96];
        char relative_path[PETSC_MAX_PATH_LEN];
        char payload_path[PETSC_MAX_PATH_LEN];
        long long expected_bytes = -1;
        char expected_bytes_text[64];
        struct stat payload_stat;

        PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                                "-checkpoint_payload_%" PetscInt_FMT "_path", payload));
        PetscCall(PetscOptionsGetString(options, NULL, option_name,
                                        relative_path, sizeof(relative_path), &found));
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
                   "Checkpoint payload %" PetscInt_FMT " has no path.", payload);
        PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                                "-checkpoint_payload_%" PetscInt_FMT "_bytes", payload));
        PetscCall(PetscOptionsGetString(options, NULL, option_name,
                                        expected_bytes_text, sizeof(expected_bytes_text), &found));
        if (found) expected_bytes = strtoll(expected_bytes_text, NULL, 10);
        PetscCheck(found && expected_bytes >= 0, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
                   "Checkpoint payload %" PetscInt_FMT " has no valid byte size.", payload);
        PetscCall(PetscSNPrintf(payload_path, sizeof(payload_path), "%s/%s",
                                checkpoint_directory, relative_path));
        PetscCheck(stat(payload_path, &payload_stat) == 0 && S_ISREG(payload_stat.st_mode),
                   PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                   "Checkpoint payload is missing: %s", payload_path);
        PetscCheck((long long)payload_stat.st_size == expected_bytes,
                   PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
                   "Checkpoint payload '%s' is %lld bytes; expected %lld.",
                   payload_path, (long long)payload_stat.st_size, expected_bytes);
    }
    PetscCall(PetscOptionsDestroy(&options));
    if (physical_time) *physical_time = saved_time;
    if (particle_count) *particle_count = saved_particle_count;
    if (particles_saved) *particles_saved = saved_particles;
    if (les_saved) *les_saved = saved_les;
    if (rans_saved) *rans_saved = saved_rans;
    PetscFunctionReturn(0);
}


// =============================================================================
//                      PUBLIC FUNCTION IMPLEMENTATIONS
// =============================================================================

/**
 * @brief Implementation of \ref TrimWhitespace().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see TrimWhitespace()
 */
void TrimWhitespace(char *str) {
    if (!str) return;
    if (str[0] == '\0') return;

    char *start = str;
    // Find the first non-whitespace character
    while (isspace((unsigned char)*start)) {
        start++;
    }

    // Find the end of the string
    char *end = str + strlen(str) - 1;
    // Move backwards from the end to find the last non-whitespace character
    while (end > start && isspace((unsigned char)*end)) {
        end--;
    }

    // Null-terminate after the last non-whitespace character
    *(end + 1) = '\0';

    // If there was leading whitespace, shift the string to the left
    if (str != start) {
        memmove(str, start, (end - start) + 2); // +2 to include the new null terminator
    }
}

/**
 * @brief Implementation of \ref ShouldWriteDataOutput().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see ShouldWriteDataOutput()
 */

PetscBool ShouldWriteDataOutput(const SimCtx *simCtx, PetscInt completed_step)
{
    if (!simCtx) {
        return PETSC_FALSE;
    }
    return (PetscBool)(simCtx->tiout > 0 && completed_step > 0 && completed_step % simCtx->tiout == 0);
}


#undef __FUNCT__
#define __FUNCT__ "ReadGridGenerationInputs"
/**
 * @brief Internal helper implementation: `ReadGridGenerationInputs()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadGridGenerationInputs(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    PetscInt       nblk = simCtx->block_number;
    PetscInt       block_index = user->_this;
    PetscBool      found;

    // Temporary arrays to hold the parsed values for ALL blocks
    PetscInt  *IMs = NULL, *JMs = NULL, *KMs = NULL, *cgrids = NULL;
    PetscReal *xMins = NULL, *xMaxs = NULL, *rxs = NULL;
    PetscReal *yMins = NULL, *yMaxs = NULL, *rys = NULL;
    PetscReal *zMins = NULL, *zMaxs = NULL, *rzs = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reading generated grid inputs for block %d.\n", simCtx->rank, block_index);

    if (block_index >= nblk) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Block index %d is out of range for nblk=%d", block_index, nblk);
    }

    // --- Allocate temporary storage for all array options ---
    ierr = PetscMalloc4(nblk, &IMs, nblk, &JMs, nblk, &KMs, nblk, &cgrids); CHKERRQ(ierr);
    ierr = PetscMalloc6(nblk, &xMins, nblk, &xMaxs, nblk, &rxs, nblk, &yMins, nblk, &yMaxs, nblk, &rys); CHKERRQ(ierr);
    ierr = PetscMalloc3(nblk, &zMins, nblk, &zMaxs, nblk, &rzs); CHKERRQ(ierr);

    // --- Set default values for the temporary arrays ---
    for (PetscInt i = 0; i < nblk; ++i) {
        IMs[i] = 10; JMs[i] = 10; KMs[i] = 10; cgrids[i] = 0;
        xMins[i] = 0.0; xMaxs[i] = 1.0; rxs[i] = 1.0;
        yMins[i] = 0.0; yMaxs[i] = 1.0; rys[i] = 1.0;
        zMins[i] = 0.0; zMaxs[i] = 1.0; rzs[i] = 1.0;
    }

    // --- Parse the array options from the command line / control file ---
    PetscInt count;
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", IMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", JMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", KMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-xMins", xMins, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-xMaxs", xMaxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-rxs", rxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-yMins", yMins, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-yMaxs", yMaxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-rys", rys, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-zMins", zMins, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-zMaxs", zMaxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-rzs", rzs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-cgrids", cgrids, &count, &found); CHKERRQ(ierr);

    // --- Assign the parsed values to the specific UserCtx struct passed in ---
    user->IM = IMs[block_index];
    user->JM = JMs[block_index];
    user->KM = KMs[block_index];
    user->Min_X = xMins[block_index];
    user->Max_X = xMaxs[block_index];
    user->rx = rxs[block_index];
    user->Min_Y = yMins[block_index];
    user->Max_Y = yMaxs[block_index];
    user->ry = rys[block_index];
    user->Min_Z = zMins[block_index];
    user->Max_Z = zMaxs[block_index];
    user->rz = rzs[block_index];
    user->cgrid = cgrids[block_index];

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d grid generation inputs set: IM=%d, JM=%d, KM=%d\n",
              simCtx->rank, block_index, user->IM, user->JM, user->KM);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d bounds: X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f]\n",
              simCtx->rank, block_index, user->Min_X, user->Max_X, user->Min_Y, user->Max_Y, user->Min_Z, user->Max_Z);

    // --- Clean up temporary storage ---
    ierr = PetscFree4(IMs, JMs, KMs, cgrids); CHKERRQ(ierr);
    ierr = PetscFree6(xMins, xMaxs, rxs, yMins, yMaxs, rys); CHKERRQ(ierr);
    ierr = PetscFree3(zMins, zMaxs, rzs); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `PopulateFinestUserGridResolutionFromOptions()`.
 * @details Local to this translation unit.
 */
PetscErrorCode PopulateFinestUserGridResolutionFromOptions(UserCtx *finest_users, PetscInt nblk)
{
    PetscErrorCode ierr;
    PetscBool      found;
    PetscInt       *IMs = NULL, *JMs = NULL, *KMs = NULL;
    SimCtx         *simCtx = NULL;

    PetscFunctionBeginUser;

    if (!finest_users) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "finest_users cannot be NULL.");
    }
    if (nblk <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "nblk must be positive. Got %d.", nblk);
    }
    simCtx = finest_users[0].simCtx;

    ierr = PetscMalloc3(nblk, &IMs, nblk, &JMs, nblk, &KMs); CHKERRQ(ierr);
    for (PetscInt i = 0; i < nblk; ++i) {
        IMs[i] = 10; JMs[i] = 10; KMs[i] = 10;
    }

    PetscInt count;
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", IMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", JMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", KMs, &count, &found); CHKERRQ(ierr);

    for (PetscInt bi = 0; bi < nblk; ++bi) {
        finest_users[bi].IM = IMs[bi];
        finest_users[bi].JM = JMs[bi];
        finest_users[bi].KM = KMs[bi];
        if (simCtx) {
            LOG_ALLOW(LOCAL, LOG_DEBUG,
                      "Rank %d: Preloaded analytical grid resolution for block %d: IM=%d, JM=%d, KM=%d\n",
                      simCtx->rank, bi, finest_users[bi].IM, finest_users[bi].JM, finest_users[bi].KM);
        }
    }

    ierr = PetscFree3(IMs, JMs, KMs); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ReadGridFile"
/**
 * @brief Internal helper implementation: `ReadGridFile()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadGridFile(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    PetscInt block_index = user->_this;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- One-Time Read and Broadcast Logic ---
    if (!g_file_has_been_read) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "First call to ReadGridFile. Reading and broadcasting grid file header from '%s'...\n", simCtx->grid_file);
        PetscMPIInt rank = simCtx->rank;
        PetscInt    nblk = simCtx->block_number;

        if (rank == 0) {
            FILE *fd = fopen(simCtx->grid_file, "r");
            if (!fd) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", simCtx->grid_file);

            // Read and validate the canonical PICGRID header.
            char firstTok[32] = {0};
            if (fscanf(fd, "%31s", firstTok) != 1)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Empty grid file: %s", simCtx->grid_file);
            if (strcmp(firstTok, "PICGRID") != 0)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                        "Grid file %s must begin with the canonical PICGRID header.", simCtx->grid_file);
            if (fscanf(fd, "%d", &g_nblk_from_file) != 1)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Expected number of blocks after \"PICGRID\" in %s", simCtx->grid_file);
            if (g_nblk_from_file != nblk) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED, "Mismatch: -nblk is %d but grid file specifies %d blocks.", nblk, g_nblk_from_file);
            }

            ierr = PetscMalloc3(nblk, &g_IMs_from_file, nblk, &g_JMs_from_file, nblk, &g_KMs_from_file); CHKERRQ(ierr);
            for (PetscInt i = 0; i < nblk; ++i) {
                if (fscanf(fd, "%d %d %d\n", &g_IMs_from_file[i], &g_JMs_from_file[i], &g_KMs_from_file[i]) != 3) {
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Expected 3 integers for block %d in %s", i, simCtx->grid_file);
                }
            }
            fclose(fd);
        }

        // Broadcast nblk to verify (optional, good practice)
        ierr = MPI_Bcast(&g_nblk_from_file, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        // Allocate on other ranks before receiving the broadcast
        if (rank != 0) {
            ierr = PetscMalloc3(nblk, &g_IMs_from_file, nblk, &g_JMs_from_file, nblk, &g_KMs_from_file); CHKERRQ(ierr);
        }
        
        // Broadcast the data arrays
        ierr = MPI_Bcast(g_IMs_from_file, nblk, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(g_JMs_from_file, nblk, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(g_KMs_from_file, nblk, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        g_file_has_been_read = PETSC_TRUE;
        LOG_ALLOW(GLOBAL, LOG_INFO, "Grid file header read and broadcast complete.\n");
    }

    // --- Per-Block Assignment Logic (runs on every call) ---
    user->IM = g_IMs_from_file[block_index];
    user->JM = g_JMs_from_file[block_index];
    user->KM = g_KMs_from_file[block_index];

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Set file inputs for Block %d: IM=%d, JM=%d, KM=%d\n",
              simCtx->rank, block_index, user->IM, user->JM, user->KM);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


//================================================================================
//
//                        PRIVATE HELPER FUNCTIONS
//
//================================================================================

/**
 * @brief Implementation of \ref FreeBC_ParamList().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see FreeBC_ParamList()
 */
void FreeBC_ParamList(BC_Param *head) {
    BC_Param *current = head;
    while (current != NULL) {
        BC_Param *next = current->next;
        PetscFree(current->key);
        PetscFree(current->value);
        PetscFree(current);
        current = next;
    }
}

/**
 * @brief Internal helper implementation: `StringToBCFace()`.
 * @details Local to this translation unit.
 */
PetscErrorCode StringToBCFace(const char* str, BCFace* face_out) {
    if      (strcasecmp(str, "-Xi")   == 0) *face_out = BC_FACE_NEG_X;
    else if (strcasecmp(str, "+Xi")   == 0) *face_out = BC_FACE_POS_X;
    else if (strcasecmp(str, "-Eta")  == 0) *face_out = BC_FACE_NEG_Y;
    else if (strcasecmp(str, "+Eta")  == 0) *face_out = BC_FACE_POS_Y;
    else if (strcasecmp(str, "-Zeta") == 0) *face_out = BC_FACE_NEG_Z;
    else if (strcasecmp(str, "+Zeta") == 0) *face_out = BC_FACE_POS_Z;
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown face specifier: %s", str);
    return 0;
}

/**
 * @brief Internal helper implementation: `StringToBCType()`.
 * @details Local to this translation unit.
 */
PetscErrorCode StringToBCType(const char* str, BCType* type_out) {
    if      (strcasecmp(str, "WALL")      == 0) *type_out = WALL;
    else if (strcasecmp(str, "SYMMETRY")  == 0) *type_out = SYMMETRY;
    else if (strcasecmp(str, "INLET")     == 0) *type_out = INLET;
    else if (strcasecmp(str, "OUTLET")    == 0) *type_out = OUTLET;
    else if (strcasecmp(str, "PERIODIC")  == 0) *type_out = PERIODIC;
    // ... add other BCTypes here ...
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown BC Type string: %s", str);
    return 0;
}

/**
 * @brief Internal helper implementation: `StringToBCHandlerType()`.
 * @details Local to this translation unit.
 */
PetscErrorCode StringToBCHandlerType(const char* str, BCHandlerType* handler_out) {
    if      (strcasecmp(str, "noslip")              == 0) *handler_out = BC_HANDLER_WALL_NOSLIP;
    else if (strcasecmp(str, "constant_velocity")   == 0) *handler_out = BC_HANDLER_INLET_CONSTANT_VELOCITY;
    else if (strcasecmp(str, "conservation")        == 0) *handler_out = BC_HANDLER_OUTLET_CONSERVATION;
    else if (strcasecmp(str, "parabolic")           == 0) *handler_out = BC_HANDLER_INLET_PARABOLIC;
    else if (strcasecmp(str, "prescribed_flow")     == 0) *handler_out = BC_HANDLER_INLET_PROFILE_FROM_FILE;
    else if (strcasecmp(str,"geometric")            == 0) *handler_out = BC_HANDLER_PERIODIC_GEOMETRIC;
    else if (strcasecmp(str,"constant_flux")        == 0) *handler_out = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    else if (strcasecmp(str,"initial_flux")         == 0) *handler_out = BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX;
    // ... add other BCHandlerTypes here ...
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown BC Handler string: %s", str);
    return 0;
}

/**
 * @brief Internal helper implementation: `ValidateBCHandlerForBCType()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ValidateBCHandlerForBCType(BCType type, BCHandlerType handler) {
    switch (type) {
        case OUTLET:
	    if(handler != BC_HANDLER_OUTLET_CONSERVATION) return PETSC_ERR_ARG_WRONG;
	    break;
        case WALL:
            if (handler != BC_HANDLER_WALL_NOSLIP && handler != BC_HANDLER_WALL_MOVING) return PETSC_ERR_ARG_WRONG;
            break;
        case INLET:
            if (handler != BC_HANDLER_INLET_CONSTANT_VELOCITY &&
                handler != BC_HANDLER_INLET_PARABOLIC &&
                handler != BC_HANDLER_INLET_PROFILE_FROM_FILE) return PETSC_ERR_ARG_WRONG;
            break;
        case PERIODIC:
            if(handler != BC_HANDLER_PERIODIC_GEOMETRIC && handler != BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX && handler != BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX) return PETSC_ERR_ARG_WRONG;
        // ... add other validation cases here ...
        default: break;
    }
    return 0; // Combination is valid
}

/**
 * @brief Internal helper implementation: `GetBCParamReal()`.
 * @details Local to this translation unit.
 */
PetscErrorCode GetBCParamReal(BC_Param *params, const char *key, PetscReal *value_out, PetscBool *found) {
    *found = PETSC_FALSE;
    *value_out = 0.0;
    if (!key) return 0; // No key to search for

    BC_Param *current = params;
    while (current) {
        if (strcasecmp(current->key, key) == 0) {
            *value_out = atof(current->value);
            *found = PETSC_TRUE;
            return 0; // Found it, we're done
        }
        current = current->next;
    }
    return 0; // It's not an error to not find the key.
}

/**
 * @brief Internal helper implementation: `GetBCParamBool()`.
 * @details Local to this translation unit.
 */
PetscErrorCode GetBCParamBool(BC_Param *params, const char *key, PetscBool *value_out, PetscBool *found) {
    *found = PETSC_FALSE;
    *value_out = PETSC_FALSE;
    if (!key) return 0; // No key to search for

    BC_Param *current = params;
    while (current) {
        if (strcasecmp(current->key, key) == 0) {
                        // Key was found.
            *found = PETSC_TRUE;
            
            // Check the value string. Default to FALSE if the value is NULL or doesn't match a "true" string.
            if (current->value && 
               (strcasecmp(current->value, "true") == 0 ||
                strcmp(current->value, "1") == 0         ||
                strcasecmp(current->value, "yes") == 0))
            {
                *value_out = PETSC_TRUE;
            } else {
                *value_out = PETSC_FALSE;
            }
            return 0; // Found it, we're done
        }
        current = current->next;
    }
    return 0; // It's not an error to not find the key.
}

#undef __FUNCT__
#define __FUNCT__ "GetDrivenSeamFluxFlag"
/**
 * @brief Implementation of \ref GetDrivenSeamFluxFlag().
 * @details The option was originally spelled `apply_trim`, which said that
 *          something was trimmed but not what or why. The canonical name is now
 *          `enforce_seam_flux`. Generated `bcs.run` files carry the canonical
 *          name, but a hand-written or archived one may still use the old
 *          spelling, so both are accepted here and the canonical name wins.
 *          The argument contract lives with the header declaration in
 *          `include/io.h`.
 * @see GetDrivenSeamFluxFlag()
 */
PetscErrorCode GetDrivenSeamFluxFlag(BC_Param *params, PetscBool *value_out, PetscBool *found)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    ierr = GetBCParamBool(params, "enforce_seam_flux", value_out, found); CHKERRQ(ierr);
    if (!*found) {
        ierr = GetBCParamBool(params, "apply_trim", value_out, found); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

//================================================================================
//
//                        PUBLIC PARSING FUNCTION
//
//================================================================================
#undef __FUNCT__
#define __FUNCT__ "ParseAllBoundaryConditions"
/**
 * @brief Internal helper implementation: `ParseAllBoundaryConditions()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ParseAllBoundaryConditions(UserCtx *user, const char *bcs_input_filename)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    // Temporary storage for rank 0 to build the configuration before broadcasting.
    BoundaryFaceConfig configs_rank0[6];
    
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        FILE *file;
        char line_buffer[1024];

        // Initialize the temporary config array with safe defaults on rank 0.
        for (int i = 0; i < 6; i++) {
            configs_rank0[i].face_id = (BCFace)i;
            configs_rank0[i].mathematical_type = WALL;
            configs_rank0[i].handler_type = BC_HANDLER_WALL_NOSLIP;
            configs_rank0[i].params = NULL;
            configs_rank0[i].handler = NULL; // Handler object is not created here.
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "Parsing BC configuration from '%s' on rank 0... \n", bcs_input_filename);
        file = fopen(bcs_input_filename, "r");
        if (!file) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Could not open BCs file '%s'.", bcs_input_filename);

        while (fgets(line_buffer, sizeof(line_buffer), file)) {
            char *current_pos = line_buffer;
            while (isspace((unsigned char)*current_pos)) current_pos++; // Skip leading whitespace
            if (*current_pos == '#' || *current_pos == '\0' || *current_pos == '\n' || *current_pos == '\r') continue;

            char *face_str = strtok(current_pos, " \t\n\r");
            char *type_str = strtok(NULL, " \t\n\r");
            char *handler_str = strtok(NULL, " \t\n\r");

            if (!face_str || !type_str || !handler_str) {
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Malformed line in bcs.dat, skipping: %s \n", line_buffer);
                continue;
            }

            BCFace      face_enum;
            BCType      type_enum;
            BCHandlerType handler_enum;
            const char* handler_name_for_log;

            // --- Convert strings to enums and validate ---
            ierr = StringToBCFace(face_str, &face_enum); CHKERRQ(ierr);
            ierr = StringToBCType(type_str, &type_enum); CHKERRQ(ierr);
            ierr = StringToBCHandlerType(handler_str, &handler_enum); CHKERRQ(ierr);
            ierr = ValidateBCHandlerForBCType(type_enum, handler_enum);
            if (ierr) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Validation failed: Handler '%s' is not valid for Type '%s' on Face '%s'.\n", handler_str, type_str, face_str);

            // Store the core types for the corresponding face
            configs_rank0[face_enum].mathematical_type = type_enum;
            configs_rank0[face_enum].handler_type = handler_enum;
            handler_name_for_log = BCHandlerTypeToString(handler_enum); // Assumes this utility exists
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "  Parsed Face '%s': Type=%s, Handler=%s \n", face_str, type_str, handler_name_for_log);

            // --- Parse optional key=value parameters for this face ---
            FreeBC_ParamList(configs_rank0[face_enum].params); // Clear any previous (default) params
            configs_rank0[face_enum].params = NULL;
            BC_Param **param_next_ptr = &configs_rank0[face_enum].params; // Pointer to the 'next' pointer to build the list

            char* token;
            while ((token = strtok(NULL, " \t\n\r")) != NULL) {
                char* equals_ptr = strchr(token, '=');
                if (!equals_ptr) {
                    LOG_ALLOW(GLOBAL, LOG_WARNING, "Malformed parameter '%s' on face '%s', skipping. \n", token, face_str);
                    continue;
                }
                
                *equals_ptr = '\0'; // Temporarily split the string at '=' to separate key and value
                char* key_str = token;
                char* value_str = equals_ptr + 1;

                BC_Param *new_param;
                ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                ierr = PetscStrallocpy(key_str, &new_param->key); CHKERRQ(ierr);
                ierr = PetscStrallocpy(value_str, &new_param->value); CHKERRQ(ierr);
                new_param->next = NULL;

                *param_next_ptr = new_param;
                param_next_ptr = &new_param->next;
                LOG_ALLOW(GLOBAL, LOG_TRACE, "    - Found param: [%s] = [%s] \n", new_param->key, new_param->value);
            }
        }
        fclose(file);
    }

    // =========================================================================
    //               BROADCASTING THE CONFIGURATION FROM RANK 0
    // =========================================================================
    // This is a critical step to ensure all processes have the same configuration.
    
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d broadcasting/receiving BC configuration.\n", rank);

    for (int i = 0; i < 6; i++) {
        // --- Broadcast simple enums ---
        if (rank == 0) {
            user->boundary_faces[i] = configs_rank0[i]; // Rank 0 populates its final struct
        }
        ierr = MPI_Bcast(&user->boundary_faces[i].mathematical_type, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(&user->boundary_faces[i].handler_type, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // --- Serialize and Broadcast the parameter linked list ---
        PetscInt n_params = 0;
        if (rank == 0) { // On rank 0, count the number of parameters to send
            for (BC_Param *p = user->boundary_faces[i].params; p; p = p->next) n_params++;
        }
        ierr = MPI_Bcast(&n_params, 1, MPI_INT, 0, PETSC_COMM_WORLD);CHKERRQ(ierr);
        
        if (rank != 0) { // Non-root ranks need to receive and build the list
            FreeBC_ParamList(user->boundary_faces[i].params); // Ensure list is empty before building
            user->boundary_faces[i].params = NULL;
        }

        BC_Param **param_next_ptr = &user->boundary_faces[i].params;

        for (int j = 0; j < n_params; j++) {
            char key_buf[256] = {0}, val_buf[256] = {0};
            if (rank == 0) {
                // On rank 0, navigate to the j-th param and copy its data to buffers
                BC_Param *p = user->boundary_faces[i].params;
                for (int k = 0; k < j; k++) p = p->next;
                strncpy(key_buf, p->key, 255);
                strncpy(val_buf, p->value, 255);
            }

            ierr = MPI_Bcast(key_buf, 256, MPI_CHAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Bcast(val_buf, 256, MPI_CHAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

            if (rank != 0) {
                // On non-root ranks, deserialize: create a new node and append it
                BC_Param *new_param;
                ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                ierr = PetscStrallocpy(key_buf, &new_param->key); CHKERRQ(ierr);
                ierr = PetscStrallocpy(val_buf, &new_param->value); CHKERRQ(ierr);
                new_param->next = NULL;
                *param_next_ptr = new_param;
                param_next_ptr = &new_param->next;
            } else {
                 // On rank 0, just advance the pointer for the next iteration
                 param_next_ptr = &((*param_next_ptr)->next);
            }
        }
        user->boundary_faces[i].face_id = (BCFace)i; // Ensure face_id is set on all ranks
    }
    
    // --- Set particle inlet lookup fields used by the particle system ---
    user->inletFaceDefined = PETSC_FALSE;
    for (int i=0; i<6; i++) {
        
        if (user->boundary_faces[i].mathematical_type == INLET && !user->inletFaceDefined) {
            user->inletFaceDefined = PETSC_TRUE;
            user->identifiedInletBCFace = (BCFace)i;
            LOG_ALLOW(GLOBAL, LOG_INFO, "Inlet face for particle initialization identified as Face %d.\n", i);
            break; // Found the first one, stop looking
        }
    }

    
    if (rank == 0) {
        // Rank 0 can now free the linked lists it created for the temporary storage.
        // As written, user->boundary_faces was populated directly on rank 0, so no extra free is needed.
        // for(int i=0; i<6; i++) FreeBC_ParamList(configs_rank0[i].params); // This would be needed if we used configs_rank0 exclusively
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//                        PRIVATE HELPER FUNCTIONS
//
//================================================================================

// ... (existing helper functions like FreeBC_ParamList, StringToBCFace, etc.) ...
#undef __FUNCT__
#define __FUNCT__ "DeterminePeriodicity"
/**
 * @brief Internal helper implementation: `DeterminePeriodicity()`.
 * @details Local to this translation unit.
 */
PetscErrorCode DeterminePeriodicity(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       periodic_flags[3] = {0, 0, 0}; // Index 0:I, 1:J, 2:K

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- Part 1: Collectively verify all BCS files exist before proceeding ---
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        const char *bcs_filename = simCtx->bcs_files[bi];
        if (!bcs_filename) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "BCS filename for block %d is not set in SimCtx.", bi);
        char desc_buf[256];
        PetscBool file_exists;
        snprintf(desc_buf, sizeof(desc_buf), "BCS file for block %d", bi);
        ierr = VerifyPathExistence(bcs_filename, PETSC_FALSE, PETSC_FALSE, desc_buf, &file_exists); CHKERRQ(ierr);
    }

    // --- Part 2: Rank 0 does the parsing, since we know all files exist ---
    if (rank == 0) {
        PetscBool global_is_periodic[3] = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
        PetscBool is_set = PETSC_FALSE;

        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            const char *bcs_filename = simCtx->bcs_files[bi];
            FILE *file = fopen(bcs_filename, "r");

            PetscBool face_is_periodic[6] = {PETSC_FALSE};
            char line_buffer[1024];

            while (fgets(line_buffer, sizeof(line_buffer), file)) {
                char *current_pos = line_buffer;
                while (isspace((unsigned char)*current_pos)) current_pos++;
                if (*current_pos == '#' || *current_pos == '\0' || *current_pos == '\n') continue;

                // --- Tokenize the line exactly like the main parser ---
                char *face_str = strtok(current_pos, " \t\n\r");
                char *type_str = strtok(NULL, " \t\n\r");

                // If the line doesn't have at least two tokens, we can't determine the type.
                if (!face_str || !type_str) continue;

                // --- Perform a direct, non-erroring check on the mathematical type string ---
                if (strcasecmp(type_str, "PERIODIC") == 0) {
                    BCFace face_enum;
                    // A malformed face string on a periodic line IS a fatal error.
                    ierr = StringToBCFace(face_str, &face_enum); CHKERRQ(ierr);
                    face_is_periodic[face_enum] = PETSC_TRUE;
                }
                // Any other type_str (e.g., "WALL", "INLET") is correctly and silently ignored.
            }
            fclose(file);

            // --- Validate consistency within this file ---
            if (face_is_periodic[BC_FACE_NEG_X] != face_is_periodic[BC_FACE_POS_X])
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Inconsistent X-periodicity in file '%s' for block %d. Both -Xi and +Xi must be periodic or neither.", bcs_filename, bi);
            if (face_is_periodic[BC_FACE_NEG_Y] != face_is_periodic[BC_FACE_POS_Y])
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Inconsistent Y-periodicity in file '%s' for block %d. Both -Eta and +Eta must be periodic or neither.", bcs_filename, bi);
            if (face_is_periodic[BC_FACE_NEG_Z] != face_is_periodic[BC_FACE_POS_Z])
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Inconsistent Z-periodicity in file '%s' for block %d. Both -Zeta and +Zeta must be periodic or neither.", bcs_filename, bi);

            PetscBool local_is_periodic[3] = {face_is_periodic[BC_FACE_NEG_X], face_is_periodic[BC_FACE_NEG_Y], face_is_periodic[BC_FACE_NEG_Z]};

            // --- Validate consistency across block files ---
            if (!is_set) {
                global_is_periodic[0] = local_is_periodic[0];
                global_is_periodic[1] = local_is_periodic[1];
                global_is_periodic[2] = local_is_periodic[2];
                is_set = PETSC_TRUE;
            } else {
                if (global_is_periodic[0] != local_is_periodic[0] || global_is_periodic[1] != local_is_periodic[1] || global_is_periodic[2] != local_is_periodic[2]) {
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP,
                            "Periodicity mismatch between blocks. Block 0 requires (I:%d, J:%d, K:%d), but block %d (file '%s') has (I:%d, J:%d, K:%d).",
                            (int)global_is_periodic[0], (int)global_is_periodic[1], (int)global_is_periodic[2],
                            bi, bcs_filename,
                            (int)local_is_periodic[0], (int)local_is_periodic[1], (int)local_is_periodic[2]);
                }
            }
        } // end loop over blocks

        periodic_flags[0] = (global_is_periodic[0]) ? 1 : 0;
        periodic_flags[1] = (global_is_periodic[1]) ? 1 : 0;
        periodic_flags[2] = (global_is_periodic[2]) ? 1 : 0;

        LOG_ALLOW(GLOBAL, LOG_INFO, "Global periodicity determined: I-periodic=%d, J-periodic=%d, K-periodic=%d\n",
                  periodic_flags[0], periodic_flags[1], periodic_flags[2]);
    }

    // --- Part 3: Broadcast the final flags from rank 0 to all other ranks ---
    ierr = MPI_Bcast(periodic_flags, 3, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    // --- All ranks now update their SimCtx ---
    simCtx->i_periodic = periodic_flags[0];
    simCtx->j_periodic = periodic_flags[1];
    simCtx->k_periodic = periodic_flags[2];

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `VerifyPathExistence()`.
 * @details Local to this translation unit.
 */
PetscErrorCode VerifyPathExistence(const char *path, PetscBool is_dir, PetscBool is_optional, const char *description, PetscBool *exists)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    MPI_Comm       comm = PETSC_COMM_WORLD;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        if (is_dir) {
            ierr = PetscTestDirectory(path, 'r', exists); CHKERRQ(ierr);
        } else {
            ierr = PetscTestFile(path, 'r', exists); CHKERRQ(ierr);
        }

        if (!(*exists)) {
            if (is_optional) {
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Optional %s not found at: %s (using defaults/ignoring).\n", description, path);
            } else {
                LOG_ALLOW(GLOBAL, LOG_ERROR, "Mandatory %s not found at: %s\n", description, path);
            }
        } else {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "Found %s: %s\n", description, path);
        }
    }

    // Broadcast the result from Rank 0
    PetscMPIInt exists_int = (rank == 0) ? (PetscMPIInt)(*exists) : 0;
    ierr = MPI_Bcast(&exists_int, 1, MPI_INT, 0, comm); CHKERRMPI(ierr);
    *exists = (PetscBool)exists_int;

    // Collective error for mandatory files
    if (!(*exists) && !is_optional) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Mandatory %s not found. Rank 0 expected it at '%s'. Check path and permissions.", description, path);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReadFieldData"
/**
 * @brief Internal helper implementation: `ReadFieldData()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadFieldData(UserCtx *user,
                             const char *field_name,
                             Vec         field_vec,
                             const char *ext)
{
   PetscErrorCode ierr;
   char           filename[PETSC_MAX_PATH_LEN];
   MPI_Comm       comm;
   PetscMPIInt    rank,size;
   SimCtx         *simCtx = user->simCtx;


   PetscFunctionBeginUser;
   PROFILE_FUNCTION_BEGIN;

    if(!simCtx->current_io_directory){
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "I/O context directory was not set before calling ReadFieldData().");
    }


   ierr = PetscObjectGetComm((PetscObject)field_vec,&comm);CHKERRQ(ierr);
   ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
   ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

   const char *source_path = NULL;
   source_path = simCtx->current_io_directory;

   if(!source_path){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "source_path was not set for the current execution mode.");
   }
   ierr = PetscSNPrintf(filename, sizeof(filename), "%s/%s.%s",
                        source_path, field_name, ext); CHKERRQ(ierr);

   LOG_ALLOW(GLOBAL,LOG_DEBUG,
             "Attempting to read <%s> on rank %d/%d\n",
             filename,(int)rank,(int)size);

   /* ======================================================================
    * 1.  SERIAL JOB  – just hand the Vec to VecLoad()
    * ==================================================================== */
   if(size==1)
   {
      PetscViewer viewer;
      PetscBool   found;
      Vec temp_vec;
      PetscInt expectedSize,loadedSize;

      ierr = PetscTestFile(filename,'r',&found);CHKERRQ(ierr);
      if(!found) SETERRQ(comm,PETSC_ERR_FILE_OPEN,
                         "Restart/Source file not found: %s",filename);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
// ---- START MODIFICATION ----
      // DO NOT load directly into field_vec, as this can resize it, which is
      // illegal for DMSwarm "view" vectors. Instead, load into a temporary vector.
      ierr = VecCreate(PETSC_COMM_SELF, &temp_vec); CHKERRQ(ierr);
      ierr = VecLoad(temp_vec,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

      // Sanity check: ensure the file size matches the expected vector size.
      ierr = VecGetSize(field_vec, &expectedSize);CHKERRQ(ierr);
      ierr = VecGetSize(temp_vec, &loadedSize);CHKERRQ(ierr);
      if (loadedSize != expectedSize) {
         SETERRQ(comm,PETSC_ERR_FILE_UNEXPECTED,
                 "File %s holds %d entries – expected %d for field '%s'",
                 filename, loadedSize, expectedSize, field_name);
      }

      // Now, safely copy the data from the temporary vector to the final destination.
      ierr = VecCopy(temp_vec, field_vec);CHKERRQ(ierr);

      // Clean up the temporary vector.
      ierr = VecDestroy(&temp_vec);CHKERRQ(ierr);

      // ---- END MODIFICATION ----
      
      /* create EMPTY sequential Vec – VecLoad() will size it correctly   */
      /*
      ierr = VecCreate(PETSC_COMM_SELF,&seq_vec);CHKERRQ(ierr);
      ierr = VecSetType(seq_vec,VECSEQ);CHKERRQ(ierr);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,
                                   FILE_MODE_READ,&viewer);CHKERRQ(ierr);

      ierr = VecLoad(field_vec,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
      */
      LOG_ALLOW(GLOBAL,LOG_INFO,
                "Loaded <%s> (serial path)\n",filename);

      PROFILE_FUNCTION_END;          
      PetscFunctionReturn(0);
   }

   /* ======================================================================
    * 2.  PARALLEL JOB
    * ==================================================================== */
   PetscInt globalSize;
   ierr = VecGetSize(field_vec,&globalSize);CHKERRQ(ierr);

   DM         dm = NULL;
   const char *dmtype = NULL;
   Vec        nat = NULL;                 /* Natural-ordered vector for DMDA */

   /* -------------------- rank-0 : read the sequential file -------------- */
   Vec            seq_vec = NULL;      /* only valid on rank-0            */
   const PetscScalar *seqArray = NULL; /* borrowed pointer on rank-0 only */

   if(rank==0)
   {
      PetscViewer viewer;
      PetscBool   found;

      ierr = PetscTestFile(filename,'r',&found);CHKERRQ(ierr);
      if(!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,
                         "Restart file not found: %s",filename);

      /* create EMPTY sequential Vec – VecLoad() will size it correctly   */
      ierr = VecCreate(PETSC_COMM_SELF,&seq_vec);CHKERRQ(ierr);
      ierr = VecSetType(seq_vec,VECSEQ);CHKERRQ(ierr);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,
                                   FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(seq_vec,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

      /* size sanity-check */
      PetscInt loaded;
      ierr = VecGetSize(seq_vec,&loaded);CHKERRQ(ierr);
      if(loaded != globalSize)
         SETERRQ(comm,PETSC_ERR_FILE_UNEXPECTED,
                 "File %s holds %d entries – expected %d",
                 filename,loaded,globalSize);

      /* borrow array for later Bcast */
      ierr = VecGetArrayRead(seq_vec,&seqArray);CHKERRQ(ierr);

      LOG_ALLOW(GLOBAL,LOG_TRACE,
                "Rank 0 successfully loaded <%s>\n",filename);
   }

   /* -------------------- Check if this is a DMDA vector ----------------- */
   ierr = VecGetDM(field_vec, &dm); CHKERRQ(ierr);
   if (dm) { ierr = DMGetType(dm, &dmtype); CHKERRQ(ierr); }

   if (dmtype && !strcmp(dmtype, DMDA)) {
      /* ==================================================================
       * DMDA PATH: File is in natural ordering, need to convert to global
       * ================================================================== */

      /* Create natural vector */
      ierr = DMDACreateNaturalVector(dm, &nat); CHKERRQ(ierr);

      /* Scatter from rank 0's seq_vec to all ranks' natural vector */
      VecScatter scatter;
      Vec nat_seq = NULL;  /* Sequential natural vector on rank 0 */

      ierr = VecScatterCreateToZero(nat, &scatter, &nat_seq); CHKERRQ(ierr);

      /* Reverse scatter: from rank 0 to all ranks */
      ierr = VecScatterBegin(scatter, (rank == 0 ? seq_vec : nat_seq), nat,
                             INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecScatterEnd(scatter, (rank == 0 ? seq_vec : nat_seq), nat,
                           INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

      /* Convert natural → global ordering */
      ierr = DMDANaturalToGlobalBegin(dm, nat, INSERT_VALUES, field_vec); CHKERRQ(ierr);
      ierr = DMDANaturalToGlobalEnd(dm, nat, INSERT_VALUES, field_vec); CHKERRQ(ierr);

      /* Cleanup */
      ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
      ierr = VecDestroy(&nat_seq); CHKERRQ(ierr);
      ierr = VecDestroy(&nat); CHKERRQ(ierr);

   } else {
      /* ==================================================================
       * NON-DMDA PATH: Use broadcast and direct copy (assumes global ordering)
       * ================================================================== */

      PetscScalar *buffer = NULL;
      if (rank == 0) {
         buffer = (PetscScalar *)seqArray;
      } else {
         ierr = PetscMalloc1(globalSize, &buffer); CHKERRQ(ierr);
      }

      ierr = MPI_Bcast(buffer, (int)globalSize, MPIU_SCALAR, 0, comm); CHKERRQ(ierr);

      /* Copy slice based on ownership range */
      PetscInt rstart, rend, loc;
      PetscScalar *locArray;

      ierr = VecGetOwnershipRange(field_vec, &rstart, &rend); CHKERRQ(ierr);
      loc = rend - rstart;

      ierr = VecGetArray(field_vec, &locArray); CHKERRQ(ierr);
      ierr = PetscMemcpy(locArray, buffer + rstart, loc * sizeof(PetscScalar)); CHKERRQ(ierr);
      ierr = VecRestoreArray(field_vec, &locArray); CHKERRQ(ierr);

      if (rank != 0) {
         ierr = PetscFree(buffer); CHKERRQ(ierr);
      }
   }

   /* -------------------- tidy up ---------------------------------------- */
   if (rank == 0) {
      ierr = VecRestoreArrayRead(seq_vec, &seqArray); CHKERRQ(ierr);
      ierr = VecDestroy(&seq_vec); CHKERRQ(ierr);
   }

   LOG_ALLOW(GLOBAL,LOG_INFO,
             "Loaded <%s> (parallel path)\n",filename);

   PROFILE_FUNCTION_END;          
   PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "RestoreDrivenFluxTarget"
/**
 * @brief Restore a latched driven-flow flux target from a checkpoint manifest.
 *
 * @details Only `initial_flux` needs this. Its target is derived from the field
 *          the run originally started with, so re-measuring it after a restart
 *          would silently retarget the controller at whatever the flux had
 *          drifted to. `constant_flux` reads its target from the bcs file on
 *          every run and is deliberately left alone, so that editing
 *          `target_flux` between segments still takes effect.
 *
 *          Checkpoints written before this metadata existed simply have no
 *          entry; the controller then falls back to re-measuring and says so.
 *
 * @param[in,out] simCtx               Simulation context receiving the target.
 * @param[in]     user                 Block context, for its boundary handler types.
 * @param[in]     checkpoint_directory Directory holding the validated manifest.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode RestoreDrivenFluxTarget(SimCtx *simCtx, UserCtx *user,
                                              const char *checkpoint_directory)
{
    char         metadata_path[PETSC_MAX_PATH_LEN];
    PetscOptions options = NULL;
    PetscBool    uses_initial_flux = PETSC_FALSE;
    PetscBool    saved_latched = PETSC_FALSE;
    PetscBool    found = PETSC_FALSE;
    PetscReal    saved_target = 0.0;

    PetscFunctionBeginUser;

    for (PetscInt face = 0; face < 6; ++face) {
        if (user->boundary_faces[face].handler_type == BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX) {
            uses_initial_flux = PETSC_TRUE;
            break;
        }
    }
    if (!uses_initial_flux) PetscFunctionReturn(0);

    PetscCall(PetscSNPrintf(metadata_path, sizeof(metadata_path),
                            "%s/checkpoint.meta", checkpoint_directory));
    PetscCall(PetscOptionsCreate(&options));
    PetscCall(PetscOptionsInsertFile(PETSC_COMM_WORLD, options, metadata_path, PETSC_TRUE));
    PetscCall(PetscOptionsGetBool(options, NULL, "-checkpoint_driven_flux_latched", &saved_latched, &found));
    if (found && saved_latched) {
        PetscCall(PetscOptionsGetReal(options, NULL, "-checkpoint_driven_flux_target", &saved_target, &found));
    } else {
        found = PETSC_FALSE;
    }
    PetscCall(PetscOptionsDestroy(&options));

    if (found) {
        simCtx->targetVolumetricFlux    = saved_target;
        simCtx->drivenFluxTargetLatched = PETSC_TRUE;
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Driven Flow: restored latched target volumetric flux %.6e from checkpoint '%s'.\n",
                  (double)saved_target, checkpoint_directory);
    } else {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Driven Flow: checkpoint '%s' records no latched flux target; the initial_flux "
                  "controller will re-measure it from the restarted field.\n", checkpoint_directory);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ReadSimulationFields()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadSimulationFields(UserCtx *user,PetscInt ti)
{
    SimCtx *simCtx = user->simCtx;
    const char *source_path = NULL;
    char checkpoint_directory[PETSC_MAX_PATH_LEN];
    PetscReal checkpoint_time = 0.0;
    PetscBool particles_saved = PETSC_FALSE;
    PetscBool les_saved = PETSC_FALSE;
    PetscBool rans_saved = PETSC_FALSE;

    PetscFunctionBeginUser;
    if(simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR){
        source_path = simCtx->pps->source_dir;
    } else if(simCtx->exec_mode == EXEC_MODE_SOLVER){
        source_path = simCtx->restart_dir;
    } else{
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Invalid execution mode for reading simulation fields.");
    }

    PetscCall(ResolveCheckpointStepDirectory(source_path, ti,
                                             checkpoint_directory, sizeof(checkpoint_directory)));
    PetscCall(ValidateCheckpointBundle(simCtx, user - user->_this,
                                       checkpoint_directory, ti, &checkpoint_time, NULL,
                                       &particles_saved, &les_saved, &rans_saved));
    PetscCheck(!(simCtx->np > 0 && !strcmp(simCtx->particleRestartMode, "load")) || particles_saved,
               PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Particle restart_mode=load was requested, but checkpoint step %" PetscInt_FMT
               " contains no particle state.", ti);
    if (simCtx->exec_mode == EXEC_MODE_SOLVER && user->_this == 0) {
        simCtx->ti = checkpoint_time;
        if (ti == simCtx->StartStep) simCtx->StartTime = checkpoint_time;
        PetscCall(RestoreDrivenFluxTarget(simCtx, user, checkpoint_directory));
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Reading Eulerian checkpoint fields for block %d from '%s'.\n",
              user->_this, checkpoint_directory);
    PetscCall(PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                            "%s/%s/block_%04" PetscInt_FMT,
                            checkpoint_directory, PICURV_EULERIAN_DIRECTORY, user->_this));
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    for (PetscInt raw_id = 0; raw_id < FIELD_ID_COUNT; ++raw_id) {
        const FieldDescriptor *descriptor = NULL;
        FieldView view;

        PetscCall(FieldGetDescriptor((FieldId)raw_id, &descriptor));
        if (!CheckpointFieldIsEnabled(simCtx, descriptor)) continue;
        if ((descriptor->availability & FIELD_AVAILABILITY_PARTICLES) &&
            (!particles_saved || strcmp(simCtx->particleRestartMode, "load"))) continue;
        if ((descriptor->availability & FIELD_AVAILABILITY_LES_DYNAMIC) && !les_saved) continue;
        if ((descriptor->availability & FIELD_AVAILABILITY_RANS) && !rans_saved) continue;
        if ((descriptor->availability & FIELD_AVAILABILITY_TURBULENCE) &&
            !((simCtx->les && les_saved) || (simCtx->rans && rans_saved))) continue;
        PetscCall(FieldGetView(user, descriptor->id, &view));
        PetscCall(ReadFieldData(user, descriptor->canonical_name, view.global_vec, "dat"));
        if (view.local_vec) PetscCall(UpdateLocalGhosts(user, descriptor->id));
    }
    if (simCtx->rans) {
        PetscCall(VecCopy(user->K_Omega, user->K_Omega_o));
        PetscCall(UpdateLocalGhosts(user, FIELD_ID_K_OMEGA_O));
    }
    simCtx->restartHistoryAvailable = PETSC_TRUE;
    simCtx->current_io_directory = NULL;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ReadStatisticsWindowState"
/** @brief Restore one window's scalar bookkeeping from a validated manifest. */
static PetscErrorCode ReadStatisticsWindowState(PetscOptions options, PetscInt window,
                                                const char *metadata_path, PetscReal checkpoint_time,
                                                PetscReal step_size, ExecutionMode exec_mode,
                                                PicurvWindow *state)
{
    char option_name[128];
    char saved_name[PICURV_WINDOW_NAME_LENGTH] = "";
    char saved_digest[65] = "";
    char saved_groups[PICURV_WINDOW_HASH_GROUP_COUNT * (PICURV_WINDOW_HASH_GROUP_LENGTH + 4)] = "";
    char saved_state[32] = "";
    char current_digest[65] = "";
    PetscBool found = PETSC_FALSE;
    PetscBool name_matches = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

#define PICURV_STATISTICS_REQUIRE(suffix, getter, target)                                        \
    do {                                                                                         \
        PetscCall(PetscSNPrintf(option_name, sizeof(option_name),                                \
                                "-checkpoint_statistics_window_%" PetscInt_FMT "_" suffix, window)); \
        PetscCall(getter(options, NULL, option_name, target, &found));                           \
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,                           \
                   "Statistics continuation requested, but '%s' records no %s for window "       \
                   "%" PetscInt_FMT ".", metadata_path, suffix, window);                         \
    } while (0)

    PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                            "-checkpoint_statistics_window_%" PetscInt_FMT "_name", window));
    PetscCall(PetscOptionsGetString(options, NULL, option_name, saved_name, sizeof(saved_name), &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Statistics continuation requested, but '%s' records no name for window %" PetscInt_FMT ".",
               metadata_path, window);
    PetscCall(PetscStrcmp(saved_name, state->definition.name, &name_matches));
    PetscCheck(name_matches, PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
               "Statistics window %" PetscInt_FMT " is '%s' in the checkpoint but '%s' in this run. "
               "Window order and names must match to continue.",
               window, saved_name, state->definition.name);

    PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                            "-checkpoint_statistics_window_%" PetscInt_FMT "_hash", window));
    PetscCall(PetscOptionsGetString(options, NULL, option_name, saved_digest, sizeof(saved_digest), &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Statistics continuation requested, but '%s' records no hash for window '%s'.",
               metadata_path, state->definition.name);
    PetscCall(PicurvWindowComputeHash(&state->definition, current_digest, NULL));
    if (strcmp(saved_digest, current_digest)) {
        const char *differing = "an unrecorded property";

        /* The saved definition itself is not stored, so the per-group digests are
         * what make it possible to name the property that changed rather than
         * reporting only that two opaque hashes differ. */
        PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                                "-checkpoint_statistics_window_%" PetscInt_FMT "_hash_groups", window));
        PetscCall(PetscOptionsGetString(options, NULL, option_name, saved_groups,
                                        sizeof(saved_groups), &found));
        if (found) {
            PetscInt group = -1;

            PetscCall(PicurvWindowFirstHashDifference(&state->definition, saved_groups, &group));
            if (group >= 0) differing = PicurvWindowHashGroupName(group);
        }
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                "Statistics window '%s' was defined differently in the checkpoint: %s changed. "
                "Continuing would merge incompatible samples; rename the window to start a new one.",
                state->definition.name, differing);
    }

    PICURV_STATISTICS_REQUIRE("sample_count", PetscOptionsGetInt, &state->sample_count);
    PICURV_STATISTICS_REQUIRE("total_weight", PetscOptionsGetReal, &state->total_weight);
    PICURV_STATISTICS_REQUIRE("represented_time", PetscOptionsGetReal, &state->represented_time);
    PICURV_STATISTICS_REQUIRE("last_accepted_time", PetscOptionsGetReal, &state->last_accepted_time);
    PICURV_STATISTICS_REQUIRE("effective_start", PetscOptionsGetReal, &state->effective_start);
    PICURV_STATISTICS_REQUIRE("effective_end", PetscOptionsGetReal, &state->effective_end);
    PICURV_STATISTICS_REQUIRE("activation_step", PetscOptionsGetInt, &state->activation_step);
    PICURV_STATISTICS_REQUIRE("last_event_step", PetscOptionsGetInt, &state->last_event_step);
    PICURV_STATISTICS_REQUIRE("next_time_target", PetscOptionsGetInt, &state->next_time_target);
    PICURV_STATISTICS_REQUIRE("restart_count", PetscOptionsGetInt, &state->restart_count);
#undef PICURV_STATISTICS_REQUIRE

    PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                            "-checkpoint_statistics_window_%" PetscInt_FMT "_state", window));
    PetscCall(PetscOptionsGetString(options, NULL, option_name, saved_state, sizeof(saved_state), &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Statistics continuation requested, but '%s' records no lifecycle state for window '%s'.",
               metadata_path, state->definition.name);
    if (!strcmp(saved_state, "complete")) state->state = PICURV_WINDOW_COMPLETE;
    else if (!strcmp(saved_state, "active")) state->state = PICURV_WINDOW_ACTIVE;
    else state->state = PICURV_WINDOW_PENDING;

    /* A window saved as complete may still be resumable: page 58 §10 permits moving
     * a bounded end forward, and the end time is deliberately outside the hash. */
    if (state->state == PICURV_WINDOW_COMPLETE && state->definition.bounded &&
        state->definition.end_time > state->effective_end) {
        /* Only across an unbroken span. Under right-rectangle weighting the first
         * sample after the former end carries the whole interval back to it, so
         * resuming across a gap would weight time the window never observed. One
         * step of slack is the closing state's own clipping, not a gap. */
        PetscCheck(checkpoint_time - state->effective_end <= PetscMax(step_size, 0.0),
                   PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                   "Statistics window '%s' ended at t=%.17g and this checkpoint is at t=%.17g, "
                   "so extending it to %.17g would weight %.17g of unobserved time. "
                   "Start a new window instead.",
                   state->definition.name, (double)state->effective_end, (double)checkpoint_time,
                   (double)state->definition.end_time,
                   (double)(checkpoint_time - state->effective_end));
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Statistics window '%s' was complete at t=%.6g and is extended to t=%.6g.\n",
                  state->definition.name, (double)state->effective_end,
                  (double)state->definition.end_time);
        state->state = PICURV_WINDOW_ACTIVE;
    }
    PetscCheck(!state->definition.bounded || state->definition.end_time >= state->effective_end,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
               "Statistics window '%s' already represents time up to %.17g, past the requested end "
               "%.17g. A window may be extended forward but never shortened; rename it to start a new one.",
               state->definition.name, (double)state->effective_end,
               (double)state->definition.end_time);
    /* A post-processing read is analysis, not a continuation, so it must not look
     * like another restart segment. Only a solver resuming the window advances the
     * lineage; otherwise deriving a series of steps would inflate it once per step. */
    if (exec_mode == EXEC_MODE_SOLVER) state->restart_count += 1;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RestoreFieldStatisticsState"
/**
 * @brief Implementation of \ref RestoreFieldStatisticsState().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see RestoreFieldStatisticsState()
 */
PetscErrorCode RestoreFieldStatisticsState(SimCtx *simCtx, PetscInt ti)
{
    UserCtx *user = NULL;
    PetscOptions options = NULL;
    char checkpoint_directory[PETSC_MAX_PATH_LEN];
    char metadata_path[PETSC_MAX_PATH_LEN];
    const char *source_path = NULL;
    PetscInt saved_window_count = 0;
    PetscReal checkpoint_time = 0.0;
    PetscBool found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");
    if (!FieldStatisticsIsActive(simCtx)) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    if (!simCtx->fieldStatisticsContinue) {
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Statistics continuation was not requested; %d window(s) start from zero.\n",
                  simCtx->fieldStatisticsWindowCount);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    PetscCheck(user != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Finest-level fields must exist before statistics state is restored.");
    source_path = (simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR) ? simCtx->pps->source_dir
                                                                 : simCtx->restart_dir;
    PetscCheck(source_path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "No checkpoint source directory is set for the current execution mode.");
    PetscCall(ResolveCheckpointStepDirectory(source_path, ti,
                                             checkpoint_directory, sizeof(checkpoint_directory)));
    PetscCall(ValidateCheckpointBundle(simCtx, user, checkpoint_directory, ti,
                                       &checkpoint_time, NULL, NULL, NULL, NULL));
    PetscCall(PetscSNPrintf(metadata_path, sizeof(metadata_path), "%s/checkpoint.meta",
                            checkpoint_directory));

    PetscCall(PetscOptionsCreate(&options));
    PetscCall(PetscOptionsInsertFile(PETSC_COMM_WORLD, options, metadata_path, PETSC_TRUE));
    PetscCall(PetscOptionsGetInt(options, NULL, "-checkpoint_statistics_window_count",
                                 &saved_window_count, &found));
    /* Missing state for a requested continuation is fatal and never silently
     * zeroed: resuming from zero would report a converged average built from a
     * fraction of the samples its metadata claims. */
    PetscCheck(found && saved_window_count == simCtx->fieldStatisticsWindowCount,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
               "Statistics continuation requested, but checkpoint '%s' holds %" PetscInt_FMT
               " window(s) and this run configures %" PetscInt_FMT ".",
               checkpoint_directory, found ? saved_window_count : 0,
               simCtx->fieldStatisticsWindowCount);

    for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
        PetscCall(ReadStatisticsWindowState(options, window, metadata_path, checkpoint_time,
                                            simCtx->dt, simCtx->exec_mode,
                                            &simCtx->fieldStatisticsWindows[window]));
    }
    PetscCall(PetscOptionsDestroy(&options));

    for (PetscInt block = 0; block < simCtx->block_number; ++block) {
        PetscCheck(user[block].fieldStatisticsStorage != NULL, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                   "Statistics accumulators were not allocated for block %" PetscInt_FMT ".", block);
        for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
            const PicurvWindowDefinition *definition = &simCtx->fieldStatisticsWindows[window].definition;
            const PicurvWindowStorage *storage = &user[block].fieldStatisticsStorage[window];
            PetscInt payload_count = 0;

            PetscCall(FormatStatisticsPath(checkpoint_directory, window, block, NULL,
                                           simCtx->_io_context_buffer,
                                           sizeof(simCtx->_io_context_buffer)));
            simCtx->current_io_directory = simCtx->_io_context_buffer;
            PetscCall(PicurvWindowStoragePayloadCount(storage, &payload_count));
            for (PetscInt index = 0; index < payload_count; ++index) {
                PicurvStatisticsPayload payload;

                PetscCall(PicurvWindowStoragePayload(&user[block], definition, storage, index, &payload));
                PetscCall(ReadFieldData(&user[block], payload.name, payload.vec, "dat"));
            }
        }
        simCtx->current_io_directory = NULL;
    }

    for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
        const PicurvWindow *state = &simCtx->fieldStatisticsWindows[window];

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Statistics window '%s' continued from '%s': state %s, %d sample(s), "
                  "total weight %.6g, represented time %.6g, restart segment %d.\n",
                  state->definition.name, checkpoint_directory,
                  PicurvWindowStateName(state->state), state->sample_count,
                  (double)state->total_weight, (double)state->represented_time,
                  state->restart_count);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ReadSwarmField()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadSwarmField(UserCtx *user, const char *field_name, const char *ext)
{
  PetscErrorCode ierr;
  DM             swarm;
  Vec            fieldVec;

  PetscFunctionBegin;

  swarm = user->swarm;

  LOG_ALLOW(GLOBAL,LOG_DEBUG," ReadSwarmField Begins \n");
 
  /* 2) Create a global vector that references the specified Swarm field. */
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," Vector created from Field \n");

  /* 3) Use the ReadFieldData() function to read data into fieldVec. */
  ierr = ReadFieldData(user, field_name, fieldVec, ext);CHKERRQ(ierr);

  /* 4) Destroy the global vector reference. */
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ReadSwarmIntField()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadSwarmIntField(UserCtx *user, const char *field_name, const char *ext)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    Vec            temp_vec;
    PetscInt       nlocal, nglobal, bs, i;
    PetscDataType  field_type;
    const PetscScalar *scalar_array; // Read-only pointer from the temp Vec
    void           *field_array_void;


    PetscFunctionBeginUser;
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Reading '%s' via temporary Vec.\n", field_name);

    // Get the properties of the swarm field to determine the expected layout
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &nglobal); CHKERRQ(ierr);
    // We get the block size but not the data pointer yet
    ierr = DMSwarmGetField(swarm, field_name, &bs, &field_type, NULL); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, field_name, &bs, NULL, NULL); CHKERRQ(ierr);
    PetscCheck(field_type == PETSC_INT || field_type == PETSC_INT64,
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Swarm field '%s' must use PETSC_INT or PETSC_INT64 for integer restart input, not %s.",
               field_name, PetscDataTypes[field_type]);

    // Create a temporary Vec with the CORRECT layout to receive the data
    ierr = VecCreate(PETSC_COMM_WORLD, &temp_vec); CHKERRQ(ierr);
    ierr = VecSetType(temp_vec, VECMPI); CHKERRQ(ierr);
    ierr = VecSetSizes(temp_vec, nlocal * bs, nglobal * bs); CHKERRQ(ierr);
    ierr = VecSetBlockSize(temp_vec, bs); CHKERRQ(ierr);
    ierr = VecSetUp(temp_vec); CHKERRQ(ierr);

    // Call your existing reader to populate the temporary Vec
    ierr = ReadFieldData(user, field_name, temp_vec, ext); CHKERRQ(ierr);

    // Get local pointers
    ierr = VecGetArrayRead(temp_vec, &scalar_array); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, NULL, NULL, &field_array_void); CHKERRQ(ierr);
    
    // Perform the cast back, using the correct loop size (nlocal * bs)
    if (field_type == PETSC_INT64) {
        PetscInt64 *int64_array = (PetscInt64 *)field_array_void;
        for (i = 0; i < nlocal * bs; i++) {
            int64_array[i] = (PetscInt64)scalar_array[i];
        }
    } else {
        PetscInt *int_array = (PetscInt *)field_array_void;
        for (i = 0; i < nlocal * bs; i++) {
            int_array[i] = (PetscInt)scalar_array[i];
        }
    }    

    // Restore access
    ierr = DMSwarmRestoreField(swarm, field_name, NULL, NULL, &field_array_void); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(temp_vec, &scalar_array); CHKERRQ(ierr);

    // 6. Clean up
    ierr = VecDestroy(&temp_vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ReadAllSwarmFields()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadAllSwarmFields(UserCtx *user, PetscInt ti)
{
  PetscInt nGlobal;
  SimCtx *simCtx = user->simCtx;
  const char *source_path = NULL;
  char checkpoint_directory[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  PetscCall(DMSwarmGetSize(user->swarm, &nGlobal));
  LOG_ALLOW(GLOBAL, LOG_INFO, "Reading DMSwarm fields for timestep %d (swarm size is %d).\n", ti, nGlobal);

  if (nGlobal == 0) {
      LOG_ALLOW(GLOBAL, LOG_INFO, "Swarm is empty for timestep %d. Nothing to read.\n", ti);
      PetscFunctionReturn(0);
  }

    // First, determine the top-level source directory based on the execution mode.
    if (simCtx->exec_mode == EXEC_MODE_SOLVER) {
        source_path = simCtx->restart_dir;
    } else if (simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR) {
        source_path = simCtx->pps->source_dir;
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Invalid execution mode for reading simulation fields.");
    }

    PetscCall(ResolveCheckpointStepDirectory(source_path, ti,
                                             checkpoint_directory, sizeof(checkpoint_directory)));
    PetscCall(PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                            "%s/%s", checkpoint_directory, PICURV_PARTICLE_DIRECTORY));
    simCtx->current_io_directory = simCtx->_io_context_buffer;

  for (PetscInt raw_id = 0; raw_id < PARTICLE_FIELD_ID_COUNT; ++raw_id) {
      const ParticleFieldDescriptor *descriptor = NULL;

      PetscCall(ParticleFieldGetDescriptor((ParticleFieldId)raw_id, &descriptor));
      if (!(descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_CHECKPOINT)) continue;
      if (descriptor->data_type == PETSC_INT || descriptor->data_type == PETSC_INT64) {
          PetscCall(ReadSwarmIntField(user, descriptor->canonical_name, "dat"));
      } else {
          PetscCall(ReadSwarmField(user, descriptor->canonical_name, "dat"));
      }
  }

  simCtx->current_io_directory = NULL;

  LOG_ALLOW(GLOBAL, LOG_INFO, "Finished reading DMSwarm fields for timestep %d.\n", ti);
  PetscFunctionReturn(0);
}

/** @brief Implementation of \ref ReadCheckpointParticleCount(). */
PetscErrorCode ReadCheckpointParticleCount(UserCtx *user, PetscInt ti, PetscInt *particle_count)
{
    SimCtx *simCtx = NULL;
    const char *source_path = NULL;
    char checkpoint_directory[PETSC_MAX_PATH_LEN];
    PetscBool particles_saved = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && particle_count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Simulation context and particle-count output are required.");
    simCtx = user->simCtx;
    if (simCtx->exec_mode == EXEC_MODE_SOLVER) source_path = simCtx->restart_dir;
    else if (simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR) source_path = simCtx->pps->source_dir;
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                 "Invalid execution mode for reading checkpoint particle metadata.");

    PetscCall(ResolveCheckpointStepDirectory(source_path, ti,
                                             checkpoint_directory, sizeof(checkpoint_directory)));
    PetscCall(ValidateCheckpointBundle(simCtx, user - user->_this,
                                       checkpoint_directory, ti, NULL, particle_count,
                                       &particles_saved, NULL, NULL));
    PetscCheck(particles_saved, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
               "Checkpoint step %" PetscInt_FMT " contains no particle state.", ti);
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "WriteFieldData" 
/**
 * @brief Internal helper implementation: `WriteFieldData()`.
 * @details Local to this translation unit.
 */
PetscErrorCode WriteFieldData(UserCtx *user,
                              const char *field_name,
                              Vec field_vec,
                              const char *ext)
{
    MPI_Comm       comm;
    PetscMPIInt    rank;
    Vec            sequential_vec = NULL;
    char           filename[PETSC_MAX_PATH_LEN];
    SimCtx         *simCtx=user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if(!simCtx->current_io_directory){
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "I/O context directory was not set before calling WriteFieldData().");
    }

    /* ------------------------------------------------------------ */
    /*                  Basic communicator information              */
    /* ------------------------------------------------------------ */
    PetscCall(PetscObjectGetComm((PetscObject)field_vec,&comm));
    PetscCallMPI(MPI_Comm_rank(comm,&rank));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "%s/%s.%s",
                            simCtx->current_io_directory, field_name, ext));

    PetscCall(GatherVectorToRankZero(field_vec, &sequential_vec));
    if (rank == 0) {
        PetscViewer viewer;
        PetscReal vmin, vmax;

        PetscCall(VecMin(sequential_vec, NULL, &vmin));
        PetscCall(VecMax(sequential_vec, NULL, &vmax));
        LOG_ALLOW(GLOBAL,LOG_DEBUG,
                  " <%s> range = [%.4e … %.4e]\n",
                  field_name,(double)vmin,(double)vmax);
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer));
        PetscCall(PetscViewerBinarySetSkipInfo(viewer, PETSC_TRUE));
        PetscCall(VecView(sequential_vec, viewer));
        PetscCall(PetscViewerDestroy(&viewer));
        LOG_ALLOW(GLOBAL, LOG_INFO, "Wrote <%s>\n", filename);
    }
    PetscCall(VecDestroy(&sequential_vec));

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref WriteSimulationFields().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see WriteSimulationFields()
 */
PetscErrorCode WriteSimulationFields(UserCtx *user, const char *checkpoint_directory)
{
    SimCtx *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PetscCheck(checkpoint_directory != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Checkpoint destination cannot be NULL.");
    PetscCall(PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                            "%s/%s/block_%04" PetscInt_FMT,
                            checkpoint_directory, PICURV_EULERIAN_DIRECTORY, user->_this));
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    if (simCtx->les) {
        // The constant model holds no coefficient field to stage.
        if (simCtx->les == DYNAMIC_SMAGORINSKY) {
            PetscCall(CopyOwnedLocalScalarToGlobal(user->da, user->lCs, user->CS));
        }
        PetscCall(CopyOwnedLocalScalarToGlobal(user->da, user->lNu_t, user->Nu_t));
    }
    for (PetscInt raw_id = 0; raw_id < FIELD_ID_COUNT; ++raw_id) {
        const FieldDescriptor *descriptor = NULL;
        FieldView view;

        PetscCall(FieldGetDescriptor((FieldId)raw_id, &descriptor));
        if (!CheckpointFieldIsEnabled(simCtx, descriptor)) continue;
        PetscCall(FieldGetView(user, descriptor->id, &view));
        PetscCall(WriteFieldData(user, descriptor->canonical_name, view.global_vec, "dat"));
    }
    simCtx->current_io_directory = NULL;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WriteStatisticsFields"
/** @brief Write every window's accumulator payloads for one block into a bundle. */
static PetscErrorCode WriteStatisticsFields(UserCtx *user, const char *checkpoint_directory)
{
    SimCtx *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(checkpoint_directory != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Checkpoint destination cannot be NULL.");
    if (!FieldStatisticsIsActive(simCtx) || !user->fieldStatisticsStorage) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }

    for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
        const PicurvWindowStorage *storage = &user->fieldStatisticsStorage[window];
        PetscInt payload_count = 0;

        PetscCall(FormatStatisticsPath(checkpoint_directory, window, user->_this, NULL,
                                       simCtx->_io_context_buffer,
                                       sizeof(simCtx->_io_context_buffer)));
        simCtx->current_io_directory = simCtx->_io_context_buffer;
        PetscCall(PicurvWindowStoragePayloadCount(storage, &payload_count));
        for (PetscInt index = 0; index < payload_count; ++index) {
            PicurvStatisticsPayload payload;

            PetscCall(PicurvWindowStoragePayload(user, &simCtx->fieldStatisticsWindows[window].definition,
                                                 storage, index, &payload));
            PetscCall(WriteFieldData(user, payload.name, payload.vec, "dat"));
        }
        LOG_ALLOW(LOCAL, LOG_DEBUG,
                  "Wrote %d statistics payload(s) for window '%s' block %d.\n",
                  (int)payload_count, simCtx->fieldStatisticsWindows[window].definition.name,
                  (int)user->_this);
    }
    simCtx->current_io_directory = NULL;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref WriteSwarmField().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see WriteSwarmField()
 */
PetscErrorCode WriteSwarmField(UserCtx *user, const char *field_name, const char *ext)
{
  PetscErrorCode ierr;
  Vec            fieldVec;
  DM     swarm;  

  PetscFunctionBeginUser; /* PETSc macro indicating start of function */

  /* 
   * 1) Retrieve the PetscSwarm from the user context.
   *    Ensure user->swarm is initialized and not NULL.
   */
  swarm = user->swarm;      

  /* 
   * 2) Create a global vector from the specified swarm field.
   *    This function is available in PETSc 3.14.4.
   *    It provides a read/write "view" of the swarm field as a global Vec.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Attempting to create global vector from field: %s\n",
            field_name);
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /*
   * 3) Use your existing WriteFieldData() to write the global vector to a file.
   *    The field name, time index, and extension are passed along for naming.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Calling WriteFieldData for field: %s\n",
            field_name);
  ierr = WriteFieldData(user, field_name, fieldVec, ext);CHKERRQ(ierr);

  /*
   * 4) Destroy the global vector once the data is successfully written.
   *    This step is crucial for avoiding memory leaks. 
   *    DMSwarmDestroyGlobalVectorFromField() is also available in PETSc 3.14.4.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Destroying the global vector for field: %s\n",
            field_name);
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /* Log and return success. */
  LOG_ALLOW(GLOBAL, LOG_INFO,
            "Successfully wrote swarm data for field: %s\n",
            field_name);

  PetscFunctionReturn(0); /* PETSc macro indicating end of function */
}

/**
 * @brief Internal helper implementation: `WriteSwarmIntField()`.
 * @details Local to this translation unit.
 */
PetscErrorCode WriteSwarmIntField(UserCtx *user, const char *field_name, const char *ext)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    Vec            temp_vec;       // Temporary Vec to hold casted data
    PetscInt       nlocal, nglobal,bs,i;
    PetscDataType  field_type;
    void           *field_array_void;
    PetscScalar    *scalar_array;  // Pointer to the temporary Vec's scalar data

    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Casting '%s' to Vec for writing.\n", field_name);

    // Get the swarm field properties
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &nglobal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, &bs, &field_type, &field_array_void); CHKERRQ(ierr);
    PetscCheck(field_type == PETSC_INT || field_type == PETSC_INT64,
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Swarm field '%s' must use PETSC_INT or PETSC_INT64 for integer output, not %s.",
               field_name, PetscDataTypes[field_type]);

    // Create Temporary parallel Vec wit the CORRECT layout
    ierr = VecCreate(PETSC_COMM_WORLD, &temp_vec); CHKERRQ(ierr);
    ierr = VecSetType(temp_vec, VECMPI); CHKERRQ(ierr);
    ierr = VecSetSizes(temp_vec, nlocal*bs, nglobal*bs); CHKERRQ(ierr);
    ierr = VecSetUp(temp_vec); CHKERRQ(ierr);

    // Defining Vector field to mandatory field 'position'
    DMSwarmVectorDefineField(swarm,ParticleFieldName(PARTICLE_FIELD_ID_POSITION));
              
    ierr = VecGetArray(temp_vec, &scalar_array); CHKERRQ(ierr);

    if (field_type == PETSC_INT64) {
        PetscInt64 *int64_array = (PetscInt64 *)field_array_void;
        // Perform the cast from PetscInt64 to PetscScalar
        for (i = 0; i < nlocal*bs; i++) {
            scalar_array[i] = (PetscScalar)int64_array[i];
        }
    }else{
        PetscInt *int_array = (PetscInt *)field_array_void;
        //Perform the cast from PetscInt to PetscScalar
        for (i = 0; i < nlocal*bs; i++) {
            scalar_array[i] = (PetscScalar)int_array[i];
        }
    }

    // Restore access to both arrays
    ierr = VecRestoreArray(temp_vec, &scalar_array); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, field_name, &bs, NULL, &field_array_void); CHKERRQ(ierr);

    // Call your existing writer with the temporary, populated Vec
    ierr = WriteFieldData(user, field_name, temp_vec, ext); CHKERRQ(ierr);

    // Clean up
    ierr = VecDestroy(&temp_vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `WriteAllSwarmFields()`.
 * @details Local to this translation unit.
 */
PetscErrorCode WriteAllSwarmFields(UserCtx *user, const char *checkpoint_directory)
{
    SimCtx         *simCtx = user->simCtx;
    
    PetscFunctionBeginUser;

    // If no swarm is configured or there are no particles, do nothing and return.
    if (!user->swarm || simCtx->np <= 0) {
        PetscFunctionReturn(0);
    }

    PetscCheck(checkpoint_directory != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Checkpoint destination cannot be NULL.");
    PetscCall(PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                            "%s/%s", checkpoint_directory, PICURV_PARTICLE_DIRECTORY));
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    for (PetscInt raw_id = 0; raw_id < PARTICLE_FIELD_ID_COUNT; ++raw_id) {
        const ParticleFieldDescriptor *descriptor = NULL;

        PetscCall(ParticleFieldGetDescriptor((ParticleFieldId)raw_id, &descriptor));
        if (!(descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_CHECKPOINT)) continue;
        if (descriptor->data_type == PETSC_INT || descriptor->data_type == PETSC_INT64) {
            PetscCall(WriteSwarmIntField(user, descriptor->canonical_name, "dat"));
        } else {
            PetscCall(WriteSwarmField(user, descriptor->canonical_name, "dat"));
        }
    }

    simCtx->current_io_directory = NULL;

    PetscFunctionReturn(0);
}

/** @brief Append one payload inventory entry to a checkpoint manifest. */
static PetscErrorCode WriteCheckpointPayloadEntry(FILE *manifest,
                                                  const char *checkpoint_directory,
                                                  PetscInt payload_index,
                                                  const char *relative_path,
                                                  const char *kind,
                                                  const char *field_name,
                                                  PetscInt block,
                                                  const char *layout,
                                                  PetscInt components,
                                                  const char *logical_type,
                                                  PetscInt global_size,
                                                  const char *encoding)
{
    char payload_path[PETSC_MAX_PATH_LEN];
    struct stat payload_stat;

    PetscFunctionBeginUser;
    PetscCheck(manifest != NULL && checkpoint_directory != NULL && relative_path != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Checkpoint payload metadata inputs are required.");
    PetscCall(PetscSNPrintf(payload_path, sizeof(payload_path), "%s/%s",
                            checkpoint_directory, relative_path));
    PetscCheck(stat(payload_path, &payload_stat) == 0 && S_ISREG(payload_stat.st_mode),
               PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
               "Expected checkpoint payload was not written: %s", payload_path);
    PetscCheck(fprintf(manifest,
                       "-checkpoint_payload_%" PetscInt_FMT "_path %s\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_kind %s\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_field %s\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_block %" PetscInt_FMT "\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_layout %s\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_components %" PetscInt_FMT "\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_logical_type %s\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_global_size %" PetscInt_FMT "\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_encoding %s\n"
                       "-checkpoint_payload_%" PetscInt_FMT "_bytes %lld\n",
                       payload_index, relative_path,
                       payload_index, kind,
                       payload_index, field_name,
                       payload_index, block,
                       payload_index, layout,
                       payload_index, components,
                       payload_index, logical_type,
                       payload_index, global_size,
                       payload_index, encoding,
                       payload_index, (long long)payload_stat.st_size) > 0,
               PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
               "Unable to append payload metadata for '%s'.", relative_path);
    PetscFunctionReturn(0);
}

/** @brief Write the complete manifest after every payload has closed successfully. */
static PetscErrorCode WriteCheckpointManifest(SimCtx *simCtx, UserCtx *user,
                                              const char *checkpoint_directory,
                                              const char *reason,
                                              const char *geometry_digest,
                                              PetscInt particle_count,
                                              char manifest_digest[65])
{
    char metadata_path[PETSC_MAX_PATH_LEN];
    PetscInt payload_count = 0;
    PetscInt payload_index = 0;
    FILE *manifest = NULL;

    PetscFunctionBeginUser;
    for (PetscInt block = 0; block < simCtx->block_number; ++block) {
        for (PetscInt raw_id = 0; raw_id < FIELD_ID_COUNT; ++raw_id) {
            const FieldDescriptor *descriptor = NULL;
            PetscCall(FieldGetDescriptor((FieldId)raw_id, &descriptor));
            if (CheckpointFieldIsEnabled(simCtx, descriptor)) ++payload_count;
        }
    }
    if (simCtx->np > 0) {
        for (PetscInt raw_id = 0; raw_id < PARTICLE_FIELD_ID_COUNT; ++raw_id) {
            const ParticleFieldDescriptor *descriptor = NULL;
            PetscCall(ParticleFieldGetDescriptor((ParticleFieldId)raw_id, &descriptor));
            if (descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_CHECKPOINT) ++payload_count;
        }
    }
    if (FieldStatisticsIsActive(simCtx)) {
        for (PetscInt block = 0; block < simCtx->block_number; ++block) {
            if (!user[block].fieldStatisticsStorage) continue;
            for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
                PetscInt storage_payloads = 0;

                PetscCall(PicurvWindowStoragePayloadCount(&user[block].fieldStatisticsStorage[window],
                                                          &storage_payloads));
                payload_count += storage_payloads;
            }
        }
    }

    if (simCtx->rank == 0) {
        PetscCall(PetscSNPrintf(metadata_path, sizeof(metadata_path),
                                "%s/checkpoint.meta", checkpoint_directory));
        manifest = fopen(metadata_path, "w");
        PetscCheck(manifest != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "Unable to create checkpoint metadata '%s'.", metadata_path);
        PetscCheck(fprintf(manifest,
                           "-checkpoint_format %s\n"
                           "-checkpoint_version %d\n"
                           "-checkpoint_picurv_release %s\n"
                           "-checkpoint_picurv_commit %s\n"
                           "-checkpoint_picurv_dirty %s\n"
                           "-checkpoint_step %" PetscInt_FMT "\n"
                           "-checkpoint_time %.17g\n"
                           "-checkpoint_dt %.17g\n"
                           "-checkpoint_reason %s\n"
                           "-checkpoint_geometry_sha256 %s\n"
                           "-checkpoint_block_count %" PetscInt_FMT "\n"
                           "-checkpoint_particles %s\n"
                           "-checkpoint_particle_count %" PetscInt_FMT "\n"
                           "-checkpoint_les %s\n"
                           "-checkpoint_rans %s\n"
                           "-checkpoint_payload_count %" PetscInt_FMT "\n",
                           PICURV_CHECKPOINT_FORMAT, PICURV_CHECKPOINT_VERSION,
                           PICURV_RELEASE_VERSION, PICURV_GIT_COMMIT, PICURV_BUILD_DIRTY,
                           simCtx->step, (double)simCtx->ti, (double)simCtx->dt,
                           reason, geometry_digest, simCtx->block_number,
                           simCtx->np > 0 ? "true" : "false",
                           particle_count,
                           simCtx->les ? "true" : "false",
                           simCtx->rans ? "true" : "false",
                           payload_count) > 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to write checkpoint metadata '%s'.", metadata_path);

        for (PetscInt block = 0; block < simCtx->block_number; ++block) {
            PetscCheck(fprintf(manifest,
                               "-checkpoint_block_%" PetscInt_FMT "_im %" PetscInt_FMT "\n"
                               "-checkpoint_block_%" PetscInt_FMT "_jm %" PetscInt_FMT "\n"
                               "-checkpoint_block_%" PetscInt_FMT "_km %" PetscInt_FMT "\n",
                               block, user[block].IM, block, user[block].JM, block, user[block].KM) > 0,
                       PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                       "Unable to write block metadata for block %" PetscInt_FMT ".", block);
        }
        PetscCheck(fprintf(manifest, "-checkpoint_periodic %d,%d,%d\n",
                           (int)simCtx->i_periodic, (int)simCtx->j_periodic,
                           (int)simCtx->k_periodic) > 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to write checkpoint periodicity metadata.");

        /* Driven-flow controller state. `initial_flux` measures its target once
         * from the starting field, so a restart must read the latched value back
         * rather than re-measure it from a field that has since drifted. */
        PetscCheck(fprintf(manifest,
                           "-checkpoint_driven_flux_latched %s\n"
                           "-checkpoint_driven_flux_target %.17g\n",
                           simCtx->drivenFluxTargetLatched ? "true" : "false",
                           (double)simCtx->targetVolumetricFlux) > 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to write driven-flow controller metadata.");

        /* Window scalars are recorded even when no window is configured, so a
         * restart can tell an absent window list from an unreadable bundle. */
        PetscCheck(fprintf(manifest, "-checkpoint_statistics_window_count %" PetscInt_FMT "\n",
                           FieldStatisticsIsActive(simCtx) ? simCtx->fieldStatisticsWindowCount : 0) > 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to write statistics window count.");
        if (FieldStatisticsIsActive(simCtx)) {
            for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
                const PicurvWindow *state = &simCtx->fieldStatisticsWindows[window];
                char digest[65] = "";
                char groups[PICURV_WINDOW_HASH_GROUP_COUNT][PICURV_WINDOW_HASH_GROUP_LENGTH];
                /* Group digests plus their separators, with headroom so a longer
                 * digest would not silently truncate into a false match. */
                char group_list[PICURV_WINDOW_HASH_GROUP_COUNT * (PICURV_WINDOW_HASH_GROUP_LENGTH + 4)] = "";
                size_t used = 0;

                PetscCall(PicurvWindowComputeHash(&state->definition, digest, groups));
                for (PetscInt group = 0; group < PICURV_WINDOW_HASH_GROUP_COUNT; ++group) {
                    PetscCall(PetscStrlen(group_list, &used));
                    PetscCall(PetscSNPrintf(group_list + used, sizeof(group_list) - used, "%s%s",
                                            group ? "," : "", groups[group]));
                }
                PetscCheck(fprintf(manifest,
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_name %s\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_hash %s\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_hash_groups %s\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_state %s\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_sample_count %" PetscInt_FMT "\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_total_weight %.17g\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_represented_time %.17g\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_last_accepted_time %.17g\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_effective_start %.17g\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_effective_end %.17g\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_activation_step %" PetscInt_FMT "\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_last_event_step %" PetscInt_FMT "\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_next_time_target %" PetscInt_FMT "\n"
                                   "-checkpoint_statistics_window_%" PetscInt_FMT "_restart_count %" PetscInt_FMT "\n",
                                   window, state->definition.name,
                                   window, digest,
                                   window, group_list,
                                   window, PicurvWindowStateName(state->state),
                                   window, state->sample_count,
                                   window, (double)state->total_weight,
                                   window, (double)state->represented_time,
                                   window, (double)state->last_accepted_time,
                                   window, (double)state->effective_start,
                                   window, (double)state->effective_end,
                                   window, state->activation_step,
                                   window, state->last_event_step,
                                   window, state->next_time_target,
                                   window, state->restart_count) > 0,
                           PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                           "Unable to write metadata for statistics window '%s'.",
                           state->definition.name);
            }
        }

        for (PetscInt block = 0; block < simCtx->block_number; ++block) {
            for (PetscInt raw_id = 0; raw_id < FIELD_ID_COUNT; ++raw_id) {
                const FieldDescriptor *descriptor = NULL;
                FieldView view;
                PetscInt global_size = 0;
                char relative_path[PETSC_MAX_PATH_LEN];

                PetscCall(FieldGetDescriptor((FieldId)raw_id, &descriptor));
                if (!CheckpointFieldIsEnabled(simCtx, descriptor)) continue;
                PetscCall(FieldGetView(&user[block], descriptor->id, &view));
                PetscCall(VecGetSize(view.global_vec, &global_size));
                PetscCall(PetscSNPrintf(relative_path, sizeof(relative_path),
                                        "%s/block_%04" PetscInt_FMT "/%s.dat",
                                        PICURV_EULERIAN_DIRECTORY, block, descriptor->canonical_name));
                PetscCall(WriteCheckpointPayloadEntry(manifest, checkpoint_directory,
                                                       payload_index++, relative_path,
                                                       "eulerian", descriptor->canonical_name,
                                                       block, FieldLayoutName(descriptor->layout),
                                                       descriptor->dof, "PetscScalar", global_size,
                                                       "petsc_vec_binary_natural"));
            }
        }
        if (simCtx->np > 0) {
            for (PetscInt raw_id = 0; raw_id < PARTICLE_FIELD_ID_COUNT; ++raw_id) {
                const ParticleFieldDescriptor *descriptor = NULL;
                char relative_path[PETSC_MAX_PATH_LEN];
                const char *encoding = "petsc_vec_binary_global";

                PetscCall(ParticleFieldGetDescriptor((ParticleFieldId)raw_id, &descriptor));
                if (!(descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_CHECKPOINT)) continue;
                if (descriptor->data_type == PETSC_INT || descriptor->data_type == PETSC_INT64) {
                    encoding = "petsc_vec_binary_scalar_cast";
                }
                PetscCall(PetscSNPrintf(relative_path, sizeof(relative_path),
                                        "%s/%s.dat", PICURV_PARTICLE_DIRECTORY,
                                        descriptor->canonical_name));
                PetscCall(WriteCheckpointPayloadEntry(manifest, checkpoint_directory,
                                                       payload_index++, relative_path,
                                                       "particle", descriptor->canonical_name,
                                                       -1, "DMSwarm", descriptor->components,
                                                       PetscDataTypes[descriptor->data_type],
                                                       particle_count * descriptor->components,
                                                       encoding));
            }
        }
        if (FieldStatisticsIsActive(simCtx)) {
            for (PetscInt block = 0; block < simCtx->block_number; ++block) {
                if (!user[block].fieldStatisticsStorage) continue;
                for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
                    const PicurvWindowDefinition *definition =
                        &simCtx->fieldStatisticsWindows[window].definition;
                    const PicurvWindowStorage *storage = &user[block].fieldStatisticsStorage[window];
                    PetscInt storage_payloads = 0;

                    PetscCall(PicurvWindowStoragePayloadCount(storage, &storage_payloads));
                    for (PetscInt index = 0; index < storage_payloads; ++index) {
                        PicurvStatisticsPayload payload;
                        PetscInt global_size = 0;
                        char relative_path[PETSC_MAX_PATH_LEN];
                        char qualified_name[PETSC_MAX_PATH_LEN];

                        PetscCall(PicurvWindowStoragePayload(&user[block], definition, storage,
                                                             index, &payload));
                        PetscCall(VecGetSize(payload.vec, &global_size));
                        PetscCall(FormatStatisticsPath(NULL, window, block, payload.name,
                                                       relative_path, sizeof(relative_path)));
                        /* The inventory field name is qualified by window so two
                         * windows accumulating the same field stay distinguishable
                         * to anything reading the manifest alone. */
                        PetscCall(PetscSNPrintf(qualified_name, sizeof(qualified_name), "%s/%s",
                                                definition->name, payload.name));
                        PetscCall(WriteCheckpointPayloadEntry(manifest, checkpoint_directory,
                                                               payload_index++, relative_path,
                                                               payload.role, qualified_name,
                                                               block, payload.layout,
                                                               payload.components, "PetscScalar",
                                                               global_size,
                                                               "petsc_vec_binary_natural"));
                    }
                }
            }
        }
        PetscCheck(fflush(manifest) == 0 && fsync(fileno(manifest)) == 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to flush checkpoint metadata '%s'.", metadata_path);
        PetscCheck(fclose(manifest) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close checkpoint metadata '%s'.", metadata_path);
        PetscCall(PicurvSHA256File(metadata_path, manifest_digest));
    }
    PetscCallMPI(MPI_Bcast(manifest_digest, 65, MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscFunctionReturn(0);
}

/** @brief Implementation of the transactional full-state checkpoint coordinator. */
PetscErrorCode WriteCheckpointBundle(SimCtx *simCtx, const char *reason)
{
    UserCtx *user = NULL;
    char checkpoints_root[PETSC_MAX_PATH_LEN];
    char final_directory[PETSC_MAX_PATH_LEN];
    char temporary_directory[PETSC_MAX_PATH_LEN];
    char nested_directory[PETSC_MAX_PATH_LEN];
    char commit_path[PETSC_MAX_PATH_LEN];
    char geometry_digest[65] = "";
    char manifest_digest[65] = "";
    PetscBool final_exists = PETSC_FALSE;
    PetscInt particle_count = 0;
    PetscMPIInt writer_pid = 0;

    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL && reason != NULL && reason[0] != '\0',
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Simulation context and checkpoint reason are required.");
    user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    PetscCheck(user != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Finest-level fields must exist before checkpoint output.");

    PetscCall(PetscSNPrintf(checkpoints_root, sizeof(checkpoints_root), "%s/%s",
                            simCtx->output_dir, PICURV_CHECKPOINTS_DIRECTORY));
    PetscCall(FormatCheckpointStepDirectory(simCtx->output_dir, simCtx->step,
                                            final_directory, sizeof(final_directory)));
    PetscCall(PetscTestDirectory(final_directory, 'r', &final_exists));
    if (final_exists) {
        PetscReal saved_time = 0.0;
        PetscCall(ValidateCheckpointBundle(simCtx, user, final_directory, simCtx->step,
                                           &saved_time, NULL, NULL, NULL, NULL));
        PetscCheck(PetscAbsReal(saved_time - simCtx->ti) <=
                   10.0 * PETSC_MACHINE_EPSILON * PetscMax(1.0, PetscAbsReal(simCtx->ti)),
                   PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
                   "Committed checkpoint step %" PetscInt_FMT " has time %.17g, current state has time %.17g.",
                   simCtx->step, (double)saved_time, (double)simCtx->ti);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Checkpoint step %d is already committed; skipping duplicate write.\n",
                  simCtx->step);
        PetscFunctionReturn(0);
    }

    if (simCtx->rank == 0) writer_pid = (PetscMPIInt)getpid();
    PetscCallMPI(MPI_Bcast(&writer_pid, 1, MPI_INT, 0, PETSC_COMM_WORLD));
    PetscCall(PetscSNPrintf(temporary_directory, sizeof(temporary_directory),
                            "%s/.step_%0*" PetscInt_FMT ".incomplete.%d",
                            checkpoints_root, PICURV_CHECKPOINT_STEP_WIDTH,
                            simCtx->step, (int)writer_pid));
    PetscCall(CreateCheckpointDirectoryCollective(simCtx, checkpoints_root));
    PetscCall(CreateCheckpointDirectoryCollective(simCtx, temporary_directory));
    PetscCall(PetscSNPrintf(nested_directory, sizeof(nested_directory), "%s/%s",
                            temporary_directory, PICURV_EULERIAN_DIRECTORY));
    PetscCall(CreateCheckpointDirectoryCollective(simCtx, nested_directory));
    for (PetscInt block = 0; block < simCtx->block_number; ++block) {
        PetscCall(PetscSNPrintf(nested_directory, sizeof(nested_directory),
                                "%s/%s/block_%04" PetscInt_FMT,
                                temporary_directory, PICURV_EULERIAN_DIRECTORY, block));
        PetscCall(CreateCheckpointDirectoryCollective(simCtx, nested_directory));
    }
    if (simCtx->np > 0) {
        PetscCall(PetscSNPrintf(nested_directory, sizeof(nested_directory), "%s/%s",
                                temporary_directory, PICURV_PARTICLE_DIRECTORY));
        PetscCall(CreateCheckpointDirectoryCollective(simCtx, nested_directory));
        PetscCall(DMSwarmGetSize(user->swarm, &particle_count));
    }
    if (FieldStatisticsIsActive(simCtx)) {
        PetscCall(FormatStatisticsPath(temporary_directory, -1, -1, NULL,
                                       nested_directory, sizeof(nested_directory)));
        PetscCall(CreateCheckpointDirectoryCollective(simCtx, nested_directory));
        for (PetscInt window = 0; window < simCtx->fieldStatisticsWindowCount; ++window) {
            PetscCall(FormatStatisticsPath(temporary_directory, window, -1, NULL,
                                           nested_directory, sizeof(nested_directory)));
            PetscCall(CreateCheckpointDirectoryCollective(simCtx, nested_directory));
            for (PetscInt block = 0; block < simCtx->block_number; ++block) {
                PetscCall(FormatStatisticsPath(temporary_directory, window, block, NULL,
                                               nested_directory, sizeof(nested_directory)));
                PetscCall(CreateCheckpointDirectoryCollective(simCtx, nested_directory));
            }
        }
    }

    PetscCall(ComputeCheckpointGeometrySHA256(simCtx, user, geometry_digest));
    for (PetscInt block = 0; block < simCtx->block_number; ++block) {
        PetscCall(WriteSimulationFields(&user[block], temporary_directory));
    }
    if (simCtx->np > 0) PetscCall(WriteAllSwarmFields(user, temporary_directory));
    for (PetscInt block = 0; block < simCtx->block_number; ++block) {
        PetscCall(WriteStatisticsFields(&user[block], temporary_directory));
    }
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));

    PetscCall(WriteCheckpointManifest(simCtx, user, temporary_directory, reason,
                                      geometry_digest, particle_count, manifest_digest));
    if (simCtx->rank == 0) {
        FILE *commit_file = NULL;

        PetscCall(PetscSNPrintf(commit_path, sizeof(commit_path), "%s/COMMITTED", temporary_directory));
        commit_file = fopen(commit_path, "w");
        PetscCheck(commit_file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "Unable to create checkpoint commit marker '%s'.", commit_path);
        PetscCheck(fprintf(commit_file, "%s\n", manifest_digest) > 0 &&
                   fflush(commit_file) == 0 && fsync(fileno(commit_file)) == 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to write checkpoint commit marker '%s'.", commit_path);
        PetscCheck(fclose(commit_file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close checkpoint commit marker '%s'.", commit_path);
        PetscCheck(rename(temporary_directory, final_directory) == 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to commit checkpoint '%s': %s", final_directory, strerror(errno));
    }
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    LOG_ALLOW(GLOBAL, LOG_INFO,
              "Committed checkpoint step %d at t=%.17g (%s): %s\n",
              simCtx->step, (double)simCtx->ti, reason, final_directory);
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `VecToArrayOnRank0()`.
 * @details Local to this translation unit.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut)
{
    MPI_Comm comm;
    PetscMPIInt rank;
    Vec sequential_vec = NULL;

    PetscFunctionBeginUser;
    PetscCheck(inVec != NULL && N != NULL && arrayOut != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Vector, size output, and array output are required.");
    PetscCall(PetscObjectGetComm((PetscObject)inVec, &comm));
    PetscCallMPI(MPI_Comm_rank(comm, &rank));
    PetscCall(VecGetSize(inVec, N));
    *arrayOut = NULL;

    PetscCall(GatherVectorToRankZero(inVec, &sequential_vec));
    if (rank == 0) {
        const PetscScalar *values = NULL;
        PetscInt local_size = 0;

        PetscCall(VecGetLocalSize(sequential_vec, &local_size));
        PetscCall(PetscMalloc1(local_size, arrayOut));
        PetscCall(VecGetArrayRead(sequential_vec, &values));
        for (PetscInt index = 0; index < local_size; ++index) {
            (*arrayOut)[index] = (double)PetscRealPart(values[index]);
        }
        PetscCall(VecRestoreArrayRead(sequential_vec, &values));
    }
    PetscCall(VecDestroy(&sequential_vec));
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `SwarmFieldToArrayOnRank0()`.
 * @details Local to this translation unit.
 */
PetscErrorCode SwarmFieldToArrayOnRank0(DM swarm, const char *field_name,
                                        PetscInt *n_total_particles, PetscInt *n_components,
                                        PetscDataType *field_type_out, void **gathered_array)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    PetscInt       nlocal, nglobal, bs;
    PetscDataType  field_type;
    void           *local_array_void;
    size_t         element_size = 0;

    PetscFunctionBeginUser;

    PetscCheck(swarm != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DMSwarm cannot be NULL.");
    PetscCheck(field_name != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Swarm field name cannot be NULL.");
    PetscCheck(n_total_particles != NULL && n_components != NULL && field_type_out != NULL && gathered_array != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Swarm gather output pointers cannot be NULL.");

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // All ranks get swarm properties to determine send/receive counts
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &nglobal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, &bs, &field_type, &local_array_void); CHKERRQ(ierr);
    
    // Determine the size of one element of the field's data type
    if (field_type == PETSC_INT64) element_size = sizeof(PetscInt64);
    else if (field_type == PETSC_INT) element_size = sizeof(PetscInt);
    else if (field_type == PETSC_REAL) element_size = sizeof(PetscReal);
#if defined(PETSC_USE_COMPLEX)
    else if (field_type == PETSC_SCALAR) element_size = sizeof(PetscScalar);
#endif
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
                 "Swarm field '%s' uses unsupported gathered data type %s.",
                 field_name, PetscDataTypes[field_type]);

    *field_type_out = field_type;
    *n_total_particles = nglobal;
    *n_components = bs;
    *gathered_array = NULL;

    if (size == 1) { // Serial case is a simple copy
        if (rank == 0) {
            ierr = PetscMalloc(nglobal * bs * element_size, gathered_array); CHKERRQ(ierr);
            ierr = PetscMemcpy(*gathered_array, local_array_void, nglobal * bs * element_size); CHKERRQ(ierr);
        }
    } else { // Parallel case: use MPI_Gatherv
        PetscInt *recvcounts = NULL, *displs = NULL;
        if (rank == 0) {
            ierr = PetscMalloc1(size, &recvcounts); CHKERRQ(ierr);
            ierr = PetscMalloc1(size, &displs); CHKERRQ(ierr);
        }
        PetscInt sendcount = nlocal * bs;
        
        // Gather the number of elements (not bytes) from each rank
        ierr = MPI_Gather(&sendcount, 1, MPIU_INT, recvcounts, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        if (rank == 0) {
            displs[0] = 0;
            // Convert counts and calculate displacements in terms of BYTES
            for (PetscMPIInt i = 0; i < size; i++) recvcounts[i] *= element_size;
            for (PetscMPIInt i = 1; i < size; i++) displs[i] = displs[i-1] + recvcounts[i-1];
            
            ierr = PetscMalloc(nglobal * bs * element_size, gathered_array); CHKERRQ(ierr);
        }
        
        // Use Gatherv with MPI_BYTE to handle any data type generically
        ierr = MPI_Gatherv(local_array_void, nlocal * bs * element_size, MPI_BYTE,
                           *gathered_array, recvcounts, displs, MPI_BYTE,
                           0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        if (rank == 0) {
            ierr = PetscFree(recvcounts); CHKERRQ(ierr);
            ierr = PetscFree(displs); CHKERRQ(ierr);
        }
    }

    ierr = DMSwarmRestoreField(swarm, field_name, &bs, NULL, &local_array_void); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Emit the rank-zero startup summary from the effective simulation context.
 * @details Reports only configuration that applies to the selected run mode. In
 *          particular, pseudo-CFL is a Dual Time Picard--Jameson RK control and
 *          is deliberately omitted for explicit and Newton--Krylov momentum solves.
 */
PetscErrorCode DisplayBanner(SimCtx *simCtx) // bboxlist is only valid on rank 0
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Cmpnts         global_min_coords, global_max_coords;
    PetscReal      StartTime;
    PetscInt       StartStep,StepsToRun,total_num_particles;
    PetscMPIInt    num_mpi_procs;
    const char    *log_level_name;
    const char    *convergence_mode_name;

    // SimCtx *simCtx = user->simCtx;
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    num_mpi_procs = simCtx->size;
    StartTime = simCtx->StartTime;
    StartStep = simCtx->StartStep;
    StepsToRun = simCtx->StepsToRun;
    total_num_particles = simCtx->np;
    BoundingBox *bboxlist_on_rank0 = simCtx->bboxlist;
    

    PetscFunctionBeginUser;

    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DisplayBanner - UserCtx pointer is NULL.");
    switch (get_log_level()) {
        case LOG_ERROR:   log_level_name = "ERROR"; break;
        case LOG_WARNING: log_level_name = "WARNING"; break;
        case LOG_INFO:    log_level_name = "INFO"; break;
        case LOG_DEBUG:   log_level_name = "DEBUG"; break;
        case LOG_TRACE:   log_level_name = "TRACE"; break;
        case LOG_VERBOSE: log_level_name = "VERBOSE"; break;
        default:          log_level_name = "UNKNOWN"; break;
    }
    switch (simCtx->solutionConvergenceMode) {
        case SOLUTION_CONVERGENCE_STEADY_DETERMINISTIC: convergence_mode_name = "STEADY_DETERMINISTIC"; break;
        case SOLUTION_CONVERGENCE_PERIODIC_DETERMINISTIC: convergence_mode_name = "PERIODIC_DETERMINISTIC"; break;
        case SOLUTION_CONVERGENCE_STATISTICAL_STEADY: convergence_mode_name = "STATISTICAL_STEADY"; break;
        case SOLUTION_CONVERGENCE_TRANSIENT: convergence_mode_name = "TRANSIENT"; break;
        default: convergence_mode_name = "UNKNOWN"; break;
    }
    global_min_coords = user->bbox.min_coords;
    global_max_coords = user->bbox.max_coords;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        // If global_domain_bbox is not pre-populated in UserCtx, compute it here from bboxlist_on_rank0
        // This assumes bboxlist_on_rank0 is valid and contains all local bounding boxes on rank 0.
        if (bboxlist_on_rank0 && num_mpi_procs > 0) {
            global_min_coords = bboxlist_on_rank0[0].min_coords;
            global_max_coords = bboxlist_on_rank0[0].max_coords;
            for (PetscMPIInt p = 1; p < num_mpi_procs; ++p) {
                global_min_coords.x = PetscMin(global_min_coords.x, bboxlist_on_rank0[p].min_coords.x);
                global_min_coords.y = PetscMin(global_min_coords.y, bboxlist_on_rank0[p].min_coords.y);
                global_min_coords.z = PetscMin(global_min_coords.z, bboxlist_on_rank0[p].min_coords.z);
                global_max_coords.x = PetscMax(global_max_coords.x, bboxlist_on_rank0[p].max_coords.x);
                global_max_coords.y = PetscMax(global_max_coords.y, bboxlist_on_rank0[p].max_coords.y);
                global_max_coords.z = PetscMax(global_max_coords.z, bboxlist_on_rank0[p].max_coords.z);
            }
            // Optionally store this in user->global_domain_bbox if it's useful elsewhere
            // user->global_domain_bbox.min_coords = global_min_coords;
            // user->global_domain_bbox.max_coords = global_max_coords;
        } else {
            // Fallback or warning if bboxlist is not available for global calculation
            LOG_ALLOW(LOCAL, LOG_WARNING, "(Rank 0) - bboxlist not provided or num_mpi_procs <=0; using user->bbox for domain bounds.\n");
	    // global_min_coords = user->bbox.min_coords; // Use local bbox of rank 0 as fallback
	    // global_max_coords = user->bbox.max_coords;
        }

        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "                          CASE SUMMARY                       \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Grid Points     : %d X %d X %d\n", user->IM, user->JM, user->KM); CHKERRQ(ierr);
	    ierr = PetscPrintf(PETSC_COMM_SELF, " Cells           : %d X %d X %d\n", user->IM - 1, user->JM - 1, user->KM - 1); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (X)    : %.6f to %.6f\n", (double)global_min_coords.x, (double)global_max_coords.x); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (Y)    : %.6f to %.6f\n", (double)global_min_coords.y, (double)global_max_coords.y); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (Z)    : %.6f to %.6f\n", (double)global_min_coords.z, (double)global_max_coords.z); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Periodic Axes (BC-derived)  : I=%s, J=%s, K=%s\n",
                           simCtx->i_periodic ? "YES" : "NO",
                           simCtx->j_periodic ? "YES" : "NO",
                           simCtx->k_periodic ? "YES" : "NO"); CHKERRQ(ierr);
        for (PetscInt axis = 0; axis < 3; axis++) {
            if (!user->periodic_translation_valid[axis]) continue;
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               " Periodic %c Translation     : (%.6e, %.6e, %.6e)\n",
                               "IJK"[axis],
                               (double)user->periodic_translation[axis].x,
                               (double)user->periodic_translation[axis].y,
                               (double)user->periodic_translation[axis].z); CHKERRQ(ierr);
        }
        if (total_num_particles > 0 &&
            (simCtx->i_periodic || simCtx->j_periodic || simCtx->k_periodic)) {
            ierr = PetscPrintf(PETSC_COMM_SELF,
                               " Particle Periodic Wrapping  : UNSUPPORTED (Eulerian periodicity only)\n"); CHKERRQ(ierr);
        }
        if(strcmp(simCtx->eulerianSource,"load")==0 || strcmp(simCtx->eulerianSource,"solve")==0){
            ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------- Boundary Conditions --------------------\n"); CHKERRQ(ierr);
            const int face_name_width = 17; // Adjusted for longer names (Zeta,Eta,Xi)
            for (PetscInt i_face = 0; i_face < 6; ++i_face) {
                BCFace current_face = (BCFace)i_face;
                // The BCFaceToString will now return the Xi, Eta, Zeta versions
                const char* face_str = BCFaceToString(current_face); 
                const char* bc_type_str = BCTypeToString(user->boundary_faces[current_face].mathematical_type);
                const char* bc_handler_type_str = BCHandlerTypeToString(user->boundary_faces[current_face].handler_type);
                if(user->boundary_faces[current_face].mathematical_type == INLET){
                    if(user->boundary_faces[current_face].handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY){
                        Cmpnts inlet_velocity = {0.0,0.0,0.0};
                        PetscBool found;
                        ierr = GetBCParamReal(user->boundary_faces[current_face].params,"vx",&inlet_velocity.x,&found); CHKERRQ(ierr);
                        ierr = GetBCParamReal(user->boundary_faces[current_face].params,"vy",&inlet_velocity.y,&found); CHKERRQ(ierr);
                        ierr = GetBCParamReal(user->boundary_faces[current_face].params,"vz",&inlet_velocity.z,&found); CHKERRQ(ierr);
                        ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s - %s - [%.4f,%.4f,%.4f]\n", 
                        face_name_width, face_str, bc_type_str, bc_handler_type_str,inlet_velocity.x,inlet_velocity.y,inlet_velocity.z); CHKERRQ(ierr);
                        } else if(user->boundary_faces[current_face].handler_type == BC_HANDLER_INLET_PARABOLIC){
                        PetscReal v_max = 0.0;
                        PetscBool found;
                        ierr = GetBCParamReal(user->boundary_faces[current_face].params,"v_max",&v_max,&found); CHKERRQ(ierr);
                        ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s - %s - v_max=%.4f\n",
                        face_name_width, face_str, bc_type_str, bc_handler_type_str, v_max); CHKERRQ(ierr);
                        } else if(user->boundary_faces[current_face].handler_type == BC_HANDLER_INLET_PROFILE_FROM_FILE){
                        const char *source_file = "(missing)";
                        for (BC_Param *param = user->boundary_faces[current_face].params; param; param = param->next) {
                            if (strcasecmp(param->key, "source_file") == 0 && param->value) {
                                source_file = param->value;
                                break;
                            }
                        }
                        ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s - %s - source_file=%s\n",
                        face_name_width, face_str, bc_type_str, bc_handler_type_str, source_file); CHKERRQ(ierr);
                        }
                    } else if(user->boundary_faces[current_face].handler_type == BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX){
                        PetscBool trimflag,foundtrimflag;
                        ierr = GetDrivenSeamFluxFlag(user->boundary_faces[current_face].params,&trimflag,&foundtrimflag); CHKERRQ(ierr);
                        ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s - %s - [from initial state] - %s\n",
                        face_name_width, face_str, bc_type_str, bc_handler_type_str,trimflag?"Enforce seam flux":"Seam flux not enforced"); CHKERRQ(ierr);
                    } else if(user->boundary_faces[current_face].handler_type == BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX){
                        PetscReal flux;
                        PetscBool trimflag,foundflux,foundtrimflag;
                        ierr = GetBCParamReal(user->boundary_faces[current_face].params,"target_flux",&flux,&foundflux); CHKERRQ(ierr);
                        ierr = GetDrivenSeamFluxFlag(user->boundary_faces[current_face].params,&trimflag,&foundtrimflag); CHKERRQ(ierr);
                        ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s - %s - [%.4f] - %s\n", 
                        face_name_width, face_str, bc_type_str, bc_handler_type_str,flux,trimflag?"Enforce seam flux":"Seam flux not enforced"); CHKERRQ(ierr);
                    } else{    
                        ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s - %s\n", 
                                face_name_width, face_str, bc_type_str,bc_handler_type_str); CHKERRQ(ierr);
                    }
            }
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------------------------------------------------\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Run Mode                   : %s\n", simCtx->OnlySetup ? "SETUP ONLY" : "Full Simulation"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Start Time                  : %.4f\n", (double)StartTime); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Timestep Size               : %.4f\n", (double)simCtx->dt); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Starting Step               : %d\n", StartStep); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Total Steps to Run          : %d\n", StepsToRun); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Ending Step                : %d\n", StartStep + StepsToRun); CHKERRQ(ierr);
        if (simCtx->tiout > 0) {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Field/Restart Cadence      : every %d step(s)\n", simCtx->tiout); CHKERRQ(ierr);
        } else {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Field/Restart Cadence      : DISABLED\n"); CHKERRQ(ierr);
        }
        /* Recorded whether or not statistics are configured, so a log says plainly
         * whether monitoring was active rather than leaving its absence ambiguous. */
        if (FieldStatisticsIsActive(simCtx) && simCtx->statisticsConsoleOutputFreq > 0) {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Statistics Console Cadence : every %d step(s), %d window(s)\n",
                               simCtx->statisticsConsoleOutputFreq, simCtx->fieldStatisticsWindowCount); CHKERRQ(ierr);
        } else if (FieldStatisticsIsActive(simCtx)) {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Statistics Console Cadence : DISABLED (%d window(s) accumulating)\n",
                               simCtx->fieldStatisticsWindowCount); CHKERRQ(ierr);
        } else {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Statistics Console Cadence : DISABLED (no window configured)\n"); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, " Immersed Boundary          : %s\n", simCtx->immersed ? "ENABLED" : "DISABLED"); CHKERRQ(ierr);
        if (simCtx->walltimeGuardEnabled) {
            ierr = PetscPrintf(
                PETSC_COMM_SELF,
                " Runtime Walltime Guard     : %s (warmup=%d, multiplier=%.2f, min=%.1f s, alpha=%.2f)\n",
                simCtx->walltimeGuardActive ? "ENABLED" : "CONFIGURED BUT INACTIVE",
                simCtx->walltimeGuardWarmupSteps,
                (double)simCtx->walltimeGuardMultiplier,
                (double)simCtx->walltimeGuardMinSeconds,
                (double)simCtx->walltimeGuardEstimatorAlpha
            ); CHKERRQ(ierr);
        } else {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Runtime Walltime Guard     : DISABLED\n"); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, " Console Log Level           : %s\n", log_level_name); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Profiling Timestep Output   : %s\n", simCtx->profilingTimestepMode); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Profiling Final Summary     : %s\n", simCtx->profilingFinalSummary ? "ENABLED" : "DISABLED"); CHKERRQ(ierr);
        if (simCtx->runtimeMemoryLogEnabled) {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Runtime Memory Log          : ENABLED (%s)\n", simCtx->runtimeMemoryLogFile); CHKERRQ(ierr);
        } else {
            ierr = PetscPrintf(PETSC_COMM_SELF, " Runtime Memory Log          : DISABLED\n"); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, " Solution Convergence Log    : %s\n",
                           simCtx->solutionConvergenceEnabled ? "ENABLED" : "DISABLED"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Number of MPI Processes     : %d\n", num_mpi_procs); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Number of Particles         : %d\n", total_num_particles); CHKERRQ(ierr);
        if (simCtx->np > 0) {
            const char *particle_init_str = ParticleInitializationToString(simCtx->ParticleInitialization);

            if (simCtx->particleConsoleOutputFreq > 0) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Console Cadence   : every %d step(s)\n", simCtx->particleConsoleOutputFreq); CHKERRQ(ierr);
            } else {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Console Cadence   : DISABLED\n"); CHKERRQ(ierr);
            }
            ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Log Row Sampling  : every %d particle(s)\n", simCtx->LoggingFrequency); CHKERRQ(ierr);
            if (simCtx->StartStep > 0) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Restart Mode      : %s\n", simCtx->particleRestartMode); CHKERRQ(ierr);
            }
            ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Initialization Mode: %s\n", particle_init_str); CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, " Interpolation Method       : %s\n",
                simCtx->interpolationMethod == INTERP_TRILINEAR ? "Trilinear (direct cell-center)" : "CornerAveraged (legacy)"); CHKERRQ(ierr);
            if (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM ||
                simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_EDGES) {
                if (user->inletFaceDefined) {
                    ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : %s (Enum Val: %d)\n", BCFaceToString(user->identifiedInletBCFace), user->identifiedInletBCFace); CHKERRQ(ierr);
                } else {
                    ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : --- (No INLET face identified)\n"); CHKERRQ(ierr);
                }
            }
        }
        if(strcmp(simCtx->eulerianSource,"solve")==0 || strcmp(simCtx->eulerianSource,"load")==0){
            ierr = PetscPrintf(PETSC_COMM_WORLD," Reynolds Number             : %le\n", simCtx->ren); CHKERRQ(ierr);
            //ierr = PetscPrintf(PETSC_COMM_WORLD," Von-Neumann Number          : %le\n", simCtx->vnn); CHKERRQ(ierr);
            if(strcmp(simCtx->eulerianSource,"solve")==0){
                //ierr = PetscPrintf(PETSC_COMM_WORLD," Stanton Number              : %le\n", simCtx->st); CHKERRQ(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD," Momentum Equation Solver    : %s\n", MomentumSolverTypeToString(simCtx->mom_solver_type)); CHKERRQ(ierr);
                if (simCtx->mom_solver_type == MOMENTUM_SOLVER_DUALTIME_PICARD_JAMESON_RK) {
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Initial Pseudo-CFL (Courant): %le\n", simCtx->pseudo_cfl); CHKERRQ(ierr);
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Pseudo-CFL Range            : [%le, %le]\n", simCtx->min_pseudo_cfl, simCtx->max_pseudo_cfl); CHKERRQ(ierr);
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Pseudo-CFL Adaptation       : growth=%le, reduction=%le, backtrack=%s\n",
                                       simCtx->pseudo_cfl_growth_factor, simCtx->pseudo_cfl_reduction_factor,
                                       simCtx->no_pseudo_cfl_backtrack ? "DISABLED" : "ENABLED"); CHKERRQ(ierr);
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Pseudo-Time Iteration Limit : %d\n", simCtx->mom_max_pseudo_steps); CHKERRQ(ierr);
                } else if (simCtx->mom_solver_type == MOMENTUM_SOLVER_NEWTON_KRYLOV) {
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Newton-Krylov PETSc Controls: SNES/KSP options (mom_nk_*)\n"); CHKERRQ(ierr);
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Newton-Krylov History Log   : %s\n", simCtx->mom_nk_monitor_history ? "ENABLED" : "DISABLED"); CHKERRQ(ierr);
                } else if (simCtx->mom_solver_type == MOMENTUM_SOLVER_EXPLICIT_RK) {
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Pseudo-Time Controller      : NOT APPLICABLE\n"); CHKERRQ(ierr);
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD," Solution Convergence Mode   : %s\n", convergence_mode_name); CHKERRQ(ierr);
                if (simCtx->solutionConvergenceMode == SOLUTION_CONVERGENCE_PERIODIC_DETERMINISTIC) {
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Convergence Period          : %d step(s)\n", simCtx->solutionConvergencePeriodSteps); CHKERRQ(ierr);
                } else if (simCtx->solutionConvergenceMode == SOLUTION_CONVERGENCE_STATISTICAL_STEADY) {
                    ierr = PetscPrintf(PETSC_COMM_WORLD," Convergence Window          : %d step(s)\n", simCtx->solutionConvergenceWindowSteps); CHKERRQ(ierr);
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD," Large Eddy Simulation Model : %s\n", LESModelToString(simCtx->les)); CHKERRQ(ierr);
                if (simCtx->les != NO_LES_MODEL) {
                    const LESConfig *les = &simCtx->les_config;

                    /* Reported with the user-facing spellings, so a line here can be
                       matched against the case file that produced it. */
                    ierr = PetscPrintf(PETSC_COMM_WORLD," LES Filter Width            : %s\n",
                                       LESFilterWidthModelToString(les->filter_width_model)); CHKERRQ(ierr);
                    if (simCtx->les == CONSTANT_SMAGORINSKY) {
                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Smagorinsky Constant    : %.4f (no coefficient field)\n",
                                           (double)les->constant_cs); CHKERRQ(ierr);
                    } else {
                        char directions[8] = "";
                        PetscInt used = 0;

                        if (les->averaging_direction[0]) directions[used++] = 'i';
                        if (les->averaging_direction[1]) directions[used++] = 'j';
                        if (les->averaging_direction[2]) directions[used++] = 'k';
                        directions[used] = '\0';

                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Test Filter             : %s (width ratio %.3f)\n",
                                           LESTestFilterKernelToString(les->test_filter_kernel),
                                           (double)les->test_filter_width_ratio); CHKERRQ(ierr);
                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Dynamic Update Cadence  : every %d step(s)\n",
                                           (int)les->dynamic_frequency); CHKERRQ(ierr);
                        /* An empty direction list under homogeneous averaging is not a
                           missing answer: the directions are derived from periodicity at
                           the first update, once the block is known. */
                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Coefficient Averaging   : %s%s%s\n",
                                           LESAveragingModeToString(les->averaging_mode),
                                           used > 0 ? " over " : "",
                                           used > 0 ? directions : ""); CHKERRQ(ierr);
                        if (les->clip_mode == LES_CLIP_CLAMP) {
                            ierr = PetscPrintf(PETSC_COMM_WORLD," LES Coefficient Limiting    : clamp (max Cs %.3f)\n",
                                               (double)les->max_cs); CHKERRQ(ierr);
                        } else {
                            ierr = PetscPrintf(PETSC_COMM_WORLD," LES Coefficient Limiting    : %s\n",
                                               LESClipModeToString(les->clip_mode)); CHKERRQ(ierr);
                        }
                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Total Viscosity Floor   : %.3f x molecular\n",
                                           (double)les->min_viscosity_ratio); CHKERRQ(ierr);
                    }
                    ierr = PetscPrintf(PETSC_COMM_WORLD," LES Gradient (Clark) Term   : %s\n",
                                       simCtx->les_gradient_model ? "ENABLED" : "DISABLED"); CHKERRQ(ierr);
                    if (les->diagnostics_enabled) {
                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Coefficient Diagnostics : ENABLED (les_coefficient.csv, every %d step(s))\n",
                                           (int)les->diagnostics_cadence); CHKERRQ(ierr);
                    } else {
                        ierr = PetscPrintf(PETSC_COMM_WORLD," LES Coefficient Diagnostics : DISABLED\n"); CHKERRQ(ierr);
                    }
                }
            }
            if (strcmp(simCtx->eulerianSource, "load") == 0) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Eulerian State Source       : load (%s)\n",
                                   simCtx->restart_dir); CHKERRQ(ierr);
            } else if (simCtx->StartStep > 0) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Eulerian State Source       : restart step %d (%s)\n",
                                   simCtx->StartStep, simCtx->restart_dir); CHKERRQ(ierr);
            } else {
                const char* field_init_str = InitialConditionModeToString(simCtx->initialConditionMode);
                ierr = PetscPrintf(PETSC_COMM_SELF, " Eulerian State Source       : initial condition (%s)\n",
                                   field_init_str); CHKERRQ(ierr);
            }
            if (strcmp(simCtx->eulerianSource, "solve") == 0 && simCtx->StartStep == 0 &&
                simCtx->initialConditionMode == IC_MODE_CONSTANT_CARTESIAN) {
                    ierr = PetscPrintf(PETSC_COMM_SELF,
                        " Constant Velocity (Cart.)   : x=%.4f  y=%.4f  z=%.4f\n",
                        (double)simCtx->InitialConstantContra.x,
                        (double)simCtx->InitialConstantContra.y,
                        (double)simCtx->InitialConstantContra.z); CHKERRQ(ierr);
            } else if (strcmp(simCtx->eulerianSource, "solve") == 0 && simCtx->StartStep == 0 &&
                       simCtx->initialConditionMode == IC_MODE_CONSTANT_STREAMWISE) {
                    ierr = PetscPrintf(PETSC_COMM_SELF,
                        " Constant Velocity (Curv.)   : speed=%.4f  direction=%s\n",
                        (double)simCtx->icVelocityPhysical,
                        FlowDirectionToString(simCtx->flowDirection)); CHKERRQ(ierr);
            } else if (strcmp(simCtx->eulerianSource, "solve") == 0 && simCtx->StartStep == 0 &&
                       simCtx->initialConditionMode == IC_MODE_POISEUILLE) {
                ierr = PetscPrintf(PETSC_COMM_SELF,
                    " Poiseuille Peak Velocity    : speed=%.4f  direction=%s\n",
                    (double)simCtx->icVelocityPhysical,
                    FlowDirectionToString(simCtx->flowDirection)); CHKERRQ(ierr);
            } else if (strcmp(simCtx->eulerianSource, "solve") == 0 && simCtx->StartStep == 0 &&
                       simCtx->initialConditionMode == IC_MODE_FILE) {
                ierr = PetscPrintf(PETSC_COMM_SELF,
                    " Initial Velocity File       : field=%s directory=%s\n",
                    simCtx->initialConditionField == IC_FIELD_UCAT ? "Ucat" : "Ucont",
                    simCtx->initialConditionDirectory); CHKERRQ(ierr);
            }
        } else if(strcmp(simCtx->eulerianSource,"analytical")==0){
            ierr = PetscPrintf(PETSC_COMM_WORLD," Analytical Solution Type : %s\n", simCtx->AnalyticalSolutionType); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ParsePostProcessingSettings"
/**
 * @brief Internal helper implementation: `ParsePostProcessingSettings()`.
 * @details Local to this translation unit.
 */
PetscErrorCode  ParsePostProcessingSettings(SimCtx *simCtx)
{
    FILE *file;
    char line[1024];
    PetscBool startTimeSet, endTimeSet, timeStepSet;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!simCtx || !simCtx->pps) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "SimCtx or its pps member is NULL in ParsePostProcessingSettings.");
    }

    char *configFile = simCtx->PostprocessingControlFile;
    PostProcessParams *pps = simCtx->pps;


    // --- 1. Set Sane Defaults First ---
    pps->startTime = 0;
    pps->endTime = 0;
    pps->timeStep = 1;
    pps->outputParticles = PETSC_FALSE;
    pps->particle_output_freq = simCtx->LoggingFrequency; // Default to logging frequency;
    strcpy(pps->process_pipeline, "");
    strcpy(pps->output_fields_instantaneous, "Ucat,P");
    strcpy(pps->output_prefix, "Field");
    strcpy(pps->particle_output_prefix,"Particle");
    strcpy(pps->particle_fields,"velocity,CellID,weight,pid");
    strcpy(pps->particle_pipeline,"");
    strncpy(pps->statistics_pipeline,      "",       MAX_PIPELINE_LENGTH - 1);
    strncpy(pps->statistics_output_prefix, "Stats",  MAX_FILENAME_LENGTH - 1);
    strcpy(pps->particleExt,"dat"); // The input file format for particles.
    strcpy(pps->eulerianExt,"dat"); // The input file format for Eulerian fields.
    /* A negative source step means "the last step this recipe covers", so a recipe
     * that does not name one still derives from a bundle that exists. */
    pps->field_statistics_windows[0] = '\0';
    strcpy(pps->field_statistics_outputs, "mean,reynolds_stress,rms,tke,flux");
    strcpy(pps->field_statistics_formats, "vtk");
    strcpy(pps->field_statistics_output_prefix, "Field");
    pps->field_statistics_source_step = -1;
    pps->reference[0] = pps->reference[1] = pps->reference[2] = 1;
    strncpy(pps->source_dir, simCtx->output_dir, sizeof(pps->source_dir) - 1);
    pps->source_dir[sizeof(pps->source_dir) - 1] = '\0'; // Ensure null-termination

    // --- 2. Parse the Configuration File (overrides defaults) ---
    file = fopen(configFile, "r");
    if (file) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Parsing post-processing config file: %s\n", configFile);
        while (fgets(line, sizeof(line), file)) {
            char *key, *value, *comment;
            comment = strchr(line, '#'); if (comment) *comment = '\0';
            TrimWhitespace(line); if (strlen(line) == 0) continue;
            key = strtok(line, "="); value = strtok(NULL, "=");
            if (key && value) {
                TrimWhitespace(key); TrimWhitespace(value);
                if (strcmp(key, "startTime") == 0) pps->startTime = atoi(value);
                else if (strcmp(key, "endTime") == 0) pps->endTime = atoi(value);
                else if (strcmp(key, "timeStep") == 0) pps->timeStep = atoi(value);
                else if (strcmp(key, "output_particles") == 0) {
                    if (strcasecmp(value, "true") == 0) pps->outputParticles = PETSC_TRUE;
                }
                else if (strcasecmp(key, "process_pipeline") == 0) {
                    strncpy(pps->process_pipeline, value, MAX_PIPELINE_LENGTH - 1);
                    pps->process_pipeline[MAX_PIPELINE_LENGTH - 1] = '\0'; // Ensure null-termination
                } else if (strcasecmp(key, "field_statistics_windows") == 0) {
                    strncpy(pps->field_statistics_windows, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->field_statistics_windows[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "field_statistics_formats") == 0) {
                    strncpy(pps->field_statistics_formats, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->field_statistics_formats[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "field_statistics_outputs") == 0) {
                    strncpy(pps->field_statistics_outputs, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->field_statistics_outputs[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "field_statistics_source_step") == 0) {
                    pps->field_statistics_source_step = atoi(value);
                } else if (strcasecmp(key, "field_statistics_output_prefix") == 0) {
                    strncpy(pps->field_statistics_output_prefix, value,
                            sizeof(pps->field_statistics_output_prefix) - 1);
                    pps->field_statistics_output_prefix[
                        sizeof(pps->field_statistics_output_prefix) - 1] = '\0';
                } else if (strcasecmp(key, "output_fields_instantaneous") == 0) {
                    strncpy(pps->output_fields_instantaneous, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->output_fields_instantaneous[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "output_prefix") == 0) {
                    strncpy(pps->output_prefix, value, MAX_FILENAME_LENGTH - 1);
                    pps->output_prefix[MAX_FILENAME_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_output_prefix") == 0) {
                    strncpy(pps->particle_output_prefix, value, MAX_FILENAME_LENGTH - 1);
                    pps->particle_output_prefix[MAX_FILENAME_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_fields_instantaneous") == 0) {
                    strncpy(pps->particle_fields, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->particle_fields[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_pipeline") == 0) {
                    strncpy(pps->particle_pipeline, value, MAX_PIPELINE_LENGTH - 1);
                    pps->particle_pipeline[MAX_PIPELINE_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_output_freq") == 0) {
                    pps->particle_output_freq = atoi(value);
                } else if (strcasecmp(key, "statistics_pipeline") == 0) {
                    strncpy(pps->statistics_pipeline, value, MAX_PIPELINE_LENGTH - 1);
                    pps->statistics_pipeline[MAX_PIPELINE_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "statistics_output_prefix") == 0) {
                    strncpy(pps->statistics_output_prefix, value, MAX_FILENAME_LENGTH - 1);
                    pps->statistics_output_prefix[MAX_FILENAME_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particleExt") == 0) {
                    strncpy(pps->particleExt, value, sizeof(pps->particleExt) - 1);
                    pps->particleExt[sizeof(pps->particleExt) - 1] = '\0';
                } else if (strcasecmp(key, "eulerianExt") == 0) {
                    strncpy(pps->eulerianExt, value, sizeof(pps->eulerianExt) - 1);
                    pps->eulerianExt[sizeof(pps->eulerianExt) - 1] = '\0';
                } else if (strcmp(key, "reference_ip") == 0) {pps->reference[0] = atoi(value);
                } else if (strcmp(key, "reference_jp") == 0) {pps->reference[1] = atoi(value);
                } else if (strcmp(key, "reference_kp") == 0) {pps->reference[2] = atoi(value);
                } else if (strcasecmp(key, "source_directory") == 0) {
                    strncpy(pps->source_dir, value, sizeof(pps->source_dir) - 1);
                    pps->source_dir[sizeof(pps->source_dir) - 1] = '\0';
                } else if (strcasecmp(key, "spectra_signature") == 0) {
                    /* Spectra are computed by the conductor's Python stage, not here. The
                       key exists so a change to the spectra recipe reaches the recipe
                       fingerprint that `--continue` compares; this executable has no use
                       for it and accepting it silently keeps the log free of a warning
                       that would appear on every post-processing run. */
                } else {
                    LOG_ALLOW(GLOBAL, LOG_WARNING, "Unknown key '%s' in post-processing config file. Ignoring.\n", key);
                }
                // Add parsing for pipeline, fields, etc. in later phases
            }
        }
        fclose(file);
    } else {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Could not open post-processing config file '%s'. Using defaults and command-line overrides.\n", configFile);
    }

    // --- 3. Parse Command-Line Options (overrides file settings and defaults) ---
    PetscOptionsGetInt(NULL, NULL, "-startTime", &pps->startTime, &startTimeSet);
    PetscOptionsGetInt(NULL, NULL, "-endTime", &pps->endTime, &endTimeSet);
    PetscOptionsGetInt(NULL, NULL, "-timeStep", &pps->timeStep, &timeStepSet);
    PetscOptionsGetBool(NULL, NULL, "-output_particles", &pps->outputParticles, NULL);

    if(pps->endTime==-1){
        pps->endTime = simCtx->StartStep + simCtx->StepsToRun; // Total steps if endTime is set to -1.
    }

    // If only startTime is given on command line, run for a single step
    if (startTimeSet && !endTimeSet) {
        pps->endTime = pps->startTime;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Post-processing configured to run from t=%d to t=%d with step %d. Particle output: %s.\n",
              pps->startTime, pps->endTime, pps->timeStep, pps->outputParticles ? "TRUE" : "FALSE");

    LOG_ALLOW(GLOBAL, LOG_INFO, "Process Pipeline: %s\n", pps->process_pipeline);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Instantaneous Output Fields: %s\n", pps->output_fields_instantaneous);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Output Prefix: %s\n", pps->output_prefix);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Output Prefix: %s\n", pps->particle_output_prefix);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Fields: %s\n", pps->particle_fields);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Pipeline: %s\n", pps->particle_pipeline);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Output Frequency: %d\n", pps->particle_output_freq);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Post input extensions: Eulerian='.%s', Particle='.%s'\n", pps->eulerianExt, pps->particleExt);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ParseScalingInformation"
/**
 * @brief Implementation of \ref ParseScalingInformation().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see ParseScalingInformation()
 */
PetscErrorCode ParseScalingInformation(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscBool      flg;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!simCtx) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx is NULL in ParseScalingInformation");

    // --- 1. Set default values to 1.0 ---
    // This represents a purely non-dimensional run if no scaling is provided.
    simCtx->scaling.L_ref   = 1.0;
    simCtx->scaling.U_ref   = 1.0;
    simCtx->scaling.rho_ref = 1.0;
    
    // --- 2. Read overrides from the command line / control file ---
    ierr = PetscOptionsGetReal(NULL, NULL, "-scaling_L_ref", &simCtx->scaling.L_ref, &flg); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-scaling_U_ref", &simCtx->scaling.U_ref, &flg); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-scaling_rho_ref", &simCtx->scaling.rho_ref, &flg); CHKERRQ(ierr);
    
    // --- 3. Calculate derived scaling factors ---
    // Check for division by zero to be safe, though U_ref should be positive.
    if (simCtx->scaling.U_ref <= 0.0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Reference velocity U_ref must be positive. Got %g", (double)simCtx->scaling.U_ref);
    }
    simCtx->scaling.P_ref = simCtx->scaling.rho_ref * simCtx->scaling.U_ref * simCtx->scaling.U_ref;

    // --- 4. Log the final, effective scales for verification ---
    LOG(GLOBAL, LOG_INFO, "---------------- Physical Scales Initialized -----------------\n");
    LOG(GLOBAL, LOG_INFO, "  L_ref: %.4f, U_ref: %.4f, rho_ref: %.4f, P_ref: %.4f\n",
        simCtx->scaling.L_ref, simCtx->scaling.U_ref, simCtx->scaling.rho_ref, simCtx->scaling.P_ref);
    LOG(GLOBAL, LOG_INFO, "--------------------------------------------------------------\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ReadDataFileToArray().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/io.h`.
 * @see ReadDataFileToArray()
 */
PetscInt ReadDataFileToArray(const char   *filename,
                        double      **data_out,
                        PetscInt          *Nout,
                        MPI_Comm      comm)
{
    /* STEP 0: Prepare local variables & log function entry */
    PetscMPIInt    rank, size;
    PetscErrorCode ierr;
    FILE  *fp = NULL;
    PetscInt    N   = 0;            /* number of lines/values read on rank 0 */
    double *array = NULL;      /* pointer to local array on each rank */
    PetscInt    fileExistsFlag = 0; /* 0 = doesn't exist, 1 = does exist */

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Start reading from file: %s\n",
              filename);

    /* Basic error checking: data_out, Nout must be non-null. */
    if (!filename || !data_out || !Nout) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Null pointer argument provided.\n");
        return 1;
    }

    /* Determine rank/size for coordinating I/O. */
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* STEP 1: On rank 0, check if file can be opened. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (fp) {
            fileExistsFlag = 1;
            fclose(fp);
        }
    }

    /* STEP 2: Broadcast file existence to all ranks. */
    // In ReadDataFileToArray:
    ierr = MPI_Bcast(&fileExistsFlag, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    if (!fileExistsFlag) {
        /* If file does not exist, log & return. */
        if (!rank) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "File '%s' not found.\n",
                      filename);
        }
        return 2;
    }

    /* STEP 3: Rank 0 re-opens and reads the file, counting lines, etc. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "File '%s' could not be opened for reading.\n",
                      filename);
            return 3;
        }

        /* (3a) Count lines first. */
        {
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                N++;
            }
        }

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "File '%s' has %d lines.\n",
                  filename, N);

        /* (3b) Allocate array on rank 0. */
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            fclose(fp);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "malloc failed for array.\n");
            return 4;
        }

        /* (3c) Rewind & read values into array. */
        rewind(fp);
        {
            PetscInt i = 0;
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                double val;
                if (sscanf(line, "%lf", &val) == 1) {
                    array[i++] = val;
                }
            }
        }
        fclose(fp);

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Successfully read %d values from '%s'.\n",
                  N, filename);
    }

    /* STEP 4: Broadcast the integer N to all ranks. */
    ierr = MPI_Bcast(&N, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    /* STEP 5: Each rank allocates an array to receive the broadcast if rank>0. */
    if (rank) {
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "malloc failed on rank %d.\n",
                      rank);
            return 5;
        }
    }

    /* STEP 6: Broadcast the actual data from rank 0 to all. */
    ierr = MPI_Bcast(array, N, MPI_DOUBLE, 0, comm); CHKERRQ(ierr);

    /* STEP 7: Assign outputs on all ranks. */
    *data_out = array;
    *Nout     = N;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Done. Provided array of length=%d to all ranks.\n",
              N);
    return 0; /* success */
}

/**
 * @brief Internal helper implementation: `ReadPositionsFromFile()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadPositionsFromFile(PetscInt timeIndex,
                                      UserCtx *user,
                                      double **coordsArray,
                                      PetscInt *Ncoords)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;
  Vec            coordsVec;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Creating coords Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

  // For example: "position" is the name of the coordinate data
  ierr = ReadFieldData(user, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), coordsVec, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "Error reading position data (ti=%d).\n",
              timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadPositions - Gathering coords Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(coordsVec, Ncoords, coordsArray);CHKERRQ(ierr);

  ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Successfully gathered coordinates. Ncoords=%d.\n", *Ncoords);
  PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `ReadFieldDataToRank0()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReadFieldDataToRank0(PetscInt timeIndex,
                                           const char *fieldName,
                                           UserCtx *user,
                                           double **scalarArray,
                                           PetscInt *Nscalars)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;
  Vec            fieldVec;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Creating field Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

  ierr = ReadFieldData(user, fieldName, fieldVec, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "Error reading field '%s' (ti=%d).\n",
              fieldName, timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Gathering field Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(fieldVec, Nscalars, scalarArray);CHKERRQ(ierr);

  ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Successfully gathered field '%s'. Nscalars=%d.\n",
            fieldName, *Nscalars);
  PetscFunctionReturn(0);
}
