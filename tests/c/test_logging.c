/**
 * @file test_logging.c
 * @brief C unit tests for runtime log-level, allow-list, conversion, and profiling helpers.
 */

#include "test_support.h"

#include "logging.h"
#include "interpolation.h"
#include "setup.h"

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
/**
 * @brief Asserts that one text file contains a required substring.
 */

static PetscErrorCode AssertFileContains(const char *path, const char *needle, const char *context)
{
    FILE *fp = NULL;
    long file_size = 0;
    char *buffer = NULL;

    PetscFunctionBeginUser;
    fp = fopen(path, "rb");
    if (!fp) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s' for assertion.", path);
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to seek '%s'.", path);
    }
    file_size = ftell(fp);
    if (file_size < 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to measure '%s'.", path);
    }
    if (fseek(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to rewind '%s'.", path);
    }

    PetscCall(PetscMalloc1((size_t)file_size + 1, &buffer));
    if (file_size > 0) {
        size_t bytes_read = fread(buffer, 1, (size_t)file_size, fp);
        if (bytes_read != (size_t)file_size) {
            fclose(fp);
            PetscCall(PetscFree(buffer));
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to read '%s'.", path);
        }
    }
    buffer[file_size] = '\0';
    fclose(fp);

    PetscCall(PicurvAssertBool((PetscBool)(strstr(buffer, needle) != NULL), context));
    PetscCall(PetscFree(buffer));
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that one text file does not contain an excluded substring.
 */

static PetscErrorCode AssertFileNotContains(const char *path, const char *needle, const char *context)
{
    FILE *fp = NULL;
    long file_size = 0;
    char *buffer = NULL;

    PetscFunctionBeginUser;
    fp = fopen(path, "rb");
    if (!fp) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s' for assertion.", path);
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to seek '%s'.", path);
    }
    file_size = ftell(fp);
    if (file_size < 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to measure '%s'.", path);
    }
    if (fseek(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to rewind '%s'.", path);
    }

    PetscCall(PetscMalloc1((size_t)file_size + 1, &buffer));
    if (file_size > 0) {
        size_t bytes_read = fread(buffer, 1, (size_t)file_size, fp);
        if (bytes_read != (size_t)file_size) {
            fclose(fp);
            PetscCall(PetscFree(buffer));
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to read '%s'.", path);
        }
    }
    buffer[file_size] = '\0';
    fclose(fp);

    PetscCall(PicurvAssertBool((PetscBool)(strstr(buffer, needle) == NULL), context));
    PetscCall(PetscFree(buffer));
    PetscFunctionReturn(0);
}

/**
 * @brief Reads the CSV header and the requested 1-based data row from a text file.
 */
static PetscErrorCode ReadCsvHeaderAndRow(const char *path,
                                          PetscInt row_index,
                                          char *header,
                                          size_t header_len,
                                          char *row,
                                          size_t row_len)
{
    FILE *fp = NULL;
    PetscInt current_row = 0;

    PetscFunctionBeginUser;
    PetscCheck(path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV path cannot be NULL.");
    PetscCheck(header != NULL && header_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV header buffer cannot be NULL or empty.");
    PetscCheck(row != NULL && row_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV row buffer cannot be NULL or empty.");
    PetscCheck(row_index >= 1, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "CSV row index must be >= 1.");

    fp = fopen(path, "r");
    PetscCheck(fp != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open CSV '%s'.", path);
    PetscCheck(fgets(header, (int)header_len, fp) != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV header missing in '%s'.", path);

    while (fgets(row, (int)row_len, fp) != NULL) {
        current_row++;
        if (current_row == row_index) {
            fclose(fp);
            PetscFunctionReturn(0);
        }
    }

    fclose(fp);
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "CSV '%s' does not contain data row %d.", path, (int)row_index);
}

/**
 * @brief Returns the zero-based column index for one CSV header field.
 */
static PetscErrorCode CsvFindColumnIndex(const char *header, const char *column_name, PetscInt *index_out)
{
    char local_header[4096];
    char *saveptr = NULL;
    char *token = NULL;
    PetscInt index = 0;

    PetscFunctionBeginUser;
    PetscCheck(header != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV header cannot be NULL.");
    PetscCheck(column_name != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV column name cannot be NULL.");
    PetscCheck(index_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV column index output cannot be NULL.");
    PetscCheck(strlen(header) < sizeof(local_header), PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "CSV header is too long for the local parser buffer.");

    PetscCall(PetscStrncpy(local_header, header, sizeof(local_header)));
    token = strtok_r(local_header, ",\r\n", &saveptr);
    while (token != NULL) {
        if (strcmp(token, column_name) == 0) {
            *index_out = index;
            PetscFunctionReturn(0);
        }
        token = strtok_r(NULL, ",\r\n", &saveptr);
        index++;
    }

    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "CSV column '%s' was not found.", column_name);
}

/**
 * @brief Extracts one CSV cell as text by header name.
 */
static PetscErrorCode CsvGetColumnText(const char *header,
                                       const char *row,
                                       const char *column_name,
                                       char *value,
                                       size_t value_len)
{
    char local_row[4096];
    char *saveptr = NULL;
    char *token = NULL;
    PetscInt target_index = -1;
    PetscInt index = 0;

    PetscFunctionBeginUser;
    PetscCheck(row != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV row cannot be NULL.");
    PetscCheck(value != NULL && value_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV value buffer cannot be NULL or empty.");
    PetscCheck(strlen(row) < sizeof(local_row), PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "CSV row is too long for the local parser buffer.");

    PetscCall(CsvFindColumnIndex(header, column_name, &target_index));
    PetscCall(PetscStrncpy(local_row, row, sizeof(local_row)));

    token = strtok_r(local_row, ",\r\n", &saveptr);
    while (token != NULL) {
        if (index == target_index) {
            PetscCall(PetscStrncpy(value, token, value_len));
            PetscFunctionReturn(0);
        }
        token = strtok_r(NULL, ",\r\n", &saveptr);
        index++;
    }

    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "CSV row is missing column '%s'.", column_name);
}

/**
 * @brief Extracts one CSV cell as an integer by header name.
 */
static PetscErrorCode CsvGetColumnInt(const char *header, const char *row, const char *column_name, PetscInt *value_out)
{
    char text[256];
    char *endptr = NULL;
    long parsed = 0;

    PetscFunctionBeginUser;
    PetscCheck(value_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV integer output cannot be NULL.");
    PetscCall(CsvGetColumnText(header, row, column_name, text, sizeof(text)));
    parsed = strtol(text, &endptr, 10);
    PetscCheck(endptr != text, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "CSV column '%s' did not contain an integer.", column_name);
    *value_out = (PetscInt)parsed;
    PetscFunctionReturn(0);
}

/**
 * @brief Extracts one CSV cell as a real by header name.
 */
static PetscErrorCode CsvGetColumnReal(const char *header, const char *row, const char *column_name, PetscReal *value_out)
{
    char text[256];
    char *endptr = NULL;
    double parsed = 0.0;

    PetscFunctionBeginUser;
    PetscCheck(value_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CSV real output cannot be NULL.");
    PetscCall(CsvGetColumnText(header, row, column_name, text, sizeof(text)));
    parsed = strtod(text, &endptr);
    PetscCheck(endptr != text, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "CSV column '%s' did not contain a real value.", column_name);
    *value_out = (PetscReal)parsed;
    PetscFunctionReturn(0);
}

/**
 * @brief Reads the column header and the requested 1-based data row from a
 *        pipe-delimited solution-convergence log file.
 *
 * Lines starting with '=' (banner) or '-' (separator) are skipped. The first
 * non-skipped line is treated as the column header; subsequent non-skipped
 * lines are data rows numbered from 1.
 */
static PetscErrorCode ReadLogHeaderAndRow(const char *path,
                                          PetscInt row_index,
                                          char *header,
                                          size_t header_len,
                                          char *row,
                                          size_t row_len)
{
    FILE      *fp = NULL;
    PetscInt   current_row = 0;
    char       line[8192];
    PetscBool  header_found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCheck(path != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log path cannot be NULL.");
    PetscCheck(header != NULL && header_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log header buffer cannot be NULL.");
    PetscCheck(row != NULL && row_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log row buffer cannot be NULL.");
    PetscCheck(row_index >= 1, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Log row index must be >= 1.");

    fp = fopen(path, "r");
    PetscCheck(fp != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open log '%s'.", path);

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '=' || line[0] == '-' || line[0] == '\n' || line[0] == '\r') continue;
        if (!header_found) {
            PetscCall(PetscStrncpy(header, line, header_len));
            header_found = PETSC_TRUE;
            continue;
        }
        current_row++;
        if (current_row == row_index) {
            PetscCall(PetscStrncpy(row, line, row_len));
            fclose(fp);
            PetscFunctionReturn(0);
        }
    }

    fclose(fp);
    if (!header_found) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Log '%s' has no column header.", path);
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Log '%s' does not contain data row %d.", path, (int)row_index);
}

/**
 * @brief Returns the zero-based column index for one pipe-delimited header field.
 */
static PetscErrorCode LogFindColumnIndex(const char *header, const char *column_name, PetscInt *index_out)
{
    char      local_header[8192];
    char     *saveptr = NULL;
    char     *token = NULL;
    PetscInt  index = 0;

    PetscFunctionBeginUser;
    PetscCheck(header != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log header cannot be NULL.");
    PetscCheck(column_name != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log column name cannot be NULL.");
    PetscCheck(index_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log column index output cannot be NULL.");

    PetscCall(PetscStrncpy(local_header, header, sizeof(local_header)));
    token = strtok_r(local_header, "|\r\n", &saveptr);
    while (token != NULL) {
        while (*token == ' ') token++;
        char *end = token + strlen(token) - 1;
        while (end > token && (*end == ' ' || *end == '\r' || *end == '\n')) end--;
        *(end + 1) = '\0';
        if (strcmp(token, column_name) == 0) {
            *index_out = index;
            PetscFunctionReturn(0);
        }
        token = strtok_r(NULL, "|\r\n", &saveptr);
        index++;
    }
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Log column '%s' was not found.", column_name);
}

/**
 * @brief Extracts one pipe-delimited log cell as text by header name.
 */
static PetscErrorCode LogGetColumnText(const char *header,
                                       const char *row,
                                       const char *column_name,
                                       char *value,
                                       size_t value_len)
{
    char      local_row[8192];
    char     *saveptr = NULL;
    char     *token = NULL;
    PetscInt  target_index = -1;
    PetscInt  index = 0;

    PetscFunctionBeginUser;
    PetscCheck(row != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log row cannot be NULL.");
    PetscCheck(value != NULL && value_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log value buffer cannot be NULL.");

    PetscCall(LogFindColumnIndex(header, column_name, &target_index));
    PetscCall(PetscStrncpy(local_row, row, sizeof(local_row)));

    token = strtok_r(local_row, "|\r\n", &saveptr);
    while (token != NULL) {
        if (index == target_index) {
            while (*token == ' ') token++;
            char *end = token + strlen(token) - 1;
            while (end > token && (*end == ' ' || *end == '\r' || *end == '\n')) end--;
            *(end + 1) = '\0';
            PetscCall(PetscStrncpy(value, token, value_len));
            PetscFunctionReturn(0);
        }
        token = strtok_r(NULL, "|\r\n", &saveptr);
        index++;
    }
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Log row is missing column '%s'.", column_name);
}

/**
 * @brief Extracts one pipe-delimited log cell as an integer by header name.
 */
static PetscErrorCode LogGetColumnInt(const char *header, const char *row, const char *column_name, PetscInt *value_out)
{
    char  text[256];
    char *endptr = NULL;
    long  parsed = 0;

    PetscFunctionBeginUser;
    PetscCheck(value_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log integer output cannot be NULL.");
    PetscCall(LogGetColumnText(header, row, column_name, text, sizeof(text)));
    parsed = strtol(text, &endptr, 10);
    PetscCheck(endptr != text, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Log column '%s' did not contain an integer.", column_name);
    *value_out = (PetscInt)parsed;
    PetscFunctionReturn(0);
}

/**
 * @brief Extracts one pipe-delimited log cell as a real by header name.
 */
static PetscErrorCode LogGetColumnReal(const char *header, const char *row, const char *column_name, PetscReal *value_out)
{
    char   text[256];
    char  *endptr = NULL;
    double parsed = 0.0;

    PetscFunctionBeginUser;
    PetscCheck(value_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Log real output cannot be NULL.");
    PetscCall(LogGetColumnText(header, row, column_name, text, sizeof(text)));
    parsed = strtod(text, &endptr);
    PetscCheck(endptr != text, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Log column '%s' did not contain a real value.", column_name);
    *value_out = (PetscReal)parsed;
    PetscFunctionReturn(0);
}

typedef PetscErrorCode (*CapturedLoggingFn)(UserCtx *user, SimCtx *simCtx, void *ctx);

typedef struct AnatomyCaptureCtx {
    const char *field_name;
    const char *stage_name;
} AnatomyCaptureCtx;

/**
 * @brief Captures stdout emitted by one logging helper into a temporary file-backed buffer.
 */
static PetscErrorCode CaptureLoggingOutput(UserCtx *user,
                                           SimCtx *simCtx,
                                           CapturedLoggingFn fn,
                                           void *ctx,
                                           char *captured,
                                           size_t captured_len)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char capture_path[PETSC_MAX_PATH_LEN];
    FILE *capture_file = NULL;
    int saved_stdout = -1;
    int capture_fd = -1;
    size_t bytes_read = 0;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PetscCheck(fn != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Capture callback cannot be NULL.");
    PetscCheck(captured != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Capture buffer cannot be NULL.");
    PetscCheck(captured_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Capture buffer must be non-empty.");

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(capture_path, sizeof(capture_path), "%s/logging.out", tmpdir));

    fflush(stdout);
    saved_stdout = dup(STDOUT_FILENO);
    PetscCheck(saved_stdout >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup(STDOUT_FILENO) failed.");
    capture_fd = open(capture_path, O_CREAT | O_TRUNC | O_WRONLY, 0600);
    PetscCheck(capture_fd >= 0, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open capture file '%s'.", capture_path);
    PetscCheck(dup2(capture_fd, STDOUT_FILENO) >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup2() failed while redirecting stdout.");
    close(capture_fd);
    capture_fd = -1;

    ierr = fn(user, simCtx, ctx);
    fflush(stdout);
    PetscCheck(dup2(saved_stdout, STDOUT_FILENO) >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup2() failed while restoring stdout.");
    close(saved_stdout);
    saved_stdout = -1;
    PetscCall(ierr);

    capture_file = fopen(capture_path, "r");
    PetscCheck(capture_file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to read capture file '%s'.", capture_path);
    bytes_read = fread(captured, 1, captured_len - 1, capture_file);
    captured[bytes_read] = '\0';
    fclose(capture_file);
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Adapts `LOG_PARTICLE_FIELDS()` to the generic stdout-capture callback shape.
 */
static PetscErrorCode InvokeParticleFieldLog(UserCtx *user, SimCtx *simCtx, void *ctx)
{
    PetscInt print_interval = *((PetscInt *)ctx);

    PetscFunctionBeginUser;
    (void)simCtx;
    PetscCall(LOG_PARTICLE_FIELDS(user, print_interval));
    PetscFunctionReturn(0);
}

/**
 * @brief Adapts `EmitParticleConsoleSnapshot()` to the generic stdout-capture callback shape.
 */
static PetscErrorCode InvokeParticleConsoleSnapshot(UserCtx *user, SimCtx *simCtx, void *ctx)
{
    PetscInt step = *((PetscInt *)ctx);

    PetscFunctionBeginUser;
    PetscCall(EmitParticleConsoleSnapshot(user, simCtx, step));
    PetscFunctionReturn(0);
}

/**
 * @brief Adapts `LOG_FIELD_ANATOMY()` to the generic stdout-capture callback shape.
 */
static PetscErrorCode InvokeFieldAnatomyLog(UserCtx *user, SimCtx *simCtx, void *ctx)
{
    const AnatomyCaptureCtx *anatomy_ctx = (const AnatomyCaptureCtx *)ctx;

    PetscFunctionBeginUser;
    (void)simCtx;
    PetscCall(LOG_FIELD_ANATOMY(user, anatomy_ctx->field_name, anatomy_ctx->stage_name));
    PetscFunctionReturn(0);
}

/**
 * @brief Creates a small particle-bearing runtime fixture used by logging tests.
 */
static PetscErrorCode SeedLoggingParticleFixture(SimCtx **simCtx_out, UserCtx **user_out)
{
    PetscReal *positions = NULL;
    PetscReal *velocities = NULL;
    PetscReal *weights = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;
    PetscReal *psi = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(simCtx_out, user_out, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(*user_out, 2, "ske"));
    (*simCtx_out)->np = 2;
    (*simCtx_out)->particleConsoleOutputFreq = 2;
    (*simCtx_out)->LoggingFrequency = 1;

    PetscCall(DMSwarmGetField((*user_out)->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "velocity", NULL, NULL, (void **)&velocities));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "Psi", NULL, NULL, (void **)&psi));

    positions[0] = 0.25; positions[1] = 0.50; positions[2] = 0.75;
    positions[3] = 0.50; positions[4] = 0.50; positions[5] = 0.50;
    velocities[0] = 1.0; velocities[1] = 2.0; velocities[2] = 3.0;
    velocities[3] = 4.0; velocities[4] = 5.0; velocities[5] = 6.0;
    weights[0] = 0.2; weights[1] = 0.3; weights[2] = 0.4;
    weights[3] = 0.5; weights[4] = 0.5; weights[5] = 0.5;
    cell_ids[0] = 0; cell_ids[1] = 0; cell_ids[2] = 0;
    cell_ids[3] = 1; cell_ids[4] = 1; cell_ids[5] = 1;
    status[0] = ACTIVE_AND_LOCATED;
    status[1] = ACTIVE_AND_LOCATED;
    psi[0] = 1.0;
    psi[1] = 3.0;

    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "velocity", NULL, NULL, (void **)&velocities));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "position", NULL, NULL, (void **)&positions));
    PetscFunctionReturn(0);
}

/**
 * @brief Fills one Cartesian velocity field with a uniform constant state.
 */
static PetscErrorCode SetUniformVelocityField(UserCtx *user, Vec field, PetscReal ux, PetscReal uy, PetscReal uz)
{
    Cmpnts ***arr = NULL;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    PetscCheck(field != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Velocity field cannot be NULL.");

    PetscCall(DMDAVecGetArray(user->fda, field, &arr));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                arr[k][j][i].x = ux;
                arr[k][j][i].y = uy;
                arr[k][j][i].z = uz;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, field, &arr));
    PetscFunctionReturn(0);
}

/**
 * @brief Fills one scalar field with a uniform constant state.
 */
static PetscErrorCode SetUniformScalarField(UserCtx *user, Vec field, PetscReal value)
{
    PetscReal ***arr = NULL;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    PetscCheck(field != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Scalar field cannot be NULL.");

    PetscCall(DMDAVecGetArray(user->da, field, &arr));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                arr[k][j][i] = value;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, field, &arr));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests string-conversion helpers for configured enums and unknown values.
 */

static PetscErrorCode TestStringConversionHelpers(void)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool(strcmp(BCFaceToString(BC_FACE_NEG_X), "-Xi (I-Min)") == 0,
                               "BCFaceToString should report the negative-x face"));
    PetscCall(PicurvAssertBool(strcmp(FieldInitializationToString(0), "Zero") == 0,
                               "FieldInitializationToString should report the zero mode"));
    PetscCall(PicurvAssertBool(strcmp(FieldInitializationToString(99), "Unknown Field Initialization") == 0,
                               "FieldInitializationToString should reject unknown selectors"));
    PetscCall(PicurvAssertBool(strcmp(ParticleInitializationToString(PARTICLE_INIT_VOLUME), "Volume") == 0,
                               "ParticleInitializationToString should report the volume mode"));
    PetscCall(PicurvAssertBool(strcmp(LESModelToString(CONSTANT_SMAGORINSKY), "Constant Smagorinsky") == 0,
                               "LESModelToString should report the constant model"));
    PetscCall(PicurvAssertBool(strcmp(MomentumSolverTypeToString(MOMENTUM_SOLVER_EXPLICIT_RK), "Explicit 4 stage Runge-Kutta ") == 0,
                               "MomentumSolverTypeToString should report the explicit solver"));
    PetscCall(PicurvAssertBool(strcmp(BCTypeToString(PERIODIC), "PERIODIC") == 0,
                               "BCTypeToString should report periodic boundaries"));
    PetscCall(PicurvAssertBool(strcmp(BCHandlerTypeToString(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX), "constant flux") == 0,
                               "BCHandlerTypeToString should report the driven periodic handler"));
    PetscCall(PicurvAssertBool(strcmp(ParticleLocationStatusToString(LOST), "LOST") == 0,
                               "ParticleLocationStatusToString should report LOST state"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that log level selection honors the environment variable.
 */

static PetscErrorCode TestGetLogLevelFromEnvironment(void)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertIntEqual(LOG_INFO, get_log_level(),
                                   "get_log_level should honor LOG_LEVEL=INFO in this test binary"));
    PetscCall(print_log_level());
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the function allow-list filter used by the logging layer.
 */

static PetscErrorCode TestAllowedFunctionsFilter(void)
{
    const char *allow_list[] = {"ComputeSpecificKE", "WriteEulerianFile"};

    PetscFunctionBeginUser;
    set_allowed_functions(allow_list, 2);
    PetscCall(PicurvAssertBool(is_function_allowed("ComputeSpecificKE"),
                               "Allowed list should include ComputeSpecificKE"));
    PetscCall(PicurvAssertBool((PetscBool)!is_function_allowed("UnlistedFunction"),
                               "Allowed list should exclude unknown function names"));

    set_allowed_functions(NULL, 0);
    PetscCall(PicurvAssertBool(is_function_allowed("AnyFunction"),
                               "Empty allow-list should permit all functions"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests periodic particle console snapshot enablement and cadence.
 */

static PetscErrorCode TestParticleConsoleSnapshotCadence(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    simCtx.np = 32;
    simCtx.particleConsoleOutputFreq = 4;

    PetscCall(PicurvAssertBool(IsParticleConsoleSnapshotEnabled(&simCtx),
                               "Particle snapshot contract should be enabled when particles and cadence are configured"));
    PetscCall(PicurvAssertBool(ShouldEmitPeriodicParticleConsoleSnapshot(&simCtx, 8),
                               "Snapshot should emit on cadence-aligned completed steps"));
    PetscCall(PicurvAssertBool((PetscBool)!ShouldEmitPeriodicParticleConsoleSnapshot(&simCtx, 7),
                               "Snapshot should not emit off-cadence"));

    simCtx.particleConsoleOutputFreq = 0;
    PetscCall(PicurvAssertBool((PetscBool)!IsParticleConsoleSnapshotEnabled(&simCtx),
                               "Zero cadence should disable periodic particle snapshots"));
    PetscCall(PicurvAssertBool((PetscBool)!ShouldEmitPeriodicParticleConsoleSnapshot(NULL, 4),
                               "NULL SimCtx should never emit periodic snapshots"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests logging-side file parsing, helper formatting, and progress utilities.
 */

static PetscErrorCode TestLoggingFileParsingAndFormattingHelpers(void)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char allow_path[PETSC_MAX_PATH_LEN];
    char dual_log_path[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;
    char **funcs = NULL;
    PetscInt nfuncs = 0;
    Cell cell;
    PetscReal distances[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    DualMonitorCtx *monctx = NULL;
    void *ctx = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&cell, sizeof(cell)));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(allow_path, sizeof(allow_path), "%s/allowed_functions.txt", tmpdir));
    PetscCall(PetscSNPrintf(dual_log_path, sizeof(dual_log_path), "%s/dual-monitor.log", tmpdir));

    file = fopen(allow_path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create allow-list file '%s'.", allow_path);
    fputs("  ComputeSpecificKE  \n", file);
    fputs("# comment-only line\n", file);
    for (PetscInt i = 0; i < 17; ++i) {
        fprintf(file, "Helper_%02d   # trailing comment\n", (int)i);
    }
    fclose(file);
    file = NULL;

    PetscCall(LoadAllowedFunctionsFromFile(allow_path, &funcs, &nfuncs));
    PetscCall(PicurvAssertIntEqual(18, nfuncs, "LoadAllowedFunctionsFromFile should trim comments and keep all identifiers"));
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(funcs[0], "ComputeSpecificKE") == 0),
                               "LoadAllowedFunctionsFromFile should trim leading and trailing whitespace"));
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(funcs[17], "Helper_16") == 0),
                               "LoadAllowedFunctionsFromFile should grow past the initial pointer capacity"));
    PetscCall(FreeAllowedFunctions(funcs, nfuncs));

    for (PetscInt i = 0; i < 8; ++i) {
        cell.vertices[i].x = (PetscReal)i;
        cell.vertices[i].y = (PetscReal)(i + 1);
        cell.vertices[i].z = (PetscReal)(i + 2);
    }
    PetscCall(LOG_CELL_VERTICES(&cell, 0));
    PetscCall(LOG_FACE_DISTANCES(distances));

    PetscCall(PetscCalloc1(1, &monctx));
    monctx->file_handle = fopen(dual_log_path, "w");
    PetscCheck(monctx->file_handle != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create dual-monitor log '%s'.", dual_log_path);
    ctx = monctx;
    PetscCall(DualMonitorDestroy(&ctx));
    PetscCall(PicurvAssertBool((PetscBool)(ctx == NULL),
                               "DualMonitorDestroy should clear the caller-owned context pointer"));

    PrintProgressBar(0, 0, 4, 0.10);
    PrintProgressBar(3, 0, 4, 0.40);
    PrintProgressBar(0, 0, 0, 0.00);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests continuity, min/max, and anatomy logging helpers on minimal runtime fixtures.
 */

static PetscErrorCode TestLoggingContinuityAndFieldDiagnostics(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char continuity_path[PETSC_MAX_PATH_LEN];
    PetscReal ***p = NULL;
    Cmpnts ***ucat = NULL;
    Cmpnts ***ucont = NULL;
    PetscErrorCode ierr_minmax = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));

    simCtx->StartStep = 0;
    simCtx->step = 1;
    simCtx->MaxDiv = 1.25;
    simCtx->MaxDivx = 1;
    simCtx->MaxDivy = 2;
    simCtx->MaxDivz = 3;
    simCtx->MaxDivFlatArg = 17;
    simCtx->summationRHS = 8.5;
    simCtx->FluxInSum = 5.0;
    simCtx->FluxOutSum = 3.25;
    PetscCall(LOG_CONTINUITY_METRICS(user));

    simCtx->step = 2;
    simCtx->MaxDiv = 0.75;
    simCtx->summationRHS = 4.5;
    simCtx->FluxInSum = 2.5;
    simCtx->FluxOutSum = 1.0;
    PetscCall(LOG_CONTINUITY_METRICS(user));

    PetscCall(PetscSNPrintf(continuity_path, sizeof(continuity_path), "%s/Continuity_Metrics.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(continuity_path, "continuity metrics log should be written"));
    PetscCall(AssertFileContains(continuity_path, "Timestep", "continuity metrics log should include the header"));
    PetscCall(AssertFileContains(continuity_path, "([3][2][1] = 17)", "continuity metrics log should include the divergence location"));
    PetscCall(AssertFileContains(continuity_path, "2          | 0", "continuity metrics log should append later timesteps"));

    PetscCall(DMDAVecGetArray(user->da, user->P, &p));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                p[k][j][i] = (PetscReal)(i + j + k);
                ucat[k][j][i].x = (PetscReal)i;
                ucat[k][j][i].y = (PetscReal)(-j);
                ucat[k][j][i].z = (PetscReal)(2 * k);
                ucont[k][j][i].x = (PetscReal)(10 + i);
                ucont[k][j][i].y = (PetscReal)(20 + j);
                ucont[k][j][i].z = (PetscReal)(30 + k);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));

    PetscCall(DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));

    PetscCall(LOG_FIELD_MIN_MAX(user, "P"));
    PetscCall(LOG_FIELD_MIN_MAX(user, "Ucat"));
    PetscCall(LOG_FIELD_MIN_MAX(user, "Coordinates"));
    PetscCall(LOG_FIELD_MIN_MAX(user, "Ucont"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_minmax = LOG_FIELD_MIN_MAX(user, "NotARealField");
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_minmax != 0),
                               "LOG_FIELD_MIN_MAX should reject unknown field names"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests interpolation-error logging against an analytically matched particle field.
 */

static PetscErrorCode TestInterpolationErrorLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal (*pos_arr)[3] = NULL;
    PetscReal (*vel_arr)[3] = NULL;
    Vec position_vec = NULL;
    Vec analytical_vec = NULL;
    const PetscScalar *analytical_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));
    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "TGV3D", sizeof(simCtx->AnalyticalSolutionType)));
    simCtx->ren = 1.0;
    simCtx->ti = 0.0;

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 0.5 * PETSC_PI; pos_arr[0][1] = 0.0;          pos_arr[0][2] = 0.0;
    pos_arr[1][0] = 0.0;           pos_arr[1][1] = 0.5 * PETSC_PI; pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(DMSwarmCreateGlobalVectorFromField(user->swarm, "position", &position_vec));
    PetscCall(VecDuplicate(position_vec, &analytical_vec));
    PetscCall(VecCopy(position_vec, analytical_vec));
    PetscCall(SetAnalyticalSolutionForParticles(analytical_vec, simCtx));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    PetscCall(VecGetArrayRead(analytical_vec, &analytical_arr));
    for (PetscInt particle = 0; particle < 2; ++particle) {
        vel_arr[particle][0] = PetscRealPart(analytical_arr[3 * particle + 0]);
        vel_arr[particle][1] = PetscRealPart(analytical_arr[3 * particle + 1]);
        vel_arr[particle][2] = PetscRealPart(analytical_arr[3 * particle + 2]);
    }
    PetscCall(VecRestoreArrayRead(analytical_vec, &analytical_arr));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(VecDestroy(&analytical_vec));
    PetscCall(DMSwarmDestroyGlobalVectorFromField(user->swarm, "position", &position_vec));
    PetscCall(LOG_INTERPOLATION_ERROR(user));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests file-backed scatter metrics logging against a fully occupied constant field.
 */

static PetscErrorCode TestScatterMetricsLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char metrics_path[PETSC_MAX_PATH_LEN];
    PetscReal *positions = NULL;
    PetscReal *psi = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;
    PetscInt particle = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 27, "ske"));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->np = 27;
    simCtx->step = 3;
    simCtx->ti = 0.3;
    simCtx->verificationScalar.enabled = PETSC_TRUE;
    PetscCall(PetscStrncpy(simCtx->verificationScalar.mode,
                           "analytical",
                           sizeof(simCtx->verificationScalar.mode)));
    PetscCall(PetscStrncpy(simCtx->verificationScalar.profile,
                           "CONSTANT",
                           sizeof(simCtx->verificationScalar.profile)));
    simCtx->verificationScalar.value = 2.0;

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));

    for (PetscInt k = 0; k < 3; ++k) {
        for (PetscInt j = 0; j < 3; ++j) {
            for (PetscInt i = 0; i < 3; ++i) {
                positions[3 * particle + 0] = (i + 0.5) / 4.0;
                positions[3 * particle + 1] = (j + 0.5) / 4.0;
                positions[3 * particle + 2] = (k + 0.5) / 4.0;
                cell_ids[3 * particle + 0] = i;
                cell_ids[3 * particle + 1] = j;
                cell_ids[3 * particle + 2] = k;
                status[particle] = ACTIVE_AND_LOCATED;
                psi[particle] = 2.0;
                ++particle;
            }
        }
    }

    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));

    PetscCall(ScatterAllParticleFieldsToEulerFields(user));
    PetscCall(LOG_SCATTER_METRICS(user));

    PetscCall(PetscSNPrintf(metrics_path, sizeof(metrics_path), "%s/scatter_metrics.csv", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(metrics_path, "LOG_SCATTER_METRICS should write scatter_metrics.csv"));
    PetscCall(AssertFileContains(metrics_path, "relative_L2_error",
                                 "Scatter metrics CSV header should include relative_L2_error"));
    PetscCall(AssertFileContains(metrics_path,
                                 "3,3.000000e-01,27,27,27,1.000000e+00,1.000000e+00,5.400000e+01,5.400000e+01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00",
                                 "Scatter metrics CSV should record the expected constant-field zero-error row"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests stdout particle-table logging on a production-like swarm fixture.
 */

static PetscErrorCode TestParticleFieldTableLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt print_interval = 1;
    char captured[8192];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(CaptureLoggingOutput(user, simCtx, InvokeParticleFieldLog, &print_interval, captured, sizeof(captured)));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Position (x,y,z)") != NULL),
                               "LOG_PARTICLE_FIELDS should print the particle table header"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Weights (a1,a2,a3)") != NULL),
                               "LOG_PARTICLE_FIELDS should print the weight-column header"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests console snapshot logging against the public periodic-snapshot helper.
 */

static PetscErrorCode TestParticleConsoleSnapshotLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt step = 4;
    char captured[8192];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(CaptureLoggingOutput(user, simCtx, InvokeParticleConsoleSnapshot, &step, captured, sizeof(captured)));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Particle states at step 4") != NULL),
                               "EmitParticleConsoleSnapshot should print the step banner"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Position (x,y,z)") != NULL),
                               "EmitParticleConsoleSnapshot should reuse the particle table output"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests file-backed particle metrics logging after derived metrics are computed.
 */

static PetscErrorCode TestParticleMetricsLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char metrics_path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->StartStep = 0;
    simCtx->step = 1;
    simCtx->particlesLostLastStep = 1;
    simCtx->particlesLostCumulative = 7;
    simCtx->particlesMigratedLastStep = 2;
    simCtx->migrationPassesLastStep = 3;

    PetscCall(CalculateParticleCountPerCell(user));
    PetscCall(CalculateAdvancedParticleMetrics(user));
    PetscCall(LOG_PARTICLE_METRICS(user, "Timestep Metrics"));

    PetscCall(PetscSNPrintf(metrics_path, sizeof(metrics_path), "%s/Particle_Metrics.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(metrics_path, "LOG_PARTICLE_METRICS should write Particle_Metrics.log"));
    PetscCall(AssertFileContains(metrics_path, "Timestep Metrics", "Particle metrics log should include the caller-provided stage label"));
    PetscCall(AssertFileContains(metrics_path, "Occupied Cells", "Particle metrics log should include the metrics table header"));
    PetscCall(AssertFileContains(metrics_path, "Lost Total", "Particle metrics log should include the cumulative-loss column"));
    PetscCall(AssertFileContains(metrics_path, "| 1          | 7          | 2", "Particle metrics log should record both per-step and cumulative loss values"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests file-backed search metrics logging with the compact CSV contract.
 */

static PetscErrorCode TestSearchMetricsLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char metrics_path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->step = 2;
    simCtx->ti = 0.2;
    simCtx->particlesLostLastStep = 1;
    simCtx->particlesLostCumulative = 7;
    simCtx->particlesMigratedLastStep = 2;
    simCtx->migrationPassesLastStep = 3;
    simCtx->particleLoadImbalance = 1.5;
    simCtx->searchMetrics.searchAttempts = 4;
    simCtx->searchMetrics.searchPopulation = 2;
    simCtx->searchMetrics.searchLocatedCount = 1;
    simCtx->searchMetrics.searchLostCount = 1;
    simCtx->searchMetrics.traversalStepsSum = 10;
    simCtx->searchMetrics.reSearchCount = 2;
    simCtx->searchMetrics.maxTraversalSteps = 6;
    simCtx->searchMetrics.maxTraversalFailCount = 1;
    simCtx->searchMetrics.tieBreakCount = 1;
    simCtx->searchMetrics.boundaryClampCount = 2;
    simCtx->searchMetrics.bboxGuessSuccessCount = 3;
    simCtx->searchMetrics.bboxGuessFallbackCount = 1;
    simCtx->searchMetrics.maxParticlePassDepth = 3;

    PetscCall(LOG_SEARCH_METRICS(user));

    PetscCall(PetscSNPrintf(metrics_path, sizeof(metrics_path), "%s/search_metrics.csv", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(metrics_path, "LOG_SEARCH_METRICS should write search_metrics.csv"));
    PetscCall(AssertFileContains(metrics_path, "search_attempts", "Search metrics CSV header should include search_attempts"));
    PetscCall(AssertFileContains(metrics_path, "search_population", "Search metrics CSV header should include search_population"));
    PetscCall(AssertFileContains(metrics_path, "max_particle_pass_depth", "Search metrics CSV header should include max_particle_pass_depth"));
    PetscCall(AssertFileContains(metrics_path, "search_work_index", "Search metrics CSV header should include search_work_index"));
    PetscCall(AssertFileContains(metrics_path, "re_search_fraction", "Search metrics CSV header should include re_search_fraction"));
    PetscCall(AssertFileContains(metrics_path, "lost_cumulative", "Search metrics CSV header should include the cumulative-loss column"));
    PetscCall(AssertFileContains(metrics_path, "2,2.000000e-01,2,1,7,2,3,4,2.500000e+00,6,1,2,3,1,3,1.500000e+00,2,1,1,10,2,1,5.000000e-01,5.000000e+00,1.000000e+00", "Search metrics CSV should record the V2 raw and derived search metrics"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests stdout field-anatomy logging on the corrected production-like DM fixture.
 */

static PetscErrorCode TestFieldAnatomyLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char captured[8192];
    AnatomyCaptureCtx anatomy_ctx = {"P", "unit-test"};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->P, 7.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(CaptureLoggingOutput(user, simCtx, InvokeFieldAnatomyLog, &anatomy_ctx, captured, sizeof(captured)));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Field Anatomy Log: [P]") != NULL),
                               "LOG_FIELD_ANATOMY should print the requested field name"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Layout: [Cell-Centered]") != NULL),
                               "LOG_FIELD_ANATOMY should report the inferred data layout"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests profiling helper lifecycle logging for timestep and final-summary outputs.
 */

static PetscErrorCode TestProfilingLifecycleHelpers(void)
{
    SimCtx simCtx;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char timestep_path[PETSC_MAX_PATH_LEN];
    char summary_path[PETSC_MAX_PATH_LEN];
    static char selected_name[] = "FlowSolver";
    char *selected_funcs[] = {selected_name};

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx.log_dir, tmpdir, sizeof(simCtx.log_dir)));
    PetscCall(PetscStrncpy(simCtx.profilingTimestepMode, "selected", sizeof(simCtx.profilingTimestepMode)));
    PetscCall(PetscStrncpy(simCtx.profilingTimestepFile, "Profiling_Timestep_Summary.csv", sizeof(simCtx.profilingTimestepFile)));
    simCtx.rank = 0;
    simCtx.exec_mode = EXEC_MODE_SOLVER;
    simCtx.StartStep = 0;
    simCtx.nProfilingSelectedFuncs = 1;
    simCtx.profilingSelectedFuncs = selected_funcs;
    simCtx.profilingFinalSummary = PETSC_TRUE;

    PetscCall(ProfilingInitialize(&simCtx));

    _ProfilingStart("FlowSolver");
    _ProfilingEnd("FlowSolver");
    _ProfilingStart("UnselectedHelper");
    _ProfilingEnd("UnselectedHelper");
    PetscCall(ProfilingLogTimestepSummary(&simCtx, 1));

    _ProfilingStart("FlowSolver");
    _ProfilingEnd("FlowSolver");
    PetscCall(ProfilingResetTimestepCounters());
    PetscCall(ProfilingLogTimestepSummary(&simCtx, 2));
    PetscCall(ProfilingFinalize(&simCtx));

    PetscCall(PetscSNPrintf(timestep_path, sizeof(timestep_path), "%s/%s", simCtx.log_dir, simCtx.profilingTimestepFile));
    PetscCall(PetscSNPrintf(summary_path, sizeof(summary_path), "%s/ProfilingSummary_Solver.log", simCtx.log_dir));
    PetscCall(PicurvAssertFileExists(timestep_path, "profiling timestep summary should be written"));
    PetscCall(PicurvAssertFileExists(summary_path, "profiling final summary should be written"));
    PetscCall(AssertFileContains(timestep_path, "step,function,calls,step_time_s",
                                 "profiling timestep summary should contain the CSV header"));
    PetscCall(AssertFileContains(timestep_path, "1,FlowSolver,1,",
                                 "profiling timestep summary should log selected functions"));
    PetscCall(AssertFileNotContains(timestep_path, "UnselectedHelper",
                                    "profiling timestep summary should omit unselected functions in selected mode"));
    PetscCall(AssertFileContains(summary_path, "FINAL PROFILING SUMMARY",
                                 "profiling final summary should include its table banner"));
    PetscCall(AssertFileContains(summary_path, "FlowSolver",
                                 "profiling final summary should include selected functions"));
    PetscCall(AssertFileContains(summary_path, "UnselectedHelper",
                                 "profiling final summary should include total-time entries for unselected functions"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests steady solution-convergence log output, IBM masking, and gauge-invariant pressure drift.
 */
static PetscErrorCode TestSolutionConvergenceSteadyLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char log_path[PETSC_MAX_PATH_LEN];
    char header[4096];
    char row1[4096];
    char row2[4096];
    char mode[128];
    Cmpnts ***ucat = NULL;
    PetscReal ***nvert = NULL;
    PetscReal mean_speed = NAN;
    PetscReal mean_speed_ref = NAN;
    PetscReal mean_speed_abs = NAN;
    PetscReal p_abs = NAN;
    PetscInt has_reference_1 = 0;
    PetscInt has_reference_2 = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->solutionConvergenceMode = SOLUTION_CONVERGENCE_STEADY_DETERMINISTIC;

    PetscCall(InitializeSolutionConvergenceState(simCtx));

    simCtx->step = 1;
    simCtx->ti = 0.1;
    PetscCall(SetUniformVelocityField(user, user->Ucat, 1.0, 0.0, 0.0));
    PetscCall(SetUniformScalarField(user, user->P, 11.0));
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    nvert[1][1][1] = 1.0;
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    ucat[1][1][1].x = 999.0;
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(LOG_SOLUTION_CONVERGENCE(simCtx));

    simCtx->step = 2;
    simCtx->ti = 0.2;
    PetscCall(SetUniformVelocityField(user, user->Ucat, 1.0, 0.0, 0.0));
    PetscCall(SetUniformVelocityField(user, user->Ucat_o, 0.5, 0.0, 0.0));
    PetscCall(SetUniformScalarField(user, user->P, 11.0));
    PetscCall(SetUniformScalarField(user, user->P_o, 7.0));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    ucat[1][1][1].x = 999.0;
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(LOG_SOLUTION_CONVERGENCE(simCtx));

    PetscCall(PetscSNPrintf(log_path, sizeof(log_path), "%s/solution_convergence.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(log_path, "steady solution-convergence logging should write the log"));
    PetscCall(ReadLogHeaderAndRow(log_path, 1, header, sizeof(header), row1, sizeof(row1)));
    PetscCall(ReadLogHeaderAndRow(log_path, 2, header, sizeof(header), row2, sizeof(row2)));
    PetscCall(LogGetColumnText(header, row1, "mode", mode, sizeof(mode)));
    PetscCall(LogGetColumnInt(header, row1, "ref", &has_reference_1));
    PetscCall(LogGetColumnInt(header, row2, "ref", &has_reference_2));
    PetscCall(LogGetColumnReal(header, row2, "mean_speed", &mean_speed));
    PetscCall(LogGetColumnReal(header, row2, "spd_ref", &mean_speed_ref));
    PetscCall(LogGetColumnReal(header, row2, "spd_abs", &mean_speed_abs));
    PetscCall(LogGetColumnReal(header, row2, "p_abs_l2", &p_abs));

    PetscCall(PicurvAssertBool((PetscBool)(strcmp(mode, "steady_deterministic") == 0),
                               "steady solution convergence row should record the mode name"));
    PetscCall(PicurvAssertIntEqual(0, has_reference_1,
                                   "the first steady solution-convergence row should be a warmup row"));
    PetscCall(PicurvAssertIntEqual(1, has_reference_2,
                                   "steady solution convergence should compare against the previous solved step"));
    PetscCall(PicurvAssertRealNear(1.0, mean_speed, 1.0e-12,
                                   "steady solution convergence should mask IBM-marked solid cells"));
    PetscCall(PicurvAssertRealNear(0.5, mean_speed_ref, 1.0e-12,
                                   "steady solution convergence should report the previous-step mean speed"));
    PetscCall(PicurvAssertRealNear(0.5, mean_speed_abs, 1.0e-12,
                                   "steady solution convergence should report mean-speed drift"));
    PetscCall(PicurvAssertRealNear(0.0, p_abs, 1.0e-12,
                                   "steady solution convergence pressure drift should be gauge invariant"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests periodic solution-convergence warmup and phase-aligned reference reuse.
 */
static PetscErrorCode TestSolutionConvergencePeriodicLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char log_path[PETSC_MAX_PATH_LEN];
    char header[4096];
    char row1[4096];
    char row2[4096];
    char row3[4096];
    PetscInt has_reference_1 = 0;
    PetscInt has_reference_2 = 0;
    PetscInt has_reference_3 = 0;
    PetscInt phase_step_1 = -1;
    PetscInt phase_step_2 = -1;
    PetscInt phase_step_3 = -1;
    PetscReal mean_speed_ref_2 = NAN;
    PetscReal mean_speed_abs_2 = NAN;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->solutionConvergenceMode = SOLUTION_CONVERGENCE_PERIODIC_DETERMINISTIC;
    simCtx->solutionConvergencePeriodSteps = 2;

    PetscCall(InitializeSolutionConvergenceState(simCtx));

    simCtx->step = 1;
    simCtx->ti = 0.1;
    PetscCall(SetUniformVelocityField(user, user->Ucat, 1.0, 0.0, 0.0));
    PetscCall(SetUniformScalarField(user, user->P, 2.0));
    PetscCall(LOG_SOLUTION_CONVERGENCE(simCtx));

    simCtx->step = 2;
    simCtx->ti = 0.2;
    PetscCall(SetUniformVelocityField(user, user->Ucat, 2.0, 0.0, 0.0));
    PetscCall(SetUniformScalarField(user, user->P, 5.0));
    PetscCall(LOG_SOLUTION_CONVERGENCE(simCtx));

    simCtx->step = 3;
    simCtx->ti = 0.3;
    PetscCall(SetUniformVelocityField(user, user->Ucat, 1.25, 0.0, 0.0));
    PetscCall(SetUniformScalarField(user, user->P, 9.0));
    PetscCall(LOG_SOLUTION_CONVERGENCE(simCtx));

    PetscCall(PetscSNPrintf(log_path, sizeof(log_path), "%s/solution_convergence.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(log_path, "periodic solution-convergence logging should write the log"));
    PetscCall(ReadLogHeaderAndRow(log_path, 1, header, sizeof(header), row1, sizeof(row1)));
    PetscCall(ReadLogHeaderAndRow(log_path, 2, header, sizeof(header), row2, sizeof(row2)));
    PetscCall(ReadLogHeaderAndRow(log_path, 3, header, sizeof(header), row3, sizeof(row3)));
    PetscCall(LogGetColumnInt(header, row1, "ref", &has_reference_1));
    PetscCall(LogGetColumnInt(header, row2, "ref", &has_reference_2));
    PetscCall(LogGetColumnInt(header, row3, "ref", &has_reference_3));
    PetscCall(LogGetColumnInt(header, row1, "ph", &phase_step_1));
    PetscCall(LogGetColumnInt(header, row2, "ph", &phase_step_2));
    PetscCall(LogGetColumnInt(header, row3, "ph", &phase_step_3));
    PetscCall(LogGetColumnReal(header, row3, "spd_ref", &mean_speed_ref_2));
    PetscCall(LogGetColumnReal(header, row3, "spd_abs", &mean_speed_abs_2));

    PetscCall(PicurvAssertIntEqual(0, has_reference_1,
                                   "the first periodic phase visit should log warmup without a reference"));
    PetscCall(PicurvAssertIntEqual(0, has_reference_2,
                                   "the first periodic cycle should fully warm up before comparisons begin"));
    PetscCall(PicurvAssertIntEqual(1, has_reference_3,
                                   "the repeated periodic phase visit should compare against the stored reference"));
    PetscCall(PicurvAssertIntEqual(1, phase_step_1,
                                   "periodic solution convergence should log the current phase slot"));
    PetscCall(PicurvAssertIntEqual(0, phase_step_2,
                                   "periodic solution convergence should log distinct phase slots during warmup"));
    PetscCall(PicurvAssertIntEqual(1, phase_step_3,
                                   "periodic solution convergence should reuse the same phase slot on later cycles"));
    PetscCall(PicurvAssertRealNear(1.0, mean_speed_ref_2, 1.0e-12,
                                   "periodic solution convergence should report the stored phase-aligned reference"));
    PetscCall(PicurvAssertRealNear(0.25, mean_speed_abs_2, 1.0e-12,
                                   "periodic solution convergence should report the phase-aligned drift"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests statistical solution-convergence sliding-window metrics.
 */
static PetscErrorCode TestSolutionConvergenceStatisticalLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char log_path[PETSC_MAX_PATH_LEN];
    char header[4096];
    char row[4096];
    PetscReal mean_speed_window = NAN;
    PetscReal mean_speed_window_prev = NAN;
    PetscReal mean_speed_window_abs = NAN;
    PetscReal mean_speed_rms_window = NAN;
    PetscReal mean_ke_window = NAN;
    PetscReal mean_ke_window_prev = NAN;
    PetscReal mean_ke_window_abs = NAN;
    PetscReal mean_ke_rms_window_abs = NAN;
    PetscInt has_reference = 0;
    const PetscReal speed_samples[4] = {1.0, 3.0, 2.0, 4.0};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->solutionConvergenceMode = SOLUTION_CONVERGENCE_STATISTICAL_STEADY;
    simCtx->solutionConvergenceWindowSteps = 2;

    PetscCall(InitializeSolutionConvergenceState(simCtx));
    for (PetscInt step = 0; step < 4; ++step) {
        simCtx->step = step + 1;
        simCtx->ti = 0.1 * (PetscReal)(step + 1);
        PetscCall(SetUniformVelocityField(user, user->Ucat, speed_samples[step], 0.0, 0.0));
        PetscCall(LOG_SOLUTION_CONVERGENCE(simCtx));
    }

    PetscCall(PetscSNPrintf(log_path, sizeof(log_path), "%s/solution_convergence.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(log_path, "statistical solution-convergence logging should write the log"));
    PetscCall(ReadLogHeaderAndRow(log_path, 4, header, sizeof(header), row, sizeof(row)));
    PetscCall(LogGetColumnInt(header, row, "ref", &has_reference));
    PetscCall(LogGetColumnReal(header, row, "spd_win", &mean_speed_window));
    PetscCall(LogGetColumnReal(header, row, "spd_win_prev", &mean_speed_window_prev));
    PetscCall(LogGetColumnReal(header, row, "spd_win_abs", &mean_speed_window_abs));
    PetscCall(LogGetColumnReal(header, row, "spd_rms_win", &mean_speed_rms_window));
    PetscCall(LogGetColumnReal(header, row, "ke_win", &mean_ke_window));
    PetscCall(LogGetColumnReal(header, row, "ke_win_prev", &mean_ke_window_prev));
    PetscCall(LogGetColumnReal(header, row, "ke_win_abs", &mean_ke_window_abs));
    PetscCall(LogGetColumnReal(header, row, "ke_rms_abs", &mean_ke_rms_window_abs));

    PetscCall(PicurvAssertIntEqual(1, has_reference,
                                   "statistical solution convergence should emit adjacent-window drift once two windows exist"));
    PetscCall(PicurvAssertRealNear(3.0, mean_speed_window, 1.0e-12,
                                   "statistical solution convergence should report the current mean-speed window"));
    PetscCall(PicurvAssertRealNear(2.0, mean_speed_window_prev, 1.0e-12,
                                   "statistical solution convergence should report the previous mean-speed window"));
    PetscCall(PicurvAssertRealNear(1.0, mean_speed_window_abs, 1.0e-12,
                                   "statistical solution convergence should report mean-speed window drift"));
    PetscCall(PicurvAssertRealNear(1.0, mean_speed_rms_window, 1.0e-12,
                                   "statistical solution convergence should report current mean-speed RMS"));
    PetscCall(PicurvAssertRealNear(5.0, mean_ke_window, 1.0e-12,
                                   "statistical solution convergence should report the current kinetic-energy window"));
    PetscCall(PicurvAssertRealNear(2.5, mean_ke_window_prev, 1.0e-12,
                                   "statistical solution convergence should report the previous kinetic-energy window"));
    PetscCall(PicurvAssertRealNear(2.5, mean_ke_window_abs, 1.0e-12,
                                   "statistical solution convergence should report kinetic-energy window drift"));
    PetscCall(PicurvAssertRealNear(1.0, mean_ke_rms_window_abs, 1.0e-12,
                                   "statistical solution convergence should report RMS kinetic-energy drift"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-logging PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"string-conversion-helpers", TestStringConversionHelpers},
        {"get-log-level-from-environment", TestGetLogLevelFromEnvironment},
        {"allowed-functions-filter", TestAllowedFunctionsFilter},
        {"particle-console-snapshot-cadence", TestParticleConsoleSnapshotCadence},
        {"logging-file-parsing-and-formatting-helpers", TestLoggingFileParsingAndFormattingHelpers},
        {"logging-continuity-and-field-diagnostics", TestLoggingContinuityAndFieldDiagnostics},
        {"interpolation-error-logging", TestInterpolationErrorLogging},
        {"scatter-metrics-logging", TestScatterMetricsLogging},
        {"particle-field-table-logging", TestParticleFieldTableLogging},
        {"particle-console-snapshot-logging", TestParticleConsoleSnapshotLogging},
        {"particle-metrics-logging", TestParticleMetricsLogging},
        {"search-metrics-logging", TestSearchMetricsLogging},
        {"field-anatomy-logging", TestFieldAnatomyLogging},
        {"profiling-lifecycle-helpers", TestProfilingLifecycleHelpers},
        {"solution-convergence-steady-logging", TestSolutionConvergenceSteadyLogging},
        {"solution-convergence-periodic-logging", TestSolutionConvergencePeriodicLogging},
        {"solution-convergence-statistical-logging", TestSolutionConvergenceStatisticalLogging},
    };

    (void)setenv("LOG_LEVEL", "INFO", 1);

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv logging tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-logging", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    set_allowed_functions(NULL, 0);

    ierr = PetscFinalize();
    return (int)ierr;
}
