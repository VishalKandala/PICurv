/**
 * @file test_io.c
 * @brief C unit tests for I/O helpers, parsers, and startup-banner output.
 */

#include "test_support.h"

#include "io.h"

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
/**
 * @brief Tests cadence-based Eulerian output triggering.
 */

static PetscErrorCode TestShouldWriteDataOutput(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    simCtx.tiout = 5;

    PetscCall(PicurvAssertBool((PetscBool)!ShouldWriteDataOutput(NULL, 5), "NULL SimCtx should never request output"));
    PetscCall(PicurvAssertBool((PetscBool)!ShouldWriteDataOutput(&simCtx, 4), "non-cadence step should not trigger output"));
    PetscCall(PicurvAssertBool(ShouldWriteDataOutput(&simCtx, 10), "cadence-aligned step should trigger output"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests filesystem existence checks for files and directories.
 */

static PetscErrorCode TestVerifyPathExistence(void)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char filepath[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(filepath, sizeof(filepath), "%s/sample.txt", tmpdir));

    file = fopen(filepath, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create temp file '%s'.", filepath);
    fputs("picurv\n", file);
    fclose(file);

    PetscCall(VerifyPathExistence(tmpdir, PETSC_TRUE, PETSC_FALSE, "temp directory", &exists));
    PetscCall(PicurvAssertBool(exists, "VerifyPathExistence should find the temp directory"));

    PetscCall(VerifyPathExistence(filepath, PETSC_FALSE, PETSC_FALSE, "temp file", &exists));
    PetscCall(PicurvAssertBool(exists, "VerifyPathExistence should find the temp file"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests writing and reloading core Eulerian field vectors.
 */

static PetscErrorCode TestWriteAndReadSimulationFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char euler_dir[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(euler_dir, sizeof(euler_dir), "%s/%s", tmpdir, simCtx->euler_subdir));
    PetscCall(PicurvEnsureDir(euler_dir));

    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(VecSet(user->P, 4.5));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(VecSet(user->Ucat, 2.0));
    PetscCall(VecSet(user->Ucont, 3.0));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    PetscCall(WriteSimulationFields(user));
    PetscCall(VecZeroEntries(user->P));
    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(VecZeroEntries(user->Ucont));

    PetscCall(ReadSimulationFields(user, simCtx->step));
    PetscCall(PicurvAssertVecConstant(user->P, 4.5, 1.0e-12, "ReadSimulationFields should restore P"));
    PetscCall(PicurvAssertVecConstant(user->Ucat, 2.0, 1.0e-12, "ReadSimulationFields should restore Ucat"));
    PetscCall(PicurvAssertVecConstant(user->Ucont, 3.0, 1.0e-12, "ReadSimulationFields should restore Ucont"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests parsing of post-processing control settings from a file.
 */

static PetscErrorCode TestParsePostProcessingSettings(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char cfg_path[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscCalloc1(1, &simCtx->pps));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(cfg_path, sizeof(cfg_path), "%s/post.run", tmpdir));
    PetscCall(PetscStrncpy(simCtx->PostprocessingControlFile, cfg_path, sizeof(simCtx->PostprocessingControlFile)));

    file = fopen(cfg_path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create temp config file '%s'.", cfg_path);
    fputs("startTime = 2\n", file);
    fputs("endTime = 6\n", file);
    fputs("timeStep = 2\n", file);
    fputs("output_particles = true\n", file);
    fputs("output_prefix = SmokeField\n", file);
    fclose(file);

    PetscCall(ParsePostProcessingSettings(simCtx));
    PetscCall(PicurvAssertIntEqual(2, simCtx->pps->startTime, "ParsePostProcessingSettings should parse startTime"));
    PetscCall(PicurvAssertIntEqual(6, simCtx->pps->endTime, "ParsePostProcessingSettings should parse endTime"));
    PetscCall(PicurvAssertIntEqual(2, simCtx->pps->timeStep, "ParsePostProcessingSettings should parse timeStep"));
    PetscCall(PicurvAssertBool(simCtx->pps->outputParticles, "ParsePostProcessingSettings should parse output_particles"));
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(simCtx->pps->output_prefix, "SmokeField") == 0),
                               "ParsePostProcessingSettings should parse output_prefix"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests trimming of leading and trailing whitespace.
 */

static PetscErrorCode TestTrimWhitespace(void)
{
    char value_a[] = "   inlet_value   ";
    char value_b[] = "   ";

    PetscFunctionBeginUser;
    TrimWhitespace(value_a);
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(value_a, "inlet_value") == 0),
                               "TrimWhitespace should remove leading and trailing whitespace"));

    TrimWhitespace(value_b);
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(value_b, "") == 0),
                               "TrimWhitespace should reduce all-whitespace strings to empty"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests boundary-condition string parsers for face, type, and handler names.
 */

static PetscErrorCode TestBoundaryConditionStringParsers(void)
{
    BCFace face = BC_FACE_NEG_X;
    BCType type = WALL;
    BCHandlerType handler = BC_HANDLER_WALL_NOSLIP;

    PetscFunctionBeginUser;
    PetscCall(StringToBCFace("+Zeta", &face));
    PetscCall(PicurvAssertIntEqual(BC_FACE_POS_Z, face, "StringToBCFace should parse +Zeta"));

    PetscCall(StringToBCType("periodic", &type));
    PetscCall(PicurvAssertIntEqual(PERIODIC, type, "StringToBCType should parse PERIODIC case-insensitively"));

    PetscCall(StringToBCHandlerType("constant_flux", &handler));
    PetscCall(PicurvAssertIntEqual(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX, handler,
                                   "StringToBCHandlerType should parse constant_flux"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests validation of boundary-type and handler compatibility.
 */

static PetscErrorCode TestValidateBCHandlerForBCType(void)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(WALL, BC_HANDLER_WALL_NOSLIP) == 0),
                               "WALL + noslip should be a valid combination"));
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(PERIODIC, BC_HANDLER_PERIODIC_GEOMETRIC) == 0),
                               "PERIODIC + geometric should be a valid combination"));
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(INLET, BC_HANDLER_WALL_NOSLIP) != 0),
                               "INLET + noslip should be rejected"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests scaling-reference parsing and derived pressure scaling.
 */

static PetscErrorCode TestParseScalingInformation(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));

    PetscCall(PetscOptionsClearValue(NULL, "-scaling_L_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_U_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_rho_ref"));

    PetscCall(ParseScalingInformation(&simCtx));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.L_ref, 1.0e-12, "Default scaling_L_ref should be 1.0"));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.U_ref, 1.0e-12, "Default scaling_U_ref should be 1.0"));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.rho_ref, 1.0e-12, "Default scaling_rho_ref should be 1.0"));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.P_ref, 1.0e-12, "Default scaling_P_ref should be 1.0"));

    PetscCall(PetscOptionsSetValue(NULL, "-scaling_L_ref", "2.5"));
    PetscCall(PetscOptionsSetValue(NULL, "-scaling_U_ref", "4.0"));
    PetscCall(PetscOptionsSetValue(NULL, "-scaling_rho_ref", "1.2"));

    PetscCall(ParseScalingInformation(&simCtx));
    PetscCall(PicurvAssertRealNear(2.5, simCtx.scaling.L_ref, 1.0e-12, "scaling_L_ref should honor options"));
    PetscCall(PicurvAssertRealNear(4.0, simCtx.scaling.U_ref, 1.0e-12, "scaling_U_ref should honor options"));
    PetscCall(PicurvAssertRealNear(1.2, simCtx.scaling.rho_ref, 1.0e-12, "scaling_rho_ref should honor options"));
    PetscCall(PicurvAssertRealNear(19.2, simCtx.scaling.P_ref, 1.0e-12, "scaling_P_ref should be rho_ref*U_ref^2"));

    PetscCall(PetscOptionsClearValue(NULL, "-scaling_L_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_U_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_rho_ref"));
    PetscFunctionReturn(0);
}
/**
 * @brief Captures the startup banner into a temporary file-backed buffer.
 */
static PetscErrorCode CaptureBannerOutput(SimCtx *simCtx, char *captured, size_t captured_len)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char capture_path[PETSC_MAX_PATH_LEN];
    FILE *capture_file = NULL;
    int saved_stdout = -1;
    int capture_fd = -1;
    size_t bytes_read = 0;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");
    PetscCheck(captured != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Capture buffer cannot be NULL.");
    PetscCheck(captured_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Capture buffer must be non-empty.");

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(capture_path, sizeof(capture_path), "%s/banner.log", tmpdir));

    fflush(stdout);
    saved_stdout = dup(STDOUT_FILENO);
    PetscCheck(saved_stdout >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup(STDOUT_FILENO) failed.");
    capture_fd = open(capture_path, O_CREAT | O_TRUNC | O_WRONLY, 0600);
    PetscCheck(capture_fd >= 0, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open capture file '%s'.", capture_path);
    PetscCheck(dup2(capture_fd, STDOUT_FILENO) >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup2() failed while redirecting stdout.");
    close(capture_fd);
    capture_fd = -1;

    ierr = DisplayBanner(simCtx);
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
 * @brief Asserts that captured banner output contains one expected substring.
 */
static PetscErrorCode AssertCapturedContains(const char *captured, const char *needle, const char *message)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, needle) != NULL), message));
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that captured banner output omits one forbidden substring.
 */
static PetscErrorCode AssertCapturedOmits(const char *captured, const char *needle, const char *message)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, needle) == NULL), message));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests conditional startup-banner fields across particle and analytical cases.
 */

static PetscErrorCode TestDisplayBannerTracksConditionalStartupFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char captured[16384];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->OnlySetup = PETSC_FALSE;
    simCtx->StepsToRun = 5;
    simCtx->immersed = PETSC_FALSE;
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "solve", sizeof(simCtx->eulerianSource)));
    simCtx->particleConsoleOutputFreq = 7;
    simCtx->LoggingFrequency = 4;
    simCtx->ParticleInitialization = PARTICLE_INIT_POINT_SOURCE;

    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Run Mode                   : Full Simulation",
                                     "DisplayBanner should include the run mode"));
    PetscCall(AssertCapturedContains(captured, "Field/Restart Cadence      : every 2 step(s)",
                                     "DisplayBanner should include field/restart cadence"));
    PetscCall(AssertCapturedContains(captured, "Immersed Boundary          : DISABLED",
                                     "DisplayBanner should include immersed-boundary state"));
    PetscCall(AssertCapturedContains(captured, "Number of Particles         : 0",
                                     "DisplayBanner should include the total particle count"));
    PetscCall(AssertCapturedOmits(captured, "Particle Console Cadence",
                                  "DisplayBanner should omit particle console cadence when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particle Log Row Sampling",
                                  "DisplayBanner should omit particle row sampling when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particle Restart Mode",
                                  "DisplayBanner should omit particle restart mode when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particle Initialization Mode",
                                  "DisplayBanner should omit particle initialization mode when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Interpolation Method",
                                  "DisplayBanner should omit interpolation method when no particles are configured"));

    simCtx->StartStep = 3;
    simCtx->np = 8;
    simCtx->particleConsoleOutputFreq = 0;
    simCtx->LoggingFrequency = 4;
    PetscCall(PetscStrncpy(simCtx->particleRestartMode, "load", sizeof(simCtx->particleRestartMode)));
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Number of Particles         : 8",
                                     "DisplayBanner should include the active particle count"));
    PetscCall(AssertCapturedContains(captured, "Particle Console Cadence   : DISABLED",
                                     "DisplayBanner should show disabled particle console cadence when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Particle Log Row Sampling  : every 4 particle(s)",
                                     "DisplayBanner should include particle row sampling when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Particle Restart Mode      : load",
                                     "DisplayBanner should include particle restart mode for restarted particle runs"));
    PetscCall(AssertCapturedContains(captured, "Particle Initialization Mode: Point Source",
                                     "DisplayBanner should include particle initialization mode when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Interpolation Method       : Trilinear (direct cell-center)",
                                     "DisplayBanner should include default interpolation method when particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particles Initialized At",
                                  "DisplayBanner should omit inlet-face placement details for point-source particle initialization"));

    simCtx->StartStep = 0;
    simCtx->particleConsoleOutputFreq = 6;
    simCtx->ParticleInitialization = PARTICLE_INIT_SURFACE_RANDOM;
    user->inletFaceDefined = PETSC_TRUE;
    user->identifiedInletBCFace = (BCFace)0;
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "analytical", sizeof(simCtx->eulerianSource)));
    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "ZERO_FLOW", sizeof(simCtx->AnalyticalSolutionType)));
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Analytical Solution Type : ZERO_FLOW",
                                     "DisplayBanner should include the analytical solution type for analytical runs"));
    PetscCall(AssertCapturedContains(captured, "Particle Console Cadence   : every 6 step(s)",
                                     "DisplayBanner should include active particle console cadence when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Particle Initialization Mode: Surface: Random",
                                     "DisplayBanner should include particle initialization mode for analytical particle runs"));
    PetscCall(AssertCapturedContains(captured, "Particles Initialized At",
                                     "DisplayBanner should include inlet-face placement details for surface particle initialization"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-io PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"should-write-data-output", TestShouldWriteDataOutput},
        {"verify-path-existence", TestVerifyPathExistence},
        {"write-and-read-simulation-fields", TestWriteAndReadSimulationFields},
        {"parse-post-processing-settings", TestParsePostProcessingSettings},
        {"trim-whitespace", TestTrimWhitespace},
        {"bc-string-parsers", TestBoundaryConditionStringParsers},
        {"validate-bc-handler-for-type", TestValidateBCHandlerForBCType},
        {"parse-scaling-information", TestParseScalingInformation},
        {"display-banner-startup-summary", TestDisplayBannerTracksConditionalStartupFields},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv I/O tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-io", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
