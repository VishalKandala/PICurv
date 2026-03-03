#include "test_support.h"

#include "io.h"

#include <stdio.h>

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
    PetscFunctionReturn(0);
}

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

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

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

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"should-write-data-output", TestShouldWriteDataOutput},
        {"verify-path-existence", TestVerifyPathExistence},
        {"write-and-read-simulation-fields", TestWriteAndReadSimulationFields},
        {"parse-post-processing-settings", TestParsePostProcessingSettings},
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
