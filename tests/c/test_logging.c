/**
 * @file test_logging.c
 * @brief C test module for PICurv.
 */

#include "test_support.h"

#include "logging.h"

#include <stdlib.h>
/**
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Test-local routine.
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
 * @brief Entry point for this unit-test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"get-log-level-from-environment", TestGetLogLevelFromEnvironment},
        {"allowed-functions-filter", TestAllowedFunctionsFilter},
        {"particle-console-snapshot-cadence", TestParticleConsoleSnapshotCadence},
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
