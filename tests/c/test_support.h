/**
 * @file test_support.h
 * @brief Shared declarations for the PICurv C test fixture and assertion layer.
 */

#ifndef PICURV_TEST_SUPPORT_H
#define PICURV_TEST_SUPPORT_H

#include <stddef.h>

#include "variables.h"

/** @brief Signature for a single C unit-test entry point. */
typedef PetscErrorCode (*PicurvTestFn)(void);

/** @brief Named test case descriptor consumed by `PicurvRunTests`. */
typedef struct PicurvTestCase {
    const char   *name;
    PicurvTestFn  fn;
} PicurvTestCase;

/** @brief Executes a table of unit tests and reports aggregated pass/fail status. */
PetscErrorCode PicurvRunTests(const char *suite_name, const PicurvTestCase *cases, size_t case_count);

/** @brief Creates minimal PETSc/solver contexts used by isolated kernel tests. */
PetscErrorCode PicurvCreateMinimalContexts(SimCtx **simCtx_out, UserCtx **user_out, PetscInt mx, PetscInt my, PetscInt mz);

/** @brief Destroys contexts previously created by `PicurvCreateMinimalContexts`. */
PetscErrorCode PicurvDestroyMinimalContexts(SimCtx **simCtx, UserCtx **user);

/** @brief Builds a tiny runtime context through the real setup path for behavior-level tests. */
PetscErrorCode PicurvBuildTinyRuntimeContext(const char *bcs_contents,
                                             PetscBool enable_particles,
                                             SimCtx **simCtx_out,
                                             UserCtx **user_out,
                                             char *tmpdir,
                                             size_t tmpdir_len);

/** @brief Finalizes and frees a runtime context built by `PicurvBuildTinyRuntimeContext`. */
PetscErrorCode PicurvDestroyRuntimeContext(SimCtx **simCtx_ptr);

/** @brief Populates cell center coordinates for a uniform grid on [0,1]^3. */
PetscErrorCode PicurvPopulateUniformCellCenters(UserCtx *user);

/** @brief Fills metric vectors with identity metrics for Cartesian-reference tests. */
PetscErrorCode PicurvPopulateIdentityMetrics(UserCtx *user);

/** @brief Creates paired runtime/post swarms with optional extra post-processing field registration. */
PetscErrorCode PicurvCreateSwarmPair(UserCtx *user, PetscInt nlocal, const char *post_field_name);

/** @brief Creates a unique temporary directory path and materializes the directory. */
PetscErrorCode PicurvMakeTempDir(char *path, size_t path_len);

/** @brief Recursively removes a temporary directory created by PicurvMakeTempDir. */
PetscErrorCode PicurvRemoveTempDir(const char *path);

/** @brief Ensures a directory exists, creating it if required. */
PetscErrorCode PicurvEnsureDir(const char *path);

/** @brief Asserts two real values are within tolerance for test validation. */
PetscErrorCode PicurvAssertRealNear(PetscReal expected, PetscReal actual, PetscReal tol, const char *context);

/** @brief Asserts two integer values are exactly equal for test validation. */
PetscErrorCode PicurvAssertIntEqual(PetscInt expected, PetscInt actual, const char *context);

/** @brief Asserts a boolean condition is true for test validation. */
PetscErrorCode PicurvAssertBool(PetscBool value, const char *context);

/** @brief Asserts a filesystem path exists. */
PetscErrorCode PicurvAssertFileExists(const char *path, const char *context);

/** @brief Asserts every entry of a vector equals a constant within tolerance. */
PetscErrorCode PicurvAssertVecConstant(Vec vec, PetscScalar expected, PetscReal tol, const char *context);

#endif
