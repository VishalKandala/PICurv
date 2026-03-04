#ifndef PICURV_TEST_SUPPORT_H
#define PICURV_TEST_SUPPORT_H

#include <stddef.h>

#include "variables.h"

typedef PetscErrorCode (*PicurvTestFn)(void);

typedef struct PicurvTestCase {
    const char   *name;
    PicurvTestFn  fn;
} PicurvTestCase;

PetscErrorCode PicurvRunTests(const char *suite_name, const PicurvTestCase *cases, size_t case_count);
PetscErrorCode PicurvCreateMinimalContexts(SimCtx **simCtx_out, UserCtx **user_out, PetscInt mx, PetscInt my, PetscInt mz);
PetscErrorCode PicurvDestroyMinimalContexts(SimCtx **simCtx, UserCtx **user);
PetscErrorCode PicurvPopulateIdentityMetrics(UserCtx *user);
PetscErrorCode PicurvCreateSwarmPair(UserCtx *user, PetscInt nlocal, const char *post_field_name);
PetscErrorCode PicurvMakeTempDir(char *path, size_t path_len);
PetscErrorCode PicurvEnsureDir(const char *path);
PetscErrorCode PicurvAssertRealNear(PetscReal expected, PetscReal actual, PetscReal tol, const char *context);
PetscErrorCode PicurvAssertIntEqual(PetscInt expected, PetscInt actual, const char *context);
PetscErrorCode PicurvAssertBool(PetscBool value, const char *context);
PetscErrorCode PicurvAssertFileExists(const char *path, const char *context);
PetscErrorCode PicurvAssertVecConstant(Vec vec, PetscScalar expected, PetscReal tol, const char *context);

#endif
