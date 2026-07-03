/**
 * @file test_momentum_newton_krylov.c
 * @brief Focused tests for the version-one matrix-free momentum solver.
 *
 * The implementation is included intentionally: its callback helpers remain
 * private in production while this translation unit can verify them directly.
 */

#include "test_support.h"
#include "initialcondition.h"

/* Rename only the included public entry point. Private callback tests use this
 * copy, while lifecycle tests below call the separately linked production object. */
#define MomentumSolver_NewtonKrylov MomentumSolver_NewtonKrylov_PrivateCopy
#include "../../src/momentum_newton_krylov.c"
#undef MomentumSolver_NewtonKrylov

static const char *geometric_periodic_bcs =
    "-Xi PERIODIC geometric\n"
    "+Xi PERIODIC geometric\n"
    "-Eta WALL noslip\n"
    "+Eta WALL noslip\n"
    "-Zeta INLET constant_velocity vx=0.0 vy=0.0 vz=1.5\n"
    "+Zeta OUTLET conservation\n";

static const char *fixed_wall_bcs =
    "-Xi WALL noslip\n"
    "+Xi WALL noslip\n"
    "-Eta WALL noslip\n"
    "+Eta WALL noslip\n"
    "-Zeta WALL noslip\n"
    "+Zeta WALL noslip\n";

static const char *parabolic_bcs =
    "-Xi WALL noslip\n"
    "+Xi WALL noslip\n"
    "-Eta WALL noslip\n"
    "+Eta WALL noslip\n"
    "-Zeta INLET parabolic v_max=1.5\n"
    "+Zeta OUTLET conservation\n";

static const char *periodic_x_bcs =
    "-Xi PERIODIC geometric\n+Xi PERIODIC geometric\n"
    "-Eta WALL noslip\n+Eta WALL noslip\n-Zeta WALL noslip\n+Zeta WALL noslip\n";

static const char *periodic_y_bcs =
    "-Xi WALL noslip\n+Xi WALL noslip\n"
    "-Eta PERIODIC geometric\n+Eta PERIODIC geometric\n-Zeta WALL noslip\n+Zeta WALL noslip\n";

static const char *periodic_z_bcs =
    "-Xi WALL noslip\n+Xi WALL noslip\n-Eta WALL noslip\n+Eta WALL noslip\n"
    "-Zeta PERIODIC geometric\n+Zeta PERIODIC geometric\n";

static const char *periodic_xy_bcs =
    "-Xi PERIODIC geometric\n+Xi PERIODIC geometric\n"
    "-Eta PERIODIC geometric\n+Eta PERIODIC geometric\n"
    "-Zeta WALL noslip\n+Zeta WALL noslip\n";

static const char *periodic_xyz_bcs =
    "-Xi PERIODIC geometric\n+Xi PERIODIC geometric\n"
    "-Eta PERIODIC geometric\n+Eta PERIODIC geometric\n"
    "-Zeta PERIODIC geometric\n+Zeta PERIODIC geometric\n";

/** @brief Checks a structured log's row count and required text after a collective solve. */
static PetscErrorCode AssertNewtonLog(const char *path, PetscInt expected_rows,
                                      const char *needle_a, const char *needle_b)
{
    FILE *file = NULL;
    char line[4096];
    PetscInt rows = 0;
    PetscBool found_a = needle_a ? PETSC_FALSE : PETSC_TRUE;
    PetscBool found_b = needle_b ? PETSC_FALSE : PETSC_TRUE;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    file = fopen(path, "r");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
               "Expected Newton log does not exist: %s", path);
    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "step:", 5) == 0) rows++;
        if (needle_a && strstr(line, needle_a)) found_a = PETSC_TRUE;
        if (needle_b && strstr(line, needle_b)) found_b = PETSC_TRUE;
    }
    fclose(file);
    if (expected_rows >= 0) {
        PetscCall(PicurvAssertIntEqual(expected_rows, rows,
                                       "Newton log must contain one nonduplicated row per solve"));
    } else {
        PetscCall(PicurvAssertBool((PetscBool)(rows >= -expected_rows),
                                   "enabled Newton history must contain the expected iteration rows"));
    }
    PetscCall(PicurvAssertBool(found_a, "Newton log is missing required structured content"));
    PetscCall(PicurvAssertBool(found_b, "Newton log is missing required structured content"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Builds and initializes a small runtime context for Newton tests.
 * @param bcs Optional boundary configuration text.
 * @param simCtx Returned simulation context.
 * @param user Returned finest-level block context.
 * @param tmpdir Returned temporary directory.
 * @param tmpdir_len Capacity of tmpdir.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode BuildNewtonFixture(const char *bcs, SimCtx **simCtx, UserCtx **user,
                                         char *tmpdir, size_t tmpdir_len)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContext(bcs, PETSC_FALSE, simCtx, user, tmpdir, tmpdir_len));
    PetscCall(InitializeEulerianState(*simCtx));
    (*simCtx)->mom_solver_type = MOMENTUM_SOLVER_NEWTON_KRYLOV;
    PetscCall(PicurvAssertBool((PetscBool)((*user)->Rhs == NULL),
                               "Newton fixture must enter with no persistent Rhs workspace"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Destroys a Newton test fixture and its temporary files.
 * @param simCtx Fixture simulation context.
 * @param tmpdir Fixture temporary directory.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode DestroyNewtonFixture(SimCtx **simCtx, char *tmpdir)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvDestroyRuntimeContext(simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Writes one static 5x5 PICSLICE profile used by the full runtime fixture.
 * @param path Output profile path.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode WriteNewtonPicSlice(const char *path)
{
    FILE *fd = NULL;

    PetscFunctionBeginUser;
    fd = fopen(path, "w");
    PetscCheck(fd != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
               "Could not create Newton test PICSLICE %s.", path);
    PetscCheck(fprintf(fd, "PICSLICE\n1\n5 5\n") >= 0,
               PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE, "Could not write PICSLICE header.");
    for (PetscInt row = 0; row < 25; ++row) {
        PetscCheck(fprintf(fd, "%.16e\n", 1.0 + 0.01 * (double)row) >= 0,
                   PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE, "Could not write PICSLICE value.");
    }
    fclose(fd);
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Checks callback repeatability, diagnostic-state independence, and X integrity.
 * @param bcs Boundary configuration text, or NULL for the standard inlet/outlet fixture.
 * @param label Configuration label used in assertion diagnostics.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode CheckResidualRepeatabilityForBC(const char *bcs, const char *label)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL, x_copy = NULL, f1 = NULL, f2 = NULL, delta = NULL;
    PetscReal norm = 0.0;
    MomentumNewtonKrylovContext ctx;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &x_copy));
    PetscCall(VecDuplicate(user->Ucont, &f1));
    PetscCall(VecDuplicate(user->Ucont, &f2));
    PetscCall(VecDuplicate(user->Ucont, &delta));
    PetscCall(VecCopy(user->Ucont, x));
    PetscCall(VecShift(x, 0.125));
    PetscCall(VecCopy(x, x_copy));
    ctx.user = user;

    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f1, &ctx));
    /* Poison every piece of hidden state a non-deterministic residual could lean
     * on: the boundary flux/area diagnostics AND the persistent Cartesian fields
     * (Ucat/lUcat). The conservation-outlet handler reads lUcat during its first
     * boundary sweep, so a callback that does not reconstruct the Cartesian state
     * from X before applying boundary conditions would produce a different F here.
     * This assertion therefore fails if the deterministic pre-boundary seed in
     * MomentumNewtonKrylov_FormResidual() is ever removed. */
    simCtx->FluxInSum = 1234.0;
    simCtx->FluxOutSum = -4321.0;
    simCtx->FarFluxInSum = 77.0;
    simCtx->FarFluxOutSum = -88.0;
    PetscCall(VecSet(user->Ucat, 7.0));
    PetscCall(VecSet(user->lUcat, 7.0));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f2, &ctx));
    PetscCall(VecWAXPY(delta, -1.0, f1, f2));
    PetscCall(VecNorm(delta, NORM_INFINITY, &norm));
    PetscCheck(norm <= 1.0e-13, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "%s residual changed between identical evaluations (inf norm=%g).", label, (double)norm);
    PetscCall(VecNorm(delta, NORM_2, &norm));
    PetscCheck(norm <= 1.0e-12, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "%s residual changed between identical evaluations (L2 norm=%g).", label, (double)norm);
    PetscCall(VecWAXPY(delta, -1.0, x_copy, x));
    PetscCall(VecNorm(delta, NORM_INFINITY, &norm));
    PetscCheck(norm == 0.0, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "%s residual callback modified X (norm=%g).", label, (double)norm);

    PetscCall(VecDestroy(&delta));
    PetscCall(VecDestroy(&f2));
    PetscCall(VecDestroy(&f1));
    PetscCall(VecDestroy(&x_copy));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Returns whether this rank owns one global DMDA grid point. */
static PetscBool OwnsStoredPoint(UserCtx *user, PetscInt i, PetscInt j, PetscInt k)
{
    return (PetscBool)(i >= user->info.xs && i < user->info.xs + user->info.xm &&
                       j >= user->info.ys && j < user->info.ys + user->info.ym &&
                       k >= user->info.zs && k < user->info.zs + user->info.zm);
}

/**
 * @brief Adds a scalar perturbation to one stored staggered component.
 * @param user Block context defining vector ownership.
 * @param vec Vector to modify.
 * @param i Global i index.
 * @param j Global j index.
 * @param k Global k index.
 * @param component Component 0=x, 1=y, 2=z.
 * @param delta Increment to apply.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode PerturbStoredValue(UserCtx *user, Vec vec, PetscInt i, PetscInt j,
                                         PetscInt k, PetscInt component, PetscScalar delta)
{
    Cmpnts ***a = NULL;

    PetscFunctionBeginUser;
    if (OwnsStoredPoint(user, i, j, k)) {
        PetscCall(DMDAVecGetArray(user->fda, vec, &a));
        if (component == 0) a[k][j][i].x += delta;
        else if (component == 1) a[k][j][i].y += delta;
        else a[k][j][i].z += delta;
        PetscCall(DMDAVecRestoreArray(user->fda, vec, &a));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Reads one globally indexed stored component on any MPI decomposition.
 * @param user Block context defining vector ownership.
 * @param vec Vector to inspect.
 * @param i Global i index.
 * @param j Global j index.
 * @param k Global k index.
 * @param component Component 0=x, 1=y, 2=z.
 * @param value Returned globally reduced scalar.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode GetStoredValue(UserCtx *user, Vec vec, PetscInt i, PetscInt j,
                                     PetscInt k, PetscInt component, PetscScalar *value)
{
    Cmpnts ***a = NULL;
    PetscScalar local = 0.0;

    PetscFunctionBeginUser;
    if (OwnsStoredPoint(user, i, j, k)) {
        PetscCall(DMDAVecGetArrayRead(user->fda, vec, &a));
        local = component == 0 ? a[k][j][i].x : (component == 1 ? a[k][j][i].y : a[k][j][i].z);
        PetscCall(DMDAVecRestoreArrayRead(user->fda, vec, &a));
    }
    PetscCallMPI(MPI_Allreduce(&local, value, 1, MPIU_SCALAR, MPIU_SUM, PETSC_COMM_WORLD));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Finite-differences one callback row with respect to one stored unknown.
 * @param user Active block context with allocated Rhs.
 * @param x Base trial vector.
 * @param row_i Residual-row i index.
 * @param row_j Residual-row j index.
 * @param row_k Residual-row k index.
 * @param row_component Residual-row component.
 * @param col_i Perturbed unknown i index.
 * @param col_j Perturbed unknown j index.
 * @param col_k Perturbed unknown k index.
 * @param col_component Perturbed unknown component.
 * @param derivative Returned finite-difference derivative of the row w.r.t. the unknown.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MeasureStoredDerivative(UserCtx *user, Vec x,
                                              PetscInt row_i, PetscInt row_j, PetscInt row_k,
                                              PetscInt row_component, PetscInt col_i, PetscInt col_j,
                                              PetscInt col_k, PetscInt col_component,
                                              PetscReal *derivative)
{
    const PetscReal epsilon = 1.0e-6;
    Vec f0 = NULL, fp = NULL, xp = NULL;
    PetscScalar base_value = 0.0, perturbed_value = 0.0;
    MomentumNewtonKrylovContext ctx = {user};

    PetscFunctionBeginUser;
    PetscCall(VecDuplicate(x, &f0));
    PetscCall(VecDuplicate(x, &fp));
    PetscCall(VecDuplicate(x, &xp));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f0, &ctx));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f0, &ctx));
    PetscCall(VecCopy(x, xp));
    PetscCall(PerturbStoredValue(user, xp, col_i, col_j, col_k, col_component, epsilon));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, xp, fp, &ctx));
    PetscCall(GetStoredValue(user, f0, row_i, row_j, row_k, row_component, &base_value));
    PetscCall(GetStoredValue(user, fp, row_i, row_j, row_k, row_component, &perturbed_value));
    *derivative = PetscRealPart((perturbed_value - base_value) / epsilon);
    PetscCall(VecDestroy(&xp));
    PetscCall(VecDestroy(&fp));
    PetscCall(VecDestroy(&f0));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Asserts one finite-differenced callback row derivative equals an expected value.
 * @param user Active block context with allocated Rhs.
 * @param x Base trial vector.
 * @param row_i Residual-row i index.
 * @param row_j Residual-row j index.
 * @param row_k Residual-row k index.
 * @param row_component Residual-row component.
 * @param col_i Perturbed unknown i index.
 * @param col_j Perturbed unknown j index.
 * @param col_k Perturbed unknown k index.
 * @param col_component Perturbed unknown component.
 * @param expected Expected derivative.
 * @param tolerance Absolute derivative tolerance.
 * @param label Assertion label.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode CheckStoredDerivative(UserCtx *user, Vec x,
                                            PetscInt row_i, PetscInt row_j, PetscInt row_k,
                                            PetscInt row_component, PetscInt col_i, PetscInt col_j,
                                            PetscInt col_k, PetscInt col_component,
                                            PetscReal expected, PetscReal tolerance, const char *label)
{
    PetscReal derivative = 0.0;

    PetscFunctionBeginUser;
    PetscCall(MeasureStoredDerivative(user, x, row_i, row_j, row_k, row_component,
                                      col_i, col_j, col_k, col_component, &derivative));
    PetscCall(PicurvAssertRealNear(expected, derivative, tolerance, label));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies repeatable callback output and read-only trial input. */
static PetscErrorCode TestResidualRepeatabilityAndInputIntegrity(void)
{
    char profile_dir[PETSC_MAX_PATH_LEN] = "";
    char profile_path[PETSC_MAX_PATH_LEN];
    char file_bcs[2 * PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(CheckResidualRepeatabilityForBC(fixed_wall_bcs, "fixed walls"));
    PetscCall(CheckResidualRepeatabilityForBC(NULL, "constant inlet/conservation outlet"));
    PetscCall(CheckResidualRepeatabilityForBC(parabolic_bcs, "parabolic inlet/conservation outlet"));
    PetscCall(CheckResidualRepeatabilityForBC(periodic_x_bcs, "x periodic"));
    PetscCall(CheckResidualRepeatabilityForBC(periodic_y_bcs, "y periodic"));
    PetscCall(CheckResidualRepeatabilityForBC(periodic_z_bcs, "z periodic"));
    PetscCall(CheckResidualRepeatabilityForBC(periodic_xy_bcs, "mixed x-y periodic"));

    PetscCall(PicurvMakeTempDir(profile_dir, sizeof(profile_dir)));
    PetscCall(PetscSNPrintf(profile_path, sizeof(profile_path), "%s/inlet.picslice", profile_dir));
    PetscCall(WriteNewtonPicSlice(profile_path));
    PetscCall(PetscSNPrintf(file_bcs, sizeof(file_bcs),
        "-Xi WALL noslip\n+Xi WALL noslip\n-Eta WALL noslip\n+Eta WALL noslip\n"
        "-Zeta INLET prescribed_flow source_file=%s\n+Zeta OUTLET conservation\n", profile_path));
    PetscCall(CheckResidualRepeatabilityForBC(file_bcs, "file inlet/conservation outlet"));
    PetscCall(PicurvRemoveTempDir(profile_dir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies fixed, periodic-duplicate, and interior residual rows. */
static PetscErrorCode TestConstraintRows(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL, f = NULL;
    Cmpnts ***xa = NULL, ***fa = NULL, ***conditioned = NULL, ***rhs = NULL;
    MomentumNewtonKrylovContext ctx;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(geometric_periodic_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &f));
    PetscCall(VecSet(x, 0.25));
    PetscCall(DMDAVecGetArray(user->fda, x, &xa));
    if (user->info.xs == 0 && 2 >= user->info.ys && 2 < user->info.ys + user->info.ym &&
        2 >= user->info.zs && 2 < user->info.zs + user->info.zm) xa[2][2][0].x = 3.0;
    if (user->info.xs + user->info.xm == user->info.mx &&
        2 >= user->info.ys && 2 < user->info.ys + user->info.ym &&
        2 >= user->info.zs && 2 < user->info.zs + user->info.zm) xa[2][2][user->info.mx - 2].x = 1.25;
    PetscCall(DMDAVecRestoreArray(user->fda, x, &xa));
    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f, &ctx));

    PetscCall(DMDAVecGetArrayRead(user->fda, x, &xa));
    PetscCall(DMDAVecGetArrayRead(user->fda, f, &fa));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &conditioned));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Rhs, &rhs));
    if (user->info.xs == 0 && 2 >= user->info.ys && 2 < user->info.ys + user->info.ym &&
        2 >= user->info.zs && 2 < user->info.zs + user->info.zm) {
        PetscCall(PicurvAssertRealNear(1.75, fa[2][2][0].x, 1.0e-12,
                                       "periodic duplicate row must be Xdup-Xrep"));
    }
    if (user->info.ys == 0 && 2 >= user->info.xs && 2 < user->info.xs + user->info.xm &&
        2 >= user->info.zs && 2 < user->info.zs + user->info.zm) {
        PetscCall(PicurvAssertRealNear(xa[2][0][2].y - conditioned[2][0][2].y,
                                       fa[2][0][2].y, 1.0e-12,
                                       "fixed wall row must be X minus conditioned boundary value"));
    }
    if (2 >= user->info.xs && 2 < user->info.xs + user->info.xm &&
        2 >= user->info.ys && 2 < user->info.ys + user->info.ym &&
        2 >= user->info.zs && 2 < user->info.zs + user->info.zm) {
        PetscCall(PicurvAssertRealNear(-rhs[2][2][2].z, fa[2][2][2].z, 1.0e-12,
                                       "unconstrained interior row must retain the physical residual"));
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Rhs, &rhs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &conditioned));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, f, &fa));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, x, &xa));

    PetscCall(VecDestroy(&f));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Proves unit derivatives for every nonperiodic stored-row category and face. */
static PetscErrorCode TestFixedConstraintDerivativesAllFaces(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL;
    const PetscInt size[3] = {7, 7, 7};

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecSet(x, 0.2));

    for (PetscInt axis = 0; axis < 3; ++axis) {
        PetscInt coord[3] = {2, 2, 2}, ri, rj, rk;
        PetscInt tangent = (axis + 1) % 3;
        MomentumNewtonKrylovRowType row;

        coord[axis] = 0;
        row = MomentumNewtonKrylov_ClassifyRow(user, coord[0], coord[1], coord[2], axis, &ri, &rj, &rk);
        PetscCall(PicurvAssertIntEqual(MOM_NK_ROW_FIXED_CONDITIONED, row,
                                       "negative face normal row classification"));
        PetscCall(CheckStoredDerivative(user, x, coord[0], coord[1], coord[2], axis,
            coord[0], coord[1], coord[2], axis, 1.0, 1.0e-8,
            "negative face normal fixed derivative"));

        row = MomentumNewtonKrylov_ClassifyRow(user, coord[0], coord[1], coord[2], tangent, &ri, &rj, &rk);
        PetscCall(PicurvAssertIntEqual(MOM_NK_ROW_FIXED_HOMOGENEOUS, row,
                                       "negative face tangential row classification"));
        PetscCall(CheckStoredDerivative(user, x, coord[0], coord[1], coord[2], tangent,
            coord[0], coord[1], coord[2], tangent, 1.0, 1.0e-8,
            "negative face tangential homogeneous derivative"));

        coord[axis] = size[axis] - 2;
        row = MomentumNewtonKrylov_ClassifyRow(user, coord[0], coord[1], coord[2], axis, &ri, &rj, &rk);
        PetscCall(PicurvAssertIntEqual(MOM_NK_ROW_FIXED_CONDITIONED, row,
                                       "positive physical normal row classification"));
        PetscCall(CheckStoredDerivative(user, x, coord[0], coord[1], coord[2], axis,
            coord[0], coord[1], coord[2], axis, 1.0, 1.0e-8,
            "positive physical normal fixed derivative"));
        row = MomentumNewtonKrylov_ClassifyRow(user, coord[0], coord[1], coord[2], tangent, &ri, &rj, &rk);
        PetscCall(PicurvAssertIntEqual(MOM_NK_ROW_PHYSICAL, row,
                                       "positive physical tangential row classification"));

        coord[axis] = size[axis] - 1;
        for (PetscInt component = 0; component < 3; ++component) {
            row = MomentumNewtonKrylov_ClassifyRow(user, coord[0], coord[1], coord[2], component, &ri, &rj, &rk);
            PetscCall(PicurvAssertIntEqual(MOM_NK_ROW_FIXED_HOMOGENEOUS, row,
                                           "positive dummy row classification"));
            PetscCall(CheckStoredDerivative(user, x, coord[0], coord[1], coord[2], component,
                coord[0], coord[1], coord[2], component, 1.0, 1.0e-8,
                "positive dummy homogeneous derivative"));
        }
    }

    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Proves admitted inlet and outlet face-normal rows have unit self derivatives. */
static PetscErrorCode TestInletOutletConstraintDerivatives(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(NULL, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecCopy(user->Ucont, x));
    PetscCall(VecShift(x, 0.1));
    PetscCall(CheckStoredDerivative(user, x, 2, 2, 0, 2, 2, 2, 0, 2,
                                     1.0, 1.0e-8, "constant inlet fixed derivative"));
    /* The constant-velocity inlet imposes a value independent of X, so its
     * conditioned row F = X - cv has an exact unit self derivative.
     *
     * The conservation outlet is different: cv is the corrected outlet flux,
     * which the deterministic residual now reconstructs from the current X
     * (Ucat is seeded from X before the first outlet pass). Perturbing the
     * outlet-normal DOF therefore changes cv, so the self derivative is
     * 1 - dcv/dX and is strictly less than one. A self derivative of exactly
     * 1.0 here was an artifact of the pre-fix residual reading a stale
     * Cartesian state, i.e. an outlet correction decoupled from X. Assert the
     * derivative is (a) deterministic across independent evaluations -- the
     * residual-purity property -- and (b) reflects real conservation coupling
     * (0 < d < 1), rather than asserting a fixture-specific magic number. */
    {
        PetscReal d0 = 0.0, d1 = 0.0;
        PetscCall(MeasureStoredDerivative(user, x, 2, 2, user->info.mz - 2, 2,
                                          2, 2, user->info.mz - 2, 2, &d0));
        PetscCall(MeasureStoredDerivative(user, x, 2, 2, user->info.mz - 2, 2,
                                          2, 2, user->info.mz - 2, 2, &d1));
        PetscCall(PicurvAssertRealNear(d0, d1, 1.0e-9,
            "conservation outlet self derivative must be deterministic"));
        PetscCheck(d0 > 1.0e-3 && d0 < 1.0 - 1.0e-3, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
            "conservation outlet self derivative must reflect X-coupling (0<d<1), got %g.",
            (double)d0);
    }
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Checks one periodic configuration's endpoint derivatives on every component.
 * @param bcs Boundary text selecting the periodic axis.
 * @param axis Periodic axis index.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode CheckSingleAxisPeriodicDerivatives(const char *bcs, PetscInt axis)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL;
    PetscInt size[3], dup[3] = {2, 2, 2}, rep[3] = {2, 2, 2}, unrelated[3] = {3, 3, 3};

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    size[0] = user->info.mx; size[1] = user->info.my; size[2] = user->info.mz;
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecSet(x, 0.15));
    for (PetscInt side = 0; side < 2; ++side) {
        dup[axis] = side == 0 ? 0 : size[axis] - 1;
        rep[axis] = side == 0 ? size[axis] - 2 : 1;
        for (PetscInt component = 0; component < 3; ++component) {
            PetscCall(CheckStoredDerivative(user, x, dup[0], dup[1], dup[2], component,
                dup[0], dup[1], dup[2], component, 1.0, 1.0e-8,
                "periodic duplicate self derivative"));
            PetscCall(CheckStoredDerivative(user, x, dup[0], dup[1], dup[2], component,
                rep[0], rep[1], rep[2], component, -1.0, 1.0e-8,
                "periodic representative derivative"));
            PetscCall(CheckStoredDerivative(user, x, dup[0], dup[1], dup[2], component,
                unrelated[0], unrelated[1], unrelated[2], (component + 1) % 3, 0.0, 1.0e-8,
                "periodic constraint unrelated derivative"));
        }
    }
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Proves single-, double-, triple-, and mixed-boundary periodic equations. */
static PetscErrorCode TestPeriodicConstraintDerivativesAndIntersections(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL, f = NULL;
    MomentumNewtonKrylovContext ctx;
    PetscScalar xdup, synced, residual;

    PetscFunctionBeginUser;
    PetscCall(CheckSingleAxisPeriodicDerivatives(periodic_x_bcs, 0));
    PetscCall(CheckSingleAxisPeriodicDerivatives(periodic_y_bcs, 1));
    PetscCall(CheckSingleAxisPeriodicDerivatives(periodic_z_bcs, 2));

    PetscCall(BuildNewtonFixture(periodic_xy_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &f));
    PetscCall(VecSet(x, 0.0));
    PetscCall(PerturbStoredValue(user, x, 0, 0, 2, 0, 3.0));
    PetscCall(PerturbStoredValue(user, x, user->info.mx - 2, user->info.my - 2, 2, 0, 1.25));
    PetscCall(VecCopy(x, user->Ucont));
    { const char *fields[] = {"Ucont"}; PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields)); }
    PetscCall(GetStoredValue(user, x, 0, 0, 2, 0, &xdup));
    PetscCall(GetStoredValue(user, user->Ucont, 0, 0, 2, 0, &synced));
    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f, &ctx));
    PetscCall(GetStoredValue(user, f, 0, 0, 2, 0, &residual));
    PetscCall(PicurvAssertRealNear(PetscRealPart(xdup - synced), PetscRealPart(residual), 1.0e-12,
                                   "doubly periodic edge must use production synchronized representative"));
    PetscCall(CheckStoredDerivative(user, x, 0, 0, 2, 0, 0, 0, 2, 0, 1.0, 1.0e-8,
                                     "doubly periodic edge self derivative"));
    PetscCall(CheckStoredDerivative(user, x, 0, 0, 2, 0,
        user->info.mx - 2, user->info.my - 2, 2, 0, -1.0, 1.0e-8,
        "doubly periodic edge representative derivative"));
    PetscCall(VecDestroy(&f)); PetscCall(VecDestroy(&x)); PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));

    tmpdir[0] = '\0';
    PetscCall(BuildNewtonFixture(periodic_xyz_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x)); PetscCall(VecSet(x, 0.1));
    PetscCall(CheckStoredDerivative(user, x, 0, 0, 0, 2, 0, 0, 0, 2, 1.0, 1.0e-8,
                                     "fully periodic corner self derivative"));
    PetscCall(CheckStoredDerivative(user, x, 0, 0, 0, 2,
        user->info.mx - 2, user->info.my - 2, user->info.mz - 2, 2,
        -1.0, 1.0e-8, "fully periodic corner representative derivative"));
    PetscCall(VecDestroy(&x)); PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));

    tmpdir[0] = '\0';
    PetscCall(BuildNewtonFixture(periodic_x_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs)); PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecSet(x, 0.1));
    PetscCall(CheckStoredDerivative(user, x, 0, 0, 2, 1, 0, 0, 2, 1, 1.0, 1.0e-8,
                                     "periodic-wall intersection self derivative"));
    PetscCall(CheckStoredDerivative(user, x, 0, 0, 2, 1,
        user->info.mx - 2, 0, 2, 1, -1.0, 1.0e-8,
        "periodic-wall intersection representative derivative"));
    PetscCall(VecDestroy(&x)); PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Compares PETSc's matrix-free action with direct differencing. */
static PetscErrorCode TestMatrixFreeDerivative(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    SNES snes = NULL;
    Mat J = NULL;
    Vec x = NULL, xp = NULL, f0 = NULL, fp = NULL, v = NULL, jv = NULL, fd = NULL;
    PetscReal h = 0.0, error = 0.0, scale = 0.0;
    MomentumNewtonKrylovContext ctx;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(NULL, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &xp));
    PetscCall(VecDuplicate(user->Ucont, &f0));
    PetscCall(VecDuplicate(user->Ucont, &fp));
    PetscCall(VecDuplicate(user->Ucont, &v));
    PetscCall(VecDuplicate(user->Ucont, &jv));
    PetscCall(VecDuplicate(user->Ucont, &fd));
    PetscCall(VecCopy(user->Ucont, x));
    PetscCall(VecShift(x, 0.05));
    PetscCall(VecSet(v, 0.5));
    ctx.user = user;

    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetDM(snes, user->fda));
    PetscCall(SNESSetFunction(snes, f0, MomentumNewtonKrylov_FormResidual, &ctx));
    PetscCall(MatCreateSNESMF(snes, &J));
    PetscCall(MomentumNewtonKrylov_FormResidual(snes, x, f0, &ctx));
    PetscCall(MatMFFDSetBase(J, x, f0));
    PetscCall(MatMult(J, v, jv));
    PetscCall(MatMFFDGetH(J, &h));
    PetscCall(VecWAXPY(xp, h, v, x));
    PetscCall(MomentumNewtonKrylov_FormResidual(snes, xp, fp, &ctx));
    PetscCall(VecWAXPY(fd, -1.0, f0, fp));
    PetscCall(VecScale(fd, 1.0 / h));
    PetscCall(VecAXPY(fd, -1.0, jv));
    PetscCall(VecNorm(fd, NORM_2, &error));
    PetscCall(VecNorm(jv, NORM_2, &scale));
    PetscCall(PicurvAssertBool((PetscBool)(error <= 1.0e-9 * PetscMax(1.0, scale)),
                               "matrix-free Jv must match direct differencing"));

    PetscCall(VecDestroy(&fd));
    PetscCall(VecDestroy(&jv));
    PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&fp));
    PetscCall(VecDestroy(&f0));
    PetscCall(VecDestroy(&xp));
    PetscCall(VecDestroy(&x));
    PetscCall(MatDestroy(&J));
    PetscCall(SNESDestroy(&snes));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Builds a compact all-wall operator fixture through real boundary handlers.
 * @param simCtx Returned simulation context.
 * @param user Returned block context.
 * @param x_periodic Whether the x faces use geometric periodicity.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode BuildMinimalWallOperatorFixture(SimCtx **simCtx, UserCtx **user,
                                                       PetscBool x_periodic)
{
    PetscFunctionBeginUser;
    /* Request the production-width (3) DMDA stencil so the complete RHS is safe
       across MPI partitions; boundary metadata below still selects physical walls. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(
        simCtx, user, 6, 6, 6, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    (*simCtx)->i_periodic = x_periodic ? 1 : 0;
    (*simCtx)->j_periodic = (*simCtx)->k_periodic = 0;
    (*simCtx)->mom_solver_type = MOMENTUM_SOLVER_NEWTON_KRYLOV;
    (*simCtx)->invicid = 1;
    (*simCtx)->dt = 0.1;
    (*simCtx)->step = 1;
    (*simCtx)->StartStep = 0;
    PetscCall(VecSet((*user)->Ucont, 0.0));
    PetscCall(VecSet((*user)->Ucont_o, 0.0));
    PetscCall(VecSet((*user)->Ucont_rm1, 0.0));
    for (PetscInt face = 0; face < 6; ++face) {
        PetscBool periodic_face = (PetscBool)(x_periodic &&
            (face == BC_FACE_NEG_X || face == BC_FACE_POS_X));
        (*user)->boundary_faces[face].face_id = (BCFace)face;
        (*user)->boundary_faces[face].mathematical_type = periodic_face ? PERIODIC : WALL;
        (*user)->boundary_faces[face].handler_type = periodic_face ?
            BC_HANDLER_PERIODIC_GEOMETRIC : BC_HANDLER_WALL_NOSLIP;
        PetscCall(BoundaryCondition_Create((*user)->boundary_faces[face].handler_type,
                                            &(*user)->boundary_faces[face].handler));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Forms the complete direct FD Jacobian, checks every row, and compares MFFD actions.
 */
static PetscErrorCode TestWholeOperatorDirectJacobian(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    SNES snes = NULL;
    Mat J = NULL;
    Vec x = NULL, xp = NULL, f0 = NULL, fp = NULL, column = NULL, square = NULL;
    Vec row_norm_sq = NULL, v[2] = {NULL, NULL}, dense_v[2] = {NULL, NULL};
    Vec mffd_v = NULL, error_vec = NULL;
    PetscInt n_global, lo, hi;
    PetscReal min_row_sq = 0.0, error = 0.0, reference = 0.0;
    const PetscReal epsilon = 1.0e-7;
    MomentumNewtonKrylovContext ctx;

    PetscFunctionBeginUser;
    PetscCall(BuildMinimalWallOperatorFixture(&simCtx, &user, PETSC_FALSE));
    PetscCall(VecDuplicate(user->Ucont, &x)); PetscCall(VecSet(x, 0.0));
    PetscCall(VecDuplicate(x, &xp)); PetscCall(VecDuplicate(x, &f0));
    PetscCall(VecDuplicate(x, &fp)); PetscCall(VecDuplicate(x, &column));
    PetscCall(VecDuplicate(x, &square)); PetscCall(VecDuplicate(x, &row_norm_sq));
    PetscCall(VecDuplicate(x, &v[0])); PetscCall(VecDuplicate(x, &v[1]));
    PetscCall(VecDuplicate(x, &dense_v[0])); PetscCall(VecDuplicate(x, &dense_v[1]));
    PetscCall(VecDuplicate(x, &mffd_v)); PetscCall(VecDuplicate(x, &error_vec));
    PetscCall(VecZeroEntries(row_norm_sq)); PetscCall(VecZeroEntries(dense_v[0]));
    PetscCall(VecZeroEntries(dense_v[1]));
    PetscCall(VecGetSize(x, &n_global));
    PetscCall(VecGetOwnershipRange(x, &lo, &hi));
    for (PetscInt which = 0; which < 2; ++which) {
        PetscScalar *a = NULL;
        PetscCall(VecGetArray(v[which], &a));
        for (PetscInt local = 0; local < hi - lo; ++local) {
            PetscInt global = lo + local;
            a[local] = which == 0 ? (PetscScalar)(1.0 + 0.05 * (global % 9))
                                  : (PetscScalar)(((global % 2) ? -1.0 : 1.0) * (0.5 + 0.03 * (global % 7)));
        }
        PetscCall(VecRestoreArray(v[which], &a));
    }

    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f0, &ctx));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f0, &ctx));
    for (PetscInt col = 0; col < n_global; ++col) {
        PetscScalar coeff[2] = {
            (PetscScalar)(1.0 + 0.05 * (col % 9)),
            (PetscScalar)(((col % 2) ? -1.0 : 1.0) * (0.5 + 0.03 * (col % 7)))
        };
        PetscCall(VecCopy(x, xp));
        if (col >= lo && col < hi) PetscCall(VecSetValue(xp, col, epsilon, ADD_VALUES));
        PetscCall(VecAssemblyBegin(xp)); PetscCall(VecAssemblyEnd(xp));
        PetscCall(MomentumNewtonKrylov_FormResidual(NULL, xp, fp, &ctx));
        PetscCall(VecWAXPY(column, -1.0, f0, fp));
        PetscCall(VecScale(column, 1.0 / epsilon));
        PetscCall(VecPointwiseMult(square, column, column));
        PetscCall(VecAXPY(row_norm_sq, 1.0, square));
        PetscCall(VecAXPY(dense_v[0], coeff[0], column));
        PetscCall(VecAXPY(dense_v[1], coeff[1], column));
    }
    PetscCall(VecMin(row_norm_sq, NULL, &min_row_sq));
    PetscCheck(min_row_sq > 0.5, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Complete Newton Jacobian contains an unexplained zero/weak row (min squared norm=%g).",
               (double)min_row_sq);
    PetscCall(CheckStoredDerivative(user, x, 2, 2, 2, 0, 2, 2, 2, 0,
                                     10.0, 1.0e-5, "interior BDF1 temporal diagonal"));

    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetDM(snes, user->fda));
    PetscCall(SNESSetFunction(snes, f0, MomentumNewtonKrylov_FormResidual, &ctx));
    PetscCall(MatCreateSNESMF(snes, &J));
    PetscCall(MatMFFDSetBase(J, x, f0));
    for (PetscInt which = 0; which < 2; ++which) {
        PetscCall(MatMult(J, v[which], mffd_v));
        PetscCall(VecWAXPY(error_vec, -1.0, dense_v[which], mffd_v));
        PetscCall(VecNorm(error_vec, NORM_2, &error));
        PetscCall(VecNorm(dense_v[which], NORM_2, &reference));
        PetscCheck(error <= 2.0e-5 * PetscMax(1.0, reference), PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "Independent dense FD Jv differs from PETSc MFFD action %d: error=%g reference=%g.",
                   which, (double)error, (double)reference);
    }

    for (PetscInt which = 0; which < 2; ++which) {
        const PetscReal steps[3] = {1.0e-4, 1.0e-6, 1.0e-8};
        PetscReal best = PETSC_MAX_REAL;
        for (PetscInt s = 0; s < 3; ++s) {
            PetscCall(VecWAXPY(xp, steps[s], v[which], x));
            PetscCall(MomentumNewtonKrylov_FormResidual(NULL, xp, fp, &ctx));
            PetscCall(VecWAXPY(column, -1.0, f0, fp));
            PetscCall(VecScale(column, 1.0 / steps[s]));
            PetscCall(VecAXPY(column, -1.0, dense_v[which]));
            PetscCall(VecNorm(column, NORM_2, &error));
            best = PetscMin(best, error);
        }
        PetscCall(VecNorm(dense_v[which], NORM_2, &reference));
        PetscCheck(best <= 2.0e-5 * PetscMax(1.0, reference), PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "Direct directional differences show no accuracy plateau for vector %d (best=%g).",
                   which, (double)best);
    }

    PetscCall(MatDestroy(&J)); PetscCall(SNESDestroy(&snes));
    PetscCall(VecDestroy(&error_vec)); PetscCall(VecDestroy(&mffd_v));
    PetscCall(VecDestroy(&dense_v[1])); PetscCall(VecDestroy(&dense_v[0]));
    PetscCall(VecDestroy(&v[1])); PetscCall(VecDestroy(&v[0]));
    PetscCall(VecDestroy(&row_norm_sq)); PetscCall(VecDestroy(&square));
    PetscCall(VecDestroy(&column)); PetscCall(VecDestroy(&fp)); PetscCall(VecDestroy(&f0));
    PetscCall(VecDestroy(&xp)); PetscCall(VecDestroy(&x));
    PetscCall(BoundarySystem_Destroy(user));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Audits every row of a complete operator containing periodic duplicates. */
static PetscErrorCode TestPeriodicOperatorHasNoZeroRows(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec x = NULL, xp = NULL, f0 = NULL, fp = NULL, column = NULL;
    Vec square = NULL, row_norm_sq = NULL;
    PetscInt n_global, lo, hi;
    PetscReal min_row_sq = 0.0;
    const PetscReal epsilon = 1.0e-7;
    MomentumNewtonKrylovContext ctx;

    PetscFunctionBeginUser;
    PetscCall(BuildMinimalWallOperatorFixture(&simCtx, &user, PETSC_TRUE));
    PetscCall(VecDuplicate(user->Ucont, &x)); PetscCall(VecSet(x, 0.0));
    PetscCall(VecDuplicate(x, &xp)); PetscCall(VecDuplicate(x, &f0));
    PetscCall(VecDuplicate(x, &fp)); PetscCall(VecDuplicate(x, &column));
    PetscCall(VecDuplicate(x, &square)); PetscCall(VecDuplicate(x, &row_norm_sq));
    PetscCall(VecZeroEntries(row_norm_sq));
    PetscCall(VecGetSize(x, &n_global)); PetscCall(VecGetOwnershipRange(x, &lo, &hi));
    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f0, &ctx));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f0, &ctx));
    for (PetscInt col = 0; col < n_global; ++col) {
        PetscCall(VecCopy(x, xp));
        if (col >= lo && col < hi) PetscCall(VecSetValue(xp, col, epsilon, ADD_VALUES));
        PetscCall(VecAssemblyBegin(xp)); PetscCall(VecAssemblyEnd(xp));
        PetscCall(MomentumNewtonKrylov_FormResidual(NULL, xp, fp, &ctx));
        PetscCall(VecWAXPY(column, -1.0, f0, fp)); PetscCall(VecScale(column, 1.0 / epsilon));
        PetscCall(VecPointwiseMult(square, column, column));
        PetscCall(VecAXPY(row_norm_sq, 1.0, square));
    }
    PetscCall(VecMin(row_norm_sq, NULL, &min_row_sq));
    PetscCheck(min_row_sq > 0.5, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Periodic Newton Jacobian contains a zero/weak row (min squared norm=%g).",
               (double)min_row_sq);
    PetscCall(VecDestroy(&row_norm_sq)); PetscCall(VecDestroy(&square));
    PetscCall(VecDestroy(&column)); PetscCall(VecDestroy(&fp)); PetscCall(VecDestroy(&f0));
    PetscCall(VecDestroy(&xp)); PetscCall(VecDestroy(&x));
    PetscCall(BoundarySystem_Destroy(user));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Exercises a converged solve, forced rollback, and per-call cleanup. */
static PetscErrorCode TestSmallSolveAndRollback(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec entry = NULL, delta = NULL;
    PetscErrorCode solve_ierr;
    PetscReal norm = 0.0;
    const char *fields[] = {"Ucont"};
    char summary_path[PETSC_MAX_PATH_LEN];
    char history_path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_rtol", "1e-4"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_max_it", "20"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_rtol", "1e-6"));
    simCtx->mom_nk_monitor_history = PETSC_TRUE;
    PetscCall(MomentumSolver_NewtonKrylov(user, NULL, NULL));
    PetscCall(PicurvAssertBool(simCtx->mom_last_converged, "small Newton solve must converge"));
    PetscCall(PicurvAssertBool((PetscBool)(user->Rhs == NULL), "successful solve must release Rhs"));
    PetscCall(PetscSNPrintf(summary_path, sizeof(summary_path),
                            "%s/Momentum_Solver_Newton_Krylov_Summary_Block_0.log", simCtx->log_dir));
    PetscCall(PetscSNPrintf(history_path, sizeof(history_path),
                            "%s/Momentum_Solver_Newton_Krylov_History_Block_0.log", simCtx->log_dir));
    PetscCall(AssertNewtonLog(summary_path, 1, "solver: Newton Krylov", "state: committed"));
    PetscCall(AssertNewtonLog(history_path, -2, "newton: 0", "nonlinear_norm:"));

    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    tmpdir[0] = '\0';
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecSet(user->Ucont, 0.2));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(VecDuplicate(user->Ucont, &entry));
    PetscCall(VecDuplicate(user->Ucont, &delta));
    PetscCall(VecCopy(user->Ucont, entry));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_max_it", "0"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_CONV_FAILED, solve_ierr,
                                   "forced nonconvergence must report PETSC_ERR_CONV_FAILED"));
    PetscCall(VecWAXPY(delta, -1.0, entry, user->Ucont));
    PetscCall(VecNorm(delta, NORM_INFINITY, &norm));
    PetscCall(PicurvAssertRealNear(0.0, norm, 1.0e-12,
                                   "failed Newton solve must restore the canonical entry state"));
    PetscCall(PicurvAssertBool((PetscBool)(user->Rhs == NULL), "failed solve must release Rhs"));
    PetscCall(PetscSNPrintf(summary_path, sizeof(summary_path),
                            "%s/Momentum_Solver_Newton_Krylov_Summary_Block_0.log", simCtx->log_dir));
    PetscCall(AssertNewtonLog(summary_path, 1, "reason_code: -", "state: rolled_back"));

    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_snes_rtol"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_snes_max_it"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_ksp_rtol"));
    PetscCall(VecDestroy(&delta));
    PetscCall(VecDestroy(&entry));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies the six-wall zero-velocity case logs zero Newton/Krylov work. */
static PetscErrorCode TestZeroIterationStructuredLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    char summary_path[PETSC_MAX_PATH_LEN];
    const char *fields[] = {"Ucont"};

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecZeroEntries(user->Ucont));
    PetscCall(VecZeroEntries(user->Ucont_o));
    PetscCall(VecZeroEntries(user->Ucont_rm1));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(MomentumSolver_NewtonKrylov(user, NULL, NULL));
    PetscCall(PetscSNPrintf(summary_path, sizeof(summary_path),
                            "%s/Momentum_Solver_Newton_Krylov_Summary_Block_0.log", simCtx->log_dir));
    PetscCall(AssertNewtonLog(summary_path, 1, "newton: 0 | evals: 1 | krylov: 0",
                             "final: 0.0000000000000000e+00 | state: committed"));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Confirms unsupported features fail before workspace allocation. */
static PetscErrorCode TestUnsupportedConfigurationFailsBeforeAllocation(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    PetscErrorCode solve_ierr;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    {
        PetscInt *unsupported_flags[] = {
            &simCtx->immersed, &simCtx->movefsi, &simCtx->rotatefsi,
            &simCtx->moveframe, &simCtx->rotateframe, &simCtx->rans,
            &simCtx->clark, &simCtx->TwoD, &simCtx->wallfunction
        };
        for (size_t flag = 0; flag < sizeof(unsupported_flags) / sizeof(unsupported_flags[0]); ++flag) {
            *unsupported_flags[flag] = 1;
            PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
            solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
            PetscCall(PetscPopErrorHandler());
            PetscCall(PicurvAssertBool((PetscBool)(solve_ierr != PETSC_SUCCESS),
                                       "unsupported Newton feature flag must fail"));
            PetscCall(PicurvAssertBool((PetscBool)(user->Rhs == NULL),
                                       "feature validation must precede workspace allocation"));
            *unsupported_flags[flag] = 0;
        }
    }
    simCtx->block_number = 2;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(solve_ierr != PETSC_SUCCESS), "multiblock must fail"));
    simCtx->block_number = 1;
    simCtx->StartStep = 1;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(solve_ierr != PETSC_SUCCESS), "restart must fail"));
    simCtx->StartStep = 0;
    PetscCall(VecSet(user->Nvert, 1.0));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(solve_ierr != PETSC_SUCCESS), "masked rows must fail"));
    PetscCall(VecSet(user->Nvert, 0.0));
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(solve_ierr != PETSC_SUCCESS),
                               "driven constant-flux controller must fail"));
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = INLET;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_INLET_INTERP_FROM_FILE;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(solve_ierr != PETSC_SUCCESS),
                               "unimplemented interpolated-file inlet must fail"));
    PetscCall(PicurvAssertBool((PetscBool)(user->Rhs == NULL),
                               "all validation failures must precede workspace allocation"));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies cleanup and rollback after an options failure following asset creation. */
static PetscErrorCode TestPostAllocationFailureCleanup(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec entry = NULL, delta = NULL;
    PetscErrorCode solve_ierr;
    PetscReal norm = 0.0;
    const char *fields[] = {"Ucont"};

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecSet(user->Ucont, 0.2));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(VecDuplicate(user->Ucont, &entry)); PetscCall(VecCopy(user->Ucont, entry));
    PetscCall(VecDuplicate(user->Ucont, &delta));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_pc_type", "jacobi"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    solve_ierr = MomentumSolver_NewtonKrylov(user, NULL, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_pc_type"));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_SUP, solve_ierr,
                                   "non-PCNONE option must fail after setup"));
    PetscCall(PicurvAssertBool((PetscBool)(user->Rhs == NULL),
                               "post-allocation failure must destroy Rhs"));
    PetscCall(VecWAXPY(delta, -1.0, entry, user->Ucont));
    PetscCall(VecNorm(delta, NORM_INFINITY, &norm));
    PetscCall(PicurvAssertRealNear(0.0, norm, 1.0e-12,
                                   "post-allocation failure must restore canonical entry"));
    PetscCall(VecDestroy(&delta)); PetscCall(VecDestroy(&entry));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Runs the focused Newton--Krylov unit suite.
 * @param argc Command-line argument count.
 * @param argv Command-line argument vector.
 * @return Process exit status.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"residual-repeatability-and-input-integrity", TestResidualRepeatabilityAndInputIntegrity},
        {"constraint-rows", TestConstraintRows},
        {"fixed-constraint-derivatives-all-faces", TestFixedConstraintDerivativesAllFaces},
        {"inlet-outlet-constraint-derivatives", TestInletOutletConstraintDerivatives},
        {"periodic-constraint-derivatives-and-intersections", TestPeriodicConstraintDerivativesAndIntersections},
        {"matrix-free-derivative", TestMatrixFreeDerivative},
        {"whole-operator-direct-jacobian", TestWholeOperatorDirectJacobian},
        {"periodic-operator-has-no-zero-rows", TestPeriodicOperatorHasNoZeroRows},
        {"zero-iteration-structured-logging", TestZeroIterationStructuredLogging},
        {"small-solve-and-rollback", TestSmallSolveAndRollback},
        {"unsupported-configuration-fails-before-allocation", TestUnsupportedConfigurationFailsBeforeAllocation},
        {"post-allocation-failure-cleanup", TestPostAllocationFailureCleanup},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv Newton Krylov tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-newton-krylov", cases, sizeof(cases) / sizeof(cases[0]));
    if (PetscFinalize()) return 1;
    return (int)ierr;
}
