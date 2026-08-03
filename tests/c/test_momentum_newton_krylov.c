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

/**
 * @brief Exercises the straight-duct BDF1 startup path used by flat_channel.
 *
 * The conservation outlet consumes lUcat during its first boundary pass.  This
 * test deliberately evaluates the callback at the initialized state before
 * installing SNES, then completes the first nonlinear solve with each shipped
 * preconditioner.  It catches a missing Ucont -> Ucat -> lUcat seed as a
 * non-finite initial residual rather than hiding it behind later MFFD work.
 */
static PetscErrorCode CheckFlatChannelStartup(PetscBool use_point_block)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    Vec x = NULL, f = NULL;
    PetscReal initial_norm = 0.0;
    MomentumNewtonKrylovContext ctx = {0};

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(NULL, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &f));
    PetscCall(VecCopy(user->Ucont, x));
    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f, &ctx));
    PetscCall(VecNorm(f, NORM_2, &initial_norm));
    PetscCheck(!PetscIsInfOrNanReal(initial_norm), PETSC_COMM_WORLD, PETSC_ERR_FP,
               "flat-channel BDF1 initial residual is non-finite (%g) with %s.",
               (double)initial_norm,
               use_point_block ? "frozen-momentum point-block" : "PCNONE");
    PetscCall(VecDestroy(&f));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_rtol", "1e-4"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_max_it", "20"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_rtol", "1e-6"));
    if (use_point_block) {
        PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model",
                                       "frozen_momentum_jacobian"));
        PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure",
                                       "point_block"));
    }
    PetscCall(MomentumSolver_NewtonKrylov(user, NULL, NULL));
    PetscCall(PicurvAssertBool(simCtx->mom_last_converged,
                               "flat-channel BDF1 Newton solve must converge"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_preconditioner_model"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_preconditioner_structure"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_snes_rtol"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_snes_max_it"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_ksp_rtol"));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Guards flat_channel's initial BDF1 residual and both shipped NK PCs. */
static PetscErrorCode TestFlatChannelStartup(void)
{
    PetscFunctionBeginUser;
    PetscCall(CheckFlatChannelStartup(PETSC_FALSE));
    PetscCall(CheckFlatChannelStartup(PETSC_TRUE));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies restarted Newton solves with both supported preconditioners. */
static PetscErrorCode TestRestartAndContinuationSolve(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    const char *fields[] = {"Ucont"};

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecSet(user->Ucont, 0.2));
    PetscCall(VecZeroEntries(user->Ucont_o));
    PetscCall(VecZeroEntries(user->Ucont_rm1));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_rtol", "1e-4"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_max_it", "20"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_rtol", "1e-6"));

    /* AdvanceSimulation() solves StartStep+1 first.  Use a nonzero checkpoint
     * state so that SNES takes a Newton/Krylov path rather than accepting a
     * trivial residual. */
    simCtx->StartStep = 7;
    simCtx->step = simCtx->StartStep + 1;
    PetscCall(MomentumSolver_NewtonKrylov(user, NULL, NULL));
    PetscCall(PicurvAssertBool(simCtx->mom_last_converged,
                               "PCNONE Newton Krylov checkpoint restart must converge"));

    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    tmpdir[0] = '\0';
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecSet(user->Ucont, 0.2));
    PetscCall(VecZeroEntries(user->Ucont_o));
    PetscCall(VecZeroEntries(user->Ucont_rm1));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model",
                                   "frozen_momentum_jacobian"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "point_block"));

    simCtx->continueMode = PETSC_TRUE;
    simCtx->StartStep = 8;
    simCtx->step = simCtx->StartStep + 1;
    PetscCall(MomentumSolver_NewtonKrylov(user, NULL, NULL));
    PetscCall(PicurvAssertBool(simCtx->mom_last_converged,
                               "frozen-momentum Newton Krylov --continue solve must converge"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_preconditioner_model"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_preconditioner_structure"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_snes_rtol"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_snes_max_it"));
    PetscCall(PetscOptionsClearValue(NULL, "-mom_nk_ksp_rtol"));
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

/** @brief Verifies finalized application-owned linearization option parsing. */
static PetscErrorCode TestLinearizationConfigParsing(void)
{
    MomentumNewtonJacobian jacobian = {0};
    MomentumPreconditionerDescription description = {0};
    PetscErrorCode config_ierr;
    const char *option_names[] = {
        "-mom_nk_jacobian_type",
        "-mom_nk_jacobian_fd_mode",
        "-mom_nk_preconditioner_model",
        "-mom_nk_preconditioner_structure"
    };

    PetscFunctionBeginUser;
    for (size_t n = 0; n < sizeof(option_names) / sizeof(option_names[0]); ++n)
        PetscCall(PetscOptionsClearValue(NULL, option_names[n]));

    PetscCall(MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description));
    PetscCall(PicurvAssertIntEqual(MOM_NK_JACOBIAN_FINITE_DIFFERENCE, jacobian.type,
                                   "default Jacobian type must be finite difference"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_FD_MODE_MATRIX_FREE,
                                   jacobian.finite_difference_mode,
                                   "default finite-difference mode must be matrix free"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_PC_MODEL_NONE, description.model,
                                   "default preconditioner model must be none"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_PC_STRUCTURE_NONE, description.structure,
                                   "default preconditioner structure must be none"));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_jacobian_type", "finite_difference"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_jacobian_fd_mode", "matrix_free"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model", "none"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "none"));
    PetscCall(MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description));
    PetscCall(PicurvAssertIntEqual(MOM_NK_JACOBIAN_FINITE_DIFFERENCE, jacobian.type,
                                   "explicit baseline Jacobian type must parse"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_FD_MODE_MATRIX_FREE,
                                   jacobian.finite_difference_mode,
                                   "explicit baseline finite-difference mode must parse"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_PC_MODEL_NONE, description.model,
                                   "explicit baseline preconditioner model must parse"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_PC_STRUCTURE_NONE, description.structure,
                                   "explicit baseline preconditioner structure must parse"));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model",
                                  "frozen_momentum_jacobian"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "point_block"));
    PetscCall(MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description));
    PetscCall(PicurvAssertIntEqual(MOM_NK_JACOBIAN_FINITE_DIFFERENCE, jacobian.type,
                                   "explicit Jacobian type must parse"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_FD_MODE_MATRIX_FREE,
                                   jacobian.finite_difference_mode,
                                   "explicit finite-difference mode must parse"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN,
                                   description.model, "frozen model must parse"));
    PetscCall(PicurvAssertIntEqual(MOM_NK_PC_STRUCTURE_POINT_BLOCK,
                                   description.structure, "point-block structure must parse"));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_jacobian_type", "analytic"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    config_ierr = MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, config_ierr,
                                   "unsupported Jacobian type must fail"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_jacobian_type", "finite_difference"));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_jacobian_fd_mode", "colored_sparse"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    config_ierr = MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, config_ierr,
                                   "unsupported finite-difference mode must fail"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_jacobian_fd_mode", "matrix_free"));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "none"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    config_ierr = MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_SUP, config_ierr,
                                   "frozen model without point block must fail"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model", "none"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "point_block"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    config_ierr = MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_SUP, config_ierr,
                                   "none model with point block must fail"));

    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model", "diagonal"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "none"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    config_ierr = MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, config_ierr,
                                   "unsupported preconditioner model must fail"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_model",
                                  "frozen_momentum_jacobian"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_preconditioner_structure", "line"));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    config_ierr = MomentumNewtonKrylov_ReadLinearizationConfig(&jacobian, &description);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, config_ierr,
                                   "unsupported preconditioner structure must fail"));

    for (size_t n = 0; n < sizeof(option_names) / sizeof(option_names[0]); ++n)
        PetscCall(PetscOptionsClearValue(NULL, option_names[n]));
    PetscFunctionReturn(PETSC_SUCCESS);
}

enum {
    ORACLE_CENTER_AJ = 1 << 0,
    ORACLE_CENTER_TRANSVERSE_METRICS = 1 << 1,
    ORACLE_CENTER_VELOCITY = 1 << 2,
    ORACLE_OMIT_A5 = 1 << 3,
    ORACLE_LEGACY_SIGN = 1 << 4
};

/** @brief Test-owned metric norm used by the independent legacy transcription. */
static PetscReal LegacyOracleMetricNormSquared(Cmpnts metric)
{
    return metric.x * metric.x + metric.y * metric.y + metric.z * metric.z;
}

/**
 * @brief Independent transcription of the audited legacy mode-2 point block.
 *
 * This intentionally shares no coefficient helper with production. Mutant flags
 * represent the historical failure modes that the nonuniform oracle must reject.
 */
static void LegacyPointBlockOracle(const SimCtx *simCtx, const Cmpnts ***u,
    const Cmpnts ***csi, const Cmpnts ***eta, const Cmpnts ***zet,
    const PetscReal ***aj, PetscInt i, PetscInt j, PetscInt k, PetscInt flags,
    PetscScalar block[9])
{
    PetscReal A[6][4] = {{0.0}};
    const PetscReal dtc = ((simCtx->step != simCtx->StartStep) && simCtx->step != 1 ? 1.5 : 1.0) /
                          simCtx->dt;
    const PetscReal AJip = flags & ORACLE_CENTER_AJ ? aj[k][j][i] :
                           0.5 * (aj[k][j][i] + aj[k][j][i + 1]);
    const PetscReal AJjp = flags & ORACLE_CENTER_AJ ? aj[k][j][i] :
                           0.5 * (aj[k][j][i] + aj[k][j + 1][i]);
    const PetscReal AJkp = flags & ORACLE_CENTER_AJ ? aj[k][j][i] :
                           0.5 * (aj[k][j][i] + aj[k + 1][j][i]);
    const PetscReal g11ip = LegacyOracleMetricNormSquared(csi[k][j][i]);
    const PetscReal g22ip = flags & ORACLE_CENTER_TRANSVERSE_METRICS ?
        LegacyOracleMetricNormSquared(eta[k][j][i]) : 0.25 * (
        LegacyOracleMetricNormSquared(eta[k][j][i]) + LegacyOracleMetricNormSquared(eta[k][j][i + 1]) +
        LegacyOracleMetricNormSquared(eta[k][j - 1][i]) + LegacyOracleMetricNormSquared(eta[k][j - 1][i + 1]));
    const PetscReal g33ip = flags & ORACLE_CENTER_TRANSVERSE_METRICS ?
        LegacyOracleMetricNormSquared(zet[k][j][i]) : 0.25 * (
        LegacyOracleMetricNormSquared(zet[k][j][i]) + LegacyOracleMetricNormSquared(zet[k][j][i + 1]) +
        LegacyOracleMetricNormSquared(zet[k - 1][j][i]) + LegacyOracleMetricNormSquared(zet[k - 1][j][i + 1]));
    const PetscReal g11jp = flags & ORACLE_CENTER_TRANSVERSE_METRICS ?
        LegacyOracleMetricNormSquared(csi[k][j][i]) : 0.25 * (
        LegacyOracleMetricNormSquared(csi[k][j][i]) + LegacyOracleMetricNormSquared(csi[k][j + 1][i]) +
        LegacyOracleMetricNormSquared(csi[k][j][i - 1]) + LegacyOracleMetricNormSquared(csi[k][j + 1][i - 1]));
    const PetscReal g22jp = LegacyOracleMetricNormSquared(eta[k][j][i]);
    const PetscReal g33jp = flags & ORACLE_CENTER_TRANSVERSE_METRICS ?
        LegacyOracleMetricNormSquared(zet[k][j][i]) : 0.25 * (
        LegacyOracleMetricNormSquared(zet[k][j][i]) + LegacyOracleMetricNormSquared(zet[k][j + 1][i]) +
        LegacyOracleMetricNormSquared(zet[k - 1][j][i]) + LegacyOracleMetricNormSquared(zet[k - 1][j + 1][i]));
    const PetscReal g11kp = flags & ORACLE_CENTER_TRANSVERSE_METRICS ?
        LegacyOracleMetricNormSquared(csi[k][j][i]) : 0.25 * (
        LegacyOracleMetricNormSquared(csi[k][j][i]) + LegacyOracleMetricNormSquared(csi[k + 1][j][i]) +
        LegacyOracleMetricNormSquared(csi[k][j][i - 1]) + LegacyOracleMetricNormSquared(csi[k + 1][j][i - 1]));
    const PetscReal g22kp = flags & ORACLE_CENTER_TRANSVERSE_METRICS ?
        LegacyOracleMetricNormSquared(eta[k][j][i]) : 0.25 * (
        LegacyOracleMetricNormSquared(eta[k][j][i]) + LegacyOracleMetricNormSquared(eta[k + 1][j][i]) +
        LegacyOracleMetricNormSquared(eta[k][j - 1][i]) + LegacyOracleMetricNormSquared(eta[k + 1][j - 1][i]));
    const PetscReal g33kp = LegacyOracleMetricNormSquared(zet[k][j][i]);
    const PetscReal U0jp = flags & ORACLE_CENTER_VELOCITY ? u[k][j][i].x : 0.25 *
        (u[k][j][i].x + u[k][j][i - 1].x + u[k][j + 1][i].x + u[k][j + 1][i - 1].x);
    const PetscReal U0kp = flags & ORACLE_CENTER_VELOCITY ? u[k][j][i].x : 0.25 *
        (u[k][j][i].x + u[k][j][i - 1].x + u[k + 1][j][i].x + u[k + 1][j][i - 1].x);
    const PetscReal U1ip = flags & ORACLE_CENTER_VELOCITY ? u[k][j][i].y : 0.25 *
        (u[k][j][i].y + u[k][j - 1][i].y + u[k][j][i + 1].y + u[k][j - 1][i + 1].y);
    const PetscReal U1kp = flags & ORACLE_CENTER_VELOCITY ? u[k][j][i].y : 0.25 *
        (u[k][j][i].y + u[k][j - 1][i].y + u[k + 1][j][i].y + u[k + 1][j - 1][i].y);
    const PetscReal U2ip = flags & ORACLE_CENTER_VELOCITY ? u[k][j][i].z : 0.25 *
        (u[k][j][i].z + u[k - 1][j][i].z + u[k][j][i + 1].z + u[k - 1][j][i + 1].z);
    const PetscReal U2jp = flags & ORACLE_CENTER_VELOCITY ? u[k][j][i].z : 0.25 *
        (u[k][j][i].z + u[k - 1][j][i].z + u[k][j + 1][i].z + u[k - 1][j + 1][i].z);
    PetscReal Su, Sv, Sw, sign = flags & ORACLE_LEGACY_SIGN ? -1.0 : 1.0;

    A[0][0] =  .125 * aj[k][j][i] * u[k][j][i].y;
    A[0][1] = -.125 * aj[k][j - 1][i] * u[k][j - 1][i].y;
    A[0][2] =  .125 * aj[k][j][i + 1] * u[k][j][i + 1].y;
    A[0][3] = -.125 * aj[k][j - 1][i + 1] * u[k][j - 1][i + 1].y;
    A[1][0] =  .125 * aj[k][j][i] * u[k][j][i].z;
    A[1][1] = -.125 * aj[k - 1][j][i] * u[k - 1][j][i].z;
    A[1][2] =  .125 * aj[k][j][i + 1] * u[k][j][i + 1].z;
    A[1][3] = -.125 * aj[k - 1][j][i + 1] * u[k - 1][j][i + 1].z;
    A[2][0] = -.125 * aj[k][j + 1][i - 1] * u[k][j + 1][i - 1].x;
    A[2][1] = -.125 * aj[k][j][i - 1] * u[k][j][i - 1].x;
    A[2][2] =  .125 * aj[k][j + 1][i] * u[k][j + 1][i].x;
    A[2][3] =  .125 * aj[k][j][i] * u[k][j][i].x;
    A[3][0] =  .125 * aj[k][j][i] * u[k][j][i].z;
    A[3][1] = -.125 * aj[k - 1][j][i] * u[k - 1][j][i].z;
    A[3][2] =  .125 * aj[k][j + 1][i] * u[k][j + 1][i].z;
    A[3][3] = -.125 * aj[k - 1][j + 1][i] * u[k - 1][j + 1][i].z;
    A[4][0] = -.125 * aj[k + 1][j][i - 1] * u[k + 1][j][i - 1].x;
    A[4][1] = -.125 * aj[k][j][i - 1] * u[k][j][i - 1].x;
    A[4][2] =  .125 * aj[k + 1][j][i] * u[k + 1][j][i].x;
    A[4][3] =  .125 * aj[k][j][i] * u[k][j][i].x;
    A[5][0] = -.125 * aj[k + 1][j - 1][i] * u[k + 1][j - 1][i].y;
    A[5][1] = -.125 * aj[k][j - 1][i] * u[k][j - 1][i].y;
    A[5][2] =  .125 * aj[k + 1][j][i] * u[k + 1][j][i].y;
    A[5][3] =  .125 * aj[k][j][i] * u[k][j][i].y;
    Su = A[0][0] + A[0][1] + A[0][2] + A[0][3] + A[1][0] + A[1][1] + A[1][2] + A[1][3];
    Sv = A[2][0] + A[2][1] + A[2][2] + A[2][3] + A[3][0] + A[3][1] + A[3][2] + A[3][3];
    Sw = A[4][0] + A[4][1] + A[4][2] + A[4][3];
    if (!(flags & ORACLE_OMIT_A5)) Sw += A[5][0] + A[5][1] + A[5][2] + A[5][3];

    block[0] = sign * (dtc + AJip * AJip * (g11ip + g22ip + g33ip) / simCtx->ren + Su);
    block[1] = sign * 0.5 * AJip * U1ip; block[2] = sign * 0.5 * AJip * U2ip;
    block[3] = sign * 0.5 * AJjp * U0jp;
    block[4] = sign * (dtc + AJjp * AJjp * (g11jp + g22jp + g33jp) / simCtx->ren + Sv);
    block[5] = sign * 0.5 * AJjp * U2jp;
    block[6] = sign * 0.5 * AJkp * U0kp; block[7] = sign * 0.5 * AJkp * U1kp;
    block[8] = sign * (dtc + AJkp * AJkp * (g11kp + g22kp + g33kp) / simCtx->ren + Sw);
}

/** @brief Seeds nonuniform, index-distinguishing coefficient fields. */
static PetscErrorCode SeedPointBlockOracleFields(UserCtx *user)
{
    Cmpnts ***u = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;
    PetscReal ***aj = NULL;
    DMDALocalInfo info = user->info;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &u));
    PetscCall(DMDAVecGetArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecGetArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecGetArray(user->fda, user->Zet, &zet));
    PetscCall(DMDAVecGetArray(user->da, user->Aj, &aj));
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                aj[k][j][i] = 0.7 + .019 * i + .043 * j + .071 * k + .003 * i * j + .002 * j * k;
                u[k][j][i] = (Cmpnts){.x = .2 + .031 * i - .017 * j + .013 * k + .004 * i * k,
                                      .y = -.3 + .011 * i + .037 * j - .019 * k + .003 * j * k,
                                      .z = .4 - .023 * i + .007 * j + .041 * k + .002 * i * j};
                csi[k][j][i] = (Cmpnts){.x = 1.1 + .029 * i + .007 * j * k,
                                        .y = .13 + .017 * j + .003 * i * k,
                                        .z = -.09 + .011 * k + .002 * i * j};
                eta[k][j][i] = (Cmpnts){.x = -.12 + .013 * i + .004 * j * k,
                                        .y = .9 + .031 * j + .003 * i * k,
                                        .z = .16 + .019 * k + .002 * i * j};
                zet[k][j][i] = (Cmpnts){.x = .08 + .023 * i + .002 * j * k,
                                        .y = -.14 + .011 * j + .005 * i * k,
                                        .z = 1.2 + .037 * k + .003 * i * j};
            }
    PetscCall(DMDAVecRestoreArray(user->da, user->Aj, &aj));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Zet, &zet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &u));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocalBegin(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscCall(DMGlobalToLocalEnd(user->da, user->Aj, INSERT_VALUES, user->lAj));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Evaluates the independent oracle on the unique owner and broadcasts it. */
static PetscErrorCode CollectiveLegacyPointBlockOracle(UserCtx *user, PetscInt i, PetscInt j,
    PetscInt k, PetscInt flags, PetscScalar block[9])
{
    Cmpnts ***u = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;
    PetscReal ***aj = NULL;
    PetscInt owns = i >= user->info.xs && i < user->info.xs + user->info.xm &&
                    j >= user->info.ys && j < user->info.ys + user->info.ym &&
                    k >= user->info.zs && k < user->info.zs + user->info.zm;
    PetscInt owners = 0;
    PetscScalar local[9] = {0.0};

    PetscFunctionBeginUser;
    if (owns) {
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &u));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, &csi));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, &eta));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, &zet));
        PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
        LegacyPointBlockOracle(user->simCtx, (const Cmpnts ***)u, (const Cmpnts ***)csi,
            (const Cmpnts ***)eta, (const Cmpnts ***)zet, (const PetscReal ***)aj,
            i, j, k, flags, local);
        PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &u));
    }
    PetscCallMPI(MPI_Allreduce(&owns, &owners, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(1, owners, "oracle point must have exactly one owner"));
    PetscCallMPI(MPI_Allreduce(local, block, 9, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Converts an in-domain or periodic-ghost DMDA stencil to PETSc ordering. */
static PetscErrorCode TestStencilToGlobal(UserCtx *user, MatStencil stencil,
                                          PetscInt *global_index)
{
    ISLocalToGlobalMapping local_to_global = NULL;
    PetscInt ghost_starts[3], ghost_sizes[3], local_index;

    PetscFunctionBeginUser;
    PetscCall(DMDAGetGhostCorners(user->fda,
                                  &ghost_starts[0], &ghost_starts[1], &ghost_starts[2],
                                  &ghost_sizes[0], &ghost_sizes[1], &ghost_sizes[2]));
    PetscCall(PicurvAssertBool((PetscBool)(
        stencil.i >= ghost_starts[0] && stencil.i < ghost_starts[0] + ghost_sizes[0] &&
        stencil.j >= ghost_starts[1] && stencil.j < ghost_starts[1] + ghost_sizes[1] &&
        stencil.k >= ghost_starts[2] && stencil.k < ghost_starts[2] + ghost_sizes[2]),
        "test stencil must lie in the local DMDA ghost region"));
    local_index = stencil.c + 3 * (
        (stencil.i - ghost_starts[0]) + ghost_sizes[0] * (
        (stencil.j - ghost_starts[1]) + ghost_sizes[1] *
        (stencil.k - ghost_starts[2])));
    PetscCall(DMGetLocalToGlobalMapping(user->fda, &local_to_global));
    PetscCall(ISLocalToGlobalMappingApply(local_to_global, 1, &local_index,
                                          global_index));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Reads one DMDA-stencil matrix entry through collective basis vectors. */
static PetscErrorCode PreconditionerMatrixStencilEntry(UserCtx *user,
                                       Mat preconditioning_matrix, MatStencil row,
                                       MatStencil col, PetscScalar *value)
{
    Vec column_basis = NULL, row_basis = NULL, product = NULL;
    AO ao = NULL;
    PetscInt mx, my, row_index, col_index;
    PetscFunctionBeginUser;
    PetscCall(DMDAGetInfo(user->fda, NULL, &mx, &my, NULL, NULL, NULL, NULL,
                          NULL, NULL, NULL, NULL, NULL, NULL));
    row_index = row.c + 3 * (row.i + mx * (row.j + my * row.k));
    col_index = col.c + 3 * (col.i + mx * (col.j + my * col.k));
    PetscCall(DMDAGetAO(user->fda, &ao));
    PetscCall(AOApplicationToPetsc(ao, 1, &row_index));
    PetscCall(AOApplicationToPetsc(ao, 1, &col_index));
    PetscCall(DMCreateGlobalVector(user->fda, &column_basis));
    PetscCall(DMCreateGlobalVector(user->fda, &row_basis));
    PetscCall(DMCreateGlobalVector(user->fda, &product));
    PetscCall(VecSet(column_basis, 0.0)); PetscCall(VecSet(row_basis, 0.0));
    PetscCall(VecSetValue(column_basis, col_index, 1.0, INSERT_VALUES));
    PetscCall(VecSetValue(row_basis, row_index, 1.0, INSERT_VALUES));
    PetscCall(VecAssemblyBegin(column_basis)); PetscCall(VecAssemblyEnd(column_basis));
    PetscCall(VecAssemblyBegin(row_basis)); PetscCall(VecAssemblyEnd(row_basis));
    PetscCall(MatMult(preconditioning_matrix, column_basis, product));
    PetscCall(VecDot(row_basis, product, value));
    PetscCall(VecDestroy(&product)); PetscCall(VecDestroy(&row_basis));
    PetscCall(VecDestroy(&column_basis));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies the exact AIJ layout and preallocation derived from row classes. */
static PetscErrorCode AssertExactPointBlockMatrixAllocation(
    UserCtx *user, Mat matrix, PetscBool require_offrank_periodic)
{
    DMDALocalInfo info;
    MatInfo matrix_info;
    PetscMPIInt comm_size = 1;
    PetscInt matrix_rows, matrix_cols, local_rows, local_cols;
    PetscInt vector_size, vector_local_size, block_size;
    PetscInt ownership_start, ownership_end;
    PetscInt expected_local = 0, expected_global = 0;
    PetscInt offrank_periodic_local = 0, offrank_periodic_global = 0;
    PetscBool is_seq_aij = PETSC_FALSE, is_mpi_aij = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(DMDAGetLocalInfo(user->fda, &info));
    PetscCall(MatGetSize(matrix, &matrix_rows, &matrix_cols));
    PetscCall(MatGetLocalSize(matrix, &local_rows, &local_cols));
    PetscCall(VecGetSize(user->Ucont, &vector_size));
    PetscCall(VecGetLocalSize(user->Ucont, &vector_local_size));
    PetscCall(VecGetOwnershipRange(user->Ucont, &ownership_start, &ownership_end));
    PetscCall(MatGetBlockSize(matrix, &block_size));
    PetscCall(PicurvAssertIntEqual(vector_size, matrix_rows,
                                   "point-block global row dimension"));
    PetscCall(PicurvAssertIntEqual(vector_size, matrix_cols,
                                   "point-block global column dimension"));
    PetscCall(PicurvAssertIntEqual(vector_local_size, local_rows,
                                   "point-block local row dimension"));
    PetscCall(PicurvAssertIntEqual(vector_local_size, local_cols,
                                   "point-block local column dimension"));
    PetscCall(PicurvAssertIntEqual(3, block_size, "point-block logical block size"));
    PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)matrix), &comm_size));
    PetscCall(PetscObjectTypeCompare((PetscObject)matrix, MATSEQAIJ, &is_seq_aij));
    PetscCall(PetscObjectTypeCompare((PetscObject)matrix, MATMPIAIJ, &is_mpi_aij));
    PetscCall(PicurvAssertBool(
        comm_size == 1 ? is_seq_aij : is_mpi_aij,
        "point-block matrix must use the expected AIJ implementation"));

    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                for (PetscInt component = 0; component < 3; ++component) {
                    PetscInt ri, rj, rk;
                    MomentumNewtonKrylovRowType type = MomentumNewtonKrylov_ClassifyRow(
                        user, i, j, k, component, &ri, &rj, &rk);
                    expected_local += type == MOM_NK_ROW_PHYSICAL ? 3 :
                                      type == MOM_NK_ROW_PERIODIC_DUPLICATE ? 2 : 1;
                    if (type == MOM_NK_ROW_PERIODIC_DUPLICATE) {
                        PetscInt representative;
                        PetscCall(TestStencilToGlobal(user,
                            (MatStencil){.i = ri, .j = rj, .k = rk, .c = component},
                            &representative));
                        if (representative < ownership_start || representative >= ownership_end)
                            ++offrank_periodic_local;
                    }
                }
            }
        }
    }
    PetscCallMPI(MPI_Allreduce(&expected_local, &expected_global, 1, MPIU_INT,
                               MPI_SUM, PetscObjectComm((PetscObject)matrix)));
    PetscCallMPI(MPI_Allreduce(&offrank_periodic_local, &offrank_periodic_global,
                               1, MPIU_INT, MPI_SUM,
                               PetscObjectComm((PetscObject)matrix)));
    PetscCall(MatGetInfo(matrix, MAT_GLOBAL_SUM, &matrix_info));
    PetscCall(PicurvAssertRealNear((PetscReal)expected_global,
        (PetscReal)matrix_info.nz_allocated, 0.0,
        "point-block matrix must allocate exactly the classified scalar pattern"));
    PetscCall(PicurvAssertRealNear((PetscReal)expected_global,
        (PetscReal)matrix_info.nz_used, 0.0,
        "point-block assembly must insert every classified structural entry"));
    PetscCall(PicurvAssertRealNear(0.0, (PetscReal)matrix_info.mallocs, 0.0,
        "point-block insertion must not reallocate matrix storage"));
    if (require_offrank_periodic && comm_size > 1)
        PetscCall(PicurvAssertBool((PetscBool)(offrank_periodic_global > 0),
            "MPI periodic fixture must exercise off-rank preallocation"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Verifies the point-block model and common preconditioning-engine wiring.
 * @return PETSc error code.
 */
static PetscErrorCode TestPointBlockPreconditionerEngine(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    MomentumNewtonKrylovContext ctx = {0};
    MomentumPreconditionerDescription description = {
        MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN,
        MOM_NK_PC_STRUCTURE_POINT_BLOCK,
        0,
        0
    };
    MomentumPreconditionerEngine engine = {0};
    Mat preconditioning_matrix = NULL;
    Vec x = NULL, f = NULL;
    PetscInt block_size = 0, velocity_dof = 0;
    PetscReal matrix_norm = 0.0, reassembled_norm = 0.0, difference_norm = 0.0;
    PetscScalar reference[9], values[9], mutant[9], legacy[9];
    const PetscInt target_i = 1, target_j = 2, target_k = 3;
    MatStencil target[3] = {
        {.i = 1, .j = 2, .k = 3, .c = 0},
        {.i = 1, .j = 2, .k = 3, .c = 1},
        {.i = 1, .j = 2, .k = 3, .c = 2}
    };
    MatStencil conditioned = {.i = 0, .j = 2, .k = 3, .c = 0};
    MatStencil homogeneous = {.i = 0, .j = 2, .k = 3, .c = 1};
    Mat saved_matrix = NULL, mffd = NULL;
    SNES snes = NULL;
    Vec direction = NULL, product = NULL, px = NULL;
    PetscScalar px_sum = 0.0;
    PetscReal px_norm = 0.0;
    KSP ksp = NULL;
    PC pc = NULL;
    Vec pc_rhs = NULL, pc_solution = NULL;
    MatInfo initial_allocation_info, repeated_allocation_info;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &f));
    PetscCall(VecCopy(user->Ucont, x));
    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f, &ctx));
    PetscCall(MomentumPreconditionerEngine_Create(user, NULL, &description, &engine));
    preconditioning_matrix = engine.preconditioning_matrix;
    PetscCall(PicurvAssertBool((PetscBool)(engine.model_ops == &frozen_momentum_point_block_ops),
                               "engine must select the frozen-momentum model callbacks"));
    PetscCall(PicurvAssertBool(engine.owns_preconditioning_matrix,
                               "point-block engine must own its separate matrix"));
    PetscCall(PicurvAssertBool((PetscBool)!engine.aliases_jacobian_operator,
                               "point-block matrix must not alias the Jacobian operator"));
    PetscCall(MomentumPreconditionerEngine_Assemble(&engine, user, x));
    PetscCall(AssertExactPointBlockMatrixAllocation(
        user, preconditioning_matrix, PETSC_FALSE));
    PetscCall(MatGetInfo(preconditioning_matrix, MAT_GLOBAL_SUM,
                         &initial_allocation_info));
    PetscCall(MatNorm(preconditioning_matrix, NORM_FROBENIUS, &matrix_norm));
    PetscCall(PicurvAssertBool((PetscBool)(matrix_norm > 0.0),
                               "model callback and common rows must insert matrix entries"));
    PetscCall(MatGetBlockSize(preconditioning_matrix, &block_size));
    PetscCall(PicurvAssertIntEqual(3, block_size, "point-block matrix block size"));
    PetscCall(DMDAGetInfo(user->fda, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           &velocity_dof, NULL, NULL, NULL, NULL, NULL));
    PetscCall(PicurvAssertIntEqual(3, velocity_dof, "Newton velocity DMDA dof"));

    /* Cartesian limit: independently evaluate and compare every entry. */
    PetscCall(CollectiveLegacyPointBlockOracle(user, target_i, target_j, target_k, 0, reference));
    for (PetscInt rr = 0; rr < 3; ++rr)
        for (PetscInt cc = 0; cc < 3; ++cc) {
            PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
                target[rr], target[cc], &values[3 * rr + cc]));
            PetscCall(PicurvAssertRealNear(PetscRealPart(reference[3 * rr + cc]),
                PetscRealPart(values[3 * rr + cc]), 1e-13,
                "Cartesian point block entry must match independent oracle"));
        }

    /* Both fixed categories are exact identity rows, with no same-cell coupling. */
    for (PetscInt cc = 0; cc < 3; ++cc) {
        MatStencil col = conditioned; col.c = cc;
        PetscScalar conditioned_value = 0.0, homogeneous_value = 0.0;
        PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
            conditioned, col, &conditioned_value));
        col = homogeneous; col.c = cc;
        PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
            homogeneous, col, &homogeneous_value));
        PetscCall(PicurvAssertRealNear(cc == conditioned.c ? 1.0 : 0.0,
            PetscRealPart(conditioned_value), 1e-14, "conditioned row must be exact identity"));
        PetscCall(PicurvAssertRealNear(cc == homogeneous.c ? 1.0 : 0.0,
            PetscRealPart(homogeneous_value), 1e-14, "homogeneous row must be exact identity"));
    }

    /* A real residual and MFFD product cannot change subsequent assembly. */
    PetscCall(MatDuplicate(preconditioning_matrix, MAT_COPY_VALUES, &saved_matrix));
    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetDM(snes, user->fda));
    PetscCall(SNESSetFunction(snes, f, MomentumNewtonKrylov_FormResidual, &ctx));
    PetscCall(MatCreateSNESMF(snes, &mffd));
    PetscCall(VecDuplicate(x, &direction)); PetscCall(VecDuplicate(x, &product));
    PetscCall(VecSet(direction, 0.375));
    PetscCall(MomentumNewtonKrylov_FormResidual(snes, x, f, &ctx));
    PetscCall(MatMFFDSetBase(mffd, x, f));
    PetscCall(MatMult(mffd, direction, product));
    PetscCall(MomentumPreconditionerEngine_Assemble(&engine, user, x));
    PetscCall(MatAXPY(saved_matrix, -1.0, preconditioning_matrix, SAME_NONZERO_PATTERN));
    PetscCall(MatNorm(saved_matrix, NORM_FROBENIUS, &difference_norm));
    PetscCall(PicurvAssertRealNear(0.0, difference_norm, 1e-12,
        "assembly must be unchanged after residual and MFFD products"));
    PetscCall(VecDestroy(&product)); PetscCall(VecDestroy(&direction));
    PetscCall(MatDestroy(&mffd)); PetscCall(SNESDestroy(&snes)); PetscCall(MatDestroy(&saved_matrix));

    /* Nonuniform oracle: i/j/k, every component, and all samples are distinct. */
    PetscCall(SeedPointBlockOracleFields(user));
    PetscCall(VecCopy(user->Ucont, x));
    simCtx->step = 1;
    PetscCall(MomentumPreconditionerEngine_Assemble(&engine, user, x));
    PetscCall(CollectiveLegacyPointBlockOracle(user, target_i, target_j, target_k, 0, reference));
    for (PetscInt rr = 0; rr < 3; ++rr)
        for (PetscInt cc = 0; cc < 3; ++cc) {
            PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
                target[rr], target[cc], &values[3 * rr + cc]));
            PetscCall(PicurvAssertRealNear(PetscRealPart(reference[3 * rr + cc]),
                PetscRealPart(values[3 * rr + cc]), 1e-13,
                "nonuniform point block entry must match independent oracle"));
        }
    PetscCall(PicurvAssertBool((PetscBool)(PetscAbsScalar(reference[1] - reference[3]) > 1e-6 &&
                                           PetscAbsScalar(reference[2] - reference[6]) > 1e-6 &&
                                           PetscAbsScalar(reference[5] - reference[7]) > 1e-6),
        "oracle must preserve nonsymmetric component ordering"));
    for (PetscInt mutant_flag = ORACLE_CENTER_AJ; mutant_flag <= ORACLE_LEGACY_SIGN;
         mutant_flag <<= 1) {
        PetscBool differs = PETSC_FALSE;
        PetscCall(CollectiveLegacyPointBlockOracle(user, target_i, target_j, target_k,
            mutant_flag, mutant));
        for (PetscInt n = 0; n < 9; ++n)
            if (PetscAbsScalar(reference[n] - mutant[n]) > 1e-8) differs = PETSC_TRUE;
        PetscCall(PicurvAssertBool(differs, "independent oracle must reject audited mutant"));
    }
    PetscCall(CollectiveLegacyPointBlockOracle(user, target_i, target_j, target_k,
        ORACLE_OMIT_A5, mutant));
    PetscCall(PicurvAssertBool((PetscBool)(PetscAbsScalar(reference[8] - mutant[8]) > 1e-8),
        "A[5] must contribute to the zeta diagonal"));
    PetscCall(CollectiveLegacyPointBlockOracle(user, target_i, target_j, target_k,
        ORACLE_LEGACY_SIGN, legacy));
    for (PetscInt n = 0; n < 9; ++n)
        PetscCall(PicurvAssertRealNear(PetscRealPart(reference[n]), -PetscRealPart(legacy[n]),
            1e-13, "modern block must be the negative legacy block"));

    /* A coordinate permutation cannot accidentally address the intended block. */
    PetscCall(CollectiveLegacyPointBlockOracle(user, target_k, target_j, target_i, 0, mutant));
    {
        PetscBool differs = PETSC_FALSE;
        for (PetscInt n = 0; n < 9; ++n)
            if (PetscAbsScalar(reference[n] - mutant[n]) > 1e-8) differs = PETSC_TRUE;
        PetscCall(PicurvAssertBool(differs, "permuted MatStencil coordinates must be detectable"));
    }
    {
        MatStencil neighbor = target[0];
        PetscScalar neighbor_value = 0.0;
        neighbor.i++;
        PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
            target[0], neighbor, &neighbor_value));
        PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(neighbor_value), 1e-14,
            "point block must not insert unintended neighbor entries"));
    }

    /* The shared time coefficient supplies BDF1 and BDF2 diagonals only. */
    simCtx->step = 2;
    PetscCall(MomentumPreconditionerEngine_Assemble(&engine, user, x));
    PetscCall(CollectiveLegacyPointBlockOracle(user, target_i, target_j, target_k, 0, mutant));
    for (PetscInt n = 0; n < 9; ++n) {
        PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
            target[n / 3], target[n % 3], &values[n]));
        PetscCall(PicurvAssertRealNear(PetscRealPart(mutant[n]), PetscRealPart(values[n]),
            1e-13, "BDF2 point block must match independent oracle"));
        if (n == 0 || n == 4 || n == 8)
            PetscCall(PicurvAssertRealNear(0.5 / simCtx->dt,
                PetscRealPart(mutant[n] - reference[n]), 1e-10, "BDF2 diagonal increment"));
        else
            PetscCall(PicurvAssertRealNear(0.0, PetscRealPart(mutant[n] - reference[n]),
                1e-13, "BDF order must not change off-diagonal entries"));
    }

    PetscCall(VecDuplicate(x, &px));
    PetscCall(MatMult(preconditioning_matrix, x, px));
    PetscCall(VecSum(px, &px_sum)); PetscCall(VecNorm(px, NORM_2, &px_norm));
    PetscCall(PicurvAssertRealNear(8.40570474622236e4, PetscRealPart(px_sum), 5e-9,
        "P*x global sum must be decomposition independent"));
    PetscCall(PicurvAssertRealNear(9.05678688399519e3, px_norm, 5e-10,
        "P*x norm must be decomposition independent"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "POINT_BLOCK_MPI_SIGNATURE sum=%.16e norm2=%.16e\n",
        (double)PetscRealPart(px_sum), (double)px_norm));
    PetscCall(VecDestroy(&px));

    {
        PetscBool assembled = PETSC_FALSE;
        PetscCall(MatAssembled(preconditioning_matrix, &assembled));
        PetscCall(PicurvAssertBool(assembled, "engine must perform final matrix assembly"));
    }
    PetscCall(MatNorm(preconditioning_matrix, NORM_FROBENIUS, &matrix_norm));
    PetscCall(MatShift(preconditioning_matrix, 7.0));
    PetscCall(MomentumPreconditionerEngine_Assemble(&engine, user, x));
    PetscCall(MatGetInfo(preconditioning_matrix, MAT_GLOBAL_SUM,
                         &repeated_allocation_info));
    PetscCall(PicurvAssertRealNear((PetscReal)initial_allocation_info.nz_allocated,
        (PetscReal)repeated_allocation_info.nz_allocated, 0.0,
        "repeated assembly must retain exact allocated storage"));
    PetscCall(PicurvAssertRealNear((PetscReal)initial_allocation_info.mallocs,
        (PetscReal)repeated_allocation_info.mallocs, 0.0,
        "repeated assembly must not add insertion reallocations"));
    PetscCall(MatNorm(preconditioning_matrix, NORM_FROBENIUS, &reassembled_norm));
    PetscCall(PicurvAssertRealNear(matrix_norm, reassembled_norm, 1e-12,
                                   "repeated engine assembly must clear old entries"));
    {
        PetscScalar value = 0.0;
        PetscCall(PreconditionerMatrixStencilEntry(user, preconditioning_matrix,
            conditioned, conditioned, &value));
        PetscCall(PicurvAssertRealNear(1.0, PetscRealPart(value), 1e-14,
            "repeated engine assembly must clear old entries"));
    }
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, preconditioning_matrix, preconditioning_matrix));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(MomentumPreconditionerEngine_ConfigurePetscPC(&engine, pc));
    PetscCall(MomentumPreconditionerEngine_ValidatePetscPC(&engine, pc));
    PetscCall(KSPSetUp(ksp));
    PetscCall(DMCreateGlobalVector(user->fda, &pc_rhs));
    PetscCall(DMCreateGlobalVector(user->fda, &pc_solution));
    PetscCall(VecSet(pc_rhs, 1.0));
    PetscCall(PCApply(pc, pc_rhs, pc_solution));
    PetscCall(VecDestroy(&pc_solution)); PetscCall(VecDestroy(&pc_rhs));
    PetscCall(KSPDestroy(&ksp));
    PetscCall(MomentumPreconditionerEngine_Destroy(&engine));
    PetscCall(PicurvAssertBool((PetscBool)(engine.preconditioning_matrix == NULL),
                               "engine destroy must clear its owned matrix"));
    PetscCall(VecDestroy(&f)); PetscCall(VecDestroy(&x)); PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Proves that one periodic matrix row has only its exact +1/-1 pair. */
static PetscErrorCode AssertPeriodicPreconditionerRow(UserCtx *user, Mat matrix,
    MatStencil row, const char *message)
{
    AO ao = NULL;
    PetscInt ri, rj, rk, global_row, global_rep, lo, hi, local_checked = 0, checked = 0;
    PetscInt ncols = 0;
    const PetscInt *cols = NULL;
    const PetscScalar *values = NULL;
    MomentumNewtonKrylovRowType type;

    PetscFunctionBeginUser;
    type = MomentumNewtonKrylov_ClassifyRow(user, row.i, row.j, row.k, row.c, &ri, &rj, &rk);
    PetscCall(PicurvAssertIntEqual(MOM_NK_ROW_PERIODIC_DUPLICATE, type, message));
    global_row = row.c + 3 * (row.i + user->info.mx * (row.j + user->info.my * row.k));
    PetscCall(DMDAGetAO(user->fda, &ao));
    PetscCall(AOApplicationToPetsc(ao, 1, &global_row));
    PetscCall(MatGetOwnershipRange(matrix, &lo, &hi));
    if (global_row >= lo && global_row < hi) {
        PetscBool found_self = PETSC_FALSE, found_rep = PETSC_FALSE;
        PetscInt nonzero_entries = 0;
        PetscCall(TestStencilToGlobal(user,
            (MatStencil){.i = ri, .j = rj, .k = rk, .c = row.c}, &global_rep));
        PetscCall(MatGetRow(matrix, global_row, &ncols, &cols, &values));
        for (PetscInt n = 0; n < ncols; ++n) {
            if (PetscAbsScalar(values[n]) <= 1e-14) continue;
            nonzero_entries++;
            if (cols[n] == global_row) {
                found_self = PETSC_TRUE;
                PetscCall(PicurvAssertRealNear(1.0, PetscRealPart(values[n]), 1e-14,
                    "periodic row self entry"));
            } else if (cols[n] == global_rep) {
                found_rep = PETSC_TRUE;
                PetscCall(PicurvAssertRealNear(-1.0, PetscRealPart(values[n]), 1e-14,
                    "periodic row representative entry"));
            } else {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
                        "Periodic row contains an unintended nonzero column.");
            }
        }
        PetscCall(PicurvAssertIntEqual(2, nonzero_entries,
            "periodic row must contain exactly two numerical entries"));
        PetscCall(PicurvAssertBool((PetscBool)(found_self && found_rep), message));
        PetscCall(MatRestoreRow(matrix, global_row, &ncols, &cols, &values));
        local_checked = 1;
    }
    PetscCallMPI(MPI_Allreduce(&local_checked, &checked, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(1, checked, "periodic row must have exactly one matrix owner"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Verifies Jacobian creation/registration and baseline alias ownership. */
static PetscErrorCode TestJacobianInterfaceAndBaselineAlias(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    MomentumNewtonKrylovContext ctx = {0};
    MomentumPreconditionerDescription description = {
        MOM_NK_PC_MODEL_NONE, MOM_NK_PC_STRUCTURE_NONE, 0, 0
    };
    SNES snes = NULL;
    Mat jacobian_operator = NULL, preconditioning_matrix = NULL;
    KSP ksp = NULL;
    PC pc = NULL;
    PetscInt rows = 0, cols = 0;
    const char *jacobian_prefix = NULL;
    const char *pc_type = NULL;
    PetscBool prefix_matches = PETSC_FALSE;
    PetscBool pc_is_none = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(fixed_wall_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
    PetscCall(SNESSetDM(snes, user->fda));
    ctx.user = user;
    ctx.jacobian.type = MOM_NK_JACOBIAN_FINITE_DIFFERENCE;
    ctx.jacobian.finite_difference_mode = MOM_NK_FD_MODE_MATRIX_FREE;
    PetscCall(SNESSetFunction(snes, NULL, MomentumNewtonKrylov_FormResidual, &ctx));
    PetscCall(MomentumNewtonJacobian_Create(snes, &ctx.jacobian));
    PetscCall(MatGetOptionsPrefix(ctx.jacobian.jacobian_operator, &jacobian_prefix));
    PetscCall(PetscStrcmp(jacobian_prefix, "mom_nk_", &prefix_matches));
    PetscCall(PicurvAssertBool(prefix_matches,
                               "Jacobian interface must apply the application prefix"));
    PetscCall(MomentumPreconditionerEngine_Create(
        user, ctx.jacobian.jacobian_operator, &description, &ctx.preconditioning_engine));
    PetscCall(MomentumNewtonJacobian_Register(
        snes, &ctx.jacobian, &ctx.preconditioning_engine, &ctx));
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(MomentumPreconditionerEngine_ConfigurePetscPC(
        &ctx.preconditioning_engine, pc));
    PetscCall(MomentumPreconditionerEngine_ValidatePetscPC(
        &ctx.preconditioning_engine, pc));
    PetscCall(PCGetType(pc, &pc_type));
    PetscCall(PetscStrcmp(pc_type, PCNONE, &pc_is_none));
    PetscCall(PicurvAssertBool(pc_is_none,
                               "baseline engine must derive PETSc PCNONE"));
    PetscCall(SNESGetJacobian(snes, &jacobian_operator, &preconditioning_matrix, NULL, NULL));
    PetscCall(PicurvAssertBool(
        (PetscBool)(jacobian_operator == ctx.jacobian.jacobian_operator),
        "SNES must receive the Jacobian-interface MFFD operator"));
    PetscCall(PicurvAssertBool((PetscBool)(preconditioning_matrix == jacobian_operator),
                               "baseline preconditioning matrix must alias the Jacobian"));
    PetscCall(PicurvAssertBool(ctx.preconditioning_engine.aliases_jacobian_operator,
                               "baseline engine must record matrix aliasing"));
    PetscCall(PicurvAssertBool((PetscBool)!ctx.preconditioning_engine.owns_preconditioning_matrix,
                               "baseline engine must not own the Jacobian alias"));
    PetscCall(MomentumPreconditionerEngine_Destroy(&ctx.preconditioning_engine));
    PetscCall(MatGetSize(ctx.jacobian.jacobian_operator, &rows, &cols));
    PetscCall(PicurvAssertBool((PetscBool)(rows > 0 && cols > 0),
                               "destroying an alias engine must preserve the Jacobian"));
    PetscCall(MomentumNewtonJacobian_Destroy(&ctx.jacobian));
    PetscCall(SNESDestroy(&snes));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(DestroyNewtonFixture(&simCtx, tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Exercises engine-owned periodic duplicate rows on every MPI layout. */
static PetscErrorCode TestPointBlockPeriodicAssembly(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    MomentumNewtonKrylovContext ctx = {0};
    MomentumPreconditionerDescription description = {
        MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN,
        MOM_NK_PC_STRUCTURE_POINT_BLOCK,
        0,
        0
    };
    Vec x = NULL, f = NULL;
    PetscReal matrix_norm = 0.0;

    PetscFunctionBeginUser;
    PetscCall(BuildNewtonFixture(periodic_xyz_bcs, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &x));
    PetscCall(VecDuplicate(user->Ucont, &f));
    PetscCall(VecCopy(user->Ucont, x));
    ctx.user = user;
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, x, f, &ctx));
    PetscCall(MomentumPreconditionerEngine_Create(
        user, NULL, &description, &ctx.preconditioning_engine));
    PetscCall(MomentumPreconditionerEngine_Assemble(&ctx.preconditioning_engine, user, x));
    PetscCall(AssertExactPointBlockMatrixAllocation(
        user, ctx.preconditioning_engine.preconditioning_matrix, PETSC_TRUE));
    PetscCall(MatNorm(ctx.preconditioning_engine.preconditioning_matrix,
                      NORM_FROBENIUS, &matrix_norm));
    PetscCall(PicurvAssertBool((PetscBool)(matrix_norm > 0.0),
                               "periodic point-block assembly must produce a nonzero matrix"));
    PetscCall(AssertPeriodicPreconditionerRow(user,
        ctx.preconditioning_engine.preconditioning_matrix,
        (MatStencil){.i = 0, .j = 2, .k = 3, .c = 0},
        "single-axis periodic row must contain exact +1/-1 entries"));
    PetscCall(AssertPeriodicPreconditionerRow(user,
        ctx.preconditioning_engine.preconditioning_matrix,
        (MatStencil){.i = 0, .j = 0, .k = 3, .c = 1},
        "periodic intersection must contain exact +1/-1 entries"));
    PetscCall(AssertPeriodicPreconditionerRow(user,
        ctx.preconditioning_engine.preconditioning_matrix,
        (MatStencil){.i = 0, .j = 0, .k = 0, .c = 2},
        "periodic origin intersection must contain exact +1/-1 entries"));
    PetscCall(MomentumPreconditionerEngine_Destroy(&ctx.preconditioning_engine));
    PetscCall(VecDestroy(&f));
    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&user->Rhs));
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
        {"flat-channel-bdf1-startup", TestFlatChannelStartup},
        {"restart-and-continuation-solve", TestRestartAndContinuationSolve},
        {"small-solve-and-rollback", TestSmallSolveAndRollback},
        {"unsupported-configuration-fails-before-allocation", TestUnsupportedConfigurationFailsBeforeAllocation},
        {"post-allocation-failure-cleanup", TestPostAllocationFailureCleanup},
        {"linearization-config-parsing", TestLinearizationConfigParsing},
        {"point-block-preconditioner-engine", TestPointBlockPreconditionerEngine},
        {"jacobian-interface-and-baseline-alias", TestJacobianInterfaceAndBaselineAlias},
        {"point-block-periodic-assembly", TestPointBlockPeriodicAssembly},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv Newton Krylov tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-newton-krylov", cases, sizeof(cases) / sizeof(cases[0]));
    if (PetscFinalize()) return 1;
    return (int)ierr;
}
