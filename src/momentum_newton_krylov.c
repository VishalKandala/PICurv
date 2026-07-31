#include "momentumsolvers.h"

typedef enum {
    MOM_NK_JACOBIAN_FINITE_DIFFERENCE
} MomentumNewtonJacobianType;

typedef enum {
    MOM_NK_FD_MODE_MATRIX_FREE
} MomentumNewtonFiniteDifferenceMode;

typedef struct {
    MomentumNewtonJacobianType type;
    MomentumNewtonFiniteDifferenceMode finite_difference_mode;
    Mat jacobian_operator;
} MomentumNewtonJacobian;

typedef enum {
    MOM_NK_PC_MODEL_NONE,
    MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN
} MomentumPreconditionerModel;

typedef enum {
    MOM_NK_PC_STRUCTURE_NONE,
    MOM_NK_PC_STRUCTURE_POINT_BLOCK
} MomentumPreconditionerStructure;

typedef struct {
    MomentumPreconditionerModel model;
    MomentumPreconditionerStructure structure;
    PetscInt block_size;
    PetscInt stencil_width;
} MomentumPreconditionerDescription;

typedef struct {
    PetscErrorCode (*Describe)(UserCtx *, MomentumPreconditionerDescription *);
    PetscErrorCode (*AssembleInterior)(UserCtx *, Vec, Mat);
} MomentumPreconditionerModelOps;

typedef struct {
    MomentumPreconditionerDescription description;
    const MomentumPreconditionerModelOps *model_ops;
    Mat preconditioning_matrix;
    PetscBool aliases_jacobian_operator;
    PetscBool owns_preconditioning_matrix;
    const char *petsc_pc_type;
} MomentumPreconditionerEngine;

typedef struct {
    UserCtx *user;
    FILE *history_file;
    PetscBool have_initial_norm;
    PetscReal initial_norm;
    /* These objects are owned by the solve context.  They are deliberately
     * kept here (rather than in UserCtx) because an SNES is per physical solve. */
    MomentumNewtonJacobian jacobian;
    MomentumPreconditionerEngine preconditioning_engine;
} MomentumNewtonKrylovContext;

typedef enum {
    MOM_NK_ROW_PHYSICAL = 0,
    MOM_NK_ROW_FIXED_CONDITIONED,
    MOM_NK_ROW_FIXED_HOMOGENEOUS,
    MOM_NK_ROW_PERIODIC_DUPLICATE
} MomentumNewtonKrylovRowType;

static PetscErrorCode MomentumNewtonKrylov_Validate(UserCtx *user);
static PetscErrorCode MomentumNewtonKrylov_FormResidual(SNES snes, Vec X, Vec F, void *ctx);
static PetscErrorCode MomentumNewtonKrylov_Monitor(SNES snes, PetscInt iteration,
                                                    PetscReal norm, void *ctx);
static void MomentumNewtonKrylov_OpenHistory(MomentumNewtonKrylovContext *ctx);
static void MomentumNewtonKrylov_WriteSummary(const MomentumNewtonKrylovContext *ctx,
                                               SNESConvergedReason reason,
                                               PetscInt nonlinear_its,
                                               PetscInt function_evals,
                                               PetscInt linear_its,
                                               PetscReal final_norm,
                                               PetscBool committed);
static PetscErrorCode MomentumNewtonKrylov_ApplyConstraints(MomentumNewtonKrylovContext *ctx,
                                                            Vec X, Vec F);
static MomentumNewtonKrylovRowType MomentumNewtonKrylov_ClassifyRow(
    UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscInt component,
    PetscInt *ri, PetscInt *rj, PetscInt *rk);
static PetscErrorCode MomentumNewtonKrylov_ReadLinearizationConfig(
    MomentumNewtonJacobian *jacobian, MomentumPreconditionerDescription *description);
static PetscErrorCode MomentumNewtonKrylov_FormJacobian(SNES snes, Vec current_solution,
    Mat jacobian_operator, Mat preconditioning_matrix, void *vctx);
static PetscErrorCode FrozenMomentumJacobian_AssemblePointBlocks(
    UserCtx *user, Vec current_solution, Mat preconditioning_matrix);

/**
 * @brief Captures SNES iteration norms and optionally writes PICurv history rows.
 * @details SNES supplies the already-computed norm, so this monitor never causes
 * an additional nonlinear residual evaluation. PETSc monitors selected through
 * `-mom_nk_snes_monitor` remain independent and may run alongside this callback.
 */
static PetscErrorCode MomentumNewtonKrylov_Monitor(SNES snes, PetscInt iteration,
                                                    PetscReal norm, void *vctx)
{
    MomentumNewtonKrylovContext *ctx = (MomentumNewtonKrylovContext *)vctx;
    SimCtx *simCtx = ctx->user->simCtx;

    (void)snes;
    PetscFunctionBeginUser;
    if (iteration == 0 && !ctx->have_initial_norm) {
        ctx->initial_norm = norm;
        ctx->have_initial_norm = PETSC_TRUE;
    }
    if (ctx->history_file) {
        (void)fprintf(ctx->history_file,
                      "step: %d | block: %d | newton: %d | nonlinear_norm: %.16e\n",
                      (int)simCtx->step, (int)ctx->user->_this, (int)iteration,
                      (double)norm);
        (void)fflush(ctx->history_file);
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Opens the optional rank-zero Newton iteration-history file. */
static void MomentumNewtonKrylov_OpenHistory(MomentumNewtonKrylovContext *ctx)
{
    SimCtx *simCtx = ctx->user->simCtx;
    char path[PETSC_MAX_PATH_LEN + 128];
    const char *mode;

    if (!simCtx->mom_nk_monitor_history || simCtx->rank != 0) return;
    if (PetscSNPrintf(path, sizeof(path),
                     "%s/Momentum_Solver_Newton_Krylov_History_Block_%d.log",
                     simCtx->log_dir, (int)ctx->user->_this)) return;
    mode = (simCtx->step == simCtx->StartStep + 1 && !simCtx->continueMode) ? "w" : "a";
    ctx->history_file = fopen(path, mode);
    if (!ctx->history_file) {
        LOG(GLOBAL, LOG_WARNING, "Could not open Newton iteration-history log '%s'.\n", path);
        return;
    }
    if (mode[0] == 'w') {
        (void)fprintf(ctx->history_file,
                      "# step | block | Newton iteration | nonlinear residual norm\n");
    } else if (simCtx->continueMode && simCtx->step == simCtx->StartStep + 1) {
        (void)fprintf(ctx->history_file, "# Continuation from step %d\n", (int)simCtx->StartStep);
    }
}

/**
 * @brief Appends one rank-zero structured Newton result for a physical step.
 * @details File failures are deliberately diagnostic-only: rollback and PETSc
 * cleanup must retain their original error behavior.
 */
static void MomentumNewtonKrylov_WriteSummary(const MomentumNewtonKrylovContext *ctx,
                                               SNESConvergedReason reason,
                                               PetscInt nonlinear_its,
                                               PetscInt function_evals,
                                               PetscInt linear_its,
                                               PetscReal final_norm,
                                               PetscBool committed)
{
    SimCtx *simCtx = ctx->user->simCtx;
    char path[PETSC_MAX_PATH_LEN + 128];
    const char *mode;
    const char *reason_name;
    FILE *file;

    if (simCtx->rank != 0) return;
    if (PetscSNPrintf(path, sizeof(path),
                     "%s/Momentum_Solver_Newton_Krylov_Summary_Block_%d.log",
                     simCtx->log_dir, (int)ctx->user->_this)) return;
    mode = (simCtx->step == simCtx->StartStep + 1 && !simCtx->continueMode) ? "w" : "a";
    file = fopen(path, mode);
    if (!file) {
        LOG(GLOBAL, LOG_WARNING, "Could not open Newton summary log '%s'.\n", path);
        return;
    }
    if (mode[0] == 'w') {
        (void)fprintf(file,
                      "# step | block | solver | Jacobian | preconditioner | SNES reason | "
                      "reason code | Newton iterations | "
                      "residual evaluations | Krylov iterations | initial nonlinear norm | "
                      "final nonlinear norm | state\n");
    } else if (simCtx->continueMode && simCtx->step == simCtx->StartStep + 1) {
        (void)fprintf(file, "# Continuation from step %d\n", (int)simCtx->StartStep);
    }
    reason_name = reason == SNES_CONVERGED_ITERATING
                    ? "SNES_CONVERGED_ITERATING" : SNESConvergedReasons[reason];
    (void)fprintf(file,
                  "step: %d | block: %d | solver: Newton Krylov | "
                  "Jacobian: finite_difference / matrix_free | Preconditioner: %s | "
                  "reason: %s | reason_code: %d | "
                  "newton: %d | evals: %d | krylov: %d | initial: ",
                  (int)simCtx->step, (int)ctx->user->_this,
                  ctx->preconditioning_engine.description.model == MOM_NK_PC_MODEL_NONE ?
                      "none" : "frozen_momentum_jacobian / point_block",
                  reason_name, (int)reason,
                  (int)nonlinear_its, (int)function_evals, (int)linear_its);
    if (ctx->have_initial_norm) (void)fprintf(file, "%.16e", (double)ctx->initial_norm);
    else                        (void)fprintf(file, "unavailable");
    (void)fprintf(file, " | final: %.16e | state: %s\n", (double)final_norm,
                  committed ? "committed" : "rolled_back");
    (void)fclose(file);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonKrylov_Validate"
/**
 * @brief Rejects configurations outside the audited version-one feature set.
 * @param user Single-block momentum context to validate.
 * @return PetscErrorCode 0 when the configuration is supported.
 */
static PetscErrorCode MomentumNewtonKrylov_Validate(UserCtx *user)
{
    SimCtx   *simCtx;
    PetscReal mask_max = 0.0;
    PetscInt  velocity_dof = 0;

    PetscFunctionBeginUser;
    PetscCheck(user != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Newton Krylov requires a non-NULL UserCtx.");
    simCtx = user->simCtx;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Newton Krylov requires UserCtx::simCtx.");
    PetscCall(DMDAGetInfo(user->fda, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           &velocity_dof, NULL, NULL, NULL, NULL, NULL));
    PetscCheck(velocity_dof == 3, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Newton Krylov requires a three-component velocity DMDA (got dof=%d).",
               velocity_dof);
    PetscCheck(simCtx->block_number == 1, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one supports exactly one block (got %d).",
               simCtx->block_number);
    PetscCheck(!simCtx->immersed, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support immersed boundaries.");
    PetscCheck(!simCtx->movefsi && !simCtx->rotatefsi, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support moving or rotating bodies/FSI.");
    PetscCheck(!simCtx->moveframe && !simCtx->rotateframe, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support moving or rotating reference frames.");
    PetscCheck(!simCtx->rans, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support RANS.");
    PetscCheck(!simCtx->clark, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support the Clark model.");
    PetscCheck(!simCtx->TwoD, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support TwoD component masking.");
    PetscCheck(!simCtx->wallfunction, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not support wall functions.");
    PetscCheck(simCtx->StartStep == 0, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one supports fresh starts only (StartStep must be zero)." );
    PetscCheck(!simCtx->continueMode, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one supports fresh starts only (--continue is not supported)." );

    for (PetscInt face = 0; face < 6; ++face) {
        const BoundaryFaceConfig *cfg = &user->boundary_faces[face];
        PetscBool supported = PETSC_FALSE;

        switch (cfg->handler_type) {
            case BC_HANDLER_WALL_NOSLIP:
                supported = (PetscBool)(cfg->mathematical_type == WALL);
                break;
            case BC_HANDLER_INLET_CONSTANT_VELOCITY:
            case BC_HANDLER_INLET_PARABOLIC:
            case BC_HANDLER_INLET_PROFILE_FROM_FILE:
                supported = (PetscBool)(cfg->mathematical_type == INLET);
                break;
            case BC_HANDLER_OUTLET_CONSERVATION:
                supported = (PetscBool)(cfg->mathematical_type == OUTLET);
                break;
            case BC_HANDLER_PERIODIC_GEOMETRIC:
                supported = (PetscBool)(cfg->mathematical_type == PERIODIC);
                break;
            default:
                supported = PETSC_FALSE;
                break;
        }
        PetscCheck(supported, PETSC_COMM_WORLD, PETSC_ERR_SUP,
                   "Newton Krylov version one does not support boundary face %d with mathematical type %d and handler %d.",
                   face, (PetscInt)cfg->mathematical_type, (PetscInt)cfg->handler_type);
    }

    PetscCheck((user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC) ==
               (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC),
               PETSC_COMM_WORLD, PETSC_ERR_SUP, "Newton Krylov requires paired x-periodic faces.");
    PetscCheck((user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC) ==
               (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC),
               PETSC_COMM_WORLD, PETSC_ERR_SUP, "Newton Krylov requires paired y-periodic faces.");
    PetscCheck((user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC) ==
               (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC),
               PETSC_COMM_WORLD, PETSC_ERR_SUP, "Newton Krylov requires paired z-periodic faces.");

    PetscCheck(user->Nvert != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Newton Krylov requires the cell-mask vector Nvert.");
    PetscCall(VecMax(user->Nvert, NULL, &mask_max));
    PetscCheck(mask_max <= 0.1, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not define equations for masked solid cells (max Nvert=%g).",
               (double)mask_max);
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonKrylov_ReadLinearizationConfig"
/**
 * @brief Reads application-owned Jacobian and preconditioner mathematics.
 */
static PetscErrorCode MomentumNewtonKrylov_ReadLinearizationConfig(
    MomentumNewtonJacobian *jacobian, MomentumPreconditionerDescription *description)
{
    char type[48] = "finite_difference";
    char finite_difference_mode[32] = "matrix_free";
    char model[48] = "none";
    char structure[32] = "none";
    PetscBool set = PETSC_FALSE, match = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCheck(jacobian != NULL && description != NULL, PETSC_COMM_SELF,
               PETSC_ERR_ARG_NULL, "Newton Krylov linearization configuration is NULL.");
    PetscCall(PetscOptionsGetString(NULL, NULL, "-mom_nk_jacobian_type",
                                    type, sizeof(type), &set));
    PetscCall(PetscStrcasecmp(type, "finite_difference", &match));
    PetscCheck(match, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "-mom_nk_jacobian_type must be 'finite_difference' (got '%s').", type);
    jacobian->type = MOM_NK_JACOBIAN_FINITE_DIFFERENCE;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-mom_nk_jacobian_fd_mode",
                                    finite_difference_mode,
                                    sizeof(finite_difference_mode), &set));
    PetscCall(PetscStrcasecmp(finite_difference_mode, "matrix_free", &match));
    PetscCheck(match, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "-mom_nk_jacobian_fd_mode must be 'matrix_free' (got '%s').",
               finite_difference_mode);
    jacobian->finite_difference_mode = MOM_NK_FD_MODE_MATRIX_FREE;

    PetscCall(PetscOptionsGetString(NULL, NULL, "-mom_nk_preconditioner_model",
                                    model, sizeof(model), &set));
    PetscCall(PetscStrcasecmp(model, "none", &match));
    if (match) description->model = MOM_NK_PC_MODEL_NONE;
    if (!match) {
        PetscCall(PetscStrcasecmp(model, "frozen_momentum_jacobian", &match));
        if (match) description->model = MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN;
    }
    PetscCheck(match, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "-mom_nk_preconditioner_model must be 'none' or "
               "'frozen_momentum_jacobian' (got '%s').", model);
    PetscCall(PetscOptionsGetString(NULL, NULL, "-mom_nk_preconditioner_structure",
                                    structure, sizeof(structure), &set));
    PetscCall(PetscStrcasecmp(structure, "none", &match));
    if (match) description->structure = MOM_NK_PC_STRUCTURE_NONE;
    if (!match) {
        PetscCall(PetscStrcasecmp(structure, "point_block", &match));
        if (match) description->structure = MOM_NK_PC_STRUCTURE_POINT_BLOCK;
    }
    PetscCheck(match, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "-mom_nk_preconditioner_structure must be 'none' or 'point_block' "
               "(got '%s').", structure);
    PetscCheck((description->model == MOM_NK_PC_MODEL_NONE &&
                description->structure == MOM_NK_PC_STRUCTURE_NONE) ||
               (description->model == MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN &&
                description->structure == MOM_NK_PC_STRUCTURE_POINT_BLOCK),
               PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Unsupported Newton Krylov preconditioner model/structure combination: "
               "model='%s', structure='%s'.", model, structure);
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Returns the preserved provisional frozen-momentum point block.
 *
 * The arithmetic in this function is carried forward unchanged for a separate
 * numerical audit. `A` retains its six existing directional face groups.
 */
static void FrozenMomentumJacobian_PointBlock(const SimCtx *simCtx,
    const Cmpnts ***ucont, const Cmpnts ***csi, const Cmpnts ***eta,
    const Cmpnts ***zet, const PetscReal ***aj, PetscInt i, PetscInt j,
    PetscInt k, PetscScalar block[9])
{
    const PetscReal dtc = MomentumBDFCoefficient((SimCtx *)simCtx) / simCtx->dt;
    const PetscReal nu = simCtx->ren > 0.0 ? 1.0 / simCtx->ren : 0.0;
    const PetscReal AJip = aj[k][j][i];
    const PetscReal AJjp = aj[k][j][i];
    const PetscReal AJkp = aj[k][j][i];
    const PetscReal g11ip = csi[k][j][i].x * csi[k][j][i].x +
                             csi[k][j][i].y * csi[k][j][i].y +
                             csi[k][j][i].z * csi[k][j][i].z;
    const PetscReal g22jp = eta[k][j][i].x * eta[k][j][i].x +
                             eta[k][j][i].y * eta[k][j][i].y +
                             eta[k][j][i].z * eta[k][j][i].z;
    const PetscReal g33kp = zet[k][j][i].x * zet[k][j][i].x +
                             zet[k][j][i].y * zet[k][j][i].y +
                             zet[k][j][i].z * zet[k][j][i].z;
    const PetscReal U1ip = ucont[k][j][i].y, U2ip = ucont[k][j][i].z;
    const PetscReal U0jp = ucont[k][j][i].x, U2jp = ucont[k][j][i].z;
    const PetscReal U0kp = ucont[k][j][i].x, U1kp = ucont[k][j][i].y;
    PetscReal A[6][4] = {{0.0}};
    PetscReal Su, Sv, Sw, nui, nuj, nuk;

    /* Preserved center-only bookkeeping: six directional face groups. */
    A[0][0] = 0.5 * AJip * ucont[k][j][i].x;
    A[1][0] = 0.5 * AJip * ucont[k][j][i - 1].x;
    A[2][1] = 0.5 * AJjp * ucont[k][j][i].y;
    A[3][1] = 0.5 * AJjp * ucont[k][j - 1][i].y;
    A[4][2] = 0.5 * AJkp * ucont[k][j][i].z;
    A[5][2] = 0.5 * AJkp * ucont[k - 1][j][i].z;
    Su = A[0][0] - A[1][0];
    Sv = A[2][1] - A[3][1];
    Sw = A[4][2] - A[5][2];
    nui = 2.0 * nu * AJip * g11ip;
    nuj = 2.0 * nu * AJjp * g22jp;
    nuk = 2.0 * nu * AJkp * g33kp;

    /* Preserve all existing signs and coefficient formulas verbatim. */
    block[0] = dtc + nui + Su; block[1] = 0.5 * AJip * U1ip; block[2] = 0.5 * AJip * U2ip;
    block[3] = 0.5 * AJjp * U0jp; block[4] = dtc + nuj + Sv; block[5] = 0.5 * AJjp * U2jp;
    block[6] = 0.5 * AJkp * U0kp; block[7] = 0.5 * AJkp * U1kp; block[8] = dtc + nuk + Sw;
}

#undef __FUNCT__
#define __FUNCT__ "FrozenMomentumJacobian_DescribePointBlock"
/** @brief Describes the provisional frozen-coefficient point-block model. */
static PetscErrorCode FrozenMomentumJacobian_DescribePointBlock(
    UserCtx *user, MomentumPreconditionerDescription *description)
{
    PetscFunctionBeginUser;
    (void)user;
    PetscCheck(description != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Preconditioner description is NULL.");
    description->model = MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN;
    description->structure = MOM_NK_PC_STRUCTURE_POINT_BLOCK;
    description->block_size = 3;
    description->stencil_width = 0;
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "FrozenMomentumJacobian_AssemblePointBlocks"
/** @brief Inserts only the existing interior frozen-momentum point blocks. */
static PetscErrorCode FrozenMomentumJacobian_AssemblePointBlocks(
    UserCtx *user, Vec current_solution, Mat preconditioning_matrix)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***ucont = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;
    PetscReal ***aj = NULL;
    PetscErrorCode ierr = PETSC_SUCCESS, cleanup_ierr;

    PetscFunctionBeginUser;
    (void)current_solution;
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont); if (ierr) goto cleanup;
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &csi); if (ierr) goto cleanup;
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &eta); if (ierr) goto cleanup;
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &zet); if (ierr) goto cleanup;
    ierr = DMDAVecGetArrayRead(user->da, user->lAj, &aj); if (ierr) goto cleanup;
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                for (PetscInt component = 0; component < 3; ++component) {
                    MatStencil row = {i, j, k, component};
                    PetscInt ri, rj, rk;
                    MomentumNewtonKrylovRowType type = MomentumNewtonKrylov_ClassifyRow(
                        user, i, j, k, component, &ri, &rj, &rk);
                    if (type == MOM_NK_ROW_PHYSICAL) {
                        PetscScalar block[9];
                        MatStencil cols[3] = {{i, j, k, 0}, {i, j, k, 1}, {i, j, k, 2}};
                        FrozenMomentumJacobian_PointBlock(user->simCtx, (const Cmpnts ***)ucont,
                            (const Cmpnts ***)csi, (const Cmpnts ***)eta, (const Cmpnts ***)zet,
                            (const PetscReal ***)aj, i, j, k, block);
                        ierr = MatSetValuesStencil(preconditioning_matrix, 1, &row, 3, cols,
                                                  &block[3 * component], INSERT_VALUES);
                        if (ierr) goto cleanup;
                    }
                }
            }
        }
    }
cleanup:
    if (aj) { cleanup_ierr = DMDAVecRestoreArrayRead(user->da, user->lAj, &aj); if (!ierr) ierr = cleanup_ierr; }
    if (zet) { cleanup_ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet); if (!ierr) ierr = cleanup_ierr; }
    if (eta) { cleanup_ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta); if (!ierr) ierr = cleanup_ierr; }
    if (csi) { cleanup_ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi); if (!ierr) ierr = cleanup_ierr; }
    if (ucont) { cleanup_ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont); if (!ierr) ierr = cleanup_ierr; }
    PetscFunctionReturn(ierr);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumPreconditionerEngine_ApplyConstraintRows"
/** @brief Inserts all common fixed, homogeneous, and periodic-duplicate rows. */
static PetscErrorCode MomentumPreconditionerEngine_ApplyConstraintRows(
    UserCtx *user, Mat preconditioning_matrix)
{
    DMDALocalInfo info = user->info;

    PetscFunctionBeginUser;
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                for (PetscInt component = 0; component < 3; ++component) {
                    PetscInt ri, rj, rk;
                    MomentumNewtonKrylovRowType type = MomentumNewtonKrylov_ClassifyRow(
                        user, i, j, k, component, &ri, &rj, &rk);
                    if (type != MOM_NK_ROW_PHYSICAL) {
                        const PetscScalar one = 1.0, minus_one = -1.0;
                        AO ao = NULL;
                        PetscInt global_row = 3 * ((k * info.my + j) * info.mx + i) + component;
                        PetscCall(DMDAGetAO(user->fda, &ao));
                        PetscCall(AOApplicationToPetsc(ao, 1, &global_row));
                        PetscCall(MatSetValue(preconditioning_matrix, global_row, global_row,
                                              one, INSERT_VALUES));
                        if (type == MOM_NK_ROW_PERIODIC_DUPLICATE) {
                            PetscInt global_col = 3 * ((rk * info.my + rj) * info.mx + ri) + component;
                            PetscCall(AOApplicationToPetsc(ao, 1, &global_col));
                            PetscCall(MatSetValue(preconditioning_matrix, global_row, global_col,
                                                  minus_one, INSERT_VALUES));
                        }
                    }
                }
            }
        }
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

static const MomentumPreconditionerModelOps frozen_momentum_point_block_ops = {
    FrozenMomentumJacobian_DescribePointBlock,
    FrozenMomentumJacobian_AssemblePointBlocks
};

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonJacobian_Create"
/** @brief Creates the selected Jacobian operator; currently PETSc MFFD only. */
static PetscErrorCode MomentumNewtonJacobian_Create(SNES snes,
    MomentumNewtonJacobian *jacobian)
{
    PetscFunctionBeginUser;
    PetscCheck(jacobian->type == MOM_NK_JACOBIAN_FINITE_DIFFERENCE &&
               jacobian->finite_difference_mode == MOM_NK_FD_MODE_MATRIX_FREE,
               PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Unsupported Newton Krylov Jacobian type/finite-difference-mode combination.");
    PetscCall(MatCreateSNESMF(snes, &jacobian->jacobian_operator));
    PetscCall(MatSetOptionsPrefix(jacobian->jacobian_operator, "mom_nk_"));
    PetscCall(PetscObjectSetName((PetscObject)jacobian->jacobian_operator,
                                 "momentum_jacobian_finite_difference_matrix_free"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Updates the matrix-free finite-difference operator base. */
static PetscErrorCode MomentumNewtonJacobian_Update(SNES snes, Vec current_solution,
    MomentumNewtonJacobian *jacobian)
{
    PetscFunctionBeginUser;
    PetscCall(MatMFFDComputeJacobian(snes, current_solution, jacobian->jacobian_operator,
                                     jacobian->jacobian_operator, NULL));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Registers the application orchestration callback and both SNES matrices. */
static PetscErrorCode MomentumNewtonJacobian_Register(SNES snes,
    MomentumNewtonJacobian *jacobian, MomentumPreconditionerEngine *engine,
    MomentumNewtonKrylovContext *ctx)
{
    PetscFunctionBeginUser;
    PetscCall(SNESSetJacobian(snes, jacobian->jacobian_operator,
                              engine->preconditioning_matrix,
                              MomentumNewtonKrylov_FormJacobian, ctx));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Destroys a partially or fully created Jacobian operator. */
static PetscErrorCode MomentumNewtonJacobian_Destroy(MomentumNewtonJacobian *jacobian)
{
    PetscFunctionBeginUser;
    PetscCall(MatDestroy(&jacobian->jacobian_operator));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumPreconditionerEngine_Create"
/** @brief Validates a model/structure and creates or aliases its matrix. */
static PetscErrorCode MomentumPreconditionerEngine_Create(UserCtx *user,
    Mat jacobian_operator, const MomentumPreconditionerDescription *requested,
    MomentumPreconditionerEngine *engine)
{
    PetscFunctionBeginUser;
    engine->description = *requested;
    if (requested->model == MOM_NK_PC_MODEL_NONE &&
        requested->structure == MOM_NK_PC_STRUCTURE_NONE) {
        engine->description.block_size = 0;
        engine->description.stencil_width = 0;
        engine->preconditioning_matrix = jacobian_operator;
        engine->aliases_jacobian_operator = PETSC_TRUE;
        engine->owns_preconditioning_matrix = PETSC_FALSE;
        engine->petsc_pc_type = PCNONE;
    } else if (requested->model == MOM_NK_PC_MODEL_FROZEN_MOMENTUM_JACOBIAN &&
               requested->structure == MOM_NK_PC_STRUCTURE_POINT_BLOCK) {
        engine->model_ops = &frozen_momentum_point_block_ops;
        PetscCall(engine->model_ops->Describe(user, &engine->description));
        PetscCall(DMCreateMatrix(user->fda, &engine->preconditioning_matrix));
        engine->owns_preconditioning_matrix = PETSC_TRUE;
        /* DMDA preallocates the point blocks; periodic duplicate equations may
         * add one wraparound column outside the local DM stencil. */
        PetscCall(MatSetOption(engine->preconditioning_matrix,
                               MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
        PetscCall(MatSetBlockSize(engine->preconditioning_matrix,
                                  engine->description.block_size));
        PetscCall(PetscObjectSetName((PetscObject)engine->preconditioning_matrix,
                                     "momentum_preconditioner_frozen_point_block"));
        engine->aliases_jacobian_operator = PETSC_FALSE;
        engine->petsc_pc_type = PCPBJACOBI;
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
                "Unsupported Newton Krylov preconditioner model/structure combination.");
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumPreconditionerEngine_Assemble"
/** @brief Runs model insertion, common row handling, and final assembly. */
static PetscErrorCode MomentumPreconditionerEngine_Assemble(
    MomentumPreconditionerEngine *engine, UserCtx *user, Vec current_solution)
{
    PetscFunctionBeginUser;
    if (engine->aliases_jacobian_operator) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCall(MatZeroEntries(engine->preconditioning_matrix));
    PetscCall(engine->model_ops->AssembleInterior(user, current_solution,
                                                   engine->preconditioning_matrix));
    PetscCall(MomentumPreconditionerEngine_ApplyConstraintRows(
        user, engine->preconditioning_matrix));
    PetscCall(MatAssemblyBegin(engine->preconditioning_matrix, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(engine->preconditioning_matrix, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Applies the validated model/structure-to-PETSc-PC mapping. */
static PetscErrorCode MomentumPreconditionerEngine_ConfigurePetscPC(
    MomentumPreconditionerEngine *engine, PC pc)
{
    PetscFunctionBeginUser;
    PetscCall(PCSetType(pc, engine->petsc_pc_type));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Rejects raw options that select an unvalidated PETSc PC backend. */
static PetscErrorCode MomentumPreconditionerEngine_ValidatePetscPC(
    MomentumPreconditionerEngine *engine, PC pc)
{
    const char *actual_type = NULL;
    PetscBool matches = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PCGetType(pc, &actual_type));
    PetscCall(PetscStrcmp(actual_type, engine->petsc_pc_type, &matches));
    PetscCheck(matches, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov preconditioner model/structure requires internal PETSc PC "
               "type '%s', but raw option processing selected '%s'.",
               engine->petsc_pc_type, actual_type ? actual_type : "(unset)");
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Destroys only a separately owned preconditioning matrix. */
static PetscErrorCode MomentumPreconditionerEngine_Destroy(
    MomentumPreconditionerEngine *engine)
{
    PetscFunctionBeginUser;
    if (!engine->owns_preconditioning_matrix) engine->preconditioning_matrix = NULL;
    PetscCall(MatDestroy(&engine->preconditioning_matrix));
    engine->aliases_jacobian_operator = PETSC_FALSE;
    engine->owns_preconditioning_matrix = PETSC_FALSE;
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonKrylov_FormJacobian"
/** @brief Updates the Jacobian and then assembles any separate preconditioning matrix. */
static PetscErrorCode MomentumNewtonKrylov_FormJacobian(SNES snes, Vec current_solution,
    Mat jacobian_operator, Mat preconditioning_matrix, void *vctx)
{
    MomentumNewtonKrylovContext *ctx = (MomentumNewtonKrylovContext *)vctx;
    PetscFunctionBeginUser;
    (void)jacobian_operator;
    (void)preconditioning_matrix;
    PetscCall(MomentumNewtonJacobian_Update(snes, current_solution, &ctx->jacobian));
    PetscCall(MomentumPreconditionerEngine_Assemble(
        &ctx->preconditioning_engine, ctx->user, current_solution));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonKrylov_ClassifyRow"
/**
 * @brief Classifies one stored staggered component row and its periodic representative.
 *
 * Negative nonperiodic planes and positive dummy planes are zeroed by the legacy
 * residual treatment. Only a face-normal row actually written by its boundary
 * handler is conditioned; the remaining zeroed rows use F=X. Periodic endpoint
 * rows use the same wrapped representatives as SynchronizePeriodicStaggeredFields().
 * At periodic/nonperiodic intersections the periodic equation owns rows that a
 * shrunken boundary-handler loop does not condition.
 *
 * @param user Block context containing boundary metadata and global dimensions.
 * @param i Global i index.
 * @param j Global j index.
 * @param k Global k index.
 * @param component Stored component, 0=x, 1=y, 2=z.
 * @param ri Returned periodic representative i index.
 * @param rj Returned periodic representative j index.
 * @param rk Returned periodic representative k index.
 * @return Exact version-one row category.
 */
static MomentumNewtonKrylovRowType MomentumNewtonKrylov_ClassifyRow(
    UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscInt component,
    PetscInt *ri, PetscInt *rj, PetscInt *rk)
{
    const PetscInt mx = user->info.mx, my = user->info.my, mz = user->info.mz;
    const PetscInt coord[3] = {i, j, k};
    const PetscInt size[3] = {mx, my, mz};
    const BCFace neg_face[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    PetscBool periodic[3], periodic_duplicate = PETSC_FALSE;
    PetscBool residual_zeroed = PETSC_FALSE, conditioned = PETSC_FALSE;

    *ri = i; *rj = j; *rk = k;
    for (PetscInt axis = 0; axis < 3; ++axis) {
        periodic[axis] = (PetscBool)(
            user->boundary_faces[neg_face[axis]].mathematical_type == PERIODIC);
        if (periodic[axis] && coord[axis] == 0) {
            periodic_duplicate = PETSC_TRUE;
            if (axis == 0) *ri = -2;
            else if (axis == 1) *rj = -2;
            else *rk = -2;
        }
        if (periodic[axis] && coord[axis] == size[axis] - 1) {
            periodic_duplicate = PETSC_TRUE;
            if (axis == 0) *ri = mx + 1;
            else if (axis == 1) *rj = my + 1;
            else *rk = mz + 1;
        }

        if (!periodic[axis] && coord[axis] == 0) residual_zeroed = PETSC_TRUE;
        if (coord[axis] == size[axis] - 1) residual_zeroed = PETSC_TRUE;
        if (!periodic[axis] && coord[axis] == size[axis] - 2 && component == axis)
            residual_zeroed = PETSC_TRUE;
    }

    if (!periodic[component] &&
        (coord[component] == 0 || coord[component] == size[component] - 2)) {
        PetscBool tangential_interior = PETSC_TRUE;
        for (PetscInt axis = 0; axis < 3; ++axis) {
            if (axis == component) continue;
            if (coord[axis] < 1 || coord[axis] > size[axis] - 2)
                tangential_interior = PETSC_FALSE;
        }
        conditioned = tangential_interior;
    }

    if (conditioned) return MOM_NK_ROW_FIXED_CONDITIONED;
    if (periodic_duplicate) return MOM_NK_ROW_PERIODIC_DUPLICATE;
    if (residual_zeroed) return MOM_NK_ROW_FIXED_HOMOGENEOUS;
    return MOM_NK_ROW_PHYSICAL;
}

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonKrylov_ApplyConstraints"
/**
 * @brief Replaces every non-independent residual row with an explicit equation.
 *
 * Conditioned face-normal rows use F=X-Uconditioned, unconditioned legacy
 * dummy/tangential rows use F=X, and periodic duplicates use Fdup=Xdup-Xrep.
 * These equations prevent the zero Jacobian rows produced by simply retaining
 * EnforceRHSBoundaryConditions() zeros in a matrix-free Newton operator. Immersed,
 * masked, TwoD, and interface rows are rejected before this callback is installed.
 * @param ctx Active solve context.
 * @param X Unconditioned PETSc trial state.
 * @param F Residual vector to update in place.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomentumNewtonKrylov_ApplyConstraints(MomentumNewtonKrylovContext *ctx,
                                                            Vec X, Vec F)
{
    UserCtx       *user = ctx->user;
    DMDALocalInfo  info = user->info;
    Vec            local_x = NULL;
    Cmpnts       ***x = NULL, ***conditioned = NULL, ***f = NULL, ***lx = NULL;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;

    PetscFunctionBeginUser;
    PetscCall(DMGetLocalVector(user->fda, &local_x));
    PetscCall(DMGlobalToLocalBegin(user->fda, X, INSERT_VALUES, local_x));
    PetscCall(DMGlobalToLocalEnd(user->fda, X, INSERT_VALUES, local_x));
    PetscCall(DMDAVecGetArrayRead(user->fda, X, &x));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &conditioned));
    PetscCall(DMDAVecGetArray(user->fda, F, &f));
    PetscCall(DMDAVecGetArrayRead(user->fda, local_x, &lx));

    for (PetscInt k = zs; k < ze; ++k) {
        for (PetscInt j = ys; j < ye; ++j) {
            for (PetscInt i = xs; i < xe; ++i) {
                PetscScalar *fv = &f[k][j][i].x;
                const PetscScalar *xv = &x[k][j][i].x;
                const PetscScalar *cv = &conditioned[k][j][i].x;

                for (PetscInt component = 0; component < 3; ++component) {
                    PetscInt ri, rj, rk;
                    MomentumNewtonKrylovRowType row = MomentumNewtonKrylov_ClassifyRow(
                        user, i, j, k, component, &ri, &rj, &rk);
                    const PetscScalar *rv = &lx[rk][rj][ri].x;

                    if (row == MOM_NK_ROW_FIXED_CONDITIONED) fv[component] = xv[component] - cv[component];
                    else if (row == MOM_NK_ROW_FIXED_HOMOGENEOUS) fv[component] = xv[component];
                    else if (row == MOM_NK_ROW_PERIODIC_DUPLICATE) fv[component] = xv[component] - rv[component];
                }
            }
        }
    }

    PetscCall(DMDAVecRestoreArrayRead(user->fda, local_x, &lx));
    PetscCall(DMDAVecRestoreArray(user->fda, F, &f));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &conditioned));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, X, &x));
    PetscCall(DMRestoreLocalVector(user->fda, &local_x));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumNewtonKrylov_FormResidual"
/**
 * @brief Adapts a PETSc trial vector to the existing momentum residual path.
 *
 * A matrix-free SNES residual must be a deterministic function of the trial
 * vector X alone: F(X) may not depend on any state left by a previous residual
 * or MFFD evaluation, or finite-difference Jacobian actions become inconsistent.
 * To honor that contract this callback fully derives the Cartesian velocity
 * state (Ucat/lUcat) from X before the first boundary sweep -- see the inline
 * comment below for why ApplyBoundaryConditions()'s own internal reconstruction
 * is not sufficient for the first outlet pass.
 *
 * State invariants:
 *   - On entry, X is the only input that determines the result; user->Ucont,
 *     user->Ucat and their local ghosts are treated as scratch and are fully
 *     overwritten from X.
 *   - Supported handlers may overwrite flux totals and other diagnostics on
 *     every call, but those values must not affect a later call at the same X;
 *     the deterministic seed guarantees this.
 *   - No histories, pressure, viscosity, or controller state advance here.
 *
 * Side effects: overwrites user->Ucont/lUcont, user->Ucat/lUcat, user->Rhs, the
 * boundary Ubcs targets, and boundary flux/area diagnostics; writes F.
 *
 * @param snes Calling nonlinear solver.
 * @param X Trial solution (read-only).
 * @param F Residual output.
 * @param vctx Pointer to MomentumNewtonKrylovContext.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MomentumNewtonKrylov_FormResidual(SNES snes, Vec X, Vec F, void *vctx)
{
    MomentumNewtonKrylovContext *ctx = (MomentumNewtonKrylovContext *)vctx;
    UserCtx *user = ctx->user;
    const char *staggered_fields[] = {"Ucont"};
    const char *cell_fields[] = {"Ucat"};

    PetscFunctionBeginUser;
    (void)snes;
    PetscCall(VecCopy(X, user->Ucont));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered_fields));

    /* Deterministic pre-boundary Cartesian seed. Establish the full Ucat/lUcat
     * state from the current X before any boundary handler runs:
     *
     *     X -> Ucont/lUcont -> Ucat -> periodic Ucat -> lUcat -> boundaries
     *
     * Why this is required, and why it is NOT redundant with the reconstruction
     * already performed inside ApplyBoundaryConditions():
     *
     *  1. A matrix-free SNES residual must be a deterministic function of X. If
     *     the Cartesian state is left over from a previous residual/MFFD call,
     *     F(X) becomes history dependent and the finite-difference Jacobian
     *     action Jv = (F(X+hv)-F(X))/h is invalidated.
     *  2. The conservation-outlet handler reads lUcat during the FIRST boundary
     *     sweep (it measures the uncorrected outflow and builds the outlet
     *     profile from the Cartesian field). Without this seed it would read the
     *     stale lUcat from the preceding evaluation.
     *  3. ApplyBoundaryConditions() does reconstruct Ucat/lUcat, but only AFTER
     *     each handler sweep (Contra2Cart runs after BoundarySystem_ExecuteStep
     *     within every pass). Those internal updates therefore prepare passes 2
     *     and 3 -- they cannot prepare the very first outlet read of pass 1.
     *  4. Contra2Cart() rebuilds the global Ucat interior from the current
     *     lUcont, but it does not by itself refresh lUcat (nor lUcont; that was
     *     done by the SynchronizePeriodicStaggeredFields call above).
     *  5. SynchronizePeriodicCellFields("Ucat") must run before the ghost
     *     scatter so periodic duplicate planes are finalized consistently (it is
     *     a no-op when no direction is periodic, as on the straight duct).
     *  6. UpdateLocalGhosts("Ucat") is required because the outlet handler reads
     *     lUcat -- the local ghosted vector -- not merely the global Ucat.
     *
     * Do NOT "simplify" this to a bare Contra2Cart(user), and do NOT delete it
     * as apparently redundant with ApplyBoundaryConditions(): the three internal
     * boundary passes remain necessary (they refresh the Cartesian state after
     * each boundary correction), but only this sequence makes pass 1's input a
     * deterministic function of X. */
    PetscCall(Contra2Cart(user));
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell_fields));
    PetscCall(UpdateLocalGhosts(user, "Ucat"));

    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(ComputeTotalResidual(user));
    PetscCall(VecCopy(user->Rhs, F));
    PetscCall(VecScale(F, -1.0));
    PetscCall(MomentumNewtonKrylov_ApplyConstraints(ctx, X, F));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumSolver_NewtonKrylov"
/**
 * @brief Runs one per-call matrix-free Newton--Krylov momentum solve.
 * @param user Single-block momentum context.
 * @param ibm Must be NULL in version one.
 * @param fsi Must be NULL in version one.
 * @return 0 on convergence or an error after restoring the canonical entry state.
 */
PetscErrorCode MomentumSolver_NewtonKrylov(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    PetscErrorCode ierr = PETSC_SUCCESS, cleanup_ierr;
    SimCtx        *simCtx;
    SNES           snes = NULL;
    Vec            solution = NULL, entry_backup = NULL;
    KSP            ksp = NULL;
    PC             pc = NULL;
    MomentumPreconditionerDescription preconditioner_description = {0};
    PetscBool      restore_entry = PETSC_FALSE;
    PetscBool      rhs_created = PETSC_FALSE;
    PetscBool      solve_started = PETSC_FALSE;
    PetscBool      committed = PETSC_FALSE;
    SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
    PetscInt       nonlinear_its = 0, function_evals = 0, linear_its = 0;
    PetscReal      final_norm = PETSC_MAX_REAL;
    MomentumNewtonKrylovContext ctx = {0};
    const char *staggered_fields[] = {"Ucont"};

    PetscFunctionBeginUser;
    PetscCall(MomentumNewtonKrylov_Validate(user));
    PetscCheck(ibm == NULL && fsi == NULL, PETSC_COMM_WORLD, PETSC_ERR_SUP,
               "Newton Krylov version one does not accept IBM or FSI objects.");
    PetscCheck(user->Rhs == NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Newton Krylov requires UserCtx::Rhs to be unallocated on entry.");
    simCtx = user->simCtx;

    ierr = VecDuplicate(user->Ucont, &solution); if (ierr) goto cleanup;
    ierr = VecDuplicate(user->Ucont, &entry_backup); if (ierr) goto cleanup;
    ierr = VecDuplicate(user->Ucont, &user->Rhs); if (ierr) goto cleanup;
    rhs_created = PETSC_TRUE;
    ierr = SNESCreate(PetscObjectComm((PetscObject)user->Ucont), &snes); if (ierr) goto cleanup;

    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); if (ierr) goto cleanup;
    ierr = ApplyBoundaryConditions(user); if (ierr) goto cleanup;
    ierr = VecCopy(user->Ucont, entry_backup); if (ierr) goto cleanup;
    restore_entry = PETSC_TRUE;
    ierr = VecCopy(user->Ucont, solution); if (ierr) goto cleanup;

    ctx.user = user;
    ierr = MomentumNewtonKrylov_ReadLinearizationConfig(
        &ctx.jacobian, &preconditioner_description); if (ierr) goto cleanup;
    ierr = SNESSetOptionsPrefix(snes, "mom_nk_"); if (ierr) goto cleanup;
    ierr = SNESSetType(snes, SNESNEWTONLS); if (ierr) goto cleanup;
    ierr = SNESSetDM(snes, user->fda); if (ierr) goto cleanup;
    ierr = SNESSetFunction(snes, NULL, MomentumNewtonKrylov_FormResidual, &ctx); if (ierr) goto cleanup;
    ierr = MomentumNewtonJacobian_Create(snes, &ctx.jacobian); if (ierr) goto cleanup;
    ierr = MomentumPreconditionerEngine_Create(user, ctx.jacobian.jacobian_operator,
                                                &preconditioner_description,
                                                &ctx.preconditioning_engine); if (ierr) goto cleanup;
    ierr = MomentumNewtonJacobian_Register(snes, &ctx.jacobian,
                                            &ctx.preconditioning_engine, &ctx); if (ierr) goto cleanup;
    ierr = SNESGetKSP(snes, &ksp); if (ierr) goto cleanup;
    ierr = KSPSetType(ksp, KSPGMRES); if (ierr) goto cleanup;
    ierr = KSPGetPC(ksp, &pc); if (ierr) goto cleanup;
    ierr = MomentumPreconditionerEngine_ConfigurePetscPC(&ctx.preconditioning_engine, pc); if (ierr) goto cleanup;
    ierr = SNESSetFromOptions(snes); if (ierr) goto cleanup;
    ierr = SNESMonitorSet(snes, MomentumNewtonKrylov_Monitor, &ctx, NULL); if (ierr) goto cleanup;
    ierr = MomentumPreconditionerEngine_ValidatePetscPC(&ctx.preconditioning_engine, pc); if (ierr) goto cleanup;

    LOG_ALLOW(GLOBAL, LOG_INFO,
              "Newton Krylov Jacobian: finite_difference / matrix_free; "
              "Preconditioner: %s; PETSc Jacobian matrix type: MATMFFD; "
              "PETSc PC type: %s.\n",
              ctx.preconditioning_engine.description.model == MOM_NK_PC_MODEL_NONE ?
                  "none" : "frozen_momentum_jacobian / point_block",
              ctx.preconditioning_engine.petsc_pc_type);

    MomentumNewtonKrylov_OpenHistory(&ctx);
    solve_started = PETSC_TRUE;
    ierr = SNESSolve(snes, NULL, solution);
    if (ierr) goto cleanup;
    ierr = SNESGetConvergedReason(snes, &reason); if (ierr) goto cleanup;
    ierr = SNESGetIterationNumber(snes, &nonlinear_its); if (ierr) goto cleanup;
    ierr = SNESGetNumberFunctionEvals(snes, &function_evals); if (ierr) goto cleanup;
    ierr = SNESGetLinearSolveIterations(snes, &linear_its); if (ierr) goto cleanup;
    ierr = SNESGetFunctionNorm(snes, &final_norm); if (ierr) goto cleanup;

    if (reason > 0) {
        ierr = VecCopy(solution, user->Ucont); if (ierr) goto cleanup;
        simCtx->mom_last_converged = PETSC_TRUE;
    } else {
        ierr = VecCopy(entry_backup, user->Ucont); if (ierr) goto cleanup;
        simCtx->mom_last_converged = PETSC_FALSE;
    }
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); if (ierr) goto cleanup;
    ierr = ApplyBoundaryConditions(user); if (ierr) goto cleanup;
    restore_entry = PETSC_FALSE;
    committed = (PetscBool)(reason > 0);

    LOG_ALLOW(GLOBAL, LOG_INFO,
              "Newton Krylov momentum solve: reason=%s (%d), Newton iterations=%d, residual evaluations=%d, Krylov iterations=%d, final norm=%.6e, state=%s.\n",
              SNESConvergedReasons[reason], (PetscInt)reason, nonlinear_its, function_evals,
              linear_its, (double)final_norm, reason > 0 ? "committed" : "rolled back");
    if (reason <= 0) ierr = PETSC_ERR_CONV_FAILED;

cleanup:
    /* A PETSc solve error can bypass the normal statistics path. Query whatever
       SNES retained without replacing the primary error so a failed attempt is
       still represented in the structured log. */
    if (solve_started && snes) {
        (void)SNESGetConvergedReason(snes, &reason);
        (void)SNESGetIterationNumber(snes, &nonlinear_its);
        (void)SNESGetNumberFunctionEvals(snes, &function_evals);
        (void)SNESGetLinearSolveIterations(snes, &linear_its);
        (void)SNESGetFunctionNorm(snes, &final_norm);
    }
    if (restore_entry && entry_backup) {
        cleanup_ierr = VecCopy(entry_backup, user->Ucont);
        if (!ierr) ierr = cleanup_ierr;
        cleanup_ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields);
        if (!ierr) ierr = cleanup_ierr;
        cleanup_ierr = ApplyBoundaryConditions(user);
        if (!ierr) ierr = cleanup_ierr;
        simCtx->mom_last_converged = PETSC_FALSE;
        committed = PETSC_FALSE;
    }
    if (ctx.history_file) {
        (void)fclose(ctx.history_file);
        ctx.history_file = NULL;
    }
    if (solve_started) {
        MomentumNewtonKrylov_WriteSummary(&ctx, reason, nonlinear_its, function_evals,
                                          linear_its, final_norm, committed);
    }
    if (rhs_created) {
        cleanup_ierr = VecDestroy(&user->Rhs);
        if (!ierr) ierr = cleanup_ierr;
    }
    cleanup_ierr = VecDestroy(&entry_backup); if (!ierr) ierr = cleanup_ierr;
    cleanup_ierr = VecDestroy(&solution); if (!ierr) ierr = cleanup_ierr;
    cleanup_ierr = MomentumPreconditionerEngine_Destroy(&ctx.preconditioning_engine); if (!ierr) ierr = cleanup_ierr;
    cleanup_ierr = MomentumNewtonJacobian_Destroy(&ctx.jacobian); if (!ierr) ierr = cleanup_ierr;
    cleanup_ierr = SNESDestroy(&snes); if (!ierr) ierr = cleanup_ierr;
    PetscFunctionReturn(ierr);
}
