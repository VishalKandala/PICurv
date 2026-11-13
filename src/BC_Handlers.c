#include "BC_Handlers.h"     // The header that declares this file's "constructor" functions
#include "expression.h"
#include "Metric.h"


//================================================================================
//              VALIDATORS
//================================================================================


/**
 * @brief Implementation of \ref Validate_DrivenFlowConfiguration().
 * @see Validate_DrivenFlowConfiguration()
 */
PetscErrorCode Validate_DrivenFlowConfiguration(UserCtx *user)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- CHECK 1: Detect if a driven flow is active. ---
    PetscBool is_driven_flow_active = PETSC_FALSE;
    char      driven_direction = ' ';
    const char* first_driven_face_name = "";

    for (int i = 0; i < 6; i++) {
        BCHandlerType handler_type = user->boundary_faces[i].handler_type;
        if (handler_type == BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX ||
            handler_type == BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX)
        {
            is_driven_flow_active = PETSC_TRUE;
            first_driven_face_name = BCFaceToString((BCFace)i);
            
            if (i <= 1)      driven_direction = 'X';
            else if (i <= 3) driven_direction = 'Y';
            else             driven_direction = 'Z';
            
            break; // Exit loop once we've confirmed it's active and found the direction.
        }
    }

    // If no driven flow handler is found, validation for this rule set is complete.
    if (!is_driven_flow_active) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "  - Driven Flow Handler detected on face %s. Applying driven flow validation rules...\n", first_driven_face_name);

    // --- CHECK 2: Ensure no conflicting BCs (Inlet/Outlet/Far-field) are present. ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "  - Checking for incompatible Inlet/Outlet/Far-field BCs...\n");
    for (int i = 0; i < 6; i++) {
        BCType math_type = user->boundary_faces[i].mathematical_type;
        if (math_type == INLET || math_type == OUTLET || math_type == FARFIELD) {
             SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                     "Configuration Error: A DRIVEN flow handler is active, which is incompatible with the %s boundary condition found on face %s.",
                     BCTypeToString(math_type), BCFaceToString((BCFace)i));
        }
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, " ... No conflicting BC types found. OK.\n");

    // --- CHECK 3: Ensure both ends of the driven direction have identical, valid setups. ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "      - Validating symmetry and mathematical types for the '%c' direction...\n", driven_direction);
    
    PetscInt neg_face_idx = 0, pos_face_idx = 0;
    if (driven_direction == 'X') {
        neg_face_idx = BC_FACE_NEG_X; pos_face_idx = BC_FACE_POS_X;
    } else if (driven_direction == 'Y') {
        neg_face_idx = BC_FACE_NEG_Y; pos_face_idx = BC_FACE_POS_Y;
    } else { // 'Z'
        neg_face_idx = BC_FACE_NEG_Z; pos_face_idx = BC_FACE_POS_Z;
    }

    BoundaryFaceConfig *neg_face_cfg = &user->boundary_faces[neg_face_idx];
    BoundaryFaceConfig *pos_face_cfg = &user->boundary_faces[pos_face_idx];

    // Rule 3a: Both faces must be PERIODIC.
    if (neg_face_cfg->mathematical_type != PERIODIC || pos_face_cfg->mathematical_type != PERIODIC) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                "Configuration Error: For a driven flow in the '%c' direction, both the %s and %s faces must be of mathematical_type PERIODIC.",
                driven_direction, BCFaceToString((BCFace)neg_face_idx), BCFaceToString((BCFace)pos_face_idx));
    }

    // Rule 3b: Both faces must use the exact same handler type.
    if (neg_face_cfg->handler_type != pos_face_cfg->handler_type) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                "Configuration Error: The DRIVEN handlers on the %s and %s faces of the '%c' direction do not match. Both must be the same type (e.g., both CONSTANT_FLUX).",
                BCFaceToString((BCFace)neg_face_idx), BCFaceToString((BCFace)pos_face_idx), driven_direction);
    }
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "          ... Symmetry and mathematical types are valid. OK.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/** @brief Cartesian boundary velocity one handler prescribes at one boundary face slot. */
typedef Cmpnts (*FaceBoundaryVelocityRule)(const BoundaryCondition *self, BCFace face_id,
                                           PetscInt i, PetscInt j, PetscInt k,
                                           Cmpnts metric, PetscReal sign, Cmpnts interior);

static PetscErrorCode ApplyFaceBoundaryVelocity(BoundaryCondition *self, BCContext *ctx,
                                                FaceBoundaryVelocityRule rule);

//================================================================================
//
//               HANDLER IMPLEMENTATION: NO-SLIP WALL
//               (Corresponds to BC_HANDLER_WALL_NOSLIP)
//
// This handler implements a stationary, impenetrable wall where the fluid
// velocity is zero (no-slip condition).
//
//================================================================================

// --- FORWARD DECLARATIONS ---
static PetscErrorCode Apply_WallNoSlip(BoundaryCondition *self, BCContext *ctx);

/**
 * @brief Implementation of \ref Create_WallNoSlip().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/BC_Handlers.h`.
 * @see Create_WallNoSlip()
 */
PetscErrorCode Create_WallNoSlip(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");

    // ✅ Set priority
    bc->priority   = BC_PRIORITY_WALL;
    
    // Assign function pointers
    bc->Initialize = NULL;
    bc->PreStep    = NULL;
    bc->Apply      = Apply_WallNoSlip;
    bc->PostStep   = NULL;
    bc->UpdateUbcs = NULL;
    bc->Destroy    = NULL;
    
    // No private data needed for this simple handler
    bc->data = NULL;

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief The no-slip wall's boundary velocity: the wall's own, which is at rest.
 */
static Cmpnts WallNoSlipVelocity(const BoundaryCondition *self, BCFace face_id,
                                 PetscInt i, PetscInt j, PetscInt k,
                                 Cmpnts metric, PetscReal sign, Cmpnts interior)
{
    const Cmpnts rest = {0.0, 0.0, 0.0};

    (void)self; (void)face_id; (void)i; (void)j; (void)k; (void)metric; (void)sign; (void)interior;
    return rest;
}

/**
 * @brief Apply no-slip velocity values to wall-adjacent cells for this boundary condition.
 * @details A zero boundary velocity also closes the face, since its flux is the
 *          velocity's product with the face metric.
 */
static PetscErrorCode Apply_WallNoSlip(BoundaryCondition *self, BCContext *ctx)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply_WallNoSlip: Applying to Face %d (%s).\n",
              ctx->face_id, BCFaceToString(ctx->face_id));
    PetscCall(ApplyFaceBoundaryVelocity(self, ctx, WallNoSlipVelocity));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//               HANDLER IMPLEMENTATION: FREE-SLIP WALL
//               (Corresponds to BC_HANDLER_WALL_FREE_SLIP)
//
// An impenetrable wall that exerts no shear: the normal velocity vanishes and the
// tangential velocity has zero normal gradient. It is also a plane of symmetry,
// the same condition read as the mirror image of the flow beyond it, which is why
// PICurv offers no separate symmetry boundary type.
//
//================================================================================

static PetscErrorCode Apply_WallFreeSlip(BoundaryCondition *self, BCContext *ctx);

/**
 * @brief Implementation of \ref Create_WallFreeSlip().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/BC_Handlers.h`.
 * @see Create_WallFreeSlip()
 */
PetscErrorCode Create_WallFreeSlip(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");
    bc->priority   = BC_PRIORITY_WALL;
    bc->Initialize = NULL;
    bc->PreStep    = NULL;
    bc->Apply      = Apply_WallFreeSlip;
    bc->PostStep   = NULL;
    /* The boundary value follows the adjacent velocity, so it is refreshed whenever
     * the velocity is: on every boundary pass and after each projection. */
    bc->UpdateUbcs = Apply_WallFreeSlip;
    bc->needs_repeat_pass = PETSC_TRUE; // Copies the adjacent Ucat, which the pass refreshes after it
    bc->Destroy    = NULL;
    bc->data       = NULL;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief The free-slip wall's boundary velocity: the adjacent velocity, less its normal part.
 * @details Removing the component along the face metric leaves a velocity tangent to the
 *          face, so the face flux vanishes; placing it as the boundary value makes the
 *          dummy cell the mirror image of the interior one, reversing the normal velocity
 *          and keeping the tangential one, which is zero normal gradient. On a curvilinear
 *          grid the metric is the face's own normal, so the construction holds there too.
 */
static inline Cmpnts WallFreeSlipVelocity(const BoundaryCondition *self, BCFace face_id,
                                          PetscInt i, PetscInt j, PetscInt k,
                                          Cmpnts metric, PetscReal sign, Cmpnts interior)
{
    (void)self; (void)face_id; (void)i; (void)j; (void)k; (void)sign;
    return RemoveFaceNormalComponent(interior, metric);
}

/**
 * @brief Apply free-slip boundary values to wall-adjacent cells for this boundary condition.
 */
static PetscErrorCode Apply_WallFreeSlip(BoundaryCondition *self, BCContext *ctx)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(ApplyFaceBoundaryVelocity(self, ctx, WallFreeSlipVelocity));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//          HANDLER IMPLEMENTATION: PRESCRIBED VELOCITY INLET
//          (Corresponds to BC_HANDLER_INLET_PRESCRIBED_VELOCITY)
//
// The velocity is a configured source (expression.h) evaluated at the centre of every
// inlet face slot, a PICSLICE file of one value per slot (a normal speed, or a Cartesian
// velocity), in one frame or a recording of frames interpolated in time, or such a file with
// expressions in its frame added to it. It is optionally rescaled so the face's mean normal
// speed or flow rate follows a target, the expressions alone when they add to a file, and
// synthetic fluctuations may be added after the rescaling. The slots are gathered at the
// first boundary pass, once the metrics exist; the velocity is re-evaluated once per step
// when it depends on t, and the rescaling is recomputed once per step when the source, the
// file or the target depends on t.
//
//================================================================================

static PetscErrorCode Initialize_InletPrescribedVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_InletPrescribedVelocity(BoundaryCondition *self, BCContext *ctx,
                                                      PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_InletVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PostStep_InletFlux(BoundaryCondition *self, BCContext *ctx,
                                         PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_WallMoving(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode Destroy_InletPrescribedVelocity(BoundaryCondition *self);

/**
 * @brief Private data of the prescribed-velocity inlet: a compiled source or a PICSLICE file,
 *        the values of the source's variables at each face slot this rank owns, and the
 *        velocity they last gave.
 */
typedef struct {
    PicurvVelocitySource source;             /* Normal frame with no expressions for a file alone. */
    char                *source_file;        /* Staged PICSLICE of solver-unit velocities, or NULL. */
    Cmpnts              *profile;            /* The file's velocity per slot when expressions add to it, or NULL. */
    FILE                *file;               /* The open PICSLICE, while later frames remain to be read. */
    fpos_t               first_frame;        /* Where the file's values begin. */
    PetscInt             frames;             /* Frames in the file; more than one varies in time. */
    PetscInt             components;         /* 1, a normal speed, or 3, a Cartesian velocity. */
    PetscInt             next_frame;         /* The frame the file's position is at. */
    PetscInt             file_n[2];          /* The file's slice dimensions. */
    PetscReal            interval;           /* Time between frames [s]. */
    PetscBool            repeat;             /* Past the last frame, return to the first. */
    PetscReal           *file_values;        /* One frame as the file stores it. */
    Cmpnts              *bracket[2];         /* This rank's slots of the frames around the step's time. */
    PetscInt             bracket_frame[2];   /* Which frames those are, or -1. */
    PicurvSourceFaceSlots slots;             /* What the source reads at each owned slot; t is refreshed. */
    Cmpnts              *velocity;           /* Solver-unit velocity per slot, the speed in x for the normal frame. */
    PetscReal            factor;             /* Bulk-velocity rescaling of `velocity`. */
    PetscInt             evaluated_step;     /* Step the cache was last refreshed for. */
    PetscBool            built;
    /* Synthetic fluctuations (digital filter, Xie & Castro 2008), added after `factor`. */
    PetscBool            fluctuate;
    PicurvVelocitySource stress[2];          /* Reynolds stresses [m^2/s^2]: (uu, vv, ww), then (uv, uw, vw). */
    PetscReal            fluct_length[2];    /* Integral lengths along the tangential axes, solver units. */
    PetscReal            fluct_time;         /* Integral time, solver units. */
    PetscInt             half[2];            /* Filter half-widths along the tangential axes, in cells. */
    PetscReal           *kernel[2];          /* Filter weights, 2 half + 1 per axis, unit sum of squares. */
    PetscReal          (*lund)[6];           /* Per slot, the stresses' Cholesky factor in solver units: a11 a21 a22 a31 a32 a33. */
    Cmpnts              *fluct;              /* Per slot, this step's fluctuation in solver units. */
    PetscBool            fluct_started;      /* The filter has a previous state to continue from. */
} InletPrescribedVelocityData;

/** @brief Index of the owned slot (i, j, k) in a prescribed-velocity inlet's per-slot arrays. */
static inline PetscInt PrescribedSlot(const InletPrescribedVelocityData *data, PetscInt i, PetscInt j, PetscInt k)
{
    const PetscInt idx[3] = {i, j, k};
    const PicurvSourceFaceSlots *slots = &data->slots;
    return (idx[slots->elements.a] - slots->elements.lo[0]) + slots->elements.n[0] * (idx[slots->elements.b] - slots->elements.lo[1]);
}

/**
 * @brief Implementation of \ref Create_InletPrescribedVelocity().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/BC_Handlers.h`.
 * @see Create_InletPrescribedVelocity()
 */
PetscErrorCode Create_InletPrescribedVelocity(BoundaryCondition *bc)
{
    InletPrescribedVelocityData *data = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");
    PetscCall(PetscNew(&data));
    bc->data       = (void*)data;
    bc->priority   = BC_PRIORITY_INLET;
    bc->Initialize = Initialize_InletPrescribedVelocity;
    bc->PreStep    = PreStep_InletPrescribedVelocity;
    bc->Apply      = Apply_InletVelocity;
    bc->PostStep   = PostStep_InletFlux;
    bc->UpdateUbcs = NULL;
    bc->Destroy    = Destroy_InletPrescribedVelocity;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref Create_WallMoving().
 * @see Create_WallMoving()
 */
PetscErrorCode Create_WallMoving(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(Create_InletPrescribedVelocity(bc));
    bc->priority = BC_PRIORITY_WALL;
    bc->Apply    = Apply_WallMoving;
    bc->PostStep = NULL;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Adds the inlet face's applied flux, the sum of its normal Ucont, to the inflow.
 */
static PetscErrorCode PostStep_InletFlux(BoundaryCondition *self, BCContext *ctx,
                                         PetscReal *local_inflow_contribution,
                                         PetscReal *local_outflow_contribution)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    PetscBool      can_service;
    
    (void)self;
    (void)local_outflow_contribution;
    
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    DMDALocalInfo *info = &user->info;
    Cmpnts ***ucont;

    PetscInt IM_nodes_global, JM_nodes_global,KM_nodes_global;

    IM_nodes_global = user->IM;
    JM_nodes_global = user->JM;
    KM_nodes_global = user->KM;
    
    ierr = CanRankServiceFace(info,IM_nodes_global,JM_nodes_global,KM_nodes_global,face_id,&can_service); CHKERRQ(ierr);
    

    if (!can_service) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    
    ierr = DMDAVecGetArrayRead(user->fda, user->Ucont, (const Cmpnts***)&ucont); CHKERRQ(ierr);

    PetscReal local_flux = 0.0;
    
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;
    PetscInt mx = info->mx, my = info->my, mz = info->mz;
    
    PetscInt lxs = xs, lxe = xe, lys = ys, lye = ye, lzs = zs, lze = ze;
    if (xs == 0) lxs = xs + 1;
    if (xe == mx) lxe = xe - 1;
    if (ys == 0) lys = ys + 1;
    if (ye == my) lye = ye - 1;
    if (zs == 0) lzs = zs + 1;
    if (ze == mz) lze = ze - 1;

    // Sum ucont components
    switch (face_id) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X: {
            PetscInt i = (face_id == BC_FACE_NEG_X) ? xs : mx - 2;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    local_flux += ucont[k][j][i].x;
                }
            }
        } break;
        
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y: {
            PetscInt j = (face_id == BC_FACE_NEG_Y) ? ys : my - 2;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    local_flux += ucont[k][j][i].y;
                }
            }
        } break;
        
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z: {
            PetscInt k = (face_id == BC_FACE_NEG_Z) ? zs : mz - 2;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    local_flux += ucont[k][j][i].z;
                }
            }
        } break;
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Ucont, (const Cmpnts***)&ucont); CHKERRQ(ierr);

    *local_inflow_contribution += local_flux;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "PostStep_InletFlux: Face %d, flux = %.6e\n",
              face_id, local_flux);

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Computes the expected PICSLICE dimensions for an inlet face.
 *
 * @param user User context containing global grid node counts.
 * @param face_id Boundary face whose tangential profile dimensions are requested.
 * @param[out] n1 First PICSLICE dimension in handler storage order.
 * @param[out] n2 Second PICSLICE dimension in handler storage order.
 * @return PetscErrorCode 0 on success, or a PETSc error for unsupported faces or invalid dimensions.
 */
static PetscErrorCode GetProfileFileExpectedDims(UserCtx *user, BCFace face_id,
                                                 PetscInt *n1, PetscInt *n2)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    switch (face_id) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X:
            *n1 = user->KM - 1;
            *n2 = user->JM - 1;
            break;
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y:
            *n1 = user->KM - 1;
            *n2 = user->IM - 1;
            break;
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z:
            *n1 = user->JM - 1;
            *n2 = user->IM - 1;
            break;
        default:
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                    "Unsupported face id %d for inlet profile dimensions.", face_id);
    }
    if (*n1 <= 0 || *n2 <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Invalid inlet profile dimensions (%" PetscInt_FMT ", %" PetscInt_FMT ") for grid (%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ").",
                *n1, *n2, user->IM, user->JM, user->KM);
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Opens a canonical PICSLICE file and reads its header.
 *
 * @details The file holds magic token `PICSLICE`, a frame count, then a line with the slice
 *          dimensions, optionally the components per slot (1, a normal speed, the default, or 3,
 *          a Cartesian velocity) and, with several frames, the time between them. Frames follow.
 *
 * @param data        Handler data; `file`, `first_frame`, `frames`, `components`, `interval` and
 *                    `file_n` are set.
 * @param expected_n1 Expected first slice dimension.
 * @param expected_n2 Expected second slice dimension.
 * @return PetscErrorCode 0 on success, or a PETSc file/validation error on malformed input.
 */
static PetscErrorCode OpenPicSlice(InletPrescribedVelocityData *data, PetscInt expected_n1, PetscInt expected_n2)
{
    const char *source_file = data->source_file;
    FILE *fd = NULL;
    char magic[32] = {0}, line[256] = {0};
    PetscInt frame_count = 0, n1 = 0, n2 = 0;
    PetscInt components = 1;
    double interval = 0.0;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    fd = fopen(source_file, "r");
    if (!fd) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                     "Cannot open PICSLICE inlet profile file: %s", source_file);

    if (fscanf(fd, "%31s", magic) != 1 || strcmp(magic, "PICSLICE") != 0) {
        fclose(fd);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                "PICSLICE inlet profile file %s must begin with PICSLICE header.", source_file);
    }
    if (fscanf(fd, "%" PetscInt_FMT " ", &frame_count) != 1) {
        fclose(fd);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                "PICSLICE inlet profile file %s missing frame count.", source_file);
    }
    if (!fgets(line, sizeof(line), fd) ||
        sscanf(line, "%" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT " %le", &n1, &n2, &components, &interval) < 2) {
        fclose(fd);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                "PICSLICE inlet profile file %s missing slice dimensions.", source_file);
    }
    if (n1 != expected_n1 || n2 != expected_n2) {
        fclose(fd);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                "PICSLICE inlet profile dimensions mismatch for %s: expected (%" PetscInt_FMT ", %" PetscInt_FMT "), found (%" PetscInt_FMT ", %" PetscInt_FMT ").",
                source_file, expected_n1, expected_n2, n1, n2);
    }
    if (frame_count < 1 || (components != 1 && components != 3) || (frame_count > 1 && !(interval > 0.0))) {
        fclose(fd);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                "PICSLICE inlet profile file %s has %" PetscInt_FMT " frames of %" PetscInt_FMT " components, %g apart; "
                "it needs at least one frame, 1 or 3 components, and a positive interval between several frames.",
                source_file, frame_count, components, interval);
    }
    if (fgetpos(fd, &data->first_frame) != 0) {
        fclose(fd);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Cannot record a position in %s.", source_file);
    }
    data->file = fd;
    data->frames = frame_count;
    data->components = components;
    data->interval = interval;
    data->file_n[0] = n1;
    data->file_n[1] = n2;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Reads the frame at the file's position into `file_values`, validating every value: finite,
 *        and a normal speed nonnegative. After the last frame the file must end.
 * @param data Handler data with an open file.
 * @return PetscErrorCode 0 on success, or a PETSc file/validation error on malformed input.
 */
static PetscErrorCode ReadPicSliceFrame(InletPrescribedVelocityData *data)
{
    const PetscInt count = data->file_n[0] * data->file_n[1] * data->components;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (PetscInt idx = 0; idx < count; idx++) {
        PetscReal value = 0.0;
        if (fscanf(data->file, "%le", &value) != 1) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                    "PICSLICE inlet profile file %s ended early in frame %" PetscInt_FMT ": expected %" PetscInt_FMT " values.",
                    data->source_file, data->next_frame, count);
        }
        if (PetscIsInfOrNanReal(value) || (data->components == 1 && value < 0.0)) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                    "PICSLICE inlet profile file %s contains invalid speed %.6e at flat index %" PetscInt_FMT ".",
                    data->source_file, (double)value, idx);
        }
        data->file_values[idx] = value;
    }
    data->next_frame++;
    if (data->next_frame == data->frames) {
        char extra[64];
        if (fscanf(data->file, "%63s", extra) == 1) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                    "PICSLICE inlet profile file %s has extra token after %" PetscInt_FMT " frames of %" PetscInt_FMT " values: %s",
                    data->source_file, data->frames, count, extra);
        }
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Prescribed velocity at one inlet face slot, as a Cartesian vector.
 * @details A normal-frame source stores its speed in `.x`; it is lifted onto the signed unit
 *          face normal. Other frames store the Cartesian velocity. Both carry the bulk factor.
 */
static inline Cmpnts EvaluateInletCartesianVelocity(const BoundaryCondition *self,
                                                     BCFace face_id,
                                                     PetscInt i, PetscInt j, PetscInt k,
                                                     Cmpnts metric, PetscReal sign, Cmpnts interior)
{
    (void)interior;
    const InletPrescribedVelocityData *data = (const InletPrescribedVelocityData*)self->data;
    Cmpnts velocity = {0.0, 0.0, 0.0};
    if (!data->built) return velocity;
    const Cmpnts v = data->velocity[PrescribedSlot(data, i, j, k)];
    const Cmpnts staged = data->profile ? data->profile[PrescribedSlot(data, i, j, k)] : velocity;
    /* Fluctuations follow the bulk factor, which rescales only the mean. */
    const Cmpnts fluct = data->fluct ? data->fluct[PrescribedSlot(data, i, j, k)] : velocity;
    if (data->source.frame != PICURV_SOURCE_FRAME_NORMAL) {
        velocity.x = data->factor * v.x + staged.x + fluct.x;
        velocity.y = data->factor * v.y + staged.y + fluct.y;
        velocity.z = data->factor * v.z + staged.z + fluct.z;
        return velocity;
    }
    const PetscReal area = sqrt(metric.x * metric.x + metric.y * metric.y + metric.z * metric.z);
    if (area <= 0.0) return velocity;
    const PetscReal speed = data->factor * v.x + staged.x;
    velocity.x = sign * speed * metric.x / area + fluct.x;
    velocity.y = sign * speed * metric.y / area + fluct.y;
    velocity.z = sign * speed * metric.z / area + fluct.z;
    return velocity;
}

/**
 * @brief Sets a handler's boundary velocity and the face flux it implies on one face.
 * @details Shared by every handler that prescribes a Cartesian velocity at its face: the
 *          inlets, the no-slip wall and the free-slip wall. For each boundary slot this
 *          rank owns it asks @p rule for the velocity, stores it as the boundary value,
 *          and sets the face flux to its product with the face metric. A rule whose
 *          velocity is tangent to the face therefore closes the face. A slot whose fluid
 *          side is blanked (`nvert > 0.1`) belongs to an immersed body, whose own boundary
 *          condition applies there, so it is left untouched.
 * @param self Handler being applied.
 * @param ctx  Boundary context naming the face.
 * @param rule Velocity the handler prescribes, given the slot, the face metric, the
 *             side (`sign`, positive on a negative face) and the adjacent interior velocity.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode ApplyFaceBoundaryVelocity(BoundaryCondition *self, BCContext *ctx,
                                                FaceBoundaryVelocityRule rule)
{
    UserCtx *user = ctx->user;
    BCFace face_id = ctx->face_id;
    PetscBool can_service;
    DMDALocalInfo *info = &user->info;
    Cmpnts ***ubcs, ***ucont, ***csi, ***eta, ***zet, ***ucat;
    PetscReal ***nvert;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    PetscCall(CanRankServiceFace(info, user->IM, user->JM, user->KM, face_id, &can_service));
    if (!can_service) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    PetscCall(DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcat, (const Cmpnts***)&ucat));

    const PetscInt xs = info->xs, xe = info->xs + info->xm;
    const PetscInt ys = info->ys, ye = info->ys + info->ym;
    const PetscInt zs = info->zs, ze = info->zs + info->zm;
    const PetscInt mx = info->mx, my = info->my, mz = info->mz;
    PetscInt lxs = xs, lxe = xe, lys = ys, lye = ye, lzs = zs, lze = ze;
    if (xs == 0) lxs++;
    if (xe == mx) lxe--;
    if (ys == 0) lys++;
    if (ye == my) lye--;
    if (zs == 0) lzs++;
    if (ze == mz) lze--;

    switch (face_id) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X: {
            const PetscReal sign = (face_id == BC_FACE_NEG_X) ? 1.0 : -1.0;
            const PetscInt i = (face_id == BC_FACE_NEG_X) ? xs : mx - 2;
            const PetscInt ib = i + (sign < 0);
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    if ((sign > 0 && nvert[k][j][i + 1] > 0.1) ||
                        (sign < 0 && nvert[k][j][i] > 0.1)) continue;
                    const Cmpnts metric = csi[k][j][i];
                    const Cmpnts velocity = rule(self, face_id, i, j, k, metric, sign,
                                                 (sign > 0) ? ucat[k][j][i + 1] : ucat[k][j][i]);
                    ubcs[k][j][ib] = velocity;
                    ucont[k][j][i].x = velocity.x * metric.x + velocity.y * metric.y + velocity.z * metric.z;
                }
            }
        } break;
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y: {
            const PetscReal sign = (face_id == BC_FACE_NEG_Y) ? 1.0 : -1.0;
            const PetscInt j = (face_id == BC_FACE_NEG_Y) ? ys : my - 2;
            const PetscInt jb = j + (sign < 0);
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if ((sign > 0 && nvert[k][j + 1][i] > 0.1) ||
                        (sign < 0 && nvert[k][j][i] > 0.1)) continue;
                    const Cmpnts metric = eta[k][j][i];
                    const Cmpnts velocity = rule(self, face_id, i, j, k, metric, sign,
                                                 (sign > 0) ? ucat[k][j + 1][i] : ucat[k][j][i]);
                    ubcs[k][jb][i] = velocity;
                    ucont[k][j][i].y = velocity.x * metric.x + velocity.y * metric.y + velocity.z * metric.z;
                }
            }
        } break;
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z: {
            const PetscReal sign = (face_id == BC_FACE_NEG_Z) ? 1.0 : -1.0;
            const PetscInt k = (face_id == BC_FACE_NEG_Z) ? zs : mz - 2;
            const PetscInt kb = k + (sign < 0);
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if ((sign > 0 && nvert[k + 1][j][i] > 0.1) ||
                        (sign < 0 && nvert[k][j][i] > 0.1)) continue;
                    const Cmpnts metric = zet[k][j][i];
                    const Cmpnts velocity = rule(self, face_id, i, j, k, metric, sign,
                                                 (sign > 0) ? ucat[k + 1][j][i] : ucat[k][j][i]);
                    ubcs[kb][j][i] = velocity;
                    ucont[k][j][i].z = velocity.x * metric.x + velocity.y * metric.y + velocity.z * metric.z;
                }
            }
        } break;
    }

    PetscCall(DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcat, (const Cmpnts***)&ucat));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief A moving wall's boundary velocity: the source's, with any part along the face removed, since
 *        a wall closes its face. PreStep refuses a source with more than a rounding of the
 *        surface's facets through it.
 */
static inline Cmpnts EvaluateMovingWallVelocity(const BoundaryCondition *self, BCFace face_id, PetscInt i, PetscInt j, PetscInt k,
                                                Cmpnts metric, PetscReal sign, Cmpnts interior)
{
    return RemoveFaceNormalComponent(EvaluateInletCartesianVelocity(self, face_id, i, j, k, metric, sign, interior), metric);
}

/** @brief Applies a moving wall's velocity through the path every prescribed face shares. */
static PetscErrorCode Apply_WallMoving(BoundaryCondition *self, BCContext *ctx)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(ApplyFaceBoundaryVelocity(self, ctx, EvaluateMovingWallVelocity));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Applies the prescribed inlet velocity through the common face-layout path.
 *
 * @details Ownership clipping, immersed-cell exclusion, metric projection, `Ubcs` placement,
 *          and the normal staggered `Ucont` write are `ApplyFaceBoundaryVelocity()`'s, shared
 *          with the walls.
 */
static PetscErrorCode Apply_InletVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(ApplyFaceBoundaryVelocity(self, ctx, EvaluateInletCartesianVelocity));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Compile the source named by the face's `source_id` parameter, with its bulk velocity or
 *        flow rate if any, record the PICSLICE file named by its `source_file` parameter, and read
 *        the synthetic-fluctuation settings `-face_source_<id>_fluct*`: the Reynolds-stress sources
 *        and the integral length and time scales.
 */
static PetscErrorCode Initialize_InletPrescribedVelocity(BoundaryCondition *self, BCContext *ctx)
{
    InletPrescribedVelocityData *data = (InletPrescribedVelocityData*)self->data;
    UserCtx   *user = ctx->user;
    PetscReal  id = -1.0;
    PetscBool  found = PETSC_FALSE;
    const char *source_file = NULL;
    char       prefix[64];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(GetBCParamReal(user->boundary_faces[ctx->face_id].params, "source_id", &id, &found));
    PetscCheck(found && id >= 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Face %s: the prescribed_velocity inlet needs a source_id parameter.", BCFaceToString(ctx->face_id));
    PetscCall(GetBCParamString(user->boundary_faces[ctx->face_id].params, "source_file", &source_file, &found));
    if (found && source_file && source_file[0] != '\0') {
        PetscCall(PetscStrallocpy(source_file, &data->source_file));
    }
    PetscCall(PetscSNPrintf(prefix, sizeof(prefix), "face_source_%d", (int)id));
    PetscCall(PicurvVelocitySourceCreate(prefix, PETSC_TRUE, (PetscBool)(data->source_file != NULL), PETSC_TRUE, &data->source));
    PetscCheck(!data->source_file || data->source.frame == PICURV_SOURCE_FRAME_NORMAL ||
               data->source.frame == PICURV_SOURCE_FRAME_CARTESIAN, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Face %s: a profile file and the expressions added to it are in the normal or cartesian frame.",
               BCFaceToString(ctx->face_id));
    PetscCall(GetBCParamBool(user->boundary_faces[ctx->face_id].params, "source_repeat", &data->repeat, &found));
    data->factor = 1.0;
    data->evaluated_step = -1;
    {
        char      option[96], stress_prefix[96];
        PetscInt  count = 2;
        PetscCall(PetscSNPrintf(option, sizeof(option), "-%s_fluct", prefix));
        PetscCall(PetscOptionsGetBool(NULL, NULL, option, &data->fluctuate, NULL));
        if (data->fluctuate) {
            PetscCheck(user->InflowFluctuation, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                       "Face %s adds fluctuations, but -inflow_fluctuations did not allocate their state.",
                       BCFaceToString(ctx->face_id));
            for (PetscInt g = 0; g < 2; ++g) {
                PetscCall(PetscSNPrintf(stress_prefix, sizeof(stress_prefix), "%s_fluct_stress_%d", prefix, (int)g));
                PetscCall(PicurvVelocitySourceCreate(stress_prefix, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE, &data->stress[g]));
                /* Three stress components, as read: no frame turns them, and no bulk velocity or
                   flow rate scales them. */
                PetscCheck(data->stress[g].frame == PICURV_SOURCE_FRAME_CARTESIAN && !data->stress[g].bulk,
                           PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                           "Face %s: the Reynolds stresses -%s_* must be in the cartesian frame, with no bulk_velocity or flow_rate.",
                           BCFaceToString(ctx->face_id), stress_prefix);
            }
            PetscCall(PetscSNPrintf(option, sizeof(option), "-%s_fluct_length", prefix));
            PetscCall(PetscOptionsGetRealArray(NULL, NULL, option, data->fluct_length, &count, &found));
            PetscCheck(found && count == 2 && data->fluct_length[0] > 0.0 && data->fluct_length[1] > 0.0,
                       PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "%s needs two positive length scales.", option);
            PetscCall(PetscSNPrintf(option, sizeof(option), "-%s_fluct_time", prefix));
            PetscCall(PetscOptionsGetReal(NULL, NULL, option, &data->fluct_time, &found));
            PetscCheck(found && data->fluct_time > 0.0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                       "%s needs a positive time scale.", option);
        }
    }
    if (data->source_file) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Face %s: prescribed velocity from %s%s.\n",
                  BCFaceToString(ctx->face_id), data->source_file, data->source.bulk ? ", rescaled to a target" : "");
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Face %s: prescribed velocity from -face_source_%d (%s frame%s%s).\n",
                  BCFaceToString(ctx->face_id), (int)id,
                  data->source.frame == PICURV_SOURCE_FRAME_NORMAL ? "normal" :
                  data->source.frame == PICURV_SOURCE_FRAME_CURVILINEAR ? "curvilinear" : "cartesian",
                  data->source.time_dependent ? ", time-dependent" : "", data->source.bulk ? ", bulk velocity" : "");
    }
    if (data->fluctuate) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Face %s: synthetic fluctuations, integral lengths %g and %g, time %g (solver units).\n",
                  BCFaceToString(ctx->face_id), (double)data->fluct_length[0], (double)data->fluct_length[1],
                  (double)data->fluct_time);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief This rank's area-weighted sums over the face, physical: the source's normal speed, the
 *        staged data's when expressions add to it, the area, and the synthetic fluctuations'
 *        inflow (zero without them). Slots blanked by an immersed body are excluded, as in Apply.
 *        The first three are what the bulk rescaling reads; the fluctuations are added after it.
 */
static PetscErrorCode PrescribedFaceSums(const InletPrescribedVelocityData *data, UserCtx *user, BCFace face_id, PetscReal local[4])
{
    const PicurvSourceFaceSlots *slots = &data->slots;
    const PetscInt normal = slots->elements.normal, a = slots->elements.a, b = slots->elements.b;
    const PetscBool negative = (PetscBool)((PetscInt)face_id % 2 == 0);
    const PetscReal sign = negative ? 1.0 : -1.0, U_ref = user->simCtx->scaling.U_ref;
    const PetscReal L_ref = user->simCtx->scaling.L_ref;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    local[0] = local[1] = local[2] = local[3] = 0.0;
    if (slots->elements.n[0] * slots->elements.n[1] > 0) {
        const PetscInt slot = slots->elements.layer;
        const PetscInt cell = negative ? slot + 1 : slot;
        const Cmpnts ***metric_arr;
        const PetscReal ***nvert;
        Vec metric_vec = (normal == 0) ? user->lCsi : (normal == 1) ? user->lEta : user->lZet;
        PetscCall(DMDAVecGetArrayRead(user->fda, metric_vec, &metric_arr));
        PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, &nvert));
        for (PetscInt q = 0; q < slots->elements.n[1]; ++q)
        for (PetscInt p = 0; p < slots->elements.n[0]; ++p) {
            PetscInt s_idx[3], c_idx[3];
            s_idx[normal] = slot; s_idx[a] = slots->elements.lo[0] + p; s_idx[b] = slots->elements.lo[1] + q;
            c_idx[normal] = cell; c_idx[a] = s_idx[a]; c_idx[b] = s_idx[b];
            if (nvert[c_idx[2]][c_idx[1]][c_idx[0]] > 0.1) continue;
            const Cmpnts m = metric_arr[s_idx[2]][s_idx[1]][s_idx[0]];
            const Cmpnts v = data->velocity[p + slots->elements.n[0] * q];
            const PetscReal area = PetscSqrtReal(m.x * m.x + m.y * m.y + m.z * m.z);
            const PetscReal physical = U_ref * L_ref * L_ref;
            local[0] += physical * ((data->source.frame == PICURV_SOURCE_FRAME_NORMAL)
                                    ? v.x * area : sign * (v.x * m.x + v.y * m.y + v.z * m.z));
            if (data->profile) {
                const Cmpnts staged = data->profile[p + slots->elements.n[0] * q];
                local[1] += physical * ((data->source.frame == PICURV_SOURCE_FRAME_NORMAL)
                                        ? staged.x * area : sign * (staged.x * m.x + staged.y * m.y + staged.z * m.z));
            }
            local[2] += L_ref * L_ref * area;
            /* Fluctuations are Cartesian in every frame (Apply adds them to the velocity). */
            if (data->fluct) {
                const Cmpnts f = data->fluct[p + slots->elements.n[0] * q];
                local[3] += physical * sign * (f.x * m.x + f.y * m.y + f.z * m.z);
            }
        }
        PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, &nvert));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, metric_vec, &metric_arr));
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Factor that makes the face's mean normal speed or flow rate equal the target. Collective.
 */
static PetscErrorCode PrescribedBulkFactor(InletPrescribedVelocityData *data, UserCtx *user, BCFace face_id,
                                           PetscReal physical_time, PetscReal *factor)
{
    PetscReal local[4];
    char      where[32];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PrescribedFaceSums(data, user, face_id, local));
    PetscCall(PetscSNPrintf(where, sizeof(where), "Face %s", BCFaceToString(face_id)));
    PetscCall(PicurvVelocitySourceBulkFactor(&data->source, physical_time, local, where, factor));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Read one frame of the face's PICSLICE into this rank's slots, moving the file there.
 *        Every rank reads the file, so a malformed one fails on all of them.
 * @param data   Handler data with an open file.
 * @param frame  Frame to read.
 * @param target Velocity per owned slot; a normal speed goes in `x`.
 */
static PetscErrorCode LoadPrescribedFrame(InletPrescribedVelocityData *data, PetscInt frame, Cmpnts *target)
{
    const PicurvSourceFaceSlots *slots = &data->slots;
    const PetscInt n2 = data->file_n[1], components = data->components;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (frame < data->next_frame) {
        PetscCheck(fsetpos(data->file, &data->first_frame) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                   "Cannot return to the first frame of %s.", data->source_file);
        data->next_frame = 0;
    }
    while (data->next_frame <= frame) PetscCall(ReadPicSliceFrame(data));
    /* A PICSLICE row runs along the face's first tangential axis; rows step along the second. */
    for (PetscInt q = 0; q < slots->elements.n[1]; ++q)
    for (PetscInt p = 0; p < slots->elements.n[0]; ++p) {
        const PetscReal *value = &data->file_values[components * ((slots->elements.lo[1] + q - 1) * n2 + (slots->elements.lo[0] + p - 1))];
        Cmpnts *v = &target[p + slots->elements.n[0] * q];
        v->x = value[0];
        v->y = (components == 3) ? value[1] : 0.0;
        v->z = (components == 3) ? value[2] : 0.0;
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Open the face's PICSLICE and read a lone frame, as the velocity itself or, under
 *        expressions, as what they add to; several frames are read as time reaches them.
 */
static PetscErrorCode OpenPrescribedProfile(InletPrescribedVelocityData *data, UserCtx *user, BCFace face_id)
{
    const PetscInt n_slots = data->slots.elements.n[0] * data->slots.elements.n[1];
    const PetscBool normal = (PetscBool)(data->source.frame == PICURV_SOURCE_FRAME_NORMAL);
    PetscInt n1 = 0, n2 = 0;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(GetProfileFileExpectedDims(user, face_id, &n1, &n2));
    PetscCall(OpenPicSlice(data, n1, n2));
    PetscCheck(data->components == (normal ? 1 : 3), PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
               "Face %s: %s holds %" PetscInt_FMT " component(s) per slot, but the source is in the %s frame.",
               BCFaceToString(face_id), data->source_file, data->components, normal ? "normal" : "cartesian");
    if (data->source.count > 0) PetscCall(PetscCalloc1(n_slots, &data->profile));
    PetscCall(PetscMalloc1(n1 * n2 * data->components, &data->file_values));
    if (data->frames == 1) {
        PetscCall(LoadPrescribedFrame(data, 0, data->profile ? data->profile : data->velocity));
        PetscCheck(fclose(data->file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Cannot close %s.", data->source_file);
        data->file = NULL;
    } else {
        /* Two allocations: the brackets swap as time advances, so neither stays the first of a pair. */
        PetscCall(PetscMalloc1(n_slots, &data->bracket[0]));
        PetscCall(PetscMalloc1(n_slots, &data->bracket[1]));
        data->bracket_frame[0] = data->bracket_frame[1] = -1;
    }
    LOG_ALLOW(GLOBAL, LOG_INFO, "Face %s: PICSLICE %s, %" PetscInt_FMT " x %" PetscInt_FMT ", %" PetscInt_FMT " component(s), %" PetscInt_FMT " frame(s).\n",
              BCFaceToString(face_id), data->source_file, n1, n2, data->components, data->frames);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Set the staged data at a time: linear between the two frames around it, the first at t = 0.
 */
static PetscErrorCode InterpolatePrescribedFrames(InletPrescribedVelocityData *data, BCFace face_id, PetscReal physical_time)
{
    const PetscInt n_slots = data->slots.elements.n[0] * data->slots.elements.n[1];
    const PetscReal position = physical_time / data->interval;
    PetscInt  frame[2] = {(PetscInt)PetscFloorReal(position), 0};
    PetscReal weight = position - PetscFloorReal(position);
    Cmpnts   *out = data->profile ? data->profile : data->velocity;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (data->repeat) {
        frame[0] = frame[0] % data->frames;
        frame[1] = (frame[0] + 1) % data->frames;
    } else {
        /* A time step landing on the last frame may round just past it. */
        if (frame[0] == data->frames - 1 && weight < 1.0e-9) { frame[0] = data->frames - 2; weight = 1.0; }
        PetscCheck(frame[0] < data->frames - 1 && frame[0] >= 0, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                   "Face %s: time %g is past the %" PetscInt_FMT " frames of %s, %g s apart. Set repeat: true to loop them.",
                   BCFaceToString(face_id), (double)physical_time, data->frames, data->source_file, (double)data->interval);
        frame[1] = frame[0] + 1;
    }
    for (PetscInt f = 0; f < 2; ++f) {
        if (data->bracket_frame[f] == frame[f]) continue;
        if (data->bracket_frame[1 - f] == frame[f]) {
            /* The later frame of the last interval is the earlier one of this. */
            Cmpnts  *kept = data->bracket[1 - f];
            data->bracket[1 - f] = data->bracket[f];
            data->bracket[f] = kept;
            data->bracket_frame[1 - f] = data->bracket_frame[f];
            data->bracket_frame[f] = frame[f];
            continue;
        }
        PetscCall(LoadPrescribedFrame(data, frame[f], data->bracket[f]));
        data->bracket_frame[f] = frame[f];
    }
    for (PetscInt s = 0; s < n_slots; ++s) {
        out[s].x = (1.0 - weight) * data->bracket[0][s].x + weight * data->bracket[1][s].x;
        out[s].y = (1.0 - weight) * data->bracket[0][s].y + weight * data->bracket[1][s].y;
        out[s].z = (1.0 - weight) * data->bracket[0][s].z + weight * data->bracket[1][s].z;
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref LogPrescribedFlowDiagnostics().
 * @see LogPrescribedFlowDiagnostics()
 */
PetscErrorCode LogPrescribedFlowDiagnostics(UserCtx *user)
{
    SimCtx *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (PetscInt face = 0; face < 6; ++face) {
        const BoundaryCondition *handler = user->boundary_faces[face].handler;
        PetscReal local[4], global[4], physical_time, flow_rate;
        FILE *file = NULL;

        if (!handler || user->boundary_faces[face].handler_type != BC_HANDLER_INLET_PRESCRIBED_VELOCITY) continue;
        const InletPrescribedVelocityData *data = (const InletPrescribedVelocityData*)handler->data;
        if (!data->built) continue;
        PetscCall(PrescribedFaceSums(data, user, (BCFace)face, local));
        PetscCallMPI(MPI_Allreduce(local, global, 4, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
        if (simCtx->rank != 0 || user->_this != 0) continue;
        PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
        flow_rate = data->factor * global[0] + global[1];
        PetscCall(PicurvOpenDiagnosticsCsv(simCtx, simCtx->analysis_dir, "inlet_flow.csv",
                                           "step,physical_time,face,area,flow_rate,mean_normal_speed,bulk_factor,"
                                           "fluctuation_flow_rate,applied_flow_rate",
                                           &file));
        /* flow_rate is the rescaled profile's, the one a bulk target controls; the synthetic
           fluctuations are added after the rescaling, so the flux the face actually admits is
           their sum. */
        fprintf(file, "%d,%.10e,%s,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n", (int)simCtx->step, (double)physical_time,
                BCFaceToString((BCFace)face), (double)global[2], (double)flow_rate,
                (double)(global[2] > 0.0 ? flow_rate / global[2] : 0.0), (double)data->factor,
                (double)global[3], (double)(flow_rate + global[3]));
        PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE, "Unable to close the inlet-flow diagnostics file.");
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/** Face cells the synthetic-inflow draw lattice extends beyond each edge; a filter half-width may not exceed it. */
#define INFLOW_LATTICE_PAD 1048576

/**
 * @brief Prepare an inlet's synthetic fluctuations once its slots exist: filter widths from the
 *        face's mean spacing, and the Reynolds stresses' Cholesky factor at each slot.
 * @details Collective over the block. The integral length L along a tangential axis spans
 *          n = L / h cells, h the face's cell width along it: the mean of its widest and
 *          narrowest cell there. The filter weights are
 *          exp(-pi |k| / n) for |k| <= 2n, scaled to a unit sum of squares (Xie & Castro 2008),
 *          so filtered unit-variance noise keeps unit variance. Every rank and every
 *          decomposition uses the same h, so the kernel is bitwise the same.
 */
static PetscErrorCode BuildInletFluctuation(InletPrescribedVelocityData *data, UserCtx *user)
{
    PicurvSourceFaceSlots *slots = &data->slots;
    const PetscInt         count = slots->elements.n[0] * slots->elements.n[1];
    const PetscInt         tangential[2] = {slots->elements.a, slots->elements.b};
    const PetscReal        U2 = user->simCtx->scaling.U_ref * user->simCtx->scaling.U_ref;
    /* Widest and narrowest cell along each axis, as {max0, max1, -min0, -min1}. Each width comes
       from its own cell's metrics and a maximum is exact in any order, so every decomposition
       builds the same kernel; a line length is summed across the ranks it crosses, and is not. */
    PetscReal              extent[4] = {0.0, 0.0, -PETSC_MAX_REAL, -PETSC_MAX_REAL};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (PetscInt s = 0; s < count; ++s) {
        for (PetscInt d = 0; d < 2; ++d) {
            const Cmpnts    e = data->slots.edges[s][tangential[d]];
            const PetscReal width = PetscSqrtReal(e.x * e.x + e.y * e.y + e.z * e.z);
            extent[d] = PetscMax(extent[d], width);
            extent[2 + d] = PetscMax(extent[2 + d], -width);
        }
    }
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, extent, 4, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    for (PetscInt d = 0; d < 2; ++d) {
        const PetscReal h = 0.5 * (extent[d] - extent[2 + d]);
        const PetscReal n = data->fluct_length[d] / h;
        PetscReal       norm = 0.0;
        data->half[d] = PetscMax(1, (PetscInt)PetscCeilReal(2.0 * n));
        PetscCheck(data->half[d] <= INFLOW_LATTICE_PAD, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                   "An inlet's fluctuation length spans %g cells, beyond the %d its draws can reach.",
                   (double)n, INFLOW_LATTICE_PAD / 2);
        PetscCall(PetscMalloc1(2 * data->half[d] + 1, &data->kernel[d]));
        for (PetscInt k = -data->half[d]; k <= data->half[d]; ++k) {
            data->kernel[d][k + data->half[d]] = PetscExpReal(-PETSC_PI * PetscAbsReal((PetscReal)k) / n);
            norm += data->kernel[d][k + data->half[d]] * data->kernel[d][k + data->half[d]];
        }
        for (PetscInt k = 0; k <= 2 * data->half[d]; ++k) data->kernel[d][k] /= PetscSqrtReal(norm);
    }
    PetscCall(PetscMalloc2(count, &data->lund, count, &data->fluct));
    for (PetscInt s = 0; s < count; ++s) {
        Cmpnts diag, cross;
        PetscCall(PicurvVelocitySourceEvaluate(&data->stress[0], data->slots.values[s], NULL, data->slots.edges[s], &diag));
        PetscCall(PicurvVelocitySourceEvaluate(&data->stress[1], data->slots.values[s], NULL, data->slots.edges[s], &cross));
        /* Lund, Wu & Squires (1998): u' = A psi with A A^T = R. A stress tensor that is not
           positive definite, such as zero at a wall, leaves the rows it cannot fill at zero. */
        PetscReal *a = data->lund[s];
        a[0] = PetscSqrtReal(PetscMax(diag.x / U2, 0.0));
        a[1] = (a[0] > 0.0) ? cross.x / U2 / a[0] : 0.0;
        a[2] = PetscSqrtReal(PetscMax(diag.y / U2 - a[1] * a[1], 0.0));
        a[3] = (a[0] > 0.0) ? cross.y / U2 / a[0] : 0.0;
        a[4] = (a[2] > 0.0) ? (cross.z / U2 - a[1] * a[3]) / a[2] : 0.0;
        a[5] = PetscSqrtReal(PetscMax(diag.z / U2 - a[3] * a[3] - a[4] * a[4], 0.0));
        data->fluct[s].x = data->fluct[s].y = data->fluct[s].z = 0.0;
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief This step's synthetic fluctuation at every owned slot of an inlet.
 * @details Keyed unit normals on the owned face cells and a halo of the filter's half-width are
 *          filtered along each tangential axis; a draw is a pure function of the seed, block,
 *          component, face cell and step, so a rank computes its halo itself and the result does
 *          not depend on the decomposition. A periodic tangential axis wraps the halo. The time
 *          filter psi <- psi exp(-pi dt / 2T) + filtered sqrt(1 - exp(-pi dt / T)) keeps unit
 *          variance; its state lives in the inlet's dummy plane of `InflowFluctuation`, which a
 *          checkpoint carries, and starts from the filtered noise when there is none. A restored
 *          state belongs to the restart step already, so the first pass after it only reads it.
 */
static PetscErrorCode StepInletFluctuation(InletPrescribedVelocityData *data, UserCtx *user, BCFace face_id)
{
    SimCtx                      *simCtx = user->simCtx;
    const PicurvSourceFaceSlots *slots = &data->slots;
    const PetscInt               n0 = slots->elements.n[0], n1 = slots->elements.n[1], h0 = data->half[0], h1 = data->half[1];
    const PetscInt               w0 = n0 + 2 * h0, w1 = n1 + 2 * h1;
    const PetscInt               nodes[3] = {user->info.mx, user->info.my, user->info.mz};
    const PetscInt               cells[2] = {nodes[slots->elements.a] - 2, nodes[slots->elements.b] - 2};
    PetscBool                    axes_periodic[3];
    GetPeriodicAxes(simCtx, axes_periodic);
    const PetscBool              periodic[2] = {axes_periodic[slots->elements.a], axes_periodic[slots->elements.b]};
    const PetscBool              negative = (PetscBool)((PetscInt)face_id % 2 == 0);
    const PetscReal              decay = PetscExpReal(-PETSC_PI * simCtx->dt / (2.0 * data->fluct_time));
    const PetscReal              gain = PetscSqrtReal(1.0 - decay * decay);
    const PetscBool              fresh = (PetscBool)(!data->fluct_started && !simCtx->inflowFluctuationsRestored);
    /* A restored state is already this step's: the checkpoint was written after its update. */
    const PetscBool              restored = (PetscBool)(!data->fluct_started && simCtx->inflowFluctuationsRestored);
    PetscReal                   *noise = NULL, *rows = NULL, *psi = NULL;
    Cmpnts                    ***state = NULL;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (n0 * n1 > 0) {
        PetscCall(PetscMalloc3(w0 * w1, &noise, n0 * w1, &rows, 3 * n0 * n1, &psi));
        for (PetscInt c = 0; c < 3 && !restored; ++c) {
            for (PetscInt qq = 0; qq < w1; ++qq) {
                for (PetscInt pp = 0; pp < w0; ++pp) {
                    PetscInt g[2] = {slots->elements.lo[0] - 1 - h0 + pp, slots->elements.lo[1] - 1 - h1 + qq};
                    for (PetscInt d = 0; d < 2; ++d) if (periodic[d]) g[d] = ((g[d] % cells[d]) + cells[d]) % cells[d];
                    /* A face cell, or one of its halo beyond a closed edge, on a lattice padded by a fixed
                       width, so a draw does not depend on the filter's half-widths. */
                    const PetscInt64 id = ((PetscInt64)face_id * (cells[1] + 2 * INFLOW_LATTICE_PAD) + (g[1] + INFLOW_LATTICE_PAD))
                                          * (cells[0] + 2 * INFLOW_LATTICE_PAD) + (g[0] + INFLOW_LATTICE_PAD);
                    noise[pp + w0 * qq] = PicurvNormalFromHash(PicurvDrawHash((PetscInt64)simCtx->sourceRandomSeed,
                                                                               PICURV_DRAW_SPACE_INFLOW,
                                                                               (uint64_t)(3 * user->_this + c), id,
                                                                               (PetscInt64)simCtx->step));
                }
            }
            for (PetscInt qq = 0; qq < w1; ++qq) {
                for (PetscInt p = 0; p < n0; ++p) {
                    PetscReal sum = 0.0;
                    for (PetscInt k = 0; k <= 2 * h0; ++k) sum += data->kernel[0][k] * noise[(p + k) + w0 * qq];
                    rows[p + n0 * qq] = sum;
                }
            }
            for (PetscInt q = 0; q < n1; ++q) {
                for (PetscInt p = 0; p < n0; ++p) {
                    PetscReal sum = 0.0;
                    for (PetscInt k = 0; k <= 2 * h1; ++k) sum += data->kernel[1][k] * rows[p + n0 * (q + k)];
                    psi[3 * (p + n0 * q) + c] = sum;
                }
            }
        }
        PetscCall(DMDAVecGetArray(user->fda, user->InflowFluctuation, &state));
        for (PetscInt q = 0; q < n1; ++q) {
            for (PetscInt p = 0; p < n0; ++p) {
                const PetscInt   s = p + n0 * q;
                PetscInt         idx[3];
                idx[slots->elements.normal] = negative ? 0 : nodes[slots->elements.normal] - 1;
                idx[slots->elements.a] = slots->elements.lo[0] + p;
                idx[slots->elements.b] = slots->elements.lo[1] + q;
                Cmpnts          *st = &state[idx[2]][idx[1]][idx[0]];
                const PetscReal *f = &psi[3 * s];
                if (fresh) { st->x = f[0]; st->y = f[1]; st->z = f[2]; }
                else if (!restored) { st->x = decay * st->x + gain * f[0]; st->y = decay * st->y + gain * f[1]; st->z = decay * st->z + gain * f[2]; }
                const PetscReal *a = data->lund[s];
                data->fluct[s].x = a[0] * st->x;
                data->fluct[s].y = a[1] * st->x + a[2] * st->y;
                data->fluct[s].z = a[3] * st->x + a[4] * st->y + a[5] * st->z;
            }
        }
        PetscCall(DMDAVecRestoreArray(user->fda, user->InflowFluctuation, &state));
        PetscCall(PetscFree3(noise, rows, psi));
    }
    data->fluct_started = PETSC_TRUE;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Refresh the inlet velocity once per step: gather the slots on the first call, then
 *        re-evaluate and rescale only what depends on t.
 */
static PetscErrorCode PreStep_InletPrescribedVelocity(BoundaryCondition *self, BCContext *ctx,
                                                      PetscReal *local_inflow_contribution,
                                                      PetscReal *local_outflow_contribution)
{
    InletPrescribedVelocityData *data = (InletPrescribedVelocityData*)self->data;
    UserCtx  *user = ctx->user;
    SimCtx   *simCtx = user->simCtx;
    PetscReal physical_time;

    (void)local_inflow_contribution;
    (void)local_outflow_contribution;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    /* Boundary passes repeat within a step (every stage and residual); the velocity is
       defined at the step's time, so it is refreshed once. */
    if (data->built && data->evaluated_step == simCtx->step) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    const PetscBool first = (PetscBool)!data->built;
    if (first) {
        Vec fraction, length;
        PetscCall(DMCreateGlobalVector(user->fda, &fraction));
        PetscCall(DMCreateGlobalVector(user->fda, &length));
        PetscCall(ComputeGridLineFractions(user, fraction, length));
        PetscCall(PicurvSourceFaceSlotsCreate(user, ctx->face_id, fraction, length, &data->slots));
        PetscCall(VecDestroy(&fraction));
        PetscCall(VecDestroy(&length));
        PetscCall(PetscMalloc1(data->slots.elements.n[0] * data->slots.elements.n[1], &data->velocity));
        if (data->source_file) PetscCall(OpenPrescribedProfile(data, user, ctx->face_id));
        if (data->fluctuate) PetscCall(BuildInletFluctuation(data, user));
    }
    const PetscBool frames_vary = (PetscBool)(data->frames > 1);
    const PetscBool evaluate = (PetscBool)(first || data->source.time_dependent);
    const PetscBool rescale = (PetscBool)(data->source.bulk && (evaluate || frames_vary || data->source.bulk_time_dependent));
    PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
    if (frames_vary) PetscCall(InterpolatePrescribedFrames(data, ctx->face_id, physical_time));
    if (evaluate && data->source.count > 0) {
        const PetscReal U_ref = simCtx->scaling.U_ref;
        const PetscInt  nodes[3] = {user->info.mx, user->info.my, user->info.mz};
        for (PetscInt s = 0; s < data->slots.elements.n[0] * data->slots.elements.n[1]; ++s) {
            const PetscInt p = s % data->slots.elements.n[0], q = s / data->slots.elements.n[0];
            const PicurvExpressionDrawKey key = PicurvFaceDrawKey((PetscInt64)simCtx->sourceRandomSeed, user->_this, nodes, ctx->face_id,
                                                                  data->slots.elements.lo[0] + p - 1, data->slots.elements.lo[1] + q - 1,
                                                                  (PetscInt64)simCtx->step);
            Cmpnts v;
            data->slots.values[s][PICURV_SOURCE_VAR_T] = physical_time;
            PetscCall(PicurvVelocitySourceEvaluate(&data->source, data->slots.values[s], &key, data->slots.edges[s], &v));
            data->velocity[s].x = v.x / U_ref;
            data->velocity[s].y = v.y / U_ref;
            data->velocity[s].z = v.z / U_ref;
            if (user->boundary_faces[ctx->face_id].handler_type == BC_HANDLER_WALL_MOVING) {
                const Cmpnts m = data->slots.edges[s][data->slots.elements.normal];
                const PetscReal speed2 = v.x * v.x + v.y * v.y + v.z * v.z;
                const PetscReal through = v.x * m.x + v.y * m.y + v.z * m.z;
                PetscCheck(through * through <= 2.5e-3 * speed2 * (m.x * m.x + m.y * m.y + m.z * m.z), PETSC_COMM_SELF, PETSC_ERR_USER_INPUT,
                           "Face %s: a moving wall's velocity must lie along the wall; (%g, %g, %g) passes through it at (%g, %g, %g).",
                           BCFaceToString(ctx->face_id), (double)v.x, (double)v.y, (double)v.z,
                           (double)data->slots.values[s][PICURV_SOURCE_VAR_X], (double)data->slots.values[s][PICURV_SOURCE_VAR_Y],
                           (double)data->slots.values[s][PICURV_SOURCE_VAR_Z]);
            }
        }
    }
    data->built = PETSC_TRUE;
    if (rescale) PetscCall(PrescribedBulkFactor(data, user, ctx->face_id, physical_time, &data->factor));
    if (data->fluctuate) PetscCall(StepInletFluctuation(data, user, ctx->face_id));
    data->evaluated_step = simCtx->step;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/** @brief Release the compiled source and the per-slot arrays. */
static PetscErrorCode Destroy_InletPrescribedVelocity(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (self && self->data) {
        InletPrescribedVelocityData *data = (InletPrescribedVelocityData*)self->data;
        PetscCall(PicurvVelocitySourceDestroy(&data->source));
        PetscCall(PetscFree(data->source_file));
        if (data->file) PetscCheck(fclose(data->file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Cannot close a PICSLICE file.");
        PetscCall(PetscFree(data->profile));
        PetscCall(PetscFree(data->file_values));
        PetscCall(PetscFree(data->bracket[0]));
        PetscCall(PetscFree(data->bracket[1]));
        if (data->built) PetscCall(PicurvSourceFaceSlotsDestroy(&data->slots));
        PetscCall(PetscFree(data->velocity));
        if (data->fluctuate) {
            PetscCall(PicurvVelocitySourceDestroy(&data->stress[0]));
            PetscCall(PicurvVelocitySourceDestroy(&data->stress[1]));
            PetscCall(PetscFree(data->kernel[0]));
            PetscCall(PetscFree(data->kernel[1]));
            PetscCall(PetscFree2(data->lund, data->fluct));
        }
        PetscCall(PetscFree(self->data));
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//          HANDLER IMPLEMENTATION: OUTLET WITH MASS CONSERVATION
//          (Corresponds to BC_HANDLER_OUTLET_CONSERVATION)
//
// This handler ensures that the total flux leaving through outlet boundaries
// balances the total flux entering through inlet and far-field boundaries.
//
// Workflow:
// 1. PreStep: Measures the *uncorrected* flux based on interior velocities.
// 2. Apply:   Calculates a global correction factor based on the flux imbalance
//             and applies it to the contravariant velocity (ucont) on the outlet face.
// 3. PostStep: Measures the *corrected* flux for verification and logging.
//
// Backflow clip. The uniform correction adds the same velocity c to every face,
// so it can leave a face with flux into the domain. Once per time step, at the
// step's first Apply, the handler asks whether F + c*A points into the domain on
// any face of this outlet (F the zero-gradient flux, signed along +xi like every
// face flux, so outward is F > 0 on a positive face and F < 0 on a negative one;
// A the face area). If not, the uniform rule
// runs unchanged. If so, the faces whose F points into the domain are clipped to
// zero flux (and zero normal velocity in Ubcs), and the rest are scaled by one factor
// psi so that this face carries the same total as the uniform rule would. The
// rule and the clipped set are frozen for the rest of the step so the residual
// stays smooth inside a Newton solve; psi is re-measured on every call, as c is,
// so the total stays exact.
//
//================================================================================

// --- 1. FORWARD DECLARATIONS ---
static PetscErrorCode PreStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                     PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution);
static PetscErrorCode Apply_OutletConservation(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PostStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                                  PetscReal *in, PetscReal *out);
static PetscErrorCode Destroy_OutletConservation(BoundaryCondition *self);

/** Per-face state of the backflow clip; see the section comment above. */
typedef struct {
    PetscInt   decisionStep; /* Step at which the rule and clipped set were chosen (-1 = never). */
    PetscBool  clipRule;     /* PETSC_TRUE: clip reversed faces and scale the rest, this step. */
    PetscInt   maskSize;     /* Owned face elements covered by clipped[]. */
    PetscBool *clipped;      /* Frozen clipped set, indexed q*n[0] + p over the owned elements. */
} OutletConservationData;

/**
 * @brief Implementation of \ref Create_OutletConservation().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/BC_Handlers.h`.
 * @see Create_OutletConservation()
 */
PetscErrorCode Create_OutletConservation(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");

    // This handler has the highest priority to ensure it runs after
    // all inflow fluxes have been calculated.
    bc->priority   = BC_PRIORITY_OUTLET;
    
    // Assign function pointers
    bc->Initialize = NULL; // No initialization needed
    bc->PreStep    = PreStep_OutletConservation;
    bc->Apply      = Apply_OutletConservation;
    bc->PostStep   = PostStep_OutletConservation;
    bc->UpdateUbcs = NULL;
    bc->needs_repeat_pass = PETSC_TRUE; // Reads the adjacent Ucat and the inlets' flux
    bc->Destroy    = Destroy_OutletConservation;

    OutletConservationData *data = NULL;
    PetscCall(PetscNew(&data));
    data->decisionStep = -1;
    bc->data = (void*)data;

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Prepare the outlet-conservation correction before advancing the solver.
 */
static PetscErrorCode PreStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                     PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    DMDALocalInfo* info = &user->info;
    PetscBool      can_service;

    // Suppress unused parameter warnings for clarity.
    (void)self;
    (void)local_inflow_contribution;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // Step 1: Use the robust utility function to determine if this MPI rank owns a computable
    // portion of the specified boundary face. If not, there is no work to do, so we exit immediately.
    const PetscInt IM_nodes_global = user->IM;
    const PetscInt JM_nodes_global = user->JM;
    const PetscInt KM_nodes_global = user->KM;
    ierr = CanRankServiceFace(info, IM_nodes_global, JM_nodes_global, KM_nodes_global, face_id, &can_service); CHKERRQ(ierr);

    if (!can_service) {
        PROFILE_HELPER_END;
        PetscFunctionReturn(0);
    }

    // Step 2: Get read-only access to the necessary PETSc arrays.
    // We use the local versions (`lUcat`, `lNvert`) which include ghost cell data,
    // ensuring we have the correct interior values adjacent to the boundary.
    Cmpnts ***ucat, ***csi, ***eta, ***zet;
    PetscReal ***nvert;
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcat, (const Cmpnts***)&ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);

    PetscReal local_flux_out = 0.0;
    const PetscInt xs=info->xs, xe=info->xs+info->xm;
    const PetscInt ys=info->ys, ye=info->ys+info->ym;
    const PetscInt zs=info->zs, ze=info->zs+info->zm;
    const PetscInt mx=info->mx, my=info->my, mz=info->mz;

    // Step 3: Replicate the legacy shrunk loop bounds to exclude corners and edges.
    PetscInt lxs = xs; if (xs == 0)    lxs = xs + 1;
    PetscInt lxe = xe; if (xe == mx)   lxe = xe - 1;
    PetscInt lys = ys; if (ys == 0)    lys = ys + 1;
    PetscInt lye = ye; if (ye == my)   lye = ye - 1;
    PetscInt lzs = zs; if (zs == 0)    lzs = zs + 1;
    PetscInt lze = ze; if (ze == mz)   lze = ze - 1;

    // Step 4: Loop over the specified face using the corrected bounds and indexing to calculate flux.
    switch (face_id) {
        case BC_FACE_NEG_X: {
            const PetscInt i_cell = xs + 1; // Index for first interior cell-centered data
            const PetscInt i_face = xs;     // Index for the -X face of that cell
            for (int k=lzs; k<lze; k++) for (int j=lys; j<lye; j++) {
                if (nvert[k][j][i_cell] < 0.1) {
                    local_flux_out += (ucat[k][j][i_cell].x * csi[k][j][i_face].x + ucat[k][j][i_cell].y * csi[k][j][i_face].y + ucat[k][j][i_cell].z * csi[k][j][i_face].z);
                }
            }
            break;
        }
        case BC_FACE_POS_X: {
            const PetscInt i_cell = xe - 2; // Index for last interior cell-centered data
            const PetscInt i_face = xe - 2; // Index for the +X face of that cell
            for (int k=lzs; k<lze; k++) for (int j=lys; j<lye; j++) {
                if (nvert[k][j][i_cell] < 0.1) {
                    local_flux_out += (ucat[k][j][i_cell].x * csi[k][j][i_face].x + ucat[k][j][i_cell].y * csi[k][j][i_face].y + ucat[k][j][i_cell].z * csi[k][j][i_face].z);
                }
            }
            break;
        }
        case BC_FACE_NEG_Y: {
            const PetscInt j_cell = ys + 1;
            const PetscInt j_face = ys;
            for (int k=lzs; k<lze; k++) for (int i=lxs; i<lxe; i++) {
                if (nvert[k][j_cell][i] < 0.1) {
                    local_flux_out += (ucat[k][j_cell][i].x * eta[k][j_face][i].x + ucat[k][j_cell][i].y * eta[k][j_face][i].y + ucat[k][j_cell][i].z * eta[k][j_face][i].z);
                }
            }
            break;
        }
        case BC_FACE_POS_Y: {
            const PetscInt j_cell = ye - 2;
            const PetscInt j_face = ye - 2;
            for (int k=lzs; k<lze; k++) for (int i=lxs; i<lxe; i++) {
                if (nvert[k][j_cell][i] < 0.1) {
                    local_flux_out += (ucat[k][j_cell][i].x * eta[k][j_face][i].x + ucat[k][j_cell][i].y * eta[k][j_face][i].y + ucat[k][j_cell][i].z * eta[k][j_face][i].z);
                }
            }
            break;
        }
        case BC_FACE_NEG_Z: {
            const PetscInt k_cell = zs + 1;
            const PetscInt k_face = zs;
            for (int j=lys; j<lye; j++) for (int i=lxs; i<lxe; i++) {
                if (nvert[k_cell][j][i] < 0.1) {
                    local_flux_out += (ucat[k_cell][j][i].x * zet[k_face][j][i].x + ucat[k_cell][j][i].y * zet[k_face][j][i].y + ucat[k_cell][j][i].z * zet[k_face][j][i].z);
                }
            }
            break;
        }
        case BC_FACE_POS_Z: {
            const PetscInt k_cell = ze - 2;
            const PetscInt k_face = ze - 2;
            for (int j=lys; j<lye; j++) for (int i=lxs; i<lxe; i++) {
                if (nvert[k_cell][j][i] < 0.1) {
                    local_flux_out += (ucat[k_cell][j][i].x * zet[k_face][j][i].x + ucat[k_cell][j][i].y * zet[k_face][j][i].y + ucat[k_cell][j][i].z * zet[k_face][j][i].z);
                }
            }
            break;
        }
    }

    // Step 5: Restore the PETSc arrays.
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcat, (const Cmpnts***)&ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);

    // Step 6: Add this face's calculated flux to the accumulator for this rank.
    *local_outflow_contribution += local_flux_out;

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/** What one OutletClipPass() over the owned outlet elements does. */
typedef enum { OUTLET_CLIP_DECIDE, OUTLET_CLIP_MEASURE, OUTLET_CLIP_WRITE } OutletClipPassMode;

/**
 * @brief One pass over this rank's elements of an outlet face for the backflow clip.
 *
 * Every fluid element's zero-gradient flux is F = Ucat(cell) . S, the flux the uniform rule
 * starts from. The pass is local; the caller reduces @p sums.
 * - DECIDE freezes the clipped set {F points into the domain} into @p data and returns
 *   sums[0] = the number of elements where F + correction*A points into the domain and
 *   sums[1] = the number clipped. Fluxes are signed along +xi, so outward is +F on a positive
 *   face and -F on a negative one.
 * - MEASURE returns sums[0] = sum F, sums[1] = sum A and sums[2] = sum F over unclipped elements.
 * - WRITE sets the outlet flux: zero on clipped elements, whose Ubcs loses its normal component,
 *   and psi*F elsewhere, with Ubcs = Ucat as in the uniform rule.
 */
static PetscErrorCode OutletClipPass(UserCtx *user, BCFace face_id, OutletConservationData *data,
                                     OutletClipPassMode mode, PetscReal correction, PetscReal psi,
                                     PetscReal sums[3])
{
    const PetscBool      negative = (PetscBool)((PetscInt)face_id % 2 == 0);
    const PetscReal      outward = negative ? -1.0 : 1.0;
    BoundaryFaceElements e;
    Cmpnts            ***ucat, ***metric, ***ubcs = NULL, ***ucont = NULL;
    PetscReal         ***nvert;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    sums[0] = sums[1] = sums[2] = 0.0;
    PetscCall(GetBoundaryFaceElements(user, face_id, &e));
    const PetscInt count = e.n[0] * e.n[1];
    if (mode == OUTLET_CLIP_DECIDE && count != data->maskSize) {
        PetscCall(PetscFree(data->clipped));
        if (count > 0) PetscCall(PetscCalloc1(count, &data->clipped));
        data->maskSize = count;
    }
    if (count == 0) { PROFILE_HELPER_END; PetscFunctionReturn(PETSC_SUCCESS); }

    const Vec      lmetric = (e.normal == 0) ? user->lCsi : (e.normal == 1) ? user->lEta : user->lZet;
    const PetscInt cell  = negative ? e.layer + 1 : e.layer;
    const PetscInt dummy = negative ? e.layer : e.layer + 1;

    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcat, (const Cmpnts***)&ucat));
    PetscCall(DMDAVecGetArrayRead(user->fda, lmetric, (const Cmpnts***)&metric));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert));
    if (mode == OUTLET_CLIP_WRITE) {
        PetscCall(DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    }
    for (PetscInt q = 0; q < e.n[1]; q++)
    for (PetscInt p = 0; p < e.n[0]; p++) {
        PetscInt c[3], f[3], d[3];

        c[e.normal] = cell; c[e.a] = e.lo[0] + p; c[e.b] = e.lo[1] + q;
        if (nvert[c[2]][c[1]][c[0]] >= 0.1) continue;
        f[0] = d[0] = c[0]; f[1] = d[1] = c[1]; f[2] = d[2] = c[2];
        f[e.normal] = e.layer;
        d[e.normal] = dummy;

        const PetscInt  idx  = q * e.n[0] + p;
        const Cmpnts    S    = metric[f[2]][f[1]][f[0]];
        const Cmpnts    u    = ucat[c[2]][c[1]][c[0]];
        const PetscReal flux = u.x * S.x + u.y * S.y + u.z * S.z;
        const PetscReal area = PetscSqrtReal(S.x * S.x + S.y * S.y + S.z * S.z);

        switch (mode) {
            case OUTLET_CLIP_DECIDE:
                data->clipped[idx] = (PetscBool)(outward * flux < 0.0);
                if (outward * (flux + correction * area) < 0.0) sums[0] += 1.0;
                if (data->clipped[idx]) sums[1] += 1.0;
                break;
            case OUTLET_CLIP_MEASURE:
                sums[0] += flux;
                sums[1] += area;
                if (!data->clipped[idx]) sums[2] += flux;
                break;
            case OUTLET_CLIP_WRITE: {
                Cmpnts    *face_flux   = &ucont[f[2]][f[1]][f[0]];
                PetscReal *normal_flux = (e.normal == 0) ? &face_flux->x : (e.normal == 1) ? &face_flux->y : &face_flux->z;

                if (data->clipped[idx]) {
                    *normal_flux = 0.0;
                    ubcs[d[2]][d[1]][d[0]] = RemoveFaceNormalComponent(u, S);
                } else {
                    *normal_flux = psi * flux;
                    ubcs[d[2]][d[1]][d[0]] = u;
                }
                break;
            }
        }
    }
    if (mode == OUTLET_CLIP_WRITE) {
        PetscCall(DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs));
        PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, lmetric, (const Cmpnts***)&metric));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcat, (const Cmpnts***)&ucat));
    PROFILE_HELPER_END;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Run the backflow clip for one Apply call; see the section comment.
 *
 * Collective: every rank calls it, whether or not it serves the face. Sets @p applied when the
 * clip rule wrote this face's fluxes; otherwise the caller applies the uniform rule unchanged.
 */
static PetscErrorCode OutletConservationClip(UserCtx *user, BCFace face_id, OutletConservationData *data,
                                             PetscReal correction, PetscBool *applied)
{
    SimCtx         *simCtx = user->simCtx;
    const PetscReal outward = ((PetscInt)face_id % 2 == 0) ? -1.0 : 1.0;
    PetscReal       local[3], global[3];
    PetscBool decided_now = PETSC_FALSE;
    PetscReal reversed = 0.0, clipped = 0.0;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    *applied = PETSC_FALSE;
    if (data->decisionStep != simCtx->step) {
        PetscCall(OutletClipPass(user, face_id, data, OUTLET_CLIP_DECIDE, correction, 0.0, local));
        PetscCallMPI(MPI_Allreduce(local, global, 2, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
        data->clipRule     = (PetscBool)(global[0] > 0.5);
        data->decisionStep = simCtx->step;
        decided_now        = PETSC_TRUE;
        reversed           = global[0];
        clipped            = global[1];
    }
    if (!data->clipRule) { PROFILE_HELPER_END; PetscFunctionReturn(PETSC_SUCCESS); }

    PetscCall(OutletClipPass(user, face_id, data, OUTLET_CLIP_MEASURE, correction, 0.0, local));
    PetscCallMPI(MPI_Allreduce(local, global, 3, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    /* The total the uniform rule would give this face: the clip only redistributes it. */
    const PetscReal target = global[0] + correction * global[1];
    if (!(outward * global[2] > 0.0 && outward * target > 0.0)) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Outlet %s, step %" PetscInt_FMT ": nothing left to scale (unclipped flux %.4e, target %.4e); "
                  "using the uniform correction for this pass.\n",
                  BCFaceToString(face_id), simCtx->step, global[2], target);
        PROFILE_HELPER_END;
        PetscFunctionReturn(PETSC_SUCCESS);
    }
    const PetscReal psi = target / global[2];
    PetscCall(OutletClipPass(user, face_id, data, OUTLET_CLIP_WRITE, correction, psi, local));
    if (decided_now) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Outlet %s, step %" PetscInt_FMT ": uniform correction would leave %.0f faces reversed; "
                  "clipped %.0f reversed faces, scaled the rest by psi = %.6f.\n",
                  BCFaceToString(face_id), simCtx->step, reversed, clipped, psi);
    }
    *applied = PETSC_TRUE;
    PROFILE_HELPER_END;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief (Handler Action) Applies mass conservation correction to the outlet face.
 *
 * This function calculates a global correction factor based on the total inflow and outflow fluxes
 * and applies it to the contravariant velocity (`ucont`) on the outlet face to ensure mass conservation.
 */
static PetscErrorCode Apply_OutletConservation(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    DMDALocalInfo* info = &user->info;
    PetscBool      can_service;
    PetscBool      clip_applied = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- STEP 1: Calculate the correction factor using pre-calculated area ---
    PetscReal total_inflow = *ctx->global_inflow_sum + *ctx->global_farfield_inflow_sum;
    PetscReal flux_imbalance = total_inflow - *ctx->global_outflow_sum;
    
    // Directly use the pre-calculated area from the simulation context.
    PetscReal velocity_correction = (PetscAbsReal(user->simCtx->AreaOutSum) > 1e-12) 
                                    ? flux_imbalance / user->simCtx->AreaOutSum 
                                    : 0.0;
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Outlet Correction on Face %d: Imbalance=%.4e, Pre-calc Area=%.4e, V_corr=%.4e\n",
              face_id, flux_imbalance, user->simCtx->AreaOutSum, velocity_correction);

    // Backflow clip (collective, so before the can_service exit). When it applies, it has
    // written this face; otherwise the uniform rule below runs unchanged.
    ierr = OutletConservationClip(user, face_id, (OutletConservationData*)self->data,
                                  velocity_correction, &clip_applied); CHKERRQ(ierr);
    if (clip_applied) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    const PetscInt IM_nodes_global = user->IM;
    const PetscInt JM_nodes_global = user->JM;
    const PetscInt KM_nodes_global = user->KM;
    ierr = CanRankServiceFace(info, IM_nodes_global, JM_nodes_global, KM_nodes_global, face_id, &can_service); CHKERRQ(ierr);
    
    if (!can_service) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // --- STEP 2: Apply the correction to ucont on the outlet face ---

    // Get read/write access to necessary arrays

    Cmpnts ***ubcs, ***ucont, ***csi, ***eta, ***zet, ***ucat;
    PetscReal ***nvert;
    ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda,user->lUcat, (const Cmpnts***)&ucat); CHKERRQ(ierr);  
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // Get local grid bounds to exclude corners/edges
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;
    PetscInt mx = info->mx, my = info->my, mz = info->mz;
    PetscInt lxs = xs, lxe = xe, lys = ys, lye = ye, lzs = zs, lze = ze;

    if (xs == 0) lxs = xs + 1;
    if (xe == mx) lxe = xe - 1;
    if (ys == 0) lys = ys + 1;
    if (ye == my) lye = ye - 1;
    if (zs == 0) lzs = zs + 1;          
    if (ze == mz) lze = ze - 1;

    // Loop over faces and apply correction    
    switch(face_id){
        case BC_FACE_NEG_X:{
            const PetscInt i_cell = xs + 1;
            const PetscInt i_face = xs;
            const PetscInt i_dummy = xs;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    if (nvert[k][j][i_cell] < 0.1) {
                        // Set ubcs
                        ubcs[k][j][i_dummy] = ucat[k][j][i_cell];

                        // Calculate Local uncorrected original flux
                        PetscReal Uncorrected_local_flux = (ubcs[k][j][i_dummy].x * csi[k][j][i_face].x) + (ubcs[k][j][i_dummy].y * csi[k][j][i_face].y) + (ubcs[k][j][i_dummy].z * csi[k][j][i_face].z);
                        
                        PetscReal Cell_Area = sqrt((csi[k][j][i_face].x*csi[k][j][i_face].x) + (csi[k][j][i_face].y*csi[k][j][i_face].y) + (csi[k][j][i_face].z*csi[k][j][i_face].z));
                                                
                        PetscReal Correction_flux = velocity_correction*Cell_Area;

                        ucont[k][j][i_face].x = Uncorrected_local_flux + Correction_flux;
                    }     
                }
            }
            break;        
        }
        case BC_FACE_POS_X:{
            const PetscInt i_cell = xe - 2;
            const PetscInt i_face = xe - 2;
            const PetscInt i_dummy = xe - 1;
            for(PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++){
                if(nvert[k][j][i_cell]<0.1){
                    // Set ubcs
                    ubcs[k][j][i_dummy] = ucat[k][j][i_cell];
                    
                    // Calculate Local uncorrected original flux
                    PetscReal Uncorrected_local_flux = (ubcs[k][j][i_dummy].x * csi[k][j][i_face].x) + (ubcs[k][j][i_dummy].y * csi[k][j][i_face].y) + (ubcs[k][j][i_dummy].z * csi[k][j][i_face].z);
                    
                    PetscReal Cell_Area = sqrt((csi[k][j][i_face].x*csi[k][j][i_face].x) + (csi[k][j][i_face].y*csi[k][j][i_face].y) + (csi[k][j][i_face].z*csi[k][j][i_face].z));
                                            
                    PetscReal Correction_flux = velocity_correction*Cell_Area;

                    ucont[k][j][i_face].x = Uncorrected_local_flux + Correction_flux;
                }
            }
            break;
        }
        case BC_FACE_NEG_Y:{
            const PetscInt j_cell = ys + 1;
            const PetscInt j_face = ys;
            const PetscInt j_dummy = ys;
            for(PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++){
                if(nvert[k][j_cell][i]<0.1){
                    // Set ubcs
                    ubcs[k][j_dummy][i] = ucat[k][j_cell][i];

                    // Calculate Local uncorrected original flux
                    PetscReal Uncorrected_local_flux = (ubcs[k][j_dummy][i].x*eta[k][j_face][i].x) + (ubcs[k][j_dummy][i].y*eta[k][j_face][i].y)  + (ubcs[k][j_dummy][i].z*eta[k][j_face][i].z);
                    
                    PetscReal Cell_Area = sqrt((eta[k][j_face][i].x*eta[k][j_face][i].x)+(eta[k][j_face][i].y*eta[k][j_face][i].y)+(eta[k][j_face][i].z*eta[k][j_face][i].z));

                    PetscReal Correction_flux = velocity_correction*Cell_Area;

                    ucont[k][j_face][i].y = Uncorrected_local_flux + Correction_flux;
                }
            }
            break;
        }
        case BC_FACE_POS_Y:{
            const PetscInt j_cell = ye - 2;
            const PetscInt j_face = ye - 2;
            const PetscInt j_dummy = ye - 1;
            for(PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++){
                if(nvert[k][j_cell][i]<0.1){
                    // Set ubcs
                    ubcs[k][j_dummy][i] = ucat[k][j_cell][i];

                    // Calculate Local uncorrected original flux
                    PetscReal Uncorrected_local_flux = (ubcs[k][j_dummy][i].x*eta[k][j_face][i].x) + (ubcs[k][j_dummy][i].y*eta[k][j_face][i].y)  + (ubcs[k][j_dummy][i].z*eta[k][j_face][i].z);
                    
                    PetscReal Cell_Area = sqrt((eta[k][j_face][i].x*eta[k][j_face][i].x)+(eta[k][j_face][i].y*eta[k][j_face][i].y)+(eta[k][j_face][i].z*eta[k][j_face][i].z));

                    PetscReal Correction_flux = velocity_correction*Cell_Area;

                    ucont[k][j_face][i].y = Uncorrected_local_flux + Correction_flux;
                }
            }
            break;
        }
        case BC_FACE_NEG_Z:{
            const PetscInt k_cell = zs + 1;
            const PetscInt k_face = zs;
            const PetscInt k_dummy = zs;
            for(PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++){
                if(nvert[k_cell][j][i]<0.1){
                    // Set ubcs
                    ubcs[k_dummy][j][i] = ucat[k_cell][j][i];

                    // Calculate Local uncorrected original flux
                    PetscReal Uncorrected_local_flux = ((ubcs[k_dummy][j][i].x*zet[k_face][j][i].x) + (ubcs[k_dummy][j][i].y*zet[k_face][j][i].y) + (ubcs[k_dummy][j][i].z*zet[k_face][j][i].z));

                    PetscReal Cell_Area = sqrt((zet[k_face][j][i].x*zet[k_face][j][i].x)+(zet[k_face][j][i].y*zet[k_face][j][i].y)+(zet[k_face][j][i].z*zet[k_face][j][i].z));

                    PetscReal Correction_flux = velocity_correction*Cell_Area;

                    ucont[k_face][j][i].z = Uncorrected_local_flux + Correction_flux;
                }
            }
            break;
        }
        case BC_FACE_POS_Z:{
            const PetscInt k_cell = ze - 2;
            const PetscInt k_face = ze - 2;
            const PetscInt k_dummy = ze - 1;
             for(PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++){
                if(nvert[k_cell][j][i]<0.1){
                    // Set ubcs
                    ubcs[k_dummy][j][i] = ucat[k_cell][j][i];

                    // Calculate Local uncorrected original flux
                    PetscReal Uncorrected_local_flux = ((ubcs[k_dummy][j][i].x*zet[k_face][j][i].x) + (ubcs[k_dummy][j][i].y*zet[k_face][j][i].y) + (ubcs[k_dummy][j][i].z*zet[k_face][j][i].z));

                    PetscReal Cell_Area = sqrt((zet[k_face][j][i].x*zet[k_face][j][i].x)+(zet[k_face][j][i].y*zet[k_face][j][i].y)+(zet[k_face][j][i].z*zet[k_face][j][i].z));

                    PetscReal Correction_flux = velocity_correction*Cell_Area;

                    ucont[k_face][j][i].z = Uncorrected_local_flux + Correction_flux;
                }
            }
            break;           
        }
    }    
    
    // Restore all arrays
    ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda,user->lUcat, (const Cmpnts***)&ucat); CHKERRQ(ierr);  
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Update outlet-conservation state after a completed solver step.
 */
static PetscErrorCode PostStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                                   PetscReal *local_inflow_contribution,
                                                   PetscReal *local_outflow_contribution)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    DMDALocalInfo* info = &user->info;
    PetscBool      can_service;
    
    (void)self;
    (void)local_inflow_contribution;
    
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    const PetscInt IM_nodes_global = user->IM;
    const PetscInt JM_nodes_global = user->JM;
    const PetscInt KM_nodes_global = user->KM;
    ierr = CanRankServiceFace(info, IM_nodes_global, JM_nodes_global, KM_nodes_global, face_id, &can_service); CHKERRQ(ierr);

    if (!can_service) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    // Get arrays (need both ucont and nvert)
    Cmpnts ***ucont;
    PetscReal ***nvert;  // ✅ ADD nvert
    
    ierr = DMDAVecGetArrayRead(user->fda, user->Ucont, (const Cmpnts***)&ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);  // ✅ ADD

    PetscReal local_flux = 0.0;
    
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;
    PetscInt mx = info->mx, my = info->my, mz = info->mz;
    
    PetscInt lxs = xs, lxe = xe, lys = ys, lye = ye, lzs = zs, lze = ze;
    if (xs == 0) lxs = xs + 1;
    if (xe == mx) lxe = xe - 1;
    if (ys == 0) lys = ys + 1;
    if (ye == my) lye = ye - 1;
    if (zs == 0) lzs = zs + 1;
    if (ze == mz) lze = ze - 1;

    // Sum ucont components, skipping solid cells (same indices as PreStep)
    switch (face_id) {
        case BC_FACE_NEG_X: {
            const PetscInt i_cell = xs + 1;  // ✅ Match PreStep
            const PetscInt i_face = xs;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    if (nvert[k][j][i_cell] < 0.1) {  // ✅ Skip solid cells
                        local_flux += ucont[k][j][i_face].x;
                    }
                }
            }
        } break;
        
        case BC_FACE_POS_X: {
            const PetscInt i_cell = xe - 2;  // ✅ Match PreStep
            const PetscInt i_face = xe - 2;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    if (nvert[k][j][i_cell] < 0.1) {  // ✅ Skip solid cells
                        local_flux += ucont[k][j][i_face].x;
                    }
                }
            }
        } break;
        
        case BC_FACE_NEG_Y: {
            const PetscInt j_cell = ys + 1;  // ✅ Match PreStep
            const PetscInt j_face = ys;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if (nvert[k][j_cell][i] < 0.1) {  // ✅ Skip solid cells
                        local_flux += ucont[k][j_face][i].y;
                    }
                }
            }
        } break;
        
        case BC_FACE_POS_Y: {
            const PetscInt j_cell = ye - 2;  // ✅ Match PreStep
            const PetscInt j_face = ye - 2;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if (nvert[k][j_cell][i] < 0.1) {  // ✅ Skip solid cells
                        local_flux += ucont[k][j_face][i].y;
                    }
                }
            }
        } break;
        
        case BC_FACE_NEG_Z: {
            const PetscInt k_cell = zs + 1;  // ✅ Match PreStep
            const PetscInt k_face = zs;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if (nvert[k_cell][j][i] < 0.1) {  // ✅ Skip solid cells
                        local_flux += ucont[k_face][j][i].z;
                    }
                }
            }
        } break;
        
        case BC_FACE_POS_Z: {
            const PetscInt k_cell = ze - 2;  // ✅ Match PreStep
            const PetscInt k_face = ze - 2;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if (nvert[k_cell][j][i] < 0.1) {  // ✅ Skip solid cells
                        local_flux += ucont[k_face][j][i].z;
                    }
                }
            }
        } break;
    }

    // Restore arrays
    ierr = DMDAVecRestoreArrayRead(user->fda, user->Ucont, (const Cmpnts***)&ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);  // ✅ ADD

    // Add to accumulator
    *local_outflow_contribution += local_flux;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "PostStep_OutletConservation: Face %d, corrected flux = %.6e\n",
              face_id, local_flux);

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Release the outlet-conservation handler's clip state.
 */
static PetscErrorCode Destroy_OutletConservation(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (self && self->data) {
        OutletConservationData *data = (OutletConservationData*)self->data;
        PetscCall(PetscFree(data->clipped));
        PetscCall(PetscFree(self->data));
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(PETSC_SUCCESS);
}


/**
 * @brief Implementation of \ref Create_PeriodicGeometric().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/BC_Handlers.h`.
 * @see Create_PeriodicGeometric()
 */

PetscErrorCode Create_PeriodicGeometric(BoundaryCondition *bc){
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");
    bc->priority = BC_PRIORITY_WALL;

    // Assign function pointers
    bc->Initialize = NULL; // No initialization needed
    bc->PreStep    = NULL;
    bc->Apply      = NULL;
    bc->PostStep   = NULL;
    bc->UpdateUbcs = NULL;
    bc->Destroy    = NULL; // No private data to destroy

    bc->data = NULL;

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Measure the two volumetric fluxes the driven-flow controller senses.
 *
 * Both periodic driven handlers steer on the same pair of measurements, so they
 * are taken here in a single sweep of `lUcont`:
 *
 *  - `*boundaryFlux` is the flux through the single periodic boundary plane.
 *    It is fast and responsive but noisy, and drives the boundary trim.
 *  - `*planarAverageFlux` is the flux averaged over every cross-sectional plane
 *    in the driven direction. It is stable and inertial, and drives the
 *    momentum source. It is also the quantity `initial_flux` latches at t=0.
 *
 * @param[in]  user              Block context supplying `lUcont`, `lNvert` and `info`.
 * @param[in]  direction         Driven direction, 'X', 'Y' or 'Z'.
 * @param[out] boundaryFlux      Globally reduced flux through the boundary plane.
 * @param[out] planarAverageFlux Globally reduced plane-averaged flux.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode MeasureDrivenFluxes(UserCtx *user, char direction,
                                          PetscReal *boundaryFlux,
                                          PetscReal *planarAverageFlux)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    PetscInt i, j, k;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // --- Get read-only access to necessary field data ---
    Cmpnts ***ucont;
    PetscReal ***nvert;
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcont, (const Cmpnts***)&ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // --- Define local loop bounds ---
    PetscInt lo[3], hi[3];
    GetOwnedInteriorCellBounds(&info, lo, hi);
    PetscInt lxs = lo[0], lys = lo[1], lzs = lo[2], lxe = hi[0], lye = hi[1], lze = hi[2];

    // --- Initialize local accumulators ---
    PetscReal localCurrentBoundaryFlux = 0.0;
    PetscReal localAveragePlanarVolumetricFluxTerm = 0.0;

    // --- Measure local contributions to the two flux types, generalized by direction ---
    switch (direction) {
        case 'X':
            if (info.xs == 0) { // Only the rank on the negative face contributes to boundary flux
                i = 0;
                for (k = lzs; k < lze; k++) for (j = lys; j < lye; j++) {
                    if (nvert[k][j][i + 1] < 0.1) localCurrentBoundaryFlux += ucont[k][j][i].x;
                }
            }
            for (i = info.xs; i < lxe; i++) {
                for (k = lzs; k < lze; k++) for (j = lys; j < lye; j++) {
                    if (nvert[k][j][i + 1] < 0.1) localAveragePlanarVolumetricFluxTerm += ucont[k][j][i].x / (PetscReal)(info.mx - 1);
                }
            }
            break;
        case 'Y':
            if (info.ys == 0) {
                j = 0;
                for (k = lzs; k < lze; k++) for (i = lxs; i < lxe; i++) {
                    if (nvert[k][j + 1][i] < 0.1) localCurrentBoundaryFlux += ucont[k][j][i].y;
                }
            }
            for (j = info.ys; j < lye; j++) {
                for (k = lzs; k < lze; k++) for (i = lxs; i < lxe; i++) {
                    if (nvert[k][j + 1][i] < 0.1) localAveragePlanarVolumetricFluxTerm += ucont[k][j][i].y / (PetscReal)(info.my - 1);
                }
            }
            break;
        case 'Z':
            if (info.zs == 0) {
                k = 0;
                for (j = lys; j < lye; j++) for (i = lxs; i < lxe; i++) {
                    if (nvert[k + 1][j][i] < 0.1) localCurrentBoundaryFlux += ucont[k][j][i].z;
                }
            }
            for (k = info.zs; k < lze; k++) {
                for (j = lys; j < lye; j++) for (i = lxs; i < lxe; i++) {
                    if (nvert[k + 1][j][i] < 0.1) localAveragePlanarVolumetricFluxTerm += ucont[k][j][i].z / (PetscReal)(info.mz - 1);
                }
            }
            break;
        default:
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                    "MeasureDrivenFluxes received an unknown driven direction '%c'.", direction);
    }

    // --- Release array access as soon as possible ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcont, (const Cmpnts***)&ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // --- Perform global reductions to get the final flux values ---
    ierr = MPI_Allreduce(&localCurrentBoundaryFlux, boundaryFlux, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&localAveragePlanarVolumetricFluxTerm, planarAverageFlux, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

// ===============================================================================
//
//      HANDLER IMPLEMENTATION: PERIODIC DRIVEN CONSTANT FLUX
//      (Corresponds to BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX)
//
// ===============================================================================

// --- 1. FORWARD DECLARATIONS & PRIVATE DATA ---

// Forward declarations for the static functions that implement this handler's behavior.
static PetscErrorCode Initialize_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx, PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode Destroy_PeriodicDrivenConstant(BoundaryCondition *self);

/**
 * @brief Private data structure shared by both periodic driven-flux handlers.
 *
 * `constant_flux` fills `targetVolumetricFlux` from the bcs file at
 * initialization; `initial_flux` latches it from the starting field at the
 * first PreStep. Everything downstream of the target is identical, so both
 * handlers reuse this struct and the PreStep/Apply/Destroy implementations
 * below.
 */
typedef struct {
    char      direction;                // 'X', 'Y', or 'Z', determined at initialization.
    PetscReal targetVolumetricFlux;     // The target flux this controller drives to.
    PetscBool isMasterController;       // Flag: PETSC_TRUE only for the handler on the negative face.
    PetscBool enforceSeamFlux;          // Flag: PETSC_TRUE to add the seam-flux correction into Ucont.
    PetscInt  lastBulkCorrectionStep;   // Physical step at which the momentum source was last set (-1 = never).
} DrivenFluxData;


// --- 2. HANDLER CONSTRUCTOR ---

/**
 * @brief Implementation of \ref Create_PeriodicDrivenConstant().
 * @see Create_PeriodicDrivenConstant()
 */
PetscErrorCode Create_PeriodicDrivenConstant(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");

    // --- Allocate the private data structure ---
    DrivenFluxData *data = NULL;
    ierr = PetscNew(&data); CHKERRQ(ierr);
    // Initialize fields to safe default values
    data->direction = ' ';
    data->targetVolumetricFlux = 0.0;
    data->isMasterController = PETSC_FALSE;
    data->enforceSeamFlux = PETSC_FALSE;
    data->lastBulkCorrectionStep = -1;
    
    // Attach the private data to the generic handler object
    bc->data = (void*)data;
    
    // --- Configure the handler's properties and methods ---
    
    // Set priority: Using BC_PRIORITY_INLET ensures this handler's PreStep runs
    // before other handlers (like outlets) that might depend on its calculations.
    // It is the caller's responsibility that there are no Inlets called along with driven periodic to avoid clash.
    bc->priority   = BC_PRIORITY_INLET;
    
    // Assign the function pointers to the implementations in this file.
    bc->Initialize = Initialize_PeriodicDrivenConstant;
    bc->PreStep    = PreStep_PeriodicDrivenConstant;
    bc->Apply      = Apply_PeriodicDrivenConstant;
    bc->PostStep   = NULL; // This handler has no action after the main solver step.
    bc->UpdateUbcs = NULL; // The boundary value is not flow-dependent (it's periodic).
    bc->Destroy    = Destroy_PeriodicDrivenConstant;
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Initialize constant forcing data for a periodically driven boundary.
 */
static PetscErrorCode Initialize_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    DrivenFluxData *data = (DrivenFluxData*)self->data;
    BCFace face_id = ctx->face_id;
    UserCtx* user = ctx->user;
    
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initializing PERIODIC_DRIVEN_CONSTANT_FLUX handler on Face %s...\n", BCFaceToString(face_id));

    // --- 1. Validation: Ensure the mathematical type is PERIODIC ---
    if (user->boundary_faces[face_id].mathematical_type != PERIODIC) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                "Configuration Error: Handler PERIODIC_DRIVEN_CONSTANT_FLUX on Face %s must be applied to a face with mathematical_type PERIODIC.",
                BCFaceToString(face_id));
    }
    
    // --- 2. Role Assignment: Determine direction and master status ---
    data->isMasterController = PETSC_FALSE;
    switch (face_id) {
        case BC_FACE_NEG_X: data->direction = 'X'; data->isMasterController = PETSC_TRUE; break;
        case BC_FACE_POS_X: data->direction = 'X'; break;
        case BC_FACE_NEG_Y: data->direction = 'Y'; data->isMasterController = PETSC_TRUE; break;
        case BC_FACE_POS_Y: data->direction = 'Y'; break;
        case BC_FACE_NEG_Z: data->direction = 'Z'; data->isMasterController = PETSC_TRUE; break;
        case BC_FACE_POS_Z: data->direction = 'Z'; break;
    }
    
    // --- 3. Parameter Parsing (Master Controller only) ---
    if (data->isMasterController) {
        PetscBool found;

        // Attempt to read the 'target_flux' parameter from the bcs.run file.
        ierr = GetBCParamReal(user->boundary_faces[face_id].params, "target_flux",
                              &data->targetVolumetricFlux, &found); CHKERRQ(ierr);                      
                
        // If the required parameter is not found, halt with an informative error.
        if (!found) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                    "Configuration Error: Handler PERIODIC_DRIVEN_CONSTANT_FLUX on Face %s requires a 'target_flux' parameter in the bcs file (e.g., target_flux=10.0).",
                    BCFaceToString(face_id));
        }
        
        LOG_ALLOW(GLOBAL, LOG_INFO, "Driven Flow (Dir %c): Constant target volumetric flux set to %le.\n",
                  data->direction, data->targetVolumetricFlux);

        // Store the target flux in the simulation context, where the driven-flow
        // controller and its diagnostics read it.
        user->simCtx->targetVolumetricFlux = data->targetVolumetricFlux;
        // The target is fixed for the run from here on; see the flag's comment in
        // variables.h for how initial_flux uses the same latch.
        user->simCtx->drivenFluxTargetLatched = PETSC_TRUE;
    }

    PetscBool trimfound;
    // Optional seam-flux enforcement; accepts the deprecated `apply_trim` spelling.
    ierr = GetDrivenSeamFluxFlag(user->boundary_faces[face_id].params,
                                 &data->enforceSeamFlux, &trimfound); CHKERRQ(ierr);

    if(!trimfound) LOG_ALLOW(GLOBAL,LOG_DEBUG,"Seam-flux enforcement not specified, defaults to %s.\n",data->enforceSeamFlux? "True":"False");
    // Enforcing the seam flux adds a trim re-measured every pass, so the pass must repeat.
    self->needs_repeat_pass = data->enforceSeamFlux;

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Prepare constant periodic-driving data before the solver step.
 */
static PetscErrorCode PreStep_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx,
                                                     PetscReal *local_inflow_contribution,
                                                     PetscReal *local_outflow_contribution)
{
    PetscErrorCode ierr;
    DrivenFluxData *data = (DrivenFluxData*)self->data;
    UserCtx* user = ctx->user;
    SimCtx* simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // --- Master Check: Only the handler on the negative face performs calculations ---
    if (!data->isMasterController) {
        PROFILE_HELPER_END;
        PetscFunctionReturn(0);
    }
    
    // The controller senses two fluxes; see MeasureDrivenFluxes() for what each
    // one is for and why the controller needs both.
    char direction = data->direction;
    PetscReal globalCurrentBoundaryFlux, globalAveragePlanarVolumetricFlux;
    ierr = MeasureDrivenFluxes(user, direction, &globalCurrentBoundaryFlux,
                               &globalAveragePlanarVolumetricFlux); CHKERRQ(ierr);

    // --- Get cross-sectional area using the dedicated geometry function ---
    Cmpnts ignored_center;
    PetscReal globalBoundaryArea;
    BCFace neg_face_id = (direction == 'X') ? BC_FACE_NEG_X : (direction == 'Y') ? BC_FACE_NEG_Y : BC_FACE_NEG_Z;
    ierr = CalculateFaceCenterAndArea(user, neg_face_id, &ignored_center, &globalBoundaryArea); CHKERRQ(ierr);

    // --- Calculate the two correction terms ---
    //
    // These are refreshed on DIFFERENT cadences, and the difference matters.
    //
    // ApplyBoundaryConditions() runs BoundarySystem_ExecuteStep() -- and hence
    // this PreStep -- three times per call, and it is itself called once per
    // Jameson RK stage under the Picard solver and once per residual evaluation
    // under Newton-Krylov. So "per PreStep" is emphatically not "per timestep".
    //
    //  - bulkVelocityCorrection scales the momentum source in ComputeRHS. Both
    //    momentum solvers are built on it being FROZEN across a timestep: the
    //    Picard shadow-Jacobian estimate treats the body force as a constant
    //    forcing with zero velocity Jacobian, and the Newton solve needs a
    //    source that does not drift between residual evaluations. It is
    //    therefore computed once per physical step, from the field at the start
    //    of that step, and held.
    //
    //  - boundaryVelocityCorrection is the tactical trim applied to the
    //    boundary fluxes in Apply(). It is deliberately re-measured on every
    //    pass: Apply() accumulates it into Ucont, and re-measuring is what makes
    //    that accumulation self-limiting as the seam converges. Freezing it
    //    would make repeated Apply() calls add the same trim over and over.
    if (globalBoundaryArea > 1.0e-12) {
        if (data->lastBulkCorrectionStep != simCtx->step) {
            simCtx->bulkVelocityCorrection = (data->targetVolumetricFlux - globalAveragePlanarVolumetricFlux) / globalBoundaryArea;
            simCtx->drivenFluxMeasured = globalAveragePlanarVolumetricFlux;
            simCtx->drivenFluxArea     = globalBoundaryArea;
            data->lastBulkCorrectionStep = simCtx->step;
        }
        simCtx->boundaryVelocityCorrection = (data->targetVolumetricFlux - globalCurrentBoundaryFlux) / globalBoundaryArea;
    } else {
        simCtx->bulkVelocityCorrection = 0.0;
        simCtx->boundaryVelocityCorrection = 0.0;
        data->lastBulkCorrectionStep = simCtx->step;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Driven Flow Controller Update (Dir %c):\n", data->direction);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Target Volumetric Flux:          %.6e\n", data->targetVolumetricFlux);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Avg Planar Volumetric Flux (Stable): %.6e\n", globalAveragePlanarVolumetricFlux);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Boundary Flux (Fast):            %.6e\n", globalCurrentBoundaryFlux);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Bulk Velocity Correction:     %.6e (For Momentum Source)\n", simCtx->bulkVelocityCorrection);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Boundary Velocity Correction: %.6e (For Boundary Trim)\n", simCtx->boundaryVelocityCorrection);
    
    // Suppress unused parameter warnings for this handler
    (void)local_inflow_contribution;
    (void)local_outflow_contribution;
    
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Apply the configured constant driving term to the periodic boundary.
 */
static PetscErrorCode Apply_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    DrivenFluxData *data = (DrivenFluxData*)self->data;
    UserCtx* user = ctx->user;
    BCFace face_id = ctx->face_id;
    PetscBool can_service;
    
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    
    // Check if this rank owns part of this boundary face
    ierr = CanRankServiceFace(&user->info, user->IM, user->JM, user->KM, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PROFILE_HELPER_END;
        PetscFunctionReturn(0);
    }
    
    // If the correction is negligible, no work is needed.
    if (PetscAbsReal(user->simCtx->boundaryVelocityCorrection) < 1e-12) {
        PROFILE_HELPER_END;
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(LOCAL, LOG_TRACE, "Apply_PeriodicDrivenConstant: Applying boundary trim on Face %s...\n", BCFaceToString(face_id));

    // --- Get read/write access to necessary arrays ---
    DMDALocalInfo info = user->info;
    Cmpnts ***ucont, ***uch, ***csi, ***eta, ***zet;
    PetscReal ***nvert;
    
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Bcs.Uch, &uch); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);
    
    PetscInt lo[3], hi[3];
    GetOwnedInteriorCellBounds(&info, lo, hi);
    PetscInt lxs = lo[0], lys = lo[1], lzs = lo[2], lxe = hi[0], lye = hi[1], lze = hi[2];

    // --- Apply correction to the appropriate face and velocity component ---
    switch (face_id) {
        case BC_FACE_NEG_X: case BC_FACE_POS_X: {
            PetscInt i_face = (face_id == BC_FACE_NEG_X) ? info.xs : info.mx - 2;
            PetscInt i_nvert = (face_id == BC_FACE_NEG_X) ? info.xs + 1 : info.mx - 2;

            for (PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++) {
                if (nvert[k][j][i_nvert] < 0.1) {
                    PetscReal faceArea = sqrt(csi[k][j][i_face].x*csi[k][j][i_face].x + csi[k][j][i_face].y*csi[k][j][i_face].y + csi[k][j][i_face].z*csi[k][j][i_face].z);
                    PetscReal fluxTrim = user->simCtx->boundaryVelocityCorrection * faceArea;
                    if(data->enforceSeamFlux) ucont[k][j][i_face].x += fluxTrim;
                    uch[k][j][i_face].x = fluxTrim; // Store correction for diagnostics
                }
            }
        } break;
        
        case BC_FACE_NEG_Y: case BC_FACE_POS_Y: {
            PetscInt j_face = (face_id == BC_FACE_NEG_Y) ? info.ys : info.my - 2;
            PetscInt j_nvert = (face_id == BC_FACE_NEG_Y) ? info.ys + 1 : info.my - 2;
            
            for (PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++) {
                if (nvert[k][j_nvert][i] < 0.1) {
                    PetscReal faceArea = sqrt(eta[k][j_face][i].x*eta[k][j_face][i].x + eta[k][j_face][i].y*eta[k][j_face][i].y + eta[k][j_face][i].z*eta[k][j_face][i].z);
                    PetscReal fluxTrim = user->simCtx->boundaryVelocityCorrection * faceArea;
                    if(data->enforceSeamFlux) ucont[k][j_face][i].y += fluxTrim;
                    uch[k][j_face][i].y = fluxTrim;
                }
            }
        } break;
            
        case BC_FACE_NEG_Z: case BC_FACE_POS_Z: {
            PetscInt k_face = (face_id == BC_FACE_NEG_Z) ? info.zs : info.mz - 2;
            PetscInt k_nvert = (face_id == BC_FACE_NEG_Z) ? info.zs + 1 : info.mz - 2;
            
            for (PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++) {
                if (nvert[k_nvert][j][i] < 0.1) {
                    PetscReal faceArea = sqrt(zet[k_face][j][i].x*zet[k_face][j][i].x + zet[k_face][j][i].y*zet[k_face][j][i].y + zet[k_face][j][i].z*zet[k_face][j][i].z);
                    PetscReal fluxTrim = user->simCtx->boundaryVelocityCorrection * faceArea;
                    if(data->enforceSeamFlux) ucont[k_face][j][i].z += fluxTrim;
                    uch[k_face][j][i].z = fluxTrim;
                }
            }
        } break;
    }

    // --- Restore arrays ---
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Uch, &uch); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);
    
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Release resources owned by the constant periodic-driving boundary.
 */
static PetscErrorCode Destroy_PeriodicDrivenConstant(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // Check that the handler object and its private data pointer are valid before trying to free.
    if (self && self->data) {
        // The private data was allocated with PetscNew(), so it must be freed with PetscFree().
        PetscFree(self->data);
        
        // It is good practice to nullify the pointer after freeing to prevent
        // any accidental use of the dangling pointer (use-after-free).
        self->data = NULL;
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "Destroy_PeriodicDrivenConstant: Private data freed successfully.\n");
    }

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}


// ===============================================================================
//
//      HANDLER IMPLEMENTATION: PERIODIC DRIVEN INITIAL FLUX
//      (Corresponds to BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX)
//
//      Identical to the CONSTANT_FLUX handler except for where the target comes
//      from: this one measures the volumetric flux of the field the run starts
//      with and then holds it, so it takes no `target_flux` parameter. Once the
//      target is latched the two handlers behave the same, so PreStep delegates
//      to the constant implementation and Apply/Destroy are shared outright.
//
// ===============================================================================

// --- 1. FORWARD DECLARATIONS ---

static PetscErrorCode Initialize_PeriodicDrivenInitial(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_PeriodicDrivenInitial(BoundaryCondition *self, BCContext *ctx, PetscReal *in, PetscReal *out);

// --- 2. HANDLER CONSTRUCTOR ---

/**
 * @brief Implementation of \ref Create_PeriodicDrivenInitial().
 * @see Create_PeriodicDrivenInitial()
 */
PetscErrorCode Create_PeriodicDrivenInitial(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    PetscCheck(bc, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");

    // --- Allocate the private data structure ---
    DrivenFluxData *data = NULL;
    ierr = PetscNew(&data); CHKERRQ(ierr);
    // Initialize fields to safe default values
    data->direction = ' ';
    data->targetVolumetricFlux = 0.0;
    data->isMasterController = PETSC_FALSE;
    data->enforceSeamFlux = PETSC_FALSE;
    data->lastBulkCorrectionStep = -1;

    // Attach the private data to the generic handler object
    bc->data = (void*)data;

    // --- Configure the handler's properties and methods ---

    // Same priority reasoning as the constant-flux handler: PreStep must run
    // before any handler that depends on the controller's corrections.
    bc->priority   = BC_PRIORITY_INLET;

    bc->Initialize = Initialize_PeriodicDrivenInitial;
    bc->PreStep    = PreStep_PeriodicDrivenInitial;
    bc->Apply      = Apply_PeriodicDrivenConstant;  // Boundary trim is target-agnostic.
    bc->PostStep   = NULL; // This handler has no action after the main solver step.
    bc->UpdateUbcs = NULL; // The boundary value is not flow-dependent (it's periodic).
    bc->Destroy    = Destroy_PeriodicDrivenConstant; // Same private data layout.

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Initialize forcing data for a periodic boundary driven to its initial flux.
 *
 * @note The target itself is NOT measured here. Boundary handlers are
 *       initialized before `InitializeEulerianState()` runs, so at this point
 *       `Ucont` still holds zeros. The measurement is deferred to the first
 *       PreStep, by which time either the initial condition has been applied or
 *       the restart target has been restored from the checkpoint manifest.
 */
static PetscErrorCode Initialize_PeriodicDrivenInitial(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    DrivenFluxData *data = (DrivenFluxData*)self->data;
    BCFace face_id = ctx->face_id;
    UserCtx* user = ctx->user;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initializing PERIODIC_DRIVEN_INITIAL_FLUX handler on Face %s...\n", BCFaceToString(face_id));

    // --- 1. Validation: Ensure the mathematical type is PERIODIC ---
    if (user->boundary_faces[face_id].mathematical_type != PERIODIC) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                "Configuration Error: Handler PERIODIC_DRIVEN_INITIAL_FLUX on Face %s must be applied to a face with mathematical_type PERIODIC.",
                BCFaceToString(face_id));
    }

    // --- 2. Role Assignment: Determine direction and master status ---
    data->isMasterController = PETSC_FALSE;
    switch (face_id) {
        case BC_FACE_NEG_X: data->direction = 'X'; data->isMasterController = PETSC_TRUE; break;
        case BC_FACE_POS_X: data->direction = 'X'; break;
        case BC_FACE_NEG_Y: data->direction = 'Y'; data->isMasterController = PETSC_TRUE; break;
        case BC_FACE_POS_Y: data->direction = 'Y'; break;
        case BC_FACE_NEG_Z: data->direction = 'Z'; data->isMasterController = PETSC_TRUE; break;
        case BC_FACE_POS_Z: data->direction = 'Z'; break;
    }

    // --- 3. Parameter Parsing (Master Controller only) ---
    if (data->isMasterController) {
        PetscReal unused_flux;
        PetscBool found;

        // This handler derives its own target, so an explicit one is a config error
        // rather than something to silently ignore. Users who want to prescribe the
        // flux should select the `constant_flux` handler instead.
        ierr = GetBCParamReal(user->boundary_faces[face_id].params, "target_flux",
                              &unused_flux, &found); CHKERRQ(ierr);
        if (found) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                    "Configuration Error: Handler PERIODIC_DRIVEN_INITIAL_FLUX on Face %s takes no 'target_flux' parameter; "
                    "it measures the flux of the initial condition and holds that. Use handler 'constant_flux' to prescribe a target.",
                    BCFaceToString(face_id));
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "Driven Flow (Dir %c): target volumetric flux will be latched from the initial state.\n",
                  data->direction);
    }

    PetscBool trimfound;
    // Optional seam-flux enforcement; accepts the deprecated `apply_trim` spelling.
    ierr = GetDrivenSeamFluxFlag(user->boundary_faces[face_id].params,
                                 &data->enforceSeamFlux, &trimfound); CHKERRQ(ierr);

    if(!trimfound) LOG_ALLOW(GLOBAL,LOG_DEBUG,"Seam-flux enforcement not specified, defaults to %s.\n",data->enforceSeamFlux? "True":"False");
    // Enforcing the seam flux adds a trim re-measured every pass, so the pass must repeat.
    self->needs_repeat_pass = data->enforceSeamFlux;

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Latch the initial-state flux once, then drive to it like a constant target.
 *
 * @details The latch is one-shot and guarded by `simCtx->drivenFluxTargetLatched`:
 *          - Fresh start: the flag is false and `Ucont` now holds the initial
 *            condition, so the plane-averaged flux is measured and stored.
 *          - Restart: `ReadSimulationFields()` restored both the target and the
 *            flag from the checkpoint manifest, so the original target survives
 *            instead of being re-measured from a drifted field.
 */
static PetscErrorCode PreStep_PeriodicDrivenInitial(BoundaryCondition *self, BCContext *ctx,
                                                    PetscReal *local_inflow_contribution,
                                                    PetscReal *local_outflow_contribution)
{
    PetscErrorCode ierr;
    DrivenFluxData *data = (DrivenFluxData*)self->data;
    SimCtx* simCtx = ctx->user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    if (data->isMasterController) {
        if (!simCtx->drivenFluxTargetLatched) {
            PetscReal boundaryFlux, planarAverageFlux;

            /* Boundary handlers are also exercised once during setup, from
             * FinalizeBlockState(), and at that point the initial condition has
             * been written to Ucont but not yet scattered into the ghosted
             * lUcont that MeasureDrivenFluxes() reads. Latching there would
             * record a target of zero. Setup runs with step == StartStep, so
             * wait for the first PreStep of the first real timestep: by then
             * lUcont holds the initial condition (or, on a restart, the state
             * the restored target already describes). Do no work at all until
             * then, so no bogus correction is derived from a zero target. */
            if (simCtx->step <= simCtx->StartStep) {
                simCtx->bulkVelocityCorrection = 0.0;
                simCtx->boundaryVelocityCorrection = 0.0;
                PROFILE_HELPER_END;
                PetscFunctionReturn(0);
            }

            ierr = MeasureDrivenFluxes(ctx->user, data->direction,
                                       &boundaryFlux, &planarAverageFlux); CHKERRQ(ierr);

            simCtx->targetVolumetricFlux    = planarAverageFlux;
            simCtx->drivenFluxTargetLatched = PETSC_TRUE;

            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "Driven Flow (Dir %c): latched initial volumetric flux %.6e as the target.\n",
                      data->direction, planarAverageFlux);
        }
        // Keep the handler's copy in step with the authoritative value in SimCtx,
        // which is also where a restart deposits the restored target.
        data->targetVolumetricFlux = simCtx->targetVolumetricFlux;
    }

    ierr = PreStep_PeriodicDrivenConstant(self, ctx,
                                          local_inflow_contribution,
                                          local_outflow_contribution); CHKERRQ(ierr);

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
