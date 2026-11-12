#include "Boundaries.h" // The main header for our project
#include <string.h>    // For strcasecmp
#include <ctype.h>     // For isspace

/**
 * @brief Implementation of \ref BoundaryFaceHasWallShear().
 * @see BoundaryFaceHasWallShear()
 */
PetscBool BoundaryFaceHasWallShear(const BoundaryFaceConfig *face)
{
    return (PetscBool)(face != NULL && face->mathematical_type == WALL &&
                       face->handler_type != BC_HANDLER_WALL_FREE_SLIP);
}


/**
 * @brief Implementation of \ref CanRankServiceFace().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see CanRankServiceFace()
 */
PetscErrorCode CanRankServiceFace(const DMDALocalInfo *info, PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
                                  BCFace face_id, PetscBool *can_service_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank_for_logging;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    rank_for_logging = PetscGlobalRank;

    *can_service_out = PETSC_FALSE; // Default to no service

    // Get the range of cells owned by this rank
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;
    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Determine the global index of the last cell (0-indexed) in each direction.
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;

    switch (face_id) {
        case BC_FACE_NEG_X:
            if (info->xs == 0 && num_owned_cells_on_rank_i > 0 &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_X:
            if (last_global_cell_idx_i >= 0 &&
                (owned_start_cell_i + num_owned_cells_on_rank_i - 1) == last_global_cell_idx_i &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Y:
            if (info->ys == 0 && num_owned_cells_on_rank_j > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Y:
            if (last_global_cell_idx_j >= 0 &&
                (owned_start_cell_j + num_owned_cells_on_rank_j - 1) == last_global_cell_idx_j &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Z:
            if (info->zs == 0 && num_owned_cells_on_rank_k > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Z:
            if (last_global_cell_idx_k >= 0 &&
                (owned_start_cell_k + num_owned_cells_on_rank_k - 1) == last_global_cell_idx_k &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        default:
             LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Unknown face enum %d. \n", rank_for_logging, face_id);
            break;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d check for face %s: Result=%s. \n",
        rank_for_logging, BCFaceToString((BCFace)face_id), (*can_service_out ? "TRUE" : "FALSE"));

    PROFILE_FUNCTION_END;
        
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref GetBoundaryFaceElements().
 * @see GetBoundaryFaceElements()
 */
PetscErrorCode GetBoundaryFaceElements(UserCtx *user, BCFace face_id, BoundaryFaceElements *elements)
{
    const DMDALocalInfo *info = &user->info;
    const PetscInt       nodes[3] = {info->mx, info->my, info->mz};
    const PetscBool      negative = (PetscBool)((PetscInt)face_id % 2 == 0);
    PetscBool            can_service = PETSC_FALSE;
    PetscInt             lo[3], hi[3];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    elements->normal = (PetscInt)face_id / 2;
    elements->a = (elements->normal == 0) ? 1 : 0;
    elements->b = (elements->normal == 2) ? 1 : 2;
    elements->layer = negative ? 0 : nodes[elements->normal] - 2;
    PetscCall(CanRankServiceFace(info, user->IM, user->JM, user->KM, face_id, &can_service));
    GetOwnedInteriorCellBounds(info, lo, hi);
    elements->lo[0] = lo[elements->a];
    elements->lo[1] = lo[elements->b];
    elements->n[0] = can_service ? hi[elements->a] - lo[elements->a] : 0;
    elements->n[1] = can_service ? hi[elements->b] - lo[elements->b] : 0;
    PROFILE_HELPER_END;
    PetscFunctionReturn(PETSC_SUCCESS);
}


/**
 * @brief Implementation of \ref ClassifyMomentumRow().
 *
 * Pure index/boundary-type arithmetic: no field reads, no communication.
 * Precedence matters. A conditioned row is reported first because its explicit
 * Dirichlet value is more specific than the homogeneous fallback, and a periodic
 * duplicate is reported before the homogeneous case because the Newton path needs
 * its representative index to build `F = X_dup - X_rep`.
 */
MomentumRowType ClassifyMomentumRow(UserCtx *user, PetscInt i, PetscInt j, PetscInt k,
                                    PetscInt component, PetscInt *ri, PetscInt *rj, PetscInt *rk)
{
    const PetscInt mx = user->info.mx, my = user->info.my, mz = user->info.mz;
    const PetscInt coord[3] = {i, j, k};
    const PetscInt size[3] = {mx, my, mz};
    PetscBool periodic[3], periodic_duplicate = PETSC_FALSE;
    PetscBool residual_zeroed = PETSC_FALSE, conditioned = PETSC_FALSE;

    *ri = i; *rj = j; *rk = k;
    GetPeriodicAxes(user->simCtx, periodic);
    for (PetscInt axis = 0; axis < 3; ++axis) {
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

    if (conditioned) return MOM_ROW_FIXED_CONDITIONED;
    if (periodic_duplicate) return MOM_ROW_PERIODIC_DUPLICATE;
    if (residual_zeroed) return MOM_ROW_FIXED_HOMOGENEOUS;
    return MOM_ROW_PHYSICAL;
}

/**
 * @brief Implementation of \ref MomentumRowIsSolidMasked().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 * @see MomentumRowIsSolidMasked()
 */
PetscBool MomentumRowIsSolidMasked(const PetscReal ***nvert, PetscInt i, PetscInt j,
                                   PetscInt k, PetscInt component)
{
    PROFILE_HELPER_BEGIN;
    /* PICURV_SOLID_THRESHOLD mirrors the 0.1 fluid test ComputeRHS() applies. */
    const PetscReal threshold = 0.1;

    if (nvert == NULL) { PROFILE_HELPER_END; return PETSC_FALSE; }
    /* A solid cell carries no momentum equation in any component. */
    if (nvert[k][j][i] > threshold) { PROFILE_HELPER_END; return PETSC_TRUE; }
    /* A staggered row also loses its equation when the cell it points into is solid. */
    switch (component) {
    case 0: PROFILE_HELPER_END; return (PetscBool)(nvert[k][j][i + 1] > threshold);
    case 1: PROFILE_HELPER_END; return (PetscBool)(nvert[k][j + 1][i] > threshold);
    case 2: PROFILE_HELPER_END; return (PetscBool)(nvert[k + 1][j][i] > threshold);
    default: PROFILE_HELPER_END; return PETSC_FALSE;
    }
    PROFILE_HELPER_END;
}

/**
 * @brief Implementation of \ref EnforceRHSBoundaryConditions().
 *
 * The sweep is deliberately expressed over every owned location rather than over
 * the six boundary slabs: restating "which indices can be non-physical" here is
 * exactly the duplication that let the periodic duplicate column go unzeroed.
 * ClassifyMomentumRow() is a handful of integer comparisons and the walk is a
 * single pass with no stencil access, which is negligible next to the several
 * ghosted stencil passes ComputeRHS() has already made over the same range.
 * MomentumNewtonKrylov_ApplyConstraints() walks the same range the same way.
 */
PetscErrorCode EnforceRHSBoundaryConditions(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info = user->info;
    Cmpnts         ***rhs;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Get a writable pointer to the local data of the global RHS vector.
    ierr = DMDAVecGetArray(user->fda, user->Rhs, &rhs); CHKERRQ(ierr);

    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                PetscScalar *row = &rhs[k][j][i].x;   /* .x/.y/.z are contiguous */
                for (PetscInt component = 0; component < 3; component++) {
                    PetscInt ri, rj, rk;
                    if (ClassifyMomentumRow(user, i, j, k, component, &ri, &rj, &rk) != MOM_ROW_PHYSICAL)
                        row[component] = 0.0;
                }
            }
        }
    }

    // --- Release the pointer to the local data ---
    ierr = DMDAVecRestoreArray(user->fda, user->Rhs, &rhs); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_TRACE, "Rank %d, Block %" PetscInt_FMT ": Finished enforcing RHS boundary conditions.\n",
              user->simCtx->rank, user->_this);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref BoundaryCondition_Create().
 * @see BoundaryCondition_Create()
 */

PetscErrorCode BoundaryCondition_Create(BCHandlerType handler_type, BoundaryCondition **new_bc_ptr)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    const char* handler_name = BCHandlerTypeToString(handler_type);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Factory called for handler type %s. \n", handler_name);

    ierr = PetscMalloc1(1, new_bc_ptr); CHKERRQ(ierr);
    BoundaryCondition *bc = *new_bc_ptr;

    bc->type        = handler_type;
    bc->priority    = -1;  // Default priority; can be overridden in specific handlers
    bc->data        = NULL;
    bc->Initialize  = NULL;
    bc->PreStep     = NULL;
    bc->Apply       = NULL;
    bc->PostStep    = NULL;
    bc->UpdateUbcs  = NULL;
    bc->Destroy     = NULL;
    bc->needs_repeat_pass = PETSC_FALSE; // A handler whose Apply() reads in-pass state sets it

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Allocated generic handler object at address %p.\n", (void*)bc);

    switch (handler_type) {

        case BC_HANDLER_OUTLET_CONSERVATION:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_OutletConservation().\n");
            ierr = Create_OutletConservation(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_WALL_NOSLIP:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_WallNoSlip().\n");
            ierr = Create_WallNoSlip(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_WALL_FREE_SLIP:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_WallFreeSlip().\n");
            ierr = Create_WallFreeSlip(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_WALL_MOVING:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_WallMoving().\n");
            ierr = Create_WallMoving(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_PERIODIC_GEOMETRIC:
            LOG_ALLOW(LOCAL,LOG_DEBUG,"Dispatching to Create_PeriodicGeometric().\n");
            ierr = Create_PeriodicGeometric(bc);
            break;
        
        case BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX:
            LOG_ALLOW(LOCAL,LOG_DEBUG,"Dispatching to Create_PeriodicDrivenConstant().\n");
            ierr = Create_PeriodicDrivenConstant(bc);
            break;
        
        case BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX:
            LOG_ALLOW(LOCAL,LOG_DEBUG,"Dispatching to Create_PeriodicDrivenInitial().\n");
            ierr = Create_PeriodicDrivenInitial(bc);
            break;
                
        case BC_HANDLER_INLET_PRESCRIBED_VELOCITY:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletPrescribedVelocity().\n");
            ierr = Create_InletPrescribedVelocity(bc); CHKERRQ(ierr);
            break;
        //Add cases for other handlers here in future phases 
        
        default:
            LOG_ALLOW(GLOBAL, LOG_ERROR, "Handler type (%s) is not recognized or implemented in the factory.\n", handler_name);
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Boundary handler type %d (%s) not recognized in factory.\n", handler_type, handler_name);
    }

    if(bc->priority < 0) {
        LOG_ALLOW(GLOBAL, LOG_ERROR, "Handler type %d (%s) did not set a valid priority during creation.\n", handler_type, handler_name);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Boundary handler type %d (%s) did not set a valid priority during creation.\n", handler_type, handler_name);
    }
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Successfully created and configured handler for %s.\n", handler_name);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref BoundarySystem_Validate().
 * @see BoundarySystem_Validate()
 */
PetscErrorCode BoundarySystem_Validate(UserCtx *user)
{
    PetscErrorCode ierr;
    const BCFace neg_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace pos_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};
    const char axis_names[3] = {'X', 'Y', 'Z'};
    DMBoundaryType bx, by, bz;
    PetscBool dm_periodic[3];
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Validating parsed boundary condition configuration...\n");
    bx = user->info.bx; by = user->info.by; bz = user->info.bz;
    dm_periodic[0] = (PetscBool)(bx == DM_BOUNDARY_PERIODIC);
    dm_periodic[1] = (PetscBool)(by == DM_BOUNDARY_PERIODIC);
    dm_periodic[2] = (PetscBool)(bz == DM_BOUNDARY_PERIODIC);

    // --- Rule Set 1: Geometric periodic faces must be paired and match the DM topology. ---
    for (PetscInt axis = 0; axis < 3; axis++) {
        const PetscBool neg_periodic =
            user->boundary_faces[neg_faces[axis]].mathematical_type == PERIODIC;
        const PetscBool pos_periodic =
            user->boundary_faces[pos_faces[axis]].mathematical_type == PERIODIC;

        PetscCheck(neg_periodic == pos_periodic, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                   "Configuration Error: Periodic boundaries in the %c direction must be paired; "
                   "%s is %s while %s is %s.",
                   axis_names[axis],
                   BCFaceToString(neg_faces[axis]), neg_periodic ? "PERIODIC" : "not periodic",
                   BCFaceToString(pos_faces[axis]), pos_periodic ? "PERIODIC" : "not periodic");
        PetscCheck(dm_periodic[axis] == neg_periodic, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                   "Configuration Error: The %c-direction DM periodic flag (%d) does not match "
                   "the paired boundary configuration (%s).",
                   axis_names[axis], (int)dm_periodic[axis], neg_periodic ? "PERIODIC" : "not periodic");
    }

    // --- Rule Set 2: Driven Flow Handler Consistency ---
    // This specialized validator will check all rules related to driven flow handlers.
    ierr = Validate_DrivenFlowConfiguration(user); CHKERRQ(ierr);

    // --- Rule Set 3: (Future Extension) Overset Interface Consistency ---
    // ierr = Validate_OversetConfiguration(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Boundary configuration is valid.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//                       PUBLIC MASTER SETUP FUNCTION
//
//================================================================================
/**
 * @brief Implementation of \ref BoundarySystem_Initialize().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see BoundarySystem_Initialize()
 */
PetscErrorCode BoundarySystem_Initialize(UserCtx *user, const char *bcs_filename)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting creation and initialization of all boundary handlers.\n");

    // =========================================================================
    // Step 0: Clear any existing boundary handlers (if re-initializing).
    // This ensures no memory leaks if this function is called multiple times.
    // =========================================================================
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        if (face_cfg->handler) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Destroying existing handler on Face %s before re-initialization.\n", BCFaceToString((BCFace)i));
            if (face_cfg->handler->Destroy) {
                ierr = face_cfg->handler->Destroy(face_cfg->handler); CHKERRQ(ierr);
            }
            ierr = PetscFree(face_cfg->handler); CHKERRQ(ierr);
            face_cfg->handler = NULL;
        }
    }
    // =========================================================================

    // Step 0.1: Initiate flux sums to zero
    user->simCtx->FluxInSum = 0.0;
    user->simCtx->FluxOutSum = 0.0;
    user->simCtx->FarFluxInSum = 0.0;
    user->simCtx->FarFluxOutSum = 0.0;
    // =========================================================================

    // Step 1: Parse the configuration file to determine user intent.
    // This function, defined in io.c, populates the configuration enums and parameter
    // lists within the user->boundary_faces array on all MPI ranks.
    ierr = ParseAllBoundaryConditions(user, bcs_filename); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Configuration file '%s' parsed successfully.\n", bcs_filename);

    // Step 1.1: Validate the parsed configuration to ensure there are no Boundary Condition conflicts
    ierr = BoundarySystem_Validate(user); CHKERRQ(ierr);

    // Step 1.2: Resolve each face's particle policy from its particle_policy parameter (remove when absent).
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        const char *policy = NULL;
        PetscBool   found;
        ierr = GetBCParamString(face_cfg->params, "particle_policy", &policy, &found); CHKERRQ(ierr);
        if (!found || strcasecmp(policy, "remove") == 0) {
            face_cfg->particle_policy = PARTICLE_BOUNDARY_REMOVE;
        } else if (strcasecmp(policy, "reflect") == 0 && face_cfg->mathematical_type != PERIODIC) {
            face_cfg->particle_policy = PARTICLE_BOUNDARY_REFLECT;
        } else if (strcasecmp(policy, "wrap") == 0 && face_cfg->mathematical_type == PERIODIC) {
            face_cfg->particle_policy = PARTICLE_BOUNDARY_WRAP;
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Face %s (%s): particle_policy '%s' is not supported; use 'remove' or '%s'.",
                    BCFaceToString((BCFace)i), BCTypeToString(face_cfg->mathematical_type), policy,
                    (face_cfg->mathematical_type == PERIODIC) ? "wrap" : "reflect");
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Face %s particle policy: %s.\n", BCFaceToString((BCFace)i),
                  ParticleBoundaryPolicyToString(face_cfg->particle_policy));
    }

    // Step 2: Create and Initialize the handler object for each of the 6 faces.
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        
        const char *face_name = BCFaceToString(face_cfg->face_id);
        const char *type_name = BCTypeToString(face_cfg->mathematical_type);
        const char *handler_name = BCHandlerTypeToString(face_cfg->handler_type);

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Creating handler for Face %s with Type %s and handler '%s'.\n", face_name, type_name,handler_name);

        // Use the private factory to construct the correct handler object based on the parsed type.
        // The factory returns a pointer to the new handler object, which we store in the config struct.
        ierr = BoundaryCondition_Create(face_cfg->handler_type, &face_cfg->handler); CHKERRQ(ierr);

        // Step 3: Call the specific Initialize() method for the newly created handler.
        // This allows the handler to perform its own setup, like reading parameters from the
        // face_cfg->params list and setting the initial field values on its face.
        if (face_cfg->handler && face_cfg->handler->Initialize) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Calling Initialize() method for handler %s(%s) on Face %s.\n",type_name,handler_name,face_name);
            
            // Prepare the context needed by the Initialize() function.
            BCContext ctx = {
                .user = user,
                .face_id = face_cfg->face_id,
                .global_inflow_sum = &user->simCtx->FluxInSum,  // Global flux sums are not relevant during initialization.
                .global_outflow_sum = &user->simCtx->FluxOutSum,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            ierr = face_cfg->handler->Initialize(face_cfg->handler, &ctx); CHKERRQ(ierr);
        } else {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Handler %s(%s) for Face %s has no Initialize() method, skipping.\n", type_name,handler_name,face_name);
        }
    }
    // =========================================================================
    // NO SYNCHRONIZATION NEEDED HERE
    // =========================================================================
    // Initialize() only reads parameters and allocates memory.
    // It does NOT modify field values (Ucat, Ucont, Ubcs).
    // Field values are set by:
    //   1. Initial conditions (before this function)
    //   2. Apply() during timestepping (after this function)
    // The first call to ApplyBoundaryConditions() will handle synchronization.
    // =========================================================================

    LOG_ALLOW(GLOBAL, LOG_INFO, "All boundary handlers created and initialized successfully.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Implementation of \ref PropagateBoundaryConfigToCoarserLevels().
 * @see PropagateBoundaryConfigToCoarserLevels()
 */
PetscErrorCode PropagateBoundaryConfigToCoarserLevels(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG *usermg = &simCtx->usermg;
    
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Propagating BC configuration from finest to coarser multigrid levels...\n");
    
    // Loop from second-finest down to coarsest
    for (PetscInt level = usermg->mglevels - 2; level >= 0; level--) {
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            UserCtx *user_coarse = &usermg->mgctx[level].user[bi];
            UserCtx *user_fine   = &usermg->mgctx[level + 1].user[bi];
            
            LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Copying BC config from level %" PetscInt_FMT " to level %" PetscInt_FMT ", block %" PetscInt_FMT "\n",
                          simCtx->rank, level + 1, level, bi);
            
            // Copy the 6 boundary face configurations
            for (int face_i = 0; face_i < 6; face_i++) {
                user_coarse->boundary_faces[face_i].face_id = user_fine->boundary_faces[face_i].face_id;
                user_coarse->boundary_faces[face_i].mathematical_type = user_fine->boundary_faces[face_i].mathematical_type;
                user_coarse->boundary_faces[face_i].handler_type = user_fine->boundary_faces[face_i].handler_type;
                user_coarse->boundary_faces[face_i].particle_policy = user_fine->boundary_faces[face_i].particle_policy;
                
                // Copy parameter list (deep copy)
                FreeBC_ParamList(user_coarse->boundary_faces[face_i].params); // Clear any existing
                user_coarse->boundary_faces[face_i].params = NULL;
                
                BC_Param **dst_next = &user_coarse->boundary_faces[face_i].params;
                for (BC_Param *src = user_fine->boundary_faces[face_i].params; src; src = src->next) {
                    BC_Param *new_param;
                    ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                    ierr = PetscStrallocpy(src->key, &new_param->key); CHKERRQ(ierr);
                    ierr = PetscStrallocpy(src->value, &new_param->value); CHKERRQ(ierr);
                    new_param->next = NULL;
                    *dst_next = new_param;
                    dst_next = &new_param->next;
                }
                
                // IMPORTANT: Do NOT create handler objects for coarser levels
                // Handlers are only needed at finest level for timestepping Apply() calls
                user_coarse->boundary_faces[face_i].handler = NULL;
            }
            
            // Propagate the particle inlet lookup fields to coarse levels as well.
            user_coarse->inletFaceDefined = user_fine->inletFaceDefined;
            user_coarse->identifiedInletBCFace = user_fine->identifiedInletBCFace;
        }
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "BC configuration propagation complete.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//                      PUBLIC MASTER TIME-STEP FUNCTION
//
//================================================================================

/**
 * @brief Implementation of \ref BoundarySystem_ExecuteStep().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see BoundarySystem_ExecuteStep()
 */
PetscErrorCode BoundarySystem_ExecuteStep(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Starting.\n");
    
    // =========================================================================
    // PRIORITY 0: INLETS
    // =========================================================================
    
        PetscReal local_inflow_pre = 0.0;
        PetscReal local_inflow_post = 0.0;
        PetscReal global_inflow_pre = 0.0;
        PetscReal global_inflow_post = 0.0;
        PetscInt  num_handlers[3] = {0,0,0};
        
        LOG_ALLOW(LOCAL, LOG_TRACE, " (INLETS): Begin.\n");
        
        // Phase 1: PreStep - Preparation (e.g., calculate profiles, read files)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, &local_inflow_pre, NULL); CHKERRQ(ierr);
        }
        
        // Optional: Global communication for PreStep (for debugging)
        if (local_inflow_pre != 0.0) {
            ierr = MPI_Allreduce(&local_inflow_pre, &global_inflow_pre, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            LOG_ALLOW(GLOBAL, LOG_TRACE, "    PreStep predicted flux: %.6e\n", global_inflow_pre);
        }
        
        // Phase 2: Apply - Set boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            if(!handler->Apply) continue; // For example Periodic BCs  
            
            num_handlers[1]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 3: PostStep - Measure actual flux
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            if (!handler->PostStep) continue;
            
            num_handlers[2]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, &local_inflow_post, NULL); CHKERRQ(ierr);
        }
        
        // Phase 4: Global communication - Sum flux for other priorities to use
        ierr = MPI_Allreduce(&local_inflow_post, &global_inflow_post, 1, MPIU_REAL,
                            MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // Store for next priority levels
        user->simCtx->FluxInSum = global_inflow_post;
        
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "  (INLETS): %" PetscInt_FMT " Prestep(s), %" PetscInt_FMT " Application(s), %" PetscInt_FMT " Poststep(s), FluxInSum = %.6e\n",
                  num_handlers[0],num_handlers[1],num_handlers[2], global_inflow_post);
    
    // =========================================================================
    // PRIORITY 1: FARFIELD
    // =========================================================================
    
        PetscReal local_farfield_in_pre = 0.0;
        PetscReal local_farfield_out_pre = 0.0;
        PetscReal local_farfield_in_post = 0.0;
        PetscReal local_farfield_out_post = 0.0;
        PetscReal global_farfield_in_pre = 0.0;
        PetscReal global_farfield_out_pre = 0.0;
        PetscReal global_farfield_in_post = 0.0;
        PetscReal global_farfield_out_post = 0.0;
        memset(num_handlers,0,sizeof(num_handlers));
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  (FARFIELD): Begin.\n");
        
        // Phase 1: PreStep - Analyze flow direction, measure initial flux
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,  // Available from Priority 0
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, &local_farfield_in_pre, &local_farfield_out_pre);
            CHKERRQ(ierr);
        }
        
        // Phase 2: Global communication (optional, for debugging)
        if (local_farfield_in_pre != 0.0 || local_farfield_out_pre != 0.0) {
            ierr = MPI_Allreduce(&local_farfield_in_pre, &global_farfield_in_pre, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Allreduce(&local_farfield_out_pre, &global_farfield_out_pre, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            
            LOG_ALLOW(GLOBAL, LOG_DEBUG, 
                      "    Farfield pre-analysis: In=%.6e, Out=%.6e\n",
                      global_farfield_in_pre, global_farfield_out_pre);
        }
        
        // Phase 3: Apply - Set farfield boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            if(!handler->Apply) continue; // For example Periodic BCs            
            
            num_handlers[1]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 4: PostStep - Measure actual farfield fluxes
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            if (!handler->PostStep) continue;
            
            num_handlers[2]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, &local_farfield_in_post, &local_farfield_out_post);
            CHKERRQ(ierr);
        }
        
        // Phase 5: Global communication - Store for outlet priority
        if (num_handlers > 0) {
            ierr = MPI_Allreduce(&local_farfield_in_post, &global_farfield_in_post, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Allreduce(&local_farfield_out_post, &global_farfield_out_post, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            
            // Store for outlet handlers to use
            user->simCtx->FarFluxInSum = global_farfield_in_post;
            user->simCtx->FarFluxOutSum = global_farfield_out_post;
            
            LOG_ALLOW(GLOBAL, LOG_INFO, 
                      "  (FARFIELD): %" PetscInt_FMT " Prestep(s), %" PetscInt_FMT " Application(s), %" PetscInt_FMT " Poststep(s) , InFlux=%.6e, OutFlux=%.6e\n",
                      num_handlers[0],num_handlers[1],num_handlers[2], global_farfield_in_post, global_farfield_out_post);
        } else {
            // No farfield handlers - zero out the fluxes
            user->simCtx->FarFluxInSum = 0.0;
            user->simCtx->FarFluxOutSum = 0.0;
        }
    
    
    // =========================================================================
    // PRIORITY 2: WALLS
    // =========================================================================
    
        memset(num_handlers,0,sizeof(num_handlers));
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  (WALLS): Begin.\n");
        
        // Phase 1: PreStep - Preparation (usually no-op for walls)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, NULL, NULL); CHKERRQ(ierr);
        }
        
        // No global communication needed for walls
        
        // Phase 2: Apply - Set boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            if(!handler->Apply) continue; // For example Periodic BCs  
            
            num_handlers[1]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 3: PostStep - Post-application processing (usually no-op for walls)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            if (!handler->PostStep) continue;
            
            num_handlers[2]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, NULL, NULL); CHKERRQ(ierr);
        }
        
        // No global communication needed for walls
        
        LOG_ALLOW(GLOBAL, LOG_INFO, "  (WALLS): %" PetscInt_FMT " Prestep(s), %" PetscInt_FMT " Application(s), %" PetscInt_FMT " Poststep(s) applied.\n",
                  num_handlers[0],num_handlers[1],num_handlers[2]);
    
    
    // =========================================================================
    // PRIORITY 3: OUTLETS
    // =========================================================================
    
        PetscReal local_outflow_pre = 0.0;
        PetscReal local_outflow_post = 0.0;
        PetscReal global_outflow_pre = 0.0;
        PetscReal global_outflow_post = 0.0;
        memset(num_handlers,0,sizeof(num_handlers));
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  (OUTLETS): Begin.\n");
        
        // Phase 1: PreStep - Measure uncorrected outflow (from ucat)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,      // From Priority 0
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, NULL, &local_outflow_pre); CHKERRQ(ierr);
        }
        
        // Phase 2: Global communication - Get uncorrected outflow sum
        ierr = MPI_Allreduce(&local_outflow_pre, &global_outflow_pre, 1, MPIU_REAL,
                            MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // Calculate total inflow (inlet + farfield inflow)
        PetscReal total_inflow = user->simCtx->FluxInSum + user->simCtx->FarFluxInSum;
        
        LOG_ALLOW(GLOBAL, LOG_DEBUG, 
                  "    Uncorrected outflow: %.6e, Total inflow: %.6e (Inlet: %.6e + Farfield: %.6e)\n",
                  global_outflow_pre, total_inflow, user->simCtx->FluxInSum, 
                  user->simCtx->FarFluxInSum);
        
        // Phase 3: Apply - Set corrected boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            if(!handler->Apply) continue; // For example Periodic BCs  
            
            num_handlers[1]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,      // From Priority 0
                .global_outflow_sum = &global_outflow_pre, // From PreStep above
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum 
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 4: PostStep - Measure corrected outflow (verification)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            if (!handler->PostStep) continue;
            
            num_handlers[2]++;

            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = &global_outflow_pre,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, NULL, &local_outflow_post); CHKERRQ(ierr);
        }
        
        // Phase 5: Global communication - Verify conservation
        ierr = MPI_Allreduce(&local_outflow_post, &global_outflow_post, 1, MPIU_REAL,
                            MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // Store for global reporting.
        user->simCtx->FluxOutSum = global_outflow_post;
        
        // Conservation check (compare total outflow vs total inflow)
        PetscReal total_outflow = global_outflow_post + user->simCtx->FarFluxOutSum;
        PetscReal flux_error = PetscAbsReal(total_outflow - total_inflow);
        PetscReal relative_error = (total_inflow > 1e-16) ? 
                                   flux_error / total_inflow : flux_error;
        
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "  (OUTLETS): %" PetscInt_FMT " Prestep(s), %" PetscInt_FMT " Application(s), %" PetscInt_FMT " Poststep(s), FluxOutSum = %.6e\n",
                  num_handlers[0],num_handlers[1],num_handlers[2], global_outflow_post);
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "    Conservation: Total In=%.6e, Total Out=%.6e, Error=%.3e (%.2e)%%)\n",
                  total_inflow, total_outflow, flux_error, relative_error * 100.0);
        
        if (relative_error > 1e-6) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, 
                     "    WARNING: Large mass conservation error (%.2e%%)!\n",
                     relative_error * 100.0);
        }
    
    
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "Complete.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

// =============================================================================
//
//                   PRIVATE "LIGHT" EXECUTION ENGINE
//
// =============================================================================

/**
 * @brief Implementation of \ref BoundarySystem_RefreshUbcs().
 * @see BoundarySystem_RefreshUbcs()
 */
PetscErrorCode BoundarySystem_RefreshUbcs(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_TRACE, "Refreshing `ubcs` targets for flow-dependent boundaries...\n");

    // Loop through all 6 faces of the domain
    for (int i = 0; i < 6; i++) {
        BoundaryCondition *handler = user->boundary_faces[i].handler;
        
        // THE FILTER:
        // This is the core logic. We only act if a handler exists for the face
        // AND that handler has explicitly implemented the `UpdateUbcs` method.
        if (handler && handler->UpdateUbcs) {
            
            const char *face_name = BCFaceToString((BCFace)i);
            LOG_ALLOW(LOCAL, LOG_TRACE, "  Calling UpdateUbcs() for handler on Face %s.\n", face_name);

            // Prepare the context. For this refresh step, we don't need to pass flux sums.
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = NULL,
                .global_farfield_outflow_sum = NULL
            };
            
            // Call the handler's specific UpdateUbcs function pointer.
            ierr = handler->UpdateUbcs(handler, &ctx); CHKERRQ(ierr);
        }
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//                         PUBLIC MASTER CLEANUP FUNCTION
//
//================================================================================
/**
 * @brief Implementation of \ref BoundarySystem_Destroy().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see BoundarySystem_Destroy()
 */
PetscErrorCode BoundarySystem_Destroy(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting destruction of all boundary handlers. \n");

    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        const char *face_name = BCFaceToString(face_cfg->face_id);

        // --- Step 1: Free the parameter linked list associated with this face ---
        if (face_cfg->params) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Freeing parameter list for Face %d (%s). \n", i, face_name);
            FreeBC_ParamList(face_cfg->params);
            face_cfg->params = NULL; // Good practice to nullify dangling pointers
        }

        // --- Step 2: Destroy the handler object itself ---
        if (face_cfg->handler) {
            const char *handler_name = BCHandlerTypeToString(face_cfg->handler->type);
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Destroying handler '%s' on Face %d (%s).\n", handler_name, i, face_name);
            
            // Call the handler's specific cleanup function first, if it exists.
            // This will free any memory stored in the handler's private `data` pointer.
            if (face_cfg->handler->Destroy) {
                ierr = face_cfg->handler->Destroy(face_cfg->handler); CHKERRQ(ierr);
            }

            // Finally, free the generic BoundaryCondition object itself.
            ierr = PetscFree(face_cfg->handler); CHKERRQ(ierr);
            face_cfg->handler = NULL;
        }
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Destruction complete.\n");
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Copies one cell field's wrapped local values onto the owned periodic duplicate plane.
 * @details Handles scalar and three-component storage for one selected logical axis.
 */
static PetscErrorCode TransferPeriodicFieldByDirection(UserCtx *user, FieldId field_id, char direction)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    FieldView field_view;
    DM        dm;
    Vec       global_vec;
    Vec       local_vec;
    PetscInt  dof;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(FieldGetView(user, field_id, &field_view));
    PetscCheck((field_view.descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_CELL_SYNC) != 0u,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Field '%s' is not registered for periodic cell synchronization.",
               field_view.descriptor->canonical_name);
    dm = field_view.dm;
    global_vec = field_view.global_vec;
    local_vec = field_view.local_vec;
    dof = field_view.descriptor->dof;

    // --- Execute the copy logic based on DoF and Direction ---
    if (dof == 1) { // --- Handle SCALAR fields (PetscReal) ---
        PetscReal ***g_array, ***l_array;
        ierr = DMDAVecGetArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr); // Use Read for safety

        switch (direction) {
            case 'i':
                if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xs] = l_array[k][j][xs-2];
                if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xe-1] = l_array[k][j][xe+1];
                break;
            case 'j':
                if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ys][i] = l_array[k][ys-2][i];
                if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ye-1][i] = l_array[k][ye+1][i];
                break;
            case 'k':
                if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[zs][j][i] = l_array[zs-2][j][i];
                if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[ze-1][j][i] = l_array[ze+1][j][i];
                break;
            default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid direction '%c'", direction);
        }
        ierr = DMDAVecRestoreArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr);

    } else if (dof == 3) { // --- Handle VECTOR fields (Cmpnts) ---
        Cmpnts ***g_array, ***l_array;
        ierr = DMDAVecGetArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr);

        switch (direction) {
            case 'i':
                if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xs] = l_array[k][j][xs-2];
                if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xe-1] = l_array[k][j][xe+1];
                break;
            case 'j':
                if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ys][i] = l_array[k][ys-2][i];
                if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ye-1][i] = l_array[k][ye+1][i];
                break;
            case 'k':
                if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[zs][j][i] = l_array[zs-2][j][i];
                if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[ze-1][j][i] = l_array[ze+1][j][i];
                break;
            default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid direction '%c'", direction);
        }
        ierr = DMDAVecRestoreArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr);
    }

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Requires a field's local vector to equal a full ghost update, when -picurv_verify_ghost_updates is set.
 * @details Rebuilds the update UpdateLocalGhosts() would make (scatter, then periodic repair)
 *          into a scratch vector and compares every local entry bitwise. The periodic
 *          synchronizers call it where they rely on an update they did not perform: a caller's
 *          locals_current claim and every plane-only refresh.
 * @param[in] user Block context owning the field.
 * @param[in] field_id Field whose local vector is checked.
 * @param[in] where Short description of the claim, for the error message.
 * @return 0 when the check is off or passes; an error naming the field otherwise.
 */
static PetscErrorCode VerifyLocalGhostsCurrent(UserCtx *user, FieldId field_id, const char *where)
{
    FieldView          view;
    Vec                full;
    const PetscScalar *expected, *actual;
    PetscInt           n;
    PetscInt64         differing = 0, total = 0;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user->simCtx->verifyGhostUpdates) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    PetscCall(FieldGetView(user, field_id, &view));
    PetscCall(VecDuplicate(view.local_vec, &full));
    PetscCall(DMGlobalToLocalBegin(view.dm, view.global_vec, INSERT_VALUES, full));
    PetscCall(DMGlobalToLocalEnd(view.dm, view.global_vec, INSERT_VALUES, full));
    PetscCall(RepairPeriodicFieldGhosts(user, &view, full));
    PetscCall(VecGetLocalSize(full, &n));
    PetscCall(VecGetArrayRead(full, &expected));
    PetscCall(VecGetArrayRead(view.local_vec, &actual));
    for (PetscInt entry = 0; entry < n; entry++) {
        if (memcmp(&expected[entry], &actual[entry], sizeof(PetscScalar))) differing++;
    }
    PetscCall(VecRestoreArrayRead(view.local_vec, &actual));
    PetscCall(VecRestoreArrayRead(full, &expected));
    PetscCall(VecDestroy(&full));
    PetscCallMPI(MPI_Allreduce(&differing, &total, 1, MPIU_INT64, MPI_SUM, PetscObjectComm((PetscObject)view.dm)));
    PetscCheck(total == 0, PetscObjectComm((PetscObject)view.dm), PETSC_ERR_PLIB,
               "Ghost-update verification failed for '%s' %s: %" PetscInt64_FMT " local entries differ from a full update.",
               view.descriptor->canonical_name, where, total);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Requires the boundary pass's fields to be finite, when -picurv_verify_ghost_updates is set.
 * @details The companion of PoisonStaleLocalGhosts(): a stale local vector it filled with NaN
 *          and something then read would leave non-finite values in these fields.
 * @param[in] user Block context after ApplyBoundaryConditions().
 * @return 0 when the check is off or passes; an error naming the field otherwise.
 */
static PetscErrorCode VerifyBoundaryPassFinite(UserCtx *user)
{
    const FieldId fields[] = {FIELD_ID_UCAT, FIELD_ID_P, FIELD_ID_UCONT};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user->simCtx->verifyGhostUpdates) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    for (size_t f = 0; f < sizeof(fields) / sizeof(fields[0]); f++) {
        FieldView view;
        PetscReal global_norm, local_norm;
        PetscCall(FieldGetView(user, fields[f], &view));
        PetscCall(VecNorm(view.global_vec, NORM_1, &global_norm));
        PetscCall(VecNorm(view.local_vec, NORM_1, &local_norm));
        PetscCheck(PetscIsNormalReal(global_norm) || global_norm == 0.0, PetscObjectComm((PetscObject)view.dm), PETSC_ERR_PLIB,
                   "Ghost-update verification failed: '%s' is not finite after the boundary pass; a stale local vector was read.",
                   view.descriptor->canonical_name);
        PetscCheck(PetscIsNormalReal(local_norm) || local_norm == 0.0, PetscObjectComm((PetscObject)view.dm), PETSC_ERR_PLIB,
                   "Ghost-update verification failed: local '%s' is not finite after the boundary pass; a stale local vector was left.",
                   view.descriptor->canonical_name);
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Fills a field's local vector with NaN, when -picurv_verify_ghost_updates is set.
 * @details Marks a local vector the boundary pass leaves out of date because a later full
 *          UpdateLocalGhosts() rewrites all of it before any read. If anything reads it in
 *          between, the NaN reaches the solution and the run's output.
 * @param[in] user Block context owning the field.
 * @param[in] field_id Field whose local vector is left stale.
 * @return 0 on success.
 */
static PetscErrorCode PoisonStaleLocalGhosts(UserCtx *user, FieldId field_id)
{
    FieldView view;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user->simCtx->verifyGhostUpdates) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    PetscCall(FieldGetView(user, field_id, &view));
    PetscCall(VecSet(view.local_vec, (PetscScalar)NAN));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief The scatter that refreshes only the local entries mirroring a DM's two end planes in one direction.
 * @details A periodic transfer in direction d writes only the global planes 0 and M_d-1. The
 *          scatter is the DMDA's global-to-local scatter restricted to the leaves whose periodic
 *          (wrapped) index along d is one of those planes: the owners' copies and every ghost
 *          copy, including those across the seam. Built on first use, collectively since every
 *          rank synchronizes the same DMs, and kept on the DM, which owns it.
 * @param[in] dm DMDA of the field being synchronized.
 * @param[in] direction 'i', 'j' or 'k'.
 * @param[out] plane The scatter, or NULL when the DM is not periodic along direction.
 * @return 0 on success.
 */
static PetscErrorCode GetPeriodicPlaneScatter(DM dm, char direction, VecScatter *plane)
{
    const PetscInt     axis = direction - 'i';
    char               key[] = "PICurvPeriodicPlaneScatter_?";
    PetscObject        cached = NULL;
    DMBoundaryType     boundary[3];
    PetscInt           size[3], dof, gs[3], gm[3], nleaves, maxleaf, nselected = 0, *selected = NULL;
    const PetscInt    *ilocal;
    VecScatter         gtol, embedded;
    PetscBool          per_point;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(axis >= 0 && axis < 3, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid direction '%c'", direction);
    *plane = NULL;
    PetscCall(DMDAGetInfo(dm, NULL, &size[0], &size[1], &size[2], NULL, NULL, NULL, &dof, NULL,
                          &boundary[0], &boundary[1], &boundary[2], NULL));
    if (boundary[axis] != DM_BOUNDARY_PERIODIC) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    key[sizeof(key) - 2] = direction;
    PetscCall(PetscObjectQuery((PetscObject)dm, key, &cached));
    if (cached) { *plane = (VecScatter)cached; PROFILE_HELPER_END; PetscFunctionReturn(0); }

    PetscCall(DMDAGetGhostCorners(dm, &gs[0], &gs[1], &gs[2], &gm[0], &gm[1], &gm[2]));
    PetscCall(DMDAGetScatter(dm, &gtol, NULL));
    PetscCall(PetscSFGetGraph(gtol, NULL, &nleaves, &ilocal, NULL));
    PetscCall(PetscSFGetLeafRange(gtol, NULL, &maxleaf));
    /* Leaves index whole points when the scatter moves dof-sized blocks, single entries otherwise. */
    per_point = (PetscBool)(dof == 1 || maxleaf < gm[0] * gm[1] * gm[2]);
    PetscCall(PetscMalloc1(nleaves, &selected));
    for (PetscInt edge = 0; edge < nleaves; edge++) {
        const PetscInt leaf  = ilocal ? ilocal[edge] : edge;
        const PetscInt point = per_point ? leaf : leaf / dof;
        const PetscInt index[3] = {gs[0] + point % gm[0], gs[1] + (point / gm[0]) % gm[1], gs[2] + point / (gm[0] * gm[1])};
        const PetscInt wrapped = ((index[axis] % size[axis]) + size[axis]) % size[axis];
        if (wrapped == 0 || wrapped == size[axis] - 1) selected[nselected++] = edge;
    }
    PetscCall(PetscSFCreateEmbeddedLeafSF(gtol, nselected, selected, &embedded));
    PetscCall(PetscFree(selected));
    PetscCall(PetscObjectCompose((PetscObject)dm, key, (PetscObject)embedded));
    PetscCall(PetscSFDestroy(&embedded));
    PetscCall(PetscObjectQuery((PetscObject)dm, key, &cached));
    *plane = (VecScatter)cached;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Brings a field's local vector up to date after a periodic transfer along one direction.
 * @details Precondition: the local vector was the ghost update of the global vector before the
 *          transfer. The transfer wrote only the end planes along direction, so rescattering the
 *          entries that mirror those planes and reapplying the periodic repair leaves exactly what
 *          UpdateLocalGhosts() would. Static so that only the periodic synchronizers below, which
 *          establish the precondition, can call it; -picurv_verify_ghost_updates checks each call.
 * @param[in,out] user Block context owning the field.
 * @param[in] field_id Field just transferred.
 * @param[in] direction Direction of that transfer.
 * @return 0 on success.
 */
static PetscErrorCode RefreshPeriodicPlaneGhosts(UserCtx *user, FieldId field_id, char direction)
{
    FieldView  view;
    VecScatter plane;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(FieldGetView(user, field_id, &view));
    PetscCall(GetPeriodicPlaneScatter(view.dm, direction, &plane));
    if (!plane) {
        PetscCall(UpdateLocalGhosts(user, field_id));
    } else {
        PetscCall(VecScatterBegin(plane, view.global_vec, view.local_vec, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(plane, view.global_vec, view.local_vec, INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(RepairPeriodicFieldGhosts(user, &view, view.local_vec));
    }
    PetscCall(VerifyLocalGhostsCurrent(user, field_id, "after a plane-only refresh"));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Body of SynchronizePeriodicCellFields(), optionally without its first or last ghost update.
 * @details UpdateLocalGhosts() rebuilds the whole local vector from the global one, so it is
 *          an exact no-op on a field whose global vector has not been written since its last
 *          update. A caller that knows this from its own call sequence passes locals_current
 *          and skips that leading exchange. A caller that updates every listed field again
 *          before anything reads its local vector passes final_update = PETSC_FALSE and skips
 *          the update after the last direction's transfers, leaving the local vectors stale
 *          until then. The transfers themselves always run.
 * @param[in,out] user Block context owning the fields.
 * @param[in] num_fields Number of cell fields to synchronize.
 * @param[in] field_ids Field identifiers.
 * @param[in] locals_current PETSC_TRUE when every listed field's local vector is already
 *            the ghost update of its current global vector.
 * @param[in] final_update PETSC_FALSE to leave the last direction's update to the caller.
 * @return 0 on success.
 */
static PetscErrorCode SynchronizePeriodicCellFieldsFrom(UserCtx *user, PetscInt num_fields,
                                                        const FieldId field_ids[],
                                                        PetscBool locals_current,
                                                        PetscBool final_update)
{
    PetscErrorCode ierr;
    PetscBool      periodic_i;
    PetscBool      periodic_j;
    PetscBool      periodic_k;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    PetscCheck(num_fields >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Number of cell fields cannot be negative.");
    if (num_fields == 0) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    PetscCheck(field_ids != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Cell field-ID array cannot be NULL.");

    {
        PetscBool periodic[3];
        GetPeriodicAxes(user->simCtx, periodic);
        periodic_i = periodic[0]; periodic_j = periodic[1]; periodic_k = periodic[2];
    }

    if (!periodic_i && !periodic_j && !periodic_k) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    for (PetscInt field = 0; field < num_fields; field++) {
        if (locals_current) { ierr = VerifyLocalGhostsCurrent(user, field_ids[field], "claimed current by the caller"); CHKERRQ(ierr); }
        else { ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr); }
    }

    /* Each transfer writes only the end planes of its direction, and the local vectors are
     * current before it, so a plane-only refresh restores exactly a full update. */
    if (periodic_i) {
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFieldByDirection(user, field_ids[field], 'i'); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; (final_update || periodic_j || periodic_k) && field < num_fields; field++) {
            ierr = RefreshPeriodicPlaneGhosts(user, field_ids[field], 'i'); CHKERRQ(ierr);
        }
    }

    if (periodic_j) {
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFieldByDirection(user, field_ids[field], 'j'); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; (final_update || periodic_k) && field < num_fields; field++) {
            ierr = RefreshPeriodicPlaneGhosts(user, field_ids[field], 'j'); CHKERRQ(ierr);
        }
    }

    if (periodic_k) {
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFieldByDirection(user, field_ids[field], 'k'); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; final_update && field < num_fields; field++) {
            ierr = RefreshPeriodicPlaneGhosts(user, field_ids[field], 'k'); CHKERRQ(ierr);
        }
    }
    for (PetscInt field = 0; !final_update && field < num_fields; field++) {
        ierr = PoisonStaleLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
    }

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SynchronizePeriodicCellFields().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 */
PetscErrorCode SynchronizePeriodicCellFields(UserCtx *user, PetscInt num_fields, const FieldId field_ids[])
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    ierr = SynchronizePeriodicCellFieldsFrom(user, num_fields, field_ids, PETSC_FALSE, PETSC_TRUE); CHKERRQ(ierr);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Resolves one registered persistent single-face-family field.
 */
static PetscErrorCode GetPersistentFaceField(UserCtx *user, FieldId field_id,
                                             char face_direction, DM *dm,
                                             Vec *global_vec, Vec *local_vec,
                                             PetscInt *dof)
{
    FieldView   field_view;
    FieldLayout expected_layout;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(face_direction == 'i' || face_direction == 'j' || face_direction == 'k',
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Invalid face direction '%c'; expected 'i', 'j', or 'k'.", face_direction);

    expected_layout = face_direction == 'i' ? FIELD_LAYOUT_I_FACE :
                      (face_direction == 'j' ? FIELD_LAYOUT_J_FACE : FIELD_LAYOUT_K_FACE);
    PetscCall(FieldGetView(user, field_id, &field_view));
    PetscCheck((field_view.descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_FACE_SYNC) != 0u &&
               field_view.descriptor->layout == expected_layout,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Field '%s' is not registered for %c-face periodic synchronization.",
               field_view.descriptor->canonical_name, face_direction);

    *dm = field_view.dm;
    *global_vec = field_view.global_vec;
    *local_vec = field_view.local_vec;
    *dof = field_view.descriptor->dof;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Returns whether a registered face field stores physical coordinates.
 */
static PetscErrorCode IsFaceCenterCoordinateField(FieldId field_id, PetscBool *is_coordinate)
{
    const FieldDescriptor *descriptor = NULL;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(is_coordinate != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Face-coordinate result cannot be NULL.");
    PetscCall(FieldGetDescriptor(field_id, &descriptor));
    *is_coordinate = (PetscBool)((descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT) != 0u);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Applies geometric translations to wrapped face-center ghost coordinates.
 */
static PetscErrorCode TranslatePeriodicFaceCenterGhosts(UserCtx *user, Vec local_vec)
{
    DMDALocalInfo info;
    Cmpnts ***array;
    PetscBool periodic[3];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    info = user->info;
    PetscCall(DMDAVecGetArray(user->fda, local_vec, &array));

    GetPeriodicAxes(user->simCtx, periodic);
    for (PetscInt axis = 0; axis < 3; axis++) {
        const PetscBool active = periodic[axis];
        Cmpnts translation;

        if (!active) continue;
        PetscCheck(user->periodic_translation_valid[axis], PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                   "Periodic face-center synchronization requires validated %c-direction geometry.",
                   "XYZ"[axis]);
        translation = user->periodic_translation[axis];

        for (PetscInt k = info.gzs; k < info.gzs + info.gzm; k++) {
            for (PetscInt j = info.gys; j < info.gys + info.gym; j++) {
                for (PetscInt i = info.gxs; i < info.gxs + info.gxm; i++) {
                    PetscReal scale = 0.0;
                    const PetscInt index = axis == 0 ? i : (axis == 1 ? j : k);
                    const PetscInt size = axis == 0 ? info.mx : (axis == 1 ? info.my : info.mz);
                    if (index < 0) scale = -1.0;
                    else if (index >= size) scale = 1.0;
                    if (scale == 0.0) continue;
                    array[k][j][i].x += scale * translation.x;
                    array[k][j][i].y += scale * translation.y;
                    array[k][j][i].z += scale * translation.z;
                }
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, local_vec, &array));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/** @brief Transfers one registered face-family field along one periodic axis. */
static PetscErrorCode TransferPeriodicFaceFieldByDirection(UserCtx *user, FieldId field_id,
                                                           char face_direction, char periodic_direction)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;
    DM dm;
    Vec global_vec, local_vec;
    PetscInt dof;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(periodic_direction == 'i' || periodic_direction == 'j' || periodic_direction == 'k',
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Invalid periodic direction '%c'; expected 'i', 'j', or 'k'.", periodic_direction);
    PetscCall(GetPersistentFaceField(user, field_id, face_direction, &dm, &global_vec, &local_vec, &dof));

    if (dof == 1) {
        PetscReal ***global_array, ***local_array;
        ierr = DMDAVecGetArray(dm, global_vec, &global_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);

        if (periodic_direction == 'i') {
            if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][0] = local_array[k][j][-2];
            if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][mx-1] = local_array[k][j][mx+1];
        } else if (periodic_direction == 'j') {
            if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][0][i] = local_array[k][-2][i];
            if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][my-1][i] = local_array[k][my+1][i];
        } else {
            if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[0][j][i] = local_array[-2][j][i];
            if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[mz-1][j][i] = local_array[mz+1][j][i];
        }
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    } else {
        Cmpnts ***global_array, ***local_array;
        ierr = DMDAVecGetArray(dm, global_vec, &global_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);

        if (periodic_direction == 'i') {
            if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][0] = local_array[k][j][-2];
            if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][mx-1] = local_array[k][j][mx+1];
        } else if (periodic_direction == 'j') {
            if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][0][i] = local_array[k][-2][i];
            if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][my-1][i] = local_array[k][my+1][i];
        } else {
            if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[0][j][i] = local_array[-2][j][i];
            if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[mz-1][j][i] = local_array[mz+1][j][i];
        }
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    }

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

// Implements SynchronizePeriodicFaceFields(); the public header owns the
// rendered API contract.
PetscErrorCode SynchronizePeriodicFaceFields(UserCtx *user, char face_direction,
                                             PetscInt num_fields, const FieldId field_ids[])
{
    PetscErrorCode ierr;
    const char periodic_directions[3] = {'i', 'j', 'k'};
    const BCFace negative_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace positive_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(num_fields >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Number of face fields cannot be negative.");
    if (num_fields == 0) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    PetscCheck(field_ids != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Face field-ID array cannot be NULL.");

    for (PetscInt field = 0; field < num_fields; field++) {
        DM dm;
        Vec global_vec, local_vec;
        PetscInt dof;
        PetscBool is_coordinate;
        PetscCall(GetPersistentFaceField(user, field_ids[field], face_direction,
                                         &dm, &global_vec, &local_vec, &dof));
        ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
        ierr = IsFaceCenterCoordinateField(field_ids[field], &is_coordinate); CHKERRQ(ierr);
        if (is_coordinate) {
            ierr = TranslatePeriodicFaceCenterGhosts(user, local_vec); CHKERRQ(ierr);
        }
    }

    for (PetscInt direction = 0; direction < 3; direction++) {
        const PetscBool active =
            user->boundary_faces[negative_faces[direction]].mathematical_type == PERIODIC ||
            user->boundary_faces[positive_faces[direction]].mathematical_type == PERIODIC;
        if (!active) continue;

        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFaceFieldByDirection(user, field_ids[field], face_direction,
                                                        periodic_directions[direction]); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; field < num_fields; field++) {
            DM dm;
            Vec global_vec, local_vec;
            PetscInt dof;
            PetscBool is_coordinate;
            ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
            ierr = IsFaceCenterCoordinateField(field_ids[field], &is_coordinate); CHKERRQ(ierr);
            if (is_coordinate) {
                PetscCall(GetPersistentFaceField(user, field_ids[field], face_direction,
                                                 &dm, &global_vec, &local_vec, &dof));
                ierr = TranslatePeriodicFaceCenterGhosts(user, local_vec); CHKERRQ(ierr);
            }
        }
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Resolves one registered persistent component-staggered field.
 */
static PetscErrorCode GetPersistentStaggeredField(UserCtx *user, FieldId field_id,
                                                  DM *dm, Vec *global_vec, Vec *local_vec)
{
    FieldView field_view;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(FieldGetView(user, field_id, &field_view));
    PetscCheck((field_view.descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_STAGGERED_SYNC) != 0u &&
               field_view.descriptor->layout == FIELD_LAYOUT_COMPONENT_STAGGERED,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Field '%s' is not registered for component-staggered periodic synchronization.",
               field_view.descriptor->canonical_name);
    *dm = field_view.dm;
    *global_vec = field_view.global_vec;
    *local_vec = field_view.local_vec;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Transfers one component-staggered field along one periodic axis.
 */
static PetscErrorCode TransferPeriodicStaggeredFieldByDirection(UserCtx *user, FieldId field_id,
                                                                char periodic_direction)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;
    DM dm;
    Vec global_vec, local_vec;
    Cmpnts ***global_array, ***local_array;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(periodic_direction == 'i' || periodic_direction == 'j' || periodic_direction == 'k',
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Invalid periodic direction '%c'; expected 'i', 'j', or 'k'.", periodic_direction);
    PetscCall(GetPersistentStaggeredField(user, field_id, &dm, &global_vec, &local_vec));

    ierr = DMDAVecGetArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);

    if (periodic_direction == 'i') {
        if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][0] = local_array[k][j][-2];
        if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][mx-1] = local_array[k][j][mx+1];
    } else if (periodic_direction == 'j') {
        if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][0][i] = local_array[k][-2][i];
        if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][my-1][i] = local_array[k][my+1][i];
    } else {
        if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
            for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[0][j][i] = local_array[-2][j][i];
        if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz)
            for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[mz-1][j][i] = local_array[mz+1][j][i];
    }

    ierr = DMDAVecRestoreArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Body of SynchronizePeriodicStaggeredFields(), optionally without its leading ghost update.
 * @details As for SynchronizePeriodicCellFieldsFrom(): locals_current skips only the leading
 *          UpdateLocalGhosts(), an exact no-op when no global vector was written since its
 *          last update; the per-direction transfers and their updates always run.
 * @param[in,out] user Block context owning the fields.
 * @param[in] num_fields Number of staggered fields to synchronize.
 * @param[in] field_ids Field identifiers.
 * @param[in] locals_current PETSC_TRUE when every listed field's local vector is already
 *            the ghost update of its current global vector.
 * @return 0 on success.
 */
static PetscErrorCode SynchronizePeriodicStaggeredFieldsFrom(UserCtx *user, PetscInt num_fields,
                                                             const FieldId field_ids[],
                                                             PetscBool locals_current)
{
    PetscErrorCode ierr;
    const char periodic_directions[3] = {'i', 'j', 'k'};
    const BCFace negative_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace positive_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(num_fields >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Number of staggered fields cannot be negative.");
    if (num_fields == 0) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    PetscCheck(field_ids != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Staggered field-ID array cannot be NULL.");

    for (PetscInt field = 0; field < num_fields; field++) {
        DM dm;
        Vec global_vec, local_vec;
        PetscCall(GetPersistentStaggeredField(user, field_ids[field], &dm, &global_vec, &local_vec));
        if (locals_current) { ierr = VerifyLocalGhostsCurrent(user, field_ids[field], "claimed current by the caller"); CHKERRQ(ierr); }
        else { ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr); }
    }

    for (PetscInt direction = 0; direction < 3; direction++) {
        const PetscBool active =
            user->boundary_faces[negative_faces[direction]].mathematical_type == PERIODIC ||
            user->boundary_faces[positive_faces[direction]].mathematical_type == PERIODIC;
        if (!active) continue;

        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicStaggeredFieldByDirection(user, field_ids[field],
                                                             periodic_directions[direction]); CHKERRQ(ierr);
        }
        /* As in SynchronizePeriodicCellFieldsFrom(): the transfer wrote only this direction's
         * end planes into current local vectors, so a plane-only refresh suffices. */
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = RefreshPeriodicPlaneGhosts(user, field_ids[field], periodic_directions[direction]); CHKERRQ(ierr);
        }
    }

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SynchronizePeriodicStaggeredFields().
 */
PetscErrorCode SynchronizePeriodicStaggeredFields(UserCtx *user, PetscInt num_fields,
                                                  const FieldId field_ids[])
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    ierr = SynchronizePeriodicStaggeredFieldsFrom(user, num_fields, field_ids, PETSC_FALSE); CHKERRQ(ierr);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PreparePeriodicQuickStencilFields().
 */
PetscErrorCode PreparePeriodicQuickStencilFields(UserCtx *user, Vec local_vector_field,
                                                 Vec local_scalar_field)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***vector_array;
    PetscReal ***scalar_array;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
    const PetscInt gxs = info.gxs, gxe = info.gxs + info.gxm;
    const PetscInt gys = info.gys, gye = info.gys + info.gym;
    const PetscInt gzs = info.gzs, gze = info.gzs + info.gzm;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(local_vector_field && local_scalar_field, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "QUICK stencil repair requires both local vector and scalar fields.");
    PetscCall(DMDAVecGetArray(user->fda, local_vector_field, &vector_array));
    PetscCall(DMDAVecGetArray(user->da, local_scalar_field, &scalar_array));

    if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt j = gys; j < gye; j++) {
            vector_array[k][j][-1] = vector_array[k][j][-3];
            scalar_array[k][j][-1] = scalar_array[k][j][-3];
        }
    }
    if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == info.mx) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt j = gys; j < gye; j++) {
            vector_array[k][j][info.mx] = vector_array[k][j][info.mx + 2];
            scalar_array[k][j][info.mx] = scalar_array[k][j][info.mx + 2];
        }
    }
    if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[k][-1][i] = vector_array[k][-3][i];
            scalar_array[k][-1][i] = scalar_array[k][-3][i];
        }
    }
    if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == info.my) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[k][info.my][i] = vector_array[k][info.my + 2][i];
            scalar_array[k][info.my][i] = scalar_array[k][info.my + 2][i];
        }
    }
    if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) {
        for (PetscInt j = gys; j < gye; j++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[-1][j][i] = vector_array[-3][j][i];
            scalar_array[-1][j][i] = scalar_array[-3][j][i];
        }
    }
    if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == info.mz) {
        for (PetscInt j = gys; j < gye; j++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[info.mz][j][i] = vector_array[info.mz + 2][j][i];
            scalar_array[info.mz][j][i] = scalar_array[info.mz + 2][j][i];
        }
    }

    PetscCall(DMDAVecRestoreArray(user->da, local_scalar_field, &scalar_array));
    PetscCall(DMDAVecRestoreArray(user->fda, local_vector_field, &vector_array));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SynchronizePeriodicLocalStaggeredField().
 */
PetscErrorCode SynchronizePeriodicLocalStaggeredField(UserCtx *user, Vec local_field)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***array;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(local_field, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Local staggered field cannot be NULL.");
    PetscCall(DMLocalToLocalBegin(user->fda, local_field, INSERT_VALUES, local_field));
    PetscCall(DMLocalToLocalEnd(user->fda, local_field, INSERT_VALUES, local_field));
    PetscCall(DMDAVecGetArray(user->fda, local_field, &array));

    if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) array[k][j][0].x = array[k][j][-2].x;
    if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == info.mx)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) array[k][j][info.mx - 1].x = array[k][j][info.mx + 1].x;
    if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) array[k][0][i].y = array[k][-2][i].y;
    if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == info.my)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) array[k][info.my - 1][i].y = array[k][info.my + 1][i].y;
    if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
        for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) array[0][j][i].z = array[-2][j][i].z;
    if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == info.mz)
        for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) array[info.mz - 1][j][i].z = array[info.mz + 1][j][i].z;

    PetscCall(DMDAVecRestoreArray(user->fda, local_field, &array));
    PetscCall(DMLocalToLocalBegin(user->fda, local_field, INSERT_VALUES, local_field));
    PetscCall(DMLocalToLocalEnd(user->fda, local_field, INSERT_VALUES, local_field));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ApplyMetricsPeriodicBCs().
 * @see ApplyMetricsPeriodicBCs()
 */
PetscErrorCode ApplyMetricsPeriodicBCs(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    const FieldId cell_fields[] = {FIELD_ID_AJ};
    const FieldId i_face_fields[] = {FIELD_ID_CENTX, FIELD_ID_CSI, FIELD_ID_ICSI,
                                     FIELD_ID_IETA, FIELD_ID_IZET, FIELD_ID_IAJ};
    const FieldId j_face_fields[] = {FIELD_ID_CENTY, FIELD_ID_ETA, FIELD_ID_JCSI,
                                     FIELD_ID_JETA, FIELD_ID_JZET, FIELD_ID_JAJ};
    const FieldId k_face_fields[] = {FIELD_ID_CENTZ, FIELD_ID_ZET, FIELD_ID_KCSI,
                                     FIELD_ID_KETA, FIELD_ID_KZET, FIELD_ID_KAJ};

    ierr = SynchronizePeriodicCellFields(user, 1, cell_fields); CHKERRQ(ierr);
    ierr = SynchronizePeriodicFaceFields(user, 'i', 6, i_face_fields); CHKERRQ(ierr);
    ierr = SynchronizePeriodicFaceFields(user, 'j', 6, j_face_fields); CHKERRQ(ierr);
    ierr = SynchronizePeriodicFaceFields(user, 'k', 6, k_face_fields); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ApplyPeriodicBCs().
 * @see ApplyPeriodicBCs()
 */
PetscErrorCode ApplyPeriodicBCs(UserCtx *user, PetscBool final_cell_update)
{
    PetscErrorCode ierr;
    PetscBool periodic[3], is_any_periodic;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    GetPeriodicAxes(user->simCtx, periodic);
    is_any_periodic = (PetscBool)(periodic[0] || periodic[1] || periodic[2]);

    if (!is_any_periodic) {
        LOG_ALLOW(GLOBAL,LOG_TRACE, "No periodic boundaries defined; skipping ApplyPeriodicBCs.\n");
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(GLOBAL, LOG_TRACE, "Applying periodic boundary conditions for all fields.\n");

    // STEP 1: Synchronize periodic cell-centered fields in deterministic direction order.
    const FieldId cell_fields[] = {FIELD_ID_UCAT, FIELD_ID_P, FIELD_ID_NVERT};
    ierr = SynchronizePeriodicCellFieldsFrom(user, 3, cell_fields, PETSC_FALSE, final_cell_update); CHKERRQ(ierr);

    // STEP 2: Synchronize persistent staggered endpoints and repair local
    // component-normal ghosts through UpdateLocalGhosts().
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref UpdateDummyCells().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 * @see UpdateDummyCells()
 */
PetscErrorCode UpdateDummyCells(UserCtx *user, FieldId field_id)
{
    DMDALocalInfo info = user->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;
    FieldView     view;
    PetscScalar ****field = NULL, ****ubcs = NULL;
    PetscInt      dof;
    PetscBool     from_boundary_value, fills[6];

    // --- Calculate shrunken loop ranges to avoid edges and corners ---
    PetscInt lxs = (xs == 0) ? xs + 1 : xs, lxe = (xe == mx) ? xe - 1 : xe;
    PetscInt lys = (ys == 0) ? ys + 1 : ys, lye = (ye == my) ? ye - 1 : ye;
    PetscInt lzs = (zs == 0) ? zs + 1 : zs, lze = (ze == mz) ? ze - 1 : ze;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(FieldGetView(user, field_id, &view));
    PetscCheck(view.descriptor->layout == FIELD_LAYOUT_CELL_CENTERED, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "UpdateDummyCells fills cell-centred fields; '%s' is not one.", view.descriptor->canonical_name);
    /* The dummy value is the quantity's boundary condition. Velocity takes the
     * boundary value the handlers set, placed so the face average equals it. Pressure
     * and the particle scalar take the adjacent cell's value: zero normal gradient. The
     * eddy viscosity mirrors only across a free-slip wall, where the flow mirrors and the
     * face average must equal the interior value; a wall that exerts shear replaces it on
     * the face, and no other face's eddy-viscosity dummy cells are written here. */
    PetscCheck(field_id == FIELD_ID_UCAT || field_id == FIELD_ID_P || field_id == FIELD_ID_PSI ||
               field_id == FIELD_ID_NU_T,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "UpdateDummyCells has no boundary rule for '%s'; one must be chosen for the quantity "
               "before its dummy cells are filled.", view.descriptor->canonical_name);
    from_boundary_value = (PetscBool)(field_id == FIELD_ID_UCAT);
    for (PetscInt face = 0; face < 6; ++face) {
        const BoundaryFaceConfig *config = &user->boundary_faces[face];

        fills[face] = (field_id == FIELD_ID_NU_T)
            ? BoundaryFaceIsFreeSlipWall(config)
            : (PetscBool)(config->mathematical_type != PERIODIC);
    }
    dof = view.descriptor->dof;

    PetscCall(DMDAVecGetArrayDOF(view.dm, view.global_vec, &field));
    if (from_boundary_value) PetscCall(DMDAVecGetArrayDOF(user->fda, user->Bcs.Ubcs, &ubcs));

    // -X Face
    if (fills[BC_FACE_NEG_X] && xs == 0) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++) for (PetscInt c = 0; c < dof; c++)
            field[k][j][xs][c] = from_boundary_value ? 2.0 * ubcs[k][j][xs][c] - field[k][j][xs + 1][c]
                                                     : field[k][j][xs + 1][c];
    }
    // +X Face
    if (fills[BC_FACE_POS_X] && xe == mx) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++) for (PetscInt c = 0; c < dof; c++)
            field[k][j][xe-1][c] = from_boundary_value ? 2.0 * ubcs[k][j][xe-1][c] - field[k][j][xe - 2][c]
                                                       : field[k][j][xe - 2][c];
    }
    // -Y Face
    if (fills[BC_FACE_NEG_Y] && ys == 0) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++) for (PetscInt c = 0; c < dof; c++)
            field[k][ys][i][c] = from_boundary_value ? 2.0 * ubcs[k][ys][i][c] - field[k][ys + 1][i][c]
                                                     : field[k][ys + 1][i][c];
    }
    // +Y Face
    if (fills[BC_FACE_POS_Y] && ye == my) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++) for (PetscInt c = 0; c < dof; c++)
            field[k][ye-1][i][c] = from_boundary_value ? 2.0 * ubcs[k][ye-1][i][c] - field[k][ye-2][i][c]
                                                       : field[k][ye-2][i][c];
    }
    // -Z Face
    if (fills[BC_FACE_NEG_Z] && zs == 0) {
        for (PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++) for (PetscInt c = 0; c < dof; c++)
            field[zs][j][i][c] = from_boundary_value ? 2.0 * ubcs[zs][j][i][c] - field[zs + 1][j][i][c]
                                                     : field[zs + 1][j][i][c];
    }
    // +Z Face
    if (fills[BC_FACE_POS_Z] && ze == mz) {
        for (PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++) for (PetscInt c = 0; c < dof; c++)
            field[ze-1][j][i][c] = from_boundary_value ? 2.0 * ubcs[ze-1][j][i][c] - field[ze-2][j][i][c]
                                                       : field[ze-2][j][i][c];
    }

    if (from_boundary_value) PetscCall(DMDAVecRestoreArrayDOF(user->fda, user->Bcs.Ubcs, &ubcs));
    PetscCall(DMDAVecRestoreArrayDOF(view.dm, view.global_vec, &field));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: the face separating an edge dummy cell from one neighbour.
 * @details Local to this translation unit. The two cells differ in exactly one index,
 *          and the target's own value of it says which side of the domain the face is on.
 */
static BCFace EdgeSeparatingFace(PetscInt tk, PetscInt tj, PetscInt ti,
                                 PetscInt nk, PetscInt nj, PetscInt ni)
{
    if (nk != tk) return (tk == 0) ? BC_FACE_NEG_Z : BC_FACE_POS_Z;
    if (nj != tj) return (tj == 0) ? BC_FACE_NEG_Y : BC_FACE_POS_Y;
    (void)ni;
    return (ti == 0) ? BC_FACE_NEG_X : BC_FACE_POS_X;
}

/**
 * @brief Internal helper: fills one edge dummy cell from its two face-adjacent dummy cells.
 * @details Local to this translation unit. Across a free-slip wall the flow is the mirror
 *          image of the flow inside, so the edge cell is the mirror image of the neighbour
 *          on the other side of that wall: a scalar is copied, and a velocity has its
 *          component along the wall's metric reversed. Elsewhere the edge cell is the mean
 *          of the two neighbours.
 */
static void FillEdgeCell(PetscScalar ****f, PetscInt dof, const PetscBool free_slip[6],
                         Cmpnts ***const metric[3],
                         PetscInt tk, PetscInt tj, PetscInt ti,
                         PetscInt ak, PetscInt aj, PetscInt ai,
                         PetscInt bk, PetscInt bj, PetscInt bi)
{
    const BCFace face_a = EdgeSeparatingFace(tk, tj, ti, ak, aj, ai);
    const BCFace face_b = EdgeSeparatingFace(tk, tj, ti, bk, bj, bi);
    PetscInt sk = ak, sj = aj, si = ai;
    BCFace mirror = face_a;

    if (!free_slip[face_a] && !free_slip[face_b]) {
        for (PetscInt c = 0; c < dof; c++) f[tk][tj][ti][c] = 0.5 * (f[ak][aj][ai][c] + f[bk][bj][bi][c]);
        return;
    }
    if (!free_slip[face_a]) { sk = bk; sj = bj; si = bi; mirror = face_b; }
    for (PetscInt c = 0; c < dof; c++) f[tk][tj][ti][c] = f[sk][sj][si][c];
    if (dof == 3) {
        const Cmpnts n = metric[mirror / 2][sk][sj][si];
        const PetscReal area_sq = n.x * n.x + n.y * n.y + n.z * n.z;
        const PetscReal normal = (area_sq > 0.0)
            ? (f[sk][sj][si][0] * n.x + f[sk][sj][si][1] * n.y + f[sk][sj][si][2] * n.z) / area_sq : 0.0;

        f[tk][tj][ti][0] -= 2.0 * normal * n.x;
        f[tk][tj][ti][1] -= 2.0 * normal * n.y;
        f[tk][tj][ti][2] -= 2.0 * normal * n.z;
    }
}

/**
 * @brief Implementation of \ref UpdateCornerNodes().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 * @see UpdateCornerNodes()
 */
PetscErrorCode UpdateCornerNodes(UserCtx *user, FieldId field_id)
{
    DMDALocalInfo info = user->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;
    FieldView     view;
    PetscScalar ****f = NULL;
    PetscInt      dof;
    PetscBool     free_slip[6], any_free_slip = PETSC_FALSE;
    Cmpnts     ***metric[3] = {NULL, NULL, NULL};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(FieldGetView(user, field_id, &view));
    PetscCheck(view.descriptor->layout == FIELD_LAYOUT_CELL_CENTERED, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "UpdateCornerNodes fills cell-centred fields; '%s' is not one.", view.descriptor->canonical_name);
    dof = view.descriptor->dof;
    for (PetscInt face = 0; face < 6; ++face) {
        const BoundaryFaceConfig *config = &user->boundary_faces[face];

        free_slip[face] = BoundaryFaceIsFreeSlipWall(config);
        any_free_slip = (PetscBool)(any_free_slip || free_slip[face]);
    }
    PetscCall(DMDAVecGetArrayDOF(view.dm, view.global_vec, &f));
    /* A velocity mirrors along the wall's own normal, which only a free-slip edge needs. */
    if (any_free_slip && dof == 3) {
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts ***)&metric[0]));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts ***)&metric[1]));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts ***)&metric[2]));
    }

#define EDGE_FILL(tk, tj, ti, ak, aj, ai, bk, bj, bi) \
    FillEdgeCell(f, dof, free_slip, metric, tk, tj, ti, ak, aj, ai, bk, bj, bi)

    // --- Update Edges and Corners ---
    // The order of these blocks ensures that corners (where 3 faces meet) are
    // computed using data from edges (where 2 faces meet), which are computed first.
    // Edges connected to the -Z face (k=zs)
    if (zs == 0) {
        if (xs == 0)  for (PetscInt j = ys; j < ye; j++) { EDGE_FILL(zs, j, xs,     zs+1, j, xs,     zs, j, xs+1); }
        if (xe == mx) for (PetscInt j = ys; j < ye; j++) { EDGE_FILL(zs, j, mx-1,   zs+1, j, mx-1,   zs, j, mx-2); }
        if (ys == 0)  for (PetscInt i = xs; i < xe; i++) { EDGE_FILL(zs, ys, i,     zs+1, ys, i,     zs, ys+1, i); }
        if (ye == my) for (PetscInt i = xs; i < xe; i++) { EDGE_FILL(zs, my-1, i,   zs+1, my-1, i,   zs, my-2, i); }
    }
    // Edges connected to the +Z face (k=ze-1)
    if (ze == mz) {
        if (xs == 0)  for (PetscInt j = ys; j < ye; j++) { EDGE_FILL(mz-1, j, xs,   mz-2, j, xs,     mz-1, j, xs+1); }
        if (xe == mx) for (PetscInt j = ys; j < ye; j++) { EDGE_FILL(mz-1, j, mx-1, mz-2, j, mx-1,   mz-1, j, mx-2); }
        if (ys == 0)  for (PetscInt i = xs; i < xe; i++) { EDGE_FILL(mz-1, ys, i,   mz-2, ys, i,     mz-1, ys+1, i); }
        if (ye == my) for (PetscInt i = xs; i < xe; i++) { EDGE_FILL(mz-1, my-1, i, mz-2, my-1, i,   mz-1, my-2, i); }
    }
    // Remaining edges on the XY plane (that are not on Z faces)
    if (ys == 0) {
        if (xs == 0)  for (PetscInt k = zs; k < ze; k++) { EDGE_FILL(k, ys, xs,     k, ys+1, xs,     k, ys, xs+1); }
        if (xe == mx) for (PetscInt k = zs; k < ze; k++) { EDGE_FILL(k, ys, mx-1,   k, ys+1, mx-1,   k, ys, mx-2); }
    }
    if (ye == my) {
        if (xs == 0)  for (PetscInt k = zs; k < ze; k++) { EDGE_FILL(k, my-1, xs,   k, my-2, xs,     k, my-1, xs+1); }
        if (xe == mx) for (PetscInt k = zs; k < ze; k++) { EDGE_FILL(k, my-1, mx-1, k, my-2, mx-1,   k, my-1, mx-2); }
    }
#undef EDGE_FILL

    if (any_free_slip && dof == 3) {
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts ***)&metric[2]));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts ***)&metric[1]));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts ***)&metric[0]));
    }
    PetscCall(DMDAVecRestoreArrayDOF(view.dm, view.global_vec, &f));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief The signed friction velocity the configured wall law gives for one wall-face element.
 *
 * Every law is evaluated from the resolved velocity of the reference cell, relative to the
 * wall, at its wall distance. The laws also return a velocity for the first cell; it is
 * discarded, because the model acts through the wall stress alone (see
 * @ref ApplyWallFunction).
 *
 * Cabot's departure from an equilibrium profile is driven by the pressure gradient in the
 * wall layer, so it is supplied from the resolved pressure at the first cell rather than
 * assumed zero; the pressure is the previous projection's, the same lag every other
 * explicit use of it carries. The other two laws ignore it.
 */
static PetscErrorCode WallModelFrictionVelocity(UserCtx *user, const PetscReal ***pressure,
                                                PetscInt i, PetscInt j, PetscInt k,
                                                PetscReal roughness_height,
                                                PetscReal distance_reference, PetscReal distance_boundary,
                                                Cmpnts velocity_wall, Cmpnts velocity_reference,
                                                const PetscReal normal[3], PetscReal *friction_velocity)
{
    Cmpnts discarded = {0.0, 0.0, 0.0};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    *friction_velocity = 0.0;
    switch ((WallFunctionModel)user->simCtx->wallfunction) {
    case WALL_FUNCTION_WERNER:
        wall_function(user, distance_reference, distance_boundary, velocity_wall, velocity_reference,
                      &discarded, friction_velocity, normal[0], normal[1], normal[2]);
        break;
    case WALL_FUNCTION_CABOT: {
        Cmpnts pressure_gradient = {0.0, 0.0, 0.0};

        if (pressure) {
            PetscCall(ComputeScalarFieldDerivatives(user, i, j, k, (PetscReal ***)pressure, &pressure_gradient));
        }
        wall_function_Cabot(user, roughness_height, distance_reference, distance_boundary,
                            velocity_wall, velocity_reference, &discarded, friction_velocity,
                            normal[0], normal[1], normal[2],
                            pressure_gradient.x, pressure_gradient.y, pressure_gradient.z, 0);
        break;
    }
    case WALL_FUNCTION_LOG_LAW:
    default:
        wall_function_loglaw(user, roughness_height, distance_reference, distance_boundary,
                             velocity_wall, velocity_reference, &discarded, friction_velocity,
                             normal[0], normal[1], normal[2]);
        break;
    }
    /* The magnitude is the friction velocity. Its sign says whether the stress opposes the
       reference flow (positive) or, for Cabot under an adverse pressure gradient strong
       enough to reverse the wall-layer balance, acts along it (negative). A non-finite
       value stops the run here rather than entering the stress. */
    PetscCheck(PetscIsInfOrNanReal(*friction_velocity) == PETSC_FALSE, PETSC_COMM_SELF, PETSC_ERR_FP,
               "Wall model (%s) returned a non-finite friction velocity at cell (%" PetscInt_FMT ", %" PetscInt_FMT
               ", %" PetscInt_FMT ").",
               WallFunctionModelToString((WallFunctionModel)user->simCtx->wallfunction), i, j, k);
    PROFILE_HELPER_END;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Implementation of \ref ApplyWallFunction().
 * @see ApplyWallFunction()
 */
PetscErrorCode ApplyWallFunction(UserCtx *user)
{
    SimCtx              *simCtx = user->simCtx;
    const DMDALocalInfo *info = &user->info;
    const PetscInt       start[3] = {info->xs, info->ys, info->zs};
    const PetscInt       end[3]   = {info->xs + info->xm, info->ys + info->ym, info->zs + info->zm};
    const PetscInt       size[3]  = {info->mx, info->my, info->mz};
    const PetscReal      molecular = 1.0 / simCtx->ren;
    const PetscReal      roughness_height = simCtx->wall_roughness_height;
    const Cmpnts      ***ucat, ***ubcs, ***metric[3];
    const PetscReal   ***nvert, ***aj, ***pressure;
    Cmpnts            ***ucont, ***stress;
    PetscReal         ***friction_velocity;
    PetscInt             wall_faces_at[3][2];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!simCtx->wallfunction) { PROFILE_HELPER_END; PetscFunctionReturn(0); }

    /* One pass, one sample set: the state describes this pass and nothing earlier. Faces the
       model does not reach must not keep a stale stress, and the friction velocity is
       accumulated per cell below, so both start from zero. */
    PetscCall(PetscMemzero(&user->wall_diagnostics, sizeof(user->wall_diagnostics)));
    PetscCall(VecSet(user->Wall_Stress, 0.0));
    PetscCall(VecSet(user->Friction_Velocity, 0.0));

    /* The resolved velocity this pass's Contra2Cart wrote, read from the global view: the
       model reads it and never writes it. */
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, (void *)&ucat));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, (void *)&ubcs));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (void *)&metric[0]));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (void *)&metric[1]));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (void *)&metric[2]));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));
    /* Cabot reads the previous projection's pressure gradient; the other laws ignore it. */
    PetscCall(DMDAVecGetArrayRead(user->da, user->lP, (void *)&pressure));
    PetscCall(DMDAVecGetArray(user->fda, user->Wall_Stress, &stress));
    PetscCall(DMDAVecGetArray(user->da, user->Friction_Velocity, &friction_velocity));

    for (PetscInt axis = 0; axis < 3; axis++) {
        wall_faces_at[axis][0] = BoundaryFaceHasWallShear(&user->boundary_faces[2 * axis]) ? 1 : 0;
        wall_faces_at[axis][1] = BoundaryFaceHasWallShear(&user->boundary_faces[2 * axis + 1]) ? 1 : 0;
    }

    /* Every wall face element is one sample of the law. A cell beside two walls, along an
       edge of the domain, has two such elements, each with its own reference cell, normal
       and stress, so the two walls never write into one another's input. */
    for (PetscInt face = 0; face < 6; face++) {
        const PetscInt  axis = face / 2;
        const PetscBool negative = (PetscBool)(face % 2 == 0);
        const PetscReal side = negative ? 1.0 : -1.0;
        BoundaryFaceElements elements;

        if (!wall_faces_at[axis][negative ? 0 : 1]) continue;
        PetscCall(GetBoundaryFaceElements(user, (BCFace)face, &elements));
        if (elements.n[0] * elements.n[1] == 0) continue;

        const PetscInt face_at   = elements.layer;
        const PetscInt cell      = negative ? face_at + 1 : face_at;
        const PetscInt reference = negative ? face_at + 2 : face_at - 1;
        const PetscInt dummy     = negative ? face_at : face_at + 1;
        /* The reference cell is read from the owned block. A rank holding fewer than three
           nodes normal to a wall would read outside it, so that layout is refused. */
        PetscCheck(reference >= start[axis] && reference < end[axis], PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                   "Wall model on face %s: this rank owns indices %" PetscInt_FMT "-%" PetscInt_FMT
                   " normal to the wall, which does not reach the reference cell %" PetscInt_FMT
                   ". Use a layout with at least three nodes per rank normal to every wall.",
                   BCFaceToString((BCFace)face), start[axis], end[axis] - 1, reference);

        for (PetscInt q = 0; q < elements.n[1]; q++)
        for (PetscInt p = 0; p < elements.n[0]; p++) {
            PetscInt  c[3];

            c[axis] = cell; c[elements.a] = elements.lo[0] + p; c[elements.b] = elements.lo[1] + q;
            const PetscInt i = c[0], j = c[1], k = c[2];
            PetscInt  c_ref[3] = {i, j, k}, c_dummy[3] = {i, j, k}, c_face[3] = {i, j, k};
            PetscReal signed_u_tau = 0.0, u_tau;

            if (nvert[k][j][i] > 0.1) continue;
            c_ref[axis] = reference;
            c_dummy[axis] = dummy;
            c_face[axis] = face_at;

            const Cmpnts    area_vector = metric[axis][c_face[2]][c_face[1]][c_face[0]];
            const PetscReal area = PetscSqrtReal(area_vector.x * area_vector.x + area_vector.y * area_vector.y +
                                                 area_vector.z * area_vector.z);
            /* Unit normal into the fluid, and the first and reference cell-centre distances
               from the wall, from each cell's volume over the wall area. */
            const PetscReal normal[3] = {side * area_vector.x / area, side * area_vector.y / area,
                                         side * area_vector.z / area};
            const PetscReal distance_boundary  = 0.5 / aj[k][j][i] / area;
            const PetscReal distance_reference = 2.0 * distance_boundary +
                                                 0.5 / aj[c_ref[2]][c_ref[1]][c_ref[0]] / area;
            const Cmpnts    velocity_wall = ubcs[c_dummy[2]][c_dummy[1]][c_dummy[0]];
            const Cmpnts    velocity_reference = ucat[c_ref[2]][c_ref[1]][c_ref[0]];
            const Cmpnts    relative_reference = {velocity_reference.x - velocity_wall.x,
                                                  velocity_reference.y - velocity_wall.y,
                                                  velocity_reference.z - velocity_wall.z};
            const Cmpnts    tangent = RemoveFaceNormalComponent(relative_reference, area_vector);
            const PetscReal tangent_speed = PetscSqrtReal(tangent.x * tangent.x + tangent.y * tangent.y +
                                                          tangent.z * tangent.z);

            PetscCall(WallModelFrictionVelocity(user, pressure, i, j, k, roughness_height,
                                                distance_reference, distance_boundary,
                                                velocity_wall, velocity_reference, normal, &signed_u_tau));
            u_tau = PetscAbsReal(signed_u_tau);

            /* The stress opposes the reference cell's motion relative to the wall, along its
               wall-parallel direction. Stored in the face's dummy cell, which no other wall
               face uses, already oriented as the viscous flux through that face, so the
               flux is the stored value times the face area (ComputeViscous). */
            if (tangent_speed > 1.0e-12) {
                const PetscReal scale = side * signed_u_tau * u_tau / tangent_speed;
                stress[c_dummy[2]][c_dummy[1]][c_dummy[0]].x = scale * tangent.x;
                stress[c_dummy[2]][c_dummy[1]][c_dummy[0]].y = scale * tangent.y;
                stress[c_dummy[2]][c_dummy[1]][c_dummy[0]].z = scale * tangent.z;
            }
            /* The wall is impermeable whatever the law. */
            if (axis == 0) ucont[c_face[2]][c_face[1]][c_face[0]].x = 0.0;
            else if (axis == 1) ucont[c_face[2]][c_face[1]][c_face[0]].y = 0.0;
            else ucont[c_face[2]][c_face[1]][c_face[0]].z = 0.0;

            /* The cell's friction velocity is the RMS over its wall faces; the sum of
               squares is accumulated here and reduced once every face has been visited. */
            friction_velocity[k][j][i] += u_tau * u_tau;

            {
                WallModelDiagnosticsState *diagnostics = &user->wall_diagnostics;
                const PetscReal            y_plus = u_tau * distance_boundary / molecular;
                const Cmpnts               first = ucat[k][j][i];
                const Cmpnts               relative_first = {first.x - velocity_wall.x, first.y - velocity_wall.y,
                                                             first.z - velocity_wall.z};
                const Cmpnts               first_tangent = RemoveFaceNormalComponent(relative_first, area_vector);
                const PetscReal            first_speed = PetscSqrtReal(first_tangent.x * first_tangent.x +
                                                                       first_tangent.y * first_tangent.y +
                                                                       first_tangent.z * first_tangent.z);
                /* The eddy viscosity a gradient closure would need on this face to carry
                   the stress across the resolved first cell: nu_eq = u_tau^2 y / |U_1|,
                   less the molecular part. A diagnostic only; the stress itself is applied. */
                if (first_speed > 1.0e-10) {
                    diagnostics->wall_viscosity_sum += PetscMax(u_tau * u_tau * distance_boundary / first_speed - molecular, 0.0);
                }
                if (diagnostics->cells == 0) {
                    diagnostics->friction_velocity_min = u_tau;
                } else {
                    diagnostics->friction_velocity_min = PetscMin(diagnostics->friction_velocity_min, u_tau);
                }
                if (diagnostics->cells == 0 || u_tau > diagnostics->friction_velocity_max) {
                    diagnostics->friction_velocity_max = u_tau;
                    diagnostics->friction_velocity_max_at[0] = face;
                    diagnostics->friction_velocity_max_at[1] = i;
                    diagnostics->friction_velocity_max_at[2] = j;
                    diagnostics->friction_velocity_max_at[3] = k;
                }
                if (diagnostics->cells == 0 || y_plus > diagnostics->y_plus_max) {
                    diagnostics->y_plus_max = y_plus;
                    diagnostics->y_plus_max_at[0] = face;
                    diagnostics->y_plus_max_at[1] = i;
                    diagnostics->y_plus_max_at[2] = j;
                    diagnostics->y_plus_max_at[3] = k;
                }
                diagnostics->friction_velocity_sum += u_tau;
                diagnostics->friction_velocity_sq  += u_tau * u_tau;
                diagnostics->wall_distance_sum     += distance_boundary;
                diagnostics->y_plus_sum            += y_plus;
                diagnostics->cells                 += 1;
            }
        }
    }

    /* Friction velocity per first cell: the RMS over the wall faces it touches. Only the
       wall-face elements are visited, and an edge cell is finalized once, by the first of its
       wall faces in face order, so the cost stays with the wall layers. */
    for (PetscInt face = 0; face < 6; face++) {
        const PetscInt axis = face / 2;
        BoundaryFaceElements elements;

        if (!wall_faces_at[axis][face % 2]) continue;
        PetscCall(GetBoundaryFaceElements(user, (BCFace)face, &elements));
        for (PetscInt q = 0; q < elements.n[1]; q++)
        for (PetscInt p = 0; p < elements.n[0]; p++) {
            PetscInt c[3];
            PetscInt faces = 0, first_face = -1;

            c[axis] = (face % 2 == 0) ? elements.layer + 1 : elements.layer;
            c[elements.a] = elements.lo[0] + p;
            c[elements.b] = elements.lo[1] + q;
            const PetscInt i = c[0], j = c[1], k = c[2];

            if (friction_velocity[k][j][i] <= 0.0) continue;
            for (PetscInt other = 0; other < 6; other++) {
                const PetscInt other_axis = other / 2;
                const PetscInt beside = (other % 2 == 0) ? 1 : size[other_axis] - 2;
                if (wall_faces_at[other_axis][other % 2] && c[other_axis] == beside) {
                    faces++;
                    if (first_face < 0) first_face = other;
                }
            }
            if (first_face != face) continue;
            friction_velocity[k][j][i] = PetscSqrtReal(friction_velocity[k][j][i] / (PetscReal)faces);
        }
    }

    PetscCall(DMDAVecRestoreArray(user->da, user->Friction_Velocity, &friction_velocity));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Wall_Stress, &stress));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lP, (void *)&pressure));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (void *)&metric[2]));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (void *)&metric[1]));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (void *)&metric[0]));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, (void *)&ubcs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, (void *)&ucat));

    /* Refresh the ghosted images here so that every caller sees consistent fields rather
       than each remembering to do it. */
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_U_TAU));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_WALL_STRESS));
    PROFILE_HELPER_END;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Implementation of \ref LogWallModelDiagnostics().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 */
PetscErrorCode LogWallModelDiagnostics(UserCtx *user)
{
    SimCtx   *simCtx = user->simCtx;
    const WallModelDiagnosticsState *state = &user->wall_diagnostics;
    MPI_Comm  comm;

    /* Index 0..3: cell count, u_tau sum, u_tau^2 sum, y+ sum. Index 4..5: wall-distance
       and wall-viscosity sums. One collective rather than six. */
    PetscReal local_sum[6], global_sum[6];
    PetscReal global_max[2];
    PetscReal local_min, global_min;
    /* Index 0: u_tau, 1: y+. MAXLOC carries the owning rank with each maximum, and that
       rank alone contributes the face and cell, so the sum below is its location. The
       lowest rank wins a tie, so the record is the same at any layout that ties. */
    struct { double value; int rank; } local_max[2], owner_max[2];
    PetscInt  local_at[8], global_at[8];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!simCtx->wallfunction) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }

    PetscCall(PetscObjectGetComm((PetscObject)user->da, &comm));

    local_sum[0] = (PetscReal)state->cells;
    local_sum[1] = state->friction_velocity_sum;
    local_sum[2] = state->friction_velocity_sq;
    local_sum[3] = state->y_plus_sum;
    local_sum[4] = state->wall_distance_sum;
    local_sum[5] = state->wall_viscosity_sum;
    /* Both are non-negative, so a rank with no wall cell offers -1 and can never own them. */
    local_max[0].value = (state->cells > 0) ? (double)state->friction_velocity_max : -1.0;
    local_max[1].value = (state->cells > 0) ? (double)state->y_plus_max : -1.0;
    local_max[0].rank  = local_max[1].rank = (int)simCtx->rank;
    /* A rank that owns no wall face must not win the minimum with a zero it never
       measured, so it contributes the identity instead. */
    local_min    = (state->cells > 0) ? state->friction_velocity_min : PETSC_MAX_REAL;

    PetscCallMPI(MPI_Allreduce(local_sum, global_sum, 6, MPIU_REAL, MPI_SUM, comm));
    PetscCallMPI(MPI_Allreduce(local_max, owner_max, 2, MPI_DOUBLE_INT, MPI_MAXLOC, comm));
    PetscCallMPI(MPI_Allreduce(&local_min, &global_min, 1, MPIU_REAL, MPI_MIN, comm));
    global_max[0] = (PetscReal)owner_max[0].value;
    global_max[1] = (PetscReal)owner_max[1].value;
    for (PetscInt n = 0; n < 4; n++) {
        local_at[n]     = (owner_max[0].rank == (int)simCtx->rank) ? state->friction_velocity_max_at[n] : 0;
        local_at[4 + n] = (owner_max[1].rank == (int)simCtx->rank) ? state->y_plus_max_at[n] : 0;
    }
    PetscCallMPI(MPI_Allreduce(local_at, global_at, 8, MPIU_INT, MPI_SUM, comm));

    if (simCtx->rank == 0) {
        const PetscReal cells = global_sum[0];
        FILE           *file  = NULL;
        PetscReal       mean, mean_square, variance;

        /* No wall face anywhere is a configuration fact, not a data point: a run with
           wall functions enabled and no WALL boundary would otherwise emit a row of
           zeros that reads like a converged answer. */
        if (cells <= 0.0) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "Wall model is enabled but no WALL face was corrected this step; "
                      "no diagnostics written.\n");
            PROFILE_FUNCTION_END;
            PetscFunctionReturn(0);
        }

        mean        = global_sum[1] / cells;
        mean_square = global_sum[2] / cells;
        variance    = PetscMax(mean_square - mean * mean, 0.0);

        PetscCall(PicurvOpenDiagnosticsCsv(simCtx, simCtx->analysis_dir, "wall_model.csv",
                                           "step,time,wall_cells,u_tau_mean,u_tau_rms,"
                                           "u_tau_min,u_tau_max,y_plus_mean,y_plus_max,"
                                           "wall_distance_mean,nu_wall_over_nu_mean,physical_time,"
                                           "u_tau_max_face,u_tau_max_i,u_tau_max_j,u_tau_max_k,"
                                           "y_plus_max_face,y_plus_max_i,y_plus_max_j,y_plus_max_k", &file));
        PetscReal physical_time = 0.0;
        PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
        fprintf(file,
                "%d,%.6e,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%d,%d,%d,%s,%d,%d,%d\n",
                (int)simCtx->step, (double)simCtx->ti, (int)cells,
                (double)mean, (double)PetscSqrtReal(variance),
                (double)global_min, (double)global_max[0],
                (double)(global_sum[3] / cells), (double)global_max[1],
                (double)(global_sum[4] / cells),
                (double)(global_sum[5] / cells * simCtx->ren), (double)physical_time,
                BCFaceToString((BCFace)global_at[0]), (int)global_at[1], (int)global_at[2], (int)global_at[3],
                BCFaceToString((BCFace)global_at[4]), (int)global_at[5], (int)global_at[6], (int)global_at[7]);
        PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close the wall-model diagnostics file.");

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "  Wall model (%s): u_tau=%.4e (max %.4e on %s at %d,%d,%d), y+ (mean)=%.2f, "
                  "y+ (max)=%.2f over %d cell(s)\n",
                  WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                  (double)mean, (double)global_max[0], BCFaceToString((BCFace)global_at[0]),
                  (int)global_at[1], (int)global_at[2], (int)global_at[3],
                  (double)(global_sum[3] / cells), (double)global_max[1], (int)cells);
    }

    /* Whether the first cell sits where the selected law is valid is a property of the
       mesh, not of the configuration, so it cannot be settled before the grid exists.
       It is checked here instead, and a run that stays outside the range is stopped
       rather than left to spend its walltime producing a wall stress the law cannot
       support. */
    {
        /* Consecutive samples tolerated outside the range. A startup transient or a
           brief excursion is not worth ending a run over; a mesh that is simply wrong
           never comes back. */
        const PetscInt  PICURV_WALL_YPLUS_GRACE_SAMPLES = 10;
        const PetscReal y_plus_mean = (global_sum[0] > 0.0) ? global_sum[3] / global_sum[0] : 0.0;
        PetscReal       lower = 0.0, upper = 300.0;
        const char     *range_reason = NULL;

        switch ((WallFunctionModel)simCtx->wallfunction) {
        case WALL_FUNCTION_LOG_LAW:
            /* Below 30 the first cell is under the logarithmic region; above 300 the
               law's own implementation reports no valid branch and the correction
               silently falls back to leaving the velocity alone. */
            lower = 30.0;
            upper = 300.0;
            range_reason = "the logarithmic region the law describes";
            break;
        case WALL_FUNCTION_WERNER:
            /* The two-layer form has a valid branch below y+ = 11.81, so only the upper
               end is a limit. */
            lower = 0.0;
            upper = 300.0;
            range_reason = "the region the power law describes";
            break;
        default:
            /* Cabot integrates across the wall layer, so it has no lower bound worth
               asserting; a first cell far outside the layer leaves its ODE nothing to
               integrate over. */
            lower = 0.0;
            upper = 1000.0;
            range_reason = "the wall layer the model integrates over";
            break;
        }

        if (global_sum[0] > 0.0 && (y_plus_mean < lower || y_plus_mean > upper)) {
            user->wall_yplus_excursions += 1;
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "Wall model (%s): first-cell y+ is %.2f, outside %s (%.0f to %.0f). "
                      "The stress it reports is not one the law supports here. Sample %d "
                      "of %d before this run is stopped.\n",
                      WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                      (double)y_plus_mean, range_reason, (double)lower, (double)upper,
                      (int)user->wall_yplus_excursions,
                      (int)PICURV_WALL_YPLUS_GRACE_SAMPLES);

            PetscCheck(user->wall_yplus_excursions < PICURV_WALL_YPLUS_GRACE_SAMPLES,
                       PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                       "Wall model (%s): first-cell y+ has been outside %s (%.0f to %.0f) "
                       "for %d consecutive samples, most recently %.2f. Refine or coarsen "
                       "the wall-normal spacing so the first cell lands in range, or "
                       "resolve the wall and disable the wall function. Stopping now "
                       "rather than spending the run producing a stress the law cannot "
                       "support.",
                       WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                       range_reason, (double)lower, (double)upper,
                       (int)user->wall_yplus_excursions, (double)y_plus_mean);
        } else {
            /* Consecutive, so a sample back in range clears the count. */
            user->wall_yplus_excursions = 0;
        }

        /* The stop above reads the mean, which a local excursion - a corner, a separation,
           an outlet - cannot move. The maximum is reported once, with where it is, so such a
           site is visible without stopping a run whose wall is otherwise in range. */
        if (global_sum[0] > 0.0 && global_max[1] > upper && !user->wall_yplus_max_reported) {
            user->wall_yplus_max_reported = PETSC_TRUE;
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "Wall model (%s): first-cell y+ reaches %.2f on face %s at cell (%d, %d, %d), above %s "
                      "(up to %.0f), although the wall-face mean %.2f is in range. Reported once; "
                      "wall_model.csv tracks the maximum and its location every step.\n",
                      WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                      (double)global_max[1], BCFaceToString((BCFace)global_at[4]), (int)global_at[5],
                      (int)global_at[6], (int)global_at[7], range_reason, (double)upper, (double)y_plus_mean);
        }
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref LogBoundaryPressureDiagnostics().
 * @see LogBoundaryPressureDiagnostics()
 */
PetscErrorCode LogBoundaryPressureDiagnostics(UserCtx *user)
{
    SimCtx              *simCtx = user->simCtx;
    /* Per face: the area-weighted pressure sum and the open area. One collective for all six. */
    PetscReal            local[12] = {0.0}, global[12];
    PetscBool            selected[6] = {PETSC_FALSE};
    PetscInt             n_selected = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    for (PetscInt face = 0; face < 6; face++) {
        const BCType type = user->boundary_faces[face].mathematical_type;
        selected[face] = (PetscBool)(type == INLET || type == OUTLET);
        n_selected += selected[face] ? 1 : 0;
    }
    if (n_selected == 0) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }

    {
        const PetscReal ***pressure, ***nvert;
        const Cmpnts    ***metric[3];

        PetscCall(DMDAVecGetArrayRead(user->da, user->P, (void *)&pressure));
        PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (void *)&metric[0]));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (void *)&metric[1]));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (void *)&metric[2]));
        for (PetscInt face = 0; face < 6; face++) {
            BoundaryFaceElements elements;

            if (!selected[face]) continue;
            /* The first interior cell beside each owned face element, and the face's own metric. */
            PetscCall(GetBoundaryFaceElements(user, (BCFace)face, &elements));
            const PetscInt axis = elements.normal;
            const PetscInt cell = (face % 2 == 0) ? elements.layer + 1 : elements.layer;
            for (PetscInt q = 0; q < elements.n[1]; q++)
            for (PetscInt p = 0; p < elements.n[0]; p++) {
                PetscInt c[3], f[3];
                c[axis] = cell; c[elements.a] = elements.lo[0] + p; c[elements.b] = elements.lo[1] + q;
                const PetscInt i = c[0], j = c[1], k = c[2];
                if (nvert[k][j][i] > 0.1) continue;
                f[0] = i; f[1] = j; f[2] = k;
                f[axis] = elements.layer;
                const Cmpnts    m    = metric[axis][f[2]][f[1]][f[0]];
                const PetscReal area = PetscSqrtReal(m.x * m.x + m.y * m.y + m.z * m.z);
                local[2 * face]     += pressure[k][j][i] * area;
                local[2 * face + 1] += area;
            }
        }
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (void *)&metric[2]));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (void *)&metric[1]));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (void *)&metric[0]));
        PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
        PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, (void *)&pressure));
    }
    PetscCallMPI(MPI_Allreduce(local, global, 12, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));

    /* One record per run, as the inlet-flow record keeps: the first block's faces. */
    if (simCtx->rank == 0 && user->_this == 0) {
        const PetscReal L_ref = simCtx->scaling.L_ref, P_ref = simCtx->scaling.P_ref;
        PetscReal       physical_time = 0.0;
        FILE           *file = NULL;

        PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
        PetscCall(PicurvOpenDiagnosticsCsv(simCtx, simCtx->analysis_dir, "boundary_pressure.csv",
                                           "step,physical_time,face,type,area,mean_pressure", &file));
        for (PetscInt face = 0; face < 6; face++) {
            if (!selected[face] || global[2 * face + 1] <= 0.0) continue;
            fprintf(file, "%d,%.10e,%s,%s,%.10e,%.10e\n", (int)simCtx->step, (double)physical_time,
                    BCFaceToString((BCFace)face), BCTypeToString(user->boundary_faces[face].mathematical_type),
                    (double)(global[2 * face + 1] * L_ref * L_ref),
                    (double)(global[2 * face] / global[2 * face + 1] * P_ref));
        }
        PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close the boundary-pressure diagnostics file.");
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref FinalizePostProjectionCellFields().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 */
PetscErrorCode FinalizePostProjectionCellFields(UserCtx *user)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Finalizing post-projection cell-centered fields.\n");

    // Ensure flow-dependent Ubcs handlers see the newly reconstructed Ucat.
    ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = BoundarySystem_RefreshUbcs(user); CHKERRQ(ierr);
    ierr = FinalizeCellFieldsFromBoundaryValues(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Post-projection cell-centered fields finalized.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref FinalizeCellFieldsFromBoundaryValues().
 * @see FinalizeCellFieldsFromBoundaryValues()
 */
PetscErrorCode FinalizeCellFieldsFromBoundaryValues(UserCtx *user)
{
    PetscErrorCode ierr;
    const FieldId cell_fields[] = {FIELD_ID_UCAT, FIELD_ID_P};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // Establish flat non-periodic faces and periodic endpoints before corners.
    ierr = UpdateDummyCells(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = UpdateDummyCells(user, FIELD_ID_P); CHKERRQ(ierr);
    ierr = SynchronizePeriodicCellFields(user, 2, cell_fields); CHKERRQ(ierr);

    // Corner averaging can overwrite periodic endpoints, so restore them after.
    ierr = UpdateCornerNodes(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = UpdateCornerNodes(user, FIELD_ID_P); CHKERRQ(ierr);
    ierr = SynchronizePeriodicCellFields(user, 2, cell_fields); CHKERRQ(ierr);

    // With a periodic axis that synchronisation ends with the local vectors current: it
    // updates them after the corner writes, and each direction's plane refresh restores what
    // its transfer changed. With none it does nothing, so they are updated here.
    {
        PetscBool periodic[3];
        GetPeriodicAxes(user->simCtx, periodic);
        if (!periodic[0] && !periodic[1] && !periodic[2]) {
            ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
            ierr = UpdateLocalGhosts(user, FIELD_ID_P); CHKERRQ(ierr);
        } else {
            ierr = VerifyLocalGhostsCurrent(user, FIELD_ID_UCAT, "after the periodic synchronisation"); CHKERRQ(ierr);
            ierr = VerifyLocalGhostsCurrent(user, FIELD_ID_P, "after the periodic synchronisation"); CHKERRQ(ierr);
        }
    }

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief One boundary pass: the handlers, Contra2Cart, the wall function, the dummy, periodic
 *        and corner updates, and the ghost synchronisations that leave the fields consistent.
 * @details Local to this translation unit; ApplyBoundaryConditions() decides how many to run.
 */
static PetscErrorCode ApplyBoundaryPass(UserCtx *user, PetscBool any_periodic)
{
    PetscErrorCode ierr;
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    // (a) Execute the boundary system. This phase calculates fluxes across
    //     the domain and then applies the physical logic for each non-periodic
    //     handler, setting the `ubcs` (boundary value) array.
    ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Boundary Condition Setup Executed.\n");

    // (b) Synchronize the updated ghost cells across all processors to ensure
    //     all ucont values are current before updating the dummy cells.
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);

    // (c) Convert updated Contravariant velocities to Cartesian velocities.
    ierr = Contra2Cart(user); CHKERRQ(ierr);

    // (d) Ucat's ghosts are not updated here: the wall function, the dummy and corner
    //     updates read only global Ucat, and (g) or (i) updates it before any local read.
    ierr = PoisonStaleLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);

    // (e) If Wall functions are enabled, apply them now to adjust near-wall velocities.
    if(user->simCtx->wallfunction){
      // Apply wall function adjustments to the boundary velocities.
      ierr = ApplyWallFunction(user); CHKERRQ(ierr);

      LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Wall Function Applied at Walls.\n");
    }

    // (f) Update the first layer of ghost cells for non-periodic faces using
    //     the newly computed `ubcs` values.
    ierr = UpdateDummyCells(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = UpdateDummyCells(user, FIELD_ID_P); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Dummy Cells/Ghost Cells Updated.\n");

    // (g) Handle all periodic boundaries. This is a parallel direct copy
    // that sets the absolute constraints for the rest of the solve.
    // There is a Ghost update happening inside this function. Its last cell-field update
    // is left to (j): (h) and the staggered transfers do not read local Ucat, P or Nvert.
    ierr = ApplyPeriodicBCs(user, PETSC_FALSE); CHKERRQ(ierr);

    // (h) Update the corner and edge ghost nodes. This routine calculates
    // values for corners/edges by averaging their neighbors, which have been
    // finalized in the steps above (both periodic and non-periodic).
    ierr = UpdateCornerNodes(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = UpdateCornerNodes(user, FIELD_ID_P); CHKERRQ(ierr);
    
    // (i) Synchronize the updated edge and corner cells across all processors to ensure
    //     consistency before the next iteration or finalization.
    //     Ucont was last written before (g)'s synchronization, or, with no periodic
    //     face, before (b)'s update unless the wall function wrote it since. With a
    //     periodic face, (j) updates P and Ucat first, so they are not updated here.
    if (!any_periodic) {
        ierr = UpdateLocalGhosts(user, FIELD_ID_P); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    }
    ierr = SynchronizePeriodicStaggeredFieldsFrom(user, 1, staggered_fields,
               (PetscBool)(any_periodic || !user->simCtx->wallfunction)); CHKERRQ(ierr);

    // (j) Ensure All the corners are synchronized with  a well defined protocol  in case of Periodic boundary conditions
    // To avoid race conditions.
    //     Its leading update also completes (g)'s, which was left to it.
    const FieldId all_fields[] = {FIELD_ID_UCAT, FIELD_ID_P, FIELD_ID_NVERT};
    ierr = SynchronizePeriodicCellFieldsFrom(user, 3, all_fields, PETSC_FALSE, PETSC_TRUE); CHKERRQ(ierr);

    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief In verify mode, checks that the passes a single-pass call skipped would change nothing.
 * @details Runs the two passes the coupled path would add and requires Ucont, Ucat, P and the
 *          boundary values to be bitwise what the first pass left. A mismatch means a handler
 *          that reads in-pass state does not set `needs_repeat_pass`. Local to this
 *          translation unit.
 */
static PetscErrorCode VerifySinglePassIsFinal(UserCtx *user, PetscBool any_periodic)
{
    Vec         fields[4] = {user->Ucont, user->Ucat, user->P, user->Bcs.Ubcs}, saved[4];
    const char *names[4] = {"Ucont", "Ucat", "P", "Ubcs"};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!user->simCtx->verifyGhostUpdates) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    for (int f = 0; f < 4; ++f) {
        PetscCall(VecDuplicate(fields[f], &saved[f]));
        PetscCall(VecCopy(fields[f], saved[f]));
    }
    for (int pass = 1; pass < 3; ++pass) PetscCall(ApplyBoundaryPass(user, any_periodic));
    for (int f = 0; f < 4; ++f) {
        PetscBool equal;
        PetscCall(VecEqual(fields[f], saved[f], &equal));
        PetscCheck(equal, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "A further boundary pass changed %s after a single-pass call: some face's handler reads "
                   "state the pass changes but does not set needs_repeat_pass.", names[f]);
        PetscCall(VecDestroy(&saved[f]));
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ApplyBoundaryConditions().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see ApplyBoundaryConditions()
 */
PetscErrorCode ApplyBoundaryConditions(UserCtx *user)
{
    PetscErrorCode ierr;
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};
    PetscBool periodic[3], any_periodic;
    PetscInt  passes = user->simCtx->wallfunction ? 3 : 1;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL,LOG_TRACE,"Boundary Condition Application begins.\n");

    GetPeriodicAxes(user->simCtx, periodic);
    any_periodic = (PetscBool)(periodic[0] || periodic[1] || periodic[2]);
    for (int face = 0; face < 6; face++) {
        if (user->boundary_faces[face].handler && user->boundary_faces[face].handler->needs_repeat_pass) passes = 3;
    }
    // STEP 1: The boundary passes. Three let information propagate between coupled boundaries,
    // like an inlet and a conserving outlet, or a free-slip wall and the interior velocity it
    // copies. Without such a face (needs_repeat_pass) and without a wall function, the first
    // pass leaves every field final, and verify mode checks that the skipped passes agree.
    for (PetscInt iter = 0; iter < passes; iter++) {
        ierr = ApplyBoundaryPass(user, any_periodic); CHKERRQ(ierr);
    }
    if (passes == 1) {
        ierr = VerifySinglePassIsFinal(user, any_periodic); CHKERRQ(ierr);
    }

    // STEP 2: Final staggered periodic synchronization. The last pass already left Ucat, P
    // and Ucont's local vectors current with their global ones, so only the transfers run: a
    // ghost update depends only on its global vector, and one that repeats an update with no
    // write to that vector since changes nothing.
    ierr = SynchronizePeriodicStaggeredFieldsFrom(user, 1, staggered_fields, PETSC_TRUE); CHKERRQ(ierr);
    ierr = VerifyBoundaryPassFinite(user); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
