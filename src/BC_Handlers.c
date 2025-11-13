#include "BC_Handlers.h"     // The header that declares this file's "constructor" functions


//================================================================================
//              VALIDATORS
//================================================================================


#undef __FUNCT__
#define __FUNCT__ "Validate_DrivenFlowConfiguration"
/**
 * @brief (Private) Validates all consistency rules for a driven flow (channel/pipe) setup.
 *
 * This function enforces a strict set of rules to ensure a driven flow simulation is
 * configured correctly. It is called by the main `BoundarySystem_Validate` dispatcher.
 *
 * The validation rules are checked in a specific order:
 * 1.  Detect if any `DRIVEN_` handler is active. If not, the function returns immediately.
 * 2.  Ensure that no `INLET`, `OUTLET`, or `FARFIELD` boundary conditions exist anywhere in the
 *     domain, as they are physically incompatible with a pressure-driven flow model.
 * 3.  Verify that both faces in the driven direction are of `mathematical_type PERIODIC`.
 * 4.  Verify that both faces in the driven direction use the exact same `DRIVEN_` handler type.
 *
 * @param user The UserCtx for a single block.
 * @return PetscErrorCode 0 on success, non-zero PETSc error code on failure.
 */
PetscErrorCode Validate_DrivenFlowConfiguration(UserCtx *user)
{
    PetscFunctionBeginUser;

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

    PetscFunctionReturn(0);
}

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

#undef __FUNCT__
#define __FUNCT__ "Create_WallNoSlip"
/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with No-Slip Wall behavior.
 *
 * A no-slip wall is simple and requires only the `Apply` method:
 *  - No special initialization needed (`Initialize` is NULL)
 *  - Does not contribute to global mass balance (`PreStep` and `PostStep` are NULL)
 *  - Enforces zero velocity at the wall (`Apply`)
 *  - Allocates no private data (`Destroy` is NULL)
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_WallNoSlip(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
                     "Input BoundaryCondition object is NULL in Create_WallNoSlip");

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

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Apply_WallNoSlip"
/**
 * @brief (Handler Action) Applies the no-slip wall condition to a specified face.
 *
 * This function enforces zero velocity at the wall by:
 *   1. Setting contravariant flux (ucont) to zero (no penetration)
 *   2. Setting boundary velocity (ubcs) to zero (no slip)
 *   
 * NOTE: Unlike the legacy code, this does NOT set ucat[ghost]. The orchestrator's
 * UpdateDummyCells function handles ghost cell extrapolation uniformly for all BC types.
 *
 * @param self The BoundaryCondition object (unused for this simple handler)
 * @param ctx  BCContext containing UserCtx and face_id
 * @return PetscErrorCode 0 on success
 */
static PetscErrorCode Apply_WallNoSlip(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    PetscBool      can_service;
    
    (void)self;  // Unused for simple handlers
    
    PetscFunctionBeginUser;
    DMDALocalInfo *info = &user->info;
    Cmpnts ***ubcs, ***ucont;
    PetscInt IM_nodes_global, JM_nodes_global,KM_nodes_global;

    IM_nodes_global = user->IM;
    JM_nodes_global = user->JM;
    KM_nodes_global = user->KM;
    
    ierr = CanRankServiceFace(info,IM_nodes_global,JM_nodes_global,KM_nodes_global,face_id,&can_service); CHKERRQ(ierr);
    // Check if this rank owns part of this boundary face
    if (!can_service) PetscFunctionReturn(0);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply_WallNoSlip: Applying to Face %d (%s).\n",
              face_id, BCFaceToString(face_id));

    // Get arrays
    
    ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);

    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;
    PetscInt mx = info->mx, my = info->my, mz = info->mz;
    
    // ✅ Use shrunken loop bounds (avoids edges/corners like inlet handler)
    PetscInt lxs = xs, lxe = xe, lys = ys, lye = ye, lzs = zs, lze = ze;
    if (xs == 0) lxs = xs + 1;
    if (xe == mx) lxe = xe - 1;
    if (ys == 0) lys = ys + 1;
    if (ye == my) lye = ye - 1;
    if (zs == 0) lzs = zs + 1;
    if (ze == mz) lze = ze - 1;

    switch (face_id) {
        case BC_FACE_NEG_X: {
            if (xs == 0){  
                PetscInt i = xs;
                for (PetscInt k = lzs; k < lze; k++) {
                    for (PetscInt j = lys; j < lye; j++) {
                        // ✅ Set contravariant flux to zero (no penetration)
                        ucont[k][j][i].x = 0.0;
                    
                        // ✅ Set boundary velocity to zero (no slip)
                        ubcs[k][j][i].x = 0.0;
                        ubcs[k][j][i].y = 0.0;
                        ubcs[k][j][i].z = 0.0;
                        }                
                    }
                }
            break;
        }

        case BC_FACE_POS_X: {
            if (xe == mx){            
                PetscInt i = xe - 1;
                for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    ucont[k][j][i-1].x = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                    }
                }    
            }
        break;
        }    
        case BC_FACE_NEG_Y: {
            if (ys == 0){
            PetscInt j = ys;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k][j][i].y = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                    }
                }
            }
        } break;

        case BC_FACE_POS_Y: {
            if (ye == my){
            PetscInt j = ye - 1;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k][j-1][i].y = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                    }
                }
            }
        } break;
            
        case BC_FACE_NEG_Z: {
            if (zs == 0){
            PetscInt k = zs;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k][j][i].z = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                    }
                }
            }
        } break;

        case BC_FACE_POS_Z: {
            if (ze == mz) {
            PetscInt k = ze - 1;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k-1][j][i].z = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                    }
                }
            }    
        } break;
    }

    // Restore arrays
    ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////

//================================================================================
//
//          HANDLER IMPLEMENTATION: CONSTANT VELOCITY INLET
//          (Corresponds to BC_HANDLER_INLET_CONSTANT_VELOCITY)
//
//================================================================================

// --- FORWARD DECLARATIONS ---
static PetscErrorCode Initialize_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx, 
                                                     PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PostStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx,
                                                      PetscReal *in, PetscReal *out);
static PetscErrorCode Destroy_InletConstantVelocity(BoundaryCondition *self);

/**
 * @brief Private data structure for the Constant Velocity Inlet handler.
 */
typedef struct{
    PetscReal normal_velocity; // The desired Cartesian velocity (vx, vy, vz)
}InletConstantData;

#undef __FUNCT__
#define __FUNCT__ "Create_InletConstantVelocity"
/*
* @brief (Handler Constructor) Populates a BoundaryCondition object with Constant Velocity Inlet behavior.
*/
PetscErrorCode Create_InletConstantVelocity(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "BoundaryCondition is NULL");

    InletConstantData *data = NULL;
    ierr = PetscMalloc1(1, &data); CHKERRQ(ierr);
    bc->data = (void*)data;
    
    bc->priority   = BC_PRIORITY_INLET;
    bc->Initialize = Initialize_InletConstantVelocity;
    bc->PreStep    = PreStep_InletConstantVelocity;
    bc->Apply      = Apply_InletConstantVelocity;
    bc->PostStep   = PostStep_InletConstantVelocity;
    bc->UpdateUbcs = NULL;
    bc->Destroy    = Destroy_InletConstantVelocity;
    
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "Initialize_InletConstantVelocity"
/**
 * @brief (Handler Action) Initializes the constant velocity inlet handler.
 *
 * Parses the appropriate velocity component ('vx', 'vy', or 'vz') from bcs.dat
 * based on the face orientation and sets the initial state on the boundary face
 * by calling the Apply function.
 */
static PetscErrorCode Initialize_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletConstantData *data = (InletConstantData*)self->data;
    PetscBool      found;
    
    PetscFunctionBeginUser;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initialize_InletConstantVelocity: Initializing handler for Face %d. \n", face_id);
    data->normal_velocity = 0.0;
    
    switch (face_id) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X:
            // For X-faces, read "vx" as normal velocity
            ierr = GetBCParamReal(user->boundary_faces[face_id].params, "vx", 
                                 &data->normal_velocity, &found); CHKERRQ(ierr);
            break;
            
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y:
            // For Y-faces, read "vy" as normal velocity
            ierr = GetBCParamReal(user->boundary_faces[face_id].params, "vy", 
                                 &data->normal_velocity, &found); CHKERRQ(ierr);
            break;
            
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z:
            // For Z-faces, read "vz" as normal velocity
            ierr = GetBCParamReal(user->boundary_faces[face_id].params, "vz", 
                                 &data->normal_velocity, &found); CHKERRQ(ierr);
            break;
    }
    
    LOG_ALLOW(LOCAL, LOG_INFO, "  Inlet Face %d: normal velocity = %.4f\n",
              face_id, data->normal_velocity);

    // Set initial boundary state
    ierr = Apply_InletConstantVelocity(self, ctx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PreStep_InletConstantVelocity"
/**
 * @brief (Handler PreStep) No preparation needed for constant velocity inlet.
 *
 * For constant velocity inlets, all parameters are parsed during Initialize
 * and stored in the handler's private data. There is nothing to prepare before
 * Apply, so this function is a no-op.
 *
 * NOTE: Other inlet types (parabolic, time-varying, file-based) DO use PreStep
 * for profile calculation, file I/O, or data preparation.
 */
static PetscErrorCode PreStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx,
                                                     PetscReal *local_inflow_contribution,
                                                     PetscReal *local_outflow_contribution)
{
    // No preparation needed for constant velocity inlet.
    // The velocity is already stored in self->data from Initialize.
    // Apply will set ucont, and PostStep will measure the actual flux.
    
    (void)self;
    (void)ctx;
    (void)local_inflow_contribution;
    (void)local_outflow_contribution;
    
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Apply_InletConstantVelocity"
/**
 * @brief (Handler Action) Applies the constant velocity inlet condition.
 *
 * This function enforces a constant normal velocity on its assigned face by:
 * 1.  Iterating over all owned nodes on the specified boundary face.
 * 2.  For each valid fluid node (not covered by a solid), setting the contravariant
 *     velocity component (Ucont) to achieve the desired normal velocity.
 * 3.  Calculating and setting the corresponding Cartesian velocity components (Ubc)
 *     in the ghost cells adjacent to the boundary face.
 */
static PetscErrorCode Apply_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletConstantData *data = (InletConstantData*)self->data;
    PetscBool      can_service;
    
    PetscFunctionBeginUser;
    
    DMDALocalInfo *info = &user->info;
    Cmpnts ***ubcs, ***ucont, ***csi, ***eta, ***zet;
    PetscReal ***nvert;
    PetscInt IM_nodes_global, JM_nodes_global,KM_nodes_global;

    IM_nodes_global = user->IM;
    JM_nodes_global = user->JM;
    KM_nodes_global = user->KM;
    
    ierr = CanRankServiceFace(info,IM_nodes_global,JM_nodes_global,KM_nodes_global,face_id,&can_service); CHKERRQ(ierr);
    
    if (!can_service) PetscFunctionReturn(0);

    // Get arrays

    ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // Get SCALAR velocity (not vector!)
    PetscReal uin_this_point = data->normal_velocity;
    
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

    switch (face_id) {
        case BC_FACE_NEG_X:
        case BC_FACE_POS_X: {
            PetscReal sign = (face_id == BC_FACE_NEG_X) ? 1.0 : -1.0;
            PetscInt i = (face_id == BC_FACE_NEG_X) ? xs : mx - 2;
            
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt j = lys; j < lye; j++) {
                    if ((sign > 0 && nvert[k][j][i+1] > 0.1) || 
                        (sign < 0 && nvert[k][j][i] > 0.1)) continue;
                    
                    PetscReal CellArea = sqrt(csi[k][j][i].x * csi[k][j][i].x + 
                                             csi[k][j][i].y * csi[k][j][i].y + 
                                             csi[k][j][i].z * csi[k][j][i].z);
                    
                    ucont[k][j][i].x = sign * uin_this_point * CellArea;
                    
                    ubcs[k][j][i + (sign < 0)].x = sign * uin_this_point * csi[k][j][i].x / CellArea;
                    ubcs[k][j][i + (sign < 0)].y = sign * uin_this_point * csi[k][j][i].y / CellArea;
                    ubcs[k][j][i + (sign < 0)].z = sign * uin_this_point * csi[k][j][i].z / CellArea;
                }
            }
        } break;
            
        case BC_FACE_NEG_Y:
        case BC_FACE_POS_Y: {
            PetscReal sign = (face_id == BC_FACE_NEG_Y) ? 1.0 : -1.0;
            PetscInt j = (face_id == BC_FACE_NEG_Y) ? ys : my - 2;
            
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if ((sign > 0 && nvert[k][j+1][i] > 0.1) || 
                        (sign < 0 && nvert[k][j][i] > 0.1)) continue;
                    
                    PetscReal CellArea = sqrt(eta[k][j][i].x * eta[k][j][i].x + 
                                             eta[k][j][i].y * eta[k][j][i].y + 
                                             eta[k][j][i].z * eta[k][j][i].z);
                    
                    ucont[k][j][i].y = sign * uin_this_point * CellArea;
                    
                    ubcs[k][j + (sign < 0)][i].x = sign * uin_this_point * eta[k][j][i].x / CellArea;
                    ubcs[k][j + (sign < 0)][i].y = sign * uin_this_point * eta[k][j][i].y / CellArea;
                    ubcs[k][j + (sign < 0)][i].z = sign * uin_this_point * eta[k][j][i].z / CellArea;
                }
            }
        } break;
            
        case BC_FACE_NEG_Z:
        case BC_FACE_POS_Z: {
            PetscReal sign = (face_id == BC_FACE_NEG_Z) ? 1.0 : -1.0;
            PetscInt k = (face_id == BC_FACE_NEG_Z) ? zs : mz - 2;
            
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    if ((sign > 0 && nvert[k+1][j][i] > 0.1) || 
                        (sign < 0 && nvert[k][j][i] > 0.1)) continue;
                    
                    PetscReal CellArea = sqrt(zet[k][j][i].x * zet[k][j][i].x + 
                                             zet[k][j][i].y * zet[k][j][i].y + 
                                             zet[k][j][i].z * zet[k][j][i].z);
                    
                    ucont[k][j][i].z = sign * uin_this_point * CellArea;
                    
                    ubcs[k + (sign < 0)][j][i].x = sign * uin_this_point * zet[k][j][i].x / CellArea;
                    ubcs[k + (sign < 0)][j][i].y = sign * uin_this_point * zet[k][j][i].y / CellArea;
                    ubcs[k + (sign < 0)][j][i].z = sign * uin_this_point * zet[k][j][i].z / CellArea;
                }
            }
        } break;
    }

    // Restore arrays
    ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PostStep_InletConstantVelocity"
/**
 * @brief (Handler PostStep) Measures actual inflow flux through the constant velocity inlet face.
 */
static PetscErrorCode PostStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx,
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

    DMDALocalInfo *info = &user->info;
    Cmpnts ***ucont;

    PetscInt IM_nodes_global, JM_nodes_global,KM_nodes_global;

    IM_nodes_global = user->IM;
    JM_nodes_global = user->JM;
    KM_nodes_global = user->KM;
    
    ierr = CanRankServiceFace(info,IM_nodes_global,JM_nodes_global,KM_nodes_global,face_id,&can_service); CHKERRQ(ierr);
    

    if (!can_service) PetscFunctionReturn(0);
    
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
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "PostStep_InletConstantVelocity: Face %d, flux = %.6e\n",
              face_id, local_flux);

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "Destroy_InletConstantVelocity"
/**
 * @brief (Handler Destructor) Frees memory for the Constant Velocity Inlet.
 */
static PetscErrorCode Destroy_InletConstantVelocity(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    if (self && self->data) {
        PetscFree(self->data);
        self->data = NULL;
    }
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
//================================================================================

// --- 1. FORWARD DECLARATIONS ---
static PetscErrorCode PreStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                     PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution);
static PetscErrorCode Apply_OutletConservation(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PostStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                                  PetscReal *in, PetscReal *out);

#undef __FUNCT__
#define __FUNCT__ "Create_OutletConservation"
/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with Outlet Conservation behavior.
 */
PetscErrorCode Create_OutletConservation(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition is NULL");

    // This handler has the highest priority to ensure it runs after
    // all inflow fluxes have been calculated.
    bc->priority   = BC_PRIORITY_OUTLET;
    
    // Assign function pointers
    bc->Initialize = NULL; // No initialization needed
    bc->PreStep    = PreStep_OutletConservation;
    bc->Apply      = Apply_OutletConservation;
    bc->PostStep   = PostStep_OutletConservation;
    bc->UpdateUbcs = NULL;
    bc->Destroy    = NULL; // No private data to destroy
    
    bc->data = NULL;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PreStep_OutletConservation"                                     
/**
 * @brief (Handler Action) Measures the current, uncorrected flux passing through a SINGLE outlet face.
 *
 * This function is called during the PreStep phase of the boundary condition cycle. It is a direct,
 * high-fidelity port of the flux measurement logic from the legacy function for
 * outlet-type boundaries.
 *
 * Its primary responsibility is to calculate the total volumetric flux that is currently passing
 * *out* of the domain through the specified outlet face, before any mass-conservation corrections
 * are applied.
 *
 * The calculation is performed by taking the dot product of the interior cell-centered Cartesian
 * velocity (`ucat`) with the corresponding face area vectors (`csi`, `eta`, or `zet`). To ensure
 * bit-for-bit identical behavior with the legacy code, this function respects the convention of
 * excluding domain corners and edges from the main face loops by using shrunk loop bounds
 * (`lxs`, `lxe`, etc.).
 *
 * The result from this rank's portion of the face is added to the `local_outflow_contribution`
 * accumulator, which is later summed across all MPI ranks to obtain the global uncorrected outflow.
 *
 * @param self A pointer to the BoundaryCondition object (unused in this specific handler).
 * @param ctx  The context for this execution, providing access to the `UserCtx` and, critically,
 *             the `face_id` that this function call should operate on.
 * @param local_inflow_contribution  Accumulator for inflow flux (unused by this handler).
 * @param[out] local_outflow_contribution Accumulator for outflow flux. This function will ADD its
 *                                        calculated flux to this value.
 * @return PetscErrorCode 0 on success.
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

    // Step 1: Use the robust utility function to determine if this MPI rank owns a computable
    // portion of the specified boundary face. If not, there is no work to do, so we exit immediately.
    const PetscInt IM_nodes_global = user->IM;
    const PetscInt JM_nodes_global = user->JM;
    const PetscInt KM_nodes_global = user->KM;
    ierr = CanRankServiceFace(info, IM_nodes_global, JM_nodes_global, KM_nodes_global, face_id, &can_service); CHKERRQ(ierr);

    if (!can_service) {
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

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Apply_OutletConservation"
/**
 * @brief (Handler Action) Applies mass conservation correction to the outlet face.
 *
 * This function calculates a global correction factor based on the total inflow and outflow fluxes
 * and applies it to the contravariant velocity (`ucont`) on the outlet face to ensure mass conservation.
 */static PetscErrorCode Apply_OutletConservation(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    DMDALocalInfo* info = &user->info;
    PetscBool      can_service;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    const PetscInt IM_nodes_global = user->IM;
    const PetscInt JM_nodes_global = user->JM;
    const PetscInt KM_nodes_global = user->KM;
    ierr = CanRankServiceFace(info, IM_nodes_global, JM_nodes_global, KM_nodes_global, face_id, &can_service); CHKERRQ(ierr);
    
    if (!can_service) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // --- STEP 1: Calculate the correction factor using pre-calculated area ---
    PetscReal total_inflow = *ctx->global_inflow_sum + *ctx->global_farfield_inflow_sum;
    PetscReal flux_imbalance = total_inflow - *ctx->global_outflow_sum;
    
    // Directly use the pre-calculated area from the simulation context.
    PetscReal velocity_correction = (PetscAbsReal(user->simCtx->AreaOutSum) > 1e-12) 
                                    ? flux_imbalance / user->simCtx->AreaOutSum 
                                    : 0.0;
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Outlet Correction on Face %d: Imbalance=%.4e, Pre-calc Area=%.4e, V_corr=%.4e\n",
              face_id, flux_imbalance, user->simCtx->AreaOutSum, velocity_correction);

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

#undef __FUNCT__
#define __FUNCT__ "PostStep_OutletConservation"
/**
 * @brief (Handler PostStep) Measures corrected outflow flux for verification.
 *
 * After Apply has set the corrected ucont values, this function measures
 * the actual flux passing through the outlet face by summing the ucont
 * components. Only fluid cells (nvert < 0.1) are included in the sum.
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
    const PetscInt IM_nodes_global = user->IM;
    const PetscInt JM_nodes_global = user->JM;
    const PetscInt KM_nodes_global = user->KM;
    ierr = CanRankServiceFace(info, IM_nodes_global, JM_nodes_global, KM_nodes_global, face_id, &can_service); CHKERRQ(ierr);

    if (!can_service) PetscFunctionReturn(0);

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

    PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///////// Geometric Periodic BC Handler
/////////////////////////////////////////////////////////////////////////////////////////////////

PetscErrorCode Create_PeriodicGeometric(BoundaryCondition *bc){
    PetscFunctionBeginUser;
    
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition is NULL");
    bc->priority = BC_PRIORITY_WALL;

    // Assign function pointers
    bc->Initialize = NULL; // No initialization needed
    bc->PreStep    = NULL;
    bc->Apply      = NULL;
    bc->PostStep   = NULL;
    bc->UpdateUbcs = NULL;
    bc->Destroy    = NULL; // No private data to destroy

    bc->data = NULL;

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

/** @brief Private data structure for the handler. */
typedef struct {
    char      direction;                // 'X', 'Y', or 'Z', determined at initialization.
    PetscReal targetVolumetricFlux;     // The constant target flux, parsed from parameters.
    PetscReal boundaryVelocityCorrection; // The "Boundary Trim" value for the Apply() step.
    PetscBool isMasterController;       // Flag: PETSC_TRUE only for the handler on the negative face.
    PetscBool applyBoundaryTrim;        // Flag: PETSC_TRUE if applying Boundary trim on ucont.
} DrivenConstantData;


// --- 2. HANDLER CONSTRUCTOR ---

#undef __FUNCT__
#define __FUNCT__ "Create_PeriodicDrivenConstant"
/**
 * @brief (Handler Constructor) Creates and configures a BoundaryCondition object for a driven periodic
 *        flow with a constant, user-defined target flux.
 *
 * This function acts as the factory entry point for the `BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX` type.
 * It performs the following steps:
 * 1.  Allocates memory for the generic `BoundaryCondition` object (done by the factory caller).
 * 2.  Allocates memory for its own private `DrivenConstantData` struct.
 * 3.  Sets the execution priority to `BC_PRIORITY_INLET` to ensure the controller's `PreStep`
 *     runs before other flux-measuring handlers.
 * 4.  Assigns the function pointers (`Initialize`, `PreStep`, `Apply`, `Destroy`) to the
 *     specific static implementations defined in this file. Other methods are set to NULL
 *     as they are not needed by this handler.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_PeriodicDrivenConstant(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition object is NULL in Create_PeriodicDrivenConstantFlux");

    // --- Allocate the private data structure ---
    DrivenConstantData *data = NULL;
    ierr = PetscNew(&data); CHKERRQ(ierr);
    // Initialize fields to safe default values
    data->direction = ' ';
    data->targetVolumetricFlux = 0.0;
    data->boundaryVelocityCorrection = 0.0;
    data->isMasterController = PETSC_FALSE;
    data->applyBoundaryTrim = PETSC_FALSE;
    
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
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Initialize_PeriodicDrivenConstant"
/**
 * @brief (Handler Initialize) Initializes the handler by validating its configuration and parsing parameters.
 *
 * This method is called once at simulation startup. It performs the following critical setup tasks:
 *
 * 1.  **Validation:** It verifies that this handler has been correctly applied to a face with
 *     `mathematical_type = PERIODIC`. If not, it halts the simulation with a clear error message.
 *
 * 2.  **Role Assignment:** It determines its operational direction ('X', 'Y', or 'Z') based on
 *     the `face_id` it is attached to. It also designates the handler on the "negative" face
 *     (e.g., BC_FACE_NEG_X) as the "master controller". This ensures that computationally
 *     expensive, domain-wide calculations in the `PreStep` method are only executed once
 *     per direction.
 *
 * 3.  **Parameter Parsing:** If this instance is the master controller, it parses the required
 *     `target_flux` parameter from the boundary condition configuration file. It will halt
 *     with an error if this parameter is missing. The parsed value is stored in the handler's
 *     private data and also copied to a shared location in the `UserCtx` for other parts of
 *     the solver (like the "Enforcer" function) to access.
 *
 * @param self The `BoundaryCondition` object containing the handler's state.
 * @param ctx  The `BCContext`, providing access to the `UserCtx` and `face_id`.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Initialize_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    DrivenConstantData *data = (DrivenConstantData*)self->data;
    BCFace face_id = ctx->face_id;
    UserCtx* user = ctx->user;
    
    PetscFunctionBeginUser;

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

        // Store the target flux in the UserCtx. This makes it globally accessible
        // to other parts of the solver, such as the `CorrectChannelFluxProfile` enforcer function.
        user->simCtx->targetVolumetricFlux = data->targetVolumetricFlux;
    }

    PetscBool trimfound;
    // Attempt  to read Trim flag.
    ierr = GetBCParamBool(user->boundary_faces[face_id].params, "apply_trim",
                        &data->applyBoundaryTrim, &trimfound); CHKERRQ(ierr);
    
    if(!trimfound) LOG_ALLOW(GLOBAL,LOG_DEBUG,"Trim Application not found,defaults to %s.\n",data->applyBoundaryTrim? "True":"False");

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PreStep_PeriodicDrivenConstant"
/**
 * @brief (Handler PreStep) Measures current fluxes and calculates the correction terms for the timestep.
 *
 * This method executes the main "Strategist" logic of the driven flow controller. It is
 * called once per timestep during the `PreStep` phase of the boundary condition cycle.
 *
 * The logic is executed *only* by the "master" handler (the one on the negative face)
 * to ensure domain-wide calculations are performed just once per direction.
 *
 * The function performs the following steps:
 * 1.  Reads the globally-set `targetVolumetricFlux` from the `SimCtx`.
 * 2.  Measures the `globalAveragePlanarVolumetricFlux` by averaging the flux across all
 *     cross-sectional planes. This provides a stable, noise-filtered signal.
 * 3.  Measures the `globalCurrentBoundaryFlux` at the single periodic boundary plane.
 *     This provides a fast, responsive signal.
 * 4.  Calculates `bulkVelocityCorrection` using the stable, averaged flux and stores it in
 *     the `SimCtx` for the `ApplyDrivenChannelFlowSource` function to use.
 * 5.  Calculates `boundaryVelocityCorrection` using the fast, boundary-specific flux and
 *     stores it in the handler's private data for its own `Apply` method to use for the
 *     "Boundary Trim".
 *
 * @param self The `BoundaryCondition` object containing the handler's state.
 * @param ctx  The `BCContext`, providing access to the `UserCtx`.
 * @param local_inflow_contribution  (Not used by this handler)
 * @param local_outflow_contribution (Not used by this handler)
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode PreStep_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx,
                                                     PetscReal *local_inflow_contribution,
                                                     PetscReal *local_outflow_contribution)
{
    PetscErrorCode ierr;
    DrivenConstantData *data = (DrivenConstantData*)self->data;
    UserCtx* user = ctx->user;
    SimCtx* simCtx = user->simCtx;

    PetscFunctionBeginUser;

    // --- Master Check: Only the handler on the negative face performs calculations ---
    if (!data->isMasterController) {
        PetscFunctionReturn(0);
    }
    
    // This function is the generalized implementation of the old InitializeChannelFlowController.
    char direction = data->direction;
    DMDALocalInfo info = user->info;
    PetscInt i, j, k;

    // --- Get read-only access to necessary field data ---
    Cmpnts ***ucont;
    PetscReal ***nvert;
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcont, (const Cmpnts***)&ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // --- Define local loop bounds ---
    PetscInt lxs = (info.xs == 0) ? 1 : info.xs;
    PetscInt lys = (info.ys == 0) ? 1 : info.ys;
    PetscInt lzs = (info.zs == 0) ? 1 : info.zs;
    PetscInt lxe = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    PetscInt lye = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    PetscInt lze = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    // ===================================================================================
    // CONTROLLER SENSORS: DUAL-FLUX MEASUREMENT STRATEGY
    //
    // To create a control system that is both stable and responsive, we measure the
    // volumetric flux in two distinct ways. One provides a fast, local signal for
    // immediate corrections, while the other provides a slow, stable signal for
    // strategic, domain-wide adjustments.
    //
    // -----------------------------------------------------------------------------------
    // ** 1. Fast/Local Sensor: `globalCurrentBoundaryFlux` **
    // -----------------------------------------------------------------------------------
    //
    // - WHAT IT IS: The total volumetric flux (m³/s) passing through the single,
    //   representative periodic boundary plane at k=0.
    //
    // - PHYSICAL MEANING: An instantaneous "snapshot" of the flow rate at the
    //   critical "seam" where the domain wraps around.
    //
    // - CHARACTERISTICS:
    //     - Fast & Responsive: It immediately reflects any local flow changes or
    //       errors occurring at the periodic interface.
    //     - Noisy: It is susceptible to local fluctuations from turbulence or
    //       numerical artifacts, causing its value to jitter from step to step.
    //
    // - CONTROLLER'S USE: This measurement is used to compute the
    //   `boundaryVelocityCorrection`. This is a TACTICAL, fast-acting "trim" that is
    //   applied immediately and directly to the boundary velocities to ensure perfect
    //   continuity at the seam.
    //
    // -----------------------------------------------------------------------------------
    // ** 2. Stable/Global Sensor: `globalAveragePlanarVolumetricFlux` **
    // -----------------------------------------------------------------------------------
    //
    // - WHAT IT IS: The average of the volumetric fluxes across ALL cross-sectional
    //   k-planes in the domain.
    //   (i.e., [Flux(k=0) + Flux(k=1) + ... + Flux(k=N)] / N)
    //
    // - PHYSICAL MEANING: It represents the overall, bulk momentum of the fluid,
    //   effectively filtering out local, transient fluctuations.
    //
    // - CHARACTERISTICS:
    //     - Stable & Robust: Local noise at one plane is averaged out by all other
    //       planes, providing a very smooth signal.
    //     - Inertial: It changes more slowly, reflecting the inertia of the entire
    //       fluid volume.
    //
    // - CONTROLLER'S USE: This measurement is used to compute the
    //   `bulkVelocityCorrection`. This is a STRATEGIC, long-term adjustment used to
    //   scale the main momentum source (body force). Using this stable signal prevents
    //   the main driving force from oscillating wildly, ensuring simulation stability.
    //
    // ===================================================================================

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
    }

    // --- Release array access as soon as possible ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcont, (const Cmpnts***)&ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);

    // --- Get cross-sectional area using the dedicated geometry function ---
    Cmpnts ignored_center;
    PetscReal globalBoundaryArea;
    BCFace neg_face_id = (direction == 'X') ? BC_FACE_NEG_X : (direction == 'Y') ? BC_FACE_NEG_Y : BC_FACE_NEG_Z;
    ierr = CalculateFaceCenterAndArea(user, neg_face_id, &ignored_center, &globalBoundaryArea); CHKERRQ(ierr);

    // --- Perform global reductions to get the final flux values ---
    PetscReal globalCurrentBoundaryFlux, globalAveragePlanarVolumetricFlux;
    ierr = MPI_Allreduce(&localCurrentBoundaryFlux, &globalCurrentBoundaryFlux, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&localAveragePlanarVolumetricFluxTerm, &globalAveragePlanarVolumetricFlux, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    // --- Calculate the two correction terms ---
    if (globalBoundaryArea > 1.0e-12) {
        // The main correction for the body force is based on the STABLE domain-averaged flux.
        simCtx->bulkVelocityCorrection = (data->targetVolumetricFlux - globalAveragePlanarVolumetricFlux) / globalBoundaryArea;
        
        // The immediate correction for the boundary trim is based on the FAST boundary-specific flux.
        data->boundaryVelocityCorrection = (data->targetVolumetricFlux - globalCurrentBoundaryFlux) / globalBoundaryArea;
    } else {
        simCtx->bulkVelocityCorrection = 0.0;
        data->boundaryVelocityCorrection = 0.0;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Driven Flow Controller Update (Dir %c):\n", data->direction);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Target Volumetric Flux:          %.6e\n", data->targetVolumetricFlux);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Avg Planar Volumetric Flux (Stable): %.6e\n", globalAveragePlanarVolumetricFlux);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Boundary Flux (Fast):            %.6e\n", globalCurrentBoundaryFlux);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Bulk Velocity Correction:     %.6e (For Momentum Source)\n", simCtx->bulkVelocityCorrection);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  - Boundary Velocity Correction: %.6e (For Boundary Trim)\n", data->boundaryVelocityCorrection);
    
    // Suppress unused parameter warnings for this handler
    (void)local_inflow_contribution;
    (void)local_outflow_contribution;
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Apply_PeriodicDrivenConstant"
/**
 * @brief (Handler Apply) Applies the immediate "Boundary Trim" velocity correction to the periodic face.
 *
 * This method executes the fast, tactical part of the control loop. It is called during the
 * `Apply` phase of the boundary condition cycle for each of the two periodic faces that have
 * this handler.
 *
 * It reads the `boundaryVelocityCorrection` value (which was computed for the entire boundary
 * by the master controller's `PreStep` method) from its private data. It then applies this
 * correction directly to the contravariant velocity (`ucont`) field on the face it manages.
 *
 * This action serves to immediately correct any flux deviation at the periodic "seam",
 * preventing errors from propagating into the domain in the subsequent timestep. The correction
 * value is also stored in the `uch` vector for diagnostic purposes.
 *
 * @param self The `BoundaryCondition` object containing the handler's state.
 * @param ctx  The `BCContext`, providing access to the `UserCtx` and `face_id`.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Apply_PeriodicDrivenConstant(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    DrivenConstantData *data = (DrivenConstantData*)self->data;
    UserCtx* user = ctx->user;
    BCFace face_id = ctx->face_id;
    PetscBool can_service;
    
    PetscFunctionBeginUser;
    
    // Check if this rank owns part of this boundary face
    ierr = CanRankServiceFace(&user->info, user->IM, user->JM, user->KM, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PetscFunctionReturn(0);
    }
    
    // If the correction is negligible, no work is needed.
    if (PetscAbsReal(data->boundaryVelocityCorrection) < 1e-12) {
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(LOCAL, LOG_TRACE, "Apply_PeriodicDrivenConstant: Applying boundary trim on Face %s...\n", BCFaceToString(face_id));

    // --- Get read/write access to necessary arrays ---
    DMDALocalInfo info = user->info;
    Cmpnts ***ucont, ***uch, ***csi, ***eta, ***zet;
    PetscReal ***nvert;
    
    ierr = DMDAVecGetArray(user->fda, user->lUcont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Bcs.Uch, &uch); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);
    
    PetscInt lxs = (info.xs == 0) ? 1 : info.xs;
    PetscInt lys = (info.ys == 0) ? 1 : info.ys;
    PetscInt lzs = (info.zs == 0) ? 1 : info.zs;
    PetscInt lxe = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    PetscInt lye = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    PetscInt lze = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    // --- Apply correction to the appropriate face and velocity component ---
    switch (face_id) {
        case BC_FACE_NEG_X: case BC_FACE_POS_X: {
            PetscInt i_face = (face_id == BC_FACE_NEG_X) ? info.xs : info.mx - 2;
            PetscInt i_nvert = (face_id == BC_FACE_NEG_X) ? info.xs + 1 : info.mx - 2;

            for (PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++) {
                if (nvert[k][j][i_nvert] < 0.1) {
                    PetscReal faceArea = sqrt(csi[k][j][i_face].x*csi[k][j][i_nvert].x + csi[k][j][i_nvert].y*csi[k][j][i_nvert].y + csi[k][j][i_face].z*csi[k][j][i_face].z);
                    PetscReal fluxTrim = data->boundaryVelocityCorrection * faceArea;
                    if(data->applyBoundaryTrim) ucont[k][j][i_face].x += fluxTrim;
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
                    PetscReal fluxTrim = data->boundaryVelocityCorrection * faceArea;
                    if(data->applyBoundaryTrim) ucont[k][j_face][i].y += fluxTrim;
                    uch[k][j_face][i].y = fluxTrim;
                }
            }
        } break;
            
        case BC_FACE_NEG_Z: case BC_FACE_POS_Z: {
            PetscInt k_face = (face_id == BC_FACE_NEG_Z) ? info.zs : info.mz - 2;
            PetscInt k_nvert = (face_id == BC_FACE_NEG_Z) ? info.zs + 1 : info.mz - 2;
            
            for (PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++) {
                if (nvert[k_nvert][j][i] < 0.1) {
                    PetscReal faceArea = sqrt(zet[k_nvert][j][i].x*zet[k_nvert][j][i].x + zet[k_nvert][j][i].y*zet[k_nvert][j][i].y + zet[k_nvert][j][i].z*zet[k_nvert][j][i].z);
                    PetscReal fluxTrim = data->boundaryVelocityCorrection * faceArea;
                    if(data->applyBoundaryTrim) ucont[k_face][j][i].z += fluxTrim;
                    uch[k_face][j][i].z = fluxTrim;
                }
            }
        } break;
    }

    // --- Restore arrays ---
    ierr = DMDAVecRestoreArray(user->fda, user->lUcont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Uch, &uch); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&nvert); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Destroy_PeriodicDrivenConstant"
/**
 * @brief (Handler Destructor) Frees the memory allocated for the handler's private data.
 *
 * This method is called once at the end of the simulation by `BoundarySystem_Destroy`.
 * Its only job is to free the `DrivenConstantData` struct that was allocated in the
 * `Create_PeriodicDrivenConstantFlux` constructor.
 *
 * This is a critical step to ensure there are no memory leaks.
 *
 * @param self The `BoundaryCondition` object containing the private data to be freed.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Destroy_PeriodicDrivenConstant(BoundaryCondition *self)
{
    PetscFunctionBeginUser;

    // Check that the handler object and its private data pointer are valid before trying to free.
    if (self && self->data) {
        // The private data was allocated with PetscNew(), so it must be freed with PetscFree().
        PetscFree(self->data);
        
        // It is good practice to nullify the pointer after freeing to prevent
        // any accidental use of the dangling pointer (use-after-free).
        self->data = NULL;
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "Destroy_PeriodicDrivenConstant: Private data freed successfully.\n");
    }

    PetscFunctionReturn(0);
}
