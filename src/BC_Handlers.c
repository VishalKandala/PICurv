#include "BC_Handlers.h"     // The header that declares this file's "constructor" functions


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
            if (xs= 0){  
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
            if (ys == 0) break;
            PetscInt j = ys;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k][j][i].y = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                }
            }
        } break;

        case BC_FACE_POS_Y: {
            if (ye != my) break;
            PetscInt j = ye - 1;
            for (PetscInt k = lzs; k < lze; k++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k][j-1][i].y = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                }
            }
        } break;
            
        case BC_FACE_NEG_Z: {
            if (zs != 0) break;
            PetscInt k = zs;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k][j][i].z = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
                }
            }
        } break;

        case BC_FACE_POS_Z: {
            if (ze != mz) break;
            PetscInt k = ze - 1;
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    ucont[k-1][j][i].z = 0.0;
                    
                    ubcs[k][j][i].x = 0.0;
                    ubcs[k][j][i].y = 0.0;
                    ubcs[k][j][i].z = 0.0;
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