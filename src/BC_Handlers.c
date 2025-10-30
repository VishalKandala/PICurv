#include "BC_Handlers.h"     // The header that declares this file's "constructor" functions

//================================================================================
//
//               HANDLER IMPLEMENTATION: NO-SLIP WALL
//               (Corresponds to BC_HANDLER_WALL_NOSLIP)
//
// This handler implements a stationary, impenetrable wall where the fluid
// velocity is zero (no-slip condition). It is one of the simplest but most
// common boundary conditions.
//
//================================================================================

/**
 * @brief (Handler Action) Applies the no-slip wall condition to a specified face.
 *
 * This function is the core implementation for the no-slip wall. It is called by the
 * BoundarySystem during the Apply phase of each time step. Its responsibilities are:
 *   1. Check if the current MPI rank owns any part of the face to be processed.
 *   2. If so, iterate over the local portion of that face.
 *   3. For each boundary cell face, set the normal contravariant velocity (flux) to zero.
 *   4. Set the Cartesian ghost cell velocity to enforce the no-slip condition
 *      (e.g., u_ghost = -u_interior), which is a common second-order accurate
 *      approximation for a wall located exactly on the cell face.
 *
 * @param self A pointer to the BoundaryCondition object. This simple handler does not
 *             use this parameter, but it is required by the standard interface.
 * @param ctx  A pointer to the BCContext, which provides access to the UserCtx (containing
 *             the vectors to be modified) and the ID of the face to be processed.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Apply_WallNoSlip(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    PetscBool      can_service;
    PetscMPIInt    rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    // This unused-variable pragma silences compiler warnings for simple handlers
    // that don't need their own 'self' pointer.
    (void)self;
    
    PetscFunctionBeginUser;

    // Step 1: Check if this MPI rank has any work to do for this face.
    // This is a critical efficiency step to avoid unnecessary work and memory access.
    // It uses a utility function that understands the parallel grid decomposition.
    ierr = CanRankServiceFace(&user->info, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PetscFunctionReturn(0); // This rank does not own this boundary face. Exit silently.
    }

    // If we proceed, this rank has work to do. Log the action.
    const char *face_name = BCFaceToString(face_id);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply_WallNoSlip: Rank %d applying condition to Face %d (%s).\n",rank, face_id, face_name);

    // Step 2: Get safe access to the local PETSc vector arrays.
    // We get the local vectors (Ucat, Ucont) because ghost cell data is only
    // guaranteed to be correct on the local representation after a ghost update.
    DMDALocalInfo  *info = &user->info;
    Cmpnts       ***ucat, ***ucont;
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);

    // Step 3: Apply the no-slip condition based on which face is being processed.
    // The loop bounds are for the *owned* nodes on this rank (xs to xe, etc.).
    // The logic inside correctly accesses the ghost nodes or face-normal fluxes.
    PetscInt i, j, k;
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;
    PetscInt mx = info->mx, my = info->my, mz = info->mz; // Global dimensions

    switch (face_id) {
        case BC_FACE_NEG_X: // -Xi face at global index i=0
            i = 0;
            for (k = zs; k < ze; k++) {
                for (j = ys; j < ye; j++) {
                    ucont[k][j][i].x = 0.0; // Set normal contravariant flux to zero.
                    // Set Cartesian ghost cell velocity for no-slip: u_ghost(i) = -u_interior(i+1)
                    ucat[k][j][i].x = -ucat[k][j][i+1].x;
                    ucat[k][j][i].y = -ucat[k][j][i+1].y;
                    ucat[k][j][i].z = -ucat[k][j][i+1].z;
                }
            }
            break;

        case BC_FACE_POS_X: // +Xi face at global index i = mx-1
            i = mx - 1;
            for (k = zs; k < ze; k++) {
                for (j = ys; j < ye; j++) {
                    ucont[k][j][i-1].x = 0.0;
                    // u_ghost(i) = -u_interior(i-1)
                    ucat[k][j][i].x = -ucat[k][j][i-1].x;
                    ucat[k][j][i].y = -ucat[k][j][i-1].y;
                    ucat[k][j][i].z = -ucat[k][j][i-1].z;
                }
            }
            break;

        case BC_FACE_NEG_Y: // -Eta face at global index j=0
            j = 0;
            for (k = zs; k < ze; k++) {
                for (i = xs; i < xe; i++) {
                    ucont[k][j][i].y = 0.0;
                    ucat[k][j][i].x = -ucat[k][j+1][i].x;
                    ucat[k][j][i].y = -ucat[k][j+1][i].y;
                    ucat[k][j][i].z = -ucat[k][j+1][i].z;
                }
            }
            break;

        case BC_FACE_POS_Y: // +Eta face at global index j = my-1
            j = my - 1;
            for (k = zs; k < ze; k++) {
                for (i = xs; i < xe; i++) {
                    ucont[k][j-1][i].y = 0.0;
                    ucat[k][j][i].x = -ucat[k][j-1][i].x;
                    ucat[k][j][i].y = -ucat[k][j-1][i].y;
                    ucat[k][j][i].z = -ucat[k][j-1][i].z;
                }
            }
            break;
            
        case BC_FACE_NEG_Z: // -Zeta face at global index k=0
            k = 0;
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    ucont[k][j][i].z = 0.0;
                    ucat[k][j][i].x = -ucat[k+1][j][i].x;
                    ucat[k][j][i].y = -ucat[k+1][j][i].y;
                    ucat[k][j][i].z = -ucat[k+1][j][i].z;
                }
            }
            break;

        case BC_FACE_POS_Z: // +Zeta face at global index k = mz-1
            k = mz - 1;
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    ucont[k-1][j][i].z = 0.0;
                    ucat[k][j][i].x = -ucat[k-1][j][i].x;
                    ucat[k][j][i].y = -ucat[k-1][j][i].y;
                    ucat[k][j][i].z = -ucat[k-1][j][i].z;
                }
            }
            break;
    }

    // Step 4: Restore safe access to the PETSc vector arrays.
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with No-Slip Wall behavior.
 *
 * This function is called by the BoundarySystem factory (`BoundaryCondition_Create` in
 * `Boundaries.c`). It "constructs" a no-slip wall handler by setting the function
 * pointers in the provided BoundaryCondition struct to point to the static functions
 * defined in this file.
 *
 * A no-slip wall is simple and requires only the `Apply` method:
 *  - It needs no special initialization (`Initialize` is NULL).
 *  - It does not contribute to the global mass balance (`PreStep` is NULL).
 *  - It has a specific `Apply` function to enforce zero velocity.
 *  - It allocates no private data, so `Destroy` is NULL.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_WallNoSlip(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition object is NULL in Create_WallNoSlip");

    // Assign the appropriate function pointers for this handler type.
    // This is the essence of the "Strategy" pattern in C.
    bc->Initialize = NULL;
    bc->PreStep    = NULL;
    bc->Apply      = Apply_WallNoSlip; // The ONLY action this handler needs to perform.
    bc->Destroy    = NULL;
    
    // No private data struct is needed for this handler, so bc->data remains NULL.

    PetscFunctionReturn(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

//================================================================================
//
//          HANDLER IMPLEMENTATION: CONSTANT VELOCITY INLET
//          (Corresponds to BC_HANDLER_INLET_CONSTANT_VELOCITY)
//
// This handler implements an inlet with a prescribed, uniform Cartesian velocity.
//
//================================================================================

// --- 1. FORWARD DECLARATIONS for this handler's static methods ---
// This tells the compiler that these functions exist before they are used.
// It resolves the "declared implicitly" warning.
static PetscErrorCode Initialize_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx, PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode Destroy_InletConstantVelocity(BoundaryCondition *self);

/**
 * @brief Private data structure for the Constant Velocity Inlet handler.
 *
 * This struct holds the specific parameters needed for this handler, which are
 * parsed from the bcs.dat file during initialization.
 */
typedef struct {
    Cmpnts specified_velocity; // The desired Cartesian velocity (vx, vy, vz)
} InletConstantData;


/**
 * @brief (Handler Action) Sets the initial state of the boundary face.
 *
 * This function is called once by BoundarySystem_Create. It performs two tasks:
 * 1. Parses the `params` list (from bcs.dat) to find the 'vx', 'vy', and 'vz'
 *    parameters and stores them in its private `data` struct.
 * 2. Sets the initial Ucont and Ucat values on the boundary face to reflect
 *    this constant velocity.
 *
 * @param self The BoundaryCondition object for this handler.
 * @param ctx  The context, providing access to UserCtx and face_id.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Initialize_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletConstantData *data = (InletConstantData*)self->data;
    
    PetscFunctionBeginUser;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initialize_InletConstantVelocity: Initializing handler for Face %d. \n", face_id);

    // 1. Parse parameters from the linked list stored in the face configuration.
    data->specified_velocity = (Cmpnts){0.0, 0.0, 0.0}; // Default to zero
    for (BC_Param *p = user->boundary_faces[face_id].params; p; p = p->next) {
        if (strcasecmp(p->key, "vx") == 0) data->specified_velocity.x = atof(p->value);
        else if (strcasecmp(p->key, "vy") == 0) data->specified_velocity.y = atof(p->value);
        else if (strcasecmp(p->key, "vz") == 0) data->specified_velocity.z = atof(p->value);
    }
    LOG_ALLOW(LOCAL, LOG_INFO, "  Inlet Face %d configured with velocity (vx,vy,vz) = (%.2f, %.2f, %.2f) \n",
              face_id, data->specified_velocity.x, data->specified_velocity.y, data->specified_velocity.z);

    // 2. Set the initial boundary state. We can simply call the Apply function to do this,
    //    as the logic is the same for the initial state and subsequent steps for this handler.
    ierr = Apply_InletConstantVelocity(self, ctx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Action) Calculates the target inflow flux for the PreStep phase.
 *
 * This function calculates the total volumetric flux that *should* be entering
 * through the portion of the inlet face owned by this MPI rank. It does this by
 * dotting the specified Cartesian velocity with the face-normal area vectors.
 *
 * @param self The BoundaryCondition object for this handler.
 * @param ctx  The context, providing UserCtx and face_id.
 * @param[out] local_inflow_contribution  The calculated inflow flux is added to this value.
 * @param[out] local_outflow_contribution This handler does not produce outflow, so this is unused.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode PreStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution)
{
    // ... (This function would calculate the flux: V_spec · Area_vector) ...
    // For Phase 2, we can leave this as a placeholder, as it's only needed
    // when we have a mass-conserving outlet that needs to know the target inflow.
    (void)self; (void)ctx; (void)local_inflow_contribution; (void)local_outflow_contribution;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

/**
 * @brief (Handler Action) Applies the constant velocity inlet condition to a specified face.
 *
 * This function is the core "workhorse" for the constant velocity inlet. It is called
 * by the BoundarySystem during the Apply phase of each time step. Its responsibility
 * is to enforce the pre-configured velocity on its assigned face. It does this by:
 *
 * 1.  Setting the ghost-cell Cartesian velocity (`Ucat`) directly to the specified
 *     velocity vector (e.g., `ucat[k][j][0] = {vx, vy, vz}`). This provides a
 *     Dirichlet condition for any calculations using Cartesian velocity.
 *
 * 2.  Calculating the corresponding contravariant velocity components (`Ucont`). This
 *     is done by projecting the specified Cartesian velocity vector onto the
 *     contravariant basis vectors at each point on the face. For example, the first
 *     contravariant component is U¹ = v_cartesian ⋅ g¹ = (vx,vy,vz) ⋅ (csi.x,csi.y,csi.z).
 *     This ensures that the volumetric flux through each face is consistent with the
 *     specified Cartesian velocity.
 *
 * @param self The BoundaryCondition object for this handler, which contains the
 *             private `data` struct holding the configured velocity.
 * @param ctx  The context, providing access to UserCtx (for PETSc vectors and grid info)
 *             and the ID of the face to be processed.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Apply_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletConstantData *data = (InletConstantData*)self->data; // Cast private data
    PetscBool      can_service;
    PetscMPIInt    rank;
    
    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Step 1: Check if this MPI rank owns any part of the face to be processed.
    ierr = CanRankServiceFace(&user->info, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PetscFunctionReturn(0); // This rank has no work to do for this face.
    }

    const char* face_name = BCFaceToString(face_id);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply_InletConstantVelocity: Rank %d applying condition to Face %d (%s).\n", rank, face_id, face_name);

    // Step 2: Get access to the necessary PETSc vector arrays.
    DMDALocalInfo  *info = &user->info;
    Cmpnts       ***ucat, ***ucont;
    Cmpnts       ***l_csi, ***l_eta, ***l_zet;
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &l_csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &l_eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &l_zet); CHKERRQ(ierr);

    const Cmpnts *v_spec = &data->specified_velocity; // Create a convenient shortcut

    // Step 3: Loop over the owned nodes on the specified face and apply the BC.
    PetscInt i, j, k;
    const PetscInt xs = info->xs, xe = info->xs + info->xm;
    const PetscInt ys = info->ys, ye = info->ys + info->ym;
    const PetscInt zs = info->zs, ze = info->zs + info->zm;
    const PetscInt mx = info->mx, my = info->my, mz = info->mz;

    // The loops iterate over the entire local ownership range. The `if` conditions
    // inside ensure that the logic is only applied to the nodes on the correct boundary face.
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {

                if (i == 0 && face_id == BC_FACE_NEG_X) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j][i].x = v_spec->x * l_csi[k][j][i].x + v_spec->y * l_csi[k][j][i].y + v_spec->z * l_csi[k][j][i].z;
                }
                else if (i == mx - 1 && face_id == BC_FACE_POS_X) {
                    ucat[k][j][i] = *v_spec;
                    // Note: Ucont.x is on the i-face, so for the max face, it's at index i-1
                    ucont[k][j][i-1].x = v_spec->x * l_csi[k][j][i-1].x + v_spec->y * l_csi[k][j][i-1].y + v_spec->z * l_csi[k][j][i-1].z;
                }
                else if (j == 0 && face_id == BC_FACE_NEG_Y) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j][i].y = v_spec->x * l_eta[k][j][i].x + v_spec->y * l_eta[k][j][i].y + v_spec->z * l_eta[k][j][i].z;
                }
                else if (j == my - 1 && face_id == BC_FACE_POS_Y) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j-1][i].y = v_spec->x * l_eta[k][j-1][i].x + v_spec->y * l_eta[k][j-1][i].y + v_spec->z * l_eta[k][j-1][i].z;
                }
                else if (k == 0 && face_id == BC_FACE_NEG_Z) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j][i].z = v_spec->x * l_zet[k][j][i].x + v_spec->y * l_zet[k][j][i].y + v_spec->z * l_zet[k][j][i].z;
                }
                else if (k == mz - 1 && face_id == BC_FACE_POS_Z) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k-1][j][i].z = v_spec->x * l_zet[k-1][j][i].x + v_spec->y * l_zet[k-1][j][i].y + v_spec->z * l_zet[k-1][j][i].z;
                }
            }
        }
    }
    
    // Step 4: Restore safe access to the PETSc vector arrays.
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &l_csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &l_eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &l_zet); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

/**
 * @brief (Handler Destructor) Frees memory allocated by the Constant Velocity Inlet handler.
 */
static PetscErrorCode Destroy_InletConstantVelocity(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    if (self && self->data) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Destroy_InletConstantVelocity: Freeing private data struct.");
        PetscFree(self->data);
        self->data = NULL;
    }
    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with Constant Velocity Inlet behavior.
 */
PetscErrorCode Create_InletConstantVelocity(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition object is NULL");

    // Allocate the private data structure for this handler
    InletConstantData *data = NULL;
    ierr = PetscMalloc1(1, &data); CHKERRQ(ierr);
    bc->data = (void*)data;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Create_InletConstantVelocity: Allocated InletConstantData (%zu bytes) at %p.\n",
              sizeof(*data), (void*)data);
    
    // Assign the function pointers for this handler type
    bc->Initialize = Initialize_InletConstantVelocity;
    bc->PreStep    = PreStep_InletConstantVelocity;
    bc->Apply      = Apply_InletConstantVelocity;
    bc->Destroy    = Destroy_InletConstantVelocity;
    
    PetscFunctionReturn(0);
}

/*****************************************************************************
 *  NOGRAD – copy-ghost handler
 *
 *  Behaviour
 *  ---------
 *    • Works only on the face(s) where it is prescribed in bcs.dat.
 *    • Copies the entire *first interior* node/face layer onto the ghost
 *      layer for both Cartesian velocity  Ucat  and contravariant velocity
 *      Ucont (normal component only – i.e. the one that lives on that face).
 *    • No Initialise() or PreStep() are required.
 *
 *  Integration
 *  -----------
 *    Called by BoundarySystem_ExecuteStep() after contra2cart().
 *****************************************************************************/

/* --- no private data is needed, keep an empty struct for future growth --- */
typedef struct { int dummy; } NgData;

/* ------------------------------------------------------------------------- */
static PetscErrorCode Apply_NogradCopyGhost(BoundaryCondition *self,
                                            BCContext         *ctx)
{
    PetscErrorCode ierr;
    UserCtx       *u   = ctx->user;
    DMDALocalInfo *inf = &u->info;

    /* Arrays */
    Cmpnts ***ucat, ***ucont;
    ierr = DMDAVecGetArray(u->fda, u->Ucat , &ucat ); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(u->fda, u->Ucont, &ucont); CHKERRQ(ierr);

    const PetscInt xs = inf->xs, xe = inf->xs + inf->xm;
    const PetscInt ys = inf->ys, ye = inf->ys + inf->ym;
    const PetscInt zs = inf->zs, ze = inf->zs + inf->zm;
    const PetscInt mx = inf->mx, my = inf->my, mz = inf->mz;

    PetscInt i,j,k;

    /*
    PetscInt lxs  = (xs == 0) ? xs + 1 : xs;
    PetscInt lxe   = (xe == mx) ? xe - 1 : xe;

    PetscInt lys  = (ys == 0) ? ys + 1 : ys;
    PetscInt lye   = (ye == my) ? ye - 1 : ye;

    PetscInt lzs  = (zs == 0) ? zs + 1 : zs;
    PetscInt lze   = (ze == mz) ? ze - 1 : ze;
    */

    PetscInt lxs  = xs;
    PetscInt lxe   = xe;

    PetscInt lys  = ys;
    PetscInt lye   = ye;

    PetscInt lzs  = zs;
    PetscInt lze   = ze;
    
    switch (ctx->face_id)
    {
    /* ------------------------------------------------------------------ */
    case BC_FACE_NEG_X:   /* i = 0 copies from i = 1 */
        if (xs == 0)
        {
            i = 0;
            for (k=lzs;k<lze;k++)
              for (j=lys;j<lye;j++)
              {
                  ucat [k][j][i]   = ucat [k][j][i+1];
                  ucont[k][j][i].x = ucont[k][j][i+1].x;    /* flux normal to face */
              }
        }
        break;

    case BC_FACE_POS_X:   /* i = mx-1 copies from i = mx-2 */
        if (xe == mx)
        {
            i = mx-1;
            for (k=lzs;k<lze;k++)
              for (j=lys;j<lye;j++)
              {
                  ucat [k][j][i]     = ucat [k][j][i-1];
		  // Add a guard to prevent out-of-bounds read on thin domains
		  if (mx >= 3) {
		    ucont[k][j][i].x = ucont[k][j][i-1].x;
		  }
              }
        }
        break;

    /* ------------------------------------------------------------------ */
    case BC_FACE_NEG_Y:   /* j = 0 from j = 1 */
        if (ys == 0)
        {
            j = 0;
            for (k=lzs;k<lze;k++)
              for (i=lxs;i<lxe;i++)
              {
                  ucat [k][j][i]   = ucat [k][j+1][i];
                  ucont[k][j][i].y = ucont[k][j+1][i].y;
              }
        }
        break;

    case BC_FACE_POS_Y:   /* j = my-1 from j = my-2 */
        if (ye == my)
        {
            j = my-1;
            for (k=zs;k<ze;k++)
              for (i=xs;i<xe;i++)
              {
                  ucat [k][j][i]     = ucat [k][j-1][i];
		  // Add a guard to prevent out-of-bounds read on thin domains
		  if (my >= 3) {
		    ucont[k][j][i].y = ucont[k][j-1][i].y;
		  }
	      }
        }
        break;

    /* ------------------------------------------------------------------ */
    case BC_FACE_NEG_Z:   /* k = 0 from k = 1 */
        if (zs == 0)
        {
            k = 0;
            for (j=lys;j<lye;j++)
              for (i=lxs;i<lxe;i++)
              {
                  ucat [k][j][i]   = ucat [k+1][j][i];
                  ucont[k][j][i].z = ucont[k+1][j][i].z;
              }
        }
        break;

    case BC_FACE_POS_Z:   /* k = mz-1 from k = mz-2 */
        if (ze == mz)
        {
            k = mz-1;
            for (j=lys;j<lye;j++)
              for (i=lxs;i<lxe;i++)
              {
                  ucat [k][j][i]     = ucat [k-1][j][i];
		  // Add a guard to prevent out-of-bounds read on thin domains
		  if (mz >= 3) {
		    ucont[k][j][i].z = ucont[k-1][j][i].z;
		  }
              }
        }
        break;
    }

    ierr = DMDAVecRestoreArray(u->fda, u->Ucat , &ucat ); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(u->fda, u->Ucont, &ucont); CHKERRQ(ierr);
    return 0;
}

/* ------------------------------------------------------------------------- */
static PetscErrorCode Destroy_NogradCopyGhost(BoundaryCondition *self)
{
    PetscFree(self->data);      /* nothing in it yet, but keep it symmetric  */
    self->data = NULL;
    return 0;
}

/* ------------------------------------------------------------------------- */
PetscErrorCode Create_NogradCopyGhost(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    ierr = PetscMalloc1(1,&bc->data); CHKERRQ(ierr);   /* allocate empty struct */

    bc->Initialize = NULL;       /* no parameters to parse */
    bc->PreStep    = NULL;       /* no flux bookkeeping    */
    bc->Apply      = Apply_NogradCopyGhost;
    bc->Destroy    = Destroy_NogradCopyGhost;
    return 0;
}

//================================================================================
//
//          HANDLER IMPLEMENTATION: PARABOLIC PROFILE INLET
//          (Corresponds to a new handler type, e.g., BC_HANDLER_INLET_PARABOLIC)
//
// This handler implements an inlet with a parabolic velocity profile.
// The flow is assumed to be aligned with the local grid direction of the face.
//
//================================================================================

// --- 1. FORWARD DECLARATIONS ---
static PetscErrorCode Initialize_InletParabolicProfile(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_InletParabolicProfile(BoundaryCondition *self, BCContext *ctx, PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_InletParabolicProfile(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode Destroy_InletParabolicProfile(BoundaryCondition *self);

/**
 * @brief Private data for the Parabolic Profile Inlet handler.
 */
typedef struct {
    PetscReal max_velocity; // The desired maximum centerline Cartesian velocity
} InletParabolicData;


/**
 * @brief (Handler Action) Initializes the parabolic inlet handler.
 *
 * Parses the 'vmax' parameter from bcs.dat and sets the initial state
 * on the boundary face by calling the Apply function.
 */
static PetscErrorCode Initialize_InletParabolicProfile(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletParabolicData *data = (InletParabolicData*)self->data;
    
    PetscFunctionBeginUser;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initialize_InletParabolicProfile: Initializing handler for Face %d. \n", face_id);

    // 1. Parse parameters to find the maximum velocity.
    data->max_velocity = 1.0; // Default to 1.0
    for (BC_Param *p = user->boundary_faces[face_id].params; p; p = p->next) {
        if (strcasecmp(p->key, "vmax") == 0) data->max_velocity = atof(p->value);
    }
    LOG_ALLOW(LOCAL, LOG_INFO, "  Inlet Face %d configured with parabolic profile, vmax = %.2f \n",
              face_id, data->max_velocity);

    // 2. Set the initial boundary state by calling the Apply function.
    ierr = Apply_InletParabolicProfile(self, ctx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Action) Calculates target inflow flux (Placeholder).
 */
static PetscErrorCode PreStep_InletParabolicProfile(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution)
{
    (void)self; (void)ctx; (void)local_inflow_contribution; (void)local_outflow_contribution;
    PetscFunctionBeginUser;
    // This could be implemented later by integrating the profile over the face area.
    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Action) Applies the parabolic velocity inlet condition.
 *
 * This function enforces the parabolic profile on its assigned face. It does this by:
 * 1.  Determining the normalized coordinates of each point on the face using grid indices.
 * 2.  Calculating the desired physical velocity magnitude using the parabolic formula.
 * 3.  Determining the local direction of flow from the face-normal area vector.
 * 4.  Constructing the full Cartesian velocity vector (magnitude * direction).
 * 5.  Setting the ghost-cell Cartesian velocity (Ucat).
 * 6.  Projecting the Cartesian velocity to find the contravariant components (Ucont).
 */
static PetscErrorCode Apply_InletParabolicProfile(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletParabolicData *data = (InletParabolicData*)self->data;
    PetscBool      can_service;
    
    PetscFunctionBeginUser;

    ierr = CanRankServiceFace(&user->info, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PetscFunctionReturn(0);
    }

    // Get necessary arrays
    DMDALocalInfo  *info = &user->info;
    Cmpnts       ***ucat, ***ucont;
    Cmpnts       ***l_csi, ***l_eta, ***l_zet;
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &l_csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &l_eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &l_zet); CHKERRQ(ierr);

    const PetscInt mx = info->mx, my = info->my, mz = info->mz;

    // Define channel geometry in "index space"
    // For an inlet on a Z-face, the cross-section is in the (i,j) plane.
    const PetscReal i_width  = (PetscReal)(mx - 2);
    const PetscReal j_width  = (PetscReal)(my - 2);
    const PetscReal i_center = 1.0 + i_width / 2.0;
    const PetscReal j_center = 1.0 + j_width / 2.0;

    // The main loop over all owned nodes. `if` conditions select the correct face.
    for (PetscInt k = info->zs; k < info->zs + info->zm; k++) {
        for (PetscInt j = info->ys; j < info->ys + info->ym; j++) {
            for (PetscInt i = info->xs; i < info->xs + info->xm; i++) {
                
                // This logic only supports Z-face inlets for simplicity, as per the bent channel problem.
                // It can be extended to X and Y faces if needed.
                if ((k == 0 && face_id == BC_FACE_NEG_Z) || (k == mz - 1 && face_id == BC_FACE_POS_Z)) {
                    
                    // --- This is the core logic, identical to the interior initialization ---
                    
                    // 1. Calculate normalized coordinates in the cross-section
                    const PetscReal i_norm = (i - i_center) / (i_width / 2.0);
                    const PetscReal j_norm = (j - j_center) / (j_width / 2.0);

                    // 2. Calculate the desired velocity magnitude using the profile
                    const PetscReal profile_i = 1.0 - i_norm * i_norm;
                    const PetscReal profile_j = 1.0 - j_norm * j_norm;
                    PetscReal v_mag = data->max_velocity * profile_i * profile_j;
                    if (v_mag < 0.0) v_mag = 0.0;
                    
                    // 3. Determine the local flow direction (unit vector) from the grid
                    const Cmpnts *area_vec = &l_zet[k][j][i];
                    PetscReal area_mag = sqrt(area_vec->x * area_vec->x + area_vec->y * area_vec->y + area_vec->z * area_vec->z);
                    Cmpnts v_dir = {0.0, 0.0, 0.0};
                    if (area_mag > 1e-9) {
                        v_dir.x = area_vec->x / area_mag;
                        v_dir.y = area_vec->y / area_mag;
                        v_dir.z = area_vec->z / area_mag;
                    }
                    
                    // On a negative face, the area vector points OUT of the domain. Flow is IN, so we flip the sign.
                    if (face_id == BC_FACE_NEG_Z) {
		      //  v_dir.x *= -1.0; v_dir.y *= -1.0; v_dir.z *= -1.0;
                    }
                    
                    // 4. Construct the final Cartesian velocity vector
                    Cmpnts v_cartesian;
                    v_cartesian.x = v_mag * v_dir.x;
                    v_cartesian.y = v_mag * v_dir.y;
                    v_cartesian.z = v_mag * v_dir.z;
                    
                    // 5. Apply the conditions
                    PetscInt ucont_k = (face_id == BC_FACE_POS_Z) ? k-1 : k; // Staggered grid correction

                    ucat[k][j][i]  = v_cartesian;
                    ucont[ucont_k][j][i].z = v_cartesian.x * l_zet[ucont_k][j][i].x +
                                             v_cartesian.y * l_zet[ucont_k][j][i].y +
                                             v_cartesian.z * l_zet[ucont_k][j][i].z;
                }
            }
        }
    }
    
    // Restore arrays
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &l_csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &l_eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &l_zet); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

/**
 * @brief (Handler Destructor) Frees memory for the Parabolic Profile Inlet.
 */
static PetscErrorCode Destroy_InletParabolicProfile(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    if (self && self->data) {
        PetscFree(self->data);
        self->data = NULL;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with Parabolic Inlet behavior.
 */
PetscErrorCode Create_InletParabolicProfile(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition object is NULL");

    InletParabolicData *data = NULL;
    ierr = PetscMalloc1(1, &data); CHKERRQ(ierr);
    bc->data = (void*)data;
    
    // Assign the function pointers for this handler type
    bc->Initialize = Initialize_InletParabolicProfile;
    bc->PreStep    = PreStep_InletParabolicProfile;
    bc->Apply      = Apply_InletParabolicProfile;
    bc->Destroy    = Destroy_InletParabolicProfile;
    
    PetscFunctionReturn(0);
}

//================================================================================
//
//          HANDLER IMPLEMENTATION: OUTLET WITH CONSERVATION OF MOMENTUM
//          (Corresponds to a new handler type, e.g., BC_HANDLER_OUTLET_CONSERVATION)
//
// This handler implements an outlet and ensures conservation of momentum by adjusting
// the contravariant velocity component at the outlet face based on inlet conditions and using a ratio.
//
//================================================================================

// --- 1. FORWARD DECLARATIONS ---
static PetscErrorCode PreStep_OutletConservation(BoundaryCondition *self, BCContext *ctx,
                                     PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution);
/**
 * @brief (Handler Action) Measures the current, uncorrected flux passing through a SINGLE outlet face.
 *
 * This function is called during the PreStep phase of the boundary condition cycle. It is a direct,
 * high-fidelity port of the flux measurement logic from the legacy `FormBCS` function for
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
    const PetscInt IM_nodes_global = user->IM + 1;
    const PetscInt JM_nodes_global = user->JM + 1;
    const PetscInt KM_nodes_global = user->KM + 1;
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

