
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