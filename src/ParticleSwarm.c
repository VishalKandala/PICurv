    // ParticleSwarm.c

#include "ParticleSwarm.h"

#define INTERPOLATION_DISTANCE_TOLERANCE 1.0e-14

#undef __FUNCT__
#define __FUNCT__ "InitializeSwarm"
/**
 * @brief Initializes the DMSwarm object within the UserCtx structure.
 *
 * This function creates the DMSwarm, sets its type and dimension, and configures basic swarm properties.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSwarm(UserCtx* user) {
    PetscErrorCode ierr;  // Error code for PETSc functions

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // Create the DMSwarm object for particle management
    ierr = DMCreate(PETSC_COMM_WORLD, &user->swarm); CHKERRQ(ierr);
    ierr = DMSetType(user->swarm, DMSWARM); CHKERRQ(ierr);
    ierr = DMSetDimension(user->swarm, 3); CHKERRQ(ierr);
    ierr = DMSwarmSetType(user->swarm, DMSWARM_BASIC); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "DMSwarm created and configured.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RegisterSwarmField"

/**
 * @brief Registers a swarm field without finalizing registration.
 *
 * This function calls DMSwarmRegisterPetscDatatypeField for the given field,
 * but does not finalize the registration. The finalization is deferred until
 * all fields have been registered.
 *
 * @param swarm      [in]  The DMSwarm object.
 * @param fieldName  [in]  Name of the field to register.
 * @param fieldDim   [in]  Dimension of the field (1 for scalar, 3 for vector, etc.).
 * @param dtype      [in]  The datatype of the swarm field being registered.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RegisterSwarmField(DM swarm, const char *fieldName, PetscInt fieldDim, PetscDataType dtype)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    ierr = DMSwarmRegisterPetscDatatypeField(swarm, fieldName, fieldDim, dtype); CHKERRQ(ierr);
    // PetscDataTypes is an extern char* [] defined in petscsystypes.h that gives string names for PetscDataType enums
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field '%s' with dimension=%d, type=%s.\n",
            fieldName, fieldDim, PetscDataTypes[dtype]);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RegisterParticleFields"

/**
 * @brief Registers necessary particle fields within the DMSwarm.
 *
 * This function registers fields such as position, velocity, CellID, and weight for each particle.
 *
 * @param[in,out] swarm   The DMSwarm object managing the particle swarm.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */

PetscErrorCode RegisterParticleFields(DM swarm)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    // Register each field using the helper function
    ierr = RegisterSwarmField(swarm, "position", 3 ,PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Registered field 'position'.\n");
    
    ierr = RegisterSwarmField(swarm, "velocity", 3, PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'velocity'.\n");
    
    ierr = RegisterSwarmField(swarm, "DMSwarm_CellID", 3, PETSC_INT); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'DMSwarm_CellID'.\n");
    
    ierr = RegisterSwarmField(swarm, "weight", 3,PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'weight'.\n");

    ierr = RegisterSwarmField(swarm,"Psi", 1,PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'Psi' - Scalar.\n");

    ierr =  RegisterSwarmField(swarm,"DMSwarm_location_status",1,PETSC_INT);CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'DMSwarm_location_status' - Status of Location of Particle(located,lost etc).\n");
    
    // Finalize the field registration after all fields have been added
    ierr = DMSwarmFinalizeFieldRegister(swarm); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO,"RegisterParticleFields - Finalized field registration.\n");
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DetermineVolumetricInitializationParameters"
/**
 * @brief Determines cell selection and intra-cell logical coordinates for volumetric initialization (Mode 1).
 *
 * This function is called when `simCtx->ParticleInitialization == 1`. It randomly selects
 * an *owned cell* on the current MPI rank and then generates random intra-cell logical
 * coordinates `[0,1)^3` within that chosen cell.
 *
 * The process involves:
 * 1. Checking if the current MPI rank owns any 3D cells.
 * 2. If it does, it randomly selects an *owned cell index* in each logical direction (i, j, k)
 *    by scaling a `[0,1)` random number with the number of owned cells in that direction.
 * 3. These local owned cell indices are then converted to the *local node indices*
 *    (`ci/cj/ck_metric_lnode_out`) corresponding to the origin of the selected cell,
 *    for use with `MetricLogicalToPhysical`. This conversion uses `xs/ys/zs_gnode`.
 * 4. All three intra-cell logical coordinates (`xi/eta/zta_metric_logic_out`) for
 *    `MetricLogicalToPhysical` are chosen randomly within the `[0,1)` range.
 * 5. A flag (`can_place_in_volume_out`) indicates if a valid placement could be determined.
 *
 * **Important Note on `DMDALocalInfo info` members (same as for surface init):**
 *   - `info->xs, info->ys, info->ks`: Global starting indices of *owned cells*.
 *   - `info->mx, info->my, info->mz`: Number of *grid points (nodes)* in each local dimension on this process.
 *     Therefore, the number of *owned cells* in a dimension is `info->mX - 1` (if `info->mX > 0`).
 *
 * @param[in]  user Pointer to `UserCtx`. (Currently not used in this specific helper, but kept for API consistency).
 * @param[in]  info Pointer to `DMDALocalInfo` for the current rank's grid portion.
 * @param[in]  xs_gnode, ys_gnode, zs_gnode Local indices (in the ghosted array) of the first *owned node*.
 * @param[in]  rand_logic_i_ptr Pointer to the RNG for i-dimension tasks [0,1).
 * @param[in]  rand_logic_j_ptr Pointer to the RNG for j-dimension tasks [0,1).
 * @param[in]  rand_logic_k_ptr Pointer to the RNG for k-dimension tasks [0,1).
 * @param[out] ci_metric_lnode_out Pointer to store the local i-node index of the selected cell's origin.
 * @param[out] cj_metric_lnode_out Pointer to store the local j-node index of the selected cell's origin.
 * @param[out] ck_metric_lnode_out Pointer to store the local k-node index of the selected cell's origin.
 * @param[out] xi_metric_logic_out Pointer to store the intra-cell logical xi-coordinate [0,1).
 * @param[out] eta_metric_logic_out Pointer to store the intra-cell logical eta-coordinate [0,1).
 * @param[out] zta_metric_logic_out Pointer to store the intra-cell logical zeta-coordinate [0,1).
 * @param[out] can_place_in_volume_out PETSC_TRUE if placement parameters were successfully determined, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, or a PETSc error code.
 */
static PetscErrorCode DetermineVolumetricInitializationParameters(
    UserCtx *user, DMDALocalInfo *info,
    PetscInt xs_gnode, PetscInt ys_gnode, PetscInt zs_gnode,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr, /* Pointers to RNGs */
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out,
    PetscBool *can_place_in_volume_out)
{
    PetscErrorCode ierr = 0;
    PetscReal      r_val; // Temporary for random numbers from [0,1) RNGs
    PetscInt       local_owned_cell_idx_i, local_owned_cell_idx_j, local_owned_cell_idx_k;
    PetscMPIInt    rank_for_logging; // For logging if needed

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_place_in_volume_out = PETSC_FALSE; // Default to: cannot place

    // Default intra-cell logicals and cell node indices (e.g. if placement fails)
    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;
    *ci_metric_lnode_out = xs_gnode; *cj_metric_lnode_out = ys_gnode; *ck_metric_lnode_out = zs_gnode;

    // Calculate number of owned cells in each direction from node counts in info
    // Assumes info->mx, info->my, info->mz are node counts on this process for each dimension.
    // Number of cells = Number of nodes - 1 (if > 0 nodes).
    PetscInt num_owned_cells_i = (info->mx > 1) ? info->mx - 1 : 0;
    PetscInt num_owned_cells_j = (info->my > 1) ? info->my - 1 : 0;
    PetscInt num_owned_cells_k = (info->mz > 1) ? info->mz - 1 : 0;

    if (num_owned_cells_i > 0 && num_owned_cells_j > 0 && num_owned_cells_k > 0) { // If rank owns any 3D cells
        *can_place_in_volume_out = PETSC_TRUE;

        // --- 1. Select a Random Owned Cell ---
        // The selected index will be a 0-based index relative to the start of this rank's owned cells.

        // Select random local owned cell index in I-direction
        ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_i = (PetscInt)(r_val * num_owned_cells_i);
        // Clamp to be safe: local_owned_cell_idx_i should be in [0, num_owned_cells_i - 1]
        local_owned_cell_idx_i = PetscMin(PetscMax(0, local_owned_cell_idx_i), num_owned_cells_i - 1);
        *ci_metric_lnode_out = xs_gnode + local_owned_cell_idx_i; // Convert to local node index for cell origin

        // Select random local owned cell index in J-direction
        ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_j = (PetscInt)(r_val * num_owned_cells_j);
        local_owned_cell_idx_j = PetscMin(PetscMax(0, local_owned_cell_idx_j), num_owned_cells_j - 1);
        *cj_metric_lnode_out = ys_gnode + local_owned_cell_idx_j;

        // Select random local owned cell index in K-direction
        ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_k = (PetscInt)(r_val * num_owned_cells_k);
        local_owned_cell_idx_k = PetscMin(PetscMax(0, local_owned_cell_idx_k), num_owned_cells_k - 1);
        *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_k;

        LOG_ALLOW(LOCAL, LOG_DEBUG, "DVP - Rank %d: Selected Cell (Owned Idx: %d,%d,%d -> LNodeStart: %d,%d,%d). OwnedCells(i,j,k): (%d,%d,%d). GhostNodeStarts(xs,ys,zs): (%d,%d,%d) \n",
                rank_for_logging, local_owned_cell_idx_i, local_owned_cell_idx_j, local_owned_cell_idx_k,
                *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
                num_owned_cells_i, num_owned_cells_j, num_owned_cells_k,
                xs_gnode, ys_gnode, zs_gnode);


        // --- 2. Generate Random Intra-Cell Logical Coordinates [0,1) for MetricLogicalToPhysical ---
        ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out);  CHKERRQ(ierr); // Re-use RNGs
        ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
        ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);

        // Ensure logical coordinates are strictly within [0,1) for robustness with MetricLogicalToPhysical
        *xi_metric_logic_out  = PetscMin(*xi_metric_logic_out,  1.0 - 1.0e-7);
        *eta_metric_logic_out = PetscMin(*eta_metric_logic_out, 1.0 - 1.0e-7);
        *zta_metric_logic_out = PetscMin(*zta_metric_logic_out, 1.0 - 1.0e-7);
        // Ensure they are not negative either (though [0,1) RNGs shouldn't produce this)
        *xi_metric_logic_out  = PetscMax(*xi_metric_logic_out,  0.0);
        *eta_metric_logic_out = PetscMax(*eta_metric_logic_out, 0.0);
        *zta_metric_logic_out = PetscMax(*zta_metric_logic_out, 0.0);

    } else {
        // This rank does not own any 3D cells (e.g., in a 1D or 2D decomposition,
        // or if the global domain itself is not 3D in terms of cells).
        // *can_place_in_volume_out remains PETSC_FALSE.
        LOG_ALLOW(LOCAL, LOG_WARNING, "DVP - Rank %d: Cannot place particle volumetrically. Rank has zero owned cells in at least one dimension (owned cells i,j,k: %d,%d,%d).\n",
                rank_for_logging, num_owned_cells_i, num_owned_cells_j, num_owned_cells_k);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitializeParticleBasicProperties"

/**
 * @brief Initializes basic properties for particles on the local process.
 *
 * This function assigns initial physical positions, Particle IDs (PIDs), and placeholder
 * cell IDs to particles. The method of position initialization depends on
 * `simCtx->ParticleInitialization`:
 * - Mode 0 (Surface): Particles are placed on a designated inlet surface if the
 *   current rank services that surface. Otherwise, they are placed at (0,0,0)
 *   to be migrated to the correct rank later.
 * - Mode 1 (Volumetric): Particles are placed randomly within a cell owned by
 *   the current rank.
 *
 * The logical coordinates for placement are generated using provided random number
 * generators. These logical coordinates are then transformed to physical coordinates
 * using `MetricLogicalToPhysical`.
 *
 * @param user Pointer to the UserCtx structure, containing simulation settings and grid information.
 * @param particlesPerProcess The number of particles to initialize on this MPI rank.
 * @param rand_logic_i Pointer to a PetscRandom generator for the xi logical coordinate.
 * @param rand_logic_j Pointer to a PetscRandom generator for the eta logical coordinate.
 * @param rand_logic_k Pointer to a PetscRandom generator for the zeta logical coordinate.
 * @param bboxlist (Unused in this function for placement) Pointer to the bounding box list;
 *                 provided for API consistency but not used for determining initial positions here.
 *                 Particle positions are determined by logical-to-physical mapping based on rank's owned cells.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
static PetscErrorCode InitializeParticleBasicProperties(UserCtx *user,
                                                PetscInt particlesPerProcess,
                                                PetscRandom *rand_logic_i,
                                                PetscRandom *rand_logic_j,
                                                PetscRandom *rand_logic_k,
                                                BoundingBox *bboxlist) // bboxlist unused for placement
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscReal      *positions_field = NULL; // Pointer to swarm field for physical positions (x,y,z)
    PetscInt64     *particleIDs = NULL;     // Pointer to swarm field for Particle IDs
    PetscInt       *cellIDs_petsc = NULL;   // Pointer to swarm field for DMSwarm_CellID (i,j,k of containing cell)
    PetscInt       *status_field  = NULL;   // Pointer to swarm field for DMSwarm_location_status(NEEDS_LOCATION etc)
    PetscMPIInt    rank,size;                    // MPI rank of the current process, and total number of ranks.
    const Cmpnts   ***coor_nodes_local_array; // Read-only access to local node coordinates (from user->da)
    Vec            Coor_local;              // Local vector for node coordinates
    DMDALocalInfo  info;                    // Local grid information (node-based) from user->da
    PetscInt       xs_gnode_rank, ys_gnode_rank, zs_gnode_rank; // Local starting node indices (incl. ghosts) of rank's DA patch
    PetscInt       IM_nodes_global, JM_nodes_global, KM_nodes_global; // Global node counts in each direction

    // Variables for surface initialization (Mode 0)
    PetscBool      can_this_rank_service_inlet = PETSC_FALSE;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    SimCtx *simCtx = user->simCtx;

    // --- 1. Input Validation and Basic Setup ---
    if (!user || !rand_logic_i || !rand_logic_j || !rand_logic_k) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null user or RNG pointer.");
    }
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);  CHKERRQ(ierr);

    // Get DMDA information for the node-centered coordinate grid (user->da)
    ierr = DMGetCoordinatesLocal(user->da, &Coor_local); CHKERRQ(ierr);
    if (!Coor_local) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "DMGetCoordinatesLocal for user->da returned NULL Coor_local.");
    ierr = DMDAVecGetArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(user->da, &xs_gnode_rank, &ys_gnode_rank, &zs_gnode_rank, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(user->da, NULL, &IM_nodes_global, &JM_nodes_global, &KM_nodes_global, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL); CHKERRQ(ierr);

    // Modification to IM_nodes_global etc. to account for 1-cell halo in each direction.
    IM_nodes_global -= 1; JM_nodes_global -= 1; KM_nodes_global -= 1; 

    const PetscInt IM_cells_global = IM_nodes_global > 0 ? IM_nodes_global - 1 : 0;
    const PetscInt JM_cells_global = JM_nodes_global > 0 ? JM_nodes_global - 1 : 0;
    const PetscInt KM_cells_global = KM_nodes_global > 0 ? KM_nodes_global - 1 : 0;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Initializing %d particles. Mode: %s.\n",
            rank, particlesPerProcess, ParticleInitializationToString(simCtx->ParticleInitialization));

    // --- 2. Pre-computation for Surface Initialization (Mode 0) ---
    if (simCtx->ParticleInitialization == 0 || simCtx->ParticleInitialization == 3) { // Surface initialization
        ierr = CanRankServiceInletFace(user, &info, IM_nodes_global, JM_nodes_global, KM_nodes_global, &can_this_rank_service_inlet); CHKERRQ(ierr);
        if (can_this_rank_service_inlet) {
            LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Will attempt to place particles on inlet face %s.\n", rank, BCFaceToString((BCFace)user->identifiedInletBCFace));
        } else {
            LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Cannot service inlet face %s. Particles will be at default (0,0,0) and rely on migration.\n", rank, BCFaceToString((BCFace)user->identifiedInletBCFace));
        }
    }

    // --- 3. Get Access to Swarm Fields ---
    ierr = DMSwarmGetField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status",NULL,NULL,(void**)&status_field); CHKERRQ(ierr);

    // --- 4. Determine Starting Global PID for this Rank ---
    PetscInt particles_per_rank_ideal = simCtx->np / size; // Assumes user->size is PETSC_COMM_WORLD size
    PetscInt remainder_particles = simCtx->np % size;
    PetscInt base_pid_for_rank = rank * particles_per_rank_ideal + PetscMin(rank, remainder_particles);
    // This calculation must match how particlesPerProcess was determined (e.g., in DistributeParticles).

    // --- 5. Loop Over Particles to Initialize ---
    for (PetscInt p = 0; p < particlesPerProcess; p++) {
        PetscInt idx = p; 
        PetscInt  ci_metric_lnode, cj_metric_lnode, ck_metric_lnode; 
        PetscReal xi_metric_logic, eta_metric_logic, zta_metric_logic; 
        Cmpnts    phys_coords = {0.0, 0.0, 0.0}; 
        PetscBool particle_placed_by_this_rank = PETSC_FALSE; 

        if (simCtx->ParticleInitialization == 0) { // --- 5.a. Surface Random Initialization ---
            if (can_this_rank_service_inlet) {
                ierr = GetRandomCellAndLogicalCoordsOnInletFace(user, &info, xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                        IM_nodes_global, JM_nodes_global, KM_nodes_global,
                                                        rand_logic_i, rand_logic_j, rand_logic_k,
                                                        &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                        &xi_metric_logic, &eta_metric_logic, &zta_metric_logic); CHKERRQ(ierr);
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                            ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                            xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                            &phys_coords); CHKERRQ(ierr);
                particle_placed_by_this_rank = PETSC_TRUE;
            }
        }else if(simCtx->ParticleInitialization == 3) { // --- 5.a1. Custom Initialization DEBUG (Mode 3) ---
            if(can_this_rank_service_inlet) {
                PetscInt64 particle_global_id = (PetscInt64)(base_pid_for_rank + p);
                ierr = GetDeterministicFaceGridLocation(user,&info,xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                        IM_cells_global, JM_cells_global, KM_cells_global,
                                                        particle_global_id,
                                                        &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                        &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                        &particle_placed_by_this_rank); CHKERRQ(ierr);
                
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                        ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                        xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                        &phys_coords); CHKERRQ(ierr);
                                        
            }
        }else if(simCtx->ParticleInitialization == 1){ // --- 5.b. Volumetric Initialization (simCtx->ParticleInitialization == 1) ---
            PetscBool can_place_volumetrically;
            ierr = DetermineVolumetricInitializationParameters(user, &info, xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                            rand_logic_i, rand_logic_j, rand_logic_k,
                                                            &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                            &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                            &can_place_volumetrically); CHKERRQ(ierr);
            if(can_place_volumetrically){
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                        ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                        xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                        &phys_coords); CHKERRQ(ierr);
                particle_placed_by_this_rank = PETSC_TRUE;
            } else {
                LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, p, 1,
                    "Rank %d: PID %lld (idx %ld) (Volumetric Mode %d) - DetermineVolumetric... returned false. Default Phys: (%.2f,%.2f,%.2f).\n",
                    rank, (long long)(base_pid_for_rank + p), (long)p, simCtx->ParticleInitialization, phys_coords.x, phys_coords.y, phys_coords.z);
            }
        }else {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Unknown ParticleInitialization mode %d.", simCtx->ParticleInitialization);
        }

        // --- 5.c. Store Particle Properties ---
        positions_field[3*p+0] = phys_coords.x; 
        positions_field[3*p+1] = phys_coords.y; 
        positions_field[3*p+2] = phys_coords.z;

        particleIDs[p]         = (PetscInt64)base_pid_for_rank + p; 
        cellIDs_petsc[3*p+0]   = -1; cellIDs_petsc[3*p+1] = -1; cellIDs_petsc[3*p+2] = -1;
        status_field[p]        = UNINITIALIZED;

        // --- 5.d. Logging for this particle ---
        if (particle_placed_by_this_rank) {
            LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, idx, user->simCtx->LoggingFrequency,//(particlesPerProcess > 20 ? particlesPerProcess/10 : 1), 
                "Rank %d: PID %lld (idx %ld) PLACED. Mode %s. Embedded Cell:(%d,%d,%d). Logical Coords: (%.2e,%.2f,%.2f).\n Final Coords: (%.6f,%.6f,%.6f).\n",
                rank, (long long)particleIDs[p], (long)p, ParticleInitializationToString(simCtx->ParticleInitialization), 
                ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                xi_metric_logic, eta_metric_logic, zta_metric_logic, 
                phys_coords.x, phys_coords.y, phys_coords.z);
            
        } else {
            LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, idx, user->simCtx->LoggingFrequency, //(particlesPerProcess > 20 ? particlesPerProcess/10 : 1), 
                "Rank %d: PID %lld (idx %ld) Mode %s NOT placed by this rank's logic. Default Coor: (%.2f,%.2f,%.2f). Relies on migration.\n",
                rank, (long long)particleIDs[p], (long)p, ParticleInitializationToString(simCtx->ParticleInitialization), 
                phys_coords.x, phys_coords.y, phys_coords.z);
        }
    } 

    // --- 6. Restore Pointers and Cleanup ---
    ierr = DMSwarmRestoreField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_field); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Completed processing for %d particles.\n",
            rank, particlesPerProcess);

    PROFILE_FUNCTION_END;
            
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitializeSwarmFieldValue"
/**
 * @brief Helper function to Initialize a given particle’s field value.
 *
 * This function performs conditional, point-level Initialization for a swarm field based on its name.
 * For example, you might want to initialize the "velocity" field to 0.0, but the "temperature"
 * field to a nonzero default (e.g., 300.0). This function can be extended for other fields.
 *
 * @param[in] fieldName  Name of the swarm field.
 * @param[in] p          Particle index.
 * @param[in] fieldDim   Dimension of the field.
 * @param[out] fieldData Pointer to the field’s data array.
 *
 * @return PetscErrorCode Returns 0 on success.
 */
static PetscErrorCode InitializeSwarmFieldValue(const char *fieldName, PetscInt p, PetscInt fieldDim, PetscReal *fieldData)
{
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    if (strcmp(fieldName, "velocity") == 0) {
        // For velocity, initialize all components to zero
        for (PetscInt d = 0; d < fieldDim; d++) {
        fieldData[fieldDim * p + d] = 0.0;
        }
    } else if (strcmp(fieldName, "temperature") == 0) {
        // For temperature, for example, initialize to a default value (e.g., 300.0)
        for (PetscInt d = 0; d < fieldDim; d++) {
        fieldData[fieldDim * p + d] = 300.0;
        }
    } else if (strcmp(fieldName, "Psi") == 0) {
        // For pressure, initialize to a default value (e.g., 101325.0)
        for (PetscInt d = 0; d < fieldDim; d++) {
        fieldData[fieldDim * p + d] = 0.0;
        }
    } else {
        // Default: initialize all components to zero
        for (PetscInt d = 0; d < fieldDim; d++) {
        fieldData[fieldDim * p + d] = 0.0;
        }
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "AssignInitialFieldToSwarm"
/**
 * @brief Initializes a generic swarm field with point-level updates.
 *
 * This field-agnostic function retrieves the specified swarm field (which may be
 * scalar or multi-component) and initializes each particle's entry using a helper
 * that performs conditional updates based on the field name.
 *
 * @param[in,out] user       Pointer to the UserCtx structure containing the swarm.
 * @param[in]     fieldName  Name of the swarm field to initialize.
 * @param[in]     fieldDim   Dimension of the field (e.g., 1 for scalar, 3 for vector).
 *
 * @return PetscErrorCode    Returns 0 on success, non-zero on failure.
 */
static PetscErrorCode AssignInitialFieldToSwarm(UserCtx *user, const char *fieldName, PetscInt fieldDim)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscReal     *fieldData = NULL;
    PetscInt       nLocal;

    PetscFunctionBeginUser;
    
    PROFILE_FUNCTION_BEGIN;

    // Get the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "%d local particles found.\n", nLocal);

    // Retrieve the swarm field pointer for the specified fieldName
    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&fieldData); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Retrieved field '%s'.\n", fieldName);

    // Loop over all particles and update the field using the helper function
    for (PetscInt p = 0; p < nLocal; p++) {
        ierr = InitializeSwarmFieldValue(fieldName, p, fieldDim, fieldData); CHKERRQ(ierr);
    PetscReal disp_data[fieldDim];
    
        for (PetscInt d = 0; d < fieldDim; d++) {
    disp_data[d] = fieldData[fieldDim* p + d];
        }
        LOG_LOOP_ALLOW(LOCAL,LOG_DEBUG,p, 100," Particle %d: %s[%d] = [%.6f, ...,%.6f].\n", p,fieldName,fieldDim,disp_data[0],disp_data[fieldDim-1]);
    }
    
    // Restore the swarm field pointer
    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&fieldData); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "Initialization of field '%s' complete.\n", fieldName);
    

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AssignInitialPropertiesToSwarm"

/**
 * @brief Initializes all particle properties in the swarm.
 *
 * This function orchestrates the initialization of particle properties.
 * It first determines the inlet face if surface initialization (Mode 0) is selected
 * by parsing "bcs.dat".
 * Then, it initializes basic particle properties (physical position, Particle ID,
 * and placeholder Cell IDs) by calling `InitializeParticleBasicProperties`. This call
 * uses the provided `rand_logic_i/j/k` RNGs, which must be pre-initialized for [0,1).
 * The `rand_phys_x/y/z` RNGs (physically bounded) are passed but may not be used by
 * `InitializeParticleBasicProperties` for position setting if all initialization paths
 * use logical-to-physical mapping.
 * Finally, it calls helper functions to initialize other registered swarm fields
 * like "velocity", "weight", and "Psi" (scalar) to default values.
 *
 * @param[in,out] user               Pointer to the `UserCtx` structure.
 * @param[in]     particlesPerProcess Number of particles assigned to this MPI process.
 * @param[in]     rand_phys_x        RNG for physical x-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_phys_y        RNG for physical y-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_phys_z        RNG for physical z-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_logic_i       RNG for i-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     rand_logic_j       RNG for j-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     rand_logic_k       RNG for k-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     bboxlist           Array of BoundingBox structures (potentially unused by IPBP).
 *
 * @return PetscErrorCode            Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialPropertiesToSwarm(UserCtx* user,
                                            PetscInt particlesPerProcess,
                                            PetscRandom *rand_phys_x, // RNG from original InitializeRandomGenerators
                                            PetscRandom *rand_phys_y, // RNG from original InitializeRandomGenerators
                                            PetscRandom *rand_phys_z, // RNG from original InitializeRandomGenerators
                                            PetscRandom *rand_logic_i, // RNG from InitializeLogicalSpaceRNGs
                                            PetscRandom *rand_logic_j, // RNG from InitializeLogicalSpaceRNGs
                                            PetscRandom *rand_logic_k, // RNG from InitializeLogicalSpaceRNGs
                                            BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    SimCtx *simCtx = user->simCtx;
    
    // --- 0. Input Validation ---
    if (!user || !bboxlist || !rand_logic_i || !rand_logic_j || !rand_logic_k || !rand_phys_x || !rand_phys_y || !rand_phys_z) {
        // Check all RNGs now as they are passed in
        LOG_ALLOW(GLOBAL, LOG_ERROR, "Null user, bboxlist, or RNG pointer.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null input detected.");
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Initializing swarm with %d particles per process. Mode: %s.\n",
            particlesPerProcess, ParticleInitializationToString(simCtx->ParticleInitialization));

    // --- 1. Parse BCS File for Inlet Information (if Mode 0) ---
    if (simCtx->ParticleInitialization == 0 || simCtx->ParticleInitialization == 3) { // Surface initialization
    if(user->inletFaceDefined == PETSC_FALSE){
    LOG_ALLOW(GLOBAL, LOG_ERROR, "Particle Initialization on inlet surface selected, but no INLET face was identified from bcs.dat. Cannot proceed.\n");
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "ParticleInitialization Mode 0 requires an INLET face to be defined in bcs.dat.");
    }else{
    LOG_ALLOW(GLOBAL, LOG_INFO, "After Parsing BCS file for Inlet, Inlet face = %s\n", BCFaceToString((BCFace)user->identifiedInletBCFace));
    }
    }

    // --- 2. Initialize Basic Particle Properties (Position, PID, Cell IDs placeholder) ---
    // The rand_logic_i/j/k are now passed directly.
    // The rand_phys_x/y/z are passed but InitializeParticleBasicProperties (refactored version)
    // will not use them for setting positions if all its paths use logical-to-physical mapping.
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Calling InitializeParticleBasicProperties.\n");
    ierr = InitializeParticleBasicProperties(user, particlesPerProcess,
                                            rand_logic_i, rand_logic_j, rand_logic_k,
                                            bboxlist); // bboxlist passed along
    CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully initialized basic particle properties.\n");

    // Note: The logical RNGs (rand_logic_i/j/k) are NOT destroyed here.
    // They were created externally (e.g., by InitializeLogicalSpaceRNGs) and
    // should be destroyed externally (e.g., in FinalizeSwarmSetup).
    // Same for rand_phys_x/y/z.

    // --- 3. Initialize Other Swarm Fields (Velocity, Weight, Pressure, etc.) ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'velocity' field.\n");
    ierr = AssignInitialFieldToSwarm(user, "velocity", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "'velocity' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'weight' field.\n");
    ierr = AssignInitialFieldToSwarm(user, "weight", 3); CHKERRQ(ierr); // Assuming weight is vec3
    LOG_ALLOW(LOCAL, LOG_INFO, "'weight' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'P' (Pressure) field.\n");
    ierr = AssignInitialFieldToSwarm(user, "Psi", 1); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'P' field initialization complete.\n");
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully completed all swarm property initialization.\n");


    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DistributeParticles"
/**
 * @brief Distributes particles evenly across MPI processes, handling any remainders.
 *
 * This function calculates the number of particles each MPI process should handle,
 * distributing the remainder particles to the first few ranks if necessary.
 *
 * @param[in]     numParticles       Total number of particles to create across all MPI processes.
 * @param[in]     rank               MPI rank of the current process.
 * @param[in]     size               Total number of MPI processes.
 * @param[out]    particlesPerProcess Number of particles assigned to the current MPI process.
 * @param[out]    remainder           Remainder particles when dividing numParticles by size.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DistributeParticles(PetscInt numParticles, PetscMPIInt rank, PetscMPIInt size, PetscInt* particlesPerProcess, PetscInt* remainder) {

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;
    // Calculate the base number of particles per process
    *particlesPerProcess = numParticles / size;
    *remainder = numParticles % size;

    // Distribute the remainder particles to the first 'remainder' ranks
    if (rank < *remainder) {
        *particlesPerProcess += 1;
        LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"Rank %d receives an extra particle. Total: %d\n", rank, *particlesPerProcess);
    } else {
        LOG_ALLOW_SYNC(GLOBAL,LOG_INFO, "Rank %d receives %d particles.\n", rank, *particlesPerProcess);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FinalizeSwarmSetup"
/**
 * @brief Finalizes the swarm setup by destroying random generators and logging completion.
 *
 * This function cleans up resources by destroying random number generators and LOG_ALLOWs the completion of swarm setup.
 *
 * @param[in]     randx             Random number generator for the x-coordinate.
 * @param[in]     randy             Random number generator for the y-coordinate.
 * @param[in]     randz             Random number generator for the z-coordinate.
 * @param[in]     rand_logic_i      Random number generator for the xi-coordinate.
 * @param[in]     rand_logic_j      Random number generator for the eta-coordinate.
 * @param[in]     rand_logic_k      Random number generator for the zeta-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSwarmSetup(PetscRandom *randx, PetscRandom *randy, PetscRandom *randz, PetscRandom *rand_logic_i, PetscRandom *rand_logic_j, PetscRandom *rand_logic_k) {
    PetscErrorCode ierr;  // Error code for PETSc functions
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // Destroy random number generators to free resources
    // Physical space
    ierr = PetscRandomDestroy(randx); CHKERRQ(ierr);
    ierr = PetscRandomDestroy(randy); CHKERRQ(ierr);
    ierr = PetscRandomDestroy(randz); CHKERRQ(ierr);
    // Logical space
    ierr = PetscRandomDestroy(rand_logic_i); CHKERRQ(ierr);
    ierr = PetscRandomDestroy(rand_logic_j); CHKERRQ(ierr);
    ierr = PetscRandomDestroy(rand_logic_k); CHKERRQ(ierr);      
    
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Destroyed all random number generators.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CreateParticleSwarm"
/**
 * @brief Creates and initializes a Particle Swarm.
 *
 * This function sets up a DMSwarm within the provided UserCtx structure, initializes
 * particle fields, and distributes particles across MPI processes. It ensures that
 * the number of particles is evenly divided among the available MPI ranks. If the total
 * number of particles isn't divisible by the number of processes, the remainder is distributed
 * to the first few ranks.
 *.
*
* @param[in,out] user          Pointer to the UserCtx structure containing the simulation context.
* @param[in]     numParticles  Total number of particles to create across all MPI processes.
* @param[in]     bboxlist      Pointer to an array of BoundingBox structures, one per rank.
*
* @param[in]     particlesPerProcess   
* @return PetscErrorCode Returns 0 on success, non-zero on failure.
*
* @note
* - Ensure that `numParticles` is a positive integer.
* - The `control.dat` file should contain necessary PETSc options.
* - The `bboxlist` array should be properly populated before calling this function.
*/
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles, PetscInt *particlesPerProcess, BoundingBox *bboxlist) {
    PetscErrorCode ierr;                      // PETSc error handling variable
    PetscMPIInt rank, size;                   // Variables to store MPI rank and size
    PetscInt remainder = 0;                   // Remainder of particles after division
    
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // Validate input parameters
    if (numParticles <= 0) {
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Number of particles must be positive. Given: %d\n", numParticles);
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    
    // Retrieve MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO," Domain dimensions: [%.2f,%.2f],[%.2f,%.2f],[%.2f,%.2f] \n", 
        user->Min_X,user->Max_X,user->Min_Y,user->Max_Y, user->Min_Z,user->Max_Z);
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "[Rank %d] Local Bounding Box: [%.2f,%.2f],[%.2f,%.2f],[%.2f,%.2f] \n", 
        rank,user->bbox.min_coords.x,user->bbox.max_coords.x,
        user->bbox.min_coords.y,user->bbox.max_coords.y,
        user->bbox.min_coords.z,user->bbox.max_coords.z);
    // Distribute particles among MPI processes
    ierr = DistributeParticles(numParticles, rank, size, particlesPerProcess, &remainder); CHKERRQ(ierr);

    // Initialize the DMSwarm - creates the swarm, sets the type and dimension
    ierr = InitializeSwarm(user); CHKERRQ(ierr);

    if (user->da) {
    ierr = DMSwarmSetCellDM(user->swarm, user->da); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO,"Associated DMSwarm with Cell DM (user->da).\n");
    } else {
    // If user->da is essential for your simulation logic with particles, this should be a fatal error.
    LOG_ALLOW(GLOBAL, LOG_WARNING, "user->da (Cell DM for Swarm) is NULL. Cell-based swarm operations might fail.\n");
    // SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "user->da (Cell DM) is NULL but required.");
    }
    
    // Register particle fields (position, velocity, CellID, weight, etc.)
    ierr = RegisterParticleFields(user->swarm); CHKERRQ(ierr);

    // Set the local number of particles for this rank and additional buffer for particle migration
    ierr = DMSwarmSetLocalSizes(user->swarm, *particlesPerProcess, numParticles); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO, "Set local swarm size: %d particles.\n", *particlesPerProcess);

    // Optionally, LOG_ALLOW detailed DM info in debug mode
    if (get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) {
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Viewing DMSwarm:\n");
        ierr = DMView(user->swarm, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL,LOG_INFO, "Particle swarm creation and initialization complete.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "UnpackSwarmFields"

/**
 * @brief Initializes a Particle struct with data from DMSwarm fields.
 *
 * This helper function populates a Particle structure using data retrieved from DMSwarm fields.
 *
 * @param[in]     i            Index of the particle in the DMSwarm.
 * @param[in]     PIDs         Pointer to the array of particle IDs.
 * @param[in]     weights      Pointer to the array of particle weights.
 * @param[in]     positions    Pointer to the array of particle positions.
 * @param[in]     cellIndices  Pointer to the array of particle cell indices.
 * @param[in]     velocities   Pointer to the array of particle velocities.
 * @param[in]     LocStatus    Pointer to the array of cell location status indicators.
 * @param[out]    particle     Pointer to the Particle struct to initialize.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode UnpackSwarmFields(PetscInt i, const PetscInt64 *PIDs, const PetscReal *weights,
                const PetscReal *positions, const PetscInt *cellIndices,
                PetscReal *velocities,PetscInt *LocStatus,Particle *particle) {
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    PetscMPIInt rank;
    PetscErrorCode ierr;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output Particle pointer is NULL. \n");
    }
    
    // logging the start of particle initialization
    LOG_ALLOW(LOCAL,LOG_INFO, "[Rank %d]Unpacking Particle [%d] with PID: %ld.\n",rank, i, PIDs[i]);
    
    // Initialize PID
    particle->PID = PIDs[i];
    LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d]Particle [%d] PID set to: %ld.\n", rank,i, particle->PID);
    
    // Initialize weights
    particle->weights.x = weights[3 * i];
    particle->weights.y = weights[3 * i + 1];
    particle->weights.z = weights[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d]Particle [%d] weights set to: (%.6f, %.6f, %.6f).\n", 
        rank,i, particle->weights.x, particle->weights.y, particle->weights.z);
    
    // Initialize locations
    particle->loc.x = positions[3 * i];
    particle->loc.y = positions[3 * i + 1];
    particle->loc.z = positions[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d]Particle [%d] location set to: (%.6f, %.6f, %.6f).\n", 
        rank,i, particle->loc.x, particle->loc.y, particle->loc.z);
    
    // Initialize velocities (assuming default zero; modify if necessary)
    particle->vel.x = velocities[3 * i];
    particle->vel.y = velocities[3 * i + 1];
    particle->vel.z = velocities[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d]Particle [%d] velocities unpacked to: [%.6f,%.6f,%.6f].\n",rank, i,particle->vel.x,particle->vel.y,particle->vel.z);
    
    // Initialize cell indices
    particle->cell[0] = cellIndices[3 * i];
    particle->cell[1] = cellIndices[3 * i + 1];
    particle->cell[2] = cellIndices[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d]Particle [%d] cell indices set to: [%d, %d, %d].\n",rank,i, particle->cell[0], particle->cell[1], particle->cell[2]);

    particle->location_status = (ParticleLocationStatus)LocStatus[i];
    LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d]Particle [%d] Status set to: %d.\n",rank, i, particle->location_status);

    // The destination_rank is only set by the location search, not read from the swarm,
    // so we initialize it to a known invalid state.
    particle->destination_rank = MPI_PROC_NULL;
    
    // logging the completion of particle initialization
    LOG_ALLOW(LOCAL,LOG_INFO,"[Rank %d]Completed initialization of Particle [%d]. \n", rank,i);
    

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "UpdateSwarmFields"

/**
 * @brief Updates DMSwarm fields with data from a Particle struct.
 *
 * This helper function writes back the modified Particle data to the corresponding DMSwarm fields.
 *
 * @param[in] i            Index of the particle in the DMSwarm.
 * @param[in] particle     Pointer to the Particle struct containing updated data.
 * @param[in,out] weights  Pointer to the array of particle weights.
 * @param[in,out] cellIndices Pointer to the array of particle cell indices.
 * @param[in,out] LocStatus   Pointer to the array of cell location status indicators.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode UpdateSwarmFields(PetscInt i, const Particle *particle,
                PetscReal *weights, PetscInt *cellIndices, PetscInt *status_field) {
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input Particle pointer is NULL.\n");
    }
    
    // logging the start of swarm fields update
    LOG_ALLOW(LOCAL,LOG_INFO,"Updating DMSwarm fields for Particle [%d].\n", i);
    
    // Update weights
    weights[3 * i]     = particle->weights.x;
    weights[3 * i + 1] = particle->weights.y;
    weights[3 * i + 2] = particle->weights.z;
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Updated weights for Particle [%d]: (%.6f, %.6f, %.6f).\n", 
        i, weights[3 * i], weights[3 * i + 1], weights[3 * i + 2]);
    
    // Update cell indices
    cellIndices[3 * i]     = particle->cell[0];
    cellIndices[3 * i + 1] = particle->cell[1];
    cellIndices[3 * i + 2] = particle->cell[2];
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Updated cell indices for Particle [%d]: [%d, %d, %d].\n", 
        i, cellIndices[3 * i], cellIndices[3 * i + 1], cellIndices[3 * i + 2]);

    status_field[i] = particle->location_status;
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Updated location status for Particle [%d]: [%d].\n", 
        i, status_field[i]);
    
    // logging the completion of swarm fields update
    LOG_ALLOW(LOCAL,LOG_INFO,"Completed updating DMSwarm fields for Particle [%d].\n", i);
    

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "IsParticleInsideBoundingBox"
/**
 * @brief Checks if a particle's location is within a specified bounding box.
 *
 * This function determines whether the given particle's location lies inside the provided bounding box.
 * It performs an axis-aligned bounding box (AABB) check by comparing the particle's coordinates to the
 * minimum and maximum coordinates of the bounding box in each dimension (x, y, z).
 *
 * logging statements are included to provide detailed information about the function's execution.
 *
 * @param[in]  bbox     Pointer to the BoundingBox structure containing minimum and maximum coordinates.
 * @param[in]  particle Pointer to the Particle structure containing the particle's location and identifier.
 *
 * @return PetscBool    Returns `PETSC_TRUE` if the particle is inside the bounding box, `PETSC_FALSE` otherwise.
 *
 * @note
 * - The function assumes that the `bbox` and `particle` pointers are valid and non-NULL.
 * - The function includes logging statements that start with the function name.
 * - The `LOG_ALLOW_SCOPE` variable is used to distinguish between `GLOBAL` and `LOCAL` LOG_ALLOW outputs.
 * - Be cautious when logging in performance-critical code sections, especially if the function is called frequently.
 */
PetscBool IsParticleInsideBoundingBox(const BoundingBox *bbox, const Particle *particle)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Validate input pointers
    if (!bbox) {
        // LOG_ALLOW error message and return PETSC_FALSE
    LOG_ALLOW(LOCAL,LOG_ERROR, "Error - 'bbox' pointer is NULL.");
        PROFILE_FUNCTION_END;
        return PETSC_FALSE;
    }
    if (!particle) {
    LOG_ALLOW(LOCAL,LOG_ERROR,"Error - 'particle' pointer is NULL.");
        PROFILE_FUNCTION_END;
        return PETSC_FALSE;
    }

    // Extract particle location and bounding box coordinates
    const Cmpnts loc = particle->loc;
    const Cmpnts min_coords = bbox->min_coords;
    const Cmpnts max_coords = bbox->max_coords;

    // LOG_ALLOW the particle location and bounding box coordinates for debugging
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Particle PID %ld location: (%.6f, %.6f, %.6f).\n",particle->PID, loc.x, loc.y, loc.z);
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "BoundingBox min_coords: (%.6f, %.6f, %.6f), max_coords: (%.6f, %.6f, %.6f).\n", 
    min_coords.x, min_coords.y, min_coords.z, max_coords.x, max_coords.y, max_coords.z);

    // Check if the particle's location is within the bounding box
    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        // Particle is inside the bounding box
        LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "Particle PID %ld is inside the bounding box.\n",particle->PID);
        PROFILE_FUNCTION_END;
        return PETSC_TRUE;
    }

    // Particle is outside the bounding box
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG,"Particle PID %ld is outside the bounding box.\n",particle->PID);
    PROFILE_FUNCTION_END;
    return PETSC_FALSE;
}


#undef __FUNCT__
#define __FUNCT__ "UpdateParticleWeights"
/**
 * @brief Updates a particle's interpolation weights based on distances to cell faces.
 *
 * This function computes interpolation weights using distances to the six
 * cell faces (`d`) and updates the `weight` field of the provided particle.
 *
 * @param[in]  d        Pointer to an array of distances to the six cell faces.
 * @param[out] particle Pointer to the Particle structure whose weights are to be updated.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode UpdateParticleWeights(PetscReal *d, Particle *particle) {

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Validate input pointers
    if (!d || !particle) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                "Null pointer argument (d or particle).");
    }


    // Validate distances
    for (PetscInt i = LEFT; i < NUM_FACES; i++) {
        if (d[i] <= INTERPOLATION_DISTANCE_TOLERANCE) {
            LOG_ALLOW(LOCAL, LOG_WARNING,
                "face distance d[%d] = %f <= %f; "
                "clamping to 1e-14 to avoid zero/negative.\n",
                i, (double)d[i], INTERPOLATION_DISTANCE_TOLERANCE);
            d[i] = INTERPOLATION_DISTANCE_TOLERANCE;
        }
    }

    // LOG_ALLOW the input distances
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "Calculating weights with distances: "
        "[LEFT=%f, RIGHT=%f, BOTTOM=%f, TOP=%f, FRONT=%f, BACK=%f].\n",
        d[LEFT], d[RIGHT], d[BOTTOM], d[TOP], d[FRONT], d[BACK]);

    // Compute and update the particle's weights
    particle->weights.x = d[LEFT] / (d[LEFT] + d[RIGHT]);
    particle->weights.y = d[BOTTOM] / (d[BOTTOM] + d[TOP]);
    particle->weights.z = d[BACK] / (d[FRONT] + d[BACK]);

    // LOG_ALLOW the updated weights
    LOG_ALLOW(LOCAL,LOG_DEBUG,
        "Updated particle weights: x=%f, y=%f, z=%f.\n",
        particle->weights.x, particle->weights.y, particle->weights.z);


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Initializes or loads the particle swarm based on the simulation context.
 *
 * This function is the central point for setting up the DMSwarm. Its behavior
 * depends on the simulation context (simCtx):
 *
 * 1.  **Fresh Start (simCtx->StartStep == 0):** A new particle population is
 *     generated according to the specified initial conditions.
 *
 * 2.  **Restart (simCtx->StartStep > 0):**
 *     - If `simCtx->particleRestartMode` is "init", a new particle population
 *       is generated, just like a fresh start. This allows injecting fresh
 *       particles into a pre-computed flow field.
 *     - If `simCtx->particleRestartMode` is "load", the particle state is loaded
 *       from restart files corresponding to the StartStep.
 *
 * @param[in,out] simCtx Pointer to the main SimulationContext, which contains all
 *                       configuration and provides access to the UserCtx.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
#undef __FUNCT__
#define __FUNCT__ "InitializeParticleSwarm"
PetscErrorCode InitializeParticleSwarm(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscInt       particlesPerProcess = 0;
    UserCtx       *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting particle swarm setup for %d particles.\n", simCtx->np);

    // --- Phase 1: Create the DMSwarm Object (Always required) ---
    // This creates the container and registers the fields. It does not add particles yet.
    ierr = CreateParticleSwarm(user, simCtx->np, &particlesPerProcess, simCtx->bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "DMSwarm object and fields created successfully.\n");


    // --- Phase 2: Decide whether to Initialize new particles or Load existing ones ---
    PetscBool should_initialize_new_particles = PETSC_FALSE;
    if(simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR){
        should_initialize_new_particles = PETSC_TRUE;
    }else{
        if (simCtx->StartStep == 0) {
            should_initialize_new_particles = PETSC_TRUE; // Standard fresh start
        } else {
        // It's a restart, so check the user's requested particle mode.
            if (strcmp(simCtx->particleRestartMode, "init") == 0) {
                should_initialize_new_particles = PETSC_TRUE; // User wants to re-initialize particles in a restarted flow.
            }
        }
    }

    // --- Phase 3: Execute the chosen particle setup path ---
    if (should_initialize_new_particles) {
        // --- PATH A: Generate a fresh population of particles ---
        LOG_ALLOW(GLOBAL, LOG_INFO, "Mode: INITIALIZE. Generating new particle population.\n");
        PetscRandom randx, randy, randz;
        PetscRandom rand_logic_i, rand_logic_j, rand_logic_k;

        ierr = InitializeRandomGenerators(user, &randx, &randy, &randz); CHKERRQ(ierr);
        ierr = InitializeLogicalSpaceRNGs(&rand_logic_i, &rand_logic_j, &rand_logic_k); CHKERRQ(ierr);
        ierr = AssignInitialPropertiesToSwarm(user, particlesPerProcess, &randx, &randy, &randz, &rand_logic_i, &rand_logic_j, &rand_logic_k, simCtx->bboxlist); CHKERRQ(ierr);
        ierr = FinalizeSwarmSetup(&randx, &randy, &randz, &rand_logic_i, &rand_logic_j, &rand_logic_k); CHKERRQ(ierr);

    } else {
        // --- PATH B: Load particle population from restart files ---
        // This path is only taken if simCtx->StartStep > 0 AND simCtx->particleRestartMode == "load"
        LOG_ALLOW(GLOBAL, LOG_INFO, "Mode: LOAD. Loading particle population from files for step %d.\n", simCtx->StartStep);

        ierr = PreCheckAndResizeSwarm(user, simCtx->StartStep, "dat"); CHKERRQ(ierr);

        ierr = ReadAllSwarmFields(user, simCtx->StartStep); CHKERRQ(ierr);
        // Note: We check for file-open errors inside ReadAllSwarmFields now.

        LOG_ALLOW(GLOBAL, LOG_INFO, "Particle data loaded. CellID and status are preserved from file.\n");
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "Particle swarm setup complete.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

