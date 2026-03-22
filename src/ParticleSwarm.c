    // ParticleSwarm.c

#include "ParticleSwarm.h"

#define INTERPOLATION_DISTANCE_TOLERANCE 1.0e-14

#undef __FUNCT__
#define __FUNCT__ "InitializeSwarm"
/**
 * @brief Implementation of \ref InitializeSwarm().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleSwarm.h`.
 * @see InitializeSwarm()
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
 * @brief Internal helper implementation: `RegisterSwarmField()`.
 * @details Local to this translation unit.
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
 * @brief Implementation of \ref RegisterParticleFields().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleSwarm.h`.
 * @see RegisterParticleFields()
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

    ierr = RegisterSwarmField(swarm,"Diffusivity", 1,PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'Diffusivity' - Scalar.\n");

    ierr = RegisterSwarmField(swarm,"DiffusivityGradient", 3,PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field 'DiffusivityGradient' - Vector.\n");

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
 * @brief Internal helper implementation: `DetermineVolumetricInitializationParameters()`.
 * @details Local to this translation unit.
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
    (void)user;
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
        // Get number of cells this rank owns in each dimension (tangential to the face mainly)
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;

    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    if (num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) { // If rank owns any 3D cells
        *can_place_in_volume_out = PETSC_TRUE;

        // --- 1. Select a Random Owned Cell ---
        // The selected index will be a 0-based index relative to the start of this rank's owned cells.

        // Select random local owned cell index in I-direction
        ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_i = (PetscInt)(r_val * num_owned_cells_on_rank_i);
        // Clamp to be safe: local_owned_cell_idx_i should be in [0, num_owned_cells_on_rank_i - 1]
        local_owned_cell_idx_i = PetscMin(PetscMax(0, local_owned_cell_idx_i), num_owned_cells_on_rank_i - 1);
        *ci_metric_lnode_out = xs_gnode + local_owned_cell_idx_i; // Convert to local node index for cell origin

        // Select random local owned cell index in J-direction
        ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_j = (PetscInt)(r_val * num_owned_cells_on_rank_j);
        local_owned_cell_idx_j = PetscMin(PetscMax(0, local_owned_cell_idx_j), num_owned_cells_on_rank_j - 1);
        *cj_metric_lnode_out = ys_gnode + local_owned_cell_idx_j;

        // Select random local owned cell index in K-direction
        ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_k = (PetscInt)(r_val * num_owned_cells_on_rank_k);
        local_owned_cell_idx_k = PetscMin(PetscMax(0, local_owned_cell_idx_k), num_owned_cells_on_rank_k - 1);
        *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_k;

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Selected Cell (Owned Idx: %d,%d,%d -> LNodeStart: %d,%d,%d). OwnedCells(i,j,k): (%d,%d,%d). GhostNodeStarts(xs,ys,zs): (%d,%d,%d) \n",
                rank_for_logging, local_owned_cell_idx_i, local_owned_cell_idx_j, local_owned_cell_idx_k,
                *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
                num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k,
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
        LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Cannot place particle volumetrically. Rank has zero owned cells in at least one dimension (owned cells i,j,k: %d,%d,%d).\n",
                rank_for_logging, num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitializeParticleBasicProperties"

/**
 * @brief Internal helper implementation: `InitializeParticleBasicProperties()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode InitializeParticleBasicProperties(UserCtx *user,
                                                PetscInt particlesPerProcess,
                                                PetscRandom *rand_logic_i,
                                                PetscRandom *rand_logic_j,
                                                PetscRandom *rand_logic_k,
                                                BoundingBox *bboxlist) // bboxlist unused for placement
{
    PetscErrorCode ierr;
    (void)bboxlist;
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
    ierr = DMDAGetCorners(user->da, &xs_gnode_rank, &ys_gnode_rank, &zs_gnode_rank, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(user->da, NULL, &IM_nodes_global, &JM_nodes_global, &KM_nodes_global, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL); CHKERRQ(ierr);

    // Modification to IM_nodes_global etc. to account for 1-cell halo in each direction.
    IM_nodes_global -= 1; JM_nodes_global -= 1; KM_nodes_global -= 1; 

    const PetscInt IM_cells_global = IM_nodes_global > 0 ? IM_nodes_global - 1 : 0;
    const PetscInt JM_cells_global = JM_nodes_global > 0 ? JM_nodes_global - 1 : 0;
    const PetscInt KM_cells_global = KM_nodes_global > 0 ? KM_nodes_global - 1 : 0;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Initializing %d particles. Mode: %s.\n",
            rank, particlesPerProcess, ParticleInitializationToString(simCtx->ParticleInitialization));

    // --- 2. Pre-computation for Surface Initialization (PARTICLE_INIT_SURFACE_RANDOM and PARTICLE_INIT_SURFACE_EDGES) ---
    if (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM ||
        simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_EDGES) { // Surface initialization
        ierr = CanRankServiceInletFace(user, &info, IM_nodes_global, JM_nodes_global, KM_nodes_global, &can_this_rank_service_inlet); CHKERRQ(ierr);
        if (can_this_rank_service_inlet) {
            LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Will attempt to place particles on inlet face %s.\n", rank, BCFaceToString((BCFace)user->identifiedInletBCFace));
        } else {
            LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Cannot service inlet face %s. Particles will be at Inlet Center (%.6f,%.6f,%.6f) and rely on migration.\n", rank, BCFaceToString((BCFace)user->identifiedInletBCFace),user->simCtx->CMx_c,user->simCtx->CMy_c,user->simCtx->CMz_c);
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

        if (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM) { // --- 5.a. Surface Random Initialization ---
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
            }else{
                // Rank cannot service inlet - place at inlet center to be migrated later
                phys_coords.x = user->simCtx->CMx_c;
                phys_coords.y = user->simCtx->CMy_c;
                phys_coords.z = user->simCtx->CMz_c;
                particle_placed_by_this_rank = PETSC_FALSE; // Relies on migration
            }
        }else if(simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_EDGES) { // --- 5.a1. Surface Edges Initialization (deterministic) ---
            if(can_this_rank_service_inlet) {
                PetscInt64 particle_global_id = (PetscInt64)(base_pid_for_rank + p);
                ierr = GetDeterministicFaceGridLocation(user,&info,xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                        IM_cells_global, JM_cells_global, KM_cells_global,
                                                        particle_global_id,
                                                        &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                        &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                        &particle_placed_by_this_rank); CHKERRQ(ierr);
                if(particle_placed_by_this_rank){
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                        ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                        xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                        &phys_coords); CHKERRQ(ierr);
                }else{
                    // Even if rank can service face, it may not own the portion of the face that a particle is placed in.
                    phys_coords.x = user->simCtx->CMx_c;
                    phys_coords.y = user->simCtx->CMy_c;
                    phys_coords.z = user->simCtx->CMz_c;
                }
            }else{
                // Rank cannot service inlet - place at inlet center to be migrated later
                phys_coords.x = user->simCtx->CMx_c;
                phys_coords.y = user->simCtx->CMy_c;
                phys_coords.z = user->simCtx->CMz_c;
                particle_placed_by_this_rank = PETSC_FALSE; // Relies on migration
            }
        }else if(simCtx->ParticleInitialization == PARTICLE_INIT_VOLUME){ // --- 5.b. Volumetric Initialization ---
            PetscBool can_place_volumetrically = PETSC_FALSE;
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
                    "Rank %d: PID %lld (idx %ld) (Volumetric Mode) - DetermineVolumetric... returned false. Default Phys: (%.2f,%.2f,%.2f).\n",
                    rank, (long long)(base_pid_for_rank + p), (long)p, phys_coords.x, phys_coords.y, phys_coords.z);
            }
        }else if(simCtx->ParticleInitialization == PARTICLE_INIT_POINT_SOURCE){ // --- 5.c. Point Source Initialization ---
            // All particles placed at the user-specified fixed point (psrc_x, psrc_y, psrc_z).
            // No random number generation or logical-to-physical conversion needed.
            phys_coords.x = simCtx->psrc_x;
            phys_coords.y = simCtx->psrc_y;
            phys_coords.z = simCtx->psrc_z;
            particle_placed_by_this_rank = PETSC_TRUE;
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
            LOG_LOOP_ALLOW(LOCAL, LOG_VERBOSE, idx, user->simCtx->LoggingFrequency,//(particlesPerProcess > 20 ? particlesPerProcess/10 : 1), 
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
 * @brief Internal helper implementation: `InitializeSwarmFieldValue()`.
 * @details Local to this translation unit.
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
    } else if (strcmp(fieldName, "Diffusivity") == 0) {
        // For Diffusivity,initialize to a default value (e.g., 1.0)
        for (PetscInt d = 0; d < fieldDim; d++) {
        fieldData[fieldDim * p + d] = 1.0;
        }
    } else if (strcmp(fieldName, "DiffusivityGradient") == 0) {
        // For DiffusivityGradient,initialize to a default value (e.g., 1.0)
        for (PetscInt d = 0; d < fieldDim; d++) {
        fieldData[fieldDim * p + d] = 1.0;
        }
    } else if (strcmp(fieldName, "Psi") == 0) {
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
 * @brief Internal helper implementation: `AssignInitialFieldToSwarm()`.
 * @details Local to this translation unit.
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
        LOG_LOOP_ALLOW(LOCAL,LOG_VERBOSE,p, 100," Particle %d: %s[%d] = [%.6f, ...,%.6f].\n", p,fieldName,fieldDim,disp_data[0],disp_data[fieldDim-1]);
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
 * @brief Internal helper implementation: `AssignInitialPropertiesToSwarm()`.
 * @details Local to this translation unit.
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

    // --- 1. Parse BCS File for Inlet Information (if surface initialization) ---
    if (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM ||
        simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_EDGES) { // Surface initialization
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

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'Diffusivity' field.\n");
    ierr = AssignInitialFieldToSwarm(user, "Diffusivity", 1); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'Diffusivity' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'DiffusivityGradient' field.\n");
    ierr = AssignInitialFieldToSwarm(user, "DiffusivityGradient", 3); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'DiffusivityGradient' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'Psi' (Scalar) field.\n");
    ierr = AssignInitialFieldToSwarm(user, "Psi", 1); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'P' field initialization complete.\n");
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully completed all swarm property initialization.\n");


    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DistributeParticles"
/**
 * @brief Implementation of \ref DistributeParticles().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleSwarm.h`.
 * @see DistributeParticles()
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
 * @brief Implementation of \ref FinalizeSwarmSetup().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleSwarm.h`.
 * @see FinalizeSwarmSetup()
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
 * @brief Internal helper implementation: `CreateParticleSwarm()`.
 * @details Local to this translation unit.
 */
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles, PetscInt *particlesPerProcess, BoundingBox *bboxlist) {
    PetscErrorCode ierr;                      // PETSc error handling variable
    (void)bboxlist;
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

// NOTE: The following two functions are helpers for unpacking and updating particle data
// between DMSwarm fields and the Particle struct used in simulation logic.
// While Swarm fields store data in arrays, the Particle struct provides a convenient
// way to manipulate individual particle properties during simulation steps.

#undef __FUNCT__
#define __FUNCT__ "UnpackSwarmFields"

/**
 * @brief Implementation of \ref UnpackSwarmFields().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleSwarm.h`.
 * @see UnpackSwarmFields()
 */
PetscErrorCode UnpackSwarmFields(PetscInt i, const PetscInt64 *PIDs, const PetscReal *weights,
                const PetscReal *positions, const PetscInt *cellIndices,
                PetscReal *velocities,PetscInt *LocStatus,PetscReal *diffusivity, Cmpnts *diffusivitygradient, PetscReal *psi, Particle *particle) {
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    PetscMPIInt rank;
    PetscErrorCode ierr;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output Particle pointer is NULL. \n");
    }
    
    // logging the start of particle initialization
    LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d]Unpacking Particle [%d] with PID: %ld.\n",rank, i, PIDs[i]);
    
    // Initialize PID
    if(PIDs == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input PIDs pointer is NULL.\n");
    }
    particle->PID = PIDs[i];
    LOG_ALLOW(LOCAL,LOG_VERBOSE, "[Rank %d]Particle [%d] PID set to: %ld.\n", rank,i, particle->PID);
    
    // Initialize weights
    if(weights == NULL){
        particle->weights.x = 1.0;
        particle->weights.y = 1.0;
        particle->weights.z = 1.0;
        LOG_ALLOW(LOCAL,LOG_WARNING, "[Rank %d]Particle [%d] weights pointer is NULL. Defaulting weights to (1.0, 1.0, 1.0).\n", rank,i);
    }else{
    particle->weights.x = weights[3 * i];
    particle->weights.y = weights[3 * i + 1];
    particle->weights.z = weights[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_VERBOSE, "[Rank %d]Particle [%d] weights set to: (%.6f, %.6f, %.6f).\n", 
        rank,i, particle->weights.x, particle->weights.y, particle->weights.z);
    }
    // Initialize locations
    if(positions == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input positions pointer is NULL.\n");
    }
    particle->loc.x = positions[3 * i];
    particle->loc.y = positions[3 * i + 1];
    particle->loc.z = positions[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_VERBOSE, "[Rank %d]Particle [%d] location set to: (%.6f, %.6f, %.6f).\n", 
        rank,i, particle->loc.x, particle->loc.y, particle->loc.z);
    
    // Initialize velocities (assuming default zero; modify if necessary)
    if(velocities == NULL){
        particle->vel.x = 0.0;
        particle->vel.y = 0.0;
        particle->vel.z = 0.0;
        LOG_ALLOW(LOCAL,LOG_WARNING, "[Rank %d]Particle [%d] velocities pointer is NULL. Defaulting velocities to (0.0, 0.0, 0.0).\n", rank,i);
    }else{
    particle->vel.x = velocities[3 * i];
    particle->vel.y = velocities[3 * i + 1];
    particle->vel.z = velocities[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_VERBOSE,"[Rank %d]Particle [%d] velocities unpacked to: [%.6f,%.6f,%.6f].\n",rank, i,particle->vel.x,particle->vel.y,particle->vel.z);
    }
    
    // Initialize diffusivity
    if(diffusivity == NULL){
        particle->diffusivity = 1.0; // Default diffusivity
    }else{
        particle->diffusivity = diffusivity[i];
    }
    LOG_ALLOW(LOCAL,LOG_VERBOSE,"[Rank %d]Particle [%d] diffusivity set to: %.6f.\n",rank,i, particle->diffusivity);
    
    // Initialize diffusivity gradient
    if(diffusivitygradient == NULL){
        particle->diffusivitygradient.x = 0.0;
        particle->diffusivitygradient.y = 0.0;
        particle->diffusivitygradient.z = 0.0;
        LOG_ALLOW(LOCAL,LOG_WARNING, "[Rank %d]Particle [%d] diffusivity gradient pointer is NULL. Defaulting to (0.0, 0.0, 0.0).\n", rank,i);
    }else{
        particle->diffusivitygradient.x = diffusivitygradient[i].x;
        particle->diffusivitygradient.y = diffusivitygradient[i].y;
        particle->diffusivitygradient.z = diffusivitygradient[i].z;
    }
    LOG_ALLOW(LOCAL,LOG_VERBOSE,"[Rank %d]Particle [%d] diffusivity gradient set to: (%.6f, %.6f, %.6f).\n", rank,i,particle->diffusivitygradient.x,particle->diffusivitygradient.y,particle->diffusivitygradient.z);

    // Initialize psi
    if(psi == NULL){
        particle->psi = 0.0; // Default psi
    }else{
    particle->psi = psi[i];
    }        
    LOG_ALLOW(LOCAL,LOG_VERBOSE,"[Rank %d]Particle [%d] psi set to: %.6f.\n",rank,i, particle->psi);
    
    // Initialize cell indices
    if(cellIndices == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input cellIndices pointer is NULL.\n");
    }
    particle->cell[0] = cellIndices[3 * i];
    particle->cell[1] = cellIndices[3 * i + 1];
    particle->cell[2] = cellIndices[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_VERBOSE,"[Rank %d]Particle [%d] cell indices set to: [%d, %d, %d].\n",rank,i, particle->cell[0], particle->cell[1], particle->cell[2]);

    if(LocStatus == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input LocStatus pointer is NULL.\n");
    }
    // Initialize location status
    particle->location_status = (ParticleLocationStatus)LocStatus[i];
    LOG_ALLOW(LOCAL,LOG_VERBOSE, "[Rank %d]Particle [%d] Status set to: %d.\n",rank, i, particle->location_status);

    // The destination_rank is only set by the location search, not read from the swarm,
    // so we initialize it to a known invalid state.
    particle->destination_rank = MPI_PROC_NULL;
    
    // logging the completion of particle initialization
    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d]Completed initialization of Particle [%d]. \n", rank,i);
    

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateSwarmFields"
/**
 * @brief Internal helper implementation: `UpdateSwarmFields()`.
 * @details Local to this translation unit.
 */
PetscErrorCode UpdateSwarmFields(PetscInt i, const Particle *particle,
                                 PetscReal *positions, 
                                 PetscReal *velocities,
                                 PetscReal *weights, 
                                 PetscInt  *cellIndices, 
                                 PetscInt  *status,
                                 PetscReal *diffusivity,
                                 Cmpnts *diffusivitygradient,
                                 PetscReal *psi) 
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    if (!particle) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input Particle pointer is NULL.\n");
    }
    
    // --- 1. Position (x, y, z) ---
    if (positions) {
        positions[3 * i + 0] = particle->loc.x;
        positions[3 * i + 1] = particle->loc.y;
        positions[3 * i + 2] = particle->loc.z;
    }

    // --- 2. Velocity (u, v, w) ---
    if (velocities) {
        velocities[3 * i + 0] = particle->vel.x;
        velocities[3 * i + 1] = particle->vel.y;
        velocities[3 * i + 2] = particle->vel.z;
    }

    // --- 3. Weights (i, j, k) ---
    if (weights) {
        weights[3 * i + 0] = particle->weights.x;
        weights[3 * i + 1] = particle->weights.y;
        weights[3 * i + 2] = particle->weights.z;
    }
    
    // --- 4. Cell Indices (i, j, k) ---
    if (cellIndices) {
        cellIndices[3 * i + 0] = particle->cell[0];
        cellIndices[3 * i + 1] = particle->cell[1];
        cellIndices[3 * i + 2] = particle->cell[2];
    }

    // --- 5. Status ---
    if (status) {
        status[i] = (PetscInt)particle->location_status;
    }

    // --- 6. Diffusivity ---
    if (diffusivity) {
        diffusivity[i] = particle->diffusivity;
    }

    if(diffusivitygradient){
        diffusivitygradient[i].x = particle->diffusivitygradient.x;
        diffusivitygradient[i].y = particle->diffusivitygradient.y;
        diffusivitygradient[i].z = particle->diffusivitygradient.z;
    }
    // --- 7. Psi ---
    if (psi) {
        psi[i] = particle->psi;
    }

    // LOG_LOOP_ALLOW(LOCAL, LOG_VERBOSE, i, 1000, "Updated fields for Particle [%d].\n", i);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IsParticleInsideBoundingBox"
/**
 * @brief Internal helper implementation: `IsParticleInsideBoundingBox()`.
 * @details Local to this translation unit.
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
    LOG_ALLOW_SYNC(LOCAL, LOG_VERBOSE, "Particle PID %ld location: (%.6f, %.6f, %.6f).\n",particle->PID, loc.x, loc.y, loc.z);
    LOG_ALLOW_SYNC(LOCAL, LOG_VERBOSE, "BoundingBox min_coords: (%.6f, %.6f, %.6f), max_coords: (%.6f, %.6f, %.6f).\n", 
    min_coords.x, min_coords.y, min_coords.z, max_coords.x, max_coords.y, max_coords.z);

    // Check if the particle's location is within the bounding box
    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        // Particle is inside the bounding box
        LOG_ALLOW_SYNC(LOCAL,LOG_VERBOSE, "Particle PID %ld is inside the bounding box.\n",particle->PID);
        PROFILE_FUNCTION_END;
        return PETSC_TRUE;
    }

    // Particle is outside the bounding box
    LOG_ALLOW_SYNC(LOCAL, LOG_VERBOSE,"Particle PID %ld is outside the bounding box.\n",particle->PID);
    PROFILE_FUNCTION_END;
    return PETSC_FALSE;
}


#undef __FUNCT__
#define __FUNCT__ "UpdateParticleWeights"
/**
 * @brief Internal helper implementation: `UpdateParticleWeights()`.
 * @details Local to this translation unit.
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
    LOG_ALLOW(LOCAL, LOG_VERBOSE,
        "Calculating weights with distances: "
        "[LEFT=%f, RIGHT=%f, BOTTOM=%f, TOP=%f, FRONT=%f, BACK=%f].\n",
        d[LEFT], d[RIGHT], d[BOTTOM], d[TOP], d[FRONT], d[BACK]);

    // Compute and update the particle's weights
    particle->weights.x = d[LEFT] / (d[LEFT] + d[RIGHT]);
    particle->weights.y = d[BOTTOM] / (d[BOTTOM] + d[TOP]);
    particle->weights.z = d[BACK] / (d[FRONT] + d[BACK]);

    // LOG_ALLOW the updated weights
    LOG_ALLOW(LOCAL,LOG_VERBOSE,
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
/**
 * @brief Implementation of \ref InitializeParticleSwarm().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleSwarm.h`.
 * @see InitializeParticleSwarm()
 */

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
