    // ParticleSwarm.c

#include "ParticleSwarm.h"

#define INTERPOLATION_DISTANCE_TOLERANCE 1.0e-14

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


/**
 * @brief Internal helper implementation: `RegisterSwarmField()`.
 * @details Local to this translation unit.
 */
PetscErrorCode RegisterSwarmField(DM swarm, const char *fieldName, PetscInt fieldDim, PetscDataType dtype)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    ierr = DMSwarmRegisterPetscDatatypeField(swarm, fieldName, fieldDim, dtype); CHKERRQ(ierr);
    // PetscDataTypes is an extern char* [] defined in petscsystypes.h that gives string names for PetscDataType enums
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Registered field '%s' with dimension=%" PetscInt_FMT ", type=%s.\n",
            fieldName, fieldDim, PetscDataTypes[dtype]);
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


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
    PROFILE_FUNCTION_BEGIN;

    for (PetscInt raw_id = 0; raw_id < PARTICLE_FIELD_ID_COUNT; ++raw_id) {
        const ParticleFieldDescriptor *descriptor = NULL;

        ierr = ParticleFieldGetDescriptor((ParticleFieldId)raw_id, &descriptor); CHKERRQ(ierr);
        if (descriptor->registration == PARTICLE_FIELD_REGISTRATION_PETSC) continue;

        ierr = RegisterSwarmField(swarm, descriptor->canonical_name,
                                  descriptor->components, descriptor->data_type); CHKERRQ(ierr);
    }
    
    // Finalize the field registration after all fields have been added
    ierr = DMSwarmFinalizeFieldRegister(swarm); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO,"RegisterParticleFields - Finalized field registration.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief The owned cells a seeding mode draws from, with each cell's share of the measure.
 * @details Volumetric seeding draws from the rank's owned cells and surface seeding from the
 *          owned cells that touch the inlet face, so both are a box of cells from `origin`
 *          with `n` cells per axis (none when the rank owns no part of the target).
 */
typedef struct {
    PetscInt   face;       /**< -1 for the volume, otherwise the BCFace the points lie on */
    PetscInt   origin[3];  /**< Origin node of the box's first cell */
    PetscInt   n[3];       /**< Cells per axis */
    PetscReal *cumulative; /**< Summed measure of cells 0..c, i-fastest (NULL when not built) */
    PetscReal *bound;      /**< Largest measure density in each cell (NULL when not built) */
    PetscReal  measure;    /**< Total owned volume or area */
} SeedingCellTable;

/**
 * @brief Measure density at a logical point: the Jacobian determinant for the volume, or the
 *        area element |dx/ds x dx/dt| across the face for a face, with the normal coordinate
 *        on the face itself.
 */
static PetscErrorCode SeedingMeasureDensity(UserCtx *user, const Cmpnts ***coor, PetscInt face,
                                            PetscInt i, PetscInt j, PetscInt k, const PetscReal logic[3],
                                            PetscReal *density)
{
    PetscReal J[3][3], detJ, l[3] = {logic[0], logic[1], logic[2]};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (face < 0) {
        PetscCall(MetricJacobian(user, coor, i, j, k, l[0], l[1], l[2], J, &detJ));
        *density = detJ;
    } else {
        const PetscInt a = face / 2, b = (a + 1) % 3, c = (a + 2) % 3;
        Cmpnts         cross;
        l[a] = (face % 2) ? 1.0 : 0.0;
        PetscCall(MetricJacobian(user, coor, i, j, k, l[0], l[1], l[2], J, NULL));
        cross.x = J[1][b] * J[2][c] - J[2][b] * J[1][c];
        cross.y = J[2][b] * J[0][c] - J[0][b] * J[2][c];
        cross.z = J[0][b] * J[1][c] - J[1][b] * J[0][c];
        *density = PetscSqrtReal(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Tabulate the measure of the owned cells a seeding mode draws from.
 * @details For `face` < 0 the cells are the rank's owned cells and the measure is volume; for
 *          a face they are the owned cells touching it, and the measure is their area on it,
 *          none when the rank does not hold part of the face. `cumulative[c]` sums cells
 *          `0..c` in i-fastest order, so a uniform draw on `[0, measure)` picks a cell in
 *          proportion to its measure. The measure integrates the density with 2x2x2 Gauss
 *          points, which is exact for the volume (the determinant has degree 2 in each
 *          coordinate) and for the area of a planar face; on a face the points coincide in
 *          pairs, which gives the 2x2 rule. `bound[c]` is the largest density on a 3x3x3
 *          lattice of the cell, the envelope used to sample a uniform point inside it.
 * @param[in]  user        Grid context.
 * @param[in]  coor        Ghosted node-coordinate array.
 * @param[in]  face        -1 for the volume, otherwise the BCFace to seed on.
 * @param[in]  build       PETSC_TRUE to allocate `cumulative` and `bound` (caller frees with
 *                         PetscFree2), PETSC_FALSE when only the total is wanted.
 * @param[out] table       The owned cells and their measure.
 * @return PetscErrorCode 0 on success; an error for a cell whose measure is not positive.
 */
static PetscErrorCode BuildSeedingCellTable(UserCtx *user, const Cmpnts ***coor, PetscInt face,
                                            PetscBool build, SeedingCellTable *table)
{
    const PetscReal gauss[2] = {0.5 - 0.5 / PetscSqrtReal(3.0), 0.5 + 0.5 / PetscSqrtReal(3.0)};
    PetscReal       sum = 0.0;
    PetscInt        start, n_cells, c = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PetscMemzero(table, sizeof(*table)));
    table->face = face;
    table->origin[0] = user->info.xs; table->origin[1] = user->info.ys; table->origin[2] = user->info.zs;
    PetscCall(GetOwnedCellRange(&user->info, 0, &start, &table->n[0]));
    PetscCall(GetOwnedCellRange(&user->info, 1, &start, &table->n[1]));
    PetscCall(GetOwnedCellRange(&user->info, 2, &start, &table->n[2]));
    if (face >= 0) {
        const PetscInt a = face / 2;
        const PetscInt nodes[3] = {user->info.mx, user->info.my, user->info.mz};
        PetscBool      holds_face = PETSC_FALSE;
        PetscCall(CanRankServiceInletFace(user, &user->info, user->info.mx - 1, user->info.my - 1, user->info.mz - 1, &holds_face));
        if (holds_face) {
            table->n[a] = 1;
            if (face % 2) table->origin[a] = nodes[a] - 3; /* origin node of the last cell */
        } else {
            table->n[0] = table->n[1] = table->n[2] = 0;
        }
    }
    n_cells = table->n[0] * table->n[1] * table->n[2];
    if (build) PetscCall(PetscMalloc2(n_cells, &table->cumulative, n_cells, &table->bound));
    for (PetscInt lk = 0; lk < table->n[2]; ++lk) {
        for (PetscInt lj = 0; lj < table->n[1]; ++lj) {
            for (PetscInt li = 0; li < table->n[0]; ++li, ++c) {
                const PetscInt i = table->origin[0] + li, j = table->origin[1] + lj, k = table->origin[2] + lk;
                PetscReal      cell_measure = 0.0, cell_bound = 0.0, density;
                for (PetscInt g = 0; g < 8; ++g) {
                    const PetscReal l[3] = {gauss[g & 1], gauss[(g >> 1) & 1], gauss[(g >> 2) & 1]};
                    PetscCall(SeedingMeasureDensity(user, coor, face, i, j, k, l, &density));
                    cell_measure += 0.125 * density;
                }
                PetscCheck(cell_measure > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                           "Cell with origin node (%d,%d,%d) has %s %g; particle seeding needs positive %s.",
                           (int)i, (int)j, (int)k, face < 0 ? "volume" : "inlet area", (double)cell_measure,
                           face < 0 ? "volumes" : "areas");
                if (build) {
                    for (PetscInt q = 0; q < 27; ++q) {
                        const PetscReal l[3] = {0.5 * (q % 3), 0.5 * ((q / 3) % 3), 0.5 * (q / 9)};
                        PetscCall(SeedingMeasureDensity(user, coor, face, i, j, k, l, &density));
                        cell_bound = PetscMax(cell_bound, density);
                    }
                    table->bound[c] = cell_bound;
                }
                sum += cell_measure;
                if (build) table->cumulative[c] = sum;
            }
        }
    }
    table->measure = sum;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Draw one seeding point uniform in the measure of a seeding table.
 * @details The cell is drawn in proportion to its measure from `table->cumulative`
 *          (see BuildSeedingCellTable); the point inside it is drawn uniform in that measure by
 *          accepting uniform logical coordinates with probability density / `bound`. A surface
 *          point is then put just inside the cell, 1e-6 of the cell from the face.
 * @param[in]  user   Grid context.
 * @param[in]  coor   Ghosted node-coordinate array.
 * @param[in]  table  A table built with `build` set; it must hold at least one cell.
 * @param[in]  rand_logic_i_ptr,rand_logic_j_ptr,rand_logic_k_ptr Logical-space RNGs.
 * @param[out] ci_metric_lnode_out,cj_metric_lnode_out,ck_metric_lnode_out Origin node of the cell.
 * @param[out] xi_metric_logic_out,eta_metric_logic_out,zta_metric_logic_out Logical point in it.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode DrawSeedingPoint(
    UserCtx *user, const Cmpnts ***coor, const SeedingCellTable *table,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out)
{
    const PetscInt n_cells = table->n[0] * table->n[1] * table->n[2];
    PetscInt       lo = 0, hi = n_cells - 1;
    PetscReal      r_val, logic[3];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // --- 1. Select a Cell in Proportion to Its Measure ---
    PetscCall(PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val));
    r_val *= table->cumulative[n_cells - 1];
    while (lo < hi) { // First cell whose cumulative measure exceeds the draw
        const PetscInt mid = lo + (hi - lo) / 2;
        if (table->cumulative[mid] > r_val) hi = mid; else lo = mid + 1;
    }
    *ci_metric_lnode_out = table->origin[0] + lo % table->n[0];
    *cj_metric_lnode_out = table->origin[1] + (lo / table->n[0]) % table->n[1];
    *ck_metric_lnode_out = table->origin[2] + lo / (table->n[0] * table->n[1]);

    // --- 2. Draw Intra-Cell Logical Coordinates Uniform in the Measure ---
    // Accept uniform logicals with probability density / bound; a cell whose density varies
    // little accepts almost every draw. The cap only guards a degenerate envelope.
    for (PetscInt attempt = 0; attempt < 1000; ++attempt) {
        PetscReal density, u;
        PetscCall(PetscRandomGetValueReal(*rand_logic_i_ptr, &logic[0]));
        PetscCall(PetscRandomGetValueReal(*rand_logic_j_ptr, &logic[1]));
        PetscCall(PetscRandomGetValueReal(*rand_logic_k_ptr, &logic[2]));
        // Ensure logical coordinates are strictly within [0,1) for robustness with MetricLogicalToPhysical
        for (PetscInt d = 0; d < 3; ++d) logic[d] = PetscMax(PetscMin(logic[d], 1.0 - 1.0e-7), 0.0);
        PetscCall(SeedingMeasureDensity(user, coor, table->face, *ci_metric_lnode_out, *cj_metric_lnode_out,
                                        *ck_metric_lnode_out, logic, &density));
        PetscCall(PetscRandomGetValueReal(*rand_logic_j_ptr, &u));
        if (u * table->bound[lo] <= density) break;
        if (attempt == 999) {
            LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: 1000 rejected draws in cell (%" PetscInt_FMT ",%" PetscInt_FMT ",%" PetscInt_FMT "); accepting the last one.\n",
                      user->simCtx->rank, *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out);
        }
    }
    if (table->face >= 0) logic[table->face / 2] = (table->face % 2) ? 1.0 - 1.0e-6 : 1.0e-6;
    *xi_metric_logic_out = logic[0]; *eta_metric_logic_out = logic[1]; *zta_metric_logic_out = logic[2];

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Initialize position-independent particle fields after a particle is created.
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
    // Owned cells and their measure, for volumetric and surface seeding
    SeedingCellTable seed_table;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    SimCtx *simCtx = user->simCtx;

    // --- 1. Input Validation and Basic Setup ---
    if (!user || !rand_logic_i || !rand_logic_j || !rand_logic_k) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null user or RNG pointer.");
    }
    rank = simCtx->rank;
    size = simCtx->size;

    // Get DMDA information for the node-centered coordinate grid (user->da)
    Coor_local = user->lCoor;
    if (!Coor_local) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "user->lCoor is NULL: grid coordinates are not set up.");
    ierr = DMDAVecGetArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
    ierr = PetscMemzero(&seed_table, sizeof(seed_table)); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Initializing %" PetscInt_FMT " particles. Mode: %s.\n",
            rank, particlesPerProcess, ParticleInitializationToString(simCtx->ParticleInitialization));

    // --- 2. Seeding Table for Volumetric and Surface Initialization ---
    // Volume draws from the owned cells, Surface from the owned cells on the inlet face; the
    // particle count was split across ranks by the same measure (DistributeParticlesByOwnedMeasure).
    const PetscBool seeds_from_table = (PetscBool)(simCtx->ParticleInitialization == PARTICLE_INIT_VOLUME ||
                                                   simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM);
    if (seeds_from_table && particlesPerProcess > 0) {
        const PetscInt face = (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM) ? (PetscInt)user->identifiedInletBCFace : -1;
        ierr = BuildSeedingCellTable(user, coor_nodes_local_array, face, PETSC_TRUE, &seed_table); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Seeding %" PetscInt_FMT " particles uniformly over %" PetscInt_FMT " owned cells of total %s %.6e.\n",
                  rank, particlesPerProcess, seed_table.n[0] * seed_table.n[1] * seed_table.n[2],
                  face < 0 ? "volume" : "inlet area", (double)seed_table.measure);
    }

    // --- 3. Get Access to Swarm Fields ---
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION),       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID),    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_LOCATION_STATUS),NULL,NULL,(void**)&status_field); CHKERRQ(ierr);

    // --- 4. Determine Starting Global PID for this Rank ---
    // PIDs are consecutive across ranks in rank order, so each rank starts where the ranks
    // before it end. Taking that from the actual local counts keeps it right however the
    // particles were split (evenly, or by owned volume or inlet area for table seeding).
    PetscInt base_pid_for_rank = 0;
    {
        PetscInt local_count = particlesPerProcess;
        ierr = MPI_Exscan(&local_count, &base_pid_for_rank, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRMPI(ierr);
        if (rank == 0) base_pid_for_rank = 0; /* MPI_Exscan leaves rank 0's result undefined */
    }

    // --- 4b. Locate the Point Source Once ---
    // Every rank places its share at the same point, so its cell is found by one search per rank;
    // the lowest rank that holds it shares the cell, and settlement then routes the particles there.
    PetscInt point_cell[3] = {-1, -1, -1};
    if (simCtx->ParticleInitialization == PARTICLE_INIT_POINT_SOURCE) {
        Particle               probe;
        ParticleLocationStatus probe_status = NEEDS_LOCATION;
        PetscMPIInt            candidate, owner = size;
        ierr = PetscMemzero(&probe, sizeof(probe)); CHKERRQ(ierr);
        probe.PID = -1;
        probe.cell[0] = -1; probe.cell[1] = -1; probe.cell[2] = -1;
        probe.loc.x = simCtx->psrc_x; probe.loc.y = simCtx->psrc_y; probe.loc.z = simCtx->psrc_z;
        ierr = LocateParticleOrFindMigrationTarget(user, NULL, &probe, &probe_status); CHKERRQ(ierr);
        candidate = (probe_status == ACTIVE_AND_LOCATED) ? rank : size;
        ierr = MPI_Allreduce(&candidate, &owner, 1, MPI_INT, MPI_MIN, PETSC_COMM_WORLD); CHKERRMPI(ierr);
        if (owner < size) {
            if (rank == owner) { point_cell[0] = probe.cell[0]; point_cell[1] = probe.cell[1]; point_cell[2] = probe.cell[2]; }
            ierr = MPI_Bcast(point_cell, 3, MPIU_INT, owner, PETSC_COMM_WORLD); CHKERRMPI(ierr);
            LOG_ALLOW(GLOBAL, LOG_INFO, "Point source (%.6e, %.6e, %.6e) lies in cell (%" PetscInt_FMT ",%" PetscInt_FMT ",%" PetscInt_FMT ") on rank %d.\n",
                      simCtx->psrc_x, simCtx->psrc_y, simCtx->psrc_z, point_cell[0], point_cell[1], point_cell[2], owner);
        } else {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Point source (%.6e, %.6e, %.6e) was not located on any rank; its particles will be searched for.\n",
                      simCtx->psrc_x, simCtx->psrc_y, simCtx->psrc_z);
        }
    }

    // --- 5. Loop Over Particles to Initialize ---
    PetscInt n_unplaced = 0; // Table-seeded particles this rank could not place
    for (PetscInt p = 0; p < particlesPerProcess; p++) {
        PetscInt  ci_metric_lnode, cj_metric_lnode, ck_metric_lnode; 
        PetscReal xi_metric_logic, eta_metric_logic, zta_metric_logic; 
        Cmpnts    phys_coords = {0.0, 0.0, 0.0}; 
        PetscBool particle_placed_by_this_rank = PETSC_FALSE; 

        if (seeds_from_table) { // --- 5.a. Volumetric or Surface Initialization ---
            if (seed_table.n[0] * seed_table.n[1] * seed_table.n[2] > 0) {
                ierr = DrawSeedingPoint(user, coor_nodes_local_array, &seed_table,
                                        rand_logic_i, rand_logic_j, rand_logic_k,
                                        &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                        &xi_metric_logic, &eta_metric_logic, &zta_metric_logic); CHKERRQ(ierr);
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                        ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                        xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                        &phys_coords); CHKERRQ(ierr);
                particle_placed_by_this_rank = PETSC_TRUE;
            } else {
                n_unplaced++;
            }
        }else if(simCtx->ParticleInitialization == PARTICLE_INIT_POINT_SOURCE){ // --- 5.c. Point Source Initialization ---
            // All particles placed at the user-specified fixed point (psrc_x, psrc_y, psrc_z).
            // No random number generation or logical-to-physical conversion needed.
            phys_coords.x = simCtx->psrc_x;
            phys_coords.y = simCtx->psrc_y;
            phys_coords.z = simCtx->psrc_z;
            ci_metric_lnode = point_cell[0]; cj_metric_lnode = point_cell[1]; ck_metric_lnode = point_cell[2];
            particle_placed_by_this_rank = (PetscBool)(point_cell[0] >= 0);
        }else {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Unknown ParticleInitialization mode %d.", simCtx->ParticleInitialization);
        }

        // --- 5.c. Store Particle Properties ---
        positions_field[3*p+0] = phys_coords.x; 
        positions_field[3*p+1] = phys_coords.y; 
        positions_field[3*p+2] = phys_coords.z;

        particleIDs[p]         = (PetscInt64)base_pid_for_rank + p; 
        // Start the first settlement's walk in the cell the particle was placed in; a particle
        // left unplaced has no cell and is located by search.
        if (particle_placed_by_this_rank) {
            cellIDs_petsc[3*p+0] = ci_metric_lnode; cellIDs_petsc[3*p+1] = cj_metric_lnode; cellIDs_petsc[3*p+2] = ck_metric_lnode;
        } else {
            cellIDs_petsc[3*p+0] = -1; cellIDs_petsc[3*p+1] = -1; cellIDs_petsc[3*p+2] = -1;
        }
        status_field[p]        = UNINITIALIZED;

        // --- 5.d. The particle's seed line ---
        if (particle_placed_by_this_rank) {
            LOG_PARTICLE(simCtx, particleIDs[p], "seed",
                         "%s placement on rank %d at (%.6e, %.6e, %.6e): cell (%" PetscInt_FMT ",%" PetscInt_FMT ",%" PetscInt_FMT ") logical (%.4f, %.4f, %.4f)\n",
                         ParticleInitializationToString(simCtx->ParticleInitialization), (int)rank,
                         phys_coords.x, phys_coords.y, phys_coords.z, ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                         xi_metric_logic, eta_metric_logic, zta_metric_logic);
        } else {
            // A point source no rank located, or a table-seeded particle on a rank with no cells.
            LOG_PARTICLE(simCtx, particleIDs[p], "seed",
                         "%s: not placed by rank %d; parked at (%.6e, %.6e, %.6e) for the settlement to locate\n",
                         ParticleInitializationToString(simCtx->ParticleInitialization), (int)rank,
                         phys_coords.x, phys_coords.y, phys_coords.z);
        }
    }
    if (n_unplaced > 0) {
        LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: %" PetscInt_FMT " %s-seeded particle(s) could not be placed in this rank's cells and are left for the settlement to locate.\n",
                  rank, n_unplaced, ParticleInitializationToString(simCtx->ParticleInitialization));
    }

    // --- 6. Restore Pointers and Cleanup ---
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION),       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID),    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_LOCATION_STATUS), NULL, NULL, (void**)&status_field); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
    ierr = PetscFree2(seed_table.cumulative, seed_table.bound); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Completed processing for %" PetscInt_FMT " particles.\n",
            rank, particlesPerProcess);

    PROFILE_FUNCTION_END;
            
    PetscFunctionReturn(0);
}

/**
 * @brief Assign one configured initial value to a swarm field entry.
 */
static PetscErrorCode InitializeSwarmFieldValue(const ParticleFieldDescriptor *descriptor,
                                                PetscInt p, PetscReal *fieldData)
{
    PetscFunctionBeginUser;

    PROFILE_HELPER_BEGIN;

    PetscCheck(descriptor != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Particle field descriptor cannot be NULL during initialization.");
    PetscCheck(fieldData != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Particle field data cannot be NULL during initialization.");
    PetscCheck(descriptor->data_type == PETSC_REAL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Default initialization currently supports only PETSC_REAL particle fields; '%s' uses %s.",
               descriptor->canonical_name, PetscDataTypes[descriptor->data_type]);
    PetscCheck((descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_DEFAULT_INITIALIZE) != 0,
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Particle field '%s' does not use catalog-default initialization.",
               descriptor->canonical_name);

    for (PetscInt d = 0; d < descriptor->components; ++d) {
        fieldData[descriptor->components * p + d] = descriptor->default_real_value;
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Apply a configured initial field specification across the local swarm.
 */
static PetscErrorCode AssignInitialFieldToSwarm(UserCtx *user, ParticleFieldId field_id)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscReal     *fieldData = NULL;
    PetscInt       nLocal;
    const ParticleFieldDescriptor *descriptor = NULL;
    const char     *fieldName = NULL;

    PetscFunctionBeginUser;
    
    PROFILE_FUNCTION_BEGIN;

    ierr = ParticleFieldGetDescriptor(field_id, &descriptor); CHKERRQ(ierr);
    fieldName = descriptor->canonical_name;

    // Get the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "%" PetscInt_FMT " local particles found.\n", nLocal);

    // Retrieve the swarm field pointer for the specified fieldName
    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&fieldData); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Retrieved field '%s'.\n", fieldName);

    // Loop over all particles and update the field using the helper function
    for (PetscInt p = 0; p < nLocal; p++) {
        ierr = InitializeSwarmFieldValue(descriptor, p, fieldData); CHKERRQ(ierr);
    PetscReal disp_data[descriptor->components];
    
        for (PetscInt d = 0; d < descriptor->components; d++) {
    disp_data[d] = fieldData[descriptor->components * p + d];
        }
        LOG_LOOP_ALLOW(LOCAL,LOG_VERBOSE,p, 100," Particle %" PetscInt_FMT ": %s[%" PetscInt_FMT "] = [%.6f, ...,%.6f].\n", p,fieldName,descriptor->components,disp_data[0],disp_data[descriptor->components-1]);
    }
    
    // Restore the swarm field pointer
    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&fieldData); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "Initialization of field '%s' complete.\n", fieldName);
    

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


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

    LOG_ALLOW(GLOBAL, LOG_INFO, "Initializing swarm with %" PetscInt_FMT " particles per process. Mode: %s.\n",
            particlesPerProcess, ParticleInitializationToString(simCtx->ParticleInitialization));

    // --- 1. Parse BCS File for Inlet Information (if surface initialization) ---
    if (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM) {
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
    ierr = AssignInitialFieldToSwarm(user, PARTICLE_FIELD_ID_VELOCITY); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "'velocity' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'weight' field.\n");
    ierr = AssignInitialFieldToSwarm(user, PARTICLE_FIELD_ID_WEIGHT); CHKERRQ(ierr); // Weight is a three-component field.
    LOG_ALLOW(LOCAL, LOG_INFO, "'weight' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'Diffusivity' field.\n");
    ierr = AssignInitialFieldToSwarm(user, PARTICLE_FIELD_ID_DIFFUSIVITY); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'Diffusivity' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'DiffusivityGradient' field.\n");
    ierr = AssignInitialFieldToSwarm(user, PARTICLE_FIELD_ID_DIFFUSIVITY_GRADIENT); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'DiffusivityGradient' field initialization complete.\n");

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Initializing 'Psi' (Scalar) field.\n");
    ierr = AssignInitialFieldToSwarm(user, PARTICLE_FIELD_ID_PSI); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "'P' field initialization complete.\n");
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully completed all swarm property initialization.\n");


    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


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
        LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"Rank %d receives an extra particle. Total: %" PetscInt_FMT "\n", rank, *particlesPerProcess);
    } else {
        LOG_ALLOW_SYNC(GLOBAL,LOG_INFO, "Rank %d receives %" PetscInt_FMT " particles.\n", rank, *particlesPerProcess);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


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

/**
 * @brief Split particles across ranks in proportion to the volume or inlet area each rank owns.
 * @details Volumetric and surface seeding place each rank's particles uniformly in its own
 *          volume or its own part of the inlet, so the split must follow that measure for the
 *          density to be uniform across rank boundaries; a rank without part of the inlet seeds
 *          none. Largest-remainder rounding keeps the total exact, and every rank computes the
 *          same answer from the gathered measures.
 * @param[in]  user          Finest-level context whose DMDA ownership defines the split.
 * @param[in]  face          -1 for the volume, otherwise the inlet BCFace.
 * @param[in]  numParticles  Global particle count.
 * @param[out] localCount    Particles this rank seeds.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode DistributeParticlesByOwnedMeasure(UserCtx *user, PetscInt face, PetscInt numParticles, PetscInt *localCount)
{
    PetscMPIInt     rank, size;
    PetscInt        assigned = 0;
    PetscInt       *count = NULL;
    PetscReal       owned, total = 0.0;
    PetscReal      *volumes = NULL, *frac = NULL;
    SeedingCellTable table;
    Vec             coor_local;
    const Cmpnts ***coor = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    rank = user->simCtx->rank;
    size = user->simCtx->size;
    coor_local = user->lCoor;
    PetscCall(DMDAVecGetArrayRead(user->fda, coor_local, (void *)&coor));
    PetscCall(BuildSeedingCellTable(user, coor, face, PETSC_FALSE, &table));
    owned = table.measure;
    PetscCall(DMDAVecRestoreArrayRead(user->fda, coor_local, (void *)&coor));
    PetscCall(PetscMalloc3(size, &volumes, size, &count, size, &frac));
    PetscCallMPI(MPI_Allgather(&owned, 1, MPIU_REAL, volumes, 1, MPIU_REAL, PETSC_COMM_WORLD));
    for (PetscMPIInt r = 0; r < size; ++r) total += volumes[r];
    PetscCheck(total > 0.0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               face < 0 ? "Volumetric particle seeding found no owned volume on any rank."
                        : "Surface particle seeding found no inlet area on any rank; it needs an INLET face.");
    for (PetscMPIInt r = 0; r < size; ++r) {
        const PetscReal share = (PetscReal)numParticles * volumes[r] / total;
        count[r] = (PetscInt)PetscFloorReal(share);
        frac[r]  = share - (PetscReal)count[r];
        assigned += count[r];
    }
    for (PetscInt extra = numParticles - assigned; extra > 0; --extra) {
        PetscMPIInt best = 0;
        for (PetscMPIInt r = 1; r < size; ++r) if (frac[r] > frac[best]) best = r;
        count[best] += 1;
        frac[best] = -1.0;
    }
    *localCount = count[rank];
    PetscCall(PetscFree3(volumes, count, frac));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

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
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Number of particles must be positive. Given: %" PetscInt_FMT "\n", numParticles);
        PROFILE_FUNCTION_END;
        return PETSC_ERR_ARG_OUTOFRANGE;
    }
    
    // Retrieve MPI rank and size
    rank = user->simCtx->rank;
    size = user->simCtx->size;
    LOG_ALLOW(GLOBAL,LOG_INFO," Domain dimensions: [%.2f,%.2f],[%.2f,%.2f],[%.2f,%.2f] \n", 
        user->Min_X,user->Max_X,user->Min_Y,user->Max_Y, user->Min_Z,user->Max_Z);
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "[Rank %d] Local Bounding Box: [%.2f,%.2f],[%.2f,%.2f],[%.2f,%.2f] \n", 
        rank,user->bbox.min_coords.x,user->bbox.max_coords.x,
        user->bbox.min_coords.y,user->bbox.max_coords.y,
        user->bbox.min_coords.z,user->bbox.max_coords.z);
    // Distribute particles among MPI processes
    if (user->simCtx->ParticleInitialization == PARTICLE_INIT_VOLUME) {
        ierr = DistributeParticlesByOwnedMeasure(user, -1, numParticles, particlesPerProcess); CHKERRQ(ierr);
    } else if (user->simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM) {
        ierr = DistributeParticlesByOwnedMeasure(user, (PetscInt)user->identifiedInletBCFace, numParticles, particlesPerProcess); CHKERRQ(ierr);
    } else {
        ierr = DistributeParticles(numParticles, rank, size, particlesPerProcess, &remainder); CHKERRQ(ierr);
    }

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

    // Set the local number of particles for this rank and spare capacity for arrivals (10% of the
    // local count; PETSc reallocates beyond it). Every resize zeroes the whole spare region, so a
    // buffer sized to the global count would cost each rank memory and work of order numParticles.
    ierr = DMSwarmSetLocalSizes(user->swarm, *particlesPerProcess, *particlesPerProcess / 10); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO, "Set local swarm size: %" PetscInt_FMT " particles.\n", *particlesPerProcess);

    // Optionally, LOG_ALLOW detailed DM info in debug mode
    if (get_log_level() >= LOG_DEBUG && is_function_allowed(__func__)) {
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

    PROFILE_HELPER_BEGIN;

    // The particle's values are reported by the stage that uses them, not here: this runs
    // several times per particle per step.
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output Particle pointer is NULL. \n");
    }
    if(PIDs == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input PIDs pointer is NULL.\n");
    }
    if(positions == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input positions pointer is NULL.\n");
    }
    if(cellIndices == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input cellIndices pointer is NULL.\n");
    }
    if(LocStatus == NULL){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input LocStatus pointer is NULL.\n");
    }

    particle->PID = PIDs[i];

    // A NULL optional field is one the caller does not need; it takes the default shown.
    if(weights == NULL){
        particle->weights.x = 1.0;
        particle->weights.y = 1.0;
        particle->weights.z = 1.0;
    }else{
        particle->weights.x = weights[3 * i];
        particle->weights.y = weights[3 * i + 1];
        particle->weights.z = weights[3 * i + 2];
    }

    particle->loc.x = positions[3 * i];
    particle->loc.y = positions[3 * i + 1];
    particle->loc.z = positions[3 * i + 2];

    if(velocities == NULL){
        particle->vel.x = 0.0;
        particle->vel.y = 0.0;
        particle->vel.z = 0.0;
    }else{
        particle->vel.x = velocities[3 * i];
        particle->vel.y = velocities[3 * i + 1];
        particle->vel.z = velocities[3 * i + 2];
    }

    particle->diffusivity = (diffusivity == NULL) ? 1.0 : diffusivity[i];

    if(diffusivitygradient == NULL){
        particle->diffusivitygradient.x = 0.0;
        particle->diffusivitygradient.y = 0.0;
        particle->diffusivitygradient.z = 0.0;
    }else{
        particle->diffusivitygradient = diffusivitygradient[i];
    }

    particle->psi = (psi == NULL) ? 0.0 : psi[i];

    particle->cell[0] = cellIndices[3 * i];
    particle->cell[1] = cellIndices[3 * i + 1];
    particle->cell[2] = cellIndices[3 * i + 2];

    particle->location_status = (ParticleLocationStatus)LocStatus[i];

    // The destination_rank is only set by the location search, not read from the swarm,
    // so we initialize it to a known invalid state.
    particle->destination_rank = MPI_PROC_NULL;

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

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
    PROFILE_HELPER_BEGIN;
    
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

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `IsParticleInsideBoundingBox()`.
 * @details Local to this translation unit.
 */
PetscBool IsParticleInsideBoundingBox(const BoundingBox *bbox, const Particle *particle)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

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

    // Check if the particle's location is within the bounding box
    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        PROFILE_FUNCTION_END;
        return PETSC_TRUE;
    }

    PROFILE_FUNCTION_END;
    return PETSC_FALSE;
}


/**
 * @brief Internal helper implementation: `UpdateParticleWeights()`.
 * @details Local to this translation unit.
 */
PetscErrorCode UpdateParticleWeights(PetscReal *d, Particle *particle) {

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    // Validate input pointers
    if (!d || !particle) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                "Null pointer argument (d or particle).");
    }


    // Validate distances
    for (PetscInt i = LEFT; i < NUM_FACES; i++) {
        if (d[i] <= INTERPOLATION_DISTANCE_TOLERANCE) {
            // A particle on a face: benign, so reported only to a VERBOSE run of this function.
            LOG_ALLOW(LOCAL, LOG_VERBOSE,
                "face distance d[%" PetscInt_FMT "] = %g <= %g; "
                "clamping to avoid zero/negative.\n",
                i, (double)d[i], INTERPOLATION_DISTANCE_TOLERANCE);
            d[i] = INTERPOLATION_DISTANCE_TOLERANCE;
        }
    }

    // Compute and update the particle's weights
    particle->weights.x = d[LEFT] / (d[LEFT] + d[RIGHT]);
    particle->weights.y = d[BOTTOM] / (d[BOTTOM] + d[TOP]);
    particle->weights.z = d[BACK] / (d[FRONT] + d[BACK]);



    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting particle swarm setup for %" PetscInt_FMT " particles.\n", simCtx->np);

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
        ierr = InitializeLogicalSpaceRNGs(user->simCtx->particleRandomSeed, &rand_logic_i, &rand_logic_j, &rand_logic_k); CHKERRQ(ierr);
        ierr = AssignInitialPropertiesToSwarm(user, particlesPerProcess, &randx, &randy, &randz, &rand_logic_i, &rand_logic_j, &rand_logic_k, simCtx->bboxlist); CHKERRQ(ierr);
        ierr = FinalizeSwarmSetup(&randx, &randy, &randz, &rand_logic_i, &rand_logic_j, &rand_logic_k); CHKERRQ(ierr);

    } else {
        // --- PATH B: Load particle population from restart files ---
        // This path is only taken if simCtx->StartStep > 0 AND simCtx->particleRestartMode == "load"
        LOG_ALLOW(GLOBAL, LOG_INFO, "Mode: LOAD. Loading particle population from files for step %" PetscInt_FMT ".\n", simCtx->StartStep);

        ierr = PreCheckAndResizeSwarm(user, simCtx->StartStep, "dat"); CHKERRQ(ierr);

        ierr = ReadAllSwarmFields(user, simCtx->StartStep); CHKERRQ(ierr);
        // Note: We check for file-open errors inside ReadAllSwarmFields now.

        LOG_ALLOW(GLOBAL, LOG_INFO, "Particle data loaded. CellID and status are preserved from file.\n");
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle swarm setup complete.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
