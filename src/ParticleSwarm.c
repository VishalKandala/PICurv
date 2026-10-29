    // ParticleSwarm.c

#include "ParticleSwarm.h"
#include "ParticleInitialConditions.h"

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
 * @brief The owned cells a seeding mode places particles in, and which particles each gets.
 * @details Volumetric seeding uses the rank's owned cells and surface seeding the owned cells
 *          that touch the inlet face, so both are a box of cells from `origin` with `n` cells
 *          per axis (none when the rank owns no part of the target). Cell `c`, in i-fastest
 *          order, seeds the particles with IDs `pid_lo[c]` to `pid_hi[c] - 1`.
 */
typedef struct {
    PetscInt    face;      /**< -1 for the volume, otherwise the BCFace the points lie on */
    PetscInt    origin[3]; /**< Origin node of the box's first cell */
    PetscInt    n[3];      /**< Cells per axis */
    PetscInt64 *pid_lo;    /**< First particle ID of each cell */
    PetscInt64 *pid_hi;    /**< One past the last particle ID of each cell */
    PetscReal  *bound;     /**< Largest measure density in each cell (NULL when not built) */
    PetscInt    count;     /**< Particles this rank seeds */
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
 * @brief Decide which particles each owned cell seeds, the same way at every rank layout.
 * @details Collective. The cells are those volumetric (`face` < 0) or surface seeding uses,
 *          and each cell's measure is its volume or its area on the face, integrated with 2x2x2
 *          Gauss points: exact for the volume (the determinant has degree 2 in each
 *          coordinate) and for the area of a planar face, where the points coincide in pairs
 *          and give the 2x2 rule.
 *
 *          Particles are stratified over the measure in global cell order (i fastest): with
 *          P(c) the measure of all cells before cell c, T the total and N the particle count,
 *          cell c seeds the IDs floor(N P(c) / T) to floor(N (P(c) + m_c) / T) - 1, so each
 *          cell holds its expected share to within one, the IDs are 0..N-1, and both depend
 *          only on the grid. For the result to be identical at any rank layout the measures
 *          are summed as integers, each cell's measure quantised against the largest cell
 *          (an exact maximum) to 2^Q units, with Q chosen so the total fits 63 bits. P(c) is
 *          the sum of the (j,k) rows before the cell's row, from one reduction over the rows,
 *          plus the cells before it in its row, from a scan over the ranks sharing the row.
 * @param[in]  user          Grid context.
 * @param[in]  coor          Ghosted node-coordinate array.
 * @param[in]  face          -1 for the volume, otherwise the BCFace to seed on.
 * @param[in]  numParticles  Global particle count.
 * @param[in]  build_bound   PETSC_TRUE to also tabulate each cell's largest measure density on
 *                           a 3x3x3 lattice, the envelope DrawSeedingPoint samples against.
 * @param[out] table         The owned cells and their particles; free with FreeSeedingCellTable.
 * @return PetscErrorCode 0 on success; an error for a cell whose measure is not positive.
 */
static PetscErrorCode BuildSeedingCellTable(UserCtx *user, const Cmpnts ***coor, PetscInt face,
                                            PetscInt numParticles, PetscBool build_bound,
                                            SeedingCellTable *table)
{
    const PetscReal gauss[2] = {0.5 - 0.5 / PetscSqrtReal(3.0), 0.5 + 0.5 / PetscSqrtReal(3.0)};
    const PetscInt  nodes[3] = {user->info.mx, user->info.my, user->info.mz};
    PetscInt        start, n_cells, n_rows_local, global_cells[3], global_origin[3] = {0, 0, 0}, quantum_bits = 62;
    PetscInt64      n_global_cells;
    PetscReal       local_max = 0.0, global_max = 0.0;
    PetscReal      *measure = NULL;
    uint64_t       *rows = NULL, *row_before = NULL, *segment = NULL, *segment_before = NULL, *quanta = NULL;
    uint64_t        total = 0;
    MPI_Comm        row_comm = MPI_COMM_NULL;
    PetscMPIInt     row_rank = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PetscMemzero(table, sizeof(*table)));
    table->face = face;
    table->origin[0] = user->info.xs; table->origin[1] = user->info.ys; table->origin[2] = user->info.zs;
    PetscCall(GetOwnedCellRange(&user->info, 0, &start, &table->n[0]));
    PetscCall(GetOwnedCellRange(&user->info, 1, &start, &table->n[1]));
    PetscCall(GetOwnedCellRange(&user->info, 2, &start, &table->n[2]));
    for (PetscInt d = 0; d < 3; ++d) global_cells[d] = nodes[d] - 2;
    if (face >= 0) {
        const PetscInt a = face / 2;
        PetscBool      holds_face = PETSC_FALSE;
        PetscCall(CanRankServiceInletFace(user, &user->info, nodes[0] - 1, nodes[1] - 1, nodes[2] - 1, &holds_face));
        global_cells[a] = 1;
        if (face % 2) global_origin[a] = nodes[a] - 3; /* origin node of the last cell */
        if (holds_face) {
            table->n[a] = 1;
            table->origin[a] = global_origin[a];
        } else {
            table->n[0] = table->n[1] = table->n[2] = 0;
        }
    }
    n_cells = table->n[0] * table->n[1] * table->n[2];
    n_rows_local = table->n[1] * table->n[2];
    n_global_cells = (PetscInt64)global_cells[0] * global_cells[1] * global_cells[2];
    while (quantum_bits > 0 && ((PetscInt64)1 << (62 - quantum_bits)) < n_global_cells) --quantum_bits;

    // --- 1. Cell measures, and the density envelope ---
    PetscCall(PetscMalloc2(n_cells, &measure, n_cells, &quanta));
    PetscCall(PetscMalloc2(n_cells, &table->pid_lo, n_cells, &table->pid_hi));
    if (build_bound) PetscCall(PetscMalloc1(n_cells, &table->bound));
    for (PetscInt c = 0; c < n_cells; ++c) {
        const PetscInt i = table->origin[0] + c % table->n[0];
        const PetscInt j = table->origin[1] + (c / table->n[0]) % table->n[1];
        const PetscInt k = table->origin[2] + c / (table->n[0] * table->n[1]);
        PetscReal      density;
        measure[c] = 0.0;
        for (PetscInt g = 0; g < 8; ++g) {
            const PetscReal l[3] = {gauss[g & 1], gauss[(g >> 1) & 1], gauss[(g >> 2) & 1]};
            PetscCall(SeedingMeasureDensity(user, coor, face, i, j, k, l, &density));
            measure[c] += 0.125 * density;
        }
        PetscCheck(measure[c] > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                   "Cell with origin node (%d,%d,%d) has %s %g; particle seeding needs positive %s.",
                   (int)i, (int)j, (int)k, face < 0 ? "volume" : "inlet area", (double)measure[c],
                   face < 0 ? "volumes" : "areas");
        local_max = PetscMax(local_max, measure[c]);
        if (build_bound) {
            table->bound[c] = 0.0;
            for (PetscInt q = 0; q < 27; ++q) {
                const PetscReal l[3] = {0.5 * (q % 3), 0.5 * ((q / 3) % 3), 0.5 * (q / 9)};
                PetscCall(SeedingMeasureDensity(user, coor, face, i, j, k, l, &density));
                table->bound[c] = PetscMax(table->bound[c], density);
            }
        }
    }
    {
        PetscInt64 local_cells = n_cells, owned_cells = 0;
        PetscCallMPI(MPI_Allreduce(&local_cells, &owned_cells, 1, MPIU_INT64, MPI_SUM, PETSC_COMM_WORLD));
        PetscCheck(owned_cells == n_global_cells, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "Seeding cells owned across ranks (%" PetscInt64_FMT ") differ from the global count (%" PetscInt64_FMT ").",
                   owned_cells, n_global_cells);
    }
    PetscCallMPI(MPI_Allreduce(&local_max, &global_max, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCheck(global_max > 0.0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               face < 0 ? "Volumetric particle seeding found no owned volume on any rank."
                        : "Surface particle seeding found no inlet area on any rank; it needs an INLET face.");

    // --- 2. Integer measures, summed over rows and over each rank's segment of a row ---
    const PetscInt n_rows = global_cells[1] * global_cells[2];
    PetscCall(PetscCalloc4(n_rows, &rows, n_rows, &row_before, n_rows_local, &segment, n_rows_local, &segment_before));
    for (PetscInt c = 0; c < n_cells; ++c) {
        const PetscInt lj = (c / table->n[0]) % table->n[1], lk = c / (table->n[0] * table->n[1]);
        const PetscInt row = (table->origin[1] + lj - global_origin[1]) + global_cells[1] * (table->origin[2] + lk - global_origin[2]);
        quanta[c] = (uint64_t)PetscMax(1.0, PetscFloorReal(ldexp((double)(measure[c] / global_max), (int)quantum_bits) + 0.5));
        rows[row] += quanta[c];
        segment[lj + table->n[1] * lk] += quanta[c];
    }
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, rows, (PetscMPIInt)n_rows, MPI_UINT64_T, MPI_SUM, PETSC_COMM_WORLD));
    for (PetscInt r = 0; r < n_rows; ++r) { row_before[r] = total; total += rows[r]; }
    // Ranks with the same j and k ranges share rows; their order along i is the order in the row.
    PetscCallMPI(MPI_Comm_split(PETSC_COMM_WORLD, n_cells > 0 ? (PetscMPIInt)(table->origin[1] + nodes[1] * table->origin[2]) : MPI_UNDEFINED,
                                (PetscMPIInt)table->origin[0], &row_comm));
    if (row_comm != MPI_COMM_NULL) {
        PetscCallMPI(MPI_Exscan(segment, segment_before, (PetscMPIInt)n_rows_local, MPI_UINT64_T, MPI_SUM, row_comm));
        PetscCallMPI(MPI_Comm_rank(row_comm, &row_rank));
        if (row_rank == 0) PetscCall(PetscArrayzero(segment_before, n_rows_local)); /* MPI_Exscan leaves rank 0's result undefined */
        PetscCallMPI(MPI_Comm_free(&row_comm));
    }

    // --- 3. Stratified particle IDs: floor(N P / T) at each cell boundary ---
    for (PetscInt lrow = 0; lrow < n_rows_local; ++lrow) {
        const PetscInt lj = lrow % table->n[1], lk = lrow / table->n[1];
        const PetscInt row = (table->origin[1] + lj - global_origin[1]) + global_cells[1] * (table->origin[2] + lk - global_origin[2]);
        uint64_t       before = row_before[row] + segment_before[lrow];
        for (PetscInt li = 0; li < table->n[0]; ++li) {
            const PetscInt c = li + table->n[0] * lrow;
            table->pid_lo[c] = (PetscInt64)(((unsigned __int128)numParticles * before) / total);
            before += quanta[c];
            table->pid_hi[c] = (PetscInt64)(((unsigned __int128)numParticles * before) / total);
            table->count += (PetscInt)(table->pid_hi[c] - table->pid_lo[c]);
        }
    }
    PetscCall(PetscFree4(rows, row_before, segment, segment_before));
    PetscCall(PetscFree2(measure, quanta));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/** @brief Free the arrays of a SeedingCellTable. */
static PetscErrorCode FreeSeedingCellTable(SeedingCellTable *table)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PetscFree2(table->pid_lo, table->pid_hi));
    PetscCall(PetscFree(table->bound));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Draw one particle's seeding point inside its cell, uniform in the cell's measure.
 * @details Uniform logical coordinates are accepted with probability density / `bound`. Every
 *          draw is keyed on the run seed, the particle ID and the attempt, so the point depends
 *          only on the particle, not on the rank that seeds it. A surface point is then put just
 *          inside the cell, 1e-6 of the cell from the face.
 * @param[in]  user   Grid context.
 * @param[in]  coor   Ghosted node-coordinate array.
 * @param[in]  table  A table built with `build_bound` set.
 * @param[in]  c      The cell, an index into the table.
 * @param[in]  pid    The particle's ID.
 * @param[out] ci_metric_lnode_out,cj_metric_lnode_out,ck_metric_lnode_out Origin node of the cell.
 * @param[out] xi_metric_logic_out,eta_metric_logic_out,zta_metric_logic_out Logical point in it.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode DrawSeedingPoint(
    UserCtx *user, const Cmpnts ***coor, const SeedingCellTable *table, PetscInt c, PetscInt64 pid,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out)
{
    PetscReal logic[3], u[4];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    *ci_metric_lnode_out = table->origin[0] + c % table->n[0];
    *cj_metric_lnode_out = table->origin[1] + (c / table->n[0]) % table->n[1];
    *ck_metric_lnode_out = table->origin[2] + c / (table->n[0] * table->n[1]);

    // Accept uniform logicals with probability density / bound; a cell whose density varies
    // little accepts almost every draw. The cap only guards a degenerate envelope.
    for (PetscInt attempt = 0; attempt < 1000; ++attempt) {
        PetscReal density;
        uint64_t  hash = PicurvParticleDrawHash(user->simCtx->particleRandomSeed, PICURV_DRAW_SPACE_PLACEMENT, 0, pid, attempt);
        for (int w = 0; w < 4; ++w) u[w] = PicurvUnitFromHash(hash = PicurvMix64(hash ^ (uint64_t)w));
        // Ensure logical coordinates are strictly within [0,1) for robustness with MetricLogicalToPhysical
        for (PetscInt d = 0; d < 3; ++d) logic[d] = PetscMax(PetscMin(u[d], 1.0 - 1.0e-7), 0.0);
        PetscCall(SeedingMeasureDensity(user, coor, table->face, *ci_metric_lnode_out, *cj_metric_lnode_out,
                                        *ck_metric_lnode_out, logic, &density));
        if (u[3] * table->bound[c] <= density) break;
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
    // Owned cells and the particles each seeds, for volumetric and surface seeding
    SeedingCellTable seed_table;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    SimCtx *simCtx = user->simCtx;

    // --- 1. Input Validation and Basic Setup ---
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null user.");
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
    // Volume seeds the owned cells, Surface the owned cells on the inlet face; the table names
    // each cell's particles, the same split CreateParticleSwarm sized the swarm with.
    const PetscBool seeds_from_table = (PetscBool)(simCtx->ParticleInitialization == PARTICLE_INIT_VOLUME ||
                                                   simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM);
    PetscInt   seed_cell = 0;  // Cell the next table-seeded particle goes in
    PetscInt64 seed_pid  = 0;  // That particle's ID
    if (seeds_from_table) {
        const PetscInt face = (simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM) ? (PetscInt)user->identifiedInletBCFace : -1;
        ierr = BuildSeedingCellTable(user, coor_nodes_local_array, face, simCtx->np, PETSC_TRUE, &seed_table); CHKERRQ(ierr);
        PetscCheck(seed_table.count == particlesPerProcess, PETSC_COMM_SELF, PETSC_ERR_PLIB,
                   "Seeding table holds %" PetscInt_FMT " particles but the swarm was sized for %" PetscInt_FMT ".",
                   seed_table.count, particlesPerProcess);
        if (seed_table.n[0] * seed_table.n[1] * seed_table.n[2] > 0) seed_pid = seed_table.pid_lo[0];
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Seeding %" PetscInt_FMT " particles stratified over %" PetscInt_FMT " owned cells by %s.\n",
                  rank, particlesPerProcess, seed_table.n[0] * seed_table.n[1] * seed_table.n[2],
                  face < 0 ? "volume" : "inlet area");
    }

    // --- 3. Get Access to Swarm Fields ---
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION),       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID),    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_LOCATION_STATUS),NULL,NULL,(void**)&status_field); CHKERRQ(ierr);

    // --- 4. Determine Starting Global PID for this Rank (point source) ---
    // Point-source PIDs are consecutive across ranks in rank order, so each rank starts where
    // the ranks before it end. Table seeding takes each particle's ID from its cell instead.
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
    for (PetscInt p = 0; p < particlesPerProcess; p++) {
        PetscInt  ci_metric_lnode, cj_metric_lnode, ck_metric_lnode; 
        PetscReal xi_metric_logic, eta_metric_logic, zta_metric_logic; 
        Cmpnts    phys_coords = {0.0, 0.0, 0.0}; 
        PetscBool particle_placed_by_this_rank = PETSC_FALSE; 
        PetscInt64 pid = (PetscInt64)base_pid_for_rank + p;

        if (seeds_from_table) { // --- 5.a. Volumetric or Surface Initialization ---
            while (seed_pid == seed_table.pid_hi[seed_cell]) seed_pid = seed_table.pid_lo[++seed_cell]; // skip filled (or empty) cells
            pid = seed_pid++;
            ierr = DrawSeedingPoint(user, coor_nodes_local_array, &seed_table, seed_cell, pid,
                                    &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                    &xi_metric_logic, &eta_metric_logic, &zta_metric_logic); CHKERRQ(ierr);
            ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                    ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                    xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                    &phys_coords); CHKERRQ(ierr);
            particle_placed_by_this_rank = PETSC_TRUE;
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

        particleIDs[p]         = pid;
        // Start the first settlement's walk in the cell the particle was placed in; a point-source
        // particle no rank located has no cell and is located by search.
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
            // A point source no rank located.
            LOG_PARTICLE(simCtx, particleIDs[p], "seed",
                         "%s: not placed by rank %d; parked at (%.6e, %.6e, %.6e) for the settlement to locate\n",
                         ParticleInitializationToString(simCtx->ParticleInitialization), (int)rank,
                         phys_coords.x, phys_coords.y, phys_coords.z);
        }
    }

    // --- 6. Restore Pointers and Cleanup ---
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION),       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID),    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_LOCATION_STATUS), NULL, NULL, (void**)&status_field); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
    ierr = FreeSeedingCellTable(&seed_table); CHKERRQ(ierr);

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
                                            BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    SimCtx *simCtx = user->simCtx;
    
    // --- 0. Input Validation ---
    if (!user || !bboxlist) {
        LOG_ALLOW(GLOBAL, LOG_ERROR, "Null user or bboxlist.\n");
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
    // Placement draws are keyed on (seed, PID), so no generator state is passed in.
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Calling InitializeParticleBasicProperties.\n");
    ierr = InitializeParticleBasicProperties(user, particlesPerProcess, bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully initialized basic particle properties.\n");

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
 * @brief Number of particles this rank seeds for volumetric or surface seeding.
 * @details Collective. The count is the rank's share of the stratified split
 *          BuildSeedingCellTable makes, so a rank seeds in proportion to the volume or inlet
 *          area it owns, none if it owns no part of the inlet, and the counts sum to
 *          `numParticles` exactly.
 * @param[in]  user          Finest-level context whose DMDA ownership defines the split.
 * @param[in]  face          -1 for the volume, otherwise the inlet BCFace.
 * @param[in]  numParticles  Global particle count.
 * @param[out] localCount    Particles this rank seeds.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode DistributeParticlesByOwnedMeasure(UserCtx *user, PetscInt face, PetscInt numParticles, PetscInt *localCount)
{
    SeedingCellTable table;
    const Cmpnts  ***coor = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCoor, (void *)&coor));
    PetscCall(BuildSeedingCellTable(user, coor, face, numParticles, PETSC_FALSE, &table));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCoor, (void *)&coor));
    *localCount = table.count;
    PetscCall(FreeSeedingCellTable(&table));
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
        ierr = AssignInitialPropertiesToSwarm(user, particlesPerProcess, simCtx->bboxlist); CHKERRQ(ierr);

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
