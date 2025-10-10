// in src/grid.c

#include "grid.h"
#include "logging.h"

#define BBOX_TOLERANCE 1e-6

#undef __FUNCT__
#define __FUNCT__ "ParseAndSetGridInputs"
/**
 * @brief Determines the grid source and calls the appropriate parsing routine.
 *
 * This function acts as a router. It checks the global `simCtx->generate_grid`
 * flag (accessed via the `user->simCtx` back-pointer) to decide whether to
 * call the parser for a programmatically generated grid or for a grid defined
 * in a file.
 *
 * @param user Pointer to the `UserCtx` for a specific block. The function will
 *             populate the geometric fields within this struct.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode ParseAndSetGridInputs(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx; // Get the global context via the back-pointer

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    if (simCtx->generate_grid) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d: Block %d is programmatically generated. Calling generation parser.\n", simCtx->rank, user->_this);
        ierr = ReadGridGenerationInputs(user); CHKERRQ(ierr);
    } else {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d: Block %d is file-based. Calling file parser.\n", simCtx->rank, user->_this);
        ierr = ReadGridFile(user); CHKERRQ(ierr);
    }

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DefineAllGridDimensions"
/**
 * @brief Orchestrates the parsing and setting of grid dimensions for all blocks.
 *
 * This function serves as the high-level entry point for defining the geometric
 * properties of each grid block in the simulation. It iterates through every
 * block defined by `simCtx->block_number`.
 *
 * For each block, it performs two key actions:
 * 1.  It explicitly sets the block's index (`_this`) in the corresponding `UserCtx`
 *     struct for the finest multigrid level. This makes the context "self-aware".
 * 2.  It calls a helper function (`ParseAndSetGridInputs`) to handle the detailed
 *     work of parsing options or files to populate the rest of the geometric
 *     properties for that specific block (e.g., `IM`, `Min_X`, `rx`).
 *
 * @param simCtx The master SimCtx, which contains the number of blocks and the
 *               UserCtx hierarchy to be configured.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode DefineAllGridDimensions(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscInt       nblk = simCtx->block_number;
    UserCtx        *finest_users;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    if (simCtx->usermg.mglevels == 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "MG levels not set. Cannot get finest_users.");
    }
    // Get the UserCtx array for the finest grid level
    finest_users = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Defining grid dimensions for %d blocks...\n", nblk);

    // Loop over each block to configure its grid dimensions and geometry.
    for (PetscInt bi = 0; bi < nblk; bi++) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d: --- Configuring Geometry for Block %d ---\n", simCtx->rank, bi);

        // Before calling any helpers, set the block index in the context.
        // This makes the UserCtx self-aware of which block it represents.
        LOG_ALLOW(GLOBAL,LOG_DEBUG,"finest_users->_this = %d, bi = %d\n",finest_users[bi]._this,bi);
        //finest_user[bi]._this = bi;

        // Call the helper function for this specific block. It can now derive
        // all necessary information from the UserCtx pointer it receives.
        ierr = ParseAndSetGridInputs(&finest_users[bi]); CHKERRQ(ierr);
    }

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitializeSingleGridDM"
/**
 * @brief Creates the DMDA objects (da and fda) for a single UserCtx.
 *
 * This function is a direct adaptation of the core logic in `MGDACreate`. It
 * creates the scalar (`da`) and vector (`fda`) DMs for a single grid level.
 *
 * If a `coarse_user` context is provided, it performs the critical processor
 * alignment calculation from the legacy code. This ensures the new (fine) DM
 * aligns with the coarse DM for multigrid efficiency. If `coarse_user` is NULL,
 * it creates the DM with a default PETSc decomposition, intended for the
 * coarsest grid level.
 *
 * @param user The UserCtx for which the DMs will be created. Its IM, JM, KM fields must be pre-populated.
 * @param coarse_user The UserCtx of the next-coarser grid level, or NULL if `user` is the coarsest level.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode InitializeSingleGridDM(UserCtx *user, UserCtx *coarse_user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;

    DMBoundaryType xperiod = (simCtx->i_periodic) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType yperiod = (simCtx->j_periodic) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    DMBoundaryType zperiod = (simCtx->k_periodic) ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
    PetscInt stencil_width = 2; // Stencil width is 2 in the legacy code

    PetscInt *lx = NULL, *ly = NULL, *lz = NULL;
    PetscInt m, n, p;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    if (coarse_user) {
        // --- This is a FINE grid; it must be aligned with the COARSE grid ---
        LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: [Aligning DM] for block %d level %d (size %dx%dx%d) with level %d\n", simCtx->rank, user->_this, user->thislevel, user->IM, user->JM, user->KM, coarse_user->thislevel);

        DMDAGetInfo(coarse_user->da, NULL, NULL, NULL, NULL, &m, &n, &p, NULL, NULL, NULL, NULL, NULL, NULL);
        LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d:   Coarse grid processor decomposition is %d x %d x %d\n", simCtx->rank, m, n, p);

        // This is the core logic from MGDACreate to ensure processor alignment.
        PetscInt *lx_contrib, *ly_contrib, *lz_contrib;
        ierr = PetscMalloc3(m, &lx_contrib, n, &ly_contrib, p, &lz_contrib); CHKERRQ(ierr);
        ierr = PetscMemzero(lx_contrib, m * sizeof(PetscInt)); CHKERRQ(ierr);
        ierr = PetscMemzero(ly_contrib, n * sizeof(PetscInt)); CHKERRQ(ierr);
        ierr = PetscMemzero(lz_contrib, p * sizeof(PetscInt)); CHKERRQ(ierr);

        DMDALocalInfo info;
        DMDAGetLocalInfo(coarse_user->da, &info);
        PetscInt xs = info.xs, xe = info.xs + info.xm, mx = info.mx;
        PetscInt ys = info.ys, ye = info.ys + info.ym, my = info.my;
        PetscInt zs = info.zs, ze = info.zs + info.zm, mz = info.mz;
        
        PetscMPIInt rank;
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
        PetscInt proc_i = rank % m;
        PetscInt proc_j = (rank / m) % n;
        PetscInt proc_k = rank / (m * n);

        // --- X-Direction Logic (Identical to MGDACreate) ---
        if (user->isc) lx_contrib[proc_i] = (xe - xs);
        else {
            if (m == 1) lx_contrib[0] = user->IM + 1;
            else if (xs == 0) lx_contrib[0] = 2 * xe - 1;
            else if (xe == mx) lx_contrib[proc_i] = user->IM + 1 - (2 * xs - 1);
            else lx_contrib[proc_i] = (xe - xs) * 2;
        }

        // --- Y-Direction Logic (Identical to MGDACreate) ---
        if (user->jsc) ly_contrib[proc_j] = (ye - ys);
        else {
            if (n == 1) ly_contrib[0] = user->JM + 1;
            else if (ys == 0) ly_contrib[0] = 2 * ye - 1;
            else if (ye == my) ly_contrib[proc_j] = user->JM + 1 - (2 * ys - 1);
            else ly_contrib[proc_j] = (ye - ys) * 2;
        }

        // --- Z-Direction Logic (Identical to MGDACreate) ---
        if (user->ksc) lz_contrib[proc_k] = (ze - zs);
        else {
            if (p == 1) lz_contrib[0] = user->KM + 1;
            else if (zs == 0) lz_contrib[0] = 2 * ze - 1;
            else if (ze == mz) lz_contrib[proc_k] = user->KM + 1 - (2 * zs - 1);
            else lz_contrib[proc_k] = (ze - zs) * 2;
        }
        LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d:   Calculated this rank's node contribution to fine grid: lx=%d, ly=%d, lz=%d\n", simCtx->rank, lx_contrib[proc_i], ly_contrib[proc_j], lz_contrib[proc_k]);

        // Allocate the final distribution arrays and Allreduce to get the global distribution
        ierr = PetscMalloc3(m, &lx, n, &ly, p, &lz); CHKERRQ(ierr);
        ierr = MPI_Allreduce(lx_contrib, lx, m, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Allreduce(ly_contrib, ly, n, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Allreduce(lz_contrib, lz, p, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        ierr = PetscFree3(lx_contrib, ly_contrib, lz_contrib); CHKERRQ(ierr);

    } else {
        // --- CASE 2: This is the COARSEST grid; use default or user-specified decomposition ---
        LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Creating coarsest DM for block %d level %d (size %dx%dx%d)\n", simCtx->rank, user->_this, user->thislevel, user->IM, user->JM, user->KM);
        m = simCtx->da_procs_x;
        n = simCtx->da_procs_y;
        p = simCtx->da_procs_z;
        // lx, ly, lz are NULL, so DMDACreate3d will use the m,n,p values.
    }

    // --- Create the DMDA for the current UserCtx ---
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d:   Calling DMDACreate3d...\n", simCtx->rank);
    ierr = DMDACreate3d(PETSC_COMM_WORLD, xperiod, yperiod, zperiod, DMDA_STENCIL_BOX,
                          user->IM + 1, user->JM + 1, user->KM + 1,
                          m, n, p,
                          1, stencil_width, lx, ly, lz, &user->da); CHKERRQ(ierr);
    
    if (coarse_user) {
        ierr = PetscFree3(lx, ly, lz); CHKERRQ(ierr);
    }
    
    // --- Standard DM setup applicable to all levels ---
    ierr = DMSetUp(user->da); CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(user->da, &user->fda); CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(user->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(user->da, &user->info); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d:   DM creation for block %d level %d complete.\n", simCtx->rank, user->_this, user->thislevel);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "InitializeAllGridDMs"
/**
 * @brief Orchestrates the creation of DMDA objects for every block and multigrid level.
 *
 * This function systematically builds the entire DMDA hierarchy. It first
 * calculates the dimensions (IM, JM, KM) for all coarse grids based on the
 * finest grid's dimensions and the semi-coarsening flags. It then iterates
 * from the coarsest to the finest level, calling a powerful helper function
 * (`InitializeSingleGridDM`) to create the DMs for each block, ensuring that
 * finer grids are properly aligned with their coarser parents for multigrid efficiency.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode InitializeAllGridDMs(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    MGCtx          *mgctx = usermg->mgctx;
    PetscInt       nblk = simCtx->block_number;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Creating DMDA objects for all levels and blocks...\n");

    // --- Part 1: Calculate Coarse Grid Dimensions & VALIDATE ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Calculating and validating coarse grid dimensions...\n");
    for (PetscInt level = usermg->mglevels - 2; level >= 0; level--) {
        for (PetscInt bi = 0; bi < nblk; bi++) {
            UserCtx *user_coarse = &mgctx[level].user[bi];
            UserCtx *user_fine   = &mgctx[level + 1].user[bi];

            user_coarse->IM = user_fine->isc ? user_fine->IM : (user_fine->IM + 1) / 2;
            user_coarse->JM = user_fine->jsc ? user_fine->JM : (user_fine->JM + 1) / 2;
            user_coarse->KM = user_fine->ksc ? user_fine->KM : (user_fine->KM + 1) / 2;

            LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Block %d, Level %d dims calculated: %d x %d x %d\n",
                      simCtx->rank, bi, level, user_coarse->IM, user_coarse->JM, user_coarse->KM);

            // Validation check from legacy MGDACreate to ensure coarsening is possible
            PetscInt check_i = user_coarse->IM * (2 - user_coarse->isc) - (user_fine->IM + 1 - user_coarse->isc);
            PetscInt check_j = user_coarse->JM * (2 - user_coarse->jsc) - (user_fine->JM + 1 - user_coarse->jsc);
            PetscInt check_k = user_coarse->KM * (2 - user_coarse->ksc) - (user_fine->KM + 1 - user_coarse->ksc);

            if (check_i + check_j + check_k != 0) {
	      //      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
	      //           "Grid at level %d, block %d cannot be coarsened from %dx%dx%d to %dx%dx%d with the given semi-coarsening flags. Check grid dimensions.",
              //          level, bi, user_fine->IM, user_fine->JM, user_fine->KM, user_coarse->IM, user_coarse->JM, user_coarse->KM);
	      LOG(GLOBAL,LOG_WARNING,"WARNING: Grid at level %d, block %d can't be consistently coarsened further.\n", level, bi);
            }
        }
    }

    // --- Part 2: Create DMs from Coarse to Fine for each Block ---
    for (PetscInt bi = 0; bi < nblk; bi++) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "--- Creating DMs for Block %d ---\n", bi);

        // Create the coarsest level DM first (passing NULL for the coarse_user)
        ierr = InitializeSingleGridDM(&mgctx[0].user[bi], NULL); CHKERRQ(ierr);

        // Create finer level DMs, passing the next-coarser context for alignment
        for (PetscInt level = 1; level < usermg->mglevels; level++) {
            ierr = InitializeSingleGridDM(&mgctx[level].user[bi], &mgctx[level-1].user[bi]); CHKERRQ(ierr);
        }
    }

    // --- Optional: View the finest DM for debugging verification ---
    if (get_log_level() >= LOG_DEBUG) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "--- Viewing Finest DMDA (Level %d, Block 0) ---\n", usermg->mglevels - 1);
        ierr = DMView(mgctx[usermg->mglevels - 1].user[0].da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "DMDA object creation complete.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

// Forward declarations for the static helper functions within this file.
static PetscErrorCode SetFinestLevelCoordinates(UserCtx *user);
static PetscErrorCode GenerateAndSetCoordinates(UserCtx *user);
static PetscErrorCode ReadAndSetCoordinates(UserCtx *user, FILE *fd);
static PetscErrorCode RestrictCoordinates(UserCtx *coarse_user, UserCtx *fine_user);

#undef __FUNCT__
#define __FUNCT__ "AssignAllGridCoordinates"
/**
 * @brief Orchestrates the assignment of physical coordinates to all DMDA objects.
 *
 * This function manages the entire process of populating the coordinate vectors
 * for every DMDA across all multigrid levels and blocks. It follows a two-part
 * strategy that is essential for multigrid methods:
 *
 * 1.  **Populate Finest Level:** It first loops through each block and calls a
 *     helper (`SetFinestLevelCoordinates`) to set the physical coordinates for
 *     the highest-resolution grid (the finest multigrid level).
 * 2.  **Restrict to Coarser Levels:** It then iterates downwards from the finest
 *     level, calling a helper (`RestrictCoordinates`) to copy the coordinate
 *     values from the fine grid nodes to their corresponding parent nodes on the
 *     coarser grids. This ensures all levels represent the exact same geometry.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode AssignAllGridCoordinates(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    PetscInt       nblk = simCtx->block_number;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Assigning physical coordinates to all grid DMs...\n");

    // --- Part 1: Populate the Finest Grid Level ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Setting coordinates for the finest grid level (%d)...\n", usermg->mglevels - 1);
    for (PetscInt bi = 0; bi < nblk; bi++) {
        UserCtx *fine_user = &usermg->mgctx[usermg->mglevels - 1].user[bi];
        ierr = SetFinestLevelCoordinates(fine_user); CHKERRQ(ierr);
        if(get_log_level()==LOG_DEBUG){
            LOG_ALLOW(GLOBAL,LOG_DEBUG,"The Finest level coordinates for block %d have been set.\n",bi);
            ierr = LOG_FIELD_MIN_MAX(fine_user,"Coordinates");
        }
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Finest level coordinates have been set for all blocks.\n");

    // --- Part 2: Restrict Coordinates to Coarser Levels ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Restricting coordinates to coarser grid levels...\n");
    for (PetscInt level = usermg->mglevels - 2; level >= 0; level--) {
        for (PetscInt bi = 0; bi < nblk; bi++) {
            UserCtx *coarse_user = &usermg->mgctx[level].user[bi];
            UserCtx *fine_user   = &usermg->mgctx[level + 1].user[bi];
            ierr = RestrictCoordinates(coarse_user, fine_user); CHKERRQ(ierr);

            if(get_log_level==LOG_DEBUG) ierr = LOG_FIELD_MIN_MAX(coarse_user,"Coordinates");
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Physical coordinates assigned to all grid levels and blocks.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SetFinestLevelCoordinates"
/**
 * @brief A router that populates the coordinates for a single finest-level DMDA.
 *
 * This function orchestrates the coordinate setting for one block. It checks the
 * global `generate_grid` flag and calls the appropriate helper for either
 * programmatic generation or reading from a file.
 *
 * After the local coordinate vector is populated by a helper, this function
 * performs the necessary DMLocalToGlobal and DMGlobalToLocal scatters to ensure
 * that the ghost node coordinate values are correctly communicated and updated
 * across all MPI ranks.
 *
 * @param user The UserCtx for a specific block on the finest level.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode SetFinestLevelCoordinates(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Setting finest level coordinates for block %d...\n", simCtx->rank, user->_this);

    if (simCtx->generate_grid) {
        ierr = GenerateAndSetCoordinates(user); CHKERRQ(ierr);
    } else {

      FILE *grid_file_handle = NULL;
      // Only Rank 0 opens the file.
      if (simCtx->rank == 0) {
	grid_file_handle = fopen(simCtx->grid_file, "r");
	if (!grid_file_handle) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", simCtx->grid_file);
        
	// Now, on Rank 0, we skip the entire header section once.
	// This is the logic from your modern code's AssignGridCoordinates.
	PetscInt headerLines = simCtx->block_number + 2; // 1 for nblk, plus one for each block's dims
	char dummy_buffer[2048];
	for (PetscInt s = 0; s < headerLines; ++s) {
	  if (!fgets(dummy_buffer, sizeof(dummy_buffer), grid_file_handle)) {
	    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Unexpected EOF while skipping grid header");
	  }
	}
	LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank 0: Skipped %d header lines, now at coordinate data.\n", headerLines);
      }

      // We now call the coordinate reader, passing the file handle.
      // It's responsible for reading its block's data and broadcasting.
      ierr = ReadAndSetCoordinates(user, grid_file_handle); CHKERRQ(ierr);
        
      // Only Rank 0, which opened the file, should close it.
      if (simCtx->rank == 0) {
	fclose(grid_file_handle);
      }
    }
    
    // After populating the local coordinate vector, we must perform a
    // Local-to-Global and then Global-to-Local scatter to correctly
    // populate the ghost node coordinates across process boundaries.
    Vec gCoor, lCoor;
    ierr = DMGetCoordinates(user->da, &gCoor); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(user->da, &lCoor); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   Scattering coordinates to update ghost nodes for block %d...\n", simCtx->rank, user->_this);
    ierr = DMLocalToGlobalBegin(user->fda, lCoor, INSERT_VALUES, gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, lCoor, INSERT_VALUES, gCoor); CHKERRQ(ierr);
    
    ierr = DMGlobalToLocalBegin(user->fda, gCoor, INSERT_VALUES, lCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, gCoor, INSERT_VALUES, lCoor); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}
/**
 * @brief Computes a stretched coordinate along one dimension.
 *
 * This function computes a coordinate based on a geometric stretching ratio.
 * If the ratio (r) is 1.0, a uniform distribution is used:
 *     x(i) = L * (i/N)
 *
 * If r != 1.0, a geometric stretching is applied:
 *     x(i) = L * [ (r^(i/N) - 1 ) / (r - 1) ]
 *
 * Here:
 * - i   : The current index along the dimension.
 * - N   : The total number of divisions along that dimension.
 * - L   : The length of the domain along that dimension.
 * - r   : The stretching ratio. r > 1.0 stretches the grid in a geometric fashion
 *         increasing spacing away from the start, whereas 0 < r < 1.0 would
 *         cluster points near the start.
 *
 * @param[in] i Current index (0 <= i <= N).
 * @param[in] N Number of segments along the dimension.
 * @param[in] L Total length of the domain.
 * @param[in] r Stretching ratio.
 *
 * @return PetscReal The computed coordinate at index i.
 *
 * @note This function does not return a PetscErrorCode because it
 *       does not allocate memory or call PETSc routines that can fail.
 *       It is just a helper calculation function.
 **/
static inline PetscReal ComputeStretchedCoord(PetscInt i, PetscInt N, PetscReal L, PetscReal r)
{
    if (N <=1) return 0.0;
    PetscReal fraction = (PetscReal)i / ((PetscReal)N - 1.0);
    if (PetscAbsReal(r - 1.0) < 1.0e-9) { // Use a tolerance for float comparison
        return L * fraction;
    } else {
        return L * (PetscPowReal(r, fraction) - 1.0) / (r - 1.0);
    }
}

#undef __FUNCT__
#define __FUNCT__ "GenerateAndSetCoordinates"
/**
 * @brief Programmatically generates and sets grid coordinates based on user parameters.
 *
 * This function populates the local coordinate vector of the provided `UserCtx`
 * using the geometric properties (`Min_X`, `Max_X`, `rx`, etc.) that were parsed
 * from command-line options. It supports non-linear grid stretching.
 *
 * @param user The UserCtx for a specific block.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode GenerateAndSetCoordinates(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Cmpnts       ***coor;
    Vec            lCoor;
    
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Generating coordinates for block %d...\n", user->simCtx->rank, user->_this);
    
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(user->da, &lCoor); CHKERRQ(ierr);

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Local Info for block %d - X range - [%d,%d], Y range - [%d,%d], Z range - [%d,%d]\n",
              user->simCtx->rank, user->_this, xs, xe, ys, ye, zs, ze);
    ierr = VecSet(lCoor, 0.0); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, lCoor, &coor); CHKERRQ(ierr);
    
    PetscReal Lx = user->Max_X - user->Min_X;
    PetscReal Ly = user->Max_Y - user->Min_Y;
    PetscReal Lz = user->Max_Z - user->Min_Z;

    // Loop over the local nodes, including ghost nodes, owned by this process.
    for (PetscInt k = zs; k < ze; k++) {
        for (PetscInt j = ys; j < ye; j++) {
            for (PetscInt i = xs; i < xe; i++) {
	            if(k < user->KM && j < user->JM && i < user->IM){
                    coor[k][j][i].x = user->Min_X + ComputeStretchedCoord(i, user->IM, Lx, user->rx);
                    coor[k][j][i].y = user->Min_Y + ComputeStretchedCoord(j, user->JM, Ly, user->ry);
                    coor[k][j][i].z = user->Min_Z + ComputeStretchedCoord(k, user->KM, Lz, user->rz);
	            }
	        }  
        }
    }

    /// DEBUG: This verifies the presence of a last "unphysical" layer of coordinates.
    /*
    PetscInt KM = user->KM;
    for (PetscInt j = ys; j < ye; j++){
        for(PetscInt i = xs; i < xe; i++){
            LOG_ALLOW(GLOBAL,LOG_DEBUG,"coor[%d][%d][%d].(x,y,z) = %le,%le,%le",KM,j,i,coor[KM][j][i].x,coor[KM][j][i].y,coor[KM][j][i].z);
        }
    }
    */



    ierr = DMDAVecRestoreArray(user->fda, lCoor, &coor); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "ReadAndSetCoordinates"
/**
 * @brief Reads physical coordinates from a file and populates the DMDA for a specific block.
 *
 * This function handles the collective read of an interleaved (X Y Z per line)
 * multi-block grid file. It assumes the file header (nblk, dimensions) has
 * already been processed by ReadGridFile.
 *
 * The process is robust for parallel execution:
 * 1.  Rank 0 opens the grid file.
 * 2.  It intelligently skips past the header section and the coordinate data
 *     for all blocks *preceding* the current block being processed (`user->_this`).
 * 3.  It then reads the entire coordinate data for the *current* block into
 *     a single contiguous buffer `gc`.
 * 4.  This global buffer `gc` is broadcast to all other MPI ranks.
 * 5.  Each rank then loops through its local (owned + ghost) node indices,
 *     calculates the corresponding index in the global `gc` array, and copies
 *     the (x,y,z) coordinates into its local PETSc coordinate vector.
 *
 * @param user The UserCtx for a specific block. Its `_this` field must be set,
 *             and its IM, JM, KM fields must be correctly pre-populated.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode ReadAndSetCoordinates(UserCtx *user, FILE *fd)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    PetscMPIInt    rank = simCtx->rank;
    PetscInt       block_index = user->_this;
    PetscInt       IM = user->IM, JM = user->JM, KM = user->KM;
    DMDALocalInfo  info;
    Cmpnts       ***coor;
    Vec            lCoor;
    PetscReal      *gc = NULL; // Global coordinate buffer, allocated on all ranks

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reading interleaved coordinates from file for block %d...\n",
              simCtx->rank, block_index);

    // 1. Allocate the buffer on ALL ranks to receive the broadcast data.
    // PetscInt n_nodes = (IM + 1) * (JM + 1) * (KM + 1);
    PetscInt n_nodes = (IM) * (JM) * (KM);
    ierr = PetscMalloc1(3 * n_nodes, &gc); CHKERRQ(ierr);

    // 2. Only Rank 0 opens the file and reads the data.
    if (rank == 0) {
        if (!fd) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Recieved a NULL file handle.\n");
                
        // Read the coordinate data for the CURRENT block.
        for (PetscInt k = 0; k < KM; k++) {
            for (PetscInt j = 0; j < JM; j++) {
                for (PetscInt i = 0; i < IM; i++) {
                    PetscInt base_index = 3 * ((k * (JM) + j) * (IM) + i);
                    if (fscanf(fd, "%le %le %le\n", &gc[base_index], &gc[base_index + 1], &gc[base_index + 2]) != 3) {
                        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Error reading coordinates for node (i,j,k)=(%d,%d,%d) in block %d", i, j, k, block_index);
                    }
                }
            }
        }
        
    }

    // 3. Broadcast the coordinate block for the current block to all other processes.
    ierr = MPI_Bcast(gc, 3 * n_nodes, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    // 4. Each rank populates its local portion of the coordinate vector.
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(user->da, &lCoor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, lCoor, &coor); CHKERRQ(ierr);

    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
	      if(k< KM && j < JM && i < IM){
		PetscInt base_idx = 3 * ((k * (JM) + j) * (IM) + i);
		coor[k][j][i].x = gc[base_idx];
                coor[k][j][i].y = gc[base_idx + 1];
                coor[k][j][i].z = gc[base_idx + 2];
	      }
	    }
        }
    }
    
    // 5. Clean up and restore.
    ierr = DMDAVecRestoreArray(user->fda, lCoor, &coor); CHKERRQ(ierr);
    ierr = PetscFree(gc); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RestrictCoordinates"
/**
 * @brief Populates coarse grid coordinates by restricting from a fine grid.
 *
 * This function is a direct adaptation of the coordinate restriction loop
 * from the legacy `MG_Initial` function. It ensures that the physical location
 * of a coarse grid node is identical to its corresponding parent node on the
 * fine grid. The mapping from coarse to fine index (`i` -> `ih`) is determined
 * by the semi-coarsening flags (`isc`, `jsc`, `ksc`) stored in the `UserCtx`.
 *
 * @param coarse_user The UserCtx for the coarse grid (destination).
 * @param fine_user The UserCtx for the fine grid (source).
 * @return PetscErrorCode
 */
static PetscErrorCode RestrictCoordinates(UserCtx *coarse_user, UserCtx *fine_user)
{
    PetscErrorCode ierr;
    Vec            c_lCoor, f_lCoor; // Coarse and Fine local coordinate vectors
    Cmpnts       ***c_coor;
    const Cmpnts ***f_coor; // Use const for read-only access
    DMDALocalInfo  c_info;
    PetscInt       ih, jh, kh; // Fine-grid indices corresponding to coarse-grid i,j,k

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Restricting coords from level %d to level %d for block %d\n",
                   fine_user->simCtx->rank, fine_user->thislevel, coarse_user->thislevel, coarse_user->_this);

    ierr = DMDAGetLocalInfo(coarse_user->da, &c_info); CHKERRQ(ierr);

    ierr = DMGetCoordinatesLocal(coarse_user->da, &c_lCoor); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(fine_user->da, &f_lCoor); CHKERRQ(ierr);
    
    ierr = DMDAVecGetArray(coarse_user->fda, c_lCoor, &c_coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fine_user->fda, f_lCoor, &f_coor); CHKERRQ(ierr);

    // Get the local owned range of the coarse grid.
    PetscInt xs = c_info.xs, xe = c_info.xs + c_info.xm;
    PetscInt ys = c_info.ys, ye = c_info.ys + c_info.ym;
    PetscInt zs = c_info.zs, ze = c_info.zs + c_info.zm;

    // Get the global dimensions of the coarse grid.
    PetscInt mx = c_info.mx, my = c_info.my, mz = c_info.mz;

    // If this process owns the maximum boundary node, contract the loop by one
    // to prevent the index doubling `2*i` from going out of bounds.
    // This is also ensuring we do not manipulate the unphysical layer of coors present in the finest level.
    if (xe == mx) xe--;
    if (ye == my) ye--;
    if (ze == mz) ze--;
    
    for (PetscInt k = zs; k < ze; k++) {
        for (PetscInt j = ys; j < ye; j++) {
            for (PetscInt i = xs; i < xe; i++) {
                // Determine the corresponding parent node index on the FINE grid,
                // respecting the semi-coarsening flags of the FINE grid's UserCtx.
                ih = coarse_user->isc ? i : 2 * i;
                jh = coarse_user->jsc ? j : 2 * j;
		kh = coarse_user->ksc ? k : 2 * k;

		//	LOG_ALLOW(GLOBAL,LOG_DEBUG," [kh][ih][jh] = %d,%d,%d  - k,j,i = %d,%d,%d.\n",kh,jh,ih,k,j,i);

                c_coor[k][j][i] = f_coor[kh][jh][ih];
            }
        }
    }

    ierr = DMDAVecRestoreArray(coarse_user->fda, c_lCoor, &c_coor); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fine_user->fda, f_lCoor, &f_coor); CHKERRQ(ierr);

    // After populating the local portion, scatter to update ghost regions.
    Vec c_gCoor;
    ierr = DMGetCoordinates(coarse_user->da, &c_gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(coarse_user->fda, c_lCoor, INSERT_VALUES, c_gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(coarse_user->fda, c_lCoor, INSERT_VALUES, c_gCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(coarse_user->fda, c_gCoor, INSERT_VALUES, c_lCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(coarse_user->fda, c_gCoor, INSERT_VALUES, c_lCoor); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeLocalBoundingBox"
/**
 * @brief Computes the local bounding box of the grid on the current process.
 *
 * This function calculates the minimum and maximum coordinates (x, y, z) of the
 * local grid points owned by the current MPI process. It iterates over the local
 * portion of the grid, examines each grid point's coordinates, and updates the
 * minimum and maximum values accordingly.
 *
 * The computed bounding box is stored in the provided `localBBox` structure,
 * and the `user->bbox` field is also updated with this bounding box for
 * consistency within the user context.
 *
 * @param[in]  user      Pointer to the user-defined context containing grid information.
 *                       This context must be properly initialized before calling this function.
 * @param[out] localBBox Pointer to the BoundingBox structure where the computed local bounding box will be stored.
 *                       The structure should be allocated by the caller.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure.
 */
PetscErrorCode ComputeLocalBoundingBox(UserCtx *user, BoundingBox *localBBox)
{
    PetscErrorCode ierr;
    PetscInt i, j, k;
    PetscMPIInt rank;
    PetscInt xs, ys, zs, xe, ye, ze;
    DMDALocalInfo info;
    Vec coordinates;
    Cmpnts ***coordArray;
    Cmpnts minCoords, maxCoords;

    PetscFunctionBeginUser;
    
    PROFILE_FUNCTION_BEGIN;

    // Start of function execution
    LOG_ALLOW(GLOBAL, LOG_INFO, "Entering the function.\n");

    // Validate input Pointers
    if (!user) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Input 'user' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!localBBox) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Output 'localBBox' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }

    // Get MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Get the local coordinates vector from the DMDA
    ierr = DMGetCoordinatesLocal(user->da, &coordinates);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error getting local coordinates vector.\n");
        return ierr;
    }

    if (!coordinates) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Coordinates vector is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }

    // Access the coordinate array for reading
    ierr = DMDAVecGetArrayRead(user->fda, coordinates, &coordArray);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error accessing coordinate array.\n");
        return ierr;
    }

    // Get the local grid information (indices and sizes)
    ierr = DMDAGetLocalInfo(user->da, &info);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error getting DMDA local info.\n");
        return ierr;
    }

    
    xs = info.gxs; xe = xs + info.gxm;
    ys = info.gys; ye = ys + info.gym;
    zs = info.gzs; ze = zs + info.gzm;
    
    /*
    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;
    */
    
    // Initialize min and max coordinates with extreme values
    minCoords.x = minCoords.y = minCoords.z = PETSC_MAX_REAL;
    maxCoords.x = maxCoords.y = maxCoords.z = PETSC_MIN_REAL;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Grid indices (Including Ghosts): xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n",rank, xs, xe, ys, ye, zs, ze);

    // Iterate over the local grid to find min and max coordinates
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {
                // Only consider nodes within the physical domain.
                if(i < user->IM && j < user->JM && k < user->KM){
                    Cmpnts coord = coordArray[k][j][i];

                    // Update min and max coordinates
                    if (coord.x < minCoords.x) minCoords.x = coord.x;
                    if (coord.y < minCoords.y) minCoords.y = coord.y;
                    if (coord.z < minCoords.z) minCoords.z = coord.z;

                    if (coord.x > maxCoords.x) maxCoords.x = coord.x;
                    if (coord.y > maxCoords.y) maxCoords.y = coord.y;
                    if (coord.z > maxCoords.z) maxCoords.z = coord.z;
                }    
            }
        }
    }


    // Add tolerance to bboxes.
    minCoords.x =  minCoords.x - BBOX_TOLERANCE;
    minCoords.y =  minCoords.y - BBOX_TOLERANCE;
    minCoords.z =  minCoords.z - BBOX_TOLERANCE;

    maxCoords.x =  maxCoords.x + BBOX_TOLERANCE;
    maxCoords.y =  maxCoords.y + BBOX_TOLERANCE;
    maxCoords.z =  maxCoords.z + BBOX_TOLERANCE;

    LOG_ALLOW(LOCAL,LOG_INFO," Tolerance added to the limits: %.8e .\n",(PetscReal)BBOX_TOLERANCE);
       
    // Log the computed min and max coordinates
     LOG_ALLOW(LOCAL, LOG_INFO,"[Rank %d] Bounding Box Ranges = X[%.6f, %.6f], Y[%.6f,%.6f], Z[%.6f, %.6f].\n",rank,minCoords.x, maxCoords.x,minCoords.y, maxCoords.y, minCoords.z, maxCoords.z);


    
    // Restore the coordinate array
    ierr = DMDAVecRestoreArrayRead(user->fda, coordinates, &coordArray);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error restoring coordinate array.\n");
        return ierr;
    }

    // Set the local bounding box
    localBBox->min_coords = minCoords;
    localBBox->max_coords = maxCoords;

    // Update the bounding box inside the UserCtx for consistency
    user->bbox = *localBBox;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Exiting the function successfully.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GatherAllBoundingBoxes"

/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * Each rank computes its local bounding box, then all ranks
 * participate in an MPI_Gather to send their BoundingBox to rank 0.
 * Rank 0 allocates the result array and returns it via allBBoxes.
 *
 * @param[in]  user       Pointer to UserCtx (must be non-NULL).
 * @param[out] allBBoxes  On rank 0, receives malloc’d array of size `size`.
 *                        On other ranks, set to NULL.
 * @return PetscErrorCode
 */
PetscErrorCode GatherAllBoundingBoxes(UserCtx *user, BoundingBox **allBBoxes)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank, size;
  BoundingBox    *bboxArray = NULL;
  BoundingBox     localBBox;

  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  /* Validate */
  if (!user || !allBBoxes) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                                    "GatherAllBoundingBoxes: NULL pointer");

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRMPI(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRMPI(ierr);

  /* Compute local bbox */
  ierr = ComputeLocalBoundingBox(user, &localBBox); CHKERRQ(ierr);

  /* Ensure everyone is synchronized before the gather */
  MPI_Barrier(PETSC_COMM_WORLD);
  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
    "Rank %d: about to MPI_Gather(localBBox)\n", rank);

  /* Allocate on root */
  if (rank == 0) {
    bboxArray = (BoundingBox*)malloc(size * sizeof(BoundingBox));
    if (!bboxArray) SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_MEM,
                                "GatherAllBoundingBoxes: malloc failed");
  }

  /* Collective: every rank must call */
  ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
                    bboxArray, sizeof(BoundingBox), MPI_BYTE,
                    0, PETSC_COMM_WORLD);
  CHKERRMPI(ierr);

  MPI_Barrier(PETSC_COMM_WORLD);
  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
    "Rank %d: completed MPI_Gather(localBBox)\n", rank);

  /* Return result */
  if (rank == 0) {
    *allBBoxes = bboxArray;
  } else {
    *allBBoxes = NULL;
  }

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BroadcastAllBoundingBoxes"

/**
 * @brief Broadcasts the bounding box list from rank 0 to all ranks.
 *
 * After GatherAllBoundingBoxes, rank 0 has an array of `size` boxes.
 * This routine makes sure every rank ends up with its own malloc’d copy.
 *
 * @param[in]  user      Pointer to UserCtx (unused here, but kept for signature).
 * @param[in,out] bboxlist  On entry: rank 0’s array; on exit: every rank’s array.
 * @return PetscErrorCode
 */
PetscErrorCode BroadcastAllBoundingBoxes(UserCtx *user, BoundingBox **bboxlist)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank, size;

  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  if (!bboxlist) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                          "BroadcastAllBoundingBoxes: NULL pointer");

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRMPI(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRMPI(ierr);

  /* Non-root ranks must allocate before the Bcast */
  if (rank != 0) {
    *bboxlist = (BoundingBox*)malloc(size * sizeof(BoundingBox));
    if (!*bboxlist) SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_MEM,
                                "BroadcastAllBoundingBoxes: malloc failed");
  }

  MPI_Barrier(PETSC_COMM_WORLD);
  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
    "Rank %d: about to MPI_Bcast(%d boxes)\n", rank, size);

  /* Collective: every rank must call */
  ierr = MPI_Bcast(*bboxlist, size * sizeof(BoundingBox), MPI_BYTE,
                   0, PETSC_COMM_WORLD);
  CHKERRMPI(ierr);

  MPI_Barrier(PETSC_COMM_WORLD);
  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
    "Rank %d: completed MPI_Bcast(%d boxes)\n", rank, size);


  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalculateInletCenter"

/**
 * @brief Calculates the geometric center of the primary inlet face.
 *
 * This function identifies the first face designated as an INLET in the boundary
 * condition configuration. It then iterates over all grid nodes on that physical
 * face across all MPI processes, calculates the average of their coordinates,
 * and stores the result in the user's SimCtx (CMx_c, CMy_c, CMz_c).
 *
 * This provides an automatic, robust way to determine the center for profiles
 * like parabolic flow, removing the need for manual user input.
 *
 * @param user The main UserCtx struct, containing BC config and the grid coordinate vector.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CalculateInletCenter(UserCtx *user)
{
    PetscErrorCode ierr;
    BCFace         inlet_face_id = -1;
    PetscBool      inlet_found = PETSC_FALSE;

    PetscReal      local_sum[3] = {0.0, 0.0, 0.0};
    PetscReal      global_sum[3] = {0.0, 0.0, 0.0};
    PetscCount     local_n_points = 0;
    PetscCount     global_n_points = 0;
    
    DM             da = user->da;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;
    Vec            lCoor;
    Cmpnts         ***coor;


    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    // 1. Identify the primary inlet face from the configuration
    for (int i = 0; i < 6; i++) {
        if (user->boundary_faces[i].mathematical_type == INLET) {
            inlet_face_id = user->boundary_faces[i].face_id;
            inlet_found = PETSC_TRUE;
            break; // Use the first inlet found
        }
    }

    if (!inlet_found) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "No INLET face found. Skipping inlet center calculation.\n");
        PetscFunctionReturn(0);
    }
    
    // 2. Get the nodal coordinates
    ierr = DMGetCoordinatesLocal(user->da,&lCoor);
    ierr = DMDAVecGetArrayRead(user->fda, lCoor, &coor); CHKERRQ(ierr);

    // 3. Loop over the identified inlet face and sum local coordinates
    switch (inlet_face_id) {
        case BC_FACE_NEG_X:
            if (xs == 0) {
                for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
                    if(j < user->JM && k < user->KM){ // Ensure within physical domain
                        local_sum[0] += coor[k][j][0].x;
                        local_sum[1] += coor[k][j][0].y;
                        local_sum[2] += coor[k][j][0].z;
                        local_n_points++;
                    }
                }
            }
            break;
        case BC_FACE_POS_X:
            if (xe == mx) { // another check could be if (xe > user->IM - 1)
                for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
                    if(j < user->JM && k < user->KM){ // Ensure within physical domain
                        local_sum[0] += coor[k][j][mx-2].x; // mx-1 is the ghost layer 
                        local_sum[1] += coor[k][j][mx-2].y; // mx-2 = IM - 1.
                        local_sum[2] += coor[k][j][mx-2].z;
                        local_n_points++;
                    }
                }
            }
            break;
        case BC_FACE_NEG_Y:
            if (ys == 0) {
                for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) {
                    if(i < user->IM && k < user->KM){ // Ensure within physical domain
                        local_sum[0] += coor[k][0][i].x;
                        local_sum[1] += coor[k][0][i].y;
                        local_sum[2] += coor[k][0][i].z;
                        local_n_points++;
                    }
                }
            }
            break;
        case BC_FACE_POS_Y:
            if (ye == my) { // another check could be if (ye > user->JM - 1)
                for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) {
                    local_sum[0] += coor[k][my-2][i].x; // my-1 is the ghost layer
                    local_sum[1] += coor[k][my-2][i].y; // my-2 = JM - 1.
                    local_sum[2] += coor[k][my-2][i].z;
                    local_n_points++;
                }
            }
            break;
        case BC_FACE_NEG_Z:
            if (zs == 0) {
                for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) {
                    if(i < user->IM && j < user->JM){ // Ensure within physical domain
                        local_sum[0] += coor[0][j][i].x;
                        local_sum[1] += coor[0][j][i].y;
                        local_sum[2] += coor[0][j][i].z;
                        local_n_points++;
                    }
                }
            }
            break;
        case BC_FACE_POS_Z:
            if (ze == mz) { // another check could be if (ze > user->KM - 1)
                for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) {
                    if(i < user->IM && j < user->JM){ // Ensure within physical domain
                        local_sum[0] += coor[mz-2][j][i].x; // mz-1 is the ghost layer
                        local_sum[1] += coor[mz-2][j][i].y; // mz-2 = KM - 1.
                        local_sum[2] += coor[mz-2][j][i].z;
                        local_n_points++;
                    }
                }
            }
            break;
    }
    
    ierr = DMDAVecRestoreArrayRead(user->fda, lCoor, &coor); CHKERRQ(ierr);

    // 4. Perform MPI Allreduce to get global sums
    ierr = MPI_Allreduce(local_sum, global_sum, 3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&local_n_points, &global_n_points, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    // 5. Calculate average and store in SimCtx
    if (global_n_points > 0) {
        user->simCtx->CMx_c = global_sum[0] / global_n_points;
        user->simCtx->CMy_c = global_sum[1] / global_n_points;
        user->simCtx->CMz_c = global_sum[2] / global_n_points;
        LOG_ALLOW(GLOBAL, LOG_INFO, "Calculated inlet center for Face %d: (x=%.4f, y=%.4f, z=%.4f)\n", 
                  inlet_face_id, user->simCtx->CMx_c, user->simCtx->CMy_c, user->simCtx->CMz_c);
    } else {
         LOG_ALLOW(GLOBAL, LOG_WARNING, "WARNING: Inlet face was identified but no grid points found on it. Center not calculated.\n");
    }


    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}