// logging.c
#include "logging.h"

/* Maximum temporary buffer size for converting numbers to strings */
#define TMP_BUF_SIZE 128

// --------------------- Static Variable for Log Level ---------------------

/**
 * @brief Static variable to cache the current logging level.
 *
 * Initialized to -1 to indicate that the log level has not been set yet.
 */
static LogLevel current_log_level = -1;

// --------------------- Static Variables for Allow-List -------------------

/**
 * @brief Global/static array of function names allowed to log.
 */
static char** gAllowedFunctions = NULL;

/**
 * @brief Number of entries in the gAllowedFunctions array.
 */
static int gNumAllowed = 0;




//-----------------------------  Logging Events definition (PetscLog) ---------------------

PetscLogEvent EVENT_Individualwalkingsearch = 0 ; //Individual walking search in (walkingsearch.c/LocateParticleInGrid()
PetscLogEvent EVENT_walkingsearch = 0 ; // Total walking search in (ParticleSwarm.c/LocateAllParticlesInGrid()
PetscLogEvent EVENT_GlobalParticleLocation = 0; // Global Particle Location routine (Search + Hand-Off)
PetscLogEvent EVENT_IndividualLocation = 0;

// --------------------- Function Implementations ---------------------

/**
 * @brief Retrieves the current logging level from the environment variable `LOG_LEVEL`.
 *
 * The function checks the `LOG_LEVEL` environment variable and sets the logging level accordingly.
 * Supported levels are "DEBUG", "INFO", "WARNING", and defaults to "ERROR" if not set or unrecognized.
 * The log level is cached after the first call to avoid repeated environment variable checks.
 *
 * @return LogLevel The current logging level.
 */
LogLevel get_log_level() {
    if (current_log_level == -1) { // Log level not set yet
        const char *env = getenv("LOG_LEVEL");
        if (!env) {
            current_log_level = LOG_ERROR; // Default level
        }
        else if (strcmp(env, "DEBUG") == 0) {
            current_log_level = LOG_DEBUG;
        }
        else if (strcmp(env, "INFO") == 0) {
            current_log_level = LOG_INFO;
        }
        else if (strcmp(env, "WARNING") == 0) {
            current_log_level = LOG_WARNING;
        }
        else if (strcmp(env, "PROFILE") == 0) {  // <-- New profile level
            current_log_level = LOG_PROFILE;
        }
        else {
            current_log_level = LOG_ERROR; // Default if unrecognized
        }
    }
    return current_log_level;
}

/**
 * @brief Prints the current logging level to the console.
 *
 * This function retrieves the log level using `get_log_level()` and prints 
 * the corresponding log level name. It helps verify the logging configuration 
 * at runtime.
 *
 * @note The log level is determined from the `LOG_LEVEL` environment variable.
 *       If `LOG_LEVEL` is not set, it defaults to `LOG_INFO`.
 *
 * @see get_log_level()
 */
PetscErrorCode print_log_level(void)
{
  PetscMPIInt     rank;
  PetscErrorCode  ierr;
  int             level;
  const char     *level_name;

  PetscFunctionBeginUser;
  /* get MPI rank */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRMPI(ierr);

  /* decide level name */
  level       = get_log_level();
  level_name = (level == LOG_ERROR)   ? "ERROR"   :
               (level == LOG_WARNING) ? "WARNING" :
               (level == LOG_INFO)    ? "INFO"    :
               (level == LOG_DEBUG)   ? "DEBUG"   :
               (level == LOG_PROFILE) ? "PROFILE" : "UNKNOWN";

  /* print it out */
  ierr = PetscPrintf(PETSC_COMM_SELF,
                     "Current log level: %s (%d) | rank: %d\n",
                     level_name, level, (int)rank);
  CHKERRMPI(ierr);

  PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Sets the global list of function names that are allowed to log.
 *
 * @param functionList An array of function name strings (e.g., {"foo", "bar"}).
 * @param count The number of entries in the array.
 *
 * The existing allow-list is cleared and replaced by the new one.
 * If you pass an empty list (count = 0), then no function will be allowed
 * unless you change it later.
 */
void set_allowed_functions(const char** functionList, int count)
{
    // 1. Free any existing entries
    if (gAllowedFunctions) {
        for (int i = 0; i < gNumAllowed; ++i) {
            free(gAllowedFunctions[i]); // each was strdup'ed
        }
        free(gAllowedFunctions);
        gAllowedFunctions = NULL;
        gNumAllowed = 0;
    }

    // 2. Allocate new array
    if (count > 0) {
        gAllowedFunctions = (char**)malloc(sizeof(char*) * count);
    }

    // 3. Copy the new entries
    for (int i = 0; i < count; ++i) {
        // strdup is a POSIX function. If not available, implement your own string copy.
        gAllowedFunctions[i] = strdup(functionList[i]);
    }
    gNumAllowed = count;
}

/**
 * @brief Checks if the given function name is in the allow-list.
 *
 * @param functionName The name of the function to check.
 * @return PETSC_TRUE if functionName is allowed, otherwise PETSC_FALSE.
 *
 * If no functions are in the list, nothing is allowed by default.
 * You can reverse this logic if you prefer to allow everything
 * unless specified otherwise.
 */
PetscBool is_function_allowed(const char* functionName)
{
    /* no list ⇒ allow all */
    if (gNumAllowed == 0) {
        return PETSC_TRUE;
    }

    /* otherwise only the listed functions are allowed */
    for (int i = 0; i < gNumAllowed; ++i) {
        if (strcmp(gAllowedFunctions[i], functionName) == 0) {
            return PETSC_TRUE;
        }
    }
    return PETSC_FALSE;
}

/**
 * @brief Prints the coordinates of a cell's vertices.
 *
 * This function iterates through the eight vertices of a given cell and prints their
 * coordinates. It is primarily used for debugging purposes to verify the correctness
 * of cell vertex assignments.
 *
 * @param[in]  cell     Pointer to a `Cell` structure representing the cell, containing its vertices.
 * @param[in]  rank     MPI rank for identification (useful in parallel environments).
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - Ensure that the `cell` pointer is not `NULL` before calling this function..
 */
PetscErrorCode LOG_CELL_VERTICES(const Cell *cell, PetscMPIInt rank)
{

    // Validate input pointers
    if (cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "LOG_CELL_VERTICES - 'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "LOG_CELL_VERTICES - Input parameter 'cell' is NULL.");
    }

      LOG_ALLOW(LOCAL,LOG_DEBUG, "LOG_CELL_VERTICES - Rank %d, Cell Vertices:\n", rank);
        for(int i = 0; i < 8; i++){ 
	  LOG_ALLOW(LOCAL,LOG_DEBUG, "  Vertex[%d]: (%.2f, %.2f, %.2f)\n", 
                       i, cell->vertices[i].x, cell->vertices[i].y, cell->vertices[i].z);
        }

    return 0; // Indicate successful execution
}


/**
 * @brief Prints the signed distances to each face of the cell.
 *
 * This function iterates through the six signed distances from a point to each face of a given cell
 * and prints their values. It is primarily used for debugging purposes to verify the correctness
 * of distance calculations.
 *
 * @param[in]  d        An array of six `PetscReal` values representing the signed distances.
 *                      The indices correspond to:
 *                      - d[LEFT]: Left Face
 *                      - d[RIGHT]: Right Face
 *                      - d[BOTTOM]: Bottom Face
 *                      - d[TOP]: Top Face
 *                      - d[FRONT]: Front Face
 *                      - d[BACK]: Back Face
 *
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - Ensure that the `d` array is correctly populated with signed distances before calling this function.
 */
PetscErrorCode LOG_FACE_DISTANCES(PetscReal* d)
{

    // Validate input array
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, " LOG_FACE_DISTANCES - 'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, " LOG_FACE_DISTANCES - Input array 'd' is NULL.");
    }

    PetscPrintf(PETSC_COMM_SELF, " LOG_FACE_DISTANCES - Face Distances:\n");
    PetscPrintf(PETSC_COMM_SELF, "  LEFT(%d):   %.15f\n", LEFT, d[LEFT]);
    PetscPrintf(PETSC_COMM_SELF, "  RIGHT(%d):  %.15f\n", RIGHT, d[RIGHT]);
    PetscPrintf(PETSC_COMM_SELF, "  BOTTOM(%d): %.15f\n", BOTTOM, d[BOTTOM]);
    PetscPrintf(PETSC_COMM_SELF, "  TOP(%d):    %.15f\n", TOP, d[TOP]);
    PetscPrintf(PETSC_COMM_SELF, "  FRONT(%d):  %.15f\n", FRONT, d[FRONT]);
    PetscPrintf(PETSC_COMM_SELF, "  BACK(%d):   %.15f\n", BACK, d[BACK]);

    return 0; // Indicate successful execution
}

/*
 * Helper function: Converts an integer (of type int) to a string.
 */
static void IntToStr(int value, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "%d", value);
}

/*
 * Helper function: Converts a 64‐bit integer to a string.
 */
static void Int64ToStr(PetscInt64 value, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "%ld", value);
}

/*
 * Helper function: Converts three integers into a formatted string "(i, j, k)".
 */
static void CellToStr(const PetscInt *cell, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "(%d, %d, %d)", cell[0], cell[1], cell[2]);
}

/*
 * Helper function: Converts three PetscReal values into a formatted string "(x, y, z)".
 */
static void TripleRealToStr(const PetscReal *arr, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "(%.4f, %.4f, %.4f)", arr[0], arr[1], arr[2]);
}

/*
 * Helper function: Computes the maximum string length for each column (across all particles).
 *
 * The function examines every particle (from 0 to nParticles-1) and converts the value to a
 * string using the helper functions above. The maximum length is stored in the pointers provided.
 *
 * @param nParticles       Number of particles.
 * @param ranks            Array of particle MPI ranks.
 * @param pids             Array of particle IDs.
 * @param cellIDs          Array of cell IDs (stored consecutively, 3 per particle).
 * @param positions        Array of positions (3 per particle).
 * @param velocities       Array of velocities (3 per particle).
 * @param weights          Array of weights (3 per particle).
 * @param wRank            [out] Maximum width for Rank column.
 * @param wPID             [out] Maximum width for PID column.
 * @param wCell            [out] Maximum width for Cell column.
 * @param wPos             [out] Maximum width for Position column.
 * @param wVel             [out] Maximum width for Velocity column.
 * @param wWt              [out] Maximum width for Weights column.
 */
static PetscErrorCode ComputeMaxColumnWidths(PetscInt nParticles,
                                               const PetscMPIInt *ranks,
                                               const PetscInt64  *pids,
                                               const PetscInt  *cellIDs,
                                               const PetscReal   *positions,
                                               const PetscReal   *velocities,
                                               const PetscReal   *weights,
                                               int *wRank, int *wPID, int *wCell,
                                               int *wPos, int *wVel, int *wWt)
{
    char tmp[TMP_BUF_SIZE];

    *wRank = strlen("Rank");  /* Start with the header label lengths */
    *wPID  = strlen("PID");
    *wCell = strlen("Cell (i,j,k)");
    *wPos  = strlen("Position (x,y,z)");
    *wVel  = strlen("Velocity (x,y,z)");
    *wWt   = strlen("Weights (a1,a2,a3)");

    for (PetscInt i = 0; i < nParticles; i++) {
        /* Rank */
        IntToStr(ranks[i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wRank) *wRank = (int)strlen(tmp);

        /* PID */
        Int64ToStr(pids[i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wPID) *wPID = (int)strlen(tmp);

        /* Cell: use the three consecutive values */
        CellToStr(&cellIDs[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wCell) *wCell = (int)strlen(tmp);

        /* Position */
        TripleRealToStr(&positions[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wPos) *wPos = (int)strlen(tmp);

        /* Velocity */
        TripleRealToStr(&velocities[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wVel) *wVel = (int)strlen(tmp);

        /* Weights */
        TripleRealToStr(&weights[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wWt) *wWt = (int)strlen(tmp);
    }
    return 0;
}

/*
 * Helper function: Builds a format string for a table row.
 *
 * The format string will include proper width specifiers for each column.
 * For example, it might create something like:
 *
 * "| %-6s | %-8s | %-20s | %-25s | %-25s | %-25s |\n"
 *
 * @param wRank   Maximum width for the Rank column.
 * @param wPID    Maximum width for the PID column.
 * @param wCell   Maximum width for the Cell column.
 * @param wPos    Maximum width for the Position column.
 * @param wVel    Maximum width for the Velocity column.
 * @param wWt     Maximum width for the Weights column.
 * @param fmtStr  Buffer in which to build the format string.
 * @param bufSize Size of fmtStr.
 */
static void BuildRowFormatString(PetscMPIInt wRank, PetscInt wPID, PetscInt wCell, PetscInt wPos, PetscInt wVel, PetscInt wWt, char *fmtStr, size_t bufSize)
{
    // Build a format string using snprintf.
    // We assume that the Rank is an int (%d), PID is a 64-bit int (%ld)
    // and the remaining columns are strings (which have been formatted already).
    snprintf(fmtStr, bufSize,
             "| %%-%dd | %%-%dd | %%-%ds | %%-%ds | %%-%ds | %%-%ds |\n",
             wRank, wPID, wCell, wPos, wVel, wWt);
}

/*
 * Helper function: Builds a header string for the table using column titles.
 */
static void BuildHeaderString(char *headerStr, size_t bufSize, PetscMPIInt wRank, PetscInt wPID, PetscInt wCell, PetscInt wPos, PetscInt wVel, PetscInt wWt)
{
    snprintf(headerStr, bufSize,
             "| %-*s | %-*s | %-*s | %-*s | %-*s | %-*s |\n",
             (int)wRank, "Rank",
             (int)wPID, "PID",
             (int)wCell, "Cell (i,j,k)",
             (int)wPos, "Position (x,y,z)",
             (int)wVel, "Velocity (x,y,z)",
             (int)wWt, "Weights (a1,a2,a3)");
}

/**
 * @brief Prints particle fields in a table that automatically adjusts its column widths.
 *
 * This function retrieves data from the particle swarm and prints a table where the
 * width of each column is determined by the maximum width needed to display the data.
 * Only every 'printInterval'-th particle is printed.
 *
 * @param[in] user           Pointer to the UserCtx structure.
 * @param[in] printInterval  Only every printInterval‑th particle is printed.
 *
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode LOG_PARTICLE_FIELDS(UserCtx* user, PetscInt printInterval)
{
    DM swarm = user->swarm;
    PetscErrorCode ierr;
    PetscInt localNumParticles;
    PetscReal *positions = NULL;
    PetscInt64 *particleIDs = NULL;
    PetscMPIInt *particleRanks = NULL;
    PetscInt *cellIDs = NULL;
    PetscReal *weights = NULL;
    PetscReal *velocities = NULL;
    PetscMPIInt rank;
    
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "PrintParticleFields - Rank %d is retrieving particle data.\n", rank);

    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"PrintParticleFields - Rank %d has %d particles.\n", rank, localNumParticles);

    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    
    /* Compute maximum column widths. */
    int wRank, wPID, wCell, wPos, wVel, wWt;
    wRank = wPID = wCell = wPos = wVel = wWt = 0;
    ierr = ComputeMaxColumnWidths(localNumParticles, particleRanks, particleIDs, cellIDs,
                                  positions, velocities, weights,
                                  &wRank, &wPID, &wCell, &wPos, &wVel, &wWt); CHKERRQ(ierr);

    /* Build a header string and a row format string. */
    char headerFmt[256];
    char rowFmt[256];
    BuildHeaderString(headerFmt, sizeof(headerFmt), wRank, wPID, wCell, wPos, wVel, wWt);
    BuildRowFormatString(wRank, wPID, wCell, wPos, wVel, wWt, rowFmt, sizeof(rowFmt));
    
    /* Print header (using synchronized printing for parallel output). */
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------------------------------------------\n"); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", headerFmt); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------------------------------------------\n"); CHKERRQ(ierr);
    
    /* Loop over particles and print every printInterval-th row. */
    char rowStr[256];
    for (PetscInt i = 0; i < localNumParticles; i++) {
        if (i % printInterval == 0) {
	  // ------- DEBUG 
	  //char cellStr[TMP_BUF_SIZE], posStr[TMP_BUF_SIZE], velStr[TMP_BUF_SIZE], wtStr[TMP_BUF_SIZE];
	  //CellToStr(&cellIDs[3*i], cellStr, TMP_BUF_SIZE);
	  //TripleRealToStr(&positions[3*i], posStr, TMP_BUF_SIZE);
	  //TripleRealToStr(&velocities[3*i], velStr, TMP_BUF_SIZE);
	  // TripleRealToStr(&weights[3*i], wtStr, TMP_BUF_SIZE);
            
	  // if (rank == 0) { // Or whatever rank is Rank 0
	    //PetscPrintf(PETSC_COMM_SELF, "[Rank 0 DEBUG LPF] Particle %lld: PID=%lld, Rank=%d\n", (long long)i, (long long)particleIDs[i], particleRanks[i]);
	    //PetscPrintf(PETSC_COMM_SELF, "[Rank 0 DEBUG LPF] Raw Pos: (%.10e, %.10e, %.10e)\n", positions[3*i+0], positions[3*i+1], positions[3*i+2]);
	    //PetscPrintf(PETSC_COMM_SELF, "[Rank 0 DEBUG LPF] Str Pos: %s\n", posStr);
	    //PetscPrintf(PETSC_COMM_SELF, "[Rank 0 DEBUG LPF] Raw Vel: (%.10e, %.10e, %.10e)\n", velocities[3*i+0], velocities[3*i+1], velocities[3*i+2]);
	    // PetscPrintf(PETSC_COMM_SELF, "[Rank 0 DEBUG LPF] Str Vel: %s\n", velStr);
	    // Add similar for cell, weights
	    // PetscPrintf(PETSC_COMM_SELF, "[Rank 0 DEBUG LPF] About to build rowStr for particle %lld\n", (long long)i);
	    //  fflush(stdout);
	    // }

	  //  snprintf(rowStr, sizeof(rowStr), rowFmt,
          //           particleRanks[i],
          //           particleIDs[i],
          //           cellStr,
          //           posStr,
          //           velStr,
	  //	     wtStr);

	     
	     //    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", rowStr); CHKERRQ(ierr);
	  
	     // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", rowStr); CHKERRQ(ierr); 
	  
	  // -------- DEBUG
            /* Format the row by converting each field to a string first.
             * We use temporary buffers and then build the row string.
             */
	  
	   char cellStr[TMP_BUF_SIZE], posStr[TMP_BUF_SIZE], velStr[TMP_BUF_SIZE], wtStr[TMP_BUF_SIZE];
            CellToStr(&cellIDs[3*i], cellStr, TMP_BUF_SIZE);
            TripleRealToStr(&positions[3*i], posStr, TMP_BUF_SIZE);
            TripleRealToStr(&velocities[3*i], velStr, TMP_BUF_SIZE);
            TripleRealToStr(&weights[3*i], wtStr, TMP_BUF_SIZE);
            
            /* Build the row string. Note that for the integer fields we can use the row format string. */
            snprintf(rowStr, sizeof(rowStr), rowFmt,
                     particleRanks[i],
                     particleIDs[i],
                     cellStr,
                     posStr,
                     velStr,
                     wtStr);
	 ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", rowStr); CHKERRQ(ierr);
        }
    }

 
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------------------------------------------\n"); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"PrintParticleFields - Completed printing on Rank %d.\n", rank);

    /* Restore fields */
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL,LOG_DEBUG, "PrintParticleFields - Restored all particle fields.\n");
    return 0;
}
 
 /**
  * @brief Logs the interpolation error between the analytical and computed solutions.
  *
  * This function creates global vectors for the "position" and "velocity" fields from the DMSwarm,
  * applies the analytical solution to the position vector, and then computes the L2 norm of the difference
  * between the analytical and computed solutions. The resulting interpolation error is logged.
  *
  * @param user The user context.
  * @return PetscErrorCode Returns 0 on success.
  */
 PetscErrorCode LOG_INTERPOLATION_ERROR(UserCtx *user)
 {
     PetscErrorCode ierr;
     DM swarm = user->swarm;
     Vec positionVec,analyticalvelocityVec,velocityVec, errorVec;
     PetscReal Interpolation_error = 0.0;
     PetscReal AnalyticalSolution_magnitude = 0.0;
     PetscReal ErrorPercentage = 0.0;
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Creating global vectors for 'position' and 'velocity'.\n");
     ierr = DMSwarmCreateGlobalVectorFromField(swarm, "position", &positionVec); CHKERRQ(ierr);
     ierr = DMSwarmCreateGlobalVectorFromField(swarm, "velocity", &velocityVec); CHKERRQ(ierr);
 
     ierr = VecDuplicate(positionVec, &analyticalvelocityVec); CHKERRQ(ierr);
     ierr = VecCopy(positionVec, analyticalvelocityVec); CHKERRQ(ierr);
     
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Applying analytical solution to position vector.\n");
     ierr = SetAnalyticalSolution(analyticalvelocityVec,user->FieldInitialization); CHKERRQ(ierr);

     ierr = VecDuplicate(analyticalvelocityVec, &errorVec); CHKERRQ(ierr);
     ierr = VecCopy(analyticalvelocityVec, errorVec); CHKERRQ(ierr);
     
     // Get the magnitude of analytical solution
     ierr = VecNorm(analyticalvelocityVec,NORM_2,&AnalyticalSolution_magnitude); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Computing difference between velocity and analytical solution.\n");
     ierr = VecAXPY(errorVec, -1.0, velocityVec); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Calculating L2 norm of the difference.\n");
     ierr = VecNorm(errorVec, NORM_2, &Interpolation_error); CHKERRQ(ierr);

     ErrorPercentage = Interpolation_error/AnalyticalSolution_magnitude;
     
     ErrorPercentage = ErrorPercentage*100;
     
     LOG_ALLOW(GLOBAL, LOG_INFO, "LOG_INTERPOLATION_ERROR - Interpolation error (%%): %g\n", ErrorPercentage);
 
     PetscPrintf(PETSC_COMM_WORLD, "Interpolation error (%%): %g\n", ErrorPercentage);

     ierr = VecDestroy(&errorVec); CHKERRQ(ierr);
     ierr = DMSwarmDestroyGlobalVectorFromField(swarm, "position", &positionVec); CHKERRQ(ierr);
     ierr = DMSwarmDestroyGlobalVectorFromField(swarm, "velocity", &velocityVec); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Completed logging interpolation error.\n");
     return 0;
 }



static void trim(char *s)
{
    if (!s) return;

    /* ---- 1. strip leading blanks ----------------------------------- */
    char *p = s;
    while (*p && isspace((unsigned char)*p))
        ++p;

    if (p != s)                      /* move the trimmed text forward   */
        memmove(s, p, strlen(p) + 1);   /* +1 to copy the final NUL     */

    /* ---- 2. strip trailing blanks ---------------------------------- */
    size_t len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1]))
        s[--len] = '\0';
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Load function names from a text file.
 *
 * The file is expected to contain **one identifier per line**.  Blank lines and
 * lines whose first non‑blank character is a <tt>#</tt> are silently skipped so
 * the file can include comments.  Example:
 *
 * @code{.txt}
 * # Allowed function list
 * main
 * InitializeSimulation
 * InterpolateAllFieldsToSwarm  # inline comments are OK, too
 * @endcode
 *
 * The routine allocates memory as needed (growing an internal buffer with
 * @c PetscRealloc()) and returns the resulting array and its length to the
 * caller.  Use FreeAllowedFunctions() to clean up when done.
 *
 * @param[in]  filename  Path of the configuration file to read.
 * @param[out] funcsOut  On success, points to a freshly‑allocated array of
 *                       <tt>char*</tt> (size @p nOut).
 * @param[out] nOut      Number of valid entries in @p funcsOut.
 *
 * @return 0 on success, or a PETSc error code on failure (e.g. I/O error, OOM).
 */
PetscErrorCode LoadAllowedFunctionsFromFile(const char   filename[],
                                            char      ***funcsOut,
                                            PetscInt   *nOut)
{
  FILE          *fp    = NULL;
  char         **funcs = NULL;
  size_t         cap   = 16;   /* initial capacity */
  size_t         n     = 0;    /* number of names  */
  char           line[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* ---------------------------------------------------------------------- */
  /* 1. Open file                                                           */
  fp = fopen(filename, "r");
  if (!fp) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                    "Cannot open %s", filename);

  /* 2. Allocate initial pointer array                                      */
  ierr = PetscMalloc1(cap, &funcs); CHKERRQ(ierr);

  /* 3. Read file line by line                                              */
  while (fgets(line, sizeof line, fp)) {
    /* Strip everything after a comment character '#'. */
    char *hash = strchr(line, '#');
    if (hash) *hash = '\0';

    trim(line);                 /* remove leading/trailing blanks */
    if (!*line) continue;       /* skip if empty                  */

    /* Grow the array if necessary */
    if (n == cap) {
      cap *= 2;
      ierr = PetscRealloc(cap * sizeof(*funcs), (void **)&funcs); CHKERRQ(ierr);
    }

    /* Deep‑copy the cleaned identifier */
    ierr = PetscStrallocpy(line, &funcs[n++]); CHKERRQ(ierr);
  }
  fclose(fp);

  /* 4. Return results to caller                                           */
  *funcsOut = funcs;
  *nOut     = (PetscInt)n;

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Free an array previously returned by LoadAllowedFunctionsFromFile().
 *
 * @param[in,out] funcs Array of strings to release (may be @c NULL).
 * @param[in]     n     Number of entries in @p funcs.  Ignored if @p funcs is
 *                      @c NULL.
 *
 * @return 0 on success or a PETSc error code.
 */
PetscErrorCode FreeAllowedFunctions(char **funcs, PetscInt n)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (funcs) {
    for (PetscInt i = 0; i < n; ++i) {
      ierr = PetscFree(funcs[i]); CHKERRQ(ierr);
    }
    ierr = PetscFree(funcs); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/**
 * @brief Helper function to convert BCFace enum to a string representation.
 * @param[in] face The BCFace enum value.
 * @return Pointer to a constant string representing the face.
 */
const char* BCFaceToString(BCFace face) {
    switch (face) {
        case BC_FACE_NEG_X: return "-Xi (I-Min)";
        case BC_FACE_POS_X: return "+Xi (I-Max)";
        case BC_FACE_NEG_Y: return "-Eta (J-Min)";
        case BC_FACE_POS_Y: return "+Eta (J-Max)";
        case BC_FACE_NEG_Z: return "-Zeta (K-Min)";
        case BC_FACE_POS_Z: return "+Zeta (K-Max)";
        default:            return "Unknown Face";
    }
}

/**
 * @brief Helper function to convert BCType enum to a string representation.
 * @param[in] type The BCType enum value.
 * @return Pointer to a constant string representing the BC type.
 */
const char* BCTypeToString(BCType type) {
    switch (type) {
      //  case DIRICHLET: return "DIRICHLET";
      //  case NEUMANN:   return "NEUMANN";
        case WALL:      return "WALL";
        case INLET:     return "INLET";
        case OUTLET:    return "OUTLET";
        case FARFIELD:  return "FARFIELD";
        case PERIODIC:  return "PERIODIC";
        case INTERFACE: return "INTERFACE";
        case NOGRAD:    return "NOGRAD";

	// case CUSTOM:    return "CUSTOM";
        default:        return "Unknown BC Type";
    }
}

/**
 * @brief Converts a BCHandlerType enum to its string representation.
 *
 * Provides a descriptive string for a specific boundary condition implementation strategy.
 * This is crucial for logging the exact behavior configured for a face.
 *
 * @param handler_type The BCHandlerType enum value (e.g., BC_HANDLER_WALL_NOSLIP).
 * @return A constant character string corresponding to the enum. Returns
 *         "UNKNOWN_HANDLER" if the enum value is not recognized.
 */
const char* BCHandlerTypeToString(BCHandlerType handler_type) {
    switch (handler_type) {
        // Wall & Symmetry Handlers
        case BC_HANDLER_WALL_NOSLIP:             return "noslip";
        case BC_HANDLER_WALL_MOVING:             return "moving";
        case BC_HANDLER_SYMMETRY_PLANE:          return "symmetry_plane";

        // Inlet Handlers
        case BC_HANDLER_INLET_CONSTANT_VELOCITY: return "constant_velocity";
        case BC_HANDLER_INLET_PULSANTILE_FLUX:   return "pulsatile_flux";
        case BC_HANDLER_INLET_PARABOLIC: return "parabolic";

        // Outlet Handlers
        case BC_HANDLER_OUTLET_CONSERVATION:     return "conservation";
        case BC_HANDLER_OUTLET_PRESSURE:         return "pressure";

        // Other Physical Handlers
        case BC_HANDLER_FARFIELD_NONREFLECTING:  return "nonreflecting";
        case BC_HANDLER_NOGRAD_COPY_GHOST: return "no_gradient";

        // Multi-Block / Interface Handlers
        case BC_HANDLER_PERIODIC:                return "periodic";
        case BC_HANDLER_INTERFACE_OVERSET:       return "overset";

        // Default case
        case BC_HANDLER_UNDEFINED:
        default:                                 return "UNKNOWN_HANDLER";
    }
}
