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
        else if (strcmp(env, "VERBOSE") == 0) {
            current_log_level = LOG_VERBOSE;
        }
        else if (strcmp(env, "TRACE") == 0) {
            current_log_level = LOG_TRACE;
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
      LOG_ALLOW(LOCAL,LOG_ERROR, "'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "LOG_CELL_VERTICES - Input parameter 'cell' is NULL.");
    }

      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Rank %d, Cell Vertices:\n", rank);
        for(int i = 0; i < 8; i++){ 
	  LOG_ALLOW(LOCAL,LOG_VERBOSE, "  Vertex[%d]: (%.2f, %.2f, %.2f)\n", 
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
      LOG_ALLOW(LOCAL,LOG_ERROR, "  'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "  Input array 'd' is NULL.");
    }

    PetscPrintf(PETSC_COMM_SELF, "  Face Distances:\n");
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
    LOG_ALLOW(LOCAL,LOG_INFO, "Rank %d is retrieving particle data.\n", rank);

    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"Rank %d has %d particles.\n", rank, localNumParticles);

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

    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"Completed printing on Rank %d.\n", rank);

    /* Restore fields */
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL,LOG_DEBUG, "Restored all particle fields.\n");
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
 * @brief Helper function to convert FieldInitialization to a string representation.
 * @param[in] PetscInt The FieldInitialization value.
 * @return Pointer to a constant string representing the FieldInitialization.
 */
const char* FieldInitializationToString(PetscInt FieldInitialization)
{
    switch(FieldInitialization){
        case 0: return "Zero";
        case 1: return "Constant Normal velocity";
        case 2: return "Poiseuille Normal velocity";
        default: return "Unknown Field Initialization";
    }
}

/**
 * @brief Helper function to convert ParticleInitialization to a string representation.
 * @param[in] PetscInt The ParticleInitialization value.
 * @return Pointer to a constant string representing the FieldInitialization.
 */
const char* ParticleInitializationToString(PetscInt ParticleInitialization)
{
    switch(ParticleInitialization){
        case 0: return "Surface: Random";
        case 1: return "Volume";
        case 3: return "Surface: At edges";
        default: return "Unknown Particle Initialization";
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
        case BC_HANDLER_INLET_PULSATILE_FLUX:   return "pulsatile_flux";
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

/**
 * @brief Destroys the DualMonitorCtx.
 *
 * This function is passed to KSPMonitorSet to ensure the viewer is
 * properly destroyed and the context memory is freed when the KSP is destroyed.
 * @param Ctx a pointer to the context pointer to be destroyed
 * @return PetscErrorCode
 */
PetscErrorCode DualMonitorDestroy(void **ctx)
{
    DualMonitorCtx *monctx = (DualMonitorCtx*)*ctx;
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    
    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
    if(!rank && monctx->file_handle){
      fclose(monctx->file_handle);
    }
    
    ierr = PetscFree(monctx); CHKERRQ(ierr);
    *ctx = NULL;
    PetscFunctionReturn(0);
}

/**
 * @brief A custom KSP monitor that logs the true residual to a file and optionally to the console.
 *
 * This function replicates the behavior of KSPMonitorTrueResidualNorm by calculating
 * the true residual norm ||b - Ax|| itself. It unconditionally logs to a file
 * viewer and conditionally logs to the console based on a flag in the context.
 *
 * @param ksp   The Krylov subspace context.
 * @param it    The current iteration number.
 * @param rnorm The preconditioned residual norm (ignored, we compute our own).
 * @param ctx   A pointer to the DualMonitorCtx structure.
 * @return      PetscErrorCode 0 on success.
 */
#undef __FUNCT__
#define __FUNCT__ "DualKSPMonitor"
PetscErrorCode DualKSPMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void *ctx)
{
    DualMonitorCtx *monctx = (DualMonitorCtx*)ctx;
    PetscErrorCode ierr;
    PetscReal      trnorm, relnorm;
    Vec            r;
    char           norm_buf[256];
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

    // 1. Calculate the true residual norm.
    ierr = KSPBuildResidual(ksp, NULL, NULL, &r); CHKERRQ(ierr);
    ierr = VecNorm(r, NORM_2, &trnorm); CHKERRQ(ierr);

    // 2. On the first iteration, compute and store the norm of the RHS vector `b`.
    if (it == 0) {
        Vec b;
        ierr = KSPGetRhs(ksp, &b); CHKERRQ(ierr);
        ierr = VecNorm(b, NORM_2, &monctx->bnorm); CHKERRQ(ierr);
    }

    if(!rank){
    // 3. Compute the relative norm and format the output string.
    if (monctx->bnorm > 1.e-15) {
      relnorm = trnorm / monctx->bnorm;
      sprintf(norm_buf, "ts: %-5d | block: %-2d | iter: %-3d | Unprecond Norm: %12.5e | True Norm: %12.5e | Rel Norm: %12.5e",(int)monctx->step, (int)monctx->block_id, (int)it, (double)rnorm, (double)trnorm, (double)relnorm);
    } else {
      sprintf(norm_buf,"ts: %-5d | block: %-2d | iter: %-3d | Unprecond Norm: %12.5e | True Norm: %12.5e",(int)monctx->step, (int)monctx->block_id, (int)it, (double)rnorm, (double)trnorm);
    }

    // 4. Log to the file viewer (unconditionally).
    if(monctx->file_handle){
      ierr = PetscFPrintf(PETSC_COMM_SELF,monctx->file_handle,"%s\n", norm_buf); CHKERRQ(ierr);
    }
    // 5. Log to the console (conditionally).
    if (monctx->log_to_console) {
      PetscFPrintf(PETSC_COMM_SELF,stdout, "%s\n", norm_buf); CHKERRQ(ierr);
    }

    } //rank

    PetscFunctionReturn(0);
}

/**
 * @brief Logs continuity metrics for a single block to a file.
 *
 * This function should be called for each block, once per timestep. It opens a
 * central log file in append mode. To ensure the header is written only once,
 * it checks if it is processing block 0 on the simulation's start step.
 *
 * @param user A pointer to the UserCtx for the specific block whose metrics
 *             are to be logged. The function accesses both global (SimCtx)
 *             and local (user->...) data.
 * @return     PetscErrorCode 0 on success.
 */
#undef __FUNCT__
#define __FUNCT__ "LOG_CONTINUITY_METRICS"
PetscErrorCode LOG_CONTINUITY_METRICS(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    SimCtx         *simCtx = user->simCtx;      // Get the shared SimCtx
    const PetscInt bi = user->_this;          // Get this block's specific ID
    const PetscInt ti = simCtx->step;         // Get the current timestep

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Only rank 0 performs file I/O.
    if (!rank) {
        FILE *f;
        char filen[128];
        sprintf(filen, "%s/Continuity_Metrics.log",simCtx->log_dir);

        // Open the log file in append mode.
        f = fopen(filen, "a");
        if (!f) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open log file: %s", filen);
        }

        // Write a header only for the very first block (bi=0) on the very
        // first timestep (ti=StartStep + 1). This ensures it's written only once.
        if (ti == simCtx->StartStep + 1 && bi == 0) {
            PetscFPrintf(PETSC_COMM_SELF, f, "%-10s | %-6s | %-18s | %-30s | %-18s | %-18s | %-18s | %-18s\n",
                         "Timestep", "Block", "Max Divergence", "Max Divergence Location ([k][j][i]=idx)", "Sum(RHS)","Total Flux In", "Total Flux Out", "Net Flux");
            PetscFPrintf(PETSC_COMM_SELF, f, "------------------------------------------------------------------------------------------------------------------------------------------\n");
        }

        // Prepare the data strings and values for the current block.
        PetscReal net_flux = simCtx->FluxInSum - simCtx->FluxOutSum;
        char location_str[64];
        sprintf(location_str, "([%d][%d][%d] = %d)", (int)simCtx->MaxDivz, (int)simCtx->MaxDivy, (int)simCtx->MaxDivx, (int)simCtx->MaxDivFlatArg);

        // Write the formatted line for the current block.
        PetscFPrintf(PETSC_COMM_SELF, f, "%-10d | %-6d | %-18.10e | %-39s | %-18.10e | %-18.10e | %-18.10e | %-18.10e\n",
                     (int)ti,
                     (int)bi,
                     (double)simCtx->MaxDiv,
                     location_str,
                     (double)simCtx->summationRHS,
		     (double)simCtx->FluxInSum,
                     (double)simCtx->FluxOutSum,
                     (double)net_flux);

        fclose(f);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief A function that outputs the name of the current level in the ParticleLocation enum.
 * @param level The ParticleLocation enum value.
 * @return A constant character string corresponding to the enum. Returns
 *        "UNKNOWN_LEVEL" if the enum value is not recognized.
 */
const char* ParticleLocationStatusToString(ParticleLocationStatus level) 
{
    switch (level) {
        case NEEDS_LOCATION:        return "NEEDS_LOCATION";
        case ACTIVE_AND_LOCATED: return "ACTIVE_AND_LOCATED";
        case MIGRATING_OUT:     return "MIGRATING_OUT";
        case LOST:               return "LOST";
        case UNINITIALIZED:      return "UNINITIALIZED";
        default:            return "UNKNOWN_LEVEL";
    }
}

///////// Profiling System /////////

// Data structure to hold profiling info for one function
typedef struct {
    const char *name;
    double     total_time;
    double     current_step_time;
    long long  total_call_count;
    long long  current_step_call_count;
    double     start_time; // Timer for the current call
    PetscBool  always_log;
} ProfiledFunction;

// Global registry for all profiled functions
static ProfiledFunction *g_profiler_registry = NULL;
static PetscInt g_profiler_count = 0;
static PetscInt g_profiler_capacity = 0;

// Internal helper to find a function in the registry or create it
static PetscErrorCode _FindOrCreateEntry(const char *func_name, PetscInt *idx)
{
    PetscFunctionBeginUser;
    // Search for existing entry
    for (PetscInt i = 0; i < g_profiler_count; ++i) {
        if (strcmp(g_profiler_registry[i].name, func_name) == 0) {
            *idx = i;
            PetscFunctionReturn(0);
        }
    }

    // Not found, create a new entry
    if (g_profiler_count >= g_profiler_capacity) {
        PetscInt new_capacity = g_profiler_capacity == 0 ? 16 : g_profiler_capacity * 2;
        PetscErrorCode ierr = PetscRealloc(sizeof(ProfiledFunction) * new_capacity, &g_profiler_registry); CHKERRQ(ierr);
        g_profiler_capacity = new_capacity;
    }

    *idx = g_profiler_count;
    g_profiler_registry[*idx].name = func_name;
    g_profiler_registry[*idx].total_time = 0.0;
    g_profiler_registry[*idx].current_step_time = 0.0;
    g_profiler_registry[*idx].total_call_count = 0;
    g_profiler_registry[*idx].current_step_call_count = 0;
    g_profiler_registry[*idx].start_time = 0.0;
    g_profiler_registry[*idx].always_log = PETSC_FALSE;
    g_profiler_count++;

    PetscFunctionReturn(0);
}

// --- Public API Implementation ---
/**
 * @brief Initializes the custom profiling system using configuration from SimCtx.
 *
 * This function sets up the internal data structures for tracking function
 * performance. It reads the list of "critical functions" from the provided
 * SimCtx and marks them for per-step logging at LOG_INFO level.
 *
 * It should be called once at the beginning of the application, after
 * CreateSimulationContext() but before the main time loop.
 *
 * @param simCtx The master simulation context, which contains the list of
 *               critical function names to always log.
 * @return PetscErrorCode
 */
PetscErrorCode ProfilingInitialize(SimCtx *simCtx)
{
    PetscFunctionBeginUser;
    if (!simCtx) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be null for ProfilingInitialize");

    // Iterate through the list of critical functions provided in SimCtx
    for (PetscInt i = 0; i < simCtx->nCriticalFuncs; ++i) {
        PetscInt idx;
        const char *func_name = simCtx->criticalFuncs[i];
        PetscErrorCode ierr = _FindOrCreateEntry(func_name, &idx); CHKERRQ(ierr);
        g_profiler_registry[idx].always_log = PETSC_TRUE;

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Marked '%s' as a critical function for profiling.\n", func_name);
    }
    PetscFunctionReturn(0);
}

void _ProfilingStart(const char *func_name)
{
    PetscInt idx;
    if (_FindOrCreateEntry(func_name, &idx) != 0) return; // Fail silently
    PetscTime(&g_profiler_registry[idx].start_time);
}

void _ProfilingEnd(const char *func_name)
{
    double end_time;
    PetscTime(&end_time);

    PetscInt idx;
    if (_FindOrCreateEntry(func_name, &idx) != 0) return; // Fail silently

    double elapsed = end_time - g_profiler_registry[idx].start_time;
    g_profiler_registry[idx].total_time += elapsed;
    g_profiler_registry[idx].current_step_time += elapsed;
    g_profiler_registry[idx].total_call_count++;
    g_profiler_registry[idx].current_step_call_count++;
}

PetscErrorCode ProfilingResetTimestepCounters(void)
{
    PetscFunctionBeginUser;
    for (PetscInt i = 0; i < g_profiler_count; ++i) {
        g_profiler_registry[i].current_step_time = 0.0;
        g_profiler_registry[i].current_step_call_count = 0;
    }
    PetscFunctionReturn(0);
}

PetscErrorCode ProfilingLogTimestepSummary(PetscInt step)
{
    LogLevel log_level = get_log_level();
    PetscBool should_print = PETSC_FALSE;

    PetscFunctionBeginUser;

    // Decide if we should print anything at all
    if (log_level >= LOG_PROFILE) {
        for (PetscInt i = 0; i < g_profiler_count; ++i) {
            if (g_profiler_registry[i].current_step_call_count > 0) {
                 if (log_level == LOG_PROFILE || g_profiler_registry[i].always_log) {
                    should_print = PETSC_TRUE;
                    break;
                 }
            }
        }
    }
    
    if (should_print) {
        PetscPrintf(PETSC_COMM_SELF, "[PROFILE] ----- Timestep %d Summary -----\n", step);
        for (PetscInt i = 0; i < g_profiler_count; ++i) {
            if (g_profiler_registry[i].current_step_call_count > 0) {
                if (log_level == LOG_PROFILE || g_profiler_registry[i].always_log) {
                    PetscPrintf(PETSC_COMM_SELF, "[PROFILE]   %-25s: %.6f s (%lld calls)\n",
                                g_profiler_registry[i].name,
                                g_profiler_registry[i].current_step_time,
                                g_profiler_registry[i].current_step_call_count);
                }
            }
        }
    }

    // Reset per-step counters for the next iteration
    for (PetscInt i = 0; i < g_profiler_count; ++i) {
        g_profiler_registry[i].current_step_time = 0.0;
        g_profiler_registry[i].current_step_call_count = 0;
    }
    PetscFunctionReturn(0);
}

// Comparison function for qsort to sort by total_time in descending order
static int _CompareProfiledFunctions(const void *a, const void *b)
{
    const ProfiledFunction *func_a = (const ProfiledFunction *)a;
    const ProfiledFunction *func_b = (const ProfiledFunction *)b;

    if (func_a->total_time < func_b->total_time) return 1;
    if (func_a->total_time > func_b->total_time) return -1;
    return 0;
}

/**
 * @brief the profiling excercise and build a profiling summary which is then printed to a log file.
 * 
 * @param simCtx  The Simulation Context Structure that can contains all the data regarding the simulation.
 * 
 * @return        PetscErrorCode 0 on success.
 */
PetscErrorCode ProfilingFinalize(SimCtx *simCtx)
{
    PetscInt rank = simCtx->rank;
    PetscFunctionBeginUser;
    if (!rank) {

        //--- Step 0: Create a file viewer for log file
        FILE *f;
        char filen[128];
        sprintf(filen, "%s/ProfilingSummary.log",simCtx->log_dir);

        // Open the log file in write mode.
        f = fopen(filen,"w");
        if(!f){
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot Open log file: %s",filen);
        }
       
        // --- Step 1: Sort the data for readability ---
        qsort(g_profiler_registry, g_profiler_count, sizeof(ProfiledFunction), _CompareProfiledFunctions);

        // --- Step 2: Dynamically determine the width for the function name column ---
        PetscInt max_name_len = strlen("Function"); // Start with the header's length
        for (PetscInt i = 0; i < g_profiler_count; ++i) {
            if (g_profiler_registry[i].total_call_count > 0) {
                PetscInt len = strlen(g_profiler_registry[i].name);
                if (len > max_name_len) {
                    max_name_len = len;
                }
            }
        }
        // Add a little padding
        max_name_len += 2; 

        // --- Step 3: Define fixed widths for numeric columns for consistent alignment ---
        const int time_width = 18;
        const int count_width = 15;
        const int avg_width = 22;

        // --- Step 4: Print the formatted table ---
        PetscFPrintf(PETSC_COMM_SELF, f, "=================================================================================================================\n");
        PetscFPrintf(PETSC_COMM_SELF, f, "                         FINAL PROFILING SUMMARY (Sorted by Total Time)\n");
        PetscFPrintf(PETSC_COMM_SELF, f, "=================================================================================================================\n");
        
        // Header Row
        PetscFPrintf(PETSC_COMM_SELF, f, "%-*s | %-*s | %-*s | %-*s\n",
                    max_name_len, "Function",
                    time_width, "Total Time (s)",
                    count_width, "Call Count",
                    avg_width, "Avg. Time/Call (ms)");
        
        // Separator Line (dynamically sized)
        for (int i = 0; i < max_name_len; i++) PetscFPrintf(PETSC_COMM_SELF, f, "-");
        PetscFPrintf(PETSC_COMM_SELF, f, "-|-");
        for (int i = 0; i < time_width; i++) PetscFPrintf(PETSC_COMM_SELF, f, "-");
        PetscFPrintf(PETSC_COMM_SELF, f, "-|-");
        for (int i = 0; i < count_width; i++) PetscFPrintf(PETSC_COMM_SELF, f, "-");
        PetscFPrintf(PETSC_COMM_SELF, f, "-|-");
        for (int i = 0; i < avg_width; i++) PetscFPrintf(PETSC_COMM_SELF, f, "-");
        PetscFPrintf(PETSC_COMM_SELF, f, "\n");

        // Data Rows
        for (PetscInt i = 0; i < g_profiler_count; ++i) {
            if (g_profiler_registry[i].total_call_count > 0) {
                double avg_time_ms = (g_profiler_registry[i].total_time / g_profiler_registry[i].total_call_count) * 1000.0;
                PetscFPrintf(PETSC_COMM_SELF, f, "%-*s | %*.*f | %*lld | %*.*f\n",
                            max_name_len, g_profiler_registry[i].name,
                            time_width, 6, g_profiler_registry[i].total_time,
                            count_width, g_profiler_registry[i].total_call_count,
                            avg_width, 6, avg_time_ms);
                PetscFPrintf(PETSC_COMM_SELF, f, "------------------------------------------------------------------------------------------------------------------\n");
            }
        }
        PetscFPrintf(PETSC_COMM_SELF, f, "==================================================================================================================\n");
    
        fclose(f);
    }

    // --- Final Cleanup ---
    PetscFree(g_profiler_registry);
    g_profiler_registry = NULL;
    g_profiler_count = 0;
    g_profiler_capacity = 0;
    PetscFunctionReturn(0);
}

/*================================================================================*
 *                          PROGRESS BAR UTILITY                                  *
 *================================================================================*/

/**
 * @brief Prints a progress bar to the console.
 *
 * This function should only be called by the root process (rank 0). It uses
 * a carriage return `\r` to overwrite the same line in the terminal, creating
 * a dynamic progress bar.
 *
 * @param step           The current step index from the loop (e.g., from 0 to N-1).
 * @param startStep      The global starting step number of the simulation.
 * @param totalSteps     The total number of steps to be run in this simulation instance.
 * @param currentTime    The current simulation time to display.
 */
void PrintProgressBar(PetscInt step, PetscInt startStep, PetscInt totalSteps, PetscReal currentTime)
{
    if (totalSteps <= 0) return;

    // --- Configuration ---
    const int barWidth = 50;

    // --- Calculation ---
    // Calculate progress as a fraction from 0.0 to 1.0
    PetscReal progress = (PetscReal)(step - startStep + 1) / totalSteps;
    // Ensure progress doesn't exceed 1.0 due to floating point inaccuracies
    if (progress > 1.0) progress = 1.0;

    int pos = (int)(barWidth * progress);

    // --- Printing ---
    // Carriage return moves cursor to the beginning of the line
    PetscPrintf(PETSC_COMM_SELF, "\rProgress: [");

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) {
            PetscPrintf(PETSC_COMM_SELF, "=");
        } else if (i == pos) {
            PetscPrintf(PETSC_COMM_SELF, ">");
        } else {
            PetscPrintf(PETSC_COMM_SELF, " ");
        }
    }

    // Print percentage, step count, and current time
    PetscPrintf(PETSC_COMM_SELF, "] %3d%% (Step %ld/%ld, t=%.4f)",
                (int)(progress * 100.0),
                step + 1,
                startStep + totalSteps,
                currentTime);

    // Flush the output buffer to ensure the bar is displayed immediately
    fflush(stdout);
}

#undef __FUNCT__
#define __FUNCT__ "LOG_FIELD_MIN_MAX"
/**
 * @brief Computes and logs the local and global min/max values of a specified field,
 *        respecting the solver's data layout architecture.
 *
 * This utility function inspects a PETSc Vec and calculates the minimum and maximum
 * values for its components, both locally and globally.
 *
 * It is "architecture-aware" and adjusts its iteration range to only include
 * physically meaningful data points:
 *
 * - **Cell-Centered Fields ("Ucat", "P"):** It uses the "Shifted Index Architecture,"
 *   iterating from index 1 to N-1 to exclude the ghost/tool values at indices 0 and N.
 * - **Node-Centered Fields ("Coordinates"):** It iterates from index 0 to N-1, covering all
 *   physical nodes.
 * - **Face-Centered Fields ("Ucont"):** It iterates from index 0 to N-1, covering all
 *   physical faces.
 *
 * The results are printed to the standard output in a formatted, easy-to-read table.
 *
 * @param[in] user      Pointer to the user-defined context.
 * @param[in] fieldName A string descriptor for the field being analyzed.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode LOG_FIELD_MIN_MAX(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscInt       i, j, k;
    DMDALocalInfo  info;
    
    Vec            fieldVec = NULL;
    DM             dm = NULL;
    PetscInt       dof;
    char           data_layout[20];

    PetscFunctionBeginUser;

    // --- 1. Map string name to PETSc objects and determine data layout ---
    if (strcasecmp(fieldName, "Ucat") == 0) {
        fieldVec = user->Ucat; dm = user->fda; dof = 3; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(fieldName, "P") == 0) {
        fieldVec = user->P; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(fieldName, "Ucont") == 0) {
        fieldVec = user->lUcont; dm = user->fda; dof = 3; strcpy(data_layout, "Face-Centered");
    } else if (strcasecmp(fieldName, "Coordinates") == 0) {
        ierr = DMGetCoordinates(user->da, &fieldVec); CHKERRQ(ierr);
        dm = user->fda; dof = 3; strcpy(data_layout, "Node-Centered");
    } else if (strcasecmp(fieldName,"Psi") == 0) {
        fieldVec = user->Psi; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered"); // Assuming Psi is cell-centered
    } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Field %s not recognized.", fieldName);
    }

    if (!fieldVec) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Vector for field '%s' is NULL.", fieldName);
    }
    if (!dm) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "DM for field '%s' is NULL.", fieldName);
    }

    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);

    // --- 2. Define Architecture-Aware Loop Bounds ---
    PetscInt i_start, i_end, j_start, j_end, k_start, k_end;

    if (strcmp(data_layout, "Cell-Centered") == 0) {
        // For cell-centered data, the physical values are stored from index 1 to N-1.
        // We find the intersection of the rank's owned range [xs, xe) with the
        // physical data range [1, IM-1).
        i_start = PetscMax(info.xs, 1);          i_end = PetscMin(info.xs + info.xm, user->IM);
        j_start = PetscMax(info.ys, 1);          j_end = PetscMin(info.ys + info.ym, user->JM);
        k_start = PetscMax(info.zs, 1);          k_end = PetscMin(info.zs + info.zm, user->KM);
    } else { // For Node- or Face-Centered data
        // The physical values are stored from index 0 to N-1.
        // We find the intersection of the rank's owned range [xs, xe) with the
        // physical data range [0, IM-1].
        i_start = PetscMax(info.xs, 0);          i_end = PetscMin(info.xs + info.xm, user->IM);
        j_start = PetscMax(info.ys, 0);          j_end = PetscMin(info.ys + info.ym, user->JM);
        k_start = PetscMax(info.zs, 0);          k_end = PetscMin(info.zs + info.zm, user->KM);
    }

    // --- 3. Barrier for clean, grouped output ---
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
    if (user->simCtx->rank == 0) {
        PetscPrintf(PETSC_COMM_SELF, "\n--- Field Ranges: [%s] (Layout: %s) ---\n", fieldName, data_layout);
    }

    // --- 4. Branch on DoF and perform calculation with correct bounds ---
    if (dof == 1) {
        PetscReal    localMin = PETSC_MAX_REAL, localMax = PETSC_MIN_REAL;
        PetscReal    globalMin, globalMax;
        const PetscScalar  ***array;

        ierr = DMDAVecGetArrayRead(dm, fieldVec, &array); CHKERRQ(ierr);
        for (k = k_start; k < k_end; k++) {
            for (j = j_start; j < j_end; j++) {
                for (i = i_start; i < i_end; i++) {
                    localMin = PetscMin(localMin, array[k][j][i]);
                    localMax = PetscMax(localMax, array[k][j][i]);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm, fieldVec, &array); CHKERRQ(ierr);

        ierr = MPI_Allreduce(&localMin, &globalMin, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Allreduce(&localMax, &globalMax, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  [Rank %d] Local Range:  [ %11.4e , %11.4e ]\n", user->simCtx->rank, localMin, localMax);
        ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
        if (user->simCtx->rank == 0) {
           PetscPrintf(PETSC_COMM_SELF, "  Global Range:         [ %11.4e , %11.4e ]\n", globalMin, globalMax);
        }

    } else if (dof == 3) {
        Cmpnts localMin = {PETSC_MAX_REAL, PETSC_MAX_REAL, PETSC_MAX_REAL};
        Cmpnts localMax = {PETSC_MIN_REAL, PETSC_MIN_REAL, PETSC_MIN_REAL};
        Cmpnts globalMin, globalMax;
        const Cmpnts ***array;

        ierr = DMDAVecGetArrayRead(dm, fieldVec, &array); CHKERRQ(ierr);
        for (k = k_start; k < k_end; k++) {
            for (j = j_start; j < j_end; j++) {
                for (i = i_start; i < i_end; i++) {
                    localMin.x = PetscMin(localMin.x, array[k][j][i].x);
                    localMin.y = PetscMin(localMin.y, array[k][j][i].y);
                    localMin.z = PetscMin(localMin.z, array[k][j][i].z);
                    localMax.x = PetscMax(localMax.x, array[k][j][i].x);
                    localMax.y = PetscMax(localMax.y, array[k][j][i].y);
                    localMax.z = PetscMax(localMax.z, array[k][j][i].z);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm, fieldVec, &array); CHKERRQ(ierr);

        ierr = MPI_Allreduce(&localMin, &globalMin, 3, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Allreduce(&localMax, &globalMax, 3, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  [Rank %d] Local X-Range: [ %11.4e , %11.4e ]\n", user->simCtx->rank, localMin.x, localMax.x);
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  [Rank %d] Local Y-Range: [ %11.4e , %11.4e ]\n", user->simCtx->rank, localMin.y, localMax.y);
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "  [Rank %d] Local Z-Range: [ %11.4e , %11.4e ]\n", user->simCtx->rank, localMin.z, localMax.z);
        ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);

        if (user->simCtx->rank == 0) {
            PetscPrintf(PETSC_COMM_SELF, "  [Global] X-Range: [ %11.4e , %11.4e ]\n", globalMin.x, globalMax.x);
            PetscPrintf(PETSC_COMM_SELF, "  [Global] Y-Range: [ %11.4e , %11.4e ]\n", globalMin.y, globalMax.y);
            PetscPrintf(PETSC_COMM_SELF, "  [Global] Z-Range: [ %11.4e , %11.4e ]\n", globalMin.z, globalMax.z);
        }

    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "LogFieldStatistics only supports fields with 1 or 3 components, but field '%s' has %D.", fieldName, dof);
    }

    // --- 5. Final barrier for clean output ordering ---
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
    if (user->simCtx->rank == 0) {
        PetscPrintf(PETSC_COMM_SELF, "--------------------------------------------\n\n");
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LOG_FIELD_ANATOMY"
/**
 * @brief Logs the anatomy of a specified field at key boundary locations,
 *        respecting the solver's specific grid and variable architecture.
 *
 * This intelligent diagnostic function inspects a PETSc Vec and prints its values
 * at critical boundary locations (-Xi/+Xi, -Eta/+Eta, -Zeta/+Zeta). It is "architecture-aware":
 *
 * - **Cell-Centered Fields ("Ucat", "P"):** It correctly applies the "Shifted Index Architecture,"
 *   where the value for geometric `Cell i` is stored at array index `i+1`. It labels
 *   the output to clearly distinguish between true physical values and ghost values.
 * - **Face-Centered Fields ("Ucont"):** It uses a direct index mapping, where the value for
 *   the face at `Node i` is stored at index `i`.
 * - **Node-Centered Fields ("Coordinates"):** It uses a direct index mapping, where the value for
 *   `Node i` is stored at index `i`.
 *
 * The output is synchronized across MPI ranks to ensure readability and focuses on a
 * slice through the center of the domain to be concise.
 *
 * @param user       A pointer to the UserCtx structure containing the DMs and Vecs.
 * @param field_name A string identifier for the field to log (e.g., "Ucat", "P", "Ucont", "Coordinates").
 * @param stage_name A string identifier for the current simulation stage (e.g., "After Advection").
 * @return           PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode LOG_FIELD_ANATOMY(UserCtx *user, const char *field_name, const char *stage_name)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscMPIInt    rank;

    Vec            vec_local = NULL;
    DM             dm = NULL;
    PetscInt       dof = 0;
    char           data_layout[20]; // To store "Cell-Centered", "Face-Centered", etc.

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- 1. Map string name to PETSc objects and determine data layout ---
    if (strcasecmp(field_name, "Ucat") == 0) {
        vec_local = user->lUcat; dm = user->fda; dof = 3; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "P") == 0) {
        vec_local = user->lP; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "Psi") == 0) {
        vec_local = user->lPsi; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "Ucont") == 0) {
        vec_local = user->lUcont; dm = user->fda; dof = 3; strcpy(data_layout, "Face-Centered");
    } else if (strcasecmp(field_name, "Coordinates") == 0) {
        ierr = DMGetCoordinatesLocal(user->da, &vec_local); CHKERRQ(ierr);
        dm = user->fda; dof = 3; strcpy(data_layout, "Node-Centered");
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown field name for LOG_FIELD_ANATOMY: %s", field_name);
    }
    
    // --- 2. Get Grid Info and Array Pointers ---
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);

    // Synchronize for clean output
    ierr = PetscBarrier(NULL);
    PetscPrintf(PETSC_COMM_WORLD, "\n--- Field Anatomy Log: [%s] | Stage: [%s] | Layout: [%s] ---\n", field_name, stage_name, data_layout);

    // We will check a slice at the center of the other two dimensions
    PetscInt im_phys = user->IM, jm_phys = user->JM, km_phys = user->KM; // Number of physical nodes
    PetscInt i_mid = im_phys / 2;
    PetscInt j_mid = jm_phys / 2;
    PetscInt k_mid = km_phys / 2;

    // --- 3. Print Boundary Information based on Data Layout ---

    // ======================================================================
    // === CASE 1: Cell-Centered Fields (Ucat, P) - USES SHIFTED INDEX    ===
    // ======================================================================
    if (strcmp(data_layout, "Cell-Centered") == 0) {
        const void *l_arr; // Use void pointer for generic array access
        ierr = DMDAVecGetArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);

        // --- I-Direction Boundaries ---
        if (info.xs == 0) { // Rank on -Xi boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost for Cell 0)      = ", rank, 0);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][0]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][0].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][0].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][0].z);

            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Value for Cell 0)      = ", rank, 1);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][1]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][1].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][1].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][1].z);
        }
        if (info.xs + info.xm == info.mx) { // Rank on +Xi boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Value for Cell %d) = ", rank, im_phys - 1, im_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][im_phys - 1]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys - 1].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys - 1].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys - 1].z);

            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost for Cell %d) = ", rank, im_phys, im_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][im_phys]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys].z);
        }

        // --- J-Direction Boundaries ---
        if (info.ys == 0) { // Rank on -Eta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Ghost for Cell 0)      = ", rank, 0);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][0][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][0][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][0][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][0][i_mid].z);
            
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Value for Cell 0)      = ", rank, 1);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][1][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][1][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][1][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][1][i_mid].z);
        }

        if (info.ys + info.ym == info.my) { // Rank on +Eta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Value for Cell %d) = ", rank, jm_phys - 1, jm_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][jm_phys - 1][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][jm_phys - 1][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][jm_phys - 1][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][jm_phys - 1][i_mid].z);

            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Ghost for Cell %d) = ", rank, jm_phys, jm_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[k_mid][jm_phys][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[k_mid][jm_phys][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][jm_phys][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][jm_phys][i_mid].z);
        }

        // --- K-Direction Boundaries ---
        if (info.zs == 0) { // Rank on -Zeta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Ghost for Cell 0)      = ", rank, 0);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[0][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[0][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[0][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[0][j_mid][i_mid].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Value for Cell 0)      = ", rank, 1);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[1][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[1][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[1][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[1][j_mid][i_mid].z);
        }
        if (info.zs + info.zm == info.mz) { // Rank on +Zeta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Value for Cell %d) = ", rank, km_phys - 1, km_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[km_phys - 1][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[km_phys - 1][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[km_phys - 1][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[km_phys - 1][j_mid][i_mid].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Ghost for Cell %d) = ", rank, km_phys, km_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e)\n", ((const PetscReal***)l_arr)[km_phys][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.3e, %.3e, %.3e)\n", ((const Cmpnts***)l_arr)[km_phys][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[km_phys][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[km_phys][j_mid][i_mid].z);
        }

        ierr = DMDAVecRestoreArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);
    }
    // ======================================================================
    // === CASE 2: Face-Centered & Node-Centered Fields - USES DIRECT INDEX      ===
    // ======================================================================
    else if (strcmp(data_layout, "Face-Centered") == 0 || strcmp(data_layout, "Node-Centered") == 0) {
        const Cmpnts ***l_arr; // Both Ucont and Coordinates are Cmpnts
        ierr = DMDAVecGetArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);
        
        // --- I-Direction Boundaries ---
        if (info.xs == 0) { // Rank on -Xi boundary
            if(strcmp(data_layout, "Face-Centered")==0){
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (First Physical Face)     = (%.3e)\n", rank, 0,
                                    l_arr[k_mid][j_mid][0].x);
            }else if(strcmp(data_layout, "Node-Centered")==0){// Node-Centered
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (First Physical Face/Node) = (%.3e, %.3e, %.3e)\n", rank, 0,
                                    l_arr[k_mid][j_mid][0].x, l_arr[k_mid][j_mid][0].y, l_arr[k_mid][j_mid][0].z);
            }    
        }
        if (info.xs + info.xm == info.mx) { // Rank on +Xi boundary
            if(strcmp(data_layout, "Face-Centered")==0){
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Last Physical Face)      = (%.3e)\n", rank, im_phys - 1,
                                    l_arr[k_mid][j_mid][im_phys - 1].x);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Unused/Ghost Location)    = (%.3e)\n", rank, im_phys,
                                    l_arr[k_mid][j_mid][im_phys].x);
            }else if(strcmp(data_layout, "Node-Centered")==0){// Node-Centered
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Last Physical Face/Node)  = (%.3e, %.3e, %.3e)\n", rank, im_phys - 1,
                                    l_arr[k_mid][j_mid][im_phys - 1].x, l_arr[k_mid][j_mid][im_phys - 1].y, l_arr[k_mid][j_mid][im_phys - 1].z);
            // Also show the value at the unused memory location for verification
             PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Unused/Ghost Location)    = (%.3e, %.3e, %.3e)\n", rank, im_phys,
                                    l_arr[k_mid][j_mid][im_phys].x, l_arr[k_mid][j_mid][im_phys].y, l_arr[k_mid][j_mid][im_phys].z);
             }
        }                 
        
        // --- J-Direction

        if (info.ys == 0) { // Rank on -Eta boundary
            if(strcmp(data_layout, "Face-Centered")==0){
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (First Physical Face)     = (%.3e)\n", rank, 0,
                                    l_arr[k_mid][0][i_mid].y);
            }else if(strcmp(data_layout, "Node-Centered")==0){// Node-Centered
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (First Physical Face/Node) = (%.3e, %.3e, %.3e)\n", rank, 0,
                                    l_arr[k_mid][0][i_mid].x, l_arr[k_mid][0][i_mid].y, l_arr[k_mid][0][i_mid].z);
            }    
        }
        if (info.ys + info.ym == info.my) { // Rank on +Eta boundary
            if(strcmp(data_layout, "Face-Centered")==0){
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Last Physical Face)      = (%.3e)\n", rank, jm_phys - 1,
                                    l_arr[k_mid][jm_phys - 1][i_mid].y);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Unused/Ghost Location)    = (%.3e)\n", rank, jm_phys,
                                    l_arr[k_mid][jm_phys][i_mid].y);
            }else if(strcmp(data_layout, "Node-Centered")==0){// Node-Centered
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Last Physical Face/Node)  = (%.3e, %.3e, %.3e)\n", rank, jm_phys - 1,
                                    l_arr[k_mid][jm_phys - 1][i_mid].x, l_arr[k_mid][jm_phys - 1][i_mid].y, l_arr[k_mid][jm_phys - 1][i_mid].z);
            // Also show the value at the unused memory location for verification
             PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Idx %2d (Unused/Ghost Location)    = (%.3e, %.3e, %.3e)\n", rank, jm_phys,
                                    l_arr[k_mid][jm_phys][i_mid].x, l_arr[k_mid][jm_phys][i_mid].y, l_arr[k_mid][jm_phys][i_mid].z);
             }
        }
        // --- K-Direction ---
        if (info.zs == 0) { // Rank on -Zeta boundary
            if(strcmp(data_layout, "Face-Centered")==0){
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (First Physical Face)     = (%.3e)\n", rank, 0,
                                    l_arr[0][j_mid][i_mid].z);
            }else if(strcmp(data_layout, "Node-Centered")==0){// Node-Centered
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (First Physical Face/Node) = (%.3e, %.3e, %.3e)\n", rank, 0,
                                    l_arr[0][j_mid][i_mid].x, l_arr[0][j_mid][i_mid].y, l_arr[0][j_mid][i_mid].z);                       
            }
        }
        if(info.zs + info.zm == info.mz) { // Rank on +Zeta boundary
            if(strcmp(data_layout, "Face-Centered")==0){
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Last Physical Face)      = (%.3e)\n", rank, km_phys - 1,
                                    l_arr[km_phys - 1][j_mid][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Unused/Ghost Location)    = (%.3e)\n", rank, km_phys,
                                    l_arr[km_phys][j_mid][i_mid].z);
            }else if(strcmp(data_layout, "Node-Centered")==0){// Node-Centered
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Last Physical Face/Node)  = (%.3e, %.3e, %.3e)\n", rank, km_phys - 1,
                                    l_arr[km_phys - 1][j_mid][i_mid].x, l_arr[km_phys - 1][j_mid][i_mid].y, l_arr[km_phys - 1][j_mid][i_mid].z);
            // Also show the value at the unused memory location for verification
             PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Unused/Ghost Location)    = (%.3e, %.3e, %.3e)\n", rank, km_phys,
                                    l_arr[km_phys][j_mid][i_mid].x, l_arr[km_phys][j_mid][i_mid].y, l_arr[km_phys][j_mid][i_mid].z);
             }
        }

        ierr = DMDAVecRestoreArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);
    }
    else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "LOG_FIELD_ANATOMY only supports fields with 1 or 3 components & certain data-layouts, but field '%s' has %D components and an unsupported data-layout %s  \n", field_name, dof,data_layout);
    }
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    ierr = PetscBarrier(NULL);
    PetscFunctionReturn(0);
}
