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
 * @brief Implementation of \ref get_log_level().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see get_log_level()
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
 * @brief Internal helper implementation: `print_log_level()`.
 * @details Local to this translation unit.
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
               (level == LOG_VERBOSE) ? "VERBOSE" :
               (level == LOG_TRACE)   ? "TRACE"   :
               "UNKNOWN";

  /* print it out */
  ierr = PetscPrintf(PETSC_COMM_SELF,
                     "Current log level: %s (%d) | rank: %d\n",
                     level_name, level, (int)rank);
  CHKERRMPI(ierr);

  PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Implementation of \ref set_allowed_functions().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see set_allowed_functions()
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
 * @brief Implementation of \ref is_function_allowed().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see is_function_allowed()
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
 * @brief Implementation of \ref LOG_CELL_VERTICES().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LOG_CELL_VERTICES()
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
 * @brief Implementation of \ref LOG_FACE_DISTANCES().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LOG_FACE_DISTANCES()
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
 * @brief Implementation of \ref LOG_PARTICLE_FIELDS().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LOG_PARTICLE_FIELDS()
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

/**
 * @brief Implementation of \ref IsParticleConsoleSnapshotEnabled().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see IsParticleConsoleSnapshotEnabled()
 */

PetscBool IsParticleConsoleSnapshotEnabled(const SimCtx *simCtx)
{
    if (!simCtx) {
        return PETSC_FALSE;
    }
    return (PetscBool)(simCtx->np > 0 &&
                       simCtx->particleConsoleOutputFreq > 0 &&
                       get_log_level() >= LOG_INFO);
}

/**
 * @brief Implementation of \ref ShouldEmitPeriodicParticleConsoleSnapshot().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see ShouldEmitPeriodicParticleConsoleSnapshot()
 */

PetscBool ShouldEmitPeriodicParticleConsoleSnapshot(const SimCtx *simCtx, PetscInt completed_step)
{
    return (PetscBool)(IsParticleConsoleSnapshotEnabled(simCtx) &&
                       completed_step > 0 &&
                       completed_step % simCtx->particleConsoleOutputFreq == 0);
}

/**
 * @brief Implementation of \ref EmitParticleConsoleSnapshot().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see EmitParticleConsoleSnapshot()
 */

PetscErrorCode EmitParticleConsoleSnapshot(UserCtx *user, SimCtx *simCtx, PetscInt step)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    LOG(GLOBAL, LOG_INFO, "Particle states at step %d:\n", step);
    ierr = LOG_PARTICLE_FIELDS(user, simCtx->LoggingFrequency); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `trim()`.
 * @details Local to this translation unit.
 */
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
 * @brief Implementation of \ref LoadAllowedFunctionsFromFile().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LoadAllowedFunctionsFromFile()
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
 * @brief Internal helper implementation: `FreeAllowedFunctions()`.
 * @details Local to this translation unit.
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
 * @brief Implementation of \ref BCFaceToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see BCFaceToString()
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
 * @brief Implementation of \ref FieldInitializationToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see FieldInitializationToString()
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
 * @brief Implementation of \ref ParticleInitializationToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see ParticleInitializationToString()
 */
const char* ParticleInitializationToString(ParticleInitializationType ParticleInitialization)
{
    switch(ParticleInitialization){
        case PARTICLE_INIT_SURFACE_RANDOM: return "Surface: Random";
        case PARTICLE_INIT_VOLUME:         return "Volume";
        case PARTICLE_INIT_POINT_SOURCE:   return "Point Source";
        case PARTICLE_INIT_SURFACE_EDGES:  return "Surface: At edges";
        default: return "Unknown Particle Initialization";
    }
}

/**
 * @brief Implementation of \ref LESModelToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LESModelToString()
 */
const char* LESModelToString(LESModelType LESFlag)
{
    switch(LESFlag){
        case NO_LES_MODEL: return "No LES";
        case CONSTANT_SMAGORINSKY: return "Constant Smagorinsky";
        case DYNAMIC_SMAGORINSKY: return "Dynamic Smagorinsky";
        default: return "Unknown LES Flag";
    }
}

/**
 * @brief Implementation of \ref MomentumSolverTypeToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see MomentumSolverTypeToString()
 */
const char* MomentumSolverTypeToString(MomentumSolverType SolverFlag)
{
    switch(SolverFlag){
        case MOMENTUM_SOLVER_EXPLICIT_RK: return "Explicit 4 stage Runge-Kutta ";
        case MOMENTUM_SOLVER_DUALTIME_PICARD_RK4: return "Dual Time Stepping with Picard Iterations and 4 stage Runge-Kutta Smoothing";
        default: return "Unknown Momentum Solver Type";
    }
}

/**
 * @brief Implementation of \ref BCTypeToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see BCTypeToString()
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

	// case CUSTOM:    return "CUSTOM";
        default:        return "Unknown BC Type";
    }
}

/**
 * @brief Internal helper implementation: `BCHandlerTypeToString()`.
 * @details Local to this translation unit.
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

        // Multi-Block / Interface Handlers
        case BC_HANDLER_PERIODIC_GEOMETRIC:      return "geometric";
        case BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX:   return "constant flux";
        case BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX:    return "initial flux";
        case BC_HANDLER_INTERFACE_OVERSET:       return "overset";

        // Default case
        case BC_HANDLER_UNDEFINED:
        default:                                 return "UNKNOWN_HANDLER";
    }
}

/**
 * @brief Implementation of \ref DualMonitorDestroy().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see DualMonitorDestroy()
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
/**
 * @brief Implementation of \ref DualKSPMonitor().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see DualKSPMonitor()
 */

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
/**
 * @brief Implementation of \ref LOG_CONTINUITY_METRICS().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LOG_CONTINUITY_METRICS()
 */

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
        char filen[PETSC_MAX_PATH_LEN + 64];
        ierr = PetscSNPrintf(filen, sizeof(filen), "%s/Continuity_Metrics.log", simCtx->log_dir); CHKERRQ(ierr);

        // Open the log file in append mode.
        f = fopen(filen, "a");
        if (!f) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open log file: %s", filen);
        }

        // Write a header only when the file is empty and it's the first block (bi=0).
        // Using ftell() instead of step comparison ensures correctness across continuations.
        if (ftell(f) == 0 && bi == 0) {
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
 * @brief Implementation of \ref ParticleLocationStatusToString().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see ParticleLocationStatusToString()
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
/**
 * @brief Internal helper implementation: `_FindOrCreateEntry()`.
 * @details Local to this translation unit.
 */
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
 * @brief Internal helper implementation: `ProfilingInitialize()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ProfilingInitialize(SimCtx *simCtx)
{
    PetscFunctionBeginUser;
    if (!simCtx) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be null for ProfilingInitialize");

    // Iterate through the list of critical functions provided in SimCtx
    for (PetscInt i = 0; i < simCtx->nProfilingSelectedFuncs; ++i) {
        PetscInt idx;
        const char *func_name = simCtx->profilingSelectedFuncs[i];
        PetscErrorCode ierr = _FindOrCreateEntry(func_name, &idx); CHKERRQ(ierr);
        g_profiler_registry[idx].always_log = PETSC_TRUE;

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Marked '%s' as a critical function for profiling.\n", func_name);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref _ProfilingStart().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see _ProfilingStart()
 */

void _ProfilingStart(const char *func_name)
{
    PetscInt idx;
    if (_FindOrCreateEntry(func_name, &idx) != 0) return; // Fail silently
    PetscTime(&g_profiler_registry[idx].start_time);
}

/**
 * @brief Implementation of \ref _ProfilingEnd().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see _ProfilingEnd()
 */

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

/**
 * @brief Implementation of \ref ProfilingResetTimestepCounters().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see ProfilingResetTimestepCounters()
 */

PetscErrorCode ProfilingResetTimestepCounters(void)
{
    PetscFunctionBeginUser;
    for (PetscInt i = 0; i < g_profiler_count; ++i) {
        g_profiler_registry[i].current_step_time = 0.0;
        g_profiler_registry[i].current_step_call_count = 0;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ProfilingLogTimestepSummary().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see ProfilingLogTimestepSummary()
 */

PetscErrorCode ProfilingLogTimestepSummary(SimCtx *simCtx, PetscInt step)
{
    PetscBool should_write = PETSC_FALSE;
    FILE *f = NULL;
    char filen[(2 * PETSC_MAX_PATH_LEN) + 16];

    PetscFunctionBeginUser;
    if (!simCtx) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be null for ProfilingLogTimestepSummary");

    if (strcmp(simCtx->profilingTimestepMode, "off") == 0) {
        for (PetscInt i = 0; i < g_profiler_count; ++i) {
            g_profiler_registry[i].current_step_time = 0.0;
            g_profiler_registry[i].current_step_call_count = 0;
        }
        PetscFunctionReturn(0);
    }

    for (PetscInt i = 0; i < g_profiler_count; ++i) {
        if (g_profiler_registry[i].current_step_call_count <= 0) {
            continue;
        }
        if (strcmp(simCtx->profilingTimestepMode, "all") == 0 || g_profiler_registry[i].always_log) {
            should_write = PETSC_TRUE;
            break;
        }
    }

    if (should_write && simCtx->rank == 0) {
        snprintf(filen, sizeof(filen), "%s/%s", simCtx->log_dir, simCtx->profilingTimestepFile);
        if (step == simCtx->StartStep + 1 && !simCtx->continueMode) {
            f = fopen(filen, "w");
            if (!f) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open profiling timestep log file: %s", filen);
            }
            PetscFPrintf(PETSC_COMM_SELF, f, "step,function,calls,step_time_s\n");
        } else {
            f = fopen(filen, "a");
            if (!f) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open profiling timestep log file: %s", filen);
            }
            if (step == simCtx->StartStep + 1 && ftell(f) == 0) {
                PetscFPrintf(PETSC_COMM_SELF, f, "step,function,calls,step_time_s\n");
            }
        }

        for (PetscInt i = 0; i < g_profiler_count; ++i) {
            if (g_profiler_registry[i].current_step_call_count <= 0) {
                continue;
            }
            if (strcmp(simCtx->profilingTimestepMode, "all") == 0 || g_profiler_registry[i].always_log) {
                PetscFPrintf(
                    PETSC_COMM_SELF,
                    f,
                    "%d,%s,%lld,%.6f\n",
                    (int)step,
                    g_profiler_registry[i].name,
                    g_profiler_registry[i].current_step_call_count,
                    g_profiler_registry[i].current_step_time
                );
            }
        }
        fclose(f);
    }

    // Reset per-step counters for the next iteration
    for (PetscInt i = 0; i < g_profiler_count; ++i) {
        g_profiler_registry[i].current_step_time = 0.0;
        g_profiler_registry[i].current_step_call_count = 0;
    }
    PetscFunctionReturn(0);
}

// Comparison function for qsort to sort by total_time in descending order
/**
 * @brief Internal helper implementation: `_CompareProfiledFunctions()`.
 * @details Local to this translation unit.
 */
static int _CompareProfiledFunctions(const void *a, const void *b)
{
    const ProfiledFunction *func_a = (const ProfiledFunction *)a;
    const ProfiledFunction *func_b = (const ProfiledFunction *)b;

    if (func_a->total_time < func_b->total_time) return 1;
    if (func_a->total_time > func_b->total_time) return -1;
    return 0;
}

/**
 * @brief Implementation of \ref ProfilingFinalize().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see ProfilingFinalize()
 */
PetscErrorCode ProfilingFinalize(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscInt rank = simCtx->rank;
    PetscFunctionBeginUser;
    if (!simCtx->profilingFinalSummary) PetscFunctionReturn(0);
    if (!rank) {
        
        char exec_mode_modifier[32] = "Unknown";
        if(simCtx->exec_mode == EXEC_MODE_SOLVER) PetscCall(PetscStrncpy(exec_mode_modifier, "Solver", sizeof(exec_mode_modifier)));
        else if(simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR) PetscCall(PetscStrncpy(exec_mode_modifier, "PostProcessor", sizeof(exec_mode_modifier)));
        //--- Step 0: Create a file viewer for log file
        FILE *f;
        char filen[PETSC_MAX_PATH_LEN + 128];
        ierr = PetscSNPrintf(filen, sizeof(filen), "%s/ProfilingSummary_%s.log",simCtx->log_dir,exec_mode_modifier); CHKERRQ(ierr);

        // Open the log file: append with section label in continue mode, truncate otherwise.
        if (simCtx->continueMode) {
            f = fopen(filen, "a");
            if (!f) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open log file: %s", filen);
            }
            fprintf(f, "\n=== Continuation from step %" PetscInt_FMT " ===\n", simCtx->StartStep);
        } else {
            f = fopen(filen, "w");
            if (!f) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open log file: %s", filen);
            }
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
 * @brief Internal helper implementation: `PrintProgressBar()`.
 * @details Local to this translation unit.
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
    PetscPrintf(PETSC_COMM_SELF, "] %3d%% (Step %" PetscInt_FMT "/%" PetscInt_FMT ", t=%.4f)",
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
 * @brief Implementation of \ref LOG_FIELD_MIN_MAX().
 * @details Full API contract is documented with the header declaration in `include/logging.h`.
 * @see LOG_FIELD_MIN_MAX()
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
    } else if (strcasecmp(fieldName, "Diffusivity") == 0) {
        fieldVec = user->Diffusivity; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(fieldName, "DiffusivityGradient") == 0) {
        fieldVec = user->DiffusivityGradient; dm = user->fda; dof = 3; strcpy(data_layout, "Cell-Centered");
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
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "LogFieldStatistics only supports fields with 1 or 3 components, but field '%s' has %" PetscInt_FMT ".", fieldName, dof);
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
 * @brief Implementation of \ref LOG_FIELD_ANATOMY().
 * @details Full API contract is documented with the header declaration in `include/logging.h`.
 * @see LOG_FIELD_ANATOMY()
 */
PetscErrorCode LOG_FIELD_ANATOMY(UserCtx *user, const char *field_name, const char *stage_name)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscMPIInt    rank;

    Vec            vec_local = NULL;
    DM             dm = NULL;
    PetscInt       dof = 0;
    char           data_layout[20];
    char           dominant_dir = '\0'; // 'x', 'y', 'z' for face-centered, 'm' for mixed (Ucont)

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- 1. Map string name to PETSc objects and determine data layout ---
    if (strcasecmp(field_name, "Ucat") == 0) {
        vec_local = user->lUcat; dm = user->fda; dof = 3; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "P") == 0) {
        vec_local = user->lP; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "Diffusivity") == 0) {
        vec_local = user->lDiffusivity; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "DiffusivityGradient") == 0) {
        vec_local = user->lDiffusivityGradient; dm = user->fda; dof = 3; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "Psi") == 0) {
        vec_local = user->lPsi; dm = user->da; dof = 1; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "Center-Coordinates") == 0) {
        vec_local = user->lCent; dm = user->fda; dof = 3; strcpy(data_layout, "Cell-Centered");
    } else if (strcasecmp(field_name, "Ucont") == 0) {
        vec_local = user->lUcont; dm = user->fda; dof = 3; strcpy(data_layout, "Face-Centered"); dominant_dir = 'm'; // Mixed
    } else if (strcasecmp(field_name, "Csi") == 0 || strcasecmp(field_name, "X-Face-Centers") == 0) {
        vec_local = (strcasecmp(field_name, "Csi") == 0) ? user->lCsi : user->Centx;
        dm = user->fda; dof = 3; strcpy(data_layout, "Face-Centered"); dominant_dir = 'x';
    } else if (strcasecmp(field_name, "Eta") == 0 || strcasecmp(field_name, "Y-Face-Centers") == 0) {
        vec_local = (strcasecmp(field_name, "Eta") == 0) ? user->lEta : user->Centy;
        dm = user->fda; dof = 3; strcpy(data_layout, "Face-Centered"); dominant_dir = 'y';
    } else if (strcasecmp(field_name, "Zet") == 0 || strcasecmp(field_name, "Z-Face-Centers") == 0) {
        vec_local = (strcasecmp(field_name, "Zet") == 0) ? user->lZet : user->Centz;
        dm = user->fda; dof = 3; strcpy(data_layout, "Face-Centered"); dominant_dir = 'z';
    } else if (strcasecmp(field_name, "Coordinates") == 0) {
        ierr = DMGetCoordinatesLocal(user->da, &vec_local); CHKERRQ(ierr);
        dm = user->fda; dof = 3; strcpy(data_layout, "Node-Centered");
    } else if  (strcasecmp(field_name, "CornerField")== 0){
        vec_local = user->lCellFieldAtCorner; strcpy(data_layout, "Node-Centered");
        PetscInt bs = 1;
        ierr = VecGetBlockSize(user->CellFieldAtCorner, &bs); CHKERRQ(ierr);
        dof = bs;
        if(dof == 1) dm = user->da;
        else         dm = user->fda;
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown field name for LOG_FIELD_ANATOMY: %s", field_name);
    }

    // --- 2. Get Grid Info and Array Pointers ---
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);

    ierr = PetscBarrier(NULL);
    PetscPrintf(PETSC_COMM_WORLD, "\n--- Field Anatomy Log: [%s] | Stage: [%s] | Layout: [%s] ---\n", field_name, stage_name, data_layout);

    // Global physical dimensions (number of cells)
    PetscInt im_phys = user->IM;
    PetscInt jm_phys = user->JM;
    PetscInt km_phys = user->KM;

    // Slice through the center of the local domain
    PetscInt i_mid = (PetscInt)(info.xs + 0.5 * info.xm) - 1;
    PetscInt j_mid = (PetscInt)(info.ys + 0.5 * info.ym) - 1;
    PetscInt k_mid = (PetscInt)(info.zs + 0.5 * info.zm) - 1;

    // --- 3. Print Boundary Information based on Data Layout ---

    // ======================================================================
    // === CASE 1: Cell-Centered Fields (Ucat, P) - USES SHIFTED INDEX    ===
    // ======================================================================
    if (strcmp(data_layout, "Cell-Centered") == 0) {
        const void *l_arr;
        ierr = DMDAVecGetArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);


        // --- I-Direction Boundaries ---
        if (info.xs == 0) { // Rank on -Xi boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost for Cell[k][j][0])      = ", rank, 0);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][0]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][0].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][0].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][0].z);

            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Value for Cell[k][j][0])      = ", rank, 1);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][1]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][1].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][1].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][1].z);
        }
        if (info.xs + info.xm == info.mx) { // Rank on +Xi boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Value for Cell[k][j][%d]) = ", rank, im_phys - 1, im_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][im_phys - 1]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys - 1].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys - 1].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys - 1].z);

            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost for Cell[k][j][%d]) = ", rank, im_phys, im_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][j_mid][im_phys]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys].x, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys].y, ((const Cmpnts***)l_arr)[k_mid][j_mid][im_phys].z);
        } 

        // --- J-Direction Boundaries ---
        if (info.ys == 0) { // Rank on -Eta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Ghost for Cell[k][0][i])      = ", rank, 0);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][0][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][0][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][0][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][0][i_mid].z);
            
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Value for Cell[k][0][i])      = ", rank, 1);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][1][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][1][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][1][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][1][i_mid].z);
        }

        if (info.ys + info.ym == info.my) { // Rank on +Eta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Value for Cell[k][%d][i]) = ", rank, jm_phys - 1, jm_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][jm_phys - 1][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][jm_phys - 1][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][jm_phys - 1][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][jm_phys - 1][i_mid].z);

            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Ghost for Cell[k][%d][i]) = ", rank, jm_phys, jm_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[k_mid][jm_phys][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[k_mid][jm_phys][i_mid].x, ((const Cmpnts***)l_arr)[k_mid][jm_phys][i_mid].y, ((const Cmpnts***)l_arr)[k_mid][jm_phys][i_mid].z);
        }

        // --- K-Direction Boundaries ---
        if (info.zs == 0) { // Rank on -Zeta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Ghost for Cell[0][j][i])      = ", rank, 0);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[0][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[0][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[0][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[0][j_mid][i_mid].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Value for Cell[0][j][i])      = ", rank, 1);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[1][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[1][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[1][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[1][j_mid][i_mid].z);
        }
        if (info.zs + info.zm == info.mz) { // Rank on +Zeta boundary
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Value for Cell[%d][j][i]) = ", rank, km_phys - 1, km_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[km_phys - 1][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[km_phys - 1][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[km_phys - 1][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[km_phys - 1][j_mid][i_mid].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Ghost for Cell[%d][j][i]) = ", rank, km_phys, km_phys - 2);
            if(dof==1) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f)\n", ((const PetscReal***)l_arr)[km_phys][j_mid][i_mid]);
            else       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "(%.5f, %.5f, %.5f)\n", ((const Cmpnts***)l_arr)[km_phys][j_mid][i_mid].x, ((const Cmpnts***)l_arr)[km_phys][j_mid][i_mid].y, ((const Cmpnts***)l_arr)[km_phys][j_mid][i_mid].z);
        }
        ierr = DMDAVecRestoreArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);
    }
    // ======================================================================
    // === CASE 2: Face-Centered Fields - NUANCED DIRECTIONAL LOGIC       ===
    // ======================================================================
    else if (strcmp(data_layout, "Face-Centered") == 0) {
        const Cmpnts ***l_arr;
        ierr = DMDAVecGetArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);

        // --- I-Direction Boundaries ---
        if (info.xs == 0) { // Rank on -Xi boundary
            if (dominant_dir == 'x') { // Node-like in I-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (First Phys. X-Face) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[k_mid][j_mid][0].x, l_arr[k_mid][j_mid][0].y, l_arr[k_mid][j_mid][0].z);
            } else if (dominant_dir == 'y' || dominant_dir == 'z') { // Cell-like in I-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost for Cell[k][j][0]) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[k_mid][j_mid][0].x, l_arr[k_mid][j_mid][0].y, l_arr[k_mid][j_mid][0].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Value for Cell[k][j][0]) = (%.5f, %.5f, %.5f)\n", rank, 1, l_arr[k_mid][j_mid][1].x, l_arr[k_mid][j_mid][1].y, l_arr[k_mid][j_mid][1].z);
            } else if (dominant_dir == 'm') { // Ucont: Mixed
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: u-comp @ Idx %2d (1st X-Face) = %.5f\n", rank, 0, l_arr[k_mid][j_mid][0].x);
            }
        }
        if (info.xs + info.xm == info.mx) { // Rank on +Xi boundary
            if (dominant_dir == 'x') { // Node-like in I-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Last Phys. X-Face)  = (%.5f, %.5f, %.5f)\n", rank, im_phys - 1, l_arr[k_mid][j_mid][im_phys - 1].x, l_arr[k_mid][j_mid][im_phys-1].y, l_arr[k_mid][j_mid][im_phys - 1].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost Location)     = (%.5f, %.5f, %.5f)\n", rank, im_phys, l_arr[k_mid][j_mid][im_phys].x, l_arr[k_mid][j_mid][im_phys].y, l_arr[k_mid][j_mid][im_phys].z);
            } else if (dominant_dir == 'y' || dominant_dir == 'z') { // Cell-like in I-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Value for Cell[k][j][%d]) = (%.5f, %.5f, %.5f)\n", rank, im_phys - 1, im_phys - 2, l_arr[k_mid][j_mid][im_phys - 1].x, l_arr[k_mid][j_mid][im_phys - 1].y, l_arr[k_mid][j_mid][im_phys-1].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Ghost for Cell[k][j][%d]) = (%.5f, %.5f, %.5f)\n", rank, im_phys, im_phys - 2, l_arr[k_mid][j_mid][im_phys].x, l_arr[k_mid][j_mid][im_phys].y, l_arr[k_mid][j_mid][im_phys].z);
            } else if (dominant_dir == 'm') { // Ucont: Mixed
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: u-comp @ Idx %2d (Last X-Face)   = %.5f\n", rank, im_phys - 1, l_arr[k_mid][j_mid][im_phys - 1].x);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: u-comp @ Idx %2d (Ghost Location)    = %.5f\n", rank, im_phys, l_arr[k_mid][j_mid][im_phys].x);
            }
        }

        // --- J-Direction Boundaries ---
        if (info.ys == 0) { // Rank on -Eta boundary
            if (dominant_dir == 'y') { // Node-like in J-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (First Phys. Y-Face) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[k_mid][0][i_mid].x, l_arr[k_mid][0][i_mid].y, l_arr[k_mid][0][i_mid].z);
            } else if (dominant_dir == 'x' || dominant_dir == 'z') { // Cell-like in J-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Ghost for Cell[k][0][i]) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[k_mid][0][i_mid].x, l_arr[k_mid][0][i_mid].y, l_arr[k_mid][0][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Value for Cell[k][0][i]) = (%.5f, %.5f, %.5f)\n", rank, 1, l_arr[k_mid][1][i_mid].x, l_arr[k_mid][1][i_mid].y, l_arr[k_mid][1][i_mid].z);
            } else if (dominant_dir == 'm') { // Ucont: Mixed
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: v-comp @ Jdx %2d (1st Y-Face) = %.5f\n", rank, 0, l_arr[k_mid][0][i_mid].y);
            }
        }
        if (info.ys + info.ym == info.my) { // Rank on +Eta boundary
             if (dominant_dir == 'y') { // Node-like in J-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Last Phys. Y-Face)  = (%.5f, %.5f, %.5f)\n", rank, jm_phys - 1, l_arr[k_mid][jm_phys - 1][i_mid].x, l_arr[k_mid][jm_phys - 1][i_mid].y, l_arr[k_mid][jm_phys - 1][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Ghost Location)     = (%.5f, %.5f, %.5f)\n", rank, jm_phys, l_arr[k_mid][jm_phys][i_mid].x, l_arr[k_mid][jm_phys][i_mid].y, l_arr[k_mid][jm_phys][i_mid].z);
            } else if (dominant_dir == 'x' || dominant_dir == 'z') { // Cell-like in J-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Value for Cell[k][%d][i]) = (%.5f, %.5f, %.5f)\n", rank, jm_phys-1, jm_phys-2, l_arr[k_mid][jm_phys - 1][i_mid].x, l_arr[k_mid][jm_phys - 1][i_mid].y, l_arr[k_mid][jm_phys - 1][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Ghost for Cell[k][%d][i]) = (%.5f, %.5f, %.5f)\n", rank, jm_phys, jm_phys-2, l_arr[k_mid][jm_phys][i_mid].x, l_arr[k_mid][jm_phys][i_mid].y, l_arr[k_mid][jm_phys][i_mid].z);
            } else if (dominant_dir == 'm') { // Ucont: Mixed
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: v-comp @ Jdx %2d (Last Y-Face)    = %.5f\n", rank, jm_phys - 1, l_arr[k_mid][jm_phys - 1][i_mid].y); 
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: v-comp @ Jdx %2d (Ghost Location)   = %.5f\n", rank, jm_phys, l_arr[k_mid][jm_phys][i_mid].y);
            }
        }

        // --- K-Direction Boundaries ---
        if (info.zs == 0) { // Rank on -Zeta boundary
            if (dominant_dir == 'z') { // Node-like in K-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (First Phys. Z-Face) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[0][j_mid][i_mid].x, l_arr[0][j_mid][i_mid].y, l_arr[0][j_mid][i_mid].z);
            } else if (dominant_dir == 'x' || dominant_dir == 'y') { // Cell-like in K-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Ghost for Cell[0][j][i]) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[0][j_mid][i_mid].x, l_arr[0][j_mid][i_mid].y, l_arr[0][j_mid][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Value for Cell[0][j][i]) = (%.5f, %.5f, %.5f)\n", rank, 1, l_arr[1][j_mid][i_mid].x, l_arr[1][j_mid][i_mid].y, l_arr[1][j_mid][i_mid].z);
            } else if (dominant_dir == 'm') { // Ucont: Mixed
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: w-comp @ Idx %2d (1st Z-Face) = %.5f\n", rank, 0, l_arr[0][j_mid][i_mid].z);
            }
        }
        if (info.zs + info.zm == info.mz) { // Rank on +Zeta boundary
             if (dominant_dir == 'z') { // Node-like in K-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Last Phys. Z-Face)  = (%.5f, %.5f, %.5f)\n", rank, km_phys - 1, l_arr[km_phys - 1][j_mid][i_mid].x, l_arr[km_phys - 1][j_mid][i_mid].y, l_arr[km_phys - 1][j_mid][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Ghost Location)     = (%.5f, %.5f, %.5f)\n", rank, km_phys, l_arr[km_phys][j_mid][i_mid].x, l_arr[km_phys][j_mid][i_mid].y, l_arr[km_phys][j_mid][i_mid].z);
            } else if (dominant_dir == 'x' || dominant_dir == 'y') { // Cell-like in K-dir
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Value for Cell[%d][j][i]) = (%.5f, %.5f, %.5f)\n", rank, km_phys-1, km_phys-2, l_arr[km_phys-1][j_mid][i_mid].x, l_arr[km_phys-1][j_mid][i_mid].y, l_arr[km_phys - 1][j_mid][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Idx %2d (Ghost for Cell[%d][j][i]) = (%.5f, %.5f, %.5f)\n", rank, km_phys, km_phys-2, l_arr[km_phys][j_mid][i_mid].x, l_arr[km_phys][j_mid][i_mid].y, l_arr[km_phys][j_mid][i_mid].z);
            } else if (dominant_dir == 'm') { // Ucont: Mixed
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: w-comp @ Idx %2d (Last Z-Face) = %.5f\n", rank, km_phys - 1, l_arr[km_phys - 1][j_mid][i_mid].z);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: w-comp @ Idx %2d (Ghost Loc.) = %.5f\n", rank, km_phys, l_arr[km_phys][j_mid][i_mid].z);

            }
        }
        ierr = DMDAVecRestoreArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);
    }
    // ======================================================================
    // === CASE 3: Node-Centered Fields - USES DIRECT INDEX               ===
    // ======================================================================
    else if (strcmp(data_layout, "Node-Centered") == 0) {
        const Cmpnts ***l_arr;
        ierr = DMDAVecGetArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);

        // --- I-Direction Boundaries ---
        if (info.xs == 0) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (First Phys. Node) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[k_mid][j_mid][0].x, l_arr[k_mid][j_mid][0].y, l_arr[k_mid][j_mid][0].z);
        }
        if (info.xs + info.xm == info.mx) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Last Phys. Node)  = (%.5f, %.5f, %.5f)\n", rank, im_phys - 1, l_arr[k_mid][j_mid][im_phys - 1].x, l_arr[k_mid][j_mid][im_phys - 1].y, l_arr[k_mid][j_mid][im_phys - 1].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, I-DIR]: Idx %2d (Unused/Ghost Loc) = (%.5f, %.5f, %.5f)\n", rank, im_phys, l_arr[k_mid][j_mid][im_phys].x, l_arr[k_mid][j_mid][im_phys].y, l_arr[k_mid][j_mid][im_phys].z);
        }
        // --- J-Direction Boundaries ---
        if (info.ys == 0) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (First Phys. Node) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[k_mid][0][i_mid].x, l_arr[k_mid][0][i_mid].y, l_arr[k_mid][0][i_mid].z);
        }
        if (info.ys + info.ym == info.my) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Last Phys. Node)  = (%.5f, %.5f, %.5f)\n", rank, jm_phys - 1, l_arr[k_mid][jm_phys - 1][i_mid].x, l_arr[k_mid][jm_phys - 1][i_mid].y, l_arr[k_mid][jm_phys - 1][i_mid].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, J-DIR]: Jdx %2d (Unused/Ghost Loc) = (%.5f, %.5f, %.5f)\n", rank, jm_phys, l_arr[k_mid][jm_phys][i_mid].x, l_arr[k_mid][jm_phys][i_mid].y, l_arr[k_mid][jm_phys][i_mid].z);
        }
        // --- K-Direction Boundaries ---
        if (info.zs == 0) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (First Phys. Node) = (%.5f, %.5f, %.5f)\n", rank, 0, l_arr[0][j_mid][i_mid].x, l_arr[0][j_mid][i_mid].y, l_arr[0][j_mid][i_mid].z);
        }
        if(info.zs + info.zm == info.mz) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Last Phys. Node)  = (%.5f, %.5f, %.5f)\n", rank, km_phys - 1, l_arr[km_phys - 1][j_mid][i_mid].x, l_arr[km_phys - 1][j_mid][i_mid].y, l_arr[km_phys - 1][j_mid][i_mid].z);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[Rank %d, K-DIR]: Kdx %2d (Unused/Ghost Loc) = (%.5f, %.5f, %.5f)\n", rank, km_phys, l_arr[km_phys][j_mid][i_mid].x, l_arr[km_phys][j_mid][i_mid].y, l_arr[km_phys][j_mid][i_mid].z);
        }
        ierr = DMDAVecRestoreArrayRead(dm, vec_local, (void*)&l_arr); CHKERRQ(ierr);
    }
    else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "LOG_FIELD_ANATOMY encountered an unknown data layout: %s", data_layout);
    }

    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    ierr = PetscBarrier(NULL);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LOG_INTERPOLATION_ERROR"
/**
 * @brief Implementation of \ref LOG_INTERPOLATION_ERROR().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LOG_INTERPOLATION_ERROR()
 */
PetscErrorCode LOG_INTERPOLATION_ERROR(UserCtx *user)
{
    SimCtx *simCtx = user->simCtx;
    PetscErrorCode ierr;
    DM swarm = user->swarm;
    Vec positionVec, analyticalvelocityVec, velocityVec, errorVec;
    PetscReal Interpolation_error = 0.0;
    PetscReal Maximum_Interpolation_error = 0.0;
    PetscReal AnalyticalSolution_magnitude = 0.0;
    PetscReal ErrorPercentage = 0.0;
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Creating global vectors.\n");
    ierr = DMSwarmCreateGlobalVectorFromField(swarm, "position", &positionVec); CHKERRQ(ierr);
    ierr = DMSwarmCreateGlobalVectorFromField(swarm, "velocity", &velocityVec); CHKERRQ(ierr);
    
    ierr = VecDuplicate(positionVec, &analyticalvelocityVec); CHKERRQ(ierr);
    ierr = VecCopy(positionVec, analyticalvelocityVec); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Computing analytical solution.\n");
    ierr = SetAnalyticalSolutionForParticles(analyticalvelocityVec, simCtx); CHKERRQ(ierr);
    
    ierr = VecDuplicate(analyticalvelocityVec, &errorVec); CHKERRQ(ierr);
    ierr = VecCopy(analyticalvelocityVec, errorVec); CHKERRQ(ierr);
    
    ierr = VecNorm(analyticalvelocityVec, NORM_2, &AnalyticalSolution_magnitude); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Computing error.\n");
    ierr = VecAXPY(errorVec, -1.0, velocityVec); CHKERRQ(ierr);
    ierr = VecNorm(errorVec, NORM_2, &Interpolation_error); CHKERRQ(ierr);
    ierr = VecNorm(errorVec,NORM_INFINITY,&Maximum_Interpolation_error); CHKERRQ(ierr);
    
    ErrorPercentage = (AnalyticalSolution_magnitude > 0) ?
                      (Interpolation_error / AnalyticalSolution_magnitude * 100.0) : 0.0;

    /* --- CSV output (always, rank 0 only) --- */
    if (simCtx->rank == 0) {
        char csv_path[PETSC_MAX_PATH_LEN];
        snprintf(csv_path, sizeof(csv_path), "%s/interpolation_error.csv", simCtx->log_dir);
        FILE *f = fopen(csv_path, "a");
        if (f) {
            if (ftell(f) == 0) {
                fprintf(f, "step,time,L2_error,Linf_error,L2_analytical,error_pct\n");
            }
            PetscReal t = (PetscReal)simCtx->ti * simCtx->dt;
            fprintf(f, "%d,%.6e,%.6e,%.6e,%.6e,%.4f\n",
                    (int)simCtx->step, t,
                    Interpolation_error, Maximum_Interpolation_error,
                    AnalyticalSolution_magnitude, ErrorPercentage);
            fclose(f);
        }
    }

    /* --- Console output (only at INFO level or above) --- */
    if (get_log_level() >= LOG_INFO) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Interpolation error (%%): %g\n", ErrorPercentage);
        PetscPrintf(PETSC_COMM_WORLD, "Interpolation error (%%): %g\n", ErrorPercentage);
        LOG_ALLOW(GLOBAL, LOG_INFO, "Maximum Interpolation error: %g\n", Maximum_Interpolation_error);
        PetscPrintf(PETSC_COMM_WORLD, "Maximum Interpolation error: %g\n", Maximum_Interpolation_error);
    }

    ierr = VecDestroy(&analyticalvelocityVec); CHKERRQ(ierr);
    ierr = VecDestroy(&errorVec); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(swarm, "position", &positionVec); CHKERRQ(ierr);
    ierr = DMSwarmDestroyGlobalVectorFromField(swarm, "velocity", &velocityVec); CHKERRQ(ierr);
    
    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "CalculateAdvancedParticleMetrics"
/**
 * @brief Internal helper implementation: `CalculateAdvancedParticleMetrics()`.
 * @details Local to this translation unit.
 */
PetscErrorCode CalculateAdvancedParticleMetrics(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    PetscMPIInt    size, rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- 1. Particle Load Imbalance ---
    PetscInt nLocal, nGlobal, nLocalMax;
    ierr = DMSwarmGetLocalSize(user->swarm, &nLocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(user->swarm, &nGlobal); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&nLocal, &nLocalMax, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscReal avg_per_rank = (size > 0) ? ((PetscReal)nGlobal / size) : 0.0;
    // Handle division by zero if there are no particles
    simCtx->particleLoadImbalance = (avg_per_rank > 1e-9) ? (nLocalMax / avg_per_rank) : 1.0;


    // --- 2. Number of Occupied Cells ---
    // This part requires access to the user->ParticleCount vector.
    PetscInt       local_occupied_cells = 0;
    PetscInt       global_occupied_cells;
    const PetscScalar *count_array;
    PetscInt       vec_local_size;

    ierr = VecGetLocalSize(user->ParticleCount, &vec_local_size); CHKERRQ(ierr);
    ierr = VecGetArrayRead(user->ParticleCount, &count_array); CHKERRQ(ierr);

    for (PetscInt i = 0; i < vec_local_size; ++i) {
        if (count_array[i] > 0.5) { // Use 0.5 to be safe with floating point
            local_occupied_cells++;
        }
    }
    ierr = VecRestoreArrayRead(user->ParticleCount, &count_array); CHKERRQ(ierr);

    ierr = MPI_Allreduce(&local_occupied_cells, &global_occupied_cells, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    simCtx->occupiedCellCount = global_occupied_cells;

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Advanced Metrics: Imbalance=%.2f, OccupiedCells=%d\n", rank, simCtx->particleLoadImbalance, simCtx->occupiedCellCount);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LOG_PARTICLE_METRICS"
/**
 * @brief Implementation of \ref LOG_PARTICLE_METRICS().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see LOG_PARTICLE_METRICS()
 */
PetscErrorCode LOG_PARTICLE_METRICS(UserCtx *user, const char *stageName)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    SimCtx         *simCtx = user->simCtx;
    const char     *stage_label = (stageName && stageName[0] != '\0') ? stageName : "N/A";

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    PetscInt totalParticles;
    ierr = DMSwarmGetSize(user->swarm, &totalParticles); CHKERRQ(ierr);

    if (!rank) {
        FILE *f;
        char filen[PETSC_MAX_PATH_LEN + 64];
        ierr = PetscSNPrintf(filen, sizeof(filen), "%s/Particle_Metrics.log", simCtx->log_dir); CHKERRQ(ierr);
        f = fopen(filen, "a");
        if (!f) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open particle log file: %s", filen);

        if (ftell(f) == 0) {
            PetscFPrintf(PETSC_COMM_SELF, f, "%-18s | %-10s | %-12s | %-10s | %-10s | %-15s | %-10s | %-10s\n",
                            "Stage", "Timestep", "Total Ptls", "Lost", "Migrated", "Occupied Cells", "Imbalance", "Mig Passes");
            PetscFPrintf(PETSC_COMM_SELF, f, "----------------------------------------------------------------------------------------------------------------------------\n");
        }

        PetscFPrintf(PETSC_COMM_SELF, f, "%-18s | %-10d | %-12d | %-10d | %-10d | %-15d | %-10.2f | %-10d\n",
                        stage_label, (int)simCtx->step, (int)totalParticles, (int)simCtx->particlesLostLastStep,
                        (int)simCtx->particlesMigratedLastStep, (int)simCtx->occupiedCellCount,
                        (double)simCtx->particleLoadImbalance, (int)simCtx->migrationPassesLastStep);
        fclose(f);
    }
    PetscFunctionReturn(0);
}
