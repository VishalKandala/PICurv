/**
 * @file io.c
 * @brief Implementation of data input/output routines, focusing on grid configuration.
 *
 * This module provides functions to parse grid geometry information, either from
 * command-line options for programmatically generated grids or by reading the
 * header of a grid definition file.
 */

#include "io.h"

// =============================================================================
//          STATIC (PRIVATE) VARIABLES FOR ONE-TIME FILE READ
// =============================================================================

/** @brief Stores the number of blocks read from the grid file. */
static PetscInt  g_nblk_from_file = 0;
/** @brief Caches the IM dimensions for all blocks read from the grid file. */
static PetscInt* g_IMs_from_file = NULL;
/** @brief Caches the JM dimensions for all blocks read from the grid file. */
static PetscInt* g_JMs_from_file = NULL;
/** @brief Caches the KM dimensions for all blocks read from the grid file. */
static PetscInt* g_KMs_from_file = NULL;
/** @brief A flag to ensure the grid file is read only once. */
static PetscBool g_file_has_been_read = PETSC_FALSE;


// =============================================================================
//                      PUBLIC FUNCTION IMPLEMENTATIONS
// =============================================================================

/**
 * @brief Trims leading and trailing whitespace from a string in-place.
 * @param str The string to be trimmed.
 */
void TrimWhitespace(char *str) {
    if (!str) return;

    char *start = str;
    // Find the first non-whitespace character
    while (isspace((unsigned char)*start)) {
        start++;
    }

    // Find the end of the string
    char *end = str + strlen(str) - 1;
    // Move backwards from the end to find the last non-whitespace character
    while (end > start && isspace((unsigned char)*end)) {
        end--;
    }

    // Null-terminate after the last non-whitespace character
    *(end + 1) = '\0';

    // If there was leading whitespace, shift the string to the left
    if (str != start) {
        memmove(str, start, (end - start) + 2); // +2 to include the new null terminator
    }
}


#undef __FUNCT__
#define __FUNCT__ "ReadGridGenerationInputs"
/**
 * @brief Parses command-line options for a programmatically generated grid for a SINGLE block.
 *
 * This function is responsible for reading all per-block array options related
 * to grid geometry, such as dimensions (`-im`), domain bounds (`-xMins`, `-xMaxs`),
 * and stretching ratios (`-rxs`). It reads the entire array for each option, then
 * uses the block index stored in `user->_this` to select the correct value and
 * populate the fields of the provided `UserCtx` struct.
 *
 * @param user Pointer to the `UserCtx` for a specific block. The function will
 *             populate the geometric fields (`IM`, `Min_X`, `rx`, etc.) within this struct.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ReadGridGenerationInputs(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    PetscInt       nblk = simCtx->block_number;
    PetscInt       block_index = user->_this;
    PetscBool      found;

    // Temporary arrays to hold the parsed values for ALL blocks
    PetscInt  *IMs = NULL, *JMs = NULL, *KMs = NULL, *cgrids = NULL;
    PetscReal *xMins = NULL, *xMaxs = NULL, *rxs = NULL;
    PetscReal *yMins = NULL, *yMaxs = NULL, *rys = NULL;
    PetscReal *zMins = NULL, *zMaxs = NULL, *rzs = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reading generated grid inputs for block %d.\n", simCtx->rank, block_index);

    if (block_index >= nblk) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Block index %d is out of range for nblk=%d", block_index, nblk);
    }

    // --- Allocate temporary storage for all array options ---
    ierr = PetscMalloc4(nblk, &IMs, nblk, &JMs, nblk, &KMs, nblk, &cgrids); CHKERRQ(ierr);
    ierr = PetscMalloc6(nblk, &xMins, nblk, &xMaxs, nblk, &rxs, nblk, &yMins, nblk, &yMaxs, nblk, &rys); CHKERRQ(ierr);
    ierr = PetscMalloc3(nblk, &zMins, nblk, &zMaxs, nblk, &rzs); CHKERRQ(ierr);

    // --- Set default values for the temporary arrays ---
    for (PetscInt i = 0; i < nblk; ++i) {
        IMs[i] = 10; JMs[i] = 10; KMs[i] = 10; cgrids[i] = 0;
        xMins[i] = 0.0; xMaxs[i] = 1.0; rxs[i] = 1.0;
        yMins[i] = 0.0; yMaxs[i] = 1.0; rys[i] = 1.0;
        zMins[i] = 0.0; zMaxs[i] = 1.0; rzs[i] = 1.0;
    }

    // --- Parse the array options from the command line / control file ---
    PetscInt count;
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", IMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", JMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", KMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-xMins", xMins, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-xMaxs", xMaxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-rxs", rxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-yMins", yMins, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-yMaxs", yMaxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-rys", rys, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-zMins", zMins, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-zMaxs", zMaxs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetRealArray(NULL, NULL, "-rzs", rzs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-cgrids", cgrids, &count, &found); CHKERRQ(ierr);

    // --- Assign the parsed values to the specific UserCtx struct passed in ---
    user->IM = IMs[block_index];
    user->JM = JMs[block_index];
    user->KM = KMs[block_index];
    user->Min_X = xMins[block_index];
    user->Max_X = xMaxs[block_index];
    user->rx = rxs[block_index];
    user->Min_Y = yMins[block_index];
    user->Max_Y = yMaxs[block_index];
    user->ry = rys[block_index];
    user->Min_Z = zMins[block_index];
    user->Max_Z = zMaxs[block_index];
    user->rz = rzs[block_index];
    user->cgrid = cgrids[block_index];

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d grid generation inputs set: IM=%d, JM=%d, KM=%d\n",
              simCtx->rank, block_index, user->IM, user->JM, user->KM);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d bounds: X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f]\n",
              simCtx->rank, block_index, user->Min_X, user->Max_X, user->Min_Y, user->Max_Y, user->Min_Z, user->Max_Z);

    // --- Clean up temporary storage ---
    ierr = PetscFree4(IMs, JMs, KMs, cgrids); CHKERRQ(ierr);
    ierr = PetscFree6(xMins, xMaxs, rxs, yMins, yMaxs, rys); CHKERRQ(ierr);
    ierr = PetscFree3(zMins, zMaxs, rzs); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ReadGridFile"
/**
 * @brief Sets grid dimensions from a file for a SINGLE block using a one-time read cache.
 *
 * This function uses a static-variable pattern to ensure the grid file header
 * is read only once, collectively, by all processes. The first time any process
 * calls this function, it triggers a collective operation where rank 0 reads the
 * file and broadcasts the dimensions for all blocks. This data is stored in
 * static "cached" arrays.
 *
 * On every call (including the first), the function retrieves the dimensions for the
 * specific block (identified by `user->_this`) from the cached arrays and populates
 * the `IM`, `JM`, and `KM` fields of the provided `UserCtx`.
 *
 * @param user Pointer to the `UserCtx` for a specific block.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ReadGridFile(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    PetscInt block_index = user->_this;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- One-Time Read and Broadcast Logic ---
    if (!g_file_has_been_read) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "First call to ReadGridFile. Reading and broadcasting grid file header from '%s'...\n", simCtx->grid_file);
        PetscMPIInt rank = simCtx->rank;
        PetscInt    nblk = simCtx->block_number;

        if (rank == 0) {
            FILE *fd = fopen(simCtx->grid_file, "r");
            if (!fd) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", simCtx->grid_file);

            fscanf(fd, "%d\n", &g_nblk_from_file);
            // ---- Read first token; decide whether it is the header or nblk ----
            char firstTok[32] = {0};
            if (fscanf(fd, "%31s", firstTok) != 1)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Empty grid file: %s", simCtx->grid_file);

            if (strcmp(firstTok, "PICGRID") == 0) {
                // Header is present – read nblk from the next line
                if (fscanf(fd, "%d", &g_nblk_from_file) != 1)
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Expected number of blocks after \"PICGRID\" in %s", simCtx->grid_file);
            } else {
                // No header – the token we just read is actually nblk
                g_nblk_from_file = (PetscInt)strtol(firstTok, NULL, 10);
            }	    
            if (g_nblk_from_file != nblk) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED, "Mismatch: -nblk is %d but grid file specifies %d blocks.", nblk, g_nblk_from_file);
            }

            ierr = PetscMalloc3(nblk, &g_IMs_from_file, nblk, &g_JMs_from_file, nblk, &g_KMs_from_file); CHKERRQ(ierr);
            for (PetscInt i = 0; i < nblk; ++i) {
                if (fscanf(fd, "%d %d %d\n", &g_IMs_from_file[i], &g_JMs_from_file[i], &g_KMs_from_file[i]) != 3) {
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Expected 3 integers for block %d in %s", i, simCtx->grid_file);
                }
            }
            fclose(fd);
        }

        // Broadcast nblk to verify (optional, good practice)
        ierr = MPI_Bcast(&g_nblk_from_file, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        // Allocate on other ranks before receiving the broadcast
        if (rank != 0) {
            ierr = PetscMalloc3(nblk, &g_IMs_from_file, nblk, &g_JMs_from_file, nblk, &g_KMs_from_file); CHKERRQ(ierr);
        }
        
        // Broadcast the data arrays
        ierr = MPI_Bcast(g_IMs_from_file, nblk, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(g_JMs_from_file, nblk, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(g_KMs_from_file, nblk, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        g_file_has_been_read = PETSC_TRUE;
        LOG_ALLOW(GLOBAL, LOG_INFO, "Grid file header read and broadcast complete.\n");
    }

    // --- Per-Block Assignment Logic (runs on every call) ---
    user->IM = g_IMs_from_file[block_index];
    user->JM = g_JMs_from_file[block_index];
    user->KM = g_KMs_from_file[block_index];

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Set file inputs for Block %d: IM=%d, JM=%d, KM=%d\n",
              simCtx->rank, block_index, user->IM, user->JM, user->KM);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


//================================================================================
//
//                        PRIVATE HELPER FUNCTIONS
//
//================================================================================

/**
 * @brief Frees the memory allocated for a linked list of BC_Param structs.
 * @param head A pointer to the head of the linked list to be freed.
 */
void FreeBC_ParamList(BC_Param *head) {
    BC_Param *current = head;
    while (current != NULL) {
        BC_Param *next = current->next;
        PetscFree(current->key);
        PetscFree(current->value);
        PetscFree(current);
        current = next;
    }
}

/**
 * @brief Converts a string representation of a face to a BCFace enum.
 * @param str The input string (e.g., "-Xi", "+Zeta"). Case-insensitive.
 * @param[out] face_out The resulting BCFace enum.
 * @return 0 on success.
 */
PetscErrorCode StringToBCFace(const char* str, BCFace* face_out) {
    if      (strcasecmp(str, "-Xi")   == 0) *face_out = BC_FACE_NEG_X;
    else if (strcasecmp(str, "+Xi")   == 0) *face_out = BC_FACE_POS_X;
    else if (strcasecmp(str, "-Eta")  == 0) *face_out = BC_FACE_NEG_Y;
    else if (strcasecmp(str, "+Eta")  == 0) *face_out = BC_FACE_POS_Y;
    else if (strcasecmp(str, "-Zeta") == 0) *face_out = BC_FACE_NEG_Z;
    else if (strcasecmp(str, "+Zeta") == 0) *face_out = BC_FACE_POS_Z;
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown face specifier: %s", str);
    return 0;
}

/**
 * @brief Converts a string representation of a BC type to a BCType enum.
 * @param str The input string (e.g., "WALL", "INLET"). Case-insensitive.
 * @param[out] type_out The resulting BCType enum.
 * @return 0 on success.
 */
PetscErrorCode StringToBCType(const char* str, BCType* type_out) {
    if      (strcasecmp(str, "WALL")      == 0) *type_out = WALL;
    else if (strcasecmp(str, "SYMMETRY")  == 0) *type_out = SYMMETRY;
    else if (strcasecmp(str, "INLET")     == 0) *type_out = INLET;
    else if (strcasecmp(str, "OUTLET")    == 0) *type_out = OUTLET;
    else if (strcasecmp(str, "PERIODIC")  == 0) *type_out = PERIODIC;
    // ... add other BCTypes here ...
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown BC Type string: %s", str);
    return 0;
}

/**
 * @brief Converts a string representation of a handler to a BCHandlerType enum.
 * @param str The input string (e.g., "noslip", "constant_velocity"). Case-insensitive.
 * @param[out] handler_out The resulting BCHandlerType enum.
 * @return 0 on success.
 */
PetscErrorCode StringToBCHandlerType(const char* str, BCHandlerType* handler_out) {
    if      (strcasecmp(str, "noslip")              == 0) *handler_out = BC_HANDLER_WALL_NOSLIP;
    else if (strcasecmp(str, "constant_velocity")   == 0) *handler_out = BC_HANDLER_INLET_CONSTANT_VELOCITY;
    else if (strcasecmp(str, "conservation")        == 0) *handler_out = BC_HANDLER_OUTLET_CONSERVATION;
    else if (strcasecmp(str, "parabolic")           == 0) *handler_out = BC_HANDLER_INLET_PARABOLIC;
    else if (strcasecmp(str,"periodic")             == 0) *handler_out = BC_HANDLER_PERIODIC;
    // ... add other BCHandlerTypes here ...
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown BC Handler string: %s", str);
    return 0;
}

/**
 * @brief Validates that a specific handler is compatible with a general BC type.
 * @param type The general BCType.
 * @param handler The specific BCHandlerType.
 * @return 0 if compatible, error code otherwise.
 */
PetscErrorCode ValidateBCHandlerForBCType(BCType type, BCHandlerType handler) {
    switch (type) {
        case OUTLET:
	    if(handler != BC_HANDLER_OUTLET_CONSERVATION) return PETSC_ERR_ARG_WRONG;
	    break;
        case WALL:
            if (handler != BC_HANDLER_WALL_NOSLIP && handler != BC_HANDLER_WALL_MOVING) return PETSC_ERR_ARG_WRONG;
            break;
        case INLET:
            if (handler != BC_HANDLER_INLET_CONSTANT_VELOCITY && handler != BC_HANDLER_INLET_PARABOLIC) return PETSC_ERR_ARG_WRONG;
            break;
        // ... add other validation cases here ...
        default: break;
    }
    return 0; // Combination is valid
}

/**
 * @brief Searches a BC_Param linked list for a key and returns its value as a double.
 * @param params The head of the BC_Param linked list.
 * @param key The key to search for (case-insensitive).
 * @param[out] value_out The found value, converted to a PetscReal.
 * @param[out] found Set to PETSC_TRUE if the key was found, PETSC_FALSE otherwise.
 * @return 0 on success.
 */
PetscErrorCode GetBCParamReal(BC_Param *params, const char *key, PetscReal *value_out, PetscBool *found) {
    *found = PETSC_FALSE;
    *value_out = 0.0;
    if (!key) return 0; // No key to search for

    BC_Param *current = params;
    while (current) {
        if (strcasecmp(current->key, key) == 0) {
            *value_out = atof(current->value);
            *found = PETSC_TRUE;
            return 0; // Found it, we're done
        }
        current = current->next;
    }
    return 0; // It's not an error to not find the key.
}

//================================================================================
//
//                        PUBLIC PARSING FUNCTION
//
//================================================================================
#undef __FUNCT__
#define __FUNCT__ "ParseAllBoundaryConditions"
/**
 * @brief Parses the boundary conditions file to configure the type, handler, and
 *        any associated parameters for all 6 global faces of the domain.
 *
 * This function performs the following steps:
 * 1.  On MPI rank 0, it reads the specified configuration file line-by-line.
 * 2.  It parses each line for `<Face> <Type> <Handler> [param=value]...` format.
 * 3.  It validates the parsed strings and stores the configuration, including a
 *     linked list of parameters, in a temporary array.
 * 4.  It then serializes this configuration and broadcasts it to all other MPI ranks.
 * 5.  All ranks (including rank 0) then deserialize the broadcasted data to populate
 *     their local `user->boundary_faces` array identically.
 * 6.  It also sets legacy fields in UserCtx for compatibility with other modules.
 *
 * @param[in,out] user               The main UserCtx struct where the final configuration
 *                                   for all ranks will be stored.
 * @param[in]     bcs_input_filename The path to the boundary conditions configuration file.
 * @return PetscErrorCode 0 on success, error code on failure.
 */
PetscErrorCode ParseAllBoundaryConditions(UserCtx *user, const char *bcs_input_filename)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    // Temporary storage for rank 0 to build the configuration before broadcasting.
    BoundaryFaceConfig configs_rank0[6];
    
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        FILE *file;
        char line_buffer[1024];

        // Initialize the temporary config array with safe defaults on rank 0.
        for (int i = 0; i < 6; i++) {
            configs_rank0[i].face_id = (BCFace)i;
            configs_rank0[i].mathematical_type = WALL;
            configs_rank0[i].handler_type = BC_HANDLER_WALL_NOSLIP;
            configs_rank0[i].params = NULL;
            configs_rank0[i].handler = NULL; // Handler object is not created here.
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "Parsing BC configuration from '%s' on rank 0... \n", bcs_input_filename);
        file = fopen(bcs_input_filename, "r");
        if (!file) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Could not open BCs file '%s'.", bcs_input_filename);

        while (fgets(line_buffer, sizeof(line_buffer), file)) {
            char *current_pos = line_buffer;
            while (isspace((unsigned char)*current_pos)) current_pos++; // Skip leading whitespace
            if (*current_pos == '#' || *current_pos == '\0' || *current_pos == '\n' || *current_pos == '\r') continue;

            char *face_str = strtok(current_pos, " \t\n\r");
            char *type_str = strtok(NULL, " \t\n\r");
            char *handler_str = strtok(NULL, " \t\n\r");

            if (!face_str || !type_str || !handler_str) {
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Malformed line in bcs.dat, skipping: %s \n", line_buffer);
                continue;
            }

            BCFace      face_enum;
            BCType      type_enum;
            BCHandlerType handler_enum;
            const char* handler_name_for_log;

            // --- Convert strings to enums and validate ---
            ierr = StringToBCFace(face_str, &face_enum); CHKERRQ(ierr);
            ierr = StringToBCType(type_str, &type_enum); CHKERRQ(ierr);
            ierr = StringToBCHandlerType(handler_str, &handler_enum); CHKERRQ(ierr);
            ierr = ValidateBCHandlerForBCType(type_enum, handler_enum);
            if (ierr) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Validation failed: Handler '%s' is not valid for Type '%s' on Face '%s'.\n", handler_str, type_str, face_str);

            // Store the core types for the corresponding face
            configs_rank0[face_enum].mathematical_type = type_enum;
            configs_rank0[face_enum].handler_type = handler_enum;
            handler_name_for_log = BCHandlerTypeToString(handler_enum); // Assumes this utility exists
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "  Parsed Face '%s': Type=%s, Handler=%s \n", face_str, type_str, handler_name_for_log);

            // --- Parse optional key=value parameters for this face ---
            FreeBC_ParamList(configs_rank0[face_enum].params); // Clear any previous (default) params
            configs_rank0[face_enum].params = NULL;
            BC_Param **param_next_ptr = &configs_rank0[face_enum].params; // Pointer to the 'next' pointer to build the list

            char* token;
            while ((token = strtok(NULL, " \t\n\r")) != NULL) {
                char* equals_ptr = strchr(token, '=');
                if (!equals_ptr) {
                    LOG_ALLOW(GLOBAL, LOG_WARNING, "Malformed parameter '%s' on face '%s', skipping. \n", token, face_str);
                    continue;
                }
                
                *equals_ptr = '\0'; // Temporarily split the string at '=' to separate key and value
                char* key_str = token;
                char* value_str = equals_ptr + 1;

                BC_Param *new_param;
                ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                ierr = PetscStrallocpy(key_str, &new_param->key); CHKERRQ(ierr);
                ierr = PetscStrallocpy(value_str, &new_param->value); CHKERRQ(ierr);
                new_param->next = NULL;

                *param_next_ptr = new_param;
                param_next_ptr = &new_param->next;
                LOG_ALLOW(GLOBAL, LOG_TRACE, "    - Found param: [%s] = [%s] \n", new_param->key, new_param->value);
            }
        }
        fclose(file);
    }

    // =========================================================================
    //               BROADCASTING THE CONFIGURATION FROM RANK 0
    // =========================================================================
    // This is a critical step to ensure all processes have the same configuration.
    
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d broadcasting/receiving BC configuration.\n", rank);

    for (int i = 0; i < 6; i++) {
        // --- Broadcast simple enums ---
        if (rank == 0) {
            user->boundary_faces[i] = configs_rank0[i]; // Rank 0 populates its final struct
        }
        ierr = MPI_Bcast(&user->boundary_faces[i].mathematical_type, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(&user->boundary_faces[i].handler_type, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // --- Serialize and Broadcast the parameter linked list ---
        PetscInt n_params = 0;
        if (rank == 0) { // On rank 0, count the number of parameters to send
            for (BC_Param *p = user->boundary_faces[i].params; p; p = p->next) n_params++;
        }
        ierr = MPI_Bcast(&n_params, 1, MPI_INT, 0, PETSC_COMM_WORLD);CHKERRQ(ierr);
        
        if (rank != 0) { // Non-root ranks need to receive and build the list
            FreeBC_ParamList(user->boundary_faces[i].params); // Ensure list is empty before building
            user->boundary_faces[i].params = NULL;
        }

        BC_Param **param_next_ptr = &user->boundary_faces[i].params;

        for (int j = 0; j < n_params; j++) {
            char key_buf[256] = {0}, val_buf[256] = {0};
            if (rank == 0) {
                // On rank 0, navigate to the j-th param and copy its data to buffers
                BC_Param *p = user->boundary_faces[i].params;
                for (int k = 0; k < j; k++) p = p->next;
                strncpy(key_buf, p->key, 255);
                strncpy(val_buf, p->value, 255);
            }

            ierr = MPI_Bcast(key_buf, 256, MPI_CHAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Bcast(val_buf, 256, MPI_CHAR, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

            if (rank != 0) {
                // On non-root ranks, deserialize: create a new node and append it
                BC_Param *new_param;
                ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                ierr = PetscStrallocpy(key_buf, &new_param->key); CHKERRQ(ierr);
                ierr = PetscStrallocpy(val_buf, &new_param->value); CHKERRQ(ierr);
                new_param->next = NULL;
                *param_next_ptr = new_param;
                param_next_ptr = &new_param->next;
            } else {
                 // On rank 0, just advance the pointer for the next iteration
                 param_next_ptr = &((*param_next_ptr)->next);
            }
        }
        user->boundary_faces[i].face_id = (BCFace)i; // Ensure face_id is set on all ranks
    }
    
    // --- Set legacy fields for compatibility with particle system ---
    user->inletFaceDefined = PETSC_FALSE;
    for (int i=0; i<6; i++) {
        
        if (user->boundary_faces[i].mathematical_type == INLET && !user->inletFaceDefined) {
            user->inletFaceDefined = PETSC_TRUE;
            user->identifiedInletBCFace = (BCFace)i;
            LOG_ALLOW(GLOBAL, LOG_INFO, "Inlet face for particle initialization identified as Face %d.\n", i);
            break; // Found the first one, stop looking
        }
    }

    
    if (rank == 0) {
        // Rank 0 can now free the linked lists it created for the temporary storage.
        // As written, user->boundary_faces was populated directly on rank 0, so no extra free is needed.
        // for(int i=0; i<6; i++) FreeBC_ParamList(configs_rank0[i].params); // This would be needed if we used configs_rank0 exclusively
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief A parallel-safe helper to verify the existence of a generic file or directory path.
 *
 * This function centralizes the logic for checking arbitrary paths. Only Rank 0 performs the
 * filesystem check, and the result is broadcast to all other processes. This ensures
 * collective and synchronized decision-making across all ranks. It is intended for
 * configuration files, source directories, etc., where the path is known completely.
 *
 * @param[in]  path         The full path to the file or directory to check.
 * @param[in]  is_dir       PETSC_TRUE if checking for a directory, PETSC_FALSE for a file.
 * @param[in]  is_optional  PETSC_TRUE if the path is optional (results in a warning),
 *                          PETSC_FALSE if mandatory (results in an error).
 * @param[in]  description  A user-friendly description of the path for logging (e.g., "Grid file").
 * @param[out] exists       The result of the check (identical on all ranks).
 *
 * @return PetscErrorCode
 */
PetscErrorCode VerifyPathExistence(const char *path, PetscBool is_dir, PetscBool is_optional, const char *description, PetscBool *exists)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    MPI_Comm       comm = PETSC_COMM_WORLD;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        if (is_dir) {
            ierr = PetscTestDirectory(path, 'r', exists); CHKERRQ(ierr);
        } else {
            ierr = PetscTestFile(path, 'r', exists); CHKERRQ(ierr);
        }

        if (!(*exists)) {
            if (is_optional) {
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Optional %s not found at: %s (using defaults/ignoring).\n", description, path);
            } else {
                LOG_ALLOW(GLOBAL, LOG_ERROR, "Mandatory %s not found at: %s\n", description, path);
            }
        } else {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "Found %s: %s\n", description, path);
        }
    }

    // Broadcast the result from Rank 0
    PetscMPIInt exists_int = (rank == 0) ? (PetscMPIInt)(*exists) : 0;
    ierr = MPI_Bcast(&exists_int, 1, MPI_INT, 0, comm); CHKERRMPI(ierr);
    *exists = (PetscBool)exists_int;

    // Collective error for mandatory files
    if (!(*exists) && !is_optional) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Mandatory %s not found. Rank 0 expected it at '%s'. Check path and permissions.", description, path);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Checks for a data file's existence in a parallel-safe manner.
 * 
 * This function's signature and behavior have been updated to use the
 * "Directory Context Pointer" pattern. The high-level calling function (e.g.,
 * ReadSimulationFields) must set the `simCtx->current_io_directory` context
 * before this function is called.
 *
 * @param[in]  user Pointer to the UserCtx structure to access SimCtx.
 * @param ti The time index of the file.
 * @param fieldName The name of the field.
 * @param ext The file extension.
 * @param fileExists [out] The result, which will be identical on all ranks.
 * @return PetscErrorCode
 */
static PetscErrorCode CheckDataFile(UserCtx *user, PetscInt ti, const char *fieldName, const char *ext, PetscBool *fileExists)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    MPI_Comm       comm = PETSC_COMM_WORLD;
    PetscInt       placeholder_int = 0;
    SimCtx         *simCtx = user->simCtx;

    PetscFunctionBeginUser;

    if(!simCtx->current_io_directory) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "I/O context directory is NULL. Ensure it is set before calling CheckDataFile().");
    }

    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        char filename[PETSC_MAX_PATH_LEN];
        // Use the same standardized, rank-independent filename format

        ierr =  PetscSNPrintf(filename, sizeof(filename), "%s/%s%05"PetscInt_FMT"_%d.%s",simCtx->current_io_directory,fieldName, ti, placeholder_int, ext); 
        ierr = PetscTestFile(filename, 'r', fileExists); CHKERRQ(ierr);
        if (!(*fileExists)) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "(Rank 0) - Optional data file '%s' not found.\n", filename);
        }
    }

    // Broadcast the result from Rank 0 to all other ranks.
    // We cast the PetscBool to a PetscMPIInt for MPI_Bcast.
    PetscMPIInt fileExists_int = (rank == 0) ? (PetscMPIInt)(*fileExists) : 0;
    ierr = MPI_Bcast(&fileExists_int, 1, MPI_INT, 0, comm); CHKERRMPI(ierr);
    *fileExists = (PetscBool)fileExists_int;

    PetscFunctionReturn(0);
}

/**
 * @brief Checks for and reads an optional Eulerian field from a file into a Vec.
 *
 * This helper function first checks if the corresponding data file exists in a
 * parallel-safe manner using CheckDataFile().
 * - If the file does not exist, it logs a warning and returns success (0),
 *   gracefully skipping the field. The target Vec is not modified.
 * - If the file exists, it calls ReadFieldData() to read the data into the
 *   provided Vec. Any error during this read operation is considered fatal
 *   (e.g., file is corrupt, permissions issue, or size mismatch), as the
 *   presence of the file implies an intent to load it.
 *
 * @param[in] user        Pointer to the UserCtx structure.
 * @param[in] field_name  Internal name of the field for the filename (e.g., "cs").
 * @param[in] field_label A user-friendly name for logging (e.g., "Smagorinsky Constant").
 * @param[in,out] field_vec The PETSc Vec to load the data into.
 * @param[in] ti          Time index for the file name.
 * @param[in] ext         File extension.
 *
 * @return PetscErrorCode Returns 0 on success or if the optional file was not found.
 *         Returns a non-zero error code on a fatal read error.
 */
static PetscErrorCode ReadOptionalField(UserCtx *user, const char *field_name, const char *field_label, Vec field_vec, PetscInt ti, const char *ext)
{
  PetscErrorCode ierr;
  PetscBool      fileExists;

  PetscFunctionBeginUser;

  /* Check if the data file for this optional field exists. */
  ierr = CheckDataFile(user,ti, field_name, ext, &fileExists); CHKERRQ(ierr);

  if (fileExists) {
    /* File exists, so we MUST be able to read it. */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "File for %s found, attempting to read...\n", field_label);
    ierr = ReadFieldData(user, field_name, field_vec, ti, ext);

    if (ierr) {
      /* Any error here is fatal. A PETSC_ERR_FILE_OPEN would mean a race condition or
         permissions issue. Other errors could be size mismatch or corruption. */
      SETERRQ(PETSC_COMM_WORLD, ierr, "Failed to read data for %s from existing file for step %d. The file may be corrupt or have an incorrect size.", field_label, ti);
    } else {
      LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully read %s field for step %d.\n", field_label, ti);
    }
  } else {
    /* File does not exist, which is acceptable for an optional field. */
    LOG_ALLOW(GLOBAL, LOG_WARNING, "Optional %s file for step %d not found. Skipping.\n", field_label, ti);
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Checks for and reads an optional DMSwarm field from a file.
 *
 * This helper function first checks if the corresponding data file exists.
 * - If the file does not exist, it logs a warning and returns success (0),
 *   gracefully skipping the field.
 * - If the file exists, it calls ReadSwarmField() to read the data. Any
 *   error during this read operation is considered fatal (e.g., file is
 *   corrupt, permissions issue, or size mismatch), as the presence of the
 *   file implies an intent to load it.
 *
 * @param[in] user        Pointer to the UserCtx structure.
 * @param[in] field_name  Internal name of the DMSwarm field (e.g., "DMSwarm_CellID").
 * @param[in] field_label A user-friendly name for logging (e.g., "Cell ID").
 * @param[in] ti          Time index for the file name.
 * @param[in] ext         File extension.
 *
 * @return PetscErrorCode Returns 0 on success or if the optional file was not found.
 *         Returns a non-zero error code on a fatal read error.
 */
static PetscErrorCode ReadOptionalSwarmField(UserCtx *user, const char *field_name, const char *field_label, PetscInt ti, const char *ext)
{
  PetscErrorCode ierr;
  PetscBool      fileExists;

  PetscFunctionBeginUser;

  /* Check if the data file for this optional field exists. */
  ierr = CheckDataFile(user,ti, field_name, ext, &fileExists); CHKERRQ(ierr);

  if (fileExists) {
    /* File exists, so we MUST be able to read it. */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "File for %s found, attempting to read...\n", field_label);
    if(strcasecmp(field_name,"DMSwarm_CellID") == 0 || strcasecmp(field_name,"DMSwarm_pid")== 0 || strcasecmp(field_name,"DMSwarm_location_status")== 0 ) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Reading integer swarm field '%s'.\n", field_name);
        ierr = ReadSwarmIntField(user,field_name,ti,ext);
    }
    else{
    ierr = ReadSwarmField(user, field_name, ti, ext);
    } 

    if (ierr) {
      /* Any error here is fatal. A PETSC_ERR_FILE_OPEN would mean a race condition or
         permissions issue. Other errors could be size mismatch or corruption. */
      SETERRQ(PETSC_COMM_WORLD, ierr, "Failed to read data for %s from existing file for step %d. The file may be corrupt or have an incorrect size.", field_label, ti);
    } else {
      LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully read %s field for step %d.\n", field_label, ti);
    }
  } else {
    /* File does not exist, which is acceptable for an optional field. */
    LOG_ALLOW(GLOBAL, LOG_WARNING, "Optional %s file for step %d not found. Skipping.\n", field_label, ti);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReadFieldData"
/************************************************************************************************
*  @file   io.c
*  @brief  Utility for (re-)loading PETSc Vec-based fields written by rank-0 in a previous run
*          – works in both serial and MPI executions.
*
*  The restart files are always written as **one** sequential Vec per field
*  ( `<path-to-file>/<field_name><step>_0.<ext>` ).  In parallel we therefore
*  have to:
*
*     • let rank-0 read that sequential file into a temporary Vec,
*     • broadcast the raw array to every rank, and
*     • copy the local slice into the distributed Vec supplied by the caller.
*
*  The function purposely avoids `VecScatterCreateToAll()` because PETSc
*  chooses an even block-cyclic layout that seldom matches an existing
*  application layout – that mismatch is the root–cause of the previous
*  “Scatter sizes of ix and iy don’t match” / “incompatible local lengths”
*  crashes.  A plain `MPI_Bcast` plus a `memcpy()` into the already
*  allocated local array is bullet-proof and costs one extra buffer only
*  on the non-root ranks.
*
*  Logging helper used below:
*      LOG_ALLOW(scope,level,fmt,…)     – your existing macro
*
*  @param[in]  user        – application context (unused here but part of original API)
*  @param[in]  field_name  – base name of the field (“position”, “ufield”, …)
*  @param[out] field_vec   – **distributed** Vec we have to fill
*  @param[in]  ti          – time-step index encoded in the filename
*  @param[in]  ext         – file-name extension (“dat”, “bin”, …)
*
*  @returns 0 on success or a PETSc error code.
************************************************************************************************/
PetscErrorCode ReadFieldData(UserCtx *user,
                             const char *field_name,
                             Vec         field_vec,
                             PetscInt    ti,
                             const char *ext)
{
   PetscErrorCode ierr;
   char           filename[PETSC_MAX_PATH_LEN];
   MPI_Comm       comm;
   PetscMPIInt    rank,size;
   SimCtx         *simCtx = user->simCtx;


   PetscFunctionBeginUser;
   PROFILE_FUNCTION_BEGIN;

    if(!simCtx->current_io_directory){
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "I/O context directory was not set before calling ReadFieldData().");
    }


   ierr = PetscObjectGetComm((PetscObject)field_vec,&comm);CHKERRQ(ierr);
   ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
   ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

   const char *source_path = NULL;
   source_path = simCtx->current_io_directory;

   if(!source_path){
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "source_path was not set for the current execution mode.");
   }
   /* ---------------------------------------------------------------------
    * Compose <path-to-file>/<field_name><step with 5 digits>_0.<ext>
    * (all restart files are written by rank-0 with that naming scheme).
    * ------------------------------------------------------------------ */
   ierr = PetscSNPrintf(filename,sizeof(filename),
                        "%s/%s%05" PetscInt_FMT "_0.%s",
                        source_path,field_name,ti,ext);CHKERRQ(ierr);

   LOG_ALLOW(GLOBAL,LOG_DEBUG,
             "Attempting to read <%s> on rank %d/%d\n",
             filename,(int)rank,(int)size);

   /* ======================================================================
    * 1.  SERIAL JOB  – just hand the Vec to VecLoad()
    * ==================================================================== */
   if(size==1)
   {
      PetscViewer viewer;
      PetscBool   found;
      Vec temp_vec;
      PetscInt expectedSize,loadedSize;

      ierr = PetscTestFile(filename,'r',&found);CHKERRQ(ierr);
      if(!found) SETERRQ(comm,PETSC_ERR_FILE_OPEN,
                         "Restart/Source file not found: %s",filename);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
// ---- START MODIFICATION ----
      // DO NOT load directly into field_vec, as this can resize it, which is
      // illegal for DMSwarm "view" vectors. Instead, load into a temporary vector.
      ierr = VecCreate(PETSC_COMM_SELF, &temp_vec); CHKERRQ(ierr);
      ierr = VecLoad(temp_vec,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

      // Sanity check: ensure the file size matches the expected vector size.
      ierr = VecGetSize(field_vec, &expectedSize);CHKERRQ(ierr);
      ierr = VecGetSize(temp_vec, &loadedSize);CHKERRQ(ierr);
      if (loadedSize != expectedSize) {
         SETERRQ(comm,PETSC_ERR_FILE_UNEXPECTED,
                 "File %s holds %d entries – expected %d for field '%s'",
                 filename, loadedSize, expectedSize, field_name);
      }

      // Now, safely copy the data from the temporary vector to the final destination.
      ierr = VecCopy(temp_vec, field_vec);CHKERRQ(ierr);

      // Clean up the temporary vector.
      ierr = VecDestroy(&temp_vec);CHKERRQ(ierr);

      // ---- END MODIFICATION ----
      
      /* create EMPTY sequential Vec – VecLoad() will size it correctly   */
      /*
      ierr = VecCreate(PETSC_COMM_SELF,&seq_vec);CHKERRQ(ierr);
      ierr = VecSetType(seq_vec,VECSEQ);CHKERRQ(ierr);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,
                                   FILE_MODE_READ,&viewer);CHKERRQ(ierr);

      ierr = VecLoad(field_vec,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
      */
      LOG_ALLOW(GLOBAL,LOG_INFO,
                "Loaded <%s> (serial path)\n",filename);

      PROFILE_FUNCTION_END;          
      PetscFunctionReturn(0);
   }

   /* ======================================================================
    * 2.  PARALLEL JOB
    * ==================================================================== */
   PetscInt globalSize;
   ierr = VecGetSize(field_vec,&globalSize);CHKERRQ(ierr);

   /* -------------------- rank-0 : read the sequential file -------------- */
   Vec            seq_vec = NULL;      /* only valid on rank-0            */
   const PetscScalar *seqArray = NULL; /* borrowed pointer on rank-0 only */

   if(rank==0)
   {
      PetscViewer viewer;
      PetscBool   found;

      ierr = PetscTestFile(filename,'r',&found);CHKERRQ(ierr);
      if(!found) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,
                         "Restart file not found: %s",filename);

      /* create EMPTY sequential Vec – VecLoad() will size it correctly   */
      ierr = VecCreate(PETSC_COMM_SELF,&seq_vec);CHKERRQ(ierr);
      ierr = VecSetType(seq_vec,VECSEQ);CHKERRQ(ierr);

      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,
                                   FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(seq_vec,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

      /* size sanity-check */
      PetscInt loaded;
      ierr = VecGetSize(seq_vec,&loaded);CHKERRQ(ierr);
      if(loaded != globalSize)
         SETERRQ(comm,PETSC_ERR_FILE_UNEXPECTED,
                 "File %s holds %d entries – expected %d",
                 filename,loaded,globalSize);

      /* borrow array for later Bcast */
      ierr = VecGetArrayRead(seq_vec,&seqArray);CHKERRQ(ierr);

      LOG_ALLOW(GLOBAL,LOG_TRACE,
                "Rank 0 successfully loaded <%s>\n",filename);
   }

   /* -------------------- everybody : broadcast raw data ----------------- */
   PetscScalar *buffer = NULL;               /* receives the full field    */
   if(rank==0)
   {
      /* shallow-copy: const-cast is safe, we do not modify the data        */
      buffer = (PetscScalar *)seqArray;
   }
   else
   { /* non-root ranks allocate a receive buffer                           */
      ierr = PetscMalloc1(globalSize,&buffer);CHKERRQ(ierr);
   }

   ierr = MPI_Bcast(buffer, (int)globalSize, MPIU_SCALAR, 0, comm);CHKERRQ(ierr);

   /* -------------------- copy my slice into field_vec ------------------- */
   PetscInt  rstart,rend,loc;
   PetscScalar *locArray;

   ierr = VecGetOwnershipRange(field_vec,&rstart,&rend);CHKERRQ(ierr);
   loc  = rend - rstart;                    /* local length               */

   ierr = VecGetArray(field_vec,&locArray);CHKERRQ(ierr);
   ierr = PetscMemcpy(locArray,
                      buffer + rstart,
                      loc*sizeof(PetscScalar));CHKERRQ(ierr);
   ierr = VecRestoreArray(field_vec,&locArray);CHKERRQ(ierr);

   /* -------------------- tidy up ---------------------------------------- */
   if(rank==0)
   {
      ierr = VecRestoreArrayRead(seq_vec,&seqArray);CHKERRQ(ierr);
      ierr = VecDestroy(&seq_vec);CHKERRQ(ierr);
   }
   else
   {
      ierr = PetscFree(buffer);CHKERRQ(ierr);
   }

   LOG_ALLOW(GLOBAL,LOG_INFO,
             "Loaded <%s> (parallel path)\n",filename);

   PROFILE_FUNCTION_END;          
   PetscFunctionReturn(0);
}


/**
 * @brief Reads simulation fields from files into their respective PETSc vectors.
 *
 * This function reads contravariant velocity, Cartesian velocity, pressure, and node state
 * fields from their respective binary files. It also conditionally reads LES, RANS, and
 * statistical fields if they are enabled.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadSimulationFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    SimCtx *simCtx = user->simCtx;
    const char *source_path = NULL;

    if(simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR){
        source_path = simCtx->pps->source_dir;
    } else if(simCtx->exec_mode == EXEC_MODE_SOLVER){
        source_path = simCtx->restart_dir;
    } else{
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Invalid execution mode for reading simulation fields.");
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to read simulation fields.\n");

    // Set the current I/O directory context
    ierr  = PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer), "%s/%s",source_path,simCtx->euler_subdir); CHKERRQ(ierr);
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    // Read Cartesian velocity field
    ierr = ReadFieldData(user, "ufield", user->Ucat, ti, "dat"); CHKERRQ(ierr);

    // Read contravariant velocity field
    ierr = ReadFieldData(user, "vfield", user->Ucont, ti, "dat"); CHKERRQ(ierr);

    // Read pressure field
    ierr = ReadFieldData(user, "pfield", user->P, ti, "dat"); CHKERRQ(ierr);

    // Read node state field (nvert)
    ierr = ReadFieldData(user, "nvfield", user->Nvert, ti, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_INFO,"Successfully read all mandatory fields. \n");

    if(simCtx->np>0){    
    // Read Particle Count field
    if(!user->ParticleCount){
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "ParticleCount Vec is NULL but np>0");
    }
    ierr = ReadOptionalField(user, "ParticleCount", "Particle Count", user->ParticleCount, ti, "dat"); CHKERRQ(ierr);
    
    ierr = ReadOptionalField(user, "psifield", "Scalar Psi Field", user->Psi, ti, "dat"); CHKERRQ(ierr);
    }
    else{
        LOG_ALLOW(GLOBAL, LOG_INFO, "No particles in simulation, skipping Particle fields reading.\n");
    }
    // Process LES fields if enabled
    if (simCtx->les) {
      ierr = ReadLESFields(user,ti); CHKERRQ(ierr);
    }

    // Process RANS fields if enabled
    if (simCtx->rans) {
      ierr = ReadRANSFields(user,ti); CHKERRQ(ierr);
    }

    // Process statistical fields if averaging is enabled
    if (simCtx->averaging) {
      ierr = ReadStatisticalFields(user,ti); CHKERRQ(ierr);
    }

    simCtx->current_io_directory = NULL; // Clear the I/O context after reading

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished reading simulation fields.\n");

    return 0;
}

/**
 * @brief Reads statistical fields for averaging purposes.
 *
 * This function reads data for fields such as Ucat_sum, Ucat_cross_sum, Ucat_square_sum,
 * and P_sum, used for statistical analysis during simulation.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadStatisticalFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to read statistical fields.\n");

    ierr = ReadOptionalField(user, "su0", "Velocity Sum",       user->Ucat_sum,        ti, "dat"); CHKERRQ(ierr);
    ierr = ReadOptionalField(user, "su1", "Velocity Cross Sum", user->Ucat_cross_sum,  ti, "dat"); CHKERRQ(ierr);
    ierr = ReadOptionalField(user, "su2", "Velocity Square Sum",user->Ucat_square_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadOptionalField(user, "sp",  "Pressure Sum",       user->P_sum,           ti, "dat"); CHKERRQ(ierr);    

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished reading statistical fields.\n");

    return 0;
}

/**
 * @brief Reads LES-related fields.
 *
 * This function reads LES-related fields such as Cs (Smagorinsky constant)
 * into their respective PETSc vectors.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadLESFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to read LES fields.\n");

    ierr = ReadOptionalField(user, "Nu_t", "Turbulent Viscosity", user->Nu_t, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadOptionalField(user, "cs", "Smagorinsky Constant (Cs)", user->CS, ti, "dat"); CHKERRQ(ierr);

    DMGlobalToLocalBegin(user->da, user->CS, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, user->CS, INSERT_VALUES, user->lCs);
    
    DMGlobalToLocalBegin(user->da, user->Nu_t, INSERT_VALUES, user->lNu_t);
    DMGlobalToLocalEnd(user->da, user->Nu_t, INSERT_VALUES, user->lNu_t);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished reading LES fields.\n");

    return 0;
}

/**
 * @brief Reads RANS-related fields.
 *
 * This function reads RANS-related fields such as K_Omega into their respective
 * PETSc vectors.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadRANSFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to read RANS fields.\n");

    ierr = ReadOptionalField(user, "kfield", "K-Omega RANS", user->K_Omega, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadOptionalField(user, "Nu_t", "Turbulent Viscosity", user->Nu_t, ti, "dat"); CHKERRQ(ierr);

    VecCopy(user->K_Omega, user->K_Omega_o);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);

    DMGlobalToLocalBegin(user->da, user->Nu_t, INSERT_VALUES, user->lNu_t);
    DMGlobalToLocalEnd(user->da, user->Nu_t, INSERT_VALUES, user->lNu_t);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished reading RANS fields.\n");

    return 0;
}


/**
 * @brief Reads data from a file into a specified field of a PETSc DMSwarm.
 *
 * This function is the counterpart to WriteSwarmField(). It creates a global PETSc vector 
 * that references the specified DMSwarm field, uses ReadFieldData() to read the data from 
 * a file, and then destroys the global vector reference.
 *
 * @param[in]  user       Pointer to the UserCtx structure (containing `user->swarm`).
 * @param[in]  field_name Name of the DMSwarm field to read into (must be previously declared/allocated).
 * @param[in]  ti         Time index used to construct the input file name.
 * @param[in]  ext        File extension (e.g., "dat" or "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note Compatible with PETSc 3.14.x.
 */
PetscErrorCode ReadSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext)
{
  PetscErrorCode ierr;
  DM             swarm;
  Vec            fieldVec;

  PetscFunctionBegin;

  swarm = user->swarm;

  LOG_ALLOW(GLOBAL,LOG_DEBUG," ReadSwarmField Begins \n");
 
  /* 2) Create a global vector that references the specified Swarm field. */
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," Vector created from Field \n");

  /* 3) Use the ReadFieldData() function to read data into fieldVec. */
  ierr = ReadFieldData(user, field_name, fieldVec, ti, ext);CHKERRQ(ierr);

  /* 4) Destroy the global vector reference. */
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Reads integer swarm data by using ReadFieldData and casting the result.
 *
 * This function is the counterpart to WriteSwarmIntField. It reads a file
 * containing floating-point data (that was originally integer) into a temporary
 * Vec and then casts it back to the integer swarm field. It works by:
 * 1. Creating a temporary parallel Vec.
 * 2. Calling the standard ReadFieldData() to populate this Vec.
 * 3. Accessing the local data of both the Vec and the swarm field.
 * 4. Populating the swarm's integer field by casting each PetscScalar back to a PetscInt.
 * 5. Destroying the temporary Vec.
 *
 * @param[in] user       Pointer to the UserCtx structure.
 * @param[in] field_name Name of the integer Swarm field to be read.
 * @param[in] ti         Time index for the input file.
 * @param[in] ext        File extension.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadSwarmIntField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    Vec            temp_vec;
    PetscInt       nlocal, nglobal, bs, i;
    const PetscScalar *scalar_array; // Read-only pointer from the temp Vec
    void           *field_array_void;


    PetscFunctionBeginUser;
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Reading '%s' via temporary Vec.\n", field_name);

    // Get the properties of the swarm field to determine the expected layout
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &nglobal); CHKERRQ(ierr);
    // We get the block size but not the data pointer yet
    ierr = DMSwarmGetField(swarm, field_name, &bs, NULL, NULL); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, field_name, &bs, NULL, NULL); CHKERRQ(ierr);

    // Create a temporary Vec with the CORRECT layout to receive the data
    ierr = VecCreate(PETSC_COMM_WORLD, &temp_vec); CHKERRQ(ierr);
    ierr = VecSetType(temp_vec, VECMPI); CHKERRQ(ierr);
    ierr = VecSetSizes(temp_vec, nlocal * bs, nglobal * bs); CHKERRQ(ierr);
    ierr = VecSetBlockSize(temp_vec, bs); CHKERRQ(ierr);
    ierr = VecSetUp(temp_vec); CHKERRQ(ierr);

    // Call your existing reader to populate the temporary Vec
    ierr = ReadFieldData(user, field_name, temp_vec, ti, ext); CHKERRQ(ierr);

    // Get local pointers
    ierr = VecGetArrayRead(temp_vec, &scalar_array); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, NULL, NULL, &field_array_void); CHKERRQ(ierr);
    
    // Perform the cast back, using the correct loop size (nlocal * bs)
    if (strcmp(field_name, "DMSwarm_pid") == 0) {
        PetscInt64 *int64_array = (PetscInt64 *)field_array_void;
        for (i = 0; i < nlocal * bs; i++) {
            int64_array[i] = (PetscInt64)scalar_array[i];
        }
    } else {
        PetscInt *int_array = (PetscInt *)field_array_void;
        for (i = 0; i < nlocal * bs; i++) {
            int_array[i] = (PetscInt)scalar_array[i];
        }
    }    

    // Restore access
    ierr = DMSwarmRestoreField(swarm, field_name, NULL, NULL, &field_array_void); CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(temp_vec, &scalar_array); CHKERRQ(ierr);

    // 6. Clean up
    ierr = VecDestroy(&temp_vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Reads multiple fields into a DMSwarm, gracefully handling missing optional files.
 *
 * This function reads all necessary and optional fields for a DMSwarm for a given
 * timestep. It assumes the swarm has already been resized to match the particle
 * count in the input files.
 *
 * The 'position' field is considered MANDATORY. If its file is missing or corrupt,
 * the function will return a fatal error.
 *
 * All other fields (velocity, CellID, weight, etc.) are OPTIONAL. If their
 * corresponding files are not found, a warning is logged, and the function
 * continues without error. If an optional file IS found but is corrupt or has a
 * size mismatch, it is treated as a fatal error.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm (user->swarm).
 * @param[in]     ti   Time index for constructing the file names.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadAllSwarmFields(UserCtx *user, PetscInt ti)
{
  PetscErrorCode ierr;
  PetscInt nGlobal;
  SimCtx *simCtx = user->simCtx;
  const char *source_path = NULL;

  PetscFunctionBeginUser;
  ierr = DMSwarmGetSize(user->swarm, &nGlobal); CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO, "Reading DMSwarm fields for timestep %d (swarm size is %d).\n", ti, nGlobal);

  if (nGlobal == 0) {
      LOG_ALLOW(GLOBAL, LOG_INFO, "Swarm is empty for timestep %d. Nothing to read.\n", ti);
      PetscFunctionReturn(0);
  }

    // First, determine the top-level source directory based on the execution mode.
    if (simCtx->exec_mode == EXEC_MODE_SOLVER) {
        source_path = simCtx->restart_dir;
    } else if (simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR) {
        source_path = simCtx->pps->source_dir;
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Invalid execution mode for reading simulation fields.");
    }

    // Set the current I/O directory context
    ierr = PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                         "%s/%s", source_path, simCtx->particle_subdir); CHKERRQ(ierr);

    simCtx->current_io_directory = simCtx->_io_context_buffer;

  /* 1) Read positions (REQUIRED) */
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Reading mandatory position field...\n");
  ierr = ReadSwarmField(user, "position", ti, "dat");
  if (ierr) {
      SETERRQ(PETSC_COMM_WORLD, ierr, "Failed to read MANDATORY 'position' field for step %d. Cannot continue.", ti);
  }
  LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully read mandatory position field for step %d.\n", ti);

  /* 2) Read all OPTIONAL fields using the helper function. */
  /* The helper will print a warning and continue if a file is not found. */
  ierr = ReadOptionalSwarmField(user, "velocity",                "Velocity",                   ti, "dat"); CHKERRQ(ierr);
  ierr = ReadOptionalSwarmField(user, "DMSwarm_pid",             "Particle ID",                ti, "dat"); CHKERRQ(ierr);
  ierr = ReadOptionalSwarmField(user, "DMSwarm_CellID",          "Cell ID",                    ti, "dat"); CHKERRQ(ierr);
  ierr = ReadOptionalSwarmField(user, "weight",                  "Particle Weight",            ti, "dat"); CHKERRQ(ierr);
  ierr = ReadOptionalSwarmField(user, "Psi",                     "Scalar Psi",                 ti, "dat"); CHKERRQ(ierr);
  ierr = ReadOptionalSwarmField(user, "DMSwarm_location_status", "Migration Status",   ti, "dat"); CHKERRQ(ierr);

  simCtx->current_io_directory = NULL; // Clear the I/O context after reading

  LOG_ALLOW(GLOBAL, LOG_INFO, "Finished reading DMSwarm fields for timestep %d.\n", ti);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "WriteFieldData" 
 /**
 * @brief Writes data from a specific PETSc vector to a single, sequential file.
 *
 * This function is now parallel-safe.
 * - In PARALLEL: All processes send their local data to Rank 0. Rank 0 assembles
 *   the data into a temporary sequential vector and writes it to a single file.
 * - In SERIAL: It performs a direct, simple write.
 *
 * This ensures the output file is always in a simple, portable format.
 *
 * @param[in] user       Pointer to the UserCtx structure.
 * @param[in] field_name Name of the field (e.g., "position").
 * @param[in] field_vec  The parallel PETSc vector containing the data to write.
 * @param[in] ti         Time index for constructing the file name.
 * @param[in] ext        File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
/**
 * @brief Dump a **distributed** PETSc Vec to the *single* sequential file
 *        format used by our restart / post-processing tools.
 *
 * The companion of ReadFieldData(): it always produces **one** file
 * (e.g. `results/ufield00006_0.dat`) regardless of how many MPI ranks
 * are running.
 *
 * **Behaviour**
 *
 * | # MPI ranks | Strategy                                                                  |
 * |-------------|---------------------------------------------------------------------------|
 * |      1      | Direct `VecView()` into the file.                                          |
 * |    > 1      | `VecScatterCreateToZero()` gathers the distributed Vec onto **rank 0**.    |
 * |             | Rank 0 writes the sequential Vec; all other ranks allocate no storage.     |
 *
 * The routine never alters or destroys the parallel Vec passed in; the
 * gather buffer is created and freed internally.
 *
 * @param[in] user        Simulation context (used only for logging).
 * @param[in] field_name  Logical field name (forms part of the filename).
 * @param[in] field_vec   Distributed PETSc Vec to write.
 * @param[in] ti          Timestep index used in the filename.
 * @param[in] ext         File extension (`"dat"` in our workflow).
 *
 * @return `0` on success or a PETSc error code.
 */
PetscErrorCode WriteFieldData(UserCtx *user,
                              const char *field_name,
                              Vec field_vec,
                              PetscInt ti,
                              const char *ext)
{
    PetscErrorCode ierr;
    MPI_Comm       comm;
    PetscMPIInt    rank, size;

    const PetscInt placeholder_int = 0;                    /* keep legacy name */
    char           filename[PETSC_MAX_PATH_LEN];
    SimCtx  *simCtx=user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if(!simCtx->output_dir){
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE, "Output directory was not set before calling WriteFieldData().");
    }

    /* ------------------------------------------------------------ */
    /*                  Basic communicator information              */
    /* ------------------------------------------------------------ */
    ierr = PetscObjectGetComm((PetscObject)field_vec,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

    ierr = PetscSNPrintf(filename,sizeof(filename),
                         "%s/%s%05" PetscInt_FMT "_%d.%s",
                         simCtx->current_io_directory,field_name,ti,placeholder_int,ext);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_DEBUG,
              " Preparing to write <%s> on rank %d/%d\n",
              filename,rank,size);

    /* ------------------------------------------------------------ */
    /*                       1.  Serial path                        */
    /* ------------------------------------------------------------ */
    if (size == 1) {
        PetscViewer viewer;

        ierr = PetscViewerBinaryOpen(comm,filename,
                                     FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = VecView(field_vec,viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_INFO,
                  " Wrote <%s> (serial path)\n",filename);
        PetscFunctionReturn(0);
    }

    /* ------------------------------------------------------------ */
    /*                      2.  Parallel path                       */
    /* ------------------------------------------------------------ */
    VecScatter scatter;
    Vec        seq_vec=NULL;               /* created by PETSc, lives only on rank 0 */

    /* 2.1  Create gather context and buffer                        */
    ierr = VecScatterCreateToZero(field_vec,&scatter,&seq_vec);CHKERRQ(ierr);

    /* 2.2  Gather distributed → sequential (on rank 0)             */
    ierr = VecScatterBegin(scatter,field_vec,seq_vec,
                           INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd  (scatter,field_vec,seq_vec,
                           INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

    /* 2.3  Rank 0 writes the file                                  */
    if (rank == 0) {
        PetscViewer viewer;

        /* (optional) value diagnostics */
        PetscReal vmin,vmax;
        ierr = VecMin(seq_vec,NULL,&vmin);CHKERRQ(ierr);
        ierr = VecMax(seq_vec,NULL,&vmax);CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL,LOG_DEBUG,
                  " <%s> range = [%.4e … %.4e]\n",
                  field_name,(double)vmin,(double)vmax);

        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,
                                     FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = VecView(seq_vec,viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_INFO,
                  " Wrote <%s> (parallel path)\n",filename);
    }

    /* 2.4  Cleanup                                                 */
    ierr = VecScatterDestroy(&scatter);CHKERRQ(ierr);
    ierr = VecDestroy(&seq_vec);CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Writes simulation fields to files.
 *
 * This function writes contravariant velocity, Cartesian velocity, pressure, and node state
 * fields to their respective binary files. It also conditionally writes LES, RANS, and
 * statistical fields if they are enabled.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteSimulationFields(UserCtx *user)
{
    PetscErrorCode ierr;

    SimCtx *simCtx = user->simCtx;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to write simulation fields.\n");

    // Set the current IO directory
    ierr = PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                         "%s/%s", simCtx->output_dir, simCtx->euler_subdir);CHKERRQ(ierr);
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    // Write contravariant velocity field
    ierr = WriteFieldData(user, "vfield", user->Ucont, simCtx->step, "dat"); CHKERRQ(ierr);

    // Write Cartesian velocity field
    ierr = WriteFieldData(user, "ufield", user->Ucat, simCtx->step, "dat"); CHKERRQ(ierr);

    // Write pressure field
    ierr = WriteFieldData(user, "pfield", user->P, simCtx->step, "dat"); CHKERRQ(ierr);

    // Write node state field (nvert)
    ierr = WriteFieldData(user, "nvfield", user->Nvert, simCtx->step, "dat"); CHKERRQ(ierr);

    // Write ParticleCountPerCell if enabled.
    if(simCtx->np>0){
    ierr = WriteFieldData(user, "ParticleCount",user->ParticleCount,simCtx->step,"dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "psifield", user->Psi, simCtx->step, "dat"); CHKERRQ(ierr);
    }
    
    // Write LES fields if enabled
    if (simCtx->les) {
        ierr = WriteLESFields(user); CHKERRQ(ierr);
    }

    // Write RANS fields if enabled
    if (simCtx->rans) {
        ierr = WriteRANSFields(user); CHKERRQ(ierr);
    }

    // Write statistical fields if averaging is enabled
    if (simCtx->averaging) {
        ierr = WriteStatisticalFields(user); CHKERRQ(ierr);
    }

    simCtx->current_io_directory = NULL;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished writing simulation fields.\n");

    return 0;
}

/**
 * @brief Writes statistical fields for averaging purposes.
 *
 * This function writes data for fields such as Ucat_sum, Ucat_cross_sum, Ucat_square_sum,
 * and P_sum to their respective binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteStatisticalFields(UserCtx *user)
{
    PetscErrorCode ierr;

    SimCtx *simCtx = user->simCtx;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to write statistical fields.\n");

    ierr = WriteFieldData(user, "su0", user->Ucat_sum, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su1", user->Ucat_cross_sum, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su2", user->Ucat_square_sum, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "sp", user->P_sum, simCtx->step, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished writing statistical fields.\n");

    return 0;
}

/**
 * @brief Writes LES-related fields.
 *
 * This function writes LES-related fields such as Cs (Smagorinsky constant)
 * to their respective binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteLESFields(UserCtx *user)
{
    PetscErrorCode ierr;

    SimCtx *simCtx = user->simCtx;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to write LES fields.\n");


    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, user->CS);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, user->CS);

    DMLocalToGlobalBegin(user->da, user->lNu_t, INSERT_VALUES, user->Nu_t);
    DMLocalToGlobalEnd(user->da, user->lNu_t, INSERT_VALUES, user->Nu_t);

    ierr = WriteFieldData(user, "Nu_t", user->Nu_t, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "cs", user->CS, simCtx->step, "dat"); CHKERRQ(ierr);
   

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished writing LES fields.\n");

    return 0;
}

/**
 * @brief Writes RANS-related fields.
 *
 * This function writes RANS-related fields such as K_Omega to their respective
 * binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteRANSFields(UserCtx *user)
{
    PetscErrorCode ierr;

    SimCtx *simCtx = user->simCtx;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to write RANS fields.\n");

    ierr = WriteFieldData(user, "kfield", user->K_Omega, simCtx->step, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished writing RANS fields.\n");

    return 0;
}

/**
 * @brief Writes data from a specific field in a PETSc Swarm to a file.
 *
 * This function retrieves the Swarm from the UserCtx (i.e., `user->swarm`) and
 * creates a global PETSc vector from the specified Swarm field. It then calls
 * the existing WriteFieldData() function to handle the actual I/O operation.
 * After writing the data, the function destroys the temporary global vector 
 * to avoid memory leaks.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context
 *                       and the PetscSwarm (as `user->swarm`).
 * @param[in] field_name Name of the Swarm field to be written (e.g., "my_field").
 * @param[in] ti         Time index used to construct the output file name.
 * @param[in] ext        File extension (e.g., "dat", "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note Compatible with PETSc 3.14.4 and newer.
 */
PetscErrorCode WriteSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext)
{
  PetscErrorCode ierr;
  Vec            fieldVec;
  DM     swarm;  

  PetscFunctionBeginUser; /* PETSc macro indicating start of function */

  /* 
   * 1) Retrieve the PetscSwarm from the user context.
   *    Ensure user->swarm is initialized and not NULL.
   */
  swarm = user->swarm;      

  /* 
   * 2) Create a global vector from the specified swarm field.
   *    This function is available in PETSc 3.14.4.
   *    It provides a read/write "view" of the swarm field as a global Vec.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Attempting to create global vector from field: %s\n",
            field_name);
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /*
   * 3) Use your existing WriteFieldData() to write the global vector to a file.
   *    The field name, time index, and extension are passed along for naming.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Calling WriteFieldData for field: %s\n",
            field_name);
  ierr = WriteFieldData(user, field_name, fieldVec, ti, ext);CHKERRQ(ierr);

  /*
   * 4) Destroy the global vector once the data is successfully written.
   *    This step is crucial for avoiding memory leaks. 
   *    DMSwarmDestroyGlobalVectorFromField() is also available in PETSc 3.14.4.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Destroying the global vector for field: %s\n",
            field_name);
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /* Log and return success. */
  LOG_ALLOW(GLOBAL, LOG_INFO,
            "Successfully wrote swarm data for field: %s\n",
            field_name);

  PetscFunctionReturn(0); /* PETSc macro indicating end of function */
}

/**
 * @brief Writes integer swarm data by casting it to a temporary Vec and using WriteFieldData.
 *
 * This function provides a bridge to write integer-based swarm fields (like DMSwarm_CellID)
 * using the existing Vec-based I/O routine (WriteFieldData). It works by:
 * 1. Creating a temporary parallel Vec with the same layout as other swarm fields.
 * 2. Accessing the local integer data from the swarm field.
 * 3. Populating the temporary Vec by casting each integer to a PetscScalar.
 * 4. Calling the standard WriteFieldData() function with the temporary Vec.
 * 5. Destroying the temporary Vec.
 *
 * @param[in] user       Pointer to the UserCtx structure.
 * @param[in] field_name Name of the integer Swarm field to be written.
 * @param[in] ti         Time index for the output file.
 * @param[in] ext        File extension.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteSwarmIntField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    Vec            temp_vec;       // Temporary Vec to hold casted data
    PetscInt       nlocal, nglobal,bs,i;
    void           *field_array_void;
    PetscScalar    *scalar_array;  // Pointer to the temporary Vec's scalar data

    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Casting '%s' to Vec for writing.\n", field_name);

    // Get the swarm field properties
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &nglobal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, &bs, NULL, &field_array_void); CHKERRQ(ierr);

    // Create Temporary parallel Vec wit the CORRECT layout
    ierr = VecCreate(PETSC_COMM_WORLD, &temp_vec); CHKERRQ(ierr);
    ierr = VecSetType(temp_vec, VECMPI); CHKERRQ(ierr);
    ierr = VecSetSizes(temp_vec, nlocal*bs, nglobal*bs); CHKERRQ(ierr);
    ierr = VecSetUp(temp_vec); CHKERRQ(ierr);

    // Defining Vector field to mandatory field 'position'
    DMSwarmVectorDefineField(swarm,"position");
              
    ierr = VecGetArray(temp_vec, &scalar_array); CHKERRQ(ierr);

    if(strcasecmp(field_name,"DMSwarm_pid") == 0){
        PetscInt64 *int64_array = (PetscInt64 *)field_array_void;
        // Perform the cast from PetscInt64 to PetscScalar
        for (i = 0; i < nlocal*bs; i++) {
            scalar_array[i] = (PetscScalar)int64_array[i];
        }
    }else{
        PetscInt *int_array = (PetscInt *)field_array_void;
        //Perform the cast from PetscInt to PetscScalar
        for (i = 0; i < nlocal*bs; i++) {
            scalar_array[i] = (PetscScalar)int_array[i];
        }
    }

    // Restore access to both arrays
    ierr = VecRestoreArray(temp_vec, &scalar_array); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, field_name, &bs, NULL, &field_array_void); CHKERRQ(ierr);

    // Call your existing writer with the temporary, populated Vec
    ierr = WriteFieldData(user, field_name, temp_vec, ti, ext); CHKERRQ(ierr);

    // Clean up
    ierr = VecDestroy(&temp_vec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Writes a predefined set of PETSc Swarm fields to files.
 *
 * This function iterates through a hardcoded list of common swarm fields 
 * (position, velocity, etc.) and calls the WriteSwarmField() helper function 
 * for each one. This provides a straightforward way to output essential particle 
 * data at a given simulation step.
 *
 * This function will only execute if particles are enabled in the simulation
 * (i.e., `user->simCtx->np > 0` and `user->swarm` is not NULL).
 *
 * @param[in] user Pointer to the UserCtx structure containing the simulation context
 *                 and the PetscSwarm.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteAllSwarmFields(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    
    PetscFunctionBeginUser;

    // If no swarm is configured or there are no particles, do nothing and return.
    if (!user->swarm || simCtx->np <= 0) {
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting to write swarm fields.\n");

    // Ensure the current IO directory is set
    ierr = PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                         "%s/%s", simCtx->output_dir, simCtx->particle_subdir);CHKERRQ(ierr);
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    // Write particle position field
    ierr = WriteSwarmField(user, "position", simCtx->step, "dat"); CHKERRQ(ierr);

    // Write particle velocity field
    ierr = WriteSwarmField(user, "velocity", simCtx->step, "dat"); CHKERRQ(ierr);

    // Write particle weight field
    ierr = WriteSwarmField(user, "weight", simCtx->step, "dat"); CHKERRQ(ierr);
    
    // Write custom particle field "Psi"
    ierr = WriteSwarmField(user, "Psi", simCtx->step, "dat"); CHKERRQ(ierr);
    
    // Integer fields require special handling

    // Write the background mesh cell ID for each particle
    ierr = WriteSwarmIntField(user, "DMSwarm_CellID", simCtx->step, "dat"); CHKERRQ(ierr);

    // Write the particle location status (e.g., inside or outside the domain)
    ierr = WriteSwarmIntField(user, "DMSwarm_location_status", simCtx->step, "dat"); CHKERRQ(ierr);

    // Write the unique particle ID
    ierr = WriteSwarmIntField(user, "DMSwarm_pid", simCtx->step, "dat"); CHKERRQ(ierr);

    simCtx->current_io_directory = NULL;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished writing swarm fields.\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Gather a (possibly distributed) PETSc Vec onto rank 0 as a contiguous C array.
 *        If the Vec has a DMDA attached, the gather is performed in DMDA "natural" ordering.
 *
 * @param[in]  inVec     The PETSc Vec (may be distributed).
 * @param[out] N         Global length of the Vec (includes dof).
 * @param[out] arrayOut  On rank 0, newly allocated buffer with the gathered values (PetscScalar-sized).
 *                       On other ranks, set to NULL.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut)
{
    PetscErrorCode ierr;
    MPI_Comm       comm;
    PetscMPIInt    rank;
    PetscInt       globalSize;
    DM             dm = NULL;
    const char    *dmtype = NULL;

    /* For DMDA path */
    Vec            nat = NULL, seqNat = NULL;
    VecScatter     scatNat = NULL;
    const PetscScalar *nar = NULL;
    PetscScalar    *buf = NULL;

    /* For generic (no DM) path */
    Vec            seq = NULL;
    VecScatter     scat = NULL;
    const PetscScalar *sar = NULL;

    PetscFunctionBeginUser;

    ierr = PetscObjectGetComm((PetscObject)inVec, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = VecGetSize(inVec, &globalSize); CHKERRQ(ierr);
    *N = globalSize;
    *arrayOut = NULL;

    ierr = VecGetDM(inVec, &dm); CHKERRQ(ierr);
    if (dm) { ierr = DMGetType(dm, &dmtype); CHKERRQ(ierr); }

    if (dmtype && !strcmp(dmtype, DMDA)) {
        /* --- DMDA path: go to NATURAL ordering, then gather to rank 0 --- */
        ierr = DMDACreateNaturalVector(dm, &nat); CHKERRQ(ierr);
        ierr = DMDAGlobalToNaturalBegin(dm, inVec, INSERT_VALUES, nat); CHKERRQ(ierr);
        ierr = DMDAGlobalToNaturalEnd  (dm, inVec, INSERT_VALUES, nat); CHKERRQ(ierr);

        ierr = VecScatterCreateToZero(nat, &scatNat, &seqNat); CHKERRQ(ierr);
        ierr = VecScatterBegin(scatNat, nat, seqNat, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd  (scatNat, nat, seqNat, INSERT_VALUES, SCATTER_FORWARD);   CHKERRQ(ierr);

        if (rank == 0) {
            PetscInt nseq;
            ierr = VecGetLocalSize(seqNat, &nseq); CHKERRQ(ierr);
            ierr = VecGetArrayRead(seqNat, &nar);  CHKERRQ(ierr);

            ierr = PetscMalloc1(nseq, &buf);       CHKERRQ(ierr);
            ierr = PetscMemcpy(buf, nar, (size_t)nseq * sizeof(PetscScalar)); CHKERRQ(ierr);

            ierr = VecRestoreArrayRead(seqNat, &nar); CHKERRQ(ierr);
            *arrayOut = (double*)buf; /* hand back as double* for drop-in compatibility */
        }

        ierr = VecScatterDestroy(&scatNat); CHKERRQ(ierr);
        ierr = VecDestroy(&seqNat);        CHKERRQ(ierr);
        ierr = VecDestroy(&nat);           CHKERRQ(ierr);
    } else {
        /* --- No DM attached: plain gather in Vec’s global (parallel) layout order --- */
        ierr = VecScatterCreateToZero(inVec, &scat, &seq); CHKERRQ(ierr);
        ierr = VecScatterBegin(scat, inVec, seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd  (scat, inVec, seq, INSERT_VALUES, SCATTER_FORWARD);   CHKERRQ(ierr);

        if (rank == 0) {
            PetscInt nseq;
            ierr = VecGetLocalSize(seq, &nseq); CHKERRQ(ierr);
            ierr = VecGetArrayRead(seq, &sar);  CHKERRQ(ierr);

            ierr = PetscMalloc1(nseq, &buf);    CHKERRQ(ierr);
            ierr = PetscMemcpy(buf, sar, (size_t)nseq * sizeof(PetscScalar)); CHKERRQ(ierr);

            ierr = VecRestoreArrayRead(seq, &sar); CHKERRQ(ierr);
            *arrayOut = (double*)buf;
        }

        ierr = VecScatterDestroy(&scat); CHKERRQ(ierr);
        ierr = VecDestroy(&seq);        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Gathers any DMSwarm field from all ranks to a single, contiguous array on rank 0.
 *
 * This is a generic, type-aware version of SwarmFieldToArrayOnRank0.
 * It is a COLLECTIVE operation.
 *
 * @param[in]  swarm             The DMSwarm to gather from.
 * @param[in]  field_name        The name of the field to gather.
 * @param[out] n_total_particles (Output on rank 0) Total number of particles in the global swarm.
 * @param[out] n_components      (Output on rank 0) Number of components for the field.
 * @param[out] gathered_array    (Output on rank 0) A newly allocated array containing the full, gathered data.
 *                               The caller is responsible for freeing this memory and for casting it to the correct type.
 * @return PetscErrorCode
 */
PetscErrorCode SwarmFieldToArrayOnRank0(DM swarm, const char *field_name, PetscInt *n_total_particles, PetscInt *n_components, void **gathered_array)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    PetscInt       nlocal, nglobal, bs;
    void           *local_array_void;
    size_t         element_size = 0;
    MPI_Datatype   mpi_type = MPI_BYTE; // We'll send raw bytes

    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // All ranks get swarm properties to determine send/receive counts
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &nglobal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, field_name, &bs, NULL, &local_array_void); CHKERRQ(ierr);
    
    // Determine the size of one element of the field's data type
    if (strcasecmp(field_name, "DMSwarm_pid") == 0) {
        element_size = sizeof(PetscInt64);
    } else if (strcasecmp(field_name, "DMSwarm_CellID") == 0 || strcasecmp(field_name, "DMSwarm_location_status") == 0) {
        element_size = sizeof(PetscInt);
    } else {
        element_size = sizeof(PetscScalar);
    }

    if (rank == 0) {
        *n_total_particles = nglobal;
        *n_components = bs;
        *gathered_array = NULL; // Initialize output
    }

    if (size == 1) { // Serial case is a simple copy
        if (rank == 0) {
            ierr = PetscMalloc(nglobal * bs * element_size, gathered_array); CHKERRQ(ierr);
            ierr = PetscMemcpy(*gathered_array, local_array_void, nglobal * bs * element_size); CHKERRQ(ierr);
        }
    } else { // Parallel case: use MPI_Gatherv
        PetscInt *recvcounts = NULL, *displs = NULL;
        if (rank == 0) {
            ierr = PetscMalloc1(size, &recvcounts); CHKERRQ(ierr);
            ierr = PetscMalloc1(size, &displs); CHKERRQ(ierr);
        }
        PetscInt sendcount = nlocal * bs;
        
        // Gather the number of elements (not bytes) from each rank
        ierr = MPI_Gather(&sendcount, 1, MPIU_INT, recvcounts, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        if (rank == 0) {
            displs[0] = 0;
            // Convert counts and calculate displacements in terms of BYTES
            for (PetscMPIInt i = 0; i < size; i++) recvcounts[i] *= element_size;
            for (PetscMPIInt i = 1; i < size; i++) displs[i] = displs[i-1] + recvcounts[i-1];
            
            ierr = PetscMalloc(nglobal * bs * element_size, gathered_array); CHKERRQ(ierr);
        }
        
        // Use Gatherv with MPI_BYTE to handle any data type generically
        ierr = MPI_Gatherv(local_array_void, nlocal * bs * element_size, MPI_BYTE,
                           *gathered_array, recvcounts, displs, MPI_BYTE,
                           0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        if (rank == 0) {
            ierr = PetscFree(recvcounts); CHKERRQ(ierr);
            ierr = PetscFree(displs); CHKERRQ(ierr);
        }
    }

    ierr = DMSwarmRestoreField(swarm, field_name, &bs, NULL, &local_array_void); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Displays a structured banner summarizing the simulation configuration.
 *
 * This function prints key simulation parameters to standard output.
 * It is intended to be called ONLY by MPI rank 0.
 * It retrieves global domain bounds from `user->global_domain_bbox` and boundary
 * conditions for all faces from `user->face_bc_types`.
 *
 * @param[in] user         Pointer to `UserCtx` structure.
 * @param[in] StartTime    Initial simulation time.
 * @param[in] StartStep    Starting timestep index.
 * @param[in] StepsToRun   Total number of timesteps to run.
 * @param[in] num_mpi_procs Total number of MPI processes.
 * @param[in] total_num_particles Total number of particles.
 * @param[in] bboxlist     (If rank 0 needed to compute global_domain_bbox here, otherwise NULL)
 *
 * @return PetscErrorCode  Returns `0` on success.
 */
PetscErrorCode DisplayBanner(SimCtx *simCtx) // bboxlist is only valid on rank 0
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Cmpnts         global_min_coords, global_max_coords;
    PetscReal      StartTime;
    PetscInt       StartStep,StepsToRun,total_num_particles;
    PetscMPIInt    num_mpi_procs;

    // SimCtx *simCtx = user->simCtx;
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    num_mpi_procs = simCtx->size;
    StartTime = simCtx->StartTime;
    StartStep = simCtx->StartStep;
    StepsToRun = simCtx->StepsToRun;
    total_num_particles = simCtx->np;
    BoundingBox *bboxlist_on_rank0 = simCtx->bboxlist;
    

    PetscFunctionBeginUser;

    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DisplayBanner - UserCtx pointer is NULL.");
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        // If global_domain_bbox is not pre-populated in UserCtx, compute it here from bboxlist_on_rank0
        // This assumes bboxlist_on_rank0 is valid and contains all local bounding boxes on rank 0.
        if (bboxlist_on_rank0 && num_mpi_procs > 0) {
            global_min_coords = bboxlist_on_rank0[0].min_coords;
            global_max_coords = bboxlist_on_rank0[0].max_coords;
            for (PetscMPIInt p = 1; p < num_mpi_procs; ++p) {
                global_min_coords.x = PetscMin(global_min_coords.x, bboxlist_on_rank0[p].min_coords.x);
                global_min_coords.y = PetscMin(global_min_coords.y, bboxlist_on_rank0[p].min_coords.y);
                global_min_coords.z = PetscMin(global_min_coords.z, bboxlist_on_rank0[p].min_coords.z);
                global_max_coords.x = PetscMax(global_max_coords.x, bboxlist_on_rank0[p].max_coords.x);
                global_max_coords.y = PetscMax(global_max_coords.y, bboxlist_on_rank0[p].max_coords.y);
                global_max_coords.z = PetscMax(global_max_coords.z, bboxlist_on_rank0[p].max_coords.z);
            }
            // Optionally store this in user->global_domain_bbox if it's useful elsewhere
            // user->global_domain_bbox.min_coords = global_min_coords;
            // user->global_domain_bbox.max_coords = global_max_coords;
        } else {
            // Fallback or warning if bboxlist is not available for global calculation
            LOG_ALLOW(PETSC_COMM_SELF, LOG_WARNING, "(Rank 0) - bboxlist not provided or num_mpi_procs <=0; using user->bbox for domain bounds.\n");
	    // global_min_coords = user->bbox.min_coords; // Use local bbox of rank 0 as fallback
	    // global_max_coords = user->bbox.max_coords;
        }


        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "                          CASE SUMMARY                       \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Grid Points     : %d X %d X %d\n", user->IM, user->JM, user->KM); CHKERRQ(ierr);
	    ierr = PetscPrintf(PETSC_COMM_SELF, " Cells           : %d X %d X %d\n", user->IM - 1, user->JM - 1, user->KM - 1); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (X)    : %.6f to %.6f\n", (double)global_min_coords.x, (double)global_max_coords.x); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (Y)    : %.6f to %.6f\n", (double)global_min_coords.y, (double)global_max_coords.y); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (Z)    : %.6f to %.6f\n", (double)global_min_coords.z, (double)global_max_coords.z); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------- Boundary Conditions --------------------\n"); CHKERRQ(ierr);
	    const int face_name_width = 17; // Adjusted for longer names (Zeta,Eta,Xi)
        const char* field_init_str = FieldInitializationToString(simCtx->FieldInitialization);
        const char* particle_init_str = ParticleInitializationToString(simCtx->ParticleInitialization);
        for (PetscInt i_face = 0; i_face < 6; ++i_face) {
	    BCFace current_face = (BCFace)i_face;
	    // The BCFaceToString will now return the Xi, Eta, Zeta versions
	    const char* face_str = BCFaceToString(current_face); 
	    const char* bc_type_str = BCTypeToString(user->boundary_faces[current_face].mathematical_type);    
	    ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s\n", 
			     face_name_width, face_str, bc_type_str); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------------------------------------------------\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Start Time                  : %.4f\n", (double)StartTime); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Timestep Size               : %.4f\n", (double)simCtx->dt); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Starting Step               : %d\n", StartStep); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Total Steps to Run          : %d\n", StepsToRun); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Number of MPI Processes     : %d\n", num_mpi_procs); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Number of Particles         : %d\n", total_num_particles); CHKERRQ(ierr);
	    ierr = PetscPrintf(PETSC_COMM_WORLD," Reynolds Number             : %le\n", simCtx->ren); CHKERRQ(ierr);
	    ierr = PetscPrintf(PETSC_COMM_WORLD," Stanton Number              : %le\n", simCtx->st); CHKERRQ(ierr);
	    ierr = PetscPrintf(PETSC_COMM_WORLD," CFL Number                  : %le\n", simCtx->cfl); CHKERRQ(ierr);
	    ierr = PetscPrintf(PETSC_COMM_WORLD," Von-Neumann Number          : %le\n", simCtx->vnn); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Initialization Mode: %s\n", particle_init_str); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Large Eddy Simulation Model : %s\n", LESFlagToString(simCtx->les)); CHKERRQ(ierr);
        if (simCtx->ParticleInitialization == 0 || simCtx->ParticleInitialization == 3) {
            if (user->inletFaceDefined) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : %s (Enum Val: %d)\n", BCFaceToString(user->identifiedInletBCFace), user->identifiedInletBCFace); CHKERRQ(ierr);
            } else {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : --- (No INLET face identified)\n"); CHKERRQ(ierr);
            }
        }

        ierr = PetscPrintf(PETSC_COMM_SELF, " Field Initialization Mode   : %s\n", field_init_str); CHKERRQ(ierr);
        if (simCtx->FieldInitialization == 1) {
	  ierr = PetscPrintf(PETSC_COMM_SELF, " Constant Velocity           : x - %.4f, y - %.4f, z - %.4f \n", (double)simCtx->InitialConstantContra.x,(double)simCtx->InitialConstantContra.y,(double)simCtx->InitialConstantContra.z ); CHKERRQ(ierr);
        }
	
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ParsePostProcessingSettings"
/**
 * @brief Initializes post-processing settings from a config file and command-line overrides.
 *
 * This function establishes the configuration for a post-processing run by:
 * 1. Setting hardcoded default values in the PostProcessParams struct.
 * 2. Reading a configuration file to override the defaults.
 * 3. Parsing command-line options (-startTime, -endTime, etc.) which can override
 *    both the defaults and the file settings.
 *
 * @param configFile The path to the configuration file to parse.
 * @param pps Pointer to the PostProcessParams struct to be populated.
 * @return PetscErrorCode
 */
PetscErrorCode  ParsePostProcessingSettings(SimCtx *simCtx)
{
    FILE *file;
    char line[1024];
    PetscBool startTimeSet, endTimeSet, timeStepSet;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!simCtx || !simCtx->pps) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "SimCtx or its pps member is NULL in ParsePostProcessingSettings.");
    }

    char *configFile = simCtx->PostprocessingControlFile;
    PostProcessParams *pps = simCtx->pps;


    // --- 1. Set Sane Defaults First ---
    pps->startTime = 0;
    pps->endTime = 0;
    pps->timeStep = 1;
    pps->outputParticles = PETSC_FALSE;
    pps->particle_output_freq = simCtx->LoggingFrequency; // Default to logging frequency;
    strcpy(pps->process_pipeline, "");
    strcpy(pps->output_fields_instantaneous, "Ucat,P");
    strcpy(pps->output_fields_averaged, "");
    strcpy(pps->output_prefix, "Field");
    strcpy(pps->particle_output_prefix,"Particle");
    strcpy(pps->particle_fields,"velocity,CellID,weight,pid");
    strcpy(pps->particle_pipeline,"");
    strcpy(pps->particleExt,"dat"); // The input file format for particles.
    strcpy(pps->eulerianExt,"dat"); // The input file format for Eulerian fields.
    pps->reference[0] = pps->reference[1] = pps->reference[2] = 1;
    strncpy(pps->source_dir, simCtx->output_dir, sizeof(pps->source_dir) - 1);
    pps->source_dir[sizeof(pps->source_dir) - 1] = '\0'; // Ensure null-termination

    // --- 2. Parse the Configuration File (overrides defaults) ---
    file = fopen(configFile, "r");
    if (file) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Parsing post-processing config file: %s\n", configFile);
        while (fgets(line, sizeof(line), file)) {
            char *key, *value, *comment;
            comment = strchr(line, '#'); if (comment) *comment = '\0';
            TrimWhitespace(line); if (strlen(line) == 0) continue;
            key = strtok(line, "="); value = strtok(NULL, "=");
            if (key && value) {
                TrimWhitespace(key); TrimWhitespace(value);
                if (strcmp(key, "startTime") == 0) pps->startTime = atoi(value);
                else if (strcmp(key, "endTime") == 0) pps->endTime = atoi(value);
                else if (strcmp(key, "timeStep") == 0) pps->timeStep = atoi(value);
                else if (strcmp(key, "output_particles") == 0) {
                    if (strcasecmp(value, "true") == 0) pps->outputParticles = PETSC_TRUE;
                }
                else if (strcasecmp(key, "process_pipeline") == 0) {
                    strncpy(pps->process_pipeline, value, MAX_PIPELINE_LENGTH - 1);
                    pps->process_pipeline[MAX_PIPELINE_LENGTH - 1] = '\0'; // Ensure null-termination
                } else if (strcasecmp(key, "output_fields_instantaneous") == 0) {
                    strncpy(pps->output_fields_instantaneous, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->output_fields_instantaneous[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "output_prefix") == 0) {
                    strncpy(pps->output_prefix, value, MAX_FILENAME_LENGTH - 1);
                    pps->output_prefix[MAX_FILENAME_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_output_prefix") == 0) {
                    strncpy(pps->particle_output_prefix, value, MAX_FILENAME_LENGTH - 1);
                    pps->particle_output_prefix[MAX_FILENAME_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_fields_instantaneous") == 0) {
                    strncpy(pps->particle_fields, value, MAX_FIELD_LIST_LENGTH - 1);
                    pps->particle_fields[MAX_FIELD_LIST_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_pipeline") == 0) {
                    strncpy(pps->particle_pipeline, value, MAX_PIPELINE_LENGTH - 1);
                    pps->particle_pipeline[MAX_PIPELINE_LENGTH - 1] = '\0';
                } else if (strcasecmp(key, "particle_output_freq") == 0) {
                    pps->particle_output_freq = atoi(value);
                } else if (strcmp(key, "reference_ip") == 0) {pps->reference[0] = atoi(value);
                } else if (strcmp(key, "reference_jp") == 0) {pps->reference[1] = atoi(value);
                } else if (strcmp(key, "reference_kp") == 0) {pps->reference[2] = atoi(value);
                } else if (strcasecmp(key, "source_directory") == 0) {
                    strncpy(pps->source_dir, value, sizeof(pps->source_dir) - 1);
                    pps->source_dir[sizeof(pps->source_dir) - 1] = '\0';
                } else {
                    LOG_ALLOW(GLOBAL, LOG_WARNING, "Unknown key '%s' in post-processing config file. Ignoring.\n", key);
                }
                // Add parsing for pipeline, fields, etc. in later phases
            }
        }
        fclose(file);
    } else {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Could not open post-processing config file '%s'. Using defaults and command-line overrides.\n", configFile);
    }

    // --- 3. Parse Command-Line Options (overrides file settings and defaults) ---
    PetscOptionsGetInt(NULL, NULL, "-startTime", &pps->startTime, &startTimeSet);
    PetscOptionsGetInt(NULL, NULL, "-endTime", &pps->endTime, &endTimeSet);
    PetscOptionsGetInt(NULL, NULL, "-timeStep", &pps->timeStep, &timeStepSet);
    PetscOptionsGetBool(NULL, NULL, "-output_particles", &pps->outputParticles, NULL);

    if(pps->endTime==-1){
        pps->endTime = simCtx->StartStep + simCtx->StepsToRun; // Total steps if endTime is set to -1.
    }

    // If only startTime is given on command line, run for a single step
    if (startTimeSet && !endTimeSet) {
        pps->endTime = pps->startTime;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Post-processing configured to run from t=%d to t=%d with step %d. Particle output: %s.\n",
              pps->startTime, pps->endTime, pps->timeStep, pps->outputParticles ? "TRUE" : "FALSE");

    LOG_ALLOW(GLOBAL, LOG_INFO, "Process Pipeline: %s\n", pps->process_pipeline);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Instantaneous Output Fields: %s\n", pps->output_fields_instantaneous);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Output Prefix: %s\n", pps->output_prefix);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Output Prefix: %s\n", pps->particle_output_prefix);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Fields: %s\n", pps->particle_fields);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Pipeline: %s\n", pps->particle_pipeline);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Output Frequency: %d\n", pps->particle_output_freq);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ParseScalingInformation"
/**
 * @brief Parses physical scaling parameters from command-line options.
 *
 * This function reads the reference length, velocity, and density from the
 * PETSc options database (provided via -scaling_L_ref, etc.). It populates
 * the simCtx->scaling struct and calculates the derived reference pressure.
 * It sets default values of 1.0 for a fully non-dimensional case if the
 * options are not provided.
 *
 * @param[in,out] simCtx The simulation context whose 'scaling' member will be populated.
 * @return PetscErrorCode
 */
PetscErrorCode ParseScalingInformation(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscBool      flg;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!simCtx) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx is NULL in ParseScalingInformation");

    // --- 1. Set default values to 1.0 ---
    // This represents a purely non-dimensional run if no scaling is provided.
    simCtx->scaling.L_ref   = 1.0;
    simCtx->scaling.U_ref   = 1.0;
    simCtx->scaling.rho_ref = 1.0;
    
    // --- 2. Read overrides from the command line / control file ---
    ierr = PetscOptionsGetReal(NULL, NULL, "-scaling_L_ref", &simCtx->scaling.L_ref, &flg); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-scaling_U_ref", &simCtx->scaling.U_ref, &flg); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-scaling_rho_ref", &simCtx->scaling.rho_ref, &flg); CHKERRQ(ierr);
    
    // --- 3. Calculate derived scaling factors ---
    // Check for division by zero to be safe, though U_ref should be positive.
    if (simCtx->scaling.U_ref <= 0.0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Reference velocity U_ref must be positive. Got %g", (double)simCtx->scaling.U_ref);
    }
    simCtx->scaling.P_ref = simCtx->scaling.rho_ref * simCtx->scaling.U_ref * simCtx->scaling.U_ref;

    // --- 4. Log the final, effective scales for verification ---
    LOG(GLOBAL, LOG_INFO, "---------------- Physical Scales Initialized -----------------\n");
    LOG(GLOBAL, LOG_INFO, "  L_ref: %.4f, U_ref: %.4f, rho_ref: %.4f, P_ref: %.4f\n",
        simCtx->scaling.L_ref, simCtx->scaling.U_ref, simCtx->scaling.rho_ref, simCtx->scaling.P_ref);
    LOG(GLOBAL, LOG_INFO, "--------------------------------------------------------------\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


//  ---------------------------------------------------------------------
//  UTILITY FUNCTIONS NOT USED IN THE MAIN CODE ( NOT  SUPPORTED  ANYMORE,INCLUDED FOR BACKWARD COMPATIBILITY)
//  ---------------------------------------------------------------------

//=====================================================================
//   1) ReadDataFileToArray

//   See the function-level comments in io.h for a summary.
//   This reads one value per line from an ASCII file. Rank 0 does I/O,
//   broadcasts the data to all ranks.
//===================================================================== */
PetscInt ReadDataFileToArray(const char   *filename,
                        double      **data_out,
                        PetscInt          *Nout,
                        MPI_Comm      comm)
{
    /* STEP 0: Prepare local variables & log function entry */
    PetscMPIInt    rank, size;
    PetscErrorCode ierr;
    FILE  *fp = NULL;
    PetscInt    N   = 0;            /* number of lines/values read on rank 0 */
    double *array = NULL;      /* pointer to local array on each rank */
    PetscInt    fileExistsFlag = 0; /* 0 = doesn't exist, 1 = does exist */

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Start reading from file: %s\n",
              filename);

    /* Basic error checking: data_out, Nout must be non-null. */
    if (!filename || !data_out || !Nout) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Null pointer argument provided.\n");
        return 1;
    }

    /* Determine rank/size for coordinating I/O. */
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* STEP 1: On rank 0, check if file can be opened. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (fp) {
            fileExistsFlag = 1;
            fclose(fp);
        }
    }

    /* STEP 2: Broadcast file existence to all ranks. */
    // In ReadDataFileToArray:
    ierr = MPI_Bcast(&fileExistsFlag, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    if (!fileExistsFlag) {
        /* If file does not exist, log & return. */
        if (!rank) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "File '%s' not found.\n",
                      filename);
        }
        return 2;
    }

    /* STEP 3: Rank 0 re-opens and reads the file, counting lines, etc. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "File '%s' could not be opened for reading.\n",
                      filename);
            return 3;
        }

        /* (3a) Count lines first. */
        {
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                N++;
            }
        }

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "File '%s' has %d lines.\n",
                  filename, N);

        /* (3b) Allocate array on rank 0. */
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            fclose(fp);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "malloc failed for array.\n");
            return 4;
        }

        /* (3c) Rewind & read values into array. */
        rewind(fp);
        {
            PetscInt i = 0;
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                double val;
                if (sscanf(line, "%lf", &val) == 1) {
                    array[i++] = val;
                }
            }
        }
        fclose(fp);

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Successfully read %d values from '%s'.\n",
                  N, filename);
    }

    /* STEP 4: Broadcast the integer N to all ranks. */
    ierr = MPI_Bcast(&N, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    /* STEP 5: Each rank allocates an array to receive the broadcast if rank>0. */
    if (rank) {
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "malloc failed on rank %d.\n",
                      rank);
            return 5;
        }
    }

    /* STEP 6: Broadcast the actual data from rank 0 to all. */
    ierr = MPI_Bcast(array, N, MPI_DOUBLE, 0, comm); CHKERRQ(ierr);

    /* STEP 7: Assign outputs on all ranks. */
    *data_out = array;
    *Nout     = N;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Done. Provided array of length=%d to all ranks.\n",
              N);
    return 0; /* success */
}

/**
 * @brief Reads coordinate data (for particles)  from file into a PETSc Vec, then gathers it to rank 0.
 *
 * This function uses \c ReadFieldData to fill a PETSc Vec with coordinate data,
 * then leverages \c VecToArrayOnRank0 to gather that data into a contiguous array
 * (valid on rank 0 only).
 *
 * @param[in]  timeIndex    The time index used to construct file names.
 * @param[in]  user         Pointer to the user context.
 * @param[out] coordsArray  On rank 0, will point to a newly allocated array holding the coordinates.
 * @param[out] Ncoords      On rank 0, the length of \p coordsArray. On other ranks, 0.
 *
 * @return PetscErrorCode  Returns 0 on success, or non-zero on failures.
 */
PetscErrorCode ReadPositionsFromFile(PetscInt timeIndex,
                                      UserCtx *user,
                                      double **coordsArray,
                                      PetscInt *Ncoords)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;
  Vec            coordsVec;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Creating coords Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

  // For example: "position" is the name of the coordinate data
  ierr = ReadFieldData(user, "position", coordsVec, timeIndex, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "Error reading position data (ti=%d).\n",
              timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadPositions - Gathering coords Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(coordsVec, Ncoords, coordsArray);CHKERRQ(ierr);

  ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Successfully gathered coordinates. Ncoords=%d.\n", *Ncoords);
  PetscFunctionReturn(0);
}


/**
 * @brief Reads a named field from file into a PETSc Vec, then gathers it to rank 0.
 *
 * This function wraps \c ReadFieldData and \c VecToArrayOnRank0 into a single step.
 * The gathered data is stored in \p scalarArray on rank 0, with its length in \p Nscalars.
 *
 * @param[in]  timeIndex     The time index used to construct file names.
 * @param[in]  fieldName     Name of the field to be read (e.g., "velocity").
 * @param[in]  user          Pointer to the user context.
 * @param[out] scalarArray   On rank 0, a newly allocated array holding the field data.
 * @param[out] Nscalars      On rank 0, length of \p scalarArray. On other ranks, 0.
 *
 * @return PetscErrorCode  Returns 0 on success, or non-zero on failures.
 */
PetscErrorCode ReadFieldDataToRank0(PetscInt timeIndex,
                                           const char *fieldName,
                                           UserCtx *user,
                                           double **scalarArray,
                                           PetscInt *Nscalars)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;
  Vec            fieldVec;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Creating field Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

  ierr = ReadFieldData(user, fieldName, fieldVec, timeIndex, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "Error reading field '%s' (ti=%d).\n",
              fieldName, timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Gathering field Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(fieldVec, Nscalars, scalarArray);CHKERRQ(ierr);

  ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "Successfully gathered field '%s'. Nscalars=%d.\n",
            fieldName, *Nscalars);
  PetscFunctionReturn(0);
}


