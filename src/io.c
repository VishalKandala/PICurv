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

            if (strcmp(firstTok, "FDFGRID") == 0) {
                // Header is present – read nblk from the next line
                if (fscanf(fd, "%d", &g_nblk_from_file) != 1)
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Expected number of blocks after \"FDFGRID\" in %s", simCtx->grid_file);
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
    else if (strcasecmp(str, "NOGRAD")    == 0) *type_out = NOGRAD;
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
    else if (strcasecmp(str, "allcopy")             == 0) *handler_out = BC_HANDLER_NOGRAD_COPY_GHOST;
    else if (strcasecmp(str, "parabolic")           == 0) *handler_out = BC_HANDLER_INLET_PARABOLIC;
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
        case NOGRAD:
	    if(handler != BC_HANDLER_NOGRAD_COPY_GHOST) return PETSC_ERR_ARG_WRONG;
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

//================================================================================
//
//                        PUBLIC PARSING FUNCTION
//
//================================================================================

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
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "    - Found param: [%s] = [%s] \n", new_param->key, new_param->value);
            }
        }
        fclose(file);
    }

    // =========================================================================
    //               BROADCASTING THE CONFIGURATION FROM RANK 0
    // =========================================================================
    // This is a critical step to ensure all processes have the same configuration.
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d broadcasting/receiving BC configuration.\n", rank);

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
        if (user->boundary_faces[i].mathematical_type == INLET) {
            user->inletFaceDefined = PETSC_TRUE;
            user->identifiedInletBCFace = (BCFace)i;
            LOG_ALLOW(GLOBAL, LOG_INFO, "Inlet face for particle initialization identified as Face %d.\n", i);
            break; // Found the first one, stop looking
        }
    }

    if (rank == 0) {
        // Rank 0 can now free the linked lists it created, as they have been broadcast.
        // Or, if user->boundary_faces was used directly, this is not needed.
    }
    
    PetscFunctionReturn(0);
}


/**
 * @brief Checks for a data file's existence in a parallel-safe manner.
 *
 * Only Rank 0 checks for the file on disk. The result (true or false) is
 * then broadcast to all other processes in the communicator. This ensures all
 * processes make a collective, synchronized decision.
 *
 * @param ti The time index of the file.
 * @param fieldName The name of the field.
 * @param ext The file extension.
 * @param fileExists [out] The result, which will be identical on all ranks.
 * @return PetscErrorCode
 */
static PetscErrorCode CheckDataFile(PetscInt ti, const char *fieldName, const char *ext, PetscBool *fileExists)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    MPI_Comm       comm = PETSC_COMM_WORLD;
    PetscInt       placeholder_int = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        char filename[PETSC_MAX_PATH_LEN];
        // Use the same standardized, rank-independent filename format

        ierr =  PetscSNPrintf(filename, sizeof(filename), "results/%s%05"PetscInt_FMT"_%d.%s", fieldName, ti, placeholder_int, ext); 
        ierr = PetscTestFile(filename, 'r', fileExists); CHKERRQ(ierr);
        if (!(*fileExists)) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "CheckDataFile (Rank 0) - Optional data file '%s' not found.\n", filename);
        }
    }

    // Broadcast the result from Rank 0 to all other ranks.
    // We cast the PetscBool to a PetscMPIInt for MPI_Bcast.
    PetscMPIInt fileExists_int = (rank == 0) ? (PetscMPIInt)(*fileExists) : 0;
    ierr = MPI_Bcast(&fileExists_int, 1, MPI_INT, 0, comm); CHKERRMPI(ierr);
    *fileExists = (PetscBool)fileExists_int;

    PetscFunctionReturn(0);
}

/************************************************************************************************
*  @file   io.c
*  @brief  Utility for (re-)loading PETSc Vec-based fields written by rank-0 in a previous run
*          – works in both serial and MPI executions.
*
*  The restart files are always written as **one** sequential Vec per field
*  ( `<results>/<field_name><step>_0.<ext>` ).  In parallel we therefore
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

   PetscFunctionBeginUser;

   ierr = PetscObjectGetComm((PetscObject)field_vec,&comm);CHKERRQ(ierr);
   ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
   ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

   /* ---------------------------------------------------------------------
    * Compose <results>/<field_name><step with 5 digits>_0.<ext>
    * (all restart files are written by rank-0 with that naming scheme).
    * ------------------------------------------------------------------ */
   ierr = PetscSNPrintf(filename,sizeof(filename),
                        "results/%s%05" PetscInt_FMT "_0.%s",
                        field_name,ti,ext);CHKERRQ(ierr);

   LOG_ALLOW(GLOBAL,LOG_DEBUG,
             "ReadFieldData - Attempting to read <%s> on rank %d/%d\n",
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
                         "Restart file not found: %s",filename);

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
                "ReadFieldData - Loaded <%s> (serial path)\n",filename);
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

      LOG_ALLOW(GLOBAL,LOG_DEBUG,
                "ReadFieldData - Rank 0 successfully loaded <%s>\n",filename);
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
             "ReadFieldData - Loaded <%s> (parallel path)\n",filename);

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadSimulationFields - Starting to read simulation fields.\n");

    // Read Cartesian velocity field
    ierr = ReadFieldData(user, "ufield", user->Ucat, ti, "dat"); CHKERRQ(ierr);

    // Read contravariant velocity field
    ierr = ReadFieldData(user, "vfield", user->Ucont, ti, "dat"); CHKERRQ(ierr);

    // Read pressure field
    ierr = ReadFieldData(user, "pfield", user->P, ti, "dat"); CHKERRQ(ierr);

    // Read node state field (nvert)
    ierr = ReadFieldData(user, "nvfield", user->Nvert_o, ti, "dat"); CHKERRQ(ierr);

    // Read Particle Count field
    ierr = ReadFieldData(user, "ParticleCount", user->ParticleCount, ti, "dat"); CHKERRQ(ierr);

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadSimulationFields - Finished reading simulation fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadStatisticalFields - Starting to read statistical fields.\n");

    ierr = ReadFieldData(user, "su0", user->Ucat_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "su1", user->Ucat_cross_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "su2", user->Ucat_square_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "sp", user->P_sum, ti, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadStatisticalFields - Finished reading statistical fields.\n");

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
    Vec Cs;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadLESFields - Starting to read LES fields.\n");

    VecDuplicate(user->P, &Cs);
    ierr = ReadFieldData(user, "cs", Cs, ti, "dat"); CHKERRQ(ierr);
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    VecDestroy(&Cs);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadLESFields - Finished reading LES fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadRANSFields - Starting to read RANS fields.\n");

    ierr = ReadFieldData(user, "kfield", user->K_Omega, ti, "dat"); CHKERRQ(ierr);
    VecCopy(user->K_Omega, user->K_Omega_o);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadRANSFields - Finished reading RANS fields.\n");

    return 0;
}

/**
 * @brief Writes data from a specific PETSc vector to a file.
 *
 * This function uses the field name to construct the file path and writes the data
 * from the provided PETSc vector to the corresponding file.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context.
 * @param[in] field_name Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[in] field_vec  PETSc vector containing the field data to write.
 * @param[in] ti         Time index for constructing the file name.
 * @param[in] ext        File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
/*
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    char filen[128];

    // Construct the file name
    snprintf(filen, sizeof(filen), "results/%s%05d_%d.%s", field_name, ti, user->_this, ext);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteFieldData - Attempting to write file: %s\n", filen);

    // Open the file for writing
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

    PetscReal vmin, vmax;
    ierr = VecMin(field_vec, NULL, &vmin); CHKERRQ(ierr);
    ierr = VecMax(field_vec, NULL, &vmax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"%s step %d  min=%.6e  max=%.6e\n",field_name, ti, (double)vmin, (double)vmax);

    
    // Write data from the vector
    ierr = VecView(field_vec, viewer); CHKERRQ(ierr);

    // Close the file viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);

    return 0;
}
*/

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

    PetscFunctionBeginUser;

    /* ------------------------------------------------------------ */
    /*                  Basic communicator information              */
    /* ------------------------------------------------------------ */
    ierr = PetscObjectGetComm((PetscObject)field_vec,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

    ierr = PetscSNPrintf(filename,sizeof(filename),
                         "results/%s%05" PetscInt_FMT "_%d.%s",
                         field_name,ti,placeholder_int,ext);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_DEBUG,
              "WriteFieldData - Preparing to write <%s> on rank %d/%d\n",
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
                  "WriteFieldData - Wrote <%s> (serial path)\n",filename);
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
                  "WriteFieldData - <%s> range = [%.4e … %.4e]\n",
                  field_name,(double)vmin,(double)vmax);

        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,
                                     FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = VecView(seq_vec,viewer);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_INFO,
                  "WriteFieldData - Wrote <%s> (parallel path)\n",filename);
    }

    /* 2.4  Cleanup                                                 */
    ierr = VecScatterDestroy(&scatter);CHKERRQ(ierr);
    ierr = VecDestroy(&seq_vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/*
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    char           filename[PETSC_MAX_PATH_LEN];
    PetscMPIInt    rank, size;
    MPI_Comm       comm;
    PetscInt       placeholder_int = 0;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)field_vec, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

    // Construct a filename that does NOT depend on the rank number.
    // This ensures a single, predictable output file.
    ierr = PetscSNPrintf(filename, sizeof(filename), "results/%s%05"PetscInt_FMT"_%d.%s", field_name, ti, placeholder_int, ext);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteFieldData - Preparing to write file: %s\n", filename);

    // Optional: Log min/max values. This is a collective operation.
    PetscReal vmin, vmax;
    ierr = VecMin(field_vec, NULL, &vmin); CHKERRQ(ierr);
    ierr = VecMax(field_vec, NULL, &vmax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"%s step %d  min=%.6e  max=%.6e\n", field_name, ti, (double)vmin, (double)vmax);

    if (size == 1) {
        // --- SERIAL CASE: Simple and direct write ---
        PetscViewer viewer;
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
        ierr = VecView(field_vec, viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    } else {
        // --- PARALLEL CASE: Gather on Rank 0 and write sequentially ---
        Vec         seq_vec = NULL; // A sequential vector on Rank 0
        VecScatter  scatter_ctx;

        // Create a scatter context to gather data from the parallel vector
        // onto a new sequential vector that will exist only on Rank 0.
        ierr = VecScatterCreateToZero(field_vec, &scatter_ctx, &seq_vec); CHKERRQ(ierr);

        // Perform the scatter (gather) operation.
        ierr = VecScatterBegin(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&scatter_ctx); CHKERRQ(ierr);

        // Now, only Rank 0 has the populated sequential vector and can write it.
        if (rank == 0) {
            PetscViewer viewer;
            ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
            ierr = VecView(seq_vec, viewer); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        }

        // All ranks must destroy the sequential vector, though it only "exists" on Rank 0.
        ierr = VecDestroy(&seq_vec); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);
    PetscFunctionReturn(0);
}
*/
/*
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    MPI_Comm       comm;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)field_vec, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

    if (size == 1) {
        // SERIAL CASE: This is simple and correct.
        PetscViewer viewer;
        char filename[PETSC_MAX_PATH_LEN];
        ierr = PetscSNPrintf(filename, sizeof(filename), "output/%s_step%d.%s", field_name, ti, ext); CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
        ierr = VecView(field_vec, viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    } else {
        // PARALLEL CASE: The robust "Gather on 0 and Write" pattern.
        Vec         seq_vec = NULL;
        VecScatter  scatter_ctx;
        IS          is_from, is_to;
        PetscInt    global_size;

        ierr = VecGetSize(field_vec, &global_size); CHKERRQ(ierr);
        ierr = ISCreateStride(comm, global_size, 0, 1, &is_from); CHKERRQ(ierr);
        ierr = ISCreateStride(PETSC_COMM_SELF, global_size, 0, 1, &is_to); CHKERRQ(ierr);

        // Rank 0 creates the destination sequential vector.
        if (rank == 0) {
            ierr = VecCreate(PETSC_COMM_SELF, &seq_vec); CHKERRQ(ierr);
            ierr = VecSetSizes(seq_vec, global_size, global_size); CHKERRQ(ierr);
            ierr = VecSetFromOptions(seq_vec); CHKERRQ(ierr);
        }

        // Create the general scatter context.
        ierr = VecScatterCreate(field_vec, is_from, seq_vec, is_to, &scatter_ctx); CHKERRQ(ierr);
        ierr = VecScatterBegin(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterEnd(scatter_ctx, field_vec, seq_vec, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

        // Only Rank 0 writes the populated sequential vector.
        if (rank == 0) {
            PetscViewer viewer;
            char filename[PETSC_MAX_PATH_LEN];
            ierr = PetscSNPrintf(filename, sizeof(filename), "output/%s_step%d.%s", field_name, ti, ext); CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
            ierr = VecView(seq_vec, viewer); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
        }

        // Clean up.
        ierr = ISDestroy(&is_from); CHKERRQ(ierr);
        ierr = ISDestroy(&is_to); CHKERRQ(ierr);
        ierr = VecScatterDestroy(&scatter_ctx); CHKERRQ(ierr);
        if (rank == 0) {
            ierr = VecDestroy(&seq_vec); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);
    PetscFunctionReturn(0);
}
*/

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
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteSimulationFields - Starting to write simulation fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteSimulationFields - Finished writing simulation fields.\n");

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
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteStatisticalFields - Starting to write statistical fields.\n");

    ierr = WriteFieldData(user, "su0", user->Ucat_sum, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su1", user->Ucat_cross_sum, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su2", user->Ucat_square_sum, simCtx->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "sp", user->P_sum, simCtx->step, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteStatisticalFields - Finished writing statistical fields.\n");

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
    Vec Cs;

    SimCtx *simCtx = user->simCtx;

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteLESFields - Starting to write LES fields.\n");

    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    ierr = WriteFieldData(user, "cs", Cs, simCtx->step, "dat"); CHKERRQ(ierr);
    VecDestroy(&Cs);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteLESFields - Finished writing LES fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteRANSFields - Starting to write RANS fields.\n");

    ierr = WriteFieldData(user, "kfield", user->K_Omega, simCtx->step, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteRANSFields - Finished writing RANS fields.\n");

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
            "WriteSwarmField - Attempting to create global vector from field: %s\n",
            field_name);
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /*
   * 3) Use your existing WriteFieldData() to write the global vector to a file.
   *    The field name, time index, and extension are passed along for naming.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "WriteSwarmField - Calling WriteFieldData for field: %s\n",
            field_name);
  ierr = WriteFieldData(user, field_name, fieldVec, ti, ext);CHKERRQ(ierr);

  /*
   * 4) Destroy the global vector once the data is successfully written.
   *    This step is crucial for avoiding memory leaks. 
   *    DMSwarmDestroyGlobalVectorFromField() is also available in PETSc 3.14.4.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "WriteSwarmField - Destroying the global vector for field: %s\n",
            field_name);
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /* Log and return success. */
  LOG_ALLOW(GLOBAL, LOG_INFO,
            "WriteSwarmField - Successfully wrote swarm data for field: %s\n",
            field_name);

  PetscFunctionReturn(0); /* PETSc macro indicating end of function */
}

/* 
   -------------- Prototypes of Internal-Use-Only (static) ---------------
   These do not appear in io.h since they are "private" to io.c:
     - WriteVTKAppendedBlock
     - WriteVTSXMLHeader, WriteVTPXMLHeader
     - WriteVTSXMLFooter, WriteVTPXMLFooter
     - WriteVTKFileHeader, WriteVTKFileFooter

static PetscInt WriteVTKAppendedBlock(FILE *fp,
                                 const void *buf,
                                 PetscInt nvals,
                                 PetscInt dsize);

static PetscInt WriteVTSXMLHeader(FILE       *fp,
                             PetscInt         mx,
                             PetscInt         my,
                             PetscInt         mz,
                             PetscInt         nnodes,
                             const char *fieldName,
                             PetscInt         numComponents,
                             PetscInt         boffset,
                             PetscInt        *boffsetOut);

static PetscInt WriteVTSXMLFooter(FILE *fp);

static PetscInt WriteVTPXMLHeader(FILE       *fp,
                             PetscInt         npoints,
                             const char *fieldName,
                             PetscInt         boffset,
                             PetscInt        *boffsetOut);

static PetscInt WriteVTPXMLFooter(FILE *fp);

static PetscInt WriteVTKFileHeader(FILE *fp,
                              const VTKMetaData *meta,
                              PetscInt boffset,
                              PetscInt *boffsetOut);

static PetscInt WriteVTKFileFooter(FILE *fp,
                              const VTKMetaData *meta);


=====================================================================
   1) ReadDataFileToArray

   See the function-level comments in io.h for a summary.
   This reads one value per line from an ASCII file. Rank 0 does I/O,
   broadcasts the data to all ranks.
===================================================================== */
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
              "ReadDataFileToArray - Start reading from file: %s\n",
              filename);

    /* Basic error checking: data_out, Nout must be non-null. */
    if (!filename || !data_out || !Nout) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ReadDataFileToArray - Null pointer argument provided.\n");
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
                      "ReadDataFileToArray - File '%s' not found.\n",
                      filename);
        }
        return 2;
    }

    /* STEP 3: Rank 0 re-opens and reads the file, counting lines, etc. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - File '%s' could not be opened for reading.\n",
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
                  "ReadDataFileToArray - File '%s' has %d lines.\n",
                  filename, N);

        /* (3b) Allocate array on rank 0. */
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            fclose(fp);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed for array.\n");
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
                  "ReadDataFileToArray - Successfully read %d values from '%s'.\n",
                  N, filename);
    }

    /* STEP 4: Broadcast the integer N to all ranks. */
    ierr = MPI_Bcast(&N, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    /* STEP 5: Each rank allocates an array to receive the broadcast if rank>0. */
    if (rank) {
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed on rank %d.\n",
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
              "ReadDataFileToArray - Done. Provided array of length=%d to all ranks.\n",
              N);
    return 0; /* success */
}

/**
 * @brief Writes a data block in appended format for a VTK file, including a 4-byte size prefix.
 *
 * This function writes the size of the data block (in bytes) as a 4-byte integer, 
 * followed immediately by the raw data bytes. It returns 0 on success, or 1 on error.
 *
 * @param fp            File pointer to write to (must be open for binary write).
 * @param data          Pointer to the data to be written.
 * @param num_elements  Number of elements in the data array.
 * @param element_size  Size (in bytes) of each element.
 *
 * @return PetscInt  Returns 0 if successful, non-zero otherwise.
 */
PetscErrorCode WriteVTKAppendedBlock(FILE *fp, const void *data, PetscInt num_elements, size_t element_size) {

    // Log the function call with parameters
  // LOG_ALLOW_SYNC(LOCAL,LOG_INFO,"WriteVTKAppendedBlock - Called with %d elements, each of size %zu bytes.\n",num_elements, element_size);

    // Calculate the block size
    PetscInt block_size = num_elements * (PetscInt)element_size;
// LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "WriteVTKAppendedBlock - Calculated block size: %d bytes.\n", block_size);

    // Write the block size as a 4-byte integer
    if (fwrite(&block_size, sizeof(PetscInt), 1, fp) != 1) {
        fprintf(stderr, "[ERROR] Failed to write block size.\n");
	//    LOG_ALLOW_SYNC(LOCAL, LOG_ERROR, "WriteVTKAppendedBlock - Error writing block size.\n");
        return 1;
    }

    // Write the actual data
    if (fwrite(data, element_size, num_elements, fp) != (size_t)num_elements) {
        fprintf(stderr, "[ERROR] Failed to write data block.\n");
	//   LOG_ALLOW_SYNC(LOCAL, LOG_ERROR, "WriteVTKAppendedBlock - Error writing data block.\n");
        return 1;
    }

    // Log success
// LOG_ALLOW_SYNC(LOCAL, LOG_INFO, "WriteVTKAppendedBlock - Successfully wrote block of %d bytes.\n", block_size);

    return 0; // Success
}

/**
 * @brief Writes the XML header portion of a .vts file, defining the structured grid extents and data arrays.
 *
 * This function prPetscInts XML tags for the grid, including whole and piece extents,
 * configures the file's byte order and floating-point precision, and sets up
 * appended data sections (for points and a single scalar field). The offsets in
 * the appended data are updated accordingly.
 *
 * @param[in]  fp         File pointer (already open) for writing the .vts header.
 * @param[in]  mx         Number of cells in the x-direction (plus ghost cells).
 * @param[in]  my         Number of cells in the y-direction (plus ghost cells).
 * @param[in]  mz         Number of cells in the z-direction (plus ghost cells).
 * @param[in]  nnodes     Total number of grid nodes.
 * @param[in]  fieldName  Name of the scalar field written to PointData.
 * @param[in]  numComponents Number of components in the fieldName array (e.g., 3 for velocity).
 * @param[in]  boffset    Current byte offset in the appended data section.
 * @param[out] boffsetOut Updated byte offset after writing the header information.
 *
 * @return 0 on success, non-zero on failure (though no explicit failure cases are handled here).
 */
static PetscInt WriteVTSXMLHeader(FILE       *fp,
                             PetscInt         mx,
                             PetscInt         my,
                             PetscInt         mz,
                             PetscInt         nnodes,
                             const char *fieldName,
                             PetscInt         numComponents,
                             PetscInt         boffset,
                             PetscInt        *boffsetOut)
{
    // Log the entry into this function and the parameters
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Writing .vts header (fieldName=%s, numComponents=%d,boffset=%d,mx=%d, my=%d, mz=%d, nnodes=%d).\n",fieldName,numComponents,boffset,
              mx, my, mz, nnodes);

    // Set XML configuration strings
    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    // Begin XML header
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",
            byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
            0, mx-1, 0, my-1, 0, mz-1);
    fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n",
            0, mx-1, 0, my-1, 0, mz-1);

    // Points section
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, boffset);
    boffset += (PetscInt)sizeof(PetscInt) + 3 * nnodes * (PetscInt)sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Updated offset to %d after Points.\n", boffset);
    fprintf(fp, "      </Points>\n");

    // PointData => "ufield,vfield,p
    fprintf(fp, "      <PointData Scalars=\"%s\">\n", fieldName);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, fieldName, numComponents,boffset);
    boffset += (PetscInt)sizeof(PetscInt) + nnodes * numComponents * (PetscInt)sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Updated offset to %d after PointData.\n", boffset);
    fprintf(fp, "      </PointData>\n");

    // Close the piece and structured grid tags
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");

    // Open the appended data section
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");

    // Return new offset
    *boffsetOut = boffset;

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Completed writing .vts header. New boffset=%d.\n", boffset);

    return 0;
}

/**
 * @brief Finalizes the .vts file by closing the appended data section and the VTKFile XML tag.
 *
 * This function writes the closing XML tags for the appended data section and the 
 * overall VTK file, ensuring a well-formed .vts file structure.
 *
 * @param[in] fp  File pointer (already open) for writing the .vts footer.
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTSXMLFooter(FILE *fp)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Closing .vts XML.\n");

    fprintf(fp, "\n </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Successfully wrote .vts footer.\n");

    return 0;
}


/**
 * @brief Writes the XML header portion of a .vtp file for a point-cloud representation.
 *
 * This function sets up the appended data sections for point coordinates, a named field
 * (e.g., velocity) with a specified number of components, and the connectivity data for
 * each point as a vertex. The byte offset for the appended data is updated after each
 * section, ensuring correct placement of raw binary data in the final output file.
 *
 * @param[in]  fp            File pointer for writing the XML header (already open).
 * @param[in]  npoints       Number of points in the point-cloud.
 * @param[in]  fieldName     Name of the point data field (e.g., "velocity").
 * @param[in]  numComponents Number of components in the fieldName array (e.g., 3 for velocity).
 * @param[in]  boffset       Current byte offset into the appended data section.
 * @param[out] boffsetOut    Updated byte offset after writing the header information.
 *
 * @return PetscInt  Returns 0 on success.
 */
/**
 * @brief Writes the XML header portion of a .vtp file for a point-cloud representation.
 *
 * This function sets up the appended data sections for point coordinates, a named field
 * (e.g., velocity) with a specified number of components, and the connectivity/offsets data for
 * each point as a vertex. The byte offset for the appended data is updated after each
 * section. The data type for connectivity and offsets ("Int32" or "Int64") is now
 * determined dynamically based on the size of PetscInt during compilation.
 *
 * @param[in]  fp            File pointer for writing the XML header (already open).
 * @param[in]  npoints       Number of points in the point-cloud.
 * @param[in]  fieldName     Name of the point data field (e.g., "velocity").
 * @param[in]  numComponents Number of components in the fieldName array (e.g., 3 for velocity).
 * @param[in]  boffset       Current byte offset into the appended data section.
 * @param[out] boffsetOut    Updated byte offset after writing the header information.
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTPXMLHeader(FILE       *fp,
                                  PetscInt    npoints,
                                  const char *fieldName,
                                  PetscInt    numComponents,
                                  PetscInt    boffset,
                                  PetscInt   *boffsetOut)
{
    // Log function entry with parameters
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Called with npoints=%d, fieldName=%s, numComponents=%d, boffset=%d.\n",
              npoints, fieldName, numComponents, boffset);

    // Determine data types based on compilation settings
    const char *precision     = "Float64"; // Assuming coordinates and field data are doubles
    const char *int_type_str  = (sizeof(PetscInt) == 8) ? "Int64" : "Int32"; // Dynamic check for PetscInt size

    // Pre-calculate byte sizes for data blocks
    PetscInt points_bytes       = npoints * 3 * sizeof(double); // 3 components, double precision
    PetscInt field_data_bytes   = npoints * numComponents * sizeof(double); // numComponents, double precision
    PetscInt connectivity_bytes = npoints * sizeof(PetscInt); // Based on actual PetscInt size
    PetscInt offsets_bytes      = npoints * sizeof(PetscInt); // Based on actual PetscInt size

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Determined types: precision=%s, int_type=%s (sizeof(PetscInt)=%zu)\n",
              precision, int_type_str, sizeof(PetscInt));

    // 1) Standard XML + <VTKFile> opening
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <PolyData>\n");

    // For simple point-cloud => #Verts = npoints
    fprintf(fp,
            "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\" NumberOfLines=\"0\" NumberOfPolys=\"0\" NumberOfStrips=\"0\">\n",
            npoints, npoints);

    // 2) Points => offset=boffset
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\"/>\n",
                precision, boffset);
    boffset += sizeof(PetscInt) + points_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after Points (%d bytes data).\n", boffset, points_bytes);
    fprintf(fp, "      </Points>\n");

    // 3) PointData => e.g. "velocity" => next offset
    fprintf(fp, "      <PointData>\n");
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" "
            "format=\"appended\" offset=\"%d\"/>\n",
            precision, fieldName, numComponents, boffset);
    boffset += sizeof(PetscInt) + field_data_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after PointData ('%s', %d bytes data).\n", boffset, fieldName, field_data_bytes);
    fprintf(fp, "      </PointData>\n");

    // 4) Verts => connectivity + offsets
    fprintf(fp, "      <Verts>\n");
    // Use the dynamically determined integer type string
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n",
            int_type_str, boffset);
    boffset += sizeof(PetscInt) + connectivity_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after connectivity (%s, %d bytes data).\n", boffset, int_type_str, connectivity_bytes);

    // Use the dynamically determined integer type string
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n",
            int_type_str, boffset);
    boffset += sizeof(PetscInt) + offsets_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after offsets (%s, %d bytes data).\n", boffset, int_type_str, offsets_bytes);
    fprintf(fp, "      </Verts>\n");

    // 5) Close out the piece + PolyData
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");

    // 6) <AppendedData encoding="raw"> + underscore
    //    => all binary must follow AFTER this line
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_"); // This underscore MUST be present before the binary data blocks

    // Assign the final calculated offset to the output parameter
    *boffsetOut = boffset;

    // Log function completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Completed writing VTP XML header. Final calculated boffset (before actual data write)=%d.\n",
              boffset);

    return 0; // Success
}


/**
 * @brief Finalizes the .vtp file by closing the appended data section and the VTKFile tag.
 *
 * This function writes the closing XML tags for the appended data section 
 * and completes the overall VTK structure, ensuring a well-formed .vtp file.
 *
 * @param[in] fp  File pointer to the .vtp file (open for writing).
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTPXMLFooter(FILE *fp)
{
    // Log the function entry
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Closing VTP XML.\n");

    // no stray prPetscInts or data -> directly close
    fprintf(fp, "\n  </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Completed writing VTP XML footer.\n");

    return 0;
}

/**
 * @brief Writes the initial XML header for a VTK file based on the provided metadata.
 *
 * This function writes the XML header for VTK output. It now supports both polydata
 * (VTK_POLYDATA) and structured grid (VTK_STRUCTURED) file types. For each type, the
 * function determines whether to write a scalar or vector field header by inspecting
 * the metadata. For polydata, it calls \c WriteVTPXMLHeader and for structured grids,
 * it calls \c WriteVTSXMLHeader.
 *
 * @param[in]  fp          File pointer (open for writing) to which the header will be written.
 * @param[in]  meta        Pointer to a \c VTKMetaData structure containing fileType, number of points,
 *                         grid dimensions, field names, and other necessary information.
 * @param[in]  boffset     Current byte offset in the appended data section.
 * @param[out] boffsetOut  Updated byte offset after writing any header-related data.
 *
 * @return PetscInt Returns 0 on success, or -1 if the fileType is not handled.
 *
 * @note The function uses PETSc-style logging (LOG_ALLOW) to trace its execution.
 *       For structured grid output, it assumes that grid dimensions (mx, my, mz) and
 *       the number of nodes (nnodes) have been set in the metadata.
 */
static PetscInt WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt boffset, PetscInt *boffsetOut)
{
    const char *fieldName = NULL;
    PetscInt numComponents = 1;

    /*-----------------------------------------------------------------------
     * Log the entry into this function with the key metadata parameters.
     *-----------------------------------------------------------------------*/
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTKFileHeader - Entered: fileType=%d, npoints=%d, scalarField=%s, vectorField=%s, "
              "numVectorFields=%d, mx=%d, my=%d, mz=%d, nnodes=%d.\n",
              meta->fileType, meta->npoints, meta->scalarFieldName, meta->vectorFieldName,
              meta->numVectorFields, meta->mx, meta->my, meta->mz, meta->nnodes);

    if (meta->fileType == VTK_POLYDATA) {
        /*===================================================================
         * Polydata Branch: VTK_POLYDATA
         *===================================================================*/
        LOG_ALLOW(GLOBAL, LOG_INFO, "WriteVTKFileHeader - Processing VTK_POLYDATA output.\n");

        /* Determine the field type for polydata:
         * If a vector field is available, use it (with 3 components);
         * otherwise, use the scalar field (with 1 component).
         */
        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components for polydata.\n",
                      fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component for polydata.\n",
                      fieldName);
        }

        /* Log the call to the polydata header-writing function */
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTPXMLHeader with boffset=%d.\n", boffset);
        return WriteVTPXMLHeader(fp,
                                 meta->npoints,
                                 fieldName,
                                 numComponents,
                                 boffset,
                                 boffsetOut);
    } else if (meta->fileType == VTK_STRUCTURED) {
        /*===================================================================
         * Structured Grid Branch: VTK_STRUCTURED
         *===================================================================*/
        LOG_ALLOW(GLOBAL, LOG_INFO, "WriteVTKFileHeader - Processing VTK_STRUCTURED output.\n");

        /* Determine the field type for structured grid:
         * If a vector field is available, use it (with 3 components);
         * otherwise, use the scalar field (with 1 component).
         */
        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components for structured grid.\n",
                      fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component for structured grid.\n",
                      fieldName);
        }

        /* Log the grid dimensions and number of nodes, then call the structured grid header-writing function */
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTSXMLHeader with boffset=%d, grid dimensions (%d, %d, %d) and nnodes=%d.\n",
                  boffset, meta->mx, meta->my, meta->mz, meta->nnodes);
        return WriteVTSXMLHeader(fp,
                                 meta->mx,
                                 meta->my,
                                 meta->mz,
                                 meta->nnodes,
                                 fieldName,
                                 numComponents,
                                 boffset,
                                 boffsetOut);
    } else {
        /*-------------------------------------------------------------------
         * Unsupported File Type: Log a warning and return an error.
         *-------------------------------------------------------------------*/
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "WriteVTKFileHeader - Unsupported file type %d encountered.\n", meta->fileType);
        return -1;
    }
}



/**
 * @brief Writes the XML footer for a VTK file based on the provided metadata.
 *
 * This function currently handles only \c VTK_POLYDATA. If the file type is 
 * \c VTK_POLYDATA, it delegates to \c WriteVTPXMLFooter. Otherwise, 
 * it logs a warning and returns -1 (indicating that other file types are not yet implemented).
 *
 * @param[in] fp    File pointer (open for writing) to which the footer will be written.
 * @param[in] meta  Pointer to a \c VTKMetaData structure containing the fileType 
 *                  (e.g., \c VTK_POLYDATA or \c VTK_STRUCTURED).
 *
 * @return PetscInt  Returns 0 or the value from \c WriteVTPXMLFooter on success, 
 *              or -1 if the file type is not handled.
 */
static PetscInt WriteVTKFileFooter(FILE *fp, const VTKMetaData *meta)
{
    // Log the entry into this function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, 
              "WriteVTKFileFooter - Called with fileType=%d.\n", meta->fileType);

    if (meta->fileType == VTK_POLYDATA) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTP XML footer.\n");
        return WriteVTPXMLFooter(fp);
    } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTS XML footer.\n");
	return WriteVTSXMLFooter(fp);
    }
}

/**
 * @brief Creates and writes a VTK file (either .vts or .vtp) based on the provided metadata.
 *
 * This function gathers the necessary information from the \c VTKMetaData structure (e.g., \c coords,
 * \c scalarField, \c vectorField, and connectivity or offsets if it is \c VTK_POLYDATA). It writes
 * the XML header, followed by the appended binary data blocks in the correct order, and finally the
 * XML footer. The I/O occurs only on \c rank 0 of the provided \c MPI_Comm.
 *
 * @param[in]  filename  The output file name (e.g., "output.vtp" or "output.vts").
 * @param[in]  meta      Pointer to a \c VTKMetaData structure containing all necessary fields.
 * @param[in]  comm      The MPI communicator used for parallel execution.
 *
 * @return PetscInt  Returns 0 on success, or 1 on file-opening failure. Always returns 0 on non-\c rank 0 processes.
 */
PetscInt CreateVTKFileFromMetadata(const char *filename,
                                   const VTKMetaData *meta,
                                   MPI_Comm comm)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Entry point. Filename: %s\n", filename);

    PetscMPIInt rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - rank=%d, size=%d.\n", rank, size);

    // 1) Only rank 0 writes the file
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "CreateVTKFileFromMetadata - Rank 0 writing file '%s'.\n", filename);

        FILE *fp = fopen(filename, "wb");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_ERROR, "CreateVTKFileFromMetadata - fopen failed for %s.\n", filename);
            return 1;
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Successfully opened file: %s\n", filename);

        PetscInt boffset = 0;
        // 2) Write the XML header => sets <DataArray ... offset="...">
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing header (initial boffset=%d).\n", boffset);
        WriteVTKFileHeader(fp, meta, boffset, &boffset);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Header written (updated boffset=%d).\n", boffset);

        // 3) Immediately write appended data blocks in EXACT order:
        //    (a) coords (3*npoints doubles)
        if (meta->coords) {
	  PetscInt ncoords = 3 * (meta->fileType == VTK_STRUCTURED ? meta->nnodes : meta->npoints);
            if (ncoords > 0) {
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing coords block: %d doubles.\n", ncoords);
                WriteVTKAppendedBlock(fp, meta->coords, ncoords, sizeof(double));
                boffset += sizeof(PetscInt) + ncoords * sizeof(double);
            }
        }

        //    (b) scalar or vector field
        if (meta->numScalarFields > 0 && meta->scalarField) {
	  PetscInt nvals = (meta->fileType == VTK_STRUCTURED ? meta->nnodes : meta->npoints);
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing scalar field block: %d doubles.\n", nvals);
            WriteVTKAppendedBlock(fp, meta->scalarField, nvals, sizeof(double));
            boffset += sizeof(PetscInt) + nvals * sizeof(double);
        }
        else if (meta->numVectorFields > 0 && meta->vectorField) {
	  PetscInt nvals = 3 * (meta->fileType == VTK_STRUCTURED ? meta->nnodes : meta->npoints); // 3 components
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing vector field block: %d doubles.\n", nvals);
            WriteVTKAppendedBlock(fp, meta->vectorField, nvals, sizeof(double));
            boffset += sizeof(PetscInt) + nvals * sizeof(double);
        }

        //    (c) connectivity (if polydata)
        if (meta->fileType == VTK_POLYDATA && meta->connectivity) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing connectivity block: %d PetscInts.\n", meta->npoints);
            WriteVTKAppendedBlock(fp, meta->connectivity, meta->npoints, sizeof(PetscInt));
            boffset += sizeof(PetscInt) + meta->npoints * sizeof(PetscInt);
        }

        //    (d) offsets (if polydata)
        if (meta->fileType == VTK_POLYDATA && meta->offsets) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing offsets block: %d PetscInts.\n", meta->npoints);
            WriteVTKAppendedBlock(fp, meta->offsets, meta->npoints, sizeof(PetscInt));
            boffset += sizeof(PetscInt) + meta->npoints * sizeof(PetscInt);
        }

        // 4) Write footer => </AppendedData> & </VTKFile>
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing footer.\n");
        WriteVTKFileFooter(fp, meta);

        fclose(fp);
        LOG_ALLOW(GLOBAL, LOG_INFO, "CreateVTKFileFromMetadata - File '%s' closed.\n", filename);
    }

    return 0;
}


/**
 * @brief Gathers a PETSc vector onto rank 0 as a contiguous array of doubles.
 *
 * This function retrieves the local portions of the input vector \p inVec from
 * all MPI ranks via \c MPI_Gatherv and assembles them into a single array on rank 0.
 * The global size of the vector is stored in \p N, and a pointer to the newly
 * allocated array is returned in \p arrayOut (valid only on rank 0).
 *
 * @param[in]  inVec      The PETSc vector to gather.
 * @param[out] N          The global size of the vector (output).
 * @param[out] arrayOut   On rank 0, points to the newly allocated array of size \p N.
 *                        On other ranks, it is set to NULL.
 *
 * @return PetscErrorCode  Returns 0 on success, or a non-zero PETSc error code.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut)
{
    PetscErrorCode    ierr;
    MPI_Comm          comm;
    PetscMPIInt       rank, size;
    PetscInt          globalSize, localSize;
    const PetscScalar *localArr = NULL;

    PetscFunctionBeginUser;
    
    // Log entry into the function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "VecToArrayOnRank0 - Start gathering vector onto rank 0.\n");

    /* Get MPI comm, rank, size */
    ierr = PetscObjectGetComm((PetscObject)inVec, &comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);

    /* Get global size (for the entire Vec) */
    ierr = VecGetSize(inVec, &globalSize);CHKERRQ(ierr);
    *N = globalSize;

    /* Get local size (portion on this rank) */
    ierr = VecGetLocalSize(inVec, &localSize);CHKERRQ(ierr);

    // Log vector sizes and process info
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "VecToArrayOnRank0 - rank=%d of %d, globalSize=%d, localSize=%d.\n",
              rank, size, globalSize, localSize);

    /* Access the local array data */
    ierr = VecGetArrayRead(inVec, &localArr);CHKERRQ(ierr);

    /*
       We'll gather the local chunks via MPI_Gatherv:
       - First, gather all local sizes into recvcounts[] on rank 0.
       - Then set up a displacement array (displs[]) to place each chunk in the correct spot.
       - Finally, gather the actual data.
    */
    PetscMPIInt *recvcounts = NULL;
    PetscMPIInt *displs     = NULL;

    if (!rank) {
        recvcounts = (PetscMPIInt *) malloc(size * sizeof(PetscMPIInt));
        displs     = (PetscMPIInt *) malloc(size * sizeof(PetscMPIInt));
    }

    /* Convert localSize (PetscInt) to PetscMPIInt for MPI calls */
    PetscMPIInt localSizeMPI = (PetscMPIInt)localSize;

    /* Gather local sizes to rank 0 */
    ierr = MPI_Gather(&localSizeMPI, 1, MPI_INT,
                      recvcounts,      1, MPI_INT,
                      0, comm);CHKERRQ(ierr);

    /* On rank 0, build displacements and allocate the big array */
    if (!rank) {
        displs[0] = 0;
        for (PetscMPIInt i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        /* Allocate a buffer for the entire (global) array */
        *arrayOut = (double *) malloc(globalSize * sizeof(double));
        if (!(*arrayOut)) SETERRQ(comm, PETSC_ERR_MEM, "Failed to allocate array on rank 0.");
    } else {
        /* On other ranks, we do not allocate anything */
        *arrayOut = NULL;
    }

    /* Gather the actual data (assuming real scalars => MPI_DOUBLE) */
    ierr = MPI_Gatherv((void *) localArr,         /* sendbuf on this rank */
                       localSizeMPI, MPI_DOUBLE,  /* how many, and type */
                       (rank == 0 ? *arrayOut : NULL),  /* recvbuf on rank 0 */
                       (rank == 0 ? recvcounts : NULL),
                       (rank == 0 ? displs    : NULL),
                       MPI_DOUBLE, 0, comm);CHKERRQ(ierr);

    /* Restore local array (cleanup) */
    ierr = VecRestoreArrayRead(inVec, &localArr);CHKERRQ(ierr);

    if (!rank) {
        free(recvcounts);
        free(displs);
    }

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_INFO, "VecToArrayOnRank0 - Successfully gathered data on rank 0.\n");

    PetscFunctionReturn(0);
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
 * @brief Reads multiple fields (positions, velocity, CellID, and weight) into a DMSwarm.
 *
 * This function is analogous to ReadSimulationFields() but targets a DMSwarm. 
 * Each Swarm field is read from a separate file using ReadSwarmField().
 * 
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm (user->swarm).
 * @param[in]     ti   Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
/**
 * @brief Reads multiple fields (positions, velocity, etc.) into a DMSwarm.
 *
 * ASSUMES the swarm has already been resized (if necessary) by a function
 * like PreCheckAndResizeSwarm to match the particle count in the input files
 * for the current timestep.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm (user->swarm).
 * @param[in]     ti   Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure. Handles
 *         missing optional files as warnings.
 */
PetscErrorCode ReadAllSwarmFields(UserCtx *user, PetscInt ti)
{
  PetscErrorCode ierr;
  PetscInt nGlobal;
  PetscBool fileExists;

  PetscFunctionBeginUser;
  ierr = DMSwarmGetSize(user->swarm, &nGlobal); CHKERRQ(ierr); // Get final size
  LOG_ALLOW(GLOBAL, LOG_INFO, "ReadAllSwarmFields - Reading DMSwarm fields for timestep %d (expected size %d).\n", ti, nGlobal);

  // If the swarm is empty after potential resize, we might still technically proceed
  // but ReadSwarmField might do nothing or error depending on implementation.
  // It's safer to just return if empty.
  if (nGlobal == 0) {
      LOG_ALLOW(GLOBAL, LOG_INFO, "ReadAllSwarmFields - Swarm is empty for timestep %d. Nothing to read.\n", ti);
      PetscFunctionReturn(0);
  }

  // 1) Read positions (required, should succeed now)
  ierr = ReadSwarmField(user, "position", ti, "dat");
  if (ierr == PETSC_ERR_FILE_OPEN) {
      // This shouldn't happen if PreCheckAndResizeSwarm ran successfully, but handle defensively.
      LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Position file missing for step %d, but PreCheck succeeded? Inconsistent state.\n", ti);
      PetscFunctionReturn(PETSC_ERR_PLIB);
  } else if (ierr) {
      LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Failed to read position field for step %d (Error %d), even though size should match.\n", ti, ierr);
      CHKERRQ(ierr); // Treat other errors as fatal
  } else {
       LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadAllSwarmFields - Successfully read position field for step %d.\n", ti);
  }

  // 2) Read velocity
  ierr = CheckDataFile(ti, "velocity", "dat", &fileExists);
  if(fileExists) ierr = ReadSwarmField(user, "velocity", ti, "dat");
   if (ierr == PETSC_ERR_FILE_OPEN) { // Handle potentially missing file as warning
        LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - Velocity file for step %d not found. Velocity data may be invalid/uninitialized.\n", ti);
   } else if (ierr) {
       LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Failed to read velocity field for step %d (Error %d). File might be corrupt or size incorrect despite pre-check.\n", ti, ierr);
       CHKERRQ(ierr); // Treat other errors as fatal
   } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadAllSwarmFields - Successfully read velocity field for step %d.\n", ti);
   }

  // 3) Read pos_phy
  ierr = CheckDataFile(ti, "pos_phy", "dat", &fileExists); 
  if(fileExists) ierr = ReadSwarmField(user, "pos_phy", ti, "dat");
   if (ierr == PETSC_ERR_FILE_OPEN) { // Handle potentially missing file as warning
        LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - Physical position file for step %d not found. physical position data may be invalid/uninitialized.\n", ti);
   } else if (ierr) {
       LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Failed to read physical position field for step %d (Error %d). File might be corrupt or size incorrect despite pre-check.\n", ti, ierr);
       CHKERRQ(ierr); // Treat other errors as fatal
   } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadAllSwarmFields - Successfully read physical position field for step %d.\n", ti);
   }
   
  // 3) Read CellID (Optional)
  // ierr = CheckDataFile(ti, "velocity", ext, &fileExists); 
  // if(fileExists)ierr = ReadSwarmField(user, "DMSwarm_CellID", ti, "dat");
  // if (ierr == PETSC_ERR_FILE_OPEN) { LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - CellID file for step %d not found.\n", ti); }
  // else CHKERRQ(ierr); // Treat other errors as fatal unless CellID read can also fail gracefully

  // 4) Read weight (Optional)
  // ierr = ReadSwarmField(user, "weight", ti, "dat");
  // if (ierr == PETSC_ERR_FILE_OPEN) { LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - Weight file for step %d not found.\n", ti); }
  // else CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_INFO, "ReadAllSwarmFields - Finished reading DMSwarm fields for timestep %d.\n", ti);
  PetscFunctionReturn(0);
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

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadPositionsFromFile - Creating coords Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

  // For example: "position" is the name of the coordinate data
  ierr = ReadFieldData(user, "position", coordsVec, timeIndex, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "ReadPositionsFromFile - Error reading position data (ti=%d).\n",
              timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadPositions - Gathering coords Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(coordsVec, Ncoords, coordsArray);CHKERRQ(ierr);

  ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "ReadPositionsFromFile - Successfully gathered coordinates. Ncoords=%d.\n", *Ncoords);
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

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadFieldDataWrapper - Creating field Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

  ierr = ReadFieldData(user, fieldName, fieldVec, timeIndex, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "ReadFieldDataWrapper - Error reading field '%s' (ti=%d).\n",
              fieldName, timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadFieldDataWrapper - Gathering field Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(fieldVec, Nscalars, scalarArray);CHKERRQ(ierr);

  ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "ReadFieldDataWrapper - Successfully gathered field '%s'. Nscalars=%d.\n",
            fieldName, *Nscalars);
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
            LOG_ALLOW(PETSC_COMM_SELF, LOG_WARNING, "DisplayBanner (Rank 0) - bboxlist not provided or num_mpi_procs <=0; using user->bbox for domain bounds.\n");
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
        ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Initialization Mode: %d\n", simCtx->ParticleInitialization); CHKERRQ(ierr);
        if (simCtx->ParticleInitialization == 0) {
            if (user->inletFaceDefined) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : %s (Enum Val: %d)\n", BCFaceToString(user->identifiedInletBCFace), user->identifiedInletBCFace); CHKERRQ(ierr);
            } else {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : --- (No INLET face identified)\n"); CHKERRQ(ierr);
            }
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, " Field Initialization Mode   : %d\n", simCtx->FieldInitialization); CHKERRQ(ierr);
        if (simCtx->FieldInitialization == 1) {
	  ierr = PetscPrintf(PETSC_COMM_SELF, " Constant Velocity           : x - %.4f, y - %.4f, z - %.4f \n", (double)simCtx->InitialConstantContra.x,(double)simCtx->InitialConstantContra.y,(double)simCtx->InitialConstantContra.z ); CHKERRQ(ierr);
        }
	
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}
