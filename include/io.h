#ifndef IO_H
#define IO_H

#include "variables.h" // Essential for SimCtx and UserCtx definitions
#include "logging.h"        // For logging macros.
#include "Boundaries.h" 
/**
 * @file io.h
 * @brief Public interface for data input/output routines.
 *
 * This header declares functions responsible for parsing grid geometry
 * information, either from command-line options for programmatically
 * generated grids or by reading the header of a grid definition file.
 */

/**
 * @brief Parses command-line options for a programmatically generated grid for a SINGLE block.
 *
 * This function reads all per-block array options related to grid geometry,
 * such as dimensions (-im), domain bounds (-xMins), and stretching ratios (-rxs).
 * It then populates the fields of the provided `UserCtx` struct using its
 * internal block index `user->_this`.
 *
 * @param user Pointer to the `UserCtx` for a specific block. The function will
 *             populate the geometric fields (`IM`, `Min_X`, `rx`, etc.) within this struct.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ReadGridGenerationInputs(UserCtx *user);

/**
 * @brief Parses grid resolution arrays (`-im`, `-jm`, `-km`) once and applies them to all finest-grid blocks.
 *
 * This helper centralizes one-time resolution ingestion for analytical grid setup.
 * It fills `IM/JM/KM` in each element of the finest-level `UserCtx` array.
 *
 * @param finest_users Pointer to the finest-level `UserCtx` array (length `nblk`).
 * @param nblk Number of blocks in the finest-level array.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode PopulateFinestUserGridResolutionFromOptions(UserCtx *finest_users, PetscInt nblk);

/**
 * @brief Sets grid dimensions from a file for a SINGLE block using a one-time read cache.
 *
 * This function uses a static-variable pattern to ensure the grid file header
 * is read only once, collectively, by all processes on the first call.
 * Subsequent calls simply retrieve the pre-loaded and broadcasted data for
 * the specified block.
 *
 * @param user Pointer to the `UserCtx` for a specific block. This function will
 *             populate the `IM`, `JM`, and `KM` fields.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ReadGridFile(UserCtx *user);

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
PetscErrorCode VerifyPathExistence(const char *path, PetscBool is_dir, PetscBool is_optional, const char *description, PetscBool *exists);

/**
 * @brief Returns whether full field/restart output should be written for the
 *
 * completed timestep.
 *
 * @param simCtx Simulation context controlling the operation.
 * @param completed_step Completed step index used by the decision helper.
 * @return PetscBool indicating the result of `ShouldWriteDataOutput()`.
 */
PetscBool ShouldWriteDataOutput(const SimCtx *simCtx, PetscInt completed_step);

/**
 * @brief Reads binary field data for velocity, pressure, and other required vectors.
 *
 * Reads contravariant velocity (`Ucont`) from `vfield`, Cartesian velocity (`Ucat`) from `ufield`,
 * pressure (`P`), node state (`Nvert_o`), and optionally statistical quantities, LES, and RANS fields
 * from binary files. Logs missing files but continues execution.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti   Timestep index used to construct restart file names.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadSimulationFields(UserCtx *user, PetscInt ti);


/**
 * @brief Reads data for a specific field from a file into the provided vector.
 *
 * This function uses the field name to construct the file path and reads the data
 * from the corresponding file into the provided PETSc vector.
 *
 * @param[in]     user        Pointer to the UserCtx structure containing simulation context.
 * @param[in]     field_name  Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[out]    field_vec   PETSc vector to store the field data.
 * @param[in]     ti          Time index for constructing the file name.
 * @param[in]     ext         File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext);


/**
 * @brief Reads statistical fields used for time-averaged simulations.
 *
 * Reads statistical quantities such as velocity sums and pressure sums. Logs missing files and initializes
 * fields to zero if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti   Timestep index used to construct statistics file names.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadStatisticalFields(UserCtx *user,PetscInt ti);

/**
 * @brief Reads LES-related fields used in turbulence modeling.
 *
 * Reads the Smagorinsky constant (`Cs`) and transfers data to local vectors. Logs missing files and
 * initializes fields to zero if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadLESFields(UserCtx *user, PetscInt ti);

/**
 * @brief Reads RANS-related fields for turbulence modeling.
 *
 * Reads `K_Omega` fields (used in RANS modeling) and initializes them if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadRANSFields(UserCtx *user, PetscInt ti);

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
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext);


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
PetscErrorCode WriteSimulationFields(UserCtx *user);


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
PetscErrorCode WriteStatisticalFields(UserCtx *user);

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
PetscErrorCode WriteLESFields(UserCtx *user);

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
PetscErrorCode WriteRANSFields(UserCtx *user);

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
 * @note Compatible with PETSc 3.14.4.
 */
PetscErrorCode WriteSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext);

/**
 * @brief Writes integer data from a specific PETSc Swarm field to a file.
 *
 * This function is designed for swarm fields that store integer data (e.g.,
 * DMSwarm_CellID), which cannot be converted to a standard PETSc Vec of
 * PetscScalars. It accesses the raw data pointer for the field on each rank
 * using DMSwarmGetField(), writes the local data to a rank-specific binary file,
 * and then restores the field access.
 *
 * @param[in] user       Pointer to the UserCtx structure containing the PetscSwarm.
 * @param[in] field_name Name of the integer Swarm field to be written.
 * @param[in] ti         Time index used to construct the output file name.
 * @param[in] ext        File extension (e.g., "dat", "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteSwarmIntField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext);

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
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteAllSwarmFields(UserCtx *user);

/**
 * @brief Reads a simple ASCII data file containing one numeric value per line.
 *
 * This helper performs rank-0 file I/O, broadcasts the parsed result to the rest
 * of the communicator, and returns a replicated array on every rank.
 *
 * @param filename Path to the input data file.
 * @param data_out Output pointer to the allocated scalar array.
 * @param Nout Output pointer storing the number of values read.
 * @param comm MPI communicator used for the coordinated read/broadcast sequence.
 * @return Integer value produced by `ReadDataFileToArray()`.
 */
PetscInt ReadDataFileToArray(const char   *filename,
                            double      **data_out,
                            PetscInt      *Nout,
                            MPI_Comm       comm);

/**
 * @brief Creates a VTK file from prepared metadata and field payloads.
 *
 * This helper dispatches to the structured-grid or polydata writer based on the
 * metadata contents and emits the assembled VTK file on the requested communicator.
 *
 * @param filename Path to the output VTK file.
 * @param meta VTK metadata describing the output geometry and field payloads.
 * @param comm MPI communicator used by the write operation.
 * @return Integer value produced by `CreateVTKFileFromMetadata()`.
 */
PetscInt CreateVTKFileFromMetadata(const char       *filename,
                                   const VTKMetaData *meta,
                                   MPI_Comm           comm);

/**
 * @brief Gathers the contents of a distributed PETSc Vec into a single array on rank 0.
 *
 * @param[in]  inVec       The input (possibly distributed) Vec.
 * @param[out] N           The global size of the vector.
 * @param[out] arrayOut    On rank 0, points to the newly allocated array holding all data.
 *                         On other ranks, it is set to NULL.
 *
 * @return PetscErrorCode  Return 0 on success, nonzero on failure.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut);

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
PetscErrorCode SwarmFieldToArrayOnRank0(DM swarm, const char *field_name, PetscInt *n_total_particles, PetscInt *n_components, void **gathered_array);

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
PetscErrorCode ReadSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext);

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
PetscErrorCode ReadSwarmIntField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext);

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
PetscErrorCode ReadAllSwarmFields(UserCtx *user, PetscInt ti);

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
				     PetscInt *Ncoords);


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
				    PetscInt *Nscalars);

/**
 * @brief Displays a structured banner summarizing the simulation configuration.
 *
 * This function prints key simulation parameters to standard output.
 * It is intended to be called ONLY by MPI rank 0.
 * It retrieves global domain bounds and block metadata from `simCtx`.
 *
 * @param[in] simCtx Pointer to the master simulation context.
 *
 * @return PetscErrorCode  Returns `0` on success.
 */
PetscErrorCode DisplayBanner(SimCtx *simCtx);


// --- Conversion and Validation Utilities ---
// These are now public and can be used by other parts of the application.

/**
 * @brief Converts a face-token string (e.g., "-Xi", "+Eta") to the internal `BCFace` enum.
 *
 * @param[in]  str      Input token from configuration.
 * @param[out] face_out Parsed enum value on success.
 * @return PetscErrorCode 0 on success, non-zero for invalid tokens or null pointers.
 */
PetscErrorCode StringToBCFace(const char* str, BCFace* face_out);

/**
 * @brief Converts a mathematical BC type string (e.g., "PERIODIC", "WALL") to `BCType`.
 *
 * @param[in]  str      Input token from configuration.
 * @param[out] type_out Parsed enum value on success.
 * @return PetscErrorCode 0 on success, non-zero for invalid tokens or null pointers.
 */
PetscErrorCode StringToBCType(const char* str, BCType* type_out);

/**
 * @brief Converts a BC handler token (implementation strategy) to `BCHandlerType`.
 *
 * @param[in]  str         Input handler token from configuration.
 * @param[out] handler_out Parsed enum value on success.
 * @return PetscErrorCode 0 on success, non-zero for invalid tokens or null pointers.
 */
PetscErrorCode StringToBCHandlerType(const char* str, BCHandlerType* handler_out);

/**
 * @brief Validates that a selected handler is compatible with a mathematical BC type.
 *
 * @param[in] type    Mathematical BC type (e.g., WALL, PERIODIC).
 * @param[in] handler Selected handler implementation enum.
 * @return PetscErrorCode 0 if compatible, non-zero if the combination is invalid.
 */
PetscErrorCode ValidateBCHandlerForBCType(BCType type, BCHandlerType handler);

// --- Memory Management ---

/**
 * @brief Frees an entire linked list of boundary-condition parameters.
 *
 * @param[in,out] head Head pointer of the `BC_Param` list to destroy.
 */
void FreeBC_ParamList(BC_Param *head);

/**
 * @brief Searches a BC_Param linked list for a key and returns its value as a double.
 * @param params The head of the BC_Param linked list.
 * @param key The key to search for (case-insensitive).
 * @param[out] value_out The found value, converted to a PetscReal.
 * @param[out] found Set to PETSC_TRUE if the key was found, PETSC_FALSE otherwise.
 * @return 0 on success.
 */
PetscErrorCode GetBCParamReal(BC_Param *params, const char *key, PetscReal *value_out, PetscBool *found);

/**
 * @brief Searches a BC_Param linked list for a key and returns its value as a bool.
 * @param params The head of the BC_Param linked list.
 * @param key The key to search for (case-insensitive).
 * @param[out] value_out The found value, converted to a PetscBool.
 * @param[out] found Set to PETSC_TRUE if the key was found, PETSC_FALSE otherwise.
 * @return 0 on success.
 */
PetscErrorCode GetBCParamBool(BC_Param *params, const char *key, PetscBool *value_out, PetscBool *found);

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
 * 6.  It also sets the particle inlet lookup fields in `UserCtx`.
 *
 * @param[in,out] user               The main UserCtx struct where the final configuration
 *                                   for all ranks will be stored.
 * @param[in]     bcs_input_filename The path to the boundary conditions configuration file.
 * @return PetscErrorCode 0 on success, error code on failure.
 */
PetscErrorCode ParseAllBoundaryConditions(UserCtx *user, const char *bcs_input_filename);

/**
 * @brief Scans all block-specific boundary condition files to determine a globally
 *        consistent periodicity for each dimension, reusing the core type parser.
 *
 * This is a lightweight pre-parser intended to be called before DMDA creation.
 * It ensures that the periodicity setting is consistent across all blocks, which is a
 * physical requirement for the domain.
 *
 * 1. It collectively verifies that the mandatory BCS file for each block exists.
 * 2. On MPI rank 0, it then iterates through the files.
 * 3. For each line, it attempts to convert the type string to a BCType enum using the
 *    standard `StringToBCType` helper.
 * 4. If the conversion is successful AND the type is PERIODIC, it flags the corresponding face.
 * 5. If the conversion fails (e.g., for "WALL", "INLET", etc.), the error is cleared
 *    and the line is simply ignored, as it's not relevant to periodicity.
 * 6. It validates consistency (e.g., -Xi and +Xi match) and ensures all block files
 *    specify the same global periodicity.
 * 7. It broadcasts the final three flags (as integers 0 or 1) to all MPI ranks.
 * 8. All ranks update the i_periodic, j_periodic, and k_periodic fields in their SimCtx.
 *
 * @param[in,out] simCtx The master SimCtx struct, containing the bcs_files list and
 *                       where the final periodicity flags will be stored.
 * @return PetscErrorCode 0 on success, error code on failure.
 */
PetscErrorCode DeterminePeriodicity(SimCtx *simCtx);

/**
 * @brief Helper function to trim leading/trailing whitespace from a string.
 * @param str The string to trim in-place.
 */
void TrimWhitespace(char *str);

/**
 * @brief Initializes post-processing settings from a config file and command-line overrides.
 *
 * This function establishes the configuration for a post-processing run by:
 * 1. Setting hardcoded default values in the PostProcessParams struct.
 * 2. Reading a configuration file to override the defaults.
 * 3. Parsing command-line options (-startTime, -endTime, etc.) which can override
 *    both the defaults and the file settings.
 *
 * @param simCtx The pointer to the simulation context that contains the postprocessing file and struct.
 * @return PetscErrorCode
 */
PetscErrorCode ParsePostProcessingSettings(SimCtx *simCtx);

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
PetscErrorCode ParseScalingInformation(SimCtx *simCtx);

#endif // IO_H
