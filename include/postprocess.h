#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "io.h"
#include "variables.h"
#include "logging.h" 
#include "ParticleSwarm.h"
#include "interpolation.h"
#include "grid.h"
#include "setup.h"
#include "Metric.h"
/* --------------------------------------------------------------------
   postprocess.h

   This header declares the interface for the post-processing executable
   or library. Typically, you'd have a function that runs your main 
   post-processing routine, or you might declare other helper functions.

   Here, we declare a single function: PostprocessMain, 
   which could be your main entry point if you want to call it from
   another place (or you might just put main() in postprocess.c).

*/


/**
 * @brief Gathers a PETSc vector or a DMSwarm field, prepares VTKMetaData, and writes to a VTK file.
 *
 * This function dynamically handles **both Eulerian (.vts) and Lagrangian (.vtp) data**:
 *  - **If the field exists in `DMSwarm`**, it **extracts it into a PETSc Vec**.
 *  - **If already a Vec**, it directly processes it.
 *  - **For `.vtp` files, it also gathers particle positions** (`DMSwarmPICField_coor`).
 *  - **For `.vts`, only the field is gathered** (structured grid).
 *
 * @param[in]  user       Pointer to user context (contains swarm and Eulerian fields).
 * @param[in]  fieldName  Name of the field (e.g., "Ucat" for vts, "Velocity" for vtp).
 * @param[in]  coorfieldName  Name of the coordinate field used (e.g., "positions" or "pos_phy".
 * @param[in]  timeIndex  Timestep index for filename (e.g., "Ucat_00010.vts").
 * @param[in]  outExt     File extension ("vts" for structured, "vtp" for particles).
 * @param[in]  prefix     File prefix(directory to store file in).
 * @param[in]  comm       MPI communicator (usually PETSC_COMM_WORLD).
 *
 * @return PetscErrorCode  0 on success, or an error code on failure.
 */
PetscErrorCode GatherAndWriteField(UserCtx *user,
                                   const char *fieldName,
				   const char *coorfieldName,
                                   PetscInt timeIndex,
                                   const char *outExt,
				   const char *prefix,
                                   MPI_Comm comm);

/**
 * @brief Writes up to four Eulerian fields (Ucat, P, Ucont, Nvert) using the
 *        GatherAndWriteField() wrapper. Each field becomes a separate file.
 *
 * @param[in]  user       Pointer to the UserCtx holding PETSc vectors.
 * @param[in]  timeIndex  Timestep index for filenames.
 * @param[in]  outExt     File extension (e.g., "vts" or "vtp").
 * @param[in]  prefix     File prefix(directory to store file in).
 *
 * @return PetscErrorCode  Returns 0 on success, error code on failure.
 */
PetscErrorCode WriteEulerianVTK(UserCtx *user, PetscInt timeIndex, const char *outExt,const char *prefix);


/**
 * @brief Writes particle (Lagrangian) data fields to separate .vtp files.
 *
 * This function loops over key DMSwarm fields (positions, velocity, etc.) and writes 
 * each as a separate .vtp file, ensuring modularity and flexibility.
 *
 * @param[in] user       Pointer to the user context (contains the particle swarm).
 * @param[in] timeIndex  Current timestep index (used in filenames).
 * @param[in] outExt     File extension ("vtp" for particles).
 * @param[in]  prefix     File prefix(directory to store file in).
 *
 * @return PetscErrorCode  0 on success, nonzero on failure.
 */
PetscErrorCode WriteParticleVTK(UserCtx *user, PetscInt timeIndex, const char *outExt,const char *prefix);

#endif /* POSTPROCESS_H */

