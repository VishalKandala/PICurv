#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H

#include "io.h"
#include "variables.h"
#include "logging.h" 
#include "ParticleSwarm.h"
#include "interpolation.h"
#include "grid.h"
#include "setup.h"
#include "Metric.h"
#include "postprocessing_kernels.h"
#include "vtk_io.h"
/* --------------------------------------------------------------------
   postprocessor.h

   This header declares the interface for the post-processing executable
   or library. Typically, you'd have a function that runs your main 
   post-processing routine, or you might declare other helper functions.

   Here, we declare a single function: PostprocessMain, 
   which could be your main entry point if you want to call it from
   another place (or you might just put main() in postprocess.c).

*/

/**
 * @brief Creates a new, dedicated DMSwarm for post-processing tasks.
 *
 * This function is called once at startup. It creates an empty DMSwarm and
 * associates it with the same grid DM as the primary swarm and registers all the required fields.
 * @param user The UserCtx where user->post_swarm will be created.
 * @param pps The PostProcessParams containing the particle_pipeline string for field registration.
 * @return PetscErrorCode
 */
PetscErrorCode SetupPostProcessSwarm(UserCtx* user, PostProcessParams* pps);

/**
 * @brief Orchestrates the writing of a combined, multi-field VTK file for a single time step.
 *
 * This function is the primary driver for generating output. It performs these steps:
 * 1. Prepares the subsampled coordinate array required for the legacy grid format.
 * 2. Parses the user-requested list of fields from the configuration.
 * 3. For each field, prepares a corresponding subsampled data array.
 * 4. Assembles all prepared arrays into a single VTKMetaData struct.
 * 5. Calls the low-level VTK writer to generate the final .vts file.
 * 6. Frees all temporary memory allocated during the preparation phase.
 *
 * @param user The UserCtx for the finest grid level.
 * @param pps  The post-processing configuration struct.
 * @param ti   The current time step index.
 * @return PetscErrorCode
 */
PetscErrorCode WriteEulerianFile(UserCtx* user, PostProcessParams* pps, PetscInt ti);

/**
 * @brief Parses the processing pipeline string and executes the requested kernels.
 * @param user The UserCtx containing the data to be transformed.
 * @param config The PostProcessConfig containing the pipeline string.
 * @return PetscErrorCode
 */
PetscErrorCode EulerianDataProcessingPipeline(UserCtx* user, PostProcessParams* pps);


/**
 * @brief Parses and executes the particle pipeline using a robust two-pass approach.
 *
 * This function ensures correctness and efficiency by separating field registration
 * from kernel execution.
 *
 * PASS 1 (Registration): The pipeline string is parsed to identify all new fields
 * that will be created. These fields are registered with the DMSwarm.
 *
 * Finalize: After Pass 1, DMSwarmFinalizeFieldRegister is called exactly once if
 * any new fields were added, preparing the swarm's memory layout.
 *
 * PASS 2 (Execution): The pipeline string is parsed again, and this time the
 * actual compute kernels are executed, filling the now-valid fields.
 *
 * @param user The UserCtx containing the DMSwarm.
 * @param pps  The PostProcessParams struct containing the particle_pipeline string.
 * @return PetscErrorCode
 */
PetscErrorCode ParticleDataProcessingPipeline(UserCtx* user, PostProcessParams* pps);

/**
 * @brief Writes particle data to a VTP file using the Prepare-Write-Cleanup pattern.
 */
PetscErrorCode WriteParticleFile(UserCtx* user, PostProcessParams* pps, PetscInt ti);

#endif /* POSTPROCESSOR_H */

