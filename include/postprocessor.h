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
PetscErrorCode RunProcessingPipeline(UserCtx* user, PostProcessParams* pps);

#endif /* POSTPROCESSOR_H */

