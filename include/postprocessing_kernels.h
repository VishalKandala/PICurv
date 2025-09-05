#ifndef POSTPROCESSING_KERNELS_H
#define POSTPROCESSING_KERNELS_H

#include "variables.h"
#include "logging.h"
#include "io.h" // For the UpdateLocalGhosts function prototype

// Function prototypes for post-processing kernels
PetscErrorCode ComputeNodalAverage(UserCtx* user, const char* in_field_name, const char* out_field_name);
PetscErrorCode ComputeQCriterion(UserCtx* user);
PetscErrorCode NormalizeRelativeField(UserCtx* user, const char* relative_field_name);

// Add more post-processing kernel prototypes as needed
#endif // POSTPROCESSING_KERNELS_H