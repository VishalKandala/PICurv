/**
 * @file ParticleInitialConditions.c
 * @brief Configured initial values of particle-carried fields.
 *
 * Each field's expressions are compiled once by the shared engine (expression.c) and
 * evaluated per particle. Draws are keyed on the particle's identity, so a value does not
 * depend on rank count, particle order, or how often it is evaluated.
 */

#include "ParticleInitialConditions.h"

#include <math.h>
#include <stdint.h>

#include "logging.h"

/* ------------------------------------------------------------------------- */
/*                                   Plans                                   */
/* ------------------------------------------------------------------------- */

/** @brief Names a particle expression may use, in the order their values are supplied. */
static const char *const kParticleVariables[] = {"x", "y", "z", "xn", "yn", "zn", "pid", "t"};
#define PARTICLE_VARIABLE_COUNT ((PetscInt)(sizeof(kParticleVariables) / sizeof(kParticleVariables[0])))
#define PARTICLE_FIELD_MAX_COMPONENTS 3
#define PARTICLE_FIELD_EXPRESSION_LENGTH 4096
#define PARTICLE_FIELD_HISTOGRAM_BINS 10

/** @brief One configured field: its identity and one expression per component. */
typedef struct {
    ParticleFieldId   field;
    PetscInt          components;
    PicurvExpression *component[PARTICLE_FIELD_MAX_COMPONENTS];
} ParticleFieldBinding;

struct ParticleFieldPlan {
    PetscInt              count;
    ParticleFieldBinding *bindings;
};

/**
 * @brief Implementation of \ref ParticleFieldPlanCreate().
 * @see ParticleFieldPlanCreate()
 */
PetscErrorCode ParticleFieldPlanCreate(ParticleFieldPlan **plan)
{
    PetscInt  count = 0;
    PetscBool found = PETSC_FALSE;
    char      option[256];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(plan, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Plan output is required.");
    *plan = NULL;
    PetscCall(PetscSNPrintf(option, sizeof(option), "-particle_fields_count"));
    PetscCall(PetscOptionsGetInt(NULL, NULL, option, &count, &found));
    if (!found || count == 0) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    PetscCheck(count > 0 && count <= PARTICLE_FIELD_ID_COUNT, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "%s must be between 0 and %d (got %" PetscInt_FMT ").", option, (int)PARTICLE_FIELD_ID_COUNT, count);

    PetscCall(PetscNew(plan));
    PetscCall(PetscCalloc1(count, &(*plan)->bindings));
    (*plan)->count = count;
    for (PetscInt i = 0; i < count; ++i) {
        ParticleFieldBinding          *binding = &(*plan)->bindings[i];
        const ParticleFieldDescriptor *descriptor = NULL;
        char                           name[64];

        PetscCall(PetscSNPrintf(option, sizeof(option), "-particle_fields_%" PetscInt_FMT "_name", i));
        PetscCall(PetscOptionsGetString(NULL, NULL, option, name, sizeof(name), &found));
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "%s is missing.", option);
        PetscCall(ParticleFieldIdFromName(name, &binding->field));
        PetscCall(ParticleFieldGetDescriptor(binding->field, &descriptor));
        PetscCheck(descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_USER_INITIALIZE, PETSC_COMM_WORLD,
                   PETSC_ERR_ARG_WRONG,
                   "Particle field '%s' cannot be given an initial value: the runtime sets it from the "
                   "Eulerian fields or from particle location.", name);
        PetscCheck(descriptor->dimension.kind == FIELD_DIMENSION_FIXED && descriptor->data_type == PETSC_REAL,
                   PETSC_COMM_WORLD, PETSC_ERR_PLIB, "Settable particle field '%s' must be a real quantity.", name);
        for (PetscInt j = 0; j < i; ++j) {
            PetscCheck((*plan)->bindings[j].field != binding->field, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                       "Particle field '%s' is configured twice.", name);
        }
        binding->components = descriptor->components;
        PetscCheck(binding->components <= PARTICLE_FIELD_MAX_COMPONENTS, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "Particle field '%s' has more components than an initial value supports.", name);
        for (PetscInt c = 0; c < binding->components; ++c) {
            char text[PARTICLE_FIELD_EXPRESSION_LENGTH];

            PetscCall(PetscSNPrintf(option, sizeof(option), "-particle_fields_%" PetscInt_FMT "_expr_%" PetscInt_FMT, i, c));
            PetscCall(PetscOptionsGetString(NULL, NULL, option, text, sizeof(text), &found));
            PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "%s is missing.", option);
            PetscCall(PicurvExpressionCompile(text, kParticleVariables, PARTICLE_VARIABLE_COUNT, PETSC_TRUE,
                                              &binding->component[c]));
        }
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/** @brief Global domain bounds, in solver units, from the replicated rank bounding boxes. */
static PetscErrorCode DomainBounds(const SimCtx *simCtx, Cmpnts *lower, Cmpnts *upper)
{
    PetscMPIInt size = 1;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(simCtx->bboxlist, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Rank bounding boxes are needed before a particle initial value is applied.");
    size = simCtx->size;
    *lower = simCtx->bboxlist[0].min_coords;
    *upper = simCtx->bboxlist[0].max_coords;
    for (PetscInt r = 1; r < (PetscInt)size * simCtx->block_number; ++r) {
        lower->x = PetscMin(lower->x, simCtx->bboxlist[r].min_coords.x);
        lower->y = PetscMin(lower->y, simCtx->bboxlist[r].min_coords.y);
        lower->z = PetscMin(lower->z, simCtx->bboxlist[r].min_coords.z);
        upper->x = PetscMax(upper->x, simCtx->bboxlist[r].max_coords.x);
        upper->y = PetscMax(upper->y, simCtx->bboxlist[r].max_coords.y);
        upper->z = PetscMax(upper->z, simCtx->bboxlist[r].max_coords.z);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/** @brief Position within the domain along one axis, in [0, 1]; 0 across a flat axis. */
static PetscReal Normalized(PetscReal value, PetscReal lower, PetscReal upper)
{
    return (upper > lower) ? (value - lower) / (upper - lower) : 0.0;
}

/**
 * @brief Implementation of \ref ParticleFieldPlanApply().
 * @see ParticleFieldPlanApply()
 */
PetscErrorCode ParticleFieldPlanApply(UserCtx *user, const ParticleFieldPlan *plan,
                                      const ParticleFieldEvent *event, PetscInt first, PetscInt end)
{
    SimCtx           *simCtx = NULL;
    const PetscReal  *positions = NULL;
    const PetscInt64 *pids = NULL;
    Cmpnts            lower, upper;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!plan || end <= first) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    PetscCheck(user && user->swarm && event, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Swarm and event are required.");
    simCtx = user->simCtx;
    PetscCall(DomainBounds(simCtx, &lower, &upper));

    PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void **)&pids));
    for (PetscInt b = 0; b < plan->count; ++b) {
        const ParticleFieldBinding    *binding = &plan->bindings[b];
        const ParticleFieldDescriptor *descriptor = NULL;
        PetscReal                     *field = NULL, scale = 1.0;

        PetscCall(ParticleFieldGetDescriptor(binding->field, &descriptor));
        /* A value is physical; the swarm stores solver units. */
        PetscCall(FieldDimensionReferenceScale(&simCtx->scaling, descriptor->dimension, &scale));
        PetscCall(DMSwarmGetField(user->swarm, descriptor->canonical_name, NULL, NULL, (void **)&field));
        for (PetscInt p = first; p < end; ++p) {
            const PetscReal *position = &positions[3 * p];
            const PetscReal  values[] = {
                position[0] * simCtx->scaling.L_ref, position[1] * simCtx->scaling.L_ref,
                position[2] * simCtx->scaling.L_ref,
                Normalized(position[0], lower.x, upper.x), Normalized(position[1], lower.y, upper.y),
                Normalized(position[2], lower.z, upper.z),
                (PetscReal)pids[p], event->physical_time,
            };
            const PicurvExpressionDrawKey key = {
                (PetscInt64)simCtx->particleRandomSeed, pids[p], event->salt, (PetscInt64)binding->field,
            };
            for (PetscInt c = 0; c < binding->components; ++c) {
                PetscReal value = 0.0;

                PetscCall(PicurvExpressionEvaluate(binding->component[c], values, &key, &value));
                PetscCheck(!PetscIsInfOrNanReal(value), PETSC_COMM_SELF, PETSC_ERR_FP,
                           "The initial value of '%s' is not finite for particle %lld at (%g, %g, %g).",
                           descriptor->canonical_name, (long long)pids[p], (double)values[0],
                           (double)values[1], (double)values[2]);
                field[binding->components * p + c] = value / scale;
            }
        }
        PetscCall(DMSwarmRestoreField(user->swarm, descriptor->canonical_name, NULL, NULL, (void **)&field));
    }
    PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void **)&pids));
    PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&positions));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ParticleFieldPlanSummarize().
 * @see ParticleFieldPlanSummarize()
 */
PetscErrorCode ParticleFieldPlanSummarize(UserCtx *user, const ParticleFieldPlan *plan)
{
    SimCtx *simCtx = NULL;
    FILE   *csv = NULL;
    PetscInt nlocal = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!plan) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    simCtx = user->simCtx;
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    if (simCtx->rank == 0) {
        char path[PETSC_MAX_PATH_LEN];

        PetscCall(PetscSNPrintf(path, sizeof(path), "%s/particle_initial_fields.csv", simCtx->analysis_dir));
        csv = fopen(path, "w");
        if (!csv) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Could not write '%s'.\n", path);
        } else {
            fprintf(csv, "field,component,count,mean,variance,min,max");
            for (int bin = 0; bin < PARTICLE_FIELD_HISTOGRAM_BINS; ++bin) fprintf(csv, ",bin_%d", bin);
            fprintf(csv, "\n");
        }
    }

    for (PetscInt b = 0; b < plan->count; ++b) {
        const ParticleFieldBinding    *binding = &plan->bindings[b];
        const ParticleFieldDescriptor *descriptor = NULL;
        const PetscReal               *field = NULL;
        PetscReal                      scale = 1.0;

        PetscCall(ParticleFieldGetDescriptor(binding->field, &descriptor));
        PetscCall(FieldDimensionReferenceScale(&simCtx->scaling, descriptor->dimension, &scale));
        PetscCall(DMSwarmGetField(user->swarm, descriptor->canonical_name, NULL, NULL, (void **)&field));
        for (PetscInt c = 0; c < binding->components; ++c) {
            /* Sums and extremes in one reduction, then a histogram over the global range. */
            PetscReal sums[3] = {0.0, 0.0, 0.0}, extremes[2] = {PETSC_MAX_REAL, PETSC_MAX_REAL};
            PetscReal bins_local[PARTICLE_FIELD_HISTOGRAM_BINS] = {0.0}, bins[PARTICLE_FIELD_HISTOGRAM_BINS];

            for (PetscInt p = 0; p < nlocal; ++p) {
                const PetscReal value = field[binding->components * p + c] * scale;
                sums[0] += 1.0;
                sums[1] += value;
                sums[2] += value * value;
                extremes[0] = PetscMin(extremes[0], value);
                extremes[1] = PetscMin(extremes[1], -value);
            }
            PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
            PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, extremes, 2, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD));
            const PetscReal count = sums[0], lowest = extremes[0], highest = -extremes[1];
            const PetscReal mean = (count > 0.0) ? sums[1] / count : 0.0;
            const PetscReal variance = (count > 0.0) ? PetscMax(sums[2] / count - mean * mean, 0.0) : 0.0;
            const PetscReal width = (highest > lowest) ? (highest - lowest) / PARTICLE_FIELD_HISTOGRAM_BINS : 1.0;

            for (PetscInt p = 0; p < nlocal; ++p) {
                const PetscReal value = field[binding->components * p + c] * scale;
                PetscInt bin = (PetscInt)((value - lowest) / width);
                bins_local[PetscMax(0, PetscMin(bin, PARTICLE_FIELD_HISTOGRAM_BINS - 1))] += 1.0;
            }
            PetscCallMPI(MPI_Allreduce(bins_local, bins, PARTICLE_FIELD_HISTOGRAM_BINS, MPIU_REAL, MPI_SUM,
                                       PETSC_COMM_WORLD));
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "[Particle IC] %s[%" PetscInt_FMT "]: n=%.0f mean=%.6g var=%.6g min=%.6g max=%.6g\n",
                      descriptor->canonical_name, c, (double)count, (double)mean, (double)variance,
                      (double)(count > 0.0 ? lowest : 0.0), (double)(count > 0.0 ? highest : 0.0));
            if (csv) {
                fprintf(csv, "%s,%d,%.0f,%.10e,%.10e,%.10e,%.10e", descriptor->canonical_name, (int)c,
                        (double)count, (double)mean, (double)variance,
                        (double)(count > 0.0 ? lowest : 0.0), (double)(count > 0.0 ? highest : 0.0));
                for (int bin = 0; bin < PARTICLE_FIELD_HISTOGRAM_BINS; ++bin) fprintf(csv, ",%.0f", (double)bins[bin]);
                fprintf(csv, "\n");
            }
        }
        PetscCall(DMSwarmRestoreField(user->swarm, descriptor->canonical_name, NULL, NULL, (void **)&field));
    }
    if (csv) fclose(csv);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ParticleFieldPlanDestroy().
 * @see ParticleFieldPlanDestroy()
 */
PetscErrorCode ParticleFieldPlanDestroy(ParticleFieldPlan **plan)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!plan || !*plan) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    for (PetscInt b = 0; b < (*plan)->count; ++b) {
        for (PetscInt c = 0; c < PARTICLE_FIELD_MAX_COMPONENTS; ++c) {
            PetscCall(PicurvExpressionDestroy(&(*plan)->bindings[b].component[c]));
        }
    }
    PetscCall(PetscFree((*plan)->bindings));
    PetscCall(PetscFree(*plan));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
