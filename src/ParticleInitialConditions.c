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

/** @brief Names a particle expression may use: those of every source (#PicurvSourceVariable), then `pid`. */
static const char *const kParticleVariables[] = {"x", "y", "z", "t", "xi", "eta", "zeta", "L_xi", "L_eta", "L_zeta", "pid"};
#define PARTICLE_VARIABLE_PID ((PetscInt)PICURV_SOURCE_VAR_CELL_COUNT)
#define PARTICLE_VARIABLE_COUNT (PARTICLE_VARIABLE_PID + 1)
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

/** @brief Whether any configured expression reads a grid line's fraction or length. */
static PetscErrorCode PlanReadsGridLines(const ParticleFieldPlan *plan, PetscBool *reads)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    *reads = PETSC_FALSE;
    for (PetscInt b = 0; b < plan->count; ++b)
    for (PetscInt c = 0; c < plan->bindings[b].components; ++c)
    for (PetscInt var = PICURV_SOURCE_VAR_XI; var <= PICURV_SOURCE_VAR_L_ZETA; ++var) {
        PetscBool uses = PETSC_FALSE;
        PetscCall(PicurvExpressionUsesVariable(plan->bindings[b].component[c], var, &uses));
        if (uses) *reads = PETSC_TRUE;
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ParticleFieldPlanApply().
 * @see ParticleFieldPlanApply()
 */
PetscErrorCode ParticleFieldPlanApply(UserCtx *user, const ParticleFieldPlan *plan,
                                      const ParticleFieldEvent *event, PetscInt first, PetscInt end)
{
    SimCtx           *simCtx = NULL;
    const PetscReal  *positions = NULL, *weights = NULL;
    const PetscInt64 *pids = NULL;
    const PetscInt   *cells = NULL;
    PetscBool         grid_lines = PETSC_FALSE;
    Vec               fraction = NULL, length = NULL, local_fraction = NULL, local_length = NULL;
    const Cmpnts   ***frac = NULL, ***len = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;
    const PetscReal ***aj = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!plan) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }
    PetscCheck(user && event, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context and event are required.");
    simCtx = user->simCtx;
    /* Before the empty-range return: the fractions are collective, and a rank may hold no particle. */
    PetscCall(PlanReadsGridLines(plan, &grid_lines));
    if (grid_lines) {
        PetscCall(DMCreateGlobalVector(user->fda, &fraction));
        PetscCall(DMCreateGlobalVector(user->fda, &length));
        PetscCall(ComputeGridLineFractions(user, fraction, length));
        PetscCall(DMCreateLocalVector(user->fda, &local_fraction));
        PetscCall(DMCreateLocalVector(user->fda, &local_length));
        PetscCall(DMGlobalToLocalBegin(user->fda, fraction, INSERT_VALUES, local_fraction));
        PetscCall(DMGlobalToLocalEnd(user->fda, fraction, INSERT_VALUES, local_fraction));
        PetscCall(DMGlobalToLocalBegin(user->fda, length, INSERT_VALUES, local_length));
        PetscCall(DMGlobalToLocalEnd(user->fda, length, INSERT_VALUES, local_length));
    }
    if (end > first) {
        PetscCheck(user->swarm, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Swarm is required.");
        PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&positions));
        PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void **)&pids));
        if (grid_lines) {
            PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void **)&cells));
            PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_WEIGHT), NULL, NULL, (void **)&weights));
            PetscCall(DMDAVecGetArrayRead(user->fda, local_fraction, &frac));
            PetscCall(DMDAVecGetArrayRead(user->fda, local_length, &len));
            PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, &csi));
            PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, &eta));
            PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, &zet));
            PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
        }
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
                const Cmpnts point = {position[0], position[1], position[2]};
                const PicurvExpressionDrawKey key = {
                    (PetscInt64)simCtx->particleRandomSeed, pids[p], event->salt, (PetscInt64)binding->field,
                };
                PetscReal values[PARTICLE_VARIABLE_COUNT] = {0.0};

                values[PICURV_SOURCE_VAR_X] = position[0] * simCtx->scaling.L_ref;
                values[PICURV_SOURCE_VAR_Y] = position[1] * simCtx->scaling.L_ref;
                values[PICURV_SOURCE_VAR_Z] = position[2] * simCtx->scaling.L_ref;
                if (grid_lines) {
                    /* The host cell in cell-centred numbering; a line's fraction is known at the
                       cell's centre, and moves with the particle by its share of the cell's edge. */
                    const PetscInt i = cells[3 * p] + 1, j = cells[3 * p + 1] + 1, k = cells[3 * p + 2] + 1;
                    const PetscReal lengths[3] = {len[k][j][i].x, len[k][j][i].y, len[k][j][i].z};
                    Cmpnts edges[3];

                    PicurvSourceSetVariables(values, point, frac[k][j][i], len[k][j][i], simCtx->scaling.L_ref);
                    PetscCall(ComputeCellEdgeVectors(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], edges));
                    for (PetscInt a = 0; a < 3; ++a) {
                        const PetscReal ds = PetscSqrtReal(edges[a].x * edges[a].x + edges[a].y * edges[a].y + edges[a].z * edges[a].z);
                        values[PICURV_SOURCE_VAR_XI + a] += (weights[3 * p + a] - 0.5) * ds / lengths[a];
                    }
                }
                values[PICURV_SOURCE_VAR_T] = event->physical_time;
                values[PARTICLE_VARIABLE_PID] = (PetscReal)pids[p];
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
        if (grid_lines) {
            PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));
            PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet));
            PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta));
            PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi));
            PetscCall(DMDAVecRestoreArrayRead(user->fda, local_length, &len));
            PetscCall(DMDAVecRestoreArrayRead(user->fda, local_fraction, &frac));
            PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_WEIGHT), NULL, NULL, (void **)&weights));
            PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void **)&cells));
        }
        PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_PID), NULL, NULL, (void **)&pids));
        PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void **)&positions));
    }
    PetscCall(VecDestroy(&local_length));
    PetscCall(VecDestroy(&local_fraction));
    PetscCall(VecDestroy(&length));
    PetscCall(VecDestroy(&fraction));
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
