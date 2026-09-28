/**
 * @file test_particle_kernels.c
 * @brief C unit tests for particle walking-search helpers and the initial-value expression engine.
 */

#include "test_support.h"

#include "walkingsearch.h"
#include "ParticleInitialConditions.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**
 * @brief Tests whether a cell index lies within the local grid ownership range.
 */

static PetscErrorCode TestCheckCellWithinLocalGrid(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscBool within = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(CheckCellWithinLocalGrid(user, 1, 1, 1, &within));
    PetscCall(PicurvAssertBool(within, "cell (1,1,1) should be within the serial local grid"));

    PetscCall(CheckCellWithinLocalGrid(user, 4, 1, 1, &within));
    PetscCall(PicurvAssertBool((PetscBool)!within, "cell (4,1,1) should fall outside the valid local cell range"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests initialization of traversal parameters for particle search.
 */

static PetscErrorCode TestInitializeTraversalParameters(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Particle particle;
    PetscInt idx = -1, idy = -1, idz = -1, steps = -1;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMemzero(&particle, sizeof(particle)));

    particle.PID = 42;
    particle.cell[0] = 1;
    particle.cell[1] = 1;
    particle.cell[2] = 1;

    PetscCall(InitializeTraversalParameters(user, &particle, &idx, &idy, &idz, &steps));
    PetscCall(PicurvAssertIntEqual(1, idx, "InitializeTraversalParameters should preserve prior i"));
    PetscCall(PicurvAssertIntEqual(1, idy, "InitializeTraversalParameters should preserve prior j"));
    PetscCall(PicurvAssertIntEqual(1, idz, "InitializeTraversalParameters should preserve prior k"));
    PetscCall(PicurvAssertIntEqual(0, steps, "InitializeTraversalParameters should reset traversal steps"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests retrieval of the current cell from logical particle coordinates.
 */

static PetscErrorCode TestRetrieveCurrentCell(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cell cell;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscMemzero(&cell, sizeof(cell)));

    PetscCall(RetrieveCurrentCell(user, 1, 1, 1, &cell));
    PetscCall(PicurvAssertRealNear(0.25, cell.vertices[0].x, 1.0e-12, "vertex 0 x coordinate"));
    PetscCall(PicurvAssertRealNear(0.25, cell.vertices[0].y, 1.0e-12, "vertex 0 y coordinate"));
    PetscCall(PicurvAssertRealNear(0.25, cell.vertices[0].z, 1.0e-12, "vertex 0 z coordinate"));
    PetscCall(PicurvAssertRealNear(0.50, cell.vertices[5].x, 1.0e-12, "vertex 5 x coordinate"));
    PetscCall(PicurvAssertRealNear(0.50, cell.vertices[5].y, 1.0e-12, "vertex 5 y coordinate"));
    PetscCall(PicurvAssertRealNear(0.50, cell.vertices[5].z, 1.0e-12, "vertex 5 z coordinate"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
static const char *const kConformanceNames[] = {"x", "y", "z"};

/**
 * @brief Tests the runtime evaluator against every shared conformance case.
 * @details `tests/fixtures/expression_conformance.txt` is also evaluated by
 *          `generators/ic.gen` in `tests/test_expression_language.py`, so the two
 *          implementations of the language cannot drift apart unnoticed.
 */
static PetscErrorCode TestExpressionConformance(void)
{
    FILE    *file = fopen("tests/fixtures/expression_conformance.txt", "r");
    char     line[1024];
    PetscInt cases = 0;

    PetscFunctionBeginUser;
    PetscCheck(file, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open the expression conformance cases.");
    while (fgets(line, sizeof(line), file)) {
        char             *fields[5];
        PetscInt          count = 0;
        PicurvExpression *expression = NULL;
        PetscReal         values[3], result = 0.0;

        if (line[0] == '#' || line[0] == '\n') continue;
        for (char *cursor = line; count < 5; ++count) {
            fields[count] = cursor;
            char *separator = strchr(cursor, ';');
            if (!separator) { ++count; break; }
            *separator = '\0';
            cursor = separator + 1;
        }
        PetscCheck(count == 5, PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED, "Malformed conformance line.");
        for (int v = 0; v < 3; ++v) values[v] = strtod(fields[v + 1], NULL);
        PetscCall(PicurvExpressionCompile(fields[0], kConformanceNames, 3, PETSC_FALSE, &expression));
        PetscCall(PicurvExpressionEvaluate(expression, values, NULL, &result));
        PetscCall(PicurvAssertRealNear(strtod(fields[4], NULL), result, 1.0e-12, fields[0]));
        PetscCall(PicurvExpressionDestroy(&expression));
        ++cases;
    }
    fclose(file);
    PetscCall(PicurvAssertBool((PetscBool)(cases >= 20), "the conformance file should hold its cases"));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the compiler refuses what the Python validator also refuses.
 */
static PetscErrorCode TestExpressionRejectsOutsideLanguage(void)
{
    const char *bad[] = {"0x1F", "1_000", "x // 2", "abs", "abs(1, 2)", "where(x, 1)", "foo",
                         "x +", "(x", "x y", "uniform(x)", "uniform(1.5)", "1.e", "not"};

    PetscFunctionBeginUser;
    for (size_t n = 0; n < sizeof(bad) / sizeof(bad[0]); ++n) {
        PicurvExpression *expression = NULL;
        PetscErrorCode    refused;

        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        refused = PicurvExpressionCompile(bad[n], kConformanceNames, 3, PETSC_TRUE, &expression);
        PetscCall(PetscPopErrorHandler());
        PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, refused, bad[n]));
    }
    {
        PicurvExpression *expression = NULL;
        PetscErrorCode    refused;

        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        refused = PicurvExpressionCompile("uniform()", kConformanceNames, 3, PETSC_FALSE, &expression);
        PetscCall(PetscPopErrorHandler());
        PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, refused, "draws need allow_random"));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that random draws are pure functions of their key.
 * @details The same particle always draws the same number; a different particle, event,
 *          or field default stream draws independently; an explicit stream is shared
 *          across fields; and the normal draw has zero mean and unit variance.
 */
static PetscErrorCode TestExpressionDrawsArePureAndKeyed(void)
{
    PicurvExpression       *uniform = NULL, *shared = NULL, *normal = NULL;
    PicurvExpressionDrawKey key = {12345, 7, 0, 6};
    PetscReal               a, b, sum = 0.0, square = 0.0;
    const PetscInt          samples = 20000;

    PetscFunctionBeginUser;
    PetscCall(PicurvExpressionCompile("uniform()", kConformanceNames, 3, PETSC_TRUE, &uniform));
    PetscCall(PicurvExpressionCompile("uniform(3)", kConformanceNames, 3, PETSC_TRUE, &shared));
    PetscCall(PicurvExpressionCompile("normal()", kConformanceNames, 3, PETSC_TRUE, &normal));

    PetscCall(PicurvExpressionEvaluate(uniform, NULL, &key, &a));
    PetscCall(PicurvExpressionEvaluate(uniform, NULL, &key, &b));
    PetscCall(PicurvAssertBool((PetscBool)(a == b && a >= 0.0 && a < 1.0), "a draw is pure and in [0, 1)"));
    key.pid = 8;
    PetscCall(PicurvExpressionEvaluate(uniform, NULL, &key, &b));
    PetscCall(PicurvAssertBool((PetscBool)(a != b), "another particle draws independently"));
    key.pid = 7;
    key.salt = 1;
    PetscCall(PicurvExpressionEvaluate(uniform, NULL, &key, &b));
    PetscCall(PicurvAssertBool((PetscBool)(a != b), "another event draws independently"));
    key.salt = 0;
    key.field = 5;
    PetscCall(PicurvExpressionEvaluate(uniform, NULL, &key, &b));
    PetscCall(PicurvAssertBool((PetscBool)(a != b), "another field's default stream is independent"));
    PetscCall(PicurvExpressionEvaluate(shared, NULL, &key, &a));
    key.field = 6;
    PetscCall(PicurvExpressionEvaluate(shared, NULL, &key, &b));
    PetscCall(PicurvAssertBool((PetscBool)(a == b), "an explicit stream is shared across fields"));

    for (PetscInt pid = 0; pid < samples; ++pid) {
        key.pid = pid;
        PetscCall(PicurvExpressionEvaluate(normal, NULL, &key, &a));
        sum += a;
        square += a * a;
    }
    PetscCall(PicurvAssertRealNear(0.0, sum / samples, 0.03, "normal draws have zero mean"));
    PetscCall(PicurvAssertRealNear(1.0, square / samples, 0.04, "normal draws have unit variance"));

    PetscCall(PicurvExpressionDestroy(&uniform));
    PetscCall(PicurvExpressionDestroy(&shared));
    PetscCall(PicurvExpressionDestroy(&normal));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the unit-particles PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"check-cell-within-local-grid", TestCheckCellWithinLocalGrid},
        {"initialize-traversal-parameters", TestInitializeTraversalParameters},
        {"retrieve-current-cell", TestRetrieveCurrentCell},
        {"expression-conformance", TestExpressionConformance},
        {"expression-rejects-outside-language", TestExpressionRejectsOutsideLanguage},
        {"expression-draws-are-pure-and-keyed", TestExpressionDrawsArePureAndKeyed},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv particle kernel tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-particles", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
