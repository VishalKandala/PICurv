/**
 * @file test_particle_kernels.c
 * @brief C unit tests for particle walking-search helpers and the initial-value expression engine.
 */

#include "test_support.h"
#include "Metric.h"

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
/**
 * @brief Fills a cell with the unit cube's vertices in walking-search order, x sheared by `shear * z`.
 */
static void BuildShearedUnitCell(PetscReal shear, Cell *cell)
{
    const PetscReal corner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                    {0, 1, 1}, {1, 1, 1}, {1, 0, 1}, {0, 0, 1}};
    for (int v = 0; v < 8; ++v) {
        cell->vertices[v].x = corner[v][0] + shear * corner[v][2];
        cell->vertices[v].y = corner[v][1];
        cell->vertices[v].z = corner[v][2];
    }
}

/**
 * @brief Tests the face mirror on a flat face, a sheared face, and a two-face corner.
 */
static PetscErrorCode TestReflectParticleAcrossCellFace(void)
{
    Cell           cell;
    CellFacePlanes planes;
    Cmpnts         p;

    PetscFunctionBeginUser;
    BuildShearedUnitCell(0.0, &cell);
    MetricBuildCellFacePlanes(&cell, &planes);
    p.x = -0.2; p.y = 0.5; p.z = 0.5;
    PetscCall(ReflectParticleAcrossCellFace(&planes, LEFT, &p));
    PetscCall(PicurvAssertRealNear(0.2, p.x, 1.0e-14, "flat face: x mirrors about x = 0"));
    PetscCall(PicurvAssertRealNear(0.5, p.y, 1.0e-14, "flat face: y is unchanged"));
    PetscCall(PicurvAssertRealNear(0.5, p.z, 1.0e-14, "flat face: z is unchanged"));

    /* LEFT face of the sheared cell is the plane x = 0.5 z; (-0.3, 0.5, 0.4) mirrors to (0.5, 0.5, 0). */
    BuildShearedUnitCell(0.5, &cell);
    MetricBuildCellFacePlanes(&cell, &planes);
    p.x = -0.3; p.y = 0.5; p.z = 0.4;
    PetscCall(ReflectParticleAcrossCellFace(&planes, LEFT, &p));
    PetscCall(PicurvAssertRealNear(0.5, p.x, 1.0e-14, "sheared face: mirrored x"));
    PetscCall(PicurvAssertRealNear(0.5, p.y, 1.0e-14, "sheared face: y is unchanged"));
    PetscCall(PicurvAssertRealNear(0.0, p.z, 1.0e-14, "sheared face: mirrored z"));

    /* Corner: reflecting about LEFT and then BOTTOM mirrors both coordinates. */
    BuildShearedUnitCell(0.0, &cell);
    MetricBuildCellFacePlanes(&cell, &planes);
    p.x = -0.1; p.y = -0.2; p.z = 0.5;
    PetscCall(ReflectParticleAcrossCellFace(&planes, LEFT, &p));
    PetscCall(ReflectParticleAcrossCellFace(&planes, BOTTOM, &p));
    PetscCall(PicurvAssertRealNear(0.1, p.x, 1.0e-14, "corner: x mirrored"));
    PetscCall(PicurvAssertRealNear(0.2, p.y, 1.0e-14, "corner: y mirrored"));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the reflection-cap fallback places a far-outside point just inside the cell.
 */
static PetscErrorCode TestPlaceParticleInsideCell(void)
{
    Cell           cell;
    CellFacePlanes planes;
    Cmpnts         p;
    PetscBool      placed = PETSC_FALSE;
    PetscReal      d[NUM_FACES];
    PetscInt       position = -1;

    PetscFunctionBeginUser;
    BuildShearedUnitCell(0.0, &cell);
    p.x = -1.0; p.y = 0.5; p.z = 0.5;
    PetscCall(PlaceParticleInsideCell(&cell, &p, &placed));
    PetscCall(PicurvAssertBool(placed, "a point outside one face of a box cell can be placed"));
    PetscCall(PicurvAssertBool((PetscBool)(p.x > 0.0 && p.x < 1.0e-5), "placed just inside the crossed face"));
    MetricBuildCellFacePlanes(&cell, &planes);
    PetscCall(GetCellCharacteristicSize(&cell, &planes.size));
    PetscCall(EvaluateParticlePositionInPlanes(&planes, d, p, &position, 1.0e-11));
    PetscCall(PicurvAssertBool((PetscBool)(position >= 0), "the placed point evaluates inside the cell"));
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
 * @brief Tests the Bessel and Kelvin functions against 30-digit references from their series,
 *        and that a non-integer Bessel order gives NaN.
 */
static PetscErrorCode TestExpressionSpecialFunctions(void)
{
    const struct { const char *text; PetscReal expected; } cases[] = {
        {"besselj(0, 1)", 0.76519768655796661},  {"besselj(1, x)", 0.4400505857449335},
        {"besselj(2, 2.5)", 0.44605905843961724}, {"besselj(0, 10)", -0.24593576445134835},
        {"bessely(0, 1)", 0.088256964215676956},  {"bessely(0, 2.5)", 0.49807035961523188},
        {"ber(1)", 0.98438178121308684},          {"bei(x)", 0.24956604003665972},
        {"ber(5)", -6.2300824786663576},          {"bei(5)", 0.11603438155020038},
        {"ber(20)", 47489.370265061763},          {"bei(20)", 114775.19736006622},
        /* Either side of the change to the asymptotic expansion at 30, and well beyond the series' reach. */
        {"ber(29.5)", -3674268.2459448087},       {"bei(29.5)", 84356546.89143799},
        {"ber(30.5)", -114878788.18171336},       {"bei(30.5)", 123131394.62887576},
        {"ber(50)", -117623968512357.44},         {"bei(50)", -50192646254462.21},
        {"ber(120)", -2.417573872141002e+35},     {"bei(120)", 9.200088423773959e+34},
        {"ber(-40)", -112596696872.07408},        {"bei(-40)", 45628130285.64126},
    };
    const PetscReal values[3] = {1.0, 0.0, 0.0};

    PetscFunctionBeginUser;
    for (size_t n = 0; n < sizeof(cases) / sizeof(cases[0]); ++n) {
        PicurvExpression *expression = NULL;
        PetscReal         result = 0.0;
        PetscCall(PicurvExpressionCompile(cases[n].text, kConformanceNames, 3, PETSC_FALSE, &expression));
        PetscCall(PicurvExpressionEvaluate(expression, values, NULL, &result));
        PetscCall(PicurvAssertRealNear(cases[n].expected, result, 1.0e-12 * PetscMax(1.0, PetscAbsReal(cases[n].expected)),
                                       cases[n].text));
        PetscCall(PicurvExpressionDestroy(&expression));
    }
    {
        PicurvExpression *expression = NULL;
        PetscReal         result = 0.0;
        PetscCall(PicurvExpressionCompile("besselj(0.5, 1)", kConformanceNames, 3, PETSC_FALSE, &expression));
        PetscCall(PicurvExpressionEvaluate(expression, values, NULL, &result));
        PetscCall(PicurvAssertBool((PetscBool)PetscIsNanReal(result), "a non-integer order is not a value"));
        PetscCall(PicurvExpressionDestroy(&expression));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Tests the documented closed-form sources. The Womersley source, written in real
 *        arithmetic with `ber` and `bei`, against the complex solution `Re{(K/i w)(1 - J0(i^1.5 a s)/J0(i^1.5 a)) e^{iwt}}`
 *        evaluated independently from its series (K = 1, w = 2 pi, a = 5, R = 0.5); the
 *        rectangular-duct series against the Python evaluator.
 */
static PetscErrorCode TestExpressionClosedFormProfiles(void)
{
    static const char *const names[] = {"x", "y", "z", "t"};
    const char *text =
        "(1.0/(2*pi))*((ber(5.0*(sqrt(x**2 + y**2)/0.5))*bei(5.0) - bei(5.0*(sqrt(x**2 + y**2)/0.5))*ber(5.0))"
        "/(ber(5.0)**2 + bei(5.0)**2)*cos(2*pi*t) + (1 - (ber(5.0*(sqrt(x**2 + y**2)/0.5))*ber(5.0)"
        " + bei(5.0*(sqrt(x**2 + y**2)/0.5))*bei(5.0))/(ber(5.0)**2 + bei(5.0)**2))*sin(2*pi*t))";
    const PetscReal points[4][5] = {{0.0, 0.0, 0.0, 0.1, 0.10894419422254727},
                                    {0.2, 0.1, 0.0, 0.3, 0.15619922873565129},
                                    {0.3, -0.2, 0.0, 0.7, -0.1332246948462772},
                                    {0.0, 0.45, 0.0, 0.9, 0.004632857508023867}};
    PicurvExpression *expression = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvExpressionCompile(text, names, 4, PETSC_FALSE, &expression));
    for (PetscInt n = 0; n < 4; ++n) {
        PetscReal result = 0.0;
        PetscCall(PicurvExpressionEvaluate(expression, points[n], NULL, &result));
        PetscCall(PicurvAssertRealNear(points[n][4], result, 1.0e-12, "Womersley velocity"));
    }
    PetscCall(PicurvExpressionDestroy(&expression));
    {
        /* The documented rectangular-duct series (unit bulk, square section, xi -> x, eta -> y),
           against the Python evaluator's value at (0.3, 0.6). */
        const char *duct =
        "(1.0)*sum(k, 0, 99, (-1)**k/(2*k+1)**3 *(1 - (exp(abs((2*k+1)*pi*(1.0)*(y - 0.5)) - (2*k+1"
        ")*pi*(1.0)/2) + exp(-abs((2*k+1)*pi*(1.0)*(y - 0.5)) - (2*k+1)*pi*(1.0)/2))/(1 + exp(-(2*k"
        "+1)*pi*(1.0)))) *cos((2*k+1)*pi*(x - 0.5))) / sum(k, 0, 99, 2/(pi*(2*k+1)**4)*(1 - 2*tanh("
        "(2*k+1)*pi*(1.0)/2)/((2*k+1)*pi*(1.0))))";
        const PetscReal at[4] = {0.3, 0.6, 0.0, 0.0};
        PetscReal       result = 0.0;
        PetscCall(PicurvExpressionCompile(duct, names, 4, PETSC_FALSE, &expression));
        PetscCall(PicurvExpressionEvaluate(expression, at, NULL, &result));
        PetscCall(PicurvAssertRealNear(1.7442020681945691, result, 1.0e-12, "duct Poiseuille series"));
        PetscCall(PicurvExpressionDestroy(&expression));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the compiler refuses what the Python validator also refuses.
 */
static PetscErrorCode TestExpressionRejectsOutsideLanguage(void)
{
    const char *bad[] = {"0x1F", "1_000", "x // 2", "abs", "abs(1, 2)", "where(x, 1)", "foo",
                         "x +", "(x", "x y", "uniform(x)", "uniform(1.5)", "1.e", "not",
                         "sum(n, 0, 1.5, n)", "sum(n, 0, y, n)", "sum(x, 0, 2, x)", "sum(n, 0, 2)",
                         "sum(n, 0, 2, sum(n, 0, 1, n))", "sum(n, 0, 2, sum(m, 0, n, m))", "sum(n, 0, 2, n) + n", "sum(n, 0, 2000000, n)",
                         "sum(pi, 0, 2, 1)", "besselj(1)", "ber(1, 2)"};

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
    key.field = 6;
    key.grid = PETSC_TRUE;
    PetscCall(PicurvExpressionEvaluate(uniform, NULL, &key, &b));
    PetscCall(PicurvAssertBool((PetscBool)(a != b), "a grid point with the same number as a particle draws independently"));
    key.grid = PETSC_FALSE;
    key.field = 5;
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
        {"reflect-particle-across-cell-face", TestReflectParticleAcrossCellFace},
        {"place-particle-inside-cell", TestPlaceParticleInsideCell},
        {"expression-conformance", TestExpressionConformance},
        {"expression-rejects-outside-language", TestExpressionRejectsOutsideLanguage},
        {"expression-special-functions", TestExpressionSpecialFunctions},
        {"expression-closed-form-profiles", TestExpressionClosedFormProfiles},
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
