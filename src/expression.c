/**
 * @file expression.c
 * @brief The expression language shared by every configured field: a parser, a postfix
 *        compiler, and a pure stack evaluator.
 *
 * Expressions are parsed into a small tree, because a chained comparison evaluates its
 * middle operand twice, and compiled to postfix code that a stack evaluator runs once per
 * point. A series `sum(k, lo, hi, body)` compiles its body once and loops over it, with the
 * index in a slot of its nesting level. Evaluation is pure: random draws are hashed from a caller-supplied key, so a value
 * does not depend on rank count, evaluation order, or how often it is evaluated.
 */

#include "expression.h"
#include "Metric.h"

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
/* ------------------------------------------------------------------------- */
/*                               Expression tree                              */
/* ------------------------------------------------------------------------- */

typedef enum {
    OP_CONST = 0, OP_VAR, OP_NEG, OP_POS, OP_NOT,
    OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_MOD, OP_POW,
    OP_EQ, OP_NE, OP_LT, OP_LE, OP_GT, OP_GE, OP_AND, OP_OR,
    OP_ABS, OP_COS, OP_ERF, OP_EXP, OP_SIN, OP_SQRT, OP_TAN, OP_TANH, OP_COSH, OP_SINH, OP_LOG, OP_FLOOR,
    OP_ATAN2, OP_MAXIMUM, OP_MINIMUM, OP_WHERE,
    OP_UNIFORM, OP_NORMAL,
    OP_BESSELJ, OP_BESSELY, OP_BER, OP_BEI,
    /* A series: START pushes the running total and sets the index, HI carries the upper bound
       and the body's length, ACC adds the body's value and loops; LOCAL reads an index. */
    OP_SUM, OP_SUM_START, OP_SUM_HI, OP_SUM_ACC, OP_LOCAL
} ExpressionOp;

/** @brief Name, arity, and opcode of every callable function. */
static const struct {
    const char  *name;
    int          min_args;
    int          max_args;
    ExpressionOp op;
    PetscBool    random;
} kFunctions[] = {
    {"abs", 1, 1, OP_ABS, PETSC_FALSE},         {"cos", 1, 1, OP_COS, PETSC_FALSE},
    {"exp", 1, 1, OP_EXP, PETSC_FALSE},         {"sin", 1, 1, OP_SIN, PETSC_FALSE},
    {"sqrt", 1, 1, OP_SQRT, PETSC_FALSE},       {"tan", 1, 1, OP_TAN, PETSC_FALSE},
    {"erf", 1, 1, OP_ERF, PETSC_FALSE},         {"tanh", 1, 1, OP_TANH, PETSC_FALSE},
    {"cosh", 1, 1, OP_COSH, PETSC_FALSE},       {"sinh", 1, 1, OP_SINH, PETSC_FALSE},
    {"log", 1, 1, OP_LOG, PETSC_FALSE},         {"floor", 1, 1, OP_FLOOR, PETSC_FALSE},
    {"atan2", 2, 2, OP_ATAN2, PETSC_FALSE},
    {"maximum", 2, 2, OP_MAXIMUM, PETSC_FALSE}, {"minimum", 2, 2, OP_MINIMUM, PETSC_FALSE},
    {"where", 3, 3, OP_WHERE, PETSC_FALSE},
    {"besselj", 2, 2, OP_BESSELJ, PETSC_FALSE}, {"bessely", 2, 2, OP_BESSELY, PETSC_FALSE},
    {"ber", 1, 1, OP_BER, PETSC_FALSE},         {"bei", 1, 1, OP_BEI, PETSC_FALSE},
    {"uniform", 0, 1, OP_UNIFORM, PETSC_TRUE},  {"normal", 0, 1, OP_NORMAL, PETSC_TRUE},
};

/** @brief Series that may nest inside one another, and terms one series may add. */
#define EXPRESSION_MAX_SERIES_DEPTH 4
#define EXPRESSION_MAX_SERIES_TERMS 1000000

/** @brief One tree node; children index the node pool. */
typedef struct {
    ExpressionOp op;
    PetscReal    value;     /* OP_CONST value; explicit stream of a draw; a series' lower bound */
    PetscReal    upper;     /* A series' upper bound */
    PetscInt     var;       /* OP_VAR index; 1 for a draw on an explicit stream; an index slot */
    PetscInt     child[3];
    PetscInt     nchild;
} ExpressionNode;

/** @brief One postfix instruction. */
typedef struct {
    ExpressionOp op;
    PetscReal    value;
    PetscInt     arg;
} ExpressionInstruction;

struct PicurvExpression {
    ExpressionInstruction *code;
    PetscInt               length;
    PetscInt               max_depth;
    PetscBool              draws;
};

/** @brief Parser state: the text, a cursor, the node pool, and the names in scope. */
typedef struct {
    const char        *text;
    size_t             pos;
    ExpressionNode    *nodes;
    PetscInt           count;
    PetscInt           capacity;
    const char *const *names;
    PetscInt           name_count;
    PetscBool          allow_random;
    char               series_index[EXPRESSION_MAX_SERIES_DEPTH][64]; /* Indices in scope, outermost first */
    PetscInt           series_depth;
    char               error[256];
} ExpressionParser;

/** @brief Record the first parse error with its position; later ones are consequences. */
static void ParserFail(ExpressionParser *parser, const char *message)
{
    PROFILE_HELPER_BEGIN;
    if (!parser->error[0]) {
        snprintf(parser->error, sizeof(parser->error), "%s at position %zu of '%s'",
                 message, parser->pos + 1, parser->text);
    }
    PROFILE_HELPER_END;
}

/** @brief Skip whitespace before the next token. */
static void SkipSpace(ExpressionParser *parser)
{
    PROFILE_FUNCTION_BEGIN;
    while (parser->text[parser->pos] && isspace((unsigned char)parser->text[parser->pos])) parser->pos++;
    PROFILE_FUNCTION_END;
}

/** @brief Consume @p token if it is next, and report whether it was. */
static PetscBool Accept(ExpressionParser *parser, const char *token)
{
    PROFILE_FUNCTION_BEGIN;
    size_t length = strlen(token);

    SkipSpace(parser);
    if (strncmp(parser->text + parser->pos, token, length) != 0) { PROFILE_FUNCTION_END; return PETSC_FALSE; }
    if (isalpha((unsigned char)token[0])) {
        const char next = parser->text[parser->pos + length];
        if (isalnum((unsigned char)next) || next == '_') { PROFILE_FUNCTION_END; return PETSC_FALSE; }
    }
    /* `*` must not match the first half of `**`, nor `<` the first half of `<=`. */
    if ((!strcmp(token, "*") && parser->text[parser->pos + 1] == '*') ||
        ((!strcmp(token, "<") || !strcmp(token, ">")) && parser->text[parser->pos + 1] == '=')) {
        PROFILE_FUNCTION_END;
        return PETSC_FALSE;
    }
    parser->pos += length;
    PROFILE_FUNCTION_END;
    return PETSC_TRUE;
}

/** @brief Append a node to the pool and return its index, or -1 once an error is set. */
static PetscInt NewNode(ExpressionParser *parser, ExpressionOp op, PetscInt nchild,
                        PetscInt a, PetscInt b, PetscInt c)
{
    PROFILE_FUNCTION_BEGIN;
    if (parser->error[0]) { PROFILE_FUNCTION_END; return -1; }
    if (parser->count == parser->capacity) {
        ParserFail(parser, "expression is too long");
        PROFILE_FUNCTION_END;
        return -1;
    }
    ExpressionNode *node = &parser->nodes[parser->count];
    memset(node, 0, sizeof(*node));
    node->op = op;
    node->nchild = nchild;
    node->child[0] = a;
    node->child[1] = b;
    node->child[2] = c;
    PROFILE_FUNCTION_END;
    return parser->count++;
}

static PetscInt ParseOr(ExpressionParser *parser);

/** @brief The value of a series bound: an integer literal, possibly negated or parenthesized. */
static PetscBool SeriesBound(const ExpressionParser *parser, PetscInt index, PetscReal *bound)
{
    PROFILE_HELPER_BEGIN;
    const ExpressionNode *node = &parser->nodes[index];
    PetscReal sign = 1.0;
    if (node->op == OP_NEG || node->op == OP_POS) {
        if (node->op == OP_NEG) sign = -1.0;
        node = &parser->nodes[node->child[0]];
    }
    *bound = sign * node->value;
    PROFILE_HELPER_END;
    return (PetscBool)(node->op == OP_CONST && node->value == PetscFloorReal(node->value));
}

/** @brief `sum(` has been read: index name, integer bounds, then the body with the index in scope. */
static PetscInt ParseSeries(ExpressionParser *parser)
{
    PROFILE_HELPER_BEGIN;
    char      name[64];
    size_t    length = 0;
    PetscReal lower = 0.0, upper = 0.0;

    SkipSpace(parser);
    const char *start = parser->text + parser->pos;
    while (isalnum((unsigned char)start[length]) || start[length] == '_') length++;
    if (!length || isdigit((unsigned char)start[0]) || length >= sizeof(name)) {
        ParserFail(parser, "a series needs an index name first");
        PROFILE_HELPER_END;
        return -1;
    }
    memcpy(name, start, length);
    name[length] = '\0';
    parser->pos += length;
    for (PetscInt v = 0; v < parser->name_count; ++v) {
        if (!strcmp(name, parser->names[v])) { ParserFail(parser, "a series index cannot reuse a variable's name"); PROFILE_HELPER_END; return -1; }
    }
    for (PetscInt level = 0; level < parser->series_depth; ++level) {
        if (!strcmp(name, parser->series_index[level])) { ParserFail(parser, "a series index cannot reuse an enclosing index"); PROFILE_HELPER_END; return -1; }
    }
    if (!strcmp(name, "pi")) { ParserFail(parser, "a series index cannot be pi"); PROFILE_HELPER_END; return -1; }
    if (parser->series_depth == EXPRESSION_MAX_SERIES_DEPTH) { ParserFail(parser, "series nest too deeply"); PROFILE_HELPER_END; return -1; }
    if (!Accept(parser, ",")) { ParserFail(parser, "expected ',' after the series index"); PROFILE_HELPER_END; return -1; }
    const PetscInt lo = ParseOr(parser);
    if (parser->error[0] || !Accept(parser, ",")) { ParserFail(parser, "expected ',' after the lower bound"); PROFILE_HELPER_END; return -1; }
    const PetscInt hi = ParseOr(parser);
    if (parser->error[0] || !Accept(parser, ",")) { ParserFail(parser, "expected ',' after the upper bound"); PROFILE_HELPER_END; return -1; }
    if (!SeriesBound(parser, lo, &lower) || !SeriesBound(parser, hi, &upper)) {
        ParserFail(parser, "series bounds must be integer numbers");
        PROFILE_HELPER_END;
        return -1;
    }
    if (upper - lower + 1.0 > (PetscReal)EXPRESSION_MAX_SERIES_TERMS) {
        ParserFail(parser, "a series may add at most 1000000 terms");
        PROFILE_HELPER_END;
        return -1;
    }
    memcpy(parser->series_index[parser->series_depth], name, length + 1);
    parser->series_depth++;
    const PetscInt body = ParseOr(parser);
    parser->series_depth--;
    if (parser->error[0]) { PROFILE_HELPER_END; return -1; }
    if (!Accept(parser, ")")) { ParserFail(parser, "expected ')' after the series body"); PROFILE_HELPER_END; return -1; }
    PetscInt node = NewNode(parser, OP_SUM, 1, body, -1, -1);
    if (node >= 0) {
        parser->nodes[node].value = lower;
        parser->nodes[node].upper = upper;
        parser->nodes[node].var = parser->series_depth;
    }
    PROFILE_HELPER_END;
    return node;
}

/** @brief number | name | name(args) | sum(index, lo, hi, body) | (expression) */
static PetscInt ParseAtom(ExpressionParser *parser)
{
    PROFILE_HELPER_BEGIN;
    SkipSpace(parser);
    const char *start = parser->text + parser->pos;

    if (isdigit((unsigned char)*start) || (*start == '.' && isdigit((unsigned char)start[1]))) {
        /* Decimal literals only, as the Python validator requires. */
        size_t length = 0;
        while (isdigit((unsigned char)start[length])) length++;
        if (start[length] == '.') { length++; while (isdigit((unsigned char)start[length])) length++; }
        if (start[length] == 'e' || start[length] == 'E') {
            size_t exponent = length + 1;
            if (start[exponent] == '+' || start[exponent] == '-') exponent++;
            if (!isdigit((unsigned char)start[exponent])) { ParserFail(parser, "malformed number"); PROFILE_HELPER_END; return -1; }
            while (isdigit((unsigned char)start[exponent])) exponent++;
            length = exponent;
        }
        if (isalnum((unsigned char)start[length]) || start[length] == '_' || start[length] == '.') {
            ParserFail(parser, "malformed number");
            PROFILE_HELPER_END;
            return -1;
        }
        const PetscReal value = strtod(start, NULL);
        parser->pos += length;
        PetscInt node = NewNode(parser, OP_CONST, 0, -1, -1, -1);
        if (node >= 0) parser->nodes[node].value = value;
        PROFILE_HELPER_END;
        return node;
    }
    if (isalpha((unsigned char)*start) || *start == '_') {
        char name[64];
        size_t length = 0;
        while (isalnum((unsigned char)start[length]) || start[length] == '_') length++;
        if (length >= sizeof(name)) { ParserFail(parser, "name is too long"); PROFILE_HELPER_END; return -1; }
        memcpy(name, start, length);
        name[length] = '\0';
        if (!strcmp(name, "and") || !strcmp(name, "or") || !strcmp(name, "not")) {
            ParserFail(parser, "operator where an operand was expected");
            PROFILE_HELPER_END;
            return -1;
        }
        parser->pos += length;
        if (!strcmp(name, "sum") && Accept(parser, "(")) { PROFILE_HELPER_END; return ParseSeries(parser); }
        if (Accept(parser, "(")) {
            for (size_t f = 0; f < sizeof(kFunctions) / sizeof(kFunctions[0]); ++f) {
                if (strcmp(name, kFunctions[f].name)) continue;
                if (kFunctions[f].random && !parser->allow_random) {
                    ParserFail(parser, "random draws are not available in this expression");
                    PROFILE_HELPER_END;
                    return -1;
                }
                PetscInt args[3] = {-1, -1, -1}, nargs = 0;
                if (!Accept(parser, ")")) {
                    do {
                        if (nargs == 3) { ParserFail(parser, "too many arguments"); PROFILE_HELPER_END; return -1; }
                        args[nargs++] = ParseOr(parser);
                        if (parser->error[0]) { PROFILE_HELPER_END; return -1; }
                    } while (Accept(parser, ","));
                    if (!Accept(parser, ")")) { ParserFail(parser, "expected ')'"); PROFILE_HELPER_END; return -1; }
                }
                if (nargs < kFunctions[f].min_args || nargs > kFunctions[f].max_args) {
                    ParserFail(parser, "wrong number of arguments");
                    PROFILE_HELPER_END;
                    return -1;
                }
                if (kFunctions[f].random) {
                    /* A draw's stream is a constant so that two expressions naming the
                       same stream draw the same number for the same particle. */
                    PetscInt node = NewNode(parser, kFunctions[f].op, 0, -1, -1, -1);
                    if (node < 0) { PROFILE_HELPER_END; return -1; }
                    if (nargs == 1) {
                        const ExpressionNode *stream = &parser->nodes[args[0]];
                        if (stream->op != OP_CONST || stream->value < 0.0 ||
                            stream->value != PetscFloorReal(stream->value)) {
                            ParserFail(parser, "a draw's stream must be a non-negative integer literal");
                            PROFILE_HELPER_END;
                            return -1;
                        }
                        parser->nodes[node].value = stream->value;
                        parser->nodes[node].var = 1;
                    }
                    PROFILE_HELPER_END;
                    return node;
                }
                PROFILE_HELPER_END;
                return NewNode(parser, kFunctions[f].op, nargs, args[0], args[1], args[2]);
            }
            ParserFail(parser, "unknown function");
            PROFILE_HELPER_END;
            return -1;
        }
        for (PetscInt level = parser->series_depth - 1; level >= 0; --level) {
            if (strcmp(name, parser->series_index[level])) continue;
            PetscInt node = NewNode(parser, OP_LOCAL, 0, -1, -1, -1);
            if (node >= 0) parser->nodes[node].var = level;
            PROFILE_HELPER_END;
            return node;
        }
        if (!strcmp(name, "pi")) {
            PetscInt node = NewNode(parser, OP_CONST, 0, -1, -1, -1);
            if (node >= 0) parser->nodes[node].value = PETSC_PI;
            PROFILE_HELPER_END;
            return node;
        }
        for (PetscInt v = 0; v < parser->name_count; ++v) {
            if (!strcmp(name, parser->names[v])) {
                PetscInt node = NewNode(parser, OP_VAR, 0, -1, -1, -1);
                if (node >= 0) parser->nodes[node].var = v;
                PROFILE_HELPER_END;
                return node;
            }
        }
        ParserFail(parser, "unknown name");
        PROFILE_HELPER_END;
        return -1;
    }
    if (Accept(parser, "(")) {
        PetscInt inner = ParseOr(parser);
        if (!Accept(parser, ")")) { ParserFail(parser, "expected ')'"); PROFILE_HELPER_END; return -1; }
        PROFILE_HELPER_END;
        return inner;
    }
    ParserFail(parser, "expected a number, name, or '('");
    PROFILE_HELPER_END;
    return -1;
}

static PetscInt ParseUnary(ExpressionParser *parser);

/**
 * @brief Parses an atom and, when a power operator follows, appends a power node.
 * @details Grammar: `atom ['**' unary]`.
 *          The exponent is parsed as a unary, so `2**-1` is valid and `a**b**c` groups to
 *          the right; a unary minus on the left is parsed by the caller, so `-2**2` is -4.
 * @param[in,out] parser Parser state; its position advances past the power.
 * @return Node index of the power (or of the bare atom), or -1 after an error.
 */
static PetscInt ParsePower(ExpressionParser *parser)
{
    PROFILE_HELPER_BEGIN;
    PetscInt base = ParseAtom(parser);
    if (Accept(parser, "**")) { PROFILE_HELPER_END; return NewNode(parser, OP_POW, 2, base, ParseUnary(parser), -1); }
    PROFILE_HELPER_END;
    return base;
}

/**
 * @brief Parses a sign prefix recursively and wraps the operand in a negate or identity node.
 * @param[in,out] parser Parser state; its position advances past the operand.
 * @return Node index of the signed operand, or -1 after an error.
 */
static PetscInt ParseUnary(ExpressionParser *parser)
{
    PROFILE_FUNCTION_BEGIN;
    if (Accept(parser, "-")) { PROFILE_FUNCTION_END; return NewNode(parser, OP_NEG, 1, ParseUnary(parser), -1, -1); }
    if (Accept(parser, "+")) { PROFILE_FUNCTION_END; return NewNode(parser, OP_POS, 1, ParseUnary(parser), -1, -1); }
    PROFILE_FUNCTION_END;
    return ParsePower(parser);
}

/**
 * @brief Parses a left-associative chain of `*`, `/`, `%`, folding each operator into a node.
 * @details `//` is refused explicitly, since the Python side does not accept floor division.
 * @param[in,out] parser Parser state; its position advances past the term.
 * @return Node index of the term, or -1 after an error.
 */
static PetscInt ParseTerm(ExpressionParser *parser)
{
    PROFILE_FUNCTION_BEGIN;
    PetscInt left = ParseUnary(parser);
    for (;;) {
        if (Accept(parser, "*")) left = NewNode(parser, OP_MUL, 2, left, ParseUnary(parser), -1);
        else if (Accept(parser, "/")) {
            if (parser->text[parser->pos] == '/') { ParserFail(parser, "'//' is not supported"); PROFILE_FUNCTION_END; return -1; }
            left = NewNode(parser, OP_DIV, 2, left, ParseUnary(parser), -1);
        } else if (Accept(parser, "%")) left = NewNode(parser, OP_MOD, 2, left, ParseUnary(parser), -1);
        else { PROFILE_FUNCTION_END; return left; }
    }
    PROFILE_FUNCTION_END;
}

/**
 * @brief Parses a left-associative chain of `+` and `-`, folding each operator into a node.
 * @param[in,out] parser Parser state; its position advances past the sum.
 * @return Node index of the sum, or -1 after an error.
 */
static PetscInt ParseSum(ExpressionParser *parser)
{
    PROFILE_FUNCTION_BEGIN;
    PetscInt left = ParseTerm(parser);
    for (;;) {
        if (Accept(parser, "+")) left = NewNode(parser, OP_ADD, 2, left, ParseTerm(parser), -1);
        else if (Accept(parser, "-")) left = NewNode(parser, OP_SUB, 2, left, ParseTerm(parser), -1);
        else { PROFILE_FUNCTION_END; return left; }
    }
    PROFILE_FUNCTION_END;
}

/** @brief sum (comparison sum)*; a chain `a < b < c` is `(a < b) and (b < c)`. */
static PetscInt ParseComparison(ExpressionParser *parser)
{
    PROFILE_FUNCTION_BEGIN;
    static const struct { const char *token; ExpressionOp op; } kComparisons[] = {
        {"==", OP_EQ}, {"!=", OP_NE}, {"<=", OP_LE}, {">=", OP_GE}, {"<", OP_LT}, {">", OP_GT},
    };
    PetscInt left = ParseSum(parser), result = -1;

    for (;;) {
        ExpressionOp op = OP_CONST;
        for (size_t c = 0; c < sizeof(kComparisons) / sizeof(kComparisons[0]); ++c) {
            if (Accept(parser, kComparisons[c].token)) { op = kComparisons[c].op; break; }
        }
        if (op == OP_CONST) { PROFILE_FUNCTION_END; return (result < 0) ? left : result; }
        const PetscInt right = ParseSum(parser);
        const PetscInt pair = NewNode(parser, op, 2, left, right, -1);
        result = (result < 0) ? pair : NewNode(parser, OP_AND, 2, result, pair, -1);
        left = right;
    }
    PROFILE_FUNCTION_END;
}

/** @brief 'not' not | comparison */
static PetscInt ParseNot(ExpressionParser *parser)
{
    PROFILE_FUNCTION_BEGIN;
    if (Accept(parser, "not")) { PROFILE_FUNCTION_END; return NewNode(parser, OP_NOT, 1, ParseNot(parser), -1, -1); }
    PROFILE_FUNCTION_END;
    return ParseComparison(parser);
}

/** @brief not ('and' not)* */
static PetscInt ParseAnd(ExpressionParser *parser)
{
    PROFILE_HELPER_BEGIN;
    PetscInt left = ParseNot(parser);
    while (Accept(parser, "and")) left = NewNode(parser, OP_AND, 2, left, ParseNot(parser), -1);
    PROFILE_HELPER_END;
    return left;
}

/** @brief and ('or' and)* */
static PetscInt ParseOr(ExpressionParser *parser)
{
    PROFILE_HELPER_BEGIN;
    PetscInt left = ParseAnd(parser);
    while (Accept(parser, "or")) left = NewNode(parser, OP_OR, 2, left, ParseAnd(parser), -1);
    PROFILE_HELPER_END;
    return left;
}

static PetscInt CountCode(const ExpressionParser *parser, PetscInt index);

/** @brief Emit postfix code for a subtree and track the stack depth it needs. */
static void EmitNode(const ExpressionParser *parser, PetscInt index, PicurvExpression *expression,
                     PetscInt *depth)
{
    PROFILE_FUNCTION_BEGIN;
    const ExpressionNode *node = &parser->nodes[index];

    if (node->op == OP_SUM) {
        const PetscInt body = CountCode(parser, node->child[0]);
        ExpressionInstruction *start = &expression->code[expression->length++];
        start->op = OP_SUM_START; start->value = node->value; start->arg = node->var;
        ExpressionInstruction *hi = &expression->code[expression->length++];
        hi->op = OP_SUM_HI; hi->value = node->upper; hi->arg = body;
        *depth += 1;  /* the running total */
        if (*depth > expression->max_depth) expression->max_depth = *depth;
        EmitNode(parser, node->child[0], expression, depth);
        ExpressionInstruction *acc = &expression->code[expression->length++];
        acc->op = OP_SUM_ACC; acc->value = (PetscReal)body; acc->arg = node->var;
        *depth -= 1;
        PROFILE_FUNCTION_END;
        return;
    }
    for (PetscInt c = 0; c < node->nchild; ++c) EmitNode(parser, node->child[c], expression, depth);
    ExpressionInstruction *instruction = &expression->code[expression->length++];
    instruction->op = node->op;
    instruction->value = node->value;
    instruction->arg = node->var;
    *depth += 1 - node->nchild;
    if (*depth > expression->max_depth) expression->max_depth = *depth;
    if (node->op == OP_UNIFORM || node->op == OP_NORMAL) expression->draws = PETSC_TRUE;
    PROFILE_FUNCTION_END;
}

/** @brief Count the instructions a subtree emits (shared subtrees count per use). */
static PetscInt CountCode(const ExpressionParser *parser, PetscInt index)
{
    PROFILE_HELPER_BEGIN;
    const ExpressionNode *node = &parser->nodes[index];
    PetscInt total = (node->op == OP_SUM) ? 3 : 1;  /* a series: START, HI and ACC around its body */

    for (PetscInt c = 0; c < node->nchild; ++c) total += CountCode(parser, node->child[c]);
    PROFILE_HELPER_END;
    return total;
}

/**
 * @brief Implementation of \ref PicurvExpressionCompile().
 * @see PicurvExpressionCompile()
 */
PetscErrorCode PicurvExpressionCompile(const char *text, const char *const *names, PetscInt name_count,
                                       PetscBool allow_random, PicurvExpression **expression)
{
    ExpressionParser parser;
    PetscInt         root, depth = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(text && expression, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Expression text and output are required.");
    PetscCall(PetscMemzero(&parser, sizeof(parser)));
    parser.text = text;
    parser.names = names;
    parser.name_count = name_count;
    parser.allow_random = allow_random;
    parser.capacity = 4 * (PetscInt)strlen(text) + 8;
    PetscCall(PetscMalloc1(parser.capacity, &parser.nodes));

    root = ParseOr(&parser);
    SkipSpace(&parser);
    if (!parser.error[0] && parser.text[parser.pos]) ParserFail(&parser, "unexpected text");
    if (parser.error[0]) {
        PetscCall(PetscFree(parser.nodes));
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid expression: %s.", parser.error);
    }

    PetscCall(PetscNew(expression));
    (*expression)->length = 0;
    PetscCall(PetscMalloc1(CountCode(&parser, root), &(*expression)->code));
    EmitNode(&parser, root, *expression, &depth);
    PetscCall(PetscFree(parser.nodes));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvMix64().
 * @see PicurvMix64()
 */
uint64_t PicurvMix64(uint64_t z)
{
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/**
 * @brief Implementation of \ref PicurvUnitFromHash().
 * @see PicurvUnitFromHash()
 */
PetscReal PicurvUnitFromHash(uint64_t hash)
{
    return (PetscReal)(hash >> 11) * (1.0 / 9007199254740992.0);
}

/**
 * @brief Implementation of \ref PicurvParticleDrawHash().
 * @see PicurvParticleDrawHash()
 */
uint64_t PicurvParticleDrawHash(PetscInt64 seed, PicurvDrawSpace space, uint64_t stream,
                                PetscInt64 pid, PetscInt64 salt)
{
    uint64_t hash = PicurvMix64((uint64_t)seed);
    hash = PicurvMix64(hash ^ (uint64_t)space);
    hash = PicurvMix64(hash ^ stream);
    hash = PicurvMix64(hash ^ (uint64_t)pid);
    return PicurvMix64(hash ^ (uint64_t)salt);
}

/** @brief A draw keyed on the particle, the event, the stream, and the function. */
static PetscReal Draw(const PicurvExpressionDrawKey *key, const ExpressionInstruction *instruction)
{
    PROFILE_HELPER_BEGIN;
    /* A default stream is the drawing field's own, so two fields draw independently; an
       explicit stream is shared by every expression that names it. */
    const PicurvDrawSpace space = instruction->arg ? PICURV_DRAW_SPACE_NAMED_STREAM : PICURV_DRAW_SPACE_FIELD;
    const uint64_t stream = instruction->arg ? (uint64_t)instruction->value : (uint64_t)key->field;
    const uint64_t hash = PicurvMix64(PicurvParticleDrawHash(key->seed, space, stream, key->pid, key->salt)
                                      ^ (uint64_t)instruction->op);
    if (instruction->op == OP_UNIFORM) { PROFILE_HELPER_END; return PicurvUnitFromHash(hash); }
    /* Box-Muller from two independent uniforms; u1 lies in (0, 1]. */
    const PetscReal u1 = 1.0 - PicurvUnitFromHash(hash);
    const PetscReal u2 = PicurvUnitFromHash(PicurvMix64(hash ^ 0xD1B54A32D192ED03ULL));
    PROFILE_HELPER_END;
    return PetscSqrtReal(-2.0 * PetscLogReal(u1)) * PetscCosReal(2.0 * PETSC_PI * u2);
}

/**
 * @brief Kelvin functions of order zero, `ber(x) + i bei(x) = J0(x e^{3 pi i / 4})`, from
 *        their power series. Terms grow to about `exp(|x| / sqrt(2))` before the series
 *        converges, so about `|x| / 5` of the sixteen digits are lost to cancellation.
 */
static PetscReal Kelvin(PetscReal x, PetscBool imaginary)
{
    PROFILE_HELPER_BEGIN;
    const PetscReal q = 0.0625 * x * x * x * x;  /* (x/2)^4 */
    PetscReal term = imaginary ? 0.25 * x * x : 1.0, total = term;
    /* ber: (-1)^m (x/2)^(4m) / ((2m)!)^2; bei: (-1)^m (x/2)^(4m+2) / ((2m+1)!)^2. */
    for (PetscInt m = 1; m < 200; ++m) {
        const PetscReal a = imaginary ? 2.0 * m : 2.0 * m - 1.0, b = a + 1.0;
        term *= -q / (a * a * b * b);
        total += term;
        if (PetscAbsReal(term) <= 1.0e-17 * PetscAbsReal(total)) break;
    }
    PROFILE_HELPER_END;
    return total;
}

/** @brief Bessel function of the first or second kind of integer order; NaN for another order. */
static PetscReal Bessel(PetscReal order, PetscReal x, PetscBool second)
{
    PROFILE_HELPER_BEGIN;
    if (order != PetscFloorReal(order) || PetscAbsReal(order) > 1.0e6) { PROFILE_HELPER_END; return NAN; }
    const int n = (int)order;
    PROFILE_HELPER_END;
    return second ? (PetscReal)yn(n, (double)x) : (PetscReal)jn(n, (double)x);
}

/** @brief Maximum or minimum that propagates NaN, as numpy does. */
static PetscReal Extreme(PetscReal a, PetscReal b, PetscBool maximum)
{
    PROFILE_HELPER_BEGIN;
    if (PetscIsNanReal(a) || PetscIsNanReal(b)) { PROFILE_HELPER_END; return NAN; }
    PROFILE_HELPER_END;
    return maximum ? PetscMax(a, b) : PetscMin(a, b);
}

/**
 * @brief Implementation of \ref PicurvExpressionEvaluate().
 * @see PicurvExpressionEvaluate()
 */
PetscErrorCode PicurvExpressionEvaluate(const PicurvExpression *expression, const PetscReal *values,
                                        const PicurvExpressionDrawKey *key, PetscReal *result)
{
    PetscReal stack[expression->max_depth + 1];
    PetscReal index[EXPRESSION_MAX_SERIES_DEPTH];
    PetscInt  top = 0;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(!expression->draws || key, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "An expression that draws random numbers needs a draw key.");
    for (PetscInt i = 0; i < expression->length; ++i) {
        const ExpressionInstruction *instruction = &expression->code[i];
        PetscReal a, b, c;

        switch (instruction->op) {
        case OP_CONST: stack[top++] = instruction->value; break;
        case OP_VAR: stack[top++] = values[instruction->arg]; break;
        case OP_UNIFORM:
        case OP_NORMAL: stack[top++] = Draw(key, instruction); break;
        case OP_NEG: stack[top - 1] = -stack[top - 1]; break;
        case OP_POS: break;
        case OP_NOT: stack[top - 1] = (stack[top - 1] == 0.0) ? 1.0 : 0.0; break;
        case OP_ABS: stack[top - 1] = PetscAbsReal(stack[top - 1]); break;
        case OP_COS: stack[top - 1] = PetscCosReal(stack[top - 1]); break;
        case OP_EXP: stack[top - 1] = PetscExpReal(stack[top - 1]); break;
        case OP_SIN: stack[top - 1] = PetscSinReal(stack[top - 1]); break;
        case OP_SQRT: stack[top - 1] = PetscSqrtReal(stack[top - 1]); break;
        case OP_TAN: stack[top - 1] = PetscTanReal(stack[top - 1]); break;
        case OP_ERF: stack[top - 1] = (PetscReal)erf((double)stack[top - 1]); break;
        case OP_TANH: stack[top - 1] = PetscTanhReal(stack[top - 1]); break;
        case OP_COSH: stack[top - 1] = PetscCoshReal(stack[top - 1]); break;
        case OP_SINH: stack[top - 1] = PetscSinhReal(stack[top - 1]); break;
        case OP_LOG: stack[top - 1] = PetscLogReal(stack[top - 1]); break;
        case OP_FLOOR: stack[top - 1] = PetscFloorReal(stack[top - 1]); break;
        case OP_BER: stack[top - 1] = Kelvin(stack[top - 1], PETSC_FALSE); break;
        case OP_BEI: stack[top - 1] = Kelvin(stack[top - 1], PETSC_TRUE); break;
        case OP_LOCAL: stack[top++] = index[instruction->arg]; break;
        case OP_SUM_START:
            stack[top++] = 0.0;
            index[instruction->arg] = instruction->value;
            /* An empty range skips the body and its ACC. */
            if (instruction->value > expression->code[i + 1].value) i += 2 + expression->code[i + 1].arg;
            else i += 1;
            break;
        case OP_SUM_ACC: {
            const PetscInt  body = (PetscInt)instruction->value;
            const PetscReal upper = expression->code[i - body - 1].value;
            const PetscReal term = stack[--top];
            stack[top - 1] += term;
            index[instruction->arg] += 1.0;
            if (index[instruction->arg] <= upper) i -= body + 1;
        } break;
        case OP_WHERE:
            c = stack[--top]; b = stack[--top]; a = stack[top - 1];
            stack[top - 1] = (a != 0.0) ? b : c;
            break;
        default:
            b = stack[--top]; a = stack[top - 1];
            switch (instruction->op) {
            case OP_ADD: a = a + b; break;
            case OP_SUB: a = a - b; break;
            case OP_MUL: a = a * b; break;
            case OP_DIV: a = a / b; break;
            /* Python's modulus takes the sign of the divisor. */
            case OP_MOD: a = (b == 0.0) ? NAN : a - b * PetscFloorReal(a / b); break;
            case OP_POW: a = PetscPowReal(a, b); break;
            case OP_EQ: a = (a == b) ? 1.0 : 0.0; break;
            case OP_NE: a = (a != b) ? 1.0 : 0.0; break;
            case OP_LT: a = (a < b) ? 1.0 : 0.0; break;
            case OP_LE: a = (a <= b) ? 1.0 : 0.0; break;
            case OP_GT: a = (a > b) ? 1.0 : 0.0; break;
            case OP_GE: a = (a >= b) ? 1.0 : 0.0; break;
            case OP_AND: a = (a != 0.0 && b != 0.0) ? 1.0 : 0.0; break;
            case OP_OR: a = (a != 0.0 || b != 0.0) ? 1.0 : 0.0; break;
            case OP_ATAN2: a = PetscAtan2Real(a, b); break;
            case OP_MAXIMUM: a = Extreme(a, b, PETSC_TRUE); break;
            case OP_MINIMUM: a = Extreme(a, b, PETSC_FALSE); break;
            case OP_BESSELJ: a = Bessel(a, b, PETSC_FALSE); break;
            case OP_BESSELY: a = Bessel(a, b, PETSC_TRUE); break;
            default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Unknown expression opcode %d.", (int)instruction->op);
            }
            stack[top - 1] = a;
        }
    }
    *result = stack[0];
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvExpressionDestroy().
 * @see PicurvExpressionDestroy()
 */
PetscErrorCode PicurvExpressionDestroy(PicurvExpression **expression)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (!expression || !*expression) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    PetscCall(PetscFree((*expression)->code));
    PetscCall(PetscFree(*expression));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Implementation of \ref PicurvExpressionUsesVariable().
 * @see PicurvExpressionUsesVariable()
 */
PetscErrorCode PicurvExpressionUsesVariable(const PicurvExpression *expression, PetscInt index, PetscBool *uses)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    *uses = PETSC_FALSE;
    for (PetscInt i = 0; i < expression->length && !*uses; ++i) {
        if (expression->code[i].op == OP_VAR && expression->code[i].arg == index) *uses = PETSC_TRUE;
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/*                         Prescribed velocity sources                        */
/* ------------------------------------------------------------------------- */

const char *const kPicurvSourceVariableNames[PICURV_SOURCE_VAR_COUNT] = {
    "x", "y", "z", "t", "xi", "eta", "zeta", "L_xi", "L_eta", "L_zeta"};

#define PICURV_SOURCE_EXPRESSION_LENGTH 4096

/**
 * @brief Implementation of \ref PicurvVelocitySourceCreate().
 * @see PicurvVelocitySourceCreate()
 */
PetscErrorCode PicurvVelocitySourceCreate(const char *prefix, PetscBool on_face, PetscBool has_data,
                                          PicurvVelocitySource *source)
{
    static const char *const kTime[] = {"t"};
    char      option[PETSC_MAX_PATH_LEN], frame[32], text[PICURV_SOURCE_EXPRESSION_LENGTH];
    PetscBool found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PetscMemzero(source, sizeof(*source)));
    PetscCall(PetscSNPrintf(option, sizeof(option), "-%s_frame", prefix));
    PetscCall(PetscOptionsGetString(NULL, NULL, option, frame, sizeof(frame), &found));
    PetscCheck(found || has_data, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Velocity source option %s is required.", option);
    if (!found) PetscCall(PetscStrncpy(frame, on_face ? "normal" : "cartesian", sizeof(frame)));
    if (!strcmp(frame, "cartesian")) source->frame = PICURV_SOURCE_FRAME_CARTESIAN;
    else if (!strcmp(frame, "curvilinear")) source->frame = PICURV_SOURCE_FRAME_CURVILINEAR;
    else if (!strcmp(frame, "normal") && on_face) source->frame = PICURV_SOURCE_FRAME_NORMAL;
    else if (!strcmp(frame, "contravariant") && !on_face) source->frame = PICURV_SOURCE_FRAME_CONTRAVARIANT;
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "%s is '%s'; expected cartesian, curvilinear, or %s.",
                 option, frame, on_face ? "normal" : "contravariant");
    source->count = (source->frame == PICURV_SOURCE_FRAME_NORMAL) ? 1 : 3;
    for (PetscInt c = 0; c < source->count; ++c) {
        PetscBool uses_t = PETSC_FALSE;
        PetscCall(PetscSNPrintf(option, sizeof(option), "-%s_expr_%d", prefix, (int)c));
        PetscCall(PetscOptionsGetString(NULL, NULL, option, text, sizeof(text), &found));
        if (!found && has_data && c == 0) { source->count = 0; break; }
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Velocity source option %s is required.", option);
        PetscCall(PicurvExpressionCompile(text, kPicurvSourceVariableNames, PICURV_SOURCE_VAR_COUNT, PETSC_FALSE,
                                          &source->component[c]));
        PetscCall(PicurvExpressionUsesVariable(source->component[c], PICURV_SOURCE_VAR_T, &uses_t));
        source->time_dependent = (PetscBool)(source->time_dependent || uses_t);
    }
    PetscCall(PetscSNPrintf(option, sizeof(option), "-%s_bulk_velocity", prefix));
    PetscCall(PetscOptionsGetString(NULL, NULL, option, text, sizeof(text), &found));
    if (found) {
        PetscCall(PicurvExpressionCompile(text, kTime, 1, PETSC_FALSE, &source->bulk));
        PetscCall(PicurvExpressionUsesVariable(source->bulk, 0, &source->bulk_time_dependent));
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvVelocitySourceBulkFactor().
 * @see PicurvVelocitySourceBulkFactor()
 */
PetscErrorCode PicurvVelocitySourceBulkFactor(const PicurvVelocitySource *source, PetscReal physical_time,
                                              const PetscReal local[3], const char *where, PetscReal *factor)
{
    PetscReal global[3], target;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCallMPI(MPI_Allreduce(local, global, 3, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    PetscCheck(global[2] > 0.0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               "%s has no open area or volume to set a bulk velocity on.", where);
    PetscCheck(global[0] != 0.0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               "%s: the source's mean speed is zero, so it cannot be scaled to a bulk velocity.", where);
    PetscCall(PicurvExpressionEvaluate(source->bulk, &physical_time, NULL, &target));
    *factor = (target * global[2] - global[1]) / global[0];
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "%s: bulk velocity %g, source mean %g, data mean %g, factor %g.\n", where,
              (double)target, (double)(global[0] / global[2]), (double)(global[1] / global[2]), (double)*factor);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvVelocitySourceEvaluate().
 * @see PicurvVelocitySourceEvaluate()
 */
PetscErrorCode PicurvVelocitySourceEvaluate(const PicurvVelocitySource *source, const PetscReal *values,
                                            const Cmpnts edges[3], Cmpnts *value)
{
    PetscReal v[3] = {0.0, 0.0, 0.0};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (PetscInt c = 0; c < source->count; ++c) {
        PetscCall(PicurvExpressionEvaluate(source->component[c], values, NULL, &v[c]));
        PetscCheck(!PetscIsInfOrNanReal(v[c]), PETSC_COMM_SELF, PETSC_ERR_FP,
                   "Velocity source value %" PetscInt_FMT " is %g at (x, y, z) = (%g, %g, %g), t = %g.",
                   c, (double)v[c], (double)values[PICURV_SOURCE_VAR_X], (double)values[PICURV_SOURCE_VAR_Y],
                   (double)values[PICURV_SOURCE_VAR_Z], (double)values[PICURV_SOURCE_VAR_T]);
    }
    if (source->frame == PICURV_SOURCE_FRAME_CURVILINEAR) {
        /* Values are components along the unit grid-line directions. */
        value->x = value->y = value->z = 0.0;
        for (PetscInt a = 0; a < 3; ++a) {
            const PetscReal norm = PetscSqrtReal(edges[a].x * edges[a].x + edges[a].y * edges[a].y + edges[a].z * edges[a].z);
            PetscCheck(norm > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "A cell edge has zero length.");
            value->x += v[a] * edges[a].x / norm;
            value->y += v[a] * edges[a].y / norm;
            value->z += v[a] * edges[a].z / norm;
        }
    } else {
        value->x = v[0];
        value->y = v[1];
        value->z = v[2];
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvSourceCellBulkFactor().
 * @see PicurvSourceCellBulkFactor()
 */
PetscErrorCode PicurvSourceCellBulkFactor(UserCtx *user, const PicurvVelocitySource *source, Vec fraction, Vec length,
                                      PetscReal physical_time, PetscReal *factor)
{
    const DMDALocalInfo *info = &user->info;
    const PetscReal L_ref = user->simCtx->scaling.L_ref, U_ref = user->simCtx->scaling.U_ref;
    const PetscInt xs = info->xs, xe = info->xs + info->xm, mx = info->mx;
    const PetscInt ys = info->ys, ye = info->ys + info->ym, my = info->my;
    const PetscInt zs = info->zs, ze = info->zs + info->zm, mz = info->mz;
    PetscReal values[PICURV_SOURCE_VAR_COUNT], local[3] = {0.0, 0.0, 0.0};
    const Cmpnts ***cent, ***frac, ***len, ***csi, ***eta, ***zet, ***ucat;
    const PetscReal ***aj;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Cent, &cent));
    PetscCall(DMDAVecGetArrayRead(user->fda, fraction, &frac));
    PetscCall(DMDAVecGetArrayRead(user->fda, length, &len));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, &csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, &eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, &zet));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    values[PICURV_SOURCE_VAR_T] = physical_time;
    for (PetscInt k = (zs == 0) ? 1 : zs; k < ((ze == mz) ? mz - 1 : ze); ++k)
    for (PetscInt j = (ys == 0) ? 1 : ys; j < ((ye == my) ? my - 1 : ye); ++j)
    for (PetscInt i = (xs == 0) ? 1 : xs; i < ((xe == mx) ? mx - 1 : xe); ++i) {
        const PetscReal volume = 1.0 / PetscAbsReal(aj[k][j][i]);
        Cmpnts edges[3], v;
        PicurvSourceSetVariables(values, cent[k][j][i], frac[k][j][i], len[k][j][i], L_ref);
        PetscCall(ComputeCellEdgeVectors(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], edges));
        PetscCall(PicurvVelocitySourceEvaluate(source, values, edges, &v));
        const PetscReal speed = PetscSqrtReal(v.x * v.x + v.y * v.y + v.z * v.z);
        local[0] += speed * volume;
        if (speed > 0.0) local[1] += U_ref * (ucat[k][j][i].x * v.x + ucat[k][j][i].y * v.y +
                                              ucat[k][j][i].z * v.z) / speed * volume;
        local[2] += volume;
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, length, &len));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, fraction, &frac));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent));
    PetscCall(PicurvVelocitySourceBulkFactor(source, physical_time, local, "The field", factor));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvSourceAddFluxes().
 * @see PicurvSourceAddFluxes()
 */
PetscErrorCode PicurvSourceAddFluxes(UserCtx *user, const PicurvVelocitySource *source, Vec fraction, Vec length,
                                     PetscReal physical_time, PetscReal factor, PetscBool boundary_faces)
{
    const DMDALocalInfo *info = &user->info;
    const PetscReal L_ref = user->simCtx->scaling.L_ref, U_ref = user->simCtx->scaling.U_ref;
    const PetscInt start[3] = {info->xs, info->ys, info->zs};
    const PetscInt end[3] = {info->xs + info->xm, info->ys + info->ym, info->zs + info->zm};
    const PetscInt nodes[3] = {info->mx, info->my, info->mz};
    const PetscInt lo[3] = {(start[0] == 0) ? 1 : start[0], (start[1] == 0) ? 1 : start[1], (start[2] == 0) ? 1 : start[2]};
    const PetscInt hi[3] = {(end[0] == nodes[0]) ? nodes[0] - 1 : end[0], (end[1] == nodes[1]) ? nodes[1] - 1 : end[1],
                            (end[2] == nodes[2]) ? nodes[2] - 1 : end[2]};
    const Cmpnts ***frac, ***len, ***csi, ***eta, ***zet, ***centx, ***centy, ***centz;
    const PetscReal ***aj;
    Cmpnts ***ucont;
    PetscReal values[PICURV_SOURCE_VAR_COUNT];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(DMDAVecGetArrayRead(user->fda, fraction, &frac));
    PetscCall(DMDAVecGetArrayRead(user->fda, length, &len));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, &csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, &eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, &zet));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCentx, &centx));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCenty, &centy));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCentz, &centz));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, &aj));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    values[PICURV_SOURCE_VAR_T] = physical_time;
    for (PetscInt k = lo[2]; k < hi[2]; ++k)
    for (PetscInt j = lo[1]; j < hi[1]; ++j)
    for (PetscInt i = lo[0]; i < hi[0]; ++i) {
        const PetscInt cell[3] = {i, j, k};
        const PetscReal lengths[3] = {len[k][j][i].x, len[k][j][i].y, len[k][j][i].z};
        Cmpnts edges[3];
        PetscCall(ComputeCellEdgeVectors(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], edges));
        for (PetscInt a = 0; a < 3; ++a) {
            const PetscReal half = 0.5 * PetscSqrtReal(edges[a].x * edges[a].x + edges[a].y * edges[a].y +
                                                       edges[a].z * edges[a].z) / lengths[a];
            /* Each cell owns its + face; the first cell of a line also owns the boundary - face. */
            const PetscInt sides = (boundary_faces && start[a] == 0 && cell[a] == 1) ? 2 : 1;
            for (PetscInt side = 0; side < sides; ++side) {
                const PetscInt offset = side ? -1 : 0;
                const PetscInt fi = i + (a == 0 ? offset : 0), fj = j + (a == 1 ? offset : 0), fk = k + (a == 2 ? offset : 0);
                const Cmpnts metric = (a == 0) ? csi[fk][fj][fi] : (a == 1) ? eta[fk][fj][fi] : zet[fk][fj][fi];
                const Cmpnts centre = (a == 0) ? centx[fk][fj][fi] : (a == 1) ? centy[fk][fj][fi] : centz[fk][fj][fi];
                PetscReal flux;
                PicurvSourceSetVariables(values, centre, frac[k][j][i], len[k][j][i], L_ref);
                values[PICURV_SOURCE_VAR_XI + a] += side ? -half : half;
                if (source->frame == PICURV_SOURCE_FRAME_CONTRAVARIANT) {
                    PetscReal value;
                    PetscCall(PicurvExpressionEvaluate(source->component[a], values, NULL, &value));
                    flux = factor * value / (U_ref * L_ref * L_ref);
                } else if (source->frame == PICURV_SOURCE_FRAME_CURVILINEAR) {
                    const PetscReal area = PetscSqrtReal(metric.x * metric.x + metric.y * metric.y + metric.z * metric.z);
                    PetscReal value;
                    PetscCall(PicurvExpressionEvaluate(source->component[a], values, NULL, &value));
                    flux = factor * (value / U_ref) * user->GridOrientation * area;
                } else {
                    Cmpnts v;
                    PetscCall(PicurvVelocitySourceEvaluate(source, values, edges, &v));
                    flux = (factor * v.x / U_ref) * metric.x + (factor * v.y / U_ref) * metric.y +
                           (factor * v.z / U_ref) * metric.z;
                }
                if (a == 0) ucont[fk][fj][fi].x += flux;
                if (a == 1) ucont[fk][fj][fi].y += flux;
                if (a == 2) ucont[fk][fj][fi].z += flux;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, &aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCentz, &centz));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCenty, &centy));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCentx, &centx));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, length, &len));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, fraction, &frac));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvVelocitySourceDestroy().
 * @see PicurvVelocitySourceDestroy()
 */
PetscErrorCode PicurvVelocitySourceDestroy(PicurvVelocitySource *source)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (PetscInt c = 0; c < 3; ++c) PetscCall(PicurvExpressionDestroy(&source->component[c]));
    PetscCall(PicurvExpressionDestroy(&source->bulk));
    source->count = 0;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
