/**
 * @file expression.c
 * @brief The expression language shared by every configured field: a parser, a postfix
 *        compiler, and a pure stack evaluator.
 *
 * Expressions are parsed into a small tree, because a chained comparison evaluates its
 * middle operand twice, and compiled to postfix code that a stack evaluator runs once per
 * point. Evaluation is pure: random draws are hashed from a caller-supplied key, so a value
 * does not depend on rank count, evaluation order, or how often it is evaluated.
 */

#include "expression.h"

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
    OP_UNIFORM, OP_NORMAL
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
    {"uniform", 0, 1, OP_UNIFORM, PETSC_TRUE},  {"normal", 0, 1, OP_NORMAL, PETSC_TRUE},
};

/** @brief One tree node; children index the node pool. */
typedef struct {
    ExpressionOp op;
    PetscReal    value;     /* OP_CONST value; explicit stream of a draw */
    PetscInt     var;       /* OP_VAR index; 1 for a draw on an explicit stream */
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

/** @brief number | name | name(args) | (expression) */
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

/** @brief Emit postfix code for a subtree and track the stack depth it needs. */
static void EmitNode(const ExpressionParser *parser, PetscInt index, PicurvExpression *expression,
                     PetscInt *depth)
{
    PROFILE_FUNCTION_BEGIN;
    const ExpressionNode *node = &parser->nodes[index];

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
    PetscInt total = 1;

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
PetscErrorCode PicurvVelocitySourceCreate(PicurvSourceTarget target, PetscInt id, PicurvVelocitySource *source)
{
    char      option[PETSC_MAX_PATH_LEN], frame[32], text[PICURV_SOURCE_EXPRESSION_LENGTH];
    PetscBool found = PETSC_FALSE;
    const PetscBool allow_normal = (PetscBool)(target == PICURV_SOURCE_TARGET_INLET);

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PetscMemzero(source, sizeof(*source)));
    if (target == PICURV_SOURCE_TARGET_INLET) PetscCall(PetscSNPrintf(option, sizeof(option), "-inlet_source_%d_frame", (int)id));
    else PetscCall(PetscSNPrintf(option, sizeof(option), "-ic_source_frame"));
    PetscCall(PetscOptionsGetString(NULL, NULL, option, frame, sizeof(frame), &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Velocity source option %s is required.", option);
    if (!strcmp(frame, "cartesian")) source->frame = PICURV_SOURCE_FRAME_CARTESIAN;
    else if (!strcmp(frame, "curvilinear")) source->frame = PICURV_SOURCE_FRAME_CURVILINEAR;
    else if (!strcmp(frame, "normal") && allow_normal) source->frame = PICURV_SOURCE_FRAME_NORMAL;
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "%s is '%s'; expected cartesian, curvilinear%s.",
                 option, frame, allow_normal ? ", or normal" : "");
    source->count = (source->frame == PICURV_SOURCE_FRAME_NORMAL) ? 1 : 3;
    for (PetscInt c = 0; c < source->count; ++c) {
        PetscBool uses_t = PETSC_FALSE;
        if (target == PICURV_SOURCE_TARGET_INLET)
            PetscCall(PetscSNPrintf(option, sizeof(option), "-inlet_source_%d_expr_%d", (int)id, (int)c));
        else PetscCall(PetscSNPrintf(option, sizeof(option), "-ic_source_expr_%d", (int)c));
        PetscCall(PetscOptionsGetString(NULL, NULL, option, text, sizeof(text), &found));
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT, "Velocity source option %s is required.", option);
        PetscCall(PicurvExpressionCompile(text, kPicurvSourceVariableNames, PICURV_SOURCE_VAR_COUNT, PETSC_FALSE,
                                          &source->component[c]));
        PetscCall(PicurvExpressionUsesVariable(source->component[c], PICURV_SOURCE_VAR_T, &uses_t));
        source->time_dependent = (PetscBool)(source->time_dependent || uses_t);
    }
    PROFILE_FUNCTION_END;
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
 * @brief Implementation of \ref PicurvVelocitySourceDestroy().
 * @see PicurvVelocitySourceDestroy()
 */
PetscErrorCode PicurvVelocitySourceDestroy(PicurvVelocitySource *source)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (PetscInt c = 0; c < 3; ++c) PetscCall(PicurvExpressionDestroy(&source->component[c]));
    source->count = 0;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
