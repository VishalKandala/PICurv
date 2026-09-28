/**
 * @file ParticleInitialConditions.c
 * @brief Configured initial values of particle-carried fields, and the expression
 *        language that defines them.
 *
 * Expressions are parsed into a small tree, because a chained comparison evaluates its
 * middle operand twice, and compiled to postfix code that a stack evaluator runs once per
 * particle. Evaluation is pure: random draws are hashed from the particle's identity, so
 * a value does not depend on rank count, particle order, or how often it is evaluated.
 */

#include "ParticleInitialConditions.h"

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
    OP_ABS, OP_COS, OP_EXP, OP_SIN, OP_SQRT, OP_TAN, OP_MAXIMUM, OP_MINIMUM, OP_WHERE,
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
    if (!parser->error[0]) {
        snprintf(parser->error, sizeof(parser->error), "%s at position %zu of '%s'",
                 message, parser->pos + 1, parser->text);
    }
}

/** @brief Skip whitespace before the next token. */
static void SkipSpace(ExpressionParser *parser)
{
    while (parser->text[parser->pos] && isspace((unsigned char)parser->text[parser->pos])) parser->pos++;
}

/** @brief Consume @p token if it is next, and report whether it was. */
static PetscBool Accept(ExpressionParser *parser, const char *token)
{
    size_t length = strlen(token);

    SkipSpace(parser);
    if (strncmp(parser->text + parser->pos, token, length) != 0) return PETSC_FALSE;
    if (isalpha((unsigned char)token[0])) {
        const char next = parser->text[parser->pos + length];
        if (isalnum((unsigned char)next) || next == '_') return PETSC_FALSE;
    }
    /* `*` must not match the first half of `**`, nor `<` the first half of `<=`. */
    if ((!strcmp(token, "*") && parser->text[parser->pos + 1] == '*') ||
        ((!strcmp(token, "<") || !strcmp(token, ">")) && parser->text[parser->pos + 1] == '=')) {
        return PETSC_FALSE;
    }
    parser->pos += length;
    return PETSC_TRUE;
}

/** @brief Append a node to the pool and return its index, or -1 once an error is set. */
static PetscInt NewNode(ExpressionParser *parser, ExpressionOp op, PetscInt nchild,
                        PetscInt a, PetscInt b, PetscInt c)
{
    if (parser->error[0]) return -1;
    if (parser->count == parser->capacity) {
        ParserFail(parser, "expression is too long");
        return -1;
    }
    ExpressionNode *node = &parser->nodes[parser->count];
    memset(node, 0, sizeof(*node));
    node->op = op;
    node->nchild = nchild;
    node->child[0] = a;
    node->child[1] = b;
    node->child[2] = c;
    return parser->count++;
}

static PetscInt ParseOr(ExpressionParser *parser);

/** @brief number | name | name(args) | (expression) */
static PetscInt ParseAtom(ExpressionParser *parser)
{
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
            if (!isdigit((unsigned char)start[exponent])) { ParserFail(parser, "malformed number"); return -1; }
            while (isdigit((unsigned char)start[exponent])) exponent++;
            length = exponent;
        }
        if (isalnum((unsigned char)start[length]) || start[length] == '_' || start[length] == '.') {
            ParserFail(parser, "malformed number");
            return -1;
        }
        const PetscReal value = strtod(start, NULL);
        parser->pos += length;
        PetscInt node = NewNode(parser, OP_CONST, 0, -1, -1, -1);
        if (node >= 0) parser->nodes[node].value = value;
        return node;
    }
    if (isalpha((unsigned char)*start) || *start == '_') {
        char name[64];
        size_t length = 0;
        while (isalnum((unsigned char)start[length]) || start[length] == '_') length++;
        if (length >= sizeof(name)) { ParserFail(parser, "name is too long"); return -1; }
        memcpy(name, start, length);
        name[length] = '\0';
        if (!strcmp(name, "and") || !strcmp(name, "or") || !strcmp(name, "not")) {
            ParserFail(parser, "operator where an operand was expected");
            return -1;
        }
        parser->pos += length;
        if (Accept(parser, "(")) {
            for (size_t f = 0; f < sizeof(kFunctions) / sizeof(kFunctions[0]); ++f) {
                if (strcmp(name, kFunctions[f].name)) continue;
                if (kFunctions[f].random && !parser->allow_random) {
                    ParserFail(parser, "random draws are not available in this expression");
                    return -1;
                }
                PetscInt args[3] = {-1, -1, -1}, nargs = 0;
                if (!Accept(parser, ")")) {
                    do {
                        if (nargs == 3) { ParserFail(parser, "too many arguments"); return -1; }
                        args[nargs++] = ParseOr(parser);
                        if (parser->error[0]) return -1;
                    } while (Accept(parser, ","));
                    if (!Accept(parser, ")")) { ParserFail(parser, "expected ')'"); return -1; }
                }
                if (nargs < kFunctions[f].min_args || nargs > kFunctions[f].max_args) {
                    ParserFail(parser, "wrong number of arguments");
                    return -1;
                }
                if (kFunctions[f].random) {
                    /* A draw's stream is a constant so that two expressions naming the
                       same stream draw the same number for the same particle. */
                    PetscInt node = NewNode(parser, kFunctions[f].op, 0, -1, -1, -1);
                    if (node < 0) return -1;
                    if (nargs == 1) {
                        const ExpressionNode *stream = &parser->nodes[args[0]];
                        if (stream->op != OP_CONST || stream->value < 0.0 ||
                            stream->value != PetscFloorReal(stream->value)) {
                            ParserFail(parser, "a draw's stream must be a non-negative integer literal");
                            return -1;
                        }
                        parser->nodes[node].value = stream->value;
                        parser->nodes[node].var = 1;
                    }
                    return node;
                }
                return NewNode(parser, kFunctions[f].op, nargs, args[0], args[1], args[2]);
            }
            ParserFail(parser, "unknown function");
            return -1;
        }
        if (!strcmp(name, "pi")) {
            PetscInt node = NewNode(parser, OP_CONST, 0, -1, -1, -1);
            if (node >= 0) parser->nodes[node].value = PETSC_PI;
            return node;
        }
        for (PetscInt v = 0; v < parser->name_count; ++v) {
            if (!strcmp(name, parser->names[v])) {
                PetscInt node = NewNode(parser, OP_VAR, 0, -1, -1, -1);
                if (node >= 0) parser->nodes[node].var = v;
                return node;
            }
        }
        ParserFail(parser, "unknown name");
        return -1;
    }
    if (Accept(parser, "(")) {
        PetscInt inner = ParseOr(parser);
        if (!Accept(parser, ")")) { ParserFail(parser, "expected ')'"); return -1; }
        return inner;
    }
    ParserFail(parser, "expected a number, name, or '('");
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
    PetscInt base = ParseAtom(parser);
    if (Accept(parser, "**")) return NewNode(parser, OP_POW, 2, base, ParseUnary(parser), -1);
    return base;
}

/**
 * @brief Parses a sign prefix recursively and wraps the operand in a negate or identity node.
 * @param[in,out] parser Parser state; its position advances past the operand.
 * @return Node index of the signed operand, or -1 after an error.
 */
static PetscInt ParseUnary(ExpressionParser *parser)
{
    if (Accept(parser, "-")) return NewNode(parser, OP_NEG, 1, ParseUnary(parser), -1, -1);
    if (Accept(parser, "+")) return NewNode(parser, OP_POS, 1, ParseUnary(parser), -1, -1);
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
    PetscInt left = ParseUnary(parser);
    for (;;) {
        if (Accept(parser, "*")) left = NewNode(parser, OP_MUL, 2, left, ParseUnary(parser), -1);
        else if (Accept(parser, "/")) {
            if (parser->text[parser->pos] == '/') { ParserFail(parser, "'//' is not supported"); return -1; }
            left = NewNode(parser, OP_DIV, 2, left, ParseUnary(parser), -1);
        } else if (Accept(parser, "%")) left = NewNode(parser, OP_MOD, 2, left, ParseUnary(parser), -1);
        else return left;
    }
}

/**
 * @brief Parses a left-associative chain of `+` and `-`, folding each operator into a node.
 * @param[in,out] parser Parser state; its position advances past the sum.
 * @return Node index of the sum, or -1 after an error.
 */
static PetscInt ParseSum(ExpressionParser *parser)
{
    PetscInt left = ParseTerm(parser);
    for (;;) {
        if (Accept(parser, "+")) left = NewNode(parser, OP_ADD, 2, left, ParseTerm(parser), -1);
        else if (Accept(parser, "-")) left = NewNode(parser, OP_SUB, 2, left, ParseTerm(parser), -1);
        else return left;
    }
}

/** @brief sum (comparison sum)*; a chain `a < b < c` is `(a < b) and (b < c)`. */
static PetscInt ParseComparison(ExpressionParser *parser)
{
    static const struct { const char *token; ExpressionOp op; } kComparisons[] = {
        {"==", OP_EQ}, {"!=", OP_NE}, {"<=", OP_LE}, {">=", OP_GE}, {"<", OP_LT}, {">", OP_GT},
    };
    PetscInt left = ParseSum(parser), result = -1;

    for (;;) {
        ExpressionOp op = OP_CONST;
        for (size_t c = 0; c < sizeof(kComparisons) / sizeof(kComparisons[0]); ++c) {
            if (Accept(parser, kComparisons[c].token)) { op = kComparisons[c].op; break; }
        }
        if (op == OP_CONST) return (result < 0) ? left : result;
        const PetscInt right = ParseSum(parser);
        const PetscInt pair = NewNode(parser, op, 2, left, right, -1);
        result = (result < 0) ? pair : NewNode(parser, OP_AND, 2, result, pair, -1);
        left = right;
    }
}

/** @brief 'not' not | comparison */
static PetscInt ParseNot(ExpressionParser *parser)
{
    if (Accept(parser, "not")) return NewNode(parser, OP_NOT, 1, ParseNot(parser), -1, -1);
    return ParseComparison(parser);
}

/** @brief not ('and' not)* */
static PetscInt ParseAnd(ExpressionParser *parser)
{
    PetscInt left = ParseNot(parser);
    while (Accept(parser, "and")) left = NewNode(parser, OP_AND, 2, left, ParseNot(parser), -1);
    return left;
}

/** @brief and ('or' and)* */
static PetscInt ParseOr(ExpressionParser *parser)
{
    PetscInt left = ParseAnd(parser);
    while (Accept(parser, "or")) left = NewNode(parser, OP_OR, 2, left, ParseAnd(parser), -1);
    return left;
}

/** @brief Emit postfix code for a subtree and track the stack depth it needs. */
static void EmitNode(const ExpressionParser *parser, PetscInt index, PicurvExpression *expression,
                     PetscInt *depth)
{
    const ExpressionNode *node = &parser->nodes[index];

    for (PetscInt c = 0; c < node->nchild; ++c) EmitNode(parser, node->child[c], expression, depth);
    ExpressionInstruction *instruction = &expression->code[expression->length++];
    instruction->op = node->op;
    instruction->value = node->value;
    instruction->arg = node->var;
    *depth += 1 - node->nchild;
    if (*depth > expression->max_depth) expression->max_depth = *depth;
    if (node->op == OP_UNIFORM || node->op == OP_NORMAL) expression->draws = PETSC_TRUE;
}

/** @brief Count the instructions a subtree emits (shared subtrees count per use). */
static PetscInt CountCode(const ExpressionParser *parser, PetscInt index)
{
    const ExpressionNode *node = &parser->nodes[index];
    PetscInt total = 1;

    for (PetscInt c = 0; c < node->nchild; ++c) total += CountCode(parser, node->child[c]);
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
    PetscFunctionReturn(0);
}

/** @brief splitmix64 finalizer: a bijective mix of 64 bits. */
static uint64_t Mix64(uint64_t z)
{
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/** @brief Uniform on [0, 1) from the top 53 bits of a hash. */
static PetscReal UnitFromHash(uint64_t hash)
{
    return (PetscReal)(hash >> 11) * (1.0 / 9007199254740992.0);
}

/** @brief A draw keyed on the particle, the event, the stream, and the function. */
static PetscReal Draw(const PicurvExpressionDrawKey *key, const ExpressionInstruction *instruction)
{
    /* A default stream is the drawing field's own, so two fields draw independently; an
       explicit stream is shared by every expression that names it. */
    const uint64_t space = instruction->arg ? 1u : 0u;
    const uint64_t stream = instruction->arg ? (uint64_t)instruction->value : (uint64_t)key->field;
    uint64_t hash = Mix64((uint64_t)key->seed);
    hash = Mix64(hash ^ space);
    hash = Mix64(hash ^ stream);
    hash = Mix64(hash ^ (uint64_t)key->pid);
    hash = Mix64(hash ^ (uint64_t)key->salt);
    hash = Mix64(hash ^ (uint64_t)instruction->op);
    if (instruction->op == OP_UNIFORM) return UnitFromHash(hash);
    /* Box-Muller from two independent uniforms; u1 lies in (0, 1]. */
    const PetscReal u1 = 1.0 - UnitFromHash(hash);
    const PetscReal u2 = UnitFromHash(Mix64(hash ^ 0xD1B54A32D192ED03ULL));
    return PetscSqrtReal(-2.0 * PetscLogReal(u1)) * PetscCosReal(2.0 * PETSC_PI * u2);
}

/** @brief Maximum or minimum that propagates NaN, as numpy does. */
static PetscReal Extreme(PetscReal a, PetscReal b, PetscBool maximum)
{
    if (PetscIsNanReal(a) || PetscIsNanReal(b)) return NAN;
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
            case OP_MAXIMUM: a = Extreme(a, b, PETSC_TRUE); break;
            case OP_MINIMUM: a = Extreme(a, b, PETSC_FALSE); break;
            default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Unknown expression opcode %d.", (int)instruction->op);
            }
            stack[top - 1] = a;
        }
    }
    *result = stack[0];
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvExpressionDestroy().
 * @see PicurvExpressionDestroy()
 */
PetscErrorCode PicurvExpressionDestroy(PicurvExpression **expression)
{
    PetscFunctionBeginUser;
    if (!expression || !*expression) PetscFunctionReturn(0);
    PetscCall(PetscFree((*expression)->code));
    PetscCall(PetscFree(*expression));
    PetscFunctionReturn(0);
}

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
    PetscCheck(plan, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Plan output is required.");
    *plan = NULL;
    PetscCall(PetscSNPrintf(option, sizeof(option), "-particle_fields_count"));
    PetscCall(PetscOptionsGetInt(NULL, NULL, option, &count, &found));
    if (!found || count == 0) PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
}

/** @brief Global domain bounds, in solver units, from the replicated rank bounding boxes. */
static PetscErrorCode DomainBounds(const SimCtx *simCtx, Cmpnts *lower, Cmpnts *upper)
{
    PetscMPIInt size = 1;

    PetscFunctionBeginUser;
    PetscCheck(simCtx->bboxlist, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Rank bounding boxes are needed before a particle initial value is applied.");
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
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
    if (!plan || end <= first) PetscFunctionReturn(0);
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
    if (!plan) PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ParticleFieldPlanDestroy().
 * @see ParticleFieldPlanDestroy()
 */
PetscErrorCode ParticleFieldPlanDestroy(ParticleFieldPlan **plan)
{
    PetscFunctionBeginUser;
    if (!plan || !*plan) PetscFunctionReturn(0);
    for (PetscInt b = 0; b < (*plan)->count; ++b) {
        for (PetscInt c = 0; c < PARTICLE_FIELD_MAX_COMPONENTS; ++c) {
            PetscCall(PicurvExpressionDestroy(&(*plan)->bindings[b].component[c]));
        }
    }
    PetscCall(PetscFree((*plan)->bindings));
    PetscCall(PetscFree(*plan));
    PetscFunctionReturn(0);
}
