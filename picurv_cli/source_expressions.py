#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""!
@file source_expressions.py
@brief The source expression language on the Python side: its functions and variables, validation,
       constant substitution, region lowering and NumPy evaluation.
@details The CLI checks every source with it before a run; `generators/ic.gen` evaluates one with
         it. The solver evaluates the same language in `src/expression.c`, and
         `tests/test_expression_language.py` holds the two to one conformance fixture.
"""

import ast
import copy
import io
import math
import re
import tokenize

import numpy as np

ALLOWED_FUNCTIONS = {
    "abs": np.abs,
    "atan2": np.arctan2,
    "cos": np.cos,
    "cosh": np.cosh,
    "erf": np.vectorize(math.erf, otypes=[float]),
    "exp": np.exp,
    "floor": np.floor,
    "log": np.log,
    "maximum": np.maximum,
    "minimum": np.minimum,
    "sin": np.sin,
    "sinh": np.sinh,
    "sqrt": np.sqrt,
    "tan": np.tan,
    "tanh": np.tanh,
    "where": np.where,
}
#: Arguments each function takes; the runtime evaluator (`src/expression.c`) accepts exactly
#: the same calls.
FUNCTION_ARITY = {"abs": 1, "cos": 1, "cosh": 1, "erf": 1, "exp": 1, "floor": 1, "log": 1, "sin": 1,
                  "sinh": 1, "sqrt": 1, "tan": 1, "tanh": 1, "atan2": 2, "maximum": 2, "minimum": 2,
                  "where": 3, "besselj": 2, "bessely": 2, "ber": 1, "bei": 1, "sum": 4}
#: Bessel functions of integer order (`besselj(n, x)`, `bessely(n, x)`) and the Kelvin
#: functions `ber`, `bei`: validated here, evaluated only by the runtime, which has them.
SOLVER_ONLY_FUNCTIONS = ("besselj", "bessely", "ber", "bei")
#: Terms one series `sum(k, lo, hi, body)` may add; the runtime holds the same limit.
MAX_SERIES_TERMS = 1000000
#: Series that may nest inside one another; the runtime holds the same limit.
MAX_SERIES_DEPTH = 4
#: Per-particle random draws: evaluated only by the runtime, from a hash of the particle's
#: identity; an optional non-negative integer literal names a stream shared across fields.
RANDOM_FUNCTIONS = ("uniform", "normal")
#: Variables an Eulerian expression may use.
EULERIAN_EXPRESSION_NAMES = frozenset({"x", "y", "z", "xi", "eta", "zeta", "i", "j", "k"})
#: Variables a particle initial-value expression may use: physical coordinates, domain-
#: normalized coordinates, particle identity, and physical time.
PARTICLE_EXPRESSION_NAMES = frozenset({"x", "y", "z", "t", "xi", "eta", "zeta", "L_xi", "L_eta", "L_zeta", "pid"})
_DECIMAL_LITERAL = re.compile(r"^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
ALLOWED_AST_NODES = (
    ast.Expression,
    ast.BinOp,
    ast.UnaryOp,
    ast.Call,
    ast.Name,
    ast.Load,
    ast.Constant,
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.Pow,
    ast.Mod,
    ast.Compare,
    ast.Eq,
    ast.NotEq,
    ast.Lt,
    ast.LtE,
    ast.Gt,
    ast.GtE,
    ast.USub,
    ast.UAdd,
    ast.BoolOp,
    ast.And,
    ast.Or,
    ast.Not,
)


def validate_expression(expression: str, names=EULERIAN_EXPRESSION_NAMES, random=False) -> None:
    """!
    @brief Reject anything outside the expression language the runtime also evaluates.
    @details Python expression syntax restricted to decimal numbers, the given names, `pi`,
             arithmetic, comparisons (chains included), `and`/`or`/`not`, and the
             documented functions at their fixed arities. `uniform(n)`/`normal(n)` are
             accepted only when @p random is true, with an optional stream `n` that must be
             a non-negative integer literal.
    @param[in] expression User expression.
    @param[in] names Variables in scope, besides `pi`.
    @param[in] random Whether per-particle random draws may be called.
    @throws ValueError naming the first construct outside the language.
    """
    if not isinstance(expression, str) or not expression.strip():
        raise ValueError("an expression must be a non-empty string.")
    if "\n" in expression or "\r" in expression:
        raise ValueError("an expression must fit on one line.")
    try:
        tree = ast.parse(expression, mode="eval")
        tokens = list(tokenize.generate_tokens(io.StringIO(expression).readline))
    except (SyntaxError, tokenize.TokenError) as exc:
        raise ValueError(f"expression '{expression}' is not valid syntax: {exc}.") from exc
    for token in tokens:
        if token.type == tokenize.NUMBER and not _DECIMAL_LITERAL.match(token.string):
            raise ValueError(f"'{token.string}' is not a decimal number; write numbers like 1, 0.5, or 2e-3.")
    _check_expression_node(tree.body, frozenset(names) | {"pi"}, (), random)


def _series_bound(node) -> int:
    """!
    @brief The integer a series bound denotes: an integer literal, possibly negated.
    @param[in] node Bound node.
    @return The bound.
    @throws ValueError for anything else.
    """
    sign = 1
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        sign, node = (-1 if isinstance(node.op, ast.USub) else 1), node.operand
    if (not isinstance(node, ast.Constant) or type(node.value) not in (int, float)
            or node.value != int(node.value)):
        raise ValueError("series bounds must be integer numbers.")
    return sign * int(node.value)


def _check_expression_node(node, names, indices, random) -> None:
    """!
    @brief Check one node of a parsed expression, with the series indices enclosing it in scope.
    @param[in] node AST node.
    @param[in] names Variables in scope, `pi` included.
    @param[in] indices Indices of the enclosing series, outermost first.
    @param[in] random Whether per-particle random draws may be called.
    @throws ValueError naming the first construct outside the language.
    """
    if not isinstance(node, ALLOWED_AST_NODES):
        raise ValueError(f"unsupported expression syntax: {type(node).__name__}.")
    if isinstance(node, ast.Constant) and type(node.value) not in (int, float):
        raise ValueError(f"only numbers may appear as constants, not {node.value!r}.")
    if isinstance(node, ast.Name):
        if node.id not in names and node.id not in indices:
            raise ValueError(f"unsupported expression name '{node.id}'.")
        return
    if isinstance(node, ast.Call):
        name = node.func.id if isinstance(node.func, ast.Name) else None
        if node.keywords:
            raise ValueError("functions take positional arguments only.")
        if name in RANDOM_FUNCTIONS and random:
            if len(node.args) > 1:
                raise ValueError(f"{name}() takes at most one stream argument.")
            if node.args:
                stream = node.args[0]
                if (not isinstance(stream, ast.Constant) or type(stream.value) not in (int, float)
                        or stream.value < 0 or stream.value != int(stream.value)):
                    raise ValueError(f"the stream of {name}() must be a non-negative integer literal.")
            return
        if name not in FUNCTION_ARITY:
            raise ValueError(f"'{name}' is not an available function.")
        if len(node.args) != FUNCTION_ARITY[name]:
            raise ValueError(f"{name}() takes {FUNCTION_ARITY[name]} argument(s).")
        if name == "sum":
            index, lower, upper, body = node.args
            if not isinstance(index, ast.Name):
                raise ValueError("a series sum(k, lo, hi, body) needs an index name first.")
            if index.id in names or index.id in indices:
                raise ValueError(f"the series index '{index.id}' reuses a name already in scope.")
            if len(indices) == MAX_SERIES_DEPTH:
                raise ValueError(f"series may nest at most {MAX_SERIES_DEPTH} deep.")
            if _series_bound(upper) - _series_bound(lower) + 1 > MAX_SERIES_TERMS:
                raise ValueError(f"a series may add at most {MAX_SERIES_TERMS} terms.")
            _check_expression_node(body, names, indices + (index.id,), random)
            return
        for argument in node.args:
            _check_expression_node(argument, names, indices, random)
        return
    for child in ast.iter_child_nodes(node):
        _check_expression_node(child, names, indices, random)


class _LogicalRewriter(ast.NodeTransformer):
    """!
    @brief Rewrite logical operators and chained comparisons into array operations.
    @details `and`, `or`, `not` and a chain `a < b < c` have no elementwise meaning for
             arrays, so they become `logical_and`, `logical_or`, `logical_not` calls, the
             same 1-or-0 semantics the runtime evaluator gives them.
    """

    @staticmethod
    def _call(name, args):
        """!
        @brief Build a call to one of the evaluation-only logical helpers.
        @param[in] name Helper name.
        @param[in] args Argument nodes.
        @return Call node.
        """
        return ast.Call(func=ast.Name(id=name, ctx=ast.Load()), args=list(args), keywords=[])

    def visit_Call(self, node):
        """!
        @brief Rewrite `sum(k, lo, hi, body)` as a loop over `lambda k: body`.
        @param[in] node Call node.
        @return Rewritten node.
        """
        self.generic_visit(node)
        if not (isinstance(node.func, ast.Name) and node.func.id == "sum"):
            return node
        index, lower, upper, body = node.args
        loop = ast.Lambda(args=ast.arguments(posonlyargs=[], args=[ast.arg(arg=index.id)], vararg=None,
                                             kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[]),
                          body=body)
        return self._call("series", [loop, lower, upper])

    def visit_BoolOp(self, node):
        """!
        @brief Fold `a and b and c` into nested logical calls.
        @param[in] node BoolOp node.
        @return Rewritten node.
        """
        self.generic_visit(node)
        name = "logical_and" if isinstance(node.op, ast.And) else "logical_or"
        result = node.values[0]
        for value in node.values[1:]:
            result = self._call(name, [result, value])
        return result

    def visit_UnaryOp(self, node):
        """!
        @brief Rewrite `not a` as `logical_not(a)`.
        @param[in] node UnaryOp node.
        @return Rewritten node.
        """
        self.generic_visit(node)
        if isinstance(node.op, ast.Not):
            return self._call("logical_not", [node.operand])
        return node

    def visit_Compare(self, node):
        """!
        @brief Rewrite `a < b < c` as `logical_and(a < b, b < c)`.
        @param[in] node Compare node.
        @return Rewritten node.
        """
        self.generic_visit(node)
        if len(node.ops) == 1:
            return node
        operands = [node.left] + node.comparators
        pairs = [ast.Compare(left=copy.deepcopy(operands[i]), ops=[op], comparators=[copy.deepcopy(operands[i + 1])])
                 for i, op in enumerate(node.ops)]
        result = pairs[0]
        for pair in pairs[1:]:
            result = self._call("logical_and", [result, pair])
        return result


def _solver_only(name: str):
    """!
    @brief A stand-in for a function only the runtime evaluates.
    @param[in] name Function name.
    @return Callable that refuses evaluation.
    """
    def refuse(*arguments):
        """!
        @brief Refuse to evaluate, naming the function and where it can be used.
        @param[in] arguments The call's arguments, unused.
        @throws ValueError always.
        """
        del arguments
        raise ValueError(f"{name}() is evaluated only by the solver; use it in a velocity source, "
                         "not in an expression picurv evaluates before launch.")
    return refuse


def evaluate_values(expression: str, values: dict):
    """!
    @brief Evaluate an already validated expression over named arrays or scalars.
    @param[in] expression Expression text.
    @param[in] values Mapping of variable name to value or array.
    @return Result as a float array, not yet checked for finiteness.
    """
    tree = ast.fix_missing_locations(_LogicalRewriter().visit(ast.parse(expression, mode="eval")))
    env = dict(ALLOWED_FUNCTIONS)
    env.update({"logical_and": np.logical_and, "logical_or": np.logical_or,
                "logical_not": np.logical_not, "pi": math.pi,
                "series": lambda body, lo, hi: sum((body(k) for k in range(int(lo), int(hi) + 1)), 0.0)})
    for name in SOLVER_ONLY_FUNCTIONS:
        env[name] = _solver_only(name)
    env.update(values)
    with np.errstate(all="ignore"):
        # Globals, not locals: a series body is a lambda, which sees only the globals.
        value = eval(compile(tree, "<ic-expression>", "eval"), dict(env, __builtins__={}))
    return np.asarray(value, dtype=float)


#: Region shapes a structured field value may paint, later regions over earlier ones.
REGION_SHAPES = ("half_space", "slab", "box", "ball", "cylinder")
#: Profiles of a smoothed region edge, each rising from 0 outside to 1 inside across `width`.
EDGE_PROFILES = {
    "linear": "({s})",
    "smoothstep": "({s})*({s})*(3 - 2*({s}))",
    "cosine": "0.5*(1 - cos(pi*({s})))",
}
_REGION_KEYS = {
    "half_space": {"axis", "at", "side"},
    "slab": {"axis", "from", "to"},
    "box": {"min", "max"},
    "ball": {"center", "radius"},
    "cylinder": {"axis", "center", "radius"},
}
_AXES = ("x", "y", "z")


def _number_text(value, params: dict, where: str) -> str:
    """!
    @brief Render one numeric region parameter, which may name a declared parameter.
    @param[in] value Number, or the name of a parameter.
    @param[in] params Declared parameters.
    @param[in] where Configuration path, for diagnostics.
    @return Parenthesized numeric text.
    """
    if isinstance(value, str) and value in params:
        return repr(float(params[value]))
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{where} must be a finite number or a declared parameter name.")
    return repr(float(value))


def _region_distance(region: dict, params: dict, where: str) -> str:
    """!
    @brief Signed distance to a region's boundary, negative inside, as expression text.
    @param[in] region Region mapping.
    @param[in] params Declared parameters.
    @param[in] where Configuration path, for diagnostics.
    @return Expression text for the signed distance.
    """
    shape = region.get("shape")
    if shape not in REGION_SHAPES:
        raise ValueError(f"{where}.shape must be one of {', '.join(REGION_SHAPES)}.")
    extra = set(region) - _REGION_KEYS[shape] - {"shape", "value", "edge"}
    missing = _REGION_KEYS[shape] - set(region) - ({"side"} if shape == "half_space" else set())
    if extra or missing:
        raise ValueError(f"{where} ({shape}) takes {sorted(_REGION_KEYS[shape])}, plus value and edge.")

    def num(key, value=None):
        """!
        @brief Render a region number (a literal or a `params` name) as expression text.
        @param[in] key Region key, for diagnostics and as the default source of the number.
        @param[in] value Number to render in place of `region[key]`, if given.
        @return Expression text of the number.
        """
        return _number_text(region[key] if value is None else value, params, f"{where}.{key}")

    def axis_of(key="axis"):
        """!
        @brief Return the region's axis name, refusing anything but x, y, z.
        @param[in] key Region key holding the axis.
        @return Axis name.
        """
        axis = region.get(key)
        if axis not in _AXES:
            raise ValueError(f"{where}.{key} must be x, y, or z.")
        return axis

    def triple(key):
        """!
        @brief Render a three-coordinate region entry as a list of expression texts.
        @param[in] key Region key holding the three coordinates.
        @return List of three expression texts.
        """
        values = region.get(key)
        if not isinstance(values, list) or len(values) != 3:
            raise ValueError(f"{where}.{key} must be a list of three coordinates.")
        return [num(key, v) for v in values]

    if shape == "half_space":
        axis, side = axis_of(), region.get("side", "above")
        if side not in ("above", "below"):
            raise ValueError(f"{where}.side must be above or below.")
        return f"({num('at')} - {axis})" if side == "above" else f"({axis} - {num('at')})"
    if shape == "slab":
        axis = axis_of()
        return f"maximum({num('from')} - {axis}, {axis} - {num('to')})"
    if shape == "box":
        lower, upper = triple("min"), triple("max")
        terms = [f"maximum({lo} - {a}, {a} - {hi})" for a, lo, hi in zip(_AXES, lower, upper)]
        return f"maximum(maximum({terms[0]}, {terms[1]}), {terms[2]})"
    if shape == "ball":
        center = triple("center")
        return (f"sqrt(({_AXES[0]} - {center[0]})**2 + ({_AXES[1]} - {center[1]})**2 + "
                f"({_AXES[2]} - {center[2]})**2) - {num('radius')}")
    axis = axis_of()
    others = [a for a in _AXES if a != axis]
    center = region.get("center")
    if not isinstance(center, list) or len(center) != 2:
        raise ValueError(f"{where}.center must list the two coordinates across the axis ({', '.join(others)}).")
    c0, c1 = (num("center", v) for v in center)
    return f"sqrt(({others[0]} - {c0})**2 + ({others[1]} - {c1})**2) - {num('radius')}"


def _substitute_params(expression: str, params: dict) -> str:
    """!
    @brief Replace declared parameter names in an expression by their values.
    @param[in] expression Expression text.
    @param[in] params Declared parameters.
    @return Expression text with parameters inlined.
    """
    pieces, cursor = [], 0
    for token in tokenize.generate_tokens(io.StringIO(expression).readline):
        if token.type == tokenize.NAME and token.string in params:
            start, end = token.start[1], token.end[1]
            pieces.append(expression[cursor:start])
            pieces.append(f"({float(params[token.string])!r})")
            cursor = end
    pieces.append(expression[cursor:])
    return "".join(pieces)


def _validate_at(expression: str, names, random: bool, where: str):
    """!
    @brief Validate an expression written in a configuration value, naming where it is.
    @param[in] expression Expression text.
    @param[in] names Variables the expression may use.
    @param[in] random Whether random draws are allowed.
    @param[in] where Configuration path of the value.
    @throws ValueError naming `where` if the expression is outside the language.
    """
    try:
        validate_expression(expression, names, random)
    except ValueError as e:
        raise ValueError(f"{where}: {e}") from None


def lower_value(spec, names=PARTICLE_EXPRESSION_NAMES, random=True, where="value") -> str:
    """!
    @brief Turn a field value - number, expression, or structured regions - into one expression.
    @details A mapping may declare `params` (named numbers usable in its expressions and
             region parameters) and gives either `value` or `background` with `regions`.
             Regions paint in order, a later region over an earlier one: a sharp region by
             `where(d <= 0, v, before)`, a smoothed one by `before*(1 - s) + v*s`, with `d`
             the signed distance to the region and `s` its edge profile over `edge.width`.
             The result is validated like any expression.
    @param[in] spec Number, expression string, or mapping.
    @param[in] names Variables in scope.
    @param[in] random Whether random draws may be called.
    @param[in] where Configuration path, for diagnostics.
    @return Validated expression text.
    """
    if isinstance(spec, bool):
        raise ValueError(f"{where} must be a number, an expression, or a mapping.")
    if isinstance(spec, (int, float)):
        if not math.isfinite(spec):
            raise ValueError(f"{where} must be finite.")
        return repr(float(spec))
    if isinstance(spec, str):
        _validate_at(spec, names, random, where)
        return spec
    if not isinstance(spec, dict):
        raise ValueError(f"{where} must be a number, an expression, or a mapping.")
    unknown = set(spec) - {"params", "value", "background", "regions"}
    if unknown:
        raise ValueError(f"{where} does not take {sorted(unknown)}; it takes params, value, background, regions.")
    params = spec.get("params") or {}
    if not isinstance(params, dict):
        raise ValueError(f"{where}.params must be a mapping of names to numbers.")
    for name, number in params.items():
        reserved = set(names) | {"pi"} | set(FUNCTION_ARITY) | set(RANDOM_FUNCTIONS) | {"and", "or", "not"}
        if not isinstance(name, str) or not name.isidentifier() or name in reserved:
            raise ValueError(f"{where}.params name '{name}' is not a free identifier.")
        _number_text(number, {}, f"{where}.params.{name}")
    scope = set(names) | set(params)

    def expression_of(value, path):
        """!
        @brief Lower one value (number, `params` name, or expression) with parameters inlined.
        @param[in] value Configured value.
        @param[in] path Configuration path, for diagnostics.
        @return Expression text.
        """
        if isinstance(value, str) and value not in params:
            _validate_at(value, scope, random, path)
            return f"({_substitute_params(value, params)})"
        return _number_text(value, params, path)

    if "value" in spec:
        if "background" in spec or "regions" in spec:
            raise ValueError(f"{where} gives either value, or background with regions.")
        lowered = expression_of(spec["value"], f"{where}.value")
    else:
        regions = spec.get("regions")
        if not isinstance(regions, list) or not regions:
            raise ValueError(f"{where}.regions must be a non-empty list.")
        lowered = expression_of(spec.get("background", 0.0), f"{where}.background")
        for index, region in enumerate(regions):
            path = f"{where}.regions[{index}]"
            if not isinstance(region, dict) or "value" not in region:
                raise ValueError(f"{path} must be a mapping with a shape and a value.")
            distance = _region_distance(region, params, path)
            value = expression_of(region["value"], f"{path}.value")
            edge = region.get("edge")
            if edge is None:
                lowered = f"where({distance} <= 0, {value}, {lowered})"
                continue
            if not isinstance(edge, dict) or set(edge) != {"profile", "width"} or edge["profile"] not in EDGE_PROFILES:
                raise ValueError(f"{path}.edge must give profile ({', '.join(EDGE_PROFILES)}) and width.")
            width = _number_text(edge["width"], params, f"{path}.edge.width")
            if float(width) <= 0.0:
                raise ValueError(f"{path}.edge.width must be positive.")
            ramp = f"maximum(0, minimum(1, 0.5 - ({distance})/{width}))"
            blend = EDGE_PROFILES[edge["profile"]].format(s=ramp)
            lowered = f"({lowered})*(1 - {blend}) + {value}*{blend}"
    validate_expression(lowered, names, random)
    return lowered
