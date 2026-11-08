"""!
@file test_expression_language.py
@brief Pytest coverage for the expression language shared by Eulerian and particle initial values.

`picurv_cli/source_expressions.py` validates and evaluates the language in Python; the runtime evaluates
it in C (`src/expression.c`). Both read
`tests/fixtures/expression_conformance.txt` and must reproduce every case.
"""

import importlib.machinery
import importlib.util
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
LANGUAGE = REPO_ROOT / "picurv_cli" / "source_expressions.py"
CONFORMANCE = REPO_ROOT / "tests" / "fixtures" / "expression_conformance.txt"


def load_ic_gen():
    """!
    @brief Load `picurv_cli/source_expressions.py` as a module.
    @return Loaded module.
    """
    loader = importlib.machinery.SourceFileLoader("expression_language_source_expressions", str(LANGUAGE))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


IC_GEN = load_ic_gen()


def conformance_cases():
    """!
    @brief Read the shared conformance cases.
    @return List of `(expression, x, y, z, expected)` tuples.
    """
    cases = []
    for line in CONFORMANCE.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        expression, x, y, z, expected = (part.strip() for part in line.split(";"))
        cases.append((expression, float(x), float(y), float(z), float(expected)))
    return cases


@pytest.mark.parametrize("expression, x, y, z, expected", conformance_cases())
def test_conformance_case(expression, x, y, z, expected):
    """!
    @brief The Python evaluator reproduces a shared conformance case.
    @param[in] expression Expression text.
    @param[in] x Value of x.
    @param[in] y Value of y.
    @param[in] z Value of z.
    @param[in] expected Expected value.
    """
    IC_GEN.validate_expression(expression, {"x", "y", "z"})
    value = float(IC_GEN.evaluate_values(expression, {"x": x, "y": y, "z": z}))
    assert value == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_logical_operators_act_elementwise_on_arrays():
    """!
    @brief `and`, `not`, and chained comparisons work on arrays, as a grid needs.
    """
    x = np.array([0.1, 0.3, 0.7])
    assert list(IC_GEN.evaluate_values("0.2 < x < 0.5", {"x": x})) == [0.0, 1.0, 0.0]
    assert list(IC_GEN.evaluate_values("not (x > 0.2 and x < 0.5)", {"x": x})) == [1.0, 0.0, 1.0]


@pytest.mark.parametrize("expression", [
    "0x1F", "1_000", "1j", "x // 2", "True", "'a'", "abs", "abs(1, 2)", "where(x, 1)",
    "foo", "sin(x=1)", "uniform()", "x\n+1", "lambda: 1", "[x]",
    "sum(n, 0, 1.5, n)", "sum(n, 0, x, n)", "sum(x, 0, 2, x)", "sum(n, 0, 2)",
    "sum(n, 0, 2, sum(n, 0, 1, n))", "sum(n, 0, 2, sum(m, 0, n, m))", "sum(n, 0, 2, n) + n",
    "sum(n, 0, 2000000, n)", "sum(pi, 0, 2, 1)", "besselj(1)", "ber(1, 2)",
])
def test_constructs_outside_the_language_are_refused(expression):
    """!
    @brief Anything the runtime evaluator would not accept is refused before launch.
    @param[in] expression Expression outside the language.
    """
    with pytest.raises(ValueError):
        IC_GEN.validate_expression(expression, {"x"})


def test_random_draws_are_particle_only_and_take_literal_streams():
    """!
    @brief Draws need the particle scope, and a stream must be a non-negative integer literal.
    """
    names = IC_GEN.PARTICLE_EXPRESSION_NAMES
    IC_GEN.validate_expression("where(uniform(3) < 0.3, 1, 0) + normal()", names, random=True)
    for bad in ("uniform(x)", "uniform(-1)", "uniform(1.5)", "normal(1, 2)"):
        with pytest.raises(ValueError):
            IC_GEN.validate_expression(bad, names, random=True)
    with pytest.raises(ValueError):
        IC_GEN.validate_expression("uniform()", names, random=False)


def test_lowered_regions_paint_later_regions_over_earlier_ones():
    """!
    @brief A slab then a ball: the ball's value wins where they overlap.
    """
    lowered = IC_GEN.lower_value({
        "params": {"r": 0.1},
        "background": 0,
        "regions": [
            {"shape": "slab", "axis": "y", "from": 0.4, "to": 0.6, "value": 1},
            {"shape": "ball", "center": [0.5, 0.5, 0.5], "radius": "r", "value": 2},
        ],
    }, names={"x", "y", "z"})
    points = {"x": np.array([0.0, 0.0, 0.5]), "y": np.array([0.0, 0.5, 0.5]), "z": np.array([0.0, 0.0, 0.5])}
    assert list(IC_GEN.evaluate_values(lowered, points)) == [0.0, 1.0, 2.0]


def test_a_smoothed_edge_is_half_way_on_the_boundary():
    """!
    @brief Each edge profile gives 0 outside, 1 inside, and one half on the boundary.
    """
    for profile in IC_GEN.EDGE_PROFILES:
        lowered = IC_GEN.lower_value({
            "background": 0,
            "regions": [{"shape": "half_space", "axis": "x", "at": 0.5, "side": "above", "value": 1,
                         "edge": {"profile": profile, "width": 0.1}}],
        }, names={"x", "y", "z"})
        values = IC_GEN.evaluate_values(lowered, {"x": np.array([0.0, 0.5, 1.0]), "y": 0.0, "z": 0.0})
        assert values == pytest.approx([0.0, 0.5, 1.0]), profile


def test_parameters_are_inlined_into_expressions():
    """!
    @brief A declared parameter may be used in an expression value.
    """
    lowered = IC_GEN.lower_value({"params": {"k": 2.0}, "value": "k * x"}, names={"x"})
    assert float(IC_GEN.evaluate_values(lowered, {"x": 3.0})) == 6.0


@pytest.mark.parametrize("spec", [
    {"value": 1, "regions": []},
    {"regions": [{"shape": "cube", "value": 1}]},
    {"regions": [{"shape": "ball", "center": [0, 0], "radius": 1, "value": 1}]},
    {"regions": [{"shape": "slab", "axis": "w", "from": 0, "to": 1, "value": 1}]},
    {"params": {"x": 1.0}, "value": "x"},
    {"regions": [{"shape": "half_space", "axis": "x", "at": 0, "value": 1,
                  "edge": {"profile": "tanh", "width": 0.1}}]},
    {"extra": 1, "value": 1},
    True,
])
def test_malformed_structured_values_are_refused(spec):
    """!
    @brief Unknown keys and shapes, bad coordinates, shadowing parameters, and bad edges fail.
    @param[in] spec Malformed structured value.
    """
    with pytest.raises(ValueError):
        IC_GEN.lower_value(spec, names={"x", "y", "z"})


def test_bessel_and_kelvin_functions_are_validated_but_left_to_the_solver():
    """!
    @brief `besselj`, `bessely`, `ber` and `bei` pass validation, and picurv refuses to
           evaluate them itself, since only the runtime has them.
    """
    expression = "besselj(0, x) + bessely(1, x) + ber(x) + bei(x)"
    IC_GEN.validate_expression(expression, {"x"})
    with pytest.raises(ValueError, match="evaluated only by the solver"):
        IC_GEN.evaluate_values(expression, {"x": 1.0})


def test_a_series_acts_elementwise_on_arrays():
    """!
    @brief A series evaluates over a grid of points, each term broadcast like any expression.
    """
    x = np.array([0.0, 0.5, 1.0])
    assert list(IC_GEN.evaluate_values("sum(n, 0, 2, x**n)", {"x": x})) == [1.0, 1.75, 3.0]


#: The rectangular-duct Poiseuille source documented on page 33, in the arc-length fractions of a
#: section of aspect ratio `beta` = height (eta) / width (xi), normalized to bulk velocity U. The
#: ratio cosh(a)/cosh(b) is written so that it cannot overflow.
DUCT_POISEUILLE = (
    "U*sum(k, 0, 99, (-1)**k/(2*k+1)**3*(1 - (exp(abs((2*k+1)*pi*beta*(eta - 0.5)) - (2*k+1)*pi*beta/2)"
    " + exp(-abs((2*k+1)*pi*beta*(eta - 0.5)) - (2*k+1)*pi*beta/2))/(1 + exp(-(2*k+1)*pi*beta)))"
    "*cos((2*k+1)*pi*(xi - 0.5))) / sum(k, 0, 99, 2/(pi*(2*k+1)**4)*(1 - 2*tanh((2*k+1)*pi*beta/2)"
    "/((2*k+1)*pi*beta)))")


def test_the_documented_duct_series_is_the_duct_solution():
    """!
    @brief On a square section the cosh-cosine series equals `profile.gen`'s double sine series,
           both at unit bulk velocity, and its own normalization gives a unit area mean.
    """
    xi, eta = np.meshgrid((np.arange(64) + 0.5) / 64, (np.arange(64) + 0.5) / 64)
    expression = IC_GEN.substitute_params(DUCT_POISEUILLE, {"U": 1.0, "beta": 1.0})
    IC_GEN.validate_expression(expression, {"xi", "eta"})
    series = IC_GEN.evaluate_values(expression, {"xi": xi, "eta": eta})
    # An independent form of the same solution: the double sine series, normalized to its own mean.
    shape = np.zeros_like(xi)
    mean = 0.0
    for m in range(1, 601, 2):
        for n in range(1, 601, 2):
            coefficient = 1.0 / (m * n * (m * m + n * n))
            shape += coefficient * np.sin(m * np.pi * xi) * np.sin(n * np.pi * eta)
            mean += coefficient * 4.0 / (np.pi * np.pi * m * n)
    assert np.max(np.abs(series - shape / mean)) < 1e-5
    for beta in (0.25, 8.0):  # wide and tall sections: finite, with the requested mean
        other = IC_GEN.evaluate_values(IC_GEN.substitute_params(DUCT_POISEUILLE, {"U": 1.0, "beta": beta}),
                                       {"xi": xi, "eta": eta})
        assert np.all(np.isfinite(other)) and float(np.mean(other)) == pytest.approx(1.0, abs=5e-3)
