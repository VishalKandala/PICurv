"""!
@file test_expression_language.py
@brief Pytest coverage for the expression language shared by Eulerian and particle initial values.

`generators/ic.gen` validates and evaluates the language in Python; the runtime evaluates
it in C (`src/ParticleInitialConditions.c`). Both read
`tests/fixtures/expression_conformance.txt` and must reproduce every case.
"""

import importlib.machinery
import importlib.util
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
IC_GENERATOR = REPO_ROOT / "generators" / "ic.gen"
CONFORMANCE = REPO_ROOT / "tests" / "fixtures" / "expression_conformance.txt"


def load_ic_gen():
    """!
    @brief Load `generators/ic.gen` as a module.
    @return Loaded module.
    """
    loader = importlib.machinery.SourceFileLoader("expression_language_ic_gen", str(IC_GENERATOR))
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
