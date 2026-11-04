"""!
@file test_velocity_sources.py
@brief Pytest coverage for prescribed velocity sources: `sources:`, the `prescribed_velocity`
       inlet, and the `source` initial condition.

The CLI resolves a source (inline or named), inlines its constants, applies the target's
constraints, and passes it to the runtime as quoted `-inlet_source_<id>_*` or `-ic_source_*`
options; `src/expression.c` compiles and evaluates it. The language itself is covered by
`test_expression_language.py`.
"""

import importlib.machinery
import importlib.util
from copy import deepcopy
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV_CORE = REPO_ROOT / "picurv_cli" / "core.py"
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "valid"


def load_core():
    """!
    @brief Load the conductor core as a module.
    @return Loaded `picurv_cli/core.py` module.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_velocity_sources_core", str(PICURV_CORE))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


CORE = load_core()

DUCT = {
    "components": "curvilinear",
    "constants": {"U": 1.5},
    "expression": {"u_zeta": "36*U*xi*(1-xi)*eta*(1-eta)"},
}


def load_configs(sources=None, inlet=None, initial_conditions=None):
    """!
    @brief Load the valid fixture configs with optional sources, inlet entry, and IC.
    @param[in] sources Mapping for `sources:`, or None.
    @param[in] inlet Replacement for the `-Zeta` boundary entry, or None.
    @param[in] initial_conditions Replacement for `properties.initial_conditions`, or None.
    @return `(case, solver, monitor)` mappings.
    """
    case = CORE.read_yaml_file(str(FIXTURES / "case.yml"))
    if sources is not None:
        case["sources"] = deepcopy(sources)
    if inlet is not None:
        case["boundary_conditions"][4] = dict({"face": "-Zeta", "type": "INLET"}, **deepcopy(inlet))
    if initial_conditions is not None:
        case["properties"]["initial_conditions"] = deepcopy(initial_conditions)
    solver = CORE.read_yaml_file(str(FIXTURES / "solver.yml"))
    monitor = CORE.read_yaml_file(str(FIXTURES / "monitor.yml"))
    return case, solver, monitor


def validate(case, solver, monitor):
    """!
    @brief Run the configuration validator on in-memory configs.
    @param[in] case Case mapping.
    @param[in] solver Solver mapping.
    @param[in] monitor Monitor mapping.
    """
    CORE.validate_simulation_configs(case, solver, monitor, "case.yml", "solver.yml", "monitor.yml")


def test_a_named_source_feeds_the_inlet_and_the_initial_condition():
    """!
    @brief One named source is resolved for both targets, with its constants inlined.
    """
    case, solver, monitor = load_configs(
        sources={"duct": DUCT},
        inlet={"handler": "prescribed_velocity", "source": "duct",
               "constraints": {"bulk_velocity": "U*(1 + 0.3*sin(2*pi*t))"}},
        initial_conditions={"source": "duct"},
    )
    validate(case, solver, monitor)
    inlet = CORE.validate_and_prepare_boundary_conditions(case)[0][4]
    assert inlet["params"]["source_id"] == 0
    assert inlet["source"]["frame"] == "curvilinear"
    assert inlet["source"]["expressions"][:2] == ["0.0", "0.0"]
    assert "(1.5)" in inlet["source"]["expressions"][2]
    assert inlet["source"]["bulk_velocity"] == "(1.5)*(1 + 0.3*sin(2*pi*t))"

    ic = CORE.resolve_initial_condition_config(
        case["properties"]["initial_conditions"], [], CORE.resolve_reference_scales(case),
        named_sources=case["sources"])
    assert ic["finit"] == 5
    assert ic["source"]["expressions"] == inlet["source"]["expressions"]


def test_sources_become_quoted_runtime_options():
    """!
    @brief Each value is one quoted option; the frame and the bulk target ride alongside.
    """
    source = CORE.apply_source_constraints(
        CORE.resolve_velocity_source({"components": "normal", "expression": "6*eta*(1-eta)"}, {},
                                     "src", allow_normal=True),
        {"scale": 2, "bulk_velocity": 1.0}, "constraints", allowed={"scale", "bulk_velocity"})
    lines = CORE.velocity_source_control_lines("-inlet_source_3", source, "src")
    assert lines == [
        "-inlet_source_3_frame normal",
        '-inlet_source_3_expr_0 "(6*eta*(1-eta))*(2.0)"',
        '-inlet_source_3_bulk_velocity "1.0"',
    ]


def test_inlet_bcs_line_carries_only_the_source_id(tmp_path):
    """!
    @brief The bcs file names the source; its expressions travel in the control file.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    case, _, _ = load_configs(inlet={"handler": "prescribed_velocity", "source": DUCT})
    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)
    generated = CORE.generate_multi_block_bcs(
        str(run_dir), "unit", case, {"Case": "case.yml", "Solver": "solver.yml", "Monitor": "monitor.yml"})
    line = next(text for text in Path(generated[0]).read_text(encoding="utf-8").splitlines()
                if text.startswith("-Zeta"))
    assert "prescribed_velocity" in line
    assert "source_id=0" in line
    assert "xi" not in line


@pytest.mark.parametrize("sources, inlet, initial_conditions, message", [
    (None, {"handler": "prescribed_velocity", "source": "missing"}, None, "no source named 'missing'"),
    (None, {"handler": "prescribed_velocity"}, None, "needs a 'source'"),
    (None, {"handler": "constant_velocity", "params": {"vx": 0, "vy": 0, "vz": 1}, "source": DUCT}, None,
     "only read by the 'prescribed_velocity' handler"),
    (None, {"handler": "prescribed_velocity", "source": {"expression": {"u_xi": "1"}}}, None,
     "not cartesian components"),
    (None, {"handler": "prescribed_velocity", "source": {"components": "normal", "expression": "q*2"}}, None,
     "unsupported expression name 'q'"),
    (None, {"handler": "prescribed_velocity", "source": {"components": "normal", "expression": "1",
                                                         "constants": {"xi": 2}}}, None,
     "cannot name a constant"),
    (None, {"handler": "prescribed_velocity", "source": DUCT, "constraints": {"project": True}}, None,
     "not constraints of this target"),
    (None, None, {"source": {"components": "normal", "expression": "1"}}, "expected one of"),
    ({"duct": DUCT}, None, {"source": "duct", "constraints": {"bulk_velocity": 1.0}},
     "not constraints of this target"),
    ({"duct": DUCT}, None, {"source": "duct", "mode": "generated"}, "takes 'source'"),
])
def test_invalid_sources_are_refused(sources, inlet, initial_conditions, message, capsys):
    """!
    @brief Unknown names, wrong frames, bad expressions and misplaced constraints fail validation.
    @param[in] sources `sources:` mapping.
    @param[in] inlet `-Zeta` entry.
    @param[in] initial_conditions IC mapping.
    @param[in] message Expected fragment of the error.
    @param[in] capsys Pytest output capture.
    """
    case, solver, monitor = load_configs(sources, inlet, initial_conditions)
    with pytest.raises((SystemExit, ValueError)) as caught:
        validate(case, solver, monitor)
        CORE.resolve_initial_condition_config(case["properties"]["initial_conditions"], [],
                                              CORE.resolve_reference_scales(case),
                                              named_sources=case.get("sources"))
    text = capsys.readouterr().err + str(caught.value)
    assert message in text
