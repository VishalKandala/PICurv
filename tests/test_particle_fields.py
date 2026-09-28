"""!
@file test_particle_fields.py
@brief Pytest coverage for configured particle initial values (`models.physics.particles.fields`).

The CLI lowers each value to one expression per component with `generators/ic.gen` and
passes it to the runtime as `-particle_fields_*` options; `src/ParticleInitialConditions.c`
compiles and applies it. The language itself is covered by `test_expression_language.py`.
"""

import importlib.machinery
import importlib.util
import re
from copy import deepcopy
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV_CORE = REPO_ROOT / "picurv_cli" / "core.py"
PARTICLE_CATALOG = REPO_ROOT / "src" / "particle_field_catalog.c"
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "valid"


def load_core():
    """!
    @brief Load the conductor core as a module.
    @return Loaded `picurv_cli/core.py` module.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_particle_fields_core", str(PICURV_CORE))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


CORE = load_core()


def load_configs(particles):
    """!
    @brief Load the valid fixture configs with the given particle block.
    @param[in] particles Mapping for `models.physics.particles`.
    @return `(case, solver, monitor)` mappings.
    """
    case = CORE.read_yaml_file(str(FIXTURES / "case.yml"))
    case.setdefault("models", {}).setdefault("physics", {})["particles"] = deepcopy(particles)
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


def test_settable_fields_match_the_runtime_catalog():
    """!
    @brief The CLI accepts exactly the fields whose catalog entry allows a user initial value.
    """
    source = PARTICLE_CATALOG.read_text(encoding="utf-8")
    settable = {}
    for entry in source.split("PARTICLE_FIELD_ENTRY(")[1:]:
        if "PARTICLE_FIELD_CAPABILITY_USER_INITIALIZE" in entry:
            match = re.search(r'"(\w+)"[^\n]*\n\s*(\d+),', entry)
            settable[match.group(1)] = int(match.group(2))
    assert settable == CORE.PARTICLE_SETTABLE_FIELDS


def test_values_become_quoted_runtime_options():
    """!
    @brief A structured value is lowered once and emitted as one quoted option per component.
    """
    case, _, _ = load_configs({
        "count": 100,
        "fields": {"Psi": {"background": 0, "regions": [
            {"shape": "half_space", "axis": "x", "at": 0.5, "side": "above", "value": 1},
        ]}},
    })
    lines = []
    CORE.parse_and_add_model_flags(case, lines)
    assert "-particle_fields_count 1" in lines
    assert "-particle_fields_0_name Psi" in lines
    expression = next(line for line in lines if line.startswith("-particle_fields_0_expr_0 "))
    text = expression.split(" ", 1)[1]
    assert text.startswith('"') and text.endswith('"')
    assert "where(" in text


def test_no_fields_emit_no_options():
    """!
    @brief Cases without `fields` keep their control file unchanged.
    """
    case, _, _ = load_configs({"count": 100})
    lines = []
    CORE.parse_and_add_model_flags(case, lines)
    assert not any(line.startswith("-particle_fields") for line in lines)


def test_a_valid_value_passes_validation():
    """!
    @brief A number, an expression, and a region block are all accepted.
    """
    for value in (0.5, "0.5 + 0.5*sin(2*pi*x)", {"params": {"k": 2}, "value": "k*xn"}):
        validate(*load_configs({"count": 100, "fields": {"Psi": value}}))


@pytest.mark.parametrize("particles, message", [
    ({"count": 100, "fields": {"Velocity": 1.0}}, "not a field a case can set"),
    ({"count": 100, "fields": {"Psi": "foo + 1"}}, "fields.Psi"),
    ({"count": 100, "fields": {}}, "non-empty mapping"),
    ({"count": 0, "fields": {"Psi": 1.0}}, "needs particles"),
])
def test_invalid_fields_are_refused(particles, message, capsys):
    """!
    @brief Unknown or re-derived fields, bad expressions, and particle-free runs fail validation.
    @param[in] particles Particle block.
    @param[in] message Expected fragment of the error.
    @param[in] capsys Pytest output capture.
    """
    with pytest.raises(SystemExit):
        validate(*load_configs(particles))
    assert message in capsys.readouterr().err


def test_a_prescribed_scalar_source_conflicts_with_an_initial_value(capsys):
    """!
    @brief A verification source that prescribes Psi every step would overwrite the value.
    @param[in] capsys Pytest output capture.
    """
    case, solver, monitor = load_configs({"count": 100, "fields": {"Psi": 1.0}})
    solver.setdefault("verification", {})["sources"] = {"scalar": {"mode": "analytical"}}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)
    assert "verification.sources.scalar" in capsys.readouterr().err


def test_a_point_source_without_draws_warns(capsys):
    """!
    @brief Particles released from one point all get the same value unless the value draws.
    @param[in] capsys Pytest output capture.
    """
    point = {"x": 0.5, "y": 0.5, "z": 0.5}
    validate(*load_configs({"count": 100, "init_mode": "PointSource", "point_source": point,
                            "fields": {"Psi": "x"}}))
    assert "has no random draw" in capsys.readouterr().err
    validate(*load_configs({"count": 100, "init_mode": "PointSource", "point_source": point,
                            "fields": {"Psi": "where(uniform() < 0.5, 0, 1)"}}))
    assert "has no random draw" not in capsys.readouterr().err
