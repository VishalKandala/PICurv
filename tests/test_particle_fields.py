"""!
@file test_particle_fields.py
@brief Pytest coverage for configured particle initial values (`models.physics.particles.fields`).

The CLI lowers each value to one expression per component with `picurv_cli/source_expressions.py` and
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
        "particles_per_cell": 1,
        "fields": {"Psi": {"expression": {"background": 0, "regions": [
            {"shape": "half_space", "axis": "x", "at": 0.5, "side": "above", "value": 1},
        ]}}},
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
    case, _, _ = load_configs({"particles_per_cell": 1})
    lines = []
    CORE.parse_and_add_model_flags(case, lines)
    assert not any(line.startswith("-particle_fields") for line in lines)


def test_a_valid_value_passes_validation():
    """!
    @brief A number, an expression with constants, a grid-line variable, and a sum are all accepted.
    """
    for value in ({"expression": 0.5}, {"expression": "0.5 + 0.5*sin(2*pi*x)"},
                  {"expression": "k*xi + L_eta", "constants": {"k": 2}},
                  {"sum": [{"expression": 0.5}, {"expression": "where(uniform() < 0.5, 0, pid)"}]}):
        validate(*load_configs({"particles_per_cell": 1, "fields": {"Psi": value}}))


def test_a_named_source_serves_a_particle_field():
    """!
    @brief A field names an entry of `sources:`, whose constants build on the case's and the fluid's.
    """
    case, _, _ = load_configs({"particles_per_cell": 1, "fields": {"Psi": "mixture"}})
    case["constants"] = {"A": 3.0}
    case["sources"] = {"mixture": {"expression": "A*y*B", "constants": {"B": 2.0}}}
    fields = CORE.normalize_particle_fields(case["models"]["physics"]["particles"]["fields"], case)
    assert fields == [{"name": "Psi", "expressions": ["(3.0)*y*(2.0)"]}]


@pytest.mark.parametrize("value, replacement", [
    (0.5, "expression"),
    ("0.5 + x", "expression"),
    ({"params": {"k": 2}, "value": "k*x"}, "retired"),
    ({"background": 0, "regions": []}, "retired"),
    ({"expression": "xn"}, "xn"),
    ({"expression": "1", "components": "normal"}, "takes expression"),
])
def test_retired_particle_forms_name_their_replacement(value, replacement, capsys):
    """!
    @brief The bare values, `params`, and domain-normalized variables are refused.
    @param[in] value Field value.
    @param[in] replacement Expected fragment of the refusal.
    @param[in] capsys Pytest output capture.
    """
    with pytest.raises(SystemExit):
        validate(*load_configs({"particles_per_cell": 1, "fields": {"Psi": value}}))
    assert replacement in capsys.readouterr().err


def test_velocity_sources_take_regions_and_refuse_draws():
    """!
    @brief Regions are part of every source's language; draws are keyed on a particle, so only it may use them.
    """
    case, _, _ = load_configs({"particles_per_cell": 1})
    regions = {"background": 1.0, "regions": [{"shape": "slab", "axis": "y", "from": 0.2, "to": 0.4, "value": 2.0}]}
    source = CORE.resolve_velocity_source({"components": "normal", "expression": regions}, case, "src", True)
    assert "where(" in source["expressions"][0]
    source = CORE.resolve_velocity_source({"expression": {"w": regions}}, case, "src", False)
    assert "where(" in source["expressions"][2]
    with pytest.raises(ValueError, match="pid"):
        CORE.resolve_velocity_source({"components": "normal", "expression": "pid"}, case, "src", True)


@pytest.mark.parametrize("particles, message", [
    ({"particles_per_cell": 1, "fields": {"Velocity": {"expression": 1.0}}}, "not a field a case can set"),
    ({"particles_per_cell": 1, "fields": {"Psi": {"expression": "foo + 1"}}}, "fields.Psi"),
    ({"particles_per_cell": 1, "fields": {"Psi": "missing"}}, "no source named"),
    ({"particles_per_cell": 1, "fields": {}}, "non-empty mapping"),
    ({"fields": {"Psi": {"expression": 1.0}}}, "needs particles"),
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
    case, solver, monitor = load_configs({"particles_per_cell": 1, "fields": {"Psi": {"expression": 1.0}}})
    solver.setdefault("verification", {})["sources"] = {"scalar": {"expression": "1.0"}}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)
    assert "verification.sources.scalar" in capsys.readouterr().err


def test_inlet_fields_conflict_with_a_prescribed_scalar_in_a_list_of_blocks(capsys):
    """!
    @brief The inlet-fields conflict is found when boundary_conditions is one list per block.
    @param[in] capsys Pytest output capture.
    """
    case, solver, monitor = load_configs({})
    inlet = next(entry for entry in case["boundary_conditions"] if entry["type"] == "INLET")
    inlet["particles"] = {"particles_per_cell": 1, "fields": {"Psi": {"expression": 1.0}}}
    case["boundary_conditions"] = [case["boundary_conditions"]]
    solver.setdefault("verification", {})["sources"] = {"scalar": {"expression": "1.0"}}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)
    assert "particles.fields cannot be combined" in capsys.readouterr().err


def test_a_point_source_without_draws_warns(capsys):
    """!
    @brief Particles released from one point all get the same value unless the value draws.
    @param[in] capsys Pytest output capture.
    """
    point = {"x": 0.5, "y": 0.5, "z": 0.5}
    validate(*load_configs({"count": 100, "init_mode": "PointSource", "point_source": point,
                            "fields": {"Psi": {"expression": "x"}}}))
    assert "has no random draw" in capsys.readouterr().err
    validate(*load_configs({"count": 100, "init_mode": "PointSource", "point_source": point,
                            "fields": {"Psi": {"expression": "where(uniform() < 0.5, 0, 1)"}}}))
    assert "has no random draw" not in capsys.readouterr().err


@pytest.mark.parametrize(
    "particles, expected",
    [
        ({"particles_per_cell": 2, "init_mode": "Volume"}, 2 * 8 * 8 * 16),
        ({"particles_per_cell": 2, "init_mode": "Surface"}, 2 * 8 * 8),
        ({"particles_per_cell": 0.5, "init_mode": "Volume"}, 512),
        ({"count": 100, "init_mode": "PointSource", "point_source": {"x": 0.5, "y": 0.5, "z": 0.5}}, 100),
        ({}, 0),
    ],
)
def test_seed_count_follows_the_cells(particles, expected):
    """!
    @brief Volume seeds particles_per_cell per cell of the grid, Surface per face cell of the
           first inlet, and a point source its count.
    @param[in] particles Mapping for `models.physics.particles`.
    @param[in] expected Seed count on the 8 x 8 x 16 fixture with its inlet on -Zeta.
    """
    case, _, _ = load_configs(particles)
    assert CORE.particle_seed_count(case, str(FIXTURES / "case.yml")) == expected
    lines = []
    CORE.parse_and_add_model_flags(case, lines, CORE.particle_seed_count(case, str(FIXTURES / "case.yml")))
    assert (f"-numParticles {expected}" in lines) == (expected > 0)


@pytest.mark.parametrize(
    "particles, message",
    [
        ({"count": 100, "init_mode": "Volume"}, "count is retired for Volume"),
        ({"count": 0}, "count is retired for Surface"),
        ({"particles_per_cell": 1, "init_mode": "PointSource", "point_source": {"x": 0.5, "y": 0.5, "z": 0.5}},
         "a PointSource has no cells"),
        ({"particles_per_cell": 0, "init_mode": "Volume"}, "must be a positive number"),
        ({"particles_per_cell": 0.001, "init_mode": "Surface"}, "seeds no particle"),
    ],
)
def test_seed_inputs_are_checked(particles, message, capsys):
    """!
    @brief Cell seeding takes particles_per_cell and a point source a count; neither takes the other.
    @param[in] particles Mapping for `models.physics.particles`.
    @param[in] message Pattern the refusal must contain.
    @param[in] capsys Pytest stream capture.
    """
    with pytest.raises(SystemExit):
        validate(*load_configs(particles))
    assert message in capsys.readouterr().err


def inlet_case(particles):
    """!
    @brief The valid fixture case with `particles` on its inlet and no seeded population.
    @param[in] particles Mapping for the inlet's `particles`.
    @return `(case, solver, monitor)` mappings.
    """
    case, solver, monitor = load_configs({"restart_mode": "init"})
    entry = next(bc for bc in case["boundary_conditions"] if bc["type"] == "INLET")
    entry["particles"] = deepcopy(particles)
    return case, solver, monitor


def test_an_inlet_becomes_runtime_options_and_reflects_by_default():
    """!
    @brief An injecting inlet's options reach the runtime; a diffusing particle stays in the domain.
    """
    case, _, _ = inlet_case({"particles_per_cell": "n0 * (1 + t)", "fields": {"Psi": {"expression": "y"}}})
    case["constants"] = {"n0": 1000.0}
    prepared = CORE.validate_and_prepare_boundary_conditions(case)
    entry = next(item for block in prepared for item in block if "particles" in item)
    assert entry["params"]["particle_policy"] == "reflect"
    assert entry["particles"]["particles_per_cell"] == "(1000.0) * (1 + t)"
    assert entry["particles"]["fields"] == [{"name": "Psi", "expressions": ["y"]}]
    assert CORE.case_injects_particles(case)

    entry = next(bc for bc in case["boundary_conditions"] if bc.get("particles"))
    entry["particle_policy"] = "remove"
    prepared = CORE.validate_and_prepare_boundary_conditions(case)
    entry = next(item for block in prepared for item in block if "particles" in item)
    assert entry["params"]["particle_policy"] == "remove"


def test_an_inlet_lets_a_case_set_fields_with_no_seeded_particles():
    """!
    @brief The case's `fields` apply to an inlet's particles, so they need no seeded population.
    """
    case, solver, monitor = inlet_case({"particles_per_cell": 1000.0})
    case["models"]["physics"]["particles"]["fields"] = {"Psi": {"expression": "x"}}
    validate(case, solver, monitor)


@pytest.mark.parametrize(
    "particles, message",
    [
        ({"fields": {"Psi": {"expression": "x"}}}, "takes particles_per_cell"),
        ({"particles_per_cell": 1.0, "rate": 2}, "takes particles_per_cell"),
        ({"particles_per_cell": 0.0}, "must be positive"),
        ({"particles_per_cell": "x + 1"}, "particles_per_cell"),
        ({"particles_per_cell": True}, "particles_per_cell"),
        ({"particles_per_cell": 1.0, "fields": {"velocity": {"expression": "1"}}}, "not a field a case can set"),
        ({"number_density": 1.0}, "number_density is retired"),
    ],
)
def test_invalid_inlets_are_refused(particles, message):
    """!
    @brief An inlet's particles per cell is a positive value in `t`, and its fields are settable ones.
    @param[in] particles The inlet's `particles` mapping.
    @param[in] message Pattern the refusal must match.
    """
    case, _, _ = inlet_case(particles)
    with pytest.raises(ValueError, match=message):
        CORE.validate_and_prepare_boundary_conditions(case)


def test_only_an_inlet_injects():
    """!
    @brief A wall or an outlet refuses `particles`.
    """
    case, _, _ = load_configs({"restart_mode": "init"})
    entry = next(bc for bc in case["boundary_conditions"] if bc["type"] == "WALL")
    entry["particles"] = {"particles_per_cell": 1.0}
    with pytest.raises(ValueError, match="only an INLET face injects"):
        CORE.validate_and_prepare_boundary_conditions(case)


def test_a_prescribed_scalar_source_conflicts_with_an_inlet_value(capsys):
    """!
    @brief A verification scalar prescribes Psi at every step, so an inlet's value would not survive.
    @param[in] capsys Pytest stream capture.
    """
    case, solver, monitor = inlet_case({"particles_per_cell": 1000.0, "fields": {"Psi": {"expression": "y"}}})
    solver.setdefault("verification", {})["sources"] = {"scalar": {"expression": "x"}}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)
    assert "particles.fields cannot be combined" in capsys.readouterr().err
