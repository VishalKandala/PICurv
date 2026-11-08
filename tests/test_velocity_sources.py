"""!
@file test_velocity_sources.py
@brief Pytest coverage for prescribed velocity sources: `sources:`, the `prescribed_velocity`
       inlet, and the `source` initial condition.

The CLI resolves a source (inline or named), inlines its constants, applies the target's
constraints, and passes it to the runtime as quoted `-face_source_<id>_*` or `-ic_source_*`
options; `src/expression.c` compiles and evaluates it. A data source on an inlet is staged as
a PICSLICE instead, which the bcs file names. The language itself is covered by
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
        case_cfg=case)
    assert ic["source"]["expressions"] == inlet["source"]["expressions"]


def test_sources_become_quoted_runtime_options():
    """!
    @brief Each value is one quoted option; the frame and the bulk target ride alongside.
    """
    source = CORE.apply_source_constraints(
        CORE.resolve_velocity_source({"components": "normal", "expression": "6*eta*(1-eta)"}, {},
                                     "src", allow_normal=True),
        {"scale": 2, "bulk_velocity": 1.0}, "constraints", allowed={"scale", "bulk_velocity"})
    lines = CORE.velocity_source_control_lines("-face_source_3", source, "src")
    assert lines == [
        "-face_source_3_frame normal",
        '-face_source_3_expr_0 "(6*eta*(1-eta))*(2.0)"',
        '-face_source_3_bulk_velocity "1.0"',
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
    (None, {"handler": "constant_velocity", "params": {"vx": 0, "vy": 0, "vz": 1}}, None,
     "is retired: use handler 'prescribed_velocity' with source: {components: normal"),
    (None, {"type": "WALL", "handler": "noslip", "source": DUCT}, None,
     "only read by the 'prescribed_velocity' and 'moving' handlers"),
    (None, None, {"source": {"type": "field_slice", "field_file": "u.dat"}}, "only an inlet reads it"),
    (None, {"handler": "prescribed_velocity", "source": {"type": "file", "path": "inlet.picslice"},
            "constraints": {"scale": -1.0}}, None, "must be positive"),
    (None, {"handler": "prescribed_velocity", "source": {"type": "generated", "generator": "square_duct_poiseuille"}},
     None, "is retired"),
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
    ({"duct": DUCT}, None, {"source": "duct", "constraints": {"target_flux": 1.0}},
     "not constraints of this target"),
    ({"duct": DUCT}, None, {"source": "duct", "mode": "generated"}, "are retired"),
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
                                              case_cfg=case)
    text = capsys.readouterr().err + str(caught.value)
    assert message in text


def test_a_data_source_is_staged_and_scaled(tmp_path):
    """!
    @brief A named file source is staged in solver units, scaled, and named in the bcs file;
           only its bulk target reaches the control file.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    profile = tmp_path / "inlet.picslice"
    case, _, _ = load_configs(
        sources={"measured": {"type": "file", "path": str(profile)}},
        inlet={"handler": "prescribed_velocity", "source": "measured",
               "constraints": {"scale": 2.0, "bulk_velocity": "1 + 0.5*t"}})
    jm, im = (case["grid"]["programmatic_settings"][key] for key in ("jm", "im"))
    speeds = [0.25 * (n % 4) for n in range(jm * im)]
    profile.write_text("PICSLICE\n1\n%d %d\n%s\n" % (jm, im, "\n".join(map(str, speeds))), encoding="utf-8")
    inlet = CORE.validate_and_prepare_boundary_conditions(case)[0][4]
    assert CORE.profile_data_source(inlet)["type"] == "file"
    assert CORE.velocity_source_control_lines("-face_source_0", inlet["source"], "src") == [
        '-face_source_0_bulk_velocity "1 + 0.5*t"']

    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)
    generated = CORE.generate_multi_block_bcs(
        str(run_dir), "unit", case, {"Case": str(tmp_path / "case.yml"), "Solver": "solver.yml",
                                     "Monitor": "monitor.yml"})
    line = next(text for text in Path(generated[0]).read_text(encoding="utf-8").splitlines()
                if text.startswith("-Zeta"))
    assert "prescribed_velocity" in line and "source_id=0" in line
    staged = Path(line.split("source_file=")[1].split()[0])
    values = [float(token) for token in staged.read_text(encoding="utf-8").split()[4:]]
    velocity_ref = CORE.resolve_reference_scales(case)["velocity_ref"]
    assert values == pytest.approx([2.0 * speed / velocity_ref for speed in speeds])


DUCT_SEED = {"generator": "spectral_velocity",
             "params": {"seed": 7, "wall_axes": ["i", "j"], "normalization": {"type": "component_rms", "target": 0.05},
                        "initial_spectra": [{"task": "line_spectrum", "axes": ["k"],
                                             "fixed_indices": {"i": 3, "j": 3}}]}}


def resolve_initial(case):
    """!
    @brief Resolve a case's initial condition against its prepared boundaries.
    @param[in] case Case mapping.
    @return Resolved initial-condition contract.
    """
    return CORE.resolve_initial_condition_config(
        case["properties"]["initial_conditions"], CORE.validate_and_prepare_boundary_conditions(case),
        CORE.resolve_reference_scales(case), case_cfg=case)


def test_a_sum_adds_a_generator_to_a_source_and_projects():
    """!
    @brief A mean source plus a fluctuation-only seed between an inlet and an outlet: the seed
           is staged as its provider stages it, added at runtime, and the sum is projected.
    """
    case, _, _ = load_configs(
        sources={"duct": DUCT},
        initial_conditions={"source": {"sum": ["duct", DUCT_SEED]}, "constraints": {"scale": 2.0}})
    ic = resolve_initial(case)
    assert ic["kind"] == "spectral_velocity"
    assert CORE.is_generated_ic_provider(ic) and "bulk_velocity" not in ic["params"]
    assert ic["source"]["expressions"][2] == "((36*(1.5)*xi*(1-xi)*eta*(1-eta)))*(2.0)"
    lines = CORE.velocity_source_control_lines("-ic_source", ic["source"], "src")
    assert lines[-2:] == ["-ic_source_field_weight 2.0", "-ic_source_project true"]


def test_source_terms_add_and_project_is_optional():
    """!
    @brief Two source terms in one frame are added; a plain source is not projected unless asked.
    """
    swirl = {"components": "curvilinear", "expression": {"u_xi": "0.1*zeta"}}
    case, _, _ = load_configs(
        sources={"duct": DUCT},
        initial_conditions={"source": {"sum": ["duct", swirl]}, "constraints": {"project": False}})
    ic = resolve_initial(case)
    assert ic["kind"] == "source" and not ic["source"]["project"]
    case, _, _ = load_configs(sources={"duct": DUCT}, initial_conditions={"source": {"sum": ["duct", swirl]}})
    assert resolve_initial(case)["source"]["project"] is True
    assert ic["source"]["expressions"][0] == "(0.0) + (0.1*zeta)"
    plain, _, _ = load_configs(initial_conditions={"source": "duct", "constraints": {"project": True}},
                               sources={"duct": DUCT})
    assert resolve_initial(plain)["source"]["project"] is True
    assert resolve_initial(load_configs(sources={"duct": DUCT},
                                        initial_conditions={"source": "duct"})[0])["source"]["project"] is False


@pytest.mark.parametrize("terms, message", [
    ([DUCT_SEED, DUCT_SEED], "at most one generator"),
    ([DUCT, {"components": "cartesian", "expression": {"w": "1"}}], "mix the"),
    ([{"generator": "spectral_velocity", "seed": 1}], "unknown keys"),
    ([], "non-empty list"),
])
def test_invalid_sums_are_refused(terms, message):
    """!
    @brief Two generators, mixed frames, misplaced generator keys and an empty sum are refused.
    @param[in] terms Sum terms.
    @param[in] message Expected fragment of the error.
    """
    case, _, _ = load_configs(initial_conditions={"source": {"sum": terms}})
    with pytest.raises(ValueError, match=message):
        resolve_initial(case)


def test_a_lone_generator_is_a_source_and_is_not_projected():
    """!
    @brief A generator alone is staged as it writes it: nothing is added and no projection runs
           unless asked for.
    """
    case, _, _ = load_configs(initial_conditions={"source": DUCT_SEED})
    ic = resolve_initial(case)
    assert ic["kind"] == "spectral_velocity"
    assert ic["source"]["expressions"] == [] and ic["source"]["project"] is False
    case, _, _ = load_configs(initial_conditions={"source": DUCT_SEED, "constraints": {"project": True}})
    assert resolve_initial(case)["source"]["project"] is True


@pytest.mark.parametrize("wrap", [lambda term: term, lambda term: {"sum": [DUCT, term]}])
def test_a_generator_inside_a_source_resolves_its_paths(wrap):
    """!
    @brief The case's relative `config_file` and `script` are anchored wherever the generator
           is written.
    @param[in] wrap Places the generator term in a source.
    """
    term = {"generator": "ic_gen", "params": {"field": "Ucat", "config_file": "ic.cfg", "script": "tools/ic.py"}}
    case = {"properties": {"initial_conditions": {"source": wrap(term)}}}
    CORE.absolutize_case_external_paths(case, "/work/case.yml")
    assert term["params"]["config_file"] == "/work/ic.cfg"
    assert term["params"]["script"] == "/work/tools/ic.py"


@pytest.mark.parametrize("source, name", [
    (DUCT_SEED, "spectral_velocity"),
    ({"sum": [DUCT, DUCT_SEED]}, "spectral_velocity"),
    (DUCT, None),
    ("duct", None),
])
def test_initial_field_term_finds_the_generator(source, name):
    """!
    @brief The generator is found alone or in a sum, and a plain source has none.
    @param[in] source Initial-condition source.
    @param[in] name Expected generator, or None.
    """
    assert CORE.initial_field_term({"source": source}).get("generator") == name


@pytest.mark.parametrize("source, provider, execution", [
    (DUCT_SEED, "spectral_velocity", "precomputable"),
    ({"sum": [DUCT, DUCT_SEED]}, "spectral_velocity", "precomputable"),
    (DUCT, "source", "runtime-c"),
])
def test_the_asset_graph_plans_a_generator_inside_a_source(source, provider, execution):
    """!
    @brief `picurv precompute` plans a generator written inside a source; a plain source runs in C.
    @param[in] source Initial-condition source.
    @param[in] provider Expected provider name.
    @param[in] execution Expected execution class.
    """
    case, _, _ = load_configs(initial_conditions={"source": source})
    graph = CORE.build_case_asset_graph(case, str(FIXTURES / "case.yml"))
    planned = next(item for item in graph["providers"] if item["kind"] == "initial-condition")
    assert (planned["provider"], planned["execution"]) == (provider, execution)


def test_every_target_takes_the_same_sum_and_bulk_velocity():
    """!
    @brief A sum of expressions and a bulk velocity mean the same on an inlet and in the interior.
    """
    terms = {"sum": ["duct", {"components": "curvilinear", "expression": {"u_zeta": "0.1*sin(2*pi*t)"}}]}
    constraints = {"bulk_velocity": 1.5}
    case, _, _ = load_configs(sources={"duct": DUCT},
                              inlet={"handler": "prescribed_velocity", "source": terms, "constraints": constraints},
                              initial_conditions={"source": terms, "constraints": constraints})
    inlet = CORE.validate_and_prepare_boundary_conditions(case)[0][4]["source"]
    initial = resolve_initial(case)["source"]
    assert inlet["expressions"] == initial["expressions"]
    assert inlet["bulk_velocity"] == initial["bulk_velocity"] == "1.5"
    assert CORE.velocity_source_control_lines("-ic_source", initial, "src")[-1] == "-ic_source_project true"
    assert '-ic_source_bulk_velocity "1.5"' in CORE.velocity_source_control_lines("-ic_source", initial, "src")


def test_an_inlet_sum_adds_expressions_to_data(tmp_path):
    """!
    @brief A pulsation added to a measured profile: the profile is staged, the expressions travel
           as options, and the bcs file names both.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    profile = tmp_path / "inlet.picslice"
    pulse = {"components": "normal", "expression": "0.1*sin(2*pi*t)"}
    case, _, _ = load_configs(inlet={"handler": "prescribed_velocity",
                                     "source": {"sum": [{"type": "file", "path": str(profile)}, pulse]}})
    jm, im = (case["grid"]["programmatic_settings"][key] for key in ("jm", "im"))
    profile.write_text("PICSLICE\n1\n%d %d\n%s\n" % (jm, im, "\n".join(["1.0"] * (jm * im))), encoding="utf-8")
    inlet = CORE.validate_and_prepare_boundary_conditions(case)[0][4]
    assert CORE.profile_data_source(inlet)["type"] == "file"
    assert inlet["source"]["expressions"] == ["(0.1*sin(2*pi*t))"]
    found = CORE.raw_inlet_data_sources(case)
    assert [term["type"] for _, _, term in found] == ["file"]


@pytest.mark.parametrize("terms, message", [
    ([{"type": "file", "path": "a.picslice"}, {"type": "file", "path": "b.picslice"}], "at most one data term"),
    ([{"type": "file", "path": "a.picslice"}, DUCT], "must be in its normal frame"),
    ([DUCT, {"components": "normal", "expression": "1"}], "mix the"),
])
def test_invalid_inlet_sums_are_refused(terms, message, capsys):
    """!
    @brief Two data terms, data with a non-normal expression, and mixed frames are refused.
    @param[in] terms Sum terms.
    @param[in] message Expected fragment of the error.
    @param[in] capsys Pytest output capture.
    """
    case, _, _ = load_configs(inlet={"handler": "prescribed_velocity", "source": {"sum": terms}})
    with pytest.raises((SystemExit, ValueError)) as caught:
        CORE.validate_and_prepare_boundary_conditions(case)
    assert message in capsys.readouterr().err + str(caught.value)


def test_a_wall_generator_has_no_mean_of_its_own():
    """!
    @brief The retired `bulk_velocity` parameter names its replacement.
    """
    seed = {"generator": "spectral_velocity",
            "params": dict(DUCT_SEED["params"], bulk_velocity=1.0)}
    case, _, _ = load_configs(initial_conditions={"source": {"sum": [DUCT, seed]}})
    with pytest.raises(ValueError, match="fluctuations only"):
        resolve_initial(case)


@pytest.mark.parametrize("retired, replacement", [
    ({"mode": "generated", "generator": "zero"}, "u: 0, v: 0, w: 0"),
    ({"mode": "generated", "generator": "constant", "params": {"w_physical": 2.0}}, "w: 2.0"),
    ({"mode": "generated", "generator": "streamwise_constant", "params": {"velocity_physical": 3.0}}, "u_zeta: 3.0"),
    ({"mode": "generated", "generator": "poiseuille", "params": {"peak_velocity_physical": 1.5}}, "16*xi*(1-xi)"),
    ({"mode": "file", "field": "Ucat", "source_file": "a.dat"}, "type: file"),
    ({"mode": "generated", "generator": "ic_gen", "params": {}}, "generator: ic_gen"),
    ({"mode": "Constant", "u_physical": 1.0}, "u: 1.0"),
])
def test_retired_initial_conditions_name_their_replacement(retired, replacement):
    """!
    @brief Every retired spelling is refused with the source that replaces it.
    @param[in] retired Retired mapping.
    @param[in] replacement Expected fragment of the refusal.
    """
    case, _, _ = load_configs(initial_conditions=retired)
    with pytest.raises(ValueError, match="are retired") as caught:
        resolve_initial(case)
    assert replacement in str(caught.value)


def test_a_file_term_stages_a_field_and_may_be_summed():
    """!
    @brief A velocity file is a term like a generator: alone, it adds nothing; in a sum, sources add.
    """
    field = {"type": "file", "field": "Ucat", "path": "velocity.dat"}
    case, _, _ = load_configs(initial_conditions={"source": field})
    ic = resolve_initial(case)
    assert ic["kind"] == "file" and ic["source_file"] == "velocity.dat" and ic["source"]["expressions"] == []
    case, _, _ = load_configs(sources={"duct": DUCT}, initial_conditions={"source": {"sum": ["duct", field]}})
    ic = resolve_initial(case)
    assert ic["kind"] == "file" and ic["source"]["expressions"][2].startswith("(36*")
    flux = dict(field, field="Ucont")
    case, _, _ = load_configs(sources={"duct": {"expression": {"w": 1.0}}},
                              initial_conditions={"source": {"sum": ["duct", flux]}})
    with pytest.raises(ValueError, match="curvilinear or contravariant"):
        resolve_initial(case)


def test_an_analytical_field_is_its_source(tmp_path):
    """!
    @brief The retired selectors name the source, and a source's pressure and leading minus reach
           the runtime as a value PETSc reads as one option.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    source = {"components": "cartesian", "constants": {"nu": 0.01},
              "expression": {"u": "sin(x)", "v": "-cos(y)*exp(-2*nu*t)"},
              "pressure": "-0.25*cos(2*x)"}
    resolved = CORE.resolve_velocity_source(source, {}, "src", allow_normal=False)
    lines = CORE.velocity_source_control_lines("-ic_source", resolved, "src")
    assert '-ic_source_expr_1 "(-cos(y)*exp(-2*(0.01)*t))"' in lines
    assert '-ic_source_pressure "(-0.25*cos(2*x))"' in lines

    case, solver, monitor = load_configs(initial_conditions={"source": source})
    solver["operation_mode"] = {"eulerian_field_source": "analytical"}
    validate(case, solver, monitor)
    for refused, message in (
            ({"source": DUCT_SEED}, "neither a generator or file term"),
            ({"source": DUCT, "constraints": {"project": True}}, "neither a generator or file term")):
        case, solver, monitor = load_configs(initial_conditions=refused)
        solver["operation_mode"] = {"eulerian_field_source": "analytical"}
        with pytest.raises(SystemExit):
            validate(case, solver, monitor)
    case, solver, monitor = load_configs(initial_conditions={"source": source})
    validate(case, solver, monitor)  # A solved case may start from a known pressure too.
    case, solver, monitor = load_configs()
    solver["operation_mode"] = {"eulerian_field_source": "analytical", "analytical_type": "TGV3D"}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)


def test_contravariant_components_are_face_fluxes():
    """!
    @brief The flux frame belongs to the cells, adds to a staged Ucont, and refuses what reads a velocity.
    """
    flux = {"components": "contravariant", "expression": {"u_zeta": "2*x"}}
    case, _, _ = load_configs(initial_conditions={"source": flux})
    assert resolve_initial(case)["source"]["frame"] == "contravariant"
    staged = {"type": "file", "field": "Ucont", "path": "flux.dat"}
    case, _, _ = load_configs(sources={"flux": flux}, initial_conditions={"source": {"sum": ["flux", staged]}})
    assert resolve_initial(case)["field_code"] == 1

    case, _, _ = load_configs(initial_conditions={"source": flux, "constraints": {"bulk_velocity": 1.0}})
    with pytest.raises(ValueError, match="contravariant components are fluxes"):
        resolve_initial(case)
    case, _, _ = load_configs(inlet={"handler": "prescribed_velocity", "source": flux})
    with pytest.raises((SystemExit, ValueError)):
        CORE.validate_and_prepare_boundary_conditions(case)
    case, solver, monitor = load_configs(initial_conditions={"source": flux})
    solver["operation_mode"] = {"eulerian_field_source": "analytical"}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)


def test_a_face_reads_r_and_a_cell_does_not():
    """! @brief `r` is the distance from a face's centre, so only a face can read it. """
    case, _, _ = load_configs()
    face = CORE.resolve_velocity_source({"components": "normal", "expression": "1 - r*r"}, case, "src", True)
    assert face["expressions"] == ["1 - r*r"]
    with pytest.raises(ValueError, match="r"):
        CORE.resolve_velocity_source({"expression": {"u": "1 - r"}}, case, "src", False)


def test_constants_come_from_the_fluid_the_case_and_the_source():
    """! @brief Each level overrides the one before it. """
    case, _, _ = load_configs()
    case["properties"]["fluid"] = {"density": 2.0, "viscosity": 0.5}
    case["constants"] = {"U": 3.0, "H": 7.0}
    source = CORE.resolve_velocity_source(
        {"components": "normal", "expression": "U*H*nu*rho", "constants": {"H": 5.0}}, case, "src", True)
    assert source["expressions"] == ["(3.0)*(5.0)*(0.25)*(2.0)"]
    case["constants"] = {"U": "three"}
    with pytest.raises(ValueError, match="constants.U must be a number"):
        CORE.resolve_velocity_source({"components": "normal", "expression": "U"}, case, "src", True)


def test_an_initial_condition_reads_the_case_constants_too():
    """! @brief The same `constants:`, `rho` and `nu` reach every target's source. """
    case, _, _ = load_configs(initial_conditions={"source": {"expression": {"w": "U*nu*rho"}}})
    case["properties"]["fluid"] = {"density": 2.0, "viscosity": 0.5}
    case["constants"] = {"U": 3.0}
    assert resolve_initial(case)["source"]["expressions"][2] == "(3.0)*(0.25)*(2.0)"


def test_a_velocity_or_pressure_source_may_draw():
    """! @brief Draws reach the runtime from an inlet and an IC under their own seed; a target that takes a value at any point or time refuses them. """
    noisy = {"components": "cartesian", "expression": {"w": "1 + 0.1*normal()", "u": "0.1*normal(2)"},
             "pressure": "0.01*uniform()"}
    case, solver, monitor = load_configs(inlet={"handler": "prescribed_velocity",
                                                "source": {"components": "normal", "expression": "1 + 0.1*normal()"}},
                                         initial_conditions={"source": noisy})
    case["source_random_seed"] = 4242
    validate(case, solver, monitor)
    lines = []
    CORE.parse_and_add_model_flags(case, lines)
    assert "-source_random_seed 4242" in lines
    resolved = resolve_initial(case)["source"]
    assert resolved["expressions"][2] == "1 + 0.1*normal()" and resolved["pressure"] == "0.01*uniform()"
    assert "-particle_random_seed" not in " ".join(lines)

    case["source_random_seed"] = -1
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)
    del case["source_random_seed"]

    solver["operation_mode"] = {"eulerian_field_source": "analytical"}
    with pytest.raises(SystemExit):
        validate(case, solver, monitor)

    with pytest.raises(ValueError, match="flow_rate"):
        source = CORE.resolve_velocity_source({"components": "normal", "expression": "1"}, case, "src", True)
        CORE.apply_source_constraints(source, {"flow_rate": "1 + uniform()"}, "c", {"flow_rate"})


def test_a_flow_rate_is_a_target_of_the_face_alone():
    """! @brief A flow rate reaches the runtime, excludes a bulk velocity, and is no initial condition's. """
    case, _, _ = load_configs()
    source = CORE.resolve_velocity_source({"components": "normal", "expression": "1"}, case, "src", True)
    constrained = CORE.apply_source_constraints(source, {"flow_rate": "2 + t"}, "c", {"scale", "bulk_velocity", "flow_rate"})
    assert CORE.velocity_source_control_lines("-face_source_0", constrained, "src")[-1] == '-face_source_0_flow_rate "2 + t"'
    with pytest.raises(ValueError, match="either bulk_velocity or flow_rate"):
        CORE.apply_source_constraints(source, {"flow_rate": 1, "bulk_velocity": 1}, "c", {"bulk_velocity", "flow_rate"})
    with pytest.raises(ValueError, match="not constraints of this target"):
        CORE.apply_source_constraints(source, {"flow_rate": 1}, "c", {"scale", "bulk_velocity"})


def test_data_may_hold_velocities_and_expressions_add_in_its_frame(tmp_path):
    """!
    @brief A cartesian file names its frame to the runtime, and a sum must be written in that frame.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    case, _, _ = load_configs()
    data = {"type": "file", "path": str(tmp_path / "v.picslice"), "components": "cartesian", "repeat": True}
    source = CORE.resolve_velocity_source(data, case, "src", True)
    assert source["frame"] == "cartesian" and source["data"]["repeat"] is True
    assert CORE.velocity_source_control_lines("-face_source_0", dict(source, bulk_velocity=None), "src") == [
        "-face_source_0_frame cartesian"]
    mixed = {"sum": [data, {"components": "normal", "expression": "1"}]}
    with pytest.raises(ValueError, match="in its cartesian frame"):
        CORE.resolve_velocity_source(mixed, case, "src", True)
    added = CORE.resolve_velocity_source({"sum": [data, {"expression": {"w": "1"}}]}, case, "src", True)
    assert added["frame"] == "cartesian" and added["expressions"] == ["(0.0)", "(0.0)", "(1)"]
    with pytest.raises(ValueError, match="components"):
        CORE.resolve_velocity_source({"type": "file", "path": "p", "components": "swirl"}, case, "src", True)


@pytest.mark.parametrize(
    "header, values, components, message",
    [
        ("3\n1 1 3 0.1", "1 2 3", "cartesian", "has 3 values, expected 9"),
        ("3\n1 1 3", "1 2 3 4 5 6 7 8 9", "cartesian", "positive time between them"),
        ("1\n1 1", "1", "cartesian", "holds 1 component"),
        ("1\n1 1 3", "-1 2 3", "normal", "holds 3 component"),
    ],
)
def test_malformed_recordings_are_refused(header, values, components, message, tmp_path):
    """!
    @brief The staging check names a recording's wrong size, missing interval, or wrong components.
    @param[in] header Text after the PICSLICE token.
    @param[in] values Value text.
    @param[in] components What the case says the file holds.
    @param[in] message Expected refusal.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    path = tmp_path / "bad.picslice"
    path.write_text(f"PICSLICE\n{header}\n{values}\n", encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        CORE.validate_and_nondimensionalize_picslice(str(path), str(tmp_path / "out"), 2.0, (1, 1), components)


def test_a_recording_is_staged_in_solver_units_with_its_interval(tmp_path):
    """!
    @brief Velocities scale, signs survive, and the interval stays in seconds.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    path, out = tmp_path / "in.picslice", tmp_path / "out.picslice"
    path.write_text("PICSLICE\n2\n1 1 3 0.25\n-2 4 6\n8 10 12\n", encoding="utf-8")
    summary = CORE.validate_and_nondimensionalize_picslice(str(path), str(out), 2.0, (1, 1), "cartesian")
    assert summary["frame_count"] == 2 and summary["value_count"] == 6
    assert out.read_text(encoding="utf-8").split()[:7] == ["PICSLICE", "2", "1", "1", "3", "0.25", "-1.00000000e+00"]


def test_a_field_slice_recording_needs_its_steps_and_interval():
    """! @brief One field file, or a recording's pattern with the time between its frames. """
    base = {"type": "field_slice", "grid_file": "g.picgrid", "velocity_scale": 1.0, "slice": {"face": "+Zeta"}}
    with pytest.raises(ValueError, match="one field_file, or field_files"):
        CORE._normalize_profile_source(base, "src")
    recording = {"pattern": "c/step_{step:012d}/Ucat.dat", "first": 10, "last": 30, "every": 10}
    normalized = CORE._normalize_profile_source(dict(base, field_files=recording, frame_interval=0.5, components="cartesian"), "src")
    assert normalized["field_files"] == recording and normalized["frame_interval"] == 0.5
    with pytest.raises(ValueError, match="frame_interval"):
        CORE._normalize_profile_source(dict(base, field_file="u.dat", frame_interval=0.5), "src")
    with pytest.raises(ValueError, match="step"):
        CORE._normalize_profile_source(dict(base, field_files=dict(recording, pattern="u.dat"), frame_interval=0.5), "src")


def test_retired_generated_profile_is_not_offered():
    """! @brief The profile source's own refusal names only the types that exist. """
    with pytest.raises(ValueError, match="type: file or field_slice"):
        CORE._normalize_profile_source("nope", "src")


@pytest.mark.parametrize(
    "entry, message",
    [
        ({}, "needs a 'source'"),
        ({"source": {"components": "normal", "expression": "1"}}, "takes expressions in the cartesian or curvilinear"),
        ({"source": {"type": "file", "path": "p"}}, "takes expressions in the cartesian or curvilinear"),
        ({"source": {"expression": {"u": "1"}}, "constraints": {"bulk_velocity": 1}}, "not constraints of this target"),
    ],
)
def test_a_moving_wall_takes_expressions_along_the_wall(entry, message):
    """!
    @brief What a moving wall refuses.
    @param[in] entry The wall's keys.
    @param[in] message Expected refusal.
    """
    case, _, _ = load_configs()
    wall = next(bc for bc in case["boundary_conditions"] if bc["type"] == "WALL")
    wall.update(handler="moving", **entry)
    with pytest.raises(ValueError, match=message):
        CORE.validate_and_prepare_boundary_conditions(case)


def test_a_moving_wall_reaches_the_runtime_and_excludes_wall_functions():
    """! @brief Its source gets a face source id, and the wall functions, which assume rest, refuse it. """
    case, _, _ = load_configs()
    wall = next(bc for bc in case["boundary_conditions"] if bc["type"] == "WALL")
    wall.update(handler="moving", source={"expression": {"u": "2.0"}})
    prepared = [bc for bc in CORE.validate_and_prepare_boundary_conditions(case)[0] if bc["handler"] == "moving"]
    assert prepared[0]["source"]["frame"] == "cartesian" and "source_id" in prepared[0]["params"]
    case.setdefault("models", {}).setdefault("physics", {})["turbulence"] = {"wall_function": {"enabled": True}}
    with pytest.raises(ValueError, match="assume a wall at rest"):
        CORE.validate_and_prepare_boundary_conditions(case)


def test_a_recordings_frames_are_part_of_its_asset_identity(tmp_path):
    """!
    @brief Changing any frame of a recording changes what a precomputed profile is keyed on.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    (tmp_path / ".picurv-workspace.yml").write_text("schema_version: 1\n", encoding="utf-8")
    (tmp_path / "config").mkdir()
    case_path = str(tmp_path / "config" / "case.yml")
    for step in (10, 20):
        (tmp_path / f"U_{step}.dat").write_bytes(b"frame")
    spec = [{"source": {"field_files": {"pattern": "U_{step}.dat", "first": 10, "last": 20, "every": 10}}}]

    before = CORE._provider_source_fingerprints(spec, case_path)
    assert len(before) == 2
    (tmp_path / "U_20.dat").write_bytes(b"changed")
    after = CORE._provider_source_fingerprints(spec, case_path)
    assert before != after
    assert [v["sha256"] for v in before.values()][0] == [v["sha256"] for v in after.values()][0]


def test_a_recording_shorter_than_a_fresh_run_is_a_problem(tmp_path):
    """!
    @brief A new run ending past the last frame is refused; one ending on it, a loop, or a restart is not.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    recording = tmp_path / "recording.picslice"
    recording.write_text("PICSLICE\n3\n2 2 1 0.5\n" + "1.0\n" * 12, encoding="utf-8")
    source = {"type": "file", "path": "recording.picslice"}
    case_path = str(tmp_path / "case.yml")
    case = {"run_control": {"start_step": 0, "total_steps": 10, "dt_physical": 0.1}}
    problem = lambda spec: CORE.recording_outlast_problem(spec, case, case_path, "-Zeta")
    assert problem(source) is None
    case["run_control"]["total_steps"] = 11
    assert "repeat: true" in problem(source)
    assert problem(dict(source, repeat=True)) is None
    assert "repeat: true" in problem({"type": "field_slice", "frame_interval": 0.5, "field_files": {"first": 10, "last": 30, "every": 10}})
    case["run_control"]["start_step"] = 5
    assert problem(source) is None
    case["run_control"]["start_step"] = 0
    (tmp_path / "steady.picslice").write_text("PICSLICE\n1\n2 2\n" + "1.0\n" * 4, encoding="utf-8")
    assert problem({"type": "file", "path": "steady.picslice"}) is None
    assert problem({"type": "file", "path": "absent.picslice"}) is None


FLUCTUATIONS = {"reynolds_stress": {"uu": "0.01*(1 - r)", "ww": 0.04, "uw": -0.005},
                "length_scales": [0.25, 0.125], "time_scale": 0.5}


def test_inlet_fluctuations_are_resolved_in_solver_units():
    """!
    @brief Stresses stay physical expressions, the runtime dividing them by U_ref^2; lengths and time
           are converted here; omitted stress components are zero.
    """
    case, solver, monitor = load_configs(
        inlet={"handler": "prescribed_velocity", "source": {"components": "normal", "expression": 1.0},
               "fluctuations": FLUCTUATIONS})
    case["properties"]["scaling"] = {"length_ref": 0.5, "velocity_ref": 2.0}
    validate(case, solver, monitor)
    entry = CORE.validate_and_prepare_boundary_conditions(case)[0][4]
    fluct = entry["fluctuations"]
    assert fluct["stress"] == [["0.01*(1 - r)", "0.0", "0.04"], ["0.0", "-0.005", "0.0"]]
    assert fluct["length_scales"] == pytest.approx([0.5, 0.25])
    assert fluct["time_scale"] == pytest.approx(0.5 * 2.0 / 0.5)


@pytest.mark.parametrize("inlet, message", [
    ({"handler": "prescribed_velocity", "source": {"components": "normal", "expression": 1.0},
      "fluctuations": dict(FLUCTUATIONS, time_scale=0)}, "positive time"),
    ({"handler": "prescribed_velocity", "source": {"components": "normal", "expression": 1.0},
      "fluctuations": dict(FLUCTUATIONS, length_scales=[1.0])}, "two positive lengths"),
    ({"handler": "prescribed_velocity", "source": {"components": "normal", "expression": 1.0},
      "fluctuations": dict(FLUCTUATIONS, reynolds_stress={"uu": "0.01*t"})}, "reynolds_stress.uu"),
    ({"handler": "prescribed_velocity", "source": {"components": "normal", "expression": 1.0},
      "fluctuations": dict(FLUCTUATIONS, reynolds_stress={"uz": 0.1})}, "maps uu, vv, ww"),
    ({"handler": "prescribed_velocity", "source": {"components": "normal", "expression": 1.0},
      "fluctuations": {"reynolds_stress": {"uu": 0.01}}}, "takes reynolds_stress, length_scales and time_scale"),
])
def test_invalid_inlet_fluctuations_are_refused(inlet, message):
    """!
    @brief Each fluctuation key is checked, and the stresses take no `t`: the runtime evaluates them once.
    @param[in] inlet The inlet entry.
    @param[in] message Pattern the refusal must match.
    """
    case, _, _ = load_configs(inlet=inlet)
    with pytest.raises(ValueError, match=message):
        CORE.validate_and_prepare_boundary_conditions(case)


def test_only_a_prescribed_velocity_inlet_fluctuates():
    """!
    @brief A wall, or an inlet of another handler, refuses fluctuations.
    """
    case, _, _ = load_configs(inlet={"handler": "prescribed_velocity",
                                     "source": {"components": "normal", "expression": 1.0}})
    wall = next(bc for bc in case["boundary_conditions"] if bc["type"] == "WALL")
    wall["fluctuations"] = FLUCTUATIONS
    with pytest.raises(ValueError, match="only a prescribed_velocity INLET"):
        CORE.validate_and_prepare_boundary_conditions(case)
