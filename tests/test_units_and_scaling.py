"""!
@file test_units_and_scaling.py
@brief Pytest coverage for the physical-units ingress rule and its dimension table.

Every configuration input is physical and is converted to solver units exactly once.
The completeness checks here run the enumeration of `tests/tooling/audit_units.py`, the
checker of the `units.nondimensionalization` contract, so a new key cannot reach the
solver without a recorded dimension; the rest check conversions at non-unit scales.
"""

import importlib.machinery
import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV_CORE = REPO_ROOT / "picurv_cli" / "core.py"
IC_GENERATOR = REPO_ROOT / "generators" / "ic.gen"
AUDIT = REPO_ROOT / "tests" / "tooling" / "audit_units.py"
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "valid"

#: Non-unit reference scales for the ingress tests: L = 2 m, U = 4 m/s, so T = 0.5 s.
SCALES = {"length_ref": 2.0, "velocity_ref": 4.0}


def load_module(name: str, path: Path):
    """!
    @brief Load a repository script as an importable module.
    @param[in] name Module name to register.
    @param[in] path Source path of the script.
    @return Loaded module object.
    """
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def core():
    """!
    @brief Load the conductor core once for this module.
    @return Loaded `picurv_cli/core.py` module.
    """
    return load_module("picurv_units_core", PICURV_CORE)


@pytest.fixture(scope="module")
def audit():
    """!
    @brief Load the units audit once for this module.
    @return Loaded `tests/tooling/audit_units.py` module.
    """
    return load_module("picurv_units_audit", AUDIT)


def test_every_solver_input_records_its_dimension(core, audit):
    """!
    @brief Every schema key, BC handler parameter, and provider parameter has an entry.
    @param[in] core Loaded conductor core.
    @param[in] audit Loaded units audit.
    """
    uncovered = audit.uncovered_inputs(core)
    assert not uncovered, "No recorded physical dimension: " + ", ".join(uncovered)


def test_the_units_audit_passes_against_page_19(audit):
    """!
    @brief Both indexes on page 19 match the code, and the audit says so.
    @param[in] audit Loaded units audit.
    """
    assert audit.main() == 0


def test_every_table_entry_names_a_known_conversion_site(core):
    """!
    @brief Each entry pairs a dimension (or none) with a declared conversion site.
    @param[in] core Loaded conductor core.
    """
    tables = [core.INPUT_QUANTITIES, core.BC_PARAM_QUANTITIES]
    tables.extend(core.IC_PARAM_QUANTITIES.values())
    for table in tables:
        for key, (dimension, site) in table.items():
            assert site in core.INPUT_CONVERSION_SITES, key
            assert dimension is None or (isinstance(dimension, tuple) and len(dimension) == 3), key
            if site == "cli":
                assert dimension is not None, f"{key} is converted by the CLI but has no dimension"


def scaled_case(core, **overrides):
    """!
    @brief Load the valid fixture case with the non-unit reference scales applied.
    @param[in] core Loaded conductor core.
    @param[in] overrides Top-level case keys to replace.
    @return Parsed case mapping.
    """
    case_cfg = core.read_yaml_file(str(FIXTURES / "case.yml"))
    case_cfg["properties"]["scaling"] = dict(SCALES)
    case_cfg.update(overrides)
    return case_cfg


def test_timestep_and_grid_lengths_use_the_reference_time_and_length(core):
    """!
    @brief `dt_physical` is divided by L/U and staged grid lengths by L.
    @param[in] core Loaded conductor core.
    """
    assert core.to_solver_units(0.1, ("case", "run_control", "dt_physical"), SCALES) == pytest.approx(0.2)
    assert core.input_reference_scale(("case", "grid", "source_file"), SCALES) == pytest.approx(2.0)
    assert core.input_reference_scale(("case", "grid", "generator"), SCALES) == pytest.approx(2.0)


def test_boundary_velocities_and_fluxes_are_converted(core):
    """!
    @brief BC velocities are divided by U and volume fluxes by U L^2.
    @param[in] core Loaded conductor core.
    """
    prepared = core.validate_and_prepare_boundary_conditions(scaled_case(core))
    inlet = next(bc for bc in prepared[0] if bc["handler"] == "constant_velocity")
    assert inlet["params"]["vz"] == pytest.approx(1.0 / 4.0)
    assert core.quantity_to_solver_units(16.0, core.BC_PARAM_QUANTITIES["target_flux"], SCALES,
                                         "target_flux") == pytest.approx(1.0)


def test_point_source_is_a_physical_position(core):
    """!
    @brief The particle point source is divided by L, like the grid it sits in.
    @param[in] core Loaded conductor core.
    """
    case_cfg = scaled_case(core)
    case_cfg["models"]["physics"]["particles"] = {
        "count": 10, "init_mode": "PointSource", "point_source": {"x": 0.5, "y": 1.0, "z": -2.0},
    }
    lines = []
    core.parse_and_add_model_flags(case_cfg, lines)
    assert "-psrc_x 0.25" in lines
    assert "-psrc_y 0.5" in lines
    assert "-psrc_z -1.0" in lines


def test_statistics_window_times_are_physical(core):
    """!
    @brief Window start, end, and cadence are divided by the reference time L/U.
    @param[in] core Loaded conductor core.
    """
    monitor_cfg = {"field_statistics": {"enabled": True, "windows": [{
        "name": "w", "start_time": 1.0, "end_time": 3.0, "weighting": "sample",
        "time_cadence": 0.25, "fields": [{"field": "Ucat", "moments": ["first"]}],
    }]}}
    case_cfg = {"models": {"physics": {"particles": {"count": 0}}},
                "properties": {"scaling": dict(SCALES)}}
    lines = core.resolve_field_statistics_flags(monitor_cfg, case_cfg)
    assert "-field_statistics_window_0_start_time 2.0" in lines
    assert "-field_statistics_window_0_end_time 6.0" in lines
    assert "-field_statistics_window_0_time_cadence 0.5" in lines


def test_spectral_provider_parameters_are_converted_before_the_provider_runs(core):
    """!
    @brief Spectral velocities are divided by U and wavenumbers multiplied by L.
    @param[in] core Loaded conductor core.
    """
    periodic = [{"face": face, "type": "PERIODIC", "handler": "geometric"}
                for face in ("-Xi", "+Xi", "-Eta", "+Eta", "-Zeta", "+Zeta")]
    resolved = core.resolve_initial_condition_config(
        {"mode": "generated", "generator": "spectral_random_velocity", "params": {
            "random": {"distribution": "gaussian", "mean": [1.0, 0.0, -2.0]},
            "spectrum": {"type": "k4_exponential", "k0": 4.0, "k_cut": 10.0},
            "normalization": {"type": "component_rms", "target": 2.0},
        }},
        [periodic], SCALES)
    params = resolved["params"]
    assert params["random"]["mean"] == pytest.approx([0.25, 0.0, -0.5])
    assert params["spectrum"]["k0"] == pytest.approx(8.0)
    assert params["spectrum"]["k_cut"] == pytest.approx(20.0)
    assert params["normalization"]["target"] == pytest.approx(0.5)


def test_expressions_see_physical_coordinates_and_write_solver_units():
    """!
    @brief An expression field is evaluated at x_phys = L x and divided by its field scale.
    """
    import numpy as np

    ic_gen = load_module("picurv_units_ic_gen_expr", IC_GENERATOR)
    axis = np.linspace(0.0, 1.0, 4)
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    nodes = np.stack([x, y, z], axis=-1)
    config = {"expressions": {"u": "x", "v": "2.0", "w": "0.0",
                              "u_xi": "1.0", "u_eta": "0.0", "u_zeta": "0.0"},
              "scale": 1.0, "zero_tolerance": 0.0, "max_magnitude": None}
    ucat = ic_gen.generate_expression_field(nodes, config, "Ucat", length_ref=2.0, velocity_ref=4.0)
    reference = ic_gen.generate_ucat(nodes, config)
    # u = x: the physical coordinate is 2x, and the value is written divided by U = 4.
    assert np.allclose(ucat[..., 0], 2.0 * reference[..., 0] / 4.0)
    assert np.allclose(ucat[..., 1], 2.0 / 4.0)
    ucont = ic_gen.generate_expression_field(nodes, config, "Ucont", length_ref=2.0, velocity_ref=4.0)
    interior = ucont[1:-1, 1:-1, 1:-1, 0]
    # A face flux is divided by U L^2 = 16.
    assert np.allclose(interior, 1.0 / 16.0)


def test_file_initial_condition_is_staged_in_solver_units(core, tmp_path):
    """!
    @brief A physical payload is divided by U; another run's field is rescaled between runs.
    @param[in] core Loaded conductor core.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    import struct

    source = tmp_path / "velocity.dat"
    source.write_bytes(struct.pack(">ii", 1211214, 3) + struct.pack(">3d", 1.0, 2.0, 3.0))
    case_path = tmp_path / "case.yml"
    case_path.write_text("{}\n", encoding="utf-8")

    def staged_values(ic, run):
        """!
        @brief Resolve and stage one file IC, returning the staged scalars.
        @param[in] ic Initial-condition mapping.
        @param[in] run Run directory name.
        @return Staged scalar values.
        """
        resolved = core.resolve_initial_condition_config(ic, [], SCALES)
        summary = core.stage_initial_condition_file(str(tmp_path / run), str(case_path), resolved)
        data = Path(summary["staged"]).read_bytes()
        return struct.unpack(">3d", data[8:])

    physical = {"mode": "file", "field": "Ucat", "source_file": str(source)}
    assert staged_values(physical, "physical") == pytest.approx((0.25, 0.5, 0.75))

    earlier_case = tmp_path / "earlier.yml"
    earlier_case.write_text("properties:\n  scaling: {length_ref: 1.0, velocity_ref: 2.0}\n", encoding="utf-8")
    from_run = dict(physical, source_case=str(earlier_case))
    # Solver units of a run at U = 2 are rescaled to this run's U = 4.
    assert staged_values(from_run, "from_run") == pytest.approx((0.5, 1.0, 1.5))

    with pytest.raises(ValueError, match="length_scale"):
        core.resolve_initial_condition_config(
            {"mode": "file", "field": "Ucont", "source_file": str(source), "velocity_scale": 2.0}, [], SCALES)


def test_transition_notice_names_changed_inputs_only_at_non_unit_scales(core):
    """!
    @brief Inputs whose reading changed are named, but only when the scales are not one.
    @param[in] core Loaded conductor core.
    """
    case_cfg = scaled_case(core)
    case_cfg["models"]["physics"]["particles"] = {"count": 1, "init_mode": "PointSource",
                                                  "point_source": {"x": 0.0, "y": 0.0, "z": 0.0}}
    notices = core.physical_units_transition_notices(case_cfg, {}, {})
    assert len(notices) == 1 and "point_source" in notices[0]
    case_cfg["properties"]["scaling"] = {"length_ref": 1.0, "velocity_ref": 1.0}
    assert core.physical_units_transition_notices(case_cfg, {}, {}) == []


@pytest.mark.parametrize(
    "dimension, expected",
    [
        ("LENGTH", 2.0),
        ("VELOCITY", 3.0),
        ("TIME", 2.0 / 3.0),
        ("WAVENUMBER", 0.5),
        ("VOLUME_FLUX", 12.0),
        ("DIFFUSIVITY", 6.0),
        ("PRESSURE", 5.0 * 9.0),
        ("DIMENSIONLESS", 1.0),
    ],
)
def test_reference_scale_of_each_dimension(core, dimension, expected):
    """!
    @brief Each dimension's reference scale is the product of powers of L, U, and rho.
    @param[in] core Loaded conductor core.
    @param[in] dimension Name of the dimension constant.
    @param[in] expected Scale at L=2, U=3, rho=5.
    """
    scales = {"length_ref": 2.0, "velocity_ref": 3.0, "density": 5.0}
    assert core.reference_scale(getattr(core, dimension), scales) == pytest.approx(expected)
