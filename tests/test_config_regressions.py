"""!
@file test_config_regressions.py
@brief Pytest coverage for ingress, post-validation, and statistics-config regressions.
"""

import csv
import importlib.machinery
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "picurv_cli" / "picurv"
PICURV_CORE = REPO_ROOT / "picurv_cli" / "core.py"
FIXTURES = REPO_ROOT / "tests" / "fixtures"


def run_picurv(args, cwd=REPO_ROOT):
    """!
    @brief Run the `picurv` CLI for configuration-regression tests.
    @param[in] args Command-line style argument list supplied to the function.
    @param[in] cwd Working directory override supplied to the function.
    @return Value returned by `run_picurv()`.
    """
    cmd = [sys.executable, str(PICURV)] + list(args)
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, timeout=60, check=False)


def load_picurv_module():
    """!
    @brief Load the conductor core as an importable module for regression tests.
    @return Value returned by `load_picurv_module()`.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_regression_module", str(PICURV_CORE))
    spec = importlib.util.spec_from_loader("picurv_regression_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_maintained_examples_use_canonical_initial_condition_contract():
    """!
    @brief Verify maintained examples and the valid fixture do not regress to legacy IC YAML.
    """
    case_paths = sorted((REPO_ROOT / "examples").rglob("*.yml"))
    case_paths.append(FIXTURES / "valid" / "case.yml")
    checked = []
    for path in case_paths:
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        initial_conditions = ((payload.get("properties") or {}).get("initial_conditions"))
        if initial_conditions is None:
            continue
        checked.append(path)
        assert initial_conditions.get("mode") in {"generated", "file"}, path
        if initial_conditions["mode"] == "generated":
            assert initial_conditions.get("generator") in {
                "zero", "constant", "streamwise_constant", "poiseuille", "ic_gen"
            }, path
    assert checked


def test_audit_ingress_manifest_matches_c_ingress():
    """!
    @brief Test that audit ingress manifest matches c ingress.
    """
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tests" / "tooling" / "audit_ingress.py")],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "[OK] Ingress manifest matches setup/io PETSc option scan." in result.stdout


def test_generate_post_recipe_supports_legacy_aliases_and_legacy_input_extensions(tmp_path):
    """!
    @brief Test that generate post recipe supports legacy aliases and legacy input extensions.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)

    post_cfg = {
        "run_control": {
            "startTime": 5,
            "endTime": 15,
            "timeStep": 2,
        },
        "source_data": {
            "directory": "<solver_output_dir>",
            "input_extensions": {
                "eulerian": ".fld",
                "particle": "prt",
            },
        },
        "statistics_pipeline": {
            "output_prefix": "stats/BrownianStats",
            "tasks": [{"task": "msd"}],
        },
        "io": {
            "output_directory": "viz",
            "output_filename_prefix": "Field",
            "particle_filename_prefix": "Particle",
            "output_particles": True,
            "eulerian_fields": ["Ucat"],
            "particle_fields": ["velocity"],
        },
    }

    recipe_path = Path(
        picurv.generate_post_recipe_file(
            str(run_dir),
            "demo_run",
            post_cfg,
            {"Case": "case.yml", "Post-Profile": "post.yml"},
        )
    )
    content = recipe_path.read_text(encoding="utf-8")

    assert "startTime = 5" in content
    assert "endTime = 15" in content
    assert "timeStep = 2" in content
    assert "eulerianExt = fld" in content
    assert "particleExt = prt" in content
    assert "statistics_pipeline = ComputeMSD" in content
    assert "statistics_output_prefix = stats/BrownianStats" in content


def test_generate_post_recipe_defaults_statistics_output_prefix_under_monitor_output_root(tmp_path):
    """!
    @brief Test that bare statistics prefixes default under the monitor output root statistics subdirectory.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)

    post_cfg = {
        "statistics_pipeline": {
            "tasks": [{"task": "msd"}],
        },
        "io": {
            "output_directory": "viz",
            "output_filename_prefix": "Field",
            "particle_filename_prefix": "Particle",
        },
    }
    monitor_cfg = {
        "io": {
            "directories": {
                "output": "results",
            },
        },
    }

    recipe_path = Path(
        picurv.generate_post_recipe_file(
            str(run_dir),
            "demo_run",
            post_cfg,
            {"Case": "case.yml", "Post-Profile": "post.yml"},
            monitor_cfg,
        )
    )
    content = recipe_path.read_text(encoding="utf-8")

    assert "statistics_pipeline = ComputeMSD" in content
    assert "statistics_output_prefix = results/statistics/Stats" in content


def test_validate_post_rejects_unsupported_eulerian_task(tmp_path):
    """!
    @brief Test that validate post rejects unsupported eulerian task.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    invalid_post = tmp_path / "post_invalid_task.yml"
    invalid_post.write_text(
        "\n".join(
            [
                "run_control:",
                "  start_step: 0",
                "  end_step: 10",
                "  step_interval: 1",
                "eulerian_pipeline:",
                "  - task: imaginary_kernel",
                "io:",
                "  output_directory: viz",
                "  output_filename_prefix: Field",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = run_picurv(["validate", "--post", str(invalid_post)])

    assert result.returncode == 1
    assert "unsupported eulerian task 'imaginary_kernel'" in result.stderr


def test_validate_post_rejects_normalize_field_for_non_pressure_field(tmp_path):
    """!
    @brief Test that validate post rejects normalize field for non pressure field.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    invalid_post = tmp_path / "post_invalid_normalize.yml"
    invalid_post.write_text(
        "\n".join(
            [
                "run_control:",
                "  start_step: 0",
                "  end_step: 10",
                "  step_interval: 1",
                "eulerian_pipeline:",
                "  - task: normalize_field",
                "    field: Ucat",
                "    reference_point: [1, 1, 1]",
                "io:",
                "  output_directory: viz",
                "  output_filename_prefix: Field",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = run_picurv(["validate", "--post", str(invalid_post)])

    assert result.returncode == 1
    assert "currently only supports 'P'" in result.stderr


def test_validate_rejects_misplaced_grid_da_processors_under_generator(tmp_path):
    """!
    @brief Test that validate rejects DMDA layout keys nested under grid.generator.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["grid"] = {
        "mode": "grid_gen",
        "generator": {
            "config_file": str(REPO_ROOT / "config" / "grids" / "coarse_square_tube_curved.cfg"),
            "grid_type": "cpipe",
            "da_processors_x": 2,
        },
    }
    case_path = tmp_path / "case_misnested_da.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(case_path),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 1
    assert "unsupported key at grid.generator: 'da_processors_x'" in result.stderr
    assert "da_processors_x" in result.stderr
    assert "This key is valid at: grid, grid.programmatic_settings." in result.stderr


def test_validate_rejects_mis_cased_grid_da_processor_key_with_suggestion(tmp_path):
    """!
    @brief Test that validate rejects mis-cased DMDA layout keys with a same-level suggestion.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["grid"]["da_processors_Z"] = 1
    case_path = tmp_path / "case_mis_cased_da.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(case_path),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 1
    assert "unsupported key at grid: 'da_processors_Z'" in result.stderr
    assert "Did you mean 'da_processors_z'?" in result.stderr


def test_validate_rejects_unknown_solver_key(tmp_path):
    """!
    @brief Test that validate rejects unknown keys in solver.yml.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["mystery_block"] = {"enabled": True}
    solver_path = tmp_path / "solver_unknown.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 1
    assert "unsupported key at <root>: 'mystery_block'" in result.stderr
    assert "mystery_block" in result.stderr


def test_validate_accepts_deprecated_rk4_solver_aliases(tmp_path):
    """!
    @brief Test validation accepts deprecated RK4 solver spellings during the compatibility window.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["strategy"]["momentum_solver"] = "Dual Time Picard RK4"
    jameson_cfg = solver_cfg["momentum_solver"].pop("dual_time_picard_jameson_rk")
    jameson_cfg["rk4_residual_noise_allowance_factor"] = 1.05
    solver_cfg["momentum_solver"]["dual_time_picard_rk4"] = jameson_cfg
    solver_path = tmp_path / "solver_rk4_compat.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 0, result.stderr


def test_dry_run_rejects_unknown_case_key_before_planning(tmp_path):
    """!
    @brief Test that run --dry-run fails on unsupported case keys.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["grid"]["generator"] = {"da_processors_x": 2}
    case_path = tmp_path / "case_unknown_nested.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(case_path),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--dry-run",
        ]
    )

    assert result.returncode == 1
    assert "unsupported key at grid.generator: 'da_processors_x'" in result.stderr
    assert "This key is valid at: grid, grid.programmatic_settings." in result.stderr
    assert "DRY-RUN PLAN" not in result.stdout


def test_validate_rejects_legacy_domain_periodic_flags(tmp_path):
    """!
    @brief Test that periodicity is configured only through paired boundary conditions.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["models"]["domain"]["i_periodic"] = True
    case_path = tmp_path / "case_legacy_periodic_flag.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(case_path),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 1
    assert "unsupported key at models.domain: 'i_periodic'" in result.stderr


def test_statistics_output_artifacts_are_relative_to_run_directory(tmp_path):
    """!
    @brief Test that statistics output artifacts are relative to run directory.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    post_cfg = {
        "statistics_pipeline": {
            "output_prefix": "stats/BrownianStats",
            "tasks": [{"task": "msd"}],
        },
        "io": {
            "output_directory": "viz",
            "output_filename_prefix": "Field",
        },
    }

    stats_paths = picurv.get_post_statistics_output_artifacts(post_cfg, str(tmp_path))

    assert stats_paths == [str((tmp_path / "stats" / "BrownianStats_msd.csv").resolve())]


def test_statistics_output_artifacts_default_under_monitor_output_statistics_dir(tmp_path):
    """!
    @brief Test that bare statistics prefixes default under the monitor output root statistics subdirectory.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    post_cfg = {
        "statistics_pipeline": {
            "tasks": [{"task": "msd"}],
        },
        "io": {
            "output_directory": "viz",
            "output_filename_prefix": "Field",
        },
    }
    monitor_cfg = {
        "io": {
            "directories": {
                "output": "results",
            },
        },
    }

    stats_paths = picurv.get_post_statistics_output_artifacts(post_cfg, str(tmp_path), monitor_cfg)

    assert stats_paths == [str((tmp_path / "results" / "statistics" / "Stats_msd.csv").resolve())]


def test_dry_run_json_reports_predicted_statistics_csv_artifact(tmp_path):
    """!
    @brief Test that dry run json reports predicted statistics csv artifact.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    post_cfg = picurv.read_yaml_file(str(valid / "post.yml"))
    post_cfg["statistics_pipeline"] = {
        "output_prefix": "stats/BrownianStats",
        "tasks": [{"task": "msd"}],
    }

    post_path = tmp_path / "post_with_stats.yml"
    picurv.write_yaml_file(str(post_path), post_cfg)

    result = run_picurv(
        [
            "run",
            "--solve",
            "--post-process",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(post_path),
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    expected_stats_path = str((Path(payload["run_dir_preview"]) / "stats" / "BrownianStats_msd.csv").resolve())
    assert expected_stats_path in payload["artifacts"]


def test_parse_solver_config_maps_uniform_flow_velocity_flags():
    """!
    @brief Test that parse_solver_config maps UNIFORM_FLOW settings into control flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "operation_mode": {
            "eulerian_field_source": "analytical",
            "analytical_type": "UNIFORM_FLOW",
            "uniform_flow": {
                "u": 0.5,
                "v": -0.25,
                "w": 0.125,
            },
        }
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-analytical_type"] == '"UNIFORM_FLOW"'
    assert flags["-analytical_uniform_u"] == 0.5
    assert flags["-analytical_uniform_v"] == -0.25
    assert flags["-analytical_uniform_w"] == 0.125


def test_parse_solver_config_maps_verification_diffusivity_flags():
    """!
    @brief Test that parse_solver_config maps verification diffusivity settings into control flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "operation_mode": {
            "eulerian_field_source": "analytical",
            "analytical_type": "ZERO_FLOW",
        },
        "verification": {
            "sources": {
                "diffusivity": {
                    "mode": "analytical",
                    "profile": "LINEAR_X",
                    "gamma0": 1.0e-3,
                    "slope_x": 2.0e-4,
                }
            }
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-verification_diffusivity_mode"] == '"analytical"'
    assert flags["-verification_diffusivity_profile"] == '"LINEAR_X"'
    assert flags["-verification_diffusivity_gamma0"] == 1.0e-3
    assert flags["-verification_diffusivity_slope_x"] == 2.0e-4


def test_parse_solver_config_maps_verification_scalar_flags():
    """!
    @brief Test that parse_solver_config maps verification scalar settings into control flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "operation_mode": {
            "eulerian_field_source": "analytical",
            "analytical_type": "ZERO_FLOW",
        },
        "verification": {
            "sources": {
                "scalar": {
                    "mode": "analytical",
                    "profile": "SIN_PRODUCT",
                    "amplitude": 2.5,
                    "kx": 3.141592653589793,
                    "ky": 1.5707963267948966,
                    "kz": 0.7853981633974483,
                }
            }
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-verification_scalar_mode"] == '"analytical"'
    assert flags["-verification_scalar_profile"] == '"SIN_PRODUCT"'
    assert flags["-verification_scalar_amplitude"] == 2.5
    assert flags["-verification_scalar_kx"] == 3.141592653589793
    assert flags["-verification_scalar_ky"] == 1.5707963267948966
    assert flags["-verification_scalar_kz"] == 0.7853981633974483


def test_parse_solver_config_maps_solution_convergence_flags():
    """!
    @brief Test that parse_solver_config maps solution-convergence settings into control flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "solution_convergence": {
            "mode": "periodic_deterministic",
            "periodic_deterministic": {
                "period_steps": 12,
            },
        }
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-solution_convergence_mode"] == '"PERIODIC_DETERMINISTIC"'
    assert flags["-solution_convergence_period_steps"] == 12


def test_parse_solver_config_maps_canonical_jameson_controls():
    """!
    @brief Test canonical Jameson solver controls map to canonical C runtime flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Dual Time Picard Jameson RK"},
        "momentum_solver": {
            "dual_time_picard_jameson_rk": {
                "jameson_residual_noise_allowance_factor": 1.07,
            },
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mom_solver_type"] == '"DUALTIME_PICARD_JAMESON_RK"'
    assert flags["-mom_dt_jameson_residual_norm_noise_allowance_factor"] == 1.07


def test_parse_solver_config_emits_newton_krylov_and_preserves_prefixed_petsc_options():
    """! @brief Newton Krylov reuses selection and raw PETSc passthrough paths. """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Newton Krylov"},
        "petsc_passthrough_options": {
            "-mom_nk_snes_rtol": 1.0e-8,
            "-mom_nk_ksp_max_it": 100,
            "-mom_nk_snes_monitor": True,
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mom_solver_type"] == '"newton_krylov"'
    assert flags["-mom_nk_snes_rtol"] == "1e-08"
    assert flags["-mom_nk_ksp_max_it"] == "100"
    assert flags["-mom_nk_snes_monitor"] == "1"


def test_parse_solver_config_rejects_version_one_newton_specific_block():
    """! @brief Version one has fixed operator/PC choices and no YAML solver block. """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Newton Krylov"},
        "momentum_solver": {"newton_krylov": {"preconditioner": "none"}},
    }

    with pytest.raises(ValueError, match="Unsupported momentum_solver"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_maps_ratio_ema_alpha():
    """!
    @brief Test ratio_ema_alpha translates to -mom_ratio_ema_alpha.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Dual Time Picard Jameson RK"},
        "momentum_solver": {
            "dual_time_picard_jameson_rk": {
                "ratio_ema_alpha": 0.5,
            },
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mom_ratio_ema_alpha"] == 0.5


def test_parse_solver_config_ratio_ema_alpha_absent_emits_no_flag():
    """!
    @brief Test ratio_ema_alpha omitted from YAML produces no flag (C default applies).
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Dual Time Picard Jameson RK"},
        "momentum_solver": {"dual_time_picard_jameson_rk": {}},
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert "-mom_ratio_ema_alpha" not in flags


def test_parse_solver_config_maps_optional_momentum_residual_tolerances():
    """!
    @brief Test optional residual tolerances map to the momentum runtime flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Dual Time Picard Jameson RK"},
        "tolerances": {
            "absolute_tol": 1.0e-8,
            "relative_tol": 1.0e-4,
            "residual_absolute_tol": 0.0,
            "residual_relative_tol": 1.0e-3,
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mom_resid_atol"] == 0.0
    assert flags["-mom_resid_rtol"] == 1.0e-3


def test_parse_solver_config_keeps_deprecated_step_tolerance_compatible():
    """!
    @brief Test deprecated step_tol remains readable during its compatibility window.
    """
    picurv = load_picurv_module()
    solver_cfg = {"tolerances": {"step_tol": 1.0e-8}}

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-imp_stol"] == 1.0e-8


@pytest.mark.parametrize(
    ("ema_alpha", "expected_error"),
    [
        (1.5, "ratio_ema_alpha must be in [0, 1]"),
        (-0.1, "ratio_ema_alpha must be in [0, 1]"),
        ("not_a_number", "ratio_ema_alpha must be numeric"),
    ],
)
def test_validate_rejects_invalid_ratio_ema_alpha(tmp_path, ema_alpha, expected_error):
    """!
    @brief Test validation rejects out-of-range or non-numeric ratio_ema_alpha values.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] ema_alpha Invalid alpha value to inject.
    @param[in] expected_error Expected validation diagnostic fragment.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["momentum_solver"]["dual_time_picard_jameson_rk"]["ratio_ema_alpha"] = ema_alpha
    solver_path = tmp_path / "solver_invalid_ema_alpha.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode != 0
    assert expected_error in result.stdout + result.stderr


@pytest.mark.parametrize(
    ("controller_update", "expected_error"),
    [
        ({"growth_factor": 0.99}, "pseudo_cfl.growth_factor must be at least 1"),
        ({"reduction_factor": 1.0}, "pseudo_cfl.reduction_factor must be in (0, 1)"),
        ({"minimum": 0.2, "initial": 0.1, "maximum": 1.0}, "pseudo_cfl requires minimum <= initial <= maximum"),
    ],
)
def test_validate_rejects_invalid_picard_jameson_controller_bounds(
    tmp_path, controller_update, expected_error
):
    """!
    @brief Test validation rejects invalid adaptive pseudo-CFL controller bounds.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] controller_update Invalid pseudo-CFL values injected into the valid fixture.
    @param[in] expected_error Expected validation diagnostic fragment.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["momentum_solver"]["dual_time_picard_jameson_rk"]["pseudo_cfl"].update(controller_update)
    solver_path = tmp_path / "solver_invalid_controller.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 1
    assert expected_error in result.stderr


def test_parse_solver_config_maps_deprecated_rk4_controls_to_jameson_flags():
    """!
    @brief Test deprecated RK4 YAML spellings remain readable but emit canonical Jameson flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Dual Time Picard RK4"},
        "momentum_solver": {
            "dual_time_picard_rk4": {
                "rk4_residual_noise_allowance_factor": 1.08,
            },
        },
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mom_solver_type"] == '"DUALTIME_PICARD_JAMESON_RK"'
    assert flags["-mom_dt_jameson_residual_norm_noise_allowance_factor"] == 1.08


def test_parse_solver_config_rejects_mixed_jameson_and_rk4_aliases():
    """!
    @brief Test canonical and deprecated Jameson solver blocks cannot be supplied together.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "momentum_solver": {
            "dual_time_picard_jameson_rk": {},
            "dual_time_picard_rk4": {},
        },
    }

    with pytest.raises(ValueError, match="do not also set its deprecated dual_time_picard_rk4 alias"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_maps_scalar_transport_flags():
    """!
    @brief Test that scalar transport properties map to C runtime flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "scalar_transport": {
            "schmidt_number": 1.0e12,
            "turbulent_schmidt_number": 0.9,
        }
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-schmidt_number"] == 1.0e12
    assert flags["-turb_schmidt_number"] == 0.9


def test_parse_model_flags_maps_structured_turbulence_options():
    """!
    @brief Test that structured turbulence settings map to C runtime flags.
    """
    picurv = load_picurv_module()
    case_cfg = {
        "models": {
            "physics": {
                "turbulence": {
                    "les": {
                        "enabled": True,
                        "model": "dynamic_smagorinsky",
                        "max_cs": 0.42,
                        "dynamic_frequency": 3,
                        "test_filter": "homogeneous_ik",
                    },
                    "rans": {"enabled": False},
                    "wall_function": {
                        "enabled": True,
                        "model": "log_law",
                        "roughness_height": 1.0e-5,
                    },
                }
            }
        }
    }
    control_lines = []

    picurv.parse_and_add_model_flags(case_cfg, control_lines)

    assert "-les 2" in control_lines
    assert "-max_cs 0.42" in control_lines
    assert "-dynamic_freq 3" in control_lines
    assert "-testfilter_ik 1" in control_lines
    assert "-rans 0" in control_lines
    assert "-wallfunction 1" in control_lines
    assert "-wall_roughness 1e-05" in control_lines


def test_parse_model_flags_preserves_legacy_les_true_constant_smagorinsky():
    """!
    @brief Test that legacy les:true still maps to constant Smagorinsky.
    """
    picurv = load_picurv_module()
    case_cfg = {
        "models": {
            "physics": {
                "turbulence": {
                    "les": True,
                }
            }
        }
    }
    control_lines = []

    picurv.parse_and_add_model_flags(case_cfg, control_lines)

    assert "-les 1" in control_lines


def test_parse_model_flags_allows_minimal_disabled_turbulence_blocks():
    """!
    @brief Test that disabled turbulence blocks do not require model selectors.
    """
    picurv = load_picurv_module()
    case_cfg = {
        "models": {
            "physics": {
                "turbulence": {
                    "les": {"enabled": False},
                    "rans": {"enabled": False},
                    "wall_function": {"enabled": False},
                }
            }
        }
    }
    control_lines = []

    picurv.parse_and_add_model_flags(case_cfg, control_lines)

    assert "-les 0" in control_lines
    assert "-rans 0" in control_lines
    assert "-wallfunction 0" in control_lines


def test_parse_solver_config_maps_structured_poisson_solver_flags():
    """!
    @brief Test that parse_solver_config maps preferred Poisson solver YAML into PETSc flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "poisson_solver": {
            "method": "fgmres",
            "absolute_tolerance": 1.0e-5,
            "relative_tolerance": 1.0e-11,
            "max_iterations": 50,
            "gmres": {"restart": 20},
            "preconditioner": {"type": "multigrid"},
            "multigrid": {
                "levels": 3,
                "pre_sweeps": 2,
                "post_sweeps": 3,
                "cycle": "v",
                "mode": "multiplicative",
                "semi_coarsening": {"i": False, "j": False, "k": True},
                "level_solvers": {
                    "level_0": {"method": "fgmres", "preconditioner": "bjacobi"},
                    "level_1": {"method": "richardson", "preconditioner": "sor"},
                },
            },
        }
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-ps_ksp_type"] == "fgmres"
    assert flags["-ps_ksp_atol"] == 1.0e-5
    assert flags["-poisson_tol"] == 1.0e-5
    assert flags["-ps_ksp_rtol"] == 1.0e-11
    assert flags["-ps_ksp_max_it"] == 50
    assert flags["-ps_ksp_gmres_restart"] == 20
    assert flags["-ps_pc_type"] == "mg"
    assert flags["-mg_level"] == 3
    assert flags["-mg_pre_it"] == 2
    assert flags["-mg_post_it"] == 3
    assert flags["-mg_i_semi"] == "0"
    assert flags["-mg_j_semi"] == "0"
    assert flags["-mg_k_semi"] == "1"
    assert flags["-ps_mg_levels_0_ksp_type"] == "fgmres"
    assert flags["-ps_mg_levels_0_pc_type"] == "bjacobi"
    assert flags["-ps_mg_levels_1_ksp_type"] == "richardson"
    assert flags["-ps_mg_levels_1_pc_type"] == "sor"


def test_parse_solver_config_keeps_legacy_pressure_solver_alias():
    """!
    @brief Test that legacy pressure_solver still maps to Poisson flags.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "pressure_solver": {
            "multigrid": {
                "levels": 4,
                "level_solvers": {
                    "level_0": {"ksp_type": "gmres", "pc_type": "ilu"},
                },
            },
        }
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mg_level"] == 4
    assert flags["-ps_mg_levels_0_ksp_type"] == "gmres"
    assert flags["-ps_mg_levels_0_pc_type"] == "ilu"


def test_parse_solver_config_rejects_non_multigrid_outer_preconditioner():
    """!
    @brief Test that unsupported outer Poisson preconditioners fail clearly.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "poisson_solver": {
            "preconditioner": {"type": "ilu"},
        }
    }

    with pytest.raises(ValueError, match="supports only 'multigrid'"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_rejects_gmres_restart_for_non_gmres_method():
    """!
    @brief Test that GMRES restart is accepted only for GMRES-family methods.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "poisson_solver": {
            "method": "cg",
            "gmres": {"restart": 20},
        }
    }

    with pytest.raises(ValueError, match="gmres.restart is valid only"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_rejects_conflicting_poisson_aliases():
    """!
    @brief Test that preferred and legacy Poisson solver blocks cannot conflict.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "poisson_solver": {"method": "fgmres"},
        "pressure_solver": {"method": "gmres"},
    }

    with pytest.raises(ValueError, match="Both 'poisson_solver' and legacy 'pressure_solver'"):
        picurv.parse_solver_config(solver_cfg)


def test_validate_rejects_verification_diffusivity_for_non_analytical_solver(tmp_path):
    """!
    @brief Test that validate rejects verification diffusivity overrides outside analytical mode.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["verification"] = {
        "sources": {
            "diffusivity": {
                "mode": "analytical",
                "profile": "LINEAR_X",
                "gamma0": 1.0e-3,
                "slope_x": 2.0e-4,
            }
        }
    }

    solver_path = tmp_path / "solver_invalid_verification.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "verification.sources.diffusivity is only valid" in result.stderr




def test_validate_rejects_verification_scalar_for_non_analytical_solver(tmp_path):
    """!
    @brief Test that validate rejects verification scalar overrides outside analytical mode.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["verification"] = {
        "sources": {
            "scalar": {
                "mode": "analytical",
                "profile": "CONSTANT",
                "value": 2.0,
            }
        }
    }

    solver_path = tmp_path / "solver_invalid_scalar_verification.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "verification.sources.scalar is only valid" in result.stderr



def test_validate_rejects_verification_scalar_missing_required_parameter(tmp_path):
    """!
    @brief Test that validate rejects scalar profiles missing required numeric parameters.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["operation_mode"] = {
        "eulerian_field_source": "analytical",
        "analytical_type": "ZERO_FLOW",
    }
    solver_cfg["verification"] = {
        "sources": {
            "scalar": {
                "mode": "analytical",
                "profile": "LINEAR_X",
                "phi0": 1.0,
            }
        }
    }

    solver_path = tmp_path / "solver_invalid_scalar_missing_parameter.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "verification.sources.scalar.slope_x is required" in result.stderr


def test_validate_rejects_solution_convergence_periodic_mode_without_period_steps(tmp_path):
    """!
    @brief Test that validate rejects periodic solution-convergence mode without period_steps.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["solution_convergence"] = {
        "enabled": True,
        "mode": "periodic_deterministic",
    }

    solver_path = tmp_path / "solver_invalid_solution_convergence_periodic.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "solution_convergence.periodic_deterministic.period_steps is required" in result.stderr


def test_validate_rejects_solution_convergence_statistical_block_under_wrong_mode(tmp_path):
    """!
    @brief Test that validate rejects statistical solution-convergence settings under the wrong mode.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["solution_convergence"] = {
        "enabled": True,
        "mode": "steady_deterministic",
        "statistical_steady": {
            "window_steps": 20,
        },
    }

    solver_path = tmp_path / "solver_invalid_solution_convergence_mismatch.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "solution_convergence.statistical_steady is only valid when mode is 'statistical_steady'" in result.stderr




def test_extract_metric_from_csv_supports_p95_and_row_ratios(tmp_path):
    """!
    @brief Test that CSV metric extraction supports p95 reductions and row-wise ratios.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    logs_dir = run_dir / "logs"
    logs_dir.mkdir(parents=True)
    (logs_dir / "search_metrics.csv").write_text(
        """step,search_work_index,migrated,search_population
1,1.0,0,10
2,2.0,5,10
3,4.0,10,20
4,8.0,20,20
""",
        encoding="utf-8",
    )

    p95 = picurv.extract_metric_from_csv(
        str(run_dir),
        {
            "file_glob": "logs/search_metrics.csv",
            "column": "search_work_index",
            "reduction": "p95",
        },
    )
    ratio_mean = picurv.extract_metric_from_csv(
        str(run_dir),
        {
            "file_glob": "logs/search_metrics.csv",
            "numerator_column": "migrated",
            "denominator_column": "search_population",
            "denominator_floor": 1.0,
            "reduction": "mean",
        },
    )

    assert abs(p95 - picurv.np.percentile([1.0, 2.0, 4.0, 8.0], 95.0)) < 1.0e-12
    assert abs(ratio_mean - 0.5) < 1.0e-12


def test_aggregate_study_metrics_supports_parameter_normalization(tmp_path):
    """!
    @brief Test that aggregated study metrics can normalize by the varied parameter space.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    logs_dir = run_dir / "logs"
    results_dir = tmp_path / "results"
    logs_dir.mkdir(parents=True)
    (logs_dir / "search_metrics.csv").write_text(
        """step,lost_cumulative,search_work_index,re_search_fraction,migrated,search_population
1,0,1.0,0.0,0,20
2,1,2.0,0.1,2,20
3,2,4.0,0.2,4,20
4,4,8.0,0.4,8,20
""",
        encoding="utf-8",
    )

    study_cfg = {
        "metrics": [
            {
                "name": "run_loss_fraction",
                "source": "statistics_csv",
                "file_glob": "logs/search_metrics.csv",
                "column": "lost_cumulative",
                "reduction": "last",
                "normalize_by_parameter": "case.models.physics.particles.count",
            },
            {
                "name": "run_swi_p95",
                "source": "statistics_csv",
                "file_glob": "logs/search_metrics.csv",
                "column": "search_work_index",
                "reduction": "p95",
            },
            {
                "name": "mean_migration_fraction",
                "source": "statistics_csv",
                "file_glob": "logs/search_metrics.csv",
                "numerator_column": "migrated",
                "denominator_column": "search_population",
                "denominator_floor": 1.0,
                "reduction": "mean",
            },
        ]
    }
    cases = [
        {
            "case_id": "case_0001",
            "run_dir": str(run_dir),
            "parameters": {
                "case.models.physics.particles.count": 20,
                "solver.operation_mode.uniform_flow.u": 0.1,
            },
        }
    ]

    out_csv = picurv.aggregate_study_metrics(study_cfg, cases, str(results_dir))
    with open(out_csv, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 1
    assert abs(float(rows[0]["run_loss_fraction"]) - 0.2) < 1.0e-12
    assert abs(float(rows[0]["run_swi_p95"]) - picurv.np.percentile([1.0, 2.0, 4.0, 8.0], 95.0)) < 1.0e-12
    assert abs(float(rows[0]["mean_migration_fraction"]) - 0.175) < 1.0e-12
