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
                "zero", "constant", "streamwise_constant", "poiseuille", "ic_gen",
                "spectral_random_velocity"
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
    assert "[OK] Ingress manifest matches the PETSc option scan of its declared sources." in result.stdout
    # Constructed names are only covered when the family check actually ran.
    assert "Scanned option families:" in result.stdout


def test_generate_post_recipe_rejects_non_checkpoint_input_extensions(tmp_path):
    """!
    @brief Test that post recipes reject extensions incompatible with committed bundles.
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

    with pytest.raises(ValueError, match="must be 'dat'"):
        picurv.generate_post_recipe_file(
            str(run_dir),
            "demo_run",
            post_cfg,
            {"Case": "case.yml", "Post-Profile": "post.yml"},
        )


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
    recipe_id = picurv.compute_post_recipe_id(post_cfg)
    assert (
        f"statistics_output_prefix = output/analysis/statistics/{recipe_id}/Stats"
        in content
    )


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

    recipe_id = picurv.compute_post_recipe_id(post_cfg)
    assert stats_paths == [str((
        tmp_path / "output" / "analysis" / "statistics" / recipe_id / "BrownianStats_msd.csv"
    ).resolve())]


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

    recipe_id = picurv.compute_post_recipe_id(post_cfg)
    assert stats_paths == [str(
        (tmp_path / "output" / "analysis" / "statistics" / recipe_id / "Stats_msd.csv").resolve()
    )]


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
    recipe_id = picurv.compute_post_recipe_id(post_cfg)
    expected_stats_path = str((
        Path(payload["run_dir_preview"])
        / "output" / "analysis" / "statistics" / recipe_id / "BrownianStats_msd.csv"
    ).resolve())
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


def test_solution_monitoring_maps_solution_convergence_flags():
    """!
    @brief Test that solution-monitoring ingress maps into the existing C flags.
    """
    picurv = load_picurv_module()
    monitor_cfg = {
        "solution_monitoring": {
            "convergence": {
                "enabled": True,
                "mode": "periodic_deterministic",
                "periodic_deterministic": {"period_steps": 12},
            }
        }
    }

    flags = picurv.resolve_solution_monitoring_flags(monitor_cfg)

    assert flags["-solution_convergence_enabled"] == "true"
    assert flags["-solution_convergence_mode"] == '"PERIODIC_DETERMINISTIC"'
    assert flags["-solution_convergence_period_steps"] == 12


def test_validate_rejects_removed_legacy_averaging_location(tmp_path):
    """!
    @brief Ensure the old case statistics switch cannot silently survive staging.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["models"]["statistics"] = {"time_averaging": True}
    case_path = tmp_path / "case_with_legacy_averaging.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    result = run_picurv([
        "validate",
        "--case", str(case_path),
        "--solver", str(valid / "solver.yml"),
        "--monitor", str(valid / "monitor.yml"),
        "--post", str(valid / "post.yml"),
    ])

    assert result.returncode == 1
    assert "models.statistics' was removed" in result.stderr
    assert "field-statistics pipeline is available" in result.stderr


def test_validate_rejects_old_solver_convergence_location(tmp_path):
    """!
    @brief Ensure physical-solution monitoring uses the monitor configuration.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    solver_cfg["solution_convergence"] = {"mode": "steady_deterministic"}
    solver_path = tmp_path / "solver_with_old_convergence.yml"
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv([
        "validate",
        "--case", str(valid / "case.yml"),
        "--solver", str(solver_path),
        "--monitor", str(valid / "monitor.yml"),
        "--post", str(valid / "post.yml"),
    ])

    assert result.returncode == 1
    assert "solution_convergence' moved" in result.stderr
    assert "solution_monitoring.convergence" in result.stderr


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
    """! @brief Newton Krylov preserves typed raw PETSc passthrough values. """
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
    assert flags["-mom_nk_snes_monitor"] is True


def _complete_newton_krylov_config():
    """!
    @brief Return the complete documented structured Newton--Krylov configuration.
    @return Complete solver configuration mapping.
    """
    return {
        "strategy": {"momentum_solver": "Newton Krylov"},
        "momentum_solver": {
            "newton_krylov": {
                "jacobian": {
                    "type": "finite_difference",
                    "finite_difference": {"mode": "matrix_free"},
                },
                "preconditioner": {"model": "none"},
                "nonlinear_solver": {
                    "method": "newtonls",
                    "absolute_tolerance": 1.0e-10,
                    "relative_tolerance": 1.0e-8,
                    "step_tolerance": 1.0e-12,
                    "max_iterations": 12,
                    "line_search": {"type": "bt"},
                },
                "linear_solver": {
                    "method": "gmres",
                    "absolute_tolerance": 1.0e-10,
                    "relative_tolerance": 1.0e-6,
                    "max_iterations": 400,
                    "gmres": {"restart": 80},
                },
            },
        },
    }


def test_parse_solver_config_maps_complete_structured_newton_krylov_controls(capsys):
    """!
    @brief Complete Newton YAML emits only Newton--Krylov, never pseudo-CFL, controls.
    @param[in] capsys Pytest output-capture fixture.
    """
    picurv = load_picurv_module()
    flags = picurv.parse_solver_config(_complete_newton_krylov_config())
    console_output = capsys.readouterr().out

    generated = {key: value for key, value in flags.items() if key.startswith("-mom_nk_")}
    assert not any("pseudo_cfl" in key for key in flags)
    assert "Momentum Solver: Newton Krylov" in console_output
    assert "pseudo_cfl" not in console_output
    assert generated == {
        "-mom_nk_jacobian_type": "finite_difference",
        "-mom_nk_jacobian_fd_mode": "matrix_free",
        "-mom_nk_preconditioner_model": "none",
        "-mom_nk_preconditioner_structure": "none",
        "-mom_nk_snes_type": "newtonls",
        "-mom_nk_snes_atol": 1.0e-10,
        "-mom_nk_snes_rtol": 1.0e-8,
        "-mom_nk_snes_stol": 1.0e-12,
        "-mom_nk_snes_max_it": 12,
        "-mom_nk_snes_linesearch_type": "bt",
        "-mom_nk_ksp_type": "gmres",
        "-mom_nk_ksp_atol": 1.0e-10,
        "-mom_nk_ksp_rtol": 1.0e-6,
        "-mom_nk_ksp_max_it": 400,
        "-mom_nk_ksp_gmres_restart": 80,
    }
    control_lines = []
    picurv.append_passthrough_flags(control_lines, generated)
    assert control_lines == [
        "-mom_nk_jacobian_type finite_difference",
        "-mom_nk_jacobian_fd_mode matrix_free",
        "-mom_nk_preconditioner_model none",
        "-mom_nk_preconditioner_structure none",
        "-mom_nk_snes_type newtonls",
        "-mom_nk_snes_atol 1e-10",
        "-mom_nk_snes_rtol 1e-08",
        "-mom_nk_snes_stol 1e-12",
        "-mom_nk_snes_max_it 12",
        "-mom_nk_snes_linesearch_type bt",
        "-mom_nk_ksp_type gmres",
        "-mom_nk_ksp_atol 1e-10",
        "-mom_nk_ksp_rtol 1e-06",
        "-mom_nk_ksp_max_it 400",
        "-mom_nk_ksp_gmres_restart 80",
    ]


@pytest.mark.parametrize(
    "newton_block",
    [
        {},
        {"nonlinear_solver": {}},
        {"linear_solver": {}},
        {"linear_solver": {"gmres": {"restart": 30}}},
    ],
)
def test_parse_solver_config_allows_omitted_newton_subfields(newton_block):
    """!
    @brief Omitted Newton controls preserve the C/PETSc defaults.
    @param[in] newton_block Parametrized partial Newton configuration.
    """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Newton Krylov"},
        "momentum_solver": {"newton_krylov": newton_block},
    }

    flags = picurv.parse_solver_config(solver_cfg)

    expected = {
        "-mom_nk_jacobian_type": "finite_difference",
        "-mom_nk_jacobian_fd_mode": "matrix_free",
        "-mom_nk_preconditioner_model": "none",
        "-mom_nk_preconditioner_structure": "none",
    }
    if newton_block.get("linear_solver", {}).get("gmres"):
        expected["-mom_nk_ksp_gmres_restart"] = 30
    assert {key: value for key, value in flags.items() if key.startswith("-mom_nk_")} == expected


def test_parse_solver_config_maps_frozen_momentum_point_block_without_petsc_pc():
    """! @brief The mathematical point-block selection omits a PETSc PC option. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    nk = solver_cfg["momentum_solver"]["newton_krylov"]
    nk["preconditioner"] = {
        "model": "frozen_momentum_jacobian",
        "structure": {"type": "point_block"},
    }
    flags = picurv.parse_solver_config(solver_cfg)
    assert flags["-mom_nk_jacobian_type"] == "finite_difference"
    assert flags["-mom_nk_jacobian_fd_mode"] == "matrix_free"
    assert flags["-mom_nk_preconditioner_model"] == "frozen_momentum_jacobian"
    assert flags["-mom_nk_preconditioner_structure"] == "point_block"
    assert "-mom_nk_pc_type" not in flags


def test_validate_newton_krylov_config_returns_final_discriminated_normalization():
    """! @brief Normalization retains mathematical discriminators, not PETSc details. """
    picurv = load_picurv_module()
    normalized = picurv.validate_newton_krylov_config({})
    assert normalized["jacobian"] == {
        "type": "finite_difference",
        "finite_difference": {"mode": "matrix_free"},
    }
    assert normalized["preconditioner"] == {
        "model": "none",
        "structure": {"type": "none"},
    }


@pytest.mark.parametrize(
    ("jacobian", "message"),
    [
        ({}, "type is required"),
        ({"type": "finite_difference"}, "finite_difference is required"),
        ({"type": "finite_difference", "finite_difference": {}}, "mode is required"),
        (
            {"type": "finite_difference", "finite_difference": {"mode": "colored_sparse"}},
            "colored_sparse.*not implemented",
        ),
        (
            {"type": "frozen_momentum_approximation"},
            "frozen_momentum_approximation.*not implemented",
        ),
        (
            {
                "type": "finite_difference",
                "finite_difference": {"mode": "matrix_free"},
                "frozen_momentum_approximation": {"structure": "diagonal"},
            },
            "unsupported key",
        ),
        (
            {"method": "residual_finite_difference", "representation": "matrix_free"},
            "unsupported key",
        ),
    ],
)
def test_parse_solver_config_rejects_incomplete_future_and_patch_only_jacobians(
        jacobian, message):
    """!
    @brief The Jacobian discriminator accepts only the implemented complete branch.
    @param[in] jacobian Jacobian configuration under test.
    @param[in] message Expected validation diagnostic.
    """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    solver_cfg["momentum_solver"]["newton_krylov"]["jacobian"] = jacobian
    with pytest.raises(ValueError, match=message):
        picurv.parse_solver_config(solver_cfg)


@pytest.mark.parametrize(
    ("model", "structure", "message"),
    [
        ("frozen_momentum_jacobian", None, "requires momentum_solver.newton_krylov.preconditioner.structure.type 'point_block'"),
        ("none", "point_block", "does not accept a matrix structure"),
        ("none", "none", "does not accept a matrix structure"),
    ],
)
def test_parse_solver_config_rejects_unsupported_preconditioner_combinations(
        model, structure, message):
    """!
    @brief Coefficient model and matrix structure are a typed, coupled configuration.
    @param[in] model Preconditioner coefficient model under test.
    @param[in] structure Optional preconditioner matrix structure under test.
    @param[in] message Expected validation diagnostic.
    """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    nk = solver_cfg["momentum_solver"]["newton_krylov"]
    nk["preconditioner"] = {"model": model}
    if structure is not None:
        nk["preconditioner"]["structure"] = {"type": structure}
    with pytest.raises(ValueError, match=message):
        picurv.parse_solver_config(solver_cfg)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (("jacobian", "type", "analytic"), "supports only 'finite_difference'"),
        (("jacobian", "finite_difference", "mode", "assembled"), "supports only 'matrix_free'"),
        (("preconditioner", "model", "diagonal"), "supports only 'none' or 'frozen_momentum_jacobian'"),
    ],
)
def test_parse_solver_config_rejects_unknown_linearization_selectors(
        mutation, message):
    """!
    @brief Unknown mathematical Jacobian and preconditioner selectors fail early.
    @param[in] mutation Configuration path and unknown selector value to apply.
    @param[in] message Expected validation diagnostic.
    """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    target = solver_cfg["momentum_solver"]["newton_krylov"]
    for segment in mutation[:-2]:
        target = target[segment]
    target[mutation[-2]] = mutation[-1]
    with pytest.raises(ValueError, match=message):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_rejects_unknown_preconditioner_structure():
    """! @brief Only the implemented point-block structure is accepted. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    solver_cfg["momentum_solver"]["newton_krylov"]["preconditioner"] = {
        "model": "frozen_momentum_jacobian",
        "structure": {"type": "line"},
    }
    with pytest.raises(ValueError, match="supports only 'none' or 'point_block'"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_accepts_deprecated_none_alias_with_notice():
    """! @brief Released PCNONE YAML remains a warning-producing compatibility alias. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    nk = solver_cfg["momentum_solver"]["newton_krylov"]
    nk.pop("preconditioner")
    nk["linear_solver"]["preconditioner"] = {"type": "none"}
    with pytest.warns(FutureWarning, match="is deprecated"):
        flags = picurv.parse_solver_config(solver_cfg)
    assert flags["-mom_nk_preconditioner_model"] == "none"
    assert flags["-mom_nk_preconditioner_structure"] == "none"
    assert "-mom_nk_pc_type" not in flags


def test_parse_solver_config_rejects_conflicting_old_and_new_preconditioner_forms():
    """! @brief A compatibility `none` cannot override a new non-none model. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    nk = solver_cfg["momentum_solver"]["newton_krylov"]
    nk["preconditioner"] = {
        "model": "frozen_momentum_jacobian",
        "structure": {"type": "point_block"},
    }
    nk["linear_solver"]["preconditioner"] = {"type": "none"}
    with pytest.raises(ValueError, match="conflicts with"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_rejects_uncommitted_low_level_operator_yaml():
    """! @brief The superseded uncommitted operator schema has no compatibility alias. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    solver_cfg["momentum_solver"]["newton_krylov"]["operator"] = {
        "preconditioning_matrix": "center_block_3x3"
    }
    with pytest.raises(ValueError, match="unsupported key"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_rejects_newton_block_for_other_solver():
    """! @brief A Newton-specific block must match the selected solver. """
    picurv = load_picurv_module()
    solver_cfg = {
        "strategy": {"momentum_solver": "Dual Time Picard Jameson RK"},
        "momentum_solver": {"newton_krylov": {}},
    }

    with pytest.raises(ValueError, match="newton_krylov is set but selected solver"):
        picurv.parse_solver_config(solver_cfg)


@pytest.mark.parametrize(
    ("mutation", "error"),
    [
        (("nonlinear_solver", "absolute_tolerance", -1.0), "must be numeric and nonnegative"),
        (("nonlinear_solver", "relative_tolerance", "bad"), "must be numeric and nonnegative"),
        (("nonlinear_solver", "method", ""), "method must be a non-empty string"),
        (("nonlinear_solver", "max_iterations", 0), "must be a positive integer"),
        (("linear_solver", "method", None), "method must be a non-empty string"),
        (("linear_solver", "max_iterations", 1.5), "must be a positive integer"),
        (("linear_solver", "preconditioner", {"type": "jacobi"}), "supports only 'none'"),
    ],
)
def test_parse_solver_config_rejects_invalid_newton_values(mutation, error):
    """!
    @brief Structured Newton tolerances, methods, counts, and PC are validated.
    @param[in] mutation Parametrized nested invalid value mutation.
    @param[in] error Expected validation message fragment.
    """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    section, key, value = mutation
    solver_cfg["momentum_solver"]["newton_krylov"][section][key] = value

    with pytest.raises(ValueError, match=error):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_rejects_newton_restart_for_non_gmres_method():
    """! @brief Newton GMRES restart is limited to the GMRES method family. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    solver_cfg["momentum_solver"]["newton_krylov"]["linear_solver"]["method"] = "cg"

    with pytest.raises(ValueError, match="gmres.restart is valid only"):
        picurv.parse_solver_config(solver_cfg)


@pytest.mark.parametrize(
    ("path", "key"),
    [
        ((), "mystery"),
        (("jacobian",), "mystery"),
        (("jacobian", "finite_difference"), "mystery"),
        (("preconditioner",), "mystery"),
        (("nonlinear_solver",), "mystery"),
        (("nonlinear_solver", "line_search"), "mystery"),
        (("linear_solver",), "mystery"),
        (("linear_solver", "gmres"), "mystery"),
        (("linear_solver", "preconditioner"), "mystery"),
    ],
)
def test_parse_solver_config_rejects_unknown_newton_keys(path, key):
    """!
    @brief Unknown keys are rejected at every Newton solver level.
    @param[in] path Nested path receiving the unknown key.
    @param[in] key Unknown key name.
    """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    newton = solver_cfg["momentum_solver"]["newton_krylov"]
    target = newton
    for segment in path:
        target = target.setdefault(segment, {})
    target[key] = 1

    with pytest.raises(ValueError, match="unsupported key"):
        picurv.parse_solver_config(solver_cfg)


def test_parse_solver_config_newton_passthrough_overrides_structured_values():
    """! @brief Advanced passthrough remains the final Newton option override layer. """
    picurv = load_picurv_module()
    solver_cfg = _complete_newton_krylov_config()
    solver_cfg["petsc_passthrough_options"] = {
        "-mom_nk_snes_rtol": 2.5e-5,
        "-mom_nk_ksp_max_it": 17,
    }

    flags = picurv.parse_solver_config(solver_cfg)

    assert flags["-mom_nk_snes_rtol"] == "2.5e-05"
    assert flags["-mom_nk_ksp_max_it"] == "17"


def test_solver_passthrough_boolean_uses_shared_bare_switch_serialization():
    """! @brief Solver Boolean passthrough matches monitor bare-switch behavior. """
    picurv = load_picurv_module()
    flags = picurv.parse_solver_config({
        "petsc_passthrough_options": {
            "-example_true_switch": True,
            "-example_false_switch": False,
        },
    })
    control_lines = []

    picurv.append_passthrough_flags(control_lines, flags)

    assert "-example_true_switch" in control_lines
    assert not any(line.startswith("-example_false_switch") for line in control_lines)
    assert "-example_true_switch 1" not in control_lines


def test_newton_pipeline_preserves_jameson_and_poisson_generated_controls():
    """! @brief Newton schema additions leave existing Jameson and Poisson mappings unchanged. """
    picurv = load_picurv_module()
    flags = picurv.parse_solver_config({
        "strategy": {"momentum_solver": "Dual Time Picard Jameson RK"},
        "momentum_solver": {
            "dual_time_picard_jameson_rk": {
                "max_pseudo_steps": 21,
                "pseudo_cfl": {"initial": 0.4},
            },
        },
        "poisson_solver": {
            "method": "fgmres",
            "absolute_tolerance": 1.0e-5,
            "gmres": {"restart": 20},
            "preconditioner": {"type": "multigrid"},
        },
    })

    assert {
        key: flags[key]
        for key in (
            "-mom_solver_type", "-mom_max_pseudo_steps", "-pseudo_cfl",
            "-ps_ksp_type", "-ps_ksp_atol", "-poisson_tol",
            "-ps_ksp_gmres_restart", "-ps_pc_type",
        )
    } == {
        "-mom_solver_type": '"DUALTIME_PICARD_JAMESON_RK"',
        "-mom_max_pseudo_steps": 21,
        "-pseudo_cfl": 0.4,
        "-ps_ksp_type": "fgmres",
        "-ps_ksp_atol": 1.0e-5,
        "-poisson_tol": 1.0e-5,
        "-ps_ksp_gmres_restart": 20,
        "-ps_pc_type": "mg",
    }


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
                        "dynamic_frequency": 3,
                        "filter_width": "max_edge",
                        "test_filter": {"kernel": "simpson_ik", "width_ratio": 2.5},
                        "averaging": {"mode": "homogeneous", "directions": ["i", "k"]},
                        "clipping": {
                            "mode": "none",
                            "max_cs": 0.42,
                            "min_viscosity_ratio": 0.25,
                        },
                        "gradient_model": {"enabled": True},
                        "diagnostics": {"enabled": True, "cadence": 5, "yoshizawa_ci": 0.11},
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
    assert "-les_dynamic_frequency 3" in control_lines
    assert "-les_filter_width 2" in control_lines
    assert "-les_test_filter_kernel 1" in control_lines
    assert "-les_test_filter_width_ratio 2.5" in control_lines
    assert "-les_averaging_mode 1" in control_lines
    assert "-les_averaging_directions ik" in control_lines
    assert "-les_clip_mode 2" in control_lines
    assert "-les_clip_max_cs 0.42" in control_lines
    assert "-les_min_viscosity_ratio 0.25" in control_lines
    assert "-les_gradient_model 1" in control_lines
    assert "-les_diagnostics true" in control_lines
    assert "-les_diagnostics_cadence 5" in control_lines
    assert "-les_yoshizawa_ci 0.11" in control_lines
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


def _les_case(les_cfg, periodic_axes="ik"):
    """!
    @brief Builds a minimal case mapping carrying one LES block and chosen periodicity.
    @param[in] les_cfg The `models.physics.turbulence.les` mapping to validate.
    @param[in] periodic_axes Subset of "ijk" whose face pairs are declared PERIODIC.
    @return A case.yml-shaped mapping.
    """
    faces = {"i": ("-Xi", "+Xi"), "j": ("-Eta", "+Eta"), "k": ("-Zeta", "+Zeta")}
    boundary = []
    for axis, pair in faces.items():
        bc_type = "PERIODIC" if axis in periodic_axes else "WALL"
        handler = "geometric" if axis in periodic_axes else "noslip"
        boundary.extend({"face": face, "type": bc_type, "handler": handler} for face in pair)
    return {
        "models": {"physics": {"turbulence": {"les": les_cfg}}},
        "boundary_conditions": boundary,
    }


def _les_validation(les_cfg, periodic_axes="ik"):
    """!
    @brief Runs LES validation over one block and returns its findings.
    @param[in] les_cfg The `models.physics.turbulence.les` mapping to validate.
    @param[in] periodic_axes Subset of "ijk" whose face pairs are declared PERIODIC.
    @return Tuple of (errors, warnings) as lists of strings.
    """
    picurv = load_picurv_module()
    case_cfg = _les_case(les_cfg, periodic_axes)
    errors, warnings = [], []
    picurv.validate_les_configuration(case_cfg, les_cfg, "case.yml", errors, warnings)
    return errors, warnings


def test_les_homogeneous_averaging_defaults_to_the_periodic_axes():
    """!
    @brief Test that homogeneous averaging needs no directions when axes are periodic.
    """
    errors, warnings = _les_validation(
        {"enabled": True, "model": "dynamic_smagorinsky", "averaging": {"mode": "homogeneous"}}
    )
    assert errors == []
    assert warnings == []


def test_les_homogeneous_averaging_without_any_periodic_axis_is_rejected():
    """!
    @brief Test that homogeneous averaging with nothing to derive from is refused.
    """
    errors, _ = _les_validation(
        {"enabled": True, "model": "dynamic_smagorinsky", "averaging": {"mode": "homogeneous"}},
        periodic_axes="",
    )
    assert any("declares none" in error for error in errors)


def test_les_averaging_direction_outside_periodic_axes_warns_but_is_allowed():
    """!
    @brief Test that claiming homogeneity the BCs do not show is the user's to make.
    """
    errors, warnings = _les_validation(
        {
            "enabled": True,
            "model": "dynamic_smagorinsky",
            "averaging": {"mode": "homogeneous", "directions": ["j"]},
        }
    )
    assert errors == []
    assert any("does not declare PERIODIC" in warning for warning in warnings)


def test_les_averaging_directions_are_rejected_outside_homogeneous_mode():
    """!
    @brief Test that local and global averaging refuse an explicit direction list.
    """
    errors, _ = _les_validation(
        {
            "enabled": True,
            "model": "dynamic_smagorinsky",
            "averaging": {"mode": "global", "directions": ["i"]},
        }
    )
    assert any("applies only to mode 'homogeneous'" in error for error in errors)


def test_les_simpson_filter_requires_declared_homogeneous_directions():
    """!
    @brief Test that the Simpson stencil is refused where its assumption does not hold.
    """
    rejected, _ = _les_validation(
        {
            "enabled": True,
            "model": "dynamic_smagorinsky",
            "test_filter": {"kernel": "simpson_ik"},
        },
        periodic_axes="i",
    )
    assert any("simpson_ik" in error for error in rejected)

    accepted, _ = _les_validation(
        {
            "enabled": True,
            "model": "dynamic_smagorinsky",
            "test_filter": {"kernel": "simpson_ik"},
        },
        periodic_axes="ik",
    )
    assert accepted == []


def test_les_test_filter_width_ratio_must_exceed_one():
    """!
    @brief Test that a test filter no wider than the grid filter is refused.
    """
    errors, _ = _les_validation(
        {
            "enabled": True,
            "model": "dynamic_smagorinsky",
            "test_filter": {"width_ratio": 1.0},
        }
    )
    assert any("width_ratio" in error for error in errors)


def test_les_dynamic_only_keys_are_rejected_for_the_constant_model():
    """!
    @brief Test that dynamic-procedure controls are refused for a prescribed coefficient.
    """
    errors, _ = _les_validation(
        {
            "enabled": True,
            "model": "constant_smagorinsky",
            "averaging": {"mode": "local"},
        }
    )
    assert any("cannot be used with model 'constant_smagorinsky'" in error for error in errors)


def test_les_max_cs_outside_clamp_mode_warns_that_it_is_ignored():
    """!
    @brief Test that a ceiling no mode applies is reported rather than silently dropped.
    """
    errors, warnings = _les_validation(
        {
            "enabled": True,
            "model": "dynamic_smagorinsky",
            "clipping": {"mode": "none", "max_cs": 0.3},
        }
    )
    assert errors == []
    assert any("only applied by" in warning for warning in warnings)


def test_les_unknown_selector_values_are_named_in_the_error():
    """!
    @brief Test that a mistyped selector reports the accepted values.
    """
    picurv = load_picurv_module()
    for normalizer, bad in (
        (picurv.normalize_les_filter_width, "cuberoot"),
        (picurv.normalize_les_averaging_mode, "planar"),
        (picurv.normalize_les_clip_mode, "clip"),
        (picurv.normalize_les_test_filter, "homogeneous_ik"),
    ):
        with pytest.raises(ValueError, match="Use one of"):
            normalizer(bad)


def test_les_averaging_directions_reject_repeats_and_unknown_axes():
    """!
    @brief Test that the direction list refuses duplicates and non-axis tokens.
    """
    picurv = load_picurv_module()
    assert picurv.normalize_les_averaging_directions(["k", "i"]) == "ik"
    with pytest.raises(ValueError, match="repeated"):
        picurv.normalize_les_averaging_directions(["i", "i"])
    with pytest.raises(ValueError, match="Unknown LES averaging direction"):
        picurv.normalize_les_averaging_directions(["x"])


def test_les_yoshizawa_constant_is_reachable_from_yaml():
    """!
    @brief Test the Yoshizawa constant is a documented key, not a passthrough-only flag.
    """
    picurv = load_picurv_module()
    control_lines = []
    picurv.append_les_parameter_flags(
        {"diagnostics": {"enabled": True, "yoshizawa_ci": 0.07}}, control_lines)
    assert "-les_yoshizawa_ci 0.07" in control_lines

    errors, _ = _les_validation(
        {"enabled": True, "model": "dynamic_smagorinsky",
         "diagnostics": {"enabled": True, "yoshizawa_ci": -1.0}}
    )
    assert any("yoshizawa_ci" in error for error in errors)


def test_les_test_filter_accepts_a_bare_kernel_name():
    """!
    @brief Test that the scalar shorthand for test_filter still emits the kernel flag.
    """
    picurv = load_picurv_module()
    control_lines = []
    picurv.append_les_parameter_flags({"test_filter": "volume_weighted_box"}, control_lines)
    assert "-les_test_filter_kernel 0" in control_lines


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
    assert flags["-ps_mg_coarse_ksp_type"] == "fgmres"
    assert flags["-ps_mg_coarse_pc_type"] == "bjacobi"
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
    assert flags["-ps_mg_coarse_ksp_type"] == "gmres"
    assert flags["-ps_mg_coarse_pc_type"] == "ilu"


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
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))
    monitor_cfg["solution_monitoring"] = {
        "convergence": {
            "enabled": True,
            "mode": "periodic_deterministic",
        }
    }

    monitor_path = tmp_path / "monitor_invalid_solution_convergence_periodic.yml"
    picurv.write_yaml_file(str(monitor_path), monitor_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(monitor_path),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "solution_monitoring.convergence.periodic_deterministic is required" in result.stderr


def test_validate_rejects_solution_convergence_statistical_block_under_wrong_mode(tmp_path):
    """!
    @brief Test that validate rejects statistical solution-convergence settings under the wrong mode.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    valid = FIXTURES / "valid"
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))
    monitor_cfg["solution_monitoring"] = {
        "convergence": {
            "enabled": True,
            "mode": "steady_deterministic",
            "statistical_steady": {"window_steps": 20},
        }
    }

    monitor_path = tmp_path / "monitor_invalid_solution_convergence_mismatch.yml"
    picurv.write_yaml_file(str(monitor_path), monitor_cfg)

    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(monitor_path),
            "--post",
            str(valid / "post.yml"),
        ]
    )

    assert result.returncode == 1
    assert "solution_monitoring.convergence.statistical_steady is only valid for statistical_steady mode" in result.stderr




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


def test_study_plot_axis_uses_grouped_physical_timestep_not_constant_first_parameter():
    """!
    @brief Timestep studies flatten grouped overrides and choose delta-t as x.
    """
    picurv = load_picurv_module()
    study_cfg = {
        "study_type": "timestep_independence",
        "parameters": {
            "case.models.physics.particles.count": [20000],
            "case.run_control": [
                {"start_step": 0, "total_steps": 400, "dt_physical": 0.005},
                {"start_step": 0, "total_steps": 200, "dt_physical": 0.01},
                {"start_step": 0, "total_steps": 100, "dt_physical": 0.02},
            ],
        },
    }
    rows = [
        {
            "case_id": f"case_{index:04d}",
            "case.models.physics.particles.count": "20000",
            "case.run_control.start_step": "0",
            "case.run_control.total_steps": str(total_steps),
            "case.run_control.dt_physical": str(dt),
            "error": str(error),
        }
        for index, (total_steps, dt, error) in enumerate(
            ((400, 0.005, 0.1), (200, 0.01, 0.2), (100, 0.02, 0.4))
        )
    ]

    assert picurv.get_study_parameter_keys(study_cfg) == [
        "case.models.physics.particles.count",
        "case.run_control.start_step",
        "case.run_control.total_steps",
        "case.run_control.dt_physical",
    ]
    axis = picurv._infer_study_plot_axis(study_cfg, rows)
    assert axis["key"] == "case.run_control.dt_physical"
    assert axis["label"] == "Physical timestep, Δt"
    assert axis["values"] == [0.005, 0.01, 0.02]
    groups = picurv._study_plot_groups(study_cfg, rows, axis, "error")
    assert len(groups) == 1
    assert len(groups[0]["points"]) == 3


def test_grid_study_plot_axis_uses_characteristic_resolution_and_coupled_groups():
    """!
    @brief Grid studies use characteristic resolution and group other varied controls.
    """
    picurv = load_picurv_module()
    grid_keys = (
        "case.grid.programmatic_settings.im",
        "case.grid.programmatic_settings.jm",
        "case.grid.programmatic_settings.km",
    )
    study_cfg = {
        "study_type": "grid_independence",
        "parameter_sets": [
            {grid_keys[0]: 8, grid_keys[1]: 16, grid_keys[2]: 32,
             "case.models.physics.particles.count": 1000},
            {grid_keys[0]: 16, grid_keys[1]: 32, grid_keys[2]: 64,
             "case.models.physics.particles.count": 2000},
        ],
    }
    rows = [
        {"case_id": "case_0000", grid_keys[0]: "8", grid_keys[1]: "16", grid_keys[2]: "32",
         "case.models.physics.particles.count": "1000", "final_error": "0.1"},
        {"case_id": "case_0001", grid_keys[0]: "16", grid_keys[1]: "32", grid_keys[2]: "64",
         "case.models.physics.particles.count": "2000", "final_error": "0.025"},
    ]

    axis = picurv._infer_study_plot_axis(study_cfg, rows)
    assert axis["key"] == "characteristic_grid_resolution"
    assert axis["values"] == pytest.approx([16.0, 32.0])
    groups = picurv._study_plot_groups(study_cfg, rows, axis, "final_error")
    assert [group["label"] for group in groups] == ["Study cases"]
    assert [point[0] for point in groups[0]["points"]] == pytest.approx([16.0, 32.0])
    assert [point[1] for point in groups[0]["points"]] == pytest.approx([0.1, 0.025])


# --- Poisson multigrid coarse-solve contract -----------------------------------
#
# `level_solvers.level_0` is the multigrid COARSE SOLVE, not a smoother, and
# maps to PETSc's `-ps_mg_coarse_*` prefix rather than `-ps_mg_levels_N_*`.
# Because multigrid is used as a preconditioner, level_0 must be a fixed linear
# operator; a Krylov method there makes the preconditioner nonlinear and
# decouples the outer FGMRES tracked residual from the true residual b-Ax, so
# the solver reports convergence on a number that no longer describes the
# constraint. See docs/pages/25_Pressure_Poisson_GMRES_Multigrid.md.

def _iter_poisson_level_zero_settings():
    """!
    @brief Yield every shipped `level_solvers.level_0` mapping with its source path.
    @return Iterator of (path, block name, level_0 mapping) tuples.
    """
    roots = [REPO_ROOT / "config" / "solvers", REPO_ROOT / "examples", REPO_ROOT / "tests" / "smoke"]
    seen = set()
    for root in roots:
        for path in sorted(root.rglob("*.yml")):
            if path in seen:
                continue
            seen.add(path)
            try:
                document = yaml.safe_load(path.read_text())
            except yaml.YAMLError:
                continue
            if not isinstance(document, dict):
                continue
            for block_name in ("poisson_solver", "pressure_solver"):
                block = document.get(block_name)
                if not isinstance(block, dict):
                    continue
                levels = (block.get("multigrid") or {}).get("level_solvers")
                if not isinstance(levels, dict):
                    continue
                level_zero = levels.get("level_0")
                if isinstance(level_zero, dict):
                    yield path, block_name, level_zero


def test_shipped_profiles_do_not_put_a_krylov_method_at_multigrid_level_zero():
    """!
    @brief Verify no shipped solver profile configures a Krylov coarse solve.
    """
    core = load_picurv_module()
    offenders = []
    checked = 0
    for path, block_name, level_zero in _iter_poisson_level_zero_settings():
        checked += 1
        method = str(level_zero.get("method", "")).strip().lower()
        if method in core.KRYLOV_KSP_TYPES:
            offenders.append(f"{path.relative_to(REPO_ROOT)} -> {block_name}: method={method!r}")
    assert checked > 0, "no shipped poisson_solver level_0 settings were found to check"
    assert not offenders, (
        "level_0 is the multigrid coarse solve, not a smoother, and must be a fixed linear "
        "operator. A Krylov method there decouples the tracked residual from b-Ax. Use "
        "{method: preonly, preconditioner: redundant}. Offenders:\n  " + "\n  ".join(offenders)
    )


def test_poisson_level_zero_maps_to_the_petsc_coarse_prefix():
    """!
    @brief Verify level_0 maps to `-ps_mg_coarse_*` and positive levels do not.
    """
    core = load_picurv_module()
    flags = core.parse_solver_config({
        "poisson_solver": {
            "preconditioner": {"type": "multigrid"},
            "multigrid": {
                "levels": 2,
                "level_solvers": {
                    "level_0": {"method": "preonly", "preconditioner": "redundant"},
                    "level_1": {"method": "richardson", "preconditioner": "bjacobi"},
                },
            },
        }
    })
    assert flags["-ps_mg_coarse_ksp_type"] == "preonly"
    assert flags["-ps_mg_coarse_pc_type"] == "redundant"
    assert flags["-ps_mg_levels_1_ksp_type"] == "richardson"
    assert "-ps_mg_levels_0_ksp_type" not in flags


def test_krylov_coarse_solver_warns_but_does_not_fail(capsys):
    """!
    @brief Verify a Krylov level_0 is warned about and still allowed.
    @param[in] capsys Pytest capture fixture used to read the emitted warning.
    @details It stays legitimate at large scale when tolerances are set against
             the true residual, so this must be a warning and not an error.
    """
    core = load_picurv_module()
    flags = core.parse_solver_config({
        "poisson_solver": {
            "preconditioner": {"type": "multigrid"},
            "multigrid": {
                "levels": 2,
                "level_solvers": {"level_0": {"method": "fgmres", "preconditioner": "bjacobi"}},
            },
        }
    })
    assert flags["-ps_mg_coarse_ksp_type"] == "fgmres"
    stderr = capsys.readouterr().err
    assert "level_0" in stderr and "COARSE SOLVE" in stderr


# --- Driven periodic flux handlers ---------------------------------------------

def test_both_driven_periodic_handlers_are_registered():
    """!
    @brief Verify `initial_flux` is a first-class handler beside `constant_flux`.
    @details `initial_flux` derives its target from the initial condition, so
             unlike `constant_flux` it must NOT require `target_flux`.
    """
    core = load_picurv_module()
    registry = core.BC_HANDLER_SPECS
    assert "initial_flux" in registry
    assert registry["initial_flux"]["types"] == {"PERIODIC"}
    assert registry["initial_flux"]["required_params"] == set()
    assert "enforce_seam_flux" in registry["initial_flux"]["optional_params"]
    # The old spelling stays accepted so existing case files keep working.
    assert "apply_trim" in registry["initial_flux"]["optional_params"]
    assert registry["constant_flux"]["required_params"] == {"target_flux"}
