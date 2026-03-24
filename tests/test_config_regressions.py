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


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "scripts" / "picurv"
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
    @brief Load `scripts/picurv` as an importable module for regression tests.
    @return Value returned by `load_picurv_module()`.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_regression_module", str(PICURV))
    spec = importlib.util.spec_from_loader("picurv_regression_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_audit_ingress_manifest_matches_c_ingress():
    """!
    @brief Test that audit ingress manifest matches c ingress.
    """
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "audit_ingress.py")],
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
