"""!
@file test_cli_smoke.py
@brief Pytest smoke and contract coverage for the `picurv` CLI.
"""

import importlib.machinery
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace
from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "scripts" / "picurv"
FIXTURES = REPO_ROOT / "tests" / "fixtures"


def run_picurv(args, cwd=REPO_ROOT, env=None):
    """!
    @brief Run the `picurv` CLI and capture the completed-process result.
    @param[in] args Command-line style argument list supplied to the function.
    @param[in] cwd Working directory override supplied to the function.
    @param[in] env Environment override mapping supplied to the function.
    @return Value returned by `run_picurv()`.
    """
    cmd = [sys.executable, str(PICURV)] + list(args)
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, timeout=60, check=False, env=merged_env)


def load_picurv_module():
    """!
    @brief Load `scripts/picurv` as an importable module for white-box CLI tests.
    @return Value returned by `load_picurv_module()`.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_module", str(PICURV))
    spec = importlib.util.spec_from_loader("picurv_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def write_legacy_1d_grid(path: Path) -> Path:
    """!
    @brief Write a minimal legacy 1D-axis grid payload for conversion tests.
    @param[in] path Filesystem path argument passed to `write_legacy_1d_grid()`.
    @return Value returned by `write_legacy_1d_grid()`.
    """
    path.write_text(
        "\n".join(
            [
                "1",
                "3 2 2",
                "0.0 0 0",
                "0.5 0 0",
                "1.0 0 0",
                "0 0.0 0",
                "0 1.0 0",
                "0 0 0.0",
                "0 0 2.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def write_canonical_picgrid(path: Path, dims=(3, 3, 3)) -> Path:
    """!
    @brief Write a minimal canonical PICGRID payload for file-grid tests.
    @param[in] path Filesystem path argument passed to `write_canonical_picgrid()`.
    @param[in] dims Argument passed to `write_canonical_picgrid()`.
    @return Value returned by `write_canonical_picgrid()`.
    """
    im, jm, km = dims
    lines = ["PICGRID", "1", f"{im} {jm} {km}"]
    for k in range(km):
        for j in range(jm):
            for i in range(im):
                lines.append(f"{float(i):.8e} {float(j):.8e} {float(k):.8e}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_executable(path: Path, content: str) -> Path:
    """!
    @brief Write a small executable script used by CLI subprocess tests.
    @param[in] path Filesystem path argument passed to `write_executable()`.
    @param[in] content Script body argument passed to `write_executable()`.
    @return Value returned by `write_executable()`.
    """
    path.write_text(content, encoding="utf-8")
    path.chmod(0o755)
    return path


def make_fake_scancel_env(tmp_path: Path):
    """!
    @brief Create a fake `scancel` executable and environment override for CLI tests.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `make_fake_scancel_env()`.
    """
    bin_dir = tmp_path / "fake-bin"
    bin_dir.mkdir()
    log_path = tmp_path / "scancel.log"
    write_executable(
        bin_dir / "scancel",
        "\n".join(
            [
                "#!/bin/sh",
                "printf '%s\\n' \"$*\" >> \"$SCANCEL_LOG_PATH\"",
                "exit \"${SCANCEL_EXIT_CODE:-0}\"",
            ]
        )
        + "\n",
    )
    return (
        {
            "PATH": str(bin_dir) + os.pathsep + os.environ.get("PATH", ""),
            "SCANCEL_LOG_PATH": str(log_path),
        },
        log_path,
    )


def make_fake_sbatch_env(tmp_path: Path, mode: str = "ok"):
    """!
    @brief Create a fake `sbatch` executable and environment override for direct module tests.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] mode Argument passed to `make_fake_sbatch_env()`.
    @return Value returned by `make_fake_sbatch_env()`.
    """
    bin_dir = tmp_path / "fake-bin"
    bin_dir.mkdir(exist_ok=True)
    log_path = tmp_path / "sbatch.log"
    write_executable(
        bin_dir / "sbatch",
        "\n".join(
            [
                "#!/bin/sh",
                "printf '%s\\n' \"$*\" >> \"$SBATCH_LOG_PATH\"",
                "if [ \"${SBATCH_MODE:-ok}\" = \"fail\" ]; then",
                "  printf '%s\\n' \"${SBATCH_ERROR_TEXT:-submit failed}\" >&2",
                "  exit \"${SBATCH_EXIT_CODE:-7}\"",
                "fi",
                "case \"$*\" in",
                "  *--dependency=afterok:*) job_id=\"${SBATCH_POST_JOB_ID:-1002}\" ;;",
                "  *) job_id=\"${SBATCH_SOLVE_JOB_ID:-1001}\" ;;",
                "esac",
                "if [ \"${SBATCH_MODE:-ok}\" = \"malformed\" ]; then",
                "  printf '%s\\n' \"${SBATCH_STDOUT_TEXT:-not a valid sbatch response}\"",
                "  exit 0",
                "fi",
                "printf 'Submitted batch job %s\\n' \"$job_id\"",
            ]
        )
        + "\n",
    )
    return (
        {
            "PATH": str(bin_dir) + os.pathsep + os.environ.get("PATH", ""),
            "SBATCH_LOG_PATH": str(log_path),
            "SBATCH_MODE": mode,
        },
        log_path,
    )


def create_summary_run_dir(
    tmp_path: Path,
    name: str = "demo_run",
    particle_console_output_frequency: int = 2,
    include_summary_logs: bool = True,
    include_particle_snapshot: bool = True,
) -> Path:
    """!
    @brief Create a run directory populated with summary-capable artifacts.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] name Argument passed to `create_summary_run_dir()`.
    @param[in] particle_console_output_frequency Argument passed to `create_summary_run_dir()`.
    @param[in] include_summary_logs Argument passed to `create_summary_run_dir()`.
    @param[in] include_particle_snapshot Argument passed to `create_summary_run_dir()`.
    @return Value returned by `create_summary_run_dir()`.
    """
    run_dir = tmp_path / "runs" / name
    config_dir = run_dir / "config"
    logs_dir = run_dir / "logs"
    scheduler_dir = run_dir / "scheduler"
    config_dir.mkdir(parents=True)
    logs_dir.mkdir()
    scheduler_dir.mkdir()

    (config_dir / "monitor.yml").write_text(
        "\n".join(
            [
                "logging:",
                "  verbosity: INFO",
                "profiling:",
                "  timestep_output:",
                "    mode: selected",
                "    functions:",
                "      - AdvanceSimulation",
                "io:",
                "  data_output_frequency: 2",
                f"  particle_console_output_frequency: {particle_console_output_frequency}",
                "  particle_log_interval: 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (config_dir / "case.yml").write_text(
        "\n".join(
            [
                "models:",
                "  physics:",
                "    particles:",
                "      count: 120",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (run_dir / "manifest.json").write_text(
        json.dumps(
            {
                "run_id": name,
                "launch_mode": "local",
                "created_at": "2026-03-10T12:00:00",
            }
        ),
        encoding="utf-8",
    )

    if include_summary_logs:
        (logs_dir / "Continuity_Metrics.log").write_text(
            "\n".join(
                [
                    "Timestep   | Block  | Max Divergence     | Max Divergence Location ([k][j][i]=idx) | Sum(RHS)           | Total Flux In      | Total Flux Out     | Net Flux",
                    "------------------------------------------------------------------------------------------------------------------------------------------",
                    "10         | 0      | 1.0000000000e-03   | ([0][0][0] = 0)                         | 2.0000000000e-04   | 1.0000000000e+00   | 9.9000000000e-01   | 1.0000000000e-02",
                    "10         | 1      | -2.5000000000e-03  | ([1][0][1] = 9)                         | 1.5000000000e-04   | 1.0000000000e+00   | 9.8500000000e-01   | 1.5000000000e-02",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        (logs_dir / "Particle_Metrics.log").write_text(
            "\n".join(
                [
                    "Stage              | Timestep   | Total Ptls   | Lost       | Migrated   | Occupied Cells  | Imbalance  | Mig Passes",
                    "----------------------------------------------------------------------------------------------------------------------------",
                    "AfterAdvection     | 10         | 120          | 2          | 5          | 40              | 1.25       | 2",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        (logs_dir / "Momentum_Solver_Convergence_History_Block_0.log").write_text(
            "Step: 10 | PseudoIter(k): 4| | Pseudo-cfl: 2.5000 |dUk|: 1.000000e-06 | |dUk|/|dUprev|: 1.500000e-01 | |Rk|: 2.000000e-05 | |Rk|/|Rprev|: 3.000000e-02\n",
            encoding="utf-8",
        )
        (logs_dir / "Poisson_Solver_Convergence_History_Block_0.log").write_text(
            "\n".join(
                [
                    "--- Convergence for Timestep 10, Block 0 ---",
                    "ts: 10    | block: 0  | iter: 7   | Unprecond Norm:  1.23456e-04 | True Norm:  7.89012e-05 | Rel Norm:  2.34567e-03",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        (logs_dir / "Profiling_Timestep_Summary.csv").write_text(
            "\n".join(
                [
                    "step,function,calls,step_time_s",
                    "10,AdvanceSimulation,1,0.400000",
                    "10,FlowSolver,1,0.300000",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    if include_particle_snapshot:
        (scheduler_dir / f"{name}_solver.log").write_text(
            "\n".join(
                [
                    "Particle states at step 8:",
                    "| Rank | PID | CellID | Position | Velocity | Weights |",
                    "| 0 | 101 | (1, 2, 2) | (0.0, 0.1, 0.2) | (0.5, 1.0, 1.5) | (0.4, 0.5, 0.6) |",
                    "| 1 | 102 | (4, 5, 5) | (0.3, 0.4, 0.5) | (0.0, 0.5, 0.0) | (0.7, 0.8, 0.9) |",
                    "",
                    "Particle states at step 10:",
                    "| Rank | PID | CellID | Position | Velocity | Weights |",
                    "| 0 | 101 | (1, 2, 3) | (0.1, 0.2, 0.3) | (1.0, 2.0, 2.0) | (0.5, 0.6, 0.7) |",
                    "| 1 | 102 | (4, 5, 6) | (0.4, 0.5, 0.6) | (0.0, 1.0, 0.0) | (0.8, 0.9, 1.0) |",
                    "| 1 | 103 | (4, 5, 6) | (0.7, 0.8, 0.9) | (0.5, 0.5, 0.5) | (1.1, 1.2, 1.3) |",
                    "",
                    "Progress: 10/100",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

    return run_dir


def test_top_level_help_smoke():
    """!
    @brief Test that top level help smoke.
    """
    result = run_picurv(["--help"])
    assert result.returncode == 0
    assert "validate" in result.stdout
    assert "Next commands:" in result.stdout


def test_run_help_smoke():
    """!
    @brief Test that run help smoke.
    """
    result = run_picurv(["run", "--help"])
    assert result.returncode == 0
    assert "--dry-run" in result.stdout
    assert "--format {text,json}" in result.stdout


def test_validate_help_smoke():
    """!
    @brief Test that validate help smoke.
    """
    result = run_picurv(["validate", "--help"])
    assert result.returncode == 0
    assert "--case" in result.stdout
    assert "--strict" in result.stdout


def test_summarize_help_smoke():
    """!
    @brief Test that summarize help smoke.
    """
    result = run_picurv(["summarize", "--help"])
    assert result.returncode == 0
    assert "--run-dir" in result.stdout
    assert "--snapshot-rows" in result.stdout


def test_submit_help_smoke():
    """!
    @brief Test that submit help smoke.
    """
    result = run_picurv(["submit", "--help"])
    assert result.returncode == 0
    assert "--run-dir" in result.stdout
    assert "--study-dir" in result.stdout
    assert "--force" in result.stdout


def test_sweep_help_smoke():
    """!
    @brief Test that sweep help smoke.
    """
    result = run_picurv(["sweep", "--help"])
    assert result.returncode == 0
    assert "--study" in result.stdout
    assert "--cluster" in result.stdout
    assert "--no-submit" in result.stdout
    assert "--study-dir" in result.stdout
    assert "--continue" in result.stdout
    assert "--reaggregate" in result.stdout


def test_cancel_help_smoke():
    """!
    @brief Test that cancel help smoke.
    """
    result = run_picurv(["cancel", "--help"])
    assert result.returncode == 0
    assert "--run-dir" in result.stdout
    assert "--stage {all,solve,post-process}" in result.stdout
    assert "--dry-run" in result.stdout


def test_removed_selector_aliases_are_rejected():
    """!
    @brief Test that removed selector aliases are rejected.
    """
    picurv = load_picurv_module()

    rejected = [
        (picurv.normalize_momentum_solver_type, "DUALTIME_PICARD_RK4"),
        (picurv.normalize_momentum_solver_type, "Dual Time NK Arnoldi"),
        (picurv.normalize_field_init_mode, "0"),
        (picurv.normalize_particle_init_mode, "2"),
        (picurv.normalize_analytical_type, "ZEROFLOW"),
        (picurv.normalize_statistics_task, "ComputeMSD"),
    ]

    for func, value in rejected:
        try:
            func(value)
        except ValueError:
            continue
        raise AssertionError(f"{func.__name__} unexpectedly accepted deprecated input {value!r}")


def test_picgrid_validation_requires_canonical_header(tmp_path):
    """!
    @brief Test that picgrid validation requires canonical header.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source = tmp_path / "legacy.grid"
    source.write_text("1\n3 3 3\n0 0 0\n", encoding="utf-8")
    dest = tmp_path / "out.picgrid"

    try:
        picurv.validate_and_nondimensionalize_picgrid(str(source), str(dest), 1.0)
    except ValueError as exc:
        assert "PICGRID header" in str(exc)
    else:
        raise AssertionError("Headerless grid file unexpectedly passed validation.")


def test_validate_valid_configs_pass():
    """!
    @brief Test that validate valid configs pass.
    """
    valid = FIXTURES / "valid"
    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
            "--cluster",
            str(valid / "cluster.yml"),
            "--study",
            str(valid / "study.yml"),
        ]
    )
    assert result.returncode == 0
    assert "[SUCCESS] Validation completed" in result.stdout


def test_validate_zero_mode_case_allows_omitting_velocity_components(tmp_path):
    """!
    @brief Test that validate zero mode case allows omitting velocity components.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_zero.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace(
        '    mode: "Constant"\n    u_physical: 0.0\n    v_physical: 0.0\n    w_physical: 1.0\n',
        '    mode: "Zero"\n',
    )
    case_path.write_text(case_text, encoding="utf-8")

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

    assert result.returncode == 0
    assert "[SUCCESS] Validation completed" in result.stdout


def test_validate_poiseuille_peak_velocity_option_passes_with_unique_inlet_axis(tmp_path):
    """!
    @brief Test that validate poiseuille peak velocity option passes with unique inlet axis.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_poiseuille_peak.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace(
        '    mode: "Constant"\n    u_physical: 0.0\n    v_physical: 0.0\n    w_physical: 1.0\n',
        '    mode: "Poiseuille"\n    peak_velocity_physical: 1.25\n',
    )
    case_path.write_text(case_text, encoding="utf-8")

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

    assert result.returncode == 0
    assert "[SUCCESS] Validation completed" in result.stdout


def test_validate_initial_condition_mode_must_be_explicit(tmp_path):
    """!
    @brief Test that validate initial condition mode must be explicit.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_missing_mode.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace('    mode: "Constant"\n', "")
    case_path.write_text(case_text, encoding="utf-8")

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
    assert "properties.initial_conditions.mode" in result.stderr


def test_validate_analytical_mode_requires_programmatic_grid(tmp_path):
    """!
    @brief Test that validate analytical mode requires programmatic grid.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_file_grid.yml"
    solver_path = tmp_path / "solver_analytical.yml"

    case_cfg = yaml.safe_load((valid / "case.yml").read_text(encoding="utf-8"))
    case_cfg["grid"] = {
        "mode": "file",
        "source_file": str(valid / "case.yml"),
    }
    case_path.write_text(yaml.safe_dump(case_cfg, sort_keys=False), encoding="utf-8")

    solver_cfg = yaml.safe_load((valid / "solver.yml").read_text(encoding="utf-8"))
    solver_cfg["operation_mode"]["eulerian_field_source"] = "analytical"
    solver_path.write_text(yaml.safe_dump(solver_cfg, sort_keys=False), encoding="utf-8")

    result = run_picurv(
        [
            "validate",
            "--case",
            str(case_path),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )

    assert result.returncode == 1
    assert "grid.mode" in result.stderr
    assert "programmatic_c" in result.stderr


def test_validate_particle_restart_mode_omission_warns_on_restart(tmp_path):
    """!
    @brief Test that validate particle restart mode omission warns on restart.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_restart_warn.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace("  start_step: 0\n", "  start_step: 5\n")
    case_text = case_text.replace("      count: 0\n", "      count: 10\n")
    case_path.write_text(case_text, encoding="utf-8")

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

    assert result.returncode == 0
    assert "default to 'load'" in result.stderr


def test_validate_invalid_case_reports_structured_error():
    """!
    @brief Test that validate invalid case reports structured error.
    """
    valid = FIXTURES / "valid"
    invalid = FIXTURES / "invalid"
    result = run_picurv(
        [
            "validate",
            "--case",
            str(invalid / "case_missing_properties.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
        ]
    )
    assert result.returncode == 1
    assert "ERROR CFG_MISSING_SECTION" in result.stderr
    assert "key=properties" in result.stderr


def test_validate_invalid_cluster_and_study_fail():
    """!
    @brief Test that validate invalid cluster and study fail.
    """
    invalid = FIXTURES / "invalid"
    cluster_result = run_picurv(["validate", "--cluster", str(invalid / "cluster_bad_scheduler.yml")])
    assert cluster_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in cluster_result.stderr

    study_result = run_picurv(["validate", "--study", str(invalid / "study_bad_parameter_key.yml")])
    assert study_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in study_result.stderr


def test_dry_run_does_not_create_run_directories():
    """!
    @brief Test that dry run does not create run directories.
    """
    valid = FIXTURES / "valid"
    runs_dir = REPO_ROOT / "runs"
    before = {p.name for p in runs_dir.iterdir()} if runs_dir.exists() else set()

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
            str(valid / "post.yml"),
            "--dry-run",
        ]
    )
    after = {p.name for p in runs_dir.iterdir()} if runs_dir.exists() else set()

    assert result.returncode == 0
    assert "DRY-RUN PLAN" in result.stdout
    assert before == after


def test_dry_run_json_output_schema():
    """!
    @brief Test that dry run json output schema.
    """
    valid = FIXTURES / "valid"
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
            str(valid / "post.yml"),
            "--dry-run",
            "--format",
            "json",
        ]
    )
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert payload["mode"] == "dry-run"
    assert payload["launch_mode"] == "local"
    assert "solve" in payload["stages"]
    assert "post-process" in payload["stages"]
    assert payload["solver_num_procs_effective"] == 1
    assert payload["post_num_procs_effective"] == 1
    assert payload["stages"]["solve"]["num_procs_effective"] == 1
    assert payload["stages"]["post-process"]["num_procs_effective"] == 1


def test_post_process_run_dir_accepts_null_source_data_mapping(tmp_path):
    """!
    @brief Test that post process run dir accepts null source data mapping.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `test_post_process_run_dir_accepts_null_source_data_mapping()`.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()

    run_dir = tmp_path / "existing_run"
    config_dir = run_dir / "config"
    logs_dir = run_dir / "logs"
    output_dir = run_dir / "output"
    config_dir.mkdir(parents=True)
    logs_dir.mkdir()
    output_dir.mkdir()
    (output_dir / "dummy.dat").write_text("0\n", encoding="utf-8")

    shutil.copy2(valid / "case.yml", config_dir / "case.yml")
    shutil.copy2(valid / "monitor.yml", config_dir / "monitor.yml")
    (config_dir / "existing_run.control").write_text("# control placeholder\n", encoding="utf-8")

    post_cfg = picurv.read_yaml_file(str(valid / "post.yml"))
    post_cfg["source_data"] = None
    post_cfg["io"]["output_directory"] = "viz/null_source_data"
    post_path = tmp_path / "post_null_source_data.yml"
    picurv.write_yaml_file(str(post_path), post_cfg)

    calls = []

    def fake_resolve_runtime_executable(name):
        """!
        @brief Return a stub runtime executable path for wrapper-stage tests.
        @param[in] name Argument passed to `fake_resolve_runtime_executable()`.
        @return Value returned by `fake_resolve_runtime_executable()`.
        """
        return f"/tmp/fake/{name}"

    def fake_execute_command(command, run_dir_arg, log_filename, monitor_cfg=None):
        """!
        @brief Record staged execute-command requests without launching a process.
        @param[in] command Argument passed to `fake_execute_command()`.
        @param[in] run_dir_arg Argument passed to `fake_execute_command()`.
        @param[in] log_filename Argument passed to `fake_execute_command()`.
        @param[in] monitor_cfg Argument passed to `fake_execute_command()`.
        """
        calls.append(
            {
                "command": command,
                "run_dir": run_dir_arg,
                "log_filename": log_filename,
            }
        )

    original_resolve = picurv.resolve_runtime_executable
    original_execute = picurv.execute_command
    picurv.resolve_runtime_executable = fake_resolve_runtime_executable
    picurv.execute_command = fake_execute_command

    try:
        args = SimpleNamespace(
            dry_run=False,
            cluster=None,
            scheduler=None,
            num_procs=1,
            solve=False,
            post_process=True,
            run_dir=str(run_dir),
            post=str(post_path),
            case=None,
            solver=None,
            monitor=None,
            no_submit=False,
        )
        picurv.run_workflow(args)
    finally:
        picurv.resolve_runtime_executable = original_resolve
        picurv.execute_command = original_execute

    assert len(calls) == 1
    assert calls[0]["command"][0].endswith("/postprocessor")
    assert calls[0]["log_filename"] == os.path.join("scheduler", "existing_run_eulerian_data.log")
    assert (config_dir / "post.run").is_file()
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stages_completed_or_submitted"] == ["post-process"]


def test_local_solve_wrapper_log_is_routed_to_scheduler_dir(tmp_path):
    """!
    @brief Test that local solver wrapper output goes to scheduler/ for new runs.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `test_local_solve_wrapper_log_is_routed_to_scheduler_dir()`.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    calls = []

    def fake_resolve_runtime_executable(name):
        """!
        @brief Return a fake runtime executable path.
        @param[in] name Argument passed to `fake_resolve_runtime_executable()`.
        @return Value returned by `fake_resolve_runtime_executable()`.
        """
        return f"/tmp/fake/{name}"

    def fake_execute_command(command, run_dir_arg, log_filename, monitor_cfg=None):
        """!
        @brief Capture execute_command inputs instead of launching.
        @param[in] command Argument passed to `fake_execute_command()`.
        @param[in] run_dir_arg Argument passed to `fake_execute_command()`.
        @param[in] log_filename Argument passed to `fake_execute_command()`.
        @param[in] monitor_cfg Argument passed to `fake_execute_command()`.
        """
        calls.append(
            {
                "command": command,
                "run_dir": run_dir_arg,
                "log_filename": log_filename,
            }
        )

    original_resolve = picurv.resolve_runtime_executable
    original_execute = picurv.execute_command
    original_cwd = os.getcwd()
    picurv.resolve_runtime_executable = fake_resolve_runtime_executable
    picurv.execute_command = fake_execute_command

    try:
        os.chdir(tmp_path)
        args = SimpleNamespace(
            dry_run=False,
            cluster=None,
            scheduler=None,
            num_procs=1,
            solve=True,
            post_process=False,
            run_dir=None,
            post=None,
            case=str(valid / "case.yml"),
            solver=str(valid / "solver.yml"),
            monitor=str(valid / "monitor.yml"),
            no_submit=False,
        )
        picurv.run_workflow(args)
    finally:
        os.chdir(original_cwd)
        picurv.resolve_runtime_executable = original_resolve
        picurv.execute_command = original_execute

    assert len(calls) == 1
    assert calls[0]["command"][0].endswith("/simulator")
    assert calls[0]["log_filename"].startswith("scheduler/")
    assert calls[0]["log_filename"].endswith("_solver.log")
    run_dir = Path(calls[0]["run_dir"])
    assert run_dir.parent == tmp_path / "runs"
    assert (run_dir / "scheduler").is_dir()


def test_execute_command_writes_scheduler_log_and_sets_log_level(tmp_path):
    """!
    @brief Test that execute_command launches a real subprocess and writes scheduler logs.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    exe = write_executable(
        tmp_path / "fake_runtime.sh",
        "\n".join(
            [
                "#!/bin/sh",
                "printf 'PWD=%s\\n' \"$PWD\"",
                "printf 'LOG_LEVEL=%s\\n' \"${LOG_LEVEL:-}\"",
                "printf 'ARGS=%s\\n' \"$*\"",
                "mkdir -p results",
                "printf 'ready\\n'",
            ]
        )
        + "\n",
    )

    picurv.execute_command(
        [str(exe), "--flag", "value"],
        str(run_dir),
        os.path.join("scheduler", "stage.log"),
        {"logging": {"verbosity": "debug"}},
    )

    log_text = (run_dir / "scheduler" / "stage.log").read_text(encoding="utf-8")
    assert f"PWD={run_dir}" in log_text
    assert "LOG_LEVEL=DEBUG" in log_text
    assert "ARGS=--flag value" in log_text
    assert (run_dir / "results").is_dir()


def test_execute_command_basename_log_uses_logs_dir_and_inherits_env(tmp_path, monkeypatch):
    """!
    @brief Test that basename log filenames are routed to `logs/` with inherited env vars.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    exe = write_executable(
        tmp_path / "fake_inherit_runtime.sh",
        "\n".join(
            [
                "#!/bin/sh",
                "printf 'LOG_LEVEL=%s\\n' \"${LOG_LEVEL:-}\"",
            ]
        )
        + "\n",
    )
    monkeypatch.setenv("LOG_LEVEL", "TRACE")

    picurv.execute_command([str(exe)], str(run_dir), "plain.log", None)

    log_text = (run_dir / "logs" / "plain.log").read_text(encoding="utf-8")
    assert "LOG_LEVEL=TRACE" in log_text


def test_execute_command_missing_binary_exits_cleanly(tmp_path):
    """!
    @brief Test that execute_command exits cleanly when the binary cannot be found.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    with pytest.raises(SystemExit) as exc_info:
        picurv.execute_command(
            [str(tmp_path / "does-not-exist")],
            str(run_dir),
            os.path.join("scheduler", "missing.log"),
            {},
        )

    assert exc_info.value.code == 1


def test_summarize_latest_json_reads_existing_runtime_artifacts(tmp_path):
    """!
    @brief Test that summarize aggregates available per-step artifacts into JSON.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_summary_run_dir(tmp_path)

    result = run_picurv(
        [
            "summarize",
            "--run-dir",
            str(run_dir),
            "--latest",
            "--snapshot-rows",
            "1",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["run_id"] == "demo_run"
    assert payload["step"] == 10
    assert payload["continuity"]["available"] is True
    assert payload["continuity"]["max_abs_divergence"] == 2.5e-03
    assert payload["momentum"]["blocks"][0]["pseudo_iterations"] == 4
    assert payload["poisson"]["blocks"][0]["iterations"] == 7
    assert payload["particles"]["total_particles"] == 120
    assert payload["profiling"]["available"] is True
    assert payload["profiling"]["functions"][0]["function"] == "AdvanceSimulation"
    assert payload["particle_snapshot"]["available"] is True
    assert payload["particle_snapshot"]["sampled_rows"] == 3
    assert len(payload["particle_snapshot"]["preview_rows"]) == 1
    assert payload["particle_snapshot"]["cadence"]["particle_console_output_frequency"] == 2
    assert payload["particle_snapshot"]["speed"]["max"] == 3.0
    assert payload["particle_snapshot"]["sampled_distribution"]["unique_cells"] == 2
    assert payload["particle_snapshot"]["sampled_distribution"]["duplicate_cells"] == 1
    assert payload["particle_snapshot"]["sampled_distribution"]["rank_counts"] == {"0": 1, "1": 2}
    assert payload["particle_snapshot"]["weights"]["component_0"]["min"] == 0.5
    assert payload["particle_snapshot"]["checks"]["duplicate_pid_count"] == 0
    assert payload["particle_snapshot"]["top_speeds"][0]["pid"] == 101
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["available"] is True
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["previous_step"] == 8
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["matched_pids"] == 2
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["rank_migrations"] == 0
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["cell_changes"] == 2
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["new_count"] == 1
    assert payload["particle_snapshot"]["delta_from_previous_snapshot"]["gone_count"] == 0


def test_summarize_text_output_renders_human_readable_sections(tmp_path):
    """!
    @brief Test that summarize text output renders the main summary sections.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_summary_run_dir(tmp_path)

    result = run_picurv(["summarize", "--run-dir", str(run_dir), "--latest"], cwd=tmp_path)

    assert result.returncode == 0, result.stderr
    assert "RUN STEP SUMMARY" in result.stdout
    assert "Step           : 10 (latest_available)" in result.stdout
    assert "Particle Snapshot (sampled):" in result.stdout
    assert "Profiling:" in result.stdout


def test_summarize_invalid_step_fails_with_available_step_hint(tmp_path):
    """!
    @brief Test that summarize rejects explicit steps missing from the available artifacts.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_summary_run_dir(tmp_path)

    result = run_picurv(
        ["summarize", "--run-dir", str(run_dir), "--step", "999", "--format", "json"],
        cwd=tmp_path,
    )

    assert result.returncode == 1
    assert "Requested step 999 is not present" in result.stderr
    assert "Available steps include" in result.stderr


def test_summarize_missing_run_artifacts_fails_cleanly(tmp_path):
    """!
    @brief Test that summarize fails when no summary-capable runtime artifacts exist.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_summary_run_dir(tmp_path, include_summary_logs=False, include_particle_snapshot=False)

    result = run_picurv(["summarize", "--run-dir", str(run_dir), "--latest"], cwd=tmp_path)

    assert result.returncode == 1
    assert "No summary-capable run artifacts were found" in result.stderr


def test_summarize_no_snapshot_reports_snapshot_unavailable(tmp_path):
    """!
    @brief Test that summarize reports unavailable particle snapshots when cadence is enabled but no stream log exists.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_summary_run_dir(tmp_path, include_particle_snapshot=False)

    result = run_picurv(
        ["summarize", "--run-dir", str(run_dir), "--latest", "--format", "json"],
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["particle_snapshot"]["available"] is False
    assert payload["profiling"]["available"] is True


def test_dry_run_restart_from_missing_run_dir_fails(tmp_path):
    """!
    @brief Test that --restart-from with a nonexistent directory fails.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()

    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    case_cfg["run_control"]["start_step"] = 3
    solver_cfg.setdefault("operation_mode", {})
    solver_cfg["operation_mode"]["eulerian_field_source"] = "load"

    case_path = tmp_path / "case_restart_missing.yml"
    solver_path = tmp_path / "solver_restart_missing.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)
    picurv.write_yaml_file(str(solver_path), solver_cfg)

    result = run_picurv(
        [
            "run",
            "--solve",
            "--restart-from",
            str(tmp_path / "does_not_exist"),
            "--case",
            str(case_path),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--dry-run",
        ],
        cwd=tmp_path,
    )

    assert result.returncode == 1
    assert "does not exist" in result.stderr


def test_post_process_run_dir_missing_config_inputs_fails(tmp_path):
    """!
    @brief Test that post process run dir missing config inputs fails.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    run_dir = tmp_path / "broken_run"
    (run_dir / "config").mkdir(parents=True)

    result = run_picurv(
        [
            "run",
            "--post-process",
            "--run-dir",
            str(run_dir),
            "--post",
            str(valid / "post.yml"),
        ],
        cwd=tmp_path,
    )

    assert result.returncode == 1
    assert "Could not automatically identify required config files" in result.stderr


def test_programmatic_grid_cell_counts_translate_to_node_counts(tmp_path):
    """!
    @brief Test that programmatic grid cell counts translate to node counts.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    run_dir = tmp_path / "programmatic_grid_translation"
    run_dir.mkdir()
    (run_dir / "config").mkdir()

    control_path = Path(
        picurv.generate_solver_control_file(
            str(run_dir),
            "programmatic_grid_translation",
            {
                "case": case_cfg,
                "solver": solver_cfg,
                "monitor": monitor_cfg,
                "case_path": str(valid / "case.yml"),
                "solver_path": str(valid / "solver.yml"),
                "monitor_path": str(valid / "monitor.yml"),
            },
            1,
            {"whitelist": "whitelist.run", "profile": "profile.run"},
        )
    )
    content = control_path.read_text(encoding="utf-8")

    assert "-im 9" in content
    assert "-jm 9" in content
    assert "-km 17" in content
    assert "-walltime_guard_enabled" not in content


def test_parse_slurm_time_limit_to_seconds_supports_common_formats():
    """!
    @brief Test that Slurm time parsing supports common scheduler formats.
    """
    picurv = load_picurv_module()

    assert picurv.parse_slurm_time_limit_to_seconds("15") == 900
    assert picurv.parse_slurm_time_limit_to_seconds("10:30") == 630
    assert picurv.parse_slurm_time_limit_to_seconds("01:02:03") == 3723
    assert picurv.parse_slurm_time_limit_to_seconds("2-03") == 183600
    assert picurv.parse_slurm_time_limit_to_seconds("2-03:04") == 183840
    assert picurv.parse_slurm_time_limit_to_seconds("2-03:04:05") == 183845


def test_resolve_walltime_guard_policy_defaults_for_slurm_cluster():
    """!
    @brief Test that omitted Slurm walltime guard config resolves to built-in defaults.
    """
    picurv = load_picurv_module()
    policy = picurv.resolve_walltime_guard_policy(
        {
            "scheduler": {"type": "slurm"},
            "resources": {"time": "00:10:00"},
            "execution": {},
        }
    )

    assert policy == {
        "enabled": True,
        "warmup_steps": 10,
        "multiplier": 2.0,
        "min_seconds": 60.0,
        "estimator_alpha": 0.35,
    }


def test_grid_gen_exports_node_counts_from_cell_inputs(tmp_path):
    """!
    @brief Test that grid gen exports node counts from cell inputs.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    output_path = tmp_path / "tiny_grid.picgrid"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "grid.gen"),
            "warp",
            "--ncells-i",
            "2",
            "--ncells-j",
            "2",
            "--ncells-k",
            "2",
            "--no-show-stats",
            "--no-write-vtk",
            "--output",
            str(output_path),
        ],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "PICGRID"
    assert lines[1] == "1"
    assert lines[2] == "3 3 3"


def test_grid_gen_legacy1d_conversion_writes_canonical_picgrid(tmp_path):
    """!
    @brief Test that grid.gen legacy1d converts headerless 1D-axis payload to canonical PICGRID.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    legacy_input = write_legacy_1d_grid(tmp_path / "legacy.grid")
    output_path = tmp_path / "converted.picgrid"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "grid.gen"),
            "legacy1d",
            "--input",
            str(legacy_input),
            "--output",
            str(output_path),
            "--no-show-stats",
            "--no-write-vtk",
        ],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "PICGRID"
    assert lines[1] == "1"
    assert lines[2] == "3 2 2"
    assert len(lines[3:]) == 12
    assert lines[3].split() == ["0.00000000e+00", "0.00000000e+00", "0.00000000e+00"]
    assert lines[-1].split() == ["1.00000000e+00", "1.00000000e+00", "2.00000000e+00"]


def test_generate_solver_control_file_applies_top_level_da_processors_for_file_grid(tmp_path):
    """!
    @brief Test that file-grid mode accepts top-level DMDA processor layout hints.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    grid_file = write_canonical_picgrid(tmp_path / "grid.picgrid")
    case_cfg["grid"] = {
        "mode": "file",
        "source_file": str(grid_file),
        "da_processors_x": 1,
        "da_processors_y": 2,
        "da_processors_z": 2,
    }
    case_path = tmp_path / "case_file_grid.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    run_dir = tmp_path / "run_file_grid"
    (run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(case_path),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(case_path),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        4,
        monitor_files,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-da_processors_x 1" in content
    assert "-da_processors_y 2" in content
    assert "-da_processors_z 2" in content


def test_generate_solver_control_file_applies_top_level_da_processors_for_grid_gen(tmp_path):
    """!
    @brief Test that grid-gen mode accepts top-level DMDA processor layout hints.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    case_cfg["grid"] = {
        "mode": "grid_gen",
        "da_processors_x": 1,
        "da_processors_y": 2,
        "da_processors_z": 2,
        "generator": {
            "config_file": str(REPO_ROOT / "config" / "grids" / "coarse_square_tube_curved.cfg"),
            "grid_type": "cpipe",
        },
    }
    case_path = tmp_path / "case_grid_gen.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    run_dir = tmp_path / "run_grid_gen"
    (run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(case_path),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(case_path),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        4,
        monitor_files,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-da_processors_x 1" in content
    assert "-da_processors_y 2" in content
    assert "-da_processors_z 2" in content


def test_generate_solver_control_file_preserves_legacy_programmatic_da_processors(tmp_path):
    """!
    @brief Test that legacy nested programmatic da_processors remain supported.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    case_cfg["grid"]["programmatic_settings"]["da_processors_x"] = 1
    case_cfg["grid"]["programmatic_settings"]["da_processors_y"] = 2
    case_cfg["grid"]["programmatic_settings"]["da_processors_z"] = 2
    case_path = tmp_path / "case_programmatic_legacy.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    run_dir = tmp_path / "run_programmatic_legacy"
    (run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(case_path),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(case_path),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        4,
        monitor_files,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-da_processors_x 1" in content
    assert "-da_processors_y 2" in content
    assert "-da_processors_z 2" in content


def test_case_local_symlinked_picurv_prefers_local_binaries(tmp_path):
    """!
    @brief Test that case local symlinked picurv prefers local binaries.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_dir = tmp_path / "case_dir"
    case_dir.mkdir()

    for name in ("case.yml", "solver.yml", "monitor.yml", "post.yml"):
        shutil.copy2(valid / name, case_dir / name)

    (case_dir / "picurv").symlink_to(PICURV)
    for exe_name in ("simulator", "postprocessor"):
        exe_path = case_dir / exe_name
        exe_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        exe_path.chmod(0o755)

    result = subprocess.run(
        [
            str(case_dir / "picurv"),
            "run",
            "--solve",
            "--post-process",
            "--case",
            "case.yml",
            "--solver",
            "solver.yml",
            "--monitor",
            "monitor.yml",
            "--post",
            "post.yml",
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=str(case_dir),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["launch_command"][0] == str(case_dir / "simulator")
    assert payload["stages"]["post-process"]["launch_command"][0] == str(case_dir / "postprocessor")


def test_case_local_copied_picurv_prefers_local_binaries(tmp_path):
    """!
    @brief Test that case local copied picurv prefers local binaries.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    case_dir = tmp_path / "copied_case_dir"
    case_dir.mkdir()

    for name in ("case.yml", "solver.yml", "monitor.yml", "post.yml"):
        shutil.copy2(valid / name, case_dir / name)

    shutil.copy2(PICURV, case_dir / "picurv")
    for exe_name in ("simulator", "postprocessor"):
        exe_path = case_dir / exe_name
        exe_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        exe_path.chmod(0o755)

    result = subprocess.run(
        [
            str(case_dir / "picurv"),
            "run",
            "--solve",
            "--post-process",
            "--case",
            "case.yml",
            "--solver",
            "solver.yml",
            "--monitor",
            "monitor.yml",
            "--post",
            "post.yml",
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=str(case_dir),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["launch_command"][0] == str(case_dir / "simulator")
    assert payload["stages"]["post-process"]["launch_command"][0] == str(case_dir / "postprocessor")


def test_init_does_not_copy_any_binaries(tmp_path):
    """!
    @brief Test that init copies only template files, no binaries.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    fake_root = tmp_path / "fake_root"
    template_dir = fake_root / "examples" / "demo_case"
    bin_dir = fake_root / "bin"
    (fake_root / "src").mkdir(parents=True)
    (fake_root / "include").mkdir()
    (fake_root / "scripts").mkdir()
    template_dir.mkdir(parents=True)
    bin_dir.mkdir(parents=True)
    (fake_root / "Makefile").write_text("all:\n\t@echo ok\n", encoding="utf-8")

    (template_dir / "case.yml").write_text("run_control:\n  start_step: 0\n", encoding="utf-8")
    for exe_name in ("picurv", "simulator", "postprocessor"):
        exe_path = bin_dir / exe_name
        exe_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        exe_path.chmod(0o755)

    original_cwd = Path.cwd()
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    try:
        os.chdir(work_dir)
        picurv.init_case(
            SimpleNamespace(
                template_name="demo_case",
                dest_name="demo_out",
                source_root=str(fake_root),
            )
        )
    finally:
        os.chdir(original_cwd)

    out_dir = work_dir / "demo_out"
    assert (out_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).is_file()
    for exe_name in ("picurv", "simulator", "postprocessor"):
        assert not (out_dir / exe_name).exists()


def test_empty_enabled_functions_omits_whitelist_and_uses_c_default(tmp_path):
    """!
    @brief Test that empty enabled functions omits whitelist and uses c default.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)

    source_files = {
        "Case": str(valid / "case.yml"),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)

    assert monitor_files["whitelist"] is None
    assert not (run_dir / "config" / "whitelist.run").exists()
    assert (run_dir / "config" / "profile.run").is_file()

    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(valid / "case.yml"),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        1,
        monitor_files,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-whitelist_config_file" not in content
    assert "-profile_config_file" in content


def test_generate_solver_control_file_converts_legacy_grid_when_enabled(tmp_path):
    """!
    @brief Test that file-grid mode can convert legacy payloads via grid.gen before staging grid.run.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    legacy_input = write_legacy_1d_grid(tmp_path / "legacy_input.grid")
    case_cfg["grid"] = {
        "mode": "file",
        "source_file": str(legacy_input),
        "legacy_conversion": {
            "enabled": True,
            "format": "legacy1d",
        },
    }
    case_path = tmp_path / "case_legacy_file.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    run_dir = tmp_path / "run_legacy_file"
    (run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(case_path),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(case_path),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        1,
        monitor_files,
    )

    grid_run_path = run_dir / "config" / "grid.run"
    converted_path = run_dir / "config" / "grid.converted.picgrid"
    assert grid_run_path.is_file()
    assert converted_path.is_file()
    lines = grid_run_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "PICGRID"
    assert lines[1] == "1"
    assert lines[2] == "3 2 2"
    assert len(lines[3:]) == 12

    content = Path(control_file).read_text(encoding="utf-8")
    assert f"-grid_file {grid_run_path}" in content


def test_particle_console_output_frequency_defaults_to_data_output_frequency(tmp_path):
    """!
    @brief Test that particle console output frequency defaults to data output frequency.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(valid / "case.yml"),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(valid / "case.yml"),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        1,
        monitor_files,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-tio 2" in content
    assert "-particle_console_output_freq 2" in content


def test_validate_rejects_negative_particle_console_output_frequency(tmp_path):
    """!
    @brief Test that validate rejects negative particle console output frequency.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))
    monitor_cfg["io"]["particle_console_output_frequency"] = -1

    monitor_path = tmp_path / "monitor_invalid.yml"
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
        ]
    )

    assert result.returncode == 1
    assert "particle_console_output_frequency" in result.stderr


def test_validate_rejects_unknown_legacy_grid_conversion_format(tmp_path):
    """!
    @brief Test that validate rejects unsupported grid.legacy_conversion.format values.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["grid"] = {
        "mode": "file",
        "source_file": str(valid / "case.yml"),
        "legacy_conversion": {
            "enabled": True,
            "format": "unknown_legacy_mode",
        },
    }
    case_path = tmp_path / "case_invalid_legacy_format.yml"
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
    assert "grid.legacy_conversion.format" in result.stderr


def test_validate_rejects_conflicting_top_level_and_legacy_da_processors(tmp_path):
    """!
    @brief Test that conflicting processor-layout definitions are rejected.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    case_cfg["grid"]["da_processors_x"] = 4
    case_cfg["grid"]["programmatic_settings"]["da_processors_x"] = 2
    case_path = tmp_path / "case_conflicting_da_layout.yml"
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
    assert "conflicts with legacy grid.programmatic_settings.da_processors_x" in result.stderr


def test_new_profiling_config_emits_explicit_timestep_flags(tmp_path):
    """!
    @brief Test that new profiling config emits explicit timestep flags.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))
    monitor_cfg["profiling"] = {
        "timestep_output": {
            "mode": "selected",
            "functions": ["AdvanceSimulation"],
            "file": "timesteps.csv",
        },
        "final_summary": {"enabled": False},
    }

    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(valid / "case.yml"),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(valid / "case.yml"),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        1,
        monitor_files,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-profiling_timestep_mode selected" in content
    assert "-profiling_timestep_file timesteps.csv" in content
    assert "-profiling_final_summary false" in content
    assert "-profile_config_file" in content


def test_restart_from_copies_checkpoint_and_sets_restart_dir(tmp_path):
    """!
    @brief Test that --restart-from copies checkpoint files and sets -restart_dir in control file.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    # Create a fake source run directory with checkpoint files at step 5
    source_run_dir = tmp_path / "old_run"
    (source_run_dir / "config").mkdir(parents=True)
    source_output = source_run_dir / "output" / "eulerian"
    source_output.mkdir(parents=True)
    source_particles = source_run_dir / "output" / "particles"
    source_particles.mkdir(parents=True)

    # Create fake checkpoint files for step 5
    (source_output / "ufield00005_0.dat").write_text("euler_data")
    (source_output / "vfield00005_0.dat").write_text("euler_data")
    (source_particles / "pfield00005_0.dat").write_text("particle_data")

    source_monitor_cfg = {
        "io": {
            "directories": {
                "output": "output",
                "restart": "restart",
            }
        }
    }
    picurv.write_yaml_file(str(source_run_dir / "config" / "monitor.yml"), source_monitor_cfg)

    case_cfg["run_control"]["start_step"] = 5
    case_path = tmp_path / "case_restart.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    # Simulate CLI args
    args = SimpleNamespace(restart_from=str(source_run_dir), continue_run=False, run_dir=None)
    new_run_dir = tmp_path / "new_run"
    new_run_dir.mkdir()

    resolved, is_continue = picurv.resolve_restart_source(
        args, case_cfg, solver_cfg, monitor_cfg, str(new_run_dir)
    )

    assert not is_continue
    assert resolved is not None
    # Verify files were copied to new_run/restart/
    assert (Path(resolved) / "eulerian" / "ufield00005_0.dat").exists()
    assert (Path(resolved) / "particles" / "pfield00005_0.dat").exists()

    # Verify control file uses the resolved path
    (new_run_dir / "config").mkdir(parents=True)
    source_files = {
        "Case": str(case_path),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(new_run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(new_run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(case_path),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        1,
        monitor_files,
        restart_source_dir=resolved,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert f"-restart_dir {resolved}" in content


def test_restart_from_load_mode_uses_direct_reference(tmp_path):
    """!
    @brief Test that --restart-from with eulerian "load" mode uses direct reference (no copy).
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    # Setup source run with all step files for "load" mode
    source_run_dir = tmp_path / "flow_run"
    (source_run_dir / "config").mkdir(parents=True)
    source_output = source_run_dir / "output" / "eulerian"
    source_output.mkdir(parents=True)

    case_cfg["run_control"]["start_step"] = 0
    case_cfg["run_control"]["total_steps"] = 3
    solver_cfg.setdefault("operation_mode", {})
    solver_cfg["operation_mode"]["eulerian_field_source"] = "load"

    # Create step files 0-3
    for step in range(4):
        (source_output / f"ufield{step:05d}_0.dat").write_text("data")

    source_monitor_cfg = {"io": {"directories": {"output": "output"}}}
    picurv.write_yaml_file(str(source_run_dir / "config" / "monitor.yml"), source_monitor_cfg)

    args = SimpleNamespace(restart_from=str(source_run_dir), continue_run=False, run_dir=None)
    new_run_dir = tmp_path / "new_run"
    new_run_dir.mkdir()

    resolved, is_continue = picurv.resolve_restart_source(
        args, case_cfg, solver_cfg, monitor_cfg, str(new_run_dir)
    )

    assert not is_continue
    # Direct reference: points to source output, not a copy
    assert resolved == str(source_run_dir / "output")


def test_restart_from_load_mode_missing_step_files_fails(tmp_path):
    """!
    @brief Test that load mode with missing step files raises error.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    source_run_dir = tmp_path / "flow_run"
    (source_run_dir / "config").mkdir(parents=True)
    source_output = source_run_dir / "output" / "eulerian"
    source_output.mkdir(parents=True)

    case_cfg["run_control"]["start_step"] = 0
    case_cfg["run_control"]["total_steps"] = 5
    solver_cfg.setdefault("operation_mode", {})
    solver_cfg["operation_mode"]["eulerian_field_source"] = "load"

    # Only create steps 0 and 1, missing 2-5
    (source_output / "ufield00000_0.dat").write_text("data")
    (source_output / "ufield00001_0.dat").write_text("data")

    source_monitor_cfg = {"io": {"directories": {"output": "output"}}}
    picurv.write_yaml_file(str(source_run_dir / "config" / "monitor.yml"), source_monitor_cfg)

    args = SimpleNamespace(restart_from=str(source_run_dir), continue_run=False, run_dir=None)

    with pytest.raises(ValueError, match="step file.*missing"):
        picurv.resolve_restart_source(
            args, case_cfg, solver_cfg, monitor_cfg, str(tmp_path / "new_run")
        )


def test_continue_mode_auto_populates_restart(tmp_path):
    """!
    @brief Test that --continue auto-populates restart/ from output/ when restart/ is empty.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    run_dir = tmp_path / "my_run"
    output_dir = run_dir / "output" / "eulerian"
    output_dir.mkdir(parents=True)
    (run_dir / "output" / "particles").mkdir()

    case_cfg["run_control"]["start_step"] = 10

    # Create checkpoint at step 10
    (output_dir / "ufield00010_0.dat").write_text("euler")
    (run_dir / "output" / "particles" / "pfield00010_0.dat").write_text("particle")

    args = SimpleNamespace(restart_from=None, continue_run=True, run_dir=str(run_dir))

    resolved, is_continue = picurv.resolve_restart_source(
        args, case_cfg, solver_cfg, monitor_cfg, str(run_dir)
    )

    assert is_continue
    assert resolved is not None
    # Auto-populated restart/ should contain the checkpoint files
    assert (Path(resolved) / "eulerian" / "ufield00010_0.dat").exists()


def test_continue_mode_prefers_curated_restart(tmp_path):
    """!
    @brief Test that --continue uses curated restart/ over output/ (warm-up-and-discard).
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    run_dir = tmp_path / "my_run"

    # Both output/ and restart/ exist, restart/ has curated files
    (run_dir / "output" / "eulerian").mkdir(parents=True)
    (run_dir / "output" / "eulerian" / "ufield00010_0.dat").write_text("output_version")
    (run_dir / "restart" / "eulerian").mkdir(parents=True)
    (run_dir / "restart" / "eulerian" / "ufield00010_0.dat").write_text("curated_version")

    case_cfg["run_control"]["start_step"] = 10

    args = SimpleNamespace(restart_from=None, continue_run=True, run_dir=str(run_dir))

    resolved, is_continue = picurv.resolve_restart_source(
        args, case_cfg, solver_cfg, monitor_cfg, str(run_dir)
    )

    assert is_continue
    # Should use restart/ (curated) not output/
    assert resolved == str(run_dir / "restart")


def test_continue_mode_sets_continue_mode_flag(tmp_path):
    """!
    @brief Test that --continue produces -continue_mode true in control file.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    run_dir = tmp_path / "my_run"
    (run_dir / "config").mkdir(parents=True)

    source_files = {
        "Case": str(valid / "case.yml"),
        "Solver": str(valid / "solver.yml"),
        "Monitor": str(valid / "monitor.yml"),
    }
    monitor_files = picurv.prepare_monitor_files(str(run_dir), "demo_run", monitor_cfg, source_files)
    control_file = picurv.generate_solver_control_file(
        str(run_dir),
        "demo_run",
        {
            "case": case_cfg,
            "case_path": str(valid / "case.yml"),
            "solver": solver_cfg,
            "solver_path": str(valid / "solver.yml"),
            "monitor": monitor_cfg,
            "monitor_path": str(valid / "monitor.yml"),
        },
        1,
        monitor_files,
        continue_mode=True,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert "-continue_mode true" in content


def test_restart_from_and_continue_mutually_exclusive(tmp_path):
    """!
    @brief Test that --restart-from and --continue together raise an error.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    case_cfg["run_control"]["start_step"] = 5

    args = SimpleNamespace(restart_from="/some/dir", continue_run=True, run_dir="/some/dir")

    with pytest.raises(ValueError, match="mutually exclusive"):
        picurv.resolve_restart_source(args, case_cfg, solver_cfg, monitor_cfg, str(tmp_path))


def test_continue_without_run_dir_fails(tmp_path):
    """!
    @brief Test that --continue without --run-dir raises an error.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    case_cfg["run_control"]["start_step"] = 5

    args = SimpleNamespace(restart_from=None, continue_run=True, run_dir=None)

    with pytest.raises(ValueError, match="--continue requires --run-dir"):
        picurv.resolve_restart_source(args, case_cfg, solver_cfg, monitor_cfg, str(tmp_path))


def test_needs_restart_source_analytical_init_no_source(tmp_path):
    """!
    @brief Test that analytical + init + start_step > 0 does not need restart source (F3).
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))

    case_cfg["run_control"]["start_step"] = 10
    solver_cfg.setdefault("operation_mode", {})
    solver_cfg["operation_mode"]["eulerian_field_source"] = "analytical"
    # Ensure particle restart mode is init (not load)
    case_cfg.setdefault("models", {}).setdefault("physics", {}).setdefault("particles", {})
    case_cfg["models"]["physics"]["particles"]["restart_mode"] = "init"

    assert not picurv.needs_restart_source(case_cfg, solver_cfg)


def test_needs_restart_source_required_when_load_mode():
    """!
    @brief Test that needs_restart_source is True when eulerian is "load".
    """
    picurv = load_picurv_module()
    case_cfg = {"run_control": {"start_step": 0}}
    solver_cfg = {"operation_mode": {"eulerian_field_source": "load"}}
    assert picurv.needs_restart_source(case_cfg, solver_cfg)


def test_needs_restart_source_no_flags_start_step_zero():
    """!
    @brief Test that needs_restart_source is False for standard fresh run.
    """
    picurv = load_picurv_module()
    case_cfg = {"run_control": {"start_step": 0}}
    solver_cfg = {"operation_mode": {"eulerian_field_source": "solve"}}
    assert not picurv.needs_restart_source(case_cfg, solver_cfg)


def test_no_restart_flags_with_restart_needed_fails(tmp_path):
    """!
    @brief Test that start_step > 0 without --restart-from or --continue raises error (F5).
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    case_cfg["run_control"]["start_step"] = 5

    args = SimpleNamespace(restart_from=None, continue_run=False, run_dir=None)

    with pytest.raises(ValueError, match="Restart data required"):
        picurv.resolve_restart_source(args, case_cfg, solver_cfg, monitor_cfg, str(tmp_path))


def test_dry_run_local_solver_vs_post_proc_counts():
    """!
    @brief Test that dry run local solver vs post proc counts.
    """
    valid = FIXTURES / "valid"
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
            str(valid / "post.yml"),
            "--num-procs",
            "8",
            "--dry-run",
            "--format",
            "json",
        ]
    )
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["num_procs_effective"] == 8
    assert payload["stages"]["post-process"]["num_procs_effective"] == 1
    assert "mpiexec" in payload["stages"]["solve"]["launch_command_string"]
    assert "-n 8" in payload["stages"]["solve"]["launch_command_string"]
    assert "mpiexec" not in payload["stages"]["post-process"]["launch_command_string"]


def test_dry_run_local_accepts_launcher_override_for_login_node_runs():
    """!
    @brief Test that local multi-rank dry-run honors launcher override env vars.
    """
    valid = FIXTURES / "valid"
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--num-procs",
            "4",
            "--dry-run",
            "--format",
            "json",
        ],
        env={"PICURV_MPI_LAUNCHER": "mpirun -mca pml ucx -mca btl '^uct,ofi' -mca mtl '^ofi'"},
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    launch_command = payload["stages"]["solve"]["launch_command"]
    assert launch_command[0] == "mpirun"
    assert launch_command[1:10] == ["-mca", "pml", "ucx", "-mca", "btl", "^uct,ofi", "-mca", "mtl", "^ofi"]
    assert launch_command[10:12] == ["-n", "4"]


def test_dry_run_local_reads_shared_runtime_execution_config(tmp_path):
    """!
    @brief Test that local multi-rank dry-run can read .picurv-execution.yml.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-execution.yml").write_text(
        "\n".join(
            [
                "default_execution:",
                "  launcher: \"mpirun\"",
                "  launcher_args:",
                "    - -mca",
                "    - pml",
                "    - ucx",
                "    - -mca",
                "    - btl",
                "    - ^uct,ofi",
                "    - -mca",
                "    - mtl",
                "    - ^ofi",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--num-procs",
            "4",
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    launch_command = payload["stages"]["solve"]["launch_command"]
    assert launch_command[0] == "mpirun"
    assert launch_command[1:10] == ["-mca", "pml", "ucx", "-mca", "btl", "^uct,ofi", "-mca", "mtl", "^ofi"]
    assert launch_command[10:12] == ["-n", "4"]


def test_env_launcher_override_wins_over_shared_runtime_execution_config(tmp_path):
    """!
    @brief Test that env launcher override takes precedence over .picurv-execution.yml.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-execution.yml").write_text(
        "\n".join(
            [
                "default_execution:",
                "  launcher: \"mpiexec\"",
                "  launcher_args:",
                "    - --from-config",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--num-procs",
            "3",
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
        env={"PICURV_MPI_LAUNCHER": "mpirun --from-env"},
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    launch_command = payload["stages"]["solve"]["launch_command"]
    assert launch_command[:4] == ["mpirun", "--from-env", "-n", "3"]


def test_dry_run_local_still_reads_legacy_local_runtime_config(tmp_path):
    """!
    @brief Test that legacy .picurv-local.yml remains supported for local runs.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-local.yml").write_text(
        "\n".join(
            [
                "local_execution:",
                "  launcher: \"mpirun\"",
                "  launcher_args:",
                "    - --legacy-local",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--num-procs",
            "3",
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["launch_command"][:4] == ["mpirun", "--legacy-local", "-n", "3"]


def test_dry_run_cluster_post_is_single_task(tmp_path):
    """!
    @brief Test that dry run cluster post is single task.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    cluster_override = tmp_path / "cluster_tmp_ntasks4.yml"
    cluster_override.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup:",
                "    - \"module purge\"",
                "  launcher: \"srun\"",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
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
            str(valid / "post.yml"),
            "--cluster",
            str(cluster_override),
            "--num-procs",
            "4",
            "--dry-run",
            "--format",
            "json",
        ]
    )
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["num_procs_effective"] == 4
    assert payload["stages"]["post-process"]["num_procs_effective"] == 1
    assert "-n 4" in payload["stages"]["solve"]["launch_command_string"]
    assert "-n 1" in payload["stages"]["post-process"]["launch_command_string"]


def test_dry_run_cluster_accepts_inline_launcher_tokens(tmp_path):
    """!
    @brief Test that cluster dry-run accepts launcher strings with inline site flags.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    cluster_override = tmp_path / "cluster_inline_launcher.yml"
    cluster_override.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup:",
                "    - \"module purge\"",
                "  launcher: \"mpirun -mca pml ucx -mca btl '^uct,ofi' -mca mtl '^ofi'\"",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--cluster",
            str(cluster_override),
            "--dry-run",
            "--format",
            "json",
        ]
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    launch_command = payload["stages"]["solve"]["launch_command"]
    assert launch_command[0] == "mpirun"
    assert launch_command[1:10] == ["-mca", "pml", "ucx", "-mca", "btl", "^uct,ofi", "-mca", "mtl", "^ofi"]
    assert launch_command[10:12] == ["-np", "4"]


def test_dry_run_cluster_reads_shared_runtime_execution_config(tmp_path):
    """!
    @brief Test that cluster dry-run falls back to .picurv-execution.yml when cluster.yml omits launcher tokens.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-execution.yml").write_text(
        "\n".join(
            [
                "default_execution:",
                "  launcher: \"mpirun\"",
                "  launcher_args:",
                "    - --from-shared-site",
                "",
            ]
        ),
        encoding="utf-8",
    )
    cluster_override = tmp_path / "cluster_shared_launcher.yml"
    cluster_override.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup:",
                "    - \"module purge\"",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--cluster",
            str(cluster_override),
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["launch_command"][:4] == ["mpirun", "--from-shared-site", "-np", "4"]


def test_cluster_yml_launcher_args_override_shared_cluster_execution(tmp_path):
    """!
    @brief Test that cluster.yml launcher_args can override shared cluster defaults without redefining launcher.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-execution.yml").write_text(
        "\n".join(
            [
                "default_execution:",
                "  launcher: \"mpiexec\"",
                "  launcher_args:",
                "    - --default-shared",
                "",
                "cluster_execution:",
                "  launcher: \"mpirun\"",
                "  launcher_args:",
                "    - --cluster-shared",
                "",
            ]
        ),
        encoding="utf-8",
    )
    cluster_override = tmp_path / "cluster_override_launcher_args.yml"
    cluster_override.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup:",
                "    - \"module purge\"",
                "  launcher_args:",
                "    - --cluster-explicit",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "run",
            "--solve",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
            "--cluster",
            str(cluster_override),
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stages"]["solve"]["launch_command"][:4] == ["mpirun", "--cluster-explicit", "-np", "4"]


def test_validate_cluster_rejects_unparseable_inline_launcher(tmp_path):
    """!
    @brief Test that cluster validation rejects malformed inline launcher quoting.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    cluster_invalid = tmp_path / "cluster_bad_launcher.yml"
    cluster_invalid.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup: []",
                "  launcher: '\"mpirun -mca pml ucx'",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(["validate", "--cluster", str(cluster_invalid)])
    assert result.returncode == 1
    assert "execution.launcher is not shell-parseable" in result.stderr


def test_validate_cluster_rejects_whitespace_packed_launcher_arg(tmp_path):
    """!
    @brief Test that cluster validation rejects launcher_args items containing embedded whitespace.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    cluster_invalid = tmp_path / "cluster_bad_launcher_args.yml"
    cluster_invalid.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup: []",
                "  launcher: \"mpirun\"",
                "  launcher_args: [\"-mca pml ucx\"]",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(["validate", "--cluster", str(cluster_invalid)])
    assert result.returncode == 1
    assert "split whitespace-separated arguments into separate list items" in result.stderr


def test_validate_cluster_warns_on_sample_placeholder_values(tmp_path):
    """!
    @brief Test that cluster validation warns when sample placeholder values remain in use.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    cluster_placeholder = tmp_path / "cluster_placeholder.yml"
    cluster_placeholder.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"my_project_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: \"user@example.edu\"",
                "  mail_type: \"END,FAIL\"",
                "",
                "execution:",
                "  module_setup: []",
                "  launcher: \"srun\"",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(["validate", "--cluster", str(cluster_placeholder)])
    assert result.returncode == 0, result.stderr
    assert "[WARN]" in result.stderr
    assert "my_project_account" in result.stderr
    assert "user@example.edu" in result.stderr


def test_validate_cluster_rejects_invalid_walltime_guard_values(tmp_path):
    """!
    @brief Test that cluster validation rejects malformed walltime guard settings.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    cluster_invalid = tmp_path / "cluster_bad_walltime_guard.yml"
    cluster_invalid.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup: []",
                "  launcher: \"srun\"",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "  walltime_guard:",
                "    warmup_steps: 0",
                "    multiplier: 6.0",
                "    min_seconds: -5",
                "    estimator_alpha: 0.0",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(["validate", "--cluster", str(cluster_invalid)])
    assert result.returncode == 1
    assert "execution.walltime_guard.warmup_steps must be a positive integer" in result.stderr
    assert "execution.walltime_guard.multiplier must be <= 5.0" in result.stderr
    assert "execution.walltime_guard.min_seconds must be a positive number" in result.stderr
    assert "execution.walltime_guard.estimator_alpha must be in (0, 1]" in result.stderr


def test_validate_rejects_invalid_shared_runtime_execution_config(tmp_path):
    """!
    @brief Test that validate reports malformed .picurv-execution.yml content.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-execution.yml").write_text(
        "\n".join(
            [
                "default_execution:",
                "  launcher: 17",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = run_picurv(
        [
            "validate",
            "--case",
            str(valid / "case.yml"),
            "--solver",
            str(valid / "solver.yml"),
            "--monitor",
            str(valid / "monitor.yml"),
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 1
    assert "runtime_execution" in result.stderr
    assert "default_execution.launcher must be a string" in result.stderr


def test_cluster_no_submit_manifest_and_scripts_use_stage_specific_counts(tmp_path):
    """!
    @brief Test that cluster no submit manifest and scripts use stage specific counts.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    cluster_override = tmp_path / "cluster_tmp_ntasks4.yml"
    cluster_override.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup:",
                "    - \"module purge\"",
                "  launcher: \"srun\"",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )

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
            str(valid / "post.yml"),
            "--cluster",
            str(cluster_override),
            "--num-procs",
            "4",
            "--no-submit",
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 0

    runs_dir = tmp_path / "runs"
    run_dirs = [p for p in runs_dir.iterdir() if p.is_dir()]
    assert len(run_dirs) == 1
    run_dir = run_dirs[0]

    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["num_procs"] == 4
    assert manifest["solver_num_procs"] == 4
    assert manifest["post_num_procs"] == 1

    solver_script = (run_dir / "scheduler" / "solver.sbatch").read_text(encoding="utf-8")
    assert "#SBATCH --ntasks-per-node=4" in solver_script
    assert f"#SBATCH --output={run_dir / 'scheduler' / 'solver_%j.out'}" in solver_script
    assert f"#SBATCH --error={run_dir / 'scheduler' / 'solver_%j.err'}" in solver_script
    assert "export PICURV_JOB_START_EPOCH=$(date +%s)" in solver_script
    assert "export PICURV_WALLTIME_LIMIT_SECONDS=600" in solver_script
    assert "srun -n 4 " in solver_script

    post_script = (run_dir / "scheduler" / "post.sbatch").read_text(encoding="utf-8")
    assert "#SBATCH --nodes=1" in post_script
    assert "#SBATCH --ntasks-per-node=1" in post_script
    assert f"#SBATCH --output={run_dir / 'scheduler' / 'post_%j.out'}" in post_script
    assert f"#SBATCH --error={run_dir / 'scheduler' / 'post_%j.err'}" in post_script
    assert "PICURV_JOB_START_EPOCH" not in post_script
    assert "srun -n 1 " in post_script

    control_file = next((run_dir / "config").glob("*.control"))
    control_text = control_file.read_text(encoding="utf-8")
    assert "-walltime_guard_enabled true" in control_text
    assert "-walltime_guard_warmup_steps 10" in control_text
    assert "-walltime_guard_multiplier 2.0" in control_text
    assert "-walltime_guard_min_seconds 60.0" in control_text
    assert "-walltime_guard_estimator_alpha 0.35" in control_text


def create_staged_run_dir(
    tmp_path: Path,
    name: str = "demo_run",
    launch_mode: str = "slurm",
    solve_meta=None,
    post_meta=None,
    include_solve_script: bool = True,
    include_post_script: bool = True,
) -> Path:
    """!
    @brief Create a staged run directory with scheduler metadata for submit tests.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] name Argument passed to `create_staged_run_dir()`.
    @param[in] launch_mode Argument passed to `create_staged_run_dir()`.
    @param[in] solve_meta Argument passed to `create_staged_run_dir()`.
    @param[in] post_meta Argument passed to `create_staged_run_dir()`.
    @param[in] include_solve_script Argument passed to `create_staged_run_dir()`.
    @param[in] include_post_script Argument passed to `create_staged_run_dir()`.
    @return Value returned by `create_staged_run_dir()`.
    """
    run_dir = tmp_path / "runs" / name
    scheduler_dir = run_dir / "scheduler"
    scheduler_dir.mkdir(parents=True)

    solver_script = scheduler_dir / "solver.sbatch"
    post_script = scheduler_dir / "post.sbatch"
    if include_solve_script:
        solver_script.write_text("#!/bin/bash\n", encoding="utf-8")
    if include_post_script:
        post_script.write_text("#!/bin/bash\n", encoding="utf-8")

    submission = {"launch_mode": launch_mode, "stages": {}}
    if solve_meta is not False:
        stage_meta = {"script": str(solver_script), "submitted": False}
        if solve_meta:
            stage_meta.update(solve_meta)
        submission["stages"]["solve"] = stage_meta
    if post_meta is not False:
        stage_meta = {"script": str(post_script), "submitted": False}
        if post_meta:
            stage_meta.update(post_meta)
        submission["stages"]["post-process"] = stage_meta

    (scheduler_dir / "submission.json").write_text(json.dumps(submission, indent=2) + "\n", encoding="utf-8")
    return run_dir


def create_staged_study_dir(
    tmp_path: Path,
    name: str = "demo_study",
    launch_mode: str = "slurm",
    solve_meta=None,
    post_meta=None,
    include_solve_script: bool = True,
    include_post_script: bool = True,
) -> Path:
    """!
    @brief Create a staged study directory with scheduler metadata for submit tests.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] name Argument passed to `create_staged_study_dir()`.
    @param[in] launch_mode Argument passed to `create_staged_study_dir()`.
    @param[in] solve_meta Argument passed to `create_staged_study_dir()`.
    @param[in] post_meta Argument passed to `create_staged_study_dir()`.
    @param[in] include_solve_script Argument passed to `create_staged_study_dir()`.
    @param[in] include_post_script Argument passed to `create_staged_study_dir()`.
    @return Value returned by `create_staged_study_dir()`.
    """
    study_dir = tmp_path / "studies" / name
    scheduler_dir = study_dir / "scheduler"
    scheduler_dir.mkdir(parents=True)

    solver_script = scheduler_dir / "solver_array.sbatch"
    post_script = scheduler_dir / "post_array.sbatch"
    if include_solve_script:
        solver_script.write_text("#!/bin/bash\n", encoding="utf-8")
    if include_post_script:
        post_script.write_text("#!/bin/bash\n", encoding="utf-8")

    submission = {
        "launch_mode": launch_mode,
        "study_id": name,
        "solver_array": {"script": str(solver_script), "submitted": False},
        "post_array": {"script": str(post_script), "submitted": False},
    }
    if solve_meta:
        submission["solver_array"].update(solve_meta)
    if post_meta:
        submission["post_array"].update(post_meta)

    (scheduler_dir / "submission.json").write_text(json.dumps(submission, indent=2) + "\n", encoding="utf-8")
    (study_dir / "study_manifest.json").write_text(
        json.dumps({"study_id": name, "submission": submission}, indent=2) + "\n",
        encoding="utf-8",
    )
    return study_dir


def test_submit_sbatch_process_boundary_parses_job_id(tmp_path, monkeypatch):
    """!
    @brief Test that submit_sbatch calls a real `sbatch` subprocess and parses its job id.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    picurv = load_picurv_module()
    script_path = tmp_path / "solver.sbatch"
    script_path.write_text("#!/bin/bash\n", encoding="utf-8")
    env, log_path = make_fake_sbatch_env(tmp_path, mode="ok")
    for key, value in env.items():
        monkeypatch.setenv(key, value)

    metadata = picurv.submit_sbatch(str(script_path), dependency="7001")

    assert metadata["job_id"] == "1002"
    assert metadata["command"] == ["sbatch", "--dependency=afterok:7001", str(script_path)]
    assert log_path.read_text(encoding="utf-8").strip() == f"--dependency=afterok:7001 {script_path}"


def test_submit_sbatch_rejects_unparseable_output(tmp_path, monkeypatch):
    """!
    @brief Test that submit_sbatch exits when sbatch output does not contain a job id.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    picurv = load_picurv_module()
    script_path = tmp_path / "solver.sbatch"
    script_path.write_text("#!/bin/bash\n", encoding="utf-8")
    env, _ = make_fake_sbatch_env(tmp_path, mode="malformed")
    for key, value in env.items():
        monkeypatch.setenv(key, value)

    with pytest.raises(SystemExit) as exc_info:
        picurv.submit_sbatch(str(script_path))

    assert exc_info.value.code == 1


def test_submit_sbatch_propagates_nonzero_exit_code(tmp_path, monkeypatch):
    """!
    @brief Test that submit_sbatch propagates sbatch subprocess failures.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    picurv = load_picurv_module()
    script_path = tmp_path / "solver.sbatch"
    script_path.write_text("#!/bin/bash\n", encoding="utf-8")
    env, _ = make_fake_sbatch_env(tmp_path, mode="fail")
    env["SBATCH_EXIT_CODE"] = "9"
    env["SBATCH_ERROR_TEXT"] = "scheduler unavailable"
    for key, value in env.items():
        monkeypatch.setenv(key, value)

    with pytest.raises(SystemExit) as exc_info:
        picurv.submit_sbatch(str(script_path))

    assert exc_info.value.code == 9


def test_submit_run_dir_dry_run_prints_planned_sbatch_calls(tmp_path):
    """!
    @brief Test that submit dry-run prints solver and dependent post sbatch commands.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(tmp_path)
    result = run_picurv(["submit", "--run-dir", str(run_dir), "--dry-run"], cwd=tmp_path)
    assert result.returncode == 0
    assert "solver.sbatch" in result.stdout
    assert "post.sbatch" in result.stdout
    assert "--dependency=afterok:<new solve job id>" in result.stdout


def test_submit_run_stage_solve_only_updates_submission_metadata(tmp_path):
    """!
    @brief Test that submit can launch only the staged solve job.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `test_submit_run_stage_solve_only_updates_submission_metadata()`.
    """
    picurv = load_picurv_module()
    run_dir = create_staged_run_dir(tmp_path)
    calls = []

    def fake_submit_sbatch(script_path, dependency=None):
        """!
        @brief Record staged submissions without calling Slurm.
        @param[in] script_path Argument passed to `fake_submit_sbatch()`.
        @param[in] dependency Argument passed to `fake_submit_sbatch()`.
        @return Value returned by `fake_submit_sbatch()`.
        """
        calls.append({"script": script_path, "dependency": dependency})
        return {
            "command": ["sbatch", script_path],
            "returncode": 0,
            "stdout": "Submitted batch job 501",
            "stderr": "",
            "script": script_path,
            "job_id": "501",
        }

    original_submit = picurv.submit_sbatch
    picurv.submit_sbatch = fake_submit_sbatch
    try:
        picurv.submit_staged_jobs(
            SimpleNamespace(run_dir=str(run_dir), study_dir=None, stage="solve", force=False, dry_run=False)
        )
    finally:
        picurv.submit_sbatch = original_submit

    assert calls == [{"script": str(run_dir / "scheduler" / "solver.sbatch"), "dependency": None}]
    submission = json.loads((run_dir / "scheduler" / "submission.json").read_text(encoding="utf-8"))
    assert submission["stages"]["solve"]["submitted"] is True
    assert submission["stages"]["solve"]["job_id"] == "501"
    assert submission["stages"]["post-process"]["submitted"] is False


def test_submit_run_stage_all_adds_post_dependency_on_new_solve_job(tmp_path):
    """!
    @brief Test that submit all submits solve first and chains post to the new solve job id.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `test_submit_run_stage_all_adds_post_dependency_on_new_solve_job()`.
    """
    picurv = load_picurv_module()
    run_dir = create_staged_run_dir(tmp_path)
    calls = []
    job_ids = iter(["601", "602"])

    def fake_submit_sbatch(script_path, dependency=None):
        """!
        @brief Record staged submissions without calling Slurm.
        @param[in] script_path Argument passed to `fake_submit_sbatch()`.
        @param[in] dependency Argument passed to `fake_submit_sbatch()`.
        @return Value returned by `fake_submit_sbatch()`.
        """
        job_id = next(job_ids)
        calls.append({"script": script_path, "dependency": dependency, "job_id": job_id})
        return {
            "command": ["sbatch", script_path],
            "returncode": 0,
            "stdout": f"Submitted batch job {job_id}",
            "stderr": "",
            "script": script_path,
            "job_id": job_id,
        }

    original_submit = picurv.submit_sbatch
    picurv.submit_sbatch = fake_submit_sbatch
    try:
        picurv.submit_staged_jobs(
            SimpleNamespace(run_dir=str(run_dir), study_dir=None, stage="all", force=False, dry_run=False)
        )
    finally:
        picurv.submit_sbatch = original_submit

    assert calls == [
        {"script": str(run_dir / "scheduler" / "solver.sbatch"), "dependency": None, "job_id": "601"},
        {"script": str(run_dir / "scheduler" / "post.sbatch"), "dependency": "601", "job_id": "602"},
    ]
    submission = json.loads((run_dir / "scheduler" / "submission.json").read_text(encoding="utf-8"))
    assert submission["stages"]["solve"]["job_id"] == "601"
    assert submission["stages"]["post-process"]["job_id"] == "602"
    assert submission["stages"]["post-process"]["dependency"] == "afterok:601"


def test_submit_post_process_requires_recorded_solve_job_id(tmp_path):
    """!
    @brief Test that post-only submit refuses when no recorded solve job id exists.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(tmp_path)
    result = run_picurv(["submit", "--run-dir", str(run_dir), "--stage", "post-process"], cwd=tmp_path)
    assert result.returncode == 1
    assert "requires a recorded solve job id" in result.stderr


def test_submit_refuses_already_submitted_stage_without_force(tmp_path):
    """!
    @brief Test that submit protects against duplicate submission unless forced.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "777"},
    )
    result = run_picurv(["submit", "--run-dir", str(run_dir), "--stage", "solve"], cwd=tmp_path)
    assert result.returncode == 1
    assert "already recorded as submitted" in result.stderr


def test_submit_force_resubmits_recorded_stage(tmp_path):
    """!
    @brief Test that --force allows re-submitting an already-submitted stage.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `test_submit_force_resubmits_recorded_stage()`.
    """
    picurv = load_picurv_module()
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "777"},
    )
    calls = []

    def fake_submit_sbatch(script_path, dependency=None):
        """!
        @brief Record forced re-submission without calling Slurm.
        @param[in] script_path Argument passed to `fake_submit_sbatch()`.
        @param[in] dependency Argument passed to `fake_submit_sbatch()`.
        @return Value returned by `fake_submit_sbatch()`.
        """
        calls.append({"script": script_path, "dependency": dependency})
        return {
            "command": ["sbatch", script_path],
            "returncode": 0,
            "stdout": "Submitted batch job 778",
            "stderr": "",
            "script": script_path,
            "job_id": "778",
        }

    original_submit = picurv.submit_sbatch
    picurv.submit_sbatch = fake_submit_sbatch
    try:
        picurv.submit_staged_jobs(
            SimpleNamespace(run_dir=str(run_dir), study_dir=None, stage="solve", force=True, dry_run=False)
        )
    finally:
        picurv.submit_sbatch = original_submit

    assert calls == [{"script": str(run_dir / "scheduler" / "solver.sbatch"), "dependency": None}]
    submission = json.loads((run_dir / "scheduler" / "submission.json").read_text(encoding="utf-8"))
    assert submission["stages"]["solve"]["job_id"] == "778"
    assert submission["stages"]["solve"]["submitted"] is True


def test_submit_study_dir_updates_submission_and_manifest(tmp_path):
    """!
    @brief Test that study submit updates both scheduler submission.json and study_manifest.json.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @return Value returned by `test_submit_study_dir_updates_submission_and_manifest()`.
    """
    picurv = load_picurv_module()
    study_dir = create_staged_study_dir(tmp_path)
    calls = []
    job_ids = iter(["801", "802"])

    def fake_submit_sbatch(script_path, dependency=None):
        """!
        @brief Record study submissions without calling Slurm.
        @param[in] script_path Argument passed to `fake_submit_sbatch()`.
        @param[in] dependency Argument passed to `fake_submit_sbatch()`.
        @return Value returned by `fake_submit_sbatch()`.
        """
        job_id = next(job_ids)
        calls.append({"script": script_path, "dependency": dependency, "job_id": job_id})
        return {
            "command": ["sbatch", script_path],
            "returncode": 0,
            "stdout": f"Submitted batch job {job_id}",
            "stderr": "",
            "script": script_path,
            "job_id": job_id,
        }

    original_submit = picurv.submit_sbatch
    picurv.submit_sbatch = fake_submit_sbatch
    try:
        picurv.submit_staged_jobs(
            SimpleNamespace(run_dir=None, study_dir=str(study_dir), stage="all", force=False, dry_run=False)
        )
    finally:
        picurv.submit_sbatch = original_submit

    assert calls == [
        {"script": str(study_dir / "scheduler" / "solver_array.sbatch"), "dependency": None, "job_id": "801"},
        {"script": str(study_dir / "scheduler" / "post_array.sbatch"), "dependency": "801", "job_id": "802"},
    ]
    submission = json.loads((study_dir / "scheduler" / "submission.json").read_text(encoding="utf-8"))
    manifest = json.loads((study_dir / "study_manifest.json").read_text(encoding="utf-8"))
    assert submission["solver_array"]["job_id"] == "801"
    assert submission["post_array"]["job_id"] == "802"
    assert submission["post_array"]["dependency"] == "afterok:801"
    assert manifest["submission"]["solver_array"]["job_id"] == "801"
    assert manifest["submission"]["post_array"]["job_id"] == "802"


def test_submit_rejects_non_slurm_or_malformed_submission_metadata(tmp_path):
    """!
    @brief Test that submit fails cleanly for local or missing submission metadata.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    local_run_dir = create_staged_run_dir(tmp_path, name="local_run", launch_mode="local")
    local_result = run_picurv(["submit", "--run-dir", str(local_run_dir), "--dry-run"], cwd=tmp_path)
    assert local_result.returncode == 1
    assert "is not Slurm" in local_result.stderr

    malformed_run_dir = tmp_path / "runs" / "missing_meta"
    (malformed_run_dir / "scheduler").mkdir(parents=True)
    malformed_result = run_picurv(["submit", "--run-dir", str(malformed_run_dir), "--dry-run"], cwd=tmp_path)
    assert malformed_result.returncode == 1
    assert "does not contain scheduler submission metadata" in malformed_result.stderr


def test_cancel_run_dir_dry_run_prints_requested_stage_commands(tmp_path):
    """!
    @brief Test that cancel dry-run prints the recorded scancel commands for requested stages.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "901"},
        post_meta={"submitted": True, "job_id": "902"},
    )

    result = run_picurv(["cancel", "--run-dir", str(run_dir), "--dry-run"], cwd=tmp_path)

    assert result.returncode == 0
    assert "Would run: scancel 901  # stage(s): solve" in result.stdout
    assert "Would run: scancel 902  # stage(s): post-process" in result.stdout
    assert "Dry-run only. No jobs were canceled." in result.stdout


def test_cancel_stage_solve_uses_scancel_stubbed_at_process_boundary(tmp_path):
    """!
    @brief Test that cancel can target the solve stage through a process-boundary scancel stub.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "911"},
        post_meta={"submitted": True, "job_id": "922"},
    )
    env, log_path = make_fake_scancel_env(tmp_path)

    result = run_picurv(
        ["cancel", "--run-dir", str(run_dir), "--stage", "solve"],
        cwd=tmp_path,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert log_path.read_text(encoding="utf-8").strip() == "911"
    assert "Canceled Slurm job 911 for stage(s): solve" in result.stdout


def test_cancel_stage_post_process_uses_scancel_stubbed_at_process_boundary(tmp_path):
    """!
    @brief Test that cancel can target the post-process stage through a process-boundary scancel stub.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "933"},
        post_meta={"submitted": True, "job_id": "944"},
    )
    env, log_path = make_fake_scancel_env(tmp_path)

    result = run_picurv(
        ["cancel", "--run-dir", str(run_dir), "--stage", "post-process"],
        cwd=tmp_path,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert log_path.read_text(encoding="utf-8").strip() == "944"
    assert "Canceled Slurm job 944 for stage(s): post-process" in result.stdout


def test_cancel_stage_all_deduplicates_recorded_job_ids(tmp_path):
    """!
    @brief Test that cancel --stage all de-duplicates shared recorded job ids.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "955"},
        post_meta={"submitted": True, "job_id": "955"},
    )
    env, log_path = make_fake_scancel_env(tmp_path)

    result = run_picurv(["cancel", "--run-dir", str(run_dir)], cwd=tmp_path, env=env)

    assert result.returncode == 0, result.stderr
    assert log_path.read_text(encoding="utf-8").splitlines() == ["955"]
    assert "Canceled Slurm job 955 for stage(s): solve, post-process" in result.stdout


def test_cancel_rejects_non_slurm_or_missing_submission_metadata(tmp_path):
    """!
    @brief Test that cancel fails cleanly for local launch metadata or missing submission metadata.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    local_run_dir = create_staged_run_dir(tmp_path, name="local_cancel", launch_mode="local")
    local_result = run_picurv(["cancel", "--run-dir", str(local_run_dir), "--dry-run"], cwd=tmp_path)
    assert local_result.returncode == 1
    assert "is not Slurm" in local_result.stderr

    malformed_run_dir = tmp_path / "runs" / "missing_cancel_meta"
    (malformed_run_dir / "scheduler").mkdir(parents=True)
    malformed_result = run_picurv(["cancel", "--run-dir", str(malformed_run_dir)], cwd=tmp_path)
    assert malformed_result.returncode == 1
    assert "does not contain scheduler submission metadata" in malformed_result.stderr


def test_cancel_fails_when_requested_stages_have_no_recorded_submitted_job_ids(tmp_path):
    """!
    @brief Test that cancel fails when the requested stages have no submitted job ids to cancel.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": False},
        post_meta={"submitted": True},
    )

    result = run_picurv(["cancel", "--run-dir", str(run_dir), "--stage", "solve"], cwd=tmp_path)

    assert result.returncode == 1
    assert "Skipping stage 'solve': job was generated but not submitted" in result.stdout
    assert "No submitted Slurm job IDs were found" in result.stderr


def test_cancel_run_jobs_module_dry_run_reports_requested_stage_commands(tmp_path, capsys):
    """!
    @brief Test that direct cancel_run_jobs dry-run reports the selected stage commands.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] capsys Pytest output-capture fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "1201"},
        post_meta={"submitted": True, "job_id": "1202"},
    )

    picurv.cancel_run_jobs(SimpleNamespace(run_dir=str(run_dir), stage="all", dry_run=True))

    captured = capsys.readouterr()
    assert "Would run: scancel 1201  # stage(s): solve" in captured.out
    assert "Would run: scancel 1202  # stage(s): post-process" in captured.out
    assert "Dry-run only. No jobs were canceled." in captured.out


def test_cancel_run_jobs_module_uses_scancel_process_boundary(tmp_path, monkeypatch, capsys):
    """!
    @brief Test that direct cancel_run_jobs uses a real `scancel` subprocess.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    @param[in] capsys Pytest output-capture fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "1301"},
        post_meta={"submitted": True, "job_id": "1302"},
    )
    env, log_path = make_fake_scancel_env(tmp_path)
    for key, value in env.items():
        monkeypatch.setenv(key, value)

    picurv.cancel_run_jobs(SimpleNamespace(run_dir=str(run_dir), stage="post-process", dry_run=False))

    captured = capsys.readouterr()
    assert log_path.read_text(encoding="utf-8").strip() == "1302"
    assert "Canceled Slurm job 1302 for stage(s): post-process" in captured.out


def test_cancel_run_jobs_module_surfaces_scancel_failures(tmp_path, monkeypatch):
    """!
    @brief Test that direct cancel_run_jobs exits when `scancel` fails.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    picurv = load_picurv_module()
    run_dir = create_staged_run_dir(
        tmp_path,
        solve_meta={"submitted": True, "job_id": "1401"},
        post_meta={"submitted": True, "job_id": "1402"},
    )
    env, _ = make_fake_scancel_env(tmp_path)
    env["SCANCEL_EXIT_CODE"] = "4"
    for key, value in env.items():
        monkeypatch.setenv(key, value)

    with pytest.raises(SystemExit) as exc_info:
        picurv.cancel_run_jobs(SimpleNamespace(run_dir=str(run_dir), stage="solve", dry_run=False))

    assert exc_info.value.code == 1


def test_sweep_workflow_module_no_submit_writes_study_artifacts(tmp_path):
    """!
    @brief Test that direct sweep_workflow writes case and scheduler artifacts in no-submit mode.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    original_cwd = os.getcwd()

    try:
        os.chdir(tmp_path)
        picurv.sweep_workflow(
            SimpleNamespace(
                study=str(valid / "study.yml"),
                cluster=str(valid / "cluster.yml"),
                no_submit=True,
            )
        )
    finally:
        os.chdir(original_cwd)

    studies_dir = tmp_path / "studies"
    study_dirs = [p for p in studies_dir.iterdir() if p.is_dir()]
    assert len(study_dirs) == 1
    study_dir = study_dirs[0]
    submission = json.loads((study_dir / "scheduler" / "submission.json").read_text(encoding="utf-8"))
    case_index_lines = (study_dir / "scheduler" / "case_index.tsv").read_text(encoding="utf-8").splitlines()

    assert submission["launch_mode"] == "slurm"
    assert submission["no_submit"] is True
    assert submission["solver_array"]["submitted"] is False
    assert submission["post_array"]["submitted"] is False
    assert len(case_index_lines) >= 1
    assert (study_dir / "scheduler" / "solver_array.sbatch").is_file()
    assert (study_dir / "scheduler" / "post_array.sbatch").is_file()
    assert (study_dir / "scheduler" / "metrics_aggregate.sbatch").is_file()


def test_sweep_no_submit_writes_array_stdout_stderr_to_scheduler_dir(tmp_path):
    """!
    @brief Test that Slurm array scripts keep scheduler stdout/stderr under scheduler/.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    valid = FIXTURES / "valid"
    cluster_override = tmp_path / "cluster_tmp_ntasks4.yml"
    cluster_override.write_text(
        "\n".join(
            [
                "scheduler:",
                "  type: slurm",
                "",
                "resources:",
                "  account: \"test_account\"",
                "  partition: \"compute\"",
                "  nodes: 1",
                "  ntasks_per_node: 4",
                "  mem: \"4G\"",
                "  time: \"00:10:00\"",
                "",
                "notifications:",
                "  mail_user: null",
                "  mail_type: null",
                "",
                "execution:",
                "  module_setup:",
                "    - \"module purge\"",
                "  launcher: \"srun\"",
                "  launcher_args: []",
                "  extra_sbatch: {}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = run_picurv(
        [
            "sweep",
            "--study",
            str(valid / "study.yml"),
            "--cluster",
            str(cluster_override),
            "--no-submit",
        ],
        cwd=tmp_path,
    )
    assert result.returncode == 0

    studies_dir = tmp_path / "studies"
    study_dirs = [p for p in studies_dir.iterdir() if p.is_dir()]
    assert len(study_dirs) == 1
    study_dir = study_dirs[0]

    solver_script = (study_dir / "scheduler" / "solver_array.sbatch").read_text(encoding="utf-8")
    assert f"#SBATCH --output={study_dir / 'scheduler' / 'solver_%A_%a.out'}" in solver_script
    assert f"#SBATCH --error={study_dir / 'scheduler' / 'solver_%A_%a.err'}" in solver_script
    assert 'export PICURV_JOB_START_EPOCH=$(date +%s)' in solver_script
    assert "export PICURV_WALLTIME_LIMIT_SECONDS=600" in solver_script

    post_script = (study_dir / "scheduler" / "post_array.sbatch").read_text(encoding="utf-8")
    assert f"#SBATCH --output={study_dir / 'scheduler' / 'post_%A_%a.out'}" in post_script
    assert f"#SBATCH --error={study_dir / 'scheduler' / 'post_%A_%a.err'}" in post_script
    assert "PICURV_JOB_START_EPOCH" not in post_script

    sample_control = next((study_dir / "cases").glob("*/config/*.control"))
    sample_control_text = sample_control.read_text(encoding="utf-8")
    assert "-walltime_guard_enabled true" in sample_control_text
    assert "-walltime_guard_warmup_steps 10" in sample_control_text

    assert not (study_dir / "logs").exists()


def test_markdown_link_checker_passes():
    """!
    @brief Test that markdown link checker passes.
    """
    checker = REPO_ROOT / "scripts" / "check_markdown_links.py"
    result = subprocess.run(
        [sys.executable, str(checker)],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
