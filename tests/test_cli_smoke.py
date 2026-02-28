import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PIC_FLOW = REPO_ROOT / "scripts" / "pic.flow"
FIXTURES = REPO_ROOT / "tests" / "fixtures"


def run_pic_flow(args, cwd=REPO_ROOT):
    cmd = [sys.executable, str(PIC_FLOW)] + list(args)
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, timeout=60, check=False)


def test_top_level_help_smoke():
    result = run_pic_flow(["--help"])
    assert result.returncode == 0
    assert "validate" in result.stdout
    assert "Next commands:" in result.stdout


def test_run_help_smoke():
    result = run_pic_flow(["run", "--help"])
    assert result.returncode == 0
    assert "--dry-run" in result.stdout
    assert "--format {text,json}" in result.stdout


def test_validate_help_smoke():
    result = run_pic_flow(["validate", "--help"])
    assert result.returncode == 0
    assert "--case" in result.stdout
    assert "--strict" in result.stdout


def test_validate_valid_configs_pass():
    valid = FIXTURES / "valid"
    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_zero.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace(
        '    mode: "Constant"\n    u_physical: 0.0\n    v_physical: 0.0\n    w_physical: 1.0\n',
        '    mode: "Zero"\n',
    )
    case_path.write_text(case_text, encoding="utf-8")

    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_poiseuille_peak.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace(
        '    mode: "Constant"\n    u_physical: 0.0\n    v_physical: 0.0\n    w_physical: 1.0\n',
        '    mode: "Poiseuille"\n    peak_velocity_physical: 1.25\n',
    )
    case_path.write_text(case_text, encoding="utf-8")

    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_missing_mode.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace('    mode: "Constant"\n', "")
    case_path.write_text(case_text, encoding="utf-8")

    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_file_grid.yml"
    solver_path = tmp_path / "solver_analytical.yml"

    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace(
        """grid:
  mode: programmatic_c
  programmatic_settings:
    im: 8
    jm: 8
    km: 16
    xMins: 0.0
    xMaxs: 1.0
    yMins: 0.0
    yMaxs: 1.0
    zMins: 0.0
    zMaxs: 4.0
    rxs: 1.0
    rys: 1.0
    rzs: 1.0
    cgrids: 0
""",
        f"""grid:\n  mode: file\n  source_file: "{valid / 'case.yml'}"\n""",
    )
    case_path.write_text(case_text, encoding="utf-8")

    solver_text = (valid / "solver.yml").read_text(encoding="utf-8")
    solver_text = solver_text.replace('  eulerian_field_source: "solve"\n', '  eulerian_field_source: "analytical"\n')
    solver_path.write_text(solver_text, encoding="utf-8")

    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    case_path = tmp_path / "case_restart_warn.yml"
    case_text = (valid / "case.yml").read_text(encoding="utf-8")
    case_text = case_text.replace("  start_step: 0\n", "  start_step: 5\n")
    case_text = case_text.replace("      count: 0\n", "      count: 10\n")
    case_path.write_text(case_text, encoding="utf-8")

    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    invalid = FIXTURES / "invalid"
    result = run_pic_flow(
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
    invalid = FIXTURES / "invalid"
    cluster_result = run_pic_flow(["validate", "--cluster", str(invalid / "cluster_bad_scheduler.yml")])
    assert cluster_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in cluster_result.stderr

    study_result = run_pic_flow(["validate", "--study", str(invalid / "study_bad_parameter_key.yml")])
    assert study_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in study_result.stderr


def test_dry_run_does_not_create_run_directories():
    valid = FIXTURES / "valid"
    runs_dir = REPO_ROOT / "runs"
    before = {p.name for p in runs_dir.iterdir()} if runs_dir.exists() else set()

    result = run_pic_flow(
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
    valid = FIXTURES / "valid"
    result = run_pic_flow(
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


def test_dry_run_local_solver_vs_post_proc_counts():
    valid = FIXTURES / "valid"
    result = run_pic_flow(
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


def test_dry_run_cluster_post_is_single_task(tmp_path):
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
    result = run_pic_flow(
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


def test_cluster_no_submit_manifest_and_scripts_use_stage_specific_counts(tmp_path):
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

    result = run_pic_flow(
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
    assert "srun -n 4 " in solver_script

    post_script = (run_dir / "scheduler" / "post.sbatch").read_text(encoding="utf-8")
    assert "#SBATCH --nodes=1" in post_script
    assert "#SBATCH --ntasks-per-node=1" in post_script
    assert "srun -n 1 " in post_script


def test_markdown_link_checker_passes():
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
