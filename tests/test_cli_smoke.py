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
