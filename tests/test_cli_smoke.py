import importlib.machinery
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "scripts" / "picurv"
FIXTURES = REPO_ROOT / "tests" / "fixtures"


def run_picurv(args, cwd=REPO_ROOT):
    cmd = [sys.executable, str(PICURV)] + list(args)
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, timeout=60, check=False)


def load_picurv_module():
    loader = importlib.machinery.SourceFileLoader("picurv_module", str(PICURV))
    spec = importlib.util.spec_from_loader("picurv_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_top_level_help_smoke():
    result = run_picurv(["--help"])
    assert result.returncode == 0
    assert "validate" in result.stdout
    assert "Next commands:" in result.stdout


def test_run_help_smoke():
    result = run_picurv(["run", "--help"])
    assert result.returncode == 0
    assert "--dry-run" in result.stdout
    assert "--format {text,json}" in result.stdout


def test_validate_help_smoke():
    result = run_picurv(["validate", "--help"])
    assert result.returncode == 0
    assert "--case" in result.stdout
    assert "--strict" in result.stdout


def test_removed_selector_aliases_are_rejected():
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
    invalid = FIXTURES / "invalid"
    cluster_result = run_picurv(["validate", "--cluster", str(invalid / "cluster_bad_scheduler.yml")])
    assert cluster_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in cluster_result.stderr

    study_result = run_picurv(["validate", "--study", str(invalid / "study_bad_parameter_key.yml")])
    assert study_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in study_result.stderr


def test_dry_run_does_not_create_run_directories():
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


def test_programmatic_grid_cell_counts_translate_to_node_counts(tmp_path):
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


def test_grid_gen_exports_node_counts_from_cell_inputs(tmp_path):
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


def test_case_local_symlinked_picurv_prefers_local_binaries(tmp_path):
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


def test_init_always_copies_full_executable_set(tmp_path):
    picurv = load_picurv_module()
    fake_root = tmp_path / "fake_root"
    template_dir = fake_root / "examples" / "demo_case"
    bin_dir = fake_root / "bin"
    template_dir.mkdir(parents=True)
    bin_dir.mkdir(parents=True)

    (template_dir / "case.yml").write_text("run_control:\n  start_step: 0\n", encoding="utf-8")
    for exe_name in ("picurv", "simulator", "postprocessor"):
        exe_path = bin_dir / exe_name
        exe_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        exe_path.chmod(0o755)

    original_project_root = picurv.PROJECT_ROOT
    original_cwd = Path.cwd()
    work_dir = tmp_path / "work"
    work_dir.mkdir()

    try:
        picurv.PROJECT_ROOT = str(fake_root)
        os.chdir(work_dir)
        picurv.init_case(SimpleNamespace(template_name="demo_case", dest_name="demo_out"))
    finally:
        os.chdir(original_cwd)
        picurv.PROJECT_ROOT = original_project_root

    out_dir = work_dir / "demo_out"
    for exe_name in ("picurv", "simulator", "postprocessor"):
        exe_path = out_dir / exe_name
        assert exe_path.is_file()
        assert not exe_path.is_symlink()


def test_empty_enabled_functions_omits_whitelist_and_uses_c_default(tmp_path):
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


def test_particle_console_output_frequency_defaults_to_data_output_frequency(tmp_path):
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


def test_new_profiling_config_emits_explicit_timestep_flags(tmp_path):
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


def test_restart_from_run_dir_resolves_previous_restart_directory(tmp_path):
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()
    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    monitor_cfg = picurv.read_yaml_file(str(valid / "monitor.yml"))

    previous_run_dir = tmp_path / "old_run"
    (previous_run_dir / "config").mkdir(parents=True)
    (previous_run_dir / "prior_results").mkdir()

    previous_monitor_cfg = {
        "io": {
            "directories": {
                "restart": "prior_results",
            }
        }
    }
    picurv.write_yaml_file(str(previous_run_dir / "config" / "monitor.yml"), previous_monitor_cfg)

    case_cfg["run_control"]["start_step"] = 5
    case_cfg["run_control"]["restart_from_run_dir"] = "old_run"
    case_path = tmp_path / "case_restart.yml"
    picurv.write_yaml_file(str(case_path), case_cfg)

    resolved = picurv.resolve_restart_source_from_run_dir(
        case_cfg,
        solver_cfg,
        monitor_cfg,
        str(case_path),
        strict=True,
    )
    assert resolved == str(previous_run_dir / "prior_results")

    run_dir = tmp_path / "new_run"
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
        restart_source_dir=resolved,
    )

    content = Path(control_file).read_text(encoding="utf-8")
    assert f"-restart_dir {resolved}" in content


def test_dry_run_local_solver_vs_post_proc_counts():
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
