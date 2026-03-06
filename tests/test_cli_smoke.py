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


def run_picurv(args, cwd=REPO_ROOT, env=None):
    """Run picurv."""
    cmd = [sys.executable, str(PICURV)] + list(args)
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    return subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, timeout=60, check=False, env=merged_env)


def load_picurv_module():
    """Load picurv module for tests."""
    loader = importlib.machinery.SourceFileLoader("picurv_module", str(PICURV))
    spec = importlib.util.spec_from_loader("picurv_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def write_legacy_1d_grid(path: Path) -> Path:
    """Write a minimal legacy 1D-axis grid payload for conversion tests."""
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


def test_top_level_help_smoke():
    """Test that top level help smoke."""
    result = run_picurv(["--help"])
    assert result.returncode == 0
    assert "validate" in result.stdout
    assert "Next commands:" in result.stdout


def test_run_help_smoke():
    """Test that run help smoke."""
    result = run_picurv(["run", "--help"])
    assert result.returncode == 0
    assert "--dry-run" in result.stdout
    assert "--format {text,json}" in result.stdout


def test_validate_help_smoke():
    """Test that validate help smoke."""
    result = run_picurv(["validate", "--help"])
    assert result.returncode == 0
    assert "--case" in result.stdout
    assert "--strict" in result.stdout


def test_removed_selector_aliases_are_rejected():
    """Test that removed selector aliases are rejected."""
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
    """Test that picgrid validation requires canonical header."""
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
    """Test that validate valid configs pass."""
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
    """Test that validate zero mode case allows omitting velocity components."""
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
    """Test that validate poiseuille peak velocity option passes with unique inlet axis."""
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
    """Test that validate initial condition mode must be explicit."""
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
    """Test that validate analytical mode requires programmatic grid."""
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
    """Test that validate particle restart mode omission warns on restart."""
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
    """Test that validate invalid case reports structured error."""
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
    """Test that validate invalid cluster and study fail."""
    invalid = FIXTURES / "invalid"
    cluster_result = run_picurv(["validate", "--cluster", str(invalid / "cluster_bad_scheduler.yml")])
    assert cluster_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in cluster_result.stderr

    study_result = run_picurv(["validate", "--study", str(invalid / "study_bad_parameter_key.yml")])
    assert study_result.returncode == 1
    assert "ERROR CFG_INVALID_VALUE" in study_result.stderr


def test_dry_run_does_not_create_run_directories():
    """Test that dry run does not create run directories."""
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
    """Test that dry run json output schema."""
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
    """Test that post process run dir accepts null source data mapping."""
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()

    run_dir = tmp_path / "existing_run"
    config_dir = run_dir / "config"
    logs_dir = run_dir / "logs"
    results_dir = run_dir / "results"
    config_dir.mkdir(parents=True)
    logs_dir.mkdir()
    results_dir.mkdir()
    (results_dir / "dummy.dat").write_text("0\n", encoding="utf-8")

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
        """Helper for fake resolve runtime executable."""
        return f"/tmp/fake/{name}"

    def fake_execute_command(command, run_dir_arg, log_filename, monitor_cfg=None):
        """Helper for fake execute command."""
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
    assert (config_dir / "post.run").is_file()
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stages_completed_or_submitted"] == ["post-process"]


def test_dry_run_restart_from_missing_run_dir_fails(tmp_path):
    """Test that dry run restart from missing run dir fails."""
    valid = FIXTURES / "valid"
    picurv = load_picurv_module()

    case_cfg = picurv.read_yaml_file(str(valid / "case.yml"))
    solver_cfg = picurv.read_yaml_file(str(valid / "solver.yml"))
    case_cfg["run_control"]["start_step"] = 3
    case_cfg["run_control"]["restart_from_run_dir"] = "does_not_exist"
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
            "--post-process",
            "--case",
            str(case_path),
            "--solver",
            str(solver_path),
            "--monitor",
            str(valid / "monitor.yml"),
            "--post",
            str(valid / "post.yml"),
            "--dry-run",
            "--format",
            "json",
        ],
        cwd=tmp_path,
    )

    assert result.returncode == 1
    assert "restart_from_run_dir does not exist" in result.stderr


def test_post_process_run_dir_missing_config_inputs_fails(tmp_path):
    """Test that post process run dir missing config inputs fails."""
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
    """Test that programmatic grid cell counts translate to node counts."""
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
    """Test that grid gen exports node counts from cell inputs."""
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
    """Test that grid.gen legacy1d converts headerless 1D-axis payload to canonical PICGRID."""
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


def test_case_local_symlinked_picurv_prefers_local_binaries(tmp_path):
    """Test that case local symlinked picurv prefers local binaries."""
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
    """Test that case local copied picurv prefers local binaries."""
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
    """Test that init always copies full executable set."""
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
        exe_path = out_dir / exe_name
        assert exe_path.is_file()
        assert not exe_path.is_symlink()


def test_empty_enabled_functions_omits_whitelist_and_uses_c_default(tmp_path):
    """Test that empty enabled functions omits whitelist and uses c default."""
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
    """Test that file-grid mode can convert legacy payloads via grid.gen before staging grid.run."""
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
    """Test that particle console output frequency defaults to data output frequency."""
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
    """Test that validate rejects negative particle console output frequency."""
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
    """Test that validate rejects unsupported grid.legacy_conversion.format values."""
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


def test_new_profiling_config_emits_explicit_timestep_flags(tmp_path):
    """Test that new profiling config emits explicit timestep flags."""
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
    """Test that restart from run dir resolves previous restart directory."""
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
    """Test that dry run local solver vs post proc counts."""
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
    """Test that local multi-rank dry-run honors launcher override env vars."""
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


def test_dry_run_local_reads_optional_local_runtime_config(tmp_path):
    """Test that local multi-rank dry-run can read .picurv-local.yml."""
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-local.yml").write_text(
        "\n".join(
            [
                "local_execution:",
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


def test_env_launcher_override_wins_over_optional_local_runtime_config(tmp_path):
    """Test that env launcher override takes precedence over .picurv-local.yml."""
    valid = FIXTURES / "valid"
    (tmp_path / ".picurv-local.yml").write_text(
        "\n".join(
            [
                "local_execution:",
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


def test_dry_run_cluster_post_is_single_task(tmp_path):
    """Test that dry run cluster post is single task."""
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
    """Test that cluster dry-run accepts launcher strings with inline site flags."""
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


def test_validate_cluster_rejects_unparseable_inline_launcher(tmp_path):
    """Test that cluster validation rejects malformed inline launcher quoting."""
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


def test_cluster_no_submit_manifest_and_scripts_use_stage_specific_counts(tmp_path):
    """Test that cluster no submit manifest and scripts use stage specific counts."""
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
    """Test that markdown link checker passes."""
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
