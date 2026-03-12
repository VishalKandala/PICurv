"""!
@file test_case_maintenance.py
@brief Pytest coverage for case initialization, sync, build, and source-maintenance workflows.
"""

import importlib.machinery
import importlib.util
import json
import os
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "scripts" / "picurv"


def load_picurv_module():
    """!
    @brief Load `scripts/picurv` as an importable module for maintenance tests.
    @return Value returned by `load_picurv_module()`.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_case_maintenance_module", str(PICURV))
    spec = importlib.util.spec_from_loader("picurv_case_maintenance_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@contextmanager
def pushd(path: Path):
    """!
    @brief Temporarily change working directory within a context manager.
    @param[in] path Filesystem path argument passed to `pushd()`.
    """
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def make_fake_source_repo(root: Path) -> Path:
    """!
    @brief Create a minimal fake PICurv source tree for case-maintenance tests.
    @param[in] root Argument passed to `make_fake_source_repo()`.
    @return Value returned by `make_fake_source_repo()`.
    """
    (root / "src").mkdir(parents=True)
    (root / "include").mkdir()
    (root / "scripts").mkdir()
    (root / "examples" / "demo" / "nested").mkdir(parents=True)
    (root / "bin").mkdir()

    (root / "Makefile").write_text("all:\n\t@echo ok\n", encoding="utf-8")
    (root / "examples" / "demo" / "case.yml").write_text("case: baseline\n", encoding="utf-8")
    (root / "examples" / "demo" / "guide.md").write_text("# demo\n", encoding="utf-8")
    (root / "examples" / "demo" / "nested" / "notes.txt").write_text("nested template file\n", encoding="utf-8")
    (root / "bin" / "simulator").write_text("simulator-v1\n", encoding="utf-8")
    (root / "bin" / "postprocessor").write_text("postprocessor-v1\n", encoding="utf-8")
    (root / "bin" / "picurv").write_text("picurv-v1\n", encoding="utf-8")
    return root


def test_init_case_writes_origin_metadata_and_copies_binaries(tmp_path):
    """!
    @brief Test that init case writes origin metadata and copies binaries.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    args = SimpleNamespace(template_name="demo", dest_name="my_case", source_root=str(source_root))
    with pushd(workspace):
        picurv.init_case(args)

    case_dir = workspace / "my_case"
    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))

    assert metadata["source_repo_root"] == str(source_root.resolve())
    assert metadata["template_name"] == "demo"
    assert metadata["template_managed_files"] == ["case.yml", "guide.md", "nested/notes.txt"]
    assert (case_dir / "case.yml").read_text(encoding="utf-8") == "case: baseline\n"
    assert (case_dir / "simulator").read_text(encoding="utf-8") == "simulator-v1\n"
    assert (case_dir / "postprocessor").read_text(encoding="utf-8") == "postprocessor-v1\n"
    assert (case_dir / "picurv").read_text(encoding="utf-8") == "picurv-v1\n"
    runtime_cfg = yaml.safe_load((case_dir / picurv.RUNTIME_EXECUTION_CONFIG_FILENAME).read_text(encoding="utf-8"))
    assert runtime_cfg == {
        "default_execution": {},
        "local_execution": {},
        "cluster_execution": {},
    }


def test_init_case_seeds_runtime_from_repo_local_config_and_omits_execution_example(tmp_path):
    """!
    @brief Test that init prefers repo-local runtime config and does not leave execution.example.yml in the case.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    (source_root / picurv.RUNTIME_EXECUTION_CONFIG_FILENAME).write_text(
        "default_execution:\n  launcher: mpirun\n  launcher_args: [--bind-to, none]\nlocal_execution: {}\ncluster_execution: {}\n",
        encoding="utf-8",
    )
    (source_root / "examples" / "demo" / picurv.RUNTIME_EXECUTION_EXAMPLE_FILENAME).write_text(
        "default_execution:\n  launcher: should_not_copy\n",
        encoding="utf-8",
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    args = SimpleNamespace(template_name="demo", dest_name="my_case", source_root=str(source_root))
    with pushd(workspace):
        picurv.init_case(args)

    case_dir = workspace / "my_case"
    runtime_cfg = yaml.safe_load((case_dir / picurv.RUNTIME_EXECUTION_CONFIG_FILENAME).read_text(encoding="utf-8"))
    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))

    assert runtime_cfg["default_execution"]["launcher"] == "mpirun"
    assert not (case_dir / picurv.RUNTIME_EXECUTION_EXAMPLE_FILENAME).exists()
    assert picurv.RUNTIME_EXECUTION_EXAMPLE_FILENAME not in metadata["template_managed_files"]


def test_init_case_prints_note_for_cluster_profile_samples(tmp_path, capsys):
    """!
    @brief Test that init warns users to edit copied cluster profile samples.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    @param[in] capsys Pytest capture fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    (source_root / "examples" / "demo" / "slurm_cluster.yml").write_text(
        "scheduler:\n  type: slurm\n",
        encoding="utf-8",
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    args = SimpleNamespace(template_name="demo", dest_name="my_case", source_root=str(source_root))
    with pushd(workspace):
        picurv.init_case(args)

    captured = capsys.readouterr()
    assert "Cluster profile sample(s) copied with this case" in captured.out
    assert "slurm_cluster.yml" in captured.out
    assert "before using --cluster" in captured.out


def test_build_project_uses_case_origin_source_root(tmp_path):
    """!
    @brief Test that build project uses case origin source root.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    picurv.write_case_origin_metadata(str(case_dir), str(source_root), template_name="demo")

    captured = {}

    def fake_execute(command, run_dir, log_filename, monitor_cfg=None):
        """!
        @brief Record build-command execution requests without launching a subprocess.
        @param[in] command Argument passed to `fake_execute()`.
        @param[in] run_dir Argument passed to `fake_execute()`.
        @param[in] log_filename Argument passed to `fake_execute()`.
        @param[in] monitor_cfg Argument passed to `fake_execute()`.
        """
        captured["command"] = command
        captured["run_dir"] = run_dir
        captured["log_filename"] = log_filename

    original_execute = picurv.execute_command
    picurv.execute_command = fake_execute
    try:
        picurv.build_project(
            SimpleNamespace(case_dir=str(case_dir), source_root=None, make_args=["clean-project"])
        )
    finally:
        picurv.execute_command = original_execute

    assert captured["command"] == ["make", "clean-project"]
    assert captured["run_dir"] == str(source_root.resolve())
    assert captured["log_filename"] == "build.log"


def test_build_project_defaults_to_all_when_make_args_have_no_goal(tmp_path):
    """!
    @brief Test that build project injects `all` when only make assignments/options are passed.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    picurv.write_case_origin_metadata(str(case_dir), str(source_root), template_name="demo")

    captured = {}

    def fake_execute(command, run_dir, log_filename, monitor_cfg=None):
        """!
        @brief Record execute-command inputs while stubbing out subprocess execution.
        @param[in] command Argument passed to `fake_execute()`.
        @param[in] run_dir Argument passed to `fake_execute()`.
        @param[in] log_filename Argument passed to `fake_execute()`.
        @param[in] monitor_cfg Argument passed to `fake_execute()`.
        """
        captured["command"] = command
        captured["run_dir"] = run_dir
        captured["log_filename"] = log_filename

    original_execute = picurv.execute_command
    picurv.execute_command = fake_execute
    try:
        picurv.build_project(
            SimpleNamespace(case_dir=str(case_dir), source_root=None, make_args=["SYSTEM=cluster"])
        )
    finally:
        picurv.execute_command = original_execute

    assert captured["command"] == ["make", "all", "SYSTEM=cluster"]
    assert captured["run_dir"] == str(source_root.resolve())
    assert captured["log_filename"] == "build.log"


def test_sync_case_binaries_command_refreshes_case_and_bootstraps_metadata(tmp_path):
    """!
    @brief Test that sync case binaries command refreshes case and bootstraps metadata.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    (case_dir / "simulator").write_text("stale\n", encoding="utf-8")

    picurv.sync_case_binaries_command(
        SimpleNamespace(case_dir=str(case_dir), source_root=str(source_root))
    )

    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))
    assert metadata["source_repo_root"] == str(source_root.resolve())
    assert (case_dir / "simulator").read_text(encoding="utf-8") == "simulator-v1\n"
    assert (case_dir / "postprocessor").read_text(encoding="utf-8") == "postprocessor-v1\n"
    assert (case_dir / "picurv").read_text(encoding="utf-8") == "picurv-v1\n"


def test_sync_config_preserves_modified_files_without_overwrite(tmp_path):
    """!
    @brief Test that sync config preserves modified files without overwrite.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    (case_dir / "case.yml").write_text("case: customized\n", encoding="utf-8")

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=False,
            prune=False,
        )
    )

    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))
    assert metadata["template_name"] == "demo"
    assert metadata["template_managed_files"] == ["case.yml", "guide.md", "nested/notes.txt"]
    assert (case_dir / "case.yml").read_text(encoding="utf-8") == "case: customized\n"
    assert (case_dir / "guide.md").read_text(encoding="utf-8") == "# demo\n"
    assert (case_dir / "nested" / "notes.txt").read_text(encoding="utf-8") == "nested template file\n"


def test_sync_config_omits_execution_example_and_bootstraps_runtime_config(tmp_path):
    """!
    @brief Test that sync-config creates .picurv-execution.yml but does not copy execution.example.yml.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    (source_root / picurv.RUNTIME_EXECUTION_CONFIG_FILENAME).write_text(
        "default_execution:\n  launcher: mpirun\n  launcher_args: [--bind-to, none]\nlocal_execution: {}\ncluster_execution: {}\n",
        encoding="utf-8",
    )
    (source_root / "examples" / "demo" / picurv.RUNTIME_EXECUTION_EXAMPLE_FILENAME).write_text(
        "default_execution:\n  launcher: should_not_copy\n",
        encoding="utf-8",
    )
    case_dir = tmp_path / "case"
    case_dir.mkdir()

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=False,
            prune=False,
        )
    )

    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))
    runtime_cfg = yaml.safe_load((case_dir / picurv.RUNTIME_EXECUTION_CONFIG_FILENAME).read_text(encoding="utf-8"))

    assert runtime_cfg["default_execution"]["launcher"] == "mpirun"
    assert not (case_dir / picurv.RUNTIME_EXECUTION_EXAMPLE_FILENAME).exists()
    assert picurv.RUNTIME_EXECUTION_EXAMPLE_FILENAME not in metadata["template_managed_files"]


def test_sync_config_overwrite_replaces_modified_files_when_requested(tmp_path):
    """!
    @brief Test that sync config overwrite replaces modified files when requested.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    (case_dir / "case.yml").write_text("case: customized\n", encoding="utf-8")

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=True,
            prune=False,
        )
    )

    assert (case_dir / "case.yml").read_text(encoding="utf-8") == "case: baseline\n"


def test_sync_config_prune_removes_only_tracked_stale_template_files(tmp_path):
    """!
    @brief Test that sync config prune removes only tracked stale template files.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=False,
            prune=False,
        )
    )
    (case_dir / "local_only.txt").write_text("keep me\n", encoding="utf-8")
    (source_root / "examples" / "demo" / "guide.md").unlink()

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=False,
            prune=True,
        )
    )

    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))
    assert metadata["template_managed_files"] == ["case.yml", "nested/notes.txt"]
    assert not (case_dir / "guide.md").exists()
    assert (case_dir / "local_only.txt").read_text(encoding="utf-8") == "keep me\n"


def test_compute_case_source_status_reports_binary_and_template_drift(tmp_path):
    """!
    @brief Test that compute case source status reports binary and template drift.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=False,
            prune=False,
        )
    )
    picurv.sync_case_binaries_command(
        SimpleNamespace(case_dir=str(case_dir), source_root=str(source_root))
    )

    (case_dir / "simulator").write_text("simulator-local\n", encoding="utf-8")
    (case_dir / "case.yml").write_text("case: customized\n", encoding="utf-8")
    (source_root / "examples" / "demo" / "guide.md").unlink()

    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))
    status = picurv.compute_case_source_status(
        str(case_dir),
        str(source_root),
        template_name="demo",
        metadata=metadata,
    )

    assert "simulator" in status["binaries"]["case_bin_different"]
    assert "case.yml" in status["config"]["case_modified_files"]
    assert "guide.md" in status["config"]["template_removed_since_last_sync"]


def test_compute_case_source_status_reports_runtime_config_seed_match(tmp_path):
    """!
    @brief Test that status-source reports case runtime config presence and repo-seed match.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    (source_root / picurv.RUNTIME_EXECUTION_CONFIG_FILENAME).write_text(
        "default_execution:\n  launcher: mpirun\n  launcher_args: [--bind-to, none]\nlocal_execution: {}\ncluster_execution: {}\n",
        encoding="utf-8",
    )
    case_dir = tmp_path / "case"
    case_dir.mkdir()

    picurv.sync_case_config_command(
        SimpleNamespace(
            case_dir=str(case_dir),
            source_root=str(source_root),
            template_name="demo",
            overwrite=False,
            prune=False,
        )
    )

    metadata = json.loads((case_dir / picurv.CASE_ORIGIN_METADATA_FILENAME).read_text(encoding="utf-8"))
    status = picurv.compute_case_source_status(
        str(case_dir),
        str(source_root),
        template_name="demo",
        metadata=metadata,
    )

    assert status["runtime_execution"]["case_config_present"] is True
    assert status["runtime_execution"]["repo_seed_present"] is True
    assert status["runtime_execution"]["case_matches_repo_seed"] is True


def test_pull_source_repo_uses_git_pull_rebase_by_default(tmp_path):
    """!
    @brief Test that pull source repo uses git pull rebase by default.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    picurv.write_case_origin_metadata(str(case_dir), str(source_root), template_name="demo")

    captured = {}

    def fake_execute(command, run_dir, log_filename, monitor_cfg=None):
        """!
        @brief Record pull-command execution requests without invoking Git.
        @param[in] command Argument passed to `fake_execute()`.
        @param[in] run_dir Argument passed to `fake_execute()`.
        @param[in] log_filename Argument passed to `fake_execute()`.
        @param[in] monitor_cfg Argument passed to `fake_execute()`.
        """
        captured["command"] = command
        captured["run_dir"] = run_dir
        captured["log_filename"] = log_filename

    original_execute = picurv.execute_command
    picurv.execute_command = fake_execute
    try:
        picurv.pull_source_repo(
            SimpleNamespace(
                case_dir=str(case_dir),
                source_root=None,
                remote=None,
                branch=None,
                no_rebase=False,
            )
        )
    finally:
        picurv.execute_command = original_execute

    assert captured["command"] == ["git", "pull", "--rebase"]
    assert captured["run_dir"] == str(source_root.resolve())
    assert captured["log_filename"] == "pull-source.log"
