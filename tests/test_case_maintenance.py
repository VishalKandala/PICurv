import importlib.machinery
import importlib.util
import json
import os
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "scripts" / "picurv"


def load_picurv_module():
    loader = importlib.machinery.SourceFileLoader("picurv_case_maintenance_module", str(PICURV))
    spec = importlib.util.spec_from_loader("picurv_case_maintenance_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@contextmanager
def pushd(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def make_fake_source_repo(root: Path) -> Path:
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


def test_build_project_uses_case_origin_source_root(tmp_path):
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    picurv.write_case_origin_metadata(str(case_dir), str(source_root), template_name="demo")

    captured = {}

    def fake_execute(command, run_dir, log_filename, monitor_cfg=None):
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


def test_sync_case_binaries_command_refreshes_case_and_bootstraps_metadata(tmp_path):
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


def test_sync_config_overwrite_replaces_modified_files_when_requested(tmp_path):
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


def test_pull_source_repo_uses_git_pull_rebase_by_default(tmp_path):
    picurv = load_picurv_module()
    source_root = make_fake_source_repo(tmp_path / "source")
    case_dir = tmp_path / "case"
    case_dir.mkdir()
    picurv.write_case_origin_metadata(str(case_dir), str(source_root), template_name="demo")

    captured = {}

    def fake_execute(command, run_dir, log_filename, monitor_cfg=None):
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
