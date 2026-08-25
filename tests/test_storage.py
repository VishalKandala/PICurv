"""!
@file test_storage.py
@brief Focused storage inventory, transport, offload, restore, and CLI integration tests.
"""

import hashlib
import json
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from picurv_cli import core
from picurv_cli import storage
from picurv_cli.cli import build_main_parser, dispatch_command


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "valid"


def _write_run(root: Path, with_checkpoint: bool = True) -> Path:
    """!
    @brief Create a small representative standalone run artifact.
    @param[in] root Value supplied through the `root` argument.
    @param[in] with_checkpoint Value supplied through the `with_checkpoint` argument.
    @return Result produced by this operation.
    """
    (root / "config").mkdir(parents=True)
    (root / "scheduler").mkdir()
    (root / "logs").mkdir()
    (root / "config" / "case.yml").write_text(
        "models: {}\nboundary_conditions: []\nrun_control:\n  start_step: 0\n  total_steps: 10\n",
        encoding="utf-8",
    )
    (root / "config" / "solver.yml").write_text("operation_mode:\n  eulerian_field_source: solve\n", encoding="utf-8")
    (root / "config" / "monitor.yml").write_text(
        "io:\n  directories:\n    output: output\n    restart: restart\n    log: logs\nlogging:\n  verbosity: INFO\n",
        encoding="utf-8",
    )
    (root / "config" / "post.yml").write_text(
        "source_data:\n  directory: <solver_output_dir>\nio:\n  output_directory: viz\n  output_filename_prefix: Field\n",
        encoding="utf-8",
    )
    (root / "config" / f"{root.name}.control").write_text(
        f"-output_dir output\n-grid_file {root / 'config' / 'grid.run'}\n",
        encoding="utf-8",
    )
    (root / "config" / "grid.run").write_bytes(b"grid")
    (root / "logs" / "solver.log").write_text("finished\n", encoding="utf-8")
    (root / "manifest.json").write_text(json.dumps({"run_id": root.name}) + "\n", encoding="utf-8")
    (root / "viz").mkdir()
    (root / "viz" / "Field_00010.vts").write_bytes(b"derived-output")
    if with_checkpoint:
        bundle = root / "output" / "checkpoints" / "step_000000000010"
        (bundle / "eulerian" / "block_0000").mkdir(parents=True)
        (bundle / "eulerian" / "block_0000" / "Ucat.dat").write_bytes(b"checkpoint-data")
        (bundle / "checkpoint.meta").write_text("test metadata\n", encoding="utf-8")
        (bundle / "COMMITTED").write_text("0" * 64 + "\n", encoding="ascii")
    return root


def _write_study(root: Path) -> Path:
    """!
    @brief Create a small study with two materialized members.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    (root / "cases").mkdir(parents=True)
    (root / "scheduler").mkdir()
    (root / "results").mkdir()
    (root / "study.yml").write_text(
        "base_configs:\n  case: base_configs/case.yml\n  solver: base_configs/solver.yml\n"
        "  monitor: base_configs/monitor.yml\n  post: base_configs/post.yml\n"
        "study_type: sensitivity\nparameters:\n  case.run_control.total_steps: [10, 20]\n",
        encoding="utf-8",
    )
    (root / "cluster.yml").write_text("scheduler:\n  type: slurm\n", encoding="utf-8")
    (root / "study_manifest.json").write_text(json.dumps({"study_id": root.name}) + "\n", encoding="utf-8")
    index_lines = []
    for index in range(2):
        case = _write_run(root / "cases" / f"case_{index:04d}")
        index_lines.append(
            f"{index}\t{case.name}\t{case}\t{case / 'config' / (case.name + '.control')}\t"
            f"{case / 'config' / 'post.run'}\tINFO\tField\t\t\n"
        )
    (root / "scheduler" / "case_index.tsv").write_text("".join(index_lines), encoding="utf-8")
    (root / "results" / "metrics_table.csv").write_text("case_id,value\ncase_0000,1\n", encoding="utf-8")
    return root


@pytest.fixture
def local_rclone(monkeypatch, tmp_path):
    """!
    @brief Replace rclone process calls with a deterministic local object store.
    @param[in] monkeypatch Value supplied through the `monkeypatch` argument.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @return Result produced by this operation.
    """
    remote_root = tmp_path / "remote"

    def resolve(value: str) -> Path:
        """!
        @brief Map an rclone-style test path into the fixture object store.
        @param[in] value Value supplied through the `value` argument.
        @return Result produced by this operation.
        """
        if str(value).startswith("fake:"):
            return remote_root / str(value)[len("fake:"):].lstrip("/")
        return Path(value)

    def completed(args, stdout="", returncode=0, stderr=""):
        """!
        @brief Build a subprocess-like result for the fake rclone boundary.
        @param[in] args Value supplied through the `args` argument.
        @param[in] stdout Value supplied through the `stdout` argument.
        @param[in] returncode Value supplied through the `returncode` argument.
        @param[in] stderr Value supplied through the `stderr` argument.
        @return Result produced by this operation.
        """
        return SimpleNamespace(args=args, stdout=stdout, stderr=stderr, returncode=returncode)

    def fake_run(arguments, check=True):
        """!
        @brief Emulate the rclone subcommands used by storage workflow tests.
        @param[in] arguments Value supplied through the `arguments` argument.
        @param[in] check Value supplied through the `check` argument.
        @return Result produced by this operation.
        """
        command = arguments[0]
        if command == "mkdir":
            resolve(arguments[1]).mkdir(parents=True, exist_ok=True)
            return completed(arguments)
        if command == "copyto":
            source = resolve(arguments[1])
            destination = resolve(arguments[2])
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            return completed(arguments)
        if command == "hashsum":
            path = resolve(arguments[2])
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            return completed(arguments, f"{digest}  {path.name}\n")
        if command == "lsf":
            root = resolve(arguments[1])
            paths = []
            if root.is_dir():
                paths = [
                    path.relative_to(root).as_posix()
                    for path in root.rglob(storage.REMOTE_MANIFEST_FILENAME)
                ]
            return completed(arguments, "".join(path + "\n" for path in sorted(paths)))
        raise AssertionError(f"Unexpected fake rclone command: {arguments}")

    def fake_read(remote_path):
        """!
        @brief Read one remote object from the fixture object store.
        @param[in] remote_path Value supplied through the `remote_path` argument.
        @return Result produced by this operation.
        """
        return resolve(remote_path).read_bytes()

    monkeypatch.setattr(storage, "_run_rclone", fake_run)
    monkeypatch.setattr(storage, "_read_remote_bytes", fake_read)
    return remote_root


def _profile(tmp_path: Path) -> dict:
    """!
    @brief Return an in-memory storage profile for direct workflow tests.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @return Result produced by this operation.
    """
    return {
        "name": "archive",
        "remote": "fake:picurv-data",
        "config_path": str(tmp_path / storage.STORAGE_CONFIG_FILENAME),
        "compression": "none",
        "chunk_size_bytes": 1024 * 1024,
    }


@pytest.mark.skipif(shutil.which("rclone") is None, reason="rclone is not installed")
def test_real_rclone_local_backend_round_trip(tmp_path):
    """!
    @brief Exercise the production rclone subprocess boundary against a local filesystem backend.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    """
    run = _write_run(tmp_path / "runs" / "real-rclone")
    remote = tmp_path / "rclone-remote"
    profile = _profile(tmp_path)
    profile["remote"] = str(remote)
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, profile, label="real transport", prune_local=True)
    shutil.rmtree(run)
    storage.restore_archive(profile, manifest["archive_id"])
    assert (run / "viz" / "Field_00010.vts").read_bytes() == b"derived-output"
    assert storage.verify_remote_archive(profile, manifest["archive_id"])["label"] == "real transport"
    assert [item["archive_id"] for item in storage.list_remote_manifests(profile)] == [manifest["archive_id"]]


def test_storage_parser_uses_existing_run_study_selector_style():
    """! @brief Storage actions parse through the main argparse command surface. """
    parser = build_main_parser()
    args = parser.parse_args([
        "storage", "offload", "--study-dir", "studies/grid", "--case-id", "case_0003", "--dry-run"
    ])
    assert args.command == "storage"
    assert args.storage_action == "offload"
    assert args.study_dir == "studies/grid"
    assert args.case_ids == ["case_0003"]


def test_storage_setup_writes_only_non_secret_profile_data(tmp_path, local_rclone):
    """!
    @brief Setup uses rclone access while persisting only remote policy, never credentials.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    config = tmp_path / storage.STORAGE_CONFIG_FILENAME
    storage.storage_setup_workflow(SimpleNamespace(
        storage_config=str(config),
        profile="archive",
        remote="fake:picurv-data",
        compression="auto",
        chunk_size_gib=2.0,
        staging_directory=None,
        dry_run=False,
    ))
    payload = yaml.safe_load(config.read_text(encoding="utf-8"))
    assert payload == {
        "default_profile": "archive",
        "profiles": {
            "archive": {
                "remote": "fake:picurv-data",
                "compression": "auto",
                "chunk_size_gib": 2.0,
            }
        },
    }
    assert (local_rclone / "picurv-data" / "objects").is_dir()


def test_inventory_recognizes_committed_checkpoint_and_components(tmp_path):
    """!
    @brief Inventory identifies canonical checkpoints without assuming every output is heavy metadata.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    """
    run = _write_run(tmp_path / "runs" / "channel_20260824-120000")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    inventory = storage.inspect_artifact(target, query_scheduler=False)
    assert inventory["checkpoint_steps"] == [10]
    assert inventory["incomplete_checkpoints"] == []
    assert inventory["total_bytes"] > 0
    assert {entry["component"] for entry in inventory["entries"]} >= {"metadata", "logs", "data", "checkpoint:10"}


def test_incomplete_checkpoint_refuses_archival(tmp_path, local_rclone):
    """!
    @brief A checkpoint still being written is never packaged or pruned.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "active")
    incomplete = run / "output" / "checkpoints" / ".step_000000000020.incomplete.123"
    incomplete.mkdir(parents=True)
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    with pytest.raises(storage.StorageError, match="incomplete checkpoint"):
        storage.archive_artifact(target, _profile(tmp_path))


def test_protect_writes_remote_complete_catalog_and_retains_local_data(tmp_path, local_rclone):
    """!
    @brief Protect is a verified copy operation and never prunes the source.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "protected")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path), label="baseline")
    state = storage.read_storage_state(str(run))
    assert state["archive_id"] == manifest["archive_id"]
    assert state["local_pruned"] is False
    assert (run / "output" / "checkpoints" / "step_000000000010" / "COMMITTED").is_file()
    assert (local_rclone / "picurv-data" / "objects" / manifest["archive_id"] / "COMPLETE").is_file()
    assert storage.verify_remote_archive(_profile(tmp_path), manifest["archive_id"])["label"] == "baseline"


def test_offload_then_remote_only_restore_recovers_deleted_run(tmp_path, local_rclone):
    """!
    @brief A remote archive ID is sufficient after the entire local tree and marker are deleted.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "cold")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path), label="cold run", prune_local=True)
    assert storage.is_artifact_cold(str(run))
    assert not (run / "viz" / "Field_00010.vts").exists()
    assert (run / "config" / "case.yml").is_file()

    shutil.rmtree(run)
    restored = tmp_path / "recovered" / "cold"
    storage.restore_archive(_profile(tmp_path), manifest["archive_id"], destination=str(restored))
    assert (restored / "viz" / "Field_00010.vts").read_bytes() == b"derived-output"
    assert (restored / "output" / "checkpoints" / "step_000000000010" / "eulerian" / "block_0000" / "Ucat.dat").read_bytes() == b"checkpoint-data"
    assert storage.storage_state_summary(str(restored))["state"] == "PROTECTED"


def test_selective_checkpoint_restore_leaves_partial_marker(tmp_path, local_rclone):
    """!
    @brief Checkpoint selection restores raw state without downloading derived output chunks.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "selective")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path), prune_local=True)
    storage.restore_archive(_profile(tmp_path), manifest["archive_id"], checkpoints=[10])
    assert (run / "output" / "checkpoints" / "step_000000000010" / "COMMITTED").is_file()
    assert not (run / "viz" / "Field_00010.vts").exists()
    assert storage.storage_state_summary(str(run))["state"] == "PARTIAL"
    storage.require_storage_payload_local(str(run), "post-processing", checkpoints=[10])
    with pytest.raises(storage.StorageError, match=r"--checkpoint 20"):
        storage.require_storage_payload_local(str(run), "post-processing", checkpoints=[10, 20])


def test_study_member_archive_restores_parent_context(tmp_path, local_rclone):
    """!
    @brief Individually archived study members carry enough control-plane context for remote-only discovery.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    study = _write_study(tmp_path / "studies" / "grid_20260824-140000")
    target = storage.resolve_local_storage_targets(None, str(study), ["case_0001"])[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path), label="64-grid", prune_local=True)
    shutil.rmtree(study)

    destination = tmp_path / "restored-study" / "cases" / "case_0001"
    storage.restore_archive(_profile(tmp_path), manifest["archive_id"], destination=str(destination))
    assert (destination / "config" / "case.yml").is_file()
    assert (destination.parent.parent / "study.yml").is_file()
    assert (destination.parent.parent / "scheduler" / "case_index.tsv").is_file()


def test_whole_study_offload_retains_each_member_control_plane(tmp_path, local_rclone):
    """!
    @brief Whole-study pruning keeps member configs/logs/results and exposes inherited cold state.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    study = _write_study(tmp_path / "studies" / "whole_20260824-150000")
    target = storage.resolve_local_storage_targets(None, str(study))[0]
    storage.archive_artifact(target, _profile(tmp_path), label="whole campaign", prune_local=True)
    for case_id in ("case_0000", "case_0001"):
        case = study / "cases" / case_id
        assert (case / "config" / "case.yml").is_file()
        assert (case / "logs" / "solver.log").is_file()
        assert (case / "viz" / "Field_00010.vts").exists() is False
        assert storage.is_artifact_cold(str(case))
    assert storage.cold_study_members(str(study)) == ["case_0000", "case_0001"]


def test_remote_catalog_search_survives_missing_local_tree(tmp_path, local_rclone):
    """!
    @brief Remote manifests remain enumerable after every local artifact is removed.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "catalogued")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path), label="64-grid pilot", tags=["campaign=grid"])
    shutil.rmtree(run)
    listed = storage.list_remote_manifests(_profile(tmp_path))
    assert [item["archive_id"] for item in listed] == [manifest["archive_id"]]
    assert listed[0]["tags"] == {"campaign": "grid"}


def test_cold_run_guard_provides_exact_restore_command(tmp_path):
    """!
    @brief Existing workflows receive an actionable refusal instead of treating cold data as absent.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    """
    run = _write_run(tmp_path / "runs" / "guarded")
    marker = {
        "storage_schema_version": 1,
        "archive_id": "a" * 32,
        "local_pruned": True,
    }
    (run / storage.STORAGE_STATE_FILENAME).write_text(json.dumps(marker) + "\n", encoding="utf-8")
    with pytest.raises(storage.StorageError, match=r"picurv storage restore --archive-id"):
        storage.require_storage_payload_local(str(run), "continuation", checkpoint=10)

    case_cfg = {"run_control": {"start_step": 10, "total_steps": 1}, "models": {"physics": {"particles": {"restart_mode": "init"}}}}
    solver_cfg = {"operation_mode": {"eulerian_field_source": "solve"}}
    args = SimpleNamespace(restart_from=str(run), continue_run=False, run_dir=None)
    with pytest.raises(ValueError, match=r"storage restore --archive-id"):
        core.resolve_restart_source(args, case_cfg, solver_cfg, {}, str(tmp_path / "new-run"))


def test_post_and_submit_refuse_cold_run_before_missing_data_is_misdiagnosed(tmp_path, capsys):
    """!
    @brief Existing CLI paths stop at the storage boundary with an exact restore instruction.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] capsys Value supplied through the `capsys` argument.
    """
    run = _write_run(tmp_path / "runs" / "cold-cli")
    marker = {"storage_schema_version": 1, "archive_id": "b" * 32, "local_pruned": True}
    (run / storage.STORAGE_STATE_FILENAME).write_text(json.dumps(marker) + "\n", encoding="utf-8")

    post_args = build_main_parser().parse_args([
        "run", "--post-process", "--run-dir", str(run), "--post", str(run / "config" / "post.yml")
    ])
    with pytest.raises(SystemExit) as post_exit:
        dispatch_command(post_args)
    assert post_exit.value.code == 1
    assert "storage restore --archive-id" in capsys.readouterr().err

    (run / "scheduler" / "submission.json").write_text(
        json.dumps({"launch_mode": "local", "stages": {}}) + "\n", encoding="utf-8"
    )
    submit_args = SimpleNamespace(run_dir=str(run), study_dir=None, stage="all", force=False, dry_run=False)
    with pytest.raises(SystemExit) as submit_exit:
        core.submit_staged_jobs(submit_args)
    assert submit_exit.value.code == 1
    assert "storage restore --archive-id" in capsys.readouterr().err


def test_sweep_continue_and_reaggregate_refuse_cold_members(tmp_path, monkeypatch, capsys):
    """!
    @brief Sweep maintenance never reclassifies a deliberately offloaded member as empty.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] monkeypatch Value supplied through the `monkeypatch` argument.
    @param[in] capsys Value supplied through the `capsys` argument.
    """
    monkeypatch.chdir(tmp_path)
    core.sweep_workflow(SimpleNamespace(
        study=str(FIXTURES / "study.yml"),
        cluster=str(FIXTURES / "cluster.yml"),
        no_submit=True,
    ))
    study = next((tmp_path / "studies").iterdir())
    case = next((study / "cases").iterdir())
    marker = {"storage_schema_version": 1, "archive_id": "c" * 32, "local_pruned": True}
    (case / storage.STORAGE_STATE_FILENAME).write_text(json.dumps(marker) + "\n", encoding="utf-8")

    with pytest.raises(SystemExit) as continue_exit:
        core.sweep_continue_workflow(SimpleNamespace(study_dir=str(study), cluster=None, no_submit=True))
    assert continue_exit.value.code == 1
    assert "cold-storage member" in capsys.readouterr().err

    with pytest.raises(SystemExit) as aggregate_exit:
        core.sweep_reaggregate_workflow(SimpleNamespace(study_dir=str(study)))
    assert aggregate_exit.value.code == 1
    assert "could replace values with blanks" in capsys.readouterr().err


def test_runtime_solver_lock_is_visible_and_cleaned(tmp_path):
    """!
    @brief Python-owned solver markers protect local execution without a C runtime change.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    """
    run = _write_run(tmp_path / "runs" / "locked", with_checkpoint=False)
    lock = run / "scheduler" / "solver.lock.json"
    with storage.runtime_stage_lock(str(run), "solver"):
        assert lock.is_file()
        target = storage.resolve_local_storage_targets(str(run), None)[0]
        inventory = storage.inspect_artifact(target, query_scheduler=False)
        assert inventory["active_locks"] == ["scheduler/solver.lock.json"]
    assert not lock.exists()


def test_relocated_restore_rebases_known_generated_paths(tmp_path, local_rclone):
    """!
    @brief Explicit relocation updates generated text artifacts while preserving archived bytes remotely.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "original")
    control = run / "config" / f"{run.name}.control"
    control.write_text(f"-restart_dir {run / 'output'}\n", encoding="utf-8")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path))
    shutil.rmtree(run)
    destination = tmp_path / "runs" / "relocated"
    storage.restore_archive(_profile(tmp_path), manifest["archive_id"], destination=str(destination))
    restored_control = destination / "config" / "original.control"
    assert str(destination / "output") in restored_control.read_text(encoding="utf-8")
    assert storage.read_storage_state(str(destination))["relocated_from"] == str(run)
