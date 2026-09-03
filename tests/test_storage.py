"""!
@file test_storage.py
@brief Focused storage inventory, transport, offload, restore, and CLI integration tests.
"""

import hashlib
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from picurv_cli import core
from picurv_cli import storage
from picurv_cli.cli import build_main_parser, dispatch_command


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "valid"


def test_core_direct_file_load_resolves_sibling_storage_module(tmp_path: Path) -> None:
    """!
    @brief Verify standalone core loading does not depend on the checkout being importable.
    @param[in] tmp_path Temporary directory outside the repository.
    @return None.
    """
    script = f"""
import importlib.machinery
import importlib.util

path = {str(REPO_ROOT / "picurv_cli" / "core.py")!r}
loader = importlib.machinery.SourceFileLoader("picurv_core_standalone_test", path)
spec = importlib.util.spec_from_loader(loader.name, loader)
module = importlib.util.module_from_spec(spec)
loader.exec_module(module)
assert module.StorageError.__name__ == "StorageError"
"""
    result = subprocess.run(
        [sys.executable, "-I", "-c", script],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


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
        "io:\n  data_output_frequency: 10\nlogging:\n  verbosity: INFO\n",
        encoding="utf-8",
    )
    (root / "config" / "post.yml").write_text(
        "source_data: {}\nio:\n  output_filename_prefix: Field\n",
        encoding="utf-8",
    )
    (root / "config" / f"{root.name}.control").write_text(
        f"-output_dir output\n-grid_file {root / 'config' / 'grid.run'}\n",
        encoding="utf-8",
    )
    (root / "config" / "grid.run").write_bytes(b"grid")
    (root / "logs" / "solver.log").write_text("finished\n", encoding="utf-8")
    (root / "manifest.json").write_text(json.dumps({"run_id": root.name}) + "\n", encoding="utf-8")
    visualization = root / "output" / "visualization" / "Field-0123456789ab"
    visualization.mkdir(parents=True)
    (visualization / "Field_00010.vts").write_bytes(b"derived-output")
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
            if not source.is_file():
                return completed(arguments, returncode=1, stderr="source not found")
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            return completed(arguments)
        if command == "hashsum":
            path = resolve(arguments[2])
            if not path.is_file():
                # Real rclone exits non-zero for a missing object rather than raising.
                return completed(arguments, returncode=1, stderr="directory not found")
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

    # The package routes every rclone call through this one boundary, so replacing it
    # here reaches all of them.
    monkeypatch.setattr(storage.transport, "_run_rclone", fake_run)
    monkeypatch.setattr(storage.transport, "_read_remote_bytes", fake_read)
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
    assert (run / "output" / "visualization" / "Field-0123456789ab" / "Field_00010.vts").read_bytes() == b"derived-output"
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
                "workers": storage.DEFAULT_STORAGE_WORKERS,
                "offload_policy": "metadata-only",
                "keep_latest_checkpoint": False,
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
    assert {entry["component"] for entry in inventory["entries"]} >= {
        "metadata", "logs", "raw-output", "checkpoint:10"
    }


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
    assert not (run / "output" / "visualization" / "Field-0123456789ab" / "Field_00010.vts").exists()
    assert (run / "config" / "case.yml").is_file()

    shutil.rmtree(run)
    restored = tmp_path / "recovered" / "cold"
    storage.restore_archive(_profile(tmp_path), manifest["archive_id"], destination=str(restored))
    assert (restored / "output" / "visualization" / "Field-0123456789ab" / "Field_00010.vts").read_bytes() == b"derived-output"
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
    assert not (run / "output" / "visualization" / "Field-0123456789ab" / "Field_00010.vts").exists()
    assert storage.storage_state_summary(str(run))["state"] == "PARTIAL"
    storage.require_storage_payload_local(str(run), "post-processing", checkpoints=[10])
    with pytest.raises(storage.StorageError, match=r"--checkpoint 20"):
        storage.require_storage_payload_local(str(run), "post-processing", checkpoints=[10, 20])
    # Only the requested step may appear. Directory entries used to be swept into the
    # metadata chunk regardless of component, so restoring one checkpoint recreated the
    # empty tree of every other one and a partial run listed steps it did not hold.
    present = sorted(p.name for p in (run / "output" / "checkpoints").iterdir())
    assert present == ["step_000000000010"], present


def test_an_unchanged_artifact_is_not_uploaded_twice(tmp_path, local_rclone):
    """!
    @brief Offloading after protecting reuses the verified archive instead of re-uploading.

    @details `protect` now, `offload` later is the documented way to keep a backup while
             still working locally. Uploading a second full copy of an unchanged artifact
             doubles remote cost for the largest runs and buys nothing.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Value supplied through the `local_rclone` argument.
    """
    run = _write_run(tmp_path / "runs" / "reused")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    protected = storage.archive_artifact(
        target, _profile(tmp_path), label="kept", prune_local=False
    )
    offloaded = storage.archive_artifact(target, _profile(tmp_path), prune_local=True)

    assert offloaded["archive_id"] == protected["archive_id"]
    objects = sorted((local_rclone / "picurv-data" / "objects").iterdir())
    assert len(objects) == 1, [p.name for p in objects]
    state = storage.storage_state_summary(str(run))
    assert state["state"] == "COLD"
    # The label from the first invocation is not lost by the second.
    assert state["label"] == "kept"

    # A real change must still produce a new archive rather than silently reusing.
    (run / "logs" / "late.log").write_text("appended after archiving", encoding="utf-8")
    changed = storage.archive_artifact(target, _profile(tmp_path), prune_local=False)
    assert changed["archive_id"] != protected["archive_id"]


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
        assert (case / "output" / "visualization" / "Field-0123456789ab" / "Field_00010.vts").exists() is False
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


def test_sweep_maintenance_handles_a_cold_member_without_losing_it(tmp_path, monkeypatch, capsys):
    """!
    @brief Sweep maintenance never reclassifies a deliberately offloaded member as empty.

    @details Continuation refuses, because it would have to rerun the member from
             scratch. Reaggregation does not: it cannot re-measure archived output, but
             the values it measured last time are still the right answer for that
             member, so they are carried forward instead of blanked.
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

    # A previous aggregation recorded a value for this member.
    results = study / "output" / "analysis"
    results.mkdir(parents=True, exist_ok=True)
    (results / "metrics_table.csv").write_text(
        f"case_id,msd_final\n{case.name},1.25\n", encoding="utf-8"
    )

    core.sweep_reaggregate_workflow(SimpleNamespace(study_dir=str(study), auto_fetch=False))
    assert "carried forward" in capsys.readouterr().out
    rows = list(csv.DictReader(
        (results / "metrics_table.csv").read_text(encoding="utf-8").splitlines()
    ))
    retained = {row["case_id"]: row for row in rows}[case.name]
    assert retained["msd_final"] == "1.25", rows


def test_sweep_auto_fetch_restores_a_cold_member_before_continuing(tmp_path, monkeypatch, local_rclone, capsys):
    """!
    @brief `--auto-fetch` restores a cold member first instead of refusing continuation.
    @details Without `--auto-fetch`, `sweep --continue` refuses a study with a cold
             member (covered by `test_sweep_maintenance_handles_a_cold_member_without_losing_it`).
             With it, the member is restored in place before the cold check is
             re-evaluated, so continuation proceeds like any other case.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] monkeypatch Value supplied through the `monkeypatch` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @param[in] capsys Captures workflow stdout.
    @return None.
    """
    monkeypatch.chdir(tmp_path)
    storage.storage_setup_workflow(SimpleNamespace(
        storage_config=str(tmp_path / storage.STORAGE_CONFIG_FILENAME),
        profile="archive", remote="fake:picurv-data", compression="none",
        chunk_size_gib=1.0, staging_directory=None, dry_run=False,
    ))
    core.sweep_workflow(SimpleNamespace(
        study=str(FIXTURES / "study.yml"),
        cluster=str(FIXTURES / "cluster.yml"),
        no_submit=True,
    ))
    study = next((tmp_path / "studies").iterdir())
    case = next((study / "cases").iterdir())

    profile = storage.load_storage_profile(None, None)
    target = storage.resolve_local_storage_targets(str(case), None)[0]
    storage.archive_artifact(target, profile, label="cold member", prune_local=True)
    assert storage.cold_study_members(str(study)) == [case.name]

    capsys.readouterr()
    core.sweep_continue_workflow(SimpleNamespace(
        study_dir=str(study), cluster=None, no_submit=True, auto_fetch=True,
    ))
    captured = capsys.readouterr()
    assert "--auto-fetch: restoring cold-storage member(s) before continuing" in captured.out
    assert "cold-storage member" not in captured.err
    assert storage.cold_study_members(str(study)) == []


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


def test_offload_policy_reports_and_retains_latest_checkpoint(tmp_path, local_rclone):
    """!
    @brief Plan reports semantic retention sizes and offload can keep only the latest checkpoint.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "retained")
    latest = run / "output" / "checkpoints" / "step_000000000020"
    shutil.copytree(run / "output" / "checkpoints" / "step_000000000010", latest)
    (run / "inputs" / "grid").mkdir(parents=True)
    (run / "inputs" / "grid" / "grid.run").write_bytes(b"grid-input")
    (run / "output" / "analysis" / "metrics").mkdir(parents=True)
    (run / "output" / "analysis" / "metrics" / "summary.json").write_text("{}\n", encoding="utf-8")
    (run / "output" / "visualization").mkdir(parents=True, exist_ok=True)
    (run / "output" / "visualization" / "field.vts").write_text("vtk\n", encoding="utf-8")
    target = storage.resolve_local_storage_targets(str(run), None)[0]

    plan = storage.build_storage_plan(
        target, _profile(tmp_path), policy="metadata-only", keep_latest_checkpoint=True,
    )
    assert plan["offload_policy"] == {
        "name": "metadata-only",
        "retained_components": ["logs", "metadata"],
        "keep_latest_checkpoint": True,
    }
    assert plan["retained_local_bytes"] > 0
    assert plan["pruned_local_bytes"] > 0
    storage.archive_artifact(
        target, _profile(tmp_path), prune_local=True,
        policy="metadata-only", keep_latest_checkpoint=True,
    )
    assert latest.is_dir()
    assert not (run / "output" / "checkpoints" / "step_000000000010").exists()
    assert not (run / "inputs" / "grid" / "grid.run").exists()
    assert (run / "config" / "case.yml").is_file()
    state = storage.read_storage_state(str(run))
    assert "checkpoint:20" in state["retained_components"]


def test_offload_policy_profiles_have_distinct_semantic_retention():
    """!
    @brief Named offload policies map to metadata, restart, and analysis use cases.
    @return None.
    """
    profile = {"offload_policy": "metadata-only", "keep_latest_checkpoint": False}
    metadata = storage._resolve_offload_policy(profile)
    restart = storage._resolve_offload_policy(profile, "restart-ready", None)
    analysis = storage._resolve_offload_policy(profile, "analysis-ready", True)
    assert metadata["retained_components"] == ["logs", "metadata"]
    assert restart["retained_components"] == ["inputs", "logs", "metadata"]
    # Nothing explicit was requested, so restart-ready's own promise to keep the
    # newest checkpoint applies.
    assert restart["keep_latest_checkpoint"] is True
    assert analysis["retained_components"] == ["analysis", "logs", "metadata", "visualization"]
    assert analysis["keep_latest_checkpoint"] is True


def test_selective_restore_can_reach_unclassified_and_workspace_components(tmp_path, local_rclone):
    """!
    @brief `--component` can select `unclassified`, `workspace-config`, and
           `workspace-inputs`, not only the six run-archive component names; a
           workspace archive's own config chunk is always included, mirroring how a
           run archive's `metadata` chunk already is.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "stray")
    (run / "NOTES.md").write_text("my own notes\n", encoding="utf-8")
    profile = _profile(tmp_path)
    manifest = storage.archive_artifact(storage.resolve_local_storage_targets(str(run), None)[0], profile)
    assert {c["component"] for c in manifest["chunks"]} >= {"metadata", "unclassified"}

    destination = tmp_path / "restored-unclassified"
    storage.restore_archive(profile, manifest["archive_id"], destination=str(destination), components=["unclassified"])
    assert (destination / "NOTES.md").read_text(encoding="utf-8") == "my own notes\n"
    # metadata is always included alongside a selected component.
    assert (destination / "config" / "case.yml").is_file()

    workspace = _write_workspace_artifact(tmp_path / "pilot64")
    workspace_manifest = storage.archive_artifact(
        storage.resolve_local_storage_targets(workspace=str(workspace))[0], profile
    )
    assert {c["component"] for c in workspace_manifest["chunks"]} >= {"workspace-config", "assets"}

    workspace_destination = tmp_path / "restored-workspace"
    storage.restore_archive(
        profile, workspace_manifest["archive_id"], destination=str(workspace_destination),
        components=["assets"],
    )
    # workspace-config is always included, the same way metadata is for run archives.
    assert (workspace_destination / storage.WORKSPACE_CONFIG_FILENAME).is_file()
    assert (workspace_destination / "assets" / "objects" / "grids" / ("a" * 64) / "asset.json").is_file()


def test_include_inputs_is_refused_without_workspace(tmp_path):
    """!
    @brief `--include-inputs` only means anything with `--workspace`; using it with
           `--run-dir`/`--study-dir` is refused rather than silently ignored.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "solo")
    with pytest.raises(storage.StorageError, match="--include-inputs is valid only with --workspace"):
        storage.resolve_local_storage_targets(str(run), None, include_inputs=True)


def test_explicit_drop_all_checkpoints_overrides_the_restart_ready_default():
    """!
    @brief `--drop-all-checkpoints` is not silently overruled by `restart-ready`.
    @details `restart-ready` keeps the newest checkpoint by default, but that default
             must yield to an explicit request, the same way an explicit
             `--keep-latest-checkpoint` already overrides `metadata-only`'s default of
             dropping it.
    @return None.
    """
    profile = {"offload_policy": "restart-ready"}
    dropped = storage._resolve_offload_policy(profile, "restart-ready", False)
    assert dropped["keep_latest_checkpoint"] is False
    kept = storage._resolve_offload_policy(profile, "metadata-only", True)
    assert kept["keep_latest_checkpoint"] is True


@pytest.mark.skipif(
    shutil.which("tar") is None or shutil.which("xz") is None,
    reason="parallel native maximum compression requires tar and xz",
)
def test_maximum_compression_uses_requested_native_cpu_workers(tmp_path):
    """!
    @brief Maximum compression passes the requested CPU count to native xz.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    root = tmp_path / "payload"
    root.mkdir()
    (root / "data.bin").write_bytes(b"repetitive-data" * 1000)
    destination = tmp_path / "chunk.tar.xz"
    compressor = storage._write_tar_chunk(
        str(root), {"entries": ["data.bin"]}, str(destination), "maximum", workers=3
    )
    assert compressor == "xz:3"
    extracted = tmp_path / "extracted"
    extracted.mkdir()
    storage._extract_chunk(str(destination), str(extracted))
    assert (extracted / "data.bin").read_bytes() == (root / "data.bin").read_bytes()


def test_python_archive_and_restore_parallelize_independent_chunks(
    tmp_path, local_rclone, monkeypatch
):
    """!
    @brief Fallback archive and restore use bounded worker pools across independent chunks.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] local_rclone Fake rclone object-store fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    observed = []
    real_executor = storage.concurrent.futures.ThreadPoolExecutor

    class RecordingExecutor(real_executor):
        """! @brief Thread pool that records requested concurrency for the test. """

        def __init__(self, max_workers=None, *args, **kwargs):
            """!
            @brief Record and delegate thread-pool construction.
            @param[in] max_workers Requested worker count.
            @param[in] args Additional positional arguments.
            @param[in] kwargs Additional keyword arguments.
            """
            observed.append(max_workers)
            super().__init__(max_workers=max_workers, *args, **kwargs)

    monkeypatch.setattr(storage.concurrent.futures, "ThreadPoolExecutor", RecordingExecutor)
    run = _write_run(tmp_path / "runs" / "parallel")
    profile = _profile(tmp_path)
    profile["chunk_size_bytes"] = 1
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    manifest = storage.archive_artifact(target, profile, workers=3)
    assert 3 in observed
    shutil.rmtree(run)
    storage.restore_archive(profile, manifest["archive_id"], destination=str(run), workers=2)
    assert 2 in observed
    assert (run / "config" / "case.yml").is_file()


def test_remote_catalog_recovers_missing_workspace_asset(tmp_path, local_rclone):
    """!
    @brief An archived run can reconstruct a deleted shared asset object by provider identity.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "asset-source")
    grid = run / "inputs" / "grid" / "grid.run"
    grid.parent.mkdir(parents=True)
    grid.write_bytes(b"canonical-grid")
    digest = hashlib.sha256(grid.read_bytes()).hexdigest()
    asset_id = "a" * 64
    spec_hash = "b" * 64
    (run / "inputs" / "assets.lock.yml").write_text(
        yaml.safe_dump({
            "assets": {
                "grid": {
                    "asset_id": asset_id,
                    "kind": "grid",
                    "provider": "file",
                    "provider_spec_sha256": spec_hash,
                    "files": [{"path": "inputs/grid/grid.run", "sha256": digest}],
                }
            }
        }, sort_keys=False),
        encoding="utf-8",
    )
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    storage.archive_artifact(target, _profile(tmp_path))

    workspace = tmp_path / "workspace"
    workspace.mkdir()
    (workspace / "assets").mkdir()
    (workspace / storage.STORAGE_CONFIG_FILENAME).write_text(
        yaml.safe_dump({
            "default_profile": "archive",
            "profiles": {
                "archive": {
                    "remote": "fake:picurv-data", "compression": "none",
                    "chunk_size_gib": 1.0, "workers": 2,
                }
            },
        }),
        encoding="utf-8",
    )
    restored = storage.restore_missing_workspace_assets(str(workspace), [spec_hash])
    assert restored[spec_hash]["asset_id"] == asset_id
    object_root = workspace / restored[spec_hash]["object"]
    assert (object_root / "asset.json").is_file()
    assert (object_root / "payload" / "inputs" / "grid" / "grid.run").read_bytes() == b"canonical-grid"


def test_unclassified_files_are_archived_but_never_pruned(tmp_path, local_rclone):
    """!
    @brief A file storage cannot classify is retained locally by every offload policy.

    @details The classifier's fallback used to be an ordinary component, which no policy
             retained, so offloading deleted any file storage did not recognise - a
             user's own notes beside a run, a tool's output, anything unregistered.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "unknown")
    stray = run / "NOTES.md"
    stray.write_text("my own notes about this run\n", encoding="utf-8")
    nested = run / "hand_analysis" / "fit.csv"
    nested.parent.mkdir()
    nested.write_text("x,y\n1,2\n", encoding="utf-8")

    target = storage.resolve_local_storage_targets(str(run), None)[0]
    inventory = storage.inspect_artifact(target, query_scheduler=False)
    unknown = {
        entry["path"] for entry in inventory["entries"]
        if entry["component"] == storage.UNCLASSIFIED_COMPONENT and entry["type"] != "directory"
    }
    assert unknown == {"NOTES.md", "hand_analysis/fit.csv"}

    for policy in storage.STORAGE_OFFLOAD_POLICIES:
        storage.archive_artifact(target, _profile(tmp_path), prune_local=True, policy=policy)
        assert stray.is_file(), f"{policy} deleted an unclassified file"
        assert nested.is_file(), f"{policy} deleted an unclassified file"
    # They are archived too, so a later restore brings them back.
    assert stray.read_text(encoding="utf-8") == "my own notes about this run\n"


def _write_workspace_artifact(root: Path) -> Path:
    """!
    @brief Create a workspace with configuration, an imported input, and one asset.
    @param[in] root Workspace directory to create.
    @return Created workspace path.
    """
    (root / "config").mkdir(parents=True)
    (root / "config" / "case.yml").write_text("title: ws\n", encoding="utf-8")
    (root / "inputs" / "grids").mkdir(parents=True)
    (root / "inputs" / "grids" / "imported.picgrid").write_bytes(b"imported-mesh")
    asset = root / "assets" / "objects" / "grids" / ("a" * 64)
    (asset / "payload").mkdir(parents=True)
    (asset / "asset.json").write_text('{"asset_id": "' + "a" * 64 + '"}\n', encoding="utf-8")
    (asset / "payload" / "grid.run").write_bytes(b"grid")
    (root / "assets" / "sets").mkdir(parents=True)
    (root / "runs").mkdir()
    (root / "studies").mkdir()
    (root / storage.WORKSPACE_CONFIG_FILENAME).write_text(
        "schema_version: 1\nworkspace:\n  id: pilot64\n", encoding="utf-8"
    )
    return root


def test_workspace_protection_covers_configuration_but_not_runs(tmp_path, local_rclone):
    """!
    @brief A workspace is its own artifact, and it does not swallow its runs.

    @details Runs and studies carry their own archives; duplicating them inside a
             workspace archive would store the same bytes twice and blur which object
             owns what. User-supplied inputs are excluded unless asked for, because
             protecting a workspace is a backup, not custody of the user's data.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    workspace = _write_workspace_artifact(tmp_path / "pilot64")
    _write_run(workspace / "runs" / "inner")

    target = storage.resolve_local_storage_targets(workspace=str(workspace))[0]
    assert target["artifact_type"] == "workspace"
    assert target["workspace_id"] == "pilot64"

    inventory = storage.inspect_artifact(target, query_scheduler=False)
    paths = {entry["path"] for entry in inventory["entries"]}
    assert not any(path.startswith("runs/") for path in paths), sorted(paths)
    assert "config/case.yml" in paths
    assert "assets/objects/grids/" + "a" * 64 + "/payload/grid.run" in paths
    assert "inputs/grids/imported.picgrid" not in paths

    with_inputs = storage.resolve_local_storage_targets(
        workspace=str(workspace), include_inputs=True
    )[0]
    inventory = storage.inspect_artifact(with_inputs, query_scheduler=False)
    assert "inputs/grids/imported.picgrid" in {e["path"] for e in inventory["entries"]}

    manifest = storage.archive_artifact(target, _profile(tmp_path), label="pilot workspace")
    assert manifest["workspace_id"] == "pilot64"
    assert manifest["workspace_assets"] == ["a" * 64]


def test_a_deleted_workspace_is_recoverable_by_its_identity(tmp_path, local_rclone):
    """!
    @brief Recovery does not depend on remembering an archive id or a directory name.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    workspace = _write_workspace_artifact(tmp_path / "pilot64")
    target = storage.resolve_local_storage_targets(workspace=str(workspace))[0]
    manifest = storage.archive_artifact(target, _profile(tmp_path), label="pilot workspace")
    shutil.rmtree(workspace)

    resolved = storage.resolve_workspace_archive_id(_profile(tmp_path), "pilot64")
    assert resolved == manifest["archive_id"]
    with pytest.raises(storage.StorageError, match="No workspace archive"):
        storage.resolve_workspace_archive_id(_profile(tmp_path), "no-such-workspace")

    recovered = tmp_path / "recovered"
    storage.restore_archive(_profile(tmp_path), resolved, destination=str(recovered))
    assert (recovered / "config" / "case.yml").read_text(encoding="utf-8") == "title: ws\n"
    assert (recovered / storage.WORKSPACE_CONFIG_FILENAME).is_file()


def test_workspace_offload_never_deletes_the_local_asset_store(tmp_path, local_rclone):
    """!
    @brief `storage offload --workspace .` must not unconditionally delete local assets.
    @details Reclaiming a shared asset is `storage prune --assets --unused-locally`'s job,
             deliberately, because that command checks whether any active local run still
             references the object first. A workspace offload that deleted assets outright
             would destroy an asset a currently-running case still needs, bypassing that
             reference-aware check entirely.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    workspace = _write_workspace_artifact(tmp_path / "pilot64")
    asset_id = "a" * 64
    run = _write_run(workspace / "runs" / "hot")
    (run / "inputs").mkdir(exist_ok=True)
    (run / "inputs" / "assets.lock.yml").write_text(
        f"assets:\n  grid:\n    asset_id: {asset_id}\n", encoding="utf-8"
    )

    asset_dir = workspace / "assets" / "objects" / "grids" / asset_id
    assert asset_dir.is_dir()
    target = storage.resolve_local_storage_targets(workspace=str(workspace))[0]
    storage.archive_artifact(target, _profile(tmp_path), prune_local=True)
    assert asset_dir.is_dir(), "workspace offload deleted an asset a local run still references"


def test_asset_pruning_is_reference_aware(tmp_path, local_rclone):
    """!
    @brief A shared asset survives locally while any active local run still needs it.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    workspace = _write_workspace_artifact(tmp_path / "pilot64")
    asset_id = "a" * 64
    for name in ("hot", "cold"):
        run = _write_run(workspace / "runs" / name)
        (run / "inputs").mkdir(exist_ok=True)
        (run / "inputs" / "assets.lock.yml").write_text(
            f"assets:\n  grid:\n    asset_id: {asset_id}\n", encoding="utf-8"
        )
    profile = _profile(tmp_path)
    storage.archive_artifact(
        storage.resolve_local_storage_targets(workspace=str(workspace))[0], profile
    )

    # Both runs local: the asset is needed and must not be removed.
    decisions = storage.prune_unused_workspace_assets(str(workspace), profile, dry_run=True)
    assert decisions[0]["active_local_runs"] == 2
    assert decisions[0]["local_removal"] == "blocked"

    # Offload both, and the local copy becomes removable while the remote keeps it.
    for name in ("hot", "cold"):
        storage.archive_artifact(
            storage.resolve_local_storage_targets(str(workspace / "runs" / name), None)[0],
            profile, prune_local=True,
        )
    decisions = storage.prune_unused_workspace_assets(str(workspace), profile, dry_run=True)
    assert decisions[0]["active_local_runs"] == 0
    assert decisions[0]["cold_runs"] == 2
    assert decisions[0]["remote_protection"] == "verified"
    assert decisions[0]["local_removal"] == "safe"
    assert (workspace / "assets" / "objects" / "grids" / asset_id).is_dir()

    storage.prune_unused_workspace_assets(str(workspace), profile)
    assert not (workspace / "assets" / "objects" / "grids" / asset_id).exists()


def test_status_completed_filters_study_members_like_other_actions(tmp_path, local_rclone, capsys):
    """!
    @brief `storage status --study-dir --completed` must select only finished members.
    @details `storage status` has its own per-member breakdown branch, built before
             `--completed` existed, that lists every `case_*` directory unconditionally.
             It must route through the same completed-member selection that `plan`,
             `protect`, and `offload` already use, not silently ignore the flag.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @param[in] capsys Captures the workflow's stdout.
    @return None.
    """
    study = _write_study(tmp_path / "grid")
    shutil.rmtree(study / "cases" / "case_0001" / "output" / "checkpoints")

    capsys.readouterr()
    storage.storage_status_workflow(SimpleNamespace(
        run_dir=None, study_dir=str(study), case_ids=[],
        workspace=None, include_inputs=False, completed=True,
        output_format="json",
    ))
    out = capsys.readouterr().out
    _info_line, payload = out.split("\n", 1)
    captured = json.loads(payload)
    assert {item["target"]["run_id"] for item in captured} == {"case_0000"}


def test_require_free_space_refuses_before_writing_anything(tmp_path):
    """!
    @brief The shared disk-space guard refuses when free space is below the estimate.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @return None.
    """
    storage._require_free_space(str(tmp_path), 1, "do something small")  # does not raise
    with pytest.raises(storage.StorageError, match="Not enough free space"):
        storage._require_free_space(str(tmp_path), 10**18, "do something huge")


def test_archive_refuses_to_stage_without_enough_free_space(tmp_path, monkeypatch, local_rclone):
    """!
    @brief `storage offload`/`protect` refuse to begin packaging a run the staging
           filesystem cannot hold, instead of discovering that mid-transfer.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] monkeypatch Value supplied through the `monkeypatch` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "big")
    target = storage.resolve_local_storage_targets(str(run), None)[0]
    monkeypatch.setattr(
        storage.operations.shutil, "disk_usage",
        lambda directory: SimpleNamespace(total=0, used=0, free=1),
    )
    with pytest.raises(storage.StorageError, match="Not enough free space"):
        storage.archive_artifact(target, _profile(tmp_path))


def test_manifest_carries_notes_and_a_captured_parameter_summary(tmp_path, local_rclone):
    """!
    @brief `--notes` and the auto-captured case parameter summary land in the manifest
           and are what `storage show` prints, so a run is identifiable months later
           without restoring anything.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "annotated")
    profile = _profile(tmp_path)
    manifest = storage.archive_artifact(
        storage.resolve_local_storage_targets(str(run), None)[0], profile,
        label="pilot 64", notes="Reference case for the pilot-64 series.",
    )
    assert manifest["notes"] == "Reference case for the pilot-64 series."
    assert manifest["parameters"]["start_step"] == 0
    assert manifest["parameters"]["total_steps"] == 10

    shown = storage.verify_remote_archive(profile, manifest["archive_id"])
    assert shown["notes"] == "Reference case for the pilot-64 series."
    assert shown["parameters"]["total_steps"] == 10


def test_storage_prune_cli_dispatch_removes_only_safe_objects(tmp_path, local_rclone, capsys):
    """!
    @brief `picurv storage prune --assets --unused-locally` reaches the real CLI parser
           and dispatch path, not only the underlying `prune_unused_workspace_assets`
           function it wraps.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @param[in] capsys Captures workflow stdout.
    @return None.
    """
    storage_config = tmp_path / storage.STORAGE_CONFIG_FILENAME
    storage.storage_setup_workflow(SimpleNamespace(
        storage_config=str(storage_config), profile="archive", remote="fake:picurv-data",
        compression="none", chunk_size_gib=1.0, staging_directory=None, dry_run=False,
    ))
    profile = storage.load_storage_profile("archive", str(storage_config))

    workspace = _write_workspace_artifact(tmp_path / "pilot64")
    asset_id = "a" * 64
    run = _write_run(workspace / "runs" / "cold")
    (run / "inputs").mkdir(exist_ok=True)
    (run / "inputs" / "assets.lock.yml").write_text(
        f"assets:\n  grid:\n    asset_id: {asset_id}\n", encoding="utf-8"
    )
    storage.archive_artifact(storage.resolve_local_storage_targets(workspace=str(workspace))[0], profile)
    storage.archive_artifact(
        storage.resolve_local_storage_targets(str(run), None)[0], profile, prune_local=True
    )

    asset_dir = workspace / "assets" / "objects" / "grids" / asset_id
    assert asset_dir.is_dir()

    parser = build_main_parser()
    args = parser.parse_args([
        "storage", "prune", "--workspace", str(workspace), "--assets", "--unused-locally",
        "--profile", "archive", "--storage-config", str(storage_config),
    ])
    dispatch_command(args)
    assert not asset_dir.exists()
    assert "Removed 1 local asset object" in capsys.readouterr().out


def test_checkpoint_range_selection_expands_to_steps():
    """!
    @brief `--checkpoints START:END[:STRIDE]` selects the steps it names.
    @return None.
    """
    assert storage._expand_checkpoint_selection(None, ["10:14"]) == [10, 11, 12, 13, 14]
    assert storage._expand_checkpoint_selection(None, ["100:500:100"]) == [100, 200, 300, 400, 500]
    assert storage._expand_checkpoint_selection([7], ["10:12"]) == [7, 10, 11, 12]
    assert storage._expand_checkpoint_selection(None, None) is None
    for bad in ("10", "10:9", "10:20:0", "a:b"):
        with pytest.raises(storage.StorageError):
            storage._expand_checkpoint_selection(None, [bad])


def test_a_marker_naming_no_archive_reports_broken(tmp_path):
    """!
    @brief A marker claiming a remote copy nobody can locate is not reported as safe.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @return None.
    """
    run = _write_run(tmp_path / "runs" / "broken")
    (run / storage.STORAGE_STATE_FILENAME).write_text(
        json.dumps({"storage_schema_version": 1, "local_pruned": True}) + "\n", encoding="utf-8"
    )
    assert storage.storage_state_summary(str(run))["state"] == "BROKEN"
    (run / storage.STORAGE_STATE_FILENAME).write_text("{ not json", encoding="utf-8")
    assert storage.storage_state_summary(str(run))["state"] == "BROKEN"


def test_identical_content_is_stored_once_across_archives(tmp_path, local_rclone):
    """!
    @brief Two archives holding the same bytes share one stored copy of them.

    @details Payload is content-addressed, so a checkpoint that appears in two runs, or
             a run archived twice under different labels, costs one copy rather than
             two. Manifests stay separate; only the payload is shared.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @return None.
    """
    profile = _profile(tmp_path)
    first = _write_run(tmp_path / "runs" / "first")
    # A faithful copy, so the payload really is identical: writing a second fixture
    # would differ by modification time alone and prove nothing about deduplication.
    second = tmp_path / "runs" / "second"
    shutil.copytree(first, second, copy_function=shutil.copy2)
    manifests = [
        storage.archive_artifact(
            storage.resolve_local_storage_targets(str(run), None)[0], profile
        )
        for run in (first, second)
    ]
    assert manifests[0]["archive_id"] != manifests[1]["archive_id"]

    blobs = local_rclone / "picurv-data" / "blobs"
    stored = {path.name for path in blobs.rglob("*") if path.is_file()}
    digests = {
        chunk["sha256"] for manifest in manifests for chunk in manifest["chunks"]
    }
    assert stored == digests
    # The two runs hold the same bytes, so their checkpoint payload is stored once.
    checkpoints = [
        chunk["sha256"] for manifest in manifests for chunk in manifest["chunks"]
        if chunk["component"].startswith("checkpoint:")
    ]
    assert len(checkpoints) == 2 and len(set(checkpoints)) == 1


def test_an_interrupted_upload_resumes_without_re_uploading(tmp_path, local_rclone, monkeypatch, capsys):
    """!
    @brief Rerunning after a failed upload skips the chunks that already landed.
    @param[in] tmp_path Value supplied through the `tmp_path` argument.
    @param[in] local_rclone Fake rclone object-store fixture.
    @param[in] monkeypatch Value supplied through the `monkeypatch` argument.
    @param[in] capsys Value supplied through the `capsys` argument.
    @return None.
    """
    profile = _profile(tmp_path)
    run = _write_run(tmp_path / "runs" / "interrupted")
    target = storage.resolve_local_storage_targets(str(run), None)[0]

    real_upload = storage.transport._upload_verified
    calls = {"count": 0}

    def failing_upload(local_path, remote_path):
        """!
        @brief Upload two chunks, then fail as a dropped transfer would.
        @param[in] local_path Chunk being uploaded.
        @param[in] remote_path Destination path.
        @return Upload result for the calls that succeed.
        """
        calls["count"] += 1
        if calls["count"] > 2:
            raise storage.StorageError("simulated transfer failure")
        return real_upload(local_path, remote_path)

    monkeypatch.setattr(storage.transport, "_upload_verified", failing_upload)
    with pytest.raises(storage.StorageError, match="simulated transfer failure"):
        storage.archive_artifact(target, profile)
    landed = {path.name for path in (local_rclone / "picurv-data" / "blobs").rglob("*") if path.is_file()}
    assert len(landed) == 2

    # Local data is untouched by the failure, and the rerun uploads only what is missing.
    assert (run / "output" / "checkpoints" / "step_000000000010" / "COMMITTED").is_file()
    monkeypatch.setattr(storage.transport, "_upload_verified", real_upload)
    capsys.readouterr()
    manifest = storage.archive_artifact(target, profile)
    assert "already stored; skipping upload" in capsys.readouterr().out
    assert {chunk["sha256"] for chunk in manifest["chunks"]} >= landed
