"""!
@file inventory.py
@brief Discovering what an artifact contains and which target was selected.
"""

import argparse
import base64
import concurrent.futures
import contextlib
import datetime
import errno
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tarfile
import tempfile
import uuid
from pathlib import Path
import yaml
from .models import (
    CHECKPOINT_DIRECTORY_PATTERN,
    INCOMPLETE_CHECKPOINT_PATTERN,
    STORAGE_LOCK_FILENAME,
    STORAGE_STATE_FILENAME,
    StorageError,
    UNCLASSIFIED_COMPONENT,
    WORKSPACE_CONFIG_FILENAME,
    WORKSPACE_EXCLUDED_ROOTS,
    storage_state_summary,
)
from .safety import (
    _assert_archive_safe,
    _lock_owner_active,
    _slurm_activity,
)


def _completed_study_members(study_root: str) -> list:
    """!
    @brief Select the study members that are finished and idle.

    @details "Finished" is deliberately conservative: a member qualifies only when it
             holds at least one committed checkpoint, is writing none, and has no
             active lock or scheduler job. Anything ambiguous is skipped rather than
             archived mid-flight.
    @param[in] study_root Study directory to scan.
    @return Sorted member ids.
    """
    cases_dir = os.path.join(study_root, "cases")
    selected = []
    for name in sorted(os.listdir(cases_dir)):
        member = os.path.join(cases_dir, name)
        if not os.path.isdir(member) or not re.fullmatch(r"case_\d+", name):
            continue
        target = {
            "artifact_type": "study-case", "root_path": member, "original_path": member,
            "run_id": name, "study_id": os.path.basename(study_root), "case_id": name,
            "study_path": study_root,
        }
        try:
            inventory = inspect_artifact(target)
            _assert_archive_safe(inventory)
        except StorageError:
            continue
        if inventory["checkpoint_steps"]:
            selected.append(name)
    return selected


def _validate_case_id(case_id: str) -> str:
    """!
    @brief Validate a canonical numbered study-member identifier.
    @param[in] case_id Value supplied through the `case_id` argument.
    @return Result produced by this operation.
    """
    if not re.fullmatch(r"case_\d+", str(case_id or "")):
        raise StorageError(f"Invalid study case ID {case_id!r}; expected a value such as case_0003.")
    return str(case_id)


def resolve_local_storage_targets(run_dir: str = None, study_dir: str = None, case_ids=None,
                                  workspace: str = None, include_inputs: bool = False,
                                  completed: bool = False) -> list:
    """!
    @brief Resolve explicit run/study/workspace selectors into artifact descriptions.
    @param[in] run_dir Value supplied through the `run_dir` argument.
    @param[in] study_dir Value supplied through the `study_dir` argument.
    @param[in] case_ids Value supplied through the `case_ids` argument.
    @param[in] workspace Workspace root to protect as its own artifact.
    @param[in] include_inputs Whether workspace protection covers user-supplied inputs.
    @param[in] completed Select every finished study member instead of naming them.
    @return Result produced by this operation.
    """
    selectors = [bool(run_dir), bool(study_dir), bool(workspace)]
    if sum(selectors) != 1:
        raise StorageError("Select exactly one of --run-dir, --study-dir, or --workspace.")
    if workspace:
        if case_ids:
            raise StorageError("--case-id is valid only with --study-dir.")
        root = os.path.abspath(workspace)
        if not os.path.isfile(os.path.join(root, WORKSPACE_CONFIG_FILENAME)):
            raise StorageError(
                f"Directory is not an initialized PICurv workspace (missing "
                f"{WORKSPACE_CONFIG_FILENAME}): {root}"
            )
        identity = {}
        try:
            with open(os.path.join(root, WORKSPACE_CONFIG_FILENAME), "r", encoding="utf-8") as stream:
                identity = (yaml.safe_load(stream) or {}).get("workspace") or {}
        except (OSError, ValueError):
            identity = {}
        return [{
            "artifact_type": "workspace",
            "root_path": root,
            "original_path": root,
            "run_id": None,
            "study_id": None,
            "case_id": None,
            "workspace_id": identity.get("id") or os.path.basename(root),
            "include_inputs": bool(include_inputs),
        }]
    if run_dir:
        if case_ids:
            raise StorageError("--case-id is valid only with --study-dir.")
        root = os.path.abspath(run_dir)
        if not os.path.isdir(root):
            raise StorageError(f"Run directory not found: {root}")
        if not os.path.isdir(os.path.join(root, "config")):
            raise StorageError(f"Directory does not look like a PICurv run (missing config/): {root}")
        return [{
            "artifact_type": "run",
            "root_path": root,
            "original_path": root,
            "run_id": os.path.basename(root),
            "study_id": None,
            "case_id": None,
        }]

    study_root = os.path.abspath(study_dir)
    if not os.path.isdir(study_root):
        raise StorageError(f"Study directory not found: {study_root}")
    if not os.path.isdir(os.path.join(study_root, "cases")):
        raise StorageError(f"Directory does not look like a PICurv study (missing cases/): {study_root}")
    if completed:
        if case_ids:
            raise StorageError("--completed selects members itself; do not also pass --case-id.")
        case_ids = _completed_study_members(study_root)
        if not case_ids:
            raise StorageError(
                f"No completed member found in {study_root}. A member is complete when it "
                "holds a committed checkpoint and no run is active on it."
            )
        print("[INFO] Selected completed member(s): " + ", ".join(case_ids))
    if case_ids:
        targets = []
        for raw_case_id in case_ids:
            case_id = _validate_case_id(raw_case_id)
            case_root = os.path.join(study_root, "cases", case_id)
            if not os.path.isdir(case_root):
                raise StorageError(f"Study member not found: {case_root}")
            targets.append({
                "artifact_type": "study-case",
                "root_path": case_root,
                "original_path": case_root,
                "run_id": case_id,
                "study_id": os.path.basename(study_root),
                "case_id": case_id,
                "study_path": study_root,
            })
        return targets
    return [{
        "artifact_type": "study",
        "root_path": study_root,
        "original_path": study_root,
        "run_id": None,
        "study_id": os.path.basename(study_root),
        "case_id": None,
    }]


def _path_is_within(root: str, candidate: str) -> bool:
    """!
    @brief Return whether an absolute candidate remains within a root directory.
    @param[in] root Value supplied through the `root` argument.
    @param[in] candidate Value supplied through the `candidate` argument.
    @return Result produced by this operation.
    """
    try:
        return os.path.commonpath([os.path.abspath(root), os.path.abspath(candidate)]) == os.path.abspath(root)
    except ValueError:
        return False


def _resolve_configured_path(root: str, value: str) -> str:
    """!
    @brief Resolve runtime-directory syntax using the run directory as working directory.
    @param[in] root Value supplied through the `root` argument.
    @param[in] value Value supplied through the `value` argument.
    @return Result produced by this operation.
    """
    return os.path.abspath(value if os.path.isabs(value) else os.path.join(root, value))


def _artifact_runtime_roots(root: str) -> list:
    """!
    @brief Return run-like roots contained by a standalone run or whole study.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    root_path = Path(os.path.abspath(root))
    case_root = root_path / "cases"
    if case_root.is_dir():
        return [
            str(path) for path in sorted(case_root.glob("case_*"))
            if path.is_dir()
        ]
    return [str(root_path)]


def _discover_external_paths(root: str) -> list:
    """!
    @brief Report explicit external-reference descriptors inside an artifact.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    external = []
    archive_root = os.path.abspath(root)
    for descriptor in Path(archive_root).glob("**/*.reference.yml"):
        try:
            with descriptor.open("r", encoding="utf-8") as stream:
                payload = yaml.safe_load(stream) or {}
            target = payload.get("picurv_external_reference")
            if isinstance(target, str) and os.path.isabs(target):
                external.append({
                    "source": os.path.relpath(descriptor, archive_root).replace(os.sep, "/"),
                    "path": os.path.abspath(target),
                })
        except (OSError, ValueError, TypeError):
            pass
    return sorted(external, key=lambda item: (item["source"], item["path"]))


def _discover_dependencies(root: str) -> list:
    """!
    @brief Discover absolute restart/source paths embedded in generated controls.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    dependencies = []
    archive_root = os.path.abspath(root)
    controls = []
    for runtime_root in _artifact_runtime_roots(root):
        controls.extend(sorted(Path(runtime_root).glob("config/*.control")))
    for control in controls:
        try:
            lines = control.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for line in lines:
            stripped = line.strip()
            if not stripped.startswith("-restart_dir "):
                continue
            try:
                tokens = __import__("shlex").split(stripped)
            except ValueError:
                continue
            if len(tokens) >= 2 and os.path.isabs(tokens[1]) and not _path_is_within(archive_root, tokens[1]):
                dependencies.append({"kind": "restart", "path": os.path.abspath(tokens[1])})
    return dependencies


def _walk_archive_entries(root: str, excluded_roots=()) -> list:
    """!
    @brief Enumerate archive entries without following symlinks.
    @param[in] root Value supplied through the `root` argument.
    @param[in] excluded_roots Top-level directory names to skip entirely.
    @return Result produced by this operation.
    """
    root_abs = os.path.abspath(root)
    excluded = set(excluded_roots or ())
    entries = []

    def visit(directory: str):
        """!
        @brief Recursively inventory directory entries without following symlinks.
        @param[in] directory Value supplied through the `directory` argument.
        """
        try:
            children = sorted(os.scandir(directory), key=lambda item: item.name)
        except OSError as exc:
            raise StorageError(f"Unable to inventory {directory}: {exc}") from exc
        for child in children:
            if child.name in {STORAGE_STATE_FILENAME, STORAGE_LOCK_FILENAME}:
                continue
            rel = os.path.relpath(child.path, root_abs).replace(os.sep, "/")
            if rel in excluded:
                continue
            try:
                stat_result = child.stat(follow_symlinks=False)
            except OSError as exc:
                raise StorageError(f"Unable to stat {child.path}: {exc}") from exc
            if child.is_symlink():
                entry_type = "symlink"
                size = 0
            elif child.is_dir(follow_symlinks=False):
                entry_type = "directory"
                size = 0
            elif child.is_file(follow_symlinks=False):
                entry_type = "file"
                size = int(stat_result.st_size)
            else:
                raise StorageError(f"Unsupported filesystem entry in artifact: {child.path}")
            entries.append({
                "path": rel,
                "type": entry_type,
                "size": size,
                "mode": int(stat_result.st_mode & 0o7777),
                "mtime_ns": int(stat_result.st_mtime_ns),
            })
            if entry_type == "directory":
                visit(child.path)

    visit(root_abs)
    return entries


def _checkpoint_component(relative_path: str):
    """!
    @brief Return a checkpoint component token for a path inside a committed step bundle.
    @param[in] relative_path Value supplied through the `relative_path` argument.
    @return Result produced by this operation.
    """
    parts = relative_path.split("/")
    for index, part in enumerate(parts):
        match = CHECKPOINT_DIRECTORY_PATTERN.fullmatch(part)
        if match and index > 0 and parts[index - 1] == "checkpoints":
            return f"checkpoint:{int(match.group(1))}"
    return None


def _classify_workspace_component(relative_path: str) -> str:
    """!
    @brief Classify one workspace-owned path for packaging and retention.
    @param[in] relative_path Workspace-relative path.
    @return Component name.
    """
    first = relative_path.split("/")[0]
    if first == "config" or relative_path.startswith("."):
        return "workspace-config"
    if first == "inputs":
        return "workspace-inputs"
    if first == "assets":
        return "assets"
    return UNCLASSIFIED_COMPONENT


def _classify_component(relative_path: str) -> str:
    """!
    @brief Classify one artifact path for packaging and local retention.
    @param[in] relative_path Value supplied through the `relative_path` argument.
    @return Result produced by this operation.
    """
    checkpoint = _checkpoint_component(relative_path)
    if checkpoint:
        return checkpoint
    parts = relative_path.split("/")
    component_parts = parts
    first = component_parts[0]
    base = os.path.basename(relative_path)
    if len(parts) >= 3 and parts[0] == "cases" and re.fullmatch(r"case_\d+", parts[1]):
        component_parts = parts[2:]
        first = component_parts[0]
    # The locks record which assets and which executables a run consumed. They are the
    # run's provenance, not its payload: pruning them would leave a cold run unable to
    # say what it was built from, and reference-aware asset pruning unable to see it.
    if first in {"config", "scheduler"} or base in {
        "manifest.json", "study_manifest.json", "study.yml", "cluster.yml",
        "assets.lock.yml", "software.lock.json",
    }:
        return "metadata"
    if first == "inputs":
        return "inputs"
    if first == "assets":
        return "assets"
    if first == "output":
        nested = component_parts[1] if len(component_parts) > 1 else ""
        if nested == "analysis":
            return "analysis"
        if nested == "visualization":
            return "visualization"
        return "raw-output"
    if first == "results":
        return "analysis"
    if first == "logs":
        return "logs"
    # Nothing recognized this path. It is archived like everything else, but it is
    # never pruned: storage must not delete a file whose purpose it cannot state.
    return UNCLASSIFIED_COMPONENT


def _checkpoint_steps(entries: list) -> list:
    """!
    @brief Return committed checkpoint steps represented in an inventory.
    @param[in] entries Value supplied through the `entries` argument.
    @return Result produced by this operation.
    """
    candidates = {}
    for entry in entries:
        component = _checkpoint_component(entry["path"])
        if component:
            step = int(component.split(":", 1)[1])
            candidates.setdefault(step, set()).add(os.path.basename(entry["path"]))
    return sorted(step for step, names in candidates.items() if {"checkpoint.meta", "COMMITTED"} <= names)


def _find_incomplete_checkpoints(entries: list) -> list:
    """!
    @brief Return incomplete checkpoint paths that make archival unsafe.
    @param[in] entries Value supplied through the `entries` argument.
    @return Result produced by this operation.
    """
    return [
        entry["path"] for entry in entries
        if any(INCOMPLETE_CHECKPOINT_PATTERN.match(part) for part in entry["path"].split("/"))
    ]


def inspect_artifact(target: dict, query_scheduler: bool = True) -> dict:
    """!
    @brief Build a read-only inventory and lifecycle assessment for one artifact.
    @param[in] target Value supplied through the `target` argument.
    @param[in] query_scheduler Value supplied through the `query_scheduler` argument.
    @return Result produced by this operation.
    """
    root = target["root_path"]
    workspace = target["artifact_type"] == "workspace"
    entries = _walk_archive_entries(
        root, WORKSPACE_EXCLUDED_ROOTS if workspace else ()
    )
    for entry in entries:
        entry["component"] = (
            _classify_workspace_component(entry["path"]) if workspace
            else _classify_component(entry["path"])
        )
    if workspace and not target.get("include_inputs"):
        # User-supplied data is theirs. It is archived only on request, because a
        # workspace protect is about the configuration and catalog, not about taking
        # custody of possibly enormous imported fields.
        entries = [
            entry for entry in entries
            if entry["component"] != "workspace-inputs" or entry["type"] == "directory"
        ]
    lock_paths = []
    for name in ("post.lock.json", "solver.lock.json"):
        for candidate in Path(root).glob(f"**/scheduler/{name}"):
            if _lock_owner_active(str(candidate)):
                lock_paths.append(os.path.relpath(candidate, root))
    slurm = _slurm_activity(root) if query_scheduler else {"job_ids": [], "active": [], "unknown": False}
    state = storage_state_summary(root)
    return {
        "target": dict(target),
        "entries": entries,
        "file_count": sum(entry["type"] in {"file", "symlink"} for entry in entries),
        "total_bytes": sum(entry["size"] for entry in entries),
        "checkpoint_steps": _checkpoint_steps(entries),
        "incomplete_checkpoints": _find_incomplete_checkpoints(entries),
        "active_locks": sorted(lock_paths),
        "slurm": slurm,
        "external_paths": _discover_external_paths(root),
        "dependencies": _discover_dependencies(root),
        "storage": state,
    }


def _inventory_fingerprint(inventory: dict) -> str:
    """!
    @brief Stable digest of the file set an archive would package.

    @details Built from each entry's run-relative path, size, and modification time -
             the same signals ordinary change detection uses - so it can be computed
             without reading file contents. Two archives of the same artifact agree on
             it exactly when nothing has been written since, which is what makes a
             re-upload detectably redundant.
    @param[in] inventory Inventory returned by `inspect_artifact()`.
    @return Hex digest over the normalized entry list.
    """
    digest = hashlib.sha256()
    for entry in sorted(inventory["entries"], key=lambda item: item["path"]):
        digest.update(
            f"{entry['path']}\0{entry.get('kind')}\0{entry.get('bytes')}\0"
            f"{entry.get('mtime_ns')}\n".encode("utf-8")
        )
    return digest.hexdigest()
