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
    read_artifact_identity,
    RUN_MANIFEST_FILENAME,
    STUDY_MANIFEST_FILENAME,
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
    study_identity = read_artifact_identity(study_root)
    selected = []
    for name in sorted(os.listdir(cases_dir)):
        member = os.path.join(cases_dir, name)
        if not os.path.isdir(member):
            continue
        # A member is a member because it carries a run manifest, not because its
        # directory is named `case_NNNN`. The name pattern is the fallback for a
        # member staged before manifests recorded study membership.
        member_identity = read_artifact_identity(member)
        if member_identity["identity_source"] != "manifest" and not re.fullmatch(r"case_\d+", name):
            continue
        target = {
            "artifact_type": "study-case", "root_path": member, "original_path": member,
            "run_id": member_identity["run_id"] or name,
            "study_id": study_identity["study_id"],
            "case_id": member_identity["case_id"] or name,
            "study_path": study_root,
            "identity_source": member_identity["identity_source"],
        }
        try:
            inventory = inspect_artifact(target)
            _assert_archive_safe(inventory)
        except StorageError:
            continue
        if inventory["checkpoint_steps"]:
            selected.append(name)
    return selected


def _select_completed_case_ids(study_root: str) -> list:
    """!
    @brief Resolve `--completed` into member ids, raising if none qualify.
    @param[in] study_root Study directory to scan.
    @return Sorted member ids; never empty.
    """
    case_ids = _completed_study_members(study_root)
    if not case_ids:
        raise StorageError(
            f"No completed member found in {study_root}. A member is complete when it "
            "holds a committed checkpoint and no run is active on it."
        )
    print("[INFO] Selected completed member(s): " + ", ".join(case_ids))
    return case_ids


def _validate_case_id(case_id: str) -> str:
    """!
    @brief Validate a canonical numbered study-member identifier.

    @details This checks a value the user typed before it becomes a path segment, so it
             stays a pattern match on purpose: it is what stops `--case-id ../..` from
             resolving outside the study. What a member *is* is still decided by its
             manifest, not by this.
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
    if include_inputs and not workspace:
        raise StorageError("--include-inputs is valid only with --workspace.")
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
        # The run's own manifest says what it is and what it is called. A directory
        # with no manifest is still accepted - runs staged by an older release have
        # none - but it must at least carry the config snapshot every run writes, and
        # the fallback identity is recorded rather than passed off as the real one.
        identity = read_artifact_identity(root)
        if identity["identity_source"] != "manifest" and not os.path.isdir(
                os.path.join(root, "config")):
            raise StorageError(
                f"Directory does not look like a PICurv run (no {RUN_MANIFEST_FILENAME} "
                f"and no config/): {root}"
            )
        if identity["artifact_type"] == "study":
            raise StorageError(f"That is a study, not a run; use --study-dir: {root}")
        return [{
            "artifact_type": "run",
            "root_path": root,
            "original_path": root,
            "run_id": identity["run_id"],
            "study_id": identity["study_id"] if identity["case_id"] else None,
            "case_id": identity["case_id"],
            "identity_source": identity["identity_source"],
        }]

    study_root = os.path.abspath(study_dir)
    if not os.path.isdir(study_root):
        raise StorageError(f"Study directory not found: {study_root}")
    study_identity = read_artifact_identity(study_root)
    if study_identity["identity_source"] != "manifest" and not os.path.isdir(
            os.path.join(study_root, "cases")):
        raise StorageError(
            f"Directory does not look like a PICurv study (no {STUDY_MANIFEST_FILENAME} "
            f"and no cases/): {study_root}"
        )
    if not os.path.isdir(os.path.join(study_root, "cases")):
        raise StorageError(f"Study has no cases/ directory: {study_root}")
    if completed:
        if case_ids:
            raise StorageError("--completed selects members itself; do not also pass --case-id.")
        case_ids = _select_completed_case_ids(study_root)
    if case_ids:
        targets = []
        for raw_case_id in case_ids:
            case_id = _validate_case_id(raw_case_id)
            case_root = os.path.join(study_root, "cases", case_id)
            if not os.path.isdir(case_root):
                raise StorageError(f"Study member not found: {case_root}")
            member_identity = read_artifact_identity(case_root)
            targets.append({
                "artifact_type": "study-case",
                "root_path": case_root,
                "original_path": case_root,
                "run_id": member_identity["run_id"] or case_id,
                "study_id": study_identity["study_id"],
                "case_id": member_identity["case_id"] or case_id,
                "study_path": study_root,
                "identity_source": member_identity["identity_source"],
            })
        return targets
    return [{
        "artifact_type": "study",
        "root_path": study_root,
        "original_path": study_root,
        "run_id": None,
        "study_id": study_identity["study_id"],
        "case_id": None,
        "identity_source": study_identity["identity_source"],
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


#: Component paths assumed when an artifact carries no manifest of its own: a run
#: staged by an older release, or a directory restored without its identity. These are
#: the fixed workspace topology the conductor writes, kept here so a manifest-less
#: artifact still classifies rather than becoming wholly "unclassified" and unprunable.
FALLBACK_COMPONENT_PATHS = {
    "config": "config",
    "scheduler": "scheduler",
    "inputs": "inputs",
    "output": "output",
    "checkpoints": "output/checkpoints",
    "analysis": "output/analysis",
    "visualization": "output/visualization",
    "logs": "logs",
}


def _artifact_component_layout(root: str, artifact_type: str) -> dict:
    """!
    @brief Read from the artifact's own manifests where each component lives.

    @details Storage classifies what it packages by asking the run what its directories
             are, not by assuming the topology from path prefixes. A run records its
             component map in `manifest.json` under `paths`, and a study's members each
             record their own, so a renamed, restored, or re-rooted artifact classifies
             the same way it did where it was written.

             Members are found by the manifest each one carries, not by matching a
             `case_NNNN` name. When a manifest is missing the fixed workspace topology
             is used and reported through `identity_source`, so a caller can tell an
             answer from a fallback.
    @param[in] root Absolute artifact root.
    @param[in] artifact_type "run", "study", "study-case", or "workspace".
    @return Mapping with `members` (member-prefix to component-path map) and
            `identity_source`.
    """
    root_abs = os.path.abspath(root)

    def component_paths(member_root: str) -> tuple:
        """!
        @brief Return one artifact's declared component map and where it came from.
        @param[in] member_root Absolute run or study-member root.
        @return Tuple of (component path mapping, identity source).
        """
        identity = read_artifact_identity(member_root)
        declared = identity.get("paths") or {}
        # A study manifest's `paths` names its scripts and tables, not the run-relative
        # component tree, so it is not a component map and must not be read as one.
        usable = {
            key: str(value).strip("/")
            for key, value in declared.items()
            if isinstance(value, str) and value and not os.path.isabs(value)
        }
        if identity["identity_source"] == "manifest" and usable:
            return usable, "manifest"
        return dict(FALLBACK_COMPONENT_PATHS), "fixed-topology"

    members = {}
    sources = set()
    if artifact_type == "study":
        cases_root = os.path.join(root_abs, "cases")
        if os.path.isdir(cases_root):
            for name in sorted(os.listdir(cases_root)):
                member_root = os.path.join(cases_root, name)
                if not os.path.isdir(member_root):
                    continue
                paths, source = component_paths(member_root)
                members[f"cases/{name}"] = paths
                sources.add(source)
    paths, source = component_paths(root_abs)
    members[""] = paths
    sources.add(source)
    return {
        "root": root_abs,
        "members": members,
        "identity_source": "manifest" if sources == {"manifest"} else "mixed"
        if "manifest" in sources else "fixed-topology",
    }


def _checkpoint_step_from_bundle(bundle_root: str):
    """!
    @brief Read a committed checkpoint's own recorded step number.

    @details The bundle writes `-checkpoint_step` into `checkpoint.meta`, which is what
             the step *is*; the directory name is a rendering of it. Reading the record
             keeps the classification on the same footing as the rest of this module,
             and stays correct for a bundle restored under a different name.
    @param[in] bundle_root Absolute path to one checkpoint bundle directory.
    @return The recorded step as an int, or None when it cannot be read.
    """
    try:
        with open(os.path.join(bundle_root, "checkpoint.meta"), "r", encoding="utf-8") as stream:
            for line in stream:
                tokens = line.split("#", 1)[0].split()
                if len(tokens) == 2 and tokens[0] == "-checkpoint_step":
                    return int(tokens[1])
    except (OSError, ValueError):
        return None
    return None


def _classify_workspace_component(relative_path: str) -> str:
    """!
    @brief Classify one workspace-owned path for packaging and retention.

    @details A workspace has no manifest of components the way a run does: its shape is
             the fixed layout `picurv init` writes, and `.picurv-workspace.yml` records
             identity rather than topology. So this stays a layout decision, and the
             three named roots are the ones the workspace contract defines.
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


def _classify_component(relative_path: str, layout: dict = None) -> str:
    """!
    @brief Classify one artifact path for packaging and local retention.
    @param[in] relative_path Artifact-relative path of the entry.
    @param[in] layout Component layout from `_artifact_component_layout()`. When absent
               the fixed workspace topology is assumed.
    @return Component name, or a `checkpoint:<step>` token.
    """
    layout = layout or {"root": "", "members": {"": dict(FALLBACK_COMPONENT_PATHS)}}
    members = layout["members"]
    # Longest member prefix wins, so a study member's own map beats the study's.
    prefix = ""
    for candidate in members:
        if not candidate:
            continue
        if relative_path == candidate or relative_path.startswith(candidate + "/"):
            if len(candidate) > len(prefix):
                prefix = candidate
    paths = members.get(prefix) or dict(FALLBACK_COMPONENT_PATHS)
    local = relative_path[len(prefix):].lstrip("/") if prefix else relative_path
    base = os.path.basename(relative_path)

    def under(component_key: str) -> bool:
        """!
        @brief Whether the entry lies at or below one declared component path.
        @param[in] component_key Component name in the artifact's own path map.
        @return True when the entry belongs to that component.
        """
        declared = paths.get(component_key)
        if not declared:
            return False
        return local == declared or local.startswith(declared + "/")

    # The locks record which assets and which executables a run consumed. They are the
    # run's provenance, not its payload: pruning them would leave a cold run unable to
    # say what it was built from, and reference-aware asset pruning unable to see it.
    if base in {
        "manifest.json", "study_manifest.json", "study.yml", "cluster.yml",
        "assets.lock.yml", "software.lock.json",
    }:
        return "metadata"
    if under("config") or under("scheduler"):
        return "metadata"
    if under("checkpoints"):
        declared = paths["checkpoints"]
        remainder = local[len(declared):].lstrip("/")
        bundle = remainder.split("/", 1)[0]
        if bundle:
            step = _checkpoint_step_from_bundle(
                os.path.join(layout["root"], *(prefix.split("/") if prefix else []),
                             *declared.split("/"), bundle)
            )
            if step is None:
                # The bundle's own record is the authority, but a checkpoint whose
                # metadata cannot be read must not thereby stop being a checkpoint:
                # that would move it into a component an ordinary policy prunes. Fall
                # back to the name it was written under.
                name_match = CHECKPOINT_DIRECTORY_PATTERN.fullmatch(bundle)
                if name_match:
                    step = int(name_match.group(1))
            if step is not None:
                return f"checkpoint:{step}"
            # Neither the record nor the name identifies a step. Leave it unclassified,
            # which is archived and never pruned, rather than guessing.
            return UNCLASSIFIED_COMPONENT
        return "raw-output"
    if under("analysis"):
        return "analysis"
    if under("visualization"):
        return "visualization"
    if under("output"):
        return "raw-output"
    if under("inputs"):
        return "inputs"
    if under("logs"):
        return "logs"
    if local.split("/")[0] == "assets":
        return "assets"
    if local.split("/")[0] == "results":
        return "analysis"
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
        component = entry.get("component") or ""
        if component.startswith("checkpoint:"):
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
    layout = None if workspace else _artifact_component_layout(root, target["artifact_type"])
    for entry in entries:
        entry["component"] = (
            _classify_workspace_component(entry["path"]) if workspace
            else _classify_component(entry["path"], layout)
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
        "component_layout_source": layout["identity_source"] if layout else "workspace-contract",
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
