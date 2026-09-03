#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""!
@file storage.py
@brief Lifecycle-aware archival, offload, verification, and restore workflows.

The storage layer deliberately remains independent of the numerical runtime.  It
packages immutable run/study artifacts, uses rclone as a transport boundary, and
leaves a small local state marker whenever payload data is pruned.
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


STORAGE_CONFIG_FILENAME = ".picurv-storage.yml"
STORAGE_STATE_FILENAME = ".picurv-storage.json"
STORAGE_LOCK_FILENAME = ".picurv-storage.lock.json"
#: Schema 2 stores chunk payloads in a content-addressed blob store shared by every
#: archive, so identical content is uploaded once and an interrupted upload resumes.
STORAGE_SCHEMA_VERSION = 2
REMOTE_OBJECTS_DIRECTORY = "objects"

#: Content-addressed payload store, shared across every archive on this remote.
REMOTE_BLOBS_DIRECTORY = "blobs"
REMOTE_MANIFEST_FILENAME = "manifest.json"
REMOTE_COMPLETE_FILENAME = "COMPLETE"
DEFAULT_PROFILE_NAME = "archive"
DEFAULT_CHUNK_SIZE_GIB = 8.0
DEFAULT_STORAGE_WORKERS = max(1, min(8, os.cpu_count() or 1))
AUTO_NO_COMPRESSION_BYTES = 256 * 1024 * 1024
AUTO_MAXIMUM_COMPRESSION_BYTES = 20 * 1024 * 1024 * 1024
CHECKPOINT_DIRECTORY_PATTERN = re.compile(r"^step_(\d{12})$")
INCOMPLETE_CHECKPOINT_PATTERN = re.compile(r"^\.step_\d{12}\.incomplete\.")
ARCHIVE_ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")
KNOWN_CHECKPOINT_VERSION = 1
#: Component assigned to any path the classifier does not recognize. Retained locally
#: by every offload policy, so an unregistered file is reported rather than removed.
UNCLASSIFIED_COMPONENT = "unclassified"

#: Workspace identity file. Storage reads it to name a workspace archive; it never
#: rewrites workspace configuration.
WORKSPACE_CONFIG_FILENAME = ".picurv-workspace.yml"

#: Directories a workspace archive never descends into. Runs and studies are their own
#: artifacts with their own archives; duplicating them inside a workspace archive would
#: store the same bytes twice and blur which object owns what.
WORKSPACE_EXCLUDED_ROOTS = ("runs", "studies")

#: Components no offload policy may prune, whatever it retains.
ALWAYS_RETAINED_COMPONENTS = frozenset(
    {UNCLASSIFIED_COMPONENT, "workspace-config", "workspace-inputs"}
)

STORAGE_RESTORE_COMPONENTS = (
    "inputs", "raw-output", "analysis", "visualization", "logs", "assets"
)
_PARALLEL_GZIP_VALUES = {"fast", "balanced"}


class StorageError(RuntimeError):
    """! @brief User-facing storage workflow failure. """


def _utc_now() -> str:
    """!
    @brief Return a stable UTC timestamp for storage metadata.
    @return Result produced by this operation.
    """
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _human_bytes(value: int) -> str:
    """!
    @brief Format a byte count for concise command output.
    @param[in] value Value supplied through the `value` argument.
    @return Result produced by this operation.
    """
    size = float(value)
    for suffix in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
        if size < 1024.0 or suffix == "PiB":
            return f"{size:.1f} {suffix}" if suffix != "B" else f"{int(size)} B"
        size /= 1024.0
    return f"{int(value)} B"


def _sha256_file(path: str) -> str:
    """!
    @brief Calculate SHA-256 without loading a potentially large file into memory.
    @param[in] path Value supplied through the `path` argument.
    @return Result produced by this operation.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while True:
            block = stream.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _atomic_write_json(path: str, payload: dict) -> None:
    """!
    @brief Atomically replace a JSON state or manifest file.
    @param[in] path Value supplied through the `path` argument.
    @param[in] payload Value supplied through the `payload` argument.
    """
    path_abs = os.path.abspath(path)
    os.makedirs(os.path.dirname(path_abs), exist_ok=True)
    temporary = f"{path_abs}.tmp.{os.getpid()}"
    with open(temporary, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path_abs)


def _read_json(path: str):
    """!
    @brief Read a JSON mapping when present, otherwise return None.
    @param[in] path Value supplied through the `path` argument.
    @return Result produced by this operation.
    """
    try:
        with open(path, "r", encoding="utf-8") as stream:
            payload = json.load(stream)
    except (OSError, ValueError):
        return None
    return payload if isinstance(payload, dict) else None


def _find_upwards(start: str, filename: str):
    """!
    @brief Find the nearest named file at or above a filesystem anchor.
    @param[in] start Value supplied through the `start` argument.
    @param[in] filename Value supplied through the `filename` argument.
    @return Result produced by this operation.
    """
    current = os.path.abspath(start)
    if os.path.isfile(current):
        current = os.path.dirname(current)
    while True:
        candidate = os.path.join(current, filename)
        if os.path.isfile(candidate):
            return candidate
        parent = os.path.dirname(current)
        if parent == current:
            return None
        current = parent


def resolve_storage_config_path(explicit_path: str = None, require: bool = True) -> str:
    """!
    @brief Resolve an explicit or nearest workspace storage configuration.
    @param[in] explicit_path Optional user-selected YAML path.
    @param[in] require Whether a missing configuration is an error.
    @return Result produced by this operation.
    """
    if explicit_path:
        result = os.path.abspath(explicit_path)
    else:
        result = _find_upwards(os.getcwd(), STORAGE_CONFIG_FILENAME)
    if result and os.path.isfile(result):
        return result
    if require:
        raise StorageError(
            "No PICurv storage configuration was found. Run "
            "'picurv storage setup --remote <rclone-remote:path>' first or pass --storage-config."
        )
    return os.path.abspath(explicit_path or os.path.join(os.getcwd(), STORAGE_CONFIG_FILENAME))


def load_storage_profile(profile_name: str = None, config_path: str = None) -> dict:
    """!
    @brief Load and validate one non-secret rclone storage profile.
    @param[in] profile_name Value supplied through the `profile_name` argument.
    @param[in] config_path Value supplied through the `config_path` argument.
    @return Result produced by this operation.
    """
    resolved_config = resolve_storage_config_path(config_path)
    with open(resolved_config, "r", encoding="utf-8") as stream:
        payload = yaml.safe_load(stream) or {}
    profiles = payload.get("profiles")
    if not isinstance(profiles, dict):
        raise StorageError(f"Storage config has no 'profiles' mapping: {resolved_config}")
    selected = profile_name or payload.get("default_profile") or DEFAULT_PROFILE_NAME
    profile = profiles.get(selected)
    if not isinstance(profile, dict):
        raise StorageError(f"Storage profile '{selected}' does not exist in {resolved_config}.")
    remote = profile.get("remote")
    if not isinstance(remote, str) or not remote.strip():
        raise StorageError(f"Storage profile '{selected}' requires a non-empty remote.")
    result = dict(profile)
    result["name"] = selected
    result["remote"] = remote.rstrip("/")
    result["config_path"] = resolved_config
    try:
        chunk_size_gib = float(result.get("chunk_size_gib", DEFAULT_CHUNK_SIZE_GIB))
    except (TypeError, ValueError) as exc:
        raise StorageError(f"Storage profile '{selected}' chunk_size_gib must be numeric.") from exc
    if chunk_size_gib <= 0.0:
        raise StorageError(f"Storage profile '{selected}' chunk_size_gib must be positive.")
    result["chunk_size_bytes"] = int(chunk_size_gib * 1024 ** 3)
    try:
        workers = int(result.get("workers", DEFAULT_STORAGE_WORKERS))
    except (TypeError, ValueError) as exc:
        raise StorageError(f"Storage profile '{selected}' workers must be an integer.") from exc
    if workers <= 0:
        raise StorageError(f"Storage profile '{selected}' workers must be positive.")
    result["workers"] = workers
    offload_policy = str(result.get("offload_policy", "metadata-only"))
    if offload_policy not in STORAGE_OFFLOAD_POLICIES:
        raise StorageError(
            f"Storage profile '{selected}' offload_policy must be one of: "
            + ", ".join(STORAGE_OFFLOAD_POLICIES)
        )
    result["offload_policy"] = offload_policy
    result["keep_latest_checkpoint"] = bool(result.get("keep_latest_checkpoint", False))
    return result


def _remote_join(remote: str, *parts: str) -> str:
    """!
    @brief Join path components without corrupting rclone remote syntax.
    @param[in] remote Value supplied through the `remote` argument.
    @param[in] parts Value supplied through the `parts` argument.
    @return Result produced by this operation.
    """
    clean_parts = [str(part).strip("/") for part in parts if str(part).strip("/")]
    suffix = "/".join(clean_parts)
    if not suffix:
        return remote
    if remote.endswith(":"):
        return remote + suffix
    return remote.rstrip("/") + "/" + suffix


def _object_remote(profile: dict, archive_id: str, *parts: str) -> str:
    """!
    @brief Return the remote path for one immutable archive object.
    @param[in] profile Value supplied through the `profile` argument.
    @param[in] archive_id Value supplied through the `archive_id` argument.
    @param[in] parts Value supplied through the `parts` argument.
    @return Result produced by this operation.
    """
    return _remote_join(profile["remote"], REMOTE_OBJECTS_DIRECTORY, archive_id, *parts)


def _blob_remote(profile: dict, digest: str) -> str:
    """!
    @brief Return the remote path of one content-addressed payload blob.
    @param[in] profile Active storage profile.
    @param[in] digest Payload SHA-256.
    @return Remote path, fanned out by digest prefix to keep directories small.
    """
    return _remote_join(profile["remote"], REMOTE_BLOBS_DIRECTORY, digest[:2], digest)


def _chunk_remote_path(profile: dict, archive_id: str, chunk: dict) -> str:
    """!
    @brief Resolve where one manifest chunk's payload lives on the remote.

    @details Schema 1 archives keep their payload beside the manifest; schema 2 and
             later resolve it out of the shared content-addressed store, so both remain
             restorable from the same catalog.
    @param[in] profile Active storage profile.
    @param[in] archive_id Owning archive id.
    @param[in] chunk Manifest chunk entry.
    @return Remote path of the chunk payload.
    """
    if chunk.get("blob"):
        return _blob_remote(profile, chunk["sha256"])
    return _object_remote(profile, archive_id, "chunks", chunk["name"])


def _remote_blob_present(profile: dict, digest: str) -> bool:
    """!
    @brief Whether a blob with this digest is already stored and intact.
    @param[in] profile Active storage profile.
    @param[in] digest Payload SHA-256.
    @return True when the remote already holds exactly this content.
    """
    try:
        return _remote_sha256(_blob_remote(profile, digest)) == digest
    except (StorageError, OSError):
        # An absent blob is the normal case, not an error; anything else that stops us
        # from confirming it is treated the same way, and the chunk is uploaded.
        return False


def _run_rclone(arguments: list, check: bool = True) -> subprocess.CompletedProcess:
    """!
    @brief Invoke rclone through the same argv-based subprocess boundary as other PICurv tools.
    @param[in] arguments Value supplied through the `arguments` argument.
    @param[in] check Value supplied through the `check` argument.
    @return Result produced by this operation.
    """
    executable = shutil.which("rclone")
    if not executable:
        raise StorageError("rclone was not found on PATH. Install/configure rclone before using PICurv storage.")
    result = subprocess.run(
        [executable] + [str(item) for item in arguments],
        text=True,
        capture_output=True,
        check=False,
    )
    if check and result.returncode != 0:
        detail = (result.stderr or result.stdout or "unknown rclone error").strip()
        raise StorageError(f"rclone {' '.join(str(item) for item in arguments[:2])} failed: {detail}")
    return result


def _remote_sha256(remote_path: str) -> str:
    """!
    @brief Ask rclone to calculate or retrieve the SHA-256 of one remote object.
    @param[in] remote_path Value supplied through the `remote_path` argument.
    @return Result produced by this operation.
    """
    result = _run_rclone(["hashsum", "SHA-256", remote_path])
    for line in result.stdout.splitlines():
        token = line.strip().split(None, 1)[0] if line.strip() else ""
        if re.fullmatch(r"[0-9a-fA-F]{64}", token):
            return token.lower()
    raise StorageError(f"rclone did not return a SHA-256 for {remote_path}.")


def _upload_verified(local_path: str, remote_path: str) -> dict:
    """!
    @brief Upload one file, then verify its remote SHA-256.
    @param[in] local_path Value supplied through the `local_path` argument.
    @param[in] remote_path Value supplied through the `remote_path` argument.
    @return Result produced by this operation.
    """
    local_digest = _sha256_file(local_path)
    _run_rclone(["copyto", local_path, remote_path])
    remote_digest = _remote_sha256(remote_path)
    if remote_digest != local_digest:
        raise StorageError(
            f"Remote checksum mismatch after upload: {remote_path} "
            f"(local {local_digest}, remote {remote_digest})."
        )
    return {"sha256": local_digest, "stored_bytes": os.path.getsize(local_path)}


def _read_remote_bytes(remote_path: str) -> bytes:
    """!
    @brief Read a small remote catalog object through rclone.
    @param[in] remote_path Value supplied through the `remote_path` argument.
    @return Result produced by this operation.
    """
    executable = shutil.which("rclone")
    if not executable:
        raise StorageError("rclone was not found on PATH.")
    result = subprocess.run(
        [executable, "cat", remote_path], capture_output=True, check=False
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or b"unknown rclone error").decode("utf-8", "replace").strip()
        raise StorageError(f"Unable to read remote object {remote_path}: {detail}")
    return result.stdout


def _load_remote_manifest(profile: dict, archive_id: str, require_complete: bool = True) -> dict:
    """!
    @brief Fetch and validate one versioned remote storage manifest.
    @param[in] profile Value supplied through the `profile` argument.
    @param[in] archive_id Value supplied through the `archive_id` argument.
    @param[in] require_complete Value supplied through the `require_complete` argument.
    @return Result produced by this operation.
    """
    if not ARCHIVE_ID_PATTERN.fullmatch(str(archive_id)):
        raise StorageError(f"Invalid archive ID: {archive_id!r}.")
    manifest_bytes = _read_remote_bytes(_object_remote(profile, archive_id, REMOTE_MANIFEST_FILENAME))
    if require_complete:
        complete = _read_remote_bytes(_object_remote(profile, archive_id, REMOTE_COMPLETE_FILENAME))
        recorded = complete.decode("ascii", "replace").strip().lower()
        actual = hashlib.sha256(manifest_bytes).hexdigest()
        if recorded != actual:
            raise StorageError(f"Archive {archive_id} has no valid completion marker.")
    try:
        manifest = json.loads(manifest_bytes.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise StorageError(f"Archive {archive_id} has an invalid manifest.") from exc
    if not isinstance(manifest, dict) or manifest.get("storage_schema_version") != STORAGE_SCHEMA_VERSION:
        raise StorageError(
            f"Archive {archive_id} uses unsupported storage schema "
            f"{manifest.get('storage_schema_version') if isinstance(manifest, dict) else 'unknown'}."
        )
    return manifest


def _state_path(root_path: str) -> str:
    """!
    @brief Return the local storage state marker path for an artifact root.
    @param[in] root_path Value supplied through the `root_path` argument.
    @return Result produced by this operation.
    """
    return os.path.join(os.path.abspath(root_path), STORAGE_STATE_FILENAME)


def read_storage_state(root_path: str):
    """!
    @brief Read the nearest applicable storage marker for a run, study, or study member.
    @param[in] root_path Value supplied through the `root_path` argument.
    @return Result produced by this operation.
    """
    root = Path(os.path.abspath(root_path))
    state = _read_json(_state_path(str(root)))
    if state:
        return state
    # A whole-study archive owns every numbered member.  Let ordinary run
    # workflows see that parent state without treating unrelated ancestors as
    # storage owners.
    if root.parent.name == "cases":
        return _read_json(_state_path(str(root.parent.parent)))
    return None


def is_artifact_cold(root_path: str) -> bool:
    """!
    @brief Return whether a local artifact marker says payload data was pruned.
    @param[in] root_path Value supplied through the `root_path` argument.
    @return Result produced by this operation.
    """
    state = read_storage_state(root_path)
    return bool(state and state.get("local_pruned"))


def workspace_asset_references(workspace_root: str) -> dict:
    """!
    @brief Count what still refers to each published workspace asset.

    @details A local copy may be removed once nothing local needs it and a verified
             remote copy exists. Runs that are themselves cold still *reference* the
             asset - that is what keeps the remote copy alive - but they do not keep
             the local one.
    @param[in] workspace_root Initialized workspace root.
    @return Mapping of asset id to its reference counts and local object path.
    """
    references = {}
    objects_root = os.path.join(workspace_root, "assets", "objects")
    if os.path.isdir(objects_root):
        for kind in sorted(os.listdir(objects_root)):
            kind_root = os.path.join(objects_root, kind)
            if not os.path.isdir(kind_root):
                continue
            for asset_id in sorted(os.listdir(kind_root)):
                object_root = os.path.join(kind_root, asset_id)
                if os.path.isdir(object_root):
                    references[asset_id] = {
                        "asset_id": asset_id, "kind": kind, "object": object_root,
                        "active_local_runs": 0, "cold_runs": 0,
                    }
    for artifacts in ("runs", "studies"):
        root = os.path.join(workspace_root, artifacts)
        if not os.path.isdir(root):
            continue
        for lock_path in Path(root).glob("**/inputs/assets.lock.yml"):
            try:
                with open(lock_path, "r", encoding="utf-8") as stream:
                    lock = yaml.safe_load(stream) or {}
            except (OSError, ValueError):
                continue
            run_root = lock_path.parent.parent
            cold = is_artifact_cold(str(run_root))
            for reference in (lock.get("assets") or {}).values():
                entry = references.get(reference.get("asset_id"))
                if entry is None:
                    continue
                entry["cold_runs" if cold else "active_local_runs"] += 1
    return references


def prune_unused_workspace_assets(workspace_root: str, profile: dict,
                                  dry_run: bool = False) -> list:
    """!
    @brief Remove local asset objects that nothing local needs and storage has verified.
    @param[in] workspace_root Initialized workspace root.
    @param[in] profile Resolved storage profile.
    @param[in] dry_run Report the decision without removing anything.
    @return Removal decisions, one per published asset.
    """
    protected = set()
    for manifest in list_remote_manifests(profile):
        if manifest.get("artifact_type") != "workspace":
            continue
        for asset_id in manifest.get("workspace_assets") or []:
            protected.add(asset_id)
    decisions = []
    for entry in workspace_asset_references(workspace_root).values():
        verified = entry["asset_id"] in protected
        removable = verified and entry["active_local_runs"] == 0
        decision = {**entry, "remote_protection": "verified" if verified else "none",
                    "local_removal": "safe" if removable else "blocked"}
        if removable and not dry_run:
            shutil.rmtree(entry["object"], ignore_errors=True)
            decision["removed"] = True
        decisions.append(decision)
    return decisions


def restore_cold_study_members(study_path: str, case_ids, profile_name: str = None,
                               storage_config: str = None) -> list:
    """!
    @brief Restore named cold study members in place from their local markers.
    @param[in] study_path Study directory owning the members.
    @param[in] case_ids Member ids to restore.
    @param[in] profile_name Optional storage profile name.
    @param[in] storage_config Optional explicit storage configuration path.
    @return Restored member ids.
    @throws StorageError when a member carries no usable archive reference.
    """
    profile = load_storage_profile(profile_name, storage_config)
    restored = []
    for case_id in case_ids:
        member = os.path.join(study_path, "cases", case_id)
        state = read_storage_state(member) or {}
        archive_id = state.get("archive_id")
        if not archive_id:
            raise StorageError(
                f"{case_id} is cold but its marker names no archive; restore it by id."
            )
        restore_archive(profile, archive_id, destination=member, force=True)
        restored.append(case_id)
    return restored


def cold_study_members(study_path: str) -> list:
    """!
    @brief Return numbered study members whose local payload was pruned.
    @param[in] study_path Value supplied through the `study_path` argument.
    @return Result produced by this operation.
    """
    cases_dir = Path(os.path.abspath(study_path)) / "cases"
    if not cases_dir.is_dir():
        return []
    return [
        child.name for child in sorted(cases_dir.iterdir())
        if child.is_dir() and is_artifact_cold(str(child))
    ]


def require_storage_payload_local(
    root_path: str,
    operation: str,
    checkpoint: int = None,
    checkpoints=None,
) -> None:
    """!
    @brief Reject a workflow that requires payload currently held in cold storage.
    @param[in] root_path Run or study-member directory checked by an existing workflow.
    @param[in] operation Human-readable consuming operation.
    @param[in] checkpoint Optional single required checkpoint step.
    @param[in] checkpoints Optional iterable of every required checkpoint step.
    """
    state = read_storage_state(root_path)
    if not state or not state.get("local_pruned"):
        return
    required_steps = set()
    if checkpoint is not None:
        required_steps.add(int(checkpoint))
    if checkpoints is not None:
        required_steps.update(int(step) for step in checkpoints)
    restored = set(state.get("restored_components") or [])
    restored.update(state.get("retained_components") or [])
    missing_steps = sorted(
        step for step in required_steps if f"checkpoint:{step}" not in restored
    )
    if required_steps and not missing_steps:
        return
    archive_id = state.get("archive_id", "<archive-id>")
    if missing_steps and len(missing_steps) <= 8:
        suffix = "".join(f" --checkpoint {step}" for step in missing_steps)
    else:
        # A full restore is clearer than printing hundreds of repeatable selectors.
        suffix = ""
    raise StorageError(
        f"{operation} requires payload archived from {os.path.abspath(root_path)}. Restore it first with:\n"
        f"  picurv storage restore --archive-id {archive_id}{suffix}"
    )


def storage_state_summary(root_path: str) -> dict:
    """!
    @brief Return a compact storage status for summarize and status commands.
    @param[in] root_path Value supplied through the `root_path` argument.
    @return Result produced by this operation.
    """
    state = read_storage_state(root_path)
    if not state:
        # A marker that exists but cannot be read is worse than none: the artifact
        # claims a remote copy nobody can locate.
        if os.path.exists(_state_path(root_path)):
            return {"state": "BROKEN", "archive_id": None, "label": None,
                    "detail": "storage marker is present but unreadable"}
        return {"state": "LOCAL", "archive_id": None, "label": None}
    if not state.get("archive_id"):
        return {"state": "BROKEN", "archive_id": None, "label": state.get("label"),
                "detail": "storage marker names no archive"}
    if state.get("local_pruned"):
        status = "PARTIAL" if state.get("restored_components") else "COLD"
    else:
        status = "PROTECTED"
    return {
        "state": status,
        "archive_id": state.get("archive_id"),
        "label": state.get("label"),
        "profile": state.get("profile"),
        "remote": state.get("remote"),
    }


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


def _lock_owner_active(metadata_path: str) -> bool:
    """!
    @brief Conservatively determine whether a solver/post/storage owner marker is active.
    @param[in] metadata_path Value supplied through the `metadata_path` argument.
    @return Result produced by this operation.
    """
    owner = _read_json(metadata_path)
    if not owner:
        return True
    host = owner.get("host")
    pid = owner.get("pid")
    if host and host != socket.gethostname():
        return True
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except (PermissionError, OSError):
        return True
    return True


def _collect_job_ids(payload) -> set:
    """!
    @brief Recursively collect submitted Slurm job IDs from scheduler metadata.
    @param[in] payload Value supplied through the `payload` argument.
    @return Result produced by this operation.
    """
    result = set()
    if isinstance(payload, dict):
        if payload.get("submitted") and payload.get("job_id") is not None:
            result.add(str(payload["job_id"]).strip())
        for value in payload.values():
            result.update(_collect_job_ids(value))
    elif isinstance(payload, list):
        for value in payload:
            result.update(_collect_job_ids(value))
    return {item for item in result if item}


def _slurm_activity(root: str) -> dict:
    """!
    @brief Query live Slurm state for every job recorded below an artifact scheduler directory.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    job_ids = set()
    scheduler_dirs = [Path(root) / "scheduler"]
    if (Path(root) / "cases").is_dir():
        scheduler_dirs.extend((Path(root) / "cases").glob("*/scheduler"))
    for scheduler in scheduler_dirs:
        if not scheduler.is_dir():
            continue
        for path in scheduler.glob("submission*.json"):
            job_ids.update(_collect_job_ids(_read_json(str(path))))
    if not job_ids:
        return {"job_ids": [], "active": [], "unknown": False}
    squeue = shutil.which("squeue")
    if not squeue:
        return {"job_ids": sorted(job_ids), "active": [], "unknown": True}
    active = []
    unknown_ids = []
    for recorded_id in sorted(job_ids):
        result = subprocess.run(
            [squeue, "-h", "-j", recorded_id, "-o", "%i|%T"],
            text=True, capture_output=True, check=False,
        )
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").lower()
            # Slurm returns nonzero for a valid historical/purged job id. Such an id
            # cannot be queued, so it is inactive rather than an ambiguous failure.
            if "invalid job id" in detail or "invalid job/step" in detail:
                continue
            unknown_ids.append(recorded_id)
            continue
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            job_id, _, state = line.partition("|")
            active.append({"job_id": job_id.strip(), "state": state.strip() or "UNKNOWN"})
    return {
        "job_ids": sorted(job_ids),
        "active": active,
        "unknown": bool(unknown_ids),
        "unknown_job_ids": unknown_ids,
    }


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


def _assert_archive_safe(inventory: dict) -> None:
    """!
    @brief Refuse to package a changing or scheduler-ambiguous artifact.
    @param[in] inventory Value supplied through the `inventory` argument.
    """
    problems = []
    if inventory["incomplete_checkpoints"]:
        problems.append("incomplete checkpoint(s): " + ", ".join(inventory["incomplete_checkpoints"][:3]))
    if inventory["active_locks"]:
        problems.append("active runtime lock(s): " + ", ".join(inventory["active_locks"]))
    if inventory["slurm"]["active"]:
        problems.append(
            "active Slurm job(s): " + ", ".join(
                f"{item['job_id']} ({item['state']})" for item in inventory["slurm"]["active"]
            )
        )
    if inventory["slurm"]["unknown"]:
        problems.append(
            "recorded Slurm job IDs could not be checked because squeue is unavailable or failed"
        )
    if problems:
        raise StorageError("Artifact is not safe to archive/offload: " + "; ".join(problems) + ".")


@contextlib.contextmanager
def storage_operation_lock(root_path: str, operation: str):
    """!
    @brief Hold an exclusive local storage-operation marker for one artifact.
    @param[in] root_path Value supplied through the `root_path` argument.
    @param[in] operation Value supplied through the `operation` argument.
    """
    root = os.path.abspath(root_path)
    lock_path = os.path.join(root, STORAGE_LOCK_FILENAME)
    if os.path.exists(lock_path):
        if _lock_owner_active(lock_path):
            raise StorageError(f"Another storage operation owns {lock_path}.")
        try:
            os.remove(lock_path)
        except OSError as exc:
            raise StorageError(f"Unable to remove stale storage lock {lock_path}: {exc}") from exc
    payload = {"operation": operation, "pid": os.getpid(), "host": socket.gethostname(), "started_at": _utc_now()}
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    try:
        descriptor = os.open(lock_path, flags, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        yield
    finally:
        try:
            os.remove(lock_path)
        except FileNotFoundError:
            pass


@contextlib.contextmanager
def runtime_stage_lock(root_path: str, stage: str):
    """!
    @brief Mark a locally executed solver/post stage as active for storage safety.
    @param[in] root_path Run directory used as the runtime working directory.
    @param[in] stage Runtime stage label; storage currently uses this for solver execution.
    """
    scheduler = os.path.join(os.path.abspath(root_path), "scheduler")
    os.makedirs(scheduler, exist_ok=True)
    lock_path = os.path.join(scheduler, f"{stage}.lock.json")
    if os.path.exists(lock_path):
        if _lock_owner_active(lock_path):
            raise StorageError(f"A {stage} runtime stage already owns {lock_path}.")
        os.remove(lock_path)
    payload = {"stage": stage, "pid": os.getpid(), "host": socket.gethostname(), "started_at": _utc_now()}
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(lock_path, flags, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    try:
        yield
    finally:
        try:
            os.remove(lock_path)
        except FileNotFoundError:
            pass


#: Output encodings shared with the rest of the CLI. Storage does not own this set.
CLI_OUTPUT_FORMATS = ("text", "json")

#: Archive compression policies `picurv storage` accepts. `auto` resolves to one of
#: the concrete policies from the payload size.
STORAGE_COMPRESSION_POLICIES = ("auto", "none", "fast", "balanced", "maximum")
STORAGE_OFFLOAD_POLICIES = ("metadata-only", "restart-ready", "analysis-ready")

#: Archive suffix each resolved compression policy produces. `auto` never appears
#: here: it resolves to one of the concrete policies before an extension is needed.
STORAGE_COMPRESSION_EXTENSIONS = {
    "none": ".tar", "fast": ".tar.gz", "balanced": ".tar.gz", "maximum": ".tar.xz",
}


def _select_compression(requested: str, total_bytes: int, profile: dict) -> str:
    """!
    @brief Resolve automatic or configured compression policy.
    @param[in] requested Value supplied through the `requested` argument.
    @param[in] total_bytes Value supplied through the `total_bytes` argument.
    @param[in] profile Value supplied through the `profile` argument.
    @return Result produced by this operation.
    """
    selected = requested or profile.get("compression", "auto")
    selected = str(selected).strip().lower()
    if selected not in STORAGE_COMPRESSION_POLICIES:
        raise StorageError("Compression must be one of: auto, none, fast, balanced, maximum.")
    if selected != "auto":
        return selected
    if total_bytes < AUTO_NO_COMPRESSION_BYTES:
        return "none"
    if total_bytes >= AUTO_MAXIMUM_COMPRESSION_BYTES:
        return "maximum"
    return "balanced"


def _chunk_extension(compression: str) -> str:
    """!
    @brief Return the archive suffix for a compression policy.
    @param[in] compression Value supplied through the `compression` argument.
    @return Result produced by this operation.
    """
    return STORAGE_COMPRESSION_EXTENSIONS[compression]


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


def _find_reusable_archive(profile: dict, target: dict, fingerprint: str) -> dict:
    """!
    @brief Find a completed archive of this artifact whose content is already current.
    @param[in] profile Resolved storage profile.
    @param[in] target Local artifact target.
    @param[in] fingerprint Inventory fingerprint of the artifact as it stands now.
    @return Matching remote manifest, or None.
    """
    if not fingerprint:
        return None
    identity = (target["artifact_type"], target.get("run_id"),
                target.get("study_id"), target.get("case_id"))
    try:
        candidates = list_remote_manifests(profile)
    except StorageError as exc:
        # Reuse is an optimization. A remote that cannot be listed - not yet created,
        # briefly unreachable - must fall through to a normal upload, never fail here.
        print(f"[INFO] Could not check for a reusable archive ({exc}); uploading.")
        return None
    newest = None
    for manifest in candidates:
        if manifest.get("inventory_sha256") != fingerprint:
            continue
        if (manifest.get("artifact_type"), manifest.get("run_id"),
                manifest.get("study_id"), manifest.get("case_id")) != identity:
            continue
        if newest is None or str(manifest.get("created_at", "")) > str(newest.get("created_at", "")):
            newest = manifest
    return newest


def _build_chunk_specs(inventory: dict, chunk_size_bytes: int) -> list:
    """!
    @brief Group archive entries into independently transferable component chunks.
    @param[in] inventory Value supplied through the `inventory` argument.
    @param[in] chunk_size_bytes Value supplied through the `chunk_size_bytes` argument.
    @return Result produced by this operation.
    """
    groups = {}
    directory_groups = {}
    for entry in inventory["entries"]:
        if entry["type"] == "directory":
            # A directory belongs to the same component as the files inside it. Sweeping
            # every directory into `metadata` instead meant restoring one checkpoint
            # recreated the empty directory tree of every other checkpoint, so a
            # partially restored run listed steps whose payload was not there.
            directory_groups.setdefault(entry["component"], []).append(entry["path"])
            continue
        groups.setdefault(entry["component"], []).append(entry)
    specs = []
    component_order = sorted(
        set(groups) | set(directory_groups),
        key=lambda name: (not name.startswith("checkpoint:"), name),
    )
    for component in component_order:
        # Directories ride in the first chunk of their own component, so extracting
        # that component creates its tree and no one else's.
        pending_directories = sorted(directory_groups.get(component, []))
        if component not in groups:
            if pending_directories:
                specs.append({
                    "component": component,
                    "entries": pending_directories,
                    "uncompressed_bytes": 0,
                })
            continue
        current = list(pending_directories)
        current_bytes = 0
        for entry in groups[component]:
            entry_size = max(1, int(entry["size"]))
            if current and current_bytes + entry_size > chunk_size_bytes:
                specs.append({"component": component, "entries": current, "uncompressed_bytes": current_bytes})
                current = []
                current_bytes = 0
            current.append(entry["path"])
            current_bytes += int(entry["size"])
        if current:
            specs.append({"component": component, "entries": current, "uncompressed_bytes": current_bytes})
    return specs


def _safe_component_name(component: str) -> str:
    """!
    @brief Convert a component token into a portable archive filename fragment.
    @param[in] component Value supplied through the `component` argument.
    @return Result produced by this operation.
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", component).strip("-") or "data"


def _write_tar_chunk(root: str, spec: dict, destination: str, compression: str,
                     workers: int = 1) -> str:
    """!
    @brief Package explicitly inventoried entries without following symlinks.
    @param[in] root Value supplied through the `root` argument.
    @param[in] spec Value supplied through the `spec` argument.
    @param[in] destination Value supplied through the `destination` argument.
    @param[in] compression Value supplied through the `compression` argument.
    @param[in] workers Maximum compressor worker count.
    @return Compressor implementation used for the archive chunk.
    """
    tar_executable = shutil.which("tar")
    compressor = None
    compressor_args = []
    if compression in _PARALLEL_GZIP_VALUES and shutil.which("pigz"):
        compressor = shutil.which("pigz")
        compressor_args = [compressor, "-c", "-p", str(workers), "-1" if compression == "fast" else "-6"]
    elif compression == "maximum" and shutil.which("xz"):
        compressor = shutil.which("xz")
        compressor_args = [compressor, "-c", "-T", str(workers), "-9e"]
    if compressor and tar_executable:
        list_file = tempfile.NamedTemporaryFile(prefix="picurv-tar-list-", delete=False)
        try:
            list_file.write(b"\0".join(path.encode("utf-8") for path in spec["entries"]) + b"\0")
            list_file.close()
            tar_command = [
                tar_executable, "-C", root, "--null", "--verbatim-files-from", "--no-recursion",
                "-T", list_file.name, "-cf", "-",
            ]
            with tempfile.TemporaryFile() as tar_stderr, open(destination, "wb") as output:
                producer = subprocess.Popen(tar_command, stdout=subprocess.PIPE, stderr=tar_stderr)
                consumer = subprocess.Popen(
                    compressor_args, stdin=producer.stdout, stdout=output, stderr=subprocess.PIPE
                )
                producer.stdout.close()
                _, compressor_stderr = consumer.communicate()
                producer_code = producer.wait()
                if producer_code != 0 or consumer.returncode != 0:
                    tar_stderr.seek(0)
                    detail = (
                        tar_stderr.read().decode("utf-8", "replace")
                        or compressor_stderr.decode("utf-8", "replace")
                        or "parallel archive command failed"
                    ).strip()
                    raise StorageError(detail)
            return f"{os.path.basename(compressor)}:{workers}"
        finally:
            try:
                os.remove(list_file.name)
            except FileNotFoundError:
                pass

    kwargs = {}
    if compression == "none":
        mode = "w"
    elif compression == "fast":
        mode = "w:gz"
        kwargs["compresslevel"] = 1
    elif compression == "balanced":
        mode = "w:gz"
        kwargs["compresslevel"] = 6
    else:
        mode = "w:xz"
        kwargs["preset"] = 9
    with tarfile.open(destination, mode, dereference=False, **kwargs) as archive:
        for relative in spec["entries"]:
            source = os.path.join(root, *relative.split("/"))
            if not os.path.lexists(source):
                raise StorageError(f"Artifact changed during packaging; entry disappeared: {source}")
            archive.add(source, arcname=relative, recursive=False)
    return "python-tarfile"


def _capture_study_context(target: dict) -> list:
    """!
    @brief Embed small study control-plane files with an individually archived member.
    @param[in] target Value supplied through the `target` argument.
    @return Result produced by this operation.
    """
    if target["artifact_type"] != "study-case":
        return []
    study_root = Path(target["study_path"])
    candidates = [
        study_root / "study.yml",
        study_root / "cluster.yml",
        study_root / "study_manifest.json",
        study_root / "scheduler" / "case_index.tsv",
        study_root / "scheduler" / "submission.json",
    ]
    captured = []
    for path in candidates:
        if not path.is_file() or path.stat().st_size > 5 * 1024 * 1024:
            continue
        captured.append({
            "path": path.relative_to(study_root).as_posix(),
            "mode": int(path.stat().st_mode & 0o7777),
            "content_base64": base64.b64encode(path.read_bytes()).decode("ascii"),
        })
    return captured


def _capture_parameter_summary(root: str) -> dict:
    """!
    @brief Capture the few case values that identify what a run actually solved.

    @details `storage show` has to answer "which run was this?" months later without
             restoring anything, and an id plus a label does not. These are read from
             the run's own configuration snapshot, so they describe the run rather than
             whatever the editable workspace says now.
    @param[in] root Artifact root directory.
    @return Parameter summary, empty when no readable case snapshot exists.
    """
    for relative in (("config", "case.yml"), ("cases", "case_0000", "config", "case.yml")):
        case_path = os.path.join(root, *relative)
        if os.path.isfile(case_path):
            break
    else:
        return {}
    try:
        with open(case_path, "r", encoding="utf-8") as stream:
            case = yaml.safe_load(stream) or {}
    except (OSError, ValueError):
        return {}
    if not isinstance(case, dict):
        return {}
    properties = case.get("properties") or {}
    scaling = properties.get("scaling") or {}
    fluid = properties.get("fluid") or {}
    grid = case.get("grid") or {}
    run_control = case.get("run_control") or {}
    summary = {
        "title": case.get("title"),
        "grid_mode": grid.get("mode"),
        "grid_dimensions": None,
        "blocks": ((case.get("models") or {}).get("domain") or {}).get("blocks"),
        "start_step": run_control.get("start_step"),
        "total_steps": run_control.get("total_steps"),
        "dt_physical": run_control.get("dt_physical"),
        "length_ref": scaling.get("length_ref"),
        "velocity_ref": scaling.get("velocity_ref"),
        "density": fluid.get("density"),
        "viscosity": fluid.get("viscosity"),
        "reynolds": None,
    }
    programmatic = grid.get("programmatic_settings")
    if isinstance(programmatic, dict):
        dims = [programmatic.get(key) for key in ("im", "jm", "km")]
        if all(isinstance(value, int) for value in dims):
            summary["grid_dimensions"] = dims
    try:
        summary["reynolds"] = (
            float(fluid["density"]) * float(scaling["velocity_ref"])
            * float(scaling["length_ref"]) / float(fluid["viscosity"])
        )
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        summary["reynolds"] = None
    return {key: value for key, value in summary.items() if value is not None}


def _capture_workspace_assets(target: dict) -> list:
    """!
    @brief List the asset ids a workspace archive carries.
    @param[in] target Local artifact target.
    @return Sorted asset ids, empty for non-workspace artifacts.
    """
    if target["artifact_type"] != "workspace":
        return []
    objects_root = os.path.join(target["root_path"], "assets", "objects")
    if not os.path.isdir(objects_root):
        return []
    found = []
    for kind in sorted(os.listdir(objects_root)):
        kind_root = os.path.join(objects_root, kind)
        if os.path.isdir(kind_root):
            found.extend(
                name for name in sorted(os.listdir(kind_root))
                if os.path.isdir(os.path.join(kind_root, name))
            )
    return sorted(found)


def _capture_run_assets(root: str) -> list:
    """!
    @brief Record reusable asset references carried by run-local input snapshots.
    @param[in] root Run or study artifact root.
    @return Portable asset references for the remote archive catalog.
    """
    captured = []
    archive_root = Path(os.path.abspath(root))
    for runtime_root_text in _artifact_runtime_roots(root):
        runtime_root = Path(runtime_root_text)
        lock_path = runtime_root / "inputs" / "assets.lock.yml"
        if not lock_path.is_file():
            continue
        try:
            payload = yaml.safe_load(lock_path.read_text(encoding="utf-8")) or {}
        except (OSError, ValueError, TypeError):
            continue
        assets = payload.get("assets") if isinstance(payload, dict) else None
        if not isinstance(assets, dict):
            continue
        runtime_prefix = runtime_root.relative_to(archive_root).as_posix()
        for kind, reference in sorted(assets.items()):
            if not isinstance(reference, dict):
                continue
            files = reference.get("files") or reference.get("exposed") or []
            normalized_files = []
            for item in files:
                if not isinstance(item, dict) or not isinstance(item.get("path"), str):
                    continue
                run_path = item["path"]
                archived_path = (
                    f"{runtime_prefix}/{run_path}" if runtime_prefix != "." else run_path
                )
                normalized_files.append({
                    "path": run_path,
                    "archived_path": archived_path,
                    "bytes": item.get("bytes"),
                    "sha256": item.get("sha256"),
                })
            captured.append({
                "kind": kind,
                "asset_id": reference.get("asset_id"),
                "provider": reference.get("provider"),
                "provider_spec_sha256": reference.get("provider_spec_sha256"),
                "files": normalized_files,
            })
    return captured


def _config_fingerprints(root: str) -> dict:
    """!
    @brief Hash small canonical YAML inputs stored under an artifact.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    result = {}
    root_path = Path(root)
    candidates = set(root_path.glob("config/*.yml"))
    candidates.update(root_path.glob("base_configs/*.yml"))
    candidates.update(root_path.glob("cases/case_*/config/*.yml"))
    for path in sorted(candidates):
        result[path.relative_to(root).as_posix()] = _sha256_file(str(path))
    for name in ("study.yml", "cluster.yml"):
        path = root_path / name
        if path.is_file():
            result[name] = _sha256_file(str(path))
    return result


def _git_provenance(root: str) -> dict:
    """!
    @brief Record best-effort current source revision and dirty state.
    @param[in] root Value supplied through the `root` argument.
    @return Result produced by this operation.
    """
    result = {"commit": None, "dirty": None}
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, capture_output=True, check=False
        )
        status = subprocess.run(
            ["git", "status", "--porcelain"], cwd=root, text=True, capture_output=True, check=False
        )
        if commit.returncode == 0:
            result["commit"] = commit.stdout.strip()
        if status.returncode == 0:
            result["dirty"] = bool(status.stdout.strip())
    except OSError:
        pass
    return result


def _parse_tags(raw_tags) -> dict:
    """!
    @brief Parse repeatable KEY=VALUE tags into deterministic metadata.
    @param[in] raw_tags Value supplied through the `raw_tags` argument.
    @return Result produced by this operation.
    """
    tags = {}
    for item in raw_tags or []:
        key, separator, value = str(item).partition("=")
        if not separator or not key.strip() or not value.strip():
            raise StorageError(f"Invalid tag {item!r}; expected KEY=VALUE.")
        tags[key.strip()] = value.strip()
    return tags


def _resolve_offload_policy(profile: dict, requested: str = None,
                            keep_latest_checkpoint=None) -> dict:
    """!
    @brief Resolve semantic local-retention behavior for an offload.
    @param[in] profile Active storage profile.
    @param[in] requested Optional command-level policy override.
    @param[in] keep_latest_checkpoint Optional checkpoint-retention override.
    @return Normalized retention policy mapping.
    """
    name = requested or profile.get("offload_policy", "metadata-only")
    if name not in STORAGE_OFFLOAD_POLICIES:
        raise StorageError("Offload policy must be one of: " + ", ".join(STORAGE_OFFLOAD_POLICIES))
    if keep_latest_checkpoint is None:
        keep_latest_checkpoint = bool(profile.get("keep_latest_checkpoint", False))
    retained = {
        "metadata-only": {"metadata", "logs"},
        "restart-ready": {"metadata", "logs", "inputs"},
        "analysis-ready": {"metadata", "logs", "analysis", "visualization"},
    }[name]
    return {
        "name": name,
        "retained_components": sorted(retained),
        "keep_latest_checkpoint": bool(keep_latest_checkpoint or name == "restart-ready"),
    }


def _entry_retained_by_policy(entry: dict, policy: dict, latest_step) -> bool:
    """!
    @brief Return whether one inventoried entry remains local after offload.
    @param[in] entry Inventoried artifact entry.
    @param[in] policy Normalized offload policy.
    @param[in] latest_step Newest committed checkpoint step, if any.
    @return True when the entry is retained locally.
    """
    component = entry["component"]
    # Never pruned, whatever the policy: a file storage cannot classify, and a
    # workspace's editable configuration and user-supplied inputs. Protecting a
    # workspace is a backup, not a handover of the user's own files.
    if component in ALWAYS_RETAINED_COMPONENTS:
        return True
    if component in set(policy["retained_components"]):
        return True
    return bool(
        policy["keep_latest_checkpoint"]
        and latest_step is not None
        and component == f"checkpoint:{latest_step}"
    )


def _compression_size_range(source_bytes: int, compression: str) -> tuple:
    """!
    @brief Return a deliberately broad planning estimate, never a promised ratio.
    @param[in] source_bytes Uncompressed payload byte count.
    @param[in] compression Selected compression policy.
    @return Estimated low and high archive sizes in bytes.
    """
    ratios = {
        "none": (1.0, 1.03),
        "fast": (0.45, 0.90),
        "balanced": (0.35, 0.85),
        "maximum": (0.25, 0.80),
    }
    low, high = ratios[compression]
    return int(source_bytes * low), int(source_bytes * high)


def build_storage_plan(target: dict, profile: dict, compression: str = None,
                       policy: str = None, keep_latest_checkpoint=None,
                       workers: int = None) -> dict:
    """!
    @brief Build the read-only plan consumed by protect and offload.
    @param[in] target Value supplied through the `target` argument.
    @param[in] profile Value supplied through the `profile` argument.
    @param[in] compression Value supplied through the `compression` argument.
    @param[in] policy Optional semantic retention policy override.
    @param[in] keep_latest_checkpoint Optional newest-checkpoint retention override.
    @param[in] workers Optional compression worker count override.
    @return Result produced by this operation.
    """
    inventory = inspect_artifact(target)
    selected_compression = _select_compression(compression, inventory["total_bytes"], profile)
    specs = _build_chunk_specs(inventory, profile["chunk_size_bytes"])
    retention = _resolve_offload_policy(profile, policy, keep_latest_checkpoint)
    latest_step = max(inventory["checkpoint_steps"], default=None)
    retained_bytes = sum(
        entry["size"] for entry in inventory["entries"]
        if entry["type"] != "directory" and _entry_retained_by_policy(entry, retention, latest_step)
    )
    worker_count = int(workers or profile.get("workers", DEFAULT_STORAGE_WORKERS))
    if worker_count <= 0:
        raise StorageError("Storage workers must be positive.")
    compressed_range = _compression_size_range(inventory["total_bytes"], selected_compression)
    return {
        "inventory": inventory,
        "compression": selected_compression,
        "chunk_count": len(specs),
        "workers": worker_count,
        "offload_policy": retention,
        "retained_local_bytes": retained_bytes,
        "pruned_local_bytes": max(0, inventory["total_bytes"] - retained_bytes),
        "estimated_stored_bytes": {"low": compressed_range[0], "high": compressed_range[1]},
        "chunks": [
            {
                "component": spec["component"],
                "file_count": len(spec["entries"]),
                "uncompressed_bytes": spec["uncompressed_bytes"],
            }
            for spec in specs
        ],
    }


def _render_plan(plan: dict) -> None:
    """!
    @brief Print a concise archive/offload plan.
    @param[in] plan Value supplied through the `plan` argument.
    """
    inventory = plan["inventory"]
    target = inventory["target"]
    print(f"[INFO] Artifact type : {target['artifact_type']}")
    print(f"[INFO] Artifact path : {target['root_path']}")
    print(f"[INFO] Local size    : {_human_bytes(inventory['total_bytes'])}")
    print(f"[INFO] Files         : {inventory['file_count']}")
    print(f"[INFO] Checkpoints   : {len(inventory['checkpoint_steps'])}")
    print(f"[INFO] Compression   : {plan['compression']}")
    print(f"[INFO] CPU workers   : {plan['workers']}")
    print(f"[INFO] Archive chunks: {plan['chunk_count']}")
    estimate = plan["estimated_stored_bytes"]
    print(
        f"[INFO] Remote estimate: {_human_bytes(estimate['low'])}–"
        f"{_human_bytes(estimate['high'])} (content-dependent)"
    )
    policy = plan["offload_policy"]
    print(f"[INFO] Offload policy: {policy['name']}")
    print(f"[INFO] Retained local: {_human_bytes(plan['retained_local_bytes'])}")
    print(f"[INFO] Pruned local  : {_human_bytes(plan['pruned_local_bytes'])}")
    if policy["keep_latest_checkpoint"]:
        latest = max(inventory["checkpoint_steps"], default=None)
        print(f"[INFO] Kept checkpoint: {latest if latest is not None else 'none available'}")
    if inventory["external_paths"]:
        print("[WARNING] External configured paths are recorded but are not followed automatically:")
        for item in inventory["external_paths"]:
            print(f"  - {item['source']}: {item['path']}")
    if inventory["dependencies"]:
        print("[WARNING] External run dependencies:")
        for item in inventory["dependencies"]:
            print(f"  - {item['kind']}: {item['path']}")
    unclassified = [
        entry for entry in inventory["entries"]
        if entry["type"] != "directory" and entry["component"] == UNCLASSIFIED_COMPONENT
    ]
    if unclassified:
        print(
            f"[WARNING] {len(unclassified)} file(s) are not part of any known component. "
            "They are archived, and retained locally by every policy, because storage "
            "does not delete files whose purpose it cannot state:"
        )
        for entry in unclassified[:10]:
            print(f"  - {entry['path']}")
        if len(unclassified) > 10:
            print(f"  - ... and {len(unclassified) - 10} more")


def archive_artifact(target: dict, profile: dict, label: str = None, tags=None,
                     compression: str = None, prune_local: bool = False,
                     policy: str = None, keep_latest_checkpoint=None,
                     workers: int = None, notes: str = None) -> dict:
    """!
    @brief Package, upload, verify, register, and optionally prune one artifact.
    @param[in] target Value supplied through the `target` argument.
    @param[in] profile Value supplied through the `profile` argument.
    @param[in] label Value supplied through the `label` argument.
    @param[in] tags Value supplied through the `tags` argument.
    @param[in] compression Value supplied through the `compression` argument.
    @param[in] prune_local Value supplied through the `prune_local` argument.
    @param[in] policy Optional semantic retention policy override.
    @param[in] keep_latest_checkpoint Optional newest-checkpoint retention override.
    @param[in] workers Optional compression worker count override.
    @param[in] notes Optional free-text note recorded with the archive.
    @return Result produced by this operation.
    """
    with storage_operation_lock(target["root_path"], "offload" if prune_local else "protect"):
        plan = build_storage_plan(
            target, profile, compression, policy=policy,
            keep_latest_checkpoint=keep_latest_checkpoint, workers=workers,
        )
        inventory = plan["inventory"]
        _assert_archive_safe(inventory)
        fingerprint = _inventory_fingerprint(inventory)
        reusable = _find_reusable_archive(profile, target, fingerprint)
        if reusable is not None:
            # `protect` then `offload` is the documented workflow for "back it up now,
            # free the space later". Re-packaging and re-uploading an unchanged artifact
            # would store a second full copy of it, which for a campaign-sized run is
            # the most expensive thing this command can do for no benefit.
            print(
                f"[INFO] Reusing verified archive {reusable['archive_id']}: "
                f"{target['root_path']} is unchanged since it was archived."
            )
            return _register_existing_archive(
                target, profile, reusable, inventory, plan, prune_local=prune_local,
                label=label, tags=tags, notes=notes,
            )
        archive_id = uuid.uuid4().hex
        specs = _build_chunk_specs(inventory, profile["chunk_size_bytes"])
        staging_parent = profile.get("staging_directory")
        if staging_parent:
            staging_parent = os.path.abspath(os.path.expanduser(str(staging_parent)))
            os.makedirs(staging_parent, exist_ok=True)
        print(f"[INFO] Creating archive {archive_id} for {target['root_path']}")
        with tempfile.TemporaryDirectory(prefix="picurv-storage-", dir=staging_parent) as staging:
            def package_and_upload(index_spec):
                """!
                @brief Package and verify one component chunk.
                @param[in] index_spec Tuple of chunk index and inventory specification.
                @return Uploaded chunk manifest entry.
                """
                index, spec = index_spec
                component = _safe_component_name(spec["component"])
                filename = f"{index:05d}_{component}{_chunk_extension(plan['compression'])}"
                local_chunk = os.path.join(staging, filename)
                print(
                    f"[INFO] Packaging chunk {index + 1}/{len(specs)}: "
                    f"{spec['component']} ({_human_bytes(spec['uncompressed_bytes'])})"
                )
                compressor = _write_tar_chunk(
                    target["root_path"], spec, local_chunk, plan["compression"],
                    workers=plan["workers"],
                )
                try:
                    digest = _sha256_file(local_chunk)
                    # Content-addressed: an identical chunk anywhere on this remote is
                    # already stored, so a rerun after a failed upload skips what
                    # succeeded, and two archives of the same checkpoint share one copy.
                    if _remote_blob_present(profile, digest):
                        print(f"[INFO] Chunk {index + 1}/{len(specs)} already stored; skipping upload.")
                        verified = {
                            "sha256": digest, "stored_bytes": os.path.getsize(local_chunk),
                        }
                    else:
                        verified = _upload_verified(local_chunk, _blob_remote(profile, digest))
                    return {
                        "name": filename,
                        "component": spec["component"],
                        "file_count": len(spec["entries"]),
                        "uncompressed_bytes": spec["uncompressed_bytes"],
                        "stored_bytes": verified["stored_bytes"],
                        "sha256": verified["sha256"],
                        "compressor": compressor,
                        "blob": True,
                    }
                finally:
                    try:
                        os.remove(local_chunk)
                    except FileNotFoundError:
                        pass

            # pigz/xz already use every requested CPU on one large chunk. The
            # Python fallback instead gains parallelism across independent chunks.
            native_parallel = (
                plan["compression"] in _PARALLEL_GZIP_VALUES and shutil.which("pigz")
            ) or (plan["compression"] == "maximum" and shutil.which("xz"))
            task_workers = 1 if native_parallel else min(plan["workers"], max(1, len(specs)))
            if task_workers == 1:
                chunks = [package_and_upload(item) for item in enumerate(specs)]
            else:
                with concurrent.futures.ThreadPoolExecutor(max_workers=task_workers) as executor:
                    chunks = list(executor.map(package_and_upload, enumerate(specs)))
            chunks.sort(key=lambda item: item["name"])

            manifest = {
                "storage_schema_version": STORAGE_SCHEMA_VERSION,
                "archive_id": archive_id,
                "created_at": _utc_now(),
                "artifact_type": target["artifact_type"],
                "run_id": target.get("run_id"),
                "study_id": target.get("study_id"),
                "case_id": target.get("case_id"),
                "workspace_id": target.get("workspace_id"),
                "label": label or os.path.basename(target["root_path"]),
                "notes": notes,
                "tags": _parse_tags(tags),
                # What was actually solved, so `storage show` answers "which run was
                # this?" without restoring anything.
                "parameters": _capture_parameter_summary(target["root_path"]),
                "original_path": target["original_path"],
                "original_study_path": target.get("study_path"),
                "profile": profile["name"],
                "remote": profile["remote"],
                "source_bytes": inventory["total_bytes"],
                "source_file_count": inventory["file_count"],
                "inventory_sha256": fingerprint,
                "compression": plan["compression"],
                "workers": plan["workers"],
                "offload_policy": plan["offload_policy"],
                "checkpoint_format_version": KNOWN_CHECKPOINT_VERSION,
                "checkpoint_steps": inventory["checkpoint_steps"],
                "chunks": chunks,
                "config_sha256": _config_fingerprints(target["root_path"]),
                "git": _git_provenance(target["root_path"]),
                "external_paths": inventory["external_paths"],
                "dependencies": inventory["dependencies"],
                "study_context": _capture_study_context(target),
                "run_assets": _capture_run_assets(target["root_path"]),
                "workspace_assets": _capture_workspace_assets(target),
                "capabilities": {
                    "restorable": True,
                    "continuable": bool(inventory["checkpoint_steps"]),
                    "reprocessable": bool(inventory["checkpoint_steps"]),
                    "exact_binary_reproduction": False,
                },
            }
            manifest_path = os.path.join(staging, REMOTE_MANIFEST_FILENAME)
            _atomic_write_json(manifest_path, manifest)
            _upload_verified(manifest_path, _object_remote(profile, archive_id, REMOTE_MANIFEST_FILENAME))
            manifest_digest = _sha256_file(manifest_path)
            complete_path = os.path.join(staging, REMOTE_COMPLETE_FILENAME)
            with open(complete_path, "w", encoding="ascii") as stream:
                stream.write(manifest_digest + "\n")
            _upload_verified(complete_path, _object_remote(profile, archive_id, REMOTE_COMPLETE_FILENAME))

        state = {
            "storage_schema_version": STORAGE_SCHEMA_VERSION,
            "archive_id": archive_id,
            "profile": profile["name"],
            "remote": profile["remote"],
            "label": manifest["label"],
            "archived_at": manifest["created_at"],
            "local_pruned": False,
            "restored_components": [],
            "retained_components": [],
        }
        _atomic_write_json(_state_path(target["root_path"]), state)
        if prune_local:
            _prune_archived_payload(target["root_path"], inventory, plan["offload_policy"])
            state["local_pruned"] = True
            state["pruned_at"] = _utc_now()
            state["offload_policy"] = plan["offload_policy"]
            latest_step = max(inventory["checkpoint_steps"], default=None)
            retained = list(plan["offload_policy"]["retained_components"])
            if plan["offload_policy"]["keep_latest_checkpoint"] and latest_step is not None:
                retained.append(f"checkpoint:{latest_step}")
            state["retained_components"] = retained
            _atomic_write_json(_state_path(target["root_path"]), state)
        print(
            f"[SUCCESS] {'Offloaded' if prune_local else 'Protected'} {target['root_path']} "
            f"as archive {archive_id}."
        )
        return manifest


def _register_existing_archive(target: dict, profile: dict, manifest: dict,
                               inventory: dict, plan: dict, *, prune_local: bool,
                               label: str = None, tags=None, notes: str = None) -> dict:
    """!
    @brief Point an artifact at an already-uploaded archive of its current content.

    @details Used when nothing has changed since a previous `protect` or `offload`. The
             remote object is immutable and is not rewritten; only this artifact's local
             marker is updated, and a label or tags supplied now are recorded there so
             the second invocation is not silently less descriptive than the first.
    @param[in] target Local artifact target.
    @param[in] profile Resolved storage profile.
    @param[in] manifest Remote manifest being reused.
    @param[in] inventory Current local inventory.
    @param[in] plan Storage plan carrying the requested offload policy.
    @param[in] prune_local Whether to prune the verified payload afterwards.
    @param[in] label Optional label supplied to this invocation.
    @param[in] tags Optional tags supplied to this invocation.
    @param[in] notes Optional free-text note supplied to this invocation.
    @return The reused remote manifest.
    """
    state = {
        "storage_schema_version": STORAGE_SCHEMA_VERSION,
        "archive_id": manifest["archive_id"],
        "profile": profile["name"],
        "remote": profile["remote"],
        "label": label or manifest.get("label") or os.path.basename(target["root_path"]),
        "archived_at": manifest.get("created_at"),
        "reused_existing_archive": True,
        "local_pruned": False,
        "restored_components": [],
        "retained_components": [],
    }
    parsed_tags = _parse_tags(tags)
    if parsed_tags:
        state["tags"] = parsed_tags
    if notes:
        state["notes"] = notes
    _atomic_write_json(_state_path(target["root_path"]), state)
    if prune_local:
        _prune_archived_payload(target["root_path"], inventory, plan["offload_policy"])
        state["local_pruned"] = True
        state["pruned_at"] = _utc_now()
        state["offload_policy"] = plan["offload_policy"]
        latest_step = max(inventory["checkpoint_steps"], default=None)
        retained = list(plan["offload_policy"]["retained_components"])
        if plan["offload_policy"]["keep_latest_checkpoint"] and latest_step is not None:
            retained.append(f"checkpoint:{latest_step}")
        state["retained_components"] = retained
        _atomic_write_json(_state_path(target["root_path"]), state)
    print(
        f"[SUCCESS] {'Offloaded' if prune_local else 'Protected'} {target['root_path']} "
        f"using existing archive {manifest['archive_id']}."
    )
    return manifest


def _prune_archived_payload(root: str, inventory: dict, policy: dict) -> None:
    """!
    @brief Remove only verified heavy payload while retaining control-plane files.
    @param[in] root Value supplied through the `root` argument.
    @param[in] inventory Value supplied through the `inventory` argument.
    @param[in] policy Normalized semantic retention policy.
    """
    latest_step = max(inventory["checkpoint_steps"], default=None)
    for entry in sorted(inventory["entries"], key=lambda item: item["path"], reverse=True):
        if entry["type"] == "directory" or _entry_retained_by_policy(entry, policy, latest_step):
            continue
        path = os.path.join(root, *entry["path"].split("/"))
        try:
            if os.path.islink(path) or os.path.isfile(path):
                os.remove(path)
        except FileNotFoundError:
            pass
    for entry in sorted(
        (item for item in inventory["entries"] if item["type"] == "directory"),
        key=lambda item: item["path"].count("/"), reverse=True,
    ):
        path = os.path.join(root, *entry["path"].split("/"))
        try:
            os.rmdir(path)
        except OSError:
            pass


def list_remote_manifests(profile: dict) -> list:
    """!
    @brief Enumerate completed archive manifests from the remote catalog.
    @param[in] profile Value supplied through the `profile` argument.
    @return Result produced by this operation.
    """
    root = _remote_join(profile["remote"], REMOTE_OBJECTS_DIRECTORY)
    result = _run_rclone([
        "lsf", root, "--recursive", "--files-only", "--include", f"*/{REMOTE_MANIFEST_FILENAME}"
    ])
    manifests = []
    for relative in sorted(line.strip() for line in result.stdout.splitlines() if line.strip()):
        archive_id = relative.split("/", 1)[0]
        if not ARCHIVE_ID_PATTERN.fullmatch(archive_id):
            continue
        try:
            manifests.append(_load_remote_manifest(profile, archive_id))
        except StorageError:
            continue
    return manifests


def verify_remote_archive(profile: dict, archive_id: str) -> dict:
    """!
    @brief Verify the completion marker and every stored chunk checksum.
    @param[in] profile Value supplied through the `profile` argument.
    @param[in] archive_id Value supplied through the `archive_id` argument.
    @return Result produced by this operation.
    """
    manifest = _load_remote_manifest(profile, archive_id)
    for chunk in manifest.get("chunks", []):
        actual = _remote_sha256(_chunk_remote_path(profile, archive_id, chunk))
        if actual != chunk.get("sha256"):
            raise StorageError(
                f"Archive {archive_id} chunk checksum mismatch: {chunk['name']} "
                f"(expected {chunk.get('sha256')}, got {actual})."
            )
    return manifest


def _validate_tar_members(archive: tarfile.TarFile) -> None:
    """!
    @brief Reject archive members that could escape the restore destination.
    @param[in] archive Value supplied through the `archive` argument.
    """
    link_paths = set()
    members = archive.getmembers()
    for member in members:
        normalized = os.path.normpath(member.name.replace("\\", "/"))
        if normalized.startswith("../") or normalized == ".." or os.path.isabs(normalized):
            raise StorageError(f"Unsafe archive member path: {member.name}")
        for link_path in link_paths:
            if normalized == link_path or normalized.startswith(link_path.rstrip("/") + "/"):
                raise StorageError(f"Archive member traverses an earlier symlink: {member.name}")
        if member.issym() or member.islnk():
            if member.islnk():
                link_target = os.path.normpath(member.linkname.replace("\\", "/"))
                if link_target.startswith("../") or link_target == ".." or os.path.isabs(link_target):
                    raise StorageError(f"Unsafe archive hardlink target: {member.linkname}")
            link_paths.add(normalized)


def _extract_chunk(path: str, destination: str) -> None:
    """!
    @brief Safely extract one verified tar chunk into a staging tree.
    @param[in] path Value supplied through the `path` argument.
    @param[in] destination Value supplied through the `destination` argument.
    """
    with tarfile.open(path, "r:*") as archive:
        _validate_tar_members(archive)
        try:
            archive.extractall(destination, filter="data")
        except TypeError:
            archive.extractall(destination)


def _restore_study_context(manifest: dict, case_destination: str) -> None:
    """!
    @brief Recreate missing study control-plane files around a restored member.
    @param[in] manifest Value supplied through the `manifest` argument.
    @param[in] case_destination Value supplied through the `case_destination` argument.
    """
    if manifest.get("artifact_type") != "study-case":
        return
    study_root = Path(case_destination).parent.parent
    for item in manifest.get("study_context", []):
        relative = item.get("path")
        if not isinstance(relative, str):
            continue
        destination = study_root.joinpath(*relative.split("/"))
        if destination.exists():
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(base64.b64decode(item.get("content_base64", "")))
        try:
            destination.chmod(int(item.get("mode", 0o644)))
        except OSError:
            pass


def _merge_tree(source: str, destination: str) -> None:
    """!
    @brief Merge a verified restore tree into a known cold artifact skeleton.
    @param[in] source Value supplied through the `source` argument.
    @param[in] destination Value supplied through the `destination` argument.
    """
    os.makedirs(destination, exist_ok=True)
    for entry in os.scandir(source):
        target = os.path.join(destination, entry.name)
        if entry.is_symlink():
            if os.path.lexists(target):
                if os.path.isdir(target) and not os.path.islink(target):
                    shutil.rmtree(target)
                else:
                    os.remove(target)
            os.symlink(os.readlink(entry.path), target)
        elif entry.is_dir(follow_symlinks=False):
            _merge_tree(entry.path, target)
        else:
            os.makedirs(os.path.dirname(target), exist_ok=True)
            shutil.copy2(entry.path, target)


def _download_archive_components(profile: dict, manifest: dict, components: set,
                                 destination: str, workers: int = None) -> None:
    """!
    @brief Download and safely merge selected semantic archive components.
    @param[in] profile Active storage profile.
    @param[in] manifest Verified remote archive manifest.
    @param[in] components Semantic components to retrieve.
    @param[in] destination Directory receiving merged content.
    @param[in] workers Optional parallel worker count.
    @return None.
    """
    chunks = [
        chunk for chunk in manifest.get("chunks", [])
        if str(chunk.get("component", "")) in components
    ]
    if not chunks:
        raise StorageError(
            f"Archive {manifest.get('archive_id')} has none of the requested components: "
            + ", ".join(sorted(components))
        )
    worker_count = int(workers or profile.get("workers", DEFAULT_STORAGE_WORKERS))
    if worker_count <= 0:
        raise StorageError("Restore workers must be positive.")
    os.makedirs(destination, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="picurv-component-restore-") as staging:
        def download_and_extract(index_chunk):
            """!
            @brief Download, verify, and extract one selected chunk.
            @param[in] index_chunk Tuple of chunk index and manifest entry.
            @return Chunk index and extraction directory.
            """
            index, chunk = index_chunk
            chunk_root = os.path.join(staging, f"extract-{index:05d}")
            os.makedirs(chunk_root)
            local_chunk = os.path.join(staging, chunk["name"])
            # Schema 1 archives keep their payload beside the manifest; schema 2 and
            # later resolve it out of the shared content-addressed store.
            _run_rclone([
                "copyto",
                _chunk_remote_path(profile, manifest["archive_id"], chunk),
                local_chunk,
            ])
            actual = _sha256_file(local_chunk)
            if actual != chunk.get("sha256"):
                raise StorageError(
                    f"Downloaded chunk checksum mismatch: {chunk['name']} "
                    f"(expected {chunk.get('sha256')}, got {actual})."
                )
            _extract_chunk(local_chunk, chunk_root)
            return index, chunk_root

        task_workers = min(worker_count, len(chunks))
        if task_workers == 1:
            extracted = [download_and_extract(item) for item in enumerate(chunks)]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=task_workers) as executor:
                extracted = list(executor.map(download_and_extract, enumerate(chunks)))
        for _, chunk_root in sorted(extracted):
            _merge_tree(chunk_root, destination)


def restore_missing_workspace_assets(workspace_root: str, provider_spec_hashes) -> dict:
    """!
    @brief Recover immutable workspace assets from archived run-local input copies.
    @param[in] workspace_root Initialized workspace receiving recovered objects.
    @param[in] provider_spec_hashes Required provider-specification hashes.
    @return Provider-specification hash to recovered asset-reference mapping.
    """
    requested = {str(value) for value in provider_spec_hashes if value}
    if not requested:
        return {}
    workspace = Path(os.path.abspath(workspace_root))
    (workspace / "assets").mkdir(parents=True, exist_ok=True)
    profile = load_storage_profile(config_path=str(workspace / STORAGE_CONFIG_FILENAME))
    matches = {}
    for manifest in sorted(
        list_remote_manifests(profile), key=lambda item: item.get("created_at", ""), reverse=True
    ):
        for reference in manifest.get("run_assets", []):
            spec_hash = reference.get("provider_spec_sha256")
            if spec_hash in requested and spec_hash not in matches:
                matches[spec_hash] = (manifest, reference)
        if requested <= set(matches):
            break
    if not matches:
        return {}

    restored = {}
    by_archive = {}
    for spec_hash, (manifest, reference) in matches.items():
        group = by_archive.setdefault(manifest["archive_id"], {"manifest": manifest, "items": []})
        group["items"].append((spec_hash, reference))
    kind_directories = {
        "grid": "grids",
        "initial-condition": "initial_conditions",
        "inlet-profiles": "inlet_profiles",
    }
    for group in by_archive.values():
        with tempfile.TemporaryDirectory(prefix=".asset-restore-", dir=workspace / "assets") as temporary:
            _download_archive_components(profile, group["manifest"], {"inputs"}, temporary)
            for spec_hash, reference in group["items"]:
                kind = reference.get("kind")
                asset_id = reference.get("asset_id")
                if kind not in kind_directories or not asset_id:
                    continue
                object_relative = f"assets/objects/{kind_directories[kind]}/{asset_id}"
                object_root = workspace / object_relative
                if not object_root.is_dir():
                    object_root.parent.mkdir(parents=True, exist_ok=True)
                    object_staging = Path(tempfile.mkdtemp(
                        prefix=f".{asset_id[:12]}-", dir=object_root.parent
                    ))
                    try:
                        files = []
                        for item in reference.get("files", []):
                            source = Path(temporary, *item["archived_path"].split("/"))
                            if not source.is_file() or (
                                item.get("sha256") and _sha256_file(str(source)) != item["sha256"]
                            ):
                                raise StorageError(
                                    f"Archived asset payload is missing or corrupt: {item['archived_path']}"
                                )
                            destination = object_staging / "payload" / item["path"]
                            destination.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy2(source, destination)
                            files.append({
                                "path": item["path"],
                                "bytes": destination.stat().st_size,
                                "sha256": _sha256_file(str(destination)),
                            })
                        _atomic_write_json(str(object_staging / "asset.json"), {
                            "schema_version": 1,
                            "asset_id": asset_id,
                            "kind": kind,
                            "provider": reference.get("provider"),
                            "provider_spec_sha256": spec_hash,
                            "recovered_from_archive": group["manifest"]["archive_id"],
                            "files": files,
                        })
                        try:
                            os.replace(object_staging, object_root)
                        except OSError as exc:
                            if exc.errno != errno.ENOTEMPTY and not object_root.is_dir():
                                raise
                    finally:
                        if object_staging.is_dir():
                            shutil.rmtree(object_staging, ignore_errors=True)
                restored[spec_hash] = {
                    "asset_id": asset_id,
                    "kind": kind,
                    "provider": reference.get("provider"),
                    "provider_spec_sha256": spec_hash,
                    "object": object_relative,
                    "files": reference.get("files", []),
                }
                print(f"[INFO] Restored shared asset {kind}: {asset_id}")
    return restored


def _rebase_restored_text_paths(root: str, replacements: list) -> list:
    """!
    @brief Rebase known generated text artifacts after an explicit relocated restore.
    @param[in] root Value supplied through the `root` argument.
    @param[in] replacements Value supplied through the `replacements` argument.
    @return Result produced by this operation.
    """
    changed = []
    allowed_suffixes = {".control", ".run", ".sbatch", ".json", ".tsv", ".yml", ".yaml"}
    for path in Path(root).rglob("*"):
        if not path.is_file() or path.is_symlink() or path.stat().st_size > 32 * 1024 * 1024:
            continue
        if path.suffix.lower() not in allowed_suffixes:
            continue
        try:
            content = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        updated = content
        for old, new in replacements:
            if old and new and old != new:
                updated = updated.replace(old, new)
        if updated != content:
            path.write_text(updated, encoding="utf-8")
            changed.append(str(path))
    return changed


def restore_archive(profile: dict, archive_id: str, destination: str = None,
                    checkpoints=None, force: bool = False,
                    workers: int = None, components=None) -> dict:
    """!
    @brief Download, verify, extract, and materialize an archive or selected checkpoints.
    @param[in] profile Value supplied through the `profile` argument.
    @param[in] archive_id Value supplied through the `archive_id` argument.
    @param[in] destination Value supplied through the `destination` argument.
    @param[in] checkpoints Value supplied through the `checkpoints` argument.
    @param[in] force Value supplied through the `force` argument.
    @param[in] workers Optional parallel download and extraction worker count.
    @param[in] components Optional semantic components to restore.
    @return Result produced by this operation.
    """
    manifest = _load_remote_manifest(profile, archive_id)
    original = os.path.abspath(manifest["original_path"])
    destination_abs = os.path.abspath(destination or original)
    selected_steps = {int(step) for step in (checkpoints or [])}
    selected_components = {str(component) for component in (components or [])}
    chunks = []
    for chunk in manifest.get("chunks", []):
        component = str(chunk.get("component", ""))
        if selected_steps or selected_components:
            if component == "metadata":
                chunks.append(chunk)
            elif component.startswith("checkpoint:") and int(component.split(":", 1)[1]) in selected_steps:
                chunks.append(chunk)
            elif component in selected_components:
                chunks.append(chunk)
        else:
            chunks.append(chunk)
    if selected_steps:
        available = set(manifest.get("checkpoint_steps", []))
        missing = selected_steps - available
        if missing:
            raise StorageError(f"Archive {archive_id} does not contain checkpoint step(s): {sorted(missing)}")
    existing_state = read_storage_state(destination_abs) if os.path.isdir(destination_abs) else None
    if os.path.exists(destination_abs) and not force:
        if not existing_state or existing_state.get("archive_id") != archive_id:
            raise StorageError(
                f"Restore destination already exists and is not the matching cold artifact: {destination_abs}. "
                "Choose --to or use --force after verifying the destination."
            )
    parent = os.path.dirname(destination_abs)
    os.makedirs(parent, exist_ok=True)
    temporary = tempfile.mkdtemp(prefix=f".picurv-restore-{archive_id[:8]}-", dir=parent)
    materialized = os.path.join(temporary, "materialized")
    os.makedirs(materialized)
    try:
        newly_restored = (
            [f"checkpoint:{step}" for step in sorted(selected_steps)]
            + sorted(selected_components)
        )
        prior_restored = (existing_state or {}).get("restored_components") or []
        worker_count = int(workers or profile.get("workers", DEFAULT_STORAGE_WORKERS))
        if worker_count <= 0:
            raise StorageError("Restore workers must be positive.")

        def download_and_extract(index_chunk):
            """!
            @brief Download, verify, and extract one archive chunk.
            @param[in] index_chunk Tuple of chunk index and manifest entry.
            @return Chunk index and extraction directory.
            """
            index, chunk = index_chunk
            chunk_root = os.path.join(temporary, f"extract-{index:05d}")
            os.makedirs(chunk_root)
            local_chunk = os.path.join(temporary, chunk["name"])
            print(f"[INFO] Restoring chunk {index + 1}/{len(chunks)}: {chunk['component']}")
            _run_rclone([
                "copyto", _chunk_remote_path(profile, archive_id, chunk), local_chunk,
            ])
            actual = _sha256_file(local_chunk)
            if actual != chunk.get("sha256"):
                raise StorageError(
                    f"Downloaded chunk checksum mismatch: {chunk['name']} "
                    f"(expected {chunk.get('sha256')}, got {actual})."
                )
            _extract_chunk(local_chunk, chunk_root)
            os.remove(local_chunk)
            return index, chunk_root

        task_workers = min(worker_count, max(1, len(chunks)))
        if task_workers == 1:
            extracted = [download_and_extract(item) for item in enumerate(chunks)]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=task_workers) as executor:
                extracted = list(executor.map(download_and_extract, enumerate(chunks)))
        for _, chunk_root in sorted(extracted):
            _merge_tree(chunk_root, materialized)

        if os.path.isdir(destination_abs):
            _merge_tree(materialized, destination_abs)
        else:
            os.replace(materialized, destination_abs)
        _restore_study_context(manifest, destination_abs)

        replacements = [(original, destination_abs)]
        original_study = manifest.get("original_study_path")
        if original_study and manifest.get("artifact_type") == "study-case":
            replacements.append((os.path.abspath(original_study), str(Path(destination_abs).parent.parent)))
        rebased = _rebase_restored_text_paths(destination_abs, replacements)
        state = {
            "storage_schema_version": STORAGE_SCHEMA_VERSION,
            "archive_id": archive_id,
            "profile": profile["name"],
            "remote": profile["remote"],
            "label": manifest.get("label"),
            "archived_at": manifest.get("created_at"),
            "restored_at": _utc_now(),
            "local_pruned": bool(selected_steps or selected_components),
            "restored_components": sorted(set(prior_restored) | set(newly_restored)),
            "retained_components": (existing_state or {}).get("retained_components") or [],
            "relocated_from": original if original != destination_abs else None,
            "rebased_files": len(rebased),
        }
        _atomic_write_json(_state_path(destination_abs), state)
    except Exception:
        if os.path.isdir(temporary):
            shutil.rmtree(temporary, ignore_errors=True)
        raise
    if os.path.isdir(temporary):
        shutil.rmtree(temporary, ignore_errors=True)
    print(f"[SUCCESS] Restored archive {archive_id} to {destination_abs}.")
    if rebased:
        print(f"[INFO] Rebased {len(rebased)} generated text artifact(s) to the restored path.")
    return manifest


def _expand_checkpoint_selection(explicit, ranges) -> list:
    """!
    @brief Combine repeated `--checkpoint` values with `START:END[:STRIDE]` ranges.
    @param[in] explicit Individually named steps, or None.
    @param[in] ranges Range expressions, or None.
    @return Sorted unique step list, or None when nothing was selected.
    """
    selected = set(explicit or ())
    for expression in ranges or ():
        parts = str(expression).split(":")
        if len(parts) not in (2, 3):
            raise StorageError(
                f"--checkpoints expects START:END or START:END:STRIDE, got {expression!r}."
            )
        try:
            numbers = [int(part) for part in parts]
        except ValueError:
            raise StorageError(
                f"--checkpoints expects integer steps, got {expression!r}."
            ) from None
        start, end = numbers[0], numbers[1]
        stride = numbers[2] if len(numbers) == 3 else 1
        if stride <= 0:
            raise StorageError(f"--checkpoints stride must be positive, got {expression!r}.")
        if end < start:
            raise StorageError(f"--checkpoints end precedes start in {expression!r}.")
        selected.update(range(start, end + 1, stride))
    return sorted(selected) or None


def resolve_workspace_archive_id(profile: dict, workspace_id: str) -> str:
    """!
    @brief Find the newest complete workspace archive carrying one workspace identity.
    @param[in] profile Resolved storage profile.
    @param[in] workspace_id Workspace identity recorded at archive time.
    @return Archive id.
    @throws StorageError when no such archive exists.
    """
    matches = [
        manifest for manifest in list_remote_manifests(profile)
        if manifest.get("artifact_type") == "workspace"
        and manifest.get("workspace_id") == workspace_id
    ]
    if not matches:
        raise StorageError(f"No workspace archive found for identity {workspace_id!r}.")
    return max(matches, key=lambda item: str(item.get("created_at", "")))["archive_id"]


def _resolve_archive_id_from_args(args) -> str:
    """!
    @brief Resolve an explicit archive ID or a local artifact marker.
    @param[in] args Value supplied through the `args` argument.
    @return Result produced by this operation.
    """
    archive_id = getattr(args, "archive_id", None)
    if archive_id:
        return archive_id
    workspace_id = getattr(args, "workspace_id", None)
    if workspace_id:
        return resolve_workspace_archive_id(
            load_storage_profile(
                getattr(args, "profile", None), getattr(args, "storage_config", None)
            ),
            workspace_id,
        )
    targets = resolve_local_storage_targets(
        getattr(args, "run_dir", None), getattr(args, "study_dir", None), getattr(args, "case_ids", None)
    )
    if len(targets) != 1:
        raise StorageError("Restore/verify by local marker requires exactly one target.")
    state = read_storage_state(targets[0]["root_path"])
    if not state or not state.get("archive_id"):
        raise StorageError(f"No storage marker with an archive ID exists under {targets[0]['root_path']}.")
    return state["archive_id"]


def _render_status(inventory: dict) -> None:
    """!
    @brief Render one artifact status row and safety details.
    @param[in] inventory Value supplied through the `inventory` argument.
    """
    target = inventory["target"]
    storage = inventory["storage"]
    activity = "BUSY" if inventory["active_locks"] or inventory["slurm"]["active"] else storage["state"]
    identity = (
        target.get("case_id") or target.get("run_id") or target.get("study_id")
        or target.get("workspace_id") or os.path.basename(target["root_path"])
    )
    print(
        f"{identity:<36} {target['artifact_type']:<12} {activity:<10} "
        f"{_human_bytes(inventory['total_bytes']):>12}  {storage.get('label') or ''}"
    )


def storage_setup_workflow(args) -> None:
    """!
    @brief Create or update a non-secret workspace storage profile.
    @param[in] args Value supplied through the `args` argument.
    """
    config_path = resolve_storage_config_path(getattr(args, "storage_config", None), require=False)
    payload = {}
    if os.path.isfile(config_path):
        with open(config_path, "r", encoding="utf-8") as stream:
            payload = yaml.safe_load(stream) or {}
    profiles = payload.setdefault("profiles", {})
    profile_name = args.profile or DEFAULT_PROFILE_NAME
    profile = {
        "remote": args.remote.rstrip("/"),
        "compression": args.compression,
        "chunk_size_gib": args.chunk_size_gib,
        "workers": getattr(args, "workers", DEFAULT_STORAGE_WORKERS),
        "offload_policy": getattr(args, "offload_policy", "metadata-only"),
        "keep_latest_checkpoint": bool(getattr(args, "keep_latest_checkpoint", False)),
    }
    if args.staging_directory:
        profile["staging_directory"] = os.path.abspath(os.path.expanduser(args.staging_directory))
    payload["default_profile"] = profile_name
    profiles[profile_name] = profile
    print(f"[INFO] Storage config : {config_path}")
    print(f"[INFO] Profile        : {profile_name}")
    print(f"[INFO] Remote         : {profile['remote']}")
    if args.dry_run:
        print("[INFO] Dry-run only. No configuration or remote directories were changed.")
        return
    _run_rclone(["mkdir", _remote_join(profile["remote"], REMOTE_OBJECTS_DIRECTORY)])
    os.makedirs(os.path.dirname(config_path), exist_ok=True)
    temporary = f"{config_path}.tmp.{os.getpid()}"
    with open(temporary, "w", encoding="utf-8") as stream:
        yaml.safe_dump(payload, stream, sort_keys=False)
    os.replace(temporary, config_path)
    print("[SUCCESS] Storage profile configured and remote access verified.")


def storage_status_workflow(args) -> None:
    """!
    @brief Print local storage and lifecycle status for selected artifacts.
    @param[in] args Value supplied through the `args` argument.
    """
    if args.study_dir and not args.case_ids:
        study_root = os.path.abspath(args.study_dir)
        case_ids = [
            path.name for path in sorted((Path(study_root) / "cases").glob("case_*")) if path.is_dir()
        ]
        targets = resolve_local_storage_targets(None, study_root, case_ids) if case_ids else resolve_local_storage_targets(None, study_root)
    else:
        targets = resolve_local_storage_targets(
        args.run_dir, args.study_dir, args.case_ids,
        workspace=getattr(args, "workspace", None),
        include_inputs=getattr(args, "include_inputs", False),
        completed=getattr(args, "completed", False),
    )
    inventories = [inspect_artifact(target) for target in targets]
    if args.output_format == "json":
        serializable = []
        for item in inventories:
            copy_item = dict(item)
            copy_item.pop("entries", None)
            serializable.append(copy_item)
        print(json.dumps(serializable, indent=2, sort_keys=True))
        return
    print(f"{'ARTIFACT':<36} {'TYPE':<12} {'STATE':<10} {'LOCAL SIZE':>12}  LABEL")
    for inventory in inventories:
        _render_status(inventory)


def storage_plan_workflow(args) -> None:
    """!
    @brief Render a read-only packaging and safety plan.
    @param[in] args Value supplied through the `args` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    targets = resolve_local_storage_targets(
        args.run_dir, args.study_dir, args.case_ids,
        workspace=getattr(args, "workspace", None),
        include_inputs=getattr(args, "include_inputs", False),
        completed=getattr(args, "completed", False),
    )
    for index, target in enumerate(targets):
        if index:
            print()
        plan = build_storage_plan(
            target, profile, args.compression,
            policy=getattr(args, "policy", None),
            keep_latest_checkpoint=getattr(args, "keep_latest_checkpoint", None),
            workers=getattr(args, "workers", None),
        )
        _render_plan(plan)
        _assert_archive_safe(plan["inventory"])


def storage_archive_workflow(args, prune_local: bool) -> None:
    """!
    @brief Execute protect or offload for one or more explicit local targets.
    @param[in] args Value supplied through the `args` argument.
    @param[in] prune_local Value supplied through the `prune_local` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    targets = resolve_local_storage_targets(
        args.run_dir, args.study_dir, args.case_ids,
        workspace=getattr(args, "workspace", None),
        include_inputs=getattr(args, "include_inputs", False),
        completed=getattr(args, "completed", False),
    )
    for target in targets:
        if args.dry_run:
            _render_plan(build_storage_plan(
                target, profile, args.compression,
                policy=getattr(args, "policy", None),
                keep_latest_checkpoint=getattr(args, "keep_latest_checkpoint", None),
                workers=getattr(args, "workers", None),
            ))
            print("[INFO] Dry-run only. No files were packaged, uploaded, or pruned.")
            continue
        archive_artifact(
            target,
            profile,
            label=args.label,
            tags=args.tags,
            compression=args.compression,
            prune_local=prune_local,
            policy=getattr(args, "policy", None),
            keep_latest_checkpoint=getattr(args, "keep_latest_checkpoint", None),
            workers=getattr(args, "workers", None),
            notes=getattr(args, "notes", None),
        )


def storage_restore_workflow(args) -> None:
    """!
    @brief Restore a remote archive by globally unique ID or local marker.
    @param[in] args Value supplied through the `args` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    archive_id = _resolve_archive_id_from_args(args)
    restore_archive(
        profile,
        archive_id,
        destination=args.destination,
        checkpoints=_expand_checkpoint_selection(
            args.checkpoints, getattr(args, "checkpoint_ranges", None)
        ),
        force=args.force,
        workers=getattr(args, "workers", None),
        components=getattr(args, "components", None),
    )


def storage_prune_workflow(args) -> None:
    """!
    @brief Remove local asset objects nothing local needs and storage has verified.
    @param[in] args Value supplied through the `args` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    workspace_root = os.path.abspath(args.workspace) if args.workspace else _find_upwards(
        os.getcwd(), WORKSPACE_CONFIG_FILENAME
    )
    if workspace_root and os.path.isfile(workspace_root):
        workspace_root = os.path.dirname(workspace_root)
    if not workspace_root:
        raise StorageError(
            "No initialized workspace found. Pass --workspace, or run inside one."
        )
    decisions = prune_unused_workspace_assets(workspace_root, profile, dry_run=args.dry_run)
    if not decisions:
        print("[INFO] The workspace asset store holds no published objects.")
        return
    for decision in decisions:
        print(
            f"{decision['asset_id'][:16]}  {decision['kind']}\n"
            f"  referenced by remote runs       {decision['cold_runs']}\n"
            f"  referenced by active local runs {decision['active_local_runs']}\n"
            f"  remote protection               {decision['remote_protection']}\n"
            f"  local removal                   {decision['local_removal']}"
        )
    removed = [item for item in decisions if item.get("removed")]
    if args.dry_run:
        safe = [item for item in decisions if item["local_removal"] == "safe"]
        print(f"[INFO] Dry-run only. {len(safe)} object(s) would be removed.")
    else:
        print(f"[SUCCESS] Removed {len(removed)} local asset object(s).")


def storage_verify_workflow(args) -> None:
    """!
    @brief Verify all remote chunks for one archive.
    @param[in] args Value supplied through the `args` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    archive_id = _resolve_archive_id_from_args(args)
    manifest = verify_remote_archive(profile, archive_id)
    print(
        f"[SUCCESS] Archive {archive_id} is complete; verified {len(manifest.get('chunks', []))} chunk(s)."
    )


def storage_list_workflow(args) -> None:
    """!
    @brief Search the remote manifest catalog without local artifact state.
    @param[in] args Value supplied through the `args` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    manifests = list_remote_manifests(profile)
    workspace_label = getattr(args, "workspace_label", None)
    if workspace_label:
        manifests = [
            item for item in manifests
            if str(item.get("workspace_id") or "") == workspace_label
        ]
    query = str(args.search or "").lower()
    if query:
        manifests = [
            item for item in manifests
            if query in " ".join(
                str(item.get(key, "")) for key in (
                    "archive_id", "label", "notes", "run_id", "study_id", "case_id",
                    "workspace_id", "tags",
                )
            ).lower()
        ]
    if args.output_format == "json":
        print(json.dumps(manifests, indent=2, sort_keys=True))
        return
    print(f"{'ARCHIVE ID':<34} {'TYPE':<12} {'IDENTITY':<32} LABEL")
    for item in manifests:
        identity = (
            item.get("case_id") or item.get("run_id") or item.get("study_id")
            or item.get("workspace_id") or "-"
        )
        print(f"{item['archive_id']:<34} {item.get('artifact_type', '-'):<12} {identity:<32} {item.get('label', '')}")


def storage_show_workflow(args) -> None:
    """!
    @brief Show one complete remote archive manifest.
    @param[in] args Value supplied through the `args` argument.
    """
    profile = load_storage_profile(args.profile, args.storage_config)
    manifest = _load_remote_manifest(profile, args.archive_id)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def storage_workflow(args) -> None:
    """!
    @brief Dispatch nested storage actions using existing PICurv workflow conventions.
    @param[in] args Value supplied through the `args` argument.
    """
    try:
        action = args.storage_action
        if action == "setup":
            storage_setup_workflow(args)
        elif action == "status":
            storage_status_workflow(args)
        elif action == "plan":
            storage_plan_workflow(args)
        elif action == "protect":
            storage_archive_workflow(args, prune_local=False)
        elif action == "offload":
            storage_archive_workflow(args, prune_local=True)
        elif action == "restore":
            storage_restore_workflow(args)
        elif action == "prune":
            storage_prune_workflow(args)
        elif action == "verify":
            storage_verify_workflow(args)
        elif action == "list":
            storage_list_workflow(args)
        elif action == "show":
            storage_show_workflow(args)
        else:
            raise StorageError(f"Unsupported storage action: {action}")
    except StorageError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        raise SystemExit(1)


def add_storage_parser(subparsers) -> argparse.ArgumentParser:
    """!
    @brief Attach the nested storage command parser to PICurv's top-level parser.
    @param[in] subparsers Value supplied through the `subparsers` argument.
    @return Result produced by this operation.
    """
    parser = subparsers.add_parser(
        "storage",
        help="Protect, offload, inspect, verify, and restore run/study artifacts.",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Manage PICurv run and study data through a configured rclone remote.\n"
            "Remote archives are checksum-verified before local payload can be pruned.\n\n"
            "Examples:\n"
            "  picurv storage setup --remote labstore:picurv-data\n"
            "  picurv storage status --run-dir runs/my_run\n"
            "  picurv storage protect --run-dir runs/my_run --label 'baseline'\n"
            "  picurv storage offload --study-dir studies/my_study --case-id case_0003\n"
            "  picurv storage list --search '64-grid'\n"
            "  picurv storage restore --archive-id <id>"
        ),
        epilog="Use `picurv storage <action> --help` for action-specific controls.",
    )
    actions = parser.add_subparsers(dest="storage_action", required=True, help="Storage action")

    setup = actions.add_parser("setup", help="Configure a non-secret rclone storage profile.")
    setup.add_argument("--remote", required=True, help="Rclone remote and base path, such as labstore:picurv-data.")
    setup.add_argument("--profile", default=DEFAULT_PROFILE_NAME, help="Profile name (default: archive).")
    setup.add_argument("--storage-config", help=f"Storage YAML path (default: ./{STORAGE_CONFIG_FILENAME}).")
    setup.add_argument("--compression", choices=list(STORAGE_COMPRESSION_POLICIES), default="auto")
    setup.add_argument("--chunk-size-gib", type=float, default=DEFAULT_CHUNK_SIZE_GIB)
    setup.add_argument("--workers", type=int, default=DEFAULT_STORAGE_WORKERS,
                       help="CPU workers for compression and restoration.")
    setup.add_argument("--offload-policy", choices=list(STORAGE_OFFLOAD_POLICIES), default="metadata-only")
    setup.add_argument("--keep-latest-checkpoint", action="store_true",
                       help="Retain the newest committed checkpoint after offload.")
    setup.add_argument("--staging-directory", help="Optional local directory for one archive chunk at a time.")
    setup.add_argument("--dry-run", action="store_true")

    def add_profile_options(action_parser):
        """!
        @brief Attach shared storage-profile selectors to one action parser.
        @param[in] action_parser Value supplied through the `action_parser` argument.
        """
        action_parser.add_argument("--profile", help="Configured storage profile name.")
        action_parser.add_argument("--storage-config", help="Explicit storage YAML path.")

    def add_local_target(action_parser, require=True, allow_workspace=True):
        """!
        @brief Attach the standard run/study/workspace target selectors to one parser.
        @param[in] action_parser Value supplied through the `action_parser` argument.
        @param[in] require Value supplied through the `require` argument.
        @param[in] allow_workspace Whether a whole workspace is a valid target here.
        """
        group = action_parser.add_mutually_exclusive_group(required=require)
        group.add_argument("--run-dir", help="Standalone run directory.")
        group.add_argument("--study-dir", help="Sweep study directory.")
        if allow_workspace:
            group.add_argument(
                "--workspace",
                help="Workspace root: its configuration, catalog, and assets, not its "
                     "runs and studies, which are their own artifacts.",
            )
        action_parser.add_argument(
            "--case-id", dest="case_ids", action="append",
            help="One numbered study member, such as case_0003; repeat to select several.",
        )
        action_parser.add_argument(
            "--completed", action="store_true",
            help="With --study-dir, select every finished member and skip the rest.",
        )

    status = actions.add_parser("status", help="Show local, protected, cold, and busy artifact state.")
    add_local_target(status)
    status.add_argument("--format", dest="output_format", choices=list(CLI_OUTPUT_FORMATS), default="text")

    plan = actions.add_parser("plan", help="Show packaging, dependencies, and safety checks without writing.")
    add_local_target(plan)
    add_profile_options(plan)
    plan.add_argument("--compression", choices=list(STORAGE_COMPRESSION_POLICIES))
    plan.add_argument("--policy", choices=list(STORAGE_OFFLOAD_POLICIES))
    plan.add_argument("--workers", type=int)
    plan_checkpoint = plan.add_mutually_exclusive_group()
    plan_checkpoint.add_argument("--keep-latest-checkpoint", dest="keep_latest_checkpoint", action="store_true")
    plan_checkpoint.add_argument("--drop-all-checkpoints", dest="keep_latest_checkpoint", action="store_false")
    plan.set_defaults(keep_latest_checkpoint=None)

    for name, help_text in (
        ("protect", "Upload and verify an archive while retaining all local files."),
        ("offload", "Upload and verify an archive, then prune heavy local payload."),
    ):
        action_parser = actions.add_parser(name, help=help_text)
        add_local_target(action_parser)
        add_profile_options(action_parser)
        action_parser.add_argument(
            "--include-inputs", action="store_true",
            help="With --workspace, also archive user-supplied files under inputs/.",
        )
        action_parser.add_argument("--label", help="Human-readable searchable label.")
        action_parser.add_argument(
            "--notes", help="Free-text note recorded with the archive and shown by `show`."
        )
        action_parser.add_argument("--tag", dest="tags", action="append", help="Repeatable KEY=VALUE catalog tag.")
        action_parser.add_argument("--compression", choices=list(STORAGE_COMPRESSION_POLICIES))
        action_parser.add_argument("--policy", choices=list(STORAGE_OFFLOAD_POLICIES))
        action_parser.add_argument("--workers", type=int)
        checkpoint_group = action_parser.add_mutually_exclusive_group()
        checkpoint_group.add_argument("--keep-latest-checkpoint", dest="keep_latest_checkpoint", action="store_true")
        checkpoint_group.add_argument("--drop-all-checkpoints", dest="keep_latest_checkpoint", action="store_false")
        action_parser.set_defaults(keep_latest_checkpoint=None)
        action_parser.add_argument("--dry-run", action="store_true")

    restore = actions.add_parser("restore", help="Restore a complete archive or selected checkpoints.")
    restore_source = restore.add_mutually_exclusive_group(required=True)
    restore_source.add_argument("--archive-id", help="Globally unique remote archive ID.")
    restore_source.add_argument(
        "--workspace-id", help="Restore a workspace archive by its recorded workspace identity."
    )
    restore_source.add_argument("--run-dir", help="Cold run containing a local storage marker.")
    restore_source.add_argument("--study-dir", help="Cold study containing a local storage marker.")
    restore.add_argument("--case-id", dest="case_ids", action="append")
    add_profile_options(restore)
    restore.add_argument("--to", dest="destination", help="Optional alternate restore destination.")
    restore.add_argument(
        "--checkpoint", dest="checkpoints", action="append", type=int,
        help="One committed step; repeat to select several.",
    )
    restore.add_argument(
        "--checkpoints", dest="checkpoint_ranges", action="append",
        help="An inclusive step range as START:END or START:END:STRIDE; repeatable.",
    )
    restore.add_argument(
        "--component", dest="components", action="append",
        choices=STORAGE_RESTORE_COMPONENTS,
        help="Restore one semantic component; repeat as needed.",
    )
    restore.add_argument("--force", action="store_true", help="Allow merge into a non-matching existing destination.")
    restore.add_argument("--workers", type=int, help="Parallel download/extraction workers.")

    prune = actions.add_parser(
        "prune", help="Remove verified local asset objects that nothing local still needs."
    )
    prune.add_argument("--workspace", help="Workspace root; defaults to discovery from the cwd.")
    prune.add_argument(
        "--assets", action="store_true", required=True,
        help="Select the workspace asset store. Required: prune removes nothing else.",
    )
    prune.add_argument(
        "--unused-locally", action="store_true", required=True,
        help="Confirm that only objects with no active local run are removed.",
    )
    prune.add_argument("--dry-run", action="store_true", help="Report the decision only.")
    add_profile_options(prune)

    verify = actions.add_parser("verify", help="Verify a remote archive completion marker and chunk checksums.")
    verify_source = verify.add_mutually_exclusive_group(required=True)
    verify_source.add_argument("--archive-id")
    verify_source.add_argument("--run-dir")
    verify_source.add_argument("--study-dir")
    verify.add_argument("--case-id", dest="case_ids", action="append")
    add_profile_options(verify)

    list_parser = actions.add_parser("list", help="List/search remote archives without local directories.")
    add_profile_options(list_parser)
    list_parser.add_argument("--search", help="Case-insensitive search across IDs, labels, identities, and tags.")
    list_parser.add_argument(
        "--workspace-label", help="Show only archives belonging to this workspace identity."
    )
    list_parser.add_argument("--format", dest="output_format", choices=list(CLI_OUTPUT_FORMATS), default="text")

    show = actions.add_parser("show", help="Print the complete manifest for one archive.")
    show.add_argument("--archive-id", required=True)
    add_profile_options(show)
    return parser
