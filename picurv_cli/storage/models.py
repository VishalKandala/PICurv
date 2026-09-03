"""!
@file models.py
@brief Storage constants, profiles, and the local artifact state marker.
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


#: Components no offload policy may prune, whatever it retains. "assets" is here
#: because no named policy lists it as retained - an ordinary offload would otherwise
#: delete the entire local asset store unconditionally, bypassing the reference-aware
#: check that `storage prune --assets --unused-locally` performs before removing an
#: object. Reclaiming assets is that command's job, deliberately, never an offload's.
ALWAYS_RETAINED_COMPONENTS = frozenset(
    {UNCLASSIFIED_COMPONENT, "workspace-config", "workspace-inputs", "assets"}
)


STORAGE_RESTORE_COMPONENTS = (
    "inputs", "raw-output", "analysis", "visualization", "logs", "assets",
    UNCLASSIFIED_COMPONENT, "workspace-config", "workspace-inputs",
)


#: Chunks a partial restore (`--checkpoint`/`--component`) always includes, whatever was
#: selected: a run/study-case archive's own identity ("metadata"), or a workspace
#: archive's own identity ("workspace-config"). Mirrors ALWAYS_RETAINED_COMPONENTS on
#: the offload side - identity and config evidence stay available either way.
ALWAYS_RESTORED_COMPONENTS = frozenset({"metadata", "workspace-config"})


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
