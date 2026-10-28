"""!
@file transport.py
@brief The rclone process boundary and remote object addressing.
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
import threading
import uuid
from pathlib import Path
import yaml
from .models import (
    ARCHIVE_ID_PATTERN,
    REMOTE_BLOBS_DIRECTORY,
    REMOTE_COMPLETE_FILENAME,
    REMOTE_MANIFEST_FILENAME,
    REMOTE_OBJECTS_DIRECTORY,
    STORAGE_SCHEMA_VERSION,
    StorageError,
    _sha256_file,
)


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


#: Remote blob directories this process has created, and the lock that serializes creating them.
_CREATED_BLOB_DIRECTORIES = set()
_BLOB_DIRECTORY_LOCK = threading.Lock()


def _ensure_blob_directory(profile: dict, digest: str) -> None:
    """!
    @brief Create the directory a blob is stored in, once, before any upload into it.
    @details Google Drive allows several folders with one name. Parallel uploads that each
             create a missing directory leave duplicates, and a later lookup by path sees
             only one of them, so a stored blob reads as missing. Creating each directory
             once, serially, gives every upload the same folder.
    @param[in] profile Active storage profile.
    @param[in] digest Payload SHA-256 whose blob directory is needed.
    """
    directory = _remote_join(profile["remote"], REMOTE_BLOBS_DIRECTORY, digest[:2])
    with _BLOB_DIRECTORY_LOCK:
        if directory not in _CREATED_BLOB_DIRECTORIES:
            _run_rclone(["mkdir", directory])
            _CREATED_BLOB_DIRECTORIES.add(directory)


def _chunk_remote_path(profile: dict, archive_id: str, chunk: dict) -> str:
    """!
    @brief Resolve where one manifest chunk's payload lives on the remote.

    @details Every chunk this schema writes is blob-addressed, resolved out of the
             shared content-addressed store. `_load_remote_manifest` already refuses
             any manifest whose schema version does not match the current one, so a
             non-blob chunk cannot reach this function from a manifest this build
             actually loads.
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


def _rclone_error_detail(text: str) -> str:
    """!
    @brief The lines of rclone output that report the failure.
    @details rclone prints advisories, such as a shared client ID's deprecation, as NOTICE
             lines ahead of the error, and reports some failures as NOTICE lines too
             ("Failed to cat"), so the lines kept are those naming a failure or error.
    @param[in] text rclone stderr or stdout.
    @return The failure lines, else the non-NOTICE lines, else the whole text.
    """
    lines = [line for line in text.strip().splitlines() if line.strip()]
    failures = [line for line in lines if re.search(r"\bERROR\b|\bFailed\b|\berror\b", line)]
    other = [line for line in lines if " NOTICE: " not in line]
    return "\n".join(failures or other or lines)


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
        detail = _rclone_error_detail(result.stderr or result.stdout or "unknown rclone error")
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


#: rclone exit codes for a missing directory or file.
_RCLONE_NOT_FOUND_CODES = (3, 4)


def _read_remote_bytes(remote_path: str) -> bytes:
    """!
    @brief Read a small remote catalog object through rclone.
    @param[in] remote_path Value supplied through the `remote_path` argument.
    @return Result produced by this operation.
    @throws FileNotFoundError when the object does not exist, so callers can name what is missing.
    """
    executable = shutil.which("rclone")
    if not executable:
        raise StorageError("rclone was not found on PATH.")
    result = subprocess.run(
        [executable, "cat", remote_path], capture_output=True, check=False
    )
    if result.returncode in _RCLONE_NOT_FOUND_CODES:
        raise FileNotFoundError(remote_path)
    if result.returncode != 0:
        detail = _rclone_error_detail(
            (result.stderr or result.stdout or b"unknown rclone error").decode("utf-8", "replace"))
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
    try:
        manifest_bytes = _read_remote_bytes(_object_remote(profile, archive_id, REMOTE_MANIFEST_FILENAME))
    except FileNotFoundError:
        raise StorageError(
            f"No archive {archive_id} on {profile['remote']}; 'picurv storage list' shows the archives there."
        ) from None
    if require_complete:
        try:
            complete = _read_remote_bytes(_object_remote(profile, archive_id, REMOTE_COMPLETE_FILENAME))
        except FileNotFoundError:
            raise StorageError(f"Archive {archive_id} has no valid completion marker.") from None
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
