"""!
@file safety.py
@brief Refusals: active jobs, runtime locks, and cold-payload requirements.
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
    STORAGE_LOCK_FILENAME,
    StorageError,
    _read_json,
    _utc_now,
    is_artifact_cold,
    read_storage_state,
)


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
