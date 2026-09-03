"""!
@file catalog.py
@brief The remote archive catalog and workspace asset references.
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
    ARCHIVE_ID_PATTERN,
    REMOTE_MANIFEST_FILENAME,
    REMOTE_OBJECTS_DIRECTORY,
    StorageError,
    is_artifact_cold,
)
from .transport import (
    _chunk_remote_path,
    _load_remote_manifest,
    _remote_join,
    _remote_sha256,
)
from . import transport as _transport


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


def list_remote_manifests(profile: dict) -> list:
    """!
    @brief Enumerate completed archive manifests from the remote catalog.
    @param[in] profile Value supplied through the `profile` argument.
    @return Result produced by this operation.
    """
    root = _remote_join(profile["remote"], REMOTE_OBJECTS_DIRECTORY)
    result = _transport._run_rclone([
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
