"""!
@file packaging.py
@brief Compression policy, component chunks, and archive extraction.
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
    STORAGE_RETENTION_COMPONENTS,
    ALWAYS_RETAINED_COMPONENTS,
    AUTO_MAXIMUM_COMPRESSION_BYTES,
    AUTO_NO_COMPRESSION_BYTES,
    STORAGE_COMPRESSION_EXTENSIONS,
    STORAGE_COMPRESSION_POLICIES,
    STORAGE_OFFLOAD_POLICIES,
    StorageError,
    _PARALLEL_GZIP_VALUES,
)


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


def _normalize_retention_selection(values, flag: str) -> set:
    """!
    @brief Validate one --retain/--drop selection into a component name set.
    @param[in] values Repeated component names from the command line, or None.
    @param[in] flag User-facing flag name, for error messages.
    @return Set of selectable component names.
    """
    selected = set()
    for raw in values or []:
        for token in str(raw).split(","):
            name = token.strip()
            if not name:
                continue
            if name not in STORAGE_RETENTION_COMPONENTS:
                raise StorageError(
                    f"{flag} {name!r} is not a selectable component. Choose from: "
                    + ", ".join(STORAGE_RETENTION_COMPONENTS)
                )
            selected.add(name)
    return selected


def _resolve_offload_policy(profile: dict, requested: str = None,
                            keep_latest_checkpoint=None, retain=None, drop=None) -> dict:
    """!
    @brief Resolve semantic local-retention behavior for an offload.

    @details A named policy is a preset, not a ceiling. `--retain`/`--drop` adjust the
             preset one component at a time, so a campaign that needs, say, the
             analysis of `analysis-ready` but also its checkpoints does not have to
             choose the policy that happens to bundle both. The preset is still what
             decides everything not named explicitly.
    @param[in] profile Active storage profile.
    @param[in] requested Optional command-level policy override.
    @param[in] keep_latest_checkpoint Optional checkpoint-retention override.
    @param[in] retain Component names to retain locally regardless of the preset.
    @param[in] drop Component names to prune locally regardless of the preset.
    @return Normalized retention policy mapping.
    """
    name = requested or profile.get("offload_policy", "metadata-only")
    if name not in STORAGE_OFFLOAD_POLICIES:
        raise StorageError("Offload policy must be one of: " + ", ".join(STORAGE_OFFLOAD_POLICIES))
    retain_set = _normalize_retention_selection(retain, "--retain")
    drop_set = _normalize_retention_selection(drop, "--drop")
    conflicting = retain_set & drop_set
    if conflicting:
        raise StorageError(
            "A component cannot be both retained and dropped: " + ", ".join(sorted(conflicting))
        )
    # The conflict is between two things the user actually asked for, so it is decided
    # on the requested value: the derived default is False far more often than
    # --drop-all-checkpoints was typed, and refusing on that would reject the ordinary
    # `--retain checkpoints` with no second flag at all.
    if "checkpoints" in retain_set and keep_latest_checkpoint is False:
        raise StorageError(
            "--retain checkpoints keeps every committed step, which contradicts "
            "--drop-all-checkpoints. Pass one or the other."
        )
    if keep_latest_checkpoint is None:
        # Nothing explicit was requested: fall back to the profile default, or to
        # `restart-ready`'s own promise to keep the newest checkpoint. An explicit
        # --keep-latest-checkpoint/--drop-all-checkpoints always overrides both,
        # rather than being silently overruled by the policy it was paired with.
        keep_latest_checkpoint = bool(profile.get("keep_latest_checkpoint", False)) or name == "restart-ready"
    retained = set({
        "metadata-only": {"metadata", "logs"},
        "restart-ready": {"metadata", "logs", "inputs"},
        "analysis-ready": {"metadata", "logs", "analysis", "visualization"},
    }[name])
    retained |= retain_set
    retained -= drop_set
    # Identity is never optional: a cold artifact that cannot say what it is, what it
    # ran, and what it consumed cannot be restored or reasoned about.
    retained.add("metadata")
    return {
        "name": name,
        "retained_components": sorted(retained),
        "explicit_retain": sorted(retain_set),
        "explicit_drop": sorted(drop_set),
        "keep_latest_checkpoint": bool(keep_latest_checkpoint),
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
    retained = set(policy["retained_components"])
    if component.startswith("checkpoint:"):
        # "checkpoints" retains every committed step; --keep-latest-checkpoint is the
        # narrower selection of just the newest one.
        if "checkpoints" in retained:
            return True
        return bool(
            policy["keep_latest_checkpoint"]
            and latest_step is not None
            and component == f"checkpoint:{latest_step}"
        )
    return component in retained


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
