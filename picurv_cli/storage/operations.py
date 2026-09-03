"""!
@file operations.py
@brief Plan, protect, offload, restore, prune, and the command surface.
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
    ALWAYS_RESTORED_COMPONENTS,
    CLI_OUTPUT_FORMATS,
    DEFAULT_CHUNK_SIZE_GIB,
    DEFAULT_PROFILE_NAME,
    DEFAULT_STORAGE_WORKERS,
    KNOWN_CHECKPOINT_VERSION,
    REMOTE_COMPLETE_FILENAME,
    REMOTE_MANIFEST_FILENAME,
    REMOTE_OBJECTS_DIRECTORY,
    STORAGE_COMPRESSION_POLICIES,
    STORAGE_CONFIG_FILENAME,
    STORAGE_OFFLOAD_POLICIES,
    STORAGE_RESTORE_COMPONENTS,
    STORAGE_SCHEMA_VERSION,
    StorageError,
    UNCLASSIFIED_COMPONENT,
    WORKSPACE_CONFIG_FILENAME,
    _PARALLEL_GZIP_VALUES,
    _atomic_write_json,
    _find_upwards,
    _human_bytes,
    _parse_tags,
    _sha256_file,
    _state_path,
    _utc_now,
    load_storage_profile,
    read_storage_state,
    resolve_storage_config_path,
)
from .transport import (
    _blob_remote,
    _chunk_remote_path,
    _load_remote_manifest,
    _object_remote,
    _remote_blob_present,
    _remote_join,
)
from .safety import (
    _assert_archive_safe,
    _require_free_space,
    storage_operation_lock,
)
from .inventory import (
    _inventory_fingerprint,
    _select_completed_case_ids,
    inspect_artifact,
    resolve_local_storage_targets,
)
from .packaging import (
    _build_chunk_specs,
    _chunk_extension,
    _compression_size_range,
    _entry_retained_by_policy,
    _extract_chunk,
    _merge_tree,
    _resolve_offload_policy,
    _safe_component_name,
    _select_compression,
    _write_tar_chunk,
)
from .compatibility import (
    _capture_parameter_summary,
    _capture_run_assets,
    _capture_study_context,
    _capture_workspace_assets,
    _config_fingerprints,
    _git_provenance,
    _rebase_restored_text_paths,
    _restore_study_context,
)
from .catalog import (
    _find_reusable_archive,
    list_remote_manifests,
    prune_unused_workspace_assets,
    resolve_workspace_archive_id,
    verify_remote_archive,
)
from . import transport as _transport


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
        # As many chunks as `plan["workers"]` can be staged (packaged, not yet
        # uploaded) at once; each is deleted right after its verified upload, so
        # staging never needs archive-sized space, only room for the largest
        # concurrent handful of uncompressed chunks.
        concurrent_chunks = sorted(
            (spec["uncompressed_bytes"] for spec in specs), reverse=True
        )[:plan["workers"]]
        _require_free_space(
            staging_parent or tempfile.gettempdir(), sum(concurrent_chunks),
            f"stage {target['root_path']}",
        )
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
                        verified = _transport._upload_verified(local_chunk, _blob_remote(profile, digest))
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
            _transport._upload_verified(manifest_path, _object_remote(profile, archive_id, REMOTE_MANIFEST_FILENAME))
            manifest_digest = _sha256_file(manifest_path)
            complete_path = os.path.join(staging, REMOTE_COMPLETE_FILENAME)
            with open(complete_path, "w", encoding="ascii") as stream:
                stream.write(manifest_digest + "\n")
            _transport._upload_verified(complete_path, _object_remote(profile, archive_id, REMOTE_COMPLETE_FILENAME))

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
    # Each chunk is extracted in full before `_merge_tree` copies it into `destination`,
    # so both copies exist on disk at once.
    extracted_bytes = sum(chunk.get("uncompressed_bytes", 0) for chunk in chunks)
    _require_free_space(
        destination, 2 * extracted_bytes, f"restore archive {manifest.get('archive_id')} components"
    )
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
            _transport._run_rclone([
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
            if component in ALWAYS_RESTORED_COMPONENTS:
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
    # Every selected chunk is extracted in full before `_merge_tree` copies it into
    # `materialized`, so both copies exist on disk at once ahead of the final move.
    extracted_bytes = sum(chunk.get("uncompressed_bytes", 0) for chunk in chunks)
    _require_free_space(parent, 2 * extracted_bytes, f"restore archive {archive_id}")
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
            _transport._run_rclone([
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
    _transport._run_rclone(["mkdir", _remote_join(profile["remote"], REMOTE_OBJECTS_DIRECTORY)])
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
        if getattr(args, "completed", False):
            case_ids = _select_completed_case_ids(study_root)
        else:
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
