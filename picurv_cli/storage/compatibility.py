"""!
@file compatibility.py
@brief Provenance capture and path rebasing across restores.
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
    _sha256_file,
)
from .inventory import (
    _artifact_runtime_roots,
)


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
