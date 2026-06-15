#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""!
@file core.py
@brief A comprehensive conductor script for the PICurv simulation platform.

This script acts as the central user interface for running simulations,
managing configurations, and orchestrating the entire end-to-end workflow.
It translates user-friendly YAML files into C-solver compatible control files,
supports full multi-block configurations, and provides live log streaming.
It features intelligent, content-based config file discovery and robustly
manages data I/O paths for the post-processor. It also supports Slurm job
generation/submission and parameter sweeps via job arrays.
"""

import yaml
import sys
import os
import argparse
import subprocess
import shutil
import glob
import csv
import json
import hashlib
import itertools
import re
import shlex
import copy
import math
import difflib
from datetime import datetime
import time
import filecmp

_NUMPY_MODULE = None
_MATPLOTLIB_PYPLOT = None


class _LazyNumpyProxy:
    """!
    @brief Module-like proxy that preserves `picurv.np` without eager import.
    """

    def __getattr__(self, name):
        """!
        @brief Resolve a NumPy attribute on first use.
        @param[in] name NumPy attribute name.
        @return Requested NumPy attribute.
        """
        return getattr(require_numpy(), name)


np = _LazyNumpyProxy()


def _prune_incompatible_python_site_paths(paths):
    """!
    @brief Remove site-package paths for a different Python major/minor version.
    @param[in] paths Candidate sys.path entries.
    @return Filtered path list.
    """
    current = (sys.version_info[0], sys.version_info[1])
    pattern = re.compile(r"python(?:-)?(\d+)\.(\d+)", re.IGNORECASE)
    filtered = []
    for path in paths:
        text = str(path)
        match = pattern.search(text)
        if match:
            path_version = (int(match.group(1)), int(match.group(2)))
            if path_version != current and ("site-packages" in text or "dist-packages" in text):
                continue
        filtered.append(path)
    return filtered


def _drop_imported_package(package_name: str):
    """!
    @brief Remove a failed/partial import package tree from sys.modules.
    @param[in] package_name Top-level package name.
    """
    prefix = package_name + "."
    for module_name in list(sys.modules):
        if module_name == package_name or module_name.startswith(prefix):
            sys.modules.pop(module_name, None)


def require_numpy():
    """!
    @brief Import NumPy only for commands that need numeric reductions.
    @return Imported NumPy module.
    """
    global _NUMPY_MODULE
    if _NUMPY_MODULE is not None:
        return _NUMPY_MODULE
    try:
        import numpy
    except Exception as exc:
        first_error = exc
        original_path = list(sys.path)
        try:
            _drop_imported_package("numpy")
            sys.path = _prune_incompatible_python_site_paths(original_path)
            import numpy
        except Exception as retry_exc:
            raise RuntimeError(
                "NumPy is required for this operation, but no compatible NumPy "
                "could be imported for this Python interpreter. PICurv ignored "
                "site-packages paths for other Python versions and retried. "
                f"First error: {first_error}. Retry error: {retry_exc}"
            ) from retry_exc
        finally:
            sys.path = original_path
    _NUMPY_MODULE = numpy
    return _NUMPY_MODULE


def optional_matplotlib_pyplot():
    """!
    @brief Import matplotlib.pyplot lazily for study plot generation.
    @return matplotlib.pyplot when available, otherwise None.
    """
    global _MATPLOTLIB_PYPLOT
    if _MATPLOTLIB_PYPLOT is not None:
        return _MATPLOTLIB_PYPLOT
    original_path = list(sys.path)
    try:
        import matplotlib.pyplot as pyplot
    except Exception:
        try:
            _drop_imported_package("matplotlib")
            sys.path = _prune_incompatible_python_site_paths(original_path)
            import matplotlib.pyplot as pyplot
        except Exception:
            return None
        finally:
            sys.path = original_path
    _MATPLOTLIB_PYPLOT = pyplot
    return _MATPLOTLIB_PYPLOT

# --- Global Path Definitions ---
# The implementation package and source-tree entrypoint live in picurv_cli/,
# while generators/ owns standalone generators.
PACKAGE_PATH = os.path.dirname(os.path.realpath(__file__))
PACKAGE_PROJECT_ROOT = os.path.dirname(PACKAGE_PATH)
INVOKED_SCRIPT_DIR = os.environ.get(
    "_PICURV_INVOKED_SCRIPT_DIR",
    PACKAGE_PATH,
)
SCRIPT_PATH = os.environ.get(
    "_PICURV_SCRIPT_PATH",
    PACKAGE_PATH,
)
PROJECT_ROOT = os.path.dirname(SCRIPT_PATH)
GENERATORS_PATH = os.path.join(PACKAGE_PROJECT_ROOT, "generators")
if os.path.basename(SCRIPT_PATH) == "bin":
    DEFAULT_BIN_DIR = SCRIPT_PATH
else:
    DEFAULT_BIN_DIR = os.path.join(PROJECT_ROOT, "bin")

PICURV_VERSION = "0.1.0"
CASE_ORIGIN_METADATA_FILENAME = ".picurv-origin.json"
RUNTIME_EXECUTION_CONFIG_FILENAME = ".picurv-execution.yml"
LEGACY_LOCAL_RUNTIME_CONFIG_FILENAME = ".picurv-local.yml"
RUNTIME_EXECUTION_EXAMPLE_FILENAME = "execution.example.yml"
RUNTIME_EXECUTION_CONFIG_FILENAMES = (
    RUNTIME_EXECUTION_CONFIG_FILENAME,
    LEGACY_LOCAL_RUNTIME_CONFIG_FILENAME,
)

DEFAULT_RUNTIME_EXECUTION_CONFIG_TEMPLATE = """# Optional shared runtime launcher overrides.
# This file is safe to leave unchanged on ordinary local machines.
# Edit it only when your site needs custom MPI launcher tokens.
#
# Precedence:
#   - local/login-node runs: local_execution -> default_execution -> built-in mpiexec
#   - generated cluster jobs: cluster.yml.execution -> cluster_execution -> default_execution -> built-in srun
#
# Example override:
# default_execution:
#   launcher: "mpirun"
#   launcher_args:
#     - --bind-to
#     - none
default_execution: {}

local_execution: {}

cluster_execution: {}
"""

CLUSTER_TEMPLATE_PLACEHOLDER_ACCOUNT = "my_project_account"
CLUSTER_TEMPLATE_PLACEHOLDER_MAIL = "user@example.edu"

DEFAULT_WALLTIME_GUARD_POLICY = {
    "enabled": True,
    "warmup_steps": 10,
    "multiplier": 2.0,
    "min_seconds": 60.0,
    "estimator_alpha": 0.35,
}
WALLTIME_GUARD_ENV_JOB_START_EPOCH = "PICURV_JOB_START_EPOCH"
WALLTIME_GUARD_ENV_LIMIT_SECONDS = "PICURV_WALLTIME_LIMIT_SECONDS"
POST_RESUME_STATE_FILENAME = "post.resume.json"
POST_LOCK_FILENAME = "post.lock"
POST_LOCK_METADATA_FILENAME = "post.lock.json"
POST_LOCK_WRAPPER_FILENAME = "post_lock_wrapper.py"
POST_RESUME_SCHEMA_VERSION = 1
POST_RECIPE_SIGNATURE_EXCLUDED_KEYS = {"startTime", "endTime"}
POST_REQUIRED_EULERIAN_SOURCE_BASENAMES = ("ufield", "vfield", "pfield", "nvfield")


def parse_slurm_time_limit_to_seconds(time_text: str) -> int:
    """!
    @brief Parse a Slurm time-limit string into total seconds.
    @param[in] time_text Argument passed to `parse_slurm_time_limit_to_seconds()`.
    @return Value returned by `parse_slurm_time_limit_to_seconds()`.
    """
    text = str(time_text).strip()
    if not text:
        raise ValueError("time limit cannot be empty")

    days = 0
    clock_text = text
    if "-" in text:
        day_text, clock_text = text.split("-", 1)
        if not day_text.isdigit():
            raise ValueError(f"invalid day field '{day_text}'")
        days = int(day_text)
        if not clock_text:
            raise ValueError("missing time portion after day field")

    parts = clock_text.split(":")
    if len(parts) > 3:
        raise ValueError(f"unsupported time format '{time_text}'")
    if any(part == "" for part in parts):
        raise ValueError(f"malformed time field '{time_text}'")
    if any(not part.isdigit() for part in parts):
        raise ValueError(f"non-numeric time field '{time_text}'")

    nums = [int(part) for part in parts]
    if days > 0:
        if len(nums) == 1:
            hours, minutes, seconds = nums[0], 0, 0
        elif len(nums) == 2:
            hours, minutes = nums
            seconds = 0
        else:
            hours, minutes, seconds = nums
    else:
        if len(nums) == 1:
            hours, minutes, seconds = 0, nums[0], 0
        elif len(nums) == 2:
            hours = 0
            minutes, seconds = nums
        else:
            hours, minutes, seconds = nums

    if minutes >= 60 or seconds >= 60:
        raise ValueError(f"minutes and seconds must be < 60 in '{time_text}'")
    if days == 0 and len(nums) == 3 and hours < 0:
        raise ValueError(f"hours must be non-negative in '{time_text}'")

    total_seconds = (((days * 24) + hours) * 60 + minutes) * 60 + seconds
    if total_seconds <= 0:
        raise ValueError("time limit must be positive")
    return total_seconds


def resolve_walltime_guard_policy(cluster_cfg: "dict | None") -> "dict | None":
    """!
    @brief Resolve the effective Slurm walltime-guard policy for generated solver jobs.
    @param[in] cluster_cfg Argument passed to `resolve_walltime_guard_policy()`.
    @return Value returned by `resolve_walltime_guard_policy()`.
    """
    if not isinstance(cluster_cfg, dict):
        return None

    scheduler = cluster_cfg.get("scheduler", {}) or {}
    if str(scheduler.get("type", "slurm")).lower() != "slurm":
        return None

    execution = cluster_cfg.get("execution", {}) or {}
    guard_cfg = execution.get("walltime_guard")
    if guard_cfg is None:
        guard_cfg = {}
    elif not isinstance(guard_cfg, dict):
        raise ValueError("execution.walltime_guard must be a mapping when provided")

    policy = copy.deepcopy(DEFAULT_WALLTIME_GUARD_POLICY)
    policy.update(guard_cfg)
    policy["enabled"] = bool(policy["enabled"])
    policy["warmup_steps"] = int(policy["warmup_steps"])
    policy["multiplier"] = float(policy["multiplier"])
    policy["min_seconds"] = float(policy["min_seconds"])
    policy["estimator_alpha"] = float(policy["estimator_alpha"])
    return policy


def build_walltime_guard_exports(cluster_cfg: "dict | None") -> dict:
    """!
    @brief Build shell-evaluated environment exports for the runtime walltime guard.
    @param[in] cluster_cfg Argument passed to `build_walltime_guard_exports()`.
    @return Value returned by `build_walltime_guard_exports()`.
    """
    policy = resolve_walltime_guard_policy(cluster_cfg)
    if not policy or not policy.get("enabled", False):
        return {}
    walltime_limit_seconds = parse_slurm_time_limit_to_seconds(cluster_cfg.get("resources", {}).get("time", ""))
    return {
        WALLTIME_GUARD_ENV_JOB_START_EPOCH: "$(date +%s)",
        WALLTIME_GUARD_ENV_LIMIT_SECONDS: str(walltime_limit_seconds),
    }

def resolve_runtime_executable(executable_name: str) -> str:
    """!
    @brief Resolve solver/post executable path, preferring local sibling binaries.
    @param[in] executable_name Argument passed to `resolve_runtime_executable()`.
    @return Value returned by `resolve_runtime_executable()`.
    """
    local_candidate = os.path.join(INVOKED_SCRIPT_DIR, executable_name)
    if os.path.isfile(local_candidate):
        return os.path.abspath(local_candidate)
    return os.path.join(DEFAULT_BIN_DIR, executable_name)

# Standardized error codes used for CLI/validation reporting.
ERROR_CODE_CLI_USAGE_INVALID = "CLI_USAGE_INVALID"
ERROR_CODE_CFG_MISSING_SECTION = "CFG_MISSING_SECTION"
ERROR_CODE_CFG_MISSING_KEY = "CFG_MISSING_KEY"
ERROR_CODE_CFG_INVALID_TYPE = "CFG_INVALID_TYPE"
ERROR_CODE_CFG_INVALID_VALUE = "CFG_INVALID_VALUE"
ERROR_CODE_CFG_FILE_NOT_FOUND = "CFG_FILE_NOT_FOUND"
ERROR_CODE_CFG_GRID_PARSE = "CFG_GRID_PARSE"
ERROR_CODE_CFG_INCONSISTENT_COMBO = "CFG_INCONSISTENT_COMBO"
ERROR_CODE_DEPENDENCY_MISSING = "DEPENDENCY_MISSING"

_ERROR_HINTS = {
    ERROR_CODE_CLI_USAGE_INVALID: "Run 'picurv <command> --help' to see valid argument combinations.",
    ERROR_CODE_CFG_MISSING_SECTION: "Add the missing section using examples/master_template/*.yml as reference.",
    ERROR_CODE_CFG_MISSING_KEY: "Add the missing key in the referenced YAML file.",
    ERROR_CODE_CFG_INVALID_TYPE: "Fix the value type to match the documented schema in docs/pages/14_Config_Contract.md.",
    ERROR_CODE_CFG_INVALID_VALUE: "Adjust the value to a supported range/enum from the config reference pages.",
    ERROR_CODE_CFG_FILE_NOT_FOUND: "Fix the path or create the missing file before running again.",
    ERROR_CODE_CFG_GRID_PARSE: "Validate grid file format and numeric payload (block count, dims, coordinates).",
    ERROR_CODE_CFG_INCONSISTENT_COMBO: "Fix conflicting options/keys so the configuration is internally consistent.",
    ERROR_CODE_DEPENDENCY_MISSING: "Install the named optional dependency for the Python interpreter used by picurv.",
}


def _sanitize_error_field(value) -> str:
    """!
    @brief Normalize error fields into a single-line string.
    @param[in] value Argument passed to `_sanitize_error_field()`.
    @return Value returned by `_sanitize_error_field()`.
    """
    if value is None:
        return "-"
    text = str(value).strip()
    if not text:
        return "-"
    return " ".join(text.splitlines())


def emit_structured_error(code: str, key: str = "-", file_path: str = "-",
                          message: str = "", hint: str = None, stream=None):
    """!
    @brief Emit one standardized error line for tooling and users.
    @param[in] code Argument passed to `emit_structured_error()`.
    @param[in] key Argument passed to `emit_structured_error()`.
    @param[in] file_path Argument passed to `emit_structured_error()`.
    @param[in] message Argument passed to `emit_structured_error()`.
    @param[in] hint Argument passed to `emit_structured_error()`.
    @param[in] stream Argument passed to `emit_structured_error()`.
    """
    if stream is None:
        stream = sys.stderr
    resolved_hint = hint if hint is not None else _ERROR_HINTS.get(code, "-")
    print(
        f"ERROR {_sanitize_error_field(code)} | "
        f"key={_sanitize_error_field(key)} | "
        f"file={_sanitize_error_field(file_path)} | "
        f"message={_sanitize_error_field(message)} | "
        f"hint={_sanitize_error_field(resolved_hint)}",
        file=stream,
    )


def fail_cli_usage(message: str, hint: str = None):
    """!
    @brief Emit a structured CLI usage error and exit with code 2.
    @param[in] message Argument passed to `fail_cli_usage()`.
    @param[in] hint Argument passed to `fail_cli_usage()`.
    """
    emit_structured_error(
        ERROR_CODE_CLI_USAGE_INVALID,
        key="-",
        file_path="-",
        message=message,
        hint=hint or _ERROR_HINTS[ERROR_CODE_CLI_USAGE_INVALID],
    )
    sys.exit(2)


def _split_error_file_and_message(raw_error: str):
    """!
    @brief Split '<file>: <message>' style validation strings when possible.
    @param[in] raw_error Argument passed to `_split_error_file_and_message()`.
    @return Value returned by `_split_error_file_and_message()`.
    """
    text = str(raw_error).strip()
    match = re.match(r"^(?P<file>[^:]+):\s*(?P<msg>.+)$", text)
    if not match:
        return "-", text
    file_candidate = match.group("file").strip()
    msg = match.group("msg").strip()
    known_suffixes = (".yml", ".yaml", ".cfg", ".picgrid", ".control", ".run", ".txt")
    if "/" in file_candidate or file_candidate.endswith(known_suffixes):
        return file_candidate, msg
    return "-", text


def _extract_key_path(message: str) -> str:
    """!
    @brief Best-effort key-path extraction from free-form validation messages.
    @param[in] message Argument passed to `_extract_key_path()`.
    @return Value returned by `_extract_key_path()`.
    """
    dotted = re.search(r"\b([A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z0-9_\[\]-]+)+)\b", message)
    if dotted:
        return dotted.group(1)

    bracketed = re.search(r"\b([A-Za-z_][A-Za-z0-9_]*\[[^\]]+\](?:\[[^\]]+\])*)\b", message)
    if bracketed:
        return bracketed.group(1)

    quoted = re.findall(r"'([A-Za-z0-9_.\[\]-]+)'", message)
    for token in quoted:
        if "." in token or "[" in token or token.isidentifier():
            return token
    return "-"


def _classify_error_code(message: str) -> str:
    """!
    @brief Map existing validation/error messages to the standardized code set.
    @param[in] message Argument passed to `_classify_error_code()`.
    @return Value returned by `_classify_error_code()`.
    """
    msg = message.lower()
    if "missing required section" in msg:
        return ERROR_CODE_CFG_MISSING_SECTION
    if "missing required key" in msg or "missing key" in msg:
        return ERROR_CODE_CFG_MISSING_KEY
    if "not found" in msg or "does not exist" in msg:
        return ERROR_CODE_CFG_FILE_NOT_FOUND
    if "invalid dimensions line" in msg or "invalid coordinate row" in msg or "grid file" in msg:
        return ERROR_CODE_CFG_GRID_PARSE
    if (
        "must both be periodic" in msg
        or "inconsistent periodicity" in msg
        or "mismatch" in msg
        or "requires --" in msg
        or "must be 1 (auto) or exactly" in msg
    ):
        return ERROR_CODE_CFG_INCONSISTENT_COMBO
    if (
        "must be a mapping" in msg
        or "must be a list" in msg
        or "must be a string" in msg
        or "must be a boolean" in msg
        or "must be either" in msg
    ):
        return ERROR_CODE_CFG_INVALID_TYPE
    if "unsupported key" in msg or "unsupported top-level section" in msg:
        return ERROR_CODE_CFG_INVALID_VALUE
    return ERROR_CODE_CFG_INVALID_VALUE

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def read_yaml_file(filepath: str) -> dict:
    """!
    @brief Safely reads a YAML file and returns its content.
    @param[in] filepath Path to the YAML file.
    @return A dictionary containing the parsed YAML content.
    @throws SystemExit if the file is not found or cannot be parsed.
    """
    if not os.path.exists(filepath):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key="-",
            file_path=filepath,
            message="Configuration file not found.",
        )
        sys.exit(1)
    try:
        with open(filepath, 'r') as f:
            return yaml.safe_load(f)
    except yaml.YAMLError as e:
        emit_structured_error(
            ERROR_CODE_CFG_INVALID_VALUE,
            key="-",
            file_path=filepath,
            message=f"YAML parse error: {e}",
            hint="Fix YAML syntax/indentation and retry validation.",
        )
        sys.exit(1)

def write_yaml_file(filepath: str, data: dict):
    """!
    @brief Write YAML with stable ordering for generated study artifacts.
    @param[in] filepath Argument passed to `write_yaml_file()`.
    @param[in] data Argument passed to `write_yaml_file()`.
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)

def write_json_file(filepath: str, payload: dict):
    """!
    @brief Write JSON metadata/manifests with a stable, readable format.
    @param[in] filepath Argument passed to `write_json_file()`.
    @param[in] payload Argument passed to `write_json_file()`.
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def write_runtime_execution_file(filepath: str, template_source_path: str = None) -> str:
    """!
    @brief Write a default runtime execution config, copying a source template when available.
    @param[in] filepath Argument passed to `write_runtime_execution_file()`.
    @param[in] template_source_path Argument passed to `write_runtime_execution_file()`.
    @return Value returned by `write_runtime_execution_file()`.
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    if template_source_path and os.path.isfile(template_source_path):
        shutil.copy2(template_source_path, filepath)
        return filepath

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(DEFAULT_RUNTIME_EXECUTION_CONFIG_TEMPLATE)
    return filepath


def _launcher_arg_contains_whitespace(token) -> bool:
    """!
    @brief Return True when a launcher arg token contains embedded whitespace and should be split.
    @param[in] token Argument passed to `_launcher_arg_contains_whitespace()`.
    @return Value returned by `_launcher_arg_contains_whitespace()`.
    """
    return isinstance(token, str) and any(ch.isspace() for ch in token.strip())


def resolve_runtime_execution_seed_source(source_project_root: str) -> "str | None":
    """!
    @brief Prefer repo-local ignored runtime config, then tracked example, then built-in defaults.
    @param[in] source_project_root Argument passed to `resolve_runtime_execution_seed_source()`.
    @return Value returned by `resolve_runtime_execution_seed_source()`.
    """
    source_root_abs = os.path.abspath(source_project_root)
    repo_local_runtime = os.path.join(source_root_abs, RUNTIME_EXECUTION_CONFIG_FILENAME)
    if os.path.isfile(repo_local_runtime):
        return repo_local_runtime

    tracked_example = os.path.join(
        source_root_abs,
        "config",
        "runtime",
        RUNTIME_EXECUTION_EXAMPLE_FILENAME,
    )
    if os.path.isfile(tracked_example):
        return tracked_example
    return None


def ensure_case_runtime_execution_config(case_dir: str, source_project_root: str, overwrite: bool = False) -> dict:
    """!
    @brief Create case-local runtime execution config if missing, seeded from repo-local config when available.
    @param[in] case_dir Argument passed to `ensure_case_runtime_execution_config()`.
    @param[in] source_project_root Argument passed to `ensure_case_runtime_execution_config()`.
    @param[in] overwrite Argument passed to `ensure_case_runtime_execution_config()`.
    @return Value returned by `ensure_case_runtime_execution_config()`.
    """
    case_dir_abs = os.path.abspath(case_dir)
    dest_path = os.path.join(case_dir_abs, RUNTIME_EXECUTION_CONFIG_FILENAME)
    if os.path.exists(dest_path) and not overwrite:
        return {
            "path": dest_path,
            "created": False,
            "seed_source": None,
        }

    seed_source = resolve_runtime_execution_seed_source(source_project_root)
    write_runtime_execution_file(dest_path, seed_source)
    return {
        "path": dest_path,
        "created": True,
        "seed_source": seed_source,
    }


def is_project_root(candidate: str) -> bool:
    """!
    @brief Return True when a directory looks like the PICurv source repository root.
    @param[in] candidate Argument passed to `is_project_root()`.
    @return Value returned by `is_project_root()`.
    """
    if not candidate:
        return False
    candidate_abs = os.path.abspath(candidate)
    return (
        os.path.isfile(os.path.join(candidate_abs, "Makefile"))
        and os.path.isdir(os.path.join(candidate_abs, "src"))
        and os.path.isdir(os.path.join(candidate_abs, "include"))
        and os.path.isdir(os.path.join(candidate_abs, "picurv_cli"))
    )


def _iter_parent_dirs(start_path: str):
    """!
    @brief Yield a path and all of its parents up to filesystem root.
    @param[in] start_path Argument passed to `_iter_parent_dirs()`.
    """
    current = os.path.abspath(start_path)
    if os.path.isfile(current):
        current = os.path.dirname(current)
    while True:
        yield current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent


def find_project_root_upwards(start_path: str):
    """!
    @brief Search upward from an anchor and return the first matching project root.
    @param[in] start_path Argument passed to `find_project_root_upwards()`.
    @return Value returned by `find_project_root_upwards()`.
    """
    if not start_path:
        return None
    for directory in _iter_parent_dirs(start_path):
        if is_project_root(directory):
            return directory
    return None


def discover_local_project_root(*extra_anchors):
    """!
    @brief Best-effort source repo discovery from runtime anchors.
    @param[in] extra_anchors Argument passed to `discover_local_project_root()`.
    @return Value returned by `discover_local_project_root()`.
    """
    anchors = list(extra_anchors) + [os.getcwd(), INVOKED_SCRIPT_DIR, SCRIPT_PATH, PROJECT_ROOT]
    seen = set()
    for anchor in anchors:
        if not anchor:
            continue
        anchor_abs = os.path.abspath(anchor)
        if anchor_abs in seen:
            continue
        seen.add(anchor_abs)
        project_root = find_project_root_upwards(anchor_abs)
        if project_root:
            return project_root
    return None


def find_case_origin_metadata_file(case_dir_hint: str = None):
    """!
    @brief Find the nearest case-origin metadata file from known runtime anchors.
    @param[in] case_dir_hint Argument passed to `find_case_origin_metadata_file()`.
    @return Value returned by `find_case_origin_metadata_file()`.
    """
    search_roots = []
    for candidate in (case_dir_hint, os.getcwd(), INVOKED_SCRIPT_DIR):
        if not candidate:
            continue
        abs_candidate = os.path.abspath(candidate)
        if abs_candidate not in search_roots:
            search_roots.append(abs_candidate)

    for root in search_roots:
        for directory in _iter_parent_dirs(root):
            metadata_path = os.path.join(directory, CASE_ORIGIN_METADATA_FILENAME)
            if os.path.isfile(metadata_path):
                return metadata_path
    return None


def load_case_origin_metadata(case_dir_hint: str = None):
    """!
    @brief Load case-origin metadata if present, returning (case_dir, metadata_path, payload).
    @param[in] case_dir_hint Argument passed to `load_case_origin_metadata()`.
    @return Value returned by `load_case_origin_metadata()`.
    """
    metadata_path = find_case_origin_metadata_file(case_dir_hint)
    if not metadata_path:
        return None, None, None
    try:
        with open(metadata_path, "r", encoding="utf-8") as f:
            payload = json.load(f)
        if not isinstance(payload, dict):
            raise ValueError("Case origin metadata must be a JSON object.")
    except Exception as exc:
        raise ValueError(f"Failed to read case origin metadata at {metadata_path}: {exc}") from exc
    return os.path.dirname(metadata_path), metadata_path, payload


def find_runtime_execution_config_file(*anchors):
    """!
    @brief Find the nearest optional execution config from runtime/case anchors.
    @param[in] anchors Argument passed to `find_runtime_execution_config_file()`.
    @return Value returned by `find_runtime_execution_config_file()`.
    """
    seen = set()
    search_roots = []
    for candidate in list(anchors) + [os.getcwd(), INVOKED_SCRIPT_DIR]:
        if not candidate:
            continue
        current = os.path.abspath(candidate)
        if os.path.isfile(current):
            current = os.path.dirname(current)
        if current in seen:
            continue
        seen.add(current)
        search_roots.append(current)

    seen_dirs = set()
    for root in search_roots:
        for directory in _iter_parent_dirs(root):
            if directory in seen_dirs:
                continue
            seen_dirs.add(directory)
            for filename in RUNTIME_EXECUTION_CONFIG_FILENAMES:
                config_path = os.path.join(directory, filename)
                if os.path.isfile(config_path):
                    return config_path
    return None


def _normalize_execution_override_section(payload: dict, section_name: str, config_path: str, config_label: str) -> dict:
    """!
    @brief Validate one execution override section while preserving missing-vs-empty semantics.
    @param[in] payload Argument passed to `_normalize_execution_override_section()`.
    @param[in] section_name Argument passed to `_normalize_execution_override_section()`.
    @param[in] config_path Argument passed to `_normalize_execution_override_section()`.
    @param[in] config_label Argument passed to `_normalize_execution_override_section()`.
    @return Value returned by `_normalize_execution_override_section()`.
    """
    section = payload.get(section_name)
    if section is None:
        return {"launcher": None, "launcher_args": None}
    if not isinstance(section, dict):
        raise ValueError(f"{config_label} at {config_path}: {section_name} must be a mapping.")

    launcher = section.get("launcher")
    if launcher is not None and not isinstance(launcher, str):
        raise ValueError(f"{config_label} at {config_path}: {section_name}.launcher must be a string.")

    launcher_args = None
    if "launcher_args" in section:
        launcher_args = section.get("launcher_args", [])
        if launcher_args is None:
            launcher_args = []
        if not isinstance(launcher_args, list):
            raise ValueError(f"{config_label} at {config_path}: {section_name}.launcher_args must be a list.")
        for i, token in enumerate(launcher_args):
            if not isinstance(token, (str, int, float)):
                raise ValueError(
                    f"{config_label} at {config_path}: {section_name}.launcher_args[{i}] "
                    "must be a scalar CLI token."
                )
            if _launcher_arg_contains_whitespace(token):
                raise ValueError(
                    f"{config_label} at {config_path}: {section_name}.launcher_args[{i}] "
                    "must be a single CLI token; split whitespace-separated arguments into separate list items."
                )
        launcher_args = [str(x) for x in launcher_args]

    return {
        "launcher": launcher,
        "launcher_args": launcher_args,
    }


def load_runtime_execution_config(config_search_anchor: str = None, extra_search_anchors=None):
    """!
    @brief Load optional shared execution launcher config from the nearest runtime config file.
    @param[in] config_search_anchor Argument passed to `load_runtime_execution_config()`.
    @param[in] extra_search_anchors Argument passed to `load_runtime_execution_config()`.
    @return Value returned by `load_runtime_execution_config()`.
    """
    anchors = []
    if config_search_anchor is not None:
        anchors.append(config_search_anchor)
    if extra_search_anchors:
        anchors.extend(extra_search_anchors)

    config_path = find_runtime_execution_config_file(*anchors)
    if not config_path:
        return None, {}

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            payload = yaml.safe_load(f) or {}
    except yaml.YAMLError as exc:
        raise ValueError(f"{os.path.basename(config_path)} YAML parse error at {config_path}: {exc}") from exc

    if not isinstance(payload, dict):
        raise ValueError(f"{os.path.basename(config_path)} at {config_path} must be a YAML mapping.")

    return config_path, {
        "default_execution": _normalize_execution_override_section(
            payload,
            "default_execution",
            config_path,
            os.path.basename(config_path),
        ),
        "local_execution": _normalize_execution_override_section(
            payload,
            "local_execution",
            config_path,
            os.path.basename(config_path),
        ),
        "cluster_execution": _normalize_execution_override_section(
            payload,
            "cluster_execution",
            config_path,
            os.path.basename(config_path),
        ),
    }


def merge_execution_overrides(base: "dict | None", override: "dict | None") -> dict:
    """!
    @brief Merge execution overrides, letting explicit override values win key-by-key.
    @param[in] base Argument passed to `merge_execution_overrides()`.
    @param[in] override Argument passed to `merge_execution_overrides()`.
    @return Value returned by `merge_execution_overrides()`.
    """
    base = base or {}
    override = override or {}

    launcher = override.get("launcher")
    if launcher is None:
        launcher = base.get("launcher")

    launcher_args = override.get("launcher_args")
    if launcher_args is None:
        launcher_args = base.get("launcher_args")

    return {
        "launcher": launcher,
        "launcher_args": None if launcher_args is None else [str(x) for x in launcher_args],
    }


def resolve_runtime_execution_context(runtime_execution_cfg: dict, context: str) -> dict:
    """!
    @brief Resolve default plus context-specific execution overrides.
    @param[in] runtime_execution_cfg Argument passed to `resolve_runtime_execution_context()`.
    @param[in] context Argument passed to `resolve_runtime_execution_context()`.
    @return Value returned by `resolve_runtime_execution_context()`.
    """
    if context not in {"local", "cluster"}:
        raise ValueError(f"Unsupported execution context '{context}'.")
    return merge_execution_overrides(
        runtime_execution_cfg.get("default_execution"),
        runtime_execution_cfg.get(f"{context}_execution"),
    )


def get_git_commit(repo_root: str = None) -> str:
    """!
    @brief Best-effort git commit lookup for run/study manifests and case metadata.
    @param[in] repo_root Argument passed to `get_git_commit()`.
    @return Value returned by `get_git_commit()`.
    """
    cwd = repo_root or PROJECT_ROOT
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=cwd,
            text=True,
            capture_output=True,
            check=False
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return None


def write_case_origin_metadata(case_dir: str, source_project_root: str, template_name: str = None,
                               existing: dict = None, template_managed_files=None):
    """!
    @brief Create or refresh case-origin metadata for repo-aware case maintenance commands.
    @param[in] case_dir Argument passed to `write_case_origin_metadata()`.
    @param[in] source_project_root Argument passed to `write_case_origin_metadata()`.
    @param[in] template_name Argument passed to `write_case_origin_metadata()`.
    @param[in] existing Argument passed to `write_case_origin_metadata()`.
    @param[in] template_managed_files Argument passed to `write_case_origin_metadata()`.
    @return Value returned by `write_case_origin_metadata()`.
    """
    payload = dict(existing or {})
    if "initialized_at" not in payload:
        payload["initialized_at"] = datetime.now().isoformat()
    payload["source_repo_root"] = os.path.abspath(source_project_root)
    if template_name:
        payload["template_name"] = template_name
    if template_managed_files is not None:
        payload["template_managed_files"] = sorted(set(str(p) for p in template_managed_files))
    payload["last_known_source_git_commit"] = get_git_commit(source_project_root)
    metadata_path = os.path.join(os.path.abspath(case_dir), CASE_ORIGIN_METADATA_FILENAME)
    write_json_file(metadata_path, payload)
    return metadata_path, payload


def make_args_include_explicit_goal(make_args: "list[str]") -> bool:
    """!
    @brief Return True when make args contain an explicit target rather than only options/assignments.
    @param[in] make_args Argument passed to `make_args_include_explicit_goal()`.
    @return Value returned by `make_args_include_explicit_goal()`.
    """
    if not make_args:
        return False

    options_with_value = {
        "-C", "-f", "-I", "-j", "-l", "-o", "-W",
        "--directory", "--file", "--makefile", "--include-dir", "--jobs",
        "--load-average", "--max-load", "--old-file", "--assume-old",
        "--what-if", "--new-file", "--assume-new",
    }
    assignment_pattern = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*[:+?]?=.*$")

    skip_next = False
    for token in make_args:
        if skip_next:
            skip_next = False
            continue
        if token in options_with_value:
            skip_next = True
            continue
        if token.startswith("-"):
            continue
        if assignment_pattern.match(token):
            continue
        return True
    return False


def resolve_case_origin_context(case_dir_hint: str = None, source_root_override: str = None, template_name_override: str = None):
    """!
    @brief Resolve case directory, source repo root, and optional template metadata.
    @param[in] case_dir_hint Argument passed to `resolve_case_origin_context()`.
    @param[in] source_root_override Argument passed to `resolve_case_origin_context()`.
    @param[in] template_name_override Argument passed to `resolve_case_origin_context()`.
    @return Value returned by `resolve_case_origin_context()`.
    """
    metadata_case_dir, metadata_path, metadata = load_case_origin_metadata(case_dir_hint)

    if metadata_case_dir:
        case_dir = metadata_case_dir
    else:
        case_dir = os.path.abspath(case_dir_hint or os.getcwd())

    source_project_root = source_root_override
    if source_project_root:
        source_project_root = os.path.abspath(source_project_root)
    elif isinstance(metadata, dict) and isinstance(metadata.get("source_repo_root"), str):
        source_project_root = os.path.abspath(metadata["source_repo_root"])
    else:
        source_project_root = discover_local_project_root(case_dir)

    template_name = template_name_override
    if not template_name and isinstance(metadata, dict):
        template_name = metadata.get("template_name")

    return {
        "case_dir": case_dir,
        "metadata_path": metadata_path,
        "metadata": metadata or {},
        "source_project_root": source_project_root,
        "template_name": template_name,
    }


def require_project_root(candidate: str, purpose: str):
    """!
    @brief Validate that a source repo root was resolved and is structurally valid.
    @param[in] candidate Argument passed to `require_project_root()`.
    @param[in] purpose Argument passed to `require_project_root()`.
    @return Value returned by `require_project_root()`.
    """
    if not candidate:
        raise ValueError(
            f"Could not determine the PICurv source repository for {purpose}. "
            "Run this command from an initialized case directory or pass --source-root."
        )
    candidate_abs = os.path.abspath(candidate)
    if not is_project_root(candidate_abs):
        raise ValueError(
            f"Resolved source repository for {purpose} is not a valid PICurv root: {candidate_abs}"
        )
    return candidate_abs


def require_existing_case_dir(case_dir: str, purpose: str, source_project_root: str = None):
    """!
    @brief Validate that a target case directory exists and is not the source repo root.
    @param[in] case_dir Argument passed to `require_existing_case_dir()`.
    @param[in] purpose Argument passed to `require_existing_case_dir()`.
    @param[in] source_project_root Argument passed to `require_existing_case_dir()`.
    @return Value returned by `require_existing_case_dir()`.
    """
    if not case_dir:
        raise ValueError(f"Could not determine the case directory for {purpose}. Pass --case-dir.")
    case_dir_abs = os.path.abspath(case_dir)
    if not os.path.isdir(case_dir_abs):
        raise ValueError(f"Case directory for {purpose} does not exist: {case_dir_abs}")
    if source_project_root and os.path.abspath(source_project_root) == case_dir_abs:
        raise ValueError(
            f"Refusing to run {purpose} against the source repository root itself: {case_dir_abs}"
        )
    return case_dir_abs


def resolve_template_directory(source_project_root: str, template_name: str):
    """!
    @brief Resolve an example template directory inside the source repository.
    @param[in] source_project_root Argument passed to `resolve_template_directory()`.
    @param[in] template_name Argument passed to `resolve_template_directory()`.
    @return Value returned by `resolve_template_directory()`.
    """
    if not template_name:
        raise ValueError(
            "Template name is required for config sync. Re-run with --template-name or from a case initialized by current picurv."
        )
    template_dir = os.path.join(source_project_root, "examples", template_name)
    if not os.path.isdir(template_dir):
        raise ValueError(f"Case template '{template_name}' not found at '{template_dir}'")
    return template_dir


def list_template_relative_files(template_dir: str, excluded_rel_paths=None):
    """!
    @brief List all files in a template directory as case-relative paths.
    @param[in] template_dir Argument passed to `list_template_relative_files()`.
    @param[in] excluded_rel_paths Argument passed to `list_template_relative_files()`.
    @return Value returned by `list_template_relative_files()`.
    """
    template_dir_abs = os.path.abspath(template_dir)
    if not os.path.isdir(template_dir_abs):
        raise ValueError(f"Template directory not found: {template_dir_abs}")
    excluded = set(excluded_rel_paths or [])
    relative_paths = []
    for root, _, files in os.walk(template_dir_abs):
        rel_root = os.path.relpath(root, template_dir_abs)
        for filename in sorted(files):
            rel_path = filename if rel_root == "." else os.path.join(rel_root, filename)
            if rel_path in excluded:
                continue
            relative_paths.append(rel_path)
    return relative_paths


def list_source_binaries(source_project_root: str):
    """!
    @brief List binary artifacts currently available in the source repo bin directory.
    @param[in] source_project_root Argument passed to `list_source_binaries()`.
    @return Value returned by `list_source_binaries()`.
    """
    source_bin_dir = os.path.join(os.path.abspath(source_project_root), "bin")
    if not os.path.isdir(source_bin_dir):
        raise ValueError(f"Source bin directory not found: {source_bin_dir}. Run 'picurv build' first.")
    binaries = sorted(
        f for f in os.listdir(source_bin_dir)
        if os.path.isfile(os.path.join(source_bin_dir, f)) and f != "picurv"
    )
    if not binaries:
        raise ValueError(f"Source bin directory contains no files: {source_bin_dir}")
    return source_bin_dir, binaries


def sync_case_binaries(case_dir: str, source_project_root: str):
    """!
    @brief Copy current source-repo binaries into a case directory for version-pinning.
    @param[in] case_dir Argument passed to `sync_case_binaries()`.
    @param[in] source_project_root Argument passed to `sync_case_binaries()`.
    @return Value returned by `sync_case_binaries()`.
    """
    case_dir_abs = os.path.abspath(case_dir)
    os.makedirs(case_dir_abs, exist_ok=True)
    source_bin_dir, binaries = list_source_binaries(source_project_root)
    copied = []
    for binary_name in binaries:
        source_path = os.path.join(source_bin_dir, binary_name)
        dest_path = os.path.join(case_dir_abs, binary_name)
        shutil.copy2(source_path, dest_path)
        copied.append(dest_path)
    return copied


def sync_case_template_files(case_dir: str, template_dir: str, overwrite: bool = False,
                             prune: bool = False, managed_rel_paths=None):
    """!
    @brief Sync template files into a case directory, preserving modified files unless overwrite is requested.
    @param[in] case_dir Argument passed to `sync_case_template_files()`.
    @param[in] template_dir Argument passed to `sync_case_template_files()`.
    @param[in] overwrite Argument passed to `sync_case_template_files()`.
    @param[in] prune Argument passed to `sync_case_template_files()`.
    @param[in] managed_rel_paths Argument passed to `sync_case_template_files()`.
    @return Value returned by `sync_case_template_files()`.
    """
    case_dir_abs = os.path.abspath(case_dir)
    template_dir_abs = os.path.abspath(template_dir)
    if not os.path.isdir(template_dir_abs):
        raise ValueError(f"Template directory not found: {template_dir_abs}")

    summary = {
        "copied": [],
        "overwritten": [],
        "skipped_modified": [],
        "unchanged": [],
        "pruned": [],
        "prune_requested_without_tracking": False,
    }
    excluded_rel_paths = {RUNTIME_EXECUTION_EXAMPLE_FILENAME}
    current_template_files = list_template_relative_files(
        template_dir_abs,
        excluded_rel_paths=excluded_rel_paths,
    )
    current_template_set = set(current_template_files)

    for root, _, files in os.walk(template_dir_abs):
        rel_root = os.path.relpath(root, template_dir_abs)
        for filename in sorted(files):
            src_path = os.path.join(root, filename)
            rel_path = filename if rel_root == "." else os.path.join(rel_root, filename)
            if rel_path in excluded_rel_paths:
                continue
            dest_path = os.path.join(case_dir_abs, rel_path)
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)

            if not os.path.exists(dest_path):
                shutil.copy2(src_path, dest_path)
                summary["copied"].append(dest_path)
                continue

            if filecmp.cmp(src_path, dest_path, shallow=False):
                summary["unchanged"].append(dest_path)
                continue

            if overwrite:
                shutil.copy2(src_path, dest_path)
                summary["overwritten"].append(dest_path)
            else:
                summary["skipped_modified"].append(dest_path)

    managed_set = set(managed_rel_paths or [])
    if prune:
        if not managed_set:
            summary["prune_requested_without_tracking"] = True
        for rel_path in sorted(managed_set - current_template_set):
            dest_path = os.path.join(case_dir_abs, rel_path)
            if os.path.isfile(dest_path):
                os.remove(dest_path)
                summary["pruned"].append(dest_path)

    summary["template_managed_files"] = current_template_files
    return summary


def compute_case_source_status(case_dir: str, source_project_root: str, template_name: str = None, metadata: dict = None):
    """!
    @brief Compute source/case drift across commits, binaries, and template-managed files.
    @param[in] case_dir Argument passed to `compute_case_source_status()`.
    @param[in] source_project_root Argument passed to `compute_case_source_status()`.
    @param[in] template_name Argument passed to `compute_case_source_status()`.
    @param[in] metadata Argument passed to `compute_case_source_status()`.
    @return Value returned by `compute_case_source_status()`.
    """
    case_dir_abs = os.path.abspath(case_dir)
    source_root_abs = os.path.abspath(source_project_root)
    metadata = metadata or {}
    status = {
        "case_dir": case_dir_abs,
        "source_repo_root": source_root_abs,
        "metadata_present": bool(metadata),
        "template_name": template_name,
        "last_known_source_git_commit": metadata.get("last_known_source_git_commit"),
        "current_source_git_commit": get_git_commit(source_root_abs),
    }
    status["source_commit_changed"] = (
        bool(status["last_known_source_git_commit"])
        and bool(status["current_source_git_commit"])
        and status["last_known_source_git_commit"] != status["current_source_git_commit"]
    )

    binary_status = {
        "source_bin_present": False,
        "source_bin_missing": [],
        "case_bin_missing": [],
        "case_bin_different": [],
        "case_bin_current": [],
    }
    try:
        source_bin_dir, binaries = list_source_binaries(source_root_abs)
        binary_status["source_bin_present"] = True
        for binary_name in binaries:
            source_path = os.path.join(source_bin_dir, binary_name)
            case_path = os.path.join(case_dir_abs, binary_name)
            if not os.path.isfile(case_path):
                binary_status["case_bin_missing"].append(binary_name)
            elif filecmp.cmp(source_path, case_path, shallow=False):
                binary_status["case_bin_current"].append(binary_name)
            else:
                binary_status["case_bin_different"].append(binary_name)
    except ValueError as exc:
        binary_status["source_bin_missing"].append(str(exc))
    status["binaries"] = binary_status

    config_status = {
        "template_available": False,
        "template_files": [],
        "case_missing_files": [],
        "case_modified_files": [],
        "case_current_files": [],
        "template_removed_since_last_sync": [],
        "tracking_available": isinstance(metadata.get("template_managed_files"), list),
    }
    if template_name:
        try:
            template_dir = resolve_template_directory(source_root_abs, template_name)
            template_files = list_template_relative_files(
                template_dir,
                excluded_rel_paths={RUNTIME_EXECUTION_EXAMPLE_FILENAME},
            )
            config_status["template_available"] = True
            config_status["template_files"] = template_files
            for rel_path in template_files:
                src_path = os.path.join(template_dir, rel_path)
                case_path = os.path.join(case_dir_abs, rel_path)
                if not os.path.isfile(case_path):
                    config_status["case_missing_files"].append(rel_path)
                elif filecmp.cmp(src_path, case_path, shallow=False):
                    config_status["case_current_files"].append(rel_path)
                else:
                    config_status["case_modified_files"].append(rel_path)
            managed_files = metadata.get("template_managed_files")
            if isinstance(managed_files, list):
                config_status["template_removed_since_last_sync"] = sorted(set(managed_files) - set(template_files))
        except ValueError:
            pass
    status["config"] = config_status

    case_runtime_cfg = os.path.join(case_dir_abs, RUNTIME_EXECUTION_CONFIG_FILENAME)
    repo_runtime_seed = os.path.join(source_root_abs, RUNTIME_EXECUTION_CONFIG_FILENAME)
    runtime_status = {
        "case_config_present": os.path.isfile(case_runtime_cfg),
        "repo_seed_present": os.path.isfile(repo_runtime_seed),
        "case_matches_repo_seed": False,
    }
    if runtime_status["case_config_present"] and runtime_status["repo_seed_present"]:
        runtime_status["case_matches_repo_seed"] = filecmp.cmp(
            case_runtime_cfg,
            repo_runtime_seed,
            shallow=False,
        )
    status["runtime_execution"] = runtime_status
    return status


def print_case_source_status(status: dict):
    """!
    @brief Render human-readable source/case drift details.
    @param[in] status Argument passed to `print_case_source_status()`.
    """
    print(f"[INFO] Case directory      : {status['case_dir']}")
    print(f"[INFO] Source repo        : {status['source_repo_root']}")
    print(f"[INFO] Template           : {status.get('template_name') or '(unknown)'}")
    if status.get("last_known_source_git_commit"):
        print(f"[INFO] Last synced commit : {status['last_known_source_git_commit']}")
    if status.get("current_source_git_commit"):
        print(f"[INFO] Current src commit : {status['current_source_git_commit']}")
    print(f"[INFO] Source changed     : {'yes' if status.get('source_commit_changed') else 'no'}")

    binaries = status["binaries"]
    if binaries["source_bin_present"]:
        print(
            f"[INFO] Binaries          : current={len(binaries['case_bin_current'])} "
            f"changed={len(binaries['case_bin_different'])} missing={len(binaries['case_bin_missing'])}"
        )
    else:
        print("[INFO] Binaries          : source bin/ unavailable")

    config = status["config"]
    if config["template_available"]:
        print(
            f"[INFO] Template files    : current={len(config['case_current_files'])} "
            f"modified={len(config['case_modified_files'])} missing={len(config['case_missing_files'])}"
        )
        if config["tracking_available"]:
            print(f"[INFO] Prune candidates  : {len(config['template_removed_since_last_sync'])}")
        else:
            print("[INFO] Prune candidates  : tracking unavailable")
    elif status.get("template_name"):
        print("[INFO] Template files    : template unavailable in source repo")

    runtime_cfg = status.get("runtime_execution", {})
    print(
        f"[INFO] Runtime config    : case={'yes' if runtime_cfg.get('case_config_present') else 'no'} "
        f"repo-seed={'yes' if runtime_cfg.get('repo_seed_present') else 'no'} "
        f"matches-repo-seed={'yes' if runtime_cfg.get('case_matches_repo_seed') else 'no'}"
    )


def status_source_command(args):
    """!
    @brief Report source/case drift for an initialized case directory.
    @param[in] args Command-line style argument list supplied to the function.
    """
    try:
        context = resolve_case_origin_context(
            case_dir_hint=getattr(args, "case_dir", None),
            source_root_override=getattr(args, "source_root", None),
            template_name_override=getattr(args, "template_name", None),
        )
        source_project_root = require_project_root(context["source_project_root"], "status-source")
        case_dir = require_existing_case_dir(context["case_dir"], "status-source", source_project_root)
        status = compute_case_source_status(
            case_dir,
            source_project_root,
            template_name=context.get("template_name"),
            metadata=context.get("metadata"),
        )
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    if getattr(args, "output_format", "text") == "json":
        print(json.dumps(status, indent=2, sort_keys=True))
        return
    print_case_source_status(status)

def resolve_path(anchor_file: str, candidate: str) -> str:
    """!
    @brief Resolve a potentially relative path against a source YAML file path.
    @param[in] anchor_file Argument passed to `resolve_path()`.
    @param[in] candidate Argument passed to `resolve_path()`.
    @return Value returned by `resolve_path()`.
    """
    if candidate is None:
        return None
    if os.path.isabs(candidate):
        return os.path.abspath(candidate)
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(anchor_file)), candidate))


POST_RUN_CONTROL_ALIASES = {
    "start_step": ("start_step", "startTime"),
    "end_step": ("end_step", "endTime"),
    "step_interval": ("step_interval", "timeStep"),
}


GRID_GENERATOR_HYPHEN_KEY_HINTS = {
    "config-file": "config_file",
    "grid-type": "grid_type",
    "cli-args": "cli_args",
    "output-file": "output_file",
    "stats-file": "stats_file",
    "vts-file": "vts_file",
}


def _mapping_value_with_aliases(mapping: dict, *keys, default=None):
    """!
    @brief Return the first defined value from a mapping across alias keys.
    @param[in] mapping Argument passed to `_mapping_value_with_aliases()`.
    @param[in] default Argument passed to `_mapping_value_with_aliases()`.
    @param[in] keys Argument passed to `_mapping_value_with_aliases()`.
    @return Value returned by `_mapping_value_with_aliases()`.
    """
    if not isinstance(mapping, dict):
        return default
    for key in keys:
        if key in mapping:
            return mapping.get(key)
    return default


def get_post_run_control_value(post_cfg: dict, canonical_key: str, default=None):
    """!
    @brief Resolve post run_control values with backwards-compatible legacy aliases.
    @param[in] post_cfg Argument passed to `get_post_run_control_value()`.
    @param[in] canonical_key Argument passed to `get_post_run_control_value()`.
    @param[in] default Argument passed to `get_post_run_control_value()`.
    @return Value returned by `get_post_run_control_value()`.
    """
    aliases = POST_RUN_CONTROL_ALIASES.get(canonical_key, (canonical_key,))
    rc = post_cfg.get("run_control", {})
    return _mapping_value_with_aliases(rc, *aliases, default=default)


def warn_on_grid_generator_hyphen_keys(generator: dict, case_path: str, warnings: list) -> None:
    """!
    @brief Warn when grid.generator uses unsupported hyphenated wrapper keys.
    @param[in] generator grid.generator mapping from case.yml.
    @param[in] case_path Case file path for diagnostics.
    @param[in,out] warnings Warning list to append to.
    """
    if not isinstance(generator, dict):
        return
    for bad_key, expected_key in GRID_GENERATOR_HYPHEN_KEY_HINTS.items():
        if bad_key in generator and bad_key != expected_key:
            warnings.append(
                f"{case_path}: grid.generator.{bad_key} is ignored; use grid.generator.{expected_key}."
            )


def get_post_source_data(post_cfg: dict):
    """!
    @brief Return source_data as a mapping when valid, else an empty mapping.
    @param[in] post_cfg Argument passed to `get_post_source_data()`.
    @return Value returned by `get_post_source_data()`.
    """
    source_cfg = post_cfg.get("source_data", {})
    if isinstance(source_cfg, dict):
        return source_cfg
    return {}


def get_post_source_directory_template(post_cfg: dict, default: str = "<solver_output_dir>") -> str:
    """!
    @brief Resolve the source directory template from source_data with a safe default.
    @param[in] post_cfg Argument passed to `get_post_source_directory_template()`.
    @param[in] default Argument passed to `get_post_source_directory_template()`.
    @return Value returned by `get_post_source_directory_template()`.
    """
    return get_post_source_data(post_cfg).get("directory", default)


def get_post_input_extensions(post_cfg: dict):
    """!
    @brief Return post input_extensions, preferring io.* and tolerating legacy source_data.* placement.
    @param[in] post_cfg Argument passed to `get_post_input_extensions()`.
    @return Value returned by `get_post_input_extensions()`.
    """
    io_cfg = post_cfg.get("io", {})
    io_ext = io_cfg.get("input_extensions") if isinstance(io_cfg, dict) else None
    if isinstance(io_ext, dict):
        return io_ext

    source_ext = get_post_source_data(post_cfg).get("input_extensions")
    if isinstance(source_ext, dict):
        return source_ext

    return {}


def get_post_statistics_task_tokens(post_cfg: dict):
    """!
    @brief Return normalized statistics pipeline tokens that will be written into post.run.
    @param[in] post_cfg Argument passed to `get_post_statistics_task_tokens()`.
    @return Value returned by `get_post_statistics_task_tokens()`.
    """
    stats_cfg = post_cfg.get("statistics_pipeline")
    stats_entries = []
    if isinstance(stats_cfg, list):
        stats_entries = stats_cfg
    elif isinstance(stats_cfg, dict):
        stats_entries = stats_cfg.get("tasks", [])

    tokens = []
    for entry in stats_entries:
        if isinstance(entry, str):
            task_name = entry
        elif isinstance(entry, dict):
            task_name = entry.get("task")
        else:
            continue
        try:
            tokens.append(normalize_statistics_task(task_name))
        except ValueError:
            continue
    return tokens


def get_monitor_output_directory(monitor_cfg: dict, default: str = "output") -> str:
    """!
    @brief Resolve the solver output root from monitor.yml, preserving the default layout.
    @param[in] monitor_cfg Argument passed to `get_monitor_output_directory()`.
    @param[in] default Argument passed to `get_monitor_output_directory()`.
    @return Value returned by `get_monitor_output_directory()`.
    """
    io_cfg = monitor_cfg.get("io")
    if isinstance(io_cfg, dict):
        directories = io_cfg.get("directories")
        if isinstance(directories, dict):
            output_dir = directories.get("output")
            if isinstance(output_dir, str) and output_dir.strip():
                return output_dir.strip()

    return default


def get_post_statistics_output_prefix(post_cfg: dict, default: str = "Stats") -> str:
    """!
    @brief Resolve the statistics CSV prefix, preserving legacy top-level override support.
    @param[in] post_cfg Argument passed to `get_post_statistics_output_prefix()`.
    @param[in] default Argument passed to `get_post_statistics_output_prefix()`.
    @return Value returned by `get_post_statistics_output_prefix()`.
    """
    stats_cfg = post_cfg.get("statistics_pipeline")
    if isinstance(stats_cfg, dict):
        prefix = stats_cfg.get("output_prefix")
        if isinstance(prefix, str) and prefix.strip():
            return prefix.strip()

    legacy_prefix = post_cfg.get("statistics_output_prefix")
    if isinstance(legacy_prefix, str) and legacy_prefix.strip():
        return legacy_prefix.strip()

    return default


def resolve_post_statistics_output_prefix(post_cfg: dict, monitor_cfg=None, default: str = "Stats") -> str:
    """!
    @brief Resolve the runtime statistics prefix, routing bare basenames under the monitor output root.
    @param[in] post_cfg Argument passed to `resolve_post_statistics_output_prefix()`.
    @param[in] monitor_cfg Optional monitor configuration used to anchor the default statistics home.
    @param[in] default Argument passed to `resolve_post_statistics_output_prefix()`.
    @return Value returned by `resolve_post_statistics_output_prefix()`.
    """
    prefix = get_post_statistics_output_prefix(post_cfg, default=default).strip()
    if os.path.isabs(prefix):
        return prefix

    if os.path.dirname(prefix):
        return prefix

    monitor_output_dir = get_monitor_output_directory(monitor_cfg or {})
    return os.path.join(monitor_output_dir, "statistics", prefix)


def get_post_statistics_output_artifacts(post_cfg: dict, run_dir: str, monitor_cfg=None):
    """!
    @brief Predict statistics CSV output paths relative to the postprocessor runtime cwd.
    @param[in] post_cfg Argument passed to `get_post_statistics_output_artifacts()`.
    @param[in] run_dir Argument passed to `get_post_statistics_output_artifacts()`.
    @param[in] monitor_cfg Optional monitor configuration used to anchor the default statistics home.
    @return Value returned by `get_post_statistics_output_artifacts()`.
    """
    token_to_suffix = {
        "ComputeMSD": "_msd.csv",
    }
    prefix = resolve_post_statistics_output_prefix(post_cfg, monitor_cfg)
    if os.path.isabs(prefix):
        base_path = os.path.abspath(prefix)
    else:
        base_path = os.path.abspath(os.path.join(run_dir, prefix))

    output_paths = []
    for token in get_post_statistics_task_tokens(post_cfg):
        suffix = token_to_suffix.get(token)
        if suffix:
            output_paths.append(base_path + suffix)

    return list(dict.fromkeys(output_paths))


def build_post_recipe_config(post_cfg: dict, monitor_cfg=None) -> dict:
    """!
    @brief Build the flat key=value mapping consumed by the C post-processor.
    @param[in] post_cfg Argument passed to `build_post_recipe_config()`.
    @param[in] monitor_cfg Optional parsed monitor YAML configuration dictionary.
    @return Value returned by `build_post_recipe_config()`.
    """
    c_config = {}

    c_config['startTime'] = get_post_run_control_value(post_cfg, 'start_step', 0)
    c_config['endTime'] = get_post_run_control_value(post_cfg, 'end_step', 0)
    c_config['timeStep'] = get_post_run_control_value(post_cfg, 'step_interval', 1)

    eulerian_pipeline_parts = []
    if post_cfg.get('global_operations', {}).get('dimensionalize', False):
        eulerian_pipeline_parts.append('DimensionalizeAllLoadedFields')

    for task in post_cfg.get('eulerian_pipeline', []):
        task_name = task.get('task')
        if task_name == 'q_criterion':
            eulerian_pipeline_parts.append('ComputeQCriterion')
        elif task_name == 'normalize_field':
            field = task.get('field', 'P')
            eulerian_pipeline_parts.append(f'NormalizeRelativeField:{field}')
            ref_point = task.get('reference_point', [1, 1, 1])
            c_config['reference_ip'] = ref_point[0]
            c_config['reference_jp'] = ref_point[1]
            c_config['reference_kp'] = ref_point[2]
        elif task_name == 'nodal_average':
            in_field = task.get('input_field')
            out_field = task.get('output_field')
            if in_field and out_field:
                eulerian_pipeline_parts.append(f'CellToNodeAverage:{in_field}>{out_field}')

    if eulerian_pipeline_parts:
        c_config['process_pipeline'] = ";".join(eulerian_pipeline_parts)

    lagrangian_pipeline_parts = []
    for task in post_cfg.get('lagrangian_pipeline', []):
        task_name = task.get('task')
        if task_name == 'specific_ke':
            in_field = task.get('input_field')
            out_field = task.get('output_field')
            if in_field and out_field:
                lagrangian_pipeline_parts.append(f'ComputeSpecificKE:{in_field}>{out_field}')
    if lagrangian_pipeline_parts:
        c_config['particle_pipeline'] = ";".join(lagrangian_pipeline_parts)

    statistics_pipeline_parts = get_post_statistics_task_tokens(post_cfg)
    statistics_output_prefix = None
    stats_cfg = post_cfg.get('statistics_pipeline')
    if isinstance(stats_cfg, dict):
        statistics_output_prefix = stats_cfg.get('output_prefix')

    if statistics_pipeline_parts:
        c_config['statistics_pipeline'] = ";".join(statistics_pipeline_parts)
        statistics_output_prefix = resolve_post_statistics_output_prefix(post_cfg, monitor_cfg)
    elif statistics_output_prefix is None:
        statistics_output_prefix = post_cfg.get('statistics_output_prefix')
    if statistics_output_prefix:
        c_config['statistics_output_prefix'] = statistics_output_prefix

    io = post_cfg.get('io', {})
    c_config['output_prefix'] = io.get('output_directory', 'viz') + '/' + io.get('output_filename_prefix', 'Field')
    c_config['particle_output_prefix'] = io.get('output_directory', 'viz') + '/' + io.get('particle_filename_prefix', 'Particle')
    c_config['output_particles'] = io.get('output_particles', False)
    c_config['particle_output_freq'] = io.get('particle_subsampling_frequency', 1)
    c_config['output_fields_instantaneous'] = ",".join(io.get('eulerian_fields', []))
    c_config['output_fields_averaged'] = ",".join(io.get('eulerian_fields_averaged', []))
    c_config['particle_fields_instantaneous'] = ",".join(io.get('particle_fields', []))
    input_extensions = get_post_input_extensions(post_cfg)
    if isinstance(input_extensions, dict):
        e_ext = input_extensions.get('eulerian')
        p_ext = input_extensions.get('particle')
        if e_ext:
            c_config['eulerianExt'] = str(e_ext).strip().lstrip('.')
        if p_ext:
            c_config['particleExt'] = str(p_ext).strip().lstrip('.')

    source_directory = get_post_source_directory_template(post_cfg, default=None)
    if source_directory is not None:
        c_config['source_directory'] = source_directory

    return c_config


def normalize_post_recipe_signature(recipe_cfg: dict) -> dict:
    """!
    @brief Normalize post recipe settings into a stable signature mapping.
    @param[in] recipe_cfg Argument passed to `normalize_post_recipe_signature()`.
    @return Value returned by `normalize_post_recipe_signature()`.
    """
    signature = {}
    for key, value in (recipe_cfg or {}).items():
        if key in POST_RECIPE_SIGNATURE_EXCLUDED_KEYS or value is None:
            continue
        if isinstance(value, bool):
            text = 'true' if value else 'false'
        else:
            text = str(value).strip()
            if text.lower() in {'true', 'false'}:
                text = text.lower()
        if text:
            signature[str(key)] = text
    return signature


def compute_post_recipe_fingerprint(recipe_cfg: dict) -> "tuple[dict, str]":
    """!
    @brief Return normalized recipe signature plus SHA-256 fingerprint.
    @param[in] recipe_cfg Argument passed to `compute_post_recipe_fingerprint()`.
    @return Value returned by `compute_post_recipe_fingerprint()`.
    """
    signature = normalize_post_recipe_signature(recipe_cfg)
    payload = json.dumps(signature, sort_keys=True, separators=(',', ':')).encode('utf-8')
    return signature, hashlib.sha256(payload).hexdigest()


def parse_post_recipe_file(post_recipe_path: str):
    """!
    @brief Parse an existing generated post.run file into a key/value mapping.
    @param[in] post_recipe_path Argument passed to `parse_post_recipe_file()`.
    @return Value returned by `parse_post_recipe_file()`.
    """
    if not post_recipe_path or not os.path.isfile(post_recipe_path):
        return None
    recipe_cfg = {}
    with open(post_recipe_path, 'r', encoding='utf-8', errors='replace') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            key, value = line.split('=', 1)
            recipe_cfg[key.strip()] = value.strip()
    return recipe_cfg


def get_post_resume_state_path(run_dir: str) -> str:
    """!
    @brief Return the JSON resume metadata path for a run directory.
    @param[in] run_dir Argument passed to `get_post_resume_state_path()`.
    @return Value returned by `get_post_resume_state_path()`.
    """
    return os.path.join(run_dir, 'config', POST_RESUME_STATE_FILENAME)


def get_post_lock_paths(run_dir: str) -> dict:
    """!
    @brief Return lock-wrapper related paths for a run directory.
    @param[in] run_dir Argument passed to `get_post_lock_paths()`.
    @return Value returned by `get_post_lock_paths()`.
    """
    scheduler_dir = os.path.join(run_dir, 'scheduler')
    return {
        'lock_file': os.path.join(scheduler_dir, POST_LOCK_FILENAME),
        'metadata_file': os.path.join(scheduler_dir, POST_LOCK_METADATA_FILENAME),
        'wrapper_path': os.path.join(scheduler_dir, POST_LOCK_WRAPPER_FILENAME),
    }


def _post_output_directory_abs(run_dir: str, post_cfg: dict) -> str:
    """!
    @brief Resolve the absolute post output directory for the current recipe.
    @param[in] run_dir Argument passed to `_post_output_directory_abs()`.
    @param[in] post_cfg Argument passed to `_post_output_directory_abs()`.
    @return Value returned by `_post_output_directory_abs()`.
    """
    io_cfg = post_cfg.get('io', {}) or {}
    return os.path.abspath(os.path.join(run_dir, io_cfg.get('output_directory', 'viz')))


def _post_requests_eulerian_output(post_cfg: dict) -> bool:
    """!
    @brief Return whether the current post recipe expects Eulerian VTK output artifacts.
    @param[in] post_cfg Argument passed to `_post_requests_eulerian_output()`.
    @return Value returned by `_post_requests_eulerian_output()`.
    """
    io_cfg = post_cfg.get('io', {}) or {}
    return bool(io_cfg.get('eulerian_fields'))


def _post_requests_particle_output(post_cfg: dict) -> bool:
    """!
    @brief Return whether the current post recipe expects particle VTP output artifacts.
    @param[in] post_cfg Argument passed to `_post_requests_particle_output()`.
    @return Value returned by `_post_requests_particle_output()`.
    """
    io_cfg = post_cfg.get('io', {}) or {}
    return bool(io_cfg.get('output_particles')) and bool(io_cfg.get('particle_fields'))


def _post_requests_statistics(post_cfg: dict) -> bool:
    """!
    @brief Return whether the current post recipe expects statistics CSV artifacts.
    @param[in] post_cfg Argument passed to `_post_requests_statistics()`.
    @return Value returned by `_post_requests_statistics()`.
    """
    return bool(get_post_statistics_task_tokens(post_cfg))


def _post_needs_particle_source(post_cfg: dict) -> bool:
    """!
    @brief Return whether the current post recipe requires particle source files to be present.
    @param[in] post_cfg Argument passed to `_post_needs_particle_source()`.
    @return Value returned by `_post_needs_particle_source()`.
    """
    io_cfg = post_cfg.get('io', {}) or {}
    return bool(io_cfg.get('output_particles')) or bool(post_cfg.get('lagrangian_pipeline')) or _post_requests_statistics(post_cfg)


def _iter_post_steps(start_step: int, end_step: int, step_interval: int):
    """!
    @brief Yield configured post-processing steps inclusively.
    @param[in] start_step Argument passed to `_iter_post_steps()`.
    @param[in] end_step Argument passed to `_iter_post_steps()`.
    @param[in] step_interval Argument passed to `_iter_post_steps()`.
    """
    if step_interval <= 0 or end_step < start_step:
        return
    step = start_step
    while step <= end_step:
        yield step
        step += step_interval


def resolve_post_requested_window(post_cfg: dict, case_cfg: dict = None) -> "tuple[int, int, int]":
    """!
    @brief Resolve post requested start/end/interval, expanding end=-1 via case.yml when available.
    @param[in] post_cfg Argument passed to `resolve_post_requested_window()`.
    @param[in] case_cfg Optional case configuration for end-step expansion.
    @return Value returned by `resolve_post_requested_window()`.
    """
    start_step = int(get_post_run_control_value(post_cfg, 'start_step', 0) or 0)
    end_step = int(get_post_run_control_value(post_cfg, 'end_step', 0) or 0)
    step_interval = int(get_post_run_control_value(post_cfg, 'step_interval', 1) or 1)
    if end_step < 0 and case_cfg:
        case_run = case_cfg.get('run_control', {}) or {}
        case_start = int(case_run.get('start_step', 0) or 0)
        case_total = int(case_run.get('total_steps', 0) or 0)
        end_step = case_start + case_total
    return start_step, end_step, step_interval


def prepare_effective_post_config(post_cfg: dict, resolved_source_dir: str, start_step: int = None, end_step: int = None) -> dict:
    """!
    @brief Return a copy of post_cfg with resolved source dir and optional effective bounds.
    @param[in] post_cfg Argument passed to `prepare_effective_post_config()`.
    @param[in] resolved_source_dir Argument passed to `prepare_effective_post_config()`.
    @param[in] start_step Argument passed to `prepare_effective_post_config()`.
    @param[in] end_step Argument passed to `prepare_effective_post_config()`.
    @return Value returned by `prepare_effective_post_config()`.
    """
    effective_cfg = copy.deepcopy(post_cfg)
    if not isinstance(effective_cfg.get('source_data'), dict):
        effective_cfg['source_data'] = {}
    effective_cfg['source_data']['directory'] = resolved_source_dir
    rc = effective_cfg.setdefault('run_control', {})
    if start_step is not None:
        rc['start_step'] = int(start_step)
    if end_step is not None:
        rc['end_step'] = int(end_step)
    return effective_cfg


def _scan_post_vtk_steps(prefix_path: str, extension: str) -> "set[int]":
    """!
    @brief Scan VTK output files matching '<prefix>_<step>.<extension>'.
    @param[in] prefix_path Argument passed to `_scan_post_vtk_steps()`.
    @param[in] extension Argument passed to `_scan_post_vtk_steps()`.
    @return Value returned by `_scan_post_vtk_steps()`.
    """
    directory = os.path.dirname(prefix_path)
    if not os.path.isdir(directory):
        return set()
    basename = os.path.basename(prefix_path)
    pattern = re.compile(rf'^{re.escape(basename)}_(\d+)\.{re.escape(extension)}$')
    steps = set()
    for name in os.listdir(directory):
        match = pattern.match(name)
        if match:
            steps.add(int(match.group(1)))
    return steps


def _scan_post_statistics_csv_steps(csv_path: str) -> "set[int]":
    """!
    @brief Scan step ids from the first CSV column of a statistics artifact.
    @param[in] csv_path Argument passed to `_scan_post_statistics_csv_steps()`.
    @return Value returned by `_scan_post_statistics_csv_steps()`.
    """
    if not os.path.isfile(csv_path):
        return set()
    steps = set()
    with open(csv_path, 'r', encoding='utf-8', errors='replace', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            head = str(row[0]).strip().lower()
            if head in {'step', 'timestep', 'time_step'}:
                continue
            step_val = _parse_int_loose(row[0])
            if step_val is not None:
                steps.add(step_val)
    return steps


def collect_post_completion_families(run_dir: str, post_cfg: dict, monitor_cfg=None) -> "list[set[int]]":
    """!
    @brief Collect per-family completed-step sets for the current post recipe.
    @param[in] run_dir Argument passed to `collect_post_completion_families()`.
    @param[in] post_cfg Argument passed to `collect_post_completion_families()`.
    @param[in] monitor_cfg Optional parsed monitor YAML configuration dictionary.
    @return Value returned by `collect_post_completion_families()`.
    """
    io_cfg = post_cfg.get('io', {}) or {}
    output_dir_abs = _post_output_directory_abs(run_dir, post_cfg)
    families = []

    if _post_requests_eulerian_output(post_cfg):
        prefix = os.path.join(output_dir_abs, io_cfg.get('output_filename_prefix', 'Field'))
        families.append(_scan_post_vtk_steps(prefix, 'vts'))

    if _post_requests_particle_output(post_cfg):
        prefix = os.path.join(output_dir_abs, io_cfg.get('particle_filename_prefix', 'Particle'))
        families.append(_scan_post_vtk_steps(prefix, 'vtp'))

    for stats_path in get_post_statistics_output_artifacts(post_cfg, run_dir, monitor_cfg):
        families.append(_scan_post_statistics_csv_steps(stats_path))

    return families


def detect_post_completed_frontier(run_dir: str, post_cfg: dict, monitor_cfg, start_step: int, end_step: int, step_interval: int) -> dict:
    """!
    @brief Detect the highest contiguous fully completed post step for the current recipe.
    @param[in] run_dir Argument passed to `detect_post_completed_frontier()`.
    @param[in] post_cfg Argument passed to `detect_post_completed_frontier()`.
    @param[in] monitor_cfg Argument passed to `detect_post_completed_frontier()`.
    @param[in] start_step Argument passed to `detect_post_completed_frontier()`.
    @param[in] end_step Argument passed to `detect_post_completed_frontier()`.
    @param[in] step_interval Argument passed to `detect_post_completed_frontier()`.
    @return Value returned by `detect_post_completed_frontier()`.
    """
    families = collect_post_completion_families(run_dir, post_cfg, monitor_cfg)
    frontier = None
    if families:
        for step in _iter_post_steps(start_step, end_step, step_interval):
            if all(step in family for family in families):
                frontier = step
            else:
                break
    return {
        'frontier_step': frontier,
        'artifact_family_count': len(families),
    }


def _nearest_step(steps: "set[int]", target: int):
    """!
    @brief Return the complete source step nearest to a target step.
    @param[in] steps Candidate step numbers.
    @param[in] target Target step number.
    @return Nearest candidate, or None when no candidates exist.
    """
    if not steps:
        return None
    return min(steps, key=lambda step: (abs(step - target), step))


def _format_optional_step(step) -> str:
    """!
    @brief Format an optional step number for user-facing diagnostics.
    @param[in] step Step number or None.
    @return Printable step text.
    """
    return 'none' if step is None else str(step)


def _scan_complete_source_steps(source_dir: str, monitor_cfg: dict, post_cfg: dict) -> "tuple[set[int], dict]":
    """!
    @brief Scan source artifacts and return steps with every file required by the recipe.
    @param[in] source_dir Source output root directory.
    @param[in] monitor_cfg Parsed monitor configuration.
    @param[in] post_cfg Parsed post-processing configuration.
    @return Tuple of complete source steps and source path metadata.
    """
    dirs = (monitor_cfg.get('io', {}) or {}).get('directories', {}) or {}
    euler_subdir = dirs.get('eulerian_subdir', 'eulerian')
    particle_subdir = dirs.get('particle_subdir', 'particles')
    euler_dir = os.path.join(source_dir, euler_subdir)
    particle_dir = os.path.join(source_dir, particle_subdir)

    input_extensions = get_post_input_extensions(post_cfg)
    euler_ext = str((input_extensions.get('eulerian') or 'dat')).strip().lstrip('.')
    particle_ext = str((input_extensions.get('particle') or 'dat')).strip().lstrip('.')

    families = []
    for basename in POST_REQUIRED_EULERIAN_SOURCE_BASENAMES:
        steps = set()
        if os.path.isdir(euler_dir):
            pattern = re.compile(rf'^{re.escape(basename)}(\d{{5}})_0\.{re.escape(euler_ext)}$')
            for name in os.listdir(euler_dir):
                match = pattern.match(name)
                if match:
                    steps.add(int(match.group(1)))
        families.append(steps)

    if _post_needs_particle_source(post_cfg):
        steps = set()
        if os.path.isdir(particle_dir):
            pattern = re.compile(rf'^position(\d{{5}})_0\.{re.escape(particle_ext)}$')
            for name in os.listdir(particle_dir):
                match = pattern.match(name)
                if match:
                    steps.add(int(match.group(1)))
        families.append(steps)

    complete_steps = set.intersection(*families) if families else set()
    return complete_steps, {
        'euler_dir': euler_dir,
        'particle_dir': particle_dir,
        'euler_ext': euler_ext,
        'particle_ext': particle_ext,
    }


def _expected_source_paths_for_step(step: int, source_scan: dict, post_cfg: dict) -> "list[str]":
    """!
    @brief Build required source file paths for a single post-processing step.
    @param[in] step Requested step number.
    @param[in] source_scan Metadata returned by `_scan_complete_source_steps()`.
    @param[in] post_cfg Parsed post-processing configuration.
    @return Required source artifact paths.
    """
    paths = [
        os.path.join(source_scan['euler_dir'], f'{basename}{step:05d}_0.{source_scan["euler_ext"]}')
        for basename in POST_REQUIRED_EULERIAN_SOURCE_BASENAMES
    ]
    if _post_needs_particle_source(post_cfg):
        paths.append(os.path.join(source_scan['particle_dir'], f'position{step:05d}_0.{source_scan["particle_ext"]}'))
    return paths


def detect_post_source_frontier(source_dir: str, monitor_cfg: dict, post_cfg: dict, start_step: int, end_step: int, step_interval: int) -> dict:
    """!
    @brief Detect the highest contiguous fully available source step for live post-processing.
    @param[in] source_dir Argument passed to `detect_post_source_frontier()`.
    @param[in] monitor_cfg Argument passed to `detect_post_source_frontier()`.
    @param[in] post_cfg Argument passed to `detect_post_source_frontier()`.
    @param[in] start_step Argument passed to `detect_post_source_frontier()`.
    @param[in] end_step Argument passed to `detect_post_source_frontier()`.
    @param[in] step_interval Argument passed to `detect_post_source_frontier()`.
    @return Value returned by `detect_post_source_frontier()`.
    """
    diagnostic = {
        'first_requested_step': start_step,
        'first_incomplete_step': None,
        'missing_files_for_first_incomplete_step': [],
        'closest_complete_step_to_start': None,
        'closest_complete_step_to_end': None,
    }
    if step_interval <= 0 or end_step < start_step or not os.path.isdir(source_dir):
        return {
            'frontier_step': None,
            'diagnostic': diagnostic,
        }

    complete_steps, source_scan = _scan_complete_source_steps(source_dir, monitor_cfg, post_cfg)
    diagnostic['closest_complete_step_to_start'] = _nearest_step(complete_steps, start_step)
    diagnostic['closest_complete_step_to_end'] = _nearest_step(complete_steps, end_step)

    frontier = None
    for step in _iter_post_steps(start_step, end_step, step_interval):
        expected_paths = _expected_source_paths_for_step(step, source_scan, post_cfg)
        if not all(os.path.isfile(path) for path in expected_paths):
            diagnostic['first_incomplete_step'] = step
            diagnostic['missing_files_for_first_incomplete_step'] = [
                os.path.relpath(path, source_dir) for path in expected_paths if not os.path.isfile(path)
            ]
            break
        frontier = step
    return {
        'frontier_step': frontier,
        'diagnostic': diagnostic,
    }


def persist_post_resume_state(run_dir: str, plan: dict, last_successful_requested_end_step=None):
    """!
    @brief Persist post resume lineage metadata for future --continue runs.
    @param[in] run_dir Argument passed to `persist_post_resume_state()`.
    @param[in] plan Argument passed to `persist_post_resume_state()`.
    @param[in] last_successful_requested_end_step Argument passed to `persist_post_resume_state()`.
    @return Value returned by `persist_post_resume_state()`.
    """
    state_path = get_post_resume_state_path(run_dir)
    payload = {
        'schema_version': POST_RESUME_SCHEMA_VERSION,
        'run_id': plan.get('run_id'),
        'recipe_fingerprint': plan.get('recipe_fingerprint'),
        'recipe_signature': plan.get('recipe_signature'),
        'requested_start_step': plan.get('requested_start_step'),
        'requested_end_step': plan.get('requested_end_step'),
        'step_interval': plan.get('step_interval'),
        'source_directory': plan.get('source_data_directory'),
        'resume_match_source': plan.get('resume_match_source'),
        'last_successful_requested_end_step': last_successful_requested_end_step,
        'updated_at': datetime.now().isoformat(),
    }
    write_json_file(state_path, payload)
    return state_path


def _build_post_lock_wrapper_source() -> str:
    """!
    @brief Return the Python wrapper used to hold an exclusive post-stage lock.
    @return Value returned by `_build_post_lock_wrapper_source()`.
    """
    return """#!/usr/bin/env python3
import argparse
import fcntl
import json
import os
import socket
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description='PICurv post-stage lock wrapper')
    parser.add_argument('--lock-file', required=True)
    parser.add_argument('--metadata-file', required=True)
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--recipe-fingerprint', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    command = list(args.command or [])
    if not command or command[0] != '--':
        parser.error("expected '-- <command ...>' after wrapper arguments")
    command = command[1:]

    os.makedirs(os.path.dirname(args.lock_file), exist_ok=True)
    fd = os.open(args.lock_file, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        owner = None
        try:
            with open(args.metadata_file, 'r', encoding='utf-8') as handle:
                owner = json.load(handle)
        except Exception:
            owner = None
        if owner:
            print(
                f"[FATAL] Post stage already active for {args.run_dir} "
                f"(pid={owner.get('pid')}, host={owner.get('host')}, started_at={owner.get('started_at')}).",
                file=sys.stderr,
            )
        else:
            print(f"[FATAL] Post stage already active for {args.run_dir}.", file=sys.stderr)
        return 2

    metadata = {
        'pid': os.getpid(),
        'host': socket.gethostname(),
        'started_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'run_dir': args.run_dir,
        'recipe_fingerprint': args.recipe_fingerprint,
        'command': command,
    }
    with open(args.metadata_file, 'w', encoding='utf-8') as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write('\\n')

    try:
        result = subprocess.run(command)
        return int(result.returncode)
    finally:
        try:
            os.remove(args.metadata_file)
        except FileNotFoundError:
            pass
        os.close(fd)


if __name__ == '__main__':
    raise SystemExit(main())
"""


def ensure_post_lock_wrapper(run_dir: str) -> str:
    """!
    @brief Ensure the lock wrapper exists for a run directory and return its path.
    @param[in] run_dir Argument passed to `ensure_post_lock_wrapper()`.
    @return Value returned by `ensure_post_lock_wrapper()`.
    """
    paths = get_post_lock_paths(run_dir)
    wrapper_path = paths['wrapper_path']
    content = _build_post_lock_wrapper_source()
    existing_content = None
    os.makedirs(os.path.dirname(wrapper_path), exist_ok=True)
    if os.path.isfile(wrapper_path):
        with open(wrapper_path, 'r', encoding='utf-8', errors='replace') as f:
            existing_content = f.read()
    if existing_content != content:
        with open(wrapper_path, 'w', encoding='utf-8') as f:
            f.write(content)
        os.chmod(wrapper_path, 0o755)
    return wrapper_path


def build_post_locked_command(run_dir: str, recipe_fingerprint: str, wrapped_command: list, create_wrapper: bool = True) -> "tuple[list, dict]":
    """!
    @brief Wrap a postprocessor command behind the run-dir-scoped lock wrapper.
    @param[in] run_dir Argument passed to `build_post_locked_command()`.
    @param[in] recipe_fingerprint Argument passed to `build_post_locked_command()`.
    @param[in] wrapped_command Argument passed to `build_post_locked_command()`.
    @param[in] create_wrapper Argument passed to `build_post_locked_command()`.
    @return Value returned by `build_post_locked_command()`.
    """
    lock_paths = get_post_lock_paths(run_dir)
    wrapper_path = ensure_post_lock_wrapper(run_dir) if create_wrapper else lock_paths['wrapper_path']
    command = [
        wrapper_path,
        '--lock-file', lock_paths['lock_file'],
        '--metadata-file', lock_paths['metadata_file'],
        '--run-dir', run_dir,
        '--recipe-fingerprint', recipe_fingerprint,
        '--',
    ] + list(wrapped_command)
    return command, lock_paths


def build_post_execution_plan(
    run_dir: str,
    run_id: str,
    case_cfg: dict,
    monitor_cfg: dict,
    post_cfg: dict,
    continue_requested: bool = False,
    allow_source_frontier_scan: bool = True,
) -> dict:
    """!
    @brief Resolve post resume/source-availability behavior into one execution plan.
    @param[in] run_dir Argument passed to `build_post_execution_plan()`.
    @param[in] run_id Argument passed to `build_post_execution_plan()`.
    @param[in] case_cfg Argument passed to `build_post_execution_plan()`.
    @param[in] monitor_cfg Argument passed to `build_post_execution_plan()`.
    @param[in] post_cfg Argument passed to `build_post_execution_plan()`.
    @param[in] continue_requested Argument passed to `build_post_execution_plan()`.
    @param[in] allow_source_frontier_scan Argument passed to `build_post_execution_plan()`.
    @return Value returned by `build_post_execution_plan()`.
    """
    requested_start_step, requested_end_step, step_interval = resolve_post_requested_window(post_cfg, case_cfg)
    resolved_source_dir = _resolve_post_source_directory_preview(run_dir, monitor_cfg, post_cfg)
    resolved_post_cfg = prepare_effective_post_config(post_cfg, resolved_source_dir)
    recipe_cfg = build_post_recipe_config(resolved_post_cfg, monitor_cfg)
    recipe_signature, recipe_fingerprint = compute_post_recipe_fingerprint(recipe_cfg)

    state_path = get_post_resume_state_path(run_dir)
    state_payload = _read_json_if_exists(state_path)
    state_match = bool(isinstance(state_payload, dict) and state_payload.get('recipe_fingerprint') == recipe_fingerprint)

    legacy_post_run_path = os.path.join(run_dir, 'config', 'post.run')
    legacy_recipe_cfg = parse_post_recipe_file(legacy_post_run_path)
    legacy_recipe_signature = normalize_post_recipe_signature(legacy_recipe_cfg or {}) if legacy_recipe_cfg else None
    legacy_match = bool(legacy_recipe_signature and legacy_recipe_signature == recipe_signature)

    resume_recipe_match = False
    resume_match_source = None
    resume_bootstrapped = False
    if continue_requested:
        if state_match:
            resume_recipe_match = True
            resume_match_source = 'state'
        elif not state_payload and legacy_match:
            resume_recipe_match = True
            resume_match_source = 'legacy_post_run'
            resume_bootstrapped = True

    completion_info = detect_post_completed_frontier(
        run_dir,
        resolved_post_cfg,
        monitor_cfg,
        requested_start_step,
        requested_end_step,
        step_interval,
    )
    completed_frontier_step = completion_info['frontier_step']
    if completion_info['artifact_family_count'] == 0 and state_match:
        completed_frontier_step = _parse_int_loose(state_payload.get('last_successful_requested_end_step'))

    if continue_requested and resume_recipe_match and completed_frontier_step is not None:
        effective_start_step = completed_frontier_step + step_interval
    else:
        effective_start_step = requested_start_step

    source_frontier_step = None
    source_frontier_diagnostic = None
    source_frontier_deferred = not allow_source_frontier_scan
    skip_reason = None
    if effective_start_step > requested_end_step:
        skip_reason = 'already-complete-window'
        effective_end_step = requested_end_step
    elif allow_source_frontier_scan:
        source_frontier_info = detect_post_source_frontier(
            resolved_source_dir,
            monitor_cfg,
            resolved_post_cfg,
            effective_start_step,
            requested_end_step,
            step_interval,
        )
        source_frontier_step = source_frontier_info['frontier_step']
        source_frontier_diagnostic = source_frontier_info['diagnostic']
        if source_frontier_step is None or source_frontier_step < effective_start_step:
            if continue_requested and resume_recipe_match and completed_frontier_step is not None:
                skip_reason = 'already-caught-up-to-current-source-frontier'
            else:
                skip_reason = 'nothing-available-yet'
            effective_end_step = None
        else:
            effective_end_step = min(requested_end_step, source_frontier_step)
    else:
        effective_end_step = requested_end_step

    effective_post_cfg = None
    if skip_reason is None:
        effective_post_cfg = prepare_effective_post_config(
            post_cfg,
            resolved_source_dir,
            start_step=effective_start_step,
            end_step=effective_end_step,
        )

    return {
        'run_id': run_id,
        'continue_requested': bool(continue_requested),
        'requested_start_step': requested_start_step,
        'requested_end_step': requested_end_step,
        'step_interval': step_interval,
        'source_data_directory': resolved_source_dir,
        'recipe_config': recipe_cfg,
        'recipe_signature': recipe_signature,
        'recipe_fingerprint': recipe_fingerprint,
        'resume_state_path': state_path,
        'resume_state_payload': state_payload,
        'resume_recipe_match': resume_recipe_match,
        'resume_match_source': resume_match_source,
        'resume_bootstrapped': resume_bootstrapped,
        'completed_frontier_step': completed_frontier_step,
        'source_frontier_step': source_frontier_step,
        'source_frontier_diagnostic': source_frontier_diagnostic,
        'source_frontier_deferred': source_frontier_deferred,
        'effective_start_step': effective_start_step,
        'effective_end_step': effective_end_step,
        'skip_reason': skip_reason,
        'resolved_post_cfg': resolved_post_cfg,
        'effective_post_cfg': effective_post_cfg,
        'lock_paths': get_post_lock_paths(run_dir),
    }


def needs_restart_source(case_cfg: dict, solver_cfg: dict) -> bool:
    """!
    @brief Return True when the solver requires restart data from disk.
    @details Correctly identifies that analytical + init + start_step > 0 does NOT
             need a restart source (C code never reads from restart_dir in that case).
    @param[in] case_cfg Parsed case YAML dictionary.
    @param[in] solver_cfg Parsed solver YAML dictionary.
    @return True if a restart source (--restart-from or --continue) is required.
    """
    try:
        start_step = int(case_cfg.get("run_control", {}).get("start_step", 0) or 0)
    except (TypeError, ValueError):
        start_step = 0
    eulerian_source = str(
        (solver_cfg.get("operation_mode", {}) or {}).get("eulerian_field_source", "solve")
    ).strip().lower()
    particle_restart_mode = str(
        (case_cfg.get("models", {}).get("physics", {}).get("particles", {}) or {}).get("restart_mode", "init")
    ).strip().lower()
    euler_needs = (eulerian_source == "load") or (eulerian_source == "solve" and start_step > 0)
    particle_needs = (particle_restart_mode == "load")
    return euler_needs or particle_needs


def resolve_run_output_dir(run_dir: str, monitor_cfg: dict) -> str:
    """!
    @brief Resolve the output data directory within a run directory.
    @param[in] run_dir Path to the run directory.
    @param[in] monitor_cfg Parsed monitor YAML dictionary.
    @return Absolute path to the output directory.
    """
    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    output_rel = dirs.get("output", "output")
    return os.path.abspath(os.path.join(run_dir, output_rel))


def resolve_run_restart_dir(run_dir: str, monitor_cfg: dict) -> str:
    """!
    @brief Resolve the restart staging directory within a run directory.
    @param[in] run_dir Path to the run directory.
    @param[in] monitor_cfg Parsed monitor YAML dictionary.
    @return Absolute path to the restart directory.
    """
    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    restart_rel = dirs.get("restart", "restart")
    return os.path.abspath(os.path.join(run_dir, restart_rel))


def populate_restart_directory(source_output: str, target_restart: str, start_step: int, monitor_cfg: dict):
    """!
    @brief Copy checkpoint files for a specific step from source output to target restart.
    @param[in] source_output Path to the source output directory containing checkpoint data.
    @param[in] target_restart Path to the target restart directory to populate.
    @param[in] start_step The step number whose checkpoint files should be copied.
    @param[in] monitor_cfg Parsed monitor YAML dictionary (for subdirectory names).
    """
    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    euler_sub = dirs.get("eulerian_subdir", "eulerian")
    particle_sub = dirs.get("particle_subdir", "particles")

    if os.path.exists(target_restart):
        shutil.rmtree(target_restart)
    os.makedirs(os.path.join(target_restart, euler_sub), exist_ok=True)
    os.makedirs(os.path.join(target_restart, particle_sub), exist_ok=True)

    step_str = f"{start_step:05d}"
    for subdir in [euler_sub, particle_sub]:
        src = os.path.join(source_output, subdir)
        dst = os.path.join(target_restart, subdir)
        if not os.path.isdir(src):
            continue
        for f_name in glob.glob(os.path.join(src, f"*{step_str}_0.*")):
            shutil.copy2(f_name, dst)

    copied = glob.glob(os.path.join(target_restart, "**", f"*{step_str}_0.*"), recursive=True)
    if not copied:
        raise ValueError(f"No checkpoint files found for step {start_step} in {source_output}")
    print(f"[INFO] Populated restart directory with {len(copied)} file(s) for step {start_step}: {target_restart}")


def detect_last_checkpoint_step(output_dir: str, euler_subdir: str = "eulerian", particle_subdir: str = "particles"):
    """!
    @brief Scan output directory for the highest step number available.
    @details Checks eulerian files first (ufield), then falls back to particle
             files (position) for analytical-mode cases that have no eulerian output.
    @param[in] output_dir Path to the output directory.
    @param[in] euler_subdir Name of the eulerian subdirectory.
    @param[in] particle_subdir Name of the particle subdirectory.
    @return The highest step number found, or None if no checkpoints exist.
    """
    import re as _re
    euler_path = os.path.join(output_dir, euler_subdir)
    if os.path.isdir(euler_path):
        pattern = _re.compile(r"ufield(\d{5})_0\.dat")
        steps = []
        for fname in os.listdir(euler_path):
            match = pattern.match(fname)
            if match:
                steps.append(int(match.group(1)))
        if steps:
            return max(steps)
    particle_path = os.path.join(output_dir, particle_subdir)
    if os.path.isdir(particle_path):
        pattern = _re.compile(r"position(\d{5})_0\.dat")
        steps = []
        for fname in os.listdir(particle_path):
            match = pattern.match(fname)
            if match:
                steps.append(int(match.group(1)))
        if steps:
            return max(steps)
    return None


def detect_case_completion_status(run_dir: str, monitor_cfg: dict, target_final_step: int) -> dict:
    """!
    @brief Determine whether a study case is complete, partially complete, or empty.
    @param[in] run_dir Path to the case run directory.
    @param[in] monitor_cfg Parsed monitor YAML dictionary.
    @param[in] target_final_step The step number the case should reach for completion.
    @return Dictionary with keys 'last_step' (int or None), 'target_step' (int),
            and 'status' ('complete', 'partial', or 'empty').
    """
    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    euler_sub = dirs.get("eulerian_subdir", "eulerian")
    particle_sub = dirs.get("particle_subdir", "particles")
    output_dir = resolve_run_output_dir(run_dir, monitor_cfg)
    last_step = detect_last_checkpoint_step(output_dir, euler_sub, particle_sub)
    if last_step is not None and last_step >= target_final_step:
        status = "complete"
    elif last_step is not None:
        status = "partial"
    else:
        status = "empty"
    return {"last_step": last_step, "target_step": target_final_step, "status": status}


def validate_load_mode_step_range(source_output: str, start_step: int, total_steps: int, monitor_cfg: dict):
    """!
    @brief Validate that all required eulerian step files exist for "load" mode.
    @details Checks that ufield files exist for every step from start_step through
             start_step + total_steps (inclusive). Reports missing steps clearly.
    @param[in] source_output Path to the output directory containing eulerian data.
    @param[in] start_step First step that will be loaded.
    @param[in] total_steps Number of steps to run.
    @param[in] monitor_cfg Parsed monitor YAML dictionary (for subdirectory names).
    """
    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    euler_sub = dirs.get("eulerian_subdir", "eulerian")
    euler_path = os.path.join(source_output, euler_sub)

    missing = []
    for step in range(start_step, start_step + total_steps + 1):
        expected = os.path.join(euler_path, f"ufield{step:05d}_0.dat")
        if not os.path.isfile(expected):
            missing.append(step)

    if missing:
        sample = missing[:3] + (["..."] if len(missing) > 6 else []) + missing[-3:]
        raise ValueError(
            f"Eulerian 'load' mode: {len(missing)} step file(s) missing in {euler_path}. "
            f"Missing steps include: {sample}"
        )


def validate_particle_checkpoint(source_dir: str, start_step: int, monitor_cfg: dict):
    """!
    @brief Validate that particle checkpoint files exist for the given step.
    @details Checks that at least a position file exists at the expected step in
             the particle subdirectory.
    @param[in] source_dir Path to the directory containing the particle subdirectory.
    @param[in] start_step The step number whose particle checkpoint is expected.
    @param[in] monitor_cfg Parsed monitor YAML dictionary (for subdirectory names).
    """
    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    particle_sub = dirs.get("particle_subdir", "particles")
    particle_path = os.path.join(source_dir, particle_sub)
    expected = os.path.join(particle_path, f"position{start_step:05d}_0.dat")
    if not os.path.isfile(expected):
        raise ValueError(
            f"Particle restart_mode='load' but checkpoint not found: {expected}"
        )


def read_monitor_from_run(run_dir: str) -> dict:
    """!
    @brief Read the monitor.yml from a run directory's config/ subdirectory.
    @param[in] run_dir Path to the run directory.
    @return Parsed monitor YAML dictionary.
    """
    monitor_path = os.path.join(run_dir, "config", "monitor.yml")
    if not os.path.isfile(monitor_path):
        raise ValueError(f"Run directory is missing config/monitor.yml: {monitor_path}")
    return read_yaml_file(monitor_path)


def resolve_restart_source(args, case_cfg: dict, solver_cfg: dict, monitor_cfg: dict, run_dir: str):
    """!
    @brief Resolve the restart source directory based on --restart-from or --continue CLI flags.
    @details Implements the full restart resolution logic including smart resolution for
             --continue (checks restart/ first for user-curated data, falls back to output/)
             and direct reference for eulerian "load" mode.
    @param[in] args Parsed CLI arguments (must have restart_from, continue_run, run_dir attrs).
    @param[in] case_cfg Parsed case YAML dictionary.
    @param[in] solver_cfg Parsed solver YAML dictionary.
    @param[in] monitor_cfg Parsed monitor YAML dictionary.
    @param[in] run_dir Path to the current run directory.
    @return Tuple of (restart_source_dir, continue_mode) where restart_source_dir is the
            resolved path (or None) and continue_mode is a boolean.
    """
    try:
        start_step = int(case_cfg.get("run_control", {}).get("start_step", 0) or 0)
    except (TypeError, ValueError):
        start_step = 0
    try:
        total_steps = int(case_cfg.get("run_control", {}).get("total_steps", 0) or 0)
    except (TypeError, ValueError):
        total_steps = 0

    eulerian_source = str(
        (solver_cfg.get("operation_mode", {}) or {}).get("eulerian_field_source", "solve")
    ).strip().lower()

    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    euler_sub = dirs.get("eulerian_subdir", "eulerian")
    particle_sub = dirs.get("particle_subdir", "particles")

    particle_restart_mode = str(
        (case_cfg.get("models", {}).get("physics", {}).get("particles", {}) or {}).get("restart_mode", "init")
    ).strip().lower()
    particle_needs = (particle_restart_mode == "load")

    requires_source = needs_restart_source(case_cfg, solver_cfg)
    restart_from = getattr(args, 'restart_from', None)
    continue_run = getattr(args, 'continue_run', False)

    if restart_from and continue_run:
        raise ValueError("--restart-from and --continue are mutually exclusive.")

    if restart_from:
        # === MODE 1: New run, restart from another run ===
        source_run = os.path.abspath(restart_from)
        if not os.path.isdir(source_run):
            raise ValueError(f"--restart-from run directory does not exist: {source_run}")
        source_monitor = read_monitor_from_run(source_run)
        source_output = resolve_run_output_dir(source_run, source_monitor)
        if not os.path.isdir(source_output):
            raise ValueError(f"Source output directory does not exist: {source_output}")

        if not requires_source:
            # R7: analytical + init — warn that --restart-from is unused
            print(
                "[WARN] --restart-from specified but no data will be read "
                "(analytical + init does not need restart data).",
                file=sys.stderr,
            )
            return None, False

        if eulerian_source == "load":
            # R3/R4/R5: Direct reference to source output
            validate_load_mode_step_range(source_output, start_step, total_steps, source_monitor)
            if particle_needs:
                validate_particle_checkpoint(source_output, start_step, source_monitor)
            return source_output, False
        else:
            # R1/R2/R6: Copy checkpoint to new run's restart/
            target_restart = resolve_run_restart_dir(run_dir, monitor_cfg)
            populate_restart_directory(source_output, target_restart, start_step, monitor_cfg)
            if particle_needs:
                validate_particle_checkpoint(target_restart, start_step, monitor_cfg)
            return target_restart, False

    elif continue_run:
        # === MODE 2: Continue in-place ===
        continue_run_dir = getattr(args, 'run_dir', None)
        if not continue_run_dir:
            raise ValueError("--continue requires --run-dir.")
        continue_run_dir = os.path.abspath(continue_run_dir)
        if not os.path.isdir(continue_run_dir):
            raise ValueError(f"--run-dir does not exist: {continue_run_dir}")

        source_output = resolve_run_output_dir(continue_run_dir, monitor_cfg)
        target_restart = resolve_run_restart_dir(continue_run_dir, monitor_cfg)

        # Warn if start_step != last checkpoint
        last_step = detect_last_checkpoint_step(source_output, euler_sub)
        if last_step is not None and last_step != start_step:
            print(
                f"[WARN] start_step={start_step} but last checkpoint in output is step {last_step}.",
                file=sys.stderr,
            )

        if eulerian_source == "load":
            # C3/C4: Direct reference to output/ (all euler steps there, can't copy)
            validate_load_mode_step_range(source_output, start_step, total_steps, monitor_cfg)
            if particle_needs:
                validate_particle_checkpoint(source_output, start_step, monitor_cfg)
            return source_output, True
        elif not requires_source:
            # C6: analytical + init — only log-append behavior, no data needed
            return None, True
        else:
            # C1/C2/C5: Smart resolution for solve/analytical
            euler_needs = (eulerian_source == "solve" and start_step > 0)

            step_str = f"{start_step:05d}"
            needed_subs = []
            if euler_needs:
                needed_subs.append(euler_sub)
            if particle_needs:
                needed_subs.append(particle_sub)

            # Check restart/ for EACH needed component
            restart_has = {}
            if os.path.isdir(target_restart):
                for sub in needed_subs:
                    restart_has[sub] = bool(
                        glob.glob(os.path.join(target_restart, sub, f"*{step_str}_0.*"))
                    )

            all_in_restart = needed_subs and all(restart_has.get(s, False) for s in needed_subs)
            some_in_restart = any(restart_has.get(s, False) for s in needed_subs)

            if all_in_restart:
                # Fully curated restart/ (warm-up-and-discard workflow)
                print(f"[INFO] Using curated restart directory: {target_restart}")
                return target_restart, True
            elif os.path.isdir(source_output):
                # Auto-populate missing components from output/ into restart/
                if some_in_restart:
                    # Partial curate: only copy components NOT already in restart/
                    missing = [s for s in needed_subs if not restart_has.get(s, False)]
                    os.makedirs(target_restart, exist_ok=True)
                    for sub in missing:
                        src = os.path.join(source_output, sub)
                        dst = os.path.join(target_restart, sub)
                        os.makedirs(dst, exist_ok=True)
                        if os.path.isdir(src):
                            for f_name in glob.glob(os.path.join(src, f"*{step_str}_0.*")):
                                shutil.copy2(f_name, dst)
                    filled = glob.glob(os.path.join(target_restart, "**", f"*{step_str}_0.*"), recursive=True)
                    print(f"[INFO] Merged curated restart/ with {len(missing)} component(s) from output/: {target_restart}")
                    if not filled:
                        raise ValueError(
                            f"After merge, no checkpoint files found for step {start_step} in {target_restart}"
                        )
                    return target_restart, True
                else:
                    # Nothing curated — full populate from output/
                    populate_restart_directory(source_output, target_restart, start_step, monitor_cfg)
                    return target_restart, True
            else:
                raise ValueError(
                    f"Neither restart/ nor output/ contain data for step {start_step}. "
                    f"Checked: {target_restart}, {source_output}"
                )

    elif requires_source:
        raise ValueError(
            "Restart data required but no source specified. Use:\n"
            "  --restart-from <run_dir>  (new run from another run's data)\n"
            "  --continue --run-dir <run_dir>  (resume in same directory)"
        )

    return None, False

def absolutize_case_external_paths(case_cfg: dict, case_anchor_path: str):
    """!
    @brief Convert external grid/generator paths in case config to absolute paths.
    @param[in] case_cfg Argument passed to `absolutize_case_external_paths()`.
    @param[in] case_anchor_path Argument passed to `absolutize_case_external_paths()`.
    """
    grid_cfg = case_cfg.get("grid", {})
    if not isinstance(grid_cfg, dict):
        return
    mode = grid_cfg.get("mode")
    if mode == "file":
        source_file = grid_cfg.get("source_file")
        if isinstance(source_file, str):
            grid_cfg["source_file"] = resolve_path(case_anchor_path, source_file)
        legacy_cfg = grid_cfg.get("legacy_conversion", {})
        if isinstance(legacy_cfg, dict):
            script_path = legacy_cfg.get("script")
            if isinstance(script_path, str):
                legacy_cfg["script"] = resolve_path(case_anchor_path, script_path)
    elif mode == "grid_gen":
        gen = grid_cfg.get("generator", {})
        if isinstance(gen, dict):
            for key in ("script", "config_file"):
                val = gen.get(key)
                if isinstance(val, str):
                    gen[key] = resolve_path(case_anchor_path, val)
    ic = (case_cfg.get("properties", {}) or {}).get("initial_conditions", {})
    if isinstance(ic, dict):
        if str(ic.get("mode", "")).strip().lower() == "file":
            source_file = ic.get("source_file")
            if isinstance(source_file, str):
                ic["source_file"] = resolve_path(case_anchor_path, source_file)
        elif str(ic.get("generator", "")).strip().lower() == "ic_gen":
            params = ic.get("params", {})
            if isinstance(params, dict):
                for key in ("script", "config_file"):
                    value = params.get(key)
                    if isinstance(value, str):
                        params[key] = resolve_path(case_anchor_path, value)
    boundary_conditions = case_cfg.get("boundary_conditions", [])
    blocks = boundary_conditions if boundary_conditions and isinstance(boundary_conditions[0], list) else [boundary_conditions]
    for block in blocks:
        if not isinstance(block, list):
            continue
        for bc in block:
            if not isinstance(bc, dict) or str(bc.get("handler", "")).strip().lower() != "prescribed_flow":
                continue
            source = ((bc.get("params") or {}).get("source") or {})
            if not isinstance(source, dict):
                continue
            source_type = str(source.get("type", "")).strip().lower()
            if source_type == "file":
                keys = ("path",)
            elif source_type == "generated":
                keys = ("script",)
            elif source_type == "field_slice":
                keys = ("script", "field_file", "grid_file", "source_case")
            else:
                keys = ()
            for key in keys:
                value = source.get(key)
                if isinstance(value, str):
                    source[key] = resolve_path(case_anchor_path, value)


def prepare_case_for_continuation(run_dir: str, case_id: str, last_step: int,
                                   target_final_step: int, cluster_cfg: dict):
    """!
    @brief Set up a partially-completed study case for continuation in-place.
    @details Updates the case config with new start_step/total_steps, sets particle
             restart_mode to 'load' if checkpoint exists, populates the restart
             directory, and regenerates the solver control file with continue_mode.
             Delegates all restart resolution to resolve_restart_source().
    @param[in] run_dir Path to the case run directory.
    @param[in] case_id The case identifier (e.g. 'case_0002').
    @param[in] last_step The last checkpoint step found in the output directory.
    @param[in] target_final_step The step number the case should reach for completion.
    @param[in] cluster_cfg Parsed cluster YAML dictionary (for num_procs, walltime guard).
    @return The absolute path to the regenerated control file.
    """
    config_dir = os.path.join(run_dir, "config")
    case_cfg = read_yaml_file(os.path.join(config_dir, "case.yml"))
    solver_cfg = read_yaml_file(os.path.join(config_dir, "solver.yml"))
    monitor_cfg = read_yaml_file(os.path.join(config_dir, "monitor.yml"))

    remaining = target_final_step - last_step
    case_cfg["run_control"]["start_step"] = last_step
    case_cfg["run_control"]["total_steps"] = remaining
    print(f"[INFO] {case_id}: updating start_step={last_step}, total_steps={remaining}")

    dirs = (monitor_cfg.get("io", {}) or {}).get("directories", {}) or {}
    particle_sub = dirs.get("particle_subdir", "particles")
    output_dir = resolve_run_output_dir(run_dir, monitor_cfg)
    particle_ckpt = os.path.join(output_dir, particle_sub, f"position{last_step:05d}_0.dat")
    particles_cfg = (case_cfg.get("models", {}).get("physics", {}) or {}).get("particles")
    if particles_cfg is not None and os.path.isfile(particle_ckpt):
        current_mode = str(particles_cfg.get("restart_mode", "init")).strip().lower()
        if current_mode != "load":
            particles_cfg["restart_mode"] = "load"
            print(f"[INFO] {case_id}: setting particle restart_mode='load' (checkpoint found)")

    write_yaml_file(os.path.join(config_dir, "case.yml"), case_cfg)

    mock_args = argparse.Namespace(restart_from=None, continue_run=True, run_dir=run_dir)
    restart_source_dir, continue_mode = resolve_restart_source(
        mock_args, case_cfg, solver_cfg, monitor_cfg, run_dir
    )

    source_files = {
        'Case': os.path.join(config_dir, "case.yml"),
        'Solver': os.path.join(config_dir, "solver.yml"),
        'Monitor': os.path.join(config_dir, "monitor.yml"),
    }
    monitor_files = prepare_monitor_files(run_dir, case_id, monitor_cfg, source_files)
    cluster_tasks = get_cluster_total_tasks(cluster_cfg)
    configs = {
        "case": case_cfg, "case_path": source_files['Case'],
        "solver": solver_cfg, "solver_path": source_files['Solver'],
        "monitor": monitor_cfg, "monitor_path": source_files['Monitor'],
        "walltime_guard_policy": resolve_walltime_guard_policy(cluster_cfg),
    }
    control_file = generate_solver_control_file(
        run_dir, case_id, configs, cluster_tasks, monitor_files,
        restart_source_dir=restart_source_dir, continue_mode=continue_mode,
    )
    print(f"[SUCCESS] {case_id}: regenerated control file for continuation")
    return control_file


def is_valid_email(email: str) -> bool:
    """!
    @brief Lightweight email validation for scheduler notifications.
    @param[in] email Argument passed to `is_valid_email()`.
    @return Value returned by `is_valid_email()`.
    """
    if not isinstance(email, str):
        return False
    pattern = r"^[^@\s]+@[^@\s]+\.[^@\s]+$"
    return re.match(pattern, email.strip()) is not None

def normalize_statistics_task(task_name: str) -> str:
    """!
    @brief Normalizes user-facing statistics task names to C pipeline keywords.
    @param[in] task_name Task name from YAML.
    @return Canonical keyword accepted by C statistics pipeline.
    @throws ValueError if task is unsupported.
    """
    # Only implemented statistics kernels belong here.
    if task_name is None:
        raise ValueError("statistics task cannot be None")
    normalized = str(task_name).strip().lower().replace("-", "_").replace(" ", "_")
    if normalized != "msd":
        raise ValueError(f"Unsupported statistics task '{task_name}'. Currently supported: 'msd'.")
    return "ComputeMSD"

def _iter_nonempty_noncomment_lines(file_obj):
    """!
    @brief Yield (lineno, stripped_line) for non-empty, non-comment lines.
    @param[in] file_obj Argument passed to `_iter_nonempty_noncomment_lines()`.
    """
    for lineno, raw in enumerate(file_obj, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield lineno, line

def validate_and_nondimensionalize_picgrid(source_grid: str, dest_grid: str, L_ref: float, expected_nblk: int = None) -> dict:
    """!
    @brief Validates PICGRID payload and writes a non-dimensionalized copy.
    @details Requires canonical PICGRID input with leading "PICGRID" token.
             Output is always written in canonical PICGRID format with header and per-block dims.
    @param[in] source_grid Input grid file path.
    @param[in] dest_grid Output grid file path.
    @param[in] L_ref Reference length for non-dimensionalization.
    @param[in] expected_nblk Optional expected block count.
    @return Summary dictionary with nblk, dims, and total_nodes.
    @throws ValueError on malformed grid.
    """
    if L_ref == 0.0:
        raise ValueError("length_ref must be non-zero when processing grid coordinates.")
    if not os.path.isfile(source_grid):
        raise ValueError(f"Grid file not found: {source_grid}")

    with open(source_grid, "r") as fin:
        line_iter = _iter_nonempty_noncomment_lines(fin)
        try:
            _, first_token = next(line_iter)
        except StopIteration:
            raise ValueError(f"Grid file '{source_grid}' is empty.")

        if first_token != "PICGRID":
            raise ValueError(
                f"Grid file '{source_grid}' must begin with the canonical PICGRID header token."
            )
        try:
            _, nblk_line = next(line_iter)
        except StopIteration:
            raise ValueError(f"Grid file '{source_grid}' missing block count after PICGRID header.")

        try:
            nblk = int(nblk_line)
        except ValueError:
            raise ValueError(f"Invalid block count '{nblk_line}' in grid file '{source_grid}'.")
        if nblk <= 0:
            raise ValueError(f"Grid file '{source_grid}' has non-positive block count: {nblk}.")
        if expected_nblk is not None and nblk != expected_nblk:
            raise ValueError(
                f"Grid file block count mismatch: case expects {expected_nblk}, grid contains {nblk}."
            )

        dims = []
        for bi in range(nblk):
            try:
                lineno, dim_line = next(line_iter)
            except StopIteration:
                raise ValueError(f"Grid file '{source_grid}' missing dimensions for block {bi}.")
            parts = dim_line.split()
            if len(parts) != 3:
                raise ValueError(
                    f"Invalid dimensions line at {source_grid}:{lineno}. Expected 3 integers, got: '{dim_line}'."
                )
            try:
                im, jm, km = (int(parts[0]), int(parts[1]), int(parts[2]))
            except ValueError:
                raise ValueError(
                    f"Invalid dimensions line at {source_grid}:{lineno}. Non-integer values: '{dim_line}'."
                )
            if im <= 0 or jm <= 0 or km <= 0:
                raise ValueError(
                    f"Invalid block dimensions at {source_grid}:{lineno}: ({im}, {jm}, {km}). Must be > 0."
                )
            dims.append((im, jm, km))

        total_nodes_expected = sum(im * jm * km for (im, jm, km) in dims)
        os.makedirs(os.path.dirname(dest_grid), exist_ok=True)
        with open(dest_grid, "w") as fout:
            fout.write("PICGRID\n")
            fout.write(f"{nblk}\n")
            for (im, jm, km) in dims:
                fout.write(f"{im} {jm} {km}\n")

            total_nodes_seen = 0
            for lineno, coord_line in line_iter:
                parts = coord_line.split()
                if len(parts) != 3:
                    raise ValueError(
                        f"Invalid coordinate row at {source_grid}:{lineno}. Expected 3 floats, got: '{coord_line}'."
                    )
                try:
                    x = float(parts[0]) / L_ref
                    y = float(parts[1]) / L_ref
                    z = float(parts[2]) / L_ref
                except ValueError:
                    raise ValueError(
                        f"Invalid coordinate row at {source_grid}:{lineno}. Non-numeric values: '{coord_line}'."
                    )
                total_nodes_seen += 1
                if total_nodes_seen > total_nodes_expected:
                    raise ValueError(
                        f"Grid file '{source_grid}' has more coordinates ({total_nodes_seen}) than expected ({total_nodes_expected})."
                    )
                fout.write(f"{x:.8e} {y:.8e} {z:.8e}\n")

        if total_nodes_seen != total_nodes_expected:
            raise ValueError(
                f"Grid file '{source_grid}' has {total_nodes_seen} coordinates, expected {total_nodes_expected} from header."
            )

    return {"nblk": nblk, "dims": dims, "total_nodes": total_nodes_expected}

def read_picgrid_header_dimensions(source_grid: str, expected_nblk: int = None) -> list:
    """!
    @brief Read only the canonical PICGRID header dimensions.
    @param[in] source_grid Input grid file path.
    @param[in] expected_nblk Optional expected block count.
    @return List of (IM, JM, KM) node-count tuples.
    @throws ValueError on malformed header.
    """
    if not os.path.isfile(source_grid):
        raise ValueError(f"Grid file not found: {source_grid}")

    with open(source_grid, "r") as fin:
        line_iter = _iter_nonempty_noncomment_lines(fin)
        try:
            _, first_token = next(line_iter)
        except StopIteration:
            raise ValueError(f"Grid file '{source_grid}' is empty.")
        if first_token != "PICGRID":
            raise ValueError(f"Grid file '{source_grid}' must begin with the canonical PICGRID header token.")

        try:
            _, nblk_line = next(line_iter)
            nblk = int(nblk_line)
        except StopIteration:
            raise ValueError(f"Grid file '{source_grid}' missing block count after PICGRID header.")
        except ValueError:
            raise ValueError(f"Invalid block count '{nblk_line}' in grid file '{source_grid}'.")
        if nblk <= 0:
            raise ValueError(f"Grid file '{source_grid}' has non-positive block count: {nblk}.")
        if expected_nblk is not None and nblk != expected_nblk:
            raise ValueError(f"Grid file block count mismatch: case expects {expected_nblk}, grid contains {nblk}.")

        dims = []
        for bi in range(nblk):
            try:
                lineno, dim_line = next(line_iter)
            except StopIteration:
                raise ValueError(f"Grid file '{source_grid}' missing dimensions for block {bi}.")
            parts = dim_line.split()
            if len(parts) != 3:
                raise ValueError(
                    f"Invalid dimensions line at {source_grid}:{lineno}. Expected 3 integers, got: '{dim_line}'."
                )
            try:
                im, jm, km = (int(parts[0]), int(parts[1]), int(parts[2]))
            except ValueError:
                raise ValueError(
                    f"Invalid dimensions line at {source_grid}:{lineno}. Non-integer values: '{dim_line}'."
                )
            if im <= 0 or jm <= 0 or km <= 0:
                raise ValueError(
                    f"Invalid block dimensions at {source_grid}:{lineno}: ({im}, {jm}, {km}). Must be > 0."
                )
            dims.append((im, jm, km))

    return dims

def validate_and_nondimensionalize_picslice(source_slice: str, dest_slice: str, U_ref: float,
                                            expected_dims: tuple = None) -> dict:
    """!
    @brief Validate a canonical PICSLICE payload and write a solver-scale copy.
    @param[in] source_slice Input PICSLICE path.
    @param[in] dest_slice Output staged PICSLICE path.
    @param[in] U_ref Reference velocity for non-dimensionalization.
    @param[in] expected_dims Optional expected (n1, n2) slice dimensions.
    @return Summary dictionary with frame_count, dims, value_count, min_speed, max_speed.
    @throws ValueError on malformed slice.
    """
    if U_ref == 0.0:
        raise ValueError("velocity_ref must be non-zero when processing PICSLICE speeds.")
    if not os.path.isfile(source_slice):
        raise ValueError(f"PICSLICE file not found: {source_slice}")

    with open(source_slice, "r") as fin:
        line_iter = _iter_nonempty_noncomment_lines(fin)
        try:
            _, first_token = next(line_iter)
        except StopIteration:
            raise ValueError(f"PICSLICE file '{source_slice}' is empty.")
        if first_token != "PICSLICE":
            raise ValueError(f"PICSLICE file '{source_slice}' must begin with the canonical PICSLICE header token.")

        try:
            _, frame_line = next(line_iter)
            frame_count = int(frame_line)
        except StopIteration:
            raise ValueError(f"PICSLICE file '{source_slice}' missing frame count after PICSLICE header.")
        except ValueError:
            raise ValueError(f"Invalid frame count '{frame_line}' in PICSLICE file '{source_slice}'.")
        if frame_count != 1:
            raise ValueError(
                f"PICSLICE file '{source_slice}' has frame count {frame_count}; Phase 1 supports exactly 1."
            )

        try:
            lineno, dim_line = next(line_iter)
        except StopIteration:
            raise ValueError(f"PICSLICE file '{source_slice}' missing slice dimensions.")
        parts = dim_line.split()
        if len(parts) != 2:
            raise ValueError(
                f"Invalid PICSLICE dimensions at {source_slice}:{lineno}. Expected 2 integers, got: '{dim_line}'."
            )
        try:
            n1, n2 = (int(parts[0]), int(parts[1]))
        except ValueError:
            raise ValueError(
                f"Invalid PICSLICE dimensions at {source_slice}:{lineno}. Non-integer values: '{dim_line}'."
            )
        if n1 <= 0 or n2 <= 0:
            raise ValueError(f"Invalid PICSLICE dimensions at {source_slice}:{lineno}: ({n1}, {n2}). Must be > 0.")
        if expected_dims is not None and (n1, n2) != tuple(expected_dims):
            raise ValueError(
                f"PICSLICE dimension mismatch for '{source_slice}': expected {tuple(expected_dims)}, found {(n1, n2)}."
            )

        values = []
        for lineno, value_line in line_iter:
            parts = value_line.split()
            if len(parts) != 1:
                raise ValueError(
                    f"Invalid PICSLICE value row at {source_slice}:{lineno}. Expected 1 float, got: '{value_line}'."
                )
            try:
                value = float(parts[0])
            except ValueError:
                raise ValueError(f"Invalid PICSLICE value at {source_slice}:{lineno}: '{value_line}'.")
            if not math.isfinite(value):
                raise ValueError(f"PICSLICE value at {source_slice}:{lineno} must be finite.")
            if value < 0.0:
                raise ValueError(f"PICSLICE value at {source_slice}:{lineno} must be nonnegative.")
            values.append(value)

    expected_count = n1 * n2
    if len(values) != expected_count:
        raise ValueError(
            f"PICSLICE file '{source_slice}' has {len(values)} values, expected {expected_count} from dimensions {(n1, n2)}."
        )

    os.makedirs(os.path.dirname(dest_slice), exist_ok=True)
    with open(dest_slice, "w") as fout:
        fout.write("PICSLICE\n")
        fout.write("1\n")
        fout.write(f"{n1} {n2}\n")
        for value in values:
            fout.write(f"{value / U_ref:.8e}\n")

    return {
        "frame_count": frame_count,
        "dims": (n1, n2),
        "value_count": len(values),
        "min_speed": min(values) if values else 0.0,
        "max_speed": max(values) if values else 0.0,
    }

def _face_artifact_token(face: str) -> str:
    """!
    @brief Convert a BC face token into a filesystem-friendly artifact token.
    @param[in] face Canonical face token such as -Zeta.
    @return Filesystem-friendly face token.
    """
    return face.replace("+", "pos").replace("-", "neg")

def _resolve_run_artifact_path(run_dir: str, configured_path: str, default_path: str,
                               default_to_config_dir: bool = False) -> str:
    """!
    @brief Resolve a run artifact path with run-dir-relative defaults.
    @param[in] run_dir Run/precompute directory root.
    @param[in] configured_path Optional user-provided artifact path.
    @param[in] default_path Default path relative to run_dir.
    @param[in] default_to_config_dir If true, bare relative names are placed under config/.
    @return Absolute artifact path.
    """
    path = configured_path if configured_path else default_path
    if not isinstance(path, str) or not path.strip():
        raise ValueError("generated profile output_file must be a non-empty path when provided.")
    path = path.strip()
    if os.path.isabs(path):
        return os.path.abspath(path)
    if default_to_config_dir and os.path.dirname(path) == "":
        path = os.path.join("config", path)
    return os.path.abspath(os.path.join(run_dir, path))

def _resolve_generator_script(configured_script: str, case_path: str, default_name: str) -> str:
    """!
    @brief Resolve an optional generator script override or repository default.
    @param[in] configured_script Optional absolute or case-relative script path.
    @param[in] case_path Current case.yml path used to anchor relative overrides.
    @param[in] default_name Repository generator filename under GENERATORS_PATH.
    @return Absolute generator script path.
    """
    if configured_script is None:
        return os.path.join(GENERATORS_PATH, default_name)
    if not isinstance(configured_script, str) or not configured_script.strip():
        raise ValueError(f"Generator script override for {default_name} must be a non-empty path.")
    script = configured_script.strip()
    if os.path.isabs(script):
        return os.path.abspath(script)
    case_dir = os.path.dirname(os.path.abspath(case_path)) if case_path else os.getcwd()
    return os.path.abspath(os.path.join(case_dir, script))

def _normalize_square_duct_poiseuille_params(params, field_name: str) -> dict:
    """!
    @brief Validate square-duct Poiseuille generator parameters.
    @param[in] params Generator params mapping.
    @param[in] field_name Human-readable YAML field name for diagnostics.
    @return Normalized params.
    """
    if params is None:
        params = {}
    if not isinstance(params, dict):
        raise ValueError(f"{field_name}.params must be a mapping when provided.")
    unknown = sorted(set(params.keys()) - {"bulk_velocity", "n_terms"})
    if unknown:
        raise ValueError(f"Unknown keys in {field_name}.params: {unknown}. Allowed: ['bulk_velocity', 'n_terms'].")
    bulk_velocity = _to_float(params.get("bulk_velocity", 1.0), f"{field_name}.params.bulk_velocity")
    if bulk_velocity <= 0.0:
        raise ValueError(f"{field_name}.params.bulk_velocity must be positive.")
    try:
        n_terms = int(params.get("n_terms", 101))
    except (TypeError, ValueError):
        raise ValueError(f"{field_name}.params.n_terms must be a positive odd integer.")
    if n_terms <= 0 or n_terms % 2 == 0:
        raise ValueError(f"{field_name}.params.n_terms must be a positive odd integer.")
    return {"bulk_velocity": bulk_velocity, "n_terms": n_terms}

GENERATED_PROFILE_GENERATORS = {"square_duct_poiseuille"}

def _normalize_field_slice_source(source, field_name: str) -> dict:
    """!
    @brief Validate a prescribed_flow field_slice source block.
    @param[in] source Source mapping from case.yml.
    @param[in] field_name Human-readable YAML path for diagnostics.
    @return Normalized source mapping.
    """
    allowed = {
        "type",
        "script",
        "field_file",
        "grid_file",
        "source_case",
        "velocity_scale",
        "source_block",
        "output_file",
        "slice",
    }
    unknown = sorted(set(source.keys()) - allowed)
    if unknown:
        raise ValueError(f"Unknown keys in {field_name}: {unknown}. Allowed: {sorted(allowed)}.")
    field_file = source.get("field_file")
    grid_file = source.get("grid_file")
    if not isinstance(field_file, str) or not field_file.strip():
        raise ValueError(f"{field_name}.field_file must be a non-empty path.")
    if not isinstance(grid_file, str) or not grid_file.strip():
        raise ValueError(f"{field_name}.grid_file must be a non-empty path.")
    if source.get("source_case") is None and source.get("velocity_scale") is None:
        raise ValueError(f"{field_name} requires source_case or velocity_scale.")

    normalized = {
        "type": "field_slice",
        "field_file": field_file.strip(),
        "grid_file": grid_file.strip(),
        "slice": _normalize_field_slice_selector(source.get("slice"), f"{field_name}.slice"),
    }
    if source.get("script") is not None:
        script = source.get("script")
        if not isinstance(script, str) or not script.strip():
            raise ValueError(f"{field_name}.script must be a non-empty path when provided.")
        normalized["script"] = script.strip()
    if source.get("source_case") is not None:
        source_case = source.get("source_case")
        if not isinstance(source_case, str) or not source_case.strip():
            raise ValueError(f"{field_name}.source_case must be a non-empty path when provided.")
        normalized["source_case"] = source_case.strip()
    if source.get("velocity_scale") is not None:
        velocity_scale = _to_float(source.get("velocity_scale"), f"{field_name}.velocity_scale")
        if velocity_scale <= 0.0:
            raise ValueError(f"{field_name}.velocity_scale must be positive.")
        normalized["velocity_scale"] = velocity_scale
    if source.get("source_block") is not None:
        try:
            source_block = int(source.get("source_block"))
        except (TypeError, ValueError):
            raise ValueError(f"{field_name}.source_block must be a non-negative integer.")
        if source_block < 0:
            raise ValueError(f"{field_name}.source_block must be a non-negative integer.")
        normalized["source_block"] = source_block
    if source.get("output_file") is not None:
        output_file = source.get("output_file")
        if not isinstance(output_file, str) or not output_file.strip():
            raise ValueError(f"{field_name}.output_file must be a non-empty path when provided.")
        normalized["output_file"] = output_file.strip()
    return normalized

def _normalize_field_slice_selector(slice_cfg, field_name: str) -> dict:
    """!
    @brief Validate the field_slice slice selector.
    @param[in] slice_cfg Slice selector mapping.
    @param[in] field_name Human-readable YAML path for diagnostics.
    @return Normalized selector mapping.
    """
    if not isinstance(slice_cfg, dict):
        raise ValueError(f"{field_name} must be a mapping.")
    orientation = str(slice_cfg.get("orientation", "opposite")).strip().lower()
    if orientation not in {"opposite", "same"}:
        raise ValueError(f"{field_name}.orientation must be 'opposite' or 'same'.")
    normal_tolerance = _to_float(slice_cfg.get("normal_tolerance", 0.99), f"{field_name}.normal_tolerance")
    if normal_tolerance <= 0.0 or normal_tolerance > 1.0:
        raise ValueError(f"{field_name}.normal_tolerance must be in the range (0, 1].")

    if slice_cfg.get("face") is not None:
        unknown = sorted(set(slice_cfg.keys()) - {"face", "orientation", "normal_tolerance"})
        if unknown:
            raise ValueError(
                f"Unknown keys in {field_name}: {unknown}. "
                "Use either face or axis/index/normal, plus orientation/normal_tolerance."
            )
        face = str(slice_cfg.get("face", "")).strip()
        if face.lower() not in BC_FACE_MAP:
            raise ValueError(f"{field_name}.face must be one of {sorted(BC_FACE_MAP.values())}.")
        return {
            "face": BC_FACE_MAP[face.lower()],
            "orientation": orientation,
            "normal_tolerance": normal_tolerance,
        }

    required = {"axis", "index", "normal"}
    missing = sorted(key for key in required if slice_cfg.get(key) is None)
    if missing:
        raise ValueError(f"{field_name} requires either face or axis/index/normal; missing {missing}.")
    unknown = sorted(set(slice_cfg.keys()) - {"axis", "index", "normal", "orientation", "normal_tolerance"})
    if unknown:
        raise ValueError(
            f"Unknown keys in {field_name}: {unknown}. "
            "Use either face or axis/index/normal, plus orientation/normal_tolerance."
        )
    axis = str(slice_cfg.get("axis", "")).strip()
    axis_map = {"xi": "Xi", "eta": "Eta", "zeta": "Zeta"}
    if axis.lower() not in axis_map:
        raise ValueError(f"{field_name}.axis must be one of Xi, Eta, Zeta.")
    normal = str(slice_cfg.get("normal", "")).strip()
    if normal.lower() not in BC_FACE_MAP:
        raise ValueError(f"{field_name}.normal must be one of {sorted(BC_FACE_MAP.values())}.")
    normal = BC_FACE_MAP[normal.lower()]
    if normal[1:].lower() != axis.lower():
        raise ValueError(f"{field_name}.normal must use the same axis as {field_name}.axis.")
    try:
        index = int(slice_cfg.get("index"))
    except (TypeError, ValueError):
        raise ValueError(f"{field_name}.index must be an integer.")
    if index < 0:
        raise ValueError(f"{field_name}.index must be non-negative.")
    return {
        "axis": axis_map[axis.lower()],
        "index": index,
        "normal": normal,
        "orientation": orientation,
        "normal_tolerance": normal_tolerance,
    }

def generate_square_duct_poiseuille_picslice(output_path: str, dims: tuple, params: dict,
                                             target_grid: str = None, target_block: int = 0,
                                             target_face: str = None, script: str = None,
                                             case_path: str = None) -> dict:
    """!
    @brief Generate a dimensional canonical PICSLICE for square-duct Poiseuille flow.
    @param[in] output_path Path to write.
    @param[in] dims PICSLICE dimensions in face storage order (n1, n2).
    @param[in] params Normalized generator params.
    @param[in] target_grid Optional canonical target PICGRID for grid-aware sampling.
    @param[in] target_block Target block index when `target_grid` is provided.
    @param[in] target_face Target inlet face when `target_grid` is provided.
    @param[in] script Optional profile.gen-compatible script override.
    @param[in] case_path Current case.yml path used to anchor relative script overrides.
    @return Summary dictionary.
    """
    n1, n2 = tuple(dims)
    profilegen_script = _resolve_generator_script(script, case_path, "profile.gen")
    if not os.path.isfile(profilegen_script):
        raise ValueError(f"profile.gen script not found: {profilegen_script}")
    cmd = [
        sys.executable,
        profilegen_script,
        "square_duct_poiseuille",
        "--output",
        output_path,
        "--dims",
        str(n1),
        str(n2),
        "--bulk-velocity",
        str(float(params["bulk_velocity"])),
        "--n-terms",
        str(int(params["n_terms"])),
    ]
    if target_grid:
        cmd.extend([
            "--target-grid",
            target_grid,
            "--target-block",
            str(int(target_block)),
            f"--target-face={target_face}",
        ])
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        details = (result.stderr or result.stdout or "").strip()
        raise ValueError(f"profile.gen failed with exit code {result.returncode}. Details:\n{details}")
    try:
        summary = json.loads((result.stdout or "").strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise ValueError(f"profile.gen did not emit valid JSON summary. Output:\n{result.stdout}") from exc
    summary["dims"] = tuple(summary["dims"])
    return summary

def generate_field_slice_picslice(output_path: str, expected_dims: tuple, source: dict,
                                  target_grid: str, target_face: str, target_block: int,
                                  case_path: str) -> dict:
    """!
    @brief Invoke profile.gen to extract a field_slice PICSLICE artifact.
    @param[in] output_path Path to write.
    @param[in] expected_dims Expected PICSLICE dimensions.
    @param[in] source Normalized field_slice source mapping.
    @param[in] target_grid Target canonical PICGRID path.
    @param[in] target_face Target inlet face token.
    @param[in] target_block Target block index.
    @param[in] case_path Path to current case.yml for relative source resolution.
    @return Summary dictionary from profile.gen.
    """
    case_dir = os.path.dirname(os.path.abspath(case_path)) if case_path else os.getcwd()
    field_file = _resolve_case_relative_path(source["field_file"], case_dir)
    source_grid = _resolve_case_relative_path(source["grid_file"], case_dir)
    velocity_scale = _resolve_field_slice_velocity_scale(source, case_dir)
    profilegen_script = _resolve_generator_script(source.get("script"), case_path, "profile.gen")
    if not os.path.isfile(profilegen_script):
        raise ValueError(f"profile.gen script not found: {profilegen_script}")
    n1, n2 = tuple(expected_dims)
    slice_cfg = source["slice"]
    cmd = [
        sys.executable,
        profilegen_script,
        "field-slice",
        "--output",
        output_path,
        "--field-file",
        field_file,
        "--source-grid",
        source_grid,
        "--target-grid",
        target_grid,
        "--source-block",
        str(int(source.get("source_block", 0))),
        "--target-block",
        str(int(target_block)),
        f"--target-face={target_face}",
        "--orientation",
        slice_cfg["orientation"],
        "--normal-tolerance",
        str(float(slice_cfg["normal_tolerance"])),
        "--velocity-scale",
        str(float(velocity_scale)),
        "--expected-dims",
        str(n1),
        str(n2),
    ]
    if "face" in slice_cfg:
        cmd.append(f"--slice-face={slice_cfg['face']}")
    else:
        cmd.extend([
            "--slice-axis",
            slice_cfg["axis"],
            "--slice-index",
            str(int(slice_cfg["index"])),
            f"--slice-normal={slice_cfg['normal']}",
        ])
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        details = (result.stderr or result.stdout or "").strip()
        raise ValueError(f"profile.gen field-slice failed with exit code {result.returncode}. Details:\n{details}")
    try:
        summary = json.loads((result.stdout or "").strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise ValueError(f"profile.gen field-slice did not emit valid JSON summary. Output:\n{result.stdout}") from exc
    summary["dims"] = tuple(summary["dims"])
    return summary

def _resolve_case_relative_path(path_value: str, case_dir: str) -> str:
    """!
    @brief Resolve a path relative to the current case directory.
    @param[in] path_value Path from case.yml.
    @param[in] case_dir Current case directory.
    @return Absolute path.
    """
    if not isinstance(path_value, str) or not path_value.strip():
        raise ValueError("path value must be a non-empty string.")
    if os.path.isabs(path_value):
        return os.path.abspath(path_value)
    return os.path.abspath(os.path.join(case_dir, path_value))

def _resolve_field_slice_velocity_scale(source: dict, case_dir: str) -> float:
    """!
    @brief Resolve field_slice dimensional velocity scale.
    @param[in] source Normalized field_slice source mapping.
    @param[in] case_dir Current case directory.
    @return Positive velocity scale.
    """
    if source.get("velocity_scale") is not None:
        return float(source["velocity_scale"])
    source_case = _resolve_case_relative_path(source["source_case"], case_dir)
    source_case_cfg = read_yaml_file(source_case)
    try:
        velocity_scale = _to_float(
            source_case_cfg.get("properties", {}).get("scaling", {}).get("velocity_ref"),
            "source_case.properties.scaling.velocity_ref",
        )
    except AttributeError as exc:
        raise ValueError("source_case must contain properties.scaling.velocity_ref.") from exc
    if velocity_scale <= 0.0:
        raise ValueError("source_case.properties.scaling.velocity_ref must be positive.")
    return velocity_scale

def resolve_target_grid_for_field_slice(case_cfg: dict, case_path: str, run_dir: str) -> str:
    """!
    @brief Resolve the target canonical PICGRID path needed for field_slice normals.
    @param[in] case_cfg Parsed current case config.
    @param[in] case_path Current case.yml path.
    @param[in] run_dir Current run/precompute directory.
    @return Absolute target PICGRID path.
    """
    grid_cfg = case_cfg.get("grid", {}) or {}
    grid_mode = grid_cfg.get("mode")
    case_dir = os.path.dirname(os.path.abspath(case_path)) if case_path else os.getcwd()
    if grid_mode == "file":
        source_grid = grid_cfg.get("source_file")
        if not isinstance(source_grid, str) or not source_grid.strip():
            raise ValueError("grid.source_file is required for field_slice target-grid normals.")
        source_grid = _resolve_case_relative_path(source_grid, case_dir)
        if isinstance(grid_cfg.get("legacy_conversion"), dict) and run_dir:
            source_grid = convert_legacy_grid_with_gridgen(case_path, run_dir, grid_cfg, source_grid)
        return source_grid
    if grid_mode == "grid_gen":
        generator = grid_cfg.get("generator", {})
        output_file = generator.get("output_file", os.path.join("config", "grid.generated.picgrid"))
        candidate = output_file if os.path.isabs(output_file) else os.path.abspath(os.path.join(run_dir, output_file))
        if os.path.isfile(candidate):
            return candidate
        staged = os.path.join(run_dir, "config", "grid.run")
        if os.path.isfile(staged):
            return staged
        raise ValueError("field_slice requires the generated target PICGRID to exist before profile extraction.")
    raise ValueError(
        f"field_slice requires grid.mode 'file' or 'grid_gen' for target-grid normals; got '{grid_mode}'."
    )

def resolve_target_grid_for_generated_profile(case_cfg: dict, case_path: str, run_dir: str) -> str:
    """!
    @brief Resolve an optional target canonical PICGRID for generated profile sampling.
    @param[in] case_cfg Parsed current case config.
    @param[in] case_path Current case.yml path.
    @param[in] run_dir Current run/precompute directory.
    @return Absolute target PICGRID path, or None when no canonical grid is available yet.
    """
    grid_mode = (case_cfg.get("grid", {}) or {}).get("mode")
    if grid_mode == "programmatic_c":
        return None
    return resolve_target_grid_for_field_slice(case_cfg, case_path, run_dir)

def write_profile_info(config_dir: str, summaries: list) -> str:
    """!
    @brief Write a profile.info summary for generated inlet profiles.
    @param[in] config_dir Run/precompute config directory.
    @param[in] summaries Generated profile summaries.
    @return Path to profile.info.
    """
    info_path = os.path.join(config_dir, "profile.info")
    os.makedirs(config_dir, exist_ok=True)
    with open(info_path, "w") as fout:
        fout.write("# PICurv generated profile summary\n")
        fout.write(f"profile_count = {len(summaries)}\n\n")
        for idx, summary in enumerate(summaries):
            dims = summary.get("dims", (0, 0))
            fout.write(f"[profile_{idx}]\n")
            fout.write(f"generator = {summary.get('generator')}\n")
            fout.write(f"block = {summary.get('block')}\n")
            fout.write(f"face = {summary.get('face')}\n")
            fout.write(f"dimensions = {dims[0]} {dims[1]}\n")
            if summary.get("bulk_velocity") is not None:
                fout.write(f"bulk_velocity = {summary.get('bulk_velocity'):.16e}\n")
            if summary.get("n_terms") is not None:
                fout.write(f"n_terms = {summary.get('n_terms')}\n")
            fout.write(f"mean_speed = {summary.get('mean_speed'):.16e}\n")
            if "area_mean_speed" in summary:
                fout.write(f"area_mean_speed = {summary.get('area_mean_speed'):.16e}\n")
            if "discrete_mean_speed" in summary:
                fout.write(f"discrete_mean_speed = {summary.get('discrete_mean_speed'):.16e}\n")
            fout.write(f"min_speed = {summary.get('min_speed'):.16e}\n")
            fout.write(f"max_speed = {summary.get('max_speed'):.16e}\n")
            if summary.get("umax_over_ubulk") is not None:
                fout.write(f"umax_over_ubulk = {summary.get('umax_over_ubulk'):.16e}\n")
            for key in (
                "normalization",
                "sampling",
                "area_weighted_mean_before_normalization",
                "area_weighted_mean_after_normalization",
                "total_inlet_area",
                "face_area_min",
                "face_area_max",
                "source_field",
                "source_grid",
                "target_grid",
                "source_block",
                "target_block",
                "target_face",
                "source_slice",
                "orientation",
                "normal_tolerance",
                "normal_dot",
                "velocity_scale",
            ):
                if key in summary:
                    fout.write(f"{key} = {summary.get(key)}\n")
            fout.write(f"output_file = {summary.get('path')}\n\n")
    return info_path

def run_grid_generator(case_path: str, run_dir: str, grid_cfg: dict) -> str:
    """!
    @brief Runs generators/grid.gen to produce a PICGRID file for this run.
    @param[in] case_path Path to case.yml (used for relative path resolution).
    @param[in] run_dir Run directory path.
    @param[in] grid_cfg The grid config section from case.yml.
    @return Absolute path to generated dimensional PICGRID file.
    @throws ValueError on invalid config or generator failure.
    """
    generator = grid_cfg.get("generator", {})
    if not isinstance(generator, dict):
        raise ValueError("grid.generator must be a mapping when grid.mode is 'grid_gen'.")

    case_dir = os.path.dirname(os.path.abspath(case_path))
    gridgen_script = generator.get("script", os.path.join(GENERATORS_PATH, "grid.gen"))
    if not os.path.isabs(gridgen_script):
        gridgen_script = os.path.abspath(os.path.join(case_dir, gridgen_script))
    if not os.path.isfile(gridgen_script):
        raise ValueError(f"grid.gen script not found: {gridgen_script}")

    config_file = generator.get("config_file")
    if not config_file:
        raise ValueError("grid.generator.config_file is required when grid.mode is 'grid_gen'.")
    if not os.path.isabs(config_file):
        config_file = os.path.abspath(os.path.join(case_dir, config_file))
    if not os.path.isfile(config_file):
        raise ValueError(f"grid.generator.config_file not found: {config_file}")

    output_file = generator.get("output_file", os.path.join("config", "grid.generated.picgrid"))
    if not os.path.isabs(output_file):
        output_file = os.path.abspath(os.path.join(run_dir, output_file))
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    grid_type = generator.get("grid_type")
    cli_args = generator.get("cli_args", [])
    if cli_args is None:
        cli_args = []
    if not isinstance(cli_args, list):
        raise ValueError("grid.generator.cli_args must be a list of CLI tokens.")

    cmd = [sys.executable, gridgen_script, "-c", config_file]
    if grid_type:
        cmd.append(str(grid_type))
    cmd.extend([str(token) for token in cli_args])
    cmd.extend(["--output", output_file])

    vts_file = generator.get("vts_file")
    if vts_file:
        if not os.path.isabs(vts_file):
            vts_file = os.path.abspath(os.path.join(run_dir, vts_file))
        os.makedirs(os.path.dirname(vts_file), exist_ok=True)
        cmd.extend(["--vts", vts_file])

    stats_file = generator.get("stats_file")
    if stats_file:
        if not os.path.isabs(stats_file):
            stats_file = os.path.abspath(os.path.join(run_dir, stats_file))
        os.makedirs(os.path.dirname(stats_file), exist_ok=True)
        cmd.extend(["--stats-file", stats_file])

    print(f"[INFO] Grid generator command: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=case_dir, text=True, capture_output=True)
    if result.returncode != 0:
        stderr = (result.stderr or "").strip()
        stdout = (result.stdout or "").strip()
        details = stderr if stderr else stdout
        raise ValueError(
            f"grid.gen failed with exit code {result.returncode}. Details:\n{details}"
        )
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())

    if not os.path.isfile(output_file):
        raise ValueError(f"grid.gen did not produce expected output file: {output_file}")

    return output_file


def convert_legacy_grid_with_gridgen(case_path: str, run_dir: str, grid_cfg: dict, source_grid: str) -> str:
    """!
    @brief Optionally convert a legacy file-grid payload to canonical PICGRID using grid.gen.
    @details Activated only when `grid.legacy_conversion.enabled` is true in case.yml.
             The converted output remains dimensional; standard nondimensionalization still
             occurs via validate_and_nondimensionalize_picgrid().
    @param[in] case_path Path to case.yml (for relative path resolution).
    @param[in] run_dir Current run directory.
    @param[in] grid_cfg Grid section from case.yml.
    @param[in] source_grid Absolute or relative path to the original grid file.
    @return Grid path that should be fed into validate_and_nondimensionalize_picgrid().
    @throws ValueError on invalid converter settings or failed conversion.
    """
    legacy_cfg = grid_cfg.get("legacy_conversion")
    if not isinstance(legacy_cfg, dict):
        return source_grid

    enabled = legacy_cfg.get("enabled", True)
    if enabled is False:
        return source_grid
    if not isinstance(enabled, bool):
        raise ValueError("grid.legacy_conversion.enabled must be a boolean.")

    raw_format = str(legacy_cfg.get("format", "legacy1d")).strip().lower()
    format_aliases = {
        "legacy1d": "legacy1d",
        "legacy_1d": "legacy1d",
        "les_flat_1d": "legacy1d",
        "les-flat-1d": "legacy1d",
    }
    command = format_aliases.get(raw_format)
    if command is None:
        raise ValueError(
            "grid.legacy_conversion.format must be one of "
            "['legacy1d', 'legacy_1d', 'les_flat_1d', 'les-flat-1d']."
        )

    case_dir = os.path.dirname(os.path.abspath(case_path))
    gridgen_script = legacy_cfg.get("script", os.path.join(GENERATORS_PATH, "grid.gen"))
    if not isinstance(gridgen_script, str) or not gridgen_script.strip():
        raise ValueError("grid.legacy_conversion.script must be a non-empty string when provided.")
    if not os.path.isabs(gridgen_script):
        gridgen_script = os.path.abspath(os.path.join(case_dir, gridgen_script))
    if not os.path.isfile(gridgen_script):
        raise ValueError(f"grid.legacy_conversion.script not found: {gridgen_script}")

    output_file = legacy_cfg.get("output_file", os.path.join("config", "grid.converted.picgrid"))
    if not isinstance(output_file, str) or not output_file.strip():
        raise ValueError("grid.legacy_conversion.output_file must be a non-empty string when provided.")
    if not os.path.isabs(output_file):
        output_file = os.path.abspath(os.path.join(run_dir, output_file))
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    axis_columns = legacy_cfg.get("axis_columns", [0, 1, 2])
    if not isinstance(axis_columns, list) or len(axis_columns) != 3:
        raise ValueError("grid.legacy_conversion.axis_columns must be a 3-item list of non-negative integers.")
    try:
        axis_columns = [int(v) for v in axis_columns]
    except (TypeError, ValueError) as exc:
        raise ValueError("grid.legacy_conversion.axis_columns must contain integers.") from exc
    if any(v < 0 for v in axis_columns):
        raise ValueError("grid.legacy_conversion.axis_columns values must be >= 0.")

    strict_trailing = legacy_cfg.get("strict_trailing", True)
    if not isinstance(strict_trailing, bool):
        raise ValueError("grid.legacy_conversion.strict_trailing must be a boolean.")

    cli_args = legacy_cfg.get("cli_args", [])
    if cli_args is None:
        cli_args = []
    if not isinstance(cli_args, list):
        raise ValueError("grid.legacy_conversion.cli_args must be a list of CLI tokens.")

    cmd = [
        sys.executable,
        gridgen_script,
        command,
        "--input",
        source_grid,
        "--output",
        output_file,
        "--axis-columns",
        str(axis_columns[0]),
        str(axis_columns[1]),
        str(axis_columns[2]),
        "--no-write-vtk",
    ]
    if strict_trailing:
        cmd.append("--strict-trailing")
    else:
        cmd.append("--allow-trailing")
    cmd.extend(str(token) for token in cli_args)

    print(f"[INFO] Legacy grid conversion command: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=case_dir, text=True, capture_output=True)
    if result.returncode != 0:
        stderr = (result.stderr or "").strip()
        stdout = (result.stdout or "").strip()
        details = stderr if stderr else stdout
        raise ValueError(
            f"legacy grid conversion failed with exit code {result.returncode}. Details:\n{details}"
        )
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())
    if not os.path.isfile(output_file):
        raise ValueError(f"legacy grid conversion did not produce expected output file: {output_file}")

    return output_file

BC_FACE_MAP = {
    "-xi": "-Xi",
    "+xi": "+Xi",
    "-eta": "-Eta",
    "+eta": "+Eta",
    "-zeta": "-Zeta",
    "+zeta": "+Zeta",
}

# Canonical BC selector maps. Add new values only with matching C parser/factory support.
BC_TYPE_MAP = {
    "wall": "WALL",
    "symmetry": "SYMMETRY",
    "inlet": "INLET",
    "outlet": "OUTLET",
    "periodic": "PERIODIC",
}

BC_HANDLER_SPECS = {
    # Only handlers that are implemented end-to-end in current C path are allowed.
    "noslip": {
        "types": {"WALL"},
        "required_params": set(),
        "optional_params": set(),
    },
    "constant_velocity": {
        "types": {"INLET"},
        "required_params": {"vx", "vy", "vz"},
        "optional_params": set(),
    },
    "conservation": {
        "types": {"OUTLET"},
        "required_params": set(),
        "optional_params": set(),
    },
    "parabolic": {
        "types": {"INLET"},
        "required_params": {"v_max"},
        "optional_params": set(),
    },
    "prescribed_flow": {
        "types": {"INLET"},
        "required_params": {"source"},
        "optional_params": set(),
    },
    "geometric": {
        "types": {"PERIODIC"},
        "required_params": set(),
        "optional_params": set(),
    },
    "constant_flux": {
        "types": {"PERIODIC"},
        "required_params": {"target_flux"},
        "optional_params": {"apply_trim"},
    },
}

_NUMERIC_BC_PARAMS = {"vx", "vy", "vz", "v_max", "target_flux"}
_BOOL_BC_PARAMS = {"apply_trim"}

def _normalize_prescribed_flow_source(source, field_name: str) -> dict:
    """!
    @brief Validate the structured source block for prescribed_flow BCs.
    @param[in] source Source mapping from case.yml.
    @param[in] field_name Human-readable YAML path for diagnostics.
    @return Normalized source mapping.
    @throws ValueError on invalid source contract.
    """
    if not isinstance(source, dict):
        raise ValueError(f"{field_name} must be a mapping with type: file, generated, or field_slice.")
    source_type = str(source.get("type", "")).strip().lower()
    if source_type == "file":
        path = source.get("path")
        if not isinstance(path, str) or not path.strip():
            raise ValueError(f"{field_name}.path must be a non-empty file path.")
        unknown = sorted(set(source.keys()) - {"type", "path"})
        if unknown:
            raise ValueError(f"Unknown keys in {field_name}: {unknown}. Allowed: ['path', 'type'].")
        return {"type": "file", "path": path.strip()}

    if source_type == "generated":
        generator = str(source.get("generator", "")).strip().lower()
        if generator not in GENERATED_PROFILE_GENERATORS:
            raise ValueError(
                f"{field_name}.generator must be one of {sorted(GENERATED_PROFILE_GENERATORS)} "
                f"(got '{source.get('generator')}')."
            )
        unknown = sorted(set(source.keys()) - {"type", "generator", "script", "output_file", "params"})
        if unknown:
            raise ValueError(
                f"Unknown keys in {field_name}: {unknown}. "
                "Allowed: ['generator', 'output_file', 'params', 'script', 'type']."
            )
        normalized = {
            "type": "generated",
            "generator": generator,
            "params": _normalize_square_duct_poiseuille_params(source.get("params", {}), field_name),
        }
        output_file = source.get("output_file")
        if output_file is not None:
            if not isinstance(output_file, str) or not output_file.strip():
                raise ValueError(f"{field_name}.output_file must be a non-empty path when provided.")
            normalized["output_file"] = output_file.strip()
        script = source.get("script")
        if script is not None:
            if not isinstance(script, str) or not script.strip():
                raise ValueError(f"{field_name}.script must be a non-empty path when provided.")
            normalized["script"] = script.strip()
        return normalized

    if source_type == "field_slice":
        return _normalize_field_slice_source(source, field_name)

    raise ValueError(f"{field_name}.type must be 'file', 'generated', or 'field_slice'.")

def _bc_profile_expected_dims(face: str, block_dims: tuple) -> tuple:
    """!
    @brief Return expected PICSLICE dimensions for a face and block node dimensions.
    @param[in] face Canonical BC face token.
    @param[in] block_dims (IM, JM, KM) node counts.
    @return (n1, n2) dimensions in profile storage order.
    """
    im, jm, km = block_dims
    if min(im, jm, km) < 2:
        raise ValueError(
            f"Block dimensions {block_dims} are too small for an inlet profile; each axis needs at least 2 nodes."
        )
    if face in {"-Xi", "+Xi"}:
        return (km - 1, jm - 1)
    if face in {"-Eta", "+Eta"}:
        return (km - 1, im - 1)
    if face in {"-Zeta", "+Zeta"}:
        return (jm - 1, im - 1)
    raise ValueError(f"Unsupported face '{face}' for prescribed_flow profile dimensions.")

def resolve_grid_block_dimensions_for_profiles(case_cfg: dict, case_path: str, run_dir: str = None) -> list:
    """!
    @brief Resolve per-block node dimensions for prescribed inlet profile validation.
    @param[in] case_cfg Parsed case.yml configuration.
    @param[in] case_path Path to case.yml for relative path resolution.
    @param[in] run_dir Current run directory, used for optional generated grid outputs.
    @return List of (IM, JM, KM) node-count tuples.
    @throws ValueError when dimensions cannot be resolved.
    """
    num_blocks = int(case_cfg.get('models', {}).get('domain', {}).get('blocks', 1))
    grid_cfg = case_cfg.get("grid", {})
    grid_mode = grid_cfg.get("mode")
    case_dir = os.path.dirname(os.path.abspath(case_path)) if case_path else os.getcwd()

    if grid_mode == "programmatic_c":
        settings = grid_cfg.get("programmatic_settings", {})
        dims_by_axis = []
        for key in ("im", "jm", "km"):
            raw = settings.get(key)
            if raw is None:
                raise ValueError(f"grid.programmatic_settings.{key} is required for prescribed_flow profiles.")
            if isinstance(raw, list):
                if len(raw) != num_blocks:
                    raise ValueError(
                        f"grid.programmatic_settings.{key} has {len(raw)} entries, expected {num_blocks} blocks."
                    )
                values = raw
            else:
                values = [raw] * num_blocks
            try:
                values = [int(v) + 1 for v in values]
            except (TypeError, ValueError):
                raise ValueError(f"grid.programmatic_settings.{key} values must be positive integer cell counts.")
            if any(v <= 1 for v in values):
                raise ValueError(f"grid.programmatic_settings.{key} values must be positive integer cell counts.")
            dims_by_axis.append(values)
        return list(zip(dims_by_axis[0], dims_by_axis[1], dims_by_axis[2]))

    if grid_mode == "file":
        source_grid = grid_cfg.get("source_file")
        if not isinstance(source_grid, str) or not source_grid.strip():
            raise ValueError("grid.source_file is required for file-grid prescribed_flow profile validation.")
        if not os.path.isabs(source_grid):
            source_grid = os.path.abspath(os.path.join(case_dir, source_grid))
        if isinstance(grid_cfg.get("legacy_conversion"), dict) and run_dir:
            source_grid = convert_legacy_grid_with_gridgen(case_path, run_dir, grid_cfg, source_grid)
        return read_picgrid_header_dimensions(source_grid, expected_nblk=num_blocks)

    if grid_mode == "grid_gen":
        generator = grid_cfg.get("generator", {})
        output_file = generator.get("output_file", os.path.join("config", "grid.generated.picgrid"))
        candidates = []
        if run_dir:
            candidates.append(output_file if os.path.isabs(output_file) else os.path.abspath(os.path.join(run_dir, output_file)))
            candidates.append(os.path.join(run_dir, "config", "grid.run"))
        for candidate in candidates:
            if os.path.isfile(candidate):
                return read_picgrid_header_dimensions(candidate, expected_nblk=num_blocks)
        raise ValueError(
            "prescribed_flow profile dimension validation for grid.mode='grid_gen' requires an existing generated "
            "PICGRID output. Run or stage the grid first, or use grid.mode='file' with the generated .picgrid."
        )

    raise ValueError(f"Unsupported grid.mode '{grid_mode}' for prescribed_flow profile validation.")

def materialize_generated_prescribed_flow_profiles(run_dir: str, case_cfg: dict, case_path: str,
                                                   profile_grid_dims: list = None) -> list:
    """!
    @brief Generate dimensional PICSLICE artifacts for generated/field_slice prescribed_flow sources.
    @param[in] run_dir Run/precompute directory root.
    @param[in] case_cfg Parsed case.yml.
    @param[in] case_path Path to case.yml for relative grid/source resolution.
    @param[in] profile_grid_dims Optional pre-resolved block node dimensions.
    @return List of generated profile summaries.
    """
    prepared_blocks = validate_and_prepare_boundary_conditions(case_cfg)
    if not any(
        bc.get("handler") == "prescribed_flow"
        and ((bc.get("params") or {}).get("source") or {}).get("type") in {"generated", "field_slice"}
        for block in prepared_blocks for bc in block
    ):
        return []
    if profile_grid_dims is None:
        profile_grid_dims = resolve_grid_block_dimensions_for_profiles(case_cfg, case_path, run_dir)

    config_dir = os.path.join(run_dir, "config")
    target_grid = None
    generated_target_grid = None
    summaries = []
    for block_idx, block in enumerate(prepared_blocks):
        for bc in block:
            if bc.get("handler") != "prescribed_flow":
                continue
            source = (bc.get("params") or {}).get("source", {})
            if source.get("type") not in {"generated", "field_slice"}:
                continue
            face = bc["face"]
            dims = _bc_profile_expected_dims(face, profile_grid_dims[block_idx])
            face_token = _face_artifact_token(face)
            suffix = "generated" if source.get("type") == "generated" else "sliced"
            default_output = os.path.join(
                "config", f"inlet_profile_block{block_idx}_{face_token}.{suffix}.picslice"
            )
            output_path = _resolve_run_artifact_path(
                run_dir,
                source.get("output_file"),
                default_output,
                default_to_config_dir=True,
            )
            if source.get("type") == "generated" and source["generator"] == "square_duct_poiseuille":
                if generated_target_grid is None:
                    generated_target_grid = resolve_target_grid_for_generated_profile(case_cfg, case_path, run_dir)
                summary = generate_square_duct_poiseuille_picslice(
                    output_path,
                    dims,
                    source["params"],
                    target_grid=generated_target_grid,
                    target_block=block_idx,
                    target_face=face,
                    script=source.get("script"),
                    case_path=case_path,
                )
            elif source.get("type") == "generated":
                raise ValueError(f"Unsupported generated profile generator '{source['generator']}'.")
            else:
                if target_grid is None:
                    target_grid = resolve_target_grid_for_field_slice(case_cfg, case_path, run_dir)
                summary = generate_field_slice_picslice(
                    output_path,
                    dims,
                    source,
                    target_grid,
                    face,
                    block_idx,
                    case_path,
                )
            summary.update({"block": block_idx, "face": face})
            summaries.append(summary)
            print(
                f"[SUCCESS] Materialized prescribed_flow profile for block {block_idx}, face {face}: "
                f"{os.path.relpath(output_path)} dims={summary['dims']}"
            )

    if summaries:
        info_path = write_profile_info(config_dir, summaries)
        print(f"[SUCCESS] Wrote generated profile summary: {os.path.relpath(info_path)}")
    return summaries

def _to_float(value, field_name: str) -> float:
    """!
    @brief Convert a YAML scalar to float with a clear error message.
    @param[in] value Argument passed to `_to_float()`.
    @param[in] field_name Argument passed to `_to_float()`.
    @return Value returned by `_to_float()`.
    """
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"'{field_name}' must be numeric (got {value!r}).")

def _to_bool(value, field_name: str) -> bool:
    """!
    @brief Convert a YAML scalar/string to bool with a clear error message.
    @param[in] value Argument passed to `_to_bool()`.
    @param[in] field_name Argument passed to `_to_bool()`.
    @return Value returned by `_to_bool()`.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        raw = value.strip().lower()
        if raw in {"true", "1", "yes"}:
            return True
        if raw in {"false", "0", "no"}:
            return False
    raise ValueError(f"'{field_name}' must be boolean (got {value!r}).")

def normalize_boundary_conditions_layout(all_blocks_bcs, num_blocks: int):
    """!
    @brief Normalize boundary_conditions to list-of-lists form and validate block count.
    @param[in] all_blocks_bcs Argument passed to `normalize_boundary_conditions_layout()`.
    @param[in] num_blocks Argument passed to `normalize_boundary_conditions_layout()`.
    @return Value returned by `normalize_boundary_conditions_layout()`.
    """
    if not all_blocks_bcs:
        raise ValueError("The 'boundary_conditions' section in case.yml is empty.")

    is_simple_list = isinstance(all_blocks_bcs[0], dict)
    if num_blocks == 1 and is_simple_list:
        all_blocks_bcs = [all_blocks_bcs]
    elif is_simple_list and num_blocks > 1:
        raise ValueError(
            f"case.yml declares {num_blocks} blocks but boundary_conditions is a single face-list. "
            "Use a list-of-lists, one inner list per block."
        )

    if len(all_blocks_bcs) != num_blocks:
        raise ValueError(
            f"Mismatch: case.yml declares {num_blocks} block(s) but found {len(all_blocks_bcs)} BC definitions."
        )
    return all_blocks_bcs

def validate_and_prepare_boundary_conditions(case_cfg: dict):
    """!
    @brief Validate BC entries against currently supported C-side handlers/types and
    @details return normalized entries ready for bcs.run generation.
    @param[in] case_cfg Argument passed to `validate_and_prepare_boundary_conditions()`.
    @return Value returned by `validate_and_prepare_boundary_conditions()`.
    """
    num_blocks = int(case_cfg.get('models', {}).get('domain', {}).get('blocks', 1))
    scales = case_cfg.get('properties', {}).get('scaling', {})
    L_ref = _to_float(scales.get('length_ref'), "properties.scaling.length_ref")
    U_ref = _to_float(scales.get('velocity_ref'), "properties.scaling.velocity_ref")
    if U_ref == 0.0:
        raise ValueError("properties.scaling.velocity_ref must be non-zero for non-dimensionalization.")
    if L_ref == 0.0:
        raise ValueError("properties.scaling.length_ref must be non-zero for non-dimensionalization.")

    all_blocks_bcs = normalize_boundary_conditions_layout(case_cfg.get('boundary_conditions', []), num_blocks)
    prepared_blocks = []

    expected_faces = {"-Xi", "+Xi", "-Eta", "+Eta", "-Zeta", "+Zeta"}
    axis_pairs = [("-Xi", "+Xi"), ("-Eta", "+Eta"), ("-Zeta", "+Zeta")]

    for bi, block_bcs in enumerate(all_blocks_bcs):
        if not isinstance(block_bcs, list):
            raise ValueError(f"boundary_conditions[{bi}] must be a list of face configs.")

        prepared_block = []
        seen_faces = {}

        for idx, bc in enumerate(block_bcs):
            if not isinstance(bc, dict):
                raise ValueError(f"boundary_conditions[{bi}][{idx}] must be a mapping.")

            for req in ("face", "type", "handler"):
                if req not in bc:
                    raise ValueError(f"boundary_conditions[{bi}][{idx}] missing required key '{req}'.")

            face_raw = str(bc["face"]).strip()
            face_key = face_raw.lower()
            face = BC_FACE_MAP.get(face_key)
            if face is None:
                raise ValueError(
                    f"Unsupported BC face '{face_raw}' at boundary_conditions[{bi}][{idx}]. "
                    f"Supported: {sorted(expected_faces)}."
                )
            if face in seen_faces:
                raise ValueError(f"Duplicate face '{face}' in boundary_conditions[{bi}] (entries {seen_faces[face]} and {idx}).")
            seen_faces[face] = idx

            bc_type_raw = str(bc["type"]).strip()
            bc_type = BC_TYPE_MAP.get(bc_type_raw.lower())
            if bc_type is None:
                raise ValueError(
                    f"Unsupported BC type '{bc_type_raw}' for face {face} in block {bi}. "
                    f"Supported: {sorted(set(BC_TYPE_MAP.values()))}."
                )

            handler = str(bc["handler"]).strip().lower()
            handler_spec = BC_HANDLER_SPECS.get(handler)
            if handler_spec is None:
                raise ValueError(
                    f"Unsupported BC handler '{bc['handler']}' for face {face} in block {bi}. "
                    f"Supported now: {sorted(BC_HANDLER_SPECS.keys())}."
                )
            if bc_type not in handler_spec["types"]:
                raise ValueError(
                    f"Invalid BC combination on block {bi}, face {face}: type '{bc_type}' cannot use handler '{handler}'."
                )

            params = bc.get("params", {})
            if params is None:
                params = {}
            if not isinstance(params, dict):
                raise ValueError(f"'params' for block {bi}, face {face} must be a mapping.")

            # Reject unsupported older structured keys explicitly.
            if "vector" in params or "velocity" in params:
                raise ValueError(
                    f"Unsupported older params key ('vector'/'velocity') found on block {bi}, face {face}. "
                    "Use scalar keys 'vx', 'vy', 'vz'."
                )

            required = handler_spec["required_params"]
            optional = handler_spec["optional_params"]
            allowed = required | optional

            missing = sorted(required - set(params.keys()))
            if missing:
                raise ValueError(
                    f"Missing required params for handler '{handler}' on block {bi}, face {face}: {missing}."
                )
            unknown = sorted(set(params.keys()) - allowed)
            if unknown:
                raise ValueError(
                    f"Unknown params for handler '{handler}' on block {bi}, face {face}: {unknown}. "
                    f"Allowed: {sorted(allowed)}."
                )

            converted_params = {}
            for key, value in params.items():
                if key in _NUMERIC_BC_PARAMS:
                    numeric = _to_float(value, f"boundary_conditions[{bi}][{idx}].params.{key}")
                    if key in {"vx", "vy", "vz", "v_max"}:
                        converted_params[key] = numeric / U_ref
                    elif key == "target_flux":
                        converted_params[key] = numeric / (U_ref * (L_ref ** 2))
                elif key in _BOOL_BC_PARAMS:
                    converted_params[key] = _to_bool(value, f"boundary_conditions[{bi}][{idx}].params.{key}")
                elif handler == "prescribed_flow" and key == "source":
                    converted_params[key] = _normalize_prescribed_flow_source(
                        value, f"boundary_conditions[{bi}][{idx}].params.source"
                    )
                else:
                    # Defensive fallback; should not happen due unknown-key gate above.
                    converted_params[key] = value

            prepared_block.append({
                "face": face,
                "type": bc_type,
                "handler": handler,
                "params": converted_params,
            })

        missing_faces = sorted(expected_faces - set(seen_faces.keys()))
        if missing_faces:
            raise ValueError(
                f"boundary_conditions[{bi}] is incomplete. Missing faces: {missing_faces}. "
                "Provide all six faces explicitly."
            )

        # Pairwise periodic consistency checks.
        face_map = {entry["face"]: entry for entry in prepared_block}
        for neg_face, pos_face in axis_pairs:
            neg = face_map[neg_face]
            pos = face_map[pos_face]
            neg_periodic = (neg["type"] == "PERIODIC")
            pos_periodic = (pos["type"] == "PERIODIC")
            if neg_periodic != pos_periodic:
                raise ValueError(
                    f"Inconsistent periodicity in block {bi}: {neg_face} and {pos_face} must both be PERIODIC or neither."
                )

            driven_handlers = {"constant_flux"}
            if (neg["handler"] in driven_handlers) or (pos["handler"] in driven_handlers):
                if neg["handler"] != pos["handler"]:
                    raise ValueError(
                        f"In block {bi}, driven periodic handlers on {neg_face}/{pos_face} must match exactly."
                    )
                if not (neg_periodic and pos_periodic):
                    raise ValueError(
                        f"In block {bi}, driven periodic handler '{neg['handler']}' requires PERIODIC type on both faces."
                    )

        prepared_blocks.append(prepared_block)

    return prepared_blocks


def _schema_path_text(path: tuple) -> str:
    """!
    @brief Render an internal schema path tuple as a user-facing YAML path.
    @param[in] path Internal path tuple.
    @return Dotted YAML path.
    """
    return ".".join(part for part in path if part != "[]") or "<root>"


def _lookup_allowed_schema_keys(schema: dict, path: tuple):
    """!
    @brief Return allowed keys for a path, honoring '*' dynamic mapping entries.
    @param[in] schema Role schema mapping.
    @param[in] path Internal path tuple.
    @return Allowed key set, None for free-form mappings, or False when path is not schema-checked.
    """
    if path in schema:
        return schema[path]
    for idx, part in enumerate(path):
        if part == "[]":
            continue
        candidate = path[:idx] + ("*",) + path[idx + 1:]
        if candidate in schema:
            return schema[candidate]
    return False


def _schema_key_hint(schema: dict, path: tuple, key: str, allowed: set) -> str:
    """!
    @brief Build a concise typo or hierarchy hint for an unsupported YAML key.
    @param[in] schema Role schema mapping.
    @param[in] path Current internal YAML path tuple.
    @param[in] key Unsupported YAML key.
    @param[in] allowed Allowed keys at the current path.
    @return Optional hint string.
    """
    hints = []
    allowed_strings = sorted(str(item) for item in allowed)
    lower_matches = [item for item in allowed_strings if item.lower() == key.lower()]
    close_matches = lower_matches or difflib.get_close_matches(key, allowed_strings, n=1, cutoff=0.80)
    if close_matches:
        hints.append(f"Did you mean '{close_matches[0]}'?")

    valid_paths = []
    for schema_path, schema_allowed in schema.items():
        if schema_path == path or not schema_allowed:
            continue
        if key in schema_allowed:
            valid_paths.append(_schema_path_text(schema_path))
    if valid_paths:
        hints.append(f"This key is valid at: {', '.join(sorted(valid_paths))}.")

    return " ".join(hints)


def _validate_yaml_schema_keys(cfg, schema: dict, file_path: str, errors: list, path: tuple = ()) -> None:
    """!
    @brief Reject unsupported YAML keys before they can be silently ignored by staging.
    @param[in] cfg Parsed YAML node.
    @param[in] schema Role schema mapping.
    @param[in] file_path Source file path for diagnostics.
    @param[in,out] errors Validation error accumulator.
    @param[in] path Current internal YAML path tuple.
    """
    if isinstance(cfg, dict):
        allowed = _lookup_allowed_schema_keys(schema, path)
        if allowed is not False and allowed is not None:
            unknown = sorted(str(key) for key in cfg.keys() if key not in allowed)
            for key in unknown:
                hint = _schema_key_hint(schema, path, key, allowed)
                hint_text = f" {hint}" if hint else ""
                errors.append(
                    f"  {file_path}: unsupported key at {_schema_path_text(path)}: '{key}'. "
                    f"Allowed keys: {sorted(allowed)}.{hint_text}"
                )
        if allowed is None:
            return
        for key, value in cfg.items():
            _validate_yaml_schema_keys(value, schema, file_path, errors, path + (key,))
    elif isinstance(cfg, list):
        for item in cfg:
            _validate_yaml_schema_keys(item, schema, file_path, errors, path + ("[]",))


_CASE_SCHEMA = {
    (): {
        "properties", "run_control", "grid", "models", "boundary_conditions", "solver_parameters",
    },
    ("run_control",): {"start_step", "total_steps", "dt_physical"},
    ("properties",): {"scaling", "fluid", "initial_conditions"},
    ("properties", "scaling"): {"length_ref", "velocity_ref"},
    ("properties", "fluid"): {"density", "viscosity"},
    ("properties", "initial_conditions"): {
        "mode", "generator", "params", "field", "source_file",
        "u_physical", "v_physical", "w_physical", "peak_velocity_physical",
        "velocity_physical", "flow_direction",
    },
    ("properties", "initial_conditions", "params"): None,
    ("grid",): {
        "mode", "source_file", "programmatic_settings", "generator", "legacy_conversion",
        "da_processors_x", "da_processors_y", "da_processors_z",
    },
    ("grid", "programmatic_settings"): {
        "im", "jm", "km", "xMins", "xMaxs", "yMins", "yMaxs", "zMins", "zMaxs",
        "rxs", "rys", "rzs", "cgrids",
        "da_processors_x", "da_processors_y", "da_processors_z",
    },
    ("grid", "generator"): {
        "script", "config_file", "grid_type", "cli_args", "output_file", "stats_file", "vts_file",
        # Retained so existing warning behavior for typo-prone hyphen keys is not bypassed.
        "config-file", "grid-type", "output-file", "stats-file", "vts-file",
    },
    ("grid", "legacy_conversion"): {
        "enabled", "format", "script", "output_file", "axis_columns", "strict_trailing", "cli_args",
    },
    ("models",): {"domain", "physics", "statistics"},
    ("models", "domain"): {"blocks"},
    ("models", "physics"): {"dimensionality", "fsi", "particles", "turbulence"},
    ("models", "physics", "fsi"): {"immersed", "moving_fsi"},
    ("models", "physics", "particles"): {"count", "init_mode", "restart_mode", "point_source"},
    ("models", "physics", "particles", "point_source"): {"x", "y", "z"},
    ("models", "physics", "turbulence"): {"les", "rans", "wall_function"},
    ("models", "physics", "turbulence", "les"): {
        "enabled", "model", "constant_cs", "max_cs", "dynamic_frequency", "test_filter",
    },
    ("models", "physics", "turbulence", "rans"): {"enabled", "model"},
    ("models", "physics", "turbulence", "wall_function"): {"enabled", "model", "roughness_height"},
    ("models", "statistics"): {"time_averaging"},
    ("boundary_conditions", "[]"): {"face", "type", "handler", "params"},
    ("boundary_conditions", "[]", "[]"): {"face", "type", "handler", "params"},
    ("boundary_conditions", "[]", "params"): None,
    ("boundary_conditions", "[]", "[]", "params"): None,
    ("solver_parameters",): None,
}


_SOLVER_SCHEMA = {
    (): {
        "operation_mode", "strategy", "tolerances", "momentum_solver", "poisson_solver",
        "pressure_solver", "interpolation", "petsc_passthrough_options", "verification",
        "scalar_transport", "solution_convergence",
    },
    ("operation_mode",): {"eulerian_field_source", "analytical_type", "uniform_flow"},
    ("operation_mode", "uniform_flow"): {"u", "v", "w"},
    ("strategy",): {"momentum_solver", "central_diff"},
    ("tolerances",): {
        "max_iterations", "absolute_tol", "relative_tol", "step_tol",
        "residual_absolute_tol", "residual_relative_tol",
    },
    ("momentum_solver",): {"type", "dual_time_picard_jameson_rk", "dual_time_picard_rk4"},
    ("momentum_solver", "dual_time_picard_jameson_rk"): {
        "max_pseudo_steps", "absolute_tol", "relative_tol", "step_tol", "pseudo_cfl",
        "jameson_residual_noise_allowance_factor", "rk4_residual_noise_allowance_factor",
    },
    ("momentum_solver", "dual_time_picard_jameson_rk", "pseudo_cfl"): {
        "initial", "minimum", "maximum", "growth_factor", "reduction_factor",
    },
    ("momentum_solver", "dual_time_picard_rk4"): {
        "max_pseudo_steps", "absolute_tol", "relative_tol", "step_tol", "pseudo_cfl",
        "jameson_residual_noise_allowance_factor", "rk4_residual_noise_allowance_factor",
    },
    ("momentum_solver", "dual_time_picard_rk4", "pseudo_cfl"): {
        "initial", "minimum", "maximum", "growth_factor", "reduction_factor",
    },
    ("poisson_solver",): {
        "method", "absolute_tolerance", "relative_tolerance", "max_iterations", "tolerance",
        "gmres", "preconditioner", "multigrid",
    },
    ("pressure_solver",): {
        "method", "absolute_tolerance", "relative_tolerance", "max_iterations", "tolerance",
        "gmres", "preconditioner", "multigrid",
    },
    ("poisson_solver", "gmres"): {"restart"},
    ("pressure_solver", "gmres"): {"restart"},
    ("poisson_solver", "preconditioner"): {"type"},
    ("pressure_solver", "preconditioner"): {"type"},
    ("poisson_solver", "multigrid"): {
        "levels", "pre_sweeps", "post_sweeps", "cycle", "mode", "semi_coarsening", "level_solvers",
    },
    ("pressure_solver", "multigrid"): {
        "levels", "pre_sweeps", "post_sweeps", "cycle", "mode", "semi_coarsening", "level_solvers",
    },
    ("poisson_solver", "multigrid", "semi_coarsening"): {"i", "j", "k"},
    ("pressure_solver", "multigrid", "semi_coarsening"): {"i", "j", "k"},
    ("poisson_solver", "multigrid", "level_solvers", "*"): {
        "method", "preconditioner", "ksp_type", "pc_type", "max_it", "rtol", "atol",
    },
    ("pressure_solver", "multigrid", "level_solvers", "*"): {
        "method", "preconditioner", "ksp_type", "pc_type", "max_it", "rtol", "atol",
    },
    ("interpolation",): {"method"},
    ("petsc_passthrough_options",): None,
    ("verification",): {"sources"},
    ("verification", "sources"): {"diffusivity", "scalar"},
    ("verification", "sources", "diffusivity"): {"mode", "profile", "gamma0", "slope_x"},
    ("verification", "sources", "scalar"): {
        "mode", "profile", "value", "phi0", "slope_x", "amplitude", "kx", "ky", "kz",
    },
    ("scalar_transport",): {"schmidt_number", "turbulent_schmidt_number"},
    ("solution_convergence",): {"enabled", "mode", "periodic_deterministic", "statistical_steady"},
    ("solution_convergence", "periodic_deterministic"): {"period_steps"},
    ("solution_convergence", "statistical_steady"): {"window_steps"},
}


_MONITOR_SCHEMA = {
    (): {"logging", "profiling", "diagnostics", "io", "solver_monitoring"},
    ("logging",): {"verbosity", "enabled_functions"},
    ("profiling",): {"timestep_output", "final_summary"},
    ("profiling", "timestep_output"): {"mode", "functions", "file"},
    ("profiling", "final_summary"): {"enabled"},
    ("diagnostics",): {"petsc", "runtime_memory_log"},
    ("diagnostics", "petsc"): {
        "malloc_debug", "malloc_test", "malloc_dump", "malloc_view", "malloc_view_threshold",
        "memory_view", "log_view", "log_view_memory", "log_all", "log_trace",
        "objects_dump", "options_left",
    },
    ("diagnostics", "runtime_memory_log"): {"enabled", "file"},
    ("io",): {
        "data_output_frequency", "particle_console_output_frequency", "particle_log_interval",
        "directories",
    },
    ("io", "directories"): {"output", "restart", "log", "eulerian_subdir", "particle_subdir"},
    ("solver_monitoring",): {"poisson", "petsc_passthrough_options"},
    ("solver_monitoring", "poisson"): {"pic_true_residual", "true_residual", "converged_reason", "view"},
    ("solver_monitoring", "petsc_passthrough_options"): None,
}


_POST_SCHEMA = {
    (): {
        "run_control", "source_data", "global_operations", "eulerian_pipeline",
        "lagrangian_pipeline", "statistics_pipeline", "statistics_output_prefix", "io",
    },
    ("run_control",): {
        "start_step", "end_step", "step_interval", "startTime", "endTime", "timeStep",
    },
    ("source_data",): {"directory", "input_extensions"},
    ("source_data", "input_extensions"): {"eulerian", "particle"},
    ("global_operations",): {"dimensionalize"},
    ("eulerian_pipeline", "[]"): {"task", "input_field", "output_field", "field", "reference_point"},
    ("lagrangian_pipeline", "[]"): {"task", "input_field", "output_field"},
    ("statistics_pipeline",): {"output_prefix", "tasks"},
    ("statistics_pipeline", "tasks", "[]"): {"task"},
    ("io",): {
        "output_directory", "output_filename_prefix", "particle_filename_prefix", "output_particles",
        "particle_subsampling_frequency", "input_extensions", "eulerian_fields_averaged",
        "eulerian_fields", "particle_fields",
    },
    ("io", "input_extensions"): {"eulerian", "particle"},
}


_CLUSTER_SCHEMA = {
    (): {"scheduler", "resources", "notifications", "execution"},
    ("scheduler",): {"type"},
    ("resources",): {"account", "partition", "nodes", "ntasks_per_node", "mem", "time"},
    ("notifications",): {"mail_user", "mail_type"},
    ("execution",): {
        "module_setup", "launcher", "launcher_args", "extra_sbatch", "walltime_guard",
    },
    ("execution", "extra_sbatch"): None,
    ("execution", "walltime_guard"): {
        "enabled", "warmup_steps", "multiplier", "min_seconds", "estimator_alpha",
    },
}


_STUDY_SCHEMA = {
    (): {
        "base_configs", "study_type", "parameters", "parameter_sets", "metrics", "plotting", "execution",
    },
    ("base_configs",): {"case", "solver", "monitor", "post"},
    ("parameters",): None,
    ("parameter_sets", "[]"): None,
    ("metrics", "[]"): {
        "name", "source", "file_glob", "column", "reduction", "normalize_by_parameter",
        "numerator_column", "denominator_column", "denominator_floor",
    },
    ("plotting",): {"enabled", "output_format"},
    ("execution",): {"max_concurrent_array_tasks"},
}


def validate_solver_configs(case_cfg: dict, solver_cfg: dict, monitor_cfg: dict,
                            case_path: str, solver_path: str, monitor_path: str):
    """!
    @brief Validates all solver input configs before any work is done.
    @details Checks for required sections, required keys, and physical sanity.
             Exits with a clear error message on the first problem found.
    @param[in] case_cfg    Parsed case YAML dictionary.
    @param[in] solver_cfg  Parsed solver YAML dictionary.
    @param[in] monitor_cfg Parsed monitor YAML dictionary.
    @param[in] case_path   Path to case file (for error messages).
    @param[in] solver_path Path to solver file (for error messages).
    @param[in] monitor_path Path to monitor file (for error messages).
    @throws SystemExit on validation failure.
    """
    errors = []
    warnings = []
    eulerian_source_mode = "solve"

    _validate_yaml_schema_keys(case_cfg, _CASE_SCHEMA, case_path, errors)
    _validate_yaml_schema_keys(solver_cfg, _SOLVER_SCHEMA, solver_path, errors)
    _validate_yaml_schema_keys(monitor_cfg, _MONITOR_SCHEMA, monitor_path, errors)

    # --- case.yml: required top-level sections ---
    required_case_sections = ['properties', 'run_control', 'grid', 'models', 'boundary_conditions']
    for section in required_case_sections:
        if section not in case_cfg:
            errors.append(f"  {case_path}: missing required section '{section}'.")

    if errors:
        _print_validation_errors(errors)

    # --- case.yml: properties sub-keys ---
    props = case_cfg.get('properties', {})
    for group, keys in [('scaling', ['length_ref', 'velocity_ref']),
                        ('fluid', ['density', 'viscosity'])]:
        sub = props.get(group, {})
        if not sub:
            errors.append(f"  {case_path}: missing 'properties.{group}' section.")
        else:
            for k in keys:
                if k not in sub:
                    errors.append(f"  {case_path}: missing key 'properties.{group}.{k}'.")

    # --- case.yml: run_control sub-keys ---
    rc = case_cfg.get('run_control', {})
    for k in ['start_step', 'total_steps', 'dt_physical']:
        if k not in rc:
            errors.append(f"  {case_path}: missing key 'run_control.{k}'.")

    # --- Physical sanity checks ---
    try:
        density = float(props.get('fluid', {}).get('density', 0))
        viscosity = float(props.get('fluid', {}).get('viscosity', 0))
        dt = float(rc.get('dt_physical', 0))
        if density <= 0:
            errors.append(f"  {case_path}: 'properties.fluid.density' must be positive (got {density}).")
        if viscosity < 0:
            errors.append(f"  {case_path}: 'properties.fluid.viscosity' must be non-negative (got {viscosity}).")
        if dt <= 0:
            errors.append(f"  {case_path}: 'run_control.dt_physical' must be positive (got {dt}).")
    except (TypeError, ValueError):
        pass  # Will be caught later during processing

    # --- case.yml: grid mode ---
    grid_cfg = case_cfg.get('grid', {})
    grid_mode = grid_cfg.get('mode')
    valid_grid_modes = ['file', 'programmatic_c', 'grid_gen']
    if grid_mode not in valid_grid_modes:
        errors.append(f"  {case_path}: 'grid.mode' must be one of {valid_grid_modes} (got '{grid_mode}').")
    elif grid_mode == 'file':
        source_file = grid_cfg.get('source_file')
        if not source_file:
            errors.append(f"  {case_path}: 'grid.source_file' is required when grid.mode is 'file'.")
        else:
            source_abs = source_file if os.path.isabs(source_file) else os.path.abspath(os.path.join(os.path.dirname(case_path), source_file))
            if not os.path.isfile(source_abs):
                errors.append(f"  {case_path}: grid.source_file does not exist: {source_abs}")

        legacy_cfg = grid_cfg.get("legacy_conversion")
        if legacy_cfg is not None:
            if not isinstance(legacy_cfg, dict):
                errors.append(f"  {case_path}: grid.legacy_conversion must be a mapping when provided.")
            else:
                enabled = legacy_cfg.get("enabled", True)
                if not isinstance(enabled, bool):
                    errors.append(f"  {case_path}: grid.legacy_conversion.enabled must be a boolean.")

                fmt = legacy_cfg.get("format")
                if fmt is not None:
                    normalized_fmt = str(fmt).strip().lower()
                    allowed_formats = {"legacy1d", "legacy_1d", "les_flat_1d", "les-flat-1d"}
                    if normalized_fmt not in allowed_formats:
                        errors.append(
                            f"  {case_path}: grid.legacy_conversion.format must be one of "
                            f"{sorted(allowed_formats)} (got '{fmt}')."
                        )

                script_path = legacy_cfg.get("script")
                if script_path is not None:
                    if not isinstance(script_path, str) or not script_path.strip():
                        errors.append(f"  {case_path}: grid.legacy_conversion.script must be a non-empty string.")
                    else:
                        script_abs = script_path if os.path.isabs(script_path) else os.path.abspath(os.path.join(os.path.dirname(case_path), script_path))
                        if not os.path.isfile(script_abs):
                            errors.append(f"  {case_path}: grid.legacy_conversion.script does not exist: {script_abs}")

                output_file = legacy_cfg.get("output_file")
                if output_file is not None and (not isinstance(output_file, str) or not output_file.strip()):
                    errors.append(f"  {case_path}: grid.legacy_conversion.output_file must be a non-empty string when provided.")

                axis_columns = legacy_cfg.get("axis_columns")
                if axis_columns is not None:
                    if not isinstance(axis_columns, list) or len(axis_columns) != 3:
                        errors.append(f"  {case_path}: grid.legacy_conversion.axis_columns must be a 3-item integer list.")
                    else:
                        for idx, value in enumerate(axis_columns):
                            if not isinstance(value, int) or value < 0:
                                errors.append(
                                    f"  {case_path}: grid.legacy_conversion.axis_columns[{idx}] must be a non-negative integer (got {value})."
                                )

                strict_trailing = legacy_cfg.get("strict_trailing")
                if strict_trailing is not None and not isinstance(strict_trailing, bool):
                    errors.append(f"  {case_path}: grid.legacy_conversion.strict_trailing must be a boolean when provided.")

                cli_args = legacy_cfg.get("cli_args")
                if cli_args is not None and not isinstance(cli_args, list):
                    errors.append(f"  {case_path}: grid.legacy_conversion.cli_args must be a list of CLI tokens.")
    elif grid_mode == 'programmatic_c':
        grid_settings = grid_cfg.get('programmatic_settings')
        if not grid_settings:
            errors.append(f"  {case_path}: 'grid.programmatic_settings' is required when grid.mode is 'programmatic_c'.")
        elif not isinstance(grid_settings, dict):
            errors.append(f"  {case_path}: 'grid.programmatic_settings' must be a mapping.")
    elif grid_mode == 'grid_gen':
        gen_cfg = grid_cfg.get('generator')
        if not isinstance(gen_cfg, dict):
            errors.append(f"  {case_path}: 'grid.generator' must be a mapping when grid.mode is 'grid_gen'.")
        else:
            warn_on_grid_generator_hyphen_keys(gen_cfg, case_path, warnings)

            config_file = gen_cfg.get('config_file')
            if not config_file:
                errors.append(f"  {case_path}: 'grid.generator.config_file' is required for grid.mode='grid_gen'.")
            else:
                config_abs = config_file if os.path.isabs(config_file) else os.path.abspath(os.path.join(os.path.dirname(case_path), config_file))
                if not os.path.isfile(config_abs):
                    errors.append(f"  {case_path}: grid.generator.config_file does not exist: {config_abs}")

            grid_type = gen_cfg.get('grid_type')
            if grid_type is not None and str(grid_type) not in {'cpipe', 'pipe', 'warp'}:
                errors.append(f"  {case_path}: grid.generator.grid_type must be one of ['cpipe','pipe','warp'] (got '{grid_type}').")

            cli_args = gen_cfg.get('cli_args', [])
            if cli_args is not None and not isinstance(cli_args, list):
                errors.append(f"  {case_path}: grid.generator.cli_args must be a list of CLI tokens.")
    try:
        resolve_grid_da_processor_layout(grid_cfg)
    except ValueError as e:
        errors.append(f"  {case_path}: {e}")

    # --- case.yml: boundary_conditions strict validation ---
    prepared_blocks = None
    try:
        prepared_blocks = validate_and_prepare_boundary_conditions(case_cfg)
    except ValueError as e:
        errors.append(f"  {case_path}: {e}")

    # --- case.yml: initial_conditions mode-aware validation ---
    ic = props.get('initial_conditions', {})
    if not ic:
        errors.append(f"  {case_path}: missing 'properties.initial_conditions' section.")
    elif not isinstance(ic, dict):
        errors.append(f"  {case_path}: 'properties.initial_conditions' must be a mapping.")
    elif 'mode' not in ic:
        errors.append(
            f"  {case_path}: missing key 'properties.initial_conditions.mode'. "
            "Specify 'generated' or 'file' explicitly."
        )
    else:
        try:
            resolve_initial_condition_config(ic, prepared_blocks, U_ref=1.0)
        except KeyError as e:
            errors.append(f"  {case_path}: missing key 'properties.initial_conditions.{e.args[0]}'.")
        except ValueError as e:
            errors.append(f"  {case_path}: {e}")

    # --- case.yml: particle initialization validation ---
    particles_cfg = case_cfg.get('models', {}).get('physics', {}).get('particles', {})
    if particles_cfg and not isinstance(particles_cfg, dict):
        errors.append(f"  {case_path}: 'models.physics.particles' must be a mapping.")
    elif isinstance(particles_cfg, dict):
        init_mode_raw = particles_cfg.get('init_mode', 'Surface')
        try:
            pinit_code = normalize_particle_init_mode(init_mode_raw)
        except ValueError as e:
            errors.append(f"  {case_path}: {e}")
            pinit_code = None

        restart_mode = particles_cfg.get('restart_mode')
        if restart_mode is not None and str(restart_mode).lower() not in {"init", "load"}:
            errors.append(
                f"  {case_path}: models.physics.particles.restart_mode must be 'init' or 'load' (got '{restart_mode}')."
            )
        elif 'restart_mode' not in particles_cfg:
            try:
                start_step = int(rc.get('start_step', 0))
                particle_count = int(particles_cfg.get('count', 0) or 0)
            except (TypeError, ValueError):
                start_step = 0
                particle_count = 0
            if start_step > 0 and particle_count > 0:
                warnings.append(
                    f"{case_path}: models.physics.particles.restart_mode is omitted for a particle restart "
                    "(run_control.start_step > 0, count > 0). C will default to 'load'."
                )

        if pinit_code == 2:
            point_cfg = particles_cfg.get('point_source', {})
            if not isinstance(point_cfg, dict):
                errors.append(f"  {case_path}: models.physics.particles.point_source must be a mapping when init_mode is PointSource.")
            else:
                for coord in ('x', 'y', 'z'):
                    if coord not in point_cfg:
                        errors.append(
                            f"  {case_path}: models.physics.particles.point_source.{coord} is required when init_mode is PointSource."
                        )

    # --- case.yml: turbulence model validation ---
    turbulence_cfg = case_cfg.get('models', {}).get('physics', {}).get('turbulence', {})
    if turbulence_cfg is not None and not isinstance(turbulence_cfg, dict):
        errors.append(f"  {case_path}: 'models.physics.turbulence' must be a mapping.")
    elif isinstance(turbulence_cfg, dict) and turbulence_cfg:
        try:
            append_turbulence_flags(case_cfg.get('models', {}), [])
        except ValueError as e:
            errors.append(f"  {case_path}: {e}")

        les_cfg = turbulence_cfg.get('les')
        rans_cfg = turbulence_cfg.get('rans')
        wall_cfg = turbulence_cfg.get('wall_function')

        if isinstance(les_cfg, dict):
            for key in ('enabled',):
                if key in les_cfg and not isinstance(les_cfg[key], bool):
                    errors.append(f"  {case_path}: models.physics.turbulence.les.{key} must be true or false.")
            for key in ('constant_cs', 'max_cs'):
                if key in les_cfg:
                    try:
                        value = float(les_cfg[key])
                        if value < 0.0:
                            errors.append(f"  {case_path}: models.physics.turbulence.les.{key} must be nonnegative.")
                    except (TypeError, ValueError):
                        errors.append(f"  {case_path}: models.physics.turbulence.les.{key} must be numeric.")
            if 'dynamic_frequency' in les_cfg:
                try:
                    value = int(les_cfg['dynamic_frequency'])
                    if value <= 0:
                        errors.append(f"  {case_path}: models.physics.turbulence.les.dynamic_frequency must be positive.")
                except (TypeError, ValueError):
                    errors.append(f"  {case_path}: models.physics.turbulence.les.dynamic_frequency must be an integer.")

        if isinstance(rans_cfg, dict):
            if 'enabled' in rans_cfg and not isinstance(rans_cfg['enabled'], bool):
                errors.append(f"  {case_path}: models.physics.turbulence.rans.enabled must be true or false.")
            try:
                rans_enabled = bool(rans_cfg.get('enabled', True)) and normalize_rans_model(rans_cfg.get('model', 'k_omega')) != 0
            except ValueError:
                rans_enabled = False
            if rans_enabled:
                warnings.append(
                    f"{case_path}: models.physics.turbulence.rans is accepted, but the k-omega runtime update is currently incomplete."
                )
        elif rans_cfg:
            warnings.append(
                f"{case_path}: models.physics.turbulence.rans is accepted, but the k-omega runtime update is currently incomplete."
            )

        if isinstance(wall_cfg, dict):
            if 'enabled' in wall_cfg and not isinstance(wall_cfg['enabled'], bool):
                errors.append(f"  {case_path}: models.physics.turbulence.wall_function.enabled must be true or false.")
            if 'roughness_height' in wall_cfg:
                try:
                    value = float(wall_cfg['roughness_height'])
                    if value < 0.0:
                        errors.append(f"  {case_path}: models.physics.turbulence.wall_function.roughness_height must be nonnegative.")
                except (TypeError, ValueError):
                    errors.append(f"  {case_path}: models.physics.turbulence.wall_function.roughness_height must be numeric.")

    # --- solver.yml: basic structure ---
    if not isinstance(solver_cfg, dict) or not solver_cfg:
        errors.append(f"  {solver_path}: solver config is empty or not a valid YAML mapping.")
    else:
        strategy_cfg = solver_cfg.get('strategy', {})
        if not isinstance(strategy_cfg, dict):
            errors.append(f"  {solver_path}: 'strategy' must be a mapping.")
        elif 'implicit' in strategy_cfg:
            errors.append(
                f"  {solver_path}: unsupported old key 'strategy.implicit' is not supported. "
                "Use 'strategy.momentum_solver' with named solver values."
            )
        if isinstance(strategy_cfg, dict) and 'momentum_solver' in strategy_cfg:
            try:
                normalize_momentum_solver_type(strategy_cfg['momentum_solver'])
            except ValueError as e:
                errors.append(f"  {solver_path}: {e}")

        op_mode_cfg = solver_cfg.get('operation_mode', {})
        if op_mode_cfg is not None and not isinstance(op_mode_cfg, dict):
            errors.append(f"  {solver_path}: 'operation_mode' must be a mapping when provided.")
        elif isinstance(op_mode_cfg, dict):
            eulerian_source_mode = None
            normalized_analytical_type = None
            if 'eulerian_field_source' in op_mode_cfg:
                try:
                    eulerian_source_mode = normalize_eulerian_field_source(op_mode_cfg.get('eulerian_field_source'))
                except ValueError as e:
                    errors.append(f"  {solver_path}: {e}")

            analytical_type = op_mode_cfg.get('analytical_type')
            if analytical_type is not None:
                try:
                    normalized_analytical_type = normalize_analytical_type(analytical_type)
                except ValueError as e:
                    errors.append(f"  {solver_path}: {e}")
                else:
                    uniform_flow_cfg = op_mode_cfg.get('uniform_flow')
                    if uniform_flow_cfg is not None and not isinstance(uniform_flow_cfg, dict):
                        errors.append(f"  {solver_path}: 'operation_mode.uniform_flow' must be a mapping when provided.")
                    elif normalized_analytical_type == "UNIFORM_FLOW":
                        if not isinstance(uniform_flow_cfg, dict):
                            errors.append(
                                f"  {solver_path}: operation_mode.uniform_flow is required when "
                                "operation_mode.analytical_type is 'UNIFORM_FLOW'."
                            )
                        else:
                            for coord in ("u", "v", "w"):
                                if coord not in uniform_flow_cfg:
                                    errors.append(
                                        f"  {solver_path}: operation_mode.uniform_flow.{coord} is required for UNIFORM_FLOW."
                                    )
                                else:
                                    try:
                                        float(uniform_flow_cfg[coord])
                                    except (TypeError, ValueError):
                                        errors.append(
                                            f"  {solver_path}: operation_mode.uniform_flow.{coord} must be numeric."
                                        )
                    elif uniform_flow_cfg is not None:
                        errors.append(
                            f"  {solver_path}: operation_mode.uniform_flow is only valid when "
                            "operation_mode.analytical_type is 'UNIFORM_FLOW'."
                        )

            if eulerian_source_mode == "analytical":
                effective_analytical_type = normalized_analytical_type or "TGV3D"
                if effective_analytical_type == "TGV3D":
                    if grid_mode != 'programmatic_c':
                        errors.append(
                            f"  {case_path}: analytical type '{effective_analytical_type}' requires grid.mode "
                            "'programmatic_c'. File-backed analytical ingestion is only supported for "
                            "ZERO_FLOW and UNIFORM_FLOW."
                        )
                    elif isinstance(grid_cfg.get('programmatic_settings'), dict):
                        missing_dims = [key for key in ('im', 'jm', 'km') if key not in grid_cfg['programmatic_settings']]
                        if missing_dims:
                            errors.append(
                                f"  {case_path}: grid.programmatic_settings must include {missing_dims} when "
                                f"operation_mode.analytical_type resolves to '{effective_analytical_type}'."
                            )
                else:
                    if grid_mode not in {'programmatic_c', 'file'}:
                        errors.append(
                            f"  {case_path}: grid.mode '{grid_mode}' is not supported when "
                            f"operation_mode.analytical_type is '{effective_analytical_type}'. "
                            "Use 'programmatic_c' or 'file'."
                        )
                    elif grid_mode == 'programmatic_c' and isinstance(grid_cfg.get('programmatic_settings'), dict):
                        missing_dims = [key for key in ('im', 'jm', 'km') if key not in grid_cfg['programmatic_settings']]
                        if missing_dims:
                            errors.append(
                                f"  {case_path}: grid.programmatic_settings must include {missing_dims} when "
                                f"operation_mode.analytical_type is '{effective_analytical_type}' and "
                                "grid.mode is 'programmatic_c'."
                            )

        verification_cfg = solver_cfg.get('verification', {})
        if verification_cfg is not None and not isinstance(verification_cfg, dict):
            errors.append(f"  {solver_path}: 'verification' must be a mapping when provided.")
        elif isinstance(verification_cfg, dict) and verification_cfg:
            sources_cfg = verification_cfg.get('sources', {})
            if sources_cfg is not None and not isinstance(sources_cfg, dict):
                errors.append(f"  {solver_path}: 'verification.sources' must be a mapping when provided.")
            elif isinstance(sources_cfg, dict) and sources_cfg:
                diff_cfg = sources_cfg.get('diffusivity')
                scalar_cfg = sources_cfg.get('scalar')

                if diff_cfg is not None:
                    if not isinstance(diff_cfg, dict):
                        errors.append(f"  {solver_path}: 'verification.sources.diffusivity' must be a mapping.")
                    else:
                        if eulerian_source_mode != "analytical":
                            errors.append(
                                f"  {solver_path}: verification.sources.diffusivity is only valid when "
                                "operation_mode.eulerian_field_source is 'analytical'."
                            )
                        mode = diff_cfg.get('mode')
                        profile = diff_cfg.get('profile')
                        if str(mode).strip().lower() != "analytical":
                            errors.append(
                                f"  {solver_path}: verification.sources.diffusivity.mode must be 'analytical'."
                            )
                        if str(profile).strip().upper() != "LINEAR_X":
                            errors.append(
                                f"  {solver_path}: verification.sources.diffusivity.profile must be 'LINEAR_X'."
                            )
                        for key in ("gamma0", "slope_x"):
                            if key not in diff_cfg:
                                errors.append(
                                    f"  {solver_path}: verification.sources.diffusivity.{key} is required."
                                )
                            else:
                                try:
                                    float(diff_cfg[key])
                                except (TypeError, ValueError):
                                    errors.append(
                                        f"  {solver_path}: verification.sources.diffusivity.{key} must be numeric."
                                    )

                if scalar_cfg is not None:
                    if not isinstance(scalar_cfg, dict):
                        errors.append(f"  {solver_path}: 'verification.sources.scalar' must be a mapping.")
                    else:
                        if eulerian_source_mode != "analytical":
                            errors.append(
                                f"  {solver_path}: verification.sources.scalar is only valid when "
                                "operation_mode.eulerian_field_source is 'analytical'."
                            )
                        mode = scalar_cfg.get('mode')
                        profile = str(scalar_cfg.get('profile', '')).strip().upper()
                        if str(mode).strip().lower() != "analytical":
                            errors.append(
                                f"  {solver_path}: verification.sources.scalar.mode must be 'analytical'."
                            )
                        if profile not in {"CONSTANT", "LINEAR_X", "SIN_PRODUCT"}:
                            errors.append(
                                f"  {solver_path}: verification.sources.scalar.profile must be one of CONSTANT, LINEAR_X, SIN_PRODUCT."
                            )
                        required_scalar_keys = {
                            "CONSTANT": ("value",),
                            "LINEAR_X": ("phi0", "slope_x"),
                            "SIN_PRODUCT": ("amplitude", "kx", "ky", "kz"),
                        }.get(profile, ())
                        for key in required_scalar_keys:
                            if key not in scalar_cfg:
                                errors.append(
                                    f"  {solver_path}: verification.sources.scalar.{key} is required for profile '{profile}'."
                                )
                            else:
                                try:
                                    float(scalar_cfg[key])
                                except (TypeError, ValueError):
                                    errors.append(
                                        f"  {solver_path}: verification.sources.scalar.{key} must be numeric."
                                    )

                unknown_source_keys = sorted(set(sources_cfg.keys()) - {"diffusivity", "scalar"})
                if unknown_source_keys:
                    errors.append(
                        f"  {solver_path}: unsupported verification.sources entries: {unknown_source_keys}. "
                        "Currently supported: 'diffusivity', 'scalar'."
                    )
            unknown_verification_keys = sorted(set(verification_cfg.keys()) - {"sources"})
            if unknown_verification_keys:
                errors.append(
                    f"  {solver_path}: unsupported verification keys: {unknown_verification_keys}. "
                    "Currently supported: 'sources'."
                )

        transport_cfg = solver_cfg.get('scalar_transport', {})
        if transport_cfg is not None and not isinstance(transport_cfg, dict):
            errors.append(f"  {solver_path}: 'scalar_transport' must be a mapping when provided.")
        elif isinstance(transport_cfg, dict):
            unknown_transport_keys = sorted(set(transport_cfg.keys()) - {"schmidt_number", "turbulent_schmidt_number"})
            if unknown_transport_keys:
                errors.append(
                    f"  {solver_path}: unsupported scalar_transport entries: {unknown_transport_keys}. "
                    "Currently supported: 'schmidt_number', 'turbulent_schmidt_number'."
                )
            for key in ("schmidt_number", "turbulent_schmidt_number"):
                if key in transport_cfg:
                    try:
                        value = float(transport_cfg[key])
                        if value <= 0.0:
                            errors.append(f"  {solver_path}: scalar_transport.{key} must be positive.")
                    except (TypeError, ValueError):
                        errors.append(f"  {solver_path}: scalar_transport.{key} must be numeric.")

        tolerances_cfg = solver_cfg.get('tolerances', {})
        if tolerances_cfg is not None and not isinstance(tolerances_cfg, dict):
            errors.append(f"  {solver_path}: 'tolerances' must be a mapping when provided.")
        elif isinstance(tolerances_cfg, dict):
            for key in ("absolute_tol", "relative_tol", "residual_absolute_tol", "residual_relative_tol"):
                if key in tolerances_cfg:
                    try:
                        float(tolerances_cfg[key])
                    except (TypeError, ValueError):
                        errors.append(f"  {solver_path}: tolerances.{key} must be numeric.")

        ms_cfg = solver_cfg.get('momentum_solver', {})
        if ms_cfg is not None and not isinstance(ms_cfg, dict):
            errors.append(f"  {solver_path}: 'momentum_solver' must be a mapping when provided.")
        elif isinstance(ms_cfg, dict):
            unsupported_flat_keys = {
                'max_pseudo_steps', 'absolute_tol', 'relative_tol', 'step_tol',
                'pseudo_cfl', 'jameson_residual_noise_allowance_factor',
                'rk4_residual_noise_allowance_factor'
            }
            present_unsupported = sorted(unsupported_flat_keys.intersection(ms_cfg.keys()))
            if present_unsupported:
                errors.append(
                    f"  {solver_path}: unsupported flat keys in 'momentum_solver' are not supported: {present_unsupported}. "
                    "Use solver-specific sub-blocks (e.g., momentum_solver.dual_time_picard_jameson_rk)."
                )

            allowed_ms_keys = {'dual_time_picard_jameson_rk', 'dual_time_picard_rk4'}
            unknown_ms_keys = sorted(set(ms_cfg.keys()) - allowed_ms_keys)
            if unknown_ms_keys:
                errors.append(
                    f"  {solver_path}: unsupported momentum_solver blocks/keys: {unknown_ms_keys}. "
                    "Currently supported: 'dual_time_picard_jameson_rk'."
                )
            if 'dual_time_picard_jameson_rk' in ms_cfg and 'dual_time_picard_rk4' in ms_cfg:
                errors.append(
                    f"  {solver_path}: use only momentum_solver.dual_time_picard_jameson_rk; "
                    "do not also set its deprecated dual_time_picard_rk4 alias."
                )

            selected_solver = None
            if isinstance(strategy_cfg, dict) and 'momentum_solver' in strategy_cfg:
                try:
                    selected_solver = normalize_momentum_solver_type(strategy_cfg['momentum_solver'])
                except ValueError:
                    pass
            if selected_solver is None:
                selected_solver = "DUALTIME_PICARD_JAMESON_RK"

            has_dualtime_block = (
                'dual_time_picard_jameson_rk' in ms_cfg or 'dual_time_picard_rk4' in ms_cfg
            )
            if selected_solver != "DUALTIME_PICARD_JAMESON_RK" and has_dualtime_block:
                errors.append(
                    f"  {solver_path}: momentum_solver.dual_time_picard_jameson_rk is set but selected solver is "
                    f"{selected_solver}. Solver-specific blocks must match the selected solver."
                )

            dt_picard_cfg = ms_cfg.get('dual_time_picard_jameson_rk', ms_cfg.get('dual_time_picard_rk4'))
            if dt_picard_cfg is not None:
                if not isinstance(dt_picard_cfg, dict):
                    errors.append(f"  {solver_path}: momentum_solver.dual_time_picard_jameson_rk must be a mapping.")
                else:
                    allowed_dt_keys = {
                        'max_pseudo_steps', 'absolute_tol', 'relative_tol', 'step_tol',
                        'pseudo_cfl', 'jameson_residual_noise_allowance_factor',
                        'rk4_residual_noise_allowance_factor'
                    }
                    unknown_dt_keys = sorted(set(dt_picard_cfg.keys()) - allowed_dt_keys)
                    if unknown_dt_keys:
                        errors.append(
                            f"  {solver_path}: unsupported keys in momentum_solver.dual_time_picard_jameson_rk: {unknown_dt_keys}."
                        )
                    if ('jameson_residual_noise_allowance_factor' in dt_picard_cfg and
                            'rk4_residual_noise_allowance_factor' in dt_picard_cfg):
                        errors.append(
                            f"  {solver_path}: use only jameson_residual_noise_allowance_factor; "
                            "do not also set its deprecated rk4_residual_noise_allowance_factor alias."
                        )
                    if 'pseudo_cfl' in dt_picard_cfg:
                        pcfl_cfg = dt_picard_cfg['pseudo_cfl']
                        if not isinstance(pcfl_cfg, dict):
                            errors.append(f"  {solver_path}: momentum_solver.dual_time_picard_jameson_rk.pseudo_cfl must be a mapping.")
                        else:
                            allowed_pcfl_keys = {'initial', 'minimum', 'maximum', 'growth_factor', 'reduction_factor'}
                            unknown_pcfl_keys = sorted(set(pcfl_cfg.keys()) - allowed_pcfl_keys)
                            if unknown_pcfl_keys:
                                errors.append(
                                    f"  {solver_path}: unsupported keys in momentum_solver.dual_time_picard_jameson_rk.pseudo_cfl: {unknown_pcfl_keys}."
                                )
                            numeric_pcfl = {}
                            for key in allowed_pcfl_keys:
                                if key in pcfl_cfg:
                                    try:
                                        numeric_pcfl[key] = float(pcfl_cfg[key])
                                    except (TypeError, ValueError):
                                        errors.append(
                                            f"  {solver_path}: momentum_solver.dual_time_picard_jameson_rk.pseudo_cfl.{key} must be numeric."
                                        )
                            if numeric_pcfl.get('minimum', 1.0) <= 0.0:
                                errors.append(f"  {solver_path}: pseudo_cfl.minimum must be positive.")
                            if numeric_pcfl.get('growth_factor', 1.0) < 1.0:
                                errors.append(f"  {solver_path}: pseudo_cfl.growth_factor must be at least 1.")
                            reduction = numeric_pcfl.get('reduction_factor', 1.0)
                            if reduction <= 0.0 or reduction >= 1.0:
                                errors.append(f"  {solver_path}: pseudo_cfl.reduction_factor must be in (0, 1).")
                            if all(key in numeric_pcfl for key in ('minimum', 'initial', 'maximum')):
                                if not numeric_pcfl['minimum'] <= numeric_pcfl['initial'] <= numeric_pcfl['maximum']:
                                    errors.append(f"  {solver_path}: pseudo_cfl requires minimum <= initial <= maximum.")
                    noise_key = (
                        'jameson_residual_noise_allowance_factor'
                        if 'jameson_residual_noise_allowance_factor' in dt_picard_cfg
                        else 'rk4_residual_noise_allowance_factor'
                    )
                    if noise_key in dt_picard_cfg:
                        try:
                            if float(dt_picard_cfg[noise_key]) < 1.0:
                                errors.append(f"  {solver_path}: {noise_key} must be at least 1.")
                        except (TypeError, ValueError):
                            errors.append(f"  {solver_path}: {noise_key} must be numeric.")

        solution_convergence_cfg = solver_cfg.get('solution_convergence', {})
        if solution_convergence_cfg is not None and not isinstance(solution_convergence_cfg, dict):
            errors.append(f"  {solver_path}: 'solution_convergence' must be a mapping when provided.")
        elif isinstance(solution_convergence_cfg, dict) and solution_convergence_cfg:
            allowed_solution_convergence_keys = {
                'enabled', 'mode', 'periodic_deterministic', 'statistical_steady'
            }
            unknown_solution_convergence_keys = sorted(set(solution_convergence_cfg.keys()) - allowed_solution_convergence_keys)
            if unknown_solution_convergence_keys:
                errors.append(
                    f"  {solver_path}: unsupported solution_convergence keys: {unknown_solution_convergence_keys}."
                )

            mode = solution_convergence_cfg.get('mode', 'steady_deterministic')
            try:
                normalized_solution_mode = normalize_solution_convergence_mode(mode)
            except ValueError as e:
                errors.append(f"  {solver_path}: {e}")
                normalized_solution_mode = None

            periodic_cfg = solution_convergence_cfg.get('periodic_deterministic')
            statistical_cfg = solution_convergence_cfg.get('statistical_steady')
            if periodic_cfg is not None and not isinstance(periodic_cfg, dict):
                errors.append(f"  {solver_path}: solution_convergence.periodic_deterministic must be a mapping when provided.")
            if statistical_cfg is not None and not isinstance(statistical_cfg, dict):
                errors.append(f"  {solver_path}: solution_convergence.statistical_steady must be a mapping when provided.")

            if isinstance(periodic_cfg, dict):
                unknown_periodic_keys = sorted(set(periodic_cfg.keys()) - {'period_steps'})
                if unknown_periodic_keys:
                    errors.append(
                        f"  {solver_path}: unsupported keys in solution_convergence.periodic_deterministic: {unknown_periodic_keys}."
                    )
                period_steps = periodic_cfg.get('period_steps')
                if period_steps is not None and (not isinstance(period_steps, int) or period_steps <= 0):
                    errors.append(f"  {solver_path}: solution_convergence.periodic_deterministic.period_steps must be a positive integer.")

            if isinstance(statistical_cfg, dict):
                unknown_statistical_keys = sorted(set(statistical_cfg.keys()) - {'window_steps'})
                if unknown_statistical_keys:
                    errors.append(
                        f"  {solver_path}: unsupported keys in solution_convergence.statistical_steady: {unknown_statistical_keys}."
                    )
                window_steps = statistical_cfg.get('window_steps')
                if window_steps is not None and (not isinstance(window_steps, int) or window_steps <= 0):
                    errors.append(f"  {solver_path}: solution_convergence.statistical_steady.window_steps must be a positive integer.")

            if normalized_solution_mode == "PERIODIC_DETERMINISTIC":
                if not isinstance(periodic_cfg, dict) or 'period_steps' not in periodic_cfg:
                    errors.append(
                        f"  {solver_path}: solution_convergence.periodic_deterministic.period_steps is required when mode is 'periodic_deterministic'."
                    )
            elif periodic_cfg is not None:
                errors.append(
                    f"  {solver_path}: solution_convergence.periodic_deterministic is only valid when mode is 'periodic_deterministic'."
                )

            if normalized_solution_mode == "STATISTICAL_STEADY":
                if not isinstance(statistical_cfg, dict) or 'window_steps' not in statistical_cfg:
                    errors.append(
                        f"  {solver_path}: solution_convergence.statistical_steady.window_steps is required when mode is 'statistical_steady'."
                    )
            elif statistical_cfg is not None:
                errors.append(
                    f"  {solver_path}: solution_convergence.statistical_steady is only valid when mode is 'statistical_steady'."
                )

    # --- solver.yml: interpolation section ---
    interp_cfg = solver_cfg.get('interpolation', {}) if isinstance(solver_cfg, dict) else {}
    if interp_cfg is not None and not isinstance(interp_cfg, dict):
        errors.append(f"  {solver_path}: 'interpolation' must be a mapping when provided.")
    elif isinstance(interp_cfg, dict) and 'method' in interp_cfg:
        try:
            normalize_interpolation_method(interp_cfg['method'])
        except ValueError as e:
            errors.append(f"  {solver_path}: {e}")

    # --- monitor.yml: basic structure ---
    if not isinstance(monitor_cfg, dict) or not monitor_cfg:
        errors.append(f"  {monitor_path}: monitor config is empty or not a valid YAML mapping.")
    else:
        io_cfg = monitor_cfg.get('io', {})
        freq = io_cfg.get('data_output_frequency')
        if freq is not None and (not isinstance(freq, int) or freq <= 0):
            errors.append(f"  {monitor_path}: 'io.data_output_frequency' must be a positive integer (got {freq}).")
        particle_console_freq = io_cfg.get('particle_console_output_frequency')
        if particle_console_freq is not None and (not isinstance(particle_console_freq, int) or particle_console_freq < 0):
            errors.append(
                f"  {monitor_path}: 'io.particle_console_output_frequency' must be a non-negative integer "
                f"(got {particle_console_freq})."
            )
        try:
            resolve_profiling_config(monitor_cfg)
        except ValueError as e:
            errors.append(f"  {monitor_path}: {e}")
        try:
            resolve_diagnostics_config(monitor_cfg)
        except ValueError as e:
            errors.append(f"  {monitor_path}: {e}")
        try:
            resolve_solver_monitoring_flags(monitor_cfg)
        except ValueError as e:
            errors.append(f"  {monitor_path}: {e}")

    if not errors:
        if needs_restart_source(case_cfg, solver_cfg):
            warnings.append(
                f"{case_path}: This configuration requires restart data (start_step > 0, "
                "eulerian_field_source='load', or particle restart_mode='load'). "
                "Use --restart-from or --continue when running."
            )

    if errors:
        _print_validation_errors(errors)
    for warning in warnings:
        print(f"[WARN] {warning}", file=sys.stderr)


def validate_post_config(post_cfg: dict, post_path: str):
    """!
    @brief Validates the post-processing config before running the post-processor.
    @param[in] post_cfg  Parsed post-processing YAML dictionary.
    @param[in] post_path Path to post file (for error messages).
    @throws SystemExit on validation failure.
    """
    errors = []

    _validate_yaml_schema_keys(post_cfg, _POST_SCHEMA, post_path, errors)

    if not isinstance(post_cfg, dict) or not post_cfg:
        errors.append(f"  {post_path}: post-processing config is empty or not a valid YAML mapping.")
        _print_validation_errors(errors)

    # --- run_control ---
    if 'run_control' not in post_cfg:
        errors.append(f"  {post_path}: missing required section 'run_control'.")
    else:
        rc = post_cfg.get('run_control', {})
        if not isinstance(rc, dict):
            errors.append(f"  {post_path}: 'run_control' must be a mapping.")
        else:
            for canonical_key, aliases in POST_RUN_CONTROL_ALIASES.items():
                if not any(alias in rc for alias in aliases):
                    alias_list = "', '".join(aliases)
                    errors.append(
                        f"  {post_path}: missing required key 'run_control.{canonical_key}' "
                        f"(accepted aliases: '{alias_list}')."
                    )
                    continue
                raw_value = _mapping_value_with_aliases(rc, *aliases)
                try:
                    int(raw_value)
                except (TypeError, ValueError):
                    alias_name = next((alias for alias in aliases if alias in rc), canonical_key)
                    errors.append(
                        f"  {post_path}: 'run_control.{alias_name}' must be an integer-compatible value."
                    )

    # --- io section ---
    io_cfg = post_cfg.get('io', {})
    source_cfg = post_cfg.get('source_data')
    if source_cfg is not None and not isinstance(source_cfg, dict):
        errors.append(f"  {post_path}: 'source_data' must be a mapping when provided.")
    global_ops = post_cfg.get('global_operations')
    if global_ops is not None:
        if not isinstance(global_ops, dict):
            errors.append(f"  {post_path}: 'global_operations' must be a mapping when provided.")
        elif 'dimensionalize' in global_ops and not isinstance(global_ops.get('dimensionalize'), bool):
            errors.append(f"  {post_path}: 'global_operations.dimensionalize' must be a boolean.")
    if not io_cfg:
        errors.append(f"  {post_path}: missing required section 'io'.")
    elif not isinstance(io_cfg, dict):
        errors.append(f"  {post_path}: 'io' must be a mapping.")
    else:
        for k in ['output_directory', 'output_filename_prefix']:
            if k not in io_cfg:
                errors.append(f"  {post_path}: missing required key 'io.{k}'.")
        for key_name in ('output_directory', 'output_filename_prefix', 'particle_filename_prefix'):
            if key_name in io_cfg and not isinstance(io_cfg.get(key_name), str):
                errors.append(f"  {post_path}: 'io.{key_name}' must be a string when provided.")
        if 'output_particles' in io_cfg and not isinstance(io_cfg.get('output_particles'), bool):
            errors.append(f"  {post_path}: 'io.output_particles' must be a boolean when provided.")
        particle_subsampling_frequency = io_cfg.get('particle_subsampling_frequency')
        if particle_subsampling_frequency is not None:
            if not isinstance(particle_subsampling_frequency, int) or particle_subsampling_frequency <= 0:
                errors.append(
                    f"  {post_path}: 'io.particle_subsampling_frequency' must be a positive integer when provided."
                )
        input_extensions = io_cfg.get('input_extensions')
        source_input_extensions = get_post_source_data(post_cfg).get('input_extensions')
        if input_extensions is not None:
            if not isinstance(input_extensions, dict):
                errors.append(f"  {post_path}: 'io.input_extensions' must be a mapping when provided.")
            else:
                for ext_key in ('eulerian', 'particle'):
                    ext_val = input_extensions.get(ext_key)
                    if ext_val is not None and not isinstance(ext_val, str):
                        errors.append(f"  {post_path}: 'io.input_extensions.{ext_key}' must be a string extension.")
        if source_input_extensions is not None:
            if not isinstance(source_input_extensions, dict):
                errors.append(f"  {post_path}: 'source_data.input_extensions' must be a mapping when provided.")
            else:
                for ext_key in ('eulerian', 'particle'):
                    ext_val = source_input_extensions.get(ext_key)
                    if ext_val is not None and not isinstance(ext_val, str):
                        errors.append(
                            f"  {post_path}: 'source_data.input_extensions.{ext_key}' must be a string extension."
                        )

        averaged_fields = io_cfg.get('eulerian_fields_averaged')
        if averaged_fields is not None and not isinstance(averaged_fields, list):
            errors.append(f"  {post_path}: 'io.eulerian_fields_averaged' must be a list when provided.")
        for list_key in ('eulerian_fields', 'particle_fields'):
            list_val = io_cfg.get(list_key)
            if list_val is not None and not isinstance(list_val, list):
                errors.append(f"  {post_path}: 'io.{list_key}' must be a list when provided.")

    # --- Check eulerian_pipeline entries have 'task' key ---
    eulerian_pipeline = post_cfg.get('eulerian_pipeline', [])
    if eulerian_pipeline is not None and not isinstance(eulerian_pipeline, list):
        errors.append(f"  {post_path}: 'eulerian_pipeline' must be a list when provided.")
        eulerian_pipeline = []
    for i, entry in enumerate(eulerian_pipeline):
        if not isinstance(entry, dict) or 'task' not in entry:
            errors.append(f"  {post_path}: 'eulerian_pipeline[{i}]' is missing the 'task' key. "
                          "Check YAML indentation (each entry needs '- task: ...' with proper spacing).")
            continue
        task_name = entry.get('task')
        if task_name == 'q_criterion':
            continue
        if task_name == 'nodal_average':
            in_field = entry.get('input_field')
            out_field = entry.get('output_field')
            if not isinstance(in_field, str) or not in_field.strip():
                errors.append(f"  {post_path}: 'eulerian_pipeline[{i}].input_field' must be a non-empty string.")
            if not isinstance(out_field, str) or not out_field.strip():
                errors.append(f"  {post_path}: 'eulerian_pipeline[{i}].output_field' must be a non-empty string.")
            if isinstance(in_field, str) and isinstance(out_field, str) and in_field == out_field:
                errors.append(
                    f"  {post_path}: 'eulerian_pipeline[{i}]' nodal_average input and output fields must differ."
                )
            continue
        if task_name == 'normalize_field':
            field = entry.get('field', 'P')
            if not isinstance(field, str) or not field.strip():
                errors.append(f"  {post_path}: 'eulerian_pipeline[{i}].field' must be a non-empty string.")
            elif field != 'P':
                errors.append(
                    f"  {post_path}: 'eulerian_pipeline[{i}].field' currently only supports 'P' "
                    f"(got '{field}')."
                )
            reference_point = entry.get('reference_point', [1, 1, 1])
            if not isinstance(reference_point, (list, tuple)) or len(reference_point) != 3:
                errors.append(
                    f"  {post_path}: 'eulerian_pipeline[{i}].reference_point' must be a 3-item list."
                )
            else:
                for rp_idx, coord in enumerate(reference_point):
                    try:
                        int(coord)
                    except (TypeError, ValueError):
                        errors.append(
                            f"  {post_path}: 'eulerian_pipeline[{i}].reference_point[{rp_idx}]' "
                            "must be integer-compatible."
                        )
            continue
        errors.append(
            f"  {post_path}: unsupported eulerian task '{task_name}' at eulerian_pipeline[{i}]."
        )

    # --- Check lagrangian_pipeline entries have 'task' key ---
    lagrangian_pipeline = post_cfg.get('lagrangian_pipeline', [])
    if lagrangian_pipeline is not None and not isinstance(lagrangian_pipeline, list):
        errors.append(f"  {post_path}: 'lagrangian_pipeline' must be a list when provided.")
        lagrangian_pipeline = []
    for i, entry in enumerate(lagrangian_pipeline):
        if not isinstance(entry, dict) or 'task' not in entry:
            errors.append(f"  {post_path}: 'lagrangian_pipeline[{i}]' is missing the 'task' key.")
            continue
        task_name = entry.get('task')
        if task_name == 'specific_ke':
            in_field = entry.get('input_field')
            out_field = entry.get('output_field')
            if not isinstance(in_field, str) or not in_field.strip():
                errors.append(f"  {post_path}: 'lagrangian_pipeline[{i}].input_field' must be a non-empty string.")
            if not isinstance(out_field, str) or not out_field.strip():
                errors.append(f"  {post_path}: 'lagrangian_pipeline[{i}].output_field' must be a non-empty string.")
            continue
        errors.append(
            f"  {post_path}: unsupported lagrangian task '{task_name}' at lagrangian_pipeline[{i}]."
        )

    # --- Check statistics pipeline entries ---
    stats_cfg = post_cfg.get('statistics_pipeline')
    stats_entries = []
    if stats_cfg is not None:
        if isinstance(stats_cfg, list):
            stats_entries = stats_cfg
        elif isinstance(stats_cfg, dict):
            stats_entries = stats_cfg.get('tasks', [])
            if not isinstance(stats_entries, list):
                errors.append(f"  {post_path}: 'statistics_pipeline.tasks' must be a list.")
            stats_output_prefix = stats_cfg.get('output_prefix')
            if stats_output_prefix is not None and not isinstance(stats_output_prefix, str):
                errors.append(f"  {post_path}: 'statistics_pipeline.output_prefix' must be a string.")
        else:
            errors.append(
                f"  {post_path}: 'statistics_pipeline' must be either a list of tasks or a mapping with a 'tasks' list."
            )
        for i, entry in enumerate(stats_entries):
            if isinstance(entry, str):
                task_name = entry
            elif isinstance(entry, dict) and 'task' in entry:
                task_name = entry.get('task')
            else:
                errors.append(
                    f"  {post_path}: statistics task entry {i} must be either a string or a mapping with key 'task'."
                )
                continue
            try:
                normalize_statistics_task(task_name)
            except ValueError as e:
                errors.append(f"  {post_path}: {e}")

    legacy_stats_output_prefix = post_cfg.get('statistics_output_prefix')
    if legacy_stats_output_prefix is not None and not isinstance(legacy_stats_output_prefix, str):
        errors.append(f"  {post_path}: 'statistics_output_prefix' must be a string when provided.")

    if errors:
        _print_validation_errors(errors)

def validate_cluster_config(cluster_cfg: dict, cluster_path: str):
    """!
    @brief Validate Slurm scheduler configuration from cluster.yml.
    @param[in] cluster_cfg Argument passed to `validate_cluster_config()`.
    @param[in] cluster_path Argument passed to `validate_cluster_config()`.
    """
    errors = []
    warnings = []
    _validate_yaml_schema_keys(cluster_cfg, _CLUSTER_SCHEMA, cluster_path, errors)
    if not isinstance(cluster_cfg, dict) or not cluster_cfg:
        errors.append(f"  {cluster_path}: cluster config is empty or not a valid YAML mapping.")
        _print_validation_errors(errors)

    scheduler = cluster_cfg.get("scheduler", {})
    if not isinstance(scheduler, dict):
        errors.append(f"  {cluster_path}: 'scheduler' must be a mapping.")
    else:
        scheduler_type = scheduler.get("type", "slurm")
        if str(scheduler_type).lower() != "slurm":
            errors.append(f"  {cluster_path}: scheduler.type must be 'slurm' in v1 (got '{scheduler_type}').")

    resources = cluster_cfg.get("resources", {})
    if not isinstance(resources, dict):
        errors.append(f"  {cluster_path}: 'resources' must be a mapping.")
    else:
        for req in ("account", "nodes", "ntasks_per_node", "mem", "time"):
            if req not in resources:
                errors.append(f"  {cluster_path}: missing required key 'resources.{req}'.")
        for int_key in ("nodes", "ntasks_per_node"):
            if int_key in resources:
                val = resources.get(int_key)
                if not isinstance(val, int) or val <= 0:
                    errors.append(f"  {cluster_path}: resources.{int_key} must be a positive integer (got {val}).")
        for str_key in ("account", "mem", "time", "partition"):
            if str_key in resources and resources.get(str_key) is not None:
                if not isinstance(resources.get(str_key), str):
                    errors.append(f"  {cluster_path}: resources.{str_key} must be a string when provided.")
        if isinstance(resources.get("time"), str):
            try:
                parse_slurm_time_limit_to_seconds(resources["time"])
            except ValueError as exc:
                errors.append(
                    f"  {cluster_path}: resources.time must be a supported finite Slurm time string ({exc})."
                )
        account = resources.get("account")
        if account == CLUSTER_TEMPLATE_PLACEHOLDER_ACCOUNT:
            warnings.append(
                f"{cluster_path}: resources.account still uses the sample placeholder "
                f"'{CLUSTER_TEMPLATE_PLACEHOLDER_ACCOUNT}'. Edit the cluster profile before submission."
            )

    notifications = cluster_cfg.get("notifications", {})
    if notifications is not None and not isinstance(notifications, dict):
        errors.append(f"  {cluster_path}: 'notifications' must be a mapping when provided.")
    elif isinstance(notifications, dict):
        mail_user = notifications.get("mail_user")
        if mail_user is not None and not is_valid_email(mail_user):
            errors.append(f"  {cluster_path}: notifications.mail_user is not a valid email '{mail_user}'.")
        if mail_user == CLUSTER_TEMPLATE_PLACEHOLDER_MAIL:
            warnings.append(
                f"{cluster_path}: notifications.mail_user still uses the sample placeholder "
                f"'{CLUSTER_TEMPLATE_PLACEHOLDER_MAIL}'. Edit the cluster profile before submission."
            )
        mail_type = notifications.get("mail_type")
        if mail_type is not None and not isinstance(mail_type, str):
            errors.append(f"  {cluster_path}: notifications.mail_type must be a string when provided.")

    execution = cluster_cfg.get("execution", {})
    if execution is not None and not isinstance(execution, dict):
        errors.append(f"  {cluster_path}: 'execution' must be a mapping when provided.")
    elif isinstance(execution, dict):
        module_setup = execution.get("module_setup", [])
        if module_setup is not None and not isinstance(module_setup, list):
            errors.append(f"  {cluster_path}: execution.module_setup must be a list of shell lines.")
        elif isinstance(module_setup, list):
            for i, line in enumerate(module_setup):
                if not isinstance(line, str):
                    errors.append(f"  {cluster_path}: execution.module_setup[{i}] must be a string.")

        launcher = execution.get("launcher")
        if launcher is not None and not isinstance(launcher, str):
            errors.append(f"  {cluster_path}: execution.launcher must be a string when provided.")
        launcher_args = execution.get("launcher_args")
        if launcher_args is not None and not isinstance(launcher_args, list):
            errors.append(f"  {cluster_path}: execution.launcher_args must be a list of CLI tokens.")
        elif isinstance(launcher_args, list):
            for i, token in enumerate(launcher_args):
                if not isinstance(token, (str, int, float)):
                    errors.append(f"  {cluster_path}: execution.launcher_args[{i}] must be a scalar CLI token.")
                elif _launcher_arg_contains_whitespace(token):
                    errors.append(
                        f"  {cluster_path}: execution.launcher_args[{i}] must be a single CLI token; "
                        "split whitespace-separated arguments into separate list items."
                    )
        if (launcher is None or isinstance(launcher, str)) and (launcher_args is None or isinstance(launcher_args, list)):
            try:
                normalize_cluster_launcher(execution)
            except ValueError as exc:
                errors.append(f"  {cluster_path}: {exc}.")

        extra_sbatch = execution.get("extra_sbatch")
        if extra_sbatch is not None and not isinstance(extra_sbatch, (dict, list)):
            errors.append(f"  {cluster_path}: execution.extra_sbatch must be a mapping or list when provided.")

        walltime_guard = execution.get("walltime_guard")
        if walltime_guard is not None and not isinstance(walltime_guard, dict):
            errors.append(f"  {cluster_path}: execution.walltime_guard must be a mapping when provided.")
        elif isinstance(walltime_guard, dict):
            enabled = walltime_guard.get("enabled")
            if enabled is not None and not isinstance(enabled, bool):
                errors.append(f"  {cluster_path}: execution.walltime_guard.enabled must be boolean when provided.")

            warmup_steps = walltime_guard.get("warmup_steps")
            if warmup_steps is not None and (not isinstance(warmup_steps, int) or isinstance(warmup_steps, bool) or warmup_steps <= 0):
                errors.append(
                    f"  {cluster_path}: execution.walltime_guard.warmup_steps must be a positive integer when provided."
                )

            multiplier = walltime_guard.get("multiplier")
            if multiplier is not None:
                if isinstance(multiplier, bool) or not isinstance(multiplier, (int, float)) or multiplier <= 0.0:
                    errors.append(
                        f"  {cluster_path}: execution.walltime_guard.multiplier must be a positive number when provided."
                    )
                elif float(multiplier) > 5.0:
                    errors.append(
                        f"  {cluster_path}: execution.walltime_guard.multiplier must be <= 5.0 (got {multiplier})."
                    )

            min_seconds = walltime_guard.get("min_seconds")
            if min_seconds is not None and (
                isinstance(min_seconds, bool) or not isinstance(min_seconds, (int, float)) or float(min_seconds) <= 0.0
            ):
                errors.append(
                    f"  {cluster_path}: execution.walltime_guard.min_seconds must be a positive number when provided."
                )

            estimator_alpha = walltime_guard.get("estimator_alpha")
            if estimator_alpha is not None:
                if isinstance(estimator_alpha, bool) or not isinstance(estimator_alpha, (int, float)):
                    errors.append(
                        f"  {cluster_path}: execution.walltime_guard.estimator_alpha must be a number in (0, 1] when provided."
                    )
                elif float(estimator_alpha) <= 0.0 or float(estimator_alpha) > 1.0:
                    errors.append(
                        f"  {cluster_path}: execution.walltime_guard.estimator_alpha must be in (0, 1] (got {estimator_alpha})."
                    )

    if warnings:
        for warning in warnings:
            print(f"[WARN] {warning}", file=sys.stderr)

    if errors:
        _print_validation_errors(errors)

def validate_study_config(study_cfg: dict, study_path: str, skip_base_file_check: bool = False):
    """!
    @brief Validate sweep/study specification from study.yml.
    @param[in] study_cfg Argument passed to `validate_study_config()`.
    @param[in] study_path Argument passed to `validate_study_config()`.
    @param[in] skip_base_file_check When True, skip file-existence check for base_configs paths.
    """
    errors = []
    _validate_yaml_schema_keys(study_cfg, _STUDY_SCHEMA, study_path, errors)
    if not isinstance(study_cfg, dict) or not study_cfg:
        errors.append(f"  {study_path}: study config is empty or not a valid YAML mapping.")
        _print_validation_errors(errors)

    base_cfgs = study_cfg.get("base_configs")
    if not isinstance(base_cfgs, dict):
        errors.append(f"  {study_path}: missing required mapping 'base_configs'.")
    else:
        for req in ("case", "solver", "monitor", "post"):
            path_val = base_cfgs.get(req)
            if not path_val or not isinstance(path_val, str):
                errors.append(f"  {study_path}: base_configs.{req} must be a path string.")
            elif not skip_base_file_check:
                resolved = resolve_path(study_path, path_val)
                if not os.path.isfile(resolved):
                    errors.append(f"  {study_path}: base_configs.{req} does not exist: {resolved}")

    study_type = study_cfg.get("study_type")
    allowed_types = {"grid_independence", "timestep_independence", "sensitivity"}
    if study_type not in allowed_types:
        errors.append(
            f"  {study_path}: study_type must be one of {sorted(allowed_types)} (got '{study_type}')."
        )

    parameters = study_cfg.get("parameters")
    parameter_sets = study_cfg.get("parameter_sets")
    allowed_roots = {"case", "solver", "monitor", "post"}
    if bool(parameters) == bool(parameter_sets):
        errors.append(f"  {study_path}: provide exactly one of 'parameters' or 'parameter_sets'.")
    elif parameter_sets:
        if not isinstance(parameter_sets, list) or not parameter_sets:
            errors.append(f"  {study_path}: 'parameter_sets' must be a non-empty list of key->value mappings.")
        else:
            for set_index, param_set in enumerate(parameter_sets):
                if not isinstance(param_set, dict) or not param_set:
                    errors.append(
                        f"  {study_path}: parameter_sets[{set_index}] must be a non-empty mapping of key->value overrides."
                    )
                    continue
                for key, value in param_set.items():
                    if not isinstance(key, str) or "." not in key:
                        errors.append(
                            f"  {study_path}: parameter_sets[{set_index}] key '{key}' must use '<target>.<yaml.path>' format."
                        )
                        continue
                    root = key.split(".", 1)[0]
                    if root not in allowed_roots:
                        errors.append(
                            f"  {study_path}: parameter_sets[{set_index}] key '{key}' must start with one of {sorted(allowed_roots)}."
                        )
                    if isinstance(value, (dict, list)):
                        errors.append(
                            f"  {study_path}: parameter_sets[{set_index}] value for '{key}' must be a scalar, not {type(value).__name__}."
                        )
    else:
        if not isinstance(parameters, dict) or not parameters:
            errors.append(f"  {study_path}: 'parameters' must be a non-empty mapping of key->list.")
        else:
            for key, values in parameters.items():
                if not isinstance(key, str) or "." not in key:
                    errors.append(
                        f"  {study_path}: parameter key '{key}' must use '<target>.<yaml.path>' format."
                    )
                    continue
                root = key.split(".", 1)[0]
                if root not in allowed_roots:
                    errors.append(
                        f"  {study_path}: parameter key '{key}' must start with one of {sorted(allowed_roots)}."
                    )
                if not isinstance(values, list) or len(values) == 0:
                    errors.append(f"  {study_path}: parameters.{key} must be a non-empty list.")

    metrics = study_cfg.get("metrics", [])
    if metrics is not None and not isinstance(metrics, list):
        errors.append(f"  {study_path}: 'metrics' must be a list when provided.")
    elif isinstance(metrics, list):
        for i, metric in enumerate(metrics):
            if isinstance(metric, str):
                continue
            if not isinstance(metric, dict):
                errors.append(
                    f"  {study_path}: metrics[{i}] must be a string or mapping."
                )
                continue
            if "name" not in metric:
                errors.append(f"  {study_path}: metrics[{i}] missing required key 'name'.")
            if "source" not in metric:
                errors.append(f"  {study_path}: metrics[{i}] missing required key 'source'.")

    plotting = study_cfg.get("plotting", {})
    if plotting is not None and not isinstance(plotting, dict):
        errors.append(f"  {study_path}: 'plotting' must be a mapping when provided.")
    elif isinstance(plotting, dict):
        enabled = plotting.get("enabled")
        if enabled is not None and not isinstance(enabled, bool):
            errors.append(f"  {study_path}: plotting.enabled must be boolean when provided.")
        output_format = plotting.get("output_format")
        if output_format is not None and output_format not in {"png", "pdf", "svg"}:
            errors.append(f"  {study_path}: plotting.output_format must be one of ['png','pdf','svg'].")

    execution = study_cfg.get("execution", {})
    if execution is not None and not isinstance(execution, dict):
        errors.append(f"  {study_path}: 'execution' must be a mapping when provided.")
    elif isinstance(execution, dict):
        max_conc = execution.get("max_concurrent_array_tasks")
        if max_conc is not None and (not isinstance(max_conc, int) or max_conc <= 0):
            errors.append(
                f"  {study_path}: execution.max_concurrent_array_tasks must be a positive integer when provided."
            )

    if errors:
        _print_validation_errors(errors)

def _deep_set(container: dict, dotted_path: str, value):
    """!
    @brief Set nested dictionary value, creating intermediate maps when needed.
    @param[in] container Argument passed to `_deep_set()`.
    @param[in] dotted_path Argument passed to `_deep_set()`.
    @param[in] value Argument passed to `_deep_set()`.
    """
    keys = dotted_path.split(".")
    current = container
    for key in keys[:-1]:
        if key not in current or not isinstance(current[key], dict):
            current[key] = {}
        current = current[key]
    current[keys[-1]] = value

def expand_parameter_matrix(parameters: dict) -> list:
    """!
    @brief Expand study parameter lists into cartesian-product combinations.
    @param[in] parameters Argument passed to `expand_parameter_matrix()`.
    @return Value returned by `expand_parameter_matrix()`.
    """
    param_keys = list(parameters.keys())
    all_values = [parameters[k] for k in param_keys]
    combos = []
    for combo in itertools.product(*all_values):
        combos.append(dict(zip(param_keys, combo)))
    return combos


def expand_study_parameter_combinations(study_cfg: dict) -> list:
    """!
    @brief Expand either cartesian-study parameters or explicit parameter sets.
    @param[in] study_cfg Argument passed to `expand_study_parameter_combinations()`.
    @return Value returned by `expand_study_parameter_combinations()`.
    """
    parameter_sets = study_cfg.get("parameter_sets")
    if parameter_sets:
        return [dict(param_set) for param_set in parameter_sets]
    return expand_parameter_matrix(study_cfg.get("parameters") or {})


def get_study_parameter_keys(study_cfg: dict) -> list:
    """!
    @brief Collect ordered parameter keys from either cross-product parameter expansions or explicit parameter sets.
    @param[in] study_cfg Argument passed to `get_study_parameter_keys()`.
    @return Value returned by `get_study_parameter_keys()`.
    """
    parameters = study_cfg.get("parameters")
    if isinstance(parameters, dict) and parameters:
        return list(parameters.keys())

    keys = []
    parameter_sets = study_cfg.get("parameter_sets") or []
    for param_set in parameter_sets:
        if not isinstance(param_set, dict):
            continue
        for key in param_set.keys():
            if key not in keys:
                keys.append(key)
    return keys


def get_cluster_total_tasks(cluster_cfg: dict) -> int:
    """!
    @brief Return cluster total tasks.
    @param[in] cluster_cfg Argument passed to `get_cluster_total_tasks()`.
    @return Value returned by `get_cluster_total_tasks()`.
    """
    resources = cluster_cfg.get("resources", {})
    return int(resources.get("nodes", 1)) * int(resources.get("ntasks_per_node", 1))

def normalize_extension(ext: str) -> str:
    """!
    @brief Normalize extension.
    @param[in] ext Argument passed to `normalize_extension()`.
    @return Value returned by `normalize_extension()`.
    """
    if ext is None:
        return None
    return str(ext).strip().lstrip(".")

def render_slurm_script(
    script_path: str,
    job_name: str,
    cluster_cfg: dict,
    command: list,
    workdir: str,
    stdout_path: str,
    stderr_path: str = None,
    env_vars: dict = None,
    shell_env_vars: dict = None,
    array_spec: str = None
):
    """!
    @brief Render a Slurm batch script for a single command.
    @param[in] script_path Argument passed to `render_slurm_script()`.
    @param[in] job_name Argument passed to `render_slurm_script()`.
    @param[in] cluster_cfg Argument passed to `render_slurm_script()`.
    @param[in] command Argument passed to `render_slurm_script()`.
    @param[in] workdir Argument passed to `render_slurm_script()`.
    @param[in] stdout_path Argument passed to `render_slurm_script()`.
    @param[in] stderr_path Argument passed to `render_slurm_script()`.
    @param[in] env_vars Argument passed to `render_slurm_script()`.
    @param[in] shell_env_vars Argument passed to `render_slurm_script()`.
    @param[in] array_spec Argument passed to `render_slurm_script()`.
    """
    resources = cluster_cfg.get("resources", {})
    notifications = cluster_cfg.get("notifications", {}) or {}
    execution = cluster_cfg.get("execution", {}) or {}
    extra_sbatch = execution.get("extra_sbatch")
    module_setup = execution.get("module_setup", []) or []

    if stderr_path is None:
        stderr_path = stdout_path.replace(".out", ".err")

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --nodes={resources['nodes']}",
        f"#SBATCH --ntasks-per-node={resources['ntasks_per_node']}",
        f"#SBATCH --mem={resources['mem']}",
        f"#SBATCH --time={resources['time']}",
        f"#SBATCH --output={stdout_path}",
        f"#SBATCH --error={stderr_path}",
        f"#SBATCH --account={resources['account']}",
    ]
    partition = resources.get("partition")
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    if array_spec:
        lines.append(f"#SBATCH --array={array_spec}")
    mail_user = notifications.get("mail_user")
    mail_type = notifications.get("mail_type")
    if mail_user:
        lines.append(f"#SBATCH --mail-user={mail_user}")
    if mail_type:
        lines.append(f"#SBATCH --mail-type={mail_type}")

    if isinstance(extra_sbatch, dict):
        for key, value in extra_sbatch.items():
            flag = str(key)
            if not flag.startswith("--"):
                flag = f"--{flag}"
            if isinstance(value, bool):
                if value:
                    lines.append(f"#SBATCH {flag}")
            elif value is not None:
                lines.append(f"#SBATCH {flag}={value}")
    elif isinstance(extra_sbatch, list):
        for token in extra_sbatch:
            lines.append(f"#SBATCH {token}")

    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"cd {shlex.quote(workdir)}",
            'echo "[$(date)] Starting job ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"',
            'echo "[$(date)] Working directory: $PWD"',
        ]
    )

    if shell_env_vars:
        for key, value in shell_env_vars.items():
            lines.append(f"export {key}={value}")

    for setup_line in module_setup:
        lines.append(str(setup_line))

    if env_vars:
        for key, value in env_vars.items():
            lines.append(f"export {key}={shlex.quote(str(value))}")

    cmd = " ".join(shlex.quote(str(tok)) for tok in command)
    lines.append(f"exec {cmd}")

    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)

def split_launcher_tokens(
    launcher: "str | None",
    launcher_args: "list | None" = None,
    label: str = "launcher",
) -> "tuple[str | None, list[str]]":
    """!
    @brief Canonicalize launcher config into executable token plus argv-style flags.
    @param[in] launcher Argument passed to `split_launcher_tokens()`.
    @param[in] launcher_args Argument passed to `split_launcher_tokens()`.
    @param[in] label Argument passed to `split_launcher_tokens()`.
    @return Value returned by `split_launcher_tokens()`.
    """
    normalized_args = [str(x) for x in (launcher_args or [])]

    if launcher is None:
        return None, normalized_args

    try:
        launcher_tokens = shlex.split(str(launcher))
    except ValueError as exc:
        raise ValueError(f"{label} is not shell-parseable: {exc}") from exc

    if not launcher_tokens:
        return None, normalized_args

    return launcher_tokens[0], launcher_tokens[1:] + normalized_args


def normalize_cluster_launcher(execution: dict) -> "tuple[str | None, list[str]]":
    """!
    @brief Canonicalize cluster launcher config into executable token plus argv-style flags.
    @param[in] execution Argument passed to `normalize_cluster_launcher()`.
    @return Value returned by `normalize_cluster_launcher()`.
    """
    return split_launcher_tokens(
        execution.get("launcher"),
        execution.get("launcher_args") or [],
        label="execution.launcher",
    )


def strip_launcher_size_flags(launcher_name: str, launcher_args: "list[str]") -> "list[str]":
    """!
    @brief Remove explicit MPI task-count flags from known launchers.
    @param[in] launcher_name Basename-normalized launcher executable.
    @param[in] launcher_args Launcher argument list.
    @return Filtered launcher arguments with explicit size flags removed.
    """
    filtered = []
    idx = 0
    while idx < len(launcher_args):
        token = str(launcher_args[idx])

        if launcher_name == "srun":
            if token in {"-n", "--ntasks"}:
                idx += 2
                continue
            if token.startswith("--ntasks="):
                idx += 1
                continue
        elif launcher_name in {"mpiexec", "mpirun"}:
            if token in {"-n", "-np"}:
                idx += 2
                continue
            if token.startswith("-n=") or token.startswith("-np="):
                idx += 1
                continue

        filtered.append(token)
        idx += 1

    return filtered


def build_serial_post_cluster_config(cluster_cfg: dict, num_procs: int = 1) -> dict:
    """!
    @brief Clone cluster config and force a single-node post stage task layout.
    @param[in] cluster_cfg Base cluster configuration.
    @param[in] num_procs Number of post tasks to request.
    @return Cluster configuration specialized for the post stage.
    """
    post_cluster_cfg = copy.deepcopy(cluster_cfg)
    post_cluster_cfg.setdefault("resources", {})
    post_cluster_cfg["resources"]["nodes"] = 1
    post_cluster_cfg["resources"]["ntasks_per_node"] = int(num_procs)
    return post_cluster_cfg


def build_local_launch_command(
    executable: str,
    executable_args: list,
    num_procs: int,
    config_search_anchor: str = None,
    allow_single_rank_launcher_override: bool = False,
    force_num_procs: "int | None" = None,
) -> list:
    """!
    @brief Build local launcher command, allowing env or shared config overrides for login-node MPI quirks.
    @param[in] executable Argument passed to `build_local_launch_command()`.
    @param[in] executable_args Argument passed to `build_local_launch_command()`.
    @param[in] num_procs Argument passed to `build_local_launch_command()`.
    @param[in] config_search_anchor Argument passed to `build_local_launch_command()`.
    @param[in] allow_single_rank_launcher_override When true, explicit launcher overrides also apply to 1-rank commands.
    @param[in] force_num_procs Optional explicit MPI rank count override applied after stripping conflicting launcher size flags.
    @return Value returned by `build_local_launch_command()`.
    """
    target_num_procs = force_num_procs if force_num_procs is not None else num_procs
    command = [executable] + executable_args
    if target_num_procs <= 1 and not allow_single_rank_launcher_override:
        return command

    launcher_override = os.environ.get("PICURV_MPI_LAUNCHER")
    if launcher_override is None:
        launcher_override = os.environ.get("MPI_LAUNCHER")

    try:
        if launcher_override is not None:
            explicit_launcher_config = True
            launcher, launcher_args = split_launcher_tokens(
                launcher_override,
                label="local MPI launcher override",
            )
        else:
            _, runtime_execution_cfg = load_runtime_execution_config(config_search_anchor)
            local_execution = resolve_runtime_execution_context(runtime_execution_cfg, "local")
            configured_launcher = local_execution.get("launcher")
            configured_args = local_execution.get("launcher_args") or []
            explicit_launcher_config = configured_launcher is not None or bool(configured_args)
            if target_num_procs <= 1 and not explicit_launcher_config:
                return command
            launcher, launcher_args = split_launcher_tokens(
                configured_launcher if configured_launcher is not None else "mpiexec",
                configured_args,
                label="local_execution.launcher",
            )
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    if not launcher:
        return command

    launcher_name = os.path.basename(launcher).lower()
    if force_num_procs is not None:
        launcher_args = strip_launcher_size_flags(launcher_name, launcher_args)
    prefix = [launcher] + launcher_args

    if launcher_name == "srun":
        has_n = any(token in {"-n", "--ntasks"} for token in launcher_args)
        if not has_n:
            prefix += ["-n", str(target_num_procs)]
    elif launcher_name in {"mpiexec", "mpirun"}:
        has_n = any(token in {"-n", "-np"} for token in launcher_args)
        if not has_n:
            prefix += ["-n", str(target_num_procs)]

    return prefix + command

def resolve_cluster_execution(cluster_cfg: dict, config_search_anchor: str = None, extra_search_anchors=None) -> dict:
    """!
    @brief Resolve cluster execution launcher settings from shared runtime config plus cluster.yml overrides.
    @param[in] cluster_cfg Argument passed to `resolve_cluster_execution()`.
    @param[in] config_search_anchor Argument passed to `resolve_cluster_execution()`.
    @param[in] extra_search_anchors Argument passed to `resolve_cluster_execution()`.
    @return Value returned by `resolve_cluster_execution()`.
    """
    _, runtime_execution_cfg = load_runtime_execution_config(config_search_anchor, extra_search_anchors=extra_search_anchors)
    shared_cluster_execution = resolve_runtime_execution_context(runtime_execution_cfg, "cluster")
    execution = cluster_cfg.get("execution", {}) or {}
    cluster_override = {
        "launcher": execution.get("launcher") if "launcher" in execution else None,
        "launcher_args": execution.get("launcher_args") if "launcher_args" in execution else None,
    }
    return merge_execution_overrides(shared_cluster_execution, cluster_override)


def build_cluster_launch_command(
    cluster_cfg: dict,
    executable: str,
    executable_args: list,
    config_search_anchor: str = None,
    extra_search_anchors=None,
    force_num_procs: "int | None" = None,
) -> list:
    """!
    @brief Build scheduler launcher command from cluster config plus optional shared execution defaults.
    @param[in] cluster_cfg Argument passed to `build_cluster_launch_command()`.
    @param[in] executable Argument passed to `build_cluster_launch_command()`.
    @param[in] executable_args Argument passed to `build_cluster_launch_command()`.
    @param[in] config_search_anchor Argument passed to `build_cluster_launch_command()`.
    @param[in] extra_search_anchors Argument passed to `build_cluster_launch_command()`.
    @param[in] force_num_procs Optional explicit MPI rank count override applied after stripping conflicting launcher size flags.
    @return Value returned by `build_cluster_launch_command()`.
    """
    try:
        execution = resolve_cluster_execution(
            cluster_cfg,
            config_search_anchor=config_search_anchor,
            extra_search_anchors=extra_search_anchors,
        )
        launcher, launcher_args = split_launcher_tokens(
            execution.get("launcher") if execution.get("launcher") is not None else "srun",
            execution.get("launcher_args") or [],
            label="cluster execution launcher",
        )
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    ntasks = int(force_num_procs) if force_num_procs is not None else get_cluster_total_tasks(cluster_cfg)
    launcher_name = launcher.lower() if launcher else ""
    if force_num_procs is not None:
        launcher_args = strip_launcher_size_flags(launcher_name, launcher_args)

    if launcher and launcher_name == "srun":
        has_n = any(token in {"-n", "--ntasks"} for token in launcher_args)
        cmd = ["srun"] + launcher_args
        if not has_n:
            cmd += ["-n", str(ntasks)]
        return cmd + [executable] + executable_args

    if launcher and launcher_name == "mpirun":
        has_np = any(token in {"-np", "-n"} for token in launcher_args)
        cmd = ["mpirun"] + launcher_args
        if not has_np:
            cmd += ["-np", str(ntasks)]
        return cmd + [executable] + executable_args

    if launcher and launcher_name == "mpiexec":
        has_np = any(token in {"-np", "-n"} for token in launcher_args)
        cmd = ["mpiexec"] + launcher_args
        if not has_np:
            cmd += ["-np", str(ntasks)]
        return cmd + [executable] + executable_args

    # Custom launcher or no launcher.
    cmd = []
    if launcher:
        cmd.append(str(launcher))
    cmd += launcher_args
    cmd += [executable] + executable_args
    return cmd

def parse_slurm_job_id(sbatch_output: str) -> str:
    """!
    @brief Extract numeric job id from standard sbatch output.
    @param[in] sbatch_output Argument passed to `parse_slurm_job_id()`.
    @return Value returned by `parse_slurm_job_id()`.
    """
    match = re.search(r"Submitted batch job\s+(\d+)", sbatch_output or "")
    return match.group(1) if match else None

def submit_sbatch(script_path: str, dependency: str = None, dependency_type: str = "afterok") -> dict:
    """!
    @brief Submit sbatch script and return submission metadata.
    @param[in] script_path Argument passed to `submit_sbatch()`.
    @param[in] dependency Argument passed to `submit_sbatch()`.
    @param[in] dependency_type Slurm dependency type (default: afterok). Common values: afterok, afterany.
    @return Value returned by `submit_sbatch()`.
    """
    cmd = ["sbatch"]
    if dependency:
        cmd.append(f"--dependency={dependency_type}:{dependency}")
    cmd.append(script_path)
    result = subprocess.run(cmd, text=True, capture_output=True, check=False)
    metadata = {
        "command": cmd,
        "returncode": result.returncode,
        "stdout": (result.stdout or "").strip(),
        "stderr": (result.stderr or "").strip(),
        "script": script_path,
    }
    if result.returncode != 0:
        print(f"[FATAL] sbatch submission failed for {script_path}\n{metadata['stderr']}", file=sys.stderr)
        sys.exit(result.returncode)
    metadata["job_id"] = parse_slurm_job_id(metadata["stdout"])
    if not metadata["job_id"]:
        print(
            f"[FATAL] Could not parse Slurm job id from sbatch output: {metadata['stdout']}",
            file=sys.stderr
        )
        sys.exit(1)
    return metadata


def _print_validation_errors(errors: list):
    """!
    @brief Prints validation errors and exits.
    @param[in] errors List of error message strings.
    """
    print(f"\n[FATAL] Configuration validation failed with {len(errors)} issue(s):", file=sys.stderr)
    for raw_error in errors:
        file_path, message = _split_error_file_and_message(raw_error)
        key_path = _extract_key_path(message)
        code = _classify_error_code(message)
        emit_structured_error(code, key=key_path, file_path=file_path, message=message)
    print(
        "\nHint: See examples/master_template/ for valid config structure and "
        "docs/pages/14_Config_Contract.md for key-level contract details.",
        file=sys.stderr,
    )
    sys.exit(1)


def generate_header(run_id: str, source_files: dict) -> str:
    """!
    @brief Creates a standard header block for all generated files.
    @param[in] run_id The unique identifier for the current simulation run.
    @param[in] source_files A dictionary of source profile files used.
    @return A formatted string containing the header.
    """
    header_parts = [
        "# ==============================================================================",
        "#                AUTO-GENERATED CONFIGURATION FILE",
        "# ------------------------------------------------------------------------------",
        f"#   Run ID:       {run_id}",
        f"#   Timestamp:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "#",
        "#   Source Configuration:"
    ]
    for name, path in source_files.items():
        header_parts.append(f"#     - {name:<12}: {os.path.basename(path)}")
    header_parts.extend([
        "#",
        "#   DO NOT EDIT THIS FILE MANUALLY. IT IS A MACHINE-READABLE ARTIFACT.",
        "# ==============================================================================\n"
    ])
    return "\n".join(header_parts)

def generate_simple_list_file(run_dir: str, run_id: str, cfg: dict, section: str, key: str, filename: str, header_sources: dict) -> str:
    """!
    @brief Generic function to create a file containing a simple list of strings.
    @param[in] run_dir The path to the main run directory.
    @param[in] run_id The unique identifier for the run.
    @param[in] cfg The dictionary containing the configuration data.
    @param[in] section The top-level key in the cfg dictionary.
    @param[in] key The second-level key whose value is the list of strings.
    @param[in] filename The name of the file to generate (e.g., 'whitelist.run').
    @param[in] header_sources A dictionary of source files for the header.
    @return The absolute path to the generated file.
    """
    print(f"[INFO] Generating {filename}...")
    config_dir = os.path.join(run_dir, "config")
    file_path = os.path.join(config_dir, filename)
    
    lines = [generate_header(run_id, header_sources)]
    items = cfg.get(section, {}).get(key, [])
    lines.extend(items)

    with open(file_path, "w") as f: f.write("\n".join(lines))
    print(f"[SUCCESS] Generated {filename}: {os.path.relpath(file_path)}")
    return os.path.abspath(file_path)


def has_explicit_monitor_whitelist(monitor_cfg: dict) -> bool:
    """!
    @brief Return True when logging.enabled_functions contains at least one entry.
    @param[in] monitor_cfg Argument passed to `has_explicit_monitor_whitelist()`.
    @return Value returned by `has_explicit_monitor_whitelist()`.
    """
    items = monitor_cfg.get("logging", {}).get("enabled_functions", [])
    return bool(items)


def resolve_profiling_config(monitor_cfg: dict) -> dict:
    """!
    @brief Resolve profiling reporting config from monitor.yml.
    @param[in] monitor_cfg Argument passed to `resolve_profiling_config()`.
    @return Value returned by `resolve_profiling_config()`.
    """
    profiling_cfg = monitor_cfg.get("profiling", {}) or {}
    timestep_cfg = profiling_cfg.get("timestep_output")
    final_cfg = profiling_cfg.get("final_summary")

    if timestep_cfg is None:
        mode = "off"
        functions = []
        timestep_file = "Profiling_Timestep_Summary.csv"
    else:
        if not isinstance(timestep_cfg, dict):
            raise ValueError("monitor.profiling.timestep_output must be a mapping when provided.")
        mode = str(timestep_cfg.get("mode", "off")).lower()
        functions = timestep_cfg.get("functions", [])
        timestep_file = str(timestep_cfg.get("file", "Profiling_Timestep_Summary.csv"))

    if mode not in {"off", "selected", "all"}:
        raise ValueError("monitor.profiling.timestep_output.mode must be one of ['off', 'selected', 'all'].")
    if functions is None:
        functions = []
    if not isinstance(functions, list):
        raise ValueError("monitor.profiling.timestep_output.functions must be a list of function names.")
    if not all(isinstance(item, str) and item.strip() for item in functions):
        raise ValueError("monitor.profiling.timestep_output.functions entries must be non-empty strings.")
    if mode == "selected" and not functions:
        raise ValueError("monitor.profiling.timestep_output.functions must be non-empty when mode is 'selected'.")
    if mode != "selected" and functions:
        raise ValueError("monitor.profiling.timestep_output.functions is only valid when mode is 'selected'.")
    if not timestep_file:
        raise ValueError("monitor.profiling.timestep_output.file must be a non-empty string.")

    if final_cfg is None:
        final_enabled = True
    elif isinstance(final_cfg, dict):
        final_enabled = bool(final_cfg.get("enabled", True))
    else:
        raise ValueError("monitor.profiling.final_summary must be a mapping when provided.")

    return {
        "mode": mode,
        "functions": functions,
        "timestep_file": timestep_file,
        "final_summary_enabled": final_enabled,
    }


DIAGNOSTICS_PETSC_KEYS = {
    "malloc_debug",
    "malloc_test",
    "malloc_dump",
    "malloc_view",
    "malloc_view_threshold",
    "memory_view",
    "log_view",
    "log_view_memory",
    "log_all",
    "log_trace",
    "objects_dump",
    "options_left",
}


def _diagnostic_bool_or_path(value, key: str):
    """!
    @brief Validate a diagnostics value that can be false, true, or a path/viewer string.
    @param[in] value Candidate value.
    @param[in] key Diagnostics key used in error messages.
    @return Normalized value.
    """
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, str) and value.strip():
        return value.strip()
    raise ValueError(f"monitor.diagnostics.petsc.{key} must be boolean, null, or a non-empty string.")


def _diagnostic_bool(value, key: str) -> bool:
    """!
    @brief Validate a diagnostics boolean value.
    @param[in] value Candidate value.
    @param[in] key Diagnostics key used in error messages.
    @return Boolean value.
    """
    if isinstance(value, bool):
        return value
    raise ValueError(f"monitor.diagnostics.petsc.{key} must be boolean.")


def _diagnostic_bool_or_all(value, key: str):
    """!
    @brief Validate a diagnostics value that can be false, true, or "all".
    @param[in] value Candidate value.
    @param[in] key Diagnostics key used in error messages.
    @return Normalized value.
    """
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, str) and value.strip().lower() == "all":
        return "all"
    raise ValueError(f"monitor.diagnostics.petsc.{key} must be boolean, null, or 'all'.")


def _diagnostic_default_file(run_dir: str, filename: str) -> str:
    """!
    @brief Return an absolute run-local diagnostics file path.
    @param[in] run_dir Run directory.
    @param[in] filename Diagnostics filename.
    @return Absolute diagnostics path under the run logs directory.
    """
    return os.path.abspath(os.path.join(run_dir, "logs", filename))


def _diagnostic_resolve_path_or_default(value, run_dir: str, default_filename: str):
    """!
    @brief Resolve true/string diagnostics values to a concrete file path.
    @param[in] value Boolean/string diagnostics value.
    @param[in] run_dir Run directory.
    @param[in] default_filename Default file name when value is true.
    @return False, or an absolute/explicit path string.
    """
    if value is True:
        return _diagnostic_default_file(run_dir, default_filename)
    if isinstance(value, str):
        if os.path.isabs(value) or value.startswith(":"):
            return value
        return os.path.abspath(os.path.join(run_dir, "logs", value))
    return False


def resolve_diagnostics_config(monitor_cfg: dict, run_dir: "str | None" = None, stage_label: str = "Solver") -> dict:
    """!
    @brief Resolve monitor diagnostics config and default run-local log paths.
    @param[in] monitor_cfg Parsed monitor.yml mapping.
    @param[in] run_dir Optional run directory for default artifact paths.
    @param[in] stage_label Solver/PostProcessor suffix used for PETSc output defaults.
    @return Normalized diagnostics config.
    """
    diagnostics_cfg = (monitor_cfg.get("diagnostics", {}) or {}) if isinstance(monitor_cfg, dict) else {}
    if not isinstance(diagnostics_cfg, dict):
        raise ValueError("monitor.diagnostics must be a mapping when provided.")

    petsc_raw = diagnostics_cfg.get("petsc", {}) or {}
    if not isinstance(petsc_raw, dict):
        raise ValueError("monitor.diagnostics.petsc must be a mapping when provided.")
    unknown = sorted(set(petsc_raw.keys()) - DIAGNOSTICS_PETSC_KEYS)
    if unknown:
        raise ValueError(f"monitor.diagnostics.petsc has unsupported key(s): {unknown}.")

    petsc = {
        "malloc_debug": _diagnostic_bool(petsc_raw.get("malloc_debug", False), "malloc_debug"),
        "malloc_test": _diagnostic_bool(petsc_raw.get("malloc_test", False), "malloc_test"),
        "malloc_dump": _diagnostic_bool(petsc_raw.get("malloc_dump", False), "malloc_dump"),
        "malloc_view": _diagnostic_bool_or_path(petsc_raw.get("malloc_view", False), "malloc_view"),
        "malloc_view_threshold": petsc_raw.get("malloc_view_threshold"),
        "memory_view": _diagnostic_bool(petsc_raw.get("memory_view", False), "memory_view"),
        "log_view": _diagnostic_bool_or_path(petsc_raw.get("log_view", False), "log_view"),
        "log_view_memory": _diagnostic_bool(petsc_raw.get("log_view_memory", False), "log_view_memory"),
        "log_all": _diagnostic_bool(petsc_raw.get("log_all", False), "log_all"),
        "log_trace": _diagnostic_bool_or_path(petsc_raw.get("log_trace", False), "log_trace"),
        "objects_dump": _diagnostic_bool_or_all(petsc_raw.get("objects_dump", False), "objects_dump"),
        "options_left": petsc_raw.get("options_left"),
    }
    if petsc["malloc_view_threshold"] is not None and not isinstance(petsc["malloc_view_threshold"], (int, float)):
        raise ValueError("monitor.diagnostics.petsc.malloc_view_threshold must be numeric or null.")
    if petsc["options_left"] is not None and not isinstance(petsc["options_left"], bool):
        raise ValueError("monitor.diagnostics.petsc.options_left must be boolean or null.")

    memory_raw = diagnostics_cfg.get("runtime_memory_log", {}) or {}
    if not isinstance(memory_raw, dict):
        raise ValueError("monitor.diagnostics.runtime_memory_log must be a mapping when provided.")
    memory_unknown = sorted(set(memory_raw.keys()) - {"enabled", "file"})
    if memory_unknown:
        raise ValueError(f"monitor.diagnostics.runtime_memory_log has unsupported key(s): {memory_unknown}.")
    memory_enabled = memory_raw.get("enabled", True)
    if not isinstance(memory_enabled, bool):
        raise ValueError("monitor.diagnostics.runtime_memory_log.enabled must be boolean.")
    memory_file = str(memory_raw.get("file", "Runtime_Memory.log")).strip()
    if not memory_file:
        raise ValueError("monitor.diagnostics.runtime_memory_log.file must be a non-empty string.")

    resolved_petsc = dict(petsc)
    artifacts = []
    if run_dir:
        suffix = "PostProcessor" if stage_label == "PostProcessor" else "Solver"
        defaults = {
            "malloc_view": f"PETSc_MallocView_{suffix}.log",
            "log_view": f"PETSc_LogView_{suffix}.log",
            "log_trace": f"PETSc_LogTrace_{suffix}.log",
        }
        for key, default_name in defaults.items():
            resolved_value = _diagnostic_resolve_path_or_default(petsc[key], run_dir, default_name)
            if key == "log_view" and resolved_value and isinstance(resolved_value, str) and not resolved_value.startswith(":"):
                resolved_value = f":{resolved_value}"
            resolved_petsc[key] = resolved_value
            if resolved_value and isinstance(resolved_value, str) and not resolved_value.startswith(":"):
                artifacts.append(resolved_value)
            elif resolved_value and isinstance(resolved_value, str) and resolved_value.startswith(":"):
                artifacts.append(resolved_value[1:])
        if memory_enabled:
            artifacts.append(os.path.abspath(os.path.join(run_dir, "logs", memory_file)))

    return {
        "petsc": resolved_petsc,
        "runtime_memory_log": {"enabled": memory_enabled, "file": memory_file},
        "artifacts": artifacts,
    }


def build_petsc_diagnostics_args(monitor_cfg: dict, run_dir: str, stage_label: str) -> list:
    """!
    @brief Build PETSc diagnostics command-line arguments for a run stage.
    @param[in] monitor_cfg Parsed monitor.yml mapping.
    @param[in] run_dir Run directory used to resolve default diagnostics files.
    @param[in] stage_label Stage label for default output names.
    @return List of executable arguments.
    """
    diagnostics = resolve_diagnostics_config(monitor_cfg, run_dir, stage_label)
    petsc = diagnostics["petsc"]
    args = []
    if petsc["malloc_debug"]:
        args.append("-malloc_debug")
    if petsc["malloc_test"]:
        args.append("-malloc_test")
    for key, flag in (
        ("malloc_dump", "-malloc_dump"),
        ("malloc_view", "-malloc_view"),
        ("memory_view", "-memory_view"),
        ("log_view", "-log_view"),
        ("log_trace", "-log_trace"),
        ("objects_dump", "-objects_dump"),
    ):
        value = petsc.get(key)
        if value is True:
            args.append(flag)
        elif value:
            args.extend([flag, str(value)])
    if petsc["malloc_view_threshold"] is not None:
        args.extend(["-malloc_view_threshold", str(petsc["malloc_view_threshold"])])
    if petsc["log_view_memory"]:
        args.append("-log_view_memory")
    if petsc["log_all"]:
        args.append("-log_all")
    if petsc["options_left"] is not None:
        args.extend(["-options_left", "true" if petsc["options_left"] else "false"])
    return args


def prepare_monitor_files(run_dir: str, run_id: str, monitor_cfg: dict, source_files: dict) -> dict:
    """!
    @brief Generate monitor sidecar files and resolve profiling reporting behavior.
    @param[in] run_dir Argument passed to `prepare_monitor_files()`.
    @param[in] run_id Argument passed to `prepare_monitor_files()`.
    @param[in] monitor_cfg Argument passed to `prepare_monitor_files()`.
    @param[in] source_files Argument passed to `prepare_monitor_files()`.
    @return Value returned by `prepare_monitor_files()`.
    """
    print("[INFO] Generating monitoring files...")

    whitelist_path = None
    if has_explicit_monitor_whitelist(monitor_cfg):
        whitelist_path = generate_simple_list_file(
            run_dir, run_id, monitor_cfg, "logging", "enabled_functions", "whitelist.run", source_files
        )
    else:
        print("[INFO] logging.enabled_functions is empty; omitting whitelist.run so the C runtime uses its default allow-list.")

    profiling_cfg = resolve_profiling_config(monitor_cfg)

    profile_path = None
    if profiling_cfg["mode"] == "selected":
        profile_path = generate_simple_list_file(
            run_dir,
            run_id,
            {"profiling": {"selected_functions": profiling_cfg["functions"]}},
            "profiling",
            "selected_functions",
            "profile.run",
            source_files,
        )
    else:
        print(f"[INFO] profiling.timestep_output.mode is '{profiling_cfg['mode']}'; no profile.run function list is needed.")

    return {"whitelist": whitelist_path, "profile": profile_path, "profiling": profiling_cfg}

def generate_multi_block_bcs(run_dir: str, run_id: str, case_cfg: dict, source_files: dict) -> list:
    """!
    @brief Parses multi-block BCs from YAML, generates a .run file for each block,
           and returns a list of their absolute paths.
    @details Handles both simple list format (for single-block cases) and a
             list-of-lists (for multi-block cases) for boundary conditions.
    @param[in] run_dir The path to the main run directory.
    @param[in] run_id The unique identifier for the run.
    @param[in] case_cfg The parsed case.yml configuration dictionary.
    @param[in] source_files A dictionary of source files for the header.
    @return A list of absolute paths to the generated BC files.
    @throws ValueError if the number of BC definitions does not match the number of blocks.
    """
    print("[INFO] Generating boundary condition files...")
    config_dir = os.path.join(run_dir, "config")
    num_blocks = int(case_cfg.get('models', {}).get('domain', {}).get('blocks', 1))
    prepared_blocks = validate_and_prepare_boundary_conditions(case_cfg)
    case_path = source_files.get("Case") if source_files else None
    profile_grid_dims = None
    scales = case_cfg.get('properties', {}).get('scaling', {})
    U_ref = _to_float(scales.get('velocity_ref'), "properties.scaling.velocity_ref")
    if U_ref == 0.0:
        raise ValueError("properties.scaling.velocity_ref must be non-zero for prescribed_flow profile staging.")

    if any(bc.get("handler") == "prescribed_flow" for block in prepared_blocks for bc in block):
        profile_grid_dims = resolve_grid_block_dimensions_for_profiles(case_cfg, case_path, run_dir)

    generated_files = []
    generated_profile_summaries = []
    generated_target_grid = None
    field_slice_target_grid = None
    for i, block_bcs_list in enumerate(prepared_blocks):
        file_name = "bcs.run" if num_blocks == 1 else f"bcs_block{i}.run"
        bcs_file_path = os.path.join(config_dir, file_name)
        bcs_lines = [generate_header(run_id, source_files)]

        for bc in block_bcs_list:
            face, bc_type, handler = bc['face'], bc['type'], bc['handler']
            params = dict(bc.get('params') or {})
            if handler == "prescribed_flow":
                source = params.pop("source")
                expected_dims = _bc_profile_expected_dims(face, profile_grid_dims[i])
                staged_name = f"inlet_profile_block{i}_{face.replace('+', 'pos').replace('-', 'neg')}.picslice"
                staged_path = os.path.join(config_dir, staged_name)
                if source["type"] == "file":
                    source_path = source["path"]
                    if case_path and not os.path.isabs(source_path):
                        source_path = os.path.abspath(os.path.join(os.path.dirname(case_path), source_path))
                elif source["type"] == "generated":
                    default_output = os.path.join(
                        "config",
                        f"inlet_profile_block{i}_{_face_artifact_token(face)}.generated.picslice",
                    )
                    source_path = _resolve_run_artifact_path(
                        run_dir,
                        source.get("output_file"),
                        default_output,
                        default_to_config_dir=True,
                    )
                    if os.path.abspath(source_path) == os.path.abspath(staged_path):
                        raise ValueError(
                            f"Generated profile output_file for block {i}, face {face} must differ from staged solver profile."
                        )
                    if source["generator"] == "square_duct_poiseuille":
                        if generated_target_grid is None:
                            generated_target_grid = resolve_target_grid_for_generated_profile(case_cfg, case_path, run_dir)
                        summary = generate_square_duct_poiseuille_picslice(
                            source_path,
                            expected_dims,
                            source["params"],
                            target_grid=generated_target_grid,
                            target_block=i,
                            target_face=face,
                            script=source.get("script"),
                            case_path=case_path,
                        )
                    else:
                        raise ValueError(f"Unsupported generated profile generator '{source['generator']}'.")
                    summary.update({"block": i, "face": face})
                    generated_profile_summaries.append(summary)
                elif source["type"] == "field_slice":
                    default_output = os.path.join(
                        "config",
                        f"inlet_profile_block{i}_{_face_artifact_token(face)}.sliced.picslice",
                    )
                    source_path = _resolve_run_artifact_path(
                        run_dir,
                        source.get("output_file"),
                        default_output,
                        default_to_config_dir=True,
                    )
                    if os.path.abspath(source_path) == os.path.abspath(staged_path):
                        raise ValueError(
                            f"field_slice output_file for block {i}, face {face} must differ from staged solver profile."
                        )
                    if field_slice_target_grid is None:
                        field_slice_target_grid = resolve_target_grid_for_field_slice(case_cfg, case_path, run_dir)
                    summary = generate_field_slice_picslice(
                        source_path,
                        expected_dims,
                        source,
                        field_slice_target_grid,
                        face,
                        i,
                        case_path,
                    )
                    summary.update({"block": i, "face": face})
                    generated_profile_summaries.append(summary)
                else:
                    raise ValueError(f"Unsupported prescribed_flow source type '{source.get('type')}'.")
                summary = validate_and_nondimensionalize_picslice(source_path, staged_path, U_ref, expected_dims)
                print(
                    f"[SUCCESS] Staged prescribed_flow profile for block {i}, face {face}: "
                    f"{os.path.relpath(staged_path)} dims={summary['dims']}"
                )
                params["source_file"] = os.path.abspath(staged_path)
            params_str = ""
            if params:
                parts = []
                for k, v in params.items():
                    if isinstance(v, bool):
                        value_str = "true" if v else "false"
                    else:
                        value_str = str(v)
                    parts.append(f"{k}={value_str}")
                params_str = " ".join(parts)
            bcs_lines.append(f"{face:<20s} {bc_type:<12s} {handler:<20s} {params_str}")
        
        with open(bcs_file_path, "w") as f: f.write("\n".join(bcs_lines))
        
        print(f"[SUCCESS] Generated BCs for Block {i}: {os.path.relpath(bcs_file_path)}")
        generated_files.append(os.path.abspath(bcs_file_path))

    if generated_profile_summaries:
        info_path = write_profile_info(config_dir, generated_profile_summaries)
        print(f"[SUCCESS] Wrote generated profile summary: {os.path.relpath(info_path)}")
        
    return generated_files

def format_flag_value(value):
    """!
    @brief Converts Python types to C-style command-line flag values.
    @param[in] value The Python object to convert (bool, list, or other).
    @return A string representation suitable for a C command-line parser.
    """
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, list):
        return ",".join(map(str, value))
    return str(value)

def translate_programmatic_grid_settings(grid_settings: dict) -> dict:
    """!
    @brief Return programmatic-grid settings translated to the C node-count contract.
    @param[in] grid_settings Argument passed to `translate_programmatic_grid_settings()`.
    @return Value returned by `translate_programmatic_grid_settings()`.
    """
    translated = dict(grid_settings)
    for dim_key in ("im", "jm", "km"):
        if dim_key in translated:
            raw_val = translated[dim_key]
            if not isinstance(raw_val, int) or raw_val <= 0:
                raise ValueError(
                    f"grid.programmatic_settings.{dim_key} must be a positive integer cell count "
                    f"(got {raw_val!r})."
                )
            translated[dim_key] = raw_val + 1
    return translated


def generate_picgrid_from_programmatic_settings(raw_settings: dict, dest_path: str, L_ref: float) -> dict:
    """!
    @brief Generate a canonical PICGRID file from programmatic Cartesian grid settings.
    @details Implements the same coordinate formula as ComputeStretchedCoord in src/grid.c.
             im/jm/km in raw_settings are cell counts; node counts are im+1, jm+1, km+1.
    @param[in] raw_settings programmatic_settings dict from case.yml.
    @param[in] dest_path Destination PICGRID file path.
    @param[in] L_ref Reference length for nondimensionalization (must be non-zero).
    @return Summary dict: nblk, dims [(IM, JM, KM)], total_nodes.
    """
    if L_ref == 0.0:
        raise ValueError("length_ref must be non-zero for programmatic grid generation.")
    IM = int(raw_settings.get("im", 0)) + 1
    JM = int(raw_settings.get("jm", 0)) + 1
    KM = int(raw_settings.get("km", 0)) + 1
    if IM < 2 or JM < 2 or KM < 2:
        raise ValueError(
            f"programmatic_settings im/jm/km must each be >= 1 "
            f"(got im={IM-1}, jm={JM-1}, km={KM-1})."
        )
    x_min = float(raw_settings.get("xMins", 0.0))
    x_max = float(raw_settings.get("xMaxs", 1.0))
    y_min = float(raw_settings.get("yMins", 0.0))
    y_max = float(raw_settings.get("yMaxs", 1.0))
    z_min = float(raw_settings.get("zMins", 0.0))
    z_max = float(raw_settings.get("zMaxs", 1.0))
    rx = float(raw_settings.get("rxs", 1.0))
    ry = float(raw_settings.get("rys", 1.0))
    rz = float(raw_settings.get("rzs", 1.0))

    def _stretched(idx, N, length, r):
        """!
        @brief Mirror of ComputeStretchedCoord from src/grid.c.
        @param[in] idx Node index along the axis.
        @param[in] N Total node count along the axis.
        @param[in] length Physical length of the axis.
        @param[in] r Geometric stretching ratio.
        @return Coordinate offset from the axis minimum.
        """
        frac = idx / (N - 1.0)
        if abs(r - 1.0) < 1.0e-9:
            return length * frac
        return length * (r ** frac - 1.0) / (r - 1.0)

    Lx, Ly, Lz = x_max - x_min, y_max - y_min, z_max - z_min
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "w") as fout:
        fout.write("PICGRID\n1\n")
        fout.write(f"{IM} {JM} {KM}\n")
        for k in range(KM):
            z = (z_min + _stretched(k, KM, Lz, rz)) / L_ref
            for j in range(JM):
                y = (y_min + _stretched(j, JM, Ly, ry)) / L_ref
                for i in range(IM):
                    x = (x_min + _stretched(i, IM, Lx, rx)) / L_ref
                    fout.write(f"{x:.8e} {y:.8e} {z:.8e}\n")
    total_nodes = IM * JM * KM
    return {"nblk": 1, "dims": [(IM, JM, KM)], "total_nodes": total_nodes}


GRID_DA_PROCESSOR_KEYS = ("da_processors_x", "da_processors_y", "da_processors_z")


def resolve_grid_da_processor_layout(grid_cfg: dict) -> dict:
    """!
    @brief Resolve optional global DMDA layout, preferring grid-level keys over legacy nested keys.
    @param[in] grid_cfg Argument passed to `resolve_grid_da_processor_layout()`.
    @return Value returned by `resolve_grid_da_processor_layout()`.
    """
    top_level = {}
    legacy = {}

    for key in GRID_DA_PROCESSOR_KEYS:
        value = grid_cfg.get(key)
        if isinstance(value, (list, tuple)):
            raise ValueError(
                f"grid.{key} must be a scalar integer. "
                "Per-block MPI decomposition is not implemented on the C side; DMDA layout is global."
            )
        if value is not None:
            if not isinstance(value, int) or value <= 0:
                raise ValueError(f"grid.{key} must be a positive integer when provided (got {value}).")
            top_level[key] = value

    legacy_settings = grid_cfg.get("programmatic_settings")
    if isinstance(legacy_settings, dict):
        for key in GRID_DA_PROCESSOR_KEYS:
            value = legacy_settings.get(key)
            if isinstance(value, (list, tuple)):
                raise ValueError(
                    f"grid.programmatic_settings.{key} must be a scalar integer. "
                    "Per-block MPI decomposition is not implemented on the C side; DMDA layout is global."
                )
            if value is not None:
                if not isinstance(value, int) or value <= 0:
                    raise ValueError(
                        f"grid.programmatic_settings.{key} must be a positive integer when provided (got {value})."
                    )
                legacy[key] = value

    resolved = {}
    for key in GRID_DA_PROCESSOR_KEYS:
        top_value = top_level.get(key)
        legacy_value = legacy.get(key)
        if top_value is not None and legacy_value is not None and top_value != legacy_value:
            raise ValueError(
                f"grid.{key} conflicts with legacy grid.programmatic_settings.{key}; "
                "define the processor layout in only one place."
            )
        if top_value is not None:
            resolved[key] = top_value
        elif legacy_value is not None:
            resolved[key] = legacy_value

    return resolved


def append_grid_da_processor_layout(control_lines: list, grid_cfg: dict, num_procs: int) -> None:
    """!
    @brief Append optional global DMDA layout flags for any grid mode.
    @param[in] control_lines Argument passed to `append_grid_da_processor_layout()`.
    @param[in] grid_cfg Argument passed to `append_grid_da_processor_layout()`.
    @param[in] num_procs Argument passed to `append_grid_da_processor_layout()`.
    """
    layout = resolve_grid_da_processor_layout(grid_cfg)
    if not layout:
        if num_procs > 1:
            print("[INFO] Letting PETSc automatically determine processor layout.")
        return

    if num_procs == 1:
        print("[INFO] Serial run, ignoring da_processors layout.")
        return

    if all(layout.get(key) is not None for key in GRID_DA_PROCESSOR_KEYS):
        total_layout = 1
        for key in GRID_DA_PROCESSOR_KEYS:
            total_layout *= layout[key]
        if total_layout != num_procs:
            raise ValueError(f"Processor layout mismatch: product ({total_layout}) != processes ({num_procs}).")
        print(f"[INFO] Applying user-defined processor layout for {num_procs} processes.")
    else:
        printable = " x ".join(str(layout.get(key, "PETSC_DECIDE")) for key in GRID_DA_PROCESSOR_KEYS)
        print(f"[INFO] Applying partial processor layout: {printable}.")

    for key in GRID_DA_PROCESSOR_KEYS:
        value = layout.get(key)
        if value is not None:
            control_lines.append(f"-{key} {value}")

def normalize_momentum_solver_type(value: str) -> str:
    """!
    @brief Maps canonical user-facing momentum solver names to C-enum CLI values.
    @param[in] value Canonical momentum solver string from YAML.
    @return Canonical value accepted by -mom_solver_type.
    @throws ValueError if the input cannot be mapped.
    """
    # Only implemented YAML values belong here. Extend only with matching C enum/parser/dispatch support.
    if value is None:
        raise ValueError("momentum solver type cannot be None")

    raw = str(value).strip()
    mapped = {
        "Explicit RK4": "EXPLICIT_RK",
        "Dual Time Picard Jameson RK": "DUALTIME_PICARD_JAMESON_RK",
        "Dual Time Picard RK4": "DUALTIME_PICARD_JAMESON_RK",
    }.get(raw)
    if mapped is None:
        raise ValueError(
            f"Unknown momentum solver '{value}'. Use one of: "
            "'Explicit RK4', 'Dual Time Picard Jameson RK'."
        )
    return mapped

def normalize_solution_convergence_mode(value: str) -> str:
    """!
    @brief Normalizes the solution-convergence mode selector to the C-side canonical string.
    @param[in] value Human-readable solution-convergence mode selector.
    @return Canonical string accepted by `-solution_convergence_mode`.
    @throws ValueError if the input cannot be mapped.
    """
    if value is None:
        raise ValueError("solution_convergence.mode cannot be None")

    normalized = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "steady_deterministic": "STEADY_DETERMINISTIC",
        "periodic_deterministic": "PERIODIC_DETERMINISTIC",
        "statistical_steady": "STATISTICAL_STEADY",
        "transient": "TRANSIENT",
    }
    mapped = aliases.get(normalized)
    if mapped is None:
        raise ValueError(
            f"Unknown solution_convergence.mode '{value}'. Use one of: "
            "'steady_deterministic', 'periodic_deterministic', 'statistical_steady', 'transient'."
        )
    return mapped

def normalize_field_init_mode(value: str) -> int:
    """!
    @brief Maps canonical field init mode names to C enum/int codes (-finit).
    @param[in] value Canonical field initialization mode.
    @return Canonical integer code accepted by -finit.
    @throws ValueError if the input cannot be mapped.
    """
    # Canonical selector map for field initialization modes.
    if value is None:
        raise ValueError("field initialization mode cannot be None")

    mapped = {
        "Zero": 0,
        "Constant": 1,
        "Poiseuille": 2,
    }.get(str(value).strip())
    if mapped is None:
        raise ValueError(
            f"Unknown initial_conditions mode '{value}'. Use one of: 'Zero', 'Constant', 'Poiseuille'."
        )
    return mapped

def normalize_initial_condition_field(value: str) -> "tuple[str, int]":
    """!
    @brief Normalize a file IC field selector to its staged basename and C enum value.
    @param[in] value User-facing Ucat or Ucont selector.
    @return Tuple of staged field basename and C enum value.
    """
    normalized = str(value or "").strip().lower()
    if normalized == "ucat":
        return "ufield", 0
    if normalized == "ucont":
        return "vfield", 1
    raise ValueError("initial_conditions.field must be 'Ucat' or 'Ucont'.")

def resolve_initial_condition_config(ic: dict, prepared_blocks, U_ref: float) -> dict:
    """!
    @brief Resolve legacy and structured initial-condition YAML into one launcher contract.
    @param[in] ic Initial-condition YAML mapping.
    @param[in] prepared_blocks Normalized boundary-condition blocks.
    @param[in] U_ref Physical reference velocity.
    @return Normalized launcher initial-condition contract.
    """
    if not isinstance(ic, dict):
        raise ValueError("properties.initial_conditions must be a mapping.")
    mode = str(ic.get("mode", "")).strip()

    # Backward-compatible legacy spelling.
    if mode in {"Zero", "Constant", "Poiseuille"}:
        finit_code = normalize_field_init_mode(mode)
        params = resolve_ic_cli_params(ic, finit_code, prepared_blocks, U_ref)
        if finit_code == 1 and params.pop("ic_coordinate_system", 0) == 1:
            finit_code = 3
        return {"finit": finit_code, "cli_params": params, "kind": "builtin", "label": mode}

    normalized_mode = mode.lower().replace("-", "_").replace(" ", "_")
    if normalized_mode == "file":
        if prepared_blocks and len(prepared_blocks) > 1:
            raise ValueError("File-backed initial conditions currently support single-block cases only.")
        source_file = ic.get("source_file")
        if not isinstance(source_file, str) or not source_file.strip():
            raise ValueError("initial_conditions.source_file is required when mode is 'file'.")
        field_name, field_code = normalize_initial_condition_field(ic.get("field"))
        return {
            "finit": 4, "cli_params": {}, "kind": "file", "label": "file",
            "source_file": source_file.strip(), "field_name": field_name, "field_code": field_code,
        }
    if normalized_mode != "generated":
        raise ValueError("initial_conditions.mode must be 'generated' or 'file'.")

    generator = str(ic.get("generator", "")).strip().lower().replace("-", "_").replace(" ", "_")
    params = ic.get("params", {})
    if not isinstance(params, dict):
        raise ValueError("initial_conditions.params must be a mapping.")
    if generator == "ic_gen":
        if prepared_blocks and len(prepared_blocks) > 1:
            raise ValueError("File-backed initial conditions currently support single-block cases only.")
        script = params.get("script")
        if script is not None and (not isinstance(script, str) or not script.strip()):
            raise ValueError("initial_conditions.params.script must be a non-empty path when provided.")
        field_name, field_code = normalize_initial_condition_field(params.get("field"))
        config_file = params.get("config_file")
        if not isinstance(config_file, str) or not config_file.strip():
            raise ValueError("initial_conditions.params.config_file is required for generator 'ic_gen'.")
        cli_args = params.get("cli_args", [])
        if cli_args is None:
            cli_args = []
        if not isinstance(cli_args, list):
            raise ValueError("initial_conditions.params.cli_args must be a list.")
        return {
            "finit": 4, "cli_params": {}, "kind": "ic_gen", "label": "ic_gen",
            "field_name": field_name, "field_code": field_code,
            "config_file": config_file.strip(),
            "script": script.strip() if script is not None else None,
            "output_file": params.get("output_file"),
            "cli_args": cli_args,
        }

    generator_modes = {
        "zero": (0, "Zero"),
        "constant": (1, "Constant"),
        "streamwise_constant": (3, "Constant"),
        "poiseuille": (2, "Poiseuille"),
    }
    if generator not in generator_modes:
        raise ValueError(
            "initial_conditions.generator must be one of: zero, constant, "
            "streamwise_constant, poiseuille, ic_gen."
        )
    finit_code, legacy_mode = generator_modes[generator]
    legacy_ic = dict(params)
    legacy_ic["mode"] = legacy_mode
    cli_params = resolve_ic_cli_params(
        legacy_ic,
        1 if finit_code == 3 else finit_code,
        prepared_blocks,
        U_ref,
    )
    cli_params.pop("ic_coordinate_system", None)
    return {"finit": finit_code, "cli_params": cli_params, "kind": "builtin", "label": generator}

def validate_petsc_vec_binary(path: str) -> dict:
    """!
    @brief Validate the basic PETSc binary VecView envelope used by ReadFieldData.
    @param[in] path PETSc binary vector path.
    @return Summary containing the absolute path and scalar count.
    """
    import struct
    with open(path, "rb") as fin:
        header = fin.read(8)
        if len(header) != 8:
            raise ValueError(f"PETSc Vec file is too short: {path}")
        class_id, scalar_count = struct.unpack(">ii", header)
        if class_id != 1211214 or scalar_count < 0:
            raise ValueError(f"Invalid PETSc Vec header in {path}.")
        payload = fin.read()
    if len(payload) != scalar_count * 8:
        raise ValueError(
            f"PETSc Vec payload size mismatch in {path}: expected {scalar_count * 8} bytes, found {len(payload)}."
        )
    return {"path": os.path.abspath(path), "scalar_count": scalar_count}

def run_initial_condition_generator(case_path: str, run_dir: str, resolved_ic: dict) -> str:
    """!
    @brief Run the repository IC generator.
    @param[in] case_path Source case YAML path.
    @param[in] run_dir Run or precompute output directory.
    @param[in] resolved_ic Normalized external-generator contract.
    @return Generated PETSc vector path.
    """
    case_dir = os.path.dirname(os.path.abspath(case_path))
    script = _resolve_generator_script(resolved_ic.get("script"), case_path, "ic.gen")
    config_file = resolved_ic["config_file"]
    config_file = config_file if os.path.isabs(config_file) else os.path.abspath(os.path.join(case_dir, config_file))
    if not os.path.isfile(script):
        raise ValueError(f"ic.gen script not found: {script}")
    if not os.path.isfile(config_file):
        raise ValueError(f"initial-condition generator config file not found: {config_file}")
    default_output = os.path.join("config", "initial_condition.generated.dat")
    output_path = _resolve_run_artifact_path(
        run_dir, resolved_ic.get("output_file"), default_output, default_to_config_dir=True
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    cmd = [sys.executable, script, "-c", config_file, "--field",
           "Ucat" if resolved_ic["field_code"] == 0 else "Ucont",
           "--output", output_path]
    staged_grid = os.path.join(run_dir, "config", "grid.run")
    if os.path.isfile(staged_grid):
        cmd.extend(["--grid", staged_grid])
    cmd.extend(str(token) for token in resolved_ic.get("cli_args", []))
    result = subprocess.run(cmd, cwd=case_dir, text=True, capture_output=True)
    if result.returncode != 0:
        details = (result.stderr or result.stdout or "").strip()
        raise ValueError(f"ic.gen failed with exit code {result.returncode}. Details:\n{details}")
    validate_petsc_vec_binary(output_path)
    return output_path

def stage_initial_condition_file(run_dir: str, case_path: str, resolved_ic: dict) -> dict:
    """!
    @brief Materialize and stage one file-backed IC in ReadFieldData's expected layout.
    @param[in] run_dir Run or precompute output directory.
    @param[in] case_path Source case YAML path.
    @param[in] resolved_ic Normalized file-backed IC contract.
    @return Source, staged path, and staging-directory summary.
    """
    if resolved_ic["kind"] == "ic_gen":
        source_path = run_initial_condition_generator(case_path, run_dir, resolved_ic)
    else:
        source_path = resolved_ic["source_file"]
        if not os.path.isabs(source_path):
            source_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(case_path)), source_path))
        if not os.path.isfile(source_path):
            raise ValueError(f"Initial-condition source file not found: {source_path}")
        validate_petsc_vec_binary(source_path)
    stage_dir = os.path.join(run_dir, "config", "initial_condition")
    os.makedirs(stage_dir, exist_ok=True)
    staged_path = os.path.join(stage_dir, f"{resolved_ic['field_name']}00000_0.dat")
    if os.path.abspath(source_path) != os.path.abspath(staged_path):
        shutil.copy2(source_path, staged_path)
    return {"source": os.path.abspath(source_path), "staged": os.path.abspath(staged_path), "directory": os.path.abspath(stage_dir)}

def normalize_flow_direction_token(value: str) -> int:
    """!
    @brief Maps a face-token flow direction string to the C FlowDirection enum integer.
    @param[in] value One of '+Xi', '-Xi', '+Eta', '-Eta', '+Zeta', '-Zeta'.
    @return Integer 0-5 matching the FlowDirection enum.
    @throws ValueError on unknown value.
    """
    mapped = {
        "+Xi": 0, "-Xi": 1,
        "+Eta": 2, "-Eta": 3,
        "+Zeta": 4, "-Zeta": 5,
    }.get(str(value).strip())
    if mapped is None:
        raise ValueError(
            f"Unknown initial_conditions.flow_direction '{value}'. "
            "Use one of: '+Xi', '-Xi', '+Eta', '-Eta', '+Zeta', '-Zeta'."
        )
    return mapped

def _ic_has_inlet(prepared_blocks) -> bool:
    """!
    @brief Return True if any prepared BC block contains an INLET face.
    @param[in] prepared_blocks List of prepared BC lists (one per domain block).
    @return True if at least one INLET entry exists across all blocks.
    """
    if not prepared_blocks:
        return False
    for block_bcs in prepared_blocks:
        for entry in block_bcs:
            if entry.get("type") == "INLET":
                return True
    return False

def resolve_ic_cli_params(ic: dict, finit_code: int, prepared_blocks, U_ref: float) -> dict:
    """!
    @brief Resolve all IC parameters and return a dict of PETSc option values.
    @param[in] ic The properties.initial_conditions mapping.
    @param[in] finit_code Normalized -finit integer code.
    @param[in] prepared_blocks Normalized BC blocks (may be None).
    @param[in] U_ref Reference velocity for non-dimensionalization.
    @return Dict with keys matching PETSc option names (without leading dash).
    @throws KeyError if a required key is absent.
    @throws ValueError on invalid combinations or values.
    """
    result = {}

    if finit_code == 0:
        return result

    if finit_code == 1:  # Constant
        has_cartesian   = any(k in ic for k in ("u_physical", "v_physical", "w_physical"))
        has_curvilinear = "velocity_physical" in ic

        if has_cartesian and has_curvilinear:
            raise ValueError(
                "initial_conditions: cannot mix u/v/w_physical (cartesian) and "
                "velocity_physical (curvilinear) — use one or the other."
            )

        if has_curvilinear:  # curvilinear: scalar speed along flow_direction axis
            cs_code = 1
            result["ic_coordinate_system"] = cs_code
            try:
                vel_phys = float(ic["velocity_physical"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid value for initial_conditions.velocity_physical: {ic['velocity_physical']!r}. "
                    "Expected a numeric value."
                ) from exc
            result["ic_velocity_physical"] = vel_phys / U_ref if U_ref != 0 else 0.0

            if "flow_direction" in ic:
                result["flow_direction"] = normalize_flow_direction_token(ic["flow_direction"])
            elif not _ic_has_inlet(prepared_blocks):
                raise ValueError(
                    "initial_conditions.flow_direction is required for curvilinear Constant IC "
                    "when no INLET face exists."
                )

        else:  # cartesian: u/v/w_physical → Cart2Contra (default when no velocity_physical)
            if "flow_direction" in ic:
                raise ValueError(
                    "initial_conditions.flow_direction is not valid for cartesian Constant IC. "
                    "Use velocity_physical + flow_direction for curvilinear mode."
                )
            cs_code = 0
            result["ic_coordinate_system"] = cs_code
            u, v, w = parse_initial_velocity_components(ic, finit_code, require_explicit=True)
            scale = 1.0 / U_ref if U_ref != 0 else 0.0
            result["ucont_x"] = u * scale
            result["ucont_y"] = v * scale
            result["ucont_z"] = w * scale

    elif finit_code == 2:  # Poiseuille
        if any(k in ic for k in ("u_physical", "v_physical", "w_physical")):
            raise ValueError(
                "For Poiseuille mode, use peak_velocity_physical, not u_physical/v_physical/w_physical."
            )
        if "velocity_physical" in ic:
            raise ValueError(
                "For Poiseuille mode, use peak_velocity_physical, not velocity_physical."
            )
        if "peak_velocity_physical" not in ic:
            raise KeyError("peak_velocity_physical")
        try:
            peak = float(ic["peak_velocity_physical"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid value for initial_conditions.peak_velocity_physical: "
                f"{ic['peak_velocity_physical']!r}. Expected a numeric value."
            ) from exc
        result["ic_velocity_physical"] = peak / U_ref if U_ref != 0 else 0.0

        if "flow_direction" in ic:
            fd_int = normalize_flow_direction_token(ic["flow_direction"])
            # Cross-check: explicit flow_direction must align with INLET face if one exists
            if _ic_has_inlet(prepared_blocks):
                inlet_axis = infer_unique_inlet_axis_from_prepared_bcs(prepared_blocks)
                fd_axis_name = {0: "x", 1: "y", 2: "z"}.get(fd_int // 2, "?")
                if inlet_axis and fd_axis_name != inlet_axis:
                    token = ic["flow_direction"]
                    raise ValueError(
                        f"initial_conditions.flow_direction '{token}' (axis '{fd_axis_name}') "
                        f"does not match INLET face axis '{inlet_axis}'."
                    )
            result["flow_direction"] = fd_int
        elif not _ic_has_inlet(prepared_blocks):
            raise ValueError(
                "initial_conditions.flow_direction is required for Poiseuille IC "
                "when no INLET face exists."
            )

    return result

def normalize_eulerian_field_source(value: str) -> str:
    """!
    @brief Normalizes the Eulerian field source selector to the C-side canonical string.
    @param[in] value Human-readable or enum-like Eulerian field source.
    @return Canonical string accepted by `-euler_field_source`.
    @throws ValueError if the input cannot be mapped.
    """
    if value is None:
        raise ValueError("eulerian_field_source cannot be None")

    normalized = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "solve": "solve",
        "load": "load",
        "analytical": "analytical",
    }
    mapped = aliases.get(normalized)
    if mapped is None:
        raise ValueError(
            f"Unknown operation_mode.eulerian_field_source '{value}'. "
            "Use one of: 'solve', 'load', 'analytical'."
        )
    return mapped

def normalize_analytical_type(value: str) -> str:
    """!
    @brief Normalizes the analytical solution selector to the C-side canonical string.
    @param[in] value Human-readable analytical solution selector.
    @return Canonical string accepted by `-analytical_type`.
    @throws ValueError if the input cannot be mapped.
    """
    # Canonical selector map for analytical solution selectors.
    if value is None:
        raise ValueError("analytical_type cannot be None")

    normalized = str(value).strip().upper().replace("-", "_").replace(" ", "_")
    if normalized not in {"TGV3D", "ZERO_FLOW", "UNIFORM_FLOW"}:
        raise ValueError(
            f"Unknown operation_mode.analytical_type '{value}'. "
            "Use one of: 'TGV3D', 'ZERO_FLOW', 'UNIFORM_FLOW'."
        )
    return normalized

def parse_initial_velocity_components(initial_conditions: dict, finit_code: int, *, require_explicit: bool) -> "tuple[float, float, float]":
    """!
    @brief Parse initial-condition velocity components with mode-aware defaults.
    @param[in] initial_conditions The `properties.initial_conditions` mapping from case.yml.
    @param[in] finit_code Normalized `-finit` integer code.
    @param[in] require_explicit If True, all three component keys must be present.
    @return Tuple `(u, v, w)` in physical units.
    @throws KeyError if a required component key is missing.
    @throws ValueError if a component cannot be converted to float.
    """
    component_keys = ("u_physical", "v_physical", "w_physical")
    components = []
    for key in component_keys:
        if key not in initial_conditions:
            if require_explicit:
                raise KeyError(key)
            raw_value = 0.0
        else:
            raw_value = initial_conditions[key]
        try:
            components.append(float(raw_value))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid value for properties.initial_conditions.{key}: {raw_value!r}. Expected a numeric value."
            ) from exc
    return tuple(components)

def infer_unique_inlet_axis_from_prepared_bcs(prepared_blocks: list) -> "str | None":
    """!
    @brief Infer the unique inlet axis across all blocks using C-side "primary inlet" ordering.
    @param[in] prepared_blocks Normalized BC blocks from `validate_and_prepare_boundary_conditions`.
    @return One of `"x"`, `"y"`, `"z"` if unique, `None` if no inlet exists.
    @throws ValueError if different blocks imply different inlet axes.
    """
    face_order = ("-Xi", "+Xi", "-Eta", "+Eta", "-Zeta", "+Zeta")
    face_axis = {
        "-Xi": "x", "+Xi": "x",
        "-Eta": "y", "+Eta": "y",
        "-Zeta": "z", "+Zeta": "z",
    }

    inlet_axes = set()
    for block_bcs in prepared_blocks:
        face_map = {entry["face"]: entry for entry in block_bcs}
        for face in face_order:
            entry = face_map.get(face)
            if entry and entry["type"] == "INLET":
                inlet_axes.add(face_axis[face])
                break

    if not inlet_axes:
        return None
    if len(inlet_axes) != 1:
        raise ValueError(
            "properties.initial_conditions.peak_velocity_physical requires all blocks to have a primary INLET "
            f"on the same axis. Found axes: {sorted(inlet_axes)}. Use u_physical/v_physical/w_physical instead."
        )
    return next(iter(inlet_axes))

def normalize_particle_init_mode(value: str) -> int:
    """!
    @brief Maps canonical particle init mode names to C enum/int codes (-pinit).
    @param[in] value Canonical particle initialization mode.
    @return Canonical integer code accepted by -pinit.
    @throws ValueError if the input cannot be mapped.
    """
    # Canonical selector map for particle initialization modes.
    if value is None:
        raise ValueError("particle init mode cannot be None")

    mapped = {
        "Surface": 0,
        "Volume": 1,
        "PointSource": 2,
        "SurfaceEdges": 3,
    }.get(str(value).strip())
    if mapped is None:
        raise ValueError(
            f"Unknown particle init_mode '{value}'. Use one of: "
            "'Surface', 'Volume', 'PointSource', 'SurfaceEdges'."
        )
    return mapped

def normalize_interpolation_method(value: str) -> int:
    """!
    @brief Maps interpolation method names to C enum/int codes (-interpolation_method).
    @param[in] value Canonical interpolation method name.
    @return Integer code accepted by -interpolation_method.
    @throws ValueError if the input cannot be mapped.
    """
    if value is None:
        raise ValueError("interpolation method cannot be None")

    mapped = {
        "Trilinear": 0,
        "CornerAveraged": 1,
    }.get(str(value).strip())
    if mapped is None:
        raise ValueError(
            f"Unknown interpolation_method '{value}'. Use one of: "
            "'Trilinear', 'CornerAveraged'."
        )
    return mapped

def normalize_les_model(value) -> int:
    """!
    @brief Maps LES model selectors to C enum/int codes (-les).
    @param[in] value LES selector name or legacy integer/bool value.
    @return Integer code accepted by -les.
    @throws ValueError if the input cannot be mapped.
    """
    if isinstance(value, bool):
        return 1 if value else 0
    if isinstance(value, int):
        if value in (0, 1, 2):
            return value
        raise ValueError("models.physics.turbulence.les must be 0, 1, 2, false/true, or a supported model block.")
    if value is None:
        raise ValueError("LES model cannot be None")

    key = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    mapped = {
        "none": 0,
        "off": 0,
        "disabled": 0,
        "no_les": 0,
        "constant": 1,
        "constant_smagorinsky": 1,
        "smagorinsky": 1,
        "dynamic": 2,
        "dynamic_smagorinsky": 2,
    }.get(key)
    if mapped is None:
        raise ValueError(
            f"Unknown LES model '{value}'. Use one of: 'none', "
            "'constant_smagorinsky', 'dynamic_smagorinsky'."
        )
    return mapped

def normalize_les_test_filter(value) -> int:
    """!
    @brief Maps LES test-filter names to the C -testfilter_ik flag.
    @param[in] value Test-filter selector name or legacy integer/bool value.
    @return 0 for volume-weighted box, 1 for homogeneous i/k Simpson filtering.
    @throws ValueError if the input cannot be mapped.
    """
    if isinstance(value, bool):
        return 1 if value else 0
    if isinstance(value, int):
        if value in (0, 1):
            return value
        raise ValueError("models.physics.turbulence.les.test_filter must be 0, 1, or a supported filter name.")
    if value is None:
        raise ValueError("LES test_filter cannot be None")

    key = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    mapped = {
        "volume_weighted_box": 0,
        "box": 0,
        "general_box": 0,
        "homogeneous_ik": 1,
        "ik_homogeneous": 1,
        "simpson_ik": 1,
    }.get(key)
    if mapped is None:
        raise ValueError(
            f"Unknown LES test_filter '{value}'. Use one of: "
            "'volume_weighted_box', 'homogeneous_ik'."
        )
    return mapped

def normalize_rans_model(value) -> int:
    """!
    @brief Maps RANS model selectors to the current C -rans switch.
    @param[in] value RANS selector name or legacy integer/bool value.
    @return Integer code accepted by -rans.
    @throws ValueError if the input cannot be mapped.
    """
    if isinstance(value, bool):
        return 1 if value else 0
    if isinstance(value, int):
        if value in (0, 1):
            return value
        raise ValueError("models.physics.turbulence.rans must be 0, 1, false/true, or a supported model block.")
    if value is None:
        raise ValueError("RANS model cannot be None")

    key = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    mapped = {
        "none": 0,
        "off": 0,
        "disabled": 0,
        "k_omega": 1,
        "komega": 1,
    }.get(key)
    if mapped is None:
        raise ValueError(f"Unknown RANS model '{value}'. Use one of: 'none', 'k_omega'.")
    return mapped

def normalize_wall_function_model(value) -> str:
    """!
    @brief Validates wall-function model selectors exposed in YAML.
    @param[in] value Wall-function selector name.
    @return Canonical wall-function model name.
    @throws ValueError if the input cannot be mapped.
    """
    if value is None:
        return "log_law"
    key = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if key in {"log_law", "loglaw"}:
        return "log_law"
    raise ValueError("Unknown wall_function model '%s'. Use: 'log_law'." % value)

def resolve_enabled_flag(cfg: dict, path: str, default: bool = True) -> bool:
    """!
    @brief Resolves a structured `enabled` flag and rejects non-boolean values.
    @param[in] cfg Mapping that may contain `enabled`.
    @param[in] path Human-readable config path for diagnostics.
    @param[in] default Value used when `enabled` is omitted.
    @return Boolean enabled state.
    @throws ValueError if `enabled` is not a YAML boolean.
    """
    if 'enabled' not in cfg:
        return default
    if not isinstance(cfg['enabled'], bool):
        raise ValueError(f"{path}.enabled must be true or false.")
    return cfg['enabled']

def append_turbulence_flags(models: dict, control_lines: list):
    """!
    @brief Appends turbulence model flags from legacy or structured case.yml blocks.
    @param[in] models Parsed case.yml `models` mapping.
    @param[out] control_lines A list of strings to which C-flags will be appended.
    """
    turbulence_cfg = models.get('physics', {}).get('turbulence', {})
    if not turbulence_cfg:
        return
    if not isinstance(turbulence_cfg, dict):
        raise ValueError("models.physics.turbulence must be a mapping.")

    les_cfg = turbulence_cfg.get('les')
    rans_cfg = turbulence_cfg.get('rans')
    wall_cfg = turbulence_cfg.get('wall_function')
    les_code = None
    rans_code = None

    if isinstance(les_cfg, dict):
        enabled = resolve_enabled_flag(les_cfg, "models.physics.turbulence.les")
        model_value = les_cfg.get('model', 'constant_smagorinsky')
        les_code = normalize_les_model(model_value) if enabled else 0
        control_lines.append(f"-les {les_code}")
        if 'constant_cs' in les_cfg:
            control_lines.append(f"-const_cs {format_flag_value(les_cfg['constant_cs'])}")
        if 'max_cs' in les_cfg:
            control_lines.append(f"-max_cs {format_flag_value(les_cfg['max_cs'])}")
        if 'dynamic_frequency' in les_cfg:
            control_lines.append(f"-dynamic_freq {format_flag_value(les_cfg['dynamic_frequency'])}")
        if 'test_filter' in les_cfg:
            control_lines.append(f"-testfilter_ik {normalize_les_test_filter(les_cfg['test_filter'])}")
    elif les_cfg is not None:
        les_code = normalize_les_model(les_cfg)
        control_lines.append(f"-les {les_code}")

    if isinstance(rans_cfg, dict):
        enabled = resolve_enabled_flag(rans_cfg, "models.physics.turbulence.rans")
        model_value = rans_cfg.get('model', 'k_omega')
        rans_code = normalize_rans_model(model_value) if enabled else 0
        control_lines.append(f"-rans {rans_code}")
    elif rans_cfg is not None:
        rans_code = normalize_rans_model(rans_cfg)
        control_lines.append(f"-rans {rans_code}")

    if les_code and rans_code:
        raise ValueError("models.physics.turbulence cannot enable both LES and RANS in the same case.")

    if isinstance(wall_cfg, dict):
        enabled = resolve_enabled_flag(wall_cfg, "models.physics.turbulence.wall_function")
        normalize_wall_function_model(wall_cfg.get('model'))
        control_lines.append(f"-wallfunction {1 if enabled else 0}")
        if 'roughness_height' in wall_cfg:
            control_lines.append(f"-wall_roughness {format_flag_value(wall_cfg['roughness_height'])}")
    elif wall_cfg is not None:
        control_lines.append(f"-wallfunction {format_flag_value(wall_cfg)}")

def append_passthrough_flags(control_lines: list, options: dict):
    """!
    @brief Appends raw CLI flags to the control list from a {flag: value} dict.
    @details Boolean `true` is emitted as a switch with no value. Boolean `false`
             is skipped. All other values are emitted as "<flag> <value>".
    @param[out] control_lines The destination list of control-file lines.
    @param[in] options Mapping of raw CLI flags to values.
    """
    if not options:
        return
    for flag, value in options.items():
        if isinstance(value, bool):
            if value:
                control_lines.append(str(flag))
            continue
        control_lines.append(f"{flag} {format_flag_value(value)}")


SOLVER_MONITORING_POISSON_FLAG_MAP = {
    "pic_true_residual": "-ps_ksp_pic_monitor_true_residual",
    "true_residual": "-ps_ksp_monitor_true_residual",
    "converged_reason": "-ps_ksp_converged_reason",
    "view": "-ps_ksp_view",
}


def resolve_solver_monitoring_flags(monitor_cfg: dict) -> dict:
    """!
    @brief Resolve human-readable solver monitoring YAML to raw control flags.
    @param[in] monitor_cfg Parsed monitor.yml mapping.
    @return Mapping of raw C/PETSc flags to values.
    """
    solver_mon_cfg = monitor_cfg.get("solver_monitoring", {}) if isinstance(monitor_cfg, dict) else {}
    if solver_mon_cfg is None:
        return {}
    if not isinstance(solver_mon_cfg, dict):
        raise ValueError("monitor.solver_monitoring must be a mapping when provided.")

    flags = {}

    poisson_cfg = solver_mon_cfg.get("poisson", {}) or {}
    if not isinstance(poisson_cfg, dict):
        raise ValueError("monitor.solver_monitoring.poisson must be a mapping when provided.")
    unknown_poisson = sorted(set(poisson_cfg.keys()) - set(SOLVER_MONITORING_POISSON_FLAG_MAP.keys()))
    if unknown_poisson:
        raise ValueError(f"monitor.solver_monitoring.poisson has unsupported key(s): {unknown_poisson}.")
    for key, flag in SOLVER_MONITORING_POISSON_FLAG_MAP.items():
        if key in poisson_cfg:
            value = poisson_cfg[key]
            if not isinstance(value, bool):
                raise ValueError(f"monitor.solver_monitoring.poisson.{key} must be boolean.")
            flags[flag] = value

    passthrough = solver_mon_cfg.get("petsc_passthrough_options", {}) or {}
    if not isinstance(passthrough, dict):
        raise ValueError("monitor.solver_monitoring.petsc_passthrough_options must be a mapping when provided.")
    flags.update(passthrough)

    legacy_raw = {
        key: value
        for key, value in solver_mon_cfg.items()
        if isinstance(key, str) and key.startswith("-")
    }
    flags.update(legacy_raw)

    unknown_top = sorted(
        key
        for key in solver_mon_cfg.keys()
        if key not in {"poisson", "petsc_passthrough_options"} and not (isinstance(key, str) and key.startswith("-"))
    )
    if unknown_top:
        raise ValueError(
            "monitor.solver_monitoring has unsupported key(s): "
            f"{unknown_top}. Use 'poisson' for structured monitors or "
            "'petsc_passthrough_options' for raw PETSc flags."
        )

    return flags


def resolve_particle_console_output_frequency(io_cfg: dict) -> "int | None":
    """!
    @brief Return the effective particle-console snapshot cadence from monitor.yml.
    @param[in] io_cfg Argument passed to `resolve_particle_console_output_frequency()`.
    @return Value returned by `resolve_particle_console_output_frequency()`.
    """
    if 'particle_console_output_frequency' in io_cfg:
        return io_cfg['particle_console_output_frequency']
    return io_cfg.get('data_output_frequency')

def parse_and_add_model_flags(case_cfg: dict, control_lines: list):
    """!
    @brief Parses the 'models' section of case.yml and adds corresponding C-solver flags.
    @param[in] case_cfg The parsed case.yml configuration dictionary.
    @param[out] control_lines A list of strings to which C-flags will be appended.
    """
    models = case_cfg.get('models', {})
    FLAG_MAP = {
        'domain': {'blocks': '-nblk'},
        'physics.fsi': {'immersed': '-imm', 'moving_fsi': '-fsi'},
        'physics.particles': {'count': '-numParticles'},
        'statistics': {'time_averaging': '-averaging'}
    }
    for section_path, flags in FLAG_MAP.items():
        current_level = models
        try:
            for key in section_path.split('.'): current_level = current_level[key]
            for yaml_key, flag in flags.items():
                if yaml_key in current_level:
                    control_lines.append(f"{flag} {format_flag_value(current_level[yaml_key])}")
        except KeyError: continue

    append_turbulence_flags(models, control_lines)

    if models.get('physics', {}).get('dimensionality') == '2D':
        control_lines.append("-TwoD 1")
    
    particles_cfg = models.get('physics', {}).get('particles', {})
    p_init_mode_str = particles_cfg.get('init_mode', 'Surface')
    pinit_code = normalize_particle_init_mode(p_init_mode_str)
    control_lines.append(f"-pinit {pinit_code}")
    print(f"  - Particle Initialization Mode: {p_init_mode_str} (Code: {pinit_code})")

    if pinit_code == 2:
        point_cfg = particles_cfg.get('point_source', {})
        if not isinstance(point_cfg, dict):
            raise ValueError("models.physics.particles.point_source must be a mapping when init_mode is PointSource.")
        try:
            psrc_x = float(point_cfg['x'])
            psrc_y = float(point_cfg['y'])
            psrc_z = float(point_cfg['z'])
        except (KeyError, TypeError, ValueError):
            raise ValueError("PointSource init_mode requires numeric point_source.{x,y,z} values.")
        control_lines.append(f"-psrc_x {psrc_x}")
        control_lines.append(f"-psrc_y {psrc_y}")
        control_lines.append(f"-psrc_z {psrc_z}")
        print(f"  - Particle Point Source: ({psrc_x}, {psrc_y}, {psrc_z})")

    p_restart_mode = particles_cfg.get('restart_mode')
    if p_restart_mode:
        p_restart_mode_normalized = str(p_restart_mode).lower()
        if p_restart_mode_normalized not in {"init", "load"}:
            raise ValueError(f"Unknown particle restart_mode '{p_restart_mode}'. Options are 'init' or 'load'.")
        control_lines.append(f"-particle_restart_mode \"{p_restart_mode}\"")

def parse_solver_config(solver_cfg: dict) -> dict:
    """!
    @brief Parses the structured solver.yml into a flat dictionary of {flag: value}.
    @param[in] solver_cfg The parsed solver.yml configuration dictionary.
    @return A dictionary where keys are C-solver flags and values are the corresponding settings.
    """
    flags = {}
    if 'operation_mode' in solver_cfg and isinstance(solver_cfg['operation_mode'], dict):
        op_mode = solver_cfg['operation_mode']
        if 'eulerian_field_source' in op_mode:
            normalized_source = normalize_eulerian_field_source(op_mode.get('eulerian_field_source'))
            flags['-euler_field_source'] = f"\"{normalized_source}\""
        if 'analytical_type' in op_mode and op_mode.get('analytical_type') is not None:
            normalized_analytical_type = normalize_analytical_type(op_mode.get('analytical_type'))
            flags['-analytical_type'] = f"\"{normalized_analytical_type}\""
            if normalized_analytical_type == "UNIFORM_FLOW":
                uniform_flow_cfg = op_mode.get('uniform_flow', {})
                if not isinstance(uniform_flow_cfg, dict):
                    raise ValueError("operation_mode.uniform_flow must be a mapping when analytical_type is 'UNIFORM_FLOW'.")
                try:
                    flags['-analytical_uniform_u'] = float(uniform_flow_cfg['u'])
                    flags['-analytical_uniform_v'] = float(uniform_flow_cfg['v'])
                    flags['-analytical_uniform_w'] = float(uniform_flow_cfg['w'])
                except KeyError as exc:
                    raise ValueError(f"operation_mode.uniform_flow.{exc.args[0]} is required when analytical_type is 'UNIFORM_FLOW'.") from exc
                except (TypeError, ValueError) as exc:
                    raise ValueError("operation_mode.uniform_flow.{u,v,w} must be numeric when analytical_type is 'UNIFORM_FLOW'.") from exc

    verification_cfg = solver_cfg.get('verification', {})
    if verification_cfg:
        if not isinstance(verification_cfg, dict):
            raise ValueError("verification must be a mapping when provided.")
        sources_cfg = verification_cfg.get('sources', {})
        if not isinstance(sources_cfg, dict):
            raise ValueError("verification.sources must be a mapping when provided.")
        diff_cfg = sources_cfg.get('diffusivity')
        if diff_cfg is not None:
            if not isinstance(diff_cfg, dict):
                raise ValueError("verification.sources.diffusivity must be a mapping.")
            try:
                flags['-verification_diffusivity_mode'] = f"\"{str(diff_cfg['mode']).strip().lower()}\""
                flags['-verification_diffusivity_profile'] = f"\"{str(diff_cfg['profile']).strip().upper()}\""
                flags['-verification_diffusivity_gamma0'] = float(diff_cfg['gamma0'])
                flags['-verification_diffusivity_slope_x'] = float(diff_cfg['slope_x'])
            except KeyError as exc:
                raise ValueError(f"verification.sources.diffusivity.{exc.args[0]} is required.") from exc
            except (TypeError, ValueError) as exc:
                raise ValueError("verification.sources.diffusivity.{gamma0,slope_x} must be numeric and mode/profile must be scalar strings.") from exc

        scalar_cfg = sources_cfg.get('scalar')
        if scalar_cfg is not None:
            if not isinstance(scalar_cfg, dict):
                raise ValueError("verification.sources.scalar must be a mapping.")
            try:
                flags['-verification_scalar_mode'] = f"\"{str(scalar_cfg['mode']).strip().lower()}\""
                flags['-verification_scalar_profile'] = f"\"{str(scalar_cfg['profile']).strip().upper()}\""
            except KeyError as exc:
                raise ValueError(f"verification.sources.scalar.{exc.args[0]} is required.") from exc

            scalar_numeric_keys = {
                'CONSTANT': ('value',),
                'LINEAR_X': ('phi0', 'slope_x'),
                'SIN_PRODUCT': ('amplitude', 'kx', 'ky', 'kz'),
            }
            profile = str(scalar_cfg.get('profile', '')).strip().upper()
            for key in scalar_numeric_keys.get(profile, ()):
                try:
                    flags[f'-verification_scalar_{key}'] = float(scalar_cfg[key])
                except KeyError as exc:
                    raise ValueError(f"verification.sources.scalar.{exc.args[0]} is required.") from exc
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"verification.sources.scalar.{key} must be numeric.") from exc

    transport_cfg = solver_cfg.get('scalar_transport', {})
    if transport_cfg:
        if not isinstance(transport_cfg, dict):
            raise ValueError("scalar_transport must be a mapping when provided.")
        transport_map = {
            'schmidt_number': '-schmidt_number',
            'turbulent_schmidt_number': '-turb_schmidt_number',
        }
        unknown_transport_keys = sorted(set(transport_cfg.keys()) - set(transport_map.keys()))
        if unknown_transport_keys:
            raise ValueError(
                f"scalar_transport has unsupported key(s): {unknown_transport_keys}. "
                "Use 'schmidt_number' or 'turbulent_schmidt_number'."
            )
        for key, flag in transport_map.items():
            if key in transport_cfg:
                try:
                    value = float(transport_cfg[key])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"scalar_transport.{key} must be numeric.") from exc
                if value <= 0.0:
                    raise ValueError(f"scalar_transport.{key} must be positive.")
                flags[flag] = value

    selected_solver = None
    if 'strategy' in solver_cfg:
        s = solver_cfg['strategy']
        if 'central_diff' in s:
            flags['-central'] = format_flag_value(s['central_diff'])
        # Preferred selector.
        if 'momentum_solver' in s:
            selected_solver = normalize_momentum_solver_type(s['momentum_solver'])
        elif 'implicit' in s:
            raise ValueError("Legacy key 'strategy.implicit' is not supported. Use 'strategy.momentum_solver'.")

    ms = solver_cfg.get('momentum_solver', {})
    if selected_solver is None:
        selected_solver = "DUALTIME_PICARD_JAMESON_RK"
    flags['-mom_solver_type'] = f"\"{selected_solver}\""

    if 'tolerances' in solver_cfg:
        t = solver_cfg['tolerances']
        tol_map = {
            'max_iterations': '-mom_max_pseudo_steps',
            'absolute_tol': '-mom_atol',
            'relative_tol': '-mom_rtol',
            'residual_absolute_tol': '-mom_resid_atol',
            'residual_relative_tol': '-mom_resid_rtol',
            'step_tol': '-imp_stol'
        }
        for key, flag in tol_map.items():
            if key in t:
                flags[flag] = t[key]

    def _append_dualtime_options(cfg: dict):
        """!
        @brief Append dualtime options.
        @param[in] cfg Argument passed to `_append_dualtime_options()`.
        """
        if 'max_pseudo_steps' in cfg:
            flags['-mom_max_pseudo_steps'] = cfg['max_pseudo_steps']
        if 'absolute_tol' in cfg:
            flags['-mom_atol'] = cfg['absolute_tol']
        if 'relative_tol' in cfg:
            flags['-mom_rtol'] = cfg['relative_tol']
        if 'step_tol' in cfg:
            flags['-imp_stol'] = cfg['step_tol']
        if 'pseudo_cfl' in cfg:
            pcfl = cfg['pseudo_cfl']
            if 'initial' in pcfl:
                flags['-pseudo_cfl'] = pcfl['initial']
            if 'minimum' in pcfl:
                flags['-min_pseudo_cfl'] = pcfl['minimum']
            if 'maximum' in pcfl:
                flags['-max_pseudo_cfl'] = pcfl['maximum']
            if 'growth_factor' in pcfl:
                flags['-pseudo_cfl_growth_factor'] = pcfl['growth_factor']
            if 'reduction_factor' in pcfl:
                flags['-pseudo_cfl_reduction_factor'] = pcfl['reduction_factor']
        if 'jameson_residual_noise_allowance_factor' in cfg:
            flags['-mom_dt_jameson_residual_norm_noise_allowance_factor'] = cfg['jameson_residual_noise_allowance_factor']
        elif 'rk4_residual_noise_allowance_factor' in cfg:
            flags['-mom_dt_jameson_residual_norm_noise_allowance_factor'] = cfg['rk4_residual_noise_allowance_factor']

    if isinstance(ms, dict):
        allowed_ms_keys = {'type', 'dual_time_picard_jameson_rk', 'dual_time_picard_rk4'}
        unknown_ms_keys = sorted(set(ms.keys()) - allowed_ms_keys)
        if unknown_ms_keys:
            raise ValueError(
                f"Unsupported momentum_solver keys/blocks: {unknown_ms_keys}. "
                "Currently supported block: 'dual_time_picard_jameson_rk'."
            )

        if 'dual_time_picard_jameson_rk' in ms and 'dual_time_picard_rk4' in ms:
            raise ValueError(
                "Use only momentum_solver.dual_time_picard_jameson_rk; "
                "do not also set its deprecated dual_time_picard_rk4 alias."
            )
        dt_picard_cfg = ms.get('dual_time_picard_jameson_rk', ms.get('dual_time_picard_rk4'))
        if dt_picard_cfg is not None:
            if selected_solver != "DUALTIME_PICARD_JAMESON_RK":
                raise ValueError(
                    f"momentum_solver.dual_time_picard_jameson_rk is set but selected solver is {selected_solver}."
                )
            if not isinstance(dt_picard_cfg, dict):
                raise ValueError("momentum_solver.dual_time_picard_jameson_rk must be a mapping.")
            if ('jameson_residual_noise_allowance_factor' in dt_picard_cfg and
                    'rk4_residual_noise_allowance_factor' in dt_picard_cfg):
                raise ValueError(
                    "Use only jameson_residual_noise_allowance_factor; "
                    "do not also set its deprecated rk4_residual_noise_allowance_factor alias."
                )
            _append_dualtime_options(dt_picard_cfg)
    solution_convergence_cfg = solver_cfg.get('solution_convergence', {})
    if solution_convergence_cfg is not None:
        if not isinstance(solution_convergence_cfg, dict):
            raise ValueError("solution_convergence must be a mapping when provided.")
        if solution_convergence_cfg and solution_convergence_cfg.get('enabled', True):
            mode = normalize_solution_convergence_mode(solution_convergence_cfg.get('mode', 'steady_deterministic'))
            flags['-solution_convergence_mode'] = f"\"{mode}\""
            if mode == "PERIODIC_DETERMINISTIC":
                periodic_cfg = solution_convergence_cfg.get('periodic_deterministic')
                if not isinstance(periodic_cfg, dict):
                    raise ValueError("solution_convergence.periodic_deterministic must be a mapping when mode is 'periodic_deterministic'.")
                flags['-solution_convergence_period_steps'] = periodic_cfg['period_steps']
            if mode == "STATISTICAL_STEADY":
                statistical_cfg = solution_convergence_cfg.get('statistical_steady')
                if not isinstance(statistical_cfg, dict):
                    raise ValueError("solution_convergence.statistical_steady must be a mapping when mode is 'statistical_steady'.")
                flags['-solution_convergence_window_steps'] = statistical_cfg['window_steps']
    def _normalize_poisson_method(value) -> str:
        """!
        @brief Normalize a user-facing Poisson linear-solver method name.
        @param[in] value Method value from the solver YAML.
        @return Lowercase PETSc KSP method token.
        """
        method = str(value).strip().lower()
        if not method:
            raise ValueError("poisson_solver.method cannot be empty.")
        return method

    def _normalize_poisson_preconditioner(value) -> str:
        """!
        @brief Normalize and validate the outer Poisson preconditioner name.
        @param[in] value Preconditioner value from the solver YAML.
        @return PETSc PC token for the supported outer preconditioner.
        """
        pc = str(value).strip().lower()
        aliases = {"mg": "multigrid", "pcmg": "multigrid"}
        pc = aliases.get(pc, pc)
        if pc != "multigrid":
            raise ValueError(
                "poisson_solver.preconditioner.type currently supports only 'multigrid'. "
                "The runtime Poisson solver still assumes PETSc PCMG setup."
            )
        return "mg"

    def _poisson_level_number(level_name) -> str:
        """!
        @brief Extract the numeric suffix from a `level_N` multigrid level key.
        @param[in] level_name YAML level key supplied by the user.
        @return Numeric level suffix as a string.
        """
        text = str(level_name).strip()
        match = re.fullmatch(r"level_(\d+)", text)
        if not match:
            raise ValueError(f"Invalid Poisson multigrid level name '{level_name}'. Expected 'level_N'.")
        return match.group(1)

    def _append_poisson_solver_flags(ps: dict, source_key: str):
        """!
        @brief Append structured Poisson solver options to the flat PETSc flag map.
        @param[in] ps The `poisson_solver` or legacy `pressure_solver` mapping.
        @param[in] source_key Name of the source YAML block, used in error messages.
        """
        if not isinstance(ps, dict):
            raise ValueError(f"{source_key} must be a mapping when provided.")

        method = None
        if 'method' in ps:
            method = _normalize_poisson_method(ps['method'])
            flags['-ps_ksp_type'] = method
        if 'absolute_tolerance' in ps:
            flags['-ps_ksp_atol'] = ps['absolute_tolerance']
            flags['-poisson_tol'] = ps['absolute_tolerance']
        if 'relative_tolerance' in ps:
            flags['-ps_ksp_rtol'] = ps['relative_tolerance']
        if 'max_iterations' in ps:
            flags['-ps_ksp_max_it'] = ps['max_iterations']
        if 'tolerance' in ps:
            flags['-poisson_tol'] = ps['tolerance']

        gmres_cfg = ps.get('gmres', {})
        if gmres_cfg is not None:
            if not isinstance(gmres_cfg, dict):
                raise ValueError(f"{source_key}.gmres must be a mapping when provided.")
            if 'restart' in gmres_cfg:
                if method is None:
                    method = _normalize_poisson_method(ps.get('method', 'fgmres'))
                    flags.setdefault('-ps_ksp_type', method)
                if method not in {"gmres", "fgmres", "lgmres"}:
                    raise ValueError(
                        f"{source_key}.gmres.restart is valid only when {source_key}.method "
                        "is one of 'gmres', 'fgmres', or 'lgmres'."
                    )
                flags['-ps_ksp_gmres_restart'] = gmres_cfg['restart']

        preconditioner_cfg = ps.get('preconditioner', {})
        if preconditioner_cfg:
            if not isinstance(preconditioner_cfg, dict):
                raise ValueError(f"{source_key}.preconditioner must be a mapping when provided.")
            if 'type' in preconditioner_cfg:
                flags['-ps_pc_type'] = _normalize_poisson_preconditioner(preconditioner_cfg['type'])

        if 'multigrid' in ps:
            mg = ps['multigrid']
            if not isinstance(mg, dict):
                raise ValueError(f"{source_key}.multigrid must be a mapping when provided.")
            mg_map = {'levels': '-mg_level', 'pre_sweeps': '-mg_pre_it', 'post_sweeps': '-mg_post_it'}
            for key, flag in mg_map.items():
                if key in mg: flags[flag] = mg[key]
            if 'cycle' in mg:
                cycle = str(mg['cycle']).strip().lower()
                if cycle not in {"v"}:
                    raise ValueError(f"{source_key}.multigrid.cycle currently supports only 'v'.")
            if 'mode' in mg:
                mode = str(mg['mode']).strip().lower()
                if mode not in {"multiplicative"}:
                    raise ValueError(f"{source_key}.multigrid.mode currently supports only 'multiplicative'.")
            if 'semi_coarsening' in mg:
                sc = mg['semi_coarsening']
                if not isinstance(sc, dict):
                    raise ValueError(f"{source_key}.multigrid.semi_coarsening must be a mapping when provided.")
                if 'i' in sc: flags['-mg_i_semi'] = format_flag_value(sc['i'])
                if 'j' in sc: flags['-mg_j_semi'] = format_flag_value(sc['j'])
                if 'k' in sc: flags['-mg_k_semi'] = format_flag_value(sc['k'])
            if 'level_solvers' in mg:
                level_solvers = mg['level_solvers']
                if not isinstance(level_solvers, dict):
                    raise ValueError(f"{source_key}.multigrid.level_solvers must be a mapping when provided.")
                for level_name, settings in level_solvers.items():
                    if not isinstance(settings, dict):
                        raise ValueError(f"{source_key}.multigrid.level_solvers.{level_name} must be a mapping.")
                    level_num = _poisson_level_number(level_name)
                    for key, value in settings.items():
                        mapped_key = {'method': 'ksp_type', 'preconditioner': 'pc_type'}.get(key, key)
                        flags[f"-ps_mg_levels_{level_num}_{mapped_key}"] = format_flag_value(value)

    if 'poisson_solver' in solver_cfg and 'pressure_solver' in solver_cfg:
        if solver_cfg['poisson_solver'] != solver_cfg['pressure_solver']:
            raise ValueError(
                "Both 'poisson_solver' and legacy 'pressure_solver' are present with different values. "
                "Use 'poisson_solver' only, or make the legacy alias identical."
            )
    poisson_cfg = solver_cfg.get('poisson_solver', solver_cfg.get('pressure_solver'))
    if poisson_cfg is not None:
        source_key = 'poisson_solver' if 'poisson_solver' in solver_cfg else 'pressure_solver'
        _append_poisson_solver_flags(poisson_cfg, source_key)
    interp_cfg = solver_cfg.get('interpolation', {})
    if isinstance(interp_cfg, dict):
        interp_method_str = interp_cfg.get('method', 'Trilinear')
    else:
        interp_method_str = 'Trilinear'
    interp_code = normalize_interpolation_method(interp_method_str)
    flags['-interpolation_method'] = interp_code
    print(f"  - Interpolation Method: {interp_method_str} (Code: {interp_code})")

    if 'petsc_passthrough_options' in solver_cfg:
        passthrough = solver_cfg['petsc_passthrough_options']
        if passthrough:
            for key, value in passthrough.items():
                flags[key] = format_flag_value(value)
    summary_bits = []
    if '-ps_ksp_type' in flags:
        summary_bits.append(f"method={flags['-ps_ksp_type']}")
    if '-ps_ksp_atol' in flags:
        summary_bits.append(f"atol={flags['-ps_ksp_atol']}")
    if '-ps_ksp_rtol' in flags:
        summary_bits.append(f"rtol={flags['-ps_ksp_rtol']}")
    if '-ps_ksp_max_it' in flags:
        summary_bits.append(f"max_it={flags['-ps_ksp_max_it']}")
    if '-mg_level' in flags:
        summary_bits.append(f"mg_levels={flags['-mg_level']}")
    if summary_bits:
        print(f"  - Poisson Solver: {', '.join(summary_bits)}")
    return flags

def generate_solver_control_file(run_dir, run_id, configs, num_procs, monitor_files, restart_source_dir=None, continue_mode=False):
    """!
    @brief Generates the main .control file for the C-solver.
    @details Orchestrates the conversion of all YAML configurations (case, solver, monitor)
    into a single, machine-readable file of command-line flags.
    @param[in] run_dir Argument passed to `generate_solver_control_file()`.
    @param[in] run_id Argument passed to `generate_solver_control_file()`.
    @param[in] configs Argument passed to `generate_solver_control_file()`.
    @param[in] num_procs Argument passed to `generate_solver_control_file()`.
    @param[in] monitor_files Argument passed to `generate_solver_control_file()`.
    @param[in] restart_source_dir Argument passed to `generate_solver_control_file()`.
    @param[in] continue_mode If True, appends -continue_mode flag for the C solver.
    @return Value returned by `generate_solver_control_file()`.
    """
    print("[INFO] Generating master solver control file...")
    case_cfg, solver_cfg, monitor_cfg = configs['case'], configs['solver'], configs['monitor']
    source_files = {'Case': configs['case_path'], 'Solver': configs['solver_path'], 'Monitor': configs['monitor_path']}
    
    control_lines = []
    try:
        props, run_ctrl = case_cfg['properties'], case_cfg['run_control']
        scales, fluid, ic = props['scaling'], props['fluid'], props['initial_conditions']
        prepared_blocks = validate_and_prepare_boundary_conditions(case_cfg)
        L_ref, U_ref, rho, mu = float(scales['length_ref']), float(scales['velocity_ref']), float(fluid['density']), float(fluid['viscosity'])
        reynolds = (rho * U_ref * L_ref) / mu if mu != 0 else float('inf')
        dt_phys = float(run_ctrl['dt_physical'])
        T_ref = L_ref / U_ref if U_ref != 0 else float('inf')
        dt_nondim = dt_phys / T_ref if T_ref != float('inf') else 0.0
        resolved_ic = resolve_initial_condition_config(ic, prepared_blocks, U_ref)
        finit_mode_str = resolved_ic["label"]
        finit_code = resolved_ic["finit"]
        ic_params = resolved_ic["cli_params"]
        print(f"  - Reynolds Number (Re) = {reynolds:.4f}")
        print(f"  - Non-Dimensional dt*  = {dt_nondim:.6f}")
        eulerian_source = normalize_eulerian_field_source(
            (solver_cfg.get("operation_mode", {}) or {}).get("eulerian_field_source", "solve")
        )
        start_step = int(run_ctrl.get("start_step", 0) or 0)
        ic_is_authoritative = eulerian_source == "solve" and start_step == 0
        ic_cli = []
        if ic_is_authoritative:
            print(f"  - Initial Condition: {finit_mode_str} (Code: {finit_code})")
            if "ucont_x" in ic_params:
                ic_cli.extend([
                    f"-ucont_x {ic_params['ucont_x']}",
                    f"-ucont_y {ic_params['ucont_y']}",
                    f"-ucont_z {ic_params['ucont_z']}",
                ])
            if "ic_velocity_physical" in ic_params:
                ic_cli.append(f"-ic_velocity_physical {ic_params['ic_velocity_physical']}")
            if "flow_direction" in ic_params:
                ic_cli.append(f"-flow_direction {ic_params['flow_direction']}")
        else:
            print(
                f"[WARN] Ignoring configured initial condition '{finit_mode_str}' because "
                f"eulerian_field_source={eulerian_source!r} and start_step={start_step} select another source.",
                file=sys.stderr,
            )
            finit_code = 0
        control_lines.extend([
            f"-start_step {run_ctrl['start_step']}", f"-totalsteps {run_ctrl['total_steps']}",
            f"-ren {reynolds}", f"-dt {dt_nondim}", f"-finit {finit_code}",
            *ic_cli,
            f"-scaling_L_ref {L_ref}", f"-scaling_U_ref {U_ref}", f"-scaling_rho_ref {rho}"
        ])
    except (KeyError, TypeError, ZeroDivisionError, ValueError) as e:
        print(f"[FATAL] Error processing case.yml: {e}", file=sys.stderr)
        sys.exit(1)
        
    # --- CORRECTED: Add paths for whitelist and profile files ---
    if monitor_files.get("whitelist"):
        control_lines.append(f"-whitelist_config_file {monitor_files['whitelist']}")
    if monitor_files.get("profile"):
        control_lines.append(f"-profile_config_file {monitor_files['profile']}")
    profiling_cfg = monitor_files.get("profiling", {})
    control_lines.append(f"-profiling_timestep_mode {profiling_cfg.get('mode', 'off')}")
    if profiling_cfg.get("mode") != "off":
        control_lines.append(f"-profiling_timestep_file {profiling_cfg.get('timestep_file', 'Profiling_Timestep_Summary.csv')}")
    control_lines.append(f"-profiling_final_summary {str(bool(profiling_cfg.get('final_summary_enabled', True))).lower()}")
    diagnostics_cfg = resolve_diagnostics_config(monitor_cfg)
    memory_log_cfg = diagnostics_cfg["runtime_memory_log"]
    control_lines.append(f"-runtime_memory_log_enabled {str(bool(memory_log_cfg.get('enabled', True))).lower()}")
    control_lines.append(f"-runtime_memory_log_file {memory_log_cfg.get('file', 'Runtime_Memory.log')}")

    walltime_guard_policy = configs.get("walltime_guard_policy")
    if walltime_guard_policy is not None:
        control_lines.extend(
            [
                f"-walltime_guard_enabled {str(bool(walltime_guard_policy.get('enabled', False))).lower()}",
                f"-walltime_guard_warmup_steps {int(walltime_guard_policy.get('warmup_steps', DEFAULT_WALLTIME_GUARD_POLICY['warmup_steps']))}",
                f"-walltime_guard_multiplier {float(walltime_guard_policy.get('multiplier', DEFAULT_WALLTIME_GUARD_POLICY['multiplier']))}",
                f"-walltime_guard_min_seconds {float(walltime_guard_policy.get('min_seconds', DEFAULT_WALLTIME_GUARD_POLICY['min_seconds']))}",
                f"-walltime_guard_estimator_alpha {float(walltime_guard_policy.get('estimator_alpha', DEFAULT_WALLTIME_GUARD_POLICY['estimator_alpha']))}",
            ]
        )
    
    grid_cfg = case_cfg.get('grid', {})
    grid_mode = grid_cfg.get('mode')
    expected_nblk = int(case_cfg.get('models', {}).get('domain', {}).get('blocks', 1))

    if grid_mode == 'file':
        print("[INFO] Grid Mode: Using external file...")
        case_file_dir = os.path.dirname(configs['case_path'])
        source_grid = grid_cfg['source_file']
        if not os.path.isabs(source_grid):
            source_grid = os.path.abspath(os.path.join(case_file_dir, source_grid))
        grid_for_validation = source_grid
        legacy_cfg = grid_cfg.get("legacy_conversion")
        if isinstance(legacy_cfg, dict):
            if legacy_cfg.get("enabled", True):
                print("[INFO] Grid file legacy conversion enabled; converting with grid.gen...")
            grid_for_validation = convert_legacy_grid_with_gridgen(
                configs['case_path'],
                run_dir,
                grid_cfg,
                source_grid,
            )
        nondim_grid_path = os.path.join(run_dir, "config", "grid.run")
        try:
            summary = validate_and_nondimensionalize_picgrid(
                grid_for_validation, nondim_grid_path, L_ref, expected_nblk=expected_nblk
            )
            print(
                f"[SUCCESS] Validated and non-dimensionalized grid: {os.path.relpath(nondim_grid_path)} "
                f"(nblk={summary['nblk']}, total_nodes={summary['total_nodes']})"
            )
            control_lines.append(f"-grid_file {nondim_grid_path}")
        except Exception as e:
            print(f"[FATAL] Failed to process grid file '{source_grid}': {e}", file=sys.stderr)
            sys.exit(1)
    elif grid_mode == 'grid_gen':
        print("[INFO] Grid Mode: Generating external grid via grid.gen...")
        nondim_grid_path = os.path.join(run_dir, "config", "grid.run")
        if continue_mode and os.path.isfile(nondim_grid_path):
            print(f"[INFO] Continue mode: reusing staged grid: {os.path.relpath(nondim_grid_path)}")
            control_lines.append(f"-grid_file {nondim_grid_path}")
        else:
            try:
                # grid.gen already converts ncells_* inputs into node-count PICGRID dims.
                generated_grid = run_grid_generator(configs['case_path'], run_dir, grid_cfg)
                summary = validate_and_nondimensionalize_picgrid(
                    generated_grid, nondim_grid_path, L_ref, expected_nblk=expected_nblk
                )
                print(
                    f"[SUCCESS] grid.gen output validated and non-dimensionalized: {os.path.relpath(nondim_grid_path)} "
                    f"(nblk={summary['nblk']}, total_nodes={summary['total_nodes']})"
                )
                control_lines.append(f"-grid_file {nondim_grid_path}")
            except Exception as e:
                print(f"[FATAL] Grid generation failed: {e}", file=sys.stderr)
                sys.exit(1)
    elif grid_mode == 'programmatic_c':
        print("[INFO] Grid Mode: Programmatic C...")
        grid_settings = translate_programmatic_grid_settings(grid_cfg.get('programmatic_settings', {}))
        control_lines.append("-grid")
        for p_key in GRID_DA_PROCESSOR_KEYS:
            grid_settings.pop(p_key, None)
        for key, value in grid_settings.items(): control_lines.append(f"-{key} {format_flag_value(value)}")
        if resolved_ic["kind"] == "ic_gen" and ic_is_authoritative:
            nondim_grid_path = os.path.join(run_dir, "config", "grid.run")
            try:
                summary = generate_picgrid_from_programmatic_settings(
                    grid_cfg.get('programmatic_settings', {}), nondim_grid_path, L_ref
                )
                print(
                    f"[INFO] Materialized grid.run for ic_gen: {os.path.relpath(nondim_grid_path)} "
                    f"(nblk={summary['nblk']}, total_nodes={summary['total_nodes']})"
                )
            except Exception as e:
                print(f"[FATAL] Failed to generate grid.run for ic_gen: {e}", file=sys.stderr)
                sys.exit(1)
    else:
        raise ValueError(f"Unknown or missing grid mode '{grid_mode}' in case.yml.")

    if resolved_ic["kind"] in {"file", "ic_gen"}:
        if ic_is_authoritative:
            try:
                staged_ic = stage_initial_condition_file(run_dir, configs["case_path"], resolved_ic)
            except Exception as e:
                print(f"[FATAL] Failed to stage initial condition: {e}", file=sys.stderr)
                sys.exit(1)
            control_lines.extend([
                f"-ic_field {resolved_ic['field_code']}",
                f"-ic_dir {staged_ic['directory']}",
            ])
            print(f"  - Staged initial condition: {os.path.relpath(staged_ic['staged'])}")

    try:
        bcs_files = generate_multi_block_bcs(run_dir, run_id, case_cfg, source_files)
    except ValueError as e:
        print(f"[FATAL] Invalid boundary_conditions in case.yml: {e}", file=sys.stderr)
        sys.exit(1)
    control_lines.append(f"-bcs_files \"{','.join(bcs_files)}\"")

    append_grid_da_processor_layout(control_lines, grid_cfg, num_procs)
    
    parse_and_add_model_flags(case_cfg, control_lines)
    
    if 'solver_parameters' in case_cfg:
        params = case_cfg['solver_parameters']
        if params:
            for key, value in params.items():
                control_lines.append(f"{key} {format_flag_value(value)}")

    try:
        solver_flags = parse_solver_config(solver_cfg)
    except ValueError as e:
        print(f"[FATAL] Invalid solver.yml settings: {e}", file=sys.stderr)
        sys.exit(1)
    for flag, value in solver_flags.items(): control_lines.append(f"{flag} {value}")

    try:
        solver_monitoring_flags = resolve_solver_monitoring_flags(monitor_cfg)
    except ValueError as e:
        print(f"[FATAL] Invalid monitor.yml solver_monitoring settings: {e}", file=sys.stderr)
        sys.exit(1)
    append_passthrough_flags(control_lines, solver_monitoring_flags)
    
    io_cfg = monitor_cfg.get('io', {})
    particle_console_output_freq = resolve_particle_console_output_frequency(io_cfg)
    if 'data_output_frequency' in io_cfg: control_lines.append(f"-tio {io_cfg['data_output_frequency']}")
    if particle_console_output_freq is not None:
        control_lines.append(f"-particle_console_output_freq {particle_console_output_freq}")
    if 'particle_log_interval' in io_cfg: control_lines.append(f"-logfreq {io_cfg['particle_log_interval']}")
    if 'directories' in io_cfg:
        dirs = io_cfg['directories']
        if 'output' in dirs: control_lines.append(f"-output_dir {dirs['output']}")
        if 'restart' in dirs and not restart_source_dir:
            control_lines.append(f"-restart_dir {dirs['restart']}")
        if 'log' in dirs: control_lines.append(f"-log_dir {dirs['log']}")
        if 'eulerian_subdir' in dirs: control_lines.append(f"-euler_subdir {dirs['eulerian_subdir']}")
        if 'particle_subdir' in dirs: control_lines.append(f"-particle_subdir {dirs['particle_subdir']}")
    if restart_source_dir:
        control_lines.append(f"-restart_dir {restart_source_dir}")
    if continue_mode:
        control_lines.append("-continue_mode true")

    final_content = generate_header(run_id, source_files) + "\n".join(control_lines)
    control_file_path = os.path.join(run_dir, "config", f"{run_id}.control")
    with open(control_file_path, "w") as f: f.write(final_content)
    print(f"[SUCCESS] Generated solver control file: {os.path.relpath(control_file_path)}")
    return os.path.abspath(control_file_path)

def generate_post_recipe_file(run_dir: str, run_id: str, post_cfg: dict, source_files: dict, monitor_cfg=None) -> str:
    """!
    @brief Generates a key=value config file (post.run) for the C post-processor.
    @details Translates the structured post-processing YAML into the specific flat
             key-value format required by the C executable, including complex,
             semicolon-separated pipeline strings.
    @param[in] run_dir The path to the main run directory.
    @param[in] run_id The unique identifier for the run.
    @param[in] post_cfg The parsed post-profile YAML configuration dictionary.
    @param[in] source_files A dictionary of source files for the header.
    @param[in] monitor_cfg Optional parsed monitor YAML configuration dictionary.
    @return The absolute path to the generated post.run recipe file.
    """
    print("[INFO] Generating post-processor recipe file (post.run)...")
    config_dir = os.path.join(run_dir, "config")
    post_recipe_path = os.path.join(config_dir, "post.run")

    lines = [generate_header(run_id, source_files)]
    c_config = build_post_recipe_config(post_cfg, monitor_cfg)

    for key, value in c_config.items():
        if value is not None and str(value) != "":
            lines.append(f"{key} = {value}")

    with open(post_recipe_path, "w") as f:
        f.write("\n".join(lines))
    print(f"[SUCCESS] Generated post-processor recipe: {os.path.relpath(post_recipe_path)}")
    return os.path.abspath(post_recipe_path)

def execute_command(command: list, run_dir: str, log_filename: str, monitor_cfg: dict = None):
    """!
    @brief Executes a command, streaming its output to the console and a log file.
    @details ...
    If None, the process inherits the parent's environment directly.
    @param[in] command Argument passed to `execute_command()`.
    @param[in] run_dir Argument passed to `execute_command()`.
    @param[in] log_filename Argument passed to `execute_command()`.
    @param[in] monitor_cfg Argument passed to `execute_command()`.
    """
    log_path = resolve_command_log_path(run_dir, log_filename)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    print(f"[INFO] Launching Command...\n  > {format_command_for_display(command)}")
    print(f"       Log file: {os.path.relpath(log_path)}")
    print("-" * 60)

    # --- Environment Handling ---
    popen_kwargs = {
        "stdout": subprocess.PIPE, "stderr": subprocess.STDOUT,
        "cwd": run_dir, "bufsize": 1, "universal_newlines": True,
        "encoding": 'utf-8', "errors": 'replace'
    }

    if monitor_cfg:
        print("[INFO] Creating custom environment to set LOG_LEVEL.")
        run_env = os.environ.copy()
        verbosity = monitor_cfg.get('logging', {}).get('verbosity', 'INFO').upper()
        run_env['LOG_LEVEL'] = verbosity
        print(f"[INFO] Setting LOG_LEVEL={verbosity} for C executable.")
        popen_kwargs['env'] = run_env
    else:
        print("[INFO] Using inherited environment for process.")

    print("-" * 60)
    try:
        # Pass the constructed keyword arguments dictionary to Popen
        process = subprocess.Popen(command, **popen_kwargs)

        with open(log_path, "w") as log_file:
            for line in process.stdout:
                sys.stdout.write(line)
                log_file.write(line)
        process.wait()
        return_code = process.returncode
        print("-" * 60)
        if return_code == 0:
            print(f"[SUCCESS] Execution finished successfully.")
        else:
            print(f"[FATAL] Execution failed with exit code {return_code}. Check log: {os.path.relpath(log_path)}", file=sys.stderr)
            sys.exit(return_code)
    except FileNotFoundError:
        print(f"[FATAL] Command not found or is not executable: '{command[0]}'", file=sys.stderr)
        print("        Please check that the path is correct and the file has execute permissions.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[FATAL] An unexpected error occurred during execution: {e}", file=sys.stderr)
        sys.exit(1)


def format_command_for_display(command: list) -> str:
    """!
    @brief Render a shell-safe command string for console and log output.
    @param[in] command Argument passed to `format_command_for_display()`.
    @return Value returned by `format_command_for_display()`.
    """
    return " ".join(shlex.quote(str(part)) for part in command)


def resolve_command_log_path(run_dir: str, log_filename: str) -> str:
    """!
    @brief Resolve a command log filename relative to the run directory.
    @param[in] run_dir Argument passed to `resolve_command_log_path()`.
    @param[in] log_filename Argument passed to `resolve_command_log_path()`.
    @return Value returned by `resolve_command_log_path()`.
    """
    if os.path.dirname(log_filename):
        return os.path.join(run_dir, log_filename)
    return os.path.join(run_dir, "logs", log_filename)


class CommandExecutionError(RuntimeError):
    """!
    @brief Raised when an external command exits unsuccessfully.
    """

    def __init__(self, command: list, returncode: int, details: str = None):
        """!
        @brief Initialize a command execution error.
        @param[in] command Argument passed to `__init__()`.
        @param[in] returncode Argument passed to `__init__()`.
        @param[in] details Argument passed to `__init__()`.
        """
        self.command = command
        self.returncode = returncode
        self.details = details
        detail_suffix = f": {details}" if details else ""
        super().__init__(
            f"Command failed with exit code {returncode}: {format_command_for_display(command)}{detail_suffix}"
        )


class PlotDependencyError(RuntimeError):
    """!
    @brief Raised when plot.gen reports a missing optional dependency.
    """


def _run_captured_command(command: list, run_dir: str) -> subprocess.CompletedProcess:
    """!
    @brief Run a command and capture combined stdout/stderr details for later inspection.
    @param[in] command Argument passed to `_run_captured_command()`.
    @param[in] run_dir Argument passed to `_run_captured_command()`.
    @return Value returned by `_run_captured_command()`.
    """
    try:
        return subprocess.run(
            command,
            cwd=run_dir,
            text=True,
            capture_output=True,
            check=False,
            encoding="utf-8",
            errors="replace",
        )
    except FileNotFoundError as exc:
        raise CommandExecutionError(command, 1, f"Command not found or is not executable: '{command[0]}'") from exc


def _require_successful_command(command: list, result: subprocess.CompletedProcess):
    """!
    @brief Raise `CommandExecutionError` when a captured command failed.
    @param[in] command Argument passed to `_require_successful_command()`.
    @param[in] result Argument passed to `_require_successful_command()`.
    """
    if result.returncode == 0:
        return
    details = (result.stderr or result.stdout).strip()
    raise CommandExecutionError(command, result.returncode, details or None)


def _capture_command_stdout(command: list, run_dir: str) -> str:
    """!
    @brief Run a command, require success, and return stripped stdout text.
    @param[in] command Argument passed to `_capture_command_stdout()`.
    @param[in] run_dir Argument passed to `_capture_command_stdout()`.
    @return Value returned by `_capture_command_stdout()`.
    """
    result = _run_captured_command(command, run_dir)
    _require_successful_command(command, result)
    return result.stdout.strip()


def _stream_command_to_console_and_log(command: list, run_dir: str, log_file):
    """!
    @brief Stream command output to stdout and an already-open log file.
    @param[in] command Argument passed to `_stream_command_to_console_and_log()`.
    @param[in] run_dir Argument passed to `_stream_command_to_console_and_log()`.
    @param[in] log_file Argument passed to `_stream_command_to_console_and_log()`.
    """
    display = format_command_for_display(command)
    print(f"[INFO] Running: {display}")
    log_file.write(f"$ {display}\n")
    log_file.flush()

    popen_kwargs = {
        "stdout": subprocess.PIPE,
        "stderr": subprocess.STDOUT,
        "cwd": run_dir,
        "bufsize": 1,
        "universal_newlines": True,
        "encoding": "utf-8",
        "errors": "replace",
    }

    try:
        process = subprocess.Popen(command, **popen_kwargs)
    except FileNotFoundError as exc:
        raise CommandExecutionError(command, 1, f"Command not found or is not executable: '{command[0]}'") from exc

    with process:
        for line in process.stdout:
            sys.stdout.write(line)
            log_file.write(line)
    return_code = process.wait()
    log_file.write("\n")
    log_file.flush()
    if return_code != 0:
        raise CommandExecutionError(command, return_code)


def _get_git_head_state(run_dir: str) -> dict:
    """!
    @brief Capture the current git HEAD branch name and commit hash.
    @param[in] run_dir Argument passed to `_get_git_head_state()`.
    @return Value returned by `_get_git_head_state()`.
    """
    head_commit = _capture_command_stdout(["git", "rev-parse", "--verify", "HEAD"], run_dir)
    branch_result = _run_captured_command(["git", "symbolic-ref", "--quiet", "--short", "HEAD"], run_dir)
    branch_name = branch_result.stdout.strip() if branch_result.returncode == 0 else None
    return {"branch": branch_name, "commit": head_commit}


def _get_local_branches_with_upstreams(run_dir: str) -> "list[tuple[str, str | None]]":
    """!
    @brief Return local branch names plus their configured upstreams.
    @param[in] run_dir Argument passed to `_get_local_branches_with_upstreams()`.
    @return Value returned by `_get_local_branches_with_upstreams()`.
    """
    output = _capture_command_stdout(
        ["git", "for-each-ref", "--sort=refname", "--format=%(refname:short)\t%(upstream:short)", "refs/heads"],
        run_dir,
    )
    branches = []
    for line in output.splitlines():
        if not line.strip():
            continue
        branch_name, _, upstream_name = line.partition("\t")
        branches.append((branch_name, upstream_name or None))
    return branches


def _working_tree_has_tracked_changes(run_dir: str) -> bool:
    """!
    @brief Return `True` when the repository has staged or unstaged tracked changes.
    @param[in] run_dir Argument passed to `_working_tree_has_tracked_changes()`.
    @return Value returned by `_working_tree_has_tracked_changes()`.
    """
    command = ["git", "status", "--porcelain", "--untracked-files=no"]
    result = _run_captured_command(command, run_dir)
    _require_successful_command(command, result)
    return bool(result.stdout.strip())


def _attempt_pull_cleanup(run_dir: str, rebase: bool, log_file):
    """!
    @brief Best-effort cleanup after a failed `git pull` so the original branch can be restored.
    @param[in] run_dir Argument passed to `_attempt_pull_cleanup()`.
    @param[in] rebase Argument passed to `_attempt_pull_cleanup()`.
    @param[in] log_file Argument passed to `_attempt_pull_cleanup()`.
    """
    cleanup_command = ["git", "rebase", "--abort"] if rebase else ["git", "merge", "--abort"]
    result = _run_captured_command(cleanup_command, run_dir)
    if result.returncode == 0:
        print(f"[INFO] Cleaned up the interrupted {'rebase' if rebase else 'merge'} state.")
        log_file.write(f"$ {format_command_for_display(cleanup_command)}\n")
        if result.stdout:
            sys.stdout.write(result.stdout)
            log_file.write(result.stdout)
        if result.stderr:
            sys.stderr.write(result.stderr)
            log_file.write(result.stderr)
        log_file.write("\n")
        log_file.flush()
        return

    details = (result.stderr or result.stdout).strip()
    if details:
        message = (
            f"[WARNING] Could not clean up a failed {'rebase' if rebase else 'merge'} automatically: {details}"
        )
        print(message, file=sys.stderr)
        log_file.write(message + "\n")
        log_file.flush()


def _restore_git_head(run_dir: str, original_head: dict, log_file):
    """!
    @brief Restore the repository back to the branch or detached commit it started on.
    @param[in] run_dir Argument passed to `_restore_git_head()`.
    @param[in] original_head Argument passed to `_restore_git_head()`.
    @param[in] log_file Argument passed to `_restore_git_head()`.
    """
    current_state = _get_git_head_state(run_dir)
    if original_head["branch"]:
        if current_state["branch"] == original_head["branch"]:
            return
        _stream_command_to_console_and_log(["git", "checkout", original_head["branch"]], run_dir, log_file)
        return

    if current_state["branch"] is None and current_state["commit"] == original_head["commit"]:
        return
    _stream_command_to_console_and_log(["git", "checkout", "--detach", original_head["commit"]], run_dir, log_file)


def pull_all_source_branches(run_dir: str, log_filename: str, rebase: bool = True):
    """!
    @brief Refresh every local tracking branch in the source repository, then restore the starting branch.
    @param[in] run_dir Argument passed to `pull_all_source_branches()`.
    @param[in] log_filename Argument passed to `pull_all_source_branches()`.
    @param[in] rebase Argument passed to `pull_all_source_branches()`.
    """
    log_path = resolve_command_log_path(run_dir, log_filename)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    print("\n" + "="*23 + " PULL SOURCE STAGE " + "="*22)
    print("[INFO] Refreshing all local source branches that track an upstream.")
    print(f"       Log file: {os.path.relpath(log_path)}")
    print("-" * 60)

    try:
        original_head = _get_git_head_state(run_dir)
        if _working_tree_has_tracked_changes(run_dir):
            raise RuntimeError(
                "Multi-branch pull requires a clean tracked working tree in the source repository. "
                "Commit or stash those changes first, or rerun with --current-branch-only."
            )
        branches = _get_local_branches_with_upstreams(run_dir)
    except (CommandExecutionError, RuntimeError) as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(getattr(exc, "returncode", 1))

    if not branches:
        print("[FATAL] No local branches were found in the source repository.", file=sys.stderr)
        sys.exit(1)

    if original_head["branch"]:
        branches = [item for item in branches if item[0] != original_head["branch"]] + [
            item for item in branches if item[0] == original_head["branch"]
        ]

    skipped_branches = []
    current_operation = None
    pull_error = None
    restore_error = None

    with open(log_path, "w", encoding="utf-8") as log_file:
        log_file.write(f"# PICurv pull-source all-branch sync\n")
        log_file.write(f"# repository: {os.path.abspath(run_dir)}\n")
        log_file.write(f"# started: {datetime.now().isoformat()}\n")
        log_file.write(
            f"# original head: {original_head['branch'] if original_head['branch'] else original_head['commit']}\n\n"
        )

        try:
            for branch_name, upstream_name in branches:
                if not upstream_name:
                    warning = f"[WARNING] Skipping branch '{branch_name}' because it has no configured upstream."
                    print(warning, file=sys.stderr)
                    log_file.write(warning + "\n")
                    skipped_branches.append(branch_name)
                    continue

                print(f"[INFO] Refreshing branch '{branch_name}' from '{upstream_name}'.")
                log_file.write(f"[INFO] Refreshing branch '{branch_name}' from '{upstream_name}'.\n")

                current_operation = f"checkout:{branch_name}"
                _stream_command_to_console_and_log(["git", "checkout", branch_name], run_dir, log_file)

                pull_command = ["git", "pull"]
                if rebase:
                    pull_command.append("--rebase")
                current_operation = f"pull:{branch_name}"
                _stream_command_to_console_and_log(pull_command, run_dir, log_file)
                current_operation = None
        except CommandExecutionError as exc:
            pull_error = exc
            if current_operation and current_operation.startswith("pull:"):
                _attempt_pull_cleanup(run_dir, rebase, log_file)
        finally:
            try:
                _restore_git_head(run_dir, original_head, log_file)
            except CommandExecutionError as exc:
                restore_error = exc

    print("-" * 60)
    if pull_error:
        if restore_error:
            print(
                f"[FATAL] Multi-branch pull failed and the original branch could not be restored. "
                f"Check log: {os.path.relpath(log_path)}",
                file=sys.stderr,
            )
            sys.exit(restore_error.returncode)
        print(
            f"[FATAL] Multi-branch pull failed. Original branch restored. "
            f"Check log: {os.path.relpath(log_path)}",
            file=sys.stderr,
        )
        sys.exit(pull_error.returncode)

    if restore_error:
        print(
            f"[FATAL] Branch updates completed, but the original branch could not be restored. "
            f"Check log: {os.path.relpath(log_path)}",
            file=sys.stderr,
        )
        sys.exit(restore_error.returncode)

    if skipped_branches:
        print(f"[WARNING] Skipped branches with no upstream: {', '.join(skipped_branches)}", file=sys.stderr)
    print("[SUCCESS] All local tracking branches are up to date.")

def auto_identify_run_inputs(config_dir: str):
    """!
    @brief Auto-detect case.yml, monitor.yml, and *.control in a run config directory.
    @param[in] config_dir Argument passed to `auto_identify_run_inputs()`.
    @return Value returned by `auto_identify_run_inputs()`.
    """
    all_yml_files = glob.glob(os.path.join(config_dir, "*.yml"))
    case_path, monitor_path = None, None
    for f_path in all_yml_files:
        try:
            content = read_yaml_file(f_path)
            if not isinstance(content, dict):
                continue
            if 'models' in content and 'boundary_conditions' in content:
                case_path = f_path
            elif 'io' in content and 'logging' in content:
                monitor_path = f_path
        except Exception as e:
            print(f"[WARNING] Could not parse or inspect '{f_path}': {e}", file=sys.stderr)
    try:
        solver_control_path = glob.glob(os.path.join(config_dir, "*.control"))[0]
    except IndexError:
        solver_control_path = None
    return case_path, monitor_path, solver_control_path

def resolve_post_source_directory(run_dir: str, monitor_cfg: dict, post_cfg: dict, strict: bool = True) -> str:
    """!
    @brief Resolve post source directory token and optionally enforce existence.
    @param[in] run_dir Argument passed to `resolve_post_source_directory()`.
    @param[in] monitor_cfg Argument passed to `resolve_post_source_directory()`.
    @param[in] post_cfg Argument passed to `resolve_post_source_directory()`.
    @param[in] strict Argument passed to `resolve_post_source_directory()`.
    @return Value returned by `resolve_post_source_directory()`.
    """
    solver_output_dir_rel = monitor_cfg.get('io', {}).get('directories', {}).get('output', 'output')
    solver_output_dir_abs = os.path.join(run_dir, solver_output_dir_rel)
    source_dir_template = get_post_source_directory_template(post_cfg)
    if source_dir_template == '<solver_output_dir>':
        resolved_source_dir = solver_output_dir_abs
        print(f"[INFO] Post-processor source data: {os.path.relpath(resolved_source_dir)}")
    else:
        resolved_source_dir = os.path.abspath(os.path.join(run_dir, source_dir_template))
        print(f"[INFO] Post-processor source data (user-defined): {os.path.relpath(resolved_source_dir)}")

    if strict and (not os.path.isdir(resolved_source_dir) or not os.listdir(resolved_source_dir)):
        print(
            f"[FATAL] Source data directory for post-processing not found or empty: {os.path.relpath(resolved_source_dir)}",
            file=sys.stderr
        )
        sys.exit(1)
    if not strict and (not os.path.isdir(resolved_source_dir) or not os.listdir(resolved_source_dir)):
        print("[WARNING] Source data directory is not available yet; keeping deferred path for scheduled post job.")
    return resolved_source_dir

def render_slurm_array_stage_script(
    script_path: str,
    job_name: str,
    cluster_cfg: dict,
    array_spec: str,
    case_index_tsv: str,
    stage: str,
    solver_exe: str,
    post_exe: str,
    stdout_path: str,
    stderr_path: str
):
    """!
    @brief Render array script that maps SLURM_ARRAY_TASK_ID to per-case run artifacts.
    @param[in] script_path Argument passed to `render_slurm_array_stage_script()`.
    @param[in] job_name Argument passed to `render_slurm_array_stage_script()`.
    @param[in] cluster_cfg Argument passed to `render_slurm_array_stage_script()`.
    @param[in] array_spec Argument passed to `render_slurm_array_stage_script()`.
    @param[in] case_index_tsv Argument passed to `render_slurm_array_stage_script()`.
    @param[in] stage Argument passed to `render_slurm_array_stage_script()`.
    @param[in] solver_exe Argument passed to `render_slurm_array_stage_script()`.
    @param[in] post_exe Argument passed to `render_slurm_array_stage_script()`.
    @param[in] stdout_path Argument passed to `render_slurm_array_stage_script()`.
    @param[in] stderr_path Argument passed to `render_slurm_array_stage_script()`.
    @return Value returned by `render_slurm_array_stage_script()`.
    """
    effective_cluster_cfg = build_serial_post_cluster_config(cluster_cfg) if stage == "post" else cluster_cfg
    resources = effective_cluster_cfg.get("resources", {})
    notifications = effective_cluster_cfg.get("notifications", {}) or {}
    execution = effective_cluster_cfg.get("execution", {}) or {}
    module_setup = execution.get("module_setup", []) or []
    extra_sbatch = execution.get("extra_sbatch")

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --nodes={resources['nodes']}",
        f"#SBATCH --ntasks-per-node={resources['ntasks_per_node']}",
        f"#SBATCH --mem={resources['mem']}",
        f"#SBATCH --time={resources['time']}",
        f"#SBATCH --output={stdout_path}",
        f"#SBATCH --error={stderr_path}",
        f"#SBATCH --account={resources['account']}",
        f"#SBATCH --array={array_spec}",
    ]
    partition = resources.get("partition")
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    mail_user = notifications.get("mail_user")
    mail_type = notifications.get("mail_type")
    if mail_user:
        lines.append(f"#SBATCH --mail-user={mail_user}")
    if mail_type:
        lines.append(f"#SBATCH --mail-type={mail_type}")
    if isinstance(extra_sbatch, dict):
        for key, value in extra_sbatch.items():
            flag = str(key)
            if not flag.startswith("--"):
                flag = f"--{flag}"
            if isinstance(value, bool):
                if value:
                    lines.append(f"#SBATCH {flag}")
            elif value is not None:
                lines.append(f"#SBATCH {flag}={value}")
    elif isinstance(extra_sbatch, list):
        for token in extra_sbatch:
            lines.append(f"#SBATCH {token}")

    lines.extend([
        "",
        "set -euo pipefail",
        "",
        f'CASE_INDEX_FILE={shlex.quote(case_index_tsv)}',
        'LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASE_INDEX_FILE")',
        'if [ -z "$LINE" ]; then',
        '  echo "No case entry for array index ${SLURM_ARRAY_TASK_ID}" >&2',
        '  exit 1',
        "fi",
        "IFS=$'\\t' read -r CASE_INDEX CASE_ID RUN_DIR CONTROL_FILE POST_RECIPE_FILE LOG_LEVEL POST_PREFIX SOLVE_DIAGNOSTIC_ARGS POST_DIAGNOSTIC_ARGS <<< \"$LINE\"",
        'cd "$RUN_DIR"',
        'echo "[$(date)] Starting case ${CASE_ID} (array index ${SLURM_ARRAY_TASK_ID})"',
    ])

    if stage == "solve":
        walltime_guard_exports = build_walltime_guard_exports(effective_cluster_cfg)
        for key, value in walltime_guard_exports.items():
            lines.append(f"export {key}={value}")

    lines.append('export LOG_LEVEL="${LOG_LEVEL}"')

    for setup_line in module_setup:
        lines.append(str(setup_line))

    if stage == "solve":
        cmd = build_cluster_launch_command(
            effective_cluster_cfg,
            solver_exe,
            ["-control_file", "$CONTROL_FILE"]
        )
    else:
        cmd = build_cluster_launch_command(
            effective_cluster_cfg,
            post_exe,
            ["-control_file", "$CONTROL_FILE", "-postprocessing_config_file", "$POST_RECIPE_FILE"],
            force_num_procs=1,
        )

    # Keep shell variables unresolved inside sbatch script.
    def _token(tok: str) -> str:
        """!
        @brief Perform token.
        @param[in] tok Argument passed to `_token()`.
        @return Value returned by `_token()`.
        """
        if tok.startswith("$"):
            return tok
        return shlex.quote(str(tok))

    diag_var = "${SOLVE_DIAGNOSTIC_ARGS}" if stage == "solve" else "${POST_DIAGNOSTIC_ARGS}"
    command_text = " ".join(_token(t) for t in cmd)
    executable_token = _token(solver_exe if stage == "solve" else post_exe)
    # Diagnostic args must be executable args, not launcher args. Insert them
    # immediately before the executable's normal control/recipe options.
    if executable_token and command_text.count(executable_token) == 1:
        command_text = command_text.replace(f"{executable_token} ", f"{executable_token} {diag_var} ", 1)
    lines.append(f"exec {command_text}")

    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)


def render_metrics_aggregate_script(
    script_path: str,
    job_name: str,
    cluster_cfg: dict,
    study_dir: str,
    picurv_path: str,
):
    """!
    @brief Generate a single-node sbatch script that runs metrics aggregation.
    @param[in] script_path Path to write the sbatch script.
    @param[in] job_name Slurm job name.
    @param[in] cluster_cfg Parsed cluster YAML dictionary.
    @param[in] study_dir Absolute path to the study directory.
    @param[in] picurv_path Absolute path to the picurv script.
    """
    resources = cluster_cfg.get("resources", {})
    notifications = cluster_cfg.get("notifications", {}) or {}
    execution = cluster_cfg.get("execution", {}) or {}
    module_setup = execution.get("module_setup", []) or []

    scheduler_dir = os.path.join(study_dir, "scheduler")
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks-per-node=1",
        "#SBATCH --mem=4G",
        "#SBATCH --time=00:10:00",
        f"#SBATCH --output={os.path.join(scheduler_dir, 'metrics_%j.out')}",
        f"#SBATCH --error={os.path.join(scheduler_dir, 'metrics_%j.err')}",
        f"#SBATCH --account={resources['account']}",
    ]
    partition = resources.get("partition")
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    mail_user = notifications.get("mail_user")
    mail_type = notifications.get("mail_type")
    if mail_user:
        lines.append(f"#SBATCH --mail-user={mail_user}")
    if mail_type:
        lines.append(f"#SBATCH --mail-type={mail_type}")

    lines.extend([
        "",
        "set -euo pipefail",
        'echo "[$(date)] Running metrics aggregation"',
        "",
    ])
    for setup_line in module_setup:
        lines.append(str(setup_line))

    lines.append(
        f"exec {shlex.quote(picurv_path)} sweep --reaggregate"
        f" --study-dir {shlex.quote(study_dir)}"
    )

    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)


def reduce_metric_values(values, reduction: str):
    """!
    @brief Reduce a metric series to one scalar according to the requested reducer.
    @param[in] values Sequence of numeric values.
    @param[in] reduction Reduction keyword.
    @return Value returned by `reduce_metric_values()`.
    """
    if not values:
        return None
    np = require_numpy()
    reduction = str(reduction).lower()
    if reduction == "mean":
        return float(np.mean(values))
    if reduction == "min":
        return float(np.min(values))
    if reduction == "max":
        return float(np.max(values))
    if reduction == "p95":
        return float(np.percentile(values, 95.0))
    return float(values[-1])


def extract_metric_from_csv(case_dir: str, spec: dict):
    """!
    @brief Extract a scalar metric from a CSV source.
    @param[in] case_dir Argument passed to `extract_metric_from_csv()`.
    @param[in] spec Argument passed to `extract_metric_from_csv()`.
    @return Value returned by `extract_metric_from_csv()`.
    """
    file_glob = spec.get("file_glob", "**/*_msd.csv")
    candidates = sorted(glob.glob(os.path.join(case_dir, file_glob), recursive=True))
    if not candidates:
        return None
    csv_path = candidates[0]
    rows = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames:
            for row in reader:
                rows.append(row)
            if not rows:
                return None
            column = spec.get("column")
            numerator_column = spec.get("numerator_column")
            denominator_column = spec.get("denominator_column")
            denominator_floor = float(spec.get("denominator_floor", 0.0) or 0.0)
            if not column and not numerator_column:
                for name in reversed(reader.fieldnames):
                    if name and name.lower() not in {"step", "time", "timestep"}:
                        column = name
                        break
            if not column and not numerator_column:
                return None
            values = []
            for row in rows:
                try:
                    if numerator_column:
                        numerator = float(row[numerator_column])
                        denominator = float(row[denominator_column])
                        denominator = max(denominator_floor, denominator)
                        if denominator == 0.0:
                            continue
                        values.append(numerator / denominator)
                    else:
                        values.append(float(row[column]))
                except Exception:
                    continue
        else:
            return None
    return reduce_metric_values(values, spec.get("reduction", "last"))


def extract_metric_from_log(case_dir: str, spec: dict):
    """!
    @brief Extract a scalar metric from a log file using regex.
    @param[in] case_dir Argument passed to `extract_metric_from_log()`.
    @param[in] spec Argument passed to `extract_metric_from_log()`.
    @return Value returned by `extract_metric_from_log()`.
    """
    file_glob = spec.get("file_glob", "logs/*.log")
    regex = spec.get("regex")
    if not regex:
        return None
    candidates = sorted(glob.glob(os.path.join(case_dir, file_glob), recursive=True))
    if not candidates:
        return None
    pattern = re.compile(regex)
    values = []
    for path in candidates:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    m = pattern.search(line)
                    if m:
                        try:
                            values.append(float(m.group(1)))
                        except Exception:
                            pass
        except OSError:
            continue
    return reduce_metric_values(values, spec.get("reduction", "last"))


def normalize_metric_spec(metric):
    """!
    @brief Normalize study metric definitions to a common dictionary form.
    @param[in] metric Argument passed to `normalize_metric_spec()`.
    @return Value returned by `normalize_metric_spec()`.
    """
    if isinstance(metric, str):
        if metric.lower() in {"msd", "msd_final"}:
            return {
                "name": "msd_final",
                "source": "statistics_csv",
                "file_glob": "**/*_msd.csv",
                "reduction": "last",
            }
        return {"name": metric, "source": "log_regex", "regex": metric}
    return dict(metric)

def aggregate_study_metrics(study_cfg: dict, cases: list, results_dir: str) -> str:
    """!
    @brief Collect metric values from generated case directories into one CSV.
    @param[in] study_cfg Argument passed to `aggregate_study_metrics()`.
    @param[in] cases Argument passed to `aggregate_study_metrics()`.
    @param[in] results_dir Argument passed to `aggregate_study_metrics()`.
    @return Value returned by `aggregate_study_metrics()`.
    """
    metrics = study_cfg.get("metrics", [])
    if not metrics:
        metrics = ["msd_final"]
    normalized_specs = [normalize_metric_spec(m) for m in metrics]

    rows = []
    for case in cases:
        row = {"case_id": case["case_id"]}
        for p_key, p_val in case["parameters"].items():
            row[p_key] = p_val
        for spec in normalized_specs:
            name = spec.get("name", "metric")
            source = str(spec.get("source", "")).lower()
            if source in {"statistics_csv", "csv"}:
                value = extract_metric_from_csv(case["run_dir"], spec)
            elif source in {"log_regex", "log"}:
                value = extract_metric_from_log(case["run_dir"], spec)
            else:
                value = None

            normalize_key = spec.get("normalize_by_parameter")
            if value is not None and normalize_key:
                denom = case.get("parameters", {}).get(normalize_key)
                try:
                    denom = float(denom)
                except Exception:
                    denom = None
                if denom not in (None, 0.0):
                    value = float(value) / denom
                else:
                    value = None

            row[name] = value
        rows.append(row)

    if not rows:
        return None

    all_keys = []
    seen = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                all_keys.append(k)

    os.makedirs(results_dir, exist_ok=True)
    out_csv = os.path.join(results_dir, "metrics_table.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[SUCCESS] Aggregated metrics table: {os.path.relpath(out_csv)}")
    return out_csv

def infer_plot_x_axis(study_cfg: dict, rows: list):
    """!
    @brief Infer x-axis key/values for study plots.
    @param[in] study_cfg Argument passed to `infer_plot_x_axis()`.
    @param[in] rows Argument passed to `infer_plot_x_axis()`.
    @return Value returned by `infer_plot_x_axis()`.
    """
    params = get_study_parameter_keys(study_cfg)
    if not params or not rows:
        return None, None

    study_type = study_cfg.get("study_type")
    if study_type == "grid_independence":
        has_im = "case.grid.programmatic_settings.im" in params
        has_jm = "case.grid.programmatic_settings.jm" in params
        has_km = "case.grid.programmatic_settings.km" in params
        if has_im and has_jm and has_km:
            xs = []
            for row in rows:
                try:
                    im = float(row["case.grid.programmatic_settings.im"])
                    jm = float(row["case.grid.programmatic_settings.jm"])
                    km = float(row["case.grid.programmatic_settings.km"])
                    xs.append((im * jm * km) ** (1.0 / 3.0))
                except Exception:
                    return None, None
            return "N^(1/3)", xs

    primary = params[0]
    xs = []
    for row in rows:
        try:
            xs.append(float(row[primary]))
        except Exception:
            return None, None
    return primary, xs

def generate_study_plots(study_cfg: dict, metrics_csv: str, plots_dir: str):
    """!
    @brief Generate metric-vs-parameter plots for completed studies.
    @param[in] study_cfg Argument passed to `generate_study_plots()`.
    @param[in] metrics_csv Argument passed to `generate_study_plots()`.
    @param[in] plots_dir Argument passed to `generate_study_plots()`.
    @return Value returned by `generate_study_plots()`.
    """
    plotting_cfg = study_cfg.get("plotting", {}) or {}
    if plotting_cfg.get("enabled", True) is False:
        print("[INFO] Plotting disabled by study.yml.")
        return []
    plt = optional_matplotlib_pyplot()
    if plt is None:
        print("[WARNING] matplotlib not available; skipping plot generation.")
        return []
    if not metrics_csv or not os.path.isfile(metrics_csv):
        return []

    with open(metrics_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        return []

    x_name, x_values = infer_plot_x_axis(study_cfg, rows)
    if not x_name or x_values is None:
        print("[WARNING] Could not infer numeric x-axis for plots; skipping.")
        return []

    metric_keys = []
    param_keys = get_study_parameter_keys(study_cfg)
    for key in rows[0].keys():
        if key in {"case_id"}:
            continue
        if key in param_keys:
            continue
        metric_keys.append(key)

    out_format = plotting_cfg.get("output_format", "png")
    os.makedirs(plots_dir, exist_ok=True)
    generated = []
    for metric in metric_keys:
        y_values = []
        ok = True
        for row in rows:
            try:
                y_values.append(float(row[metric]))
            except Exception:
                ok = False
                break
        if not ok:
            continue
        plt.figure(figsize=(7.0, 4.2))
        plt.plot(x_values, y_values, marker="o", linewidth=1.5)
        plt.xlabel(x_name)
        plt.ylabel(metric)
        plt.title(f"{metric} vs {x_name}")
        plt.grid(True, alpha=0.3)
        out_path = os.path.join(plots_dir, f"{metric}_vs_{x_name.replace('/', '_')}.{out_format}")
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()
        generated.append(out_path)
    if generated:
        print(f"[SUCCESS] Generated {len(generated)} plot(s) in {os.path.relpath(plots_dir)}")
    return generated


def _command_to_string(command_tokens: list) -> str:
    """!
    @brief Render a command list as a shell-safe display string.
    @param[in] command_tokens Argument passed to `_command_to_string()`.
    @return Value returned by `_command_to_string()`.
    """
    return " ".join(shlex.quote(str(tok)) for tok in command_tokens)


def _resolve_post_source_directory_preview(run_dir: str, monitor_cfg: dict, post_cfg: dict) -> str:
    """!
    @brief Resolve post source directory without side effects or stdout/stderr output.
    @param[in] run_dir Argument passed to `_resolve_post_source_directory_preview()`.
    @param[in] monitor_cfg Argument passed to `_resolve_post_source_directory_preview()`.
    @param[in] post_cfg Argument passed to `_resolve_post_source_directory_preview()`.
    @return Value returned by `_resolve_post_source_directory_preview()`.
    """
    solver_output_dir_rel = monitor_cfg.get('io', {}).get('directories', {}).get('output', 'output')
    solver_output_dir_abs = os.path.join(run_dir, solver_output_dir_rel)
    source_dir_template = get_post_source_directory_template(post_cfg)
    if source_dir_template == '<solver_output_dir>':
        return solver_output_dir_abs
    return os.path.abspath(os.path.join(run_dir, source_dir_template))


def build_run_dry_plan(args) -> dict:
    """!
    @brief Build a no-write execution plan for `run --dry-run`.
    @param[in] args Command-line style argument list supplied to the function.
    @return Value returned by `build_run_dry_plan()`.
    """
    plan = {
        "mode": "dry-run",
        "created_at": datetime.now().isoformat(),
        "warnings": [],
        "inputs": {},
        "stages": {},
        "artifacts": [],
    }

    if args.dry_run and args.no_submit:
        plan["warnings"].append("--dry-run takes precedence over --no-submit; no files will be written.")

    cluster_mode = bool(getattr(args, "cluster", None))
    cluster_cfg = None
    cluster_path = None
    solver_num_procs_effective = args.num_procs
    post_num_procs_effective = 1
    run_id = None
    run_dir = None
    solver_control_path = None
    loaded_case_cfg = None
    loaded_monitor_cfg = None
    resolved_restart_source_dir = None

    if cluster_mode:
        cluster_path = os.path.abspath(args.cluster)
        cluster_cfg = read_yaml_file(cluster_path)
        validate_cluster_config(cluster_cfg, cluster_path)
        scheduler_type = str(cluster_cfg.get("scheduler", {}).get("type", "slurm")).lower()
        if args.scheduler and args.scheduler.lower() != scheduler_type:
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key="scheduler.type",
                file_path=cluster_path,
                message=f"--scheduler={args.scheduler} does not match cluster.yml scheduler.type={scheduler_type}.",
            )
            sys.exit(1)
        if scheduler_type != "slurm":
            emit_structured_error(
                ERROR_CODE_CFG_INVALID_VALUE,
                key="scheduler.type",
                file_path=cluster_path,
                message=f"Unsupported scheduler '{scheduler_type}'. Only Slurm is supported in v1.",
            )
            sys.exit(1)
        cluster_tasks = get_cluster_total_tasks(cluster_cfg)
        if args.solve and args.num_procs not in (1, cluster_tasks):
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key="resources.ntasks_per_node",
                file_path=cluster_path,
                message=(
                    "--num-procs applies to the solver stage and must be 1 (auto) or "
                    f"exactly nodes*ntasks_per_node ({cluster_tasks}) in cluster mode."
                ),
            )
            sys.exit(1)
        if args.solve:
            solver_num_procs_effective = cluster_tasks
        plan["launch_mode"] = "slurm"
        plan["inputs"]["cluster"] = cluster_path
    else:
        if getattr(args, "scheduler", None):
            fail_cli_usage("--scheduler requires --cluster in this version.")
        plan["launch_mode"] = "local"

    # --- Guard: restart flags without --solve ---
    if not args.solve:
        if getattr(args, 'restart_from', None):
            print("[WARNING] --restart-from has no effect without --solve and will be ignored.", file=sys.stderr)
        if getattr(args, 'continue_run', False) and not args.post_process:
            print("[WARNING] --continue has no effect without --solve or --post-process and will be ignored.", file=sys.stderr)

    if args.solve:
        case_path = os.path.abspath(args.case)
        solver_path = os.path.abspath(args.solver)
        monitor_path = os.path.abspath(args.monitor)
        loaded_case_cfg = read_yaml_file(case_path)
        solver_cfg = read_yaml_file(solver_path)
        loaded_monitor_cfg = read_yaml_file(monitor_path)
        validate_solver_configs(loaded_case_cfg, solver_cfg, loaded_monitor_cfg, case_path, solver_path, monitor_path)

        continue_mode = getattr(args, 'continue_run', False)

        if continue_mode:
            # --continue reuses existing run directory
            if not args.run_dir:
                fail_cli_usage("--continue requires --run-dir.")
            run_dir = os.path.abspath(args.run_dir)
            if not os.path.isdir(run_dir):
                emit_structured_error(
                    ERROR_CODE_CFG_FILE_NOT_FOUND,
                    key="run-dir",
                    file_path=run_dir,
                    message="Specified run directory not found.",
                )
                sys.exit(1)
            run_id = os.path.basename(run_dir)
        else:
            case_name = os.path.splitext(os.path.basename(case_path))[0]
            timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            run_id = f"{case_name}_{timestamp}"
            run_dir = os.path.abspath(os.path.join("runs", run_id))

        try:
            resolved_restart_source_dir, is_continue = resolve_restart_source(
                args, loaded_case_cfg, solver_cfg, loaded_monitor_cfg, run_dir
            )
        except ValueError as e:
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key="restart",
                file_path=case_path,
                message=str(e),
            )
            sys.exit(1)

        config_dir = os.path.join(run_dir, "config")
        scheduler_dir = os.path.join(run_dir, "scheduler")
        logs_dir = os.path.join(run_dir, "logs")
        solver_control_path = os.path.join(config_dir, f"{run_id}.control")
        profile_path = os.path.join(config_dir, "profile.run")
        profiling_preview = resolve_profiling_config(loaded_monitor_cfg)

        plan["run_id_preview"] = run_id
        plan["run_dir_preview"] = run_dir
        plan["inputs"].update({"case": case_path, "solver": solver_path, "monitor": monitor_path})
        plan["artifacts"].extend(
            [
                run_dir,
                config_dir,
                logs_dir,
                os.path.join(run_dir, "output"),
                scheduler_dir,
                os.path.join(config_dir, "case.yml"),
                os.path.join(config_dir, "solver.yml"),
                os.path.join(config_dir, "monitor.yml"),
                solver_control_path,
                os.path.join(run_dir, "manifest.json"),
            ]
        )
        add_planned_grid_artifacts(plan, loaded_case_cfg, run_dir)
        add_planned_profile_artifacts(plan, loaded_case_cfg, run_dir)
        add_planned_initial_condition_artifacts(plan, loaded_case_cfg, solver_cfg, run_dir)
        if has_explicit_monitor_whitelist(loaded_monitor_cfg):
            plan["artifacts"].append(os.path.join(config_dir, "whitelist.run"))
        if profiling_preview["mode"] == "selected":
            plan["artifacts"].append(profile_path)
        solve_diagnostics = resolve_diagnostics_config(loaded_monitor_cfg, run_dir, "Solver")
        plan["artifacts"].extend(solve_diagnostics["artifacts"])
        if cluster_mode:
            plan["artifacts"].append(os.path.join(config_dir, "cluster.yml"))
            plan["artifacts"].append(os.path.join(scheduler_dir, "submission.json"))

        solver_exe = resolve_runtime_executable("simulator")
        solver_args = build_petsc_diagnostics_args(loaded_monitor_cfg, run_dir, "Solver") + ["-control_file", solver_control_path]
        if cluster_mode:
            solver_script = os.path.join(scheduler_dir, "solver.sbatch")
            solver_cmd = build_cluster_launch_command(
                cluster_cfg,
                solver_exe,
                solver_args,
                config_search_anchor=case_path,
                extra_search_anchors=[cluster_path],
            )
            plan["artifacts"].append(solver_script)
            plan["stages"]["solve"] = {
                "mode": "slurm",
                "script": solver_script,
                "num_procs_effective": solver_num_procs_effective,
                "launch_command": solver_cmd,
                "launch_command_string": _command_to_string(solver_cmd),
            }
        else:
            solver_cmd = build_local_launch_command(
                solver_exe,
                solver_args,
                solver_num_procs_effective,
                config_search_anchor=case_path,
            )
            solver_stream_log = os.path.join(scheduler_dir, f"{run_id}_solver.log")
            plan["artifacts"].append(solver_stream_log)
            plan["stages"]["solve"] = {
                "mode": "local",
                "num_procs_effective": solver_num_procs_effective,
                "stream_log": solver_stream_log,
                "launch_command": solver_cmd,
                "launch_command_string": _command_to_string(solver_cmd),
            }
        if resolved_restart_source_dir:
            plan["stages"]["solve"]["restart_source_directory"] = resolved_restart_source_dir
        if is_continue:
            plan["stages"]["solve"]["continue_mode"] = True

    if args.post_process:
        post_path = os.path.abspath(args.post)
        plan["inputs"]["post"] = post_path
        post_cfg = read_yaml_file(post_path)
        validate_post_config(post_cfg, post_path)

        if args.run_dir:
            run_dir = os.path.abspath(args.run_dir)
            if not os.path.isdir(run_dir):
                emit_structured_error(
                    ERROR_CODE_CFG_FILE_NOT_FOUND,
                    key="run-dir",
                    file_path=run_dir,
                    message="Specified run directory not found.",
                )
                sys.exit(1)
            run_id = os.path.basename(run_dir)
        elif not args.solve:
            fail_cli_usage("--post-process requires --run-dir when not used with --solve.")

        if args.run_dir:
            config_dir = os.path.join(run_dir, "config")
            case_path, monitor_path, solver_control_path = auto_identify_run_inputs(config_dir)
            if not all([case_path, monitor_path, solver_control_path]):
                emit_structured_error(
                    ERROR_CODE_CFG_MISSING_KEY,
                    key="run_dir.config",
                    file_path=config_dir,
                    message=(
                        "Could not auto-identify required run inputs "
                        "(case.yml/monitor.yml/*.control) in run config directory."
                    ),
                )
                sys.exit(1)
            loaded_case_cfg = read_yaml_file(case_path)
            loaded_monitor_cfg = read_yaml_file(monitor_path)
        else:
            config_dir = os.path.join(run_dir, "config")
            case_path = os.path.join(config_dir, "case.yml")
            monitor_path = os.path.join(config_dir, "monitor.yml")
            if solver_control_path is None:
                solver_control_path = os.path.join(config_dir, f"{run_id}.control")

        allow_source_frontier_scan = not args.solve
        post_plan = build_post_execution_plan(
            run_dir,
            run_id,
            loaded_case_cfg,
            loaded_monitor_cfg,
            post_cfg,
            continue_requested=getattr(args, 'continue_run', False),
            allow_source_frontier_scan=allow_source_frontier_scan,
        )

        post_recipe_path = os.path.join(config_dir, "post.run")
        output_dir_rel = post_cfg.get("io", {}).get("output_directory")
        output_prefix = post_cfg.get("io", {}).get("output_filename_prefix")
        if not output_dir_rel or not output_prefix:
            emit_structured_error(
                ERROR_CODE_CFG_MISSING_KEY,
                key="io.output_directory/io.output_filename_prefix",
                file_path=post_path,
                message="Missing required post IO keys.",
            )
            sys.exit(1)
        output_dir_abs = os.path.abspath(os.path.join(run_dir, output_dir_rel))
        statistics_output_paths = get_post_statistics_output_artifacts(post_cfg, run_dir, loaded_monitor_cfg)
        post_exe = resolve_runtime_executable("postprocessor")
        post_diagnostics = resolve_diagnostics_config(loaded_monitor_cfg, run_dir, "PostProcessor")
        plan["artifacts"].extend(post_diagnostics["artifacts"])
        post_args = build_petsc_diagnostics_args(loaded_monitor_cfg, run_dir, "PostProcessor") + [
            "-control_file",
            solver_control_path,
            "-postprocessing_config_file",
            post_recipe_path,
        ]
        plan["artifacts"].extend([
            post_recipe_path,
            output_dir_abs,
            post_plan["resume_state_path"],
            post_plan["lock_paths"]["wrapper_path"],
            post_plan["lock_paths"]["lock_file"],
            post_plan["lock_paths"]["metadata_file"],
        ])
        plan["artifacts"].extend(statistics_output_paths)

        stage_meta = {
            "source_data_directory": post_plan["source_data_directory"],
            "requested_start_step": post_plan["requested_start_step"],
            "requested_end_step": post_plan["requested_end_step"],
            "step_interval": post_plan["step_interval"],
            "resume_applied": bool(post_plan["continue_requested"] and post_plan["resume_recipe_match"]),
            "resume_recipe_match": post_plan["resume_recipe_match"],
            "resume_bootstrapped": post_plan["resume_bootstrapped"],
            "resume_match_source": post_plan["resume_match_source"],
            "completed_frontier_step": post_plan["completed_frontier_step"],
            "source_frontier_step": post_plan["source_frontier_step"],
            "source_frontier_diagnostic": post_plan["source_frontier_diagnostic"],
            "source_frontier_deferred": post_plan["source_frontier_deferred"],
            "effective_start_step": post_plan["effective_start_step"],
            "effective_end_step": post_plan["effective_end_step"],
            "skip_reason": post_plan["skip_reason"],
            "post_skipped_as_complete": post_plan["skip_reason"] == "already-complete-window",
            "recipe_fingerprint": post_plan["recipe_fingerprint"],
            "num_procs_effective": post_num_procs_effective,
        }

        if post_plan["skip_reason"] is None:
            if cluster_mode:
                scheduler_dir = os.path.join(run_dir, "scheduler")
                post_script = os.path.join(scheduler_dir, "post.sbatch")
                post_cluster_cfg = build_serial_post_cluster_config(cluster_cfg, post_num_procs_effective)
                raw_post_cmd = build_cluster_launch_command(
                    post_cluster_cfg,
                    post_exe,
                    post_args,
                    config_search_anchor=case_path,
                    extra_search_anchors=[cluster_path],
                    force_num_procs=post_num_procs_effective,
                )
                post_cmd, _ = build_post_locked_command(
                    run_dir,
                    post_plan["recipe_fingerprint"],
                    raw_post_cmd,
                    create_wrapper=False,
                )
                plan["artifacts"].append(post_script)
                stage_meta.update({
                    "mode": "slurm",
                    "script": post_script,
                    "launch_command": post_cmd,
                    "launch_command_string": _command_to_string(post_cmd),
                })
            else:
                raw_post_cmd = build_local_launch_command(
                    post_exe,
                    post_args,
                    post_num_procs_effective,
                    config_search_anchor=case_path,
                    allow_single_rank_launcher_override=True,
                    force_num_procs=post_num_procs_effective,
                )
                post_cmd, _ = build_post_locked_command(
                    run_dir,
                    post_plan["recipe_fingerprint"],
                    raw_post_cmd,
                    create_wrapper=False,
                )
                post_stream_log = os.path.join(run_dir, "scheduler", f"{run_id}_{output_prefix}.log")
                plan["artifacts"].append(post_stream_log)
                stage_meta.update({
                    "mode": "local",
                    "stream_log": post_stream_log,
                    "launch_command": post_cmd,
                    "launch_command_string": _command_to_string(post_cmd),
                })
        else:
            stage_meta.update({
                "mode": "slurm" if cluster_mode else "local",
                "launch_command": [],
                "launch_command_string": "",
            })

        plan["stages"]["post-process"] = stage_meta

    # Preserve insertion order while removing duplicates.
    deduped = []
    seen = set()
    for item in plan["artifacts"]:
        if item not in seen:
            seen.add(item)
            deduped.append(item)
    plan["artifacts"] = deduped
    if run_id and "run_id_preview" not in plan:
        plan["run_id_preview"] = run_id
    if run_dir and "run_dir_preview" not in plan:
        plan["run_dir_preview"] = run_dir
    plan["solver_num_procs_effective"] = solver_num_procs_effective
    plan["post_num_procs_effective"] = post_num_procs_effective
    plan["num_procs_effective"] = solver_num_procs_effective
    return plan


def add_planned_grid_artifacts(plan: dict, case_cfg: dict, run_dir: str) -> None:
    """!
    @brief Add grid-mode-specific staged artifacts to a dry-run plan.
    @param[in,out] plan Dry-run plan to update.
    @param[in] case_cfg Parsed case configuration.
    @param[in] run_dir Preview run directory for relative artifact resolution.
    """
    grid_cfg = case_cfg.get("grid", {})
    if not isinstance(grid_cfg, dict):
        return

    mode = grid_cfg.get("mode")
    config_dir = os.path.join(run_dir, "config")

    if mode == "file":
        plan["artifacts"].append(os.path.join(config_dir, "grid.run"))
        legacy_cfg = grid_cfg.get("legacy_conversion")
        if isinstance(legacy_cfg, dict) and legacy_cfg.get("enabled", True):
            output_file = legacy_cfg.get("output_file", os.path.join("config", "grid.converted.picgrid"))
            if isinstance(output_file, str) and output_file.strip():
                if not os.path.isabs(output_file):
                    output_file = os.path.abspath(os.path.join(run_dir, output_file))
                plan["artifacts"].append(output_file)
    elif mode == "grid_gen":
        generator = grid_cfg.get("generator", {})
        if not isinstance(generator, dict):
            return
        plan["artifacts"].append(os.path.join(config_dir, "grid.run"))
        for key, default in (
            ("output_file", os.path.join("config", "grid.generated.picgrid")),
            ("stats_file", None),
            ("vts_file", None),
        ):
            artifact_path = generator.get(key, default)
            if isinstance(artifact_path, str) and artifact_path.strip():
                if not os.path.isabs(artifact_path):
                    artifact_path = os.path.abspath(os.path.join(run_dir, artifact_path))
                plan["artifacts"].append(artifact_path)

def add_planned_profile_artifacts(plan: dict, case_cfg: dict, run_dir: str) -> None:
    """!
    @brief Add generated prescribed-flow profile artifacts to a dry-run plan.
    @param[in,out] plan Dry-run plan to update.
    @param[in] case_cfg Parsed case configuration.
    @param[in] run_dir Preview run directory.
    """
    try:
        prepared_blocks = validate_and_prepare_boundary_conditions(case_cfg)
    except ValueError:
        return
    config_dir = os.path.join(run_dir, "config")
    has_generated = False
    for block_idx, block in enumerate(prepared_blocks):
        for bc in block:
            if bc.get("handler") != "prescribed_flow":
                continue
            source = (bc.get("params") or {}).get("source", {})
            if source.get("type") not in {"generated", "field_slice"}:
                continue
            has_generated = True
            face_token = _face_artifact_token(bc["face"])
            suffix = "generated" if source.get("type") == "generated" else "sliced"
            default_output = os.path.join(
                "config", f"inlet_profile_block{block_idx}_{face_token}.{suffix}.picslice"
            )
            generated_path = _resolve_run_artifact_path(
                run_dir,
                source.get("output_file"),
                default_output,
                default_to_config_dir=True,
            )
            staged_path = os.path.join(config_dir, f"inlet_profile_block{block_idx}_{face_token}.picslice")
            plan["artifacts"].append(generated_path)
            plan["artifacts"].append(staged_path)
    if has_generated:
        plan["artifacts"].append(os.path.join(config_dir, "profile.info"))

def add_planned_initial_condition_artifacts(plan: dict, case_cfg: dict, solver_cfg: dict, run_dir: str) -> None:
    """!
    @brief Add authoritative file-backed initial-condition artifacts to a dry-run plan.
    @param[in,out] plan Dry-run plan receiving artifact paths.
    @param[in] case_cfg Parsed case configuration.
    @param[in] solver_cfg Parsed solver configuration.
    @param[in] run_dir Planned run directory.
    """
    source = normalize_eulerian_field_source(
        (solver_cfg.get("operation_mode", {}) or {}).get("eulerian_field_source", "solve")
    )
    start_step = int((case_cfg.get("run_control", {}) or {}).get("start_step", 0) or 0)
    if source != "solve" or start_step != 0:
        return
    try:
        resolved = resolve_initial_condition_config(
            (case_cfg.get("properties", {}) or {}).get("initial_conditions", {}),
            validate_and_prepare_boundary_conditions(case_cfg),
            U_ref=1.0,
        )
    except (KeyError, ValueError):
        return
    if resolved["kind"] not in {"file", "ic_gen"}:
        return
    config_dir = os.path.join(run_dir, "config")
    plan["artifacts"].append(
        os.path.join(config_dir, "initial_condition", f"{resolved['field_name']}00000_0.dat")
    )
    if resolved["kind"] == "ic_gen":
        plan["artifacts"].append(_resolve_run_artifact_path(
            run_dir, resolved.get("output_file"), os.path.join("config", "initial_condition.generated.dat"),
            default_to_config_dir=True,
        ))


def render_run_dry_plan(plan: dict, output_format: str = "text"):
    """!
    @brief Render dry-run plan in human or JSON format.
    @param[in] plan Argument passed to `render_run_dry_plan()`.
    @param[in] output_format Argument passed to `render_run_dry_plan()`.
    """
    if output_format == "json":
        print(json.dumps(plan, indent=2, sort_keys=True))
        return

    print("\n" + "=" * 60)
    print("                      DRY-RUN PLAN")
    print("=" * 60)
    print(f"  Launch mode    : {plan.get('launch_mode')}")
    print(f"  Created at     : {plan.get('created_at')}")
    if plan.get("run_id_preview"):
        print(f"  Run ID preview : {plan.get('run_id_preview')}")
    if plan.get("run_dir_preview"):
        print(f"  Run dir preview: {plan.get('run_dir_preview')}")
    print(f"  Solver MPI procs: {plan.get('solver_num_procs_effective')}")
    print(f"  Post MPI procs  : {plan.get('post_num_procs_effective')}")
    if plan.get("warnings"):
        print("  Warnings       :")
        for warning in plan["warnings"]:
            print(f"    - {warning}")

    if plan.get("inputs"):
        print("\n  Inputs:")
        for key, value in plan["inputs"].items():
            print(f"    - {key}: {value}")

    if plan.get("stages"):
        print("\n  Planned stage commands:")
        for stage, details in plan["stages"].items():
            print(f"    - {stage} ({details.get('mode')}):")
            if details.get('skip_reason'):
                print(f"      skipped: {details.get('skip_reason')}")
            else:
                print(f"      {details.get('launch_command_string')}")

    diagnostics_artifacts = [item for item in plan.get("artifacts", []) if "PETSc_" in os.path.basename(str(item)) or os.path.basename(str(item)) == "Runtime_Memory.log"]
    if diagnostics_artifacts:
        print("\n  Diagnostics artifacts:")
        for artifact in diagnostics_artifacts:
            print(f"    - {artifact}")

    print("\n  Planned artifacts (no files created in dry-run):")
    for artifact in plan.get("artifacts", []):
        print(f"    - {artifact}")
    print("=" * 60)


def validate_workflow(args):
    """!
    @brief Implements `picurv validate` without launching solver/post workflows.
    @param[in] args Command-line style argument list supplied to the function.
    """
    checked = []
    solver_group_selected = any([args.case, args.solver, args.monitor])
    any_group_selected = solver_group_selected or any([args.post, args.cluster, args.study])
    case_path = None
    cluster_path = None

    if not any_group_selected:
        fail_cli_usage(
            "validate requires at least one config group. Provide solver trio and/or --post/--cluster/--study.",
            hint="Example: picurv validate --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml",
        )

    if solver_group_selected and not all([args.case, args.solver, args.monitor]):
        fail_cli_usage("When solver validation is requested, --case, --solver, and --monitor are all required.")

    # --- Guard: restart flags without solver group ---
    restart_from = getattr(args, 'restart_from', None)
    continue_run = getattr(args, 'continue_run', False)
    run_dir_val = getattr(args, 'run_dir', None)
    if not solver_group_selected:
        if restart_from:
            print("[WARNING] --restart-from has no effect without --case/--solver/--monitor and will be ignored.", file=sys.stderr)
        if continue_run and not args.post:
            print("[WARNING] --continue has no effect without solver configs or --post and will be ignored.", file=sys.stderr)
    if continue_run and not run_dir_val:
        fail_cli_usage("--continue requires --run-dir.")

    if solver_group_selected:
        case_path = os.path.abspath(args.case)
        solver_path = os.path.abspath(args.solver)
        monitor_path = os.path.abspath(args.monitor)
        case_cfg = read_yaml_file(case_path)
        solver_cfg = read_yaml_file(solver_path)
        monitor_cfg = read_yaml_file(monitor_path)
        validate_solver_configs(case_cfg, solver_cfg, monitor_cfg, case_path, solver_path, monitor_path)
        checked.extend([case_path, solver_path, monitor_path])

        # Validate restart flags if provided
        if restart_from or continue_run:
            target_run_dir = os.path.abspath(run_dir_val) if run_dir_val else os.path.abspath("runs/_validate_dummy")
            try:
                resolve_restart_source(args, case_cfg, solver_cfg, monitor_cfg, target_run_dir)
                print("[SUCCESS] Restart source validation passed.")
            except ValueError as e:
                print(f"[ERROR] Restart validation failed: {e}", file=sys.stderr)
                sys.exit(1)

    post_cfg = None
    if args.post:
        post_path = os.path.abspath(args.post)
        post_cfg = read_yaml_file(post_path)
        validate_post_config(post_cfg, post_path)
        checked.append(post_path)

    cluster_cfg = None
    if args.cluster:
        cluster_path = os.path.abspath(args.cluster)
        cluster_cfg = read_yaml_file(cluster_path)
        validate_cluster_config(cluster_cfg, cluster_path)
        checked.append(cluster_path)

    try:
        runtime_execution_path, _ = load_runtime_execution_config(
            case_path,
            extra_search_anchors=[cluster_path] if cluster_path else None,
        )
    except ValueError as exc:
        emit_structured_error(
            ERROR_CODE_CFG_INVALID_VALUE,
            key="runtime_execution",
            file_path=case_path or cluster_path or os.getcwd(),
            message=str(exc),
        )
        sys.exit(1)
    if runtime_execution_path:
        checked.append(runtime_execution_path)

    study_cfg = None
    if args.study:
        study_path = os.path.abspath(args.study)
        study_cfg = read_yaml_file(study_path)
        validate_study_config(study_cfg, study_path)
        checked.append(study_path)

    if post_cfg is not None and run_dir_val:
        post_path = os.path.abspath(args.post)
        validate_run_dir = os.path.abspath(run_dir_val)
        if os.path.isdir(validate_run_dir):
            monitor_for_post = monitor_cfg if solver_group_selected else None
            if monitor_for_post is None:
                config_dir_candidate = os.path.join(validate_run_dir, "config")
                monitor_candidate = os.path.join(config_dir_candidate, "monitor.yml")
                if os.path.isfile(monitor_candidate):
                    monitor_for_post = read_yaml_file(monitor_candidate)
            if monitor_for_post is not None:
                resolved_source = _resolve_post_source_directory_preview(validate_run_dir, monitor_for_post, post_cfg)
                if os.path.isdir(resolved_source) and os.listdir(resolved_source):
                    print(f"[SUCCESS] Post-processor source data directory exists: {resolved_source}")
                else:
                    print(f"[WARNING] Post-processor source data directory is missing or empty: {resolved_source}", file=sys.stderr)

    if args.strict and post_cfg is not None:
        post_path = os.path.abspath(args.post)
        source_dir = post_cfg.get("source_data", {}).get("directory")
        if source_dir and source_dir != "<solver_output_dir>":
            resolved = resolve_path(post_path, source_dir)
            if not os.path.isdir(resolved):
                emit_structured_error(
                    ERROR_CODE_CFG_FILE_NOT_FOUND,
                    key="source_data.directory",
                    file_path=post_path,
                    message=f"strict mode: source_data.directory resolves to missing directory '{resolved}'.",
                )
                sys.exit(1)

    if args.strict and study_cfg is not None:
        study_path = os.path.abspath(args.study)
        base_cfgs = study_cfg.get("base_configs", {})
        if isinstance(base_cfgs, dict):
            base_case_path = resolve_path(study_path, base_cfgs.get("case"))
            base_solver_path = resolve_path(study_path, base_cfgs.get("solver"))
            base_monitor_path = resolve_path(study_path, base_cfgs.get("monitor"))
            base_post_path = resolve_path(study_path, base_cfgs.get("post"))
            if all([base_case_path, base_solver_path, base_monitor_path]):
                validate_solver_configs(
                    read_yaml_file(base_case_path),
                    read_yaml_file(base_solver_path),
                    read_yaml_file(base_monitor_path),
                    base_case_path,
                    base_solver_path,
                    base_monitor_path,
                )
            if base_post_path:
                validate_post_config(read_yaml_file(base_post_path), base_post_path)

    print(f"[SUCCESS] Validation completed for {len(checked)} file(s).")
    for path in checked:
        print(f"  - {path}")

def precompute_workflow(args):
    """!
    @brief Generate deterministic case artifacts without launching solver/post stages.
    @param[in] args Parsed precompute command arguments.
    """
    case_path = os.path.abspath(args.case)
    case_cfg = read_yaml_file(case_path)
    case_name = os.path.splitext(os.path.basename(case_path))[0]
    output_dir = args.output_dir or os.path.join("precomputed", case_name)
    output_dir = os.path.abspath(output_dir)
    config_dir = os.path.join(output_dir, "config")
    os.makedirs(config_dir, exist_ok=True)

    print(f"[INFO] Precomputing deterministic artifacts for case: {case_path}")
    print(f"[INFO] Output directory: {output_dir}")
    validate_and_prepare_boundary_conditions(case_cfg)

    artifacts = []
    grid_cfg = case_cfg.get("grid", {}) or {}
    grid_mode = grid_cfg.get("mode")
    scaling = (case_cfg.get("properties", {}) or {}).get("scaling", {}) or {}
    length_ref = float(scaling.get("length_ref", 1.0))
    expected_nblk = int((case_cfg.get("models", {}) or {}).get("domain", {}).get("blocks", 1))
    staged_grid = os.path.join(config_dir, "grid.run")
    if grid_mode == "grid_gen":
        print("[INFO] Precomputing grid via grid.gen...")
        generated_grid = run_grid_generator(case_path, output_dir, grid_cfg)
        artifacts.append(os.path.abspath(generated_grid))
        validate_and_nondimensionalize_picgrid(generated_grid, staged_grid, length_ref, expected_nblk=expected_nblk)
        artifacts.append(os.path.abspath(staged_grid))
    elif grid_mode == "file":
        source_grid = _resolve_case_relative_path(grid_cfg.get("source_file"), os.path.dirname(case_path))
        if isinstance(grid_cfg.get("legacy_conversion"), dict):
            source_grid = convert_legacy_grid_with_gridgen(case_path, output_dir, grid_cfg, source_grid)
        validate_and_nondimensionalize_picgrid(source_grid, staged_grid, length_ref, expected_nblk=expected_nblk)
        artifacts.append(os.path.abspath(staged_grid))
        print(f"[INFO] Precomputed validated file grid: {staged_grid}")
    elif grid_mode == "programmatic_c":
        print(f"[INFO] Grid mode '{grid_mode}' does not require precomputed grid generation.")
    else:
        raise ValueError(f"Unsupported grid.mode '{grid_mode}' for precompute.")

    profile_summaries = materialize_generated_prescribed_flow_profiles(output_dir, case_cfg, case_path)
    artifacts.extend(summary["path"] for summary in profile_summaries)
    if profile_summaries:
        artifacts.append(os.path.join(config_dir, "profile.info"))

    initial_condition = None
    resolved_ic = resolve_initial_condition_config(
        (case_cfg.get("properties", {}) or {}).get("initial_conditions", {}),
        validate_and_prepare_boundary_conditions(case_cfg),
        U_ref=float((case_cfg.get("properties", {}).get("scaling", {}) or {}).get("velocity_ref", 1.0)),
    )
    if resolved_ic["kind"] == "ic_gen":
        if grid_mode == "programmatic_c":
            try:
                summary = generate_picgrid_from_programmatic_settings(
                    grid_cfg.get('programmatic_settings', {}), staged_grid, length_ref
                )
                artifacts.append(os.path.abspath(staged_grid))
                print(
                    f"[INFO] Materialized programmatic grid.run for ic_gen: {staged_grid} "
                    f"(nblk={summary['nblk']}, total_nodes={summary['total_nodes']})"
                )
            except Exception as e:
                raise RuntimeError(f"Failed to generate grid.run for ic_gen: {e}") from e
        initial_condition = stage_initial_condition_file(output_dir, case_path, resolved_ic)
        artifacts.extend([initial_condition["source"], initial_condition["staged"]])
    elif resolved_ic["kind"] == "file":
        print("[INFO] File initial condition does not require generated precompute output.")
    else:
        print(f"[INFO] Built-in initial-condition generator '{resolved_ic['label']}' runs in the C solver.")

    manifest = {
        "case": case_path,
        "output_dir": output_dir,
        "grid_mode": grid_mode,
        "artifacts": artifacts,
        "profiles": profile_summaries,
        "initial_condition": initial_condition,
    }
    manifest_path = os.path.join(config_dir, "precompute.manifest.json")
    write_json_file(manifest_path, manifest)
    print(f"[SUCCESS] Wrote precompute manifest: {os.path.relpath(manifest_path)}")
    print(f"[SUCCESS] Precompute completed with {len(artifacts)} artifact(s).")

def run_workflow(args):
    """!
    @brief Main orchestrator for the 'run' command (local and Slurm modes).
    @param[in] args Command-line style argument list supplied to the function.
    """
    if getattr(args, "dry_run", False):
        plan = build_run_dry_plan(args)
        render_run_dry_plan(plan, output_format=getattr(args, "output_format", "text"))
        return

    run_dir = None
    run_id = None
    output_dir_abs = None
    statistics_output_paths = []
    workflow_start = time.time()
    stages_completed = []
    configs = None
    submission_meta = {"launch_mode": "local", "no_submit": bool(args.no_submit), "stages": {}}

    cluster_mode = bool(getattr(args, "cluster", None))
    cluster_cfg = None
    cluster_path = None
    solver_num_procs_effective = args.num_procs
    post_num_procs_effective = 1

    if cluster_mode:
        cluster_path = os.path.abspath(args.cluster)
        cluster_cfg = read_yaml_file(cluster_path)
        validate_cluster_config(cluster_cfg, cluster_path)
        scheduler_type = str(cluster_cfg.get("scheduler", {}).get("type", "slurm")).lower()
        if args.scheduler and args.scheduler.lower() != scheduler_type:
            print(
                f"[FATAL] --scheduler={args.scheduler} does not match cluster.yml scheduler.type={scheduler_type}.",
                file=sys.stderr
            )
            sys.exit(1)
        if scheduler_type != "slurm":
            print(f"[FATAL] Unsupported scheduler '{scheduler_type}'. Only Slurm is supported in v1.", file=sys.stderr)
            sys.exit(1)
        cluster_tasks = get_cluster_total_tasks(cluster_cfg)
        if args.solve and args.num_procs not in (1, cluster_tasks):
            print(
                "[FATAL] In cluster mode, --num-procs applies to the solver stage and must be "
                f"1 (auto) or exactly nodes*ntasks_per_node ({cluster_tasks}).",
                file=sys.stderr
            )
            sys.exit(1)
        if args.solve:
            solver_num_procs_effective = cluster_tasks
        submission_meta["launch_mode"] = "slurm"
        submission_meta["cluster_config"] = cluster_path
        submission_meta["no_submit"] = bool(args.no_submit)
        if args.solve:
            print(
                f"[INFO] Cluster mode enabled (Slurm). Solver uses {solver_num_procs_effective} MPI tasks "
                f"from cluster.yml; post stage defaults to {post_num_procs_effective} task."
            )
        else:
            print(f"[INFO] Cluster mode enabled (Slurm). Post stage defaults to {post_num_procs_effective} task.")
    elif getattr(args, "scheduler", None):
        print("[FATAL] --scheduler requires --cluster in this version.", file=sys.stderr)
        sys.exit(1)

    # --- Guard: restart flags without --solve ---
    if not args.solve:
        if getattr(args, 'restart_from', None):
            print("[WARNING] --restart-from has no effect without --solve and will be ignored.", file=sys.stderr)
        if getattr(args, 'continue_run', False) and not args.post_process:
            print("[WARNING] --continue has no effect without --solve or --post-process and will be ignored.", file=sys.stderr)

    # --- Stage 1: Solver (if requested) ---
    if args.solve:
        walltime_guard_policy = resolve_walltime_guard_policy(cluster_cfg) if cluster_mode else None
        configs = {
            'case': read_yaml_file(args.case), 'case_path': os.path.abspath(args.case),
            'solver': read_yaml_file(args.solver), 'solver_path': os.path.abspath(args.solver),
            'monitor': read_yaml_file(args.monitor), 'monitor_path': os.path.abspath(args.monitor),
            'walltime_guard_policy': walltime_guard_policy,
        }

        print("\n[INFO] Validating configuration files...")
        validate_solver_configs(
            configs['case'], configs['solver'], configs['monitor'],
            args.case, args.solver, args.monitor
        )
        print("[SUCCESS] All configuration files passed validation.\n")

        continue_mode = getattr(args, 'continue_run', False)

        if continue_mode:
            if not args.run_dir:
                fail_cli_usage("--continue requires --run-dir.")
            run_dir = os.path.abspath(args.run_dir)
            if not os.path.isdir(run_dir):
                emit_structured_error(
                    ERROR_CODE_CFG_FILE_NOT_FOUND,
                    key="run-dir",
                    file_path=run_dir,
                    message="Specified run directory not found.",
                )
                sys.exit(1)
            run_id = os.path.basename(run_dir)
        else:
            case_name = os.path.splitext(os.path.basename(args.case))[0]
            timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            run_id = f"{case_name}_{timestamp}"
            run_dir = os.path.abspath(os.path.join("runs", run_id))

        try:
            resolved_restart_source_dir, is_continue = resolve_restart_source(
                args, configs["case"], configs["solver"], configs["monitor"], run_dir
            )
        except ValueError as e:
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key="restart",
                file_path=args.case,
                message=str(e),
            )
            sys.exit(1)

        config_dir = os.path.join(run_dir, "config")
        if not continue_mode:
            for d in [config_dir, os.path.join(run_dir, "scheduler")]:
                os.makedirs(d, exist_ok=True)
        else:
            os.makedirs(config_dir, exist_ok=True)
        if continue_mode:
            print(f"[INFO] Continuing in existing run directory: {os.path.relpath(run_dir)}")
        else:
            print(f"[INFO] Created new self-contained run directory: {os.path.relpath(run_dir)}")

        shutil.copy(args.case, os.path.join(config_dir, "case.yml"))
        shutil.copy(args.solver, os.path.join(config_dir, "solver.yml"))
        shutil.copy(args.monitor, os.path.join(config_dir, "monitor.yml"))
        if cluster_mode:
            shutil.copy(cluster_path, os.path.join(config_dir, "cluster.yml"))

        print("\n" + "="*25 + " SOLVER STAGE " + "="*25)
        source_files = {'Case': args.case, 'Solver': args.solver, 'Monitor': args.monitor}
        monitor_files = prepare_monitor_files(run_dir, run_id, configs['monitor'], source_files)
        if resolved_restart_source_dir:
            print(f"[INFO] Restart source: {resolved_restart_source_dir}")
        if is_continue:
            print("[INFO] Continue mode: logs will be appended, not overwritten.")
        control_file = generate_solver_control_file(
            run_dir,
            run_id,
            configs,
            solver_num_procs_effective,
            monitor_files,
            restart_source_dir=resolved_restart_source_dir,
            continue_mode=is_continue,
        )

        solver_exe = resolve_runtime_executable("simulator")
        solver_args = build_petsc_diagnostics_args(configs["monitor"], run_dir, "Solver") + ["-control_file", control_file]
        if cluster_mode:
            scheduler_dir = os.path.join(run_dir, "scheduler")
            solver_script = os.path.join(scheduler_dir, "solver.sbatch")
            solver_log = os.path.join(scheduler_dir, "solver_%j.out")
            solver_err = os.path.join(scheduler_dir, "solver_%j.err")
            solver_cmd = build_cluster_launch_command(
                cluster_cfg,
                solver_exe,
                solver_args,
                config_search_anchor=args.case,
                extra_search_anchors=[cluster_path],
            )
            render_slurm_script(
                solver_script,
                f"{run_id}_solve",
                cluster_cfg,
                solver_cmd,
                run_dir,
                solver_log,
                solver_err,
                env_vars={"LOG_LEVEL": configs['monitor'].get('logging', {}).get('verbosity', 'INFO').upper()},
                shell_env_vars=build_walltime_guard_exports(cluster_cfg),
            )
            submission_meta["stages"]["solve"] = {
                "script": solver_script,
                "submitted": False,
                "num_procs_effective": solver_num_procs_effective,
            }
            print(f"[SUCCESS] Generated solver Slurm script: {os.path.relpath(solver_script)}")
            if not args.no_submit:
                submit_info = submit_sbatch(solver_script)
                submission_meta["stages"]["solve"].update(submit_info)
                submission_meta["stages"]["solve"]["submitted"] = True
                print(f"[SUCCESS] Submitted solver job: {submit_info['job_id']}")
            stages_completed.append('solve')
        else:
            command = build_local_launch_command(
                solver_exe,
                solver_args,
                solver_num_procs_effective,
                config_search_anchor=configs["case_path"],
            )
            solver_log = os.path.join("scheduler", f"{run_id}_solver.log")
            submission_meta["stages"]["solve"] = {
                "command": command,
                "command_string": format_command_for_display(command),
                "log_file": solver_log,
                "submitted": False,
                "num_procs_effective": solver_num_procs_effective,
            }
            if args.no_submit:
                print(f"[SUCCESS] Staged local solver command: {solver_log}")
            else:
                execute_command(command, run_dir, solver_log, configs['monitor'])
                submission_meta["stages"]["solve"]["submitted"] = True
                submission_meta["stages"]["solve"]["executed"] = True
                submission_meta["stages"]["solve"]["completed_at"] = datetime.now().isoformat()
            stages_completed.append('solve')

    # --- Stage 2: Post-Processing (if requested) ---
    if args.post_process:
        if args.run_dir:
            run_dir = os.path.abspath(args.run_dir)
            if not os.path.isdir(run_dir):
                print(f"[FATAL] Specified run directory not found: {run_dir}", file=sys.stderr)
                sys.exit(1)
            print(f"[INFO] Operating on existing run directory: {os.path.relpath(run_dir)}")
            run_id = os.path.basename(run_dir)
        elif not args.solve:
            print("[FATAL] --post-process requires --run-dir when not used with --solve.", file=sys.stderr)
            sys.exit(1)

        print("\n" + "="*20 + " POST-PROCESSING STAGE " + "="*20)
        config_dir = os.path.join(run_dir, "config")
        case_path, monitor_path, solver_control_path = auto_identify_run_inputs(config_dir)

        if not all([case_path, monitor_path, solver_control_path]):
            print(f"[FATAL] Could not automatically identify required config files in {config_dir}", file=sys.stderr)
            if not case_path:
                print("         - No 'case' file found (expected 'models' + 'boundary_conditions').", file=sys.stderr)
            if not monitor_path:
                print("         - No 'monitor' file found (expected 'io' + 'logging').", file=sys.stderr)
            if not solver_control_path:
                print("         - No '.control' file found.", file=sys.stderr)
            sys.exit(1)

        print(f"[INFO] Auto-identified Case file:    {os.path.basename(case_path)}")
        print(f"[INFO] Auto-identified Monitor file: {os.path.basename(monitor_path)}")

        case_cfg = read_yaml_file(case_path)
        monitor_cfg = read_yaml_file(monitor_path)
        post_cfg = read_yaml_file(args.post)

        print("[INFO] Validating post-processing configuration...")
        validate_post_config(post_cfg, args.post)
        print("[SUCCESS] Post-processing configuration passed validation.\n")

        solver_sources_deferred = bool(args.solve and (cluster_mode or args.no_submit))
        allow_source_frontier_scan = not solver_sources_deferred
        post_plan = build_post_execution_plan(
            run_dir,
            run_id,
            case_cfg,
            monitor_cfg,
            post_cfg,
            continue_requested=getattr(args, 'continue_run', False),
            allow_source_frontier_scan=allow_source_frontier_scan,
        )

        source_template = get_post_source_directory_template(post_cfg)
        if source_template == '<solver_output_dir>':
            print(f"[INFO] Post-processor source data: {os.path.relpath(post_plan['source_data_directory'])}")
        else:
            print(f"[INFO] Post-processor source data (user-defined): {os.path.relpath(post_plan['source_data_directory'])}")

        if getattr(args, 'continue_run', False):
            if post_plan['resume_recipe_match']:
                print(f"[INFO] Post resume recipe match: yes ({post_plan['resume_match_source']}).")
            else:
                print("[INFO] Post resume recipe match: no. Using the configured start_step for this recipe.")
        if post_plan['completed_frontier_step'] is not None:
            print(f"[INFO] Completed post frontier: step {post_plan['completed_frontier_step']}")
        else:
            print("[INFO] Completed post frontier: none")
        if post_plan['source_frontier_deferred']:
            print("[INFO] Source availability frontier: deferred because the solver stage will populate the requested window before post starts.")
        elif post_plan['source_frontier_step'] is not None:
            print(f"[INFO] Current source availability frontier: step {post_plan['source_frontier_step']}")
        else:
            print("[INFO] Current source availability frontier: none")

        persist_post_resume_state(run_dir, post_plan, last_successful_requested_end_step=post_plan['completed_frontier_step'])

        if post_plan['skip_reason'] == 'already-complete-window':
            print("[INFO] Requested post window is already complete; skipping postprocessor launch.")
            persist_post_resume_state(run_dir, post_plan, last_successful_requested_end_step=post_plan['requested_end_step'])
        elif post_plan['skip_reason'] == 'already-caught-up-to-current-source-frontier':
            print("[INFO] Post outputs are already caught up to the current fully available source frontier; nothing new to launch right now.")
            diagnostic = post_plan.get('source_frontier_diagnostic') or {}
            first_incomplete = diagnostic.get('first_incomplete_step')
            if first_incomplete is not None:
                print(f"[INFO] First incomplete requested source step: {first_incomplete}")
                print(
                    "[INFO] Closest complete source steps: "
                    f"near start={_format_optional_step(diagnostic.get('closest_complete_step_to_start'))}, "
                    f"near end={_format_optional_step(diagnostic.get('closest_complete_step_to_end'))}"
                )
        elif post_plan['skip_reason'] == 'nothing-available-yet':
            diagnostic = post_plan.get('source_frontier_diagnostic') or {}
            first_incomplete = diagnostic.get('first_incomplete_step')
            if first_incomplete is not None:
                print(
                    f"[INFO] First requested source step {first_incomplete} is incomplete; "
                    "skipping postprocessor launch for now."
                )
                print(
                    "[INFO] Closest complete source steps: "
                    f"near start={_format_optional_step(diagnostic.get('closest_complete_step_to_start'))}, "
                    f"near end={_format_optional_step(diagnostic.get('closest_complete_step_to_end'))}"
                )
                missing_files = diagnostic.get('missing_files_for_first_incomplete_step') or []
                if missing_files:
                    print(f"[INFO] Missing files for step {first_incomplete}: {', '.join(missing_files[:4])}")
            else:
                print("[INFO] No fully available source steps exist yet in the requested window; skipping postprocessor launch for now.")
        else:
            print(
                f"[INFO] Effective post window: {post_plan['effective_start_step']}..{post_plan['effective_end_step']} "
                f"(stride {post_plan['step_interval']})"
            )

            post_effective_cfg = post_plan['effective_post_cfg']
            post_io_cfg = post_effective_cfg.get('io', {})
            try:
                output_dir_rel = post_io_cfg['output_directory']
                output_prefix = post_io_cfg['output_filename_prefix']
            except KeyError as e:
                print(f"[FATAL] Missing required key '{e.args[0]}' in the 'io' section of {args.post}", file=sys.stderr)
                sys.exit(1)

            output_dir_abs = os.path.abspath(os.path.join(run_dir, output_dir_rel))
            os.makedirs(output_dir_abs, exist_ok=True)
            print(f"[INFO] Post-processor output directory: {os.path.relpath(output_dir_abs)}")
            statistics_output_paths = get_post_statistics_output_artifacts(post_effective_cfg, run_dir, monitor_cfg)
            for stats_path in statistics_output_paths:
                print(f"[INFO] Statistics CSV output: {os.path.relpath(stats_path)}")

            source_files_post = {'Case': case_path, 'Post-Profile': args.post}
            post_recipe_file = generate_post_recipe_file(run_dir, run_id, post_effective_cfg, source_files_post, monitor_cfg)

            post_exe = resolve_runtime_executable("postprocessor")
            post_args = build_petsc_diagnostics_args(monitor_cfg, run_dir, "PostProcessor") + [
                "-control_file",
                solver_control_path,
                "-postprocessing_config_file",
                post_recipe_file,
            ]
            if cluster_mode:
                scheduler_dir = os.path.join(run_dir, "scheduler")
                os.makedirs(scheduler_dir, exist_ok=True)
                post_script = os.path.join(scheduler_dir, "post.sbatch")
                post_log = os.path.join(scheduler_dir, "post_%j.out")
                post_err = os.path.join(scheduler_dir, "post_%j.err")
                post_cluster_cfg = build_serial_post_cluster_config(cluster_cfg, post_num_procs_effective)
                raw_post_cmd = build_cluster_launch_command(
                    post_cluster_cfg,
                    post_exe,
                    post_args,
                    config_search_anchor=case_path,
                    extra_search_anchors=[cluster_path],
                    force_num_procs=post_num_procs_effective,
                )
                post_cmd, _ = build_post_locked_command(
                    run_dir,
                    post_plan['recipe_fingerprint'],
                    raw_post_cmd,
                    create_wrapper=True,
                )
                render_slurm_script(
                    post_script,
                    f"{run_id}_post",
                    post_cluster_cfg,
                    post_cmd,
                    run_dir,
                    post_log,
                    post_err,
                    env_vars={"LOG_LEVEL": monitor_cfg.get('logging', {}).get('verbosity', 'INFO').upper()},
                )
                submission_meta["stages"]["post-process"] = {
                    "script": post_script,
                    "submitted": False,
                    "num_procs_effective": post_num_procs_effective,
                    "resume_recipe_match": post_plan['resume_recipe_match'],
                    "resume_bootstrapped": post_plan['resume_bootstrapped'],
                    "resume_match_source": post_plan['resume_match_source'],
                    "effective_start_step": post_plan['effective_start_step'],
                    "effective_end_step": post_plan['effective_end_step'],
                    "completed_frontier_step": post_plan['completed_frontier_step'],
                    "source_frontier_step": post_plan['source_frontier_step'],
                    "source_frontier_deferred": post_plan['source_frontier_deferred'],
                    "recipe_fingerprint": post_plan['recipe_fingerprint'],
                }
                print(f"[SUCCESS] Generated post Slurm script: {os.path.relpath(post_script)}")

                if not args.no_submit:
                    dependency_job = None
                    if args.solve:
                        dependency_job = submission_meta.get("stages", {}).get("solve", {}).get("job_id")
                    submit_info = submit_sbatch(post_script, dependency=dependency_job)
                    submission_meta["stages"]["post-process"].update(submit_info)
                    submission_meta["stages"]["post-process"]["submitted"] = True
                    if dependency_job:
                        submission_meta["stages"]["post-process"]["dependency"] = f"afterok:{dependency_job}"
                    print(f"[SUCCESS] Submitted post job: {submit_info['job_id']}")
                stages_completed.append('post-process')
            else:
                raw_command = build_local_launch_command(
                    post_exe,
                    post_args,
                    post_num_procs_effective,
                    config_search_anchor=case_path,
                    allow_single_rank_launcher_override=True,
                    force_num_procs=post_num_procs_effective,
                )
                command, _ = build_post_locked_command(
                    run_dir,
                    post_plan['recipe_fingerprint'],
                    raw_command,
                    create_wrapper=True,
                )
                post_log = os.path.join("scheduler", f"{run_id}_{output_prefix}.log")
                submission_meta["stages"]["post-process"] = {
                    "command": command,
                    "command_string": format_command_for_display(command),
                    "log_file": post_log,
                    "submitted": False,
                    "num_procs_effective": post_num_procs_effective,
                    "resume_recipe_match": post_plan['resume_recipe_match'],
                    "resume_bootstrapped": post_plan['resume_bootstrapped'],
                    "resume_match_source": post_plan['resume_match_source'],
                    "effective_start_step": post_plan['effective_start_step'],
                    "effective_end_step": post_plan['effective_end_step'],
                    "completed_frontier_step": post_plan['completed_frontier_step'],
                    "source_frontier_step": post_plan['source_frontier_step'],
                    "source_frontier_deferred": post_plan['source_frontier_deferred'],
                    "recipe_fingerprint": post_plan['recipe_fingerprint'],
                }
                if args.no_submit:
                    print(f"[SUCCESS] Staged local post command: {post_log}")
                else:
                    execute_command(command, run_dir, post_log, monitor_cfg)
                    persist_post_resume_state(run_dir, post_plan, last_successful_requested_end_step=post_plan['effective_end_step'])
                    submission_meta["stages"]["post-process"]["submitted"] = True
                    submission_meta["stages"]["post-process"]["executed"] = True
                    submission_meta["stages"]["post-process"]["completed_at"] = datetime.now().isoformat()
                stages_completed.append('post-process')

    if run_dir:
        manifest = {
            "run_id": run_id,
            "created_at": datetime.now().isoformat(),
            "launch_mode": "slurm" if cluster_mode else "local",
            "git_commit": get_git_commit(),
            "num_procs": solver_num_procs_effective,
            "solver_num_procs": solver_num_procs_effective,
            "post_num_procs": post_num_procs_effective,
            "stages_requested": {"solve": bool(args.solve), "post_process": bool(args.post_process)},
            "stages_completed_or_submitted": stages_completed,
            "inputs": {},
        }
        if args.solve:
            manifest["inputs"]["case"] = os.path.abspath(args.case)
            manifest["inputs"]["solver"] = os.path.abspath(args.solver)
            manifest["inputs"]["monitor"] = os.path.abspath(args.monitor)
        if args.post_process:
            manifest["inputs"]["post"] = os.path.abspath(args.post)
        if cluster_mode:
            manifest["inputs"]["cluster"] = cluster_path
        if submission_meta.get("stages"):
            write_json_file(os.path.join(run_dir, "scheduler", "submission.json"), submission_meta)
        write_json_file(os.path.join(run_dir, "manifest.json"), manifest)

    if stages_completed:
        elapsed = time.time() - workflow_start
        mins, secs = divmod(int(elapsed), 60)
        hrs, mins = divmod(mins, 60)
        if hrs > 0:
            time_str = f"{hrs}h {mins}m {secs}s"
        elif mins > 0:
            time_str = f"{mins}m {secs}s"
        else:
            time_str = f"{secs}s"

        print("\n" + "=" * 60)
        print("                       RUN SUMMARY")
        print("=" * 60)
        print(f"  Run ID         : {run_id}")
        print(f"  Run directory  : {os.path.relpath(run_dir)}")
        print(f"  Wall-clock     : {time_str}")
        print(f"  Stages         : {', '.join(stages_completed)}")
        print(f"  Launch mode    : {'slurm' if cluster_mode else 'local'}")
        if args.solve:
            print(f"  Solver MPI procs: {solver_num_procs_effective}")
        if args.post_process:
            print(f"  Post MPI procs  : {post_num_procs_effective}")
        if args.solve and configs:
            total_steps = configs['case'].get('run_control', {}).get('total_steps', '?')
            result_dir = os.path.join(run_dir, configs['monitor'].get('io', {}).get('directories', {}).get('output', 'output'))
            print(f"  Steps run      : {total_steps}")
            print(f"  Solver output  : {os.path.relpath(result_dir)}")
        if 'post-process' in stages_completed and output_dir_abs:
            print(f"  Post output    : {os.path.relpath(output_dir_abs)}")
        for stats_path in statistics_output_paths:
            print(f"  Stats output   : {os.path.relpath(stats_path)}")
        print(f"  Logs           : {os.path.relpath(os.path.join(run_dir, 'logs'))}")
        if cluster_mode or submission_meta.get("stages"):
            submission_file = os.path.join(run_dir, "scheduler", "submission.json")
            print(f"  Submission meta: {os.path.relpath(submission_file)}")
        print("=" * 60)


def parse_case_index_tsv(tsv_path: str) -> list:
    """!
    @brief Parse a case_index.tsv file back into a list of case entry dicts.
    @param[in] tsv_path Path to the case_index.tsv file.
    @return List of dicts with keys: index, case_id, run_dir, control_file,
            post_recipe_file, log_level, post_prefix.
    """
    entries = []
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            entries.append({
                "index": int(parts[0]),
                "case_id": parts[1],
                "run_dir": parts[2],
                "control_file": parts[3],
                "post_recipe_file": parts[4],
                "log_level": parts[5],
                "post_prefix": parts[6],
            })
    return entries


def sweep_workflow(args):
    """!
    @brief Study/sweep orchestration using Slurm job arrays.
    @param[in] args Command-line style argument list supplied to the function.
    """
    study_path = os.path.abspath(args.study)
    cluster_path = os.path.abspath(args.cluster)

    study_cfg = read_yaml_file(study_path)
    cluster_cfg = read_yaml_file(cluster_path)
    validate_study_config(study_cfg, study_path)
    validate_cluster_config(cluster_cfg, cluster_path)

    study_name = os.path.splitext(os.path.basename(study_path))[0]
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    study_id = f"{study_name}_{timestamp}"
    study_dir = os.path.abspath(os.path.join("studies", study_id))
    cases_dir = os.path.join(study_dir, "cases")
    scheduler_dir = os.path.join(study_dir, "scheduler")
    results_dir = os.path.join(study_dir, "results")
    for path in [cases_dir, scheduler_dir, results_dir]:
        os.makedirs(path, exist_ok=True)

    print(f"[INFO] Creating study directory: {os.path.relpath(study_dir)}")
    shutil.copy(study_path, os.path.join(study_dir, "study.yml"))
    shutil.copy(cluster_path, os.path.join(study_dir, "cluster.yml"))

    base_cfgs = study_cfg["base_configs"]
    base_paths = {k: resolve_path(study_path, v) for k, v in base_cfgs.items()}
    base_case = read_yaml_file(base_paths["case"])
    base_solver = read_yaml_file(base_paths["solver"])
    base_monitor = read_yaml_file(base_paths["monitor"])
    base_post = read_yaml_file(base_paths["post"])
    validate_solver_configs(base_case, base_solver, base_monitor, base_paths["case"], base_paths["solver"], base_paths["monitor"])
    validate_post_config(base_post, base_paths["post"])

    combinations = expand_study_parameter_combinations(study_cfg)
    if not combinations:
        print("[FATAL] Study parameter matrix expanded to zero cases.", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] Expanded sweep matrix to {len(combinations)} case(s).")

    cluster_tasks = get_cluster_total_tasks(cluster_cfg)
    case_entries = []
    case_index_file = os.path.join(scheduler_dir, "case_index.tsv")

    for idx, combo in enumerate(combinations):
        case_id = f"case_{idx:04d}"
        run_dir = os.path.join(cases_dir, case_id)
        config_dir = os.path.join(run_dir, "config")
        os.makedirs(config_dir, exist_ok=True)
        os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(run_dir, "output"), exist_ok=True)

        case_cfg = copy.deepcopy(base_case)
        solver_cfg = copy.deepcopy(base_solver)
        monitor_cfg = copy.deepcopy(base_monitor)
        post_cfg = copy.deepcopy(base_post)
        target_map = {"case": case_cfg, "solver": solver_cfg, "monitor": monitor_cfg, "post": post_cfg}
        for full_key, value in combo.items():
            root, nested = full_key.split(".", 1)
            _deep_set(target_map[root], nested, value)

        # Preserve file-based/grid-gen workflows when study cases are materialized
        # into new directories by rewriting external paths as absolute.
        absolutize_case_external_paths(case_cfg, base_paths["case"])

        case_path = os.path.join(config_dir, "case.yml")
        solver_path = os.path.join(config_dir, "solver.yml")
        monitor_path = os.path.join(config_dir, "monitor.yml")
        post_path = os.path.join(config_dir, "post.yml")
        write_yaml_file(case_path, case_cfg)
        write_yaml_file(solver_path, solver_cfg)
        write_yaml_file(monitor_path, monitor_cfg)
        write_yaml_file(post_path, post_cfg)

        validate_solver_configs(case_cfg, solver_cfg, monitor_cfg, case_path, solver_path, monitor_path)
        validate_post_config(post_cfg, post_path)

        source_files = {'Case': case_path, 'Solver': solver_path, 'Monitor': monitor_path}
        monitor_files = prepare_monitor_files(run_dir, case_id, monitor_cfg, source_files)
        configs = {
            "case": case_cfg, "case_path": case_path,
            "solver": solver_cfg, "solver_path": solver_path,
            "monitor": monitor_cfg, "monitor_path": monitor_path,
            "walltime_guard_policy": resolve_walltime_guard_policy(cluster_cfg),
        }
        control_file = generate_solver_control_file(run_dir, case_id, configs, cluster_tasks, monitor_files)

        source_dir = resolve_post_source_directory(run_dir, monitor_cfg, post_cfg, strict=False)
        if not isinstance(post_cfg.get('source_data'), dict):
            post_cfg['source_data'] = {}
        post_cfg['source_data']['directory'] = source_dir
        output_prefix = post_cfg.get("io", {}).get("output_filename_prefix", "post")
        post_recipe = generate_post_recipe_file(run_dir, case_id, post_cfg, {'Case': case_path, 'Post-Profile': post_path}, monitor_cfg)

        case_entries.append({
            "index": idx,
            "case_id": case_id,
            "run_dir": os.path.abspath(run_dir),
            "control_file": control_file,
            "post_recipe_file": post_recipe,
            "log_level": str(monitor_cfg.get("logging", {}).get("verbosity", "INFO")).upper(),
            "post_prefix": output_prefix,
            "solve_diagnostic_args": shlex.join(build_petsc_diagnostics_args(monitor_cfg, run_dir, "Solver")),
            "post_diagnostic_args": shlex.join(build_petsc_diagnostics_args(monitor_cfg, run_dir, "PostProcessor")),
            "parameters": combo,
        })

    with open(case_index_file, "w") as f:
        for entry in case_entries:
            f.write(
                "\t".join(
                    [
                        str(entry["index"]),
                        entry["case_id"],
                        entry["run_dir"],
                        entry["control_file"],
                        entry["post_recipe_file"],
                        entry["log_level"],
                        entry["post_prefix"],
                        entry["solve_diagnostic_args"],
                        entry["post_diagnostic_args"],
                    ]
                ) + "\n"
            )
    print(f"[SUCCESS] Wrote sweep case index: {os.path.relpath(case_index_file)}")

    max_idx = len(case_entries) - 1
    max_conc = study_cfg.get("execution", {}).get("max_concurrent_array_tasks")
    array_spec = f"0-{max_idx}"
    if max_conc:
        array_spec = f"{array_spec}%{max_conc}"

    solver_exe = resolve_runtime_executable("simulator")
    post_exe = resolve_runtime_executable("postprocessor")
    solver_array_script = os.path.join(scheduler_dir, "solver_array.sbatch")
    post_array_script = os.path.join(scheduler_dir, "post_array.sbatch")
    render_slurm_array_stage_script(
        solver_array_script,
        f"{study_id}_solve",
        cluster_cfg,
        array_spec,
        case_index_file,
        "solve",
        solver_exe,
        post_exe,
        os.path.join(scheduler_dir, "solver_%A_%a.out"),
        os.path.join(scheduler_dir, "solver_%A_%a.err")
    )
    render_slurm_array_stage_script(
        post_array_script,
        f"{study_id}_post",
        cluster_cfg,
        array_spec,
        case_index_file,
        "post",
        solver_exe,
        post_exe,
        os.path.join(scheduler_dir, "post_%A_%a.out"),
        os.path.join(scheduler_dir, "post_%A_%a.err")
    )
    print(f"[SUCCESS] Generated Slurm array scripts in {os.path.relpath(scheduler_dir)}")

    picurv_path = os.path.abspath(os.path.join(INVOKED_SCRIPT_DIR, "picurv"))
    metrics_aggregate_script = os.path.join(scheduler_dir, "metrics_aggregate.sbatch")
    render_metrics_aggregate_script(
        metrics_aggregate_script,
        f"{study_id}_metrics",
        cluster_cfg,
        study_dir,
        picurv_path,
    )
    print(f"[SUCCESS] Generated metrics aggregation script: {os.path.relpath(metrics_aggregate_script)}")

    submission = {
        "launch_mode": "slurm",
        "study_id": study_id,
        "solver_array": {"script": solver_array_script, "submitted": False},
        "post_array": {"script": post_array_script, "submitted": False},
        "metrics_aggregate": {"script": metrics_aggregate_script, "submitted": False},
        "no_submit": bool(args.no_submit),
    }
    if not args.no_submit:
        solver_submit = submit_sbatch(solver_array_script)
        submission["solver_array"].update(solver_submit)
        submission["solver_array"]["submitted"] = True
        post_submit = submit_sbatch(post_array_script, dependency=solver_submit["job_id"])
        submission["post_array"].update(post_submit)
        submission["post_array"]["submitted"] = True
        submission["post_array"]["dependency"] = f"afterok:{solver_submit['job_id']}"
        metrics_submit = submit_sbatch(metrics_aggregate_script, dependency=post_submit["job_id"], dependency_type="afterany")
        submission["metrics_aggregate"].update(metrics_submit)
        submission["metrics_aggregate"]["submitted"] = True
        submission["metrics_aggregate"]["dependency"] = f"afterany:{post_submit['job_id']}"
        print(f"[SUCCESS] Submitted solver array job:  {solver_submit['job_id']}")
        print(f"[SUCCESS] Submitted post array job:    {post_submit['job_id']}")
        print(f"[SUCCESS] Submitted metrics agg. job:  {metrics_submit['job_id']}")

    metrics_csv = aggregate_study_metrics(study_cfg, case_entries, results_dir)
    plots = generate_study_plots(study_cfg, metrics_csv, os.path.join(results_dir, "plots"))

    summary = {
        "study_id": study_id,
        "created_at": datetime.now().isoformat(),
        "git_commit": get_git_commit(),
        "study_type": study_cfg.get("study_type"),
        "num_cases": len(case_entries),
        "paths": {
            "study_dir": study_dir,
            "case_index": case_index_file,
            "solver_array_script": solver_array_script,
            "post_array_script": post_array_script,
            "metrics_table": metrics_csv,
            "plots_dir": os.path.join(results_dir, "plots"),
        },
        "submission": submission,
    }
    write_json_file(os.path.join(scheduler_dir, "submission.json"), submission)
    write_json_file(os.path.join(study_dir, "study_manifest.json"), summary)
    write_json_file(os.path.join(results_dir, "summary.json"), {"study_id": study_id, "metrics_csv": metrics_csv, "plots": plots})

    print("\n" + "=" * 60)
    print("                     STUDY SUMMARY")
    print("=" * 60)
    print(f"  Study ID        : {study_id}")
    print(f"  Study directory : {os.path.relpath(study_dir)}")
    print(f"  Cases generated : {len(case_entries)}")
    print(f"  Array spec      : {array_spec}")
    print(f"  Solver script   : {os.path.relpath(solver_array_script)}")
    print(f"  Post script     : {os.path.relpath(post_array_script)}")
    if metrics_csv:
        print(f"  Metrics table   : {os.path.relpath(metrics_csv)}")
    if plots:
        print(f"  Plots           : {os.path.relpath(os.path.join(results_dir, 'plots'))}")
    print("=" * 60)


def sweep_continue_workflow(args):
    """!
    @brief Continue a partially-completed Slurm parameter sweep study.
    @details Detects incomplete cases, prepares them for continuation (updating
             start_step, populating restart directories, regenerating control files),
             and submits new solver/post/metrics Slurm jobs. If all cases are already
             complete, performs metrics aggregation automatically.
    @param[in] args Parsed CLI arguments with study_dir and optional cluster override.
    """
    study_dir = os.path.abspath(args.study_dir)
    manifest_path = os.path.join(study_dir, "study_manifest.json")
    if not os.path.isfile(manifest_path):
        print(f"[FATAL] Study manifest not found: {manifest_path}", file=sys.stderr)
        sys.exit(1)
    manifest = _read_json_if_exists(manifest_path)
    study_id = manifest["study_id"]

    study_path = os.path.join(study_dir, "study.yml")
    cluster_path = os.path.abspath(args.cluster) if args.cluster else os.path.join(study_dir, "cluster.yml")
    study_cfg = read_yaml_file(study_path)
    cluster_cfg = read_yaml_file(cluster_path)
    validate_study_config(study_cfg, study_path, skip_base_file_check=True)
    validate_cluster_config(cluster_cfg, cluster_path)

    if args.cluster:
        shutil.copy(os.path.abspath(args.cluster), os.path.join(study_dir, "cluster.yml"))
        print(f"[INFO] Updated study cluster config from: {os.path.relpath(args.cluster)}")

    scheduler_dir = os.path.join(study_dir, "scheduler")
    cases_dir = os.path.join(study_dir, "cases")
    results_dir = os.path.join(study_dir, "results")
    case_index_file = os.path.join(scheduler_dir, "case_index.tsv")
    if not os.path.isfile(case_index_file):
        print(f"[FATAL] Case index not found: {case_index_file}", file=sys.stderr)
        sys.exit(1)

    parsed_entries = parse_case_index_tsv(case_index_file)

    base_cfgs = study_cfg["base_configs"]
    base_paths = {k: resolve_path(study_path, v) for k, v in base_cfgs.items()}
    base_case = read_yaml_file(base_paths["case"])

    combinations = expand_study_parameter_combinations(study_cfg)
    if len(combinations) != len(parsed_entries):
        print(
            f"[FATAL] Parameter matrix ({len(combinations)} cases) does not match "
            f"case_index.tsv ({len(parsed_entries)} entries).",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"\n[INFO] Study: {study_id}")
    print(f"[INFO] Scanning {len(combinations)} case(s) for completion status...")

    incomplete_indices = []
    all_case_entries = []
    for idx, combo in enumerate(combinations):
        case_id = f"case_{idx:04d}"
        entry = parsed_entries[idx]
        entry["parameters"] = combo
        run_dir = entry["run_dir"]

        effective_case = copy.deepcopy(base_case)
        for full_key, value in combo.items():
            root, nested = full_key.split(".", 1)
            if root == "case":
                _deep_set(effective_case, nested, value)
        try:
            eff_start = int(effective_case.get("run_control", {}).get("start_step", 0) or 0)
        except (TypeError, ValueError):
            eff_start = 0
        eff_total = int(effective_case["run_control"]["total_steps"])
        target = eff_start + eff_total

        monitor_cfg = read_yaml_file(os.path.join(run_dir, "config", "monitor.yml"))
        status = detect_case_completion_status(run_dir, monitor_cfg, target)
        entry["_status"] = status

        if status["status"] == "complete":
            print(f"  {case_id}: complete (step {status['last_step']}/{target})")
        elif status["status"] == "partial":
            print(f"  {case_id}: incomplete (step {status['last_step']}/{target}) — will continue")
            incomplete_indices.append(idx)
        else:
            print(f"  {case_id}: no checkpoint — will re-run from scratch")
            incomplete_indices.append(idx)

        all_case_entries.append(entry)

    if not incomplete_indices:
        print("\n[INFO] All cases are complete. Running metrics aggregation...")
        metrics_csv = aggregate_study_metrics(study_cfg, all_case_entries, results_dir)
        plots = generate_study_plots(study_cfg, metrics_csv, os.path.join(results_dir, "plots"))
        print("\n" + "=" * 60)
        print("              STUDY CONTINUATION SUMMARY")
        print("=" * 60)
        print(f"  Study ID        : {study_id}")
        print(f"  Status          : ALL COMPLETE")
        if metrics_csv:
            print(f"  Metrics table   : {os.path.relpath(metrics_csv)}")
        if plots:
            print(f"  Plots           : {os.path.relpath(os.path.join(results_dir, 'plots'))}")
        print("=" * 60)
        return

    print(f"\n[INFO] {len(incomplete_indices)} incomplete case(s) to continue/re-run.")

    skipped = []
    for idx in incomplete_indices:
        entry = all_case_entries[idx]
        status = entry["_status"]
        if status["status"] == "partial":
            prepare_case_for_continuation(
                entry["run_dir"], entry["case_id"],
                status["last_step"], status["target_step"],
                cluster_cfg,
            )
        elif status["status"] == "empty":
            print(f"[INFO] {entry['case_id']}: re-running from scratch (no control file changes)")

    solver_array_spec = ",".join(str(i) for i in incomplete_indices)
    max_conc = study_cfg.get("execution", {}).get("max_concurrent_array_tasks")
    if max_conc:
        solver_array_spec = f"{solver_array_spec}%{max_conc}"

    max_idx = len(combinations) - 1
    post_array_spec = f"0-{max_idx}"
    if max_conc:
        post_array_spec = f"{post_array_spec}%{max_conc}"

    solver_exe = resolve_runtime_executable("simulator")
    post_exe = resolve_runtime_executable("postprocessor")

    solver_continue_script = os.path.join(scheduler_dir, "solver_continue_array.sbatch")
    post_continue_script = os.path.join(scheduler_dir, "post_continue_array.sbatch")
    render_slurm_array_stage_script(
        solver_continue_script,
        f"{study_id}_solve_cont",
        cluster_cfg,
        solver_array_spec,
        case_index_file,
        "solve",
        solver_exe, post_exe,
        os.path.join(scheduler_dir, "solver_cont_%A_%a.out"),
        os.path.join(scheduler_dir, "solver_cont_%A_%a.err"),
    )
    render_slurm_array_stage_script(
        post_continue_script,
        f"{study_id}_post_cont",
        cluster_cfg,
        post_array_spec,
        case_index_file,
        "post",
        solver_exe, post_exe,
        os.path.join(scheduler_dir, "post_cont_%A_%a.out"),
        os.path.join(scheduler_dir, "post_cont_%A_%a.err"),
    )

    picurv_path = os.path.abspath(os.path.join(INVOKED_SCRIPT_DIR, "picurv"))
    metrics_aggregate_script = os.path.join(scheduler_dir, "metrics_continue_aggregate.sbatch")
    render_metrics_aggregate_script(
        metrics_aggregate_script,
        f"{study_id}_metrics_cont",
        cluster_cfg,
        study_dir,
        picurv_path,
    )
    print(f"[SUCCESS] Generated continuation scripts in {os.path.relpath(scheduler_dir)}")

    submission = {
        "launch_mode": "slurm",
        "study_id": study_id,
        "continuation": True,
        "incomplete_cases": [all_case_entries[i]["case_id"] for i in incomplete_indices],
        "solver_continue_array": {"script": solver_continue_script, "submitted": False},
        "post_continue_array": {"script": post_continue_script, "submitted": False},
        "metrics_aggregate": {"script": metrics_aggregate_script, "submitted": False},
        "no_submit": bool(args.no_submit),
    }
    if not args.no_submit:
        solver_submit = submit_sbatch(solver_continue_script)
        submission["solver_continue_array"].update(solver_submit)
        submission["solver_continue_array"]["submitted"] = True
        post_submit = submit_sbatch(post_continue_script, dependency=solver_submit["job_id"])
        submission["post_continue_array"].update(post_submit)
        submission["post_continue_array"]["submitted"] = True
        submission["post_continue_array"]["dependency"] = f"afterok:{solver_submit['job_id']}"
        metrics_submit = submit_sbatch(metrics_aggregate_script, dependency=post_submit["job_id"], dependency_type="afterany")
        submission["metrics_aggregate"].update(metrics_submit)
        submission["metrics_aggregate"]["submitted"] = True
        submission["metrics_aggregate"]["dependency"] = f"afterany:{post_submit['job_id']}"
        print(f"[SUCCESS] Submitted continuation solver array: {solver_submit['job_id']}")
        print(f"[SUCCESS] Submitted continuation post array:   {post_submit['job_id']}")
        print(f"[SUCCESS] Submitted metrics aggregation job:   {metrics_submit['job_id']}")

    write_json_file(os.path.join(scheduler_dir, "submission_continue.json"), submission)

    manifest["continuation"] = {
        "continued_at": datetime.now().isoformat(),
        "incomplete_cases": [all_case_entries[i]["case_id"] for i in incomplete_indices],
        "submission": submission,
    }
    write_json_file(manifest_path, manifest)

    print("\n" + "=" * 60)
    print("              STUDY CONTINUATION SUMMARY")
    print("=" * 60)
    print(f"  Study ID          : {study_id}")
    print(f"  Incomplete cases  : {len(incomplete_indices)}/{len(combinations)}")
    print(f"  Solver array spec : {solver_array_spec}")
    print(f"  Post array spec   : {post_array_spec}")
    print(f"  Solver script     : {os.path.relpath(solver_continue_script)}")
    print(f"  Post script       : {os.path.relpath(post_continue_script)}")
    print(f"  Metrics script    : {os.path.relpath(metrics_aggregate_script)}")
    if not args.no_submit:
        print(f"  [Metrics aggregation will run automatically after post-processing]")
    else:
        print(f"  [--no-submit] Scripts generated but not submitted.")
        print(f"  After manual submission and completion, run:")
        print(f"    picurv sweep --reaggregate --study-dir {os.path.relpath(study_dir)}")
    print("=" * 60)


def sweep_reaggregate_workflow(args):
    """!
    @brief Re-run metrics aggregation and plot generation for an existing study.
    @param[in] args Parsed CLI arguments with study_dir.
    """
    study_dir = os.path.abspath(args.study_dir)
    study_path = os.path.join(study_dir, "study.yml")
    if not os.path.isfile(study_path):
        print(f"[FATAL] Study config not found: {study_path}", file=sys.stderr)
        sys.exit(1)
    study_cfg = read_yaml_file(study_path)
    validate_study_config(study_cfg, study_path, skip_base_file_check=True)

    case_index_file = os.path.join(study_dir, "scheduler", "case_index.tsv")
    if not os.path.isfile(case_index_file):
        print(f"[FATAL] Case index not found: {case_index_file}", file=sys.stderr)
        sys.exit(1)

    parsed_entries = parse_case_index_tsv(case_index_file)
    combinations = expand_study_parameter_combinations(study_cfg)
    if len(combinations) != len(parsed_entries):
        print(
            f"[FATAL] Parameter matrix ({len(combinations)} cases) does not match "
            f"case_index.tsv ({len(parsed_entries)} entries).",
            file=sys.stderr,
        )
        sys.exit(1)

    case_entries = []
    for idx, combo in enumerate(combinations):
        entry = parsed_entries[idx]
        entry["parameters"] = combo
        case_entries.append(entry)

    results_dir = os.path.join(study_dir, "results")
    metrics_csv = aggregate_study_metrics(study_cfg, case_entries, results_dir)
    plots = generate_study_plots(study_cfg, metrics_csv, os.path.join(results_dir, "plots"))

    print("\n" + "=" * 60)
    print("              REAGGREGATION SUMMARY")
    print("=" * 60)
    if metrics_csv:
        print(f"  Metrics table   : {os.path.relpath(metrics_csv)}")
    if plots:
        print(f"  Plots generated : {len(plots)}")
    print("=" * 60)


_SUMMARY_NUMERIC_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _read_yaml_if_exists(filepath: str):
    """!
    @brief Read YAML when present, otherwise return None.
    @param[in] filepath Argument passed to `_read_yaml_if_exists()`.
    @return Value returned by `_read_yaml_if_exists()`.
    """
    if not filepath or not os.path.isfile(filepath):
        return None
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except yaml.YAMLError:
        return None


def _read_json_if_exists(filepath: str):
    """!
    @brief Read JSON when present, otherwise return None.
    @param[in] filepath Argument passed to `_read_json_if_exists()`.
    @return Value returned by `_read_json_if_exists()`.
    """
    if not filepath or not os.path.isfile(filepath):
        return None
    with open(filepath, "r", encoding="utf-8") as f:
        return json.load(f)


def _parse_int_loose(value):
    """!
    @brief Best-effort integer parsing for summary extraction.
    @param[in] value Argument passed to `_parse_int_loose()`.
    @return Value returned by `_parse_int_loose()`.
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        try:
            return int(float(text))
        except ValueError:
            return None


def _parse_float_loose(value):
    """!
    @brief Best-effort float parsing for summary extraction.
    @param[in] value Argument passed to `_parse_float_loose()`.
    @return Value returned by `_parse_float_loose()`.
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _extract_numeric_tuple(text: str):
    """!
    @brief Extract a numeric tuple from a string like '(1, 2, 3)'.
    @param[in] text Argument passed to `_extract_numeric_tuple()`.
    @return Value returned by `_extract_numeric_tuple()`.
    """
    if not text:
        return []
    return [float(token) for token in _SUMMARY_NUMERIC_RE.findall(text)]


def _build_summary_context(run_dir: str) -> dict:
    """!
    @brief Resolve run-local config and artifact paths for summarize.
    @param[in] run_dir Argument passed to `_build_summary_context()`.
    @return Value returned by `_build_summary_context()`.
    """
    run_dir = os.path.abspath(run_dir)
    if not os.path.isdir(run_dir):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key="run_dir",
            file_path=run_dir,
            message="Run directory not found.",
        )
        sys.exit(1)

    config_dir = os.path.join(run_dir, "config")
    config_paths = {
        "case": os.path.join(config_dir, "case.yml"),
        "solver": os.path.join(config_dir, "solver.yml"),
        "monitor": os.path.join(config_dir, "monitor.yml"),
    }
    monitor_cfg = _read_yaml_if_exists(config_paths["monitor"]) or {}
    case_cfg = _read_yaml_if_exists(config_paths["case"]) or {}
    solver_cfg = _read_yaml_if_exists(config_paths["solver"]) or {}
    manifest = _read_json_if_exists(os.path.join(run_dir, "manifest.json")) or {}

    io_cfg = monitor_cfg.get("io", {}) if isinstance(monitor_cfg, dict) else {}
    io_dirs = io_cfg.get("directories", {}) if isinstance(io_cfg, dict) else {}
    log_dir_name = io_dirs.get("log", "logs")
    scheduler_dir = os.path.join(run_dir, "scheduler")

    profiling_cfg = {"mode": "off", "functions": [], "timestep_file": "Profiling_Timestep_Summary.csv", "final_summary_enabled": True}
    if monitor_cfg:
        profiling_cfg = resolve_profiling_config(monitor_cfg)

    particle_console_output_freq = None
    particle_log_interval = None
    if monitor_cfg:
        particle_console_output_freq = resolve_particle_console_output_frequency(io_cfg)
        particle_log_interval = io_cfg.get("particle_log_interval")

    particle_count_cfg = None
    if case_cfg:
        particle_count_cfg = (
            case_cfg.get("models", {})
            .get("physics", {})
            .get("particles", {})
            .get("count")
        )

    return {
        "run_dir": run_dir,
        "config_dir": config_dir,
        "log_dir": os.path.join(run_dir, log_dir_name),
        "scheduler_dir": scheduler_dir,
        "monitor_cfg": monitor_cfg,
        "case_cfg": case_cfg,
        "solver_cfg": solver_cfg,
        "config_paths": config_paths,
        "manifest": manifest,
        "profiling_cfg": profiling_cfg,
        "particle_console_output_freq": particle_console_output_freq,
        "particle_log_interval": particle_log_interval,
        "particle_count_cfg": particle_count_cfg,
    }


def _require_summary_config(context: dict, name: str) -> dict:
    """!
    @brief Return one explicitly requested copied config or fail with a structured error.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @param[in] name Config selector name.
    @return Parsed config mapping.
    """
    path = context["config_paths"][name]
    cfg = context.get(f"{name}_cfg")
    if not os.path.isfile(path):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key=name,
            file_path=path,
            message=f"Copied run config '{name}.yml' was not found.",
            hint="Use a staged run directory containing the requested copied config.",
        )
        sys.exit(1)
    if not isinstance(cfg, dict) or not cfg:
        emit_structured_error(
            ERROR_CODE_CFG_INVALID_VALUE,
            key=name,
            file_path=path,
            message=f"Copied run config '{name}.yml' is empty or is not a YAML mapping.",
        )
        sys.exit(1)
    return cfg


def _build_run_overview(context: dict) -> dict:
    """!
    @brief Build timestep-independent run metadata for summarize.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @return Curated run metadata mapping.
    """
    manifest = context["manifest"]
    return {
        "run_id": manifest.get("run_id", os.path.basename(context["run_dir"])),
        "run_dir": context["run_dir"],
        "created_at": manifest.get("created_at"),
        "launch_mode": manifest.get("launch_mode"),
        "git_commit": manifest.get("git_commit"),
        "solver_num_procs": manifest.get("solver_num_procs", manifest.get("num_procs")),
        "post_num_procs": manifest.get("post_num_procs"),
        "stages_requested": manifest.get("stages_requested"),
        "stages_completed_or_submitted": manifest.get("stages_completed_or_submitted"),
    }


def _summarize_turbulence(turbulence_cfg: dict) -> dict:
    """!
    @brief Build compact turbulence and wall-model selections.
    @param[in] turbulence_cfg Case turbulence configuration mapping.
    @return Curated turbulence and wall-model mapping.
    """
    result = {}
    for key in ("les", "rans", "wall_function"):
        value = turbulence_cfg.get(key)
        if isinstance(value, dict):
            result[key] = {
                "enabled": value.get("enabled", True),
                "model": value.get("model"),
                **{k: v for k, v in value.items() if k not in {"enabled", "model"}},
            }
        elif value is not None:
            result[key] = value
    return result


def _build_case_overview(context: dict) -> dict:
    """!
    @brief Build a curated case.yml summary with useful derived quantities.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @return Curated case configuration mapping.
    """
    cfg = _require_summary_config(context, "case")
    props = cfg.get("properties", {})
    scaling = props.get("scaling", {})
    fluid = props.get("fluid", {})
    run = cfg.get("run_control", {})
    grid = cfg.get("grid", {})
    models = cfg.get("models", {})
    domain = models.get("domain", {})
    physics = models.get("physics", {})
    particles = physics.get("particles", {})
    start = int(run.get("start_step", 0))
    total = int(run.get("total_steps", 0))
    dt = float(run.get("dt_physical", 0.0))
    length_ref = float(scaling.get("length_ref"))
    velocity_ref = float(scaling.get("velocity_ref"))
    density = float(fluid.get("density"))
    viscosity = float(fluid.get("viscosity"))
    prepared_bcs = validate_and_prepare_boundary_conditions(cfg)
    first_block_faces = {row["face"]: row for row in prepared_bcs[0]}
    periodic_axes = {
        "i": first_block_faces["-Xi"]["type"] == "PERIODIC",
        "j": first_block_faces["-Eta"]["type"] == "PERIODIC",
        "k": first_block_faces["-Zeta"]["type"] == "PERIODIC",
    }
    bc_blocks = []
    for block_idx, block in enumerate(prepared_bcs):
        bc_blocks.append(
            {
                "block": block_idx,
                "faces": [
                    {"face": row["face"], "type": row["type"], "handler": row["handler"]}
                    for row in block
                ],
            }
        )
    return {
        "run_control": {
            "start_step": start,
            "total_steps": total,
            "end_step": start + total,
            "dt_physical": dt,
            "duration_physical": total * dt,
            "dt_nondimensional": dt * velocity_ref / length_ref,
        },
        "properties": {
            "length_ref": length_ref,
            "velocity_ref": velocity_ref,
            "density": density,
            "viscosity": viscosity,
            "reynolds_number": density * velocity_ref * length_ref / viscosity if viscosity else None,
            "initial_conditions": props.get("initial_conditions", {}),
        },
        "grid": {
            "mode": grid.get("mode"),
            "processor_layout": resolve_grid_da_processor_layout(grid),
            "programmatic_settings": grid.get("programmatic_settings") if grid.get("mode") == "programmatic_c" else None,
            "source_file": grid.get("source_file"),
        },
        "domain": {
            "blocks": domain.get("blocks", 1),
            "dimensionality": physics.get("dimensionality", "3D"),
            "periodic": periodic_axes,
        },
        "physics": {
            "fsi": physics.get("fsi", {}),
            "particles": particles,
            "turbulence": _summarize_turbulence(physics.get("turbulence", {}) or {}),
            "statistics": models.get("statistics", {}),
        },
        "boundary_conditions": bc_blocks,
    }


def _build_solver_overview(context: dict) -> dict:
    """!
    @brief Build a curated solver.yml summary with normalized selections.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @return Curated solver configuration mapping.
    """
    cfg = _require_summary_config(context, "solver")
    strategy = cfg.get("strategy", {}) or {}
    selected = normalize_momentum_solver_type(strategy.get("momentum_solver", "Dual Time Picard Jameson RK"))
    momentum_cfg = cfg.get("momentum_solver", {}) or {}
    dualtime = momentum_cfg.get("dual_time_picard_jameson_rk", momentum_cfg.get("dual_time_picard_rk4", {})) or {}
    poisson = cfg.get("poisson_solver", cfg.get("pressure_solver", {})) or {}
    convergence = cfg.get("solution_convergence", {}) or {}
    convergence_mode = normalize_solution_convergence_mode(convergence.get("mode", "steady_deterministic"))
    operation_mode = cfg.get("operation_mode", {}) or {}
    operation_mode = {
        **operation_mode,
        "eulerian_field_source": normalize_eulerian_field_source(operation_mode.get("eulerian_field_source", "solve")),
    }
    if operation_mode.get("analytical_type") is not None:
        operation_mode["analytical_type"] = normalize_analytical_type(operation_mode["analytical_type"])
    passthrough = cfg.get("petsc_passthrough_options", {}) or {}
    return {
        "operation_mode": operation_mode,
        "momentum": {
            "type": selected,
            "central_diff": bool(strategy.get("central_diff", False)),
            "tolerances": cfg.get("tolerances", {}),
            "controls": dualtime,
        },
        "poisson": poisson,
        "interpolation": cfg.get("interpolation", {"method": "Trilinear"}),
        "solution_convergence": {**convergence, "mode": convergence_mode},
        "scalar_transport": cfg.get("scalar_transport", {}),
        "verification": cfg.get("verification", {}),
        "petsc_passthrough": {"count": len(passthrough), "options": sorted(passthrough.keys())},
    }


def _build_monitor_overview(context: dict) -> dict:
    """!
    @brief Build a curated monitor.yml summary with resolved defaults.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @return Curated monitor configuration mapping.
    """
    cfg = _require_summary_config(context, "monitor")
    logging_cfg = cfg.get("logging", {}) or {}
    io_cfg = cfg.get("io", {}) or {}
    diagnostics = resolve_diagnostics_config(cfg, context["run_dir"], "Solver")
    monitoring_flags = resolve_solver_monitoring_flags(cfg)
    enabled_petsc = sorted(key for key, value in diagnostics["petsc"].items() if value not in (False, None))
    return {
        "logging": {
            "verbosity": logging_cfg.get("verbosity", "WARNING"),
            "enabled_functions": logging_cfg.get("enabled_functions", []),
        },
        "profiling": resolve_profiling_config(cfg),
        "diagnostics": {
            "enabled_petsc": enabled_petsc,
            "petsc": diagnostics["petsc"],
            "runtime_memory_log": diagnostics["runtime_memory_log"],
        },
        "io": {
            "data_output_frequency": io_cfg.get("data_output_frequency"),
            "particle_console_output_frequency": resolve_particle_console_output_frequency(io_cfg),
            "particle_log_interval": io_cfg.get("particle_log_interval"),
            "directories": io_cfg.get("directories", {}),
        },
        "solver_monitoring": {
            "enabled_flags": sorted(flag for flag, value in monitoring_flags.items() if value not in (False, None)),
            "flags": monitoring_flags,
        },
    }


def _parse_continuity_metrics_log(filepath: str) -> "tuple[dict, list[int]]":
    """!
    @brief Parse Continuity_Metrics.log into latest rows by step plus observed order.
    @param[in] filepath Argument passed to `_parse_continuity_metrics_log()`.
    @return Value returned by `_parse_continuity_metrics_log()`.
    """
    rows_by_step = {}
    step_order = []
    active_step = None
    if not os.path.isfile(filepath):
        return rows_by_step, step_order

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("-") or line.startswith("Timestep"):
                continue
            parts = [part.strip() for part in raw_line.split("|")]
            if len(parts) < 8:
                continue
            step = _parse_int_loose(parts[0])
            block = _parse_int_loose(parts[1])
            max_div = _parse_float_loose(parts[2])
            rhs_sum = _parse_float_loose(parts[4])
            flux_in = _parse_float_loose(parts[5])
            flux_out = _parse_float_loose(parts[6])
            net_flux = _parse_float_loose(parts[7])
            if step is None or block is None:
                continue
            if step != active_step:
                active_step = step
                step_order.append(step)
                rows_by_step[step] = {}
            rows_by_step.setdefault(step, {})[block] = {
                "block": block,
                "max_divergence": max_div,
                "max_divergence_location": parts[3],
                "rhs_sum": rhs_sum,
                "flux_in": flux_in,
                "flux_out": flux_out,
                "net_flux": net_flux,
            }
    return {step: list(block_rows.values()) for step, block_rows in rows_by_step.items()}, step_order


def _parse_particle_metrics_log(filepath: str) -> "tuple[dict, list[int]]":
    """!
    @brief Parse Particle_Metrics.log into latest rows by step plus observed order.
    @param[in] filepath Argument passed to `_parse_particle_metrics_log()`.
    @return Value returned by `_parse_particle_metrics_log()`.
    """
    rows_by_step = {}
    step_order = []
    if not os.path.isfile(filepath):
        return rows_by_step, step_order

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("-") or line.startswith("Stage"):
                continue
            parts = [part.strip() for part in raw_line.split("|")]
            if len(parts) < 8:
                continue
            step = _parse_int_loose(parts[1])
            if step is None:
                continue
            row = {
                "stage": parts[0],
                "total_particles": _parse_int_loose(parts[2]),
                "lost_particles": _parse_int_loose(parts[3]),
                "lost_particles_cumulative": None,
                "migrated_particles": None,
                "occupied_cells": None,
                "load_imbalance": None,
                "migration_passes": None,
            }
            if len(parts) >= 9:
                row.update(
                    {
                        "lost_particles_cumulative": _parse_int_loose(parts[4]),
                        "migrated_particles": _parse_int_loose(parts[5]),
                        "occupied_cells": _parse_int_loose(parts[6]),
                        "load_imbalance": _parse_float_loose(parts[7]),
                        "migration_passes": _parse_int_loose(parts[8]),
                    }
                )
            else:
                row.update(
                    {
                        "migrated_particles": _parse_int_loose(parts[4]),
                        "occupied_cells": _parse_int_loose(parts[5]),
                        "load_imbalance": _parse_float_loose(parts[6]),
                        "migration_passes": _parse_int_loose(parts[7]),
                    }
                )
            rows_by_step[step] = row
            step_order.append(step)
    return rows_by_step, step_order


def _parse_momentum_convergence_logs(log_dir: str) -> "tuple[dict, dict, list[int]]":
    """!
    @brief Parse per-block momentum convergence logs.
    @param[in] log_dir Argument passed to `_parse_momentum_convergence_logs()`.
    @return Value returned by `_parse_momentum_convergence_logs()`.
    """
    rows_by_step = {}
    sources = {}
    step_order = []
    pattern = os.path.join(log_dir, "Momentum_Solver_Convergence_History_Block_*.log")
    regex = re.compile(
        r"Step:\s*(?P<step>\d+)\s*\|\s*PseudoIter\(k\):\s*(?P<pseudo_iter>\d+)\s*\|"
        r"\s*Pseudo-cfl:\s*(?P<pseudo_cfl>[-+0-9.eE]+)\s*\|\s*\|dUk\|:\s*(?P<delta>[-+0-9.eE]+)\s*\|"
        r"\s*\|dUk\|/\|dU0\|:\s*(?P<delta_rel>[-+0-9.eE]+)\s*\|\s*\|Rk\|:\s*(?P<resid>[-+0-9.eE]+)\s*\|"
        r"\s*\|Rk\|/\|R0\|:\s*(?P<resid_rel>[-+0-9.eE]+)"
        r"(?:\s*\|\s*trial_ratio:\s*(?P<trial_ratio>[-+0-9.eE]+)\s*\|\s*status:\s*(?P<status>\w+)\s*\|\s*cfl_after:\s*(?P<cfl_after>[-+0-9.eE]+))?"
    )

    # Per (step, block): track accepted/rejected counts and last committed state.
    _state = {}

    for path in sorted(glob.glob(pattern)):
        block_match = re.search(r"Block_(\d+)\.log$", path)
        if not block_match:
            continue
        block = int(block_match.group(1))
        sources[block] = path
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                match = regex.search(raw_line)
                if not match:
                    continue
                step = int(match.group("step"))
                step_order.append(step)
                key = (step, block)
                if key not in _state:
                    _state[key] = {"accepted_count": 0, "rejected_count": 0, "last_accepted": None, "last_rejected": None}
                entry = _state[key]
                status = match.group("status")
                row = {
                    "block": block,
                    "pseudo_iterations": int(match.group("pseudo_iter")),
                    "pseudo_cfl": float(match.group("pseudo_cfl")),
                    "cfl_after": _parse_float_loose(match.group("cfl_after")),
                    "delta_norm": float(match.group("delta")),
                    "delta_rel": float(match.group("delta_rel")),
                    "residual_norm": float(match.group("resid")),
                    "residual_rel": float(match.group("resid_rel")),
                    "trial_ratio": _parse_float_loose(match.group("trial_ratio")),
                    "status": status,
                }
                if status == "accepted":
                    entry["accepted_count"] += 1
                    entry["last_accepted"] = row
                else:
                    entry["rejected_count"] += 1
                    entry["last_rejected"] = row

    for (step, block), entry in _state.items():
        display_row = entry["last_accepted"] or entry["last_rejected"]
        if display_row:
            rows_by_step.setdefault(step, {})[block] = {
                **display_row,
                "accepted_count": entry["accepted_count"],
                "rejected_count": entry["rejected_count"],
            }

    return rows_by_step, sources, step_order


def _parse_poisson_convergence_logs(log_dir: str) -> "tuple[dict, dict, list[int]]":
    """!
    @brief Parse per-block Poisson convergence logs.
    @param[in] log_dir Argument passed to `_parse_poisson_convergence_logs()`.
    @return Value returned by `_parse_poisson_convergence_logs()`.
    """
    rows_by_step = {}
    sources = {}
    step_order = []
    pattern = os.path.join(log_dir, "Poisson_Solver_Convergence_History_Block_*.log")
    regex = re.compile(
        r"ts:\s*(?P<step>\d+)\s*\|\s*block:\s*(?P<block>\d+)\s*\|\s*iter:\s*(?P<iter>\d+)\s*\|"
        r"\s*Unprecond Norm:\s*(?P<unpre>[-+0-9.eE]+)\s*\|\s*True Norm:\s*(?P<true>[-+0-9.eE]+)"
        r"(?:\s*\|\s*Rel Norm:\s*(?P<rel>[-+0-9.eE]+))?"
    )

    for path in sorted(glob.glob(pattern)):
        block_match = re.search(r"Block_(\d+)\.log$", path)
        if not block_match:
            continue
        block = int(block_match.group(1))
        sources[block] = path
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                match = regex.search(raw_line)
                if not match:
                    continue
                step = int(match.group("step"))
                step_order.append(step)
                rows_by_step.setdefault(step, {})[block] = {
                    "block": block,
                    "iterations": int(match.group("iter")),
                    "unpreconditioned_norm": float(match.group("unpre")),
                    "true_norm": float(match.group("true")),
                    "relative_norm": _parse_float_loose(match.group("rel")),
                }
    return rows_by_step, sources, step_order


def _parse_profiling_timestep_csv(filepath: str) -> "tuple[dict, list[int]]":
    """!
    @brief Parse profiling timestep CSV into latest rows by step plus observed order.
    @param[in] filepath Argument passed to `_parse_profiling_timestep_csv()`.
    @return Value returned by `_parse_profiling_timestep_csv()`.
    """
    rows_by_step = {}
    step_order = []
    active_step = None
    if not os.path.isfile(filepath):
        return rows_by_step, step_order

    with open(filepath, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            step = _parse_int_loose(row.get("step"))
            if step is None:
                continue
            if step != active_step:
                active_step = step
                step_order.append(step)
                rows_by_step[step] = []
            rows_by_step.setdefault(step, []).append(
                {
                    "function": row.get("function"),
                    "calls": _parse_int_loose(row.get("calls")),
                    "step_time_s": _parse_float_loose(row.get("step_time_s")),
                }
            )
    return rows_by_step, step_order


def _parse_runtime_memory_log(filepath: str) -> "tuple[dict, list[int], dict]":
    """!
    @brief Parse Runtime_Memory.log into latest rows by step and final status.
    @param[in] filepath Runtime memory log path.
    @return Tuple of rows by step, observed step order, and final/shutdown metadata.
    """
    rows_by_step = {}
    step_order = []
    final_row = None
    latest_sample_row = None
    max_process_change_mb = None
    if not os.path.isfile(filepath):
        return rows_by_step, step_order, {"available": False}

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("Step"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            step = _parse_int_loose(parts[0])
            if step is None:
                continue
            row = {
                "step": step,
                "event": parts[1],
                "process_current_mb_max": _parse_float_loose(parts[2]),
                "process_peak_mb_max": _parse_float_loose(parts[3]),
                "petsc_allocated_mb_max": _parse_float_loose(parts[4]),
                "petsc_peak_allocated_mb_max": _parse_float_loose(parts[5]),
                "process_change_mb_max": _parse_float_loose(parts[6]),
                "reason": parts[7],
            }
            if row["process_change_mb_max"] is not None:
                max_process_change_mb = (
                    row["process_change_mb_max"]
                    if max_process_change_mb is None
                    else max(max_process_change_mb, row["process_change_mb_max"])
                )
            if row["event"] in {"Step", "Post"}:
                rows_by_step[step] = row
                step_order.append(step)
                latest_sample_row = row
            elif row["event"] in {"Shutdown", "Final"}:
                final_row = row

    meta = {
        "available": bool(rows_by_step or final_row),
        "source": filepath,
        "final_event": final_row.get("event") if final_row else None,
        "final_reason": final_row.get("reason") if final_row else None,
        "max_process_change_mb": max_process_change_mb,
        "latest_sample_row": latest_sample_row,
        "final_row": final_row,
    }
    return rows_by_step, step_order, meta


def _parse_solution_convergence_log(filepath: str) -> "tuple[dict, list[int]]":
    """!
    @brief Parse solution_convergence.log into latest rows by step plus observed order.

    The log format uses pipe-delimited aligned columns. The first line of the
    file is a banner (starts with '=') containing the mode tag; the second line
    is the column header; the third line is a separator (starts with '-').
    Subsequent lines are one data row per timestep.

    @param[in] filepath Path to solution_convergence.log.
    @return Mapping of step number to a dict of column values.
    """
    rows_by_step = {}
    step_order = []
    if not os.path.isfile(filepath):
        return rows_by_step, step_order

    mode = None
    col_names = None

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("="):
                m = re.search(r"\[mode:\s*([\w_]+)", line)
                if m:
                    mode = m.group(1)
                continue
            if line.startswith("-"):
                continue
            if col_names is None:
                col_names = [p.strip() for p in raw_line.split("|")]
                continue
            parts = [p.strip() for p in raw_line.split("|")]
            if len(parts) < 4 or col_names is None:
                continue
            step = _parse_int_loose(parts[0]) if col_names[0] == "step" else None
            if step is None:
                continue
            step_order.append(step)
            row = {"mode": mode}
            for name, val in zip(col_names, parts):
                if not name:
                    continue
                float_val = _parse_float_loose(val)
                int_val = _parse_int_loose(val)
                if float_val is not None and ("." in val or "e" in val.lower()):
                    row[name] = float_val
                elif int_val is not None:
                    row[name] = int_val
                else:
                    row[name] = val
            rows_by_step[step] = row
    return rows_by_step, step_order


def _find_solver_stream_log_candidates(run_dir: str, log_dir: str) -> "list[str]":
    """!
    @brief Return plausible solver stream logs for local and Slurm runs.
    @param[in] run_dir Argument passed to `_find_solver_stream_log_candidates()`.
    @param[in] log_dir Argument passed to `_find_solver_stream_log_candidates()`.
    @return Value returned by `_find_solver_stream_log_candidates()`.
    """
    patterns = [
        os.path.join(run_dir, "scheduler", "*_solver.log"),
        os.path.join(run_dir, "scheduler", "solver_*.out"),
        os.path.join(log_dir, "*_solver.log"),
    ]
    found = []
    seen = set()
    for pattern in patterns:
        for path in sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True):
            if path not in seen:
                seen.add(path)
                found.append(path)
    return found


def _parse_particle_snapshot_file(filepath: str) -> dict:
    """!
    @brief Parse sampled particle snapshots from a solver stream log.
    @param[in] filepath Argument passed to `_parse_particle_snapshot_file()`.
    @return Value returned by `_parse_particle_snapshot_file()`.
    """
    snapshots = {}
    if not os.path.isfile(filepath):
        return snapshots

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    idx = 0
    while idx < len(lines):
        match = re.search(r"Particle states at step\s+(\d+):", lines[idx])
        if not match:
            idx += 1
            continue
        step = int(match.group(1))
        idx += 1
        rows = []
        while idx < len(lines):
            stripped = lines[idx].strip()
            if re.search(r"Particle states at step\s+\d+:", lines[idx]):
                break
            if stripped.startswith("|"):
                parts = [part.strip() for part in stripped.split("|")[1:-1]]
                if len(parts) >= 6 and parts[0] != "Rank":
                    velocity = _extract_numeric_tuple(parts[4])
                    rows.append(
                        {
                            "rank": _parse_int_loose(parts[0]),
                            "pid": _parse_int_loose(parts[1]),
                            "cell": [int(value) for value in _extract_numeric_tuple(parts[2])],
                            "position": _extract_numeric_tuple(parts[3]),
                            "velocity": velocity,
                            "weights": _extract_numeric_tuple(parts[5]),
                            "sample_speed": math.sqrt(sum(component * component for component in velocity)) if len(velocity) == 3 else None,
                        }
                    )
            elif rows and (not stripped or stripped.startswith("Progress:")):
                break
            idx += 1
        if rows:
            snapshots[step] = rows
        else:
            snapshots.setdefault(step, [])
    return snapshots


def _find_previous_snapshot_step(snapshot_steps: "list[int]", step: int) -> "int | None":
    """!
    @brief Return the nearest earlier snapshot step when available.
    @param[in] snapshot_steps Argument passed to `_find_previous_snapshot_step()`.
    @param[in] step Argument passed to `_find_previous_snapshot_step()`.
    @return Value returned by `_find_previous_snapshot_step()`.
    """
    earlier_steps = [candidate for candidate in snapshot_steps if candidate < step]
    if not earlier_steps:
        return None
    return max(earlier_steps)


def _compute_particle_snapshot_delta(current_rows: "list[dict]", previous_rows: "list[dict]") -> dict:
    """!
    @brief Compute sampled deltas between two particle snapshot samples.
    @param[in] current_rows Argument passed to `_compute_particle_snapshot_delta()`.
    @param[in] previous_rows Argument passed to `_compute_particle_snapshot_delta()`.
    @return Value returned by `_compute_particle_snapshot_delta()`.
    """
    np = require_numpy()
    previous_by_pid = {
        row.get("pid"): row
        for row in previous_rows
        if row.get("pid") is not None
    }
    current_by_pid = {
        row.get("pid"): row
        for row in current_rows
        if row.get("pid") is not None
    }
    matched_pids = sorted(set(previous_by_pid) & set(current_by_pid))
    if not matched_pids:
        return {"available": False}

    displacements = []
    rank_migrations = 0
    cell_changes = 0
    speed_changes = []
    for pid in matched_pids:
        current_row = current_by_pid[pid]
        previous_row = previous_by_pid[pid]
        current_pos = current_row.get("position") or []
        previous_pos = previous_row.get("position") or []
        if len(current_pos) == len(previous_pos) and current_pos:
            displacements.append(
                math.sqrt(
                    sum(
                        (float(current_pos[idx]) - float(previous_pos[idx])) ** 2
                        for idx in range(len(current_pos))
                    )
                )
            )
        current_speed = current_row.get("sample_speed")
        previous_speed = previous_row.get("sample_speed")
        if current_speed is not None and previous_speed is not None:
            speed_changes.append(float(current_speed) - float(previous_speed))
        if current_row.get("rank") is not None and previous_row.get("rank") is not None:
            if current_row["rank"] != previous_row["rank"]:
                rank_migrations += 1
        if current_row.get("cell") and previous_row.get("cell"):
            if current_row["cell"] != previous_row["cell"]:
                cell_changes += 1

    payload = {
        "available": True,
        "matched_pids": len(matched_pids),
        "new_count": len(set(current_by_pid) - set(previous_by_pid)),
        "gone_count": len(set(previous_by_pid) - set(current_by_pid)),
        "rank_migrations": rank_migrations,
        "cell_changes": cell_changes,
    }
    if displacements:
        payload["mean_displacement"] = float(np.mean(displacements))
        payload["max_displacement"] = float(np.max(displacements))
    if speed_changes:
        payload["mean_speed_change"] = float(np.mean(speed_changes))
        payload["max_abs_speed_change"] = float(np.max(np.abs(speed_changes)))
    return payload


def _build_particle_snapshot_summary(
    source: str,
    step: int,
    rows: "list[dict]",
    preview_rows: int,
    particle_console_output_freq,
    particle_log_interval,
    previous_step: "int | None" = None,
    previous_rows: "list[dict] | None" = None,
) -> dict:
    """!
    @brief Build sampled diagnostics for one particle console snapshot.
    @param[in] source Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] step Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] rows Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] preview_rows Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] particle_console_output_freq Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] particle_log_interval Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] previous_step Argument passed to `_build_particle_snapshot_summary()`.
    @param[in] previous_rows Argument passed to `_build_particle_snapshot_summary()`.
    @return Value returned by `_build_particle_snapshot_summary()`.
    """
    np = require_numpy()
    payload = {
        "available": True,
        "sampled": True,
        "source": source,
        "step": step,
        "sampled_rows": len(rows),
        "preview_rows": rows[:preview_rows],
        "cadence": {
            "particle_console_output_frequency": particle_console_output_freq,
            "particle_log_interval": particle_log_interval,
        },
    }
    if not rows:
        return payload

    rank_counts = {}
    duplicate_pid_count = 0
    duplicate_cell_count = 0
    nan_count = 0
    inf_count = 0
    zero_weight_count = 0
    negative_weight_count = 0
    unique_pid_count = 0

    seen_pids = set()
    cell_counter = {}
    position_components = [[], [], []]
    weight_components = {}
    speeds = []

    for row in rows:
        pid = row.get("pid")
        if pid is not None:
            if pid in seen_pids:
                duplicate_pid_count += 1
            else:
                seen_pids.add(pid)
        rank = row.get("rank")
        if rank is not None:
            rank_counts[str(rank)] = rank_counts.get(str(rank), 0) + 1

        cell = row.get("cell") or []
        if cell:
            key = tuple(cell)
            cell_counter[key] = cell_counter.get(key, 0) + 1

        position = row.get("position") or []
        for idx, value in enumerate(position[:3]):
            if not np.isfinite(value):
                if np.isnan(value):
                    nan_count += 1
                else:
                    inf_count += 1
                continue
            position_components[idx].append(float(value))

        velocity = row.get("velocity") or []
        if any(not np.isfinite(value) for value in velocity):
            for value in velocity:
                if not np.isfinite(value):
                    if np.isnan(value):
                        nan_count += 1
                    else:
                        inf_count += 1
        speed = row.get("sample_speed")
        if speed is not None and np.isfinite(speed):
            speeds.append(float(speed))

        weights = row.get("weights") or []
        for idx, value in enumerate(weights):
            if not np.isfinite(value):
                if np.isnan(value):
                    nan_count += 1
                else:
                    inf_count += 1
                continue
            numeric = float(value)
            weight_components.setdefault(idx, []).append(numeric)
            if abs(numeric) <= 1.0e-15:
                zero_weight_count += 1
            if numeric < 0.0:
                negative_weight_count += 1

    unique_pid_count = len(seen_pids)
    duplicate_cell_count = sum(1 for count in cell_counter.values() if count > 1)

    payload["sampled_distribution"] = {
        "unique_cells": len(cell_counter),
        "duplicate_cells": duplicate_cell_count,
        "rank_counts": rank_counts,
        "unique_pids": unique_pid_count,
    }
    payload["checks"] = {
        "duplicate_pid_count": duplicate_pid_count,
        "nan_count": nan_count,
        "inf_count": inf_count,
        "zero_weight_count": zero_weight_count,
        "negative_weight_count": negative_weight_count,
    }

    if speeds:
        payload["speed"] = {
            "min": float(np.min(speeds)),
            "mean": float(np.mean(speeds)),
            "max": float(np.max(speeds)),
            "std": float(np.std(speeds)),
            "stagnant_count": sum(1 for speed in speeds if abs(speed) < 1.0e-6),
        }

        fastest_rows = sorted(
            [row for row in rows if row.get("sample_speed") is not None],
            key=lambda row: row["sample_speed"],
            reverse=True,
        )[:3]
        payload["top_speeds"] = [
            {
                "pid": row.get("pid"),
                "rank": row.get("rank"),
                "speed": row.get("sample_speed"),
                "cell": row.get("cell"),
            }
            for row in fastest_rows
        ]

    if any(position_components):
        axes = ["x", "y", "z"]
        payload["position_bounds"] = {}
        centroid = []
        for idx, axis in enumerate(axes):
            values = position_components[idx]
            if values:
                payload["position_bounds"][axis] = [float(np.min(values)), float(np.max(values))]
                centroid.append(float(np.mean(values)))
            else:
                centroid.append(None)
        payload["position_centroid"] = centroid

    if weight_components:
        payload["weights"] = {}
        for idx, values in sorted(weight_components.items()):
            payload["weights"][f"component_{idx}"] = {
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }

    delta_summary = {"available": False}
    if previous_step is not None and previous_rows:
        delta_summary = _compute_particle_snapshot_delta(rows, previous_rows)
        if delta_summary.get("available"):
            delta_summary["previous_step"] = previous_step
    payload["delta_from_previous_snapshot"] = delta_summary
    return payload


def _find_particle_snapshot_for_step(
    run_dir: str,
    log_dir: str,
    step: int,
    preview_rows: int,
    particle_console_output_freq,
    particle_log_interval,
) -> dict:
    """!
    @brief Locate and summarize a particle console snapshot for one step.
    @param[in] run_dir Argument passed to `_find_particle_snapshot_for_step()`.
    @param[in] log_dir Argument passed to `_find_particle_snapshot_for_step()`.
    @param[in] step Argument passed to `_find_particle_snapshot_for_step()`.
    @param[in] preview_rows Argument passed to `_find_particle_snapshot_for_step()`.
    @param[in] particle_console_output_freq Argument passed to `_find_particle_snapshot_for_step()`.
    @param[in] particle_log_interval Argument passed to `_find_particle_snapshot_for_step()`.
    @return Value returned by `_find_particle_snapshot_for_step()`.
    """
    best = None
    for path in _find_solver_stream_log_candidates(run_dir, log_dir):
        snapshots = _parse_particle_snapshot_file(path)
        rows = snapshots.get(step)
        if not rows:
            continue
        if best is None or len(rows) > len(best["rows"]):
            best = {"source": path, "rows": rows, "snapshots": snapshots}

    if not best:
        return {"available": False}

    previous_step = _find_previous_snapshot_step(list(best["snapshots"].keys()), step)
    previous_rows = best["snapshots"].get(previous_step, []) if previous_step is not None else None
    return _build_particle_snapshot_summary(
        best["source"],
        step,
        best["rows"],
        preview_rows,
        particle_console_output_freq=particle_console_output_freq,
        particle_log_interval=particle_log_interval,
        previous_step=previous_step,
        previous_rows=previous_rows,
    )


def _resolve_summary_step(
    requested_step,
    continuity_rows,
    particle_rows,
    momentum_rows,
    poisson_rows,
    profiling_rows,
    memory_rows=None,
    convergence_rows=None,
    step_orders=None,
    selection_mode: str = "latest",
):
    """!
    @brief Select a step to summarize from available metric artifacts.
    @param[in] requested_step Argument passed to `_resolve_summary_step()`.
    @param[in] continuity_rows Argument passed to `_resolve_summary_step()`.
    @param[in] particle_rows Argument passed to `_resolve_summary_step()`.
    @param[in] momentum_rows Argument passed to `_resolve_summary_step()`.
    @param[in] poisson_rows Argument passed to `_resolve_summary_step()`.
    @param[in] profiling_rows Argument passed to `_resolve_summary_step()`.
    @param[in] memory_rows Argument passed to `_resolve_summary_step()`.
    @param[in] convergence_rows Argument passed to `_resolve_summary_step()`.
    @param[in] step_orders Argument passed to `_resolve_summary_step()`.
    @param[in] selection_mode Argument passed to `_resolve_summary_step()`.
    @return Value returned by `_resolve_summary_step()`.
    """
    if memory_rows is None:
        memory_rows = {}
    if convergence_rows is None:
        convergence_rows = {}
    if step_orders is None:
        step_orders = []
    available_steps = (
        set(continuity_rows) | set(particle_rows) | set(momentum_rows)
        | set(poisson_rows) | set(profiling_rows) | set(memory_rows) | set(convergence_rows)
    )
    if not available_steps:
        return None, []

    if requested_step is not None:
        return requested_step, sorted(available_steps)

    if selection_mode == "max_step":
        return max(available_steps), sorted(available_steps)

    for order in step_orders:
        if order:
            return order[-1], sorted(available_steps)
    return max(available_steps), sorted(available_steps)


def _format_summary_float(value, spec: str = ".6e", missing: str = "n/a") -> str:
    """!
    @brief Format optional numeric values for summary text output.
    @param[in] value Argument passed to `_format_summary_float()`.
    @param[in] spec Argument passed to `_format_summary_float()`.
    @param[in] missing Argument passed to `_format_summary_float()`.
    @return Value returned by `_format_summary_float()`.
    """
    if value is None:
        return missing
    return format(value, spec)


def _summary_source_mtime(paths) -> float:
    """!
    @brief Return the newest modification time among one or more summary sources.
    @param[in] paths Path string, iterable of paths, or mapping of paths.
    @return Newest modification time, or -1.0 when no source exists.
    """
    if isinstance(paths, dict):
        paths = paths.values()
    elif isinstance(paths, str):
        paths = [paths]
    newest = -1.0
    for path in paths or []:
        if path and os.path.isfile(path):
            newest = max(newest, os.path.getmtime(path))
    return newest


def _order_summary_step_orders(sources: "list[tuple[list[int], object]]") -> "list[list[int]]":
    """!
    @brief Order observed step sequences by the recency of their source files.
    @param[in] sources Pairs of observed steps and filesystem source path(s).
    @return Step-order lists sorted so active append sources are considered first.
    """
    ranked = []
    for priority, (order, paths) in enumerate(sources):
        if order:
            ranked.append((_summary_source_mtime(paths), priority, order))
    ranked.sort(key=lambda item: (-item[0], item[1]))
    return [order for _, _, order in ranked]


def build_run_summary_payload(run_dir: str, step: "int | None" = None, snapshot_rows: int = 5, selection_mode: str = "latest") -> dict:
    """!
    @brief Build a read-only run-step summary from existing PICurv artifacts.
    @param[in] run_dir Argument passed to `build_run_summary_payload()`.
    @param[in] step Argument passed to `build_run_summary_payload()`.
    @param[in] snapshot_rows Argument passed to `build_run_summary_payload()`.
    @param[in] selection_mode Argument passed to `build_run_summary_payload()`.
    @return Value returned by `build_run_summary_payload()`.
    """
    context = _build_summary_context(run_dir)
    log_dir = context["log_dir"]
    continuity_path = os.path.join(log_dir, "Continuity_Metrics.log")
    particle_metrics_path = os.path.join(log_dir, "Particle_Metrics.log")

    continuity_rows, continuity_order = _parse_continuity_metrics_log(continuity_path)
    particle_rows, particle_order = _parse_particle_metrics_log(particle_metrics_path)
    momentum_rows, momentum_sources, momentum_order = _parse_momentum_convergence_logs(log_dir)
    poisson_rows, poisson_sources, poisson_order = _parse_poisson_convergence_logs(log_dir)

    profiling_rows = {}
    profiling_order = []
    profiling_path = os.path.join(log_dir, context["profiling_cfg"].get("timestep_file", "Profiling_Timestep_Summary.csv"))
    if context["profiling_cfg"].get("mode") != "off":
        profiling_rows, profiling_order = _parse_profiling_timestep_csv(profiling_path)

    diagnostics_cfg = resolve_diagnostics_config(context["monitor_cfg"])
    memory_log_file = diagnostics_cfg["runtime_memory_log"].get("file", "Runtime_Memory.log")
    memory_path = os.path.join(log_dir, memory_log_file)
    memory_rows, memory_order, memory_meta = _parse_runtime_memory_log(memory_path)

    convergence_log_path = os.path.join(log_dir, "solution_convergence.log")
    convergence_rows, convergence_order = _parse_solution_convergence_log(convergence_log_path)
    step_orders = _order_summary_step_orders(
        [
            (continuity_order, continuity_path),
            (particle_order, particle_metrics_path),
            (convergence_order, convergence_log_path),
            (profiling_order, profiling_path),
            (memory_order, memory_path),
            (momentum_order, momentum_sources),
            (poisson_order, poisson_sources),
        ]
    )

    resolved_step, available_steps = _resolve_summary_step(
        step,
        continuity_rows,
        particle_rows,
        momentum_rows,
        poisson_rows,
        profiling_rows,
        convergence_rows,
        memory_rows,
        step_orders=step_orders,
        selection_mode=selection_mode,
    )
    if resolved_step is None:
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key="summary",
            file_path=log_dir,
            message="No summary-capable run artifacts were found under the run log directory.",
            hint="Run the solver first, then retry summarize on a run directory that contains continuity or solver convergence logs.",
        )
        sys.exit(1)

    if step is not None and step not in set(available_steps):
        emit_structured_error(
            ERROR_CODE_CFG_INVALID_VALUE,
            key="step",
            file_path=context["run_dir"],
            message=f"Requested step {step} is not present in the available summary artifacts.",
            hint=f"Available steps include: {available_steps[:10]}{'...' if len(available_steps) > 10 else ''}",
        )
        sys.exit(1)

    continuity_step_rows = sorted(continuity_rows.get(resolved_step, []), key=lambda row: row["block"])
    continuity_summary = {"available": bool(continuity_step_rows), "blocks": continuity_step_rows}
    if continuity_step_rows:
        divergence_values = [
            abs(row["max_divergence"])
            for row in continuity_step_rows
            if row["max_divergence"] is not None
        ]
        continuity_summary["max_abs_divergence"] = max(divergence_values) if divergence_values else None
        continuity_summary["net_flux"] = continuity_step_rows[0].get("net_flux")
        continuity_summary["flux_in"] = continuity_step_rows[0].get("flux_in")
        continuity_summary["flux_out"] = continuity_step_rows[0].get("flux_out")

    momentum_step_rows = [row for _, row in sorted(momentum_rows.get(resolved_step, {}).items())]
    momentum_summary = {"available": bool(momentum_step_rows), "blocks": momentum_step_rows}

    poisson_step_rows = [row for _, row in sorted(poisson_rows.get(resolved_step, {}).items())]
    poisson_summary = {"available": bool(poisson_step_rows), "blocks": poisson_step_rows}

    particle_summary = {"available": resolved_step in particle_rows}
    if resolved_step in particle_rows:
        particle_summary.update(particle_rows[resolved_step])

    profiling_summary = {"available": resolved_step in profiling_rows}
    if resolved_step in profiling_rows:
        functions = sorted(
            profiling_rows[resolved_step],
            key=lambda row: (row.get("step_time_s") or 0.0),
            reverse=True,
        )
        profiling_summary["functions"] = functions
        profiling_summary["total_logged_step_time_s"] = sum(
            row.get("step_time_s") or 0.0 for row in functions
        )

    memory_summary = {"available": resolved_step in memory_rows}
    if resolved_step in memory_rows:
        memory_summary.update(memory_rows[resolved_step])
        memory_summary["source"] = memory_path
        memory_summary["max_process_change_mb"] = memory_meta.get("max_process_change_mb")
        memory_summary["final_event"] = memory_meta.get("final_event")
        memory_summary["final_reason"] = memory_meta.get("final_reason")
        memory_summary["selected_step"] = resolved_step
        memory_summary["step_match"] = True
    elif memory_meta.get("available"):
        latest_sample_row = memory_meta.get("latest_sample_row")
        if latest_sample_row:
            memory_summary.update(latest_sample_row)
        memory_summary.update(memory_meta)
        memory_summary["selected_step"] = resolved_step
        memory_summary["step_match"] = False

    snapshot_summary = {"available": False}
    if context["particle_console_output_freq"] and context["particle_console_output_freq"] > 0:
        snapshot_summary = _find_particle_snapshot_for_step(
            context["run_dir"],
            log_dir,
            resolved_step,
            preview_rows=max(1, snapshot_rows),
            particle_console_output_freq=context["particle_console_output_freq"],
            particle_log_interval=context["particle_log_interval"],
        )

    monitor_info = {
        "profiling_timestep_mode": context["profiling_cfg"].get("mode"),
        "profiling_timestep_file": context["profiling_cfg"].get("timestep_file"),
        "particle_console_output_frequency": context["particle_console_output_freq"],
        "particle_log_interval": context["particle_log_interval"],
    }

    return {
        "run_id": context["manifest"].get("run_id", os.path.basename(context["run_dir"])),
        "run_dir": context["run_dir"],
        "step": resolved_step,
        "selected_via": "explicit" if step is not None else ("max_step" if selection_mode == "max_step" else "latest_available"),
        "available_steps": available_steps,
        "launch_mode": context["manifest"].get("launch_mode"),
        "created_at": context["manifest"].get("created_at"),
        "monitor": monitor_info,
        "particles_configured": context["particle_count_cfg"],
        "sources": {
            "continuity_log": continuity_path if os.path.isfile(continuity_path) else None,
            "particle_metrics_log": particle_metrics_path if os.path.isfile(particle_metrics_path) else None,
            "momentum_logs": momentum_sources,
            "poisson_logs": poisson_sources,
            "profiling_timestep_csv": profiling_path if os.path.isfile(profiling_path) else None,
            "solution_convergence_log": convergence_log_path if os.path.isfile(convergence_log_path) else None,
            "runtime_memory_log": memory_path if os.path.isfile(memory_path) else None,
        },
        "continuity": continuity_summary,
        "momentum": momentum_summary,
        "poisson": poisson_summary,
        "particles": particle_summary,
        "particle_snapshot": snapshot_summary,
        "profiling": profiling_summary,
        "memory": memory_summary,
        "convergence": convergence_rows.get(resolved_step) if convergence_rows else None,
    }


def render_run_summary(payload: dict, output_format: str = "text"):
    """!
    @brief Render a run-step summary in human or JSON form.
    @param[in] payload Argument passed to `render_run_summary()`.
    @param[in] output_format Argument passed to `render_run_summary()`.
    """
    if output_format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
        return

    print("\n" + "=" * 60)
    print("                     RUN STEP SUMMARY")
    print("=" * 60)
    print(f"  Run ID         : {payload.get('run_id')}")
    print(f"  Run directory  : {os.path.relpath(payload.get('run_dir'))}")
    print(f"  Step           : {payload.get('step')} ({payload.get('selected_via')})")
    if payload.get("launch_mode"):
        print(f"  Launch mode    : {payload.get('launch_mode')}")
    if payload.get("created_at"):
        print(f"  Created at     : {payload.get('created_at')}")

    continuity = payload.get("continuity", {})
    print("\n  Continuity:")
    if continuity.get("available"):
        if continuity.get("max_abs_divergence") is not None:
            print(f"    max |div|    : {continuity['max_abs_divergence']:.6e}")
        if continuity.get("net_flux") is not None:
            print(f"    net flux     : {continuity['net_flux']:.6e}")
        for row in continuity.get("blocks", []):
            print(
                "    "
                f"block {row['block']}: div={_format_summary_float(row.get('max_divergence'))} "
                f"rhs={_format_summary_float(row.get('rhs_sum'))} location={row['max_divergence_location']}"
            )
    else:
        print("    unavailable")

    momentum = payload.get("momentum", {})
    print("\n  Momentum:")
    if momentum.get("available"):
        for row in momentum.get("blocks", []):
            status = row.get("status") or "unknown"
            accepted = row.get("accepted_count")
            rejected = row.get("rejected_count")
            counts_str = f" accepted={accepted} rejected={rejected}" if accepted is not None else ""
            cfl_in = row.get("pseudo_cfl")
            cfl_out = row.get("cfl_after")
            if cfl_in is not None and cfl_out is not None:
                cfl_str = f"cfl {cfl_in:.4f}->{cfl_out:.4f}"
            elif cfl_in is not None:
                cfl_str = f"cfl={_format_summary_float(cfl_in, '.4f')}"
            else:
                cfl_str = "cfl=n/a"
            ratio = row.get("trial_ratio")
            ratio_str = f"  ratio={_format_summary_float(ratio)}" if ratio is not None else ""
            print(f"    block {row['block']} [{status}]:{counts_str}  {cfl_str}{ratio_str}")
            print(
                f"      resid={_format_summary_float(row.get('residual_norm'))}"
                f"  delta={_format_summary_float(row.get('delta_norm'))}"
            )
    else:
        print("    unavailable")

    poisson = payload.get("poisson", {})
    print("\n  Poisson:")
    if poisson.get("available"):
        for row in poisson.get("blocks", []):
            print(
                "    "
                f"block {row['block']}: iter={row['iterations']} "
                f"true={_format_summary_float(row.get('true_norm'))} "
                f"rel={_format_summary_float(row.get('relative_norm'))}"
            )
    else:
        print("    unavailable")

    convergence = payload.get("convergence")
    print("\n  Solution Convergence:")
    if convergence is not None:
        mode = convergence.get("mode", "unknown")
        ref = convergence.get("ref")
        print(f"    mode           : {mode}  (ref={'yes' if ref else 'no'})")
        if mode in ("steady_deterministic", "transient"):
            print(f"    u_abs_l2       : {_format_summary_float(convergence.get('u_abs_l2'))}")
            print(f"    mean_speed     : {_format_summary_float(convergence.get('mean_speed'))}  drift={_format_summary_float(convergence.get('spd_abs'))}")
            print(f"    mean_ke        : {_format_summary_float(convergence.get('mean_ke'))}  drift={_format_summary_float(convergence.get('ke_abs'))}")
        elif mode == "periodic_deterministic":
            ph = convergence.get("ph")
            per = convergence.get("per")
            print(f"    phase          : {ph}/{per}")
            print(f"    u_abs_l2       : {_format_summary_float(convergence.get('u_abs_l2'))}")
            print(f"    mean_speed     : {_format_summary_float(convergence.get('mean_speed'))}  drift={_format_summary_float(convergence.get('spd_abs'))}")
            print(f"    mean_ke        : {_format_summary_float(convergence.get('mean_ke'))}  drift={_format_summary_float(convergence.get('ke_abs'))}")
        elif mode == "statistical_steady":
            print(f"    mean_speed     : {_format_summary_float(convergence.get('mean_speed'))}  win={_format_summary_float(convergence.get('spd_win'))}  win_drift={_format_summary_float(convergence.get('spd_win_abs'))}")
            print(f"    mean_ke        : {_format_summary_float(convergence.get('mean_ke'))}  win={_format_summary_float(convergence.get('ke_win'))}  win_drift={_format_summary_float(convergence.get('ke_win_abs'))}")
    else:
        print("    unavailable")

    particles = payload.get("particles", {})
    print("\n  Particles:")
    if particles.get("available"):
        loss_summary = f"lost={particles.get('lost_particles')}"
        if particles.get("lost_particles_cumulative") is not None:
            loss_summary = (
                f"lost(step/total)={particles.get('lost_particles')}/"
                f"{particles.get('lost_particles_cumulative')}"
            )
        print(
            "    "
            f"total={particles.get('total_particles')} {loss_summary} "
            f"migrated={particles.get('migrated_particles')} occupied={particles.get('occupied_cells')} "
            f"imbalance={_format_summary_float(particles.get('load_imbalance'), '.2f')}"
        )
    else:
        print("    unavailable")

    memory = payload.get("memory", {})
    print("\n  Runtime Memory:")
    if memory.get("available"):
        if memory.get("source"):
            print(f"    source        : {os.path.relpath(memory.get('source'))}")
        if memory.get("step") is not None and not memory.get("step_match", True):
            print(f"    memory step   : {memory.get('step')} (latest memory row; selected step {memory.get('selected_step')} has no row yet)")
        if memory.get("event"):
            print(f"    event         : {memory.get('event')}  reason={memory.get('reason', '-')}")
        print(f"    process max   : {_format_summary_float(memory.get('process_current_mb_max'), '.3f')} MB current, {_format_summary_float(memory.get('process_peak_mb_max'), '.3f')} MB peak")
        print(f"    PETSc max     : {_format_summary_float(memory.get('petsc_allocated_mb_max'), '.3f')} MB allocated, {_format_summary_float(memory.get('petsc_peak_allocated_mb_max'), '.3f')} MB peak")
        print(f"    max change    : {_format_summary_float(memory.get('max_process_change_mb'), '.3f')} MB")
        if memory.get("final_reason"):
            print(f"    final reason  : {memory.get('final_reason')}")
    else:
        print("    unavailable")

    snapshot = payload.get("particle_snapshot", {})
    if snapshot.get("available"):
        print("\n  Particle Snapshot (sampled):")
        print(f"    source        : {os.path.relpath(snapshot.get('source'))}")
        cadence = snapshot.get("cadence", {})
        print(
            "    "
            f"cadence       : every {cadence.get('particle_console_output_frequency', 'n/a')} steps, "
            f"row interval {cadence.get('particle_log_interval', 'n/a')}"
        )
        print(f"    sampled rows  : {snapshot.get('sampled_rows')}")
        speed = snapshot.get("speed", {})
        if speed:
            print(
                "    "
                f"sampled speeds: min={_format_summary_float(speed.get('min'))} "
                f"mean={_format_summary_float(speed.get('mean'))} "
                f"max={_format_summary_float(speed.get('max'))} "
                f"std={_format_summary_float(speed.get('std'))} "
                f"stagnant(<1e-6)={speed.get('stagnant_count', 0)}"
            )
        bounds = snapshot.get("position_bounds", {})
        centroid = snapshot.get("position_centroid")
        if bounds:
            bound_parts = []
            for axis in ("x", "y", "z"):
                if axis in bounds:
                    bound_parts.append(
                        f"{axis}=[{_format_summary_float(bounds[axis][0])}, { _format_summary_float(bounds[axis][1])}]"
                    )
            print(f"    sampled bounds: {' '.join(bound_parts)}")
        if centroid:
            print(
                "    "
                f"sampled center: ({_format_summary_float(centroid[0])}, "
                f"{_format_summary_float(centroid[1])}, {_format_summary_float(centroid[2])})"
            )
        distribution = snapshot.get("sampled_distribution", {})
        if distribution:
            print(
                "    "
                f"sampled spread: unique_cells={distribution.get('unique_cells', 'n/a')} "
                f"duplicate_cells={distribution.get('duplicate_cells', 'n/a')} "
                f"unique_pids={distribution.get('unique_pids', 'n/a')} "
                f"ranks={distribution.get('rank_counts', {})}"
            )
        weights = snapshot.get("weights", {})
        if weights:
            weight_parts = []
            for component, summary in sorted(weights.items()):
                weight_parts.append(
                    f"{component}[min/max]=[{_format_summary_float(summary.get('min'))}, { _format_summary_float(summary.get('max'))}]"
                )
            print(f"    sampled weights: {' '.join(weight_parts)}")
        checks = snapshot.get("checks", {})
        if checks:
            print(
                "    "
                f"checks        : duplicate_pid={checks.get('duplicate_pid_count', 0)} "
                f"nan={checks.get('nan_count', 0)} inf={checks.get('inf_count', 0)} "
                f"zero_weight={checks.get('zero_weight_count', 0)} "
                f"negative_weight={checks.get('negative_weight_count', 0)}"
            )
        top_speeds = snapshot.get("top_speeds", [])
        if top_speeds:
            summary = ", ".join(
                f"pid={row.get('pid')} {_format_summary_float(row.get('speed'))}"
                for row in top_speeds
            )
            print(f"    top speeds    : {summary}")
        delta_summary = snapshot.get("delta_from_previous_snapshot", {})
        if delta_summary.get("available"):
            print(
                "    "
                f"vs prev snap  : step={delta_summary.get('previous_step')} "
                f"matched_pids={delta_summary.get('matched_pids')} "
                f"mean_disp={_format_summary_float(delta_summary.get('mean_displacement'))} "
                f"max_disp={_format_summary_float(delta_summary.get('max_displacement'))} "
                f"rank_moves={delta_summary.get('rank_migrations')} "
                f"cell_changes={delta_summary.get('cell_changes')} "
                f"new={delta_summary.get('new_count')} gone={delta_summary.get('gone_count')}"
            )
        print("    preview rows  :")
        for row in snapshot.get("preview_rows", []):
            print(
                "    "
                f"pid={row.get('pid')} rank={row.get('rank')} "
                f"cell={row.get('cell')} pos={row.get('position')} vel={row.get('velocity')}"
            )

    profiling = payload.get("profiling", {})
    print("\n  Profiling:")
    if profiling.get("available"):
        print(f"    total logged step time: {profiling.get('total_logged_step_time_s', 0.0):.6f}s")
        for row in profiling.get("functions", [])[:5]:
            print(
                "    "
                f"{row.get('function')}: calls={row.get('calls')} "
                f"time={_format_summary_float(row.get('step_time_s'), '.6f', '0.000000')}s"
            )
    else:
        print("    unavailable")
    print("=" * 60)


_CONFIG_SUMMARY_WIDTH = 78


def _summary_display_value(value) -> str:
    """!
    @brief Format one configuration-summary value for compact text output.
    @param[in] value Value to format.
    @return Compact human-readable value.
    """
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "enabled" if value else "disabled"
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (list, tuple)):
        return ", ".join(_summary_display_value(item) for item in value) if value else "none"
    if isinstance(value, dict):
        if not value:
            return "none"
        return ", ".join(f"{key}={_summary_display_value(item)}" for key, item in value.items())
    return str(value)


def _print_config_header(title: str, subtitle: "str | None" = None):
    """!
    @brief Print a strong dashboard-style configuration summary header.
    @param[in] title Section title.
    @param[in] subtitle Optional one-line section subtitle.
    """
    print("\n" + "=" * _CONFIG_SUMMARY_WIDTH)
    print(f"{title:^78}")
    if subtitle:
        print(f"{subtitle:^78}")
    print("=" * _CONFIG_SUMMARY_WIDTH)


def _print_config_group(title: str, rows: list):
    """!
    @brief Print an aligned configuration-summary field group.
    @param[in] title Group title.
    @param[in] rows Sequence of `(label, value)` pairs.
    """
    visible_rows = [(label, value) for label, value in rows if value is not None]
    if not visible_rows:
        return
    print(f"\n  {title}")
    print(f"  {'-' * (len(title) + 1)}")
    for label, value in visible_rows:
        print(f"    {label:<32} {_summary_display_value(value)}")


def _flatten_summary_mapping(mapping: dict, prefix: str = "") -> list:
    """!
    @brief Flatten nested summary mappings into readable dotted field rows.
    @param[in] mapping Mapping to flatten.
    @param[in] prefix Optional parent-field prefix.
    @return Sequence of `(field, value)` pairs.
    """
    rows = []
    for key, value in mapping.items():
        label = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict) and value:
            rows.extend(_flatten_summary_mapping(value, label))
        else:
            rows.append((label, value))
    return rows


def _render_run_overview_text(summary: dict):
    """!
    @brief Render run metadata as a compact dashboard.
    @param[in] summary Curated run overview mapping.
    """
    _print_config_header("RUN OVERVIEW", summary.get("run_id"))
    _print_config_group(
        "Identity",
        [
            ("Run directory", os.path.relpath(summary.get("run_dir")) if summary.get("run_dir") else None),
            ("Created", summary.get("created_at")),
            ("Launch mode", summary.get("launch_mode")),
            ("Git commit", summary.get("git_commit")),
        ],
    )
    _print_config_group(
        "Execution",
        [
            ("Solver MPI processes", summary.get("solver_num_procs")),
            ("Post MPI processes", summary.get("post_num_procs")),
            ("Stages requested", summary.get("stages_requested")),
            ("Stages ready/completed", summary.get("stages_completed_or_submitted")),
        ],
    )


def _render_case_summary_text(summary: dict):
    """!
    @brief Render the case summary as a glanceable simulation dashboard.
    @param[in] summary Curated case configuration mapping.
    """
    run = summary.get("run_control", {})
    props = summary.get("properties", {})
    grid = summary.get("grid", {})
    domain = summary.get("domain", {})
    physics = summary.get("physics", {})
    subtitle = (
        f"{domain.get('dimensionality', '-')} | {domain.get('blocks', '-')} block(s) | "
        f"Re={_summary_display_value(props.get('reynolds_number'))}"
    )
    _print_config_header("CASE SUMMARY", subtitle)
    _print_config_group(
        "Simulation",
        [
            ("Step range", f"{run.get('start_step')} -> {run.get('end_step')} ({run.get('total_steps')} steps)"),
            ("Physical timestep", run.get("dt_physical")),
            ("Nondimensional timestep", run.get("dt_nondimensional")),
            ("Physical duration", run.get("duration_physical")),
            ("Initial conditions", props.get("initial_conditions")),
        ],
    )
    _print_config_group(
        "Fluid And Scaling",
        [
            ("Reynolds number", props.get("reynolds_number")),
            ("Reference length", props.get("length_ref")),
            ("Reference velocity", props.get("velocity_ref")),
            ("Density", props.get("density")),
            ("Viscosity", props.get("viscosity")),
        ],
    )
    _print_config_group(
        "Domain And Grid",
        [
            ("Grid mode", grid.get("mode")),
            ("Blocks", domain.get("blocks")),
            ("Dimensionality", domain.get("dimensionality")),
            ("Periodic axes", domain.get("periodic")),
            ("MPI grid layout", grid.get("processor_layout")),
            ("Grid source", grid.get("source_file")),
        ],
    )
    if grid.get("programmatic_settings"):
        _print_config_group("Programmatic Grid", _flatten_summary_mapping(grid["programmatic_settings"]))
    _print_config_group(
        "Physics",
        [
            ("Particles", physics.get("particles")),
            ("FSI", physics.get("fsi")),
            ("Turbulence", physics.get("turbulence")),
            ("Statistics", physics.get("statistics")),
        ],
    )
    boundary_blocks = summary.get("boundary_conditions", [])
    if boundary_blocks:
        print("\n  Boundary Conditions")
        print("  --------------------")
        print(f"    {'Block':<7} {'Face':<8} {'Type':<12} Handler")
        print(f"    {'-' * 7} {'-' * 8} {'-' * 12} {'-' * 20}")
        for block in boundary_blocks:
            for face in block.get("faces", []):
                print(
                    f"    {block.get('block', '-')!s:<7} {face.get('face', '-'):<8} "
                    f"{face.get('type', '-'):<12} {face.get('handler', '-')}"
                )


def _render_solver_summary_text(summary: dict):
    """!
    @brief Render the solver summary as a glanceable numerical-method dashboard.
    @param[in] summary Curated solver configuration mapping.
    """
    momentum = summary.get("momentum", {})
    poisson = summary.get("poisson", {})
    operation = summary.get("operation_mode", {})
    subtitle = (
        f"Field: {operation.get('eulerian_field_source', '-')} | "
        f"Momentum: {momentum.get('type', '-')} | Poisson: {poisson.get('method', '-')}"
    )
    _print_config_header("SOLVER SUMMARY", subtitle)
    _print_config_group("Operation", _flatten_summary_mapping(operation))
    _print_config_group(
        "Primary Methods",
        [
            ("Momentum solver", momentum.get("type")),
            ("Central differencing", momentum.get("central_diff")),
            ("Poisson method", poisson.get("method")),
            ("Interpolation", summary.get("interpolation")),
            ("Convergence mode", summary.get("solution_convergence", {}).get("mode")),
        ],
    )
    _print_config_group("Momentum Tolerances", _flatten_summary_mapping(momentum.get("tolerances", {})))
    _print_config_group("Momentum Controls", _flatten_summary_mapping(momentum.get("controls", {})))
    _print_config_group("Poisson Configuration", _flatten_summary_mapping(poisson))
    _print_config_group("Solution Convergence", _flatten_summary_mapping(summary.get("solution_convergence", {})))
    _print_config_group("Scalar Transport", _flatten_summary_mapping(summary.get("scalar_transport", {})))
    _print_config_group("Verification Sources", _flatten_summary_mapping(summary.get("verification", {})))
    passthrough = summary.get("petsc_passthrough", {})
    _print_config_group(
        "Advanced PETSc Options",
        [("Option count", passthrough.get("count")), ("Option names", passthrough.get("options"))],
    )


def _render_monitor_summary_text(summary: dict):
    """!
    @brief Render the monitor summary as a glanceable observability dashboard.
    @param[in] summary Curated monitor configuration mapping.
    """
    logging_cfg = summary.get("logging", {})
    profiling = summary.get("profiling", {})
    diagnostics = summary.get("diagnostics", {})
    io_cfg = summary.get("io", {})
    memory_log = diagnostics.get("runtime_memory_log", {})
    subtitle = (
        f"Verbosity: {logging_cfg.get('verbosity', '-')} | Profiling: {profiling.get('mode', '-')} | "
        f"Output every {_summary_display_value(io_cfg.get('data_output_frequency'))} steps"
    )
    _print_config_header("MONITOR SUMMARY", subtitle)
    _print_config_group(
        "Logging",
        [
            ("Verbosity", logging_cfg.get("verbosity")),
            ("Enabled functions", logging_cfg.get("enabled_functions")),
        ],
    )
    _print_config_group("Profiling", _flatten_summary_mapping(profiling))
    _print_config_group(
        "Output Cadence",
        [
            ("Field output", io_cfg.get("data_output_frequency")),
            ("Particle snapshots", io_cfg.get("particle_console_output_frequency")),
            ("Particle row interval", io_cfg.get("particle_log_interval")),
        ],
    )
    _print_config_group("Output Directories", _flatten_summary_mapping(io_cfg.get("directories", {})))
    _print_config_group(
        "Diagnostics",
        [
            ("Enabled PETSc diagnostics", diagnostics.get("enabled_petsc")),
            ("Runtime memory log", memory_log.get("enabled")),
            ("Runtime memory file", memory_log.get("file")),
        ],
    )
    _print_config_group("PETSc Diagnostic Settings", _flatten_summary_mapping(diagnostics.get("petsc", {})))
    solver_monitoring = summary.get("solver_monitoring", {})
    _print_config_group(
        "Solver Monitoring",
        [
            ("Enabled flags", solver_monitoring.get("enabled_flags")),
            ("All flags", solver_monitoring.get("flags")),
        ],
    )


def render_selected_summary(payload: dict, output_format: str = "text"):
    """!
    @brief Render selected timestep-independent config views and optional health.
    @param[in] payload Combined selected summary payload.
    @param[in] output_format Output format.
    """
    if output_format == "json":
        json_payload = {key: value for key, value in payload.items() if key != "_health_requested"}
        print(json.dumps(json_payload, indent=2, sort_keys=True))
        return

    if payload.get("run_overview") is not None:
        _render_run_overview_text(payload["run_overview"])
    renderers = {
        "case": _render_case_summary_text,
        "solver": _render_solver_summary_text,
        "monitor": _render_monitor_summary_text,
    }
    for key in ("case", "solver", "monitor"):
        if key in payload.get("configuration", {}):
            renderers[key](payload["configuration"][key])
    if payload.get("_health_requested"):
        health_payload = {key: value for key, value in payload.items() if key not in {"run_overview", "configuration", "_health_requested"}}
        render_run_summary(health_payload, output_format="text")


_SUMMARY_PLOT_LOG_SCALE_FIELDS = {
    "delta_norm", "delta_rel", "residual_norm", "residual_rel",
    "unpreconditioned_norm", "true_norm", "relative_norm",
    "u_abs_l2", "u_rel_l2", "p_abs_l2", "p_rel_l2",
}


def _append_summary_plot_record(records: list, source: str, step, line: str, values: dict, source_path: str, segment: int = 0):
    """!
    @brief Append one numeric append-ordered record for summarize plotting.
    @param[out] records Destination record list.
    @param[in] source Qualified source prefix.
    @param[in] step Logged timestep.
    @param[in] line Human-readable line identity.
    @param[in] values Candidate field mapping.
    @param[in] source_path Source artifact path.
    @param[in] segment Zero-based continuation segment within the source artifact.
    """
    numeric = {
        key: value
        for key, value in values.items()
        if isinstance(value, (int, float)) and not isinstance(value, bool)
    }
    if step is not None and numeric:
        records.append({
            "source": source,
            "step": int(step),
            "line": line,
            "values": numeric,
            "source_path": source_path,
            "segment": int(segment),
        })


def _is_summary_plot_continuation_marker(line: str) -> bool:
    """!
    @brief Return whether a log line starts a new continuation segment.
    @param[in] line Candidate raw or stripped log line.
    @return True for the shared continuation marker syntax.
    """
    return bool(re.match(r"^\s*#?\s*=*\s*Continuation from step\s+\d+", line, re.IGNORECASE))


def _collect_summary_plot_records(context: dict) -> list:
    """!
    @brief Collect append-ordered numeric records from summarize-supported scalar logs.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @return Append-ordered plot record list.
    """
    records = []
    log_dir = context["log_dir"]

    continuity_path = os.path.join(log_dir, "Continuity_Metrics.log")
    if os.path.isfile(continuity_path):
        segment = 0
        with open(continuity_path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                if _is_summary_plot_continuation_marker(raw_line):
                    segment += 1
                    continue
                parts = [part.strip() for part in raw_line.split("|")]
                if len(parts) < 8:
                    continue
                step, block = _parse_int_loose(parts[0]), _parse_int_loose(parts[1])
                _append_summary_plot_record(
                    records, "continuity", step, f"block {block}",
                    {
                        "max_divergence": _parse_float_loose(parts[2]),
                        "rhs_sum": _parse_float_loose(parts[4]),
                        "flux_in": _parse_float_loose(parts[5]),
                        "flux_out": _parse_float_loose(parts[6]),
                        "net_flux": _parse_float_loose(parts[7]),
                    },
                    continuity_path, segment,
                )

    particle_path = os.path.join(log_dir, "Particle_Metrics.log")
    if os.path.isfile(particle_path):
        segment = 0
        with open(particle_path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                if _is_summary_plot_continuation_marker(raw_line):
                    segment += 1
                    continue
                parts = [part.strip() for part in raw_line.split("|")]
                if len(parts) < 8:
                    continue
                step = _parse_int_loose(parts[1])
                offset = 1 if len(parts) >= 9 else 0
                _append_summary_plot_record(
                    records, "particles", step, "particles",
                    {
                        "total_particles": _parse_int_loose(parts[2]),
                        "lost_particles": _parse_int_loose(parts[3]),
                        "lost_particles_cumulative": _parse_int_loose(parts[4]) if offset else None,
                        "migrated_particles": _parse_int_loose(parts[4 + offset]),
                        "occupied_cells": _parse_int_loose(parts[5 + offset]),
                        "load_imbalance": _parse_float_loose(parts[6 + offset]),
                        "migration_passes": _parse_int_loose(parts[7 + offset]),
                    },
                    particle_path, segment,
                )

    momentum_regex = re.compile(
        r"Step:\s*(?P<step>\d+)\s*\|\s*PseudoIter\(k\):\s*(?P<pseudo_iter>\d+)\s*\|"
        r"\s*Pseudo-cfl:\s*(?P<pseudo_cfl>[-+0-9.eE]+)\s*\|\s*\|dUk\|:\s*(?P<delta>[-+0-9.eE]+)\s*\|"
        r"\s*\|dUk\|/\|dU0\|:\s*(?P<delta_rel>[-+0-9.eE]+)\s*\|\s*\|Rk\|:\s*(?P<resid>[-+0-9.eE]+)\s*\|"
        r"\s*\|Rk\|/\|R0\|:\s*(?P<resid_rel>[-+0-9.eE]+)"
        r"(?:\s*\|\s*trial_ratio:\s*(?P<trial_ratio>[-+0-9.eE]+)\s*\|\s*status:\s*(?P<status>\w+)\s*\|\s*cfl_after:\s*(?P<cfl_after>[-+0-9.eE]+))?"
    )
    for path in sorted(glob.glob(os.path.join(log_dir, "Momentum_Solver_Convergence_History_Block_*.log"))):
        block_match = re.search(r"Block_(\d+)\.log$", path)
        if not block_match:
            continue
        segment = 0
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                if _is_summary_plot_continuation_marker(raw_line):
                    segment += 1
                    continue
                match = momentum_regex.search(raw_line)
                if match:
                    _append_summary_plot_record(
                        records, "momentum", int(match.group("step")), f"block {block_match.group(1)}",
                        {
                            "pseudo_iterations": int(match.group("pseudo_iter")),
                            "pseudo_cfl": float(match.group("pseudo_cfl")),
                            "delta_norm": float(match.group("delta")),
                            "delta_rel": float(match.group("delta_rel")),
                            "residual_norm": float(match.group("resid")),
                            "residual_rel": float(match.group("resid_rel")),
                            "trial_ratio": _parse_float_loose(match.group("trial_ratio")),
                            "cfl_after": _parse_float_loose(match.group("cfl_after")),
                        },
                        path, segment,
                    )

    poisson_regex = re.compile(
        r"ts:\s*(?P<step>\d+)\s*\|\s*block:\s*(?P<block>\d+)\s*\|\s*iter:\s*(?P<iter>\d+)\s*\|"
        r"\s*Unprecond Norm:\s*(?P<unpre>[-+0-9.eE]+)\s*\|\s*True Norm:\s*(?P<true>[-+0-9.eE]+)"
        r"(?:\s*\|\s*Rel Norm:\s*(?P<rel>[-+0-9.eE]+))?"
    )
    for path in sorted(glob.glob(os.path.join(log_dir, "Poisson_Solver_Convergence_History_Block_*.log"))):
        segment = 0
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                if _is_summary_plot_continuation_marker(raw_line):
                    segment += 1
                    continue
                match = poisson_regex.search(raw_line)
                if match:
                    _append_summary_plot_record(
                        records, "poisson", int(match.group("step")), f"block {match.group('block')}",
                        {
                            "iterations": int(match.group("iter")),
                            "unpreconditioned_norm": float(match.group("unpre")),
                            "true_norm": float(match.group("true")),
                            "relative_norm": _parse_float_loose(match.group("rel")),
                        },
                        path, segment,
                    )

    profiling_path = os.path.join(log_dir, context["profiling_cfg"].get("timestep_file", "Profiling_Timestep_Summary.csv"))
    if os.path.isfile(profiling_path):
        segment = 0
        columns = None
        with open(profiling_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            for raw_line in f:
                if _is_summary_plot_continuation_marker(raw_line):
                    segment += 1
                    continue
                values = next(csv.reader([raw_line]))
                if not values:
                    continue
                if columns is None:
                    columns = values
                    continue
                row = dict(zip(columns, values))
                _append_summary_plot_record(
                    records, "profiling", _parse_int_loose(row.get("step")), row.get("function") or "unknown",
                    {"calls": _parse_int_loose(row.get("calls")), "step_time_s": _parse_float_loose(row.get("step_time_s"))},
                    profiling_path, segment,
                )

    diagnostics = resolve_diagnostics_config(context["monitor_cfg"])
    memory_path = os.path.join(log_dir, diagnostics["runtime_memory_log"].get("file", "Runtime_Memory.log"))
    if os.path.isfile(memory_path):
        segment = 0
        with open(memory_path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                if _is_summary_plot_continuation_marker(raw_line):
                    segment += 1
                    continue
                parts = raw_line.split()
                if len(parts) >= 8 and parts[1] in {"Step", "Post"}:
                    _append_summary_plot_record(
                        records, "memory", _parse_int_loose(parts[0]), "memory",
                        {
                            "process_current_mb_max": _parse_float_loose(parts[2]),
                            "process_peak_mb_max": _parse_float_loose(parts[3]),
                            "petsc_allocated_mb_max": _parse_float_loose(parts[4]),
                            "petsc_peak_allocated_mb_max": _parse_float_loose(parts[5]),
                            "process_change_mb_max": _parse_float_loose(parts[6]),
                        },
                        memory_path, segment,
                    )

    convergence_path = os.path.join(log_dir, "solution_convergence.log")
    if os.path.isfile(convergence_path):
        columns = None
        segment = 0
        with open(convergence_path, "r", encoding="utf-8", errors="replace") as f:
            for raw_line in f:
                line = raw_line.strip()
                if _is_summary_plot_continuation_marker(line):
                    segment += 1
                    continue
                if not line or line.startswith(("=", "-")):
                    continue
                if columns is None:
                    columns = [part.strip() for part in raw_line.split("|")]
                    continue
                parts = [part.strip() for part in raw_line.split("|")]
                step = _parse_int_loose(parts[0])
                values = {name: _parse_float_loose(value) for name, value in zip(columns, parts) if name not in {"step", "mode", "ref"}}
                _append_summary_plot_record(records, "convergence", step, "convergence", values, convergence_path, segment)
    return records


def _build_summary_plot_catalog(records: list) -> list:
    """!
    @brief Build available qualified-series metadata from plot records.
    @param[in] records Append-ordered plot record list.
    @return Available series catalog.
    """
    catalog = {}
    for record in records:
        for field in record["values"]:
            name = f"{record['source']}.{field}"
            item = catalog.setdefault(name, {"series": name, "lines": {}, "source_paths": set(), "sample_count": 0})
            item["lines"][record["line"]] = item["lines"].get(record["line"], 0) + 1
            item["source_paths"].add(record["source_path"])
            item["sample_count"] += 1
    return [
        {
            **item,
            "lines": [{"label": label, "sample_count": count} for label, count in sorted(item["lines"].items())],
            "source_paths": sorted(item["source_paths"]),
        }
        for _, item in sorted(catalog.items())
    ]


def _build_summary_plot_request(context: dict, records: list, series: str, last_n: "int | None", linear_y: bool, output_path: "str | None") -> dict:
    """!
    @brief Build one normalized plot.gen request from collected summarize records.
    @param[in] context Summary context returned by `_build_summary_context()`.
    @param[in] records Append-ordered plot record list.
    @param[in] series Qualified series name.
    @param[in] last_n Optional last-N records per plotted line.
    @param[in] linear_y Whether to force linear scaling.
    @param[in] output_path Optional explicit output path.
    @return Versioned normalized plot request.
    """
    source, separator, field = series.partition(".")
    if not separator:
        raise ValueError("plot series must be qualified as '<source>.<field>'")
    matching = [record for record in records if record["source"] == source and field in record["values"]]
    if not matching:
        raise ValueError(f"Plot series '{series}' is unavailable. Use --list-plot-series to inspect available series.")
    latest_segments = {}
    for record in matching:
        source_path = record["source_path"]
        latest_segments[source_path] = max(latest_segments.get(source_path, 0), record.get("segment", 0))
    matching = [
        record for record in matching
        if record.get("segment", 0) == latest_segments[record["source_path"]]
    ]
    grouped = {}
    for record in matching:
        grouped.setdefault(record["line"], []).append([record["step"], record["values"][field]])
    if last_n is not None:
        grouped = {label: points[-last_n:] for label, points in grouped.items()}
    all_values = [point[1] for points in grouped.values() for point in points]
    use_log = not linear_y and field in _SUMMARY_PLOT_LOG_SCALE_FIELDS and all(value > 0 for value in all_values)
    window_token = f"last-{last_n}" if last_n is not None else "full"
    safe_series = re.sub(r"[^A-Za-z0-9_.-]+", "_", series)
    fallback = os.path.join(context["run_dir"], "summary", "plots", f"{safe_series}_{window_token}.png")
    return {
        "schema_version": 1,
        "plot_type": "time_history",
        "series": series,
        "title": f"{series} time history",
        "x_label": "Timestep",
        "y_label": series,
        "y_scale": "log" if use_log else "linear",
        "window": {"mode": "last" if last_n is not None else "full", "last": last_n},
        "lines": [{"label": label, "points": points} for label, points in sorted(grouped.items())],
        "output_path": os.path.abspath(output_path) if output_path else None,
        "fallback_output_path": fallback,
    }


def _render_summary_plot_catalog(catalog: list, output_format: str):
    """!
    @brief Render available summarize plot-series metadata.
    @param[in] catalog Available series catalog.
    @param[in] output_format Text or JSON output format.
    """
    if output_format == "json":
        print(json.dumps({"available_series": catalog}, indent=2, sort_keys=True))
        return
    print("\nAVAILABLE TIME-HISTORY SERIES")
    print("=" * 78)
    for item in catalog:
        labels = ", ".join(line["label"] for line in item["lines"])
        print(f"  {item['series']:<42} samples={item['sample_count']:<5} lines={labels}")
        print(f"    source: {', '.join(os.path.relpath(path) for path in item['source_paths'])}")


def _invoke_plot_gen(request: dict):
    """!
    @brief Invoke standalone plot.gen with one normalized request over stdin.
    @param[in] request Versioned normalized plot request.
    """
    plotgen_path = os.path.join(GENERATORS_PATH, "plot.gen")
    if not os.path.isfile(plotgen_path):
        raise ValueError(f"plot.gen script not found: {plotgen_path}")
    result = subprocess.run(
        [sys.executable, plotgen_path, "--input", "-"],
        input=json.dumps(request),
        text=True,
        capture_output=True,
        check=False,
    )
    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0:
        details = (result.stderr or result.stdout or "unknown plotting error").strip()
        if result.returncode == 3:
            raise PlotDependencyError(details)
        raise ValueError(f"plot.gen failed with exit code {result.returncode}: {details}")


def summarize_workflow(args):
    """!
    @brief Build and render a read-only health summary for a run step.
    @param[in] args Command-line style argument list supplied to the function.
    """
    if args.step is not None and args.step < 0:
        fail_cli_usage("--step must be a non-negative integer.")
    if args.snapshot_rows < 1:
        fail_cli_usage("--snapshot-rows must be at least 1.")
    plot_series = getattr(args, "plot_series", None)
    list_plot_series = bool(getattr(args, "list_plot_series", False))
    last_n = getattr(args, "last_n", None)
    plot_output = getattr(args, "plot_output", None)
    linear_y = bool(getattr(args, "linear_y", False))
    plot_mode = bool(plot_series or list_plot_series)
    existing_selectors = any(
        [
            getattr(args, "overview", False),
            getattr(args, "case", False),
            getattr(args, "solver", False),
            getattr(args, "monitor", False),
            args.step is not None,
            getattr(args, "latest", False),
            getattr(args, "max_step", False),
        ]
    )
    if plot_mode and existing_selectors:
        fail_cli_usage("Plot discovery and --plot cannot be combined with config or selected-step selectors.")
    if not plot_series and (last_n is not None or plot_output or linear_y):
        fail_cli_usage("--last, --plot-output, and --linear-y require --plot.")
    if last_n is not None and last_n < 1:
        fail_cli_usage("--last must be a positive integer.")
    if plot_series and args.output_format == "json":
        fail_cli_usage("--plot does not support --format json; use --list-plot-series --format json for structured discovery.")
    if plot_mode:
        context = _build_summary_context(args.run_dir)
        try:
            records = _collect_summary_plot_records(context)
            catalog = _build_summary_plot_catalog(records)
            if list_plot_series:
                if not catalog:
                    raise ValueError("No plottable scalar histories were found in the run logs.")
                _render_summary_plot_catalog(catalog, args.output_format)
                return
            request = _build_summary_plot_request(context, records, plot_series, last_n, linear_y, plot_output)
            _invoke_plot_gen(request)
            return
        except PlotDependencyError as exc:
            emit_structured_error(
                ERROR_CODE_DEPENDENCY_MISSING,
                key="plotting",
                file_path=sys.executable,
                message=str(exc),
            )
            sys.exit(1)
        except ValueError as exc:
            emit_structured_error(
                ERROR_CODE_CFG_INVALID_VALUE,
                key="plot",
                file_path=context["log_dir"],
                message=str(exc),
            )
            sys.exit(1)

    selected_configs = {
        name
        for name in ("case", "solver", "monitor")
        if bool(getattr(args, name, False))
    }
    if getattr(args, "overview", False):
        selected_configs.update({"case", "solver", "monitor"})
    explicit_health = args.step is not None or bool(getattr(args, "latest", False)) or bool(getattr(args, "max_step", False))
    health_requested = explicit_health or (not selected_configs and not getattr(args, "overview", False))

    context = None
    combined = {}
    if selected_configs or getattr(args, "overview", False):
        context = _build_summary_context(args.run_dir)
        if getattr(args, "overview", False):
            combined["run_overview"] = _build_run_overview(context)
        combined["configuration"] = {}
        builders = {
            "case": _build_case_overview,
            "solver": _build_solver_overview,
            "monitor": _build_monitor_overview,
        }
        for name in ("case", "solver", "monitor"):
            if name in selected_configs:
                try:
                    combined["configuration"][name] = builders[name](context)
                except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
                    emit_structured_error(
                        ERROR_CODE_CFG_INVALID_VALUE,
                        key=name,
                        file_path=context["config_paths"][name],
                        message=f"Could not summarize copied {name}.yml: {exc}",
                    )
                    sys.exit(1)

    if not health_requested:
        combined["_health_requested"] = False
        render_selected_summary(combined, output_format=args.output_format)
        return

    requested_step = args.step
    if requested_step is None and getattr(args, "latest", False):
        requested_step = None
    selection_mode = "max_step" if getattr(args, "max_step", False) else "latest"
    health_payload = build_run_summary_payload(
        args.run_dir,
        step=requested_step,
        snapshot_rows=args.snapshot_rows,
        selection_mode=selection_mode,
    )
    if not combined:
        render_run_summary(health_payload, output_format=args.output_format)
        return
    combined = {**health_payload, **combined, "_health_requested": True}
    render_selected_summary(combined, output_format=args.output_format)


def _resolve_submission_target(run_dir: str = None, study_dir: str = None) -> dict:
    """!
    @brief Resolve a run/study submission target from explicit directory flags.
    @param[in] run_dir Argument passed to `_resolve_submission_target()`.
    @param[in] study_dir Argument passed to `_resolve_submission_target()`.
    @return Value returned by `_resolve_submission_target()`.
    """
    has_run_dir = bool(run_dir)
    has_study_dir = bool(study_dir)
    if has_run_dir == has_study_dir:
        fail_cli_usage("submit requires exactly one of --run-dir or --study-dir.")

    target_kind = "run" if has_run_dir else "study"
    target_key = "run_dir" if target_kind == "run" else "study_dir"
    root_dir = os.path.abspath(run_dir if has_run_dir else study_dir)
    if not os.path.isdir(root_dir):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key=target_key,
            file_path=root_dir,
            message=f"{'Run' if target_kind == 'run' else 'Study'} directory not found.",
        )
        sys.exit(1)

    scheduler_dir = os.path.join(root_dir, "scheduler")
    submission_path = os.path.join(scheduler_dir, "submission.json")
    submission_meta = _read_json_if_exists(submission_path)
    if not isinstance(submission_meta, dict):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key="scheduler.submission",
            file_path=submission_path,
            message="Target directory does not contain scheduler submission metadata.",
            hint="Use a Slurm-staged run/study directory with scheduler/submission.json, or submit the script manually.",
        )
        sys.exit(1)

    launch_mode = str(submission_meta.get("launch_mode", "")).lower()
    if launch_mode == "local" and target_kind != "run":
        emit_structured_error(
            ERROR_CODE_CFG_INCONSISTENT_COMBO,
            key="scheduler.launch_mode",
            file_path=submission_path,
            message="Local staged submission is supported for run directories only.",
            hint="Use --run-dir for local staged execution; study submit remains Slurm-only.",
        )
        sys.exit(1)
    if launch_mode not in {"slurm", "local"}:
        emit_structured_error(
            ERROR_CODE_CFG_INCONSISTENT_COMBO,
            key="scheduler.launch_mode",
            file_path=submission_path,
            message=f"Target launch_mode={launch_mode or 'unknown'} is not supported.",
            hint="Use a staged run/study directory with launch_mode 'slurm' or a run directory with launch_mode 'local'.",
        )
        sys.exit(1)

    if target_kind == "run":
        script_map = {
            "solve": os.path.join(scheduler_dir, "solver.sbatch"),
            "post-process": os.path.join(scheduler_dir, "post.sbatch"),
        }
        display_label = "Run directory"
        manifest_path = None
    else:
        script_map = {
            "solve": os.path.join(scheduler_dir, "solver_array.sbatch"),
            "post-process": os.path.join(scheduler_dir, "post_array.sbatch"),
        }
        display_label = "Study directory"
        manifest_path = os.path.join(root_dir, "study_manifest.json")

    return {
        "target_kind": target_kind,
        "target_key": target_key,
        "root_dir": root_dir,
        "scheduler_dir": scheduler_dir,
        "submission_path": submission_path,
        "submission_meta": submission_meta,
        "launch_mode": launch_mode,
        "script_map": script_map,
        "display_label": display_label,
        "manifest_path": manifest_path,
    }


def _get_submission_stage_metadata(target_context: dict, stage_name: str) -> dict:
    """!
    @brief Return stored metadata for one staged submission target.
    @param[in] target_context Argument passed to `_get_submission_stage_metadata()`.
    @param[in] stage_name Argument passed to `_get_submission_stage_metadata()`.
    @return Value returned by `_get_submission_stage_metadata()`.
    """
    submission_meta = target_context["submission_meta"]
    if target_context["target_kind"] == "run":
        stages = submission_meta.get("stages", {})
        if not isinstance(stages, dict):
            return {}
        stage_meta = stages.get(stage_name)
        return copy.deepcopy(stage_meta) if isinstance(stage_meta, dict) else {}

    key = "solver_array" if stage_name == "solve" else "post_array"
    stage_meta = submission_meta.get(key)
    return copy.deepcopy(stage_meta) if isinstance(stage_meta, dict) else {}


def _get_recorded_submission_stages(target_context: dict) -> list:
    """!
    @brief Return stage names explicitly recorded in scheduler submission metadata.
    @param[in] target_context Argument passed to `_get_recorded_submission_stages()`.
    @return Value returned by `_get_recorded_submission_stages()`.
    """
    submission_meta = target_context["submission_meta"]
    recorded = []
    if target_context["target_kind"] == "run":
        stages = submission_meta.get("stages", {})
        if isinstance(stages, dict):
            for stage_name in ["solve", "post-process"]:
                if isinstance(stages.get(stage_name), dict):
                    recorded.append(stage_name)
        return recorded

    if isinstance(submission_meta.get("solver_array"), dict):
        recorded.append("solve")
    if isinstance(submission_meta.get("post_array"), dict):
        recorded.append("post-process")
    return recorded


def _format_stage_list(stage_names: list) -> str:
    """!
    @brief Format a human-readable stage list for submit diagnostics.
    @param[in] stage_names Argument passed to `_format_stage_list()`.
    @return Value returned by `_format_stage_list()`.
    """
    return ", ".join(stage_names) if stage_names else "none"


def _build_submit_missing_stage_hint(target_context: dict, requested_stage: str, selected_stages: list) -> str:
    """!
    @brief Build an actionable hint for requested submit stages missing from metadata.
    @param[in] target_context Argument passed to `_build_submit_missing_stage_hint()`.
    @param[in] requested_stage Argument passed to `_build_submit_missing_stage_hint()`.
    @param[in] selected_stages Argument passed to `_build_submit_missing_stage_hint()`.
    @return Value returned by `_build_submit_missing_stage_hint()`.
    """
    recorded_stages = _get_recorded_submission_stages(target_context)
    recorded_set = set(recorded_stages)
    selected_set = set(selected_stages)
    target_flag = "--run-dir" if target_context["target_kind"] == "run" else "--study-dir"
    target_path = os.path.relpath(target_context["root_dir"])
    submit_prefix = f"picurv submit {target_flag} {target_path}"
    solve_stage_command = (
        "picurv run --solve ... --no-submit"
        if target_context["target_kind"] == "run"
        else "picurv sweep --cluster <cluster.yml> ... --no-submit"
    )
    post_stage_command = (
        "picurv run --post-process --post <post.yml> ... --no-submit"
        if target_context["target_kind"] == "run"
        else "picurv sweep --cluster <cluster.yml> ... --no-submit"
    )
    solve_post_command = (
        "picurv run --solve --post-process --post <post.yml> ... --no-submit"
        if target_context["target_kind"] == "run"
        else "picurv sweep --cluster <cluster.yml> ... --no-submit"
    )

    if requested_stage == "all":
        if recorded_set == {"solve"}:
            return (
                "--stage all requests solve and post-process, but this target records only solve. "
                f"Use `{submit_prefix} --stage solve`, or re-stage with post-processing enabled "
                f"(`{solve_post_command}`)."
            )
        if recorded_set == {"post-process"}:
            return (
                "--stage all requests solve and post-process, but this target records only post-process. "
                f"Use `{submit_prefix} --stage post-process`, or re-stage including the solve stage "
                f"(`{solve_stage_command}`)."
            )
        missing = [stage for stage in selected_stages if stage not in recorded_set]
        if missing:
            return (
                "--stage all requests solve and post-process, but submission metadata records "
                f"{_format_stage_list(recorded_stages)}. Re-stage the missing stage(s): "
                f"{_format_stage_list(missing)}."
            )

    if selected_set == {"solve"} and "solve" not in recorded_set:
        return (
            "The solve stage was requested, but submission metadata does not record a staged solve command/script. "
            f"Re-stage with `{solve_stage_command}`."
        )
    if selected_set == {"post-process"} and "post-process" not in recorded_set:
        return (
            "The post-process stage was requested, but submission metadata does not record a staged post-process command/script. "
            f"Re-stage with post-processing enabled (`{post_stage_command}`, or `{solve_post_command}`)."
        )

    return "Re-stage the requested stage(s) with picurv run/sweep --no-submit before calling picurv submit."


def _set_submission_stage_metadata(target_context: dict, stage_name: str, stage_meta: dict):
    """!
    @brief Persist one stage's metadata back into the submission payload.
    @param[in] target_context Argument passed to `_set_submission_stage_metadata()`.
    @param[in] stage_name Argument passed to `_set_submission_stage_metadata()`.
    @param[in] stage_meta Argument passed to `_set_submission_stage_metadata()`.
    """
    submission_meta = target_context["submission_meta"]
    if target_context["target_kind"] == "run":
        stages = submission_meta.get("stages")
        if not isinstance(stages, dict):
            stages = {}
            submission_meta["stages"] = stages
        stages[stage_name] = stage_meta
        return

    key = "solver_array" if stage_name == "solve" else "post_array"
    submission_meta[key] = stage_meta


def _write_submission_target_metadata(target_context: dict):
    """!
    @brief Write updated submission metadata back to disk.
    @param[in] target_context Argument passed to `_write_submission_target_metadata()`.
    """
    write_json_file(target_context["submission_path"], target_context["submission_meta"])

    manifest_path = target_context.get("manifest_path")
    if manifest_path and os.path.isfile(manifest_path):
        manifest_payload = _read_json_if_exists(manifest_path)
        if isinstance(manifest_payload, dict):
            manifest_payload["submission"] = target_context["submission_meta"]
            write_json_file(manifest_path, manifest_payload)


def submit_staged_jobs(args):
    """!
    @brief Submit previously staged Slurm artifacts from an existing run/study directory.
    @param[in] args Command-line style argument list supplied to the function.
    """
    target_context = _resolve_submission_target(
        run_dir=getattr(args, "run_dir", None),
        study_dir=getattr(args, "study_dir", None),
    )
    stage_order = ["solve", "post-process"]
    requested_stage = args.stage
    selected_stages = stage_order if requested_stage == "all" else [requested_stage]

    print(f"[INFO] {target_context['display_label']:<20}: {os.path.relpath(target_context['root_dir'])}")
    print(f"[INFO] Submission metadata  : {os.path.relpath(target_context['submission_path'])}")
    print(f"[INFO] Requested stages    : {', '.join(selected_stages)}")

    if target_context.get("launch_mode") == "local":
        submit_staged_local_run(args, target_context, selected_stages)
        return

    stage_plans = []
    solve_existing_meta = _get_submission_stage_metadata(target_context, "solve")
    solve_existing_job_id = str(solve_existing_meta.get("job_id", "")).strip()

    for stage_name in selected_stages:
        existing_meta = _get_submission_stage_metadata(target_context, stage_name)
        script_path = target_context["script_map"][stage_name]
        missing_stage_hint = _build_submit_missing_stage_hint(target_context, requested_stage, selected_stages)
        if not existing_meta:
            emit_structured_error(
                ERROR_CODE_CFG_MISSING_KEY,
                key=f"scheduler.{stage_name}.metadata",
                file_path=target_context["submission_path"],
                message=f"Submission metadata does not record stage '{stage_name}'.",
                hint=missing_stage_hint,
            )
            sys.exit(1)

        if not os.path.isfile(script_path):
            emit_structured_error(
                ERROR_CODE_CFG_FILE_NOT_FOUND,
                key=f"scheduler.{stage_name}.script",
                file_path=script_path,
                message=f"Required {stage_name} sbatch artifact is missing.",
                hint=missing_stage_hint,
            )
            sys.exit(1)

        if existing_meta.get("submitted") and not args.force:
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key=f"scheduler.{stage_name}.submitted",
                file_path=target_context["submission_path"],
                message=f"Stage '{stage_name}' is already recorded as submitted.",
                hint="Use --force to resubmit this stage intentionally.",
            )
            sys.exit(1)

        dependency = None
        if stage_name == "post-process":
            if "solve" in selected_stages:
                dependency = "__NEW_SOLVE_JOB_ID__"
            else:
                if not (solve_existing_meta.get("submitted") and solve_existing_job_id):
                    emit_structured_error(
                        ERROR_CODE_CFG_INCONSISTENT_COMBO,
                        key="scheduler.post-process.dependency",
                        file_path=target_context["submission_path"],
                        message="Post-process submission requires a recorded solve job id when solve is not being submitted in the same command.",
                        hint="Submit --stage solve or --stage all first, or use --force only after solve metadata exists.",
                    )
                    sys.exit(1)
                dependency = solve_existing_job_id

        stage_plans.append(
            {
                "stage": stage_name,
                "script": script_path,
                "dependency": dependency,
                "existing_meta": existing_meta,
            }
        )

    if args.dry_run:
        for plan in stage_plans:
            cmd = ["sbatch"]
            dependency = plan["dependency"]
            if dependency == "__NEW_SOLVE_JOB_ID__":
                cmd.append("--dependency=afterok:<new solve job id>")
            elif dependency:
                cmd.append(f"--dependency=afterok:{dependency}")
            cmd.append(plan["script"])
            print(f"[DRY-RUN] Would run: {' '.join(cmd)}")
        print("[INFO] Dry-run only. No jobs were submitted.")
        return

    latest_solve_job_id = None
    for plan in stage_plans:
        dependency = plan["dependency"]
        if dependency == "__NEW_SOLVE_JOB_ID__":
            dependency = latest_solve_job_id

        submit_info = submit_sbatch(plan["script"], dependency=dependency)
        stage_meta = copy.deepcopy(plan["existing_meta"])
        stage_meta.update(submit_info)
        stage_meta["script"] = plan["script"]
        stage_meta["submitted"] = True
        if dependency:
            stage_meta["dependency"] = f"afterok:{dependency}"
        else:
            stage_meta.pop("dependency", None)

        _set_submission_stage_metadata(target_context, plan["stage"], stage_meta)
        print(f"[SUCCESS] Submitted {plan['stage']} job: {submit_info['job_id']}")

        if plan["stage"] == "solve":
            latest_solve_job_id = submit_info["job_id"]

    _write_submission_target_metadata(target_context)


def submit_staged_local_run(args, target_context: dict, selected_stages: list):
    """!
    @brief Execute previously staged local run commands from scheduler/submission.json.
    @param[in] args Command-line style argument list supplied to the function.
    @param[in] target_context Resolved submission target context.
    @param[in] selected_stages Ordered stage names selected by the user.
    """
    if target_context["target_kind"] != "run":
        emit_structured_error(
            ERROR_CODE_CFG_INCONSISTENT_COMBO,
            key="scheduler.launch_mode",
            file_path=target_context["submission_path"],
            message="Local staged execution is supported for run directories only.",
            hint="Use --run-dir for local staged execution.",
        )
        sys.exit(1)

    stage_plans = []
    solve_existing_meta = _get_submission_stage_metadata(target_context, "solve")
    solve_already_done = bool(solve_existing_meta.get("submitted") or solve_existing_meta.get("executed"))

    for stage_name in selected_stages:
        existing_meta = _get_submission_stage_metadata(target_context, stage_name)
        command = existing_meta.get("command")
        if not isinstance(command, list) or not command:
            if existing_meta:
                hint = "Re-stage the run with picurv run --no-submit before calling picurv submit."
            else:
                hint = _build_submit_missing_stage_hint(target_context, args.stage, selected_stages)
            emit_structured_error(
                ERROR_CODE_CFG_MISSING_KEY,
                key=f"scheduler.{stage_name}.command",
                file_path=target_context["submission_path"],
                message=f"Required local command metadata for stage '{stage_name}' is missing.",
                hint=hint,
            )
            sys.exit(1)

        if existing_meta.get("submitted") and not args.force:
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key=f"scheduler.{stage_name}.submitted",
                file_path=target_context["submission_path"],
                message=f"Stage '{stage_name}' is already recorded as submitted.",
                hint="Use --force to execute this stage again intentionally.",
            )
            sys.exit(1)

        if stage_name == "post-process" and "solve" not in selected_stages and not args.force and not solve_already_done:
            emit_structured_error(
                ERROR_CODE_CFG_INCONSISTENT_COMBO,
                key="scheduler.post-process.dependency",
                file_path=target_context["submission_path"],
                message="Post-process local execution requires a recorded completed solve stage when solve is not being executed in the same command.",
                hint="Submit --stage solve or --stage all first, or use --force after confirming source data exists.",
            )
            sys.exit(1)

        log_file = existing_meta.get("log_file")
        if not isinstance(log_file, str) or not log_file.strip():
            log_file = os.path.join("scheduler", f"{os.path.basename(target_context['root_dir'])}_{stage_name}.log")

        stage_plans.append(
            {
                "stage": stage_name,
                "command": [str(token) for token in command],
                "log_file": log_file,
                "existing_meta": existing_meta,
            }
        )

    if args.dry_run:
        for plan in stage_plans:
            print(f"[DRY-RUN] Would run: {format_command_for_display(plan['command'])}")
            print(f"[DRY-RUN] Log file : {plan['log_file']}")
        print("[INFO] Dry-run only. No local commands were executed.")
        return

    monitor_cfg = None
    monitor_path = os.path.join(target_context["root_dir"], "config", "monitor.yml")
    if os.path.isfile(monitor_path):
        monitor_cfg = read_yaml_file(monitor_path)

    for plan in stage_plans:
        execute_command(plan["command"], target_context["root_dir"], plan["log_file"], monitor_cfg)
        stage_meta = copy.deepcopy(plan["existing_meta"])
        stage_meta["command"] = plan["command"]
        stage_meta["command_string"] = format_command_for_display(plan["command"])
        stage_meta["log_file"] = plan["log_file"]
        stage_meta["submitted"] = True
        stage_meta["executed"] = True
        stage_meta["completed_at"] = datetime.now().isoformat()
        _set_submission_stage_metadata(target_context, plan["stage"], stage_meta)
        print(f"[SUCCESS] Executed local {plan['stage']} stage.")

    _write_submission_target_metadata(target_context)


def cancel_run_jobs(args):
    """!
    @brief Cancel Slurm-submitted jobs for an existing run directory.
    @param[in] args Command-line style argument list supplied to the function.
    """
    run_dir = os.path.abspath(args.run_dir)
    if not os.path.isdir(run_dir):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key="run_dir",
            file_path=run_dir,
            message="Run directory not found.",
        )
        sys.exit(1)

    submission_path = os.path.join(run_dir, "scheduler", "submission.json")
    submission_meta = _read_json_if_exists(submission_path)
    if not isinstance(submission_meta, dict):
        emit_structured_error(
            ERROR_CODE_CFG_FILE_NOT_FOUND,
            key="scheduler.submission",
            file_path=submission_path,
            message="Run directory does not contain scheduler submission metadata.",
            hint="Use a Slurm-submitted run directory with scheduler/submission.json, or cancel the job manually.",
        )
        sys.exit(1)

    launch_mode = str(submission_meta.get("launch_mode", "")).lower()
    if launch_mode != "slurm":
        emit_structured_error(
            ERROR_CODE_CFG_INCONSISTENT_COMBO,
            key="scheduler.launch_mode",
            file_path=submission_path,
            message=f"Run directory launch_mode={launch_mode or 'unknown'} is not Slurm.",
            hint="picurv cancel currently supports Slurm-submitted runs only.",
        )
        sys.exit(1)

    stage_order = ["solve", "post-process"]
    requested_stage = args.stage
    selected_stages = stage_order if requested_stage == "all" else [requested_stage]
    recorded_stages = submission_meta.get("stages", {})
    if not isinstance(recorded_stages, dict):
        recorded_stages = {}

    job_to_stages = {}
    skipped = []
    for stage_name in selected_stages:
        stage_meta = recorded_stages.get(stage_name)
        if not isinstance(stage_meta, dict):
            skipped.append((stage_name, "no stage metadata recorded"))
            continue

        job_id = str(stage_meta.get("job_id", "")).strip()
        if not stage_meta.get("submitted"):
            skipped.append((stage_name, "job was generated but not submitted"))
            continue
        if not job_id:
            skipped.append((stage_name, "submitted stage is missing a recorded job id"))
            continue

        job_to_stages.setdefault(job_id, []).append(stage_name)

    if not job_to_stages:
        print(f"[INFO] Run directory       : {os.path.relpath(run_dir)}")
        print(f"[INFO] Submission metadata: {os.path.relpath(submission_path)}")
        for stage_name, reason in skipped:
            print(f"[INFO] Skipping stage '{stage_name}': {reason}")
        print("[FATAL] No submitted Slurm job IDs were found for the requested stage selection.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Run directory       : {os.path.relpath(run_dir)}")
    print(f"[INFO] Submission metadata: {os.path.relpath(submission_path)}")
    print(f"[INFO] Requested stages   : {', '.join(selected_stages)}")

    if skipped:
        for stage_name, reason in skipped:
            print(f"[INFO] Skipping stage '{stage_name}': {reason}")

    graceful = bool(getattr(args, "graceful", False))
    failures = []
    for job_id, stage_names in job_to_stages.items():
        joined_stage_names = ", ".join(stage_names)
        use_graceful_signal = graceful and "solve" in stage_names
        scancel_cmd = ["scancel"]
        if use_graceful_signal:
            scancel_cmd.append("--signal=USR1")
        scancel_cmd.append(job_id)

        if args.dry_run:
            print(f"[DRY-RUN] Would run: {' '.join(scancel_cmd)}  # stage(s): {joined_stage_names}")
            continue

        result = subprocess.run(scancel_cmd, text=True, capture_output=True, check=False)
        stderr_text = (result.stderr or "").strip()
        stdout_text = (result.stdout or "").strip()
        if result.returncode == 0:
            if use_graceful_signal:
                print(
                    f"[SUCCESS] Requested graceful shutdown for Slurm job {job_id} for stage(s): {joined_stage_names}. "
                    "Solver jobs trap SIGUSR1 and write the latest safe off-cadence step at the next checkpoint."
                )
            else:
                print(f"[SUCCESS] Canceled Slurm job {job_id} for stage(s): {joined_stage_names}")
            continue

        detail = stderr_text or stdout_text or "unknown scancel failure"
        failures.append((job_id, joined_stage_names, detail, result.returncode))
        print(
            f"[ERROR] Failed to cancel Slurm job {job_id} for stage(s) {joined_stage_names}: {detail}",
            file=sys.stderr,
        )

    if args.dry_run:
        print("[INFO] Dry-run only. No jobs were canceled.")
        return

    if failures:
        sys.exit(1)


def init_case(args):
    """!
    @brief Implements the 'init' command.
    @details Creates a new case study directory by copying a template.
             Runtime binaries are resolved from the project bin/ directory
             via PATH; use 'sync-binaries' to pin specific versions locally.
    @param[in] args The command-line arguments parsed by argparse.
    """
    context = resolve_case_origin_context(source_root_override=getattr(args, "source_root", None))
    try:
        source_project_root = require_project_root(context["source_project_root"], "init")
        template_path = resolve_template_directory(source_project_root, args.template_name)
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    # The destination path is relative to the current working directory.
    dest_path = os.path.abspath(os.path.join(os.getcwd(), args.dest_name if args.dest_name else args.template_name))

    if os.path.exists(dest_path):
        print(f"[FATAL] Destination directory '{dest_path}' already exists.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Initializing new case '{os.path.basename(dest_path)}' from template '{args.template_name}'...")
    
    shutil.copytree(template_path, dest_path)
    print(f"[SUCCESS] Copied template files to: {dest_path}")

    copied_runtime_example = os.path.join(dest_path, RUNTIME_EXECUTION_EXAMPLE_FILENAME)
    if os.path.isfile(copied_runtime_example):
        os.remove(copied_runtime_example)

    try:
        runtime_result = ensure_case_runtime_execution_config(dest_path, source_project_root, overwrite=True)
        print(f"[INFO] Wrote optional runtime launcher config: {os.path.relpath(runtime_result['path'])}")
        if runtime_result["seed_source"] and os.path.basename(runtime_result["seed_source"]) == RUNTIME_EXECUTION_CONFIG_FILENAME:
            print("       Seeded from repo-local '.picurv-execution.yml'.")
        print("       Leave it unchanged for ordinary local runs; edit it only if your site needs custom MPI launcher tokens.")
    except Exception as e:
        print(f"[ERROR] Failed to write runtime execution config: {e}", file=sys.stderr)

    try:
        metadata_path, _ = write_case_origin_metadata(
            dest_path,
            source_project_root,
            template_name=args.template_name,
            template_managed_files=list_template_relative_files(
                template_path,
                excluded_rel_paths={RUNTIME_EXECUTION_EXAMPLE_FILENAME},
            ),
        )
        print(f"[INFO] Wrote case origin metadata: {os.path.relpath(metadata_path)}")
    except Exception as e:
        print(f"[ERROR] Failed to write case origin metadata: {e}", file=sys.stderr)

    cluster_profile_candidates = sorted(
        {
            os.path.basename(path)
            for pattern in ("*cluster*.yml", "*cluster*.yaml")
            for path in glob.glob(os.path.join(dest_path, pattern))
        }
    )
    if cluster_profile_candidates:
        print("[INFO] Cluster profile sample(s) copied with this case:")
        for profile_name in cluster_profile_candidates:
            print(f"  - {profile_name}")
        print("       Edit account/partition/module_setup and any batch-specific launcher overrides before using --cluster.")

    if getattr(args, "pin_binaries", False):
        print("[INFO] Pinning runtime binaries into case directory...")
        try:
            copied_binaries = sync_case_binaries(dest_path, source_project_root)
            for dest_file_path in copied_binaries:
                print(f"  - Pinned '{os.path.basename(dest_file_path)}'")
            print("[SUCCESS] Case directory is ready with pinned binaries.")
            print("          These local copies will be used instead of bin/ originals.")
        except ValueError as exc:
            print(f"[WARNING] {exc}", file=sys.stderr)
            print("          No binaries were pinned. Run 'picurv build' first.", file=sys.stderr)
    else:
        print("[SUCCESS] Case directory is ready.")
        print("          Runtime binaries (simulator, postprocessor) are resolved from bin/ automatically.")
        print("          To pin specific binary versions, re-run with --pin-binaries or use: picurv sync-binaries")
    print("          Ensure 'picurv' is on your PATH (source etc/picurv.sh) to run from any directory.")


def sync_case_binaries_command(args):
    """!
    @brief Refresh case-local executables from the source repository bin directory.
    @param[in] args Command-line style argument list supplied to the function.
    """
    try:
        context = resolve_case_origin_context(
            case_dir_hint=getattr(args, "case_dir", None),
            source_root_override=getattr(args, "source_root", None),
        )
        source_project_root = require_project_root(context["source_project_root"], "sync-binaries")
        case_dir = require_existing_case_dir(context["case_dir"], "sync-binaries", source_project_root)
        copied = sync_case_binaries(case_dir, source_project_root)
        metadata_path, metadata = write_case_origin_metadata(
            case_dir,
            source_project_root,
            template_name=context.get("template_name"),
            existing=context.get("metadata"),
        )
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"[SUCCESS] Refreshed {len(copied)} binaries in: {case_dir}")
    for dest_path in copied:
        print(f"  - {os.path.basename(dest_path)}")
    print(f"[INFO] Case origin metadata refreshed: {os.path.relpath(metadata_path)}")
    if metadata.get("last_known_source_git_commit"):
        print(f"[INFO] Source commit recorded: {metadata['last_known_source_git_commit']}")


def sync_case_config_command(args):
    """!
    @brief Refresh template-managed config/docs files in a case directory.
    @param[in] args Command-line style argument list supplied to the function.
    """
    try:
        context = resolve_case_origin_context(
            case_dir_hint=getattr(args, "case_dir", None),
            source_root_override=getattr(args, "source_root", None),
            template_name_override=getattr(args, "template_name", None),
        )
        source_project_root = require_project_root(context["source_project_root"], "sync-config")
        case_dir = require_existing_case_dir(context["case_dir"], "sync-config", source_project_root)
        template_name = context.get("template_name")
        template_dir = resolve_template_directory(source_project_root, template_name)
        existing_managed = context.get("metadata", {}).get("template_managed_files")
        if not isinstance(existing_managed, list):
            existing_managed = None
        summary = sync_case_template_files(
            case_dir,
            template_dir,
            overwrite=getattr(args, "overwrite", False),
            prune=getattr(args, "prune", False),
            managed_rel_paths=existing_managed,
        )
        metadata_path, _ = write_case_origin_metadata(
            case_dir,
            source_project_root,
            template_name=template_name,
            existing=context.get("metadata"),
            template_managed_files=summary["template_managed_files"],
        )
        runtime_result = ensure_case_runtime_execution_config(case_dir, source_project_root, overwrite=False)
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"[SUCCESS] Synced template files from '{template_name}' into: {case_dir}")
    print(f"[INFO] Copied new files     : {len(summary['copied'])}")
    print(f"[INFO] Overwritten files   : {len(summary['overwritten'])}")
    print(f"[INFO] Skipped modified    : {len(summary['skipped_modified'])}")
    print(f"[INFO] Already unchanged   : {len(summary['unchanged'])}")
    print(f"[INFO] Pruned stale files  : {len(summary['pruned'])}")
    if runtime_result["created"]:
        print(f"[INFO] Created runtime launcher config: {os.path.relpath(runtime_result['path'])}")
        if runtime_result["seed_source"] and os.path.basename(runtime_result["seed_source"]) == RUNTIME_EXECUTION_CONFIG_FILENAME:
            print("[INFO] Seed source              : repo-local .picurv-execution.yml")
    if summary.get("prune_requested_without_tracking"):
        print("[WARNING] Prune tracking unavailable for this case; no removed template files were deleted.", file=sys.stderr)
    print(f"[INFO] Case origin metadata refreshed: {os.path.relpath(metadata_path)}")


def pull_source_repo(args):
    """!
    @brief Refresh source branches in the repository resolved from a case directory.
    @param[in] args Command-line style argument list supplied to the function.
    """
    try:
        context = resolve_case_origin_context(
            case_dir_hint=getattr(args, "case_dir", None),
            source_root_override=getattr(args, "source_root", None),
        )
        source_project_root = require_project_root(context["source_project_root"], "pull-source")
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    rebase = not getattr(args, "no_rebase", False)
    remote = getattr(args, "remote", None)
    branch = getattr(args, "branch", None)
    current_branch_only = (
        getattr(args, "current_branch_only", False)
        or remote is not None
        or branch is not None
    )

    if not current_branch_only:
        pull_all_source_branches(source_project_root, "pull-source.log", rebase=rebase)
        return

    command = ["git", "pull"]
    if rebase:
        command.append("--rebase")
    if remote:
        command.append(remote)
        if branch:
            command.append(branch)
    elif branch:
        command.extend(["origin", branch])

    execute_command(command, source_project_root, "pull-source.log", {})

def build_project(args):
    """!
    @brief Implements the 'build' command.
    @details Executes the top-level Makefile directly, passing through any
             additional arguments to `make`. This allows for building,
             cleaning, and other Makefile targets via the orchestrator
             without maintaining a separate build wrapper script.
    @param[in] args The command-line arguments parsed by argparse.
    """

    print("\n" + "="*27 + " BUILD STAGE " + "="*27)
    try:
        context = resolve_case_origin_context(
            case_dir_hint=getattr(args, "case_dir", None),
            source_root_override=getattr(args, "source_root", None),
        )
        source_project_root = require_project_root(context["source_project_root"], "build")
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        sys.exit(1)

    makefile_path = os.path.join(source_project_root, "Makefile")

    if not os.path.isfile(makefile_path):
        print(f"[FATAL] Makefile not found at expected location: {makefile_path}", file=sys.stderr)
        print("        Please ensure the project root contains a valid Makefile.", file=sys.stderr)
        sys.exit(1)

    make_args = list(args.make_args or [])
    if make_args_include_explicit_goal(make_args):
        command = ["make"] + make_args
    else:
        command = ["make", "all"] + make_args
        print("[INFO] No explicit make target supplied; defaulting to 'all'.")
        print("       Use 'picurv build clean-project ...' or another target when you want a non-build make action.")
    # For the build process, we don't have a monitor.yml, so we pass an empty
    # dict to execute_command. The command should be run in the project root.
    execute_command(command, source_project_root, "build.log", {})
   



# ==============================================================================
