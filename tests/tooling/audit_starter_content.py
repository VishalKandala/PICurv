#!/usr/bin/env python3
"""Certify that shipped starter templates and reusable configuration profiles stay usable."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "tests" / "tooling" / "starter_content_contract.json"
PICURV = REPO_ROOT / "picurv_cli" / "picurv"
RUNTIME_EXECUTION_EXAMPLE = "execution.example.yml"


def run_cli(args: list[str], cwd: Path = REPO_ROOT) -> subprocess.CompletedProcess[str]:
    """!
    @brief Run one CLI command used to validate a shipped starter artifact.
    @param[in] args CLI argument list.
    @param[in] cwd Working directory for the command.
    @return Completed CLI process.
    """
    return subprocess.run(
        [sys.executable, str(PICURV), *args], cwd=cwd, text=True, capture_output=True, timeout=90, check=False
    )


def fail(context: str, result: subprocess.CompletedProcess[str]) -> None:
    """!
    @brief Raise a compact error that preserves the relevant CLI output.
    @param[in] context Human-readable operation description.
    @param[in] result Completed failing CLI process.
    @return None.
    """
    raise RuntimeError(f"{context} failed:\n{result.stdout}\n{result.stderr}")


def validate_bundle(bundle: dict[str, str], label: str) -> None:
    """!
    @brief Validate one declared case/solver/monitor/post or cluster/study composition.
    @param[in] bundle Declared role-to-path mapping.
    @param[in] label Human-readable bundle description.
    @return None.
    """
    args = ["validate"]
    for role, path in bundle.items():
        args.extend([f"--{role}", str(REPO_ROOT / path)])
    result = run_cli(args)
    if result.returncode:
        fail(label, result)


def _yaml_role(path: Path) -> Optional[str]:
    """!
    @brief Classify a starter YAML by the same stable top-level shapes used by init.
    @param[in] path Source YAML path.
    @return Canonical configuration role, or None for unclassified YAML.
    """
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(payload, dict):
        return None
    keys = set(payload)
    if {"grid", "properties", "run_control"} <= keys:
        return "case"
    if "base_configs" in keys and ("study_type" in keys or "parameters" in keys or "parameter_sets" in keys):
        return "study"
    if "scheduler" in keys and "resources" in keys:
        return "cluster"
    if "source_data" in keys or "eulerian_pipeline" in keys or "lagrangian_pipeline" in keys:
        return "post"
    if "io" in keys and ("logging" in keys or "profiling" in keys or "diagnostics" in keys):
        return "monitor"
    if "momentum_solver" in keys or "poisson_solver" in keys or "operation_mode" in keys:
        return "solver"
    return None


def audit_template_copy(template_name: str, temporary_root: Path) -> None:
    """!
    @brief Initialize one declared template and verify canonical relocation preserves its content.
    @param[in] template_name Top-level example directory name.
    @param[in] temporary_root Temporary parent directory for initialized cases.
    @return None.
    """
    source = REPO_ROOT / "examples" / template_name
    destination = temporary_root / template_name
    result = run_cli(["init", template_name, "--dest", str(destination)])
    if result.returncode:
        fail(f"picurv init {template_name}", result)
    required_dirs = (
        "config", "config/studies", "inputs", "assets/objects", "assets/sets", "runs", "studies"
    )
    for relative in required_dirs:
        if not (destination / relative).is_dir():
            raise RuntimeError(f"picurv init {template_name} omitted workspace directory {relative}")
    if not (destination / ".picurv-workspace.yml").is_file():
        raise RuntimeError(f"picurv init {template_name} omitted .picurv-workspace.yml")

    source_yamls = [path for path in source.rglob("*.yml") if path.name != RUNTIME_EXECUTION_EXAMPLE]
    initialized_yamls = list((destination / "config").rglob("*.yml"))
    if len(initialized_yamls) < len(source_yamls):
        raise RuntimeError(f"picurv init {template_name} lost one or more YAML configurations")

    for source_path in source.rglob("*"):
        if not source_path.is_file():
            continue
        relative = source_path.relative_to(source)
        if relative.as_posix() == RUNTIME_EXECUTION_EXAMPLE:
            continue
        if source_path.suffix.lower() in {".yml", ".yaml"}:
            role = _yaml_role(source_path)
            stem = source_path.stem
            candidates = [
                destination / "config" / source_path.name,
                destination / "config" / f"{role}.yml" if role else destination / source_path.name,
                destination / "config" / f"{role}-{stem}{source_path.suffix}" if role else destination / source_path.name,
                destination / "config" / "studies" / source_path.name,
            ]
            if not any(candidate.is_file() for candidate in candidates):
                raise RuntimeError(f"picurv init {template_name} did not relocate YAML {relative}")
            continue
        matches = [path for path in destination.rglob(source_path.name) if path.is_file()]
        if not any(path.read_bytes() == source_path.read_bytes() for path in matches):
            raise RuntimeError(f"picurv init {template_name} did not preserve {relative}")
    if (destination / RUNTIME_EXECUTION_EXAMPLE).exists():
        raise RuntimeError(f"picurv init {template_name} copied site-specific {RUNTIME_EXECUTION_EXAMPLE}")


def main() -> int:
    """!
    @brief Run the starter-content inventory, composition, and initializer audit.
    @return Process status code.
    """
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    actual_templates = sorted(path.name for path in (REPO_ROOT / "examples").iterdir() if path.is_dir())
    if actual_templates != sorted(contract["template_directories"]):
        raise RuntimeError(f"Template inventory differs from contract: actual={actual_templates}")
    unknown_reference_templates = set(contract["reference_only_templates"]) - set(actual_templates)
    if unknown_reference_templates:
        raise RuntimeError(f"Reference-only templates are not in the template inventory: {sorted(unknown_reference_templates)}")

    actual_config_assets = sorted(str(path.relative_to(REPO_ROOT)) for path in (REPO_ROOT / "config").rglob("*") if path.is_file())
    if actual_config_assets != sorted(contract["config_assets"]):
        raise RuntimeError("Configuration asset inventory differs from the starter-content contract.")

    declared_example_yamls = set(contract["auxiliary_example_yamls"])
    for bundle in contract["case_bundles"] + contract["study_bundles"]:
        declared_example_yamls.update(path for path in bundle.values() if path.startswith("examples/"))
    actual_example_yamls = {str(path.relative_to(REPO_ROOT)) for path in (REPO_ROOT / "examples").rglob("*.yml")}
    if actual_example_yamls != declared_example_yamls:
        raise RuntimeError("Example YAML inventory differs from the declared runnable/reference compositions.")

    for bundle in contract["case_bundles"]:
        validate_bundle(bundle, f"case bundle {bundle['case']}")
    for bundle in contract["config_role_bundles"]:
        validate_bundle(bundle, f"config role bundle {bundle['solver']}")
    for bundle in contract["study_bundles"]:
        validate_bundle(bundle, f"study bundle {bundle['study']}")

    with tempfile.TemporaryDirectory(prefix="picurv-starter-content-") as temporary_directory:
        temporary_root = Path(temporary_directory)
        for template_name in contract["template_directories"]:
            audit_template_copy(template_name, temporary_root)

    print("Starter template, example, and configuration audit passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"Starter-content audit failed: {error}", file=sys.stderr)
        raise SystemExit(1)
