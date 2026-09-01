#!/usr/bin/env python3
"""Extract a normalized run-artifact topology snapshot from the CLI's own planner."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "tests" / "tooling" / "artifact_topology.json"
SNAPSHOT_PATH = REPO_ROOT / "docs" / "generated" / "artifact_topology_snapshot.json"
PICURV = REPO_ROOT / "picurv_cli" / "picurv"


def normalize(paths: list, workspace: Path, run_id: str) -> list:
    """!
    @brief Replace unstable path components with logical tokens.

    The workspace root, the generated run id, and any embedded timestamp differ on
    every invocation. Fingerprinting the raw plan would therefore report drift on
    every run; fingerprinting the normalized shape reports drift only when the
    topology actually changes.
    @param[in] paths Planned artifact paths.
    @param[in] workspace Temporary workspace root used for the plan.
    @param[in] run_id Generated run identifier.
    @return Sorted, normalized artifact tokens.
    """
    normalized = set()
    for raw in paths:
        text = str(raw)
        text = text.replace(str(workspace), "<workspace>")
        text = text.replace(run_id, "<run_id>")
        text = re.sub(r"\d{8}-\d{6}", "<timestamp>", text)
        normalized.add(text)
    return sorted(normalized)


def run_plan(workspace: Path, case_dir: Path, monitor: str, extra: list,
             post: str = "standard_analysis.yml") -> dict:
    """!
    @brief Invoke the dry-run planner for one scenario and return its raw plan.
    @param[in] workspace Scenario workspace.
    @param[in] case_dir Directory holding the initialized case.
    @param[in] monitor Monitor file name inside the case directory.
    @param[in] extra Additional CLI arguments for the scenario.
    @param[in] post Post recipe file name inside the case directory.
    @return Parsed plan mapping.
    @throws RuntimeError when the planner fails.
    """
    result = subprocess.run(
        [
            sys.executable, str(PICURV), "run", "--dry-run", "--format", "json",
            "--case", str(case_dir / "flat_channel.yml"),
            "--solver", str(case_dir / "Imp-MG-Standard.yml"),
            "--monitor", str(case_dir / monitor),
            "--post", str(case_dir / post),
            *extra,
        ],
        cwd=workspace, capture_output=True, text=True, check=False, timeout=180,
    )
    if result.returncode != 0:
        raise RuntimeError(f"dry-run failed:\n{result.stdout}\n{result.stderr}")
    return json.loads(result.stdout)


def init_case(workspace: Path) -> Path:
    """!
    @brief Initialize a shipped case inside a scenario workspace.
    @param[in] workspace Scenario workspace.
    @return Path to the created case directory.
    @throws RuntimeError when initialization fails.
    """
    case_dir = workspace / "case"
    result = subprocess.run(
        [sys.executable, str(PICURV), "init", "flat_channel", "--dest", str(case_dir)],
        cwd=workspace, capture_output=True, text=True, check=False, timeout=120,
    )
    if result.returncode != 0:
        raise RuntimeError(f"picurv init failed:\n{result.stdout}\n{result.stderr}")
    return case_dir


def custom_monitor(case_dir: Path) -> str:
    """!
    @brief Write a monitor variant with non-default output and log directory names.

    Custom directory names must appear in the plan. A planner that reported defaults
    regardless would hide the directory the runtime deletes on a fresh solve.
    @param[in] case_dir Case directory to write into.
    @return File name of the variant monitor.
    """
    source = (case_dir / "Standard_Output.yml").read_text(encoding="utf-8")
    source = re.sub(r'^(\s*)output:\s*"output"', r'\1output: "custom_out"', source, flags=re.M)
    source = re.sub(r'^(\s*)log:\s*"logs"', r'\1log: "custom_logs"', source, flags=re.M)
    (case_dir / "Custom_Dirs.yml").write_text(source, encoding="utf-8")
    return "Custom_Dirs.yml"


def flat_post_recipe(case_dir: Path) -> str:
    """!
    @brief Write a post recipe using the flat default output directory.

    @details Shipped recipes use both a flat `viz` and a nested
             `visualization/<recipe>` layout, and `post.yml -> io.output_directory` is
             resolved directly against the run directory. Exercising only the nested
             form let the snapshot look as though `visualization/<recipe>` were the
             fixed shape, which is how that identity came to be recorded as
             non-configurable.
    @param[in] case_dir Case directory to write into.
    @return File name of the variant recipe.
    """
    source = (case_dir / "standard_analysis.yml").read_text(encoding="utf-8")
    source = re.sub(r'^(\s*)output_directory:\s*.*$', r'\1output_directory: "viz"',
                    source, flags=re.M)
    (case_dir / "Flat_Viz_Analysis.yml").write_text(source, encoding="utf-8")
    return "Flat_Viz_Analysis.yml"


def map_to_identities(artifacts: list, contract: dict, configured: dict = None) -> dict:
    """!
    @brief Map each normalized artifact token onto a declared logical identity.

    An artifact the contract does not name is reported rather than ignored: an
    unmapped path is a documented-layout gap, which is exactly what a topology
    contract exists to surface.
    @param[in] artifacts Normalized artifact tokens.
    @param[in] contract Parsed topology contract.
    @param[in] configured Actual directory names this scenario used, so the two
               configurable identities resolve to distinct paths rather than sharing an
               ambiguous single-segment wildcard.
    @return Mapping of identity id to matched tokens, plus an `unmapped` list.
    """
    configured = configured or {}
    concrete = {
        "run.runtime_logs": configured.get("log", "logs"),
        "run.solver_output": configured.get("output", "output"),
        # The post output directory may be several segments deep, so the scenario's
        # actual configured value is substituted rather than a single-segment wildcard.
        "run.visualization": configured.get("post_output", "viz"),
    }
    rules = []
    for record in contract["artifacts"]:
        rule = record["path_rule"]
        pattern = re.escape(rule)
        pattern = pattern.replace(re.escape("<workspace>"), r"<workspace>")
        pattern = pattern.replace(re.escape("<run.root>"), r"<workspace>/runs/<run_id>")
        pattern = pattern.replace(re.escape("<run.config>"), r"<workspace>/runs/<run_id>/config")
        pattern = pattern.replace(re.escape("<run.solver_output>"), r"<workspace>/runs/<run_id>/[^/]+")
        pattern = pattern.replace(re.escape("<run.runtime_logs>"), r"<workspace>/runs/<run_id>/[^/]+")
        pattern = pattern.replace(re.escape("<run.scheduler>"), r"<workspace>/runs/<run_id>/scheduler")
        pattern = pattern.replace(
            re.escape('<configured log dir, default "logs">'),
            re.escape(concrete["run.runtime_logs"]),
        )
        pattern = pattern.replace(
            re.escape('<configured output dir, default "output">'),
            re.escape(concrete["run.solver_output"]),
        )
        pattern = pattern.replace(
            re.escape('<configured post output dir, default "viz">'),
            re.escape(concrete["run.visualization"]),
        )
        pattern = pattern.replace(
            re.escape("<run.runtime_logs>"),
            r"<workspace>/runs/<run_id>/" + re.escape(concrete["run.runtime_logs"]),
        )
        pattern = pattern.replace(
            re.escape("<run.solver_output>"),
            r"<workspace>/runs/<run_id>/" + re.escape(concrete["run.solver_output"]),
        )
        pattern = pattern.replace(re.escape("<role>"), r"[^/]+")
        pattern = pattern.replace(re.escape("<ext>"), r"[^/]+")
        pattern = pattern.replace(re.escape("<name>"), r"[^/]+")
        pattern = pattern.replace(re.escape("<recipe>"), r"[^/]+")
        pattern = pattern.replace(re.escape("<run_id>"), r"<run_id>")
        pattern = pattern.replace(re.escape("<n>"), r"[^/]+")
        rules.append((record["id"], re.compile("^" + pattern + "$")))

    # Specific rules must win over wildcard ones. A configurable directory such as
    # `<run.root>/<configured log dir>` matches any single segment, so first-match
    # order let it swallow manifest.json, output, and scheduler.
    def specificity(item) -> tuple:
        """!
        @brief Rank a rule so literal path rules are tried before wildcard ones.
        @param[in] item Tuple of identity id and compiled pattern.
        @return Sort key placing more specific rules first.
        """
        pattern = item[1].pattern
        return (pattern.count("[^/]+"), -len(pattern))

    ordered = sorted(rules, key=specificity)
    mapped: dict = {rid: [] for rid, _ in rules}
    unmapped = []
    for token in artifacts:
        for rid, pattern in ordered:
            if pattern.match(token):
                mapped[rid].append(token)
                break
        else:
            unmapped.append(token)
    return {"mapped": {k: sorted(v) for k, v in mapped.items() if v}, "unmapped": sorted(unmapped)}


def build_snapshot() -> dict:
    """!
    @brief Build the normalized topology snapshot across scenarios.
    @return Snapshot mapping.
    """
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    scenarios = []
    workspace = Path(tempfile.mkdtemp(prefix="picurv-topology-"))
    try:
        case_dir = init_case(workspace)
        variant = custom_monitor(case_dir)
        flat_post = flat_post_recipe(case_dir)
        cases = [
            ("fresh_local_solve_and_post", "Standard_Output.yml",
             ["--solve", "--post-process"], "standard_analysis.yml"),
            ("fresh_local_solve_only", "Standard_Output.yml", ["--solve"],
             "standard_analysis.yml"),
            ("custom_directories", variant, ["--solve", "--post-process"],
             "standard_analysis.yml"),
            # A second post layout, so the snapshot records that the post output
            # directory is configured rather than fixed.
            ("flat_post_output_directory", "Standard_Output.yml",
             ["--solve", "--post-process"], flat_post),
        ]
        for name, monitor, extra, post in cases:
            dirs = {"log": "custom_logs", "output": "custom_out"} if monitor == variant else {}
            dirs["post_output"] = ("viz" if post == flat_post
                                   else "visualization/standard_analysis")
            plan = run_plan(workspace, case_dir, monitor, extra, post)
            artifacts = normalize(plan.get("artifacts", []), workspace, plan["run_id_preview"])
            scenarios.append({
                "scenario": name,
                "launch_mode": plan.get("launch_mode"),
                "artifacts": artifacts,
                "identities": map_to_identities(artifacts, contract, dirs),
            })
    finally:
        shutil.rmtree(workspace, ignore_errors=True)
    return {
        "default_layout": contract["default_layout"],
        "isolation_enforced": False,
        "logical_artifacts": [a["id"] for a in contract["artifacts"]],
        "scenarios": scenarios,
    }


def main() -> int:
    """!
    @brief Write or verify the artifact topology snapshot.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(description="Extract the run artifact topology snapshot.")
    parser.add_argument("--check", action="store_true", help="Fail if the snapshot is stale.")
    args = parser.parse_args()

    try:
        snapshot = build_snapshot()
    except (RuntimeError, subprocess.TimeoutExpired, json.JSONDecodeError) as error:
        print(f"Artifact topology extraction failed: {error}", file=sys.stderr)
        return 1

    content = json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not SNAPSHOT_PATH.is_file() or SNAPSHOT_PATH.read_text(encoding="utf-8") != content:
            print(
                "Run artifact topology has changed. The planned artifact set no longer matches "
                "the recorded snapshot.\n"
                "  Review the pages that document run layout, then refresh with:\n"
                "    make docs-topology",
                file=sys.stderr,
            )
            return 1
        print(f"Artifact topology snapshot is current ({len(snapshot['scenarios'])} scenario(s)).")
        return 0

    SNAPSHOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    SNAPSHOT_PATH.write_text(content, encoding="utf-8")
    total = sum(len(s["artifacts"]) for s in snapshot["scenarios"])
    print(f"Wrote artifact topology snapshot: {len(snapshot['scenarios'])} scenario(s), {total} artifacts.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
