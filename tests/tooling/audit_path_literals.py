#!/usr/bin/env python3
"""Reject unmanaged run-path literals so a layout change cannot leave prose stale."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from repo_files import enumerate_repository_files  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
TOPOLOGY = REPO_ROOT / "tests" / "tooling" / "artifact_topology.json"
SKIP_DIRS = {".git", "docs_build", "obj", "bin", "stubs", "runs", "studies", "__pycache__", ".pytest_cache"}

# Historical records describe layouts as they were; rewriting them would falsify history.
HISTORICAL = {"docs/CHANGELOG.md"}

# Directories a run owns. Their names are now fixed rather than configurable, but a
# bare `logs/...` in prose is still wrong for a different reason: three different
# owners ship a directory by that kind of name - the repository's own build logs and
# `config/` library, a workspace's editable `config/` and `inputs/`, and a run's
# generated tree. A bare prefix does not say which one is meant, and a reader
# following it lands in the wrong place. Logical identities name the owner as well as
# the path, and artifact_topology.json is where the mapping is enumerated, so a layout
# change is a contract change rather than a prose hunt.
RUN_OWNED = "config|logs|output|scheduler|visualization|checkpoints"

# Concrete run-relative directories that a layout change would invalidate.
MANAGED = re.compile(
    r"(?<![\w./<])runs/(?:<[^>]+>|\$\{[^}]+\}|[A-Za-z0-9_.*-]+)/(" + RUN_OWNED + r")\b")

# A run-owned directory named against an unresolved run root - `<run_dir>/logs/...`.
# The `<run_dir>` placeholder says "somewhere in the run" but still hardcodes which
# subdirectory, so it drifts exactly like the concrete form.
PLACEHOLDER_ROOT = re.compile(
    r"<(?:run_dir|run|RUN_DIR)>/(" + RUN_OWNED + r")\b")

# A bare run-owned prefix at the start of a path - `logs/Runtime_Memory.log`. This is
# the form the audit used to miss entirely, and the one pages reach for most often.
# `<repo>/logs/...` is the escape hatch for the repository's own build and test logs,
# which are a different directory that no run configuration moves.
# The trailing part is optional: `logs/` on its own is exactly as wrong as
# `logs/Runtime_Memory.log`, and it is the form pages reach for when naming a
# directory rather than a file. Requiring a character after the slash let every
# terminal directory reference through.
BARE_PREFIX = re.compile(
    r"(?<![\w./<>-])(" + RUN_OWNED + r")/(?![A-Za-z0-9_./*-]*[<>])[A-Za-z0-9_./*-]*")

# A repository-owned path, written with an explicit `<repo>/` prefix. That prefix is
# what distinguishes the repository's build log from a run's runtime logs, and its
# shipped `config/` library from a run's config directory - all of which otherwise
# read identically in prose while naming directories with different owners.
REPO_OWNED = re.compile(r"<repo>/[A-Za-z0-9_][A-Za-z0-9_./*<>-]*")

# A workspace-owned path, written with an explicit `<workspace>/` prefix. An
# initialized workspace owns `config/`, `inputs/`, and `assets/` directories whose
# names match a run's, so the prefix is what tells a reader which of the two a
# sentence means.
WORKSPACE_OWNED = re.compile(r"<workspace>/[A-Za-z0-9_][A-Za-z0-9_./*<>-]*")

# Bare prefixes that are not run-owned paths at all: repository directories that
# happen to share a name. `config/` is a real top-level source directory.
# The repository ships a top-level `config/` library of reusable profiles. Its
# subdirectories are source, not run output, and are named from their real paths.
REPO_CONFIG_SUBDIRS = ("guide.md", "build", "grids", "initial_conditions", "monitors",
                       "postprocessors", "profiles", "runtime", "schedulers", "solvers",
                       "studies")
NOT_RUN_OWNED_PREFIXES = re.compile(
    r"(?<![\w./<>-])config/(?:" + "|".join(
        name.replace(".", r"\.") for name in REPO_CONFIG_SUBDIRS) + r")\b")

# Occurrences that legitimately show a concrete path.
ALLOWED_CONTEXTS = (
    "```",          # runnable command examples
    "@verbinclude", # embedded executable templates
    "@code",
)


def tracked_markdown() -> list:
    """!
    @brief Markdown files the current commit carries.
    @return Sorted Markdown paths.
    """
    found = enumerate_repository_files(REPO_ROOT, ".md", frozenset(SKIP_DIRS))
    if found is None:
        return sorted(path for path in REPO_ROOT.rglob("*.md") if path.is_file())
    return found


def in_code_block(lines: list, index: int) -> bool:
    """!
    @brief Whether a line sits inside a fenced code block.
    @param[in] lines All lines of the file.
    @param[in] index Zero-based line index.
    @return True when the line is inside a fence.
    """
    fences = sum(1 for line in lines[:index] if line.lstrip().startswith("```"))
    return fences % 2 == 1


def main() -> int:
    """!
    @brief Fail when narrative prose hardcodes a run-relative path.
    @return Process status code.
    """
    contract = json.loads(TOPOLOGY.read_text(encoding="utf-8"))
    identities = [a["id"] for a in contract["artifacts"]]
    violations: list = []
    scanned = 0
    for path in tracked_markdown():
        if str(path.relative_to(REPO_ROOT)) in HISTORICAL:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        scanned += 1
        for number, line in enumerate(lines):
            if in_code_block(lines, number) or any(token in line for token in ALLOWED_CONTEXTS):
                continue
            probe = NOT_RUN_OWNED_PREFIXES.sub("", line)
            # `<repo>/logs/...` is the repository's own build and test log directory,
            # which no run configuration moves. Remove the whole reference before
            # looking for bare prefixes, so the distinct notation actually buys the
            # author something rather than leaving a bare `logs/` behind.
            probe = REPO_OWNED.sub("", probe)
            probe = WORKSPACE_OWNED.sub("", probe)
            if MANAGED.search(line):
                found, advice = "unmanaged run-path literal", (
                    "Use logical notation such as `<run.config>`, or move the concrete "
                    "path into a runnable command block.")
            elif PLACEHOLDER_ROOT.search(line):
                found, advice = "run-owned directory named under an unresolved run root", (
                    "`<run_dir>/logs/...` names a subdirectory without naming the "
                    "contract that fixes it. Use the logical identity - "
                    "`<run.runtime_logs>/...` - which artifact_topology.json maps.")
            elif BARE_PREFIX.search(probe):
                match = BARE_PREFIX.search(probe)
                found, advice = f"bare run-owned prefix `{match.group(0)}`", (
                    "A bare `logs/...` or `config/...` does not say which owner is "
                    "meant. Use the run's logical identity - `<run.runtime_logs>/...`, "
                    "`<run.solver_output>/...` - or name the other owner explicitly "
                    "with `<repo>/logs/...` or `<workspace>/config/...`.")
            else:
                continue
            violations.append(
                f"{path.relative_to(REPO_ROOT)}:{number + 1}: {found}\n"
                f"      {line.strip()[:96]}\n"
                f"      {advice}"
            )
    if violations:
        print("Unmanaged run-path literals in narrative prose:", file=sys.stderr)
        for violation in violations:
            print(f"  {violation}", file=sys.stderr)
        print(
            f"\nLogical identities are declared in tests/tooling/artifact_topology.json "
            f"({len(identities)} identities). Narrative pages should refer to those, so a layout "
            f"change is a contract change rather than a prose hunt.",
            file=sys.stderr,
        )
        return 1
    print(f"Path-literal audit passed: {scanned} Markdown files, no unmanaged run-path literals in prose.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
