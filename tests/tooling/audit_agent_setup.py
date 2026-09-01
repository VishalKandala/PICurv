#!/usr/bin/env python3
"""Audit and synchronize the repository-portable agent instruction setup."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_ROOT = REPO_ROOT / ".agents" / "skills"
CLAUDE_ROOT = REPO_ROOT / ".claude" / "skills"
EXPECTED_SKILLS = {
    "picurv-capability-change",
    "picurv-solver-debugging",
    "picurv-spatial-kernel-change",
}
LOCAL_SETTINGS = ".claude/settings.local.json"
LOCAL_SETTINGS_PATTERN = "/.claude/settings.local.json"


def frontmatter(path: Path) -> dict[str, str]:
    """!
    @brief Parse the simple YAML front matter used by a skill file.
    @param[in] path Skill file to inspect.
    @return Parsed scalar front-matter fields.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "---":
        return {}
    fields: dict[str, str] = {}
    for line in lines[1:]:
        if line == "---":
            return fields
        if ":" in line:
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
    return {}


def skill_entries(root: Path) -> dict[str, Path]:
    """!
    @brief Return skill directories directly below one discovery root.
    @param[in] root Skill discovery directory.
    @return Mapping from directory name to path, including symlinked directories.
    """

    if not root.is_dir():
        return {}
    return {
        child.name: child
        for child in root.iterdir()
        if child.is_dir() or child.is_symlink()
    }


def entry_names(root: Path) -> set[str]:
    """!
    @brief Return every direct child name under a skill discovery root.
    @param[in] root Skill discovery directory.
    @return Child-name set, including unexpected files and links.
    """

    return {child.name for child in root.iterdir()} if root.is_dir() else set()


def remove_entry(path: Path) -> None:
    """!
    @brief Remove one skill-tree entry without following a directory symlink.
    @param[in] path Exact materialized or symlinked skill path.
    @return None.
    """

    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def synchronize() -> None:
    """!
    @brief Materialize byte-identical Claude skill copies from the canonical tree.
    @return None.
    @throws RuntimeError when the canonical skill set is incomplete or unexpected.
    """

    canonical = skill_entries(CANONICAL_ROOT)
    canonical_names = entry_names(CANONICAL_ROOT)
    if canonical_names != EXPECTED_SKILLS:
        missing = sorted(EXPECTED_SKILLS - canonical_names)
        extra = sorted(canonical_names - EXPECTED_SKILLS)
        raise RuntimeError(f"canonical skill set differs: missing={missing}, extra={extra}")

    CLAUDE_ROOT.mkdir(parents=True, exist_ok=True)
    for path in list(CLAUDE_ROOT.iterdir()):
        remove_entry(path)
    for name in sorted(EXPECTED_SKILLS):
        shutil.copytree(canonical[name], CLAUDE_ROOT / name, symlinks=False)


def git(*args: str) -> subprocess.CompletedProcess[str]:
    """!
    @brief Run Git with machine-global excludes disabled.
    @param[in] args Arguments following `git`.
    @return Completed process with captured text output.
    """

    return subprocess.run(
        ["git", "-c", "core.excludesFile=/dev/null", "-C", str(REPO_ROOT), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def audit() -> list[str]:
    """!
    @brief Validate shared instructions, skill discovery copies, and local-settings hygiene.
    @return Human-readable violations; an empty list means the audit passed.
    """

    errors: list[str] = []
    agents = REPO_ROOT / "AGENTS.md"
    claude = REPO_ROOT / "CLAUDE.md"

    if not agents.is_file() or agents.is_symlink():
        errors.append("AGENTS.md must be a tracked regular file, not a symlink")
    if not claude.is_file() or claude.is_symlink():
        errors.append("CLAUDE.md must be a regular file, not a symlink")
    elif "@AGENTS.md" not in {line.strip() for line in claude.read_text(encoding="utf-8").splitlines()}:
        errors.append("CLAUDE.md must import the canonical @AGENTS.md instructions")

    canonical = skill_entries(CANONICAL_ROOT)
    copies = skill_entries(CLAUDE_ROOT)
    canonical_names = entry_names(CANONICAL_ROOT)
    copy_names = entry_names(CLAUDE_ROOT)
    if canonical_names != EXPECTED_SKILLS:
        errors.append(
            "canonical skill set mismatch: "
            f"missing={sorted(EXPECTED_SKILLS - canonical_names)}, "
            f"extra={sorted(canonical_names - EXPECTED_SKILLS)}"
        )
    if copy_names != EXPECTED_SKILLS:
        errors.append(
            "Claude skill set mismatch: "
            f"missing={sorted(EXPECTED_SKILLS - copy_names)}, "
            f"extra={sorted(copy_names - EXPECTED_SKILLS)}; run make sync-agent-skills"
        )

    for name in sorted(EXPECTED_SKILLS & set(canonical)):
        source_dir = canonical[name]
        source = source_dir / "SKILL.md"
        if source_dir.is_symlink() or not source.is_file() or source.is_symlink():
            errors.append(f"canonical skill {name} must contain a regular SKILL.md")
            continue
        metadata = frontmatter(source)
        if metadata.get("name") != name:
            errors.append(f"{source.relative_to(REPO_ROOT)} front-matter name must be {name!r}")
        if not metadata.get("description"):
            errors.append(f"{source.relative_to(REPO_ROOT)} requires a non-empty description")

        target_dir = copies.get(name)
        target = target_dir / "SKILL.md" if target_dir else None
        if target_dir is None:
            continue
        if target_dir.is_symlink() or target is None or not target.is_file() or target.is_symlink():
            errors.append(f"Claude skill {name} must be a materialized regular-file copy")
            continue
        source_files = sorted(
            path.relative_to(source_dir) for path in source_dir.rglob("*") if path.is_file()
        )
        target_files = sorted(
            path.relative_to(target_dir) for path in target_dir.rglob("*") if path.is_file()
        )
        if source_files != target_files:
            errors.append(f"Claude skill {name} file set differs; run make sync-agent-skills")
            continue
        for relative in source_files:
            if (source_dir / relative).read_bytes() != (target_dir / relative).read_bytes():
                errors.append(
                    f"Claude skill {name}/{relative} differs byte-for-byte; "
                    "run make sync-agent-skills"
                )

    tracked = git("ls-files", "--error-unmatch", "--", LOCAL_SETTINGS)
    if tracked.returncode == 0:
        errors.append(f"{LOCAL_SETTINGS} is machine-local and must not be tracked")
    ignored = git("check-ignore", "-v", "--", LOCAL_SETTINGS)
    if ignored.returncode != 0:
        errors.append(f"{LOCAL_SETTINGS} must be ignored by the repository .gitignore")
    else:
        fields = ignored.stdout.rstrip().split("\t", 1)[0].rsplit(":", 2)
        source = fields[0] if len(fields) == 3 else ""
        pattern = fields[2] if len(fields) == 3 else ""
        if Path(source).name != ".gitignore" or pattern != LOCAL_SETTINGS_PATTERN:
            errors.append(
                f"{LOCAL_SETTINGS} must use the exact repository ignore pattern "
                f"{LOCAL_SETTINGS_PATTERN!r}, not a user-global or incidental rule"
            )

    return errors


def parse_args(argv: list[str]) -> argparse.Namespace:
    """!
    @brief Parse audit command-line arguments.
    @param[in] argv Arguments excluding the executable name.
    @return Parsed command-line namespace.
    """

    parser = argparse.ArgumentParser(
        description="Audit or synchronize PICurv's portable shared agent setup."
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="replace Claude skill entries with materialized copies of canonical .agents skills",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """!
    @brief Synchronize when requested, then fail closed on portability drift.
    @param[in] argv Optional arguments excluding the executable name.
    @return Zero on success and one on any setup violation.
    """

    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        if args.sync:
            synchronize()
    except (OSError, RuntimeError) as error:
        print(f"Agent setup synchronization failed: {error}", file=sys.stderr)
        return 1

    errors = audit()
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    action = "synchronized and audited" if args.sync else "audited"
    print(f"Agent setup {action}: shared instructions and {len(EXPECTED_SKILLS)} skills are portable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
