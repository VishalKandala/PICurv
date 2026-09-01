#!/usr/bin/env python3
"""One git-backed enumeration of the repository's reviewable files.

Every documentation scanner needs the same set: files this commit carries, plus
non-ignored untracked ones, minus generated scratch. Four scanners each had their own
copy, and every copy shared one defect - `git ls-files` reports a path the index still
tracks even after the worktree file is deleted, so a scanner would enumerate it and
then fail on `read_text()`. Deleting a page is ordinary work, and the constituent
audits must keep running against a dirty tree even though `certify-docs` refuses one.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


def enumerate_repository_files(repo_root: Path, suffix: str = "",
                               skip_dirs: frozenset = frozenset()) -> list:
    """!
    @brief List tracked plus non-ignored untracked files that exist on disk.

    @details Paths the index tracks but the worktree no longer holds are filtered out:
             a deletion staged or unstaged is still a deletion, and reading it would
             raise rather than report a finding.
    @param[in] repo_root Repository root directory.
    @param[in] suffix Restrict to paths with this suffix, or empty for all files.
    @param[in] skip_dirs Path components whose presence excludes a file.
    @return Sorted existing paths, or None when git enumeration is unavailable.
    """
    entries: list = []
    for args in (["ls-files", "-z"], ["ls-files", "-z", "--others", "--exclude-standard"]):
        result = subprocess.run(
            ["git", "-C", str(repo_root), *args], capture_output=True, text=True, check=False
        )
        if result.returncode != 0:
            return None
        entries.extend(
            entry for entry in result.stdout.split("\0")
            if entry and (not suffix or entry.endswith(suffix))
        )
    return sorted({
        repo_root / entry for entry in entries
        if not skip_dirs.intersection(Path(entry).parts)
        and (repo_root / entry).is_file()
    })
