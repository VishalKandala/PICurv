"""Regressions for git-backed file enumeration.

`git ls-files` reports a path the index still tracks after the worktree file is gone.
Every documentation scanner enumerates that set, so a deleted page used to crash them
with FileNotFoundError. Deleting a page is ordinary work, and the constituent audits
must keep working against a dirty tree.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLING = REPO_ROOT / "tests" / "tooling"


def _load(name):
    """!
    @brief Import a tooling script by path.
    @param[in] name Script stem.
    @return The imported module.
    """
    spec = importlib.util.spec_from_file_location(name, TOOLING / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


repo_files = _load("repo_files")


@pytest.fixture()
def repo_with_deleted_tracked_file(tmp_path):
    """!
    @brief A git repository whose index tracks a Markdown file the worktree lacks.
    @param[in] tmp_path Temporary directory fixture.
    @return Tuple of (repository root, deleted path, surviving path).
    """
    root = tmp_path / "repo"
    root.mkdir()
    run = lambda *args: subprocess.run(["git", "-C", str(root), *args], check=True,
                                       capture_output=True)
    run("init", "-q")
    run("config", "user.email", "t@example.com")
    run("config", "user.name", "t")
    doomed = root / "doomed.md"
    survivor = root / "survivor.md"
    doomed.write_text("# doomed\n", encoding="utf-8")
    survivor.write_text("# survivor\n", encoding="utf-8")
    run("add", "-A")
    run("commit", "-qm", "seed")
    doomed.unlink()
    return root, doomed, survivor


def test_deleted_tracked_file_is_not_enumerated(repo_with_deleted_tracked_file):
    """!
    @brief A tracked path missing from the worktree must not be returned.
    @param[in] repo_with_deleted_tracked_file Repository fixture.
    @return None.
    """
    root, doomed, survivor = repo_with_deleted_tracked_file
    found = repo_files.enumerate_repository_files(root, ".md")
    assert survivor in found
    assert doomed not in found


def test_enumeration_survives_a_staged_deletion(repo_with_deleted_tracked_file):
    """!
    @brief Staging the deletion must not change the answer.
    @param[in] repo_with_deleted_tracked_file Repository fixture.
    @return None.
    """
    root, doomed, survivor = repo_with_deleted_tracked_file
    subprocess.run(["git", "-C", str(root), "add", "-A"], check=True, capture_output=True)
    found = repo_files.enumerate_repository_files(root, ".md")
    assert found == [survivor]


def test_every_enumerated_path_can_be_read(repo_with_deleted_tracked_file):
    """!
    @brief The contract callers rely on: what is returned can be opened.
    @param[in] repo_with_deleted_tracked_file Repository fixture.
    @return None.
    """
    root, _, _ = repo_with_deleted_tracked_file
    for path in repo_files.enumerate_repository_files(root, ".md"):
        path.read_text(encoding="utf-8")


def test_untracked_files_are_included(repo_with_deleted_tracked_file):
    """!
    @brief Non-ignored untracked Markdown is part of the reviewable set.
    @param[in] repo_with_deleted_tracked_file Repository fixture.
    @return None.
    """
    root, _, _ = repo_with_deleted_tracked_file
    fresh = root / "fresh.md"
    fresh.write_text("# fresh\n", encoding="utf-8")
    assert fresh in repo_files.enumerate_repository_files(root, ".md")


def test_skip_dirs_are_honoured(repo_with_deleted_tracked_file):
    """!
    @brief A skipped component excludes files beneath it.
    @param[in] repo_with_deleted_tracked_file Repository fixture.
    @return None.
    """
    root, _, _ = repo_with_deleted_tracked_file
    nested = root / "scratch"
    nested.mkdir()
    (nested / "note.md").write_text("# note\n", encoding="utf-8")
    found = repo_files.enumerate_repository_files(root, ".md", frozenset({"scratch"}))
    assert all("scratch" not in path.parts for path in found)


def test_non_git_directory_reports_unavailable(tmp_path):
    """!
    @brief Outside a repository the helper reports that git enumeration failed.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    assert repo_files.enumerate_repository_files(tmp_path / "nowhere", ".md") is None


@pytest.mark.parametrize("script", [
    "audit_docs_site.py", "audit_path_literals.py", "check_markdown_links.py",
])
def test_scanners_run_against_the_current_dirty_tree(script):
    """!
    @brief Constituent audits must pass while the tree carries a deleted tracked page.
    @param[in] script Scanner to run.
    @return None.
    """
    result = subprocess.run([sys.executable, str(TOOLING / script)], cwd=REPO_ROOT,
                            capture_output=True, text=True)
    assert result.returncode == 0, f"{script} failed:\n{result.stdout}\n{result.stderr}"
