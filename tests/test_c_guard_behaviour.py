"""Behavioural regressions for the C-side log-directory guard.

Source-text assertions pass just as happily on dead code. These drive the built
simulator on a disposable run tree and assert what it actually refuses.

Every case uses a temporary directory. None of them targets the filesystem root: `/` is
covered by classification only, in tests/test_run_directory_containment.py.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SIMULATOR = REPO_ROOT / "bin" / "simulator"
EXAMPLE = REPO_ROOT / "examples" / "flat_channel"

pytestmark = pytest.mark.skipif(
    not SIMULATOR.is_file(),
    reason="bin/simulator is not built; run `make simulator` to exercise the C guard",
)


def _staged_run(tmp_path):
    """!
    @brief Stage a real run whose control file the guard will read.
    @param[in] tmp_path Temporary directory fixture.
    @return Tuple of (run directory, base control file).
    """
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    for name in ("quickstart_flat_channel.yml", "quickstart_Standard_Output.yml",
                 "quickstart_standard_analysis.yml", "Imp-MG-Standard.yml"):
        shutil.copy(EXAMPLE / name, workspace / name)
    result = subprocess.run(
        [str(REPO_ROOT / "bin" / "picurv"), "run", "--solve", "-n", "1",
         "--case", "quickstart_flat_channel.yml",
         "--solver", "Imp-MG-Standard.yml",
         "--monitor", "quickstart_Standard_Output.yml"],
        cwd=workspace, capture_output=True, text=True, timeout=900,
    )
    if result.returncode != 0:
        pytest.skip(f"could not stage a run in this environment:\n{result.stdout}\n{result.stderr}")
    runs = sorted((workspace / "runs").iterdir())
    assert runs, "staging produced no run directory"
    run_dir = runs[-1]
    control = next(iter(sorted((run_dir / "config").glob("*.control"))))
    return run_dir, control


def _probe(run_dir, control, log_dir, authorized):
    """!
    @brief Run the simulator with one log-directory value and report the refusal.
    @param[in] run_dir Working directory for the run.
    @param[in] control Base control file to derive from.
    @param[in] log_dir Value to place in `-log_dir`.
    @param[in] authorized Whether to supply `-allow_unsafe_log_dir true`.
    @return Combined stdout/stderr text.
    """
    lines = [
        line for line in control.read_text(encoding="utf-8").splitlines()
        if not line.startswith("-log_dir") and not line.startswith("-allow_unsafe_log_dir")
    ]
    if authorized:
        lines.append("-allow_unsafe_log_dir true")
    lines.append(f"-log_dir {log_dir}")
    probe = control.parent / "probe.control"
    probe.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = subprocess.run([str(SIMULATOR), "-control_file", str(probe)],
                            cwd=run_dir, capture_output=True, text=True, timeout=600)
    return result.stdout + result.stderr


@pytest.fixture(scope="module")
def staged(tmp_path_factory):
    """!
    @brief One staged run reused by every probe in this module.
    @param[in] tmp_path_factory Temporary directory factory fixture.
    @return Tuple of (run directory, base control file).
    """
    return _staged_run(tmp_path_factory.mktemp("guard"))


def test_authorized_run_root_is_refused(staged):
    """!
    @brief An authorized absolute path naming the run itself must not be deleted.
    @param[in] staged Staged run fixture.
    @return None.
    """
    run_dir, control = staged
    sentinel = run_dir / "sentinel.txt"
    sentinel.write_text("keep", encoding="utf-8")
    output = _probe(run_dir, control, str(run_dir.resolve()), authorized=True)
    assert "Refusing to delete log directory" in output
    assert "working directory itself" in output
    assert sentinel.read_text(encoding="utf-8") == "keep"


def test_authorized_ancestor_is_refused(staged):
    """!
    @brief An authorized absolute ancestor would take the run and its siblings with it.
    @param[in] staged Staged run fixture.
    @return None.
    """
    run_dir, control = staged
    for ancestor in (run_dir.parent, run_dir.parent.parent):
        output = _probe(run_dir, control, str(ancestor.resolve()), authorized=True)
        assert "Refusing to delete log directory" in output, ancestor
        assert "contains the working directory" in output, ancestor
        assert ancestor.is_dir()


def test_authorized_relative_symlink_escape_is_refused(staged, tmp_path):
    """!
    @brief A contained name that resolves out through a symlink stays refused.
    @param[in] staged Staged run fixture.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run_dir, control = staged
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "keep.txt").write_text("keep", encoding="utf-8")
    link = run_dir / "sneaky"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(outside)
    try:
        output = _probe(run_dir, control, "sneaky", authorized=True)
        assert "Refusing to delete log directory" in output
        assert "symlink" in output
        assert (outside / "keep.txt").read_text(encoding="utf-8") == "keep"
    finally:
        link.unlink()


@pytest.mark.parametrize("value,fragment", [
    ("..", "relative traversal"),
    ("config", "reserved run directory"),
    ("output", "overlaps the solver output directory"),
    ("~/logs", "nothing expands"),
    ("~", "nothing expands"),
])
def test_authorization_never_waives_these(staged, value, fragment):
    """!
    @brief Reserved names, overlap, and traversal stay fatal with authorization set.
    @param[in] staged Staged run fixture.
    @param[in] value Log-directory value to probe.
    @param[in] fragment Expected phrase in the refusal.
    @return None.
    """
    run_dir, control = staged
    output = _probe(run_dir, control, value, authorized=True)
    assert "Refusing to delete log directory" in output
    assert fragment in output


def test_contained_directory_is_accepted(staged):
    """!
    @brief The ordinary case must still run without an authorization.
    @param[in] staged Staged run fixture.
    @return None.
    """
    run_dir, control = staged
    output = _probe(run_dir, control, "logs", authorized=False)
    assert "Refusing to delete log directory" not in output


def test_external_absolute_is_accepted_only_with_authorization(staged, tmp_path):
    """!
    @brief The one waivable verdict: refused by default, permitted when authorized.
    @param[in] staged Staged run fixture.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run_dir, control = staged
    external = tmp_path / "archive"
    external.mkdir()
    refused = _probe(run_dir, control, str(external), authorized=False)
    assert "Refusing to delete log directory" in refused
    assert "no -allow_unsafe_log_dir authorization" in refused
    allowed = _probe(run_dir, control, str(external), authorized=True)
    assert "Refusing to delete log directory" not in allowed
