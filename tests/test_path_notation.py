"""Regressions for the run-path notation contract.

The audit used to match only `runs/<id>/logs/...`. It therefore passed while pages
hardcoded `logs/Runtime_Memory.log` - the default location of a **configurable**
directory, and so exactly the thing that breaks when someone relocates their storage.
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


literals = _load("audit_path_literals")


def _violations(markdown, tmp_path, monkeypatch):
    """!
    @brief Run the audit over a single synthetic page.
    @param[in] markdown Page text.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return Captured stderr from the audit.
    """
    page = tmp_path / "page.md"
    page.write_text(markdown, encoding="utf-8")
    monkeypatch.setattr(literals, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(literals, "TOPOLOGY", TOOLING / "artifact_topology.json")
    monkeypatch.setattr(literals, "tracked_markdown", lambda: [page])
    import io
    import contextlib
    stream = io.StringIO()
    with contextlib.redirect_stderr(stream):
        code = literals.main()
    return code, stream.getvalue()


@pytest.mark.parametrize("line", [
    "The solver writes logs/Runtime_Memory.log every step.",
    "Look in logs/Continuity_Metrics.log for the divergence history.",
    "Checkpoints land in output/checkpoints/step_000000000100.",
    # Terminal directory references: the form the pattern used to let through
    # entirely, because it required a character after the slash.
    "Normalized runtime config artifacts live under `config/`.",
    "Solver and post logs are written under `logs/`.",
    "Solver field output lands in `output/`.",
    "Scheduler scripts and submission metadata are under `scheduler/`.",
    "Rendered frames are written to `visualization/`.",
])
def test_a_bare_run_owned_prefix_is_rejected(line, tmp_path, monkeypatch):
    """!
    @brief Prose naming a configurable run directory concretely must fail.
    @param[in] line The offending prose.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(line + "\n", tmp_path, monkeypatch)
    assert code == 1
    assert "bare run-owned prefix" in output


def test_a_placeholder_run_root_is_rejected(tmp_path, monkeypatch):
    """!
    @brief `<run_dir>/logs/...` still hardcodes the subdirectory.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "Errors are written to <run_dir>/logs/interpolation_error.csv.\n",
        tmp_path, monkeypatch)
    assert code == 1
    assert "unresolved run root" in output


def test_the_repository_log_directory_stays_valid(tmp_path, monkeypatch):
    """!
    @brief `<repo>/logs/...` must survive: it is a different directory.

    @details The repository's build log and a run's runtime log both read as `logs/`
             in prose. One moves with configuration and one does not, so the notation
             has to tell them apart - and the escape hatch has to actually work.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "`picurv build` writes `<repo>/logs/build.log` and `<repo>/logs/doxygen.warnings`.\n",
        tmp_path, monkeypatch)
    assert code == 0, output


def test_logical_identities_are_accepted(tmp_path, monkeypatch):
    """!
    @brief The prescribed notation must pass.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "Memory diagnostics land in `<run.runtime_logs>/Runtime_Memory.log`, and "
        "checkpoints in `<run.solver_output>/checkpoints/`.\n", tmp_path, monkeypatch)
    assert code == 0, output


def test_runnable_command_blocks_keep_concrete_paths(tmp_path, monkeypatch):
    """!
    @brief A reader has to be able to type something.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "Run it:\n\n```bash\ntail -f logs/Runtime_Memory.log\n```\n",
        tmp_path, monkeypatch)
    assert code == 0, output


def test_the_repository_config_library_is_not_a_run_directory(tmp_path, monkeypatch):
    """!
    @brief `config/solvers` is shipped source, not a run's config directory.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "Reusable profiles live in `config/solvers` and `config/monitors`.\n",
        tmp_path, monkeypatch)
    assert code == 0, output


def test_the_repository_passes_the_strengthened_audit():
    """!
    @brief Every page in the repository already satisfies the notation contract.
    @return None.
    """
    result = subprocess.run([sys.executable, str(TOOLING / "audit_path_literals.py")],
                            cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_terminal_directory_literals_are_caught(tmp_path, monkeypatch):
    """!
    @brief A directory named with no file after it drifts exactly the same way.

    @details `BARE_PREFIX` required a character after the slash, so `logs/` on its own
             passed while `logs/x.log` failed - and a bare directory is the form pages
             reach for when describing a layout, which is the thing a relocation breaks.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    for directory in ("config", "logs", "output", "scheduler", "visualization",
                      "checkpoints"):
        code, output = _violations(
            f"Artifacts are written under `{directory}/` during a run.\n",
            tmp_path, monkeypatch)
        assert code == 1, f"terminal `{directory}/` was not caught"
        assert directory in output


def test_repository_owned_directories_survive_with_the_prefix(tmp_path, monkeypatch):
    """!
    @brief `<repo>/config/` and `<repo>/logs/` name the shipped tree, not a run.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "Reusable profiles live in `<repo>/config/`, and build output in `<repo>/logs/`.\n",
        tmp_path, monkeypatch)
    assert code == 0, output


def test_a_source_subdirectory_needs_no_prefix(tmp_path, monkeypatch):
    """!
    @brief `config/solvers/` is unambiguous already - the subdirectory identifies it.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned tree.
    @return None.
    """
    code, output = _violations(
        "Profiles live in `config/solvers/` and `config/monitors/`.\n",
        tmp_path, monkeypatch)
    assert code == 0, output


def test_the_notation_contract_is_documented():
    """!
    @brief Page 71 must state the rule the audit enforces.
    @return None.
    """
    text = (REPO_ROOT / "docs" / "pages" / "71_Invariant_Contracts.md").read_text(
        encoding="utf-8")
    assert "<run.runtime_logs>" in text
    assert "<repo>/logs/build.log" in text
    assert "p71_notation_sub" in text
