"""Regressions for the named-choice-set contract.

The census discovers named constants. This audit is the other half: it finds inline
literals the census cannot see, so "exhaustive" means something. These tests check that
it actually catches a new one.
"""

import importlib.util
import json
import subprocess
import sys
import textwrap
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


inline = _load("audit_inline_choices")


@pytest.fixture(scope="module")
def waivers():
    """!
    @brief The committed waiver document.
    @return Parsed document.
    """
    return json.loads((TOOLING / "inline_choice_waivers.json").read_text(encoding="utf-8"))


def test_the_repository_passes_the_audit():
    """!
    @brief Every inline choice set in the CLI package is classified.
    @return None.
    """
    result = subprocess.run([sys.executable, str(TOOLING / "audit_inline_choices.py")],
                            cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_an_unclassified_inline_choice_is_rejected(tmp_path, monkeypatch):
    """!
    @brief A public closed choice added inline must fail the audit.

    @details This is the gap the audit exists to close: the census could call itself
             exhaustive while a set like this sat in a literal nobody enumerated.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned package.
    @return None.
    """
    package = tmp_path / "picurv_cli"
    package.mkdir()
    (package / "core.py").write_text(textwrap.dedent('''
        def validate_wall_model(value):
            """A new public choice, written inline."""
            if value not in {"algebraic", "ode", "equilibrium"}:
                raise ValueError("bad wall model")
            return value
    '''), encoding="utf-8")
    monkeypatch.setattr(inline, "REPO_ROOT", tmp_path)
    found = inline.find_inline_choices()
    key = inline.literal_key("picurv_cli.core", ["algebraic", "ode", "equilibrium"])
    assert key in found


def test_an_inline_argparse_choices_list_is_detected(tmp_path, monkeypatch):
    """!
    @brief An unnamed argparse `choices=` literal is detected too.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned package.
    @return None.
    """
    package = tmp_path / "picurv_cli"
    package.mkdir()
    (package / "cli.py").write_text(textwrap.dedent('''
        def build(parser):
            """Add a flag with unnamed choices."""
            parser.add_argument("--flavour", choices=["vanilla", "chocolate"])
    '''), encoding="utf-8")
    monkeypatch.setattr(inline, "REPO_ROOT", tmp_path)
    found = inline.find_inline_choices()
    assert inline.literal_key("picurv_cli.cli", ["vanilla", "chocolate"]) in found


def test_naming_the_set_removes_it_from_the_findings(tmp_path, monkeypatch):
    """!
    @brief Giving the set a name is the fix, and the audit recognises it.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned package.
    @return None.
    """
    package = tmp_path / "picurv_cli"
    package.mkdir()
    (package / "core.py").write_text(textwrap.dedent('''
        WALL_MODEL_TYPES = ("algebraic", "ode", "equilibrium")


        def validate_wall_model(value):
            """A new public choice, now named."""
            if value not in WALL_MODEL_TYPES:
                raise ValueError("bad wall model")
            return value
    '''), encoding="utf-8")
    monkeypatch.setattr(inline, "REPO_ROOT", tmp_path)
    assert inline.find_inline_choices() == {}


def test_a_named_set_becomes_visible_to_the_census(tmp_path, monkeypatch):
    """!
    @brief The two audits meet: what one stops being blind to, the other must classify.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned package.
    @return None.
    """
    census = _load("audit_family_census")
    package = tmp_path / "picurv_cli"
    package.mkdir()
    (package / "core.py").write_text(
        'WALL_MODEL_TYPES = ("algebraic", "ode", "equilibrium")\n', encoding="utf-8")
    (package / "cli.py").write_text("", encoding="utf-8")
    (package / "storage.py").write_text("", encoding="utf-8")
    monkeypatch.setattr(census, "REPO_ROOT", tmp_path)
    assert "picurv_cli.core::WALL_MODEL_TYPES" in census.candidate_surfaces()


def test_single_member_sets_are_not_treated_as_choices(tmp_path, monkeypatch):
    """!
    @brief A one-element literal is a comparison, not a menu.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the scanned package.
    @return None.
    """
    package = tmp_path / "picurv_cli"
    package.mkdir()
    (package / "core.py").write_text(
        'def f(v):\n    """Doc."""\n    return v in {"only"}\n', encoding="utf-8")
    monkeypatch.setattr(inline, "REPO_ROOT", tmp_path)
    assert inline.find_inline_choices() == {}


def test_every_waiver_is_typed_and_reasoned(waivers):
    """!
    @brief A waiver without a real reason is not a classification.
    @param[in] waivers The committed waiver document.
    @return None.
    """
    for key, entry in waivers["inline_choices"].items():
        assert entry["classification"] in inline.VALID_CLASSIFICATIONS, key
        assert len(entry["reason"].strip()) >= inline.MIN_REASON_CHARS, key


def test_no_waiver_claims_a_pending_family(waivers):
    """!
    @brief An inline literal cannot be a pending family; it has no surface to name.
    @param[in] waivers The committed waiver document.
    @return None.
    """
    assert "public_pending" not in inline.VALID_CLASSIFICATIONS
    for entry in waivers["inline_choices"].values():
        assert entry["classification"] != "public_pending"


def test_the_cli_choice_policy_is_stated(waivers):
    """!
    @brief The intended policy for argparse commands and choices must be written down.
    @param[in] waivers The committed waiver document.
    @return None.
    """
    policy = waivers["cli_choices_policy"]
    assert "cli.reference" in policy
    assert "when to use" in policy
