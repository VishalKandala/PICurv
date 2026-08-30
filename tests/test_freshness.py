"""Tests for the tiered freshness mechanism.

The mechanism is only useful if hard suspicion blocks, soft suspicion does not, and a
broken manifest is neither of those but an integrity failure.
"""

import copy
import importlib.util
import json
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


freshness = _load("audit_freshness")


@pytest.fixture(scope="module")
def manifest():
    """!
    @brief The committed freshness manifest.
    @return Parsed manifest document.
    """
    return json.loads((TOOLING / "freshness_manifest.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def published():
    """!
    @brief Ids of the published pages.
    @return Set of page ids.
    """
    return set(
        json.loads((TOOLING / "page_types.json").read_text(encoding="utf-8"))["assignments"]
    )


@pytest.fixture()
def surfaces(manifest):
    """!
    @brief A mutable copy of the manifest surfaces.
    @param[in] manifest The committed freshness manifest.
    @return List of surface records.
    """
    return copy.deepcopy(manifest["surfaces"])


def _hard(surfaces):
    """!
    @brief The first hard-tier surface.
    @param[in] surfaces Manifest surfaces.
    @return One surface record.
    """
    return next(s for s in surfaces if s["tier"] == "hard")


def _soft(surfaces):
    """!
    @brief The first advisory soft-tier surface.
    @param[in] surfaces Manifest surfaces.
    @return One surface record.
    """
    return next(s for s in surfaces if s["tier"] == "soft"
                and freshness.enforcement_of(s) == "report")


def test_committed_manifest_is_structurally_sound(surfaces, published):
    """!
    @brief The committed manifest has no integrity failures.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    assert freshness.validate_manifest(surfaces, published) == []


def test_committed_manifest_has_no_blocking_suspicion(surfaces, published):
    """!
    @brief The repository is currently free of blocking staleness.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    report = freshness.evaluate(surfaces, published)
    assert report["integrity_failures"] == []
    assert freshness.blocking_entries(report) == []


def test_hard_drift_blocks(surfaces, published):
    """!
    @brief A changed normalized artifact fails the build.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _hard(surfaces)["attested_digest"] = "sha256:" + "0" * 64
    report = freshness.evaluate(surfaces, published)
    assert report["hard_suspect"]
    assert freshness.blocking_entries(report)


def test_soft_drift_does_not_block(surfaces, published):
    """!
    @brief A raw-source change is advisory, not a failure.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    surface = _soft(surfaces)
    surface["attested_digest"] = "sha256:" + "0" * 64
    report = freshness.evaluate(surfaces, published)
    assert any(e["id"] == surface["id"] for e in report["soft_suspect"])
    # Only this surface's verdict is asserted: another surface may legitimately be
    # blocking in the committed manifest, and that is a different test's business.
    assert surface["id"] not in {e["id"] for e in freshness.blocking_entries(report)}


def test_promoted_soft_surface_blocks(surfaces, published):
    """!
    @brief A soft surface may be promoted to blocking when it guards something critical.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    surface = _soft(surfaces)
    surface["enforcement"] = "blocking"
    surface["promotion_reason"] = "This surface guards something destructive."
    surface["attested_digest"] = "sha256:" + "0" * 64
    report = freshness.evaluate(surfaces, published)
    assert freshness.blocking_entries(report)


def test_promotion_without_a_reason_is_an_integrity_failure(surfaces, published):
    """!
    @brief Promotion to blocking must state why.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    surface = _soft(surfaces)
    surface["enforcement"] = "blocking"
    problems = freshness.validate_manifest(surfaces, published)
    assert any("must state why" in p for p in problems)


def test_missing_watched_path_is_an_integrity_failure(surfaces, published):
    """!
    @brief A vanished watched path is broken coverage, not staleness.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _soft(surfaces)["watched_paths"] = ["src/this_file_was_deleted.c"]
    problems = freshness.validate_manifest(surfaces, published)
    assert any("coverage is broken" in p for p in problems)


def test_integrity_failure_suppresses_the_ordinary_report(surfaces, published):
    """!
    @brief With coverage broken, the absence of suspicion proves nothing.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _soft(surfaces)["watched_paths"] = ["src/this_file_was_deleted.c"]
    report = freshness.evaluate(surfaces, published)
    assert report["integrity_failures"]
    assert report["hard_current"] == [] and report["soft_current"] == []


def test_a_surface_with_no_owning_page_is_rejected(surfaces, published):
    """!
    @brief A fingerprint must have somewhere to route its suspicion.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _hard(surfaces)["owning_pages"] = []
    problems = freshness.validate_manifest(surfaces, published)
    assert any("nowhere to route its suspicion" in p for p in problems)


def test_unpublished_owning_page_is_rejected(surfaces, published):
    """!
    @brief Review cannot be routed to a page the site does not publish.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _hard(surfaces)["owning_pages"] = ["99_No_Such_Page"]
    problems = freshness.validate_manifest(surfaces, published)
    assert any("is not published" in p for p in problems)


def test_hard_surface_must_name_its_regeneration_command(surfaces, published):
    """!
    @brief A hard surface must tell the reader how to make it current.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _hard(surfaces).pop("regenerate")
    problems = freshness.validate_manifest(surfaces, published)
    assert any("regenerates it" in p for p in problems)


def test_unattested_hard_surface_blocks(surfaces, published):
    """!
    @brief A hard surface nobody has ever reviewed fails the build.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    _hard(surfaces)["attested_digest"] = freshness.UNATTESTED
    report = freshness.evaluate(surfaces, published)
    assert report["unattested"]
    assert freshness.blocking_entries(report)


def test_unattested_soft_surface_does_not_block(surfaces, published):
    """!
    @brief A never-reviewed soft surface is reported honestly, not enforced.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    surface = _soft(surfaces)
    surface["attested_digest"] = freshness.UNATTESTED
    report = freshness.evaluate(surfaces, published)
    assert any(e["id"] == surface["id"] for e in report["unattested"])
    assert surface["id"] not in {e["id"] for e in freshness.blocking_entries(report)}


def test_digest_covers_the_path_not_only_the_bytes(tmp_path, monkeypatch):
    """!
    @brief Identical content under a different name is a different fingerprint.
    @param[in] tmp_path Fixture.
    @param[in] monkeypatch Fixture.
    @return None.
    """
    (tmp_path / "first.c").write_bytes(b"same content")
    (tmp_path / "second.c").write_bytes(b"same content")
    monkeypatch.setattr(freshness, "REPO_ROOT", tmp_path)
    assert freshness.digest_of(["first.c"]) != freshness.digest_of(["second.c"])


def test_digest_is_order_independent_and_deterministic(tmp_path, monkeypatch):
    """!
    @brief Declaration order must not change a surface's fingerprint.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the repository root.
    @return None.
    """
    (tmp_path / "a.c").write_bytes(b"alpha")
    (tmp_path / "b.c").write_bytes(b"beta")
    monkeypatch.setattr(freshness, "REPO_ROOT", tmp_path)
    assert freshness.digest_of(["a.c", "b.c"]) == freshness.digest_of(["b.c", "a.c"])


def test_moving_content_between_watched_files_changes_the_digest(tmp_path, monkeypatch):
    """!
    @brief Relocating content between watched files is itself a change.
    @param[in] tmp_path Temporary directory fixture.
    @param[in] monkeypatch Fixture used to relocate the repository root.
    @return None.
    """
    monkeypatch.setattr(freshness, "REPO_ROOT", tmp_path)
    (tmp_path / "a.c").write_bytes(b"body")
    (tmp_path / "b.c").write_bytes(b"")
    before = freshness.digest_of(["a.c", "b.c"])
    (tmp_path / "a.c").write_bytes(b"")
    (tmp_path / "b.c").write_bytes(b"body")
    assert freshness.digest_of(["a.c", "b.c"]) != before


def test_attestation_records_the_current_digest(surfaces):
    """!
    @brief Re-attesting stores the digest as it is now.
    @param[in] surfaces Fixture.
    @return None.
    """
    surface = _hard(surfaces)
    surface["attested_digest"] = "sha256:" + "0" * 64
    changed = freshness.attest(surfaces, [surface["id"]])
    assert changed == [surface["id"]]
    assert surface["attested_digest"] == freshness.digest_of(freshness.surface_paths(surface))


def test_attesting_all_touches_every_surface(surfaces):
    """!
    @brief The 'all' form re-attests the whole manifest.
    @param[in] surfaces Fixture.
    @return None.
    """
    for surface in surfaces:
        surface["attested_digest"] = "sha256:" + "0" * 64
    changed = freshness.attest(surfaces, ["all"])
    assert len(changed) == len(surfaces)


def test_no_page_stores_a_commit_sha_for_attestation():
    """!
    @brief Attestation state lives in the manifest, never inside a page.
    @return None.
    """
    for page in sorted((REPO_ROOT / "docs" / "pages").glob("*.md")):
        text = page.read_text(encoding="utf-8")
        assert "reviewed_at:" not in text, page.name
        assert "attested_digest" not in text, page.name


def test_every_owning_page_appears_in_the_routed_review_list(surfaces, published):
    """!
    @brief A suspicion names the pages that owe review.
    @param[in] surfaces Fixture.
    @param[in] published Fixture.
    @return None.
    """
    surface = _hard(surfaces)
    surface["attested_digest"] = "sha256:" + "0" * 64
    report = freshness.evaluate(surfaces, published)
    entry = next(e for e in report["hard_suspect"] if e["id"] == surface["id"])
    for page in surface["owning_pages"]:
        assert page in entry["pages"]
