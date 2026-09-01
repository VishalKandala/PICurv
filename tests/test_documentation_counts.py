"""Regressions for counts and statuses that prose restates from a registry.

A number written into prose drifts the moment the registry behind it changes, and the
freshness fingerprint routes review but does not do arithmetic. These checks compare
the two directly, so a contradiction fails rather than waits to be noticed.
"""

import json
import re
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PAGES = REPO_ROOT / "docs" / "pages"
TOOLING = REPO_ROOT / "tests" / "tooling"


def _read(path):
    """!
    @brief Read a UTF-8 text file.
    @param[in] path File to read.
    @return File contents.
    """
    return path.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def registry():
    """!
    @brief The invariant contract registry.
    @return Parsed registry document.
    """
    return json.loads(_read(TOOLING / "contract_registry.json"))


@pytest.fixture(scope="module")
def page71():
    """!
    @brief The invariant contracts page.
    @return Page text.
    """
    return _read(PAGES / "71_Invariant_Contracts.md")


def test_page_metadata_matches_the_contract_registry(registry, page71):
    """!
    @brief Page 71's metadata line must agree with the registry it describes.
    @param[in] registry Contract registry.
    @param[in] page71 Page text.
    @return None.
    """
    enforced = sum(1 for c in registry["contracts"] if c["status"] == "enforced")
    total = len(registry["contracts"])
    match = re.search(r"@pagemeta\{[^}]*?(\d+) enforced of (\d+) registered", page71)
    assert match, "page 71 metadata no longer states an enforced/registered count"
    assert (int(match.group(1)), int(match.group(2))) == (enforced, total)


def test_page_body_matches_the_contract_registry(registry, page71):
    """!
    @brief The register section must agree with the registry too.
    @param[in] registry Contract registry.
    @param[in] page71 Page text.
    @return None.
    """
    counts = {"enforced": 0, "tracked": 0, "planned": 0, "retired": 0}
    for contract in registry["contracts"]:
        counts[contract["status"]] += 1
    match = re.search(r"\*\*(\d+) enforced, (\d+) tracked,\s*\n?(\d+) planned\*\*", page71)
    assert match, "page 71 no longer states an enforced/tracked/planned breakdown"
    assert [int(g) for g in match.groups()] == [
        counts["enforced"], counts["tracked"], counts["planned"]
    ]


def test_page_metadata_and_body_agree_with_each_other(page71):
    """!
    @brief The two statements on the same page must not contradict.
    @param[in] page71 Page text.
    @return None.
    """
    meta = re.search(r"@pagemeta\{[^}]*?(\d+) enforced of", page71)
    body = re.search(r"\*\*(\d+) enforced,", page71)
    assert meta and body
    assert meta.group(1) == body.group(1)


def test_logical_identity_count_matches_the_snapshot(page71):
    """!
    @brief The identity count in prose must match the generated snapshot.
    @param[in] page71 Page text.
    @return None.
    """
    snapshot = json.loads(_read(REPO_ROOT / "docs" / "generated" /
                                "artifact_topology_snapshot.json"))
    identities = snapshot["logical_artifacts"]
    match = re.search(r"names (\w+)\s*\n?\s*logical identities", page71)
    assert match, "page 71 no longer states a logical identity count"
    assert match.group(1) == str(len(identities))
    for identity in identities:
        assert f"`{identity}`" in page71, identity


def test_enumerated_containment_points_match_the_stated_number(page71):
    """!
    @brief "Checked at N points" must be followed by exactly N enumerated points.
    @param[in] page71 Page text.
    @return None.
    """
    words = {"three": 3, "four": 4, "five": 5, "six": 6, "seven": 7}
    match = re.search(r"Containment is checked at (\w+) points:", page71)
    assert match, "page 71 no longer enumerates containment points"
    stated = words[match.group(1)]
    tail = page71[match.end():]
    # Count only the consecutive run 1..N that immediately follows, so a later
    # numbered section elsewhere on the page cannot inflate the total.
    enumerated = 0
    while re.search(rf"^{enumerated + 1}\. ", tail, re.M):
        enumerated += 1
    assert enumerated == stated, f"stated {stated}, enumerated {enumerated}"


def test_no_family_is_described_as_awaiting_backfill():
    """!
    @brief The registry must not claim a backfill that has completed.
    @return None.
    """
    document = json.loads(_read(TOOLING / "capability_families.json"))
    advisory = [f["id"] for f in document["families"] if not f.get("coverage_enforced")]
    prose = " ".join(str(v) for v in document.values() if isinstance(v, str))
    prose += " ".join(str(f.get("coverage_note", "")) for f in document["families"])
    if not advisory:
        # "nothing is awaiting backfill" is a statement that the backfill is done, not
        # a claim that one is pending; the failing phrasing is the claim of a family
        # still waiting for one.
        assert "still awaiting backfill" not in prose
        assert "Families still awaiting" not in prose
        assert "Advisory until" not in prose
        assert "backfill pending" not in prose


def test_census_status_line_matches_its_own_classifications():
    """!
    @brief The census summary must agree with the classifications beneath it.
    @return None.
    """
    census = json.loads(_read(TOOLING / "family_census.json"))
    pending = [k for k, v in census["acknowledged_uncovered"].items()
               if isinstance(v, dict) and v.get("classification") == "public_pending"]
    if not pending:
        assert "pending backfill" not in census["status"]


def test_capability_status_page_does_not_claim_to_be_unenforced():
    """!
    @brief Page 62 must not describe an enforced contract as a target.
    @return None.
    """
    text = _read(PAGES / "62_Capability_Status_Vocabulary.md")
    assert "not yet enforced" not in text
    assert "Contract status: enforced." in text
