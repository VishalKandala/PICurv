"""Negative tests for the subsystem lifecycle gate.

The gate is only worth having if it rejects things. Each test here constructs a
record that a careless author could plausibly write and asserts the audit refuses it.
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


lifecycle = _load("audit_subsystem_lifecycle")


@pytest.fixture(scope="module")
def pages():
    """!
    @brief The real page/anchor index, so references resolve as they do in the audit.
    @return Mapping of page id to anchors.
    """
    return lifecycle.page_index()


@pytest.fixture(scope="module")
def families():
    """!
    @brief The real capability families, keyed by id.
    @return Mapping of family id to family record.
    """
    document = json.loads((TOOLING / "capability_families.json").read_text(encoding="utf-8"))
    return {family["id"]: family for family in document["families"]}


@pytest.fixture(scope="module")
def records():
    """!
    @brief The committed subsystem records.
    @return List of records.
    """
    document = json.loads((TOOLING / "subsystem_records.json").read_text(encoding="utf-8"))
    return document["subsystems"]


@pytest.fixture()
def supported(records):
    """!
    @brief A deep copy of a real `supported` record, for mutation.
    @param[in] records The committed subsystem records.
    @return One supported subsystem record.
    """
    return copy.deepcopy(next(r for r in records if r["status"] == "supported"))


@pytest.fixture(scope="module")
def published():
    """!
    @brief The ids of the pages the site publishes.
    @return Set of page ids.
    """
    return set(
        json.loads((TOOLING / "page_types.json").read_text(encoding="utf-8"))["assignments"]
    )


def _violations(record, pages, families, published):
    """!
    @brief Validate a single record and return the violation strings.
    @param[in] record One subsystem record.
    @param[in] pages Page/anchor index.
    @param[in] families Capability families by id.
    @param[in] published Ids of published pages.
    @return List of violations.
    """
    return lifecycle.validate([record], pages, families, published)


def test_committed_records_pass(records, pages, families, published):
    """!
    @brief The records in the repository satisfy the contract as committed.
    @param[in] records Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    assert lifecycle.validate(records, pages, families, published) == []


def test_supported_owes_the_full_ladder(supported, pages, families, published):
    """!
    @brief A supported subsystem cannot omit an obligation its own rung requires.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"].pop("troubleshooting")
    problems = _violations(supported, pages, families, published)
    assert any("owes 'troubleshooting'" in p for p in problems)


def test_supported_inherits_lower_rung_obligations(supported, pages, families, published):
    """!
    @brief Claiming a rung owes every obligation below it as well.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"].pop("scope_boundary")
    problems = _violations(supported, pages, families, published)
    assert any("owes 'scope_boundary'" in p for p in problems)


def test_planned_does_not_owe_supported_obligations(records, pages, families, published):
    """!
    @brief A planned subsystem is not blocked by documentation it cannot honestly write.
    @param[in] records Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    planned = copy.deepcopy(next(r for r in records if r["status"] == "planned"))
    assert _violations(planned, pages, families, published) == []


def test_claiming_supported_from_planned_documentation_fails(records, pages, families, published):
    """!
    @brief Raising the claimed status raises what must already be written.
    @param[in] records Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    planned = copy.deepcopy(next(r for r in records if r["status"] == "planned"))
    planned["status"] = "supported"
    problems = _violations(planned, pages, families, published)
    assert any("owes 'evidence'" in p for p in problems)


def test_skipping_deprecation_on_the_way_to_removed_is_invalid(supported, pages, families, published):
    """!
    @brief Nothing reaches removed without first being deprecated or declared defective.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["status"] = "removed"
    supported["peak_status"] = "supported"
    problems = _violations(supported, pages, families, published)
    assert any("not a valid lifecycle transition" in p for p in problems)


def test_reviving_a_removed_subsystem_is_invalid(supported, pages, families, published):
    """!
    @brief A removed subsystem is terminal; reintroduction is a new record.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["previous_status"] = "removed"
    problems = _violations(supported, pages, families, published)
    assert any("not a valid lifecycle transition" in p for p in problems)


def test_a_brand_new_record_cannot_start_removed(supported, pages, families, published):
    """!
    @brief A subsystem that never existed cannot be documented as removed.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["previous_status"] = None
    supported["status"] = "removed"
    supported["peak_status"] = "supported"
    problems = _violations(supported, pages, families, published)
    assert any("not a valid lifecycle transition" in p for p in problems)


def test_defective_status_requires_a_peak(records, pages, families, published):
    """!
    @brief Off-ladder statuses must say how far up the ladder they got.
    @param[in] records Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    defective = copy.deepcopy(next(r for r in records if r["status"] == "known-defective"))
    defective.pop("peak_status")
    problems = _violations(defective, pages, families, published)
    assert any("requires peak_status" in p for p in problems)


def test_bare_not_applicable_is_rejected(supported, pages, families, published):
    """!
    @brief A reasoned N/A is a real answer; a bare one is an evasion.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"]["troubleshooting"] = {"not_applicable": "N/A"}
    problems = _violations(supported, pages, families, published)
    assert any("without a stated reason" in p for p in problems)


def test_reasoned_not_applicable_is_accepted(supported, pages, families, published):
    """!
    @brief A stated reason satisfies an obligation without inventing prose.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"]["troubleshooting"] = {
        "not_applicable": "This subsystem has no failure mode a user can act on; every "
        "rejection is a validation error covered by the fatal-error page."
    }
    assert _violations(supported, pages, families, published) == []


def test_dangling_page_reference_is_rejected(supported, pages, families, published):
    """!
    @brief An obligation pointing at a page that does not exist is unsatisfied.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"]["operations"] = {"page": "99_Does_Not_Exist"}
    problems = _violations(supported, pages, families, published)
    assert any("does not exist" in p for p in problems)


def test_dangling_anchor_is_rejected(supported, pages, families, published):
    """!
    @brief An anchor must resolve, not merely name a real page.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    reference = supported["obligations"]["operations"]
    supported["obligations"]["operations"] = {
        "page": reference["page"], "anchor": "p99_not_a_real_anchor_sec"
    }
    problems = _violations(supported, pages, families, published)
    assert any("which does not define it" in p for p in problems)


def test_empty_obligation_object_is_rejected(supported, pages, families, published):
    """!
    @brief An obligation declared but not answered fails.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"]["examples"] = {}
    problems = _violations(supported, pages, families, published)
    assert any("is unsatisfied" in p for p in problems)


def test_internal_status_cannot_be_publicly_visible(supported, pages, families, published):
    """!
    @brief Internal and publicly reachable is a contradiction.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["status"] = "internal"
    supported["previous_status"] = "planned"
    problems = _violations(supported, pages, families, published)
    assert any("contradicts visibility" in p for p in problems)


def test_unknown_concern_is_rejected(supported, pages, families, published):
    """!
    @brief The concern vocabulary is closed.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["concerns"]["quantum_entanglement"] = {"page": "13_Code_Architecture"}
    problems = _violations(supported, pages, families, published)
    assert any("not in the concern vocabulary" in p for p in problems)


def test_planned_subsystem_may_not_own_supported_values(records, pages, families, published):
    """!
    @brief A planned subsystem must not already appear as working behavior.
    @param[in] records Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    planned = copy.deepcopy(next(r for r in records if r["status"] == "planned"))
    planned["capability_families"] = ["momentum.solver"]
    problems = _violations(planned, pages, families, published)
    assert any("must not appear as supported behavior" in p for p in problems)


def test_unregistered_capability_family_is_rejected(supported, pages, families, published):
    """!
    @brief A record cannot cite a family the registry does not know.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["capability_families"] = ["nonexistent.family"]
    problems = _violations(supported, pages, families, published)
    assert any("is not registered" in p for p in problems)


def test_obligation_a_status_does_not_owe_is_rejected(supported, pages, families, published):
    """!
    @brief Declaring an obligation of a higher status does not confer that status.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["obligations"]["migration_path"] = {"page": "13_Code_Architecture"}
    problems = _violations(supported, pages, families, published)
    assert any("does not\nowe" in p or "does not owe" in p for p in problems)


def test_unknown_record_field_is_rejected(supported, pages, families, published):
    """!
    @brief The record schema fails closed on unknown fields.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    supported["mood"] = "confident"
    problems = _violations(supported, pages, families, published)
    assert any("unknown field" in p for p in problems)


def test_duplicate_ids_are_rejected(supported, pages, families, published):
    """!
    @brief Two records cannot claim the same subsystem id.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    problems = lifecycle.validate(
        [supported, copy.deepcopy(supported)], pages, families, published
    )
    assert any("declared more than once" in p for p in problems)


def test_public_subsystem_may_not_cite_an_unpublished_page(supported, pages, families,
                                                           published):
    """!
    @brief A public subsystem cannot route its reader to an unpublished page.
    @param[in] supported Fixture.
    @param[in] pages Fixture.
    @param[in] families Fixture.
    @param[in] published Fixture.
    @return None.
    """
    trimmed = set(published) - {supported["obligations"]["operations"]["page"]}
    problems = lifecycle.validate([supported], pages, families, trimmed)
    assert any("not a published page" in p for p in problems)


def test_every_status_in_the_vocabulary_has_an_obligation_set():
    """!
    @brief No status in the vocabulary owes nothing.
    @return None.
    """
    for status in lifecycle.VALID_STATUSES:
        record = {"status": status, "peak_status": "experimental"}
        assert lifecycle.required_obligations(record), status


def test_cancelling_a_planned_subsystem_is_allowed(records, pages, families, published):
    """!
    @brief A design that was never built may be cancelled straight to removed.
    @param[in] records The committed subsystem records.
    @param[in] pages Page/anchor index.
    @param[in] families Capability families by id.
    @param[in] published Ids of published pages.
    @return None.
    """
    assert "removed" in lifecycle.TRANSITIONS["planned"]


def test_every_other_route_to_removed_passes_through_deprecation():
    """!
    @brief Only `planned` and the two terminal statuses may precede `removed`.
    @return None.
    """
    allowed_from = {status for status, targets in lifecycle.TRANSITIONS.items()
                    if status is not None and "removed" in targets}
    assert allowed_from == {"planned", "known-defective", "deprecated", "removed"}


def test_a_promotion_must_state_what_the_owner_is_deciding(supported, pages, families,
                                                           published):
    """!
    @brief A proposed promotion without a rationale is noise.
    @param[in] supported A supported subsystem record.
    @param[in] pages Page/anchor index.
    @param[in] families Capability families by id.
    @param[in] published Ids of published pages.
    @return None.
    """
    supported["status"] = "experimental"
    supported["proposed_status"] = "supported"
    problems = _violations(supported, pages, families, published)
    assert any("no promotion_rationale" in p for p in problems)


def test_a_proposal_must_be_a_promotion(supported, pages, families, published):
    """!
    @brief Proposing a status at or below the claimed one is rejected.
    @param[in] supported A supported subsystem record.
    @param[in] pages Page/anchor index.
    @param[in] families Capability families by id.
    @param[in] published Ids of published pages.
    @return None.
    """
    supported["proposed_status"] = "experimental"
    supported["promotion_rationale"] = "A rationale that should not save this."
    problems = _violations(supported, pages, families, published)
    assert any("not above the claimed status" in p for p in problems)


def test_obligations_are_checked_against_the_proposal(records, pages, families, published):
    """!
    @brief A deferred promotion is still held to the higher documentation bar.
    @param[in] records The committed subsystem records.
    @param[in] pages Page/anchor index.
    @param[in] families Capability families by id.
    @param[in] published Ids of published pages.
    @return None.
    """
    proposals = [r for r in records if r.get("proposed_status")]
    assert proposals, "this branch defers at least one promotion"
    for record in proposals:
        owed = lifecycle.required_obligations(record)
        assert "evidence" in owed and "troubleshooting" in owed, record["id"]
        broken = copy.deepcopy(record)
        broken["obligations"].pop("evidence")
        assert any("owes 'evidence'" in p
                   for p in _violations(broken, pages, families, published)), record["id"]


def test_no_promotion_is_recorded_as_already_taken(records):
    """!
    @brief Every promotion this branch introduced is still pending, not asserted.
    @param[in] records The committed subsystem records.
    @return None.
    """
    for record in records:
        if record.get("proposed_status"):
            assert record["status"] != record["proposed_status"], record["id"]
