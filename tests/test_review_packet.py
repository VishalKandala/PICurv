"""!
@file test_review_packet.py
@brief Pytest coverage for page, registry, xref, and changed-set review-packet modes.
"""

import contextlib
import importlib.machinery
import importlib.util
import io
import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
REVIEW_PACKET = REPO_ROOT / "tests" / "tooling" / "review_packet.py"
CONTRACT_REGISTRY = REPO_ROOT / "tests" / "tooling" / "contract_registry.json"
CAPABILITY_REGISTRY = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
SUBSYSTEM_REGISTRY = REPO_ROOT / "tests" / "tooling" / "subsystem_records.json"
FRESHNESS_REGISTRY = REPO_ROOT / "tests" / "tooling" / "freshness_manifest.json"
PAGES_DIR = REPO_ROOT / "docs" / "pages"

PAGE_FRESHNESS_HEADING = "FRESHNESS SURFACES ROUTING REVIEW HERE"
CONTRACT_FRESHNESS_HEADING = "FRESHNESS SURFACES ON THIS CONTRACT'S PAGES"
ENTRY_ID = re.compile(r"^\s*\[[^\]]+\]\s+(\S+)")


def load_review_packet_module():
    """!
    @brief Load the review-packet assembler as an importable module.
    @return Value returned by `load_review_packet_module()`.
    """
    loader = importlib.machinery.SourceFileLoader("review_packet_module", str(REVIEW_PACKET))
    spec = importlib.util.spec_from_loader("review_packet_module", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def render(module, argv):
    """!
    @brief Render one review packet in-process and capture its output.
    @param[in] module Loaded review-packet module supplied to the function.
    @param[in] argv Argument list following the program name supplied to the function.
    @return Tuple of process status code and captured stdout text.
    """
    buffer = io.StringIO()
    saved = sys.argv
    sys.argv = ["review_packet.py"] + list(argv)
    try:
        with contextlib.redirect_stdout(buffer):
            status = module.main()
    finally:
        sys.argv = saved
    return status, buffer.getvalue()


def section_entry_ids(text, heading):
    """!
    @brief Collect the identifiers listed under one review-packet section heading.
    @param[in] text Captured review-packet output supplied to the function.
    @param[in] heading Section heading supplied to the function.
    @return Set of identifiers listed directly under the heading.
    """
    lines = text.splitlines()
    if heading not in lines:
        return set()
    found = set()
    for line in lines[lines.index(heading) + 1:]:
        if not line.strip():
            break
        match = ENTRY_ID.match(line)
        if match:
            found.add(match.group(1))
    return found


def without_section(text, heading):
    """!
    @brief Drop one section heading and its body from review-packet output.
    @param[in] text Captured review-packet output supplied to the function.
    @param[in] heading Section heading supplied to the function.
    @return Remaining output lines with that section removed.
    """
    lines = text.splitlines()
    if heading not in lines:
        return lines
    start = lines.index(heading)
    end = start + 1
    while end < len(lines) and lines[end].strip():
        end += 1
    return lines[:start] + lines[end:]


def published_page_stems():
    """!
    @brief List the published documentation page stems.
    @return Sorted page stems excluding the directory guide.
    """
    return sorted(p.stem for p in PAGES_DIR.glob("*.md") if p.stem != "guide")


def registered_contract_ids():
    """!
    @brief List every invariant contract identifier in the registry.
    @return Sorted contract identifiers.
    """
    registry = json.loads(CONTRACT_REGISTRY.read_text(encoding="utf-8"))
    return sorted(contract["id"] for contract in registry["contracts"])


def registry_ids(path, key):
    """!
    @brief Return sorted identifiers from one review-packet routing registry.
    @param[in] path Registry JSON path supplied to the function.
    @param[in] key Top-level record-list key supplied to the function.
    @return Sorted declared identifiers.
    """
    return sorted(record["id"] for record in json.loads(path.read_text(encoding="utf-8"))[key])


def test_page_mode_renders_every_published_page():
    """!
    @brief Verify page mode succeeds for every published page, including pages with no contracts.
    """
    module = load_review_packet_module()
    stems = published_page_stems()
    assert stems, "no published pages discovered"
    for stem in stems:
        status, text = render(module, [stem])
        assert status == 0, f"page mode failed for {stem}"
        assert text.startswith("REVIEW PACKET: docs/pages/"), f"unexpected header for {stem}"


def test_page_mode_covers_pages_with_zero_one_and_many_contracts():
    """!
    @brief Verify the three contract-count classes all render, and that all three classes exist.
    """
    module = load_review_packet_module()
    counts = {stem: len(module.contracts_for(stem)) for stem in published_page_stems()}
    zero = [stem for stem, count in counts.items() if count == 0]
    single = [stem for stem, count in counts.items() if count == 1]
    many = [stem for stem, count in counts.items() if count > 1]
    assert zero and single and many, "expected pages owning zero, one, and several contracts"
    for stem in (zero[0], single[0], many[0]):
        status, text = render(module, [stem])
        assert status == 0, f"page mode failed for {stem}"
        assert PAGE_FRESHNESS_HEADING in text


def test_representative_pages_render():
    """!
    @brief Verify the pages that exercised the page-mode contract-scope defect still render.
    """
    module = load_review_packet_module()
    for stem in ("16_Config_Extension_Playbook",
                 "20_Grid_Cell_Architecture_Guide",
                 "52_Run_Lifecycle_Guide"):
        assert (PAGES_DIR / f"{stem}.md").is_file(), f"{stem} is no longer published"
        status, _ = render(module, [stem])
        assert status == 0, f"page mode failed for {stem}"


def test_page_mode_reports_only_freshness_surfaces_routing_to_that_page():
    """!
    @brief Verify page mode lists exactly the page's own surfaces and no contract-scoped section.
    """
    module = load_review_packet_module()
    for stem in published_page_stems():
        _, text = render(module, [stem])
        expected = {surface["id"] for _, surface in module.freshness_for(stem)}
        assert section_entry_ids(text, PAGE_FRESHNESS_HEADING) == expected, stem
        assert CONTRACT_FRESHNESS_HEADING not in text, (
            f"{stem} reports freshness scoped to a contract rather than to the page")


def test_page_mode_output_is_independent_of_contract_iteration_order(monkeypatch):
    """!
    @brief Verify reversing the contract registry order does not change page-mode output.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    module = load_review_packet_module()
    multi = [stem for stem in published_page_stems() if len(module.contracts_for(stem)) > 1]
    assert multi, "expected at least one page owning several contracts"
    ordered = module.contracts_for
    for stem in multi:
        _, forward = render(module, [stem])
        monkeypatch.setattr(module, "contracts_for", lambda name: list(reversed(ordered(name))))
        _, backward = render(module, [stem])
        monkeypatch.undo()
        assert (without_section(forward, "INVARIANT CONTRACTS")
                == without_section(backward, "INVARIANT CONTRACTS")), (
            f"{stem} renders differently when its contracts are iterated in another order")


def test_contract_mode_renders_every_registered_contract():
    """!
    @brief Verify contract mode still renders each registered contract and reports its own surfaces.
    """
    module = load_review_packet_module()
    contract_ids = registered_contract_ids()
    assert contract_ids, "no contracts registered"
    scoped = 0
    for contract_id in contract_ids:
        status, text = render(module, ["--contract", contract_id])
        assert status == 0, f"contract mode failed for {contract_id}"
        assert text.startswith(f"REVIEW PACKET: contract {contract_id}\n")
        if CONTRACT_FRESHNESS_HEADING in text:
            scoped += 1
            assert section_entry_ids(text, CONTRACT_FRESHNESS_HEADING)
    assert scoped, "contract mode no longer reports contract-scoped freshness surfaces"


def test_contract_mode_rejects_an_unknown_contract():
    """!
    @brief Verify an unknown contract id exits non-zero rather than rendering a packet.
    """
    module = load_review_packet_module()
    status, text = render(module, ["--contract", "no.such.contract"])
    assert status == 2
    assert "REVIEW PACKET" not in text


def test_registry_modes_render_every_declared_identifier():
    """!
    @brief Verify capability, subsystem, and freshness joins never produce an empty success.
    """
    module = load_review_packet_module()
    modes = (
        ("--capability", CAPABILITY_REGISTRY, "families", "capability"),
        ("--subsystem", SUBSYSTEM_REGISTRY, "subsystems", "subsystem"),
        ("--surface", FRESHNESS_REGISTRY, "surfaces", "surface"),
    )
    for option, path, key, label in modes:
        identifiers = registry_ids(path, key)
        assert identifiers, f"no {label} identifiers discovered"
        for identifier in identifiers:
            status, text = render(module, [option, identifier])
            assert status == 0, f"{label} mode failed for {identifier}\n{text}"
            assert text.startswith(f"REVIEW PACKET: {label} {identifier}\n")
            assert module.COMPLETE in text
            assert "OPTIONAL SOURCE CROSS-REFERENCES" in text


def test_registry_modes_reject_unknown_identifiers():
    """!
    @brief Verify each registry selector differentiates invalid input from incomplete routing.
    """
    module = load_review_packet_module()
    for option in ("--capability", "--subsystem", "--surface"):
        status, text = render(module, [option, "no.such.identifier"])
        assert status == 2
        assert "REVIEW PACKET" not in text


def test_source_xref_stamp_distinguishes_missing_current_and_stale(tmp_path, monkeypatch):
    """!
    @brief Verify optional xrefs are consumed only when schema and dirty-byte stamps match.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    @param[in] monkeypatch Pytest monkeypatch fixture supplied to the function.
    """
    module = load_review_packet_module()
    index = tmp_path / "xref.json"
    monkeypatch.setattr(module, "XREF", index)
    assert module.load_xref() == ("missing", None)

    paths = module.xref_source_files(REPO_ROOT)
    payload = {
        "schema_version": module.XREF_SCHEMA_VERSION,
        "source_digest": module.xref_digest_files(REPO_ROOT, paths),
        "source_files": [path.as_posix() for path in paths],
        "doxyfile_digest": module.xref_file_digest(REPO_ROOT / "docs" / "Doxyfile"),
        "symbols": {},
    }
    index.write_text(json.dumps(payload), encoding="utf-8")
    state, loaded = module.load_xref()
    assert state == "current"
    assert loaded == payload

    payload["source_digest"] = "sha256:stale"
    index.write_text(json.dumps(payload), encoding="utf-8")
    assert module.load_xref() == ("stale", None)


def test_changed_path_classification_and_declared_routing():
    """!
    @brief Verify changed-source coverage distinguishes routed production from tooling and gaps.
    """
    module = load_review_packet_module()
    assert module.classify_path("src/Boundaries.c") == "production"
    routes = module.routes_for_path("src/Boundaries.c")
    assert "boundary.handler" in routes["capability"]
    assert "boundary.system" in routes["surface"]
    assert module.classify_path("tests/tooling/review_packet.py") == (
        "documentation/configuration/test/tooling"
    )
    assert module.classify_path("docs/generated/capability_inventory.json") == "generated"
    assert module.routes_for_path("src/not_declared_anywhere.c") == {}


def test_changed_mode_reports_coverage_before_changed_paths():
    """!
    @brief Verify working-tree routing is advisory, nonempty, and reports its route state first.
    """
    module = load_review_packet_module()
    status, text = render(module, ["--changed", "working-tree"])
    assert status in (0, 3)
    assert text.startswith("REVIEW PACKET: changed working-tree\nROUTE:")
    assert "CHANGED PATHS" in text
    assert render(module, ["--changed", "not-a-set"])[0] == 2
