#!/usr/bin/env python3
"""Regression tests for the capability inventory extractor and coverage audit.

These cover the failure modes the gates exist to catch. They are deliberately
permanent tests rather than manual probes: a gate that stops detecting its own
failure mode is worse than no gate, because it reports a false assurance.
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    """!
    @brief Load one tooling module by path.
    @param[in] name Module file stem under tests/tooling.
    @return Imported module object.
    """
    path = REPO_ROOT / "tests" / "tooling" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GEN = load("generate_capability_inventory")
AUDIT = load("audit_capability_coverage")


MINIMUM_PYTHON = (3, 8)


def test_tooling_uses_no_syntax_newer_than_minimum_python():
    """!
    @brief The tooling must be *parseable* by the minimum supported Python.

    This checks syntax compatibility only. It compiles each tooling module against the
    minimum feature version, which catches syntax the declared floor cannot parse -
    the real risk when the local interpreter and CI differ. It does NOT establish
    runtime compatibility: a module can parse under 3.8 and still call a standard
    library API that only exists later. Full runtime compatibility is established by
    actually running the suite on the minimum interpreter.
    @return None.
    """
    import ast

    for name in ("generate_capability_inventory", "audit_capability_coverage"):
        source = (REPO_ROOT / "tests" / "tooling" / f"{name}.py").read_text(encoding="utf-8")
        ast.parse(source, feature_version=MINIMUM_PYTHON)


def test_tooling_runs_on_this_interpreter():
    """!
    @brief Record which interpreter actually exercised the tooling in this run.
    @return None.
    """
    assert sys.version_info[:2] >= MINIMUM_PYTHON, (
        f"tooling requires Python {MINIMUM_PYTHON}; running {sys.version_info[:2]}"
    )


# --------------------------------------------------------------------------
# Extractor
# --------------------------------------------------------------------------

def write_c(tmp_path: Path, body: str) -> Path:
    """!
    @brief Write a throwaway C source file for extractor tests.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] body File contents.
    @return Path to the written file.
    """
    path = tmp_path / "sample.c"
    path.write_text(body, encoding="utf-8")
    return path


def test_function_body_is_scoped_to_one_function(tmp_path, monkeypatch):
    """!
    @brief Extraction must not leak across function boundaries.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    write_c(tmp_path, """
PetscErrorCode Wanted(const char* str) {
    if (strcasecmp(str, "alpha") == 0) *handler_out = ENUM_ALPHA;
}
PetscErrorCode Other(const char* str) {
    if (strcasecmp(str, "beta") == 0) *handler_out = ENUM_BETA;
}
""")
    monkeypatch.setattr(GEN, "REPO_ROOT", tmp_path)
    found = GEN.c_string_map_values("sample.c", "Wanted")
    assert found == {"alpha": "ENUM_ALPHA"}, "extraction leaked into a neighbouring function"


def test_token_map_handles_shared_alias_arms(tmp_path, monkeypatch):
    """!
    @brief Two tokens sharing one assignment must both be captured.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    write_c(tmp_path, """
PetscErrorCode Setup(void) {
    if (strcmp(buf, "CANON") == 0 || strcmp(buf, "LEGACY") == 0) {
        ctx->field = ENUM_CANON;
    } else if (strcmp(buf, "OTHER") == 0) {
        ctx->field = ENUM_OTHER;
    }
}
""")
    monkeypatch.setattr(GEN, "REPO_ROOT", tmp_path)
    found = GEN.c_token_map_values("sample.c", "Setup", "buf", "field")
    assert found == {"CANON": "ENUM_CANON", "LEGACY": "ENUM_CANON", "OTHER": "ENUM_OTHER"}


def test_enum_extraction_does_not_span_enums(tmp_path, monkeypatch):
    """!
    @brief An enum regex must not swallow preceding enum blocks.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    write_c(tmp_path, """
typedef enum { OTHER_A = 0, OTHER_B = 1 } OtherType;
typedef enum { WANTED_A = 0, WANTED_B = 1 } WantedType;
""")
    monkeypatch.setattr(GEN, "REPO_ROOT", tmp_path)
    assert GEN.c_enum_values("sample.c", "WantedType") == {"WANTED_A", "WANTED_B"}


def test_unknown_parity_kind_is_rejected(monkeypatch):
    """!
    @brief An unrecognized parity source must fail loudly, never silently pass.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    family = {
        "id": "x",
        "title": "X",
        "selector": "s",
        "family_page": "p",
        "public_surface": {"kind": "python_dict", "module": "picurv_cli.core", "symbol": "BC_TYPE_MAP"},
        "parity_sources": [{"kind": "invented_kind"}],
    }
    with pytest.raises(RuntimeError, match="unknown parity source kind"):
        GEN.collect(family)


# --------------------------------------------------------------------------
# Audit
# --------------------------------------------------------------------------

def family_record(**overrides) -> dict:
    """!
    @brief Build a minimal inventory record for audit tests.
    @param[in] overrides Fields replacing the defaults.
    @return Inventory record.
    """
    record = {
        "id": "demo.family",
        "public_values": {"alpha": {"maps_to": "TOKEN_A"}},
        "parity": [
            {
                "source": {"kind": "c_token_map", "path": "src/x.c"},
                "values": ["TOKEN_A"],
                "mapping": {"TOKEN_A": "ENUM_A"},
            },
            {"source": {"kind": "c_dispatch", "path": "src/y.c"}, "values": ["ENUM_A"]},
        ],
    }
    record.update(overrides)
    return record


def test_parity_detects_token_the_parser_rejects():
    """!
    @brief A validator token absent from the C parser must be reported.
    @return None.
    """
    record = family_record(public_values={"alpha": {"maps_to": "MISSING_TOKEN"}})
    problems = AUDIT.check_parity(record)
    assert any("MISSING_TOKEN" in p for p in problems)


def test_parity_detects_enum_the_runtime_never_dispatches():
    """!
    @brief A token that parses but is never dispatched must be reported.
    @return None.
    """
    record = family_record()
    record["parity"][1]["values"] = []
    problems = AUDIT.check_parity(record)
    assert any("does not handle" in p for p in problems)


def test_parity_reports_unverified_source_kind():
    """!
    @brief A registered but unhandled parity kind must be reported, not ignored.
    @return None.
    """
    record = family_record()
    record["parity"].append({"source": {"kind": "brand_new_kind"}, "values": []})
    problems = AUDIT.check_parity(record)
    assert any("not verified by this audit" in p for p in problems)


def test_metadata_detects_stale_and_undeclared_values():
    """!
    @brief Metadata must match the values the sources expose, in both directions.
    @return None.
    """
    record = family_record(public_values={"alpha": {"maps_to": "T"}, "beta": {"maps_to": "T2"}})
    entry = {"value_metadata": {"alpha": {"canonical": True}, "removed_value": {"canonical": True}}}
    problems = AUDIT.check_metadata(record, entry)
    assert any("removed_value" in p and "no longer a public value" in p for p in problems)
    assert any("beta" in p and "no metadata entry" in p for p in problems)


def test_alias_must_name_its_target():
    """!
    @brief A non-canonical value without alias_of is a modelling error.
    @return None.
    """
    record = family_record()
    entry = {"value_metadata": {"alpha": {"canonical": False, "status": "deprecated"}}}
    problems = AUDIT.check_metadata(record, entry)
    assert any("names no alias_of or spelling_of target" in p for p in problems)


def coverage_page(tmp_path: Path, body: str) -> dict:
    """!
    @brief Write a family page and return the registry entry pointing at it.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] body Page contents.
    @return Registry entry for the audit.
    """
    pages = tmp_path / "docs" / "pages"
    pages.mkdir(parents=True, exist_ok=True)
    (pages / "FAM.md").write_text(body, encoding="utf-8")
    return {
        "family_page": "FAM",
        "entry_anchor_prefix": "cap_",
        "coverage_enforced": True,
        "value_metadata": {"alpha": {"canonical": True}},
    }


FULL_ENTRY = """
@anchor cap_alpha

**Identity.** x
**What it does.** x
**When to choose it.** x
**Parameters it owns.** x
**Interactions.** x
**Diagnostics.** x
**Evidence.** x
**Limitations.** x
"""


def test_coverage_accepts_a_complete_entry(tmp_path, monkeypatch):
    """!
    @brief A complete entry passes coverage.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, FULL_ENTRY)
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(family_record(), entry)
    assert blocking == []


def test_coverage_rejects_missing_required_field(tmp_path, monkeypatch):
    """!
    @brief An entry missing a contract part must be reported.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, FULL_ENTRY.replace("**Evidence.** x\n", ""))
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(family_record(), entry)
    assert any("missing the **Evidence** part" in p for p in blocking)


def test_coverage_ignores_anchor_inside_code_fence(tmp_path, monkeypatch):
    """!
    @brief An anchor shown inside a code example must not satisfy coverage.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, "```\n@anchor cap_alpha\n```\n")
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(family_record(), entry)
    assert any("no Tier-2 entry" in p for p in blocking)


def test_coverage_reports_stale_entry(tmp_path, monkeypatch):
    """!
    @brief An entry for a value that no longer exists must be reported.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, FULL_ENTRY + "\n@anchor cap_deleted_value\n")
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(family_record(), entry)
    assert any("stale entry" in p and "cap_deleted_value" in p for p in blocking)


def test_coverage_detects_slug_collision(tmp_path, monkeypatch):
    """!
    @brief Two values whose slugs collide would share one entry.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, FULL_ENTRY)
    entry["value_metadata"]["alpha-x"] = {"canonical": True}
    record = family_record(public_values={"alpha_x": {"maps_to": "T"}, "alpha-x": {"maps_to": "T"}})
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(record, entry)
    assert any("same anchor slug" in p for p in blocking)


def test_coverage_skips_latent_values(tmp_path, monkeypatch):
    """!
    @brief A latent value is not selectable and owes no entry.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, FULL_ENTRY)
    record = family_record(
        public_values={"alpha": {"maps_to": "T"}, "ghost": {"maps_to": "G", "reachability": "latent"}}
    )
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(record, entry)
    assert not any("ghost" in p for p in blocking)


def test_alias_entry_requires_only_the_stub_shape(tmp_path, monkeypatch):
    """!
    @brief A deprecated alias needs Identity/Status/Migration, not the full contract.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    body = "\n@anchor cap_alpha\n\n**Identity.** x\n**Status.** x\n**Migration.** x\n"
    entry = coverage_page(tmp_path, body)
    entry["value_metadata"] = {"alpha": {"canonical": False, "alias_of": "beta"}}
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, _ = AUDIT.check_coverage(family_record(), entry)
    assert blocking == []


def test_advisory_family_does_not_block(tmp_path, monkeypatch):
    """!
    @brief A family still awaiting backfill reports without failing the gate.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    entry = coverage_page(tmp_path, "no entries here\n")
    entry["coverage_enforced"] = False
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    blocking, advisory = AUDIT.check_coverage(family_record(), entry)
    assert blocking == [] and advisory


# --------------------------------------------------------------------------
# End to end against the real repository
# --------------------------------------------------------------------------

def test_repository_inventory_is_current_and_consistent():
    """!
    @brief The committed inventory must match the sources and pass every check.
    @return None.
    """
    inventory = json.loads((REPO_ROOT / "docs" / "generated" / "capability_inventory.json").read_text())
    registry = {
        f["id"]: f
        for f in json.loads((REPO_ROOT / "tests" / "tooling" / "capability_families.json").read_text())["families"]
    }
    assert AUDIT.check_generated_current() == []
    for family in inventory:
        entry = registry[family["id"]]
        assert AUDIT.check_parity(family) == []
        assert AUDIT.check_metadata(family, entry) == []
        blocking, _ = AUDIT.check_coverage(family, entry, list(registry.values()))
        assert blocking == []


# --------------------------------------------------------------------------
# Invariant contract registry
# --------------------------------------------------------------------------

CONTRACTS = load("audit_contracts")


def _record(**overrides) -> dict:
    """!
    @brief Build a minimal valid contract record for registry validation tests.
    @param[in] overrides Fields replacing the defaults.
    @return Contract record.
    """
    record = {
        "id": "x", "title": "X", "kind": "units", "status": "tracked",
        "authoritative_sources": ["Makefile"], "checker": None,
        "canonical_documentation": [], "dependent_pages": [], "enforcement": "report-only",
    }
    record.update(overrides)
    return record


def test_contract_registry_records_are_valid():
    """!
    @brief Every registered contract must point at sources and pages that exist.
    @return None.
    """
    contracts = CONTRACTS.load_registry()["contracts"]
    assert CONTRACTS.validate_records(contracts, CONTRACTS.load_registry()) == []


def test_contract_registry_detects_missing_source():
    """!
    @brief A contract naming a nonexistent source must be reported.
    @return None.
    """
    bad = [_record(authoritative_sources=["src/does_not_exist.c"])]
    problems = CONTRACTS.validate_records(bad, CONTRACTS.load_registry())
    assert any("does not exist" in p for p in problems)


def test_contract_registry_detects_missing_page():
    """!
    @brief A contract naming a nonexistent page must be reported.
    @return None.
    """
    bad = [_record(canonical_documentation=["99_No_Such_Page"])]
    problems = CONTRACTS.validate_records(bad, CONTRACTS.load_registry())
    assert any("99_No_Such_Page" in p for p in problems)


def test_enforced_contract_requires_a_checker():
    """!
    @brief A contract cannot claim enforced status without a checker.
    @return None.
    """
    bad = [_record(status="enforced", enforcement="blocking", checker=None)]
    problems = CONTRACTS.validate_records(bad, CONTRACTS.load_registry())
    assert any("no checker is registered" in p for p in problems)


def test_contract_ids_are_unique():
    """!
    @brief Duplicate contract ids must be reported.
    @return None.
    """
    dup = _record()
    problems = CONTRACTS.validate_records([dup, dict(dup)], CONTRACTS.load_registry())
    assert any("duplicate" in p for p in problems)


def test_artifact_topology_snapshot_is_current():
    """!
    @brief The recorded run topology must match what the CLI planner produces.
    @return None.
    """
    snapshot = REPO_ROOT / "docs" / "generated" / "artifact_topology_snapshot.json"
    assert snapshot.is_file(), "topology snapshot missing; run 'make docs-topology'"
    recorded = json.loads(snapshot.read_text(encoding="utf-8"))
    assert recorded["scenarios"], "snapshot records no scenarios"
    artifacts = recorded["scenarios"][0]["artifacts"]
    # Normalization must have removed every unstable component.
    assert all("<run_id>" in a or a == "<workspace>/runs" or "<workspace>" in a for a in artifacts)
    assert not any(re.search(r"\d{8}-\d{6}", a) for a in artifacts), "timestamps not normalized"


def test_topology_contract_declares_destructive_scope():
    """!
    @brief Every artifact record must state what destructive operations can reach it.
    @return None.
    """
    contract = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "artifact_topology.json").read_text(encoding="utf-8")
    )
    for artifact in contract["artifacts"]:
        assert "destructive_scope" in artifact, f"{artifact['id']} declares no destructive scope"
        assert artifact["isolation_model"] in contract["isolation_models"]


def test_registry_rejects_misspelled_status():
    """!
    @brief A status typo must fail loudly rather than dropping the contract silently.
    @return None.
    """
    registry = CONTRACTS.load_registry()
    bad = [dict(registry["contracts"][0], status="enfroced")]
    problems = CONTRACTS.validate_records(bad, registry)
    assert any("closed vocabulary" in p for p in problems)


def test_registry_rejects_status_enforcement_mismatch():
    """!
    @brief Enforcement must agree with status in both directions.
    @return None.
    """
    registry = CONTRACTS.load_registry()
    bad = [dict(registry["contracts"][0], status="tracked", enforcement="blocking")]
    problems = CONTRACTS.validate_records(bad, registry)
    assert any("requires status 'enforced'" in p for p in problems)


def test_registry_rejects_missing_required_field():
    """!
    @brief Every required field must be present.
    @return None.
    """
    registry = CONTRACTS.load_registry()
    bad = dict(registry["contracts"][0])
    bad.pop("title")
    problems = CONTRACTS.validate_records([bad], registry)
    assert any("missing required field 'title'" in p for p in problems)


def test_registry_rejects_bad_checker_path():
    """!
    @brief A checker naming a nonexistent script must be reported.
    @return None.
    """
    registry = CONTRACTS.load_registry()
    bad = [dict(registry["contracts"][0], checker=["tests/tooling/not_a_checker.py"])]
    problems = CONTRACTS.validate_records(bad, registry)
    assert any("does not exist" in p for p in problems)


def test_built_docs_prerequisite_is_modelled():
    """!
    @brief The site checker must declare that it needs the built HTML.

    Without this the runner would execute it on a clean checkout and either fail
    confusingly or pass on a stale build.
    @return None.
    """
    registry = CONTRACTS.load_registry()
    site = next(c for c in registry["contracts"] if c["id"] == "docs.published_site")
    assert site.get("requires_built_docs") is True


def test_every_planned_artifact_maps_to_an_identity():
    """!
    @brief No planned artifact may be left without a declared logical identity.
    @return None.
    """
    snapshot = json.loads(
        (REPO_ROOT / "docs" / "generated" / "artifact_topology_snapshot.json").read_text(encoding="utf-8")
    )
    for scenario in snapshot["scenarios"]:
        assert scenario["identities"]["unmapped"] == [], (
            f"{scenario['scenario']} has unmapped artifacts: {scenario['identities']['unmapped']}"
        )


def test_topology_snapshot_covers_custom_directories():
    """!
    @brief A custom output/log directory must appear in the plan, not be hidden by defaults.
    @return None.
    """
    snapshot = json.loads(
        (REPO_ROOT / "docs" / "generated" / "artifact_topology_snapshot.json").read_text(encoding="utf-8")
    )
    names = {s["scenario"] for s in snapshot["scenarios"]}
    assert "custom_directories" in names
    custom = next(s for s in snapshot["scenarios"] if s["scenario"] == "custom_directories")
    assert any("custom_logs" in a for a in custom["artifacts"]), "configured log dir not visible in plan"
    assert any("custom_out" in a for a in custom["artifacts"]), "configured output dir not visible in plan"


def test_isolation_is_not_claimed_as_enforced():
    """!
    @brief The topology contract must not claim an enforced isolation model.
    @return None.
    """
    contract = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "artifact_topology.json").read_text(encoding="utf-8")
    )
    assert "current_isolation_model" not in contract
    assert contract["default_layout"]
    # The destructive-path defect is closed. It lives under resolved history rather
    # than under a present-tense "known defect" key, because a reader scanning the
    # contract should not have to read a paragraph to learn that it is fixed.
    assert "known_safety_defect" not in contract, (
        "a resolved defect must not sit under a present-tense key"
    )
    history = contract["resolved_safety_history"]["entries"]
    assert history, "the enforcement is only intelligible next to what motivated it"
    for entry in history:
        assert entry["state"] == "resolved"
        assert entry["resolution"]["residual_risk"], (
            "a resolution must state what it does not cover"
        )
        assert "not_changed" not in entry["resolution"], (
            "a field that became false is removed, not corrected in place"
        )
    snapshot = json.loads(
        (REPO_ROOT / "docs" / "generated" / "artifact_topology_snapshot.json").read_text(encoding="utf-8")
    )
    assert snapshot["isolation_enforced"] is False


# --------------------------------------------------------------------------
# Scaffolding and lifecycle enforcement
# --------------------------------------------------------------------------

SCAFFOLD = load("scaffold_documentation")


def test_scaffolded_capability_entry_satisfies_the_contract():
    """!
    @brief A generated entry must already carry every required contract part.

    A scaffold that produces something the audit rejects is worse than none: it teaches
    the wrong shape.
    @return None.
    """
    body = SCAFFOLD.capability_entry("boundary.handler", "example_value")
    for field in AUDIT.CANONICAL_FIELDS:
        assert re.search(rf"\*\*{re.escape(field)}", body), f"scaffold omits **{field}**"


def test_scaffolded_alias_stub_satisfies_the_alias_contract():
    """!
    @brief A generated alias stub owes three parts, not eight.
    @return None.
    """
    body = SCAFFOLD.alias_stub("momentum.solver", "Old Name", "Explicit RK4")
    for field in AUDIT.ALIAS_FIELDS:
        assert re.search(rf"\*\*{re.escape(field)}", body)
    assert "**Evidence.**" not in body, "an alias stub must not claim the full contract"


def test_scaffolded_capability_anchor_matches_audit_expectation():
    """!
    @brief The scaffold's anchor must be the one coverage looks for.
    @return None.
    """
    body = SCAFFOLD.capability_entry("boundary.handler", "Example Value")
    expected = "p44_cap_" + AUDIT.slug("Example Value")
    assert f"@anchor {expected}" in body


def test_scaffolded_page_declares_a_valid_page_type():
    """!
    @brief Every page scaffold must carry a three-field pagemeta and a declared type.
    @return None.
    """
    for kind in SCAFFOLD.PAGE_TYPES:
        text = SCAFFOLD.page(kind, "Sample Page", 99)
        match = re.search(r"@pagemeta\{([^}]*)\}", text)
        assert match, f"{kind} scaffold has no pagemeta"
        assert len(match.group(1).split(",")) == 3, f"{kind} pagemeta must have exactly 3 fields"
        assert SCAFFOLD.PAGE_TYPES[kind] in match.group(1)


def test_scaffolded_page_has_a_page_directive_and_anchor():
    """!
    @brief A page scaffold must be adoptable into the navigation hierarchy.
    @return None.
    """
    text = SCAFFOLD.page("reference", "Sample Page", 99)
    assert re.search(r"^@page\s+99_Sample_Page", text, re.M)
    assert "@anchor _Sample_Page" in text


def test_scaffold_rejects_unknown_family():
    """!
    @brief Scaffolding against a nonexistent family must fail loudly.
    @return None.
    """
    with pytest.raises(SystemExit):
        SCAFFOLD.capability_entry("no.such.family", "x")


def test_subsystem_charter_is_scoped_to_planned_status():
    """!
    @brief A new subsystem charter must ask only what `planned` owes.
    @return None.
    """
    text = SCAFFOLD.subsystem_charter("ibm", "Immersed Boundaries")
    assert "planned" in text
    assert "Explicitly out of scope" in text
    assert "Concern modules" in text


def test_lifecycle_requirements_cover_each_status_that_owes_something():
    """!
    @brief The lifecycle table must demand evidence, limitations, and migration.
    @return None.
    """
    requirements = AUDIT.LIFECYCLE_REQUIREMENTS
    assert requirements["supported"] == ("Evidence",)
    assert "Limitations" in requirements["experimental"]
    assert "Limitations" in requirements["known-defective"]
    assert "Migration" in requirements["deprecated"]


def test_repository_satisfies_lifecycle_requirements():
    """!
    @brief Every declared value must already meet what its status owes.
    @return None.
    """
    inventory = json.loads(
        (REPO_ROOT / "docs" / "generated" / "capability_inventory.json").read_text(encoding="utf-8")
    )
    registry = {
        f["id"]: f
        for f in json.loads(
            (REPO_ROOT / "tests" / "tooling" / "capability_families.json").read_text(encoding="utf-8")
        )["families"]
    }
    for family in inventory:
        assert AUDIT.check_lifecycle_requirements(family, registry[family["id"]]) == []


# --------------------------------------------------------------------------
# Lifecycle status must fail closed
# --------------------------------------------------------------------------

def _family_pair(status=None, extra=None, value="Newton Krylov"):
    """!
    @brief Build an (inventory record, registry entry) pair for status tests.
    @param[in] status Status to declare, or None to omit the key entirely.
    @param[in] extra Additional metadata fields.
    @param[in] value Value name to mutate.
    @return Tuple of (inventory record, registry entry).
    """
    inventory = json.loads(
        (REPO_ROOT / "docs" / "generated" / "capability_inventory.json").read_text(encoding="utf-8")
    )
    registry = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "capability_families.json").read_text(encoding="utf-8")
    )
    record = next(f for f in inventory if f["id"] == "momentum.solver")
    entry = next(f for f in registry["families"] if f["id"] == "momentum.solver")
    entry = json.loads(json.dumps(entry))
    meta = entry["value_metadata"][value]
    meta.pop("status", None)
    if status is not None:
        meta["status"] = status
    if extra:
        meta.update(extra)
    return record, entry


def test_missing_status_is_rejected():
    """!
    @brief A canonical value with no status must not pass silently.
    @return None.
    """
    record, entry = _family_pair(status=None)
    problems = AUDIT.check_metadata(record, entry)
    assert any("declares no status" in p for p in problems)


@pytest.mark.parametrize("typo", ["suported", "SUPPORTED", "stable", "", "deprecated ", None])
def test_status_typos_are_rejected(typo):
    """!
    @brief Only the closed vocabulary is accepted; a typo must not lose the obligation.
    @param[in] typo Invalid status value.
    @return None.
    """
    record, entry = _family_pair(status=typo)
    assert AUDIT.check_metadata(record, entry) != []


@pytest.mark.parametrize("status", AUDIT.VALID_STATUSES)
def test_every_valid_status_is_accepted_by_the_vocabulary(status):
    """!
    @brief Each declared status must be recognized, whatever else it then requires.
    @param[in] status Status from the closed vocabulary.
    @return None.
    """
    record, entry = _family_pair(status=status)
    problems = AUDIT.check_metadata(record, entry)
    assert not any("closed vocabulary" in p for p in problems)


@pytest.mark.parametrize("status", AUDIT.NON_SELECTABLE_STATUSES)
def test_non_selectable_status_on_a_reachable_value_is_rejected(status):
    """!
    @brief A planned, internal, or removed capability must not be publicly selectable.
    @param[in] status Non-selectable status.
    @return None.
    """
    record, entry = _family_pair(status=status)
    problems = AUDIT.check_metadata(record, entry)
    assert any("must not be reachable" in p for p in problems)


def test_alias_must_declare_deprecated_status():
    """!
    @brief An alias is by definition deprecated; any other status is incoherent.
    @return None.
    """
    record, entry = _family_pair(status="supported", extra={"alias_of": "Explicit RK4"})
    problems = AUDIT.check_metadata(record, entry)
    assert any("an alias is by definition deprecated" in p for p in problems)


def test_spelling_must_not_declare_its_own_status():
    """!
    @brief A spelling inherits status; declaring one would let it diverge.
    @return None.
    """
    record, entry = _family_pair(status="supported", extra={"spelling_of": "Explicit RK4"})
    problems = AUDIT.check_metadata(record, entry)
    assert any("must not\ndeclare its own status" in p or "must not declare its own status" in p
               for p in problems)


def test_alias_keeps_its_own_status_in_the_inventory():
    """!
    @brief A deprecated alias must render as deprecated, not inherit the canonical status.
    @return None.
    """
    inventory = json.loads(
        (REPO_ROOT / "docs" / "generated" / "capability_inventory.json").read_text(encoding="utf-8")
    )
    solvers = next(f for f in inventory if f["id"] == "momentum.solver")
    alias = solvers["public_values"]["Dual Time Picard RK4"]
    assert alias["status"] == "deprecated"
    assert alias["alias_of"] == "Dual Time Picard Jameson RK"


def test_latent_is_reachability_not_a_lifecycle_status():
    """!
    @brief `latent` must never appear as a declared lifecycle status.
    @return None.
    """
    registry = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "capability_families.json").read_text(encoding="utf-8")
    )
    assert "latent" not in AUDIT.VALID_STATUSES
    for family in registry["families"]:
        for name, meta in family.get("value_metadata", {}).items():
            assert meta.get("status") != "latent", f"{family['id']}::{name} uses latent as a status"


# --------------------------------------------------------------------------
# Scaffolds must not certify while incomplete
# --------------------------------------------------------------------------

def test_every_scaffold_carries_the_incomplete_marker():
    """!
    @brief Structural correctness and completion are separate conditions.
    @return None.
    """
    marker = SCAFFOLD.INCOMPLETE_MARKER
    assert marker in SCAFFOLD.capability_entry("boundary.handler", "x")
    assert marker in SCAFFOLD.alias_stub("momentum.solver", "Old", "Explicit RK4")
    assert marker in SCAFFOLD.subsystem_charter("ibm", "Immersed Boundaries")
    for kind in SCAFFOLD.PAGE_TYPES:
        assert marker in SCAFFOLD.page(kind, "Sample", 99)


def test_incomplete_marker_is_a_forbidden_signature():
    """!
    @brief The publication gate must reject an untouched scaffold.
    @return None.
    """
    contract = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "generic_expansion_contract.json").read_text(encoding="utf-8")
    )
    assert SCAFFOLD.INCOMPLETE_MARKER in contract["forbidden_markers"]


# --------------------------------------------------------------------------
# Identifier sanitization
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "title,expected",
    [("A/B: Test?", "A_B_Test"), ("Running On A Cluster", "Running_On_A_Cluster"),
     ("123 Start", "p_123_Start"), ("a---b", "a_b")],
)
def test_identifier_slug_is_safe(title, expected):
    """!
    @brief Free-form titles must produce legal identifiers.
    @param[in] title Input title.
    @param[in] expected Expected identifier.
    @return None.
    """
    assert SCAFFOLD.identifier_slug(title) == expected


@pytest.mark.parametrize("title", ["   ", "///", "???"])
def test_empty_slug_is_rejected(title):
    """!
    @brief A title with nothing usable must fail rather than emit a broken identifier.
    @param[in] title Input title.
    @return None.
    """
    with pytest.raises(SystemExit):
        SCAFFOLD.identifier_slug(title)


def test_scaffolded_page_keeps_the_display_title():
    """!
    @brief Sanitization applies to the identifier, not to what the reader sees.
    @return None.
    """
    text = SCAFFOLD.page("reference", "A/B: Test?", 72)
    assert "@page 72_A_B_Test A/B: Test?" in text
    assert "@anchor _A_B_Test" in text


def test_alias_scaffold_validates_family_and_canonical_target():
    """!
    @brief Aliases get the same friendly validation as capabilities.
    @return None.
    """
    with pytest.raises(SystemExit, match="Unknown family"):
        SCAFFOLD.alias_stub("no.such.family", "x", "y")
    with pytest.raises(SystemExit, match="not a declared value"):
        SCAFFOLD.alias_stub("momentum.solver", "x", "No Such Canonical")


# --------------------------------------------------------------------------
# Page templates carry their type's contract
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "kind,required",
    [("tutorial", ["Outcome", "Prerequisites", "Confirm It Worked"]),
     ("how-to", ["Goal", "Procedure", "Verify"]),
     ("reference", ["Scope", "Constraints"]),
     ("explanation", ["Concepts", "Alternatives Considered", "Limitations"]),
     ("hub", ["Routes"])],
)
def test_page_templates_are_type_specific(kind, required):
    """!
    @brief Each page type must scaffold the sections it actually owes.
    @param[in] kind Page type key.
    @param[in] required Section headings that must be present.
    @return None.
    """
    text = SCAFFOLD.page(kind, "Sample Page", 99)
    for heading in required:
        assert heading in text, f"{kind} scaffold omits '{heading}'"


def test_page_templates_differ_between_types():
    """!
    @brief Types must not collapse into one generic outline.
    @return None.
    """
    shapes = {
        kind: tuple(h for h, _ in SCAFFOLD.PAGE_SECTIONS[kind]) for kind in SCAFFOLD.PAGE_TYPES
    }
    assert len(set(shapes.values())) == len(shapes), "two page types share an identical shape"


# --------------------------------------------------------------------------
# Markdown fragment-link extraction
# --------------------------------------------------------------------------

SITE = load("audit_docs_site")


@pytest.mark.parametrize(
    "text,expected",
    [("[a](#simple)", [(None, "simple")]),
     ("[a](#hyphenated-anchor)", [(None, "hyphenated-anchor")]),
     ('[a](#titled "some title")', [(None, "titled")]),
     ("[a](other.md#cross)", [("other.md", "cross")]),
     ("[a](#dotted.name)", [(None, "dotted.name")])],
)
def test_fragment_links_are_extracted(text, expected):
    """!
    @brief Every legal fragment form must be seen, not only word characters.
    @param[in] text Markdown source.
    @param[in] expected Expected extraction.
    @return None.
    """
    assert SITE.markdown_fragment_links(text) == expected


@pytest.mark.parametrize(
    "text",
    ["```\n[a](#in_fence)\n```",
     "~~~\n[a](#in_tilde_fence)\n~~~",
     "<!-- [a](#in_comment) -->",
     "`[a](#in_inline_code)`"],
)
def test_fragment_links_inside_non_prose_are_ignored(text):
    """!
    @brief A link shown inside an example is documentation of a link, not a link.
    @param[in] text Markdown source.
    @return None.
    """
    assert SITE.markdown_fragment_links(text) == []


def test_strip_non_prose_preserves_line_count():
    """!
    @brief Blanking rather than deleting keeps line numbers meaningful.
    @return None.
    """
    text = "one\n```\ntwo\nthree\n```\nfour\n"
    assert SITE.strip_non_prose(text).count("\n") == text.count("\n")


def test_repository_fragment_links_all_resolve():
    """!
    @brief Every hand-written fragment link in the docs must resolve.
    @return None.
    """
    html_dir = REPO_ROOT / "docs_build" / "html"
    if not (html_dir / "index.html").is_file():
        pytest.skip("site not built; run 'make build-docs'")
    assert SITE.check_page_cross_references(html_dir) == []


# --------------------------------------------------------------------------
# Fail-open closures found in review
# --------------------------------------------------------------------------

def test_absent_metadata_block_reports_every_value():
    """!
    @brief A family exposing values with no metadata block must not pass silently.
    @return None.
    """
    inventory = json.loads(
        (REPO_ROOT / "docs" / "generated" / "capability_inventory.json").read_text(encoding="utf-8")
    )
    record = next(f for f in inventory if f["id"] == "boundary.type")
    problems = AUDIT.check_metadata(record, {})
    assert len(problems) == len(record["public_values"])
    assert all("has no metadata entry" in p for p in problems)


def test_reachability_never_overwrites_lifecycle_status():
    """!
    @brief The two axes must both survive into the generated inventory.

    Reachability previously overwrote `status` with `public`/`latent`, which are not
    lifecycle values, silently disabling lifecycle enforcement for the family.
    @return None.
    """
    inventory = json.loads(
        (REPO_ROOT / "docs" / "generated" / "capability_inventory.json").read_text(encoding="utf-8")
    )
    types = next(f for f in inventory if f["id"] == "boundary.type")
    for name, spec in types["public_values"].items():
        assert spec["status"] in AUDIT.VALID_STATUSES, f"{name} status is not a lifecycle value"
        assert spec["reachability"] in ("reachable", "latent")
    assert types["public_values"]["symmetry"]["status"] == "internal"
    assert types["public_values"]["symmetry"]["reachability"] == "latent"
    assert types["public_values"]["inlet"]["reachability"] == "reachable"


def test_experimental_value_in_a_reachable_family_still_owes_limitations(tmp_path, monkeypatch):
    """!
    @brief Lifecycle enforcement must survive in families that compute reachability.

    Reachability previously overwrote `status`, so the audit read "reachable" instead of
    the declared lifecycle and found no obligation at all.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    pages = tmp_path / "docs" / "pages"
    pages.mkdir(parents=True)
    # An entry with every part except Limitations.
    (pages / "FAM.md").write_text(
        "\n@anchor cap_alpha\n\n"
        "**Identity.** x\n**What it does.** x\n**When to choose it.** x\n"
        "**Parameters it owns.** x\n**Interactions.** x\n**Diagnostics.** x\n**Evidence.** x\n",
        encoding="utf-8",
    )
    record = {
        "id": "demo.family",
        "public_values": {
            "alpha": {"maps_to": "T", "status": "experimental", "reachability": "reachable"}
        },
        "parity": [],
    }
    entry = {
        "family_page": "FAM",
        "entry_anchor_prefix": "cap_",
        "value_metadata": {"alpha": {"status": "experimental", "canonical": True, "evidence": {}}},
    }
    monkeypatch.setattr(AUDIT, "REPO_ROOT", tmp_path)
    problems = AUDIT.check_lifecycle_requirements(record, entry)
    assert any("experimental" in p and "Limitations" in p for p in problems)


def test_scaffold_cli_prints_a_safe_filename():
    """!
    @brief The suggested filename must use the sanitized identifier, not the raw title.
    @return None.
    """
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tests" / "tooling" / "scaffold_documentation.py"),
         "page", "--type", "reference", "--title", "A/B: Test?"],
        capture_output=True, text=True, check=False, cwd=str(REPO_ROOT),
    )
    assert result.returncode == 0
    hint = result.stderr
    assert "72_A_B_Test.md" in hint or re.search(r"docs/pages/\d+_A_B_Test\.md", hint)
    assert "/B:" not in hint and "?" not in hint, "raw title leaked into the filename hint"


def test_fragment_gate_covers_markdown_outside_docs_pages():
    """!
    @brief README and guide.md files must be inside fragment validation.
    @return None.
    """
    scanned = {p.name for p in SITE.tracked_markdown()}
    assert "README.md" in scanned
    assert "guide.md" in scanned


def test_fragment_gate_does_not_conflate_duplicate_basenames():
    """!
    @brief Many files share the name guide.md; targets must resolve relative to the source.
    @return None.
    """
    guides = [p for p in SITE.tracked_markdown() if p.name == "guide.md"]
    assert len(guides) > 1, "expected several guide.md files in this repository"
    assert len({p.resolve() for p in guides}) == len(guides)


def test_markdown_heading_anchors_are_available_for_unrendered_files():
    """!
    @brief A plain Markdown file offers heading anchors, not nothing.
    @return None.
    """
    assert SITE.heading_anchor("## Where work goes") == "where-work-goes"
    assert SITE.heading_anchor("### A `code` heading!") == "a-code-heading"
    anchors = SITE.markdown_anchors(REPO_ROOT / "CLAUDE.md")
    assert anchors, "CLAUDE.md should expose heading anchors"


def test_repository_fragment_links_resolve_repo_wide():
    """!
    @brief Every fragment link in every tracked Markdown file must resolve.
    @return None.
    """
    html_dir = REPO_ROOT / "docs_build" / "html"
    if not (html_dir / "index.html").is_file():
        pytest.skip("site not built; run 'make build-docs'")
    assert SITE.check_page_cross_references(html_dir) == []


def test_solver_output_is_recorded_as_containment_validated():
    """!
    @brief The lifecycle record must not claim configured solver output is unvalidated.

    @details It said "configured overrides are not containment-validated" long after
             `io.directories.output` joined the shared Python rule set. A reader
             checking whether their configured output directory is checked would have
             been told no.
    @return None.
    """
    contract = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "artifact_topology.json").read_text(encoding="utf-8")
    )
    records = {a["id"]: a for a in contract["artifacts"]}
    output = records["run.solver_output"]
    assert "not containment-validated" not in output["containment"]
    assert "validated by picurv" in output["containment"]
    # And it must not overclaim in the other direction: the independent C guard stands
    # in front of PetscRMTree, which only ever deletes the log directory.
    assert "no independent C guard" in output["containment"]
    note = contract["default_layout_note"]
    assert "BOTH" in note and "only the log directory" in note


def test_post_output_is_recorded_as_configurable():
    """!
    @brief Post-processing output is configured, not a fixed `visualization/<recipe>`.

    @details `post.yml -> io.output_directory` defaults to `viz` and is resolved
             directly against the run directory, and shipped recipes use both a flat
             and a nested layout. Recording it as fixed made the contract describe one
             example rather than the rule.
    @return None.
    """
    contract = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "artifact_topology.json").read_text(encoding="utf-8")
    )
    visualization = {a["id"]: a for a in contract["artifacts"]}["run.visualization"]
    assert visualization["configurable"] is True
    assert "configured post output dir" in visualization["path_rule"]
    assert "NOT containment-validated" in visualization["containment"]


def test_the_snapshot_exercises_more_than_one_post_layout():
    """!
    @brief The extractor must not let the snapshot collapse to one post layout.
    @return None.
    """
    snapshot = json.loads(
        (REPO_ROOT / "docs" / "generated" / "artifact_topology_snapshot.json").read_text(
            encoding="utf-8")
    )
    layouts = set()
    for scenario in snapshot["scenarios"]:
        for token in scenario["identities"]["mapped"].get("run.visualization", []):
            layouts.add(token)
    assert len(layouts) >= 2, (
        "only one post output layout is exercised; the snapshot would record the "
        "quickstart's nested path as though it were the fixed shape"
    )
    assert any(token.endswith("/viz") for token in layouts)
    assert not any(scenario["identities"]["unmapped"] for scenario in snapshot["scenarios"])
