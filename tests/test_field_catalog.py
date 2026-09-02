"""!
@file test_field_catalog.py
@brief Pytest coverage for the field-catalog identity and layout audit.
"""

import importlib.machinery
import importlib.util
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
AUDIT = REPO_ROOT / "tests" / "tooling" / "audit_field_catalog.py"
PAGE = REPO_ROOT / "docs" / "pages" / "56_Field_Identity_and_Layout_Catalog.md"


def load_audit_module():
    """!
    @brief Load the field-catalog audit as an importable module.
    @return Value returned by `load_audit_module()`.
    """
    loader = importlib.machinery.SourceFileLoader("field_catalog_audit", str(AUDIT))
    spec = importlib.util.spec_from_loader("field_catalog_audit", loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def audit_against(tmp_path, page_text, module):
    """!
    @brief Run the audit against a modified copy of the catalog page.
    @param[in] tmp_path Pytest temporary directory supplied to the function.
    @param[in] page_text Replacement page content supplied to the function.
    @param[in] module Loaded audit module supplied to the function.
    @return Process status code the audit returned.
    """
    copy = tmp_path / PAGE.name
    copy.write_text(page_text, encoding="utf-8")
    module.PAGE = copy
    return module.main()


def test_the_repository_passes_the_audit():
    """!
    @brief Verify the committed page and compiled catalogs agree.
    """
    assert load_audit_module().main() == 0


def test_every_compiled_identity_is_documented():
    """!
    @brief Verify the audit reads a non-trivial catalog rather than an empty one.
    """
    module = load_audit_module()
    eulerian, particle = module.compiled_eulerian(), module.compiled_particle()
    assert len(eulerian) > 20 and len(particle) > 5
    assert set(eulerian.values()) <= module.compiled_layouts()


def test_a_field_missing_from_the_inventory_is_rejected(tmp_path):
    """!
    @brief Verify a compiled field absent from section 6 fails the audit.
    @param[in] tmp_path Pytest temporary directory supplied to the function.
    """
    module = load_audit_module()
    text = PAGE.read_text(encoding="utf-8")
    name = sorted(module.compiled_eulerian())[0]
    assert audit_against(tmp_path, text.replace(f"`{name}`, ", "", 1), module) == 1


def test_a_field_documented_under_the_wrong_layout_is_rejected(tmp_path):
    """!
    @brief Verify a layout-group mismatch fails the audit.
    @param[in] tmp_path Pytest temporary directory supplied to the function.
    """
    module = load_audit_module()
    text = PAGE.read_text(encoding="utf-8")
    moved = text.replace("- component-staggered: `Ucont`, ", "- component-staggered: ", 1)
    moved = moved.replace("- I-face family: `Csi`,", "- I-face family: `Ucont`, `Csi`,", 1)
    assert moved != text
    assert audit_against(tmp_path, moved, module) == 1


def test_an_undocumented_layout_value_is_rejected(tmp_path):
    """!
    @brief Verify a layout defined in the header but untabulated fails the audit.
    @param[in] tmp_path Pytest temporary directory supplied to the function.
    """
    module = load_audit_module()
    text = PAGE.read_text(encoding="utf-8")
    stripped = "\n".join(
        line for line in text.splitlines()
        if not line.startswith("| `FIELD_LAYOUT_J_FACE`")
    )
    assert stripped != text
    assert audit_against(tmp_path, stripped, module) == 1


def test_a_particle_field_missing_from_the_inventory_is_rejected(tmp_path):
    """!
    @brief Verify a compiled particle field absent from section 6 fails the audit.
    @param[in] tmp_path Pytest temporary directory supplied to the function.
    """
    module = load_audit_module()
    text = PAGE.read_text(encoding="utf-8")
    assert audit_against(tmp_path, text.replace("`position`, ", "", 1), module) == 1


def test_an_unrecognized_layout_group_is_rejected(tmp_path):
    """!
    @brief Verify a renamed inventory group fails rather than silently dropping fields.
    @param[in] tmp_path Pytest temporary directory supplied to the function.
    """
    module = load_audit_module()
    text = PAGE.read_text(encoding="utf-8")
    renamed = text.replace("- I-face family:", "- I-face group:", 1)
    assert renamed != text
    assert audit_against(tmp_path, renamed, module) == 1


def test_the_contract_registry_points_at_this_audit():
    """!
    @brief Verify the invariant contract registers this script as its blocking checker.
    """
    registry = json.loads(
        (REPO_ROOT / "tests" / "tooling" / "contract_registry.json").read_text(encoding="utf-8")
    )
    record = next(c for c in registry["contracts"] if c["id"] == "field.identity_and_layout")
    assert record["checker"] == ["tests/tooling/audit_field_catalog.py"]
    assert record["status"] == "enforced"
    assert record["enforcement"] == "blocking"
