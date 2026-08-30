"""The census classifications are claims, and these check them.

A classification such as "documented as a parameter of X" is only worth writing if X
exists and actually documents it. These tests hold the census to what it says, and
fail when a surface outgrows the classification it was given.
"""

import ast
import json
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLING = REPO_ROOT / "tests" / "tooling"
PAGES = REPO_ROOT / "docs" / "pages"

sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(scope="module")
def census():
    """!
    @brief The committed census document.
    @return Parsed census.
    """
    return json.loads((TOOLING / "family_census.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def families():
    """!
    @brief Registered capability families, keyed by id.
    @return Mapping of family id to record.
    """
    document = json.loads((TOOLING / "capability_families.json").read_text(encoding="utf-8"))
    return {f["id"]: f for f in document["families"]}


def _constant(module_symbol):
    """!
    @brief Read a module-level constant named in `module::SYMBOL` form.
    @param[in] module_symbol Census key.
    @return The literal value, or None when it is not a module constant.
    """
    module, _, symbol = module_symbol.partition("::")
    path = REPO_ROOT / Path(*module.split(".")).with_suffix(".py")
    if not path.is_file():
        return None
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == symbol for t in node.targets):
            try:
                return ast.literal_eval(node.value)
            except (ValueError, TypeError):
                return None
    return None


def test_no_family_is_public_pending(census):
    """!
    @brief The completion criterion: nothing public is left undocumented.
    @param[in] census The census document.
    @return None.
    """
    pending = [k for k, v in census["acknowledged_uncovered"].items()
               if isinstance(v, dict) and v.get("classification") == "public_pending"]
    assert pending == []


def test_every_entry_is_typed_and_reasoned(census):
    """!
    @brief No entry may be a bare string or an unreasoned classification.
    @param[in] census The census document.
    @return None.
    """
    for key, entry in census["acknowledged_uncovered"].items():
        assert isinstance(entry, dict), f"{key} is a legacy free-text entry"
        assert entry.get("classification"), key
        assert len(entry.get("reason", "").strip()) > 40, f"{key}: reason is too thin"


def test_parameter_and_alias_entries_name_a_real_owner(census, families):
    """!
    @brief An owner pointer must resolve to a family, a value, or a page anchor.
    @param[in] census The census document.
    @param[in] families Registered families.
    @return None.
    """
    for key, entry in census["acknowledged_uncovered"].items():
        if entry.get("classification") not in ("parameter_of_entry", "spelling_alias"):
            continue
        owner = entry.get("owner")
        assert owner, key
        head, _, tail = owner.partition("::")
        if head in families:
            if tail:
                assert tail in families[head]["value_metadata"], f"{key}: {owner}"
            continue
        page = PAGES / f"{head}.md"
        if page.is_file():
            assert not tail or f"@section {tail}" in page.read_text(encoding="utf-8"), \
                f"{key}: {owner} names no such anchor"
            continue
        # Otherwise it must name a module constant that exists.
        assert _constant(owner) is not None, f"{key}: owner '{owner}' resolves to nothing"


def test_single_valued_parameters_stay_single_valued(census):
    """!
    @brief A one-value parameter that grows a second value is a real choice.

    @details `poisson_solver.preconditioner.type` is classified as a parameter because
             there is nothing to choose between. The moment a second canonical value
             appears, that reasoning stops holding and it owes a capability family.
    @param[in] census The census document.
    @return None.
    """
    key = "picurv_cli.core::POISSON_PRECONDITIONER_TYPES"
    entry = census["acknowledged_uncovered"].get(key)
    assert entry and entry["classification"] == "parameter_of_entry"
    values = _constant(key)
    assert values is not None and len(values) == 1, (
        "POISSON_PRECONDITIONER_TYPES now offers a choice; register a capability "
        "family for it and reclassify the census entry"
    )


def test_structural_entries_are_not_closed_picurv_sets(census):
    """!
    @brief A `structural` classification claims PICurv does not own the vocabulary.
    @param[in] census The census document.
    @return None.
    """
    for key, entry in census["acknowledged_uncovered"].items():
        if entry.get("classification") != "structural":
            continue
        reason = entry["reason"].lower()
        assert "petsc" in reason or "free-form" in reason or "not a closed" in reason, \
            f"{key}: a structural claim must say who owns the vocabulary"


def test_discovery_covers_every_authoritative_constant():
    """!
    @brief Every named choice set is discovered, so none can hide as an inline literal.
    @return None.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "afc", TOOLING / "audit_family_census.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["afc"] = module
    spec.loader.exec_module(module)
    found = module.candidate_surfaces()
    for expected in (
        "picurv_cli.core::GRID_MODES",
        "picurv_cli.core::GRID_GENERATOR_TYPES",
        "picurv_cli.core::PARTICLE_RESTART_MODES",
        "picurv_cli.core::POST_EULERIAN_PIPELINE_TASKS",
        "picurv_cli.core::POST_LAGRANGIAN_PIPELINE_TASKS",
        "picurv_cli.core::POST_SPECTRA_TASKS",
        "picurv_cli.core::STUDY_TYPES",
        "picurv_cli.core::NEWTON_KRYLOV_PRECONDITIONER_MODELS",
        "picurv_cli.core::STATISTICS_WEIGHTING_MODES",
        "picurv_cli.storage::STORAGE_COMPRESSION_POLICIES",
    ):
        assert expected in found, f"{expected} is no longer discovered by the census"


def test_choice_sets_named_by_the_user_review_are_all_accounted_for(census, families):
    """!
    @brief Every surface the closeout review named must be covered or classified.
    @param[in] census The census document.
    @param[in] families Registered families.
    @return None.
    """
    covered_symbols = {
        f["public_surface"].get("symbol") for f in families.values()
    }
    acknowledged = {k.partition("::")[2] for k in census["acknowledged_uncovered"]}
    for symbol in (
        "GRID_MODES", "GRID_GENERATOR_TYPES", "PARTICLE_RESTART_MODES",
        "POST_EULERIAN_PIPELINE_TASKS", "POST_LAGRANGIAN_PIPELINE_TASKS",
        "POST_SPECTRA_TASKS", "POST_FIELD_STATISTICS_OUTPUTS",
        "POST_FIELD_STATISTICS_FORMATS", "STUDY_TYPES",
        "NEWTON_KRYLOV_PRECONDITIONER_MODELS",
        "NEWTON_KRYLOV_PRECONDITIONER_STRUCTURES", "KRYLOV_KSP_TYPES",
        "POISSON_PRECONDITIONER_TYPES", "STORAGE_COMPRESSION_POLICIES",
    ):
        assert symbol in covered_symbols or symbol in acknowledged, symbol
