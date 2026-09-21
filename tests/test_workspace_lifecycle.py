"""!
@file test_workspace_lifecycle.py
@brief Workspace topology, reusable asset, version, input, and recipe lifecycle tests.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from picurv_cli import core
from picurv_cli.cli import build_main_parser


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "valid"


def _write_workspace(root: Path) -> Path:
    """!
    @brief Create a minimal initialized workspace for lifecycle tests.
    @param[in] root Workspace directory to create.
    @return Created workspace path.
    """
    root.mkdir()
    core.initialize_workspace_root(str(root), "test")
    return root


def _write_file_grid_case(workspace: Path):
    """!
    @brief Write a file-grid case whose grid can become a reusable asset.
    @param[in] workspace Initialized workspace path.
    @return Parsed case mapping and canonical case path.
    """
    grid = workspace / "inputs" / "grids" / "grid.picgrid"
    shutil.copy2(REPO_ROOT / "examples" / "search_robustness" / "bent_channel_coarse.picgrid", grid)
    case = yaml.safe_load((FIXTURES / "case.yml").read_text(encoding="utf-8"))
    case["title"] = "asset-case"
    case["grid"] = {"mode": "file", "source_file": "inputs/grids/grid.picgrid"}
    case_path = workspace / "config" / "case.yml"
    core.write_yaml_file(str(case_path), case)
    return case, case_path


def test_init_creates_uniform_workspace_and_uses_case_title(tmp_path):
    """!
    @brief Initialized cases expose the full stable workspace control surface.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = tmp_path / "channel"
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "picurv_cli" / "picurv"),
         "init", "flat_channel", "--dest", str(workspace)],
        cwd=tmp_path, text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    for relative in core.WORKSPACE_DIRECTORY_LAYOUT:
        assert (workspace / relative).is_dir(), relative
    marker = yaml.safe_load((workspace / core.WORKSPACE_CONFIG_FILENAME).read_text(encoding="utf-8"))
    assert marker["software"] == {}
    case_path = workspace / "config" / "case.yml"
    case = yaml.safe_load(case_path.read_text(encoding="utf-8"))
    assert case["title"] == "flat_channel"
    assert core.case_run_label(case, str(case_path)) == "flat_channel"


def test_precompute_publishes_and_run_reuses_content_addressed_grid(tmp_path):
    """!
    @brief Precompute publishes one immutable grid object that a run locks and reuses.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "workspace")
    case, case_path = _write_file_grid_case(workspace)

    published = core.precompute_case_assets(
        str(workspace), case, str(case_path), requested=["grid"]
    )
    reference = published["assets"]["grid"]
    object_root = workspace / reference["object"]
    assert (object_root / "asset.json").is_file()
    assert (object_root / "payload" / "inputs" / "grid" / "grid.run").is_file()
    assert core.plan_run_assets(case, str(case_path))["actions"][0]["action"] == "reuse"

    run_dir = workspace / "runs" / "asset-case_20260902-120000"
    core.ensure_run_layout(str(run_dir))
    lock = core.materialize_run_assets(
        str(run_dir), case, str(case_path), require_precomputed=True
    )
    assert lock["assets"]["grid"]["asset_id"] == reference["asset_id"]
    assert (run_dir / "inputs" / "grid" / "grid.run").is_file()
    assert yaml.safe_load((run_dir / "inputs" / "assets.lock.yml").read_text(encoding="utf-8"))[
        "assets"
    ]["grid"]["provider_spec_sha256"] == reference["provider_spec_sha256"]


def test_precompute_runtime_c_dependency_fails_without_partial_publication(tmp_path):
    """!
    @brief A Python provider depending on a C-only grid fails before publishing assets.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "workspace")
    case = yaml.safe_load((FIXTURES / "case.yml").read_text(encoding="utf-8"))
    case["properties"]["initial_conditions"] = {
        "mode": "generated", "generator": "ic_gen",
        "params": {"field": "Ucat", "config_file": "config/initial_conditions/expressions.cfg"},
    }
    case_path = workspace / "config" / "case.yml"
    core.write_yaml_file(str(case_path), case)

    with pytest.raises(ValueError, match="requires C generation"):
        core.precompute_case_assets(str(workspace), case, str(case_path))
    assert not list((workspace / "assets" / "objects").rglob("asset.json"))
    assert not list((workspace / "assets" / "sets").glob("*.yml"))


def test_workspace_version_pin_is_optional_but_enforced_when_present(tmp_path):
    """!
    @brief An empty software block follows the active release and an incompatible pin fails loudly.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "workspace")
    assert core.enforce_workspace_version(str(workspace))["release_version"] == core.PICURV_RELEASE_VERSION
    marker_path = workspace / core.WORKSPACE_CONFIG_FILENAME
    marker = yaml.safe_load(marker_path.read_text(encoding="utf-8"))
    marker["software"]["picurv"] = ">=99"
    core.write_yaml_file(str(marker_path), marker)
    with pytest.raises(ValueError, match="versions activate"):
        core.enforce_workspace_version(str(workspace))


def test_workspace_config_rejects_unsupported_keys(tmp_path):
    """!
    @brief An unrecognized top-level or nested key in `.picurv-workspace.yml` fails loudly.
    @details A silently ignored typo (`reproducibility.require_clena_release` for
             `require_clean_release`, say) would look configured but enforce nothing,
             the same failure mode the case/solver/monitor/post/cluster/study schemas
             already guard against.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "workspace")
    marker_path = workspace / core.WORKSPACE_CONFIG_FILENAME
    marker = yaml.safe_load(marker_path.read_text(encoding="utf-8"))
    marker["reproducibility"] = {"require_clean_release": True}
    core.write_yaml_file(str(marker_path), marker)
    assert core.load_workspace_config(str(workspace))["reproducibility"]["require_clean_release"] is True

    marker["typo_top_level_field"] = True
    core.write_yaml_file(str(marker_path), marker)
    with pytest.raises(ValueError, match="typo_top_level_field"):
        core.load_workspace_config(str(workspace))

    del marker["typo_top_level_field"]
    marker["reproducibility"] = {"require_clena_release": True}
    core.write_yaml_file(str(marker_path), marker)
    with pytest.raises(ValueError, match="require_clena_release"):
        core.load_workspace_config(str(workspace))


def test_input_import_modes_are_explicit_and_catalogued(tmp_path):
    """!
    @brief Copy and external-reference input modes create distinct durable catalog records.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "workspace")
    source = tmp_path / "source.picgrid"
    source.write_text("PICGRID\n", encoding="utf-8")
    copied = core.import_workspace_input(str(workspace), "grid", str(source), mode="copy")
    referenced = core.import_workspace_input(
        str(workspace), "reference-field", str(source), name="source.dat", mode="reference"
    )
    assert copied["path"] == "inputs/grids/source.picgrid"
    assert (workspace / copied["path"]).read_text(encoding="utf-8") == "PICGRID\n"
    assert referenced["path"].endswith(".reference.yml")
    reference = yaml.safe_load((workspace / referenced["path"]).read_text(encoding="utf-8"))
    assert reference["picurv_external_reference"] == str(source)
    catalog = yaml.safe_load((workspace / "inputs" / "catalog.yml").read_text(encoding="utf-8"))
    assert {item["mode"] for item in catalog["inputs"].values()} == {"copy", "reference"}


def test_restart_statistics_modes_and_post_recipes_have_stable_control_surfaces(tmp_path):
    """!
    @brief Restart statistics choices parse and distinct post recipes coexist under one run.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    parser = build_main_parser()
    for mode in ("reset", "carry"):
        args = parser.parse_args([
            "run", "--solve", "--statistics-state", mode,
            "--case", "case.yml", "--solver", "solver.yml", "--monitor", "monitor.yml",
        ])
        assert args.statistics_state == mode

    run_dir = tmp_path / "run"
    core.ensure_run_layout(str(run_dir))
    base = yaml.safe_load((FIXTURES / "post.yml").read_text(encoding="utf-8"))
    alternative = json.loads(json.dumps(base))
    alternative["io"]["eulerian_fields"] = ["Qcrit"]
    first = Path(core.generate_post_recipe_file(
        str(run_dir), "run", base, {"Post-Profile": str(FIXTURES / "post.yml")}
    ))
    second = Path(core.generate_post_recipe_file(
        str(run_dir), "run", alternative, {"Post-Profile": str(FIXTURES / "post.yml")}
    ))
    assert first != second
    assert first.is_file() and second.is_file()
    assert first.parent.parent == run_dir / "config" / "post-recipes"
    assert "output/visualization/" in first.read_text(encoding="utf-8")


def test_restart_statistics_carry_emits_native_continue_flag(tmp_path):
    """!
    @brief Carry mode crosses the Python/native boundary in the generated solver control.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "workspace")
    case, case_path = _write_file_grid_case(workspace)
    solver_path = FIXTURES / "solver.yml"
    monitor_path = FIXTURES / "monitor.yml"
    solver = core.read_yaml_file(str(solver_path))
    monitor = core.read_yaml_file(str(monitor_path))
    run_dir = workspace / "runs" / "carry"
    core.ensure_run_layout(str(run_dir))
    sources = {
        "Case": str(case_path), "Solver": str(solver_path), "Monitor": str(monitor_path)
    }
    monitor_files = core.prepare_monitor_files(str(run_dir), "carry", monitor, sources)
    control = core.generate_solver_control_file(
        str(run_dir), "carry",
        {
            "case": case, "case_path": str(case_path),
            "solver": solver, "solver_path": str(solver_path),
            "monitor": monitor, "monitor_path": str(monitor_path),
            "statistics_state": "carry",
        },
        1, monitor_files,
    )
    assert "-field_statistics_continue true" in Path(control).read_text(encoding="utf-8")


def _fake_binary(path: Path, output: str) -> Path:
    """!
    @brief Write an executable stub that prints one fixed --version line.
    @param[in] path Executable path to create.
    @param[in] output Line the stub prints on stdout.
    @return Created executable path.
    """
    path.write_text(f'#!/bin/sh\necho "{output}"\n', encoding="utf-8")
    path.chmod(0o755)
    return path


def test_binary_build_identity_is_read_back_from_the_executable(tmp_path):
    """!
    @brief The stamped identity is parsed from what the executable actually reports.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    clean = _fake_binary(tmp_path / "simulator", "simulator 0.2.0+gabcdef123456")
    dirty = _fake_binary(tmp_path / "postprocessor", "postprocessor 0.2.0+gabcdef123456.dirty")

    clean_identity = core.read_binary_build_identity(str(clean))
    assert clean_identity["available"] is True
    assert clean_identity["release_version"] == "0.2.0"
    assert clean_identity["git_commit"] == "abcdef123456"
    assert clean_identity["dirty"] is False
    assert core.read_binary_build_identity(str(dirty))["dirty"] is True

    # A binary predating the identity flag, or any other output, is recorded as
    # unavailable rather than guessed at: a wrong provenance claim is worse than none.
    legacy = _fake_binary(tmp_path / "legacy", "PICurv simulator")
    assert core.read_binary_build_identity(str(legacy))["available"] is False
    missing = core.read_binary_build_identity(str(tmp_path / "absent"))
    assert missing["available"] is False


def test_stale_runtime_binaries_are_reported_against_the_active_source(tmp_path, capsys, monkeypatch):
    """!
    @brief A binary built from another commit is named, not silently accepted.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] capsys Pytest capture fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    monkeypatch.setitem(core.PICURV_BUILD, "git_commit", "abcdef1234567890abcdef")
    monkeypatch.setitem(core.PICURV_BUILD, "dirty", False)
    monkeypatch.setitem(core.PICURV_BUILD, "build_id", "0.2.0+gabcdef123456")
    _fake_binary(tmp_path / "simulator", "simulator 0.2.0+gabcdef123456")
    _fake_binary(tmp_path / "postprocessor", "postprocessor 0.2.0+g999999999999")
    monkeypatch.setattr(core, "INVOKED_SCRIPT_DIR", str(tmp_path))

    identities = core.runtime_build_identities()
    assert identities["simulator"]["matches_source"] is True
    assert identities["postprocessor"]["matches_source"] is False

    stale = core.warn_on_stale_runtime_binaries(identities)
    assert stale == ["postprocessor"]
    assert "postprocessor was built from" in capsys.readouterr().err


def test_asset_identity_covers_case_values_the_build_reads(tmp_path):
    """!
    @brief A change that alters an asset's bytes must not be reported as reusable.

    @details The grid build nondimensionalizes by properties.scaling.length_ref, which
             lives outside the `grid:` subtree. When identity covered only that subtree,
             editing length_ref left the staleness check reporting `reuse` and the solver
             silently received geometry scaled by the wrong reference length.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "ws")
    case, case_path = _write_file_grid_case(workspace)
    core.precompute_case_assets(str(workspace), case, str(case_path), requested=["grid"])

    def actions(cfg):
        """!
        @brief Planned action per provider kind.
        @param[in] cfg Case mapping.
        @return Mapping of provider kind to planned action.
        """
        plan = core.plan_run_assets(cfg, str(case_path))
        return {item["kind"]: item["action"] for item in plan["actions"]}

    assert actions(case)["grid"] == "reuse"

    # An edit the build never reads must not invalidate a published object.
    unrelated = yaml.safe_load(case_path.read_text(encoding="utf-8"))
    unrelated.setdefault("run_control", {})["total_steps"] = 999
    assert actions(unrelated)["grid"] == "reuse"

    # An edit the build does read must.
    rescaled = yaml.safe_load(case_path.read_text(encoding="utf-8"))
    rescaled["properties"]["scaling"]["length_ref"] = float(
        rescaled["properties"]["scaling"]["length_ref"]
    ) * 2.0
    assert actions(rescaled)["grid"] == "build"


def test_asset_identity_follows_the_dependencies_it_declares(tmp_path):
    """!
    @brief A provider re-identifies when the asset it is built on top of changes.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    workspace = _write_workspace(tmp_path / "ws")
    case, case_path = _write_file_grid_case(workspace)
    graph = core.build_case_asset_graph(case, str(case_path))
    dependent = [p for p in graph["providers"] if p.get("dependencies")]
    if not dependent:
        pytest.skip("this case declares no dependent providers")
    before = {p["kind"]: p["spec_sha256"] for p in graph["providers"]}

    rescaled = yaml.safe_load(case_path.read_text(encoding="utf-8"))
    rescaled["properties"]["scaling"]["length_ref"] = float(
        rescaled["properties"]["scaling"]["length_ref"]
    ) * 2.0
    after = {p["kind"]: p["spec_sha256"]
             for p in core.build_case_asset_graph(rescaled, str(case_path))["providers"]}

    assert before["grid"] != after["grid"]
    for provider in dependent:
        assert before[provider["kind"]] != after[provider["kind"]], (
            f"{provider['kind']} declares a dependency on {provider['dependencies']} "
            "but did not re-identify when it changed"
        )


def test_a_refused_run_leaves_no_empty_run_directory(tmp_path):
    """!
    @brief A run refused during staging must not leave a skeleton behind.

    @details The layout is created before restart resolution can refuse, so every
             rejected attempt used to leave an empty directory in runs/ that is
             indistinguishable from a real run to anyone listing the workspace.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    run_dir = tmp_path / "runs" / "case_20260101-000000"
    core.ensure_run_layout(str(run_dir))
    assert run_dir.is_dir()

    assert core.discard_unused_run_directory(str(run_dir), created=True) is True
    assert not run_dir.exists()

    # A directory that already holds output is never removed by a staging failure,
    # and neither is one this invocation did not create.
    populated = tmp_path / "runs" / "populated"
    core.ensure_run_layout(str(populated))
    (populated / "logs" / "solver.log").write_text("output", encoding="utf-8")
    assert core.discard_unused_run_directory(str(populated), created=True) is False
    assert populated.is_dir()

    existing = tmp_path / "runs" / "existing"
    core.ensure_run_layout(str(existing))
    assert core.discard_unused_run_directory(str(existing), created=False) is False
    assert existing.is_dir()


def test_software_lock_pins_the_executables_a_run_used(tmp_path, monkeypatch):
    """!
    @brief The run records the bytes it executed, not only the source it was staged from.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    for name in ("simulator", "postprocessor"):
        path = fake_bin / name
        path.write_text(f'#!/bin/sh\necho "{name} 0.2.0+gabcdef123456"\n', encoding="utf-8")
        path.chmod(0o755)
    monkeypatch.setattr(core, "INVOKED_SCRIPT_DIR", str(fake_bin))

    run_dir = tmp_path / "runs" / "locked"
    core.ensure_run_layout(str(run_dir))
    core.write_software_lock(str(run_dir))

    lock = json.loads((run_dir / "inputs" / "software.lock.json").read_text(encoding="utf-8"))
    assert lock["picurv_version"] == core.PICURV_RELEASE_VERSION
    for name in ("simulator", "postprocessor"):
        entry = lock["executables"][name]
        assert entry["sha256"], f"{name} was not hashed"
        assert entry["build_id"] == "0.2.0+gabcdef123456"
    assert lock["python_conductor_sha256"]
    assert lock["generators"], "no generator was hashed"


def test_reproducibility_policy_refuses_a_development_build(tmp_path, monkeypatch):
    """!
    @brief An opt-in workspace policy refuses staging from a modified or untagged tree.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "ws")
    config = workspace / core.WORKSPACE_CONFIG_FILENAME
    payload = yaml.safe_load(config.read_text(encoding="utf-8"))
    payload["reproducibility"] = {"require_clean_release": True}
    core.write_yaml_file(str(config), payload)

    monkeypatch.setitem(core.PICURV_BUILD, "dirty", True)
    monkeypatch.setitem(core.PICURV_BUILD, "dev_distance", 12)
    monkeypatch.setitem(core.PICURV_BUILD, "build_id", "0.2.0.dev12+gabc.dirty")
    with pytest.raises(ValueError, match="require_clean_release"):
        core.enforce_reproducibility_policy(str(workspace))

    # A clean release satisfies it, and no policy at all is never enforced.
    monkeypatch.setitem(core.PICURV_BUILD, "dirty", False)
    monkeypatch.setitem(core.PICURV_BUILD, "dev_distance", 0)
    assert core.enforce_reproducibility_policy(str(workspace))["require_clean_release"] is True
    payload.pop("reproducibility")
    core.write_yaml_file(str(config), payload)
    assert core.enforce_reproducibility_policy(str(workspace)) == {}


def test_published_assets_carry_inspection_material(tmp_path):
    """!
    @brief A published asset is something the user can look at, not only solver input.

    @details Precompute exists so a grid, field, or profile can be checked and changed
             before a solve is committed. Objects that carried only opaque payload
             could not serve that purpose.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return None.
    """
    workspace = _write_workspace(tmp_path / "ws")
    case, case_path = _write_file_grid_case(workspace)
    published = core.precompute_case_assets(
        str(workspace), case, str(case_path), requested=["grid"]
    )
    reference = published["assets"]["grid"]
    obj = workspace / reference["object"]

    assert "validation.json" in reference["inspection"]
    validation = json.loads((obj / "validation.json").read_text(encoding="utf-8"))
    assert validation["asset_kind"] == "grid"
    geometry = validation["geometry"]
    assert geometry["total_nodes"] > 0
    assert len(geometry["bounds_min"]) == 3 and len(geometry["bounds_max"]) == 3
    assert all(
        high >= low for low, high in zip(geometry["bounds_min"], geometry["bounds_max"])
    )

    preview = obj / "preview.vts"
    assert preview.is_file(), sorted(p.name for p in obj.iterdir())
    assert preview.read_text(encoding="utf-8").lstrip().startswith("<?xml")

    # Inspection material describes the payload; it must not change which object a run
    # resolves to, or every asset would re-identify whenever a preview format changed.
    again = core.precompute_case_assets(
        str(workspace), case, str(case_path), requested=["grid"]
    )
    assert again["assets"]["grid"]["asset_id"] == reference["asset_id"]


def test_build_identity_problems_reports_stale_and_workspace_mismatch(tmp_path, monkeypatch):
    """!
    @brief Coherence problems are stated, covering executables and the workspace pin.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    monkeypatch.setitem(core.PICURV_BUILD, "git_commit", "abcdef1234567890abcdef")
    monkeypatch.setitem(core.PICURV_BUILD, "dirty", False)
    monkeypatch.setitem(core.PICURV_BUILD, "build_id", "0.2.0+gabcdef123456")
    _fake_binary(tmp_path / "simulator", "simulator 0.2.0+gabcdef123456")
    _fake_binary(tmp_path / "postprocessor", "postprocessor 0.2.0+g999999999999")
    monkeypatch.setattr(core, "INVOKED_SCRIPT_DIR", str(tmp_path))
    identities = core.runtime_build_identities()

    problems = core.build_identity_problems(identities)
    assert len(problems) == 1
    assert problems[0].startswith("postprocessor: built from")

    unsatisfiable = f">{core.PICURV_RELEASE_VERSION}"
    problems = core.build_identity_problems(identities, unsatisfiable)
    assert any(problem.startswith("workspace: requires") for problem in problems)
    assert any(problem.startswith("workspace: ") and "not a valid" in problem
               for problem in core.build_identity_problems(identities, "not-a-version"))


def test_version_status_exits_non_zero_when_the_build_is_incoherent(tmp_path):
    """!
    @brief `version status` validates and fails; bare `version` reports and succeeds.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    picurv = str(REPO_ROOT / "picurv_cli" / "picurv")
    reported = subprocess.run(
        [sys.executable, picurv, "version"], cwd=tmp_path,
        text=True, capture_output=True, check=False,
    )
    assert reported.returncode == 0
    assert "PICurv release" in reported.stdout

    validated = subprocess.run(
        [sys.executable, picurv, "version", "status", "--format", "json"], cwd=tmp_path,
        text=True, capture_output=True, check=False,
    )
    payload = json.loads(validated.stdout)
    assert payload["coherent"] is (validated.returncode == 0)
    assert payload["coherent"] == (not payload["problems"])


def test_a_branched_run_records_the_parent_it_started_from(tmp_path):
    """!
    @brief The manifest names the parent run, its checkpoint, and the statistics decision.

    @details The parent's identity is read from its own manifest, so a parent that was
             renamed is still named correctly rather than reported by directory.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    workspace = _write_workspace(tmp_path / "ws")
    parent = workspace / "runs" / "parent_20260101-000000"
    core.ensure_run_layout(str(parent))
    (parent / "manifest.json").write_text(
        json.dumps({"artifact_type": "run", "run_id": "parent_20260101-000000"}) + "\n",
        encoding="utf-8",
    )
    renamed = parent.parent / "moved_aside"
    parent.rename(renamed)

    lineage = core.build_run_lineage(
        str(renamed), 250, workspace_root=str(workspace),
        statistics_state="carry", requested_source="latest",
    )
    assert lineage["relationship"] == "branch"
    assert lineage["parent_run_id"] == "parent_20260101-000000"
    assert lineage["parent_identity_source"] == "manifest"
    assert lineage["checkpoint_step"] == 250
    assert lineage["statistics_state"] == "carry"
    assert lineage["requested_source"] == "latest"

    child = workspace / "runs" / "child_20260102-000000"
    core.ensure_run_layout(str(child))
    manifest = core.build_run_manifest(
        str(child), "child_20260102-000000", workspace_root=str(workspace), lineage=lineage
    )
    assert manifest["lineage"] == lineage

    # A fresh run says so, rather than leaving the reader to infer it from a missing key.
    assert core.build_run_manifest(str(child), "child_20260102-000000")["lineage"] == {
        "relationship": "root"
    }


def test_a_rebuilt_manifest_keeps_the_lineage_the_run_was_created_with(tmp_path):
    """!
    @brief Resuming a branched run does not erase what it branched from.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    workspace = _write_workspace(tmp_path / "ws")
    run = workspace / "runs" / "child_20260102-000000"
    core.ensure_run_layout(str(run))
    lineage = {"relationship": "branch", "parent_run_id": "parent_20260101-000000",
               "checkpoint_step": 250}
    (run / "manifest.json").write_text(
        json.dumps(core.build_run_manifest(
            str(run), "child_20260102-000000", workspace_root=str(workspace), lineage=lineage
        )) + "\n",
        encoding="utf-8",
    )

    rebuilt = core.build_run_manifest(
        str(run), "child_20260102-000000", workspace_root=str(workspace)
    )
    assert rebuilt["lineage"] == lineage


def test_a_run_root_refuses_a_directory_outside_the_layout(tmp_path):
    """!
    @brief An unrouted peer directory is an error; an unexpected file is only reported.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    workspace = _write_workspace(tmp_path / "ws")
    run = workspace / "runs" / "run_20260101-000000"
    core.ensure_run_layout(str(run))
    assert core.validate_run_directory_structure(str(run)) == ([], [])

    (run / "diagnostics").mkdir()
    errors, warnings = core.validate_run_directory_structure(str(run))
    assert len(errors) == 1
    assert "'diagnostics/' is not part of the run layout" in errors[0]
    assert warnings == []

    (run / "diagnostics").rmdir()
    (run / "notes.txt").write_text("kept\n", encoding="utf-8")
    errors, warnings = core.validate_run_directory_structure(str(run))
    assert errors == []
    assert len(warnings) == 1 and "notes.txt" in warnings[0]

    # The files a run legitimately carries at its root are not reported as strays.
    (run / "notes.txt").unlink()
    (run / "manifest.json").write_text("{}\n", encoding="utf-8")
    (run / core.STORAGE_STATE_FILENAME).write_text("{}\n", encoding="utf-8")
    assert core.validate_run_directory_structure(str(run)) == ([], [])


def test_a_version_pin_failure_names_the_installation_it_would_change(tmp_path, monkeypatch):
    """!
    @brief The refusal states that one shared installation serves every workspace.

    @details PICurv installs once and `versions activate` rewrites that checkout in
             place, so satisfying one workspace's pin re-points every other. A message
             that only said "activate a matching installation" would read as though
             several could coexist.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    workspace = _write_workspace(tmp_path / "ws")
    config_path = workspace / core.WORKSPACE_CONFIG_FILENAME
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    payload["software"] = {"picurv": f">{core.PICURV_RELEASE_VERSION}"}
    config_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError) as failure:
        core.enforce_workspace_version(str(workspace))

    message = str(failure.value)
    assert "single shared installation" in message
    assert core.PACKAGE_PROJECT_ROOT in message
    assert "--pin-executables" in message


def _git_repo_with_tag(root: Path, tag: str) -> Path:
    """!
    @brief Create a one-commit Git repository carrying one tag.
    @param[in] root Repository directory to create.
    @param[in] tag Tag name to attach to the commit.
    @return Repository root.
    """
    root.mkdir(parents=True)
    for command in (
        ["git", "init", "-q"],
        ["git", "config", "user.email", "test@example.invalid"],
        ["git", "config", "user.name", "test"],
        ["git", "commit", "-q", "--allow-empty", "-m", "initial"],
        ["git", "tag", tag],
    ):
        subprocess.run(command, cwd=root, check=True, capture_output=True)
    return root


def test_a_bare_release_resolves_to_its_v_prefixed_tag(tmp_path, monkeypatch):
    """!
    @brief `VERSION` and workspace pins say `0.1.0`; the release tag is `v0.1.0`.

    @details `versions activate` with no argument reads the bare release from the
             workspace pin. Checking that out verbatim fails, because no ref carries
             the unprefixed name.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    repo = _git_repo_with_tag(tmp_path / "source", "v1.2.3")
    monkeypatch.setattr(core, "PACKAGE_PROJECT_ROOT", str(repo))

    assert core._resolve_version_ref("1.2.3") == "v1.2.3"
    assert core._resolve_version_ref("v1.2.3") == "v1.2.3"
    with pytest.raises(ValueError, match="'9.9.9' or 'v9.9.9'"):
        core._resolve_version_ref("9.9.9")


def _stub_version_install(monkeypatch, identities_by_call):
    """!
    @brief Replace the Git, build, and identity owners that `versions install` drives.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @param[in] identities_by_call Executable identities returned after the build.
    @return Mapping recording the make arguments and checked-out ref.
    """
    recorded = {}
    monkeypatch.setattr(core, "_require_clean_source_checkout", lambda action: None)
    monkeypatch.setattr(core, "_git_source_command",
                        lambda arguments: recorded.setdefault("git", []).append(arguments))
    monkeypatch.setattr(core, "_resolve_version_ref", lambda version: f"v{version}")
    monkeypatch.setattr(core, "run_project_make",
                        lambda root, make_args: recorded.__setitem__("make_args", make_args))
    fresh = {"git_commit": "fedcba9876543210", "dirty": False, "build_id": "1.2.3+gfedcba987654"}
    monkeypatch.setattr(core, "_source_build_identity", lambda release: dict(fresh))
    monkeypatch.setattr(core, "runtime_build_identities", lambda source_identity=None: identities_by_call)
    return recorded


def test_versions_install_passes_site_build_settings_to_make(tmp_path, monkeypatch):
    """!
    @brief A cluster install reaches `make` with its SYSTEM, as `picurv build` does.

    @details The install used to run a bare `make all`, so a cluster build silently used
             the local configuration unless SYSTEM happened to be exported in the shell.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    matching = {name: {"available": True, "matches_source": True, "build_id": "1.2.3+gfedcba987654"}
                for name in ("simulator", "postprocessor")}
    recorded = _stub_version_install(monkeypatch, matching)
    args = build_main_parser().parse_args(["versions", "install", "1.2.3", "SYSTEM=cluster", "-j4"])

    core.versions_workflow(args)

    assert recorded["make_args"] == ["SYSTEM=cluster", "-j4"]
    assert ["checkout", "--detach", "v1.2.3"] in recorded["git"]


def test_versions_activate_reads_a_leading_assignment_as_a_make_argument(tmp_path, monkeypatch):
    """!
    @brief `versions activate SYSTEM=cluster` keeps the workspace version and the build setting.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    matching = {name: {"available": True, "matches_source": True, "build_id": "1.2.3+gfedcba987654"}
                for name in ("simulator", "postprocessor")}
    recorded = _stub_version_install(monkeypatch, matching)
    monkeypatch.setattr(core, "_workspace_requested_version", lambda root: "1.2.3")
    args = build_main_parser().parse_args(
        ["versions", "activate", "--workspace", str(tmp_path), "SYSTEM=cluster"]
    )

    core.versions_workflow(args)

    assert recorded["make_args"] == ["SYSTEM=cluster"]
    assert ["checkout", "--detach", "v1.2.3"] in recorded["git"]


def test_versions_install_refuses_a_make_target(tmp_path, monkeypatch):
    """!
    @brief A make target would replace the build the install is supposed to verify.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    _stub_version_install(monkeypatch, {})
    args = build_main_parser().parse_args(["versions", "install", "1.2.3", "clean-project"])
    with pytest.raises(ValueError, match="not a target"):
        core.versions_workflow(args)


def test_versions_install_fails_when_the_build_does_not_carry_the_new_identity(tmp_path, monkeypatch):
    """!
    @brief Success is reported only after the executables are checked against the new source.

    @details The conductor's own PICURV_BUILD still describes the commit it started on,
             so the check must use the identity read after the checkout moved.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    stale = {
        "simulator": {"available": True, "matches_source": False, "build_id": "1.2.2+g0123456789ab"},
        "postprocessor": {"available": True, "matches_source": True, "build_id": "1.2.3+gfedcba987654"},
    }
    _stub_version_install(monkeypatch, stale)
    args = build_main_parser().parse_args(["versions", "install", "1.2.3"])
    with pytest.raises(ValueError, match=r"simulator: built from 1\.2\.2\+g0123456789ab.*1\.2\.3\+gfedcba987654"):
        core.versions_workflow(args)


def test_a_binary_that_cannot_start_is_not_reported_as_predating_identity(tmp_path):
    """!
    @brief A failing `--version` records the exit status and first error line.

    @details A run manifest recorded "no build identity reported" for a binary that does
             report one; a binary that cannot load its libraries in the staging shell was
             indistinguishable from one built before the identity flag existed.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    broken = tmp_path / "simulator"
    broken.write_text('#!/bin/sh\necho "error while loading shared libraries: libpetsc.so" >&2\nexit 127\n',
                      encoding="utf-8")
    broken.chmod(0o755)

    identity = core.read_binary_build_identity(str(broken))

    assert identity["available"] is False
    assert identity["reason"].startswith("--version exited 127")
    assert "libpetsc.so" in identity["reason"]


def test_staging_warns_when_a_build_identity_cannot_be_read(tmp_path, capsys):
    """!
    @brief An unreadable identity is said at staging, not only recorded in the manifest.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] capsys Pytest capture fixture.
    """
    identities = {
        "simulator": {"available": False, "reason": "--version exited 127", "path": "/x/simulator"},
    }
    assert core.warn_on_stale_runtime_binaries(identities) == []
    err = capsys.readouterr().err
    assert "build identity of simulator could not be read" in err
    assert "--version exited 127" in err


def _installation_with_executables(tmp_path: Path, monkeypatch, commit: str = "0123456789ab") -> Path:
    """!
    @brief Point executable resolution at stub simulator and postprocessor binaries.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @param[in] commit Commit the stubs report.
    @return Directory holding the stubs.
    """
    bin_dir = tmp_path / "installation_bin"
    bin_dir.mkdir(exist_ok=True)
    for name in core.RUN_EXECUTABLE_NAMES:
        _fake_binary(bin_dir / name, f"{name} 0.1.0+g{commit}")
    monkeypatch.setattr(core, "INVOKED_SCRIPT_DIR", str(bin_dir))
    return bin_dir


def _run_with_initial_config(root: Path) -> Path:
    """!
    @brief Create a run directory holding an initial configuration record.
    @param[in] root Run directory to create.
    @return Run directory.
    """
    (root / "config").mkdir(parents=True)
    core.write_json_file(str(root / "config" / "active.json"),
                         {"schema_version": 1, "revision": "initial", "files": {}})
    return root


def test_pinned_executables_survive_a_rebuild_of_the_installation(tmp_path, monkeypatch):
    """!
    @brief A run launches the copy made at staging, not what `bin/` holds later.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    bin_dir = _installation_with_executables(tmp_path, monkeypatch)
    run = _run_with_initial_config(tmp_path / "run")

    recorded = core.pin_run_executables(str(run), "initial")
    _fake_binary(bin_dir / "simulator", "simulator 0.1.0+gffffffffffff")

    assert recorded["simulator"]["path"] == "config/bin/simulator"
    assert recorded["simulator"]["version_line"] == "simulator 0.1.0+g0123456789ab"
    pinned = core.run_executable_path(str(run), "simulator")
    assert pinned == str(run / "config" / "bin" / "simulator")
    assert core.read_binary_build_identity(pinned)["build_id"] == "0.1.0+g0123456789ab"
    assert core.runtime_build_identities(run_dir=str(run))["simulator"]["build_id"] == "0.1.0+g0123456789ab"
    assert core.build_software_lock(str(run))["executables"]["simulator"]["path"] == pinned


def test_a_continuation_keeps_the_pin_unless_it_re_pins(tmp_path, monkeypatch):
    """!
    @brief A continuation revision carries the pin forward; an explicit re-pin records a new one.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    bin_dir = _installation_with_executables(tmp_path, monkeypatch)
    run = tmp_path / "run"
    (run / "config").mkdir(parents=True)
    source = tmp_path / "case.yml"
    source.write_text("title: x\n", encoding="utf-8")
    core.snapshot_run_configuration(str(run), {"case": str(source)})
    core.pin_run_executables(str(run), "initial")
    _fake_binary(bin_dir / "simulator", "simulator 0.1.0+gffffffffffff")

    continued = core.snapshot_run_configuration(str(run), {"case": str(source)}, continuation=True)
    kept = core.apply_run_executable_pinning(
        str(run), continued["revision"], SimpleNamespace(pin_executables=None), None
    )
    assert kept["simulator"]["path"] == str(run / "config" / "bin" / "simulator")

    with pytest.raises(ValueError, match="--no-pin-executables would launch"):
        core.apply_run_executable_pinning(
            str(run), continued["revision"], SimpleNamespace(pin_executables=False), None
        )

    moved = core.apply_run_executable_pinning(
        str(run), continued["revision"], SimpleNamespace(pin_executables=True), None
    )
    assert moved["simulator"]["path"] == str(
        run / "config" / "history" / continued["revision"] / "bin" / "simulator"
    )
    assert moved["simulator"]["build_id"] == "0.1.0+gffffffffffff"
    assert core.read_binary_build_identity(str(run / "config" / "bin" / "simulator"))["build_id"] \
        == "0.1.0+g0123456789ab"


def test_pinning_is_opt_in_through_the_flag_or_the_workspace(tmp_path, monkeypatch):
    """!
    @brief Nothing is pinned unless a switch or the workspace asks, and a switch wins.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    unset = SimpleNamespace(pin_executables=None)
    assert core.resolve_executable_pinning(unset, None) is False
    assert core.resolve_executable_pinning(SimpleNamespace(pin_executables=True), None) is True

    workspace = _write_workspace(tmp_path / "ws")
    assert core.resolve_executable_pinning(unset, str(workspace)) is False
    config_path = workspace / core.WORKSPACE_CONFIG_FILENAME
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    payload["reproducibility"] = {"pin_executables": True}
    config_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    assert core.resolve_executable_pinning(unset, str(workspace)) is True
    assert core.resolve_executable_pinning(SimpleNamespace(pin_executables=False), str(workspace)) is False


def test_an_unbuilt_executable_is_left_unpinned_with_a_warning(tmp_path, monkeypatch, capsys):
    """!
    @brief Staging without a built postprocessor still pins the simulator and says so.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @param[in] capsys Pytest capture fixture.
    """
    bin_dir = _installation_with_executables(tmp_path, monkeypatch)
    (bin_dir / "postprocessor").unlink()
    monkeypatch.setattr(core, "DEFAULT_BIN_DIR", str(bin_dir))
    run = _run_with_initial_config(tmp_path / "run")

    recorded = core.pin_run_executables(str(run), "initial")

    assert sorted(recorded) == ["simulator"]
    assert "postprocessor is not built" in capsys.readouterr().err
    assert core.run_executable_path(str(run), "postprocessor") == str(bin_dir / "postprocessor")


def _run_identity_check(tmp_path: Path, executable: Path, expected) -> subprocess.CompletedProcess:
    """!
    @brief Execute the generated job-start identity check in a real shell.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] executable Executable the check probes.
    @param[in] expected Version line captured at staging, or None.
    @return Completed shell process.
    """
    script = tmp_path / "check.sh"
    lines = ["#!/bin/bash", "set -euo pipefail",
             *core.executable_identity_check_lines(str(executable), expected), "echo LAUNCHED"]
    script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return subprocess.run(["bash", str(script)], capture_output=True, text=True)


def test_a_job_refuses_to_launch_an_executable_rebuilt_after_staging(tmp_path):
    """!
    @brief The job-start check passes the staged build and stops a different one before launch.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    exe = _fake_binary(tmp_path / "simulator", "simulator 0.1.0+g0123456789ab")

    same = _run_identity_check(tmp_path, exe, "simulator 0.1.0+g0123456789ab")
    assert same.returncode == 0 and "LAUNCHED" in same.stdout

    _fake_binary(exe, "simulator 0.1.0+gffffffffffff")
    rebuilt = _run_identity_check(tmp_path, exe, "simulator 0.1.0+g0123456789ab")
    assert rebuilt.returncode == 1
    assert "LAUNCHED" not in rebuilt.stdout
    assert "rebuilt after staging" in rebuilt.stderr

    exe.write_text("#!/bin/sh\necho 'cannot load libpetsc.so' >&2\nexit 127\n", encoding="utf-8")
    broken = _run_identity_check(tmp_path, exe, "simulator 0.1.0+g0123456789ab")
    assert broken.returncode == 1 and "--version failed" in broken.stderr

    _fake_binary(exe, "simulator 0.1.0+gffffffffffff")
    unknown = _run_identity_check(tmp_path, exe, None)
    assert unknown.returncode == 0 and "LAUNCHED" in unknown.stdout


def test_slurm_script_checks_identity_after_modules_and_before_launch(tmp_path):
    """!
    @brief The check needs the modules' shared libraries and must precede the launcher.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    script = tmp_path / "solver.sbatch"
    cluster_cfg = {
        "resources": {"nodes": 1, "ntasks_per_node": 2, "mem": "1G", "time": "00:10:00", "account": "a"},
        "execution": {"module_setup": ["module restore petsc-prod"]},
    }
    core.render_slurm_script(
        str(script), "job", cluster_cfg, ["mpirun", "-np", "2", "/x/config/bin/simulator"],
        str(tmp_path), str(tmp_path / "o.out"),
        identity_check=("/x/config/bin/simulator", "simulator 0.1.0+g0123456789ab"),
    )
    text = script.read_text(encoding="utf-8")
    assert text.index("module restore petsc-prod") < text.index("--version") < text.index("exec mpirun")


def test_study_members_launch_their_own_pins_through_one_array_script(tmp_path, monkeypatch):
    """!
    @brief Pinned members share a run-relative path; a half-pinned study is refused.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    """
    _installation_with_executables(tmp_path, monkeypatch)
    members = [_run_with_initial_config(tmp_path / f"case_{i}") for i in range(2)]
    for member in members:
        core.pin_run_executables(str(member), "initial")

    chosen = core.resolve_sweep_stage_executables([str(m) for m in members])
    assert chosen["simulator"] == ("$RUN_DIR/config/bin/simulator", "simulator 0.1.0+g0123456789ab")

    unpinned = _run_with_initial_config(tmp_path / "case_2")
    with pytest.raises(ValueError, match="disagree"):
        core.resolve_sweep_stage_executables([str(m) for m in members] + [str(unpinned)])


def _git(arguments, cwd):
    """!
    @brief Run a git command for a version-workflow fixture and fail loudly if it fails.
    @param[in] arguments Git arguments excluding the executable.
    @param[in] cwd Working directory.
    @return Captured standard output.
    """
    result = subprocess.run(["git", *arguments], cwd=str(cwd), text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    return result.stdout


def _source_checkout_with_origin(tmp_path: Path):
    """!
    @brief Build a real origin repository, a clone of it, and a second clone that can push.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @return Tuple of (checkout the conductor manages, publisher clone).
    """
    seed = tmp_path / "seed"
    seed.mkdir()
    (seed / "VERSION").write_text("1.0.0\n", encoding="utf-8")
    for arguments in (["init", "-q"], ["config", "user.email", "tests@example.com"],
                      ["config", "user.name", "PICurv Tests"], ["add", "."],
                      ["commit", "-q", "-m", "first"], ["tag", "v1.0.0"]):
        _git(arguments, seed)
    origin = tmp_path / "origin.git"
    _git(["clone", "-q", "--bare", str(seed), str(origin)], tmp_path)
    checkout = tmp_path / "checkout"
    publisher = tmp_path / "publisher"
    for clone in (checkout, publisher):
        _git(["clone", "-q", str(origin), str(clone)], tmp_path)
        _git(["config", "user.email", "tests@example.com"], clone)
        _git(["config", "user.name", "PICurv Tests"], clone)
    return checkout, publisher


def test_source_update_fetches_new_tags_without_moving_the_checkout(tmp_path, monkeypatch, capsys):
    """!
    @brief `source update` makes a new release visible and leaves the running code alone.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @param[in] capsys Pytest output-capture fixture.
    @return None.
    """
    checkout, publisher = _source_checkout_with_origin(tmp_path)
    head_before = _git(["rev-parse", "HEAD"], checkout).strip()
    (publisher / "VERSION").write_text("1.1.0\n", encoding="utf-8")
    _git(["commit", "-q", "-am", "second"], publisher)
    _git(["tag", "v1.1.0"], publisher)
    _git(["push", "-q", "origin", "HEAD", "--tags"], publisher)
    monkeypatch.setattr(core, "PACKAGE_PROJECT_ROOT", str(checkout))

    core.source_workflow(build_main_parser().parse_args(["source", "update"]))

    assert "v1.1.0" in _git(["tag", "--list"], checkout).split()
    assert _git(["rev-parse", "HEAD"], checkout).strip() == head_before
    assert (checkout / "VERSION").read_text(encoding="utf-8") == "1.0.0\n"
    assert "active checkout was not changed" in capsys.readouterr().out


def test_source_update_reports_an_unreachable_remote(tmp_path, monkeypatch):
    """!
    @brief A fetch failure is an error, not a silent no-op.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    checkout, _ = _source_checkout_with_origin(tmp_path)
    monkeypatch.setattr(core, "PACKAGE_PROJECT_ROOT", str(checkout))
    args = build_main_parser().parse_args(["source", "update", "--remote", "nowhere"])
    with pytest.raises(ValueError):
        core.source_workflow(args)


def test_versions_list_orders_tags_by_version_not_text(tmp_path, monkeypatch, capsys):
    """!
    @brief `versions list` names the active build and every tag, newest release first.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @param[in] capsys Pytest output-capture fixture.
    @return None.
    """
    checkout, _ = _source_checkout_with_origin(tmp_path)
    for tag in ("v1.2.0", "v1.10.0", "v1.9.0"):
        _git(["tag", tag], checkout)
    monkeypatch.setattr(core, "PACKAGE_PROJECT_ROOT", str(checkout))

    core.versions_workflow(build_main_parser().parse_args(["versions", "list"]))

    out = capsys.readouterr().out
    assert out.startswith(f"Active: {core.PICURV_BUILD['build_id']}")
    listed = [line.strip() for line in out.splitlines() if line.startswith("  ")]
    assert listed == ["v1.10.0", "v1.9.0", "v1.2.0", "v1.0.0"]


def test_versions_activate_reads_a_leading_option_as_a_make_argument(tmp_path, monkeypatch):
    """!
    @brief `versions activate -- -j8` keeps the workspace version and passes -j8 to make.
    @details argparse binds '-j8' to the optional version positional; git then received
             it as a ref and printed its usage.
    @param[in] tmp_path Pytest temporary-directory fixture.
    @param[in] monkeypatch Pytest monkeypatch fixture.
    @return None.
    """
    matching = {name: {"available": True, "matches_source": True, "build_id": "1.2.3+gfedcba987654"}
                for name in ("simulator", "postprocessor")}
    recorded = _stub_version_install(monkeypatch, matching)
    monkeypatch.setattr(core, "_workspace_requested_version", lambda root: "1.2.3")
    args = build_main_parser().parse_args(
        ["versions", "activate", "--workspace", str(tmp_path), "--", "-j8"]
    )

    core.versions_workflow(args)

    assert recorded["make_args"] == ["-j8"]
    assert ["checkout", "--detach", "v1.2.3"] in recorded["git"]
