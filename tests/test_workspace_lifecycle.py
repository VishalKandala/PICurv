"""!
@file test_workspace_lifecycle.py
@brief Workspace topology, reusable asset, version, input, and recipe lifecycle tests.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

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
    shutil.copy2(REPO_ROOT / "examples" / "bent_channel" / "bent_channel_coarse.picgrid", grid)
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
