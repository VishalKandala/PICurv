"""Contract tests for the spectra block of post.yml."""

import importlib.machinery
import importlib.util
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_core():
    """! @brief Load the conductor core module by path. @return Module. """
    path = ROOT / "picurv_cli" / "core.py"
    loader = importlib.machinery.SourceFileLoader("picurv_core_spectra_config_tests", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


CORE = load_core()


def recipe(**overrides):
    """! @brief Build a post.yml spectra block. @param[in] overrides Task overrides. @return Mapping. """
    task = {"task": "shell_spectrum"}
    task.update(overrides)
    return {"spectra": {"tasks": [task]}}


FACES = ["-Xi", "+Xi", "-Eta", "+Eta", "-Zeta", "+Zeta"]


def _case(block_bcs, blocks=1):
    """! @brief Build a case configuration around a boundary-condition layout. @param[in] block_bcs Layout. @param[in] blocks Block count. @return Mapping. """
    return {
        "properties": {"scaling": {"length_ref": 1.0, "velocity_ref": 1.0}},
        "models": {"domain": {"blocks": blocks}},
        "boundary_conditions": block_bcs,
    }


def periodic_case(blocks=1):
    """! @brief Build a triply periodic case configuration. @param[in] blocks Block count. @return Mapping. """
    one = [{"face": face, "type": "PERIODIC", "handler": "geometric"} for face in FACES]
    # A multi-block case must present one face list per block.
    return _case(one if blocks == 1 else [list(one) for _ in range(blocks)], blocks)


def walled_case():
    """! @brief Build a wall-bounded case configuration. @return Mapping. """
    return _case([{"face": face, "type": "WALL", "handler": "noslip"} for face in FACES])


def test_absent_block_normalizes_to_no_tasks():
    """! @brief A recipe without spectra must request nothing. """
    result = CORE.normalize_post_spectra_config({})
    assert result["tasks"] == []
    assert result["output_prefix"] == "Spectrum"


def test_defaults_are_filled_in():
    """! @brief Optional keys must resolve to the documented defaults. """
    result = CORE.normalize_post_spectra_config(recipe())
    assert result["tasks"] == [{
        "task": "shell_spectrum", "field": "Ucat", "block": 0,
        "symbol": "continuum", "subtract_mean": "none",
        "mean_source_step": None, "reference": None,
    }]


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"task": "imaginary_spectrum"}, "unknown task"),
        ({"field": "P"}, "is not supported"),
        ({"symbol": "fourier"}, "unknown symbol"),
        ({"subtract_mean": "median"}, "unknown subtract_mean"),
        ({"subtract_mean": "window:"}, "needs a window name"),
        ({"block": -1}, "non-negative integer"),
        ({"reference": "  "}, "must be a non-empty string"),
    ],
)
def test_malformed_tasks_are_rejected(overrides, message):
    """! @brief Each malformed key must be named in its own error. @param[in] overrides Task overrides. @param[in] message Expected text. """
    with pytest.raises(ValueError, match=message):
        CORE.normalize_post_spectra_config(recipe(**overrides))


def test_window_mean_form_is_accepted():
    """! @brief A window-scoped mean must survive normalization verbatim. """
    result = CORE.normalize_post_spectra_config(recipe(subtract_mean="window:early_decay"))
    assert result["tasks"][0]["subtract_mean"] == "window:early_decay"


def test_duplicate_task_identity_is_rejected():
    """! @brief Two identical tasks would overwrite one file, so they are refused. """
    cfg = {"spectra": {"tasks": [{"task": "shell_spectrum"}, {"task": "shell_spectrum"}]}}
    with pytest.raises(ValueError, match="more than once"):
        CORE.normalize_post_spectra_config(cfg)


def test_same_task_with_a_different_symbol_is_allowed():
    """! @brief Distinct binning abscissae write distinct files and must coexist. """
    cfg = {"spectra": {"tasks": [
        {"task": "shell_spectrum", "symbol": "continuum"},
        {"task": "shell_spectrum", "symbol": "discrete"},
    ]}}
    result = CORE.normalize_post_spectra_config(cfg)
    assert [entry["symbol"] for entry in result["tasks"]] == ["continuum", "discrete"]


def test_empty_task_list_is_rejected():
    """! @brief An empty task list is a mistake rather than a no-op. """
    with pytest.raises(ValueError, match="non-empty list"):
        CORE.normalize_post_spectra_config({"spectra": {"tasks": []}})


def test_periodic_box_satisfies_the_preconditions():
    """! @brief A triply periodic single-block case must be accepted. """
    normalized = CORE.normalize_post_spectra_config(recipe())
    assert CORE.validate_post_spectra_preconditions(normalized, periodic_case(), "post.yml") == []


def test_wall_bounded_case_is_refused_with_every_offending_face():
    """! @brief A wall-bounded case must be refused and its faces named. """
    normalized = CORE.normalize_post_spectra_config(recipe())
    errors = CORE.validate_post_spectra_preconditions(normalized, walled_case(), "post.yml")
    assert len(errors) == 1
    for face in ("-Xi", "+Xi", "-Eta", "+Eta", "-Zeta", "+Zeta"):
        assert face in errors[0]


def test_multi_block_domain_is_refused():
    """! @brief A shell spectrum is defined for one block only. """
    normalized = CORE.normalize_post_spectra_config(recipe())
    errors = CORE.validate_post_spectra_preconditions(normalized, periodic_case(blocks=2), "post.yml")
    assert any("single-block domain" in message for message in errors)


def test_block_outside_the_domain_is_refused():
    """! @brief A task naming a block the case lacks must be refused. """
    normalized = CORE.normalize_post_spectra_config(recipe(block=3))
    errors = CORE.validate_post_spectra_preconditions(normalized, periodic_case(), "post.yml")
    assert any("targets block 3" in message for message in errors)


def test_task_basename_separates_every_distinguishing_choice():
    """! @brief Two tasks differing in any identity key must not share a filename. """
    normalized = CORE.normalize_post_spectra_config({"spectra": {"tasks": [
        {"task": "shell_spectrum", "symbol": "continuum"},
        {"task": "shell_spectrum", "symbol": "discrete"},
    ]}})
    names = {CORE.post_spectra_task_basename(entry, "Spectrum") for entry in normalized["tasks"]}
    assert len(names) == 2
    assert "Spectrum_shell_spectrum_Ucat_block0000_continuum" in names


def test_stage_selection_defaults_to_every_stage():
    """! @brief An absent selector must run the whole post pipeline. """
    assert CORE.resolve_post_stage_selection(None) == set(CORE.POST_STAGE_NAMES)
    assert CORE.resolve_post_stage_selection("") == set(CORE.POST_STAGE_NAMES)


def test_stage_selection_accepts_a_subset():
    """! @brief A named subset must select exactly those stages. """
    assert CORE.resolve_post_stage_selection("spectra") == {"spectra"}
    assert CORE.resolve_post_stage_selection("spectra,fields") == {"spectra", "fields"}


@pytest.mark.parametrize("selector", ["curl", "spectra,curl", "spectra,"])
def test_stage_selection_rejects_unknown_or_empty_names(selector):
    """! @brief An unknown or empty stage must be refused. @param[in] selector Selector text. """
    with pytest.raises(ValueError):
        CORE.resolve_post_stage_selection(selector)


def test_signature_tracks_every_recipe_choice():
    """! @brief Any change that changes the output must change the signature. """
    base = CORE.normalize_post_spectra_config(recipe())
    baseline = CORE.compute_post_spectra_signature(base)
    assert CORE.compute_post_spectra_signature(CORE.normalize_post_spectra_config(recipe())) == baseline
    for overrides in ({"symbol": "discrete"}, {"subtract_mean": "domain"}, {"block": 1}):
        changed = CORE.normalize_post_spectra_config(recipe(**overrides))
        assert CORE.compute_post_spectra_signature(changed) != baseline


def test_recipe_config_carries_the_signature_only_when_spectra_are_requested():
    """! @brief The fingerprint must see spectra edits without inventing a key otherwise. """
    post = {"run_control": {"start_step": 0, "end_step": 10, "step_interval": 5}, "io": {}}
    assert "spectra_signature" not in CORE.build_post_recipe_config(dict(post))
    with_spectra = dict(post)
    with_spectra.update(recipe())
    built = CORE.build_post_recipe_config(with_spectra)
    assert built["spectra_signature"] == CORE.compute_post_spectra_signature(
        CORE.normalize_post_spectra_config(with_spectra)
    )


def test_predicted_artifacts_match_the_task_basenames(tmp_path):
    """! @brief Predicted CSV paths must be one per task under the spectra directory. @param[in] tmp_path Temp dir. """
    cfg = {"spectra": {"tasks": [
        {"task": "shell_spectrum", "symbol": "continuum"},
        {"task": "shell_spectrum", "symbol": "discrete"},
    ]}}
    paths = CORE.get_post_spectra_output_artifacts(cfg, str(tmp_path), {})
    assert len(paths) == 2
    for path in paths:
        assert os.path.dirname(path).endswith(CORE.CANONICAL_RUN_PATHS["spectra"])
    assert CORE.get_post_spectra_output_artifacts({}, str(tmp_path), {}) == []


def test_window_mean_resolves_to_the_checkpoint_payloads():
    """! @brief A window-scoped mean must resolve to the stored mean and count payloads. """
    bundle = {
        "bundle": "/runs/x/output/checkpoints/step_000000000150",
        "payloads": [
            {"kind": "mean", "field": "early_decay/Ucat_mean", "block": "0",
             "path": "statistics/window_0000/block_0000/Ucat_mean.dat"},
            {"kind": "occupancy", "field": "early_decay/count", "block": "0",
             "path": "statistics/window_0000/block_0000/count.dat"},
        ],
    }
    task = CORE.normalize_post_spectra_config(recipe(subtract_mean="window:early_decay"))["tasks"][0]
    args = CORE._spectra_mean_arguments(task, bundle)
    assert args[:2] == ["--subtract-mean", "field"]
    assert args[3].endswith("statistics/window_0000/block_0000/Ucat_mean.dat")
    assert args[5].endswith("statistics/window_0000/block_0000/count.dat")


def test_missing_window_payload_is_reported_by_name():
    """! @brief A window the bundle lacks must be named in the error. """
    bundle = {"bundle": "/runs/x/step", "payloads": []}
    task = CORE.normalize_post_spectra_config(recipe(subtract_mean="window:absent"))["tasks"][0]
    with pytest.raises(ValueError, match="absent/Ucat_mean"):
        CORE._spectra_mean_arguments(task, bundle)


def test_plain_mean_modes_pass_straight_through():
    """! @brief Modes needing no payload must not touch the bundle. """
    for mode in ("none", "domain"):
        task = CORE.normalize_post_spectra_config(recipe(subtract_mean=mode))["tasks"][0]
        assert CORE._spectra_mean_arguments(task, {"bundle": "", "payloads": []}) == \
            ["--subtract-mean", mode]


def test_mean_source_step_requires_a_window_mean():
    """! @brief Pinning a mean bundle is meaningless without a window mean. """
    with pytest.raises(ValueError, match="only applies to"):
        CORE.normalize_post_spectra_config(recipe(mean_source_step=400))


def test_mean_source_step_is_accepted_with_a_window_mean():
    """! @brief A pinned bundle must survive normalization. """
    result = CORE.normalize_post_spectra_config(
        recipe(subtract_mean="window:early_decay", mean_source_step=400)
    )
    assert result["tasks"][0]["mean_source_step"] == 400


@pytest.mark.parametrize("value", [-1, 1.5, True])
def test_malformed_mean_source_step_is_rejected(value):
    """! @brief A non-integer or negative pinned step must be refused. @param[in] value Candidate. """
    with pytest.raises(ValueError, match="mean_source_step"):
        CORE.normalize_post_spectra_config(
            recipe(subtract_mean="window:early_decay", mean_source_step=value)
        )


def checkpoint_bundle(tag):
    """! @brief Build a bundle carrying one window's mean payloads. @param[in] tag Bundle directory name. @return Mapping. """
    return {"bundle": f"/runs/x/{tag}", "payloads": [
            {"kind": "mean", "field": "early_decay/Ucat_mean", "block": "0",
             "path": "statistics/window_0000/block_0000/Ucat_mean.dat"},
        {"kind": "occupancy", "field": "early_decay/count", "block": "0",
         "path": "statistics/window_0000/block_0000/count.dat"},
    ]}


def test_pinned_mean_bundle_overrides_the_step_bundle():
    """! @brief A pinned bundle must supply the mean instead of the processed step's. """
    task = CORE.normalize_post_spectra_config(
        recipe(subtract_mean="window:early_decay", mean_source_step=400)
    )["tasks"][0]
    args = CORE._spectra_mean_arguments(
        task, checkpoint_bundle("step_150"), checkpoint_bundle("step_400")
    )
    assert "/runs/x/step_400/" in args[3]
    assert "/runs/x/step_400/" in args[5]


def test_pinned_mean_changes_the_recipe_signature():
    """! @brief Pinning a different mean bundle must change the fingerprint. """
    base = CORE.normalize_post_spectra_config(recipe(subtract_mean="window:w"))
    pinned = CORE.normalize_post_spectra_config(
        recipe(subtract_mean="window:w", mean_source_step=400)
    )
    assert CORE.compute_post_spectra_signature(base) != CORE.compute_post_spectra_signature(pinned)


def test_stage_measures_only_committed_steps(tmp_path, monkeypatch):
    """! @brief A window reaching past the last checkpoint must measure what exists. @param[in] tmp_path Temp dir. @param[in] monkeypatch Fixture. """
    calls = []
    monkeypatch.setattr(CORE, "_scan_committed_checkpoint_steps", lambda *a, **k: {0, 50})
    monkeypatch.setattr(CORE, "validate_committed_checkpoint",
                        lambda src, step: calls.append(step) or {
                            "bundle": str(tmp_path), "metadata": {"checkpoint_time": "0.0"},
                            "payloads": [{"kind": "eulerian", "field": "Ucat", "block": "0",
                                          "path": "eulerian/block_0000/Ucat.dat"}]})
    monkeypatch.setattr(CORE.os.path, "isfile", lambda p: True)

    class Result:
        returncode = 0
        stdout = '{"shell_spectrum": [], ' + ", ".join(
            f'"{name}": 0.0' for name in CORE.POST_SPECTRA_SCALAR_COLUMNS) + "}"
        stderr = ""
    monkeypatch.setattr(CORE.subprocess, "run", lambda *a, **k: Result())

    summary = CORE.run_post_spectra_stage(
        str(tmp_path), recipe(), {}, str(tmp_path), [0, 50, 100, 150], quiet=True
    )
    assert summary["steps"] == [0, 50]
    assert calls == [0, 50]


def test_stage_is_a_no_op_when_nothing_is_committed(tmp_path, monkeypatch):
    """! @brief A window with no committed checkpoint must return empty, not raise. @param[in] tmp_path Temp dir. @param[in] monkeypatch Fixture. """
    monkeypatch.setattr(CORE, "_scan_committed_checkpoint_steps", lambda *a, **k: set())
    monkeypatch.setattr(CORE.os.path, "isfile", lambda p: True)
    summary = CORE.run_post_spectra_stage(
        str(tmp_path), recipe(), {}, str(tmp_path), [0, 50], quiet=True
    )
    assert summary == {"tasks": [], "steps": [], "artifacts": []}


def test_follow_command_is_bare_and_targets_the_run(tmp_path):
    """! @brief The batch spectra step must not be wrapped in the MPI launcher. @param[in] tmp_path Temp dir. """
    cmd = CORE.build_spectra_follow_command(str(tmp_path), "post.yml", recipe())
    assert cmd, "a recipe requesting spectra must produce a follow command"
    assert "--only" in cmd and cmd[cmd.index("--only") + 1] == "spectra"
    assert "--run-dir" in cmd and "--post" in cmd
    # Running under srun/mpirun would launch one identical copy per task.
    assert not any(token in ("srun", "mpirun", "mpiexec") for token in cmd)


def test_no_follow_command_without_a_spectra_recipe(tmp_path):
    """! @brief A recipe requesting no spectra must add no batch step. @param[in] tmp_path Temp dir. """
    assert CORE.build_spectra_follow_command(str(tmp_path), "post.yml", {}) == []


def test_follow_commands_replace_exec_in_the_batch_script(tmp_path):
    """! @brief Trailing commands require dropping exec, which would end the script. @param[in] tmp_path Temp dir. """
    script = tmp_path / "job.sbatch"
    cluster = {"scheduler": {"type": "slurm"},
               "resources": {"account": "A", "nodes": 1, "ntasks_per_node": 1,
                             "mem": "1G", "time": "00:10:00"},
               "execution": {}}
    CORE.render_slurm_script(str(script), "job", cluster, ["main", "arg"], str(tmp_path),
                             str(tmp_path / "out.log"),
                             follow_commands=[["after", "one"]])
    body = script.read_text(encoding="utf-8")
    assert "exec main" not in body
    assert "\nmain arg\n" in body
    assert body.rstrip().endswith("after one")

    CORE.render_slurm_script(str(script), "job", cluster, ["main", "arg"], str(tmp_path),
                             str(tmp_path / "out.log"))
    assert "exec main arg" in script.read_text(encoding="utf-8")
