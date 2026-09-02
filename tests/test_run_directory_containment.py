#!/usr/bin/env python3
"""Regression tests for run-directory containment.

These guard a data-loss path: `io.directories.log` is passed through to `-log_dir`,
and the C runtime recursively deletes that directory at the start of a fresh solve.
An escaping value must not be silently accepted.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from picurv_cli.core import (  # noqa: E402
    classify_run_directory_value,
    validate_run_directory_containment,
)


@pytest.mark.parametrize("value", ["logs", "output", "nested/logs", "./logs", "a/../b"])
def test_relative_paths_inside_the_run_are_contained(value):
    """!
    @brief Ordinary relative directory names must remain acceptable.
    @param[in] value Configured directory value.
    @return None.
    """
    assert classify_run_directory_value(value) == "contained"


# `~/logs` moved out of this list when it gained its own verdict: it is not an escape
# to somewhere, it is a value no layer agrees on. See the tilde tests below.
@pytest.mark.parametrize(
    "value", ["/tmp/shared_logs", "/", "../logs", "../../shared", "a/../../escape"]
)
def test_escaping_paths_are_detected(value):
    """!
    @brief Absolute paths and parent traversal must be classified as escaping.
    @param[in] value Configured directory value.
    @return None.
    """
    assert classify_run_directory_value(value) == "escaping"


@pytest.mark.parametrize("value", ["", "   ", None, 5, []])
def test_malformed_values_are_detected(value):
    """!
    @brief Non-string or empty directory values must be rejected as invalid.
    @param[in] value Configured directory value.
    @return None.
    """
    assert classify_run_directory_value(value) == "invalid"


def monitor(**dirs) -> dict:
    """!
    @brief Build a monitor config carrying the given io.directories mapping.
    @param[in] dirs Directory entries.
    @return Monitor configuration mapping.
    """
    return {"io": {"directories": dirs}}


def test_escaping_log_directory_is_an_error():
    """!
    @brief The destructive log path must be rejected, naming the consequence.
    @return None.
    """
    errors, warnings = validate_run_directory_containment(monitor(log="/tmp/victim"), "m.yml")
    assert len(errors) == 1 and not warnings
    assert "RECURSIVELY DELETES" in errors[0]


def test_escaping_output_directory_is_an_error():
    """!
    @brief Run output must stay inside the run tree.
    @return None.
    """
    errors, _ = validate_run_directory_containment(monitor(output="/tmp/elsewhere"), "m.yml")
    assert len(errors) == 1
    assert "archived and restored" in errors[0]


def test_default_directories_pass():
    """!
    @brief The shipped defaults must remain valid.
    @return None.
    """
    errors, warnings = validate_run_directory_containment(
        monitor(output="output", restart="restart", log="logs"), "m.yml"
    )
    assert errors == [] and warnings == []


def test_restart_directory_is_not_containment_checked():
    """!
    @brief `restart` is a read path and may legitimately point at another run.
    @return None.
    """
    errors, warnings = validate_run_directory_containment(monitor(restart="/other/run/output"), "m.yml")
    assert errors == [] and warnings == []


def test_override_downgrades_error_to_warning():
    """!
    @brief The explicit override must permit an escaping path, loudly.
    @return None.
    """
    errors, warnings = validate_run_directory_containment(
        monitor(log="/tmp/victim", allow_unsafe_paths=True), "m.yml"
    )
    assert errors == []
    assert warnings and any("Allowed only because" in w for w in warnings)


def test_override_must_be_explicitly_true():
    """!
    @brief A falsy override must not weaken the check.
    @return None.
    """
    errors, _ = validate_run_directory_containment(
        monitor(log="/tmp/victim", allow_unsafe_paths=False), "m.yml"
    )
    assert len(errors) == 1


def test_absent_directories_block_is_accepted():
    """!
    @brief A monitor without an io.directories block must not error.
    @return None.
    """
    errors, warnings = validate_run_directory_containment({"io": {}}, "m.yml")
    assert errors == [] and warnings == []


# --------------------------------------------------------------------------
# Bypasses found in review
# --------------------------------------------------------------------------

from picurv_cli.core import (  # noqa: E402
    RESERVED_DIRECTORY_FLAGS,
    append_passthrough_flags,
    paths_overlap,
    preflight_staged_run_directories,
    resolve_unsafe_paths_override,
    validate_reserved_directory_flags,
)


@pytest.mark.parametrize("raw", ["false", "true", "no", "0", 1, 0, "", [], {}, None])
def test_non_boolean_override_is_rejected(raw):
    """!
    @brief Only a real YAML boolean may enable the unsafe override.
    @param[in] raw Override value under test.
    @return None.
    """
    enabled, errors = resolve_unsafe_paths_override({"allow_unsafe_paths": raw}, "m.yml")
    assert enabled is False
    if raw is not False:
        assert errors, f"{raw!r} must be rejected rather than silently interpreted"


def test_boolean_override_values_behave():
    """!
    @brief Genuine booleans enable and disable the override without error.
    @return None.
    """
    assert resolve_unsafe_paths_override({"allow_unsafe_paths": True}, "m.yml") == (True, [])
    assert resolve_unsafe_paths_override({"allow_unsafe_paths": False}, "m.yml") == (False, [])


@pytest.mark.parametrize("value", [".", "a/..", "./.", "x/y/../.."])
def test_run_root_values_are_rejected(value):
    """!
    @brief A directory resolving to the run root would delete the run itself.
    @param[in] value Configured directory value.
    @return None.
    """
    assert classify_run_directory_value(value) == "run_root"
    errors, _ = validate_run_directory_containment(monitor(log=value), "m.yml")
    assert len(errors) == 1 and "run directory itself" in errors[0]


@pytest.mark.parametrize("reserved", ["config", "scheduler", "checkpoints", "visualization"])
def test_reserved_directory_collisions_are_rejected(reserved):
    """!
    @brief A run-owned directory must not collide with a reserved run-tree directory.
    @param[in] reserved Reserved directory name.
    @return None.
    """
    errors, _ = validate_run_directory_containment(monitor(log=reserved), "m.yml")
    assert len(errors) == 1 and "reserved run directory" in errors[0]


@pytest.mark.parametrize(
    "log,output",
    [("logs", "logs"), ("logs", "logs/sub"), ("logs/sub", "logs"), ("a/b", "a/b")],
)
def test_overlapping_log_and_output_are_rejected(log, output):
    """!
    @brief Overlapping log and output would delete solver output.
    @param[in] log Configured log directory.
    @param[in] output Configured output directory.
    @return None.
    """
    errors, _ = validate_run_directory_containment(monitor(log=log, output=output), "m.yml")
    assert any("overlap" in e for e in errors)


def test_non_overlapping_paths_are_accepted():
    """!
    @brief Distinct sibling directories remain valid.
    @return None.
    """
    errors, warnings = validate_run_directory_containment(monitor(log="logs", output="output"), "m.yml")
    assert errors == [] and warnings == []


def test_paths_overlap_helper():
    """!
    @brief Overlap detection must not fire on merely similar prefixes.
    @return None.
    """
    assert paths_overlap("logs", "logs")
    assert paths_overlap("logs", "logs/sub")
    assert not paths_overlap("logs", "logs_extra")


@pytest.mark.parametrize("flag", RESERVED_DIRECTORY_FLAGS)
def test_passthrough_cannot_set_reserved_flags(flag):
    """!
    @brief Raw passthrough must not be able to set a run-directory flag.
    @param[in] flag Reserved runtime flag.
    @return None.
    """
    config = {"petsc_passthrough_options": {flag: "/tmp/victim"}}
    problems = validate_reserved_directory_flags(config, "s.yml", "solver")
    assert len(problems) == 1 and flag in problems[0]


def test_passthrough_scan_reaches_nested_surfaces():
    """!
    @brief The scan must find reserved flags nested anywhere in the mapping.
    @return None.
    """
    config = {"solver_monitoring": {"petsc_passthrough_options": {"-log_dir": "/tmp/x"}}}
    assert validate_reserved_directory_flags(config, "m.yml", "monitor") != []


def test_emitter_refuses_reserved_flags():
    """!
    @brief The control-file emitter is the last line of defence.
    @return None.
    """
    with pytest.raises(ValueError, match="reserved directory flag"):
        append_passthrough_flags([], {"-log_dir": "/tmp/victim"})


def test_emitter_still_writes_ordinary_flags():
    """!
    @brief Ordinary passthrough options must be unaffected.
    @return None.
    """
    lines: list = []
    append_passthrough_flags(lines, {"-ksp_view": True, "-ksp_rtol": 1e-6})
    assert "-ksp_view" in lines
    assert any("-ksp_rtol" in line for line in lines)


def test_submission_preflight_rejects_escaping_control(tmp_path):
    """!
    @brief A staged run whose control file escapes must not be submittable.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    config = tmp_path / "config"
    config.mkdir()
    (config / "r.control").write_text(
        "-log_dir /tmp/victim\n-output_dir output\n-tio 100\n", encoding="utf-8"
    )
    errors, warnings = preflight_staged_run_directories(str(tmp_path))
    assert errors and warnings == []
    assert any("RECURSIVELY DELETES" in e for e in errors)


def test_submission_preflight_rejects_reserved_collision(tmp_path):
    """!
    @brief A staged control colliding with a reserved directory must be caught.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    config = tmp_path / "config"
    config.mkdir()
    (config / "r.control").write_text("-log_dir config\n-output_dir output\n", encoding="utf-8")
    errors, _ = preflight_staged_run_directories(str(tmp_path))
    # Two independent refusals: the value is not the canonical one, and it also names
    # a reserved directory the runtime would delete.
    assert len(errors) == 2
    assert any("the canonical value is 'logs'" in message for message in errors)
    assert any("reserved run directory" in message for message in errors)


def test_submission_preflight_accepts_contained_control(tmp_path):
    """!
    @brief A well-formed staged run must pass preflight.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    config = tmp_path / "config"
    config.mkdir()
    (config / "r.control").write_text("-log_dir logs\n-output_dir output\n", encoding="utf-8")
    assert preflight_staged_run_directories(str(tmp_path)) == ([], [])


# --------------------------------------------------------------------------
# Effective defaults and shared preflight logic
# --------------------------------------------------------------------------

from picurv_cli.core import (  # noqa: E402
    RUN_DIRECTORY_DEFAULTS,
    effective_run_directories,
    evaluate_run_directories,
    staged_control_directories,
)


def test_effective_directories_fill_in_defaults():
    """!
    @brief An omitted key takes its runtime default, not absence.
    @return None.
    """
    assert effective_run_directories({}) == RUN_DIRECTORY_DEFAULTS
    assert effective_run_directories({"log": "custom"})["output"] == RUN_DIRECTORY_DEFAULTS["output"]
    assert effective_run_directories({"log": "custom"})["log"] == "custom"


@pytest.mark.parametrize(
    "configured,label",
    [({"log": "output"}, "log collides with default output"),
     ({"output": "logs"}, "output collides with default log")],
)
def test_one_sided_default_collisions_are_rejected(configured, label):
    """!
    @brief A collision against a defaulted directory is as destructive as an explicit one.
    @param[in] configured Partial directory mapping.
    @param[in] label Human-readable case description.
    @return None.
    """
    errors, _ = validate_run_directory_containment(monitor(**configured), "m.yml")
    assert len(errors) == 1, label
    assert "overlap" in errors[0] and "(default)" in errors[0]


def test_explicit_defaults_still_pass():
    """!
    @brief Writing the defaults out explicitly must remain valid.
    @return None.
    """
    errors, warnings = validate_run_directory_containment(
        monitor(log="logs", output="output"), "m.yml"
    )
    assert errors == [] and warnings == []


def test_shared_evaluation_is_used_by_both_callers():
    """!
    @brief Config validation and preflight must agree, because they share one rule set.
    @return None.
    """
    effective = effective_run_directories({"log": "output"})
    shared_errors, _ = evaluate_run_directories(effective, override=False)
    config_errors, _ = validate_run_directory_containment(monitor(log="output"), "m.yml")
    assert len(shared_errors) == len(config_errors) == 1


def test_staged_control_parsing():
    """!
    @brief Directory flags must be read out of a staged control file.
    @return None.
    """
    import tempfile

    path = Path(tempfile.mkdtemp()) / "r.control"
    path.write_text("-tio 100\n-log_dir custom_logs\n-output_dir custom_out\n", encoding="utf-8")
    values, parse_errors = staged_control_directories(str(path))
    assert values == {"log": "custom_logs", "output": "custom_out"}
    assert parse_errors == []


def staged(tmp_path: Path, control: str, monitor_yaml: str = None) -> str:
    """!
    @brief Build a staged run directory for preflight tests.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] control Contents of the staged control file.
    @param[in] monitor_yaml Optional staged monitor YAML.
    @return Path to the staged run directory.
    """
    config = tmp_path / "config"
    config.mkdir(exist_ok=True)
    (config / "r.control").write_text(control, encoding="utf-8")
    if monitor_yaml is not None:
        (config / "monitor.yml").write_text(monitor_yaml, encoding="utf-8")
    return str(tmp_path)


def test_preflight_rejects_explicit_staged_overlap(tmp_path):
    """!
    @brief A staged control with log and output equal must be refused.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, _ = preflight_staged_run_directories(
        staged(tmp_path, "-log_dir output\n-output_dir output\n")
    )
    assert len(errors) == 2
    assert any("the canonical value is 'logs'" in message for message in errors)
    assert any("overlap" in message for message in errors)


def test_preflight_rejects_overlap_against_missing_flag(tmp_path):
    """!
    @brief A staged control omitting one flag still collides through its default.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, _ = preflight_staged_run_directories(staged(tmp_path, "-log_dir output\n"))
    assert len(errors) == 2
    assert any("the canonical value is 'logs'" in message for message in errors)
    assert any("overlap" in message and "(default)" in message for message in errors)


def test_preflight_refuses_staged_unsafe_override_under_fixed_topology(tmp_path):
    """!
    @brief The legacy unsafe-path override no longer waives a non-canonical staged path.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, warnings = preflight_staged_run_directories(
        staged(
            tmp_path,
            "-log_dir /tmp/victim\n-output_dir output\n",
            "io:\n  directories:\n    allow_unsafe_paths: true\n",
        )
    )
    # Run-owned directories are fixed, so there is nothing left to waive: the escape is
    # refused outright rather than downgraded to a warning by a staged override.
    assert warnings == []
    assert any("the canonical value is 'logs'" in message for message in errors)
    assert any("escapes the run directory" in message for message in errors)


def test_preflight_ignores_non_boolean_staged_override(tmp_path):
    """!
    @brief A quoted override in a staged monitor must not grant an unsafe path.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, _ = preflight_staged_run_directories(
        staged(
            tmp_path,
            "-log_dir /tmp/victim\n-output_dir output\n",
            'io:\n  directories:\n    allow_unsafe_paths: "true"\n',
        )
    )
    assert errors


def test_preflight_without_directory_flags_is_silent(tmp_path):
    """!
    @brief A control file carrying no directory flags must not be judged.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    assert preflight_staged_run_directories(staged(tmp_path, "-tio 100\n")) == ([], [])


# --------------------------------------------------------------------------
# PETSc options-file parsing contract
# --------------------------------------------------------------------------

from picurv_cli.core import directory_value_charset_problem  # noqa: E402


def test_quoted_absolute_path_with_spaces_is_parsed_as_one_token(tmp_path):
    """!
    @brief PETSc treats a quoted span as one token; naive splitting hid an absolute path.

    Splitting on whitespace read `-log_dir "/tmp/VICTIM DIR"` as `"/tmp/VICTIM`, which
    looks relative and contained while PETSc would use the absolute path.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    control = tmp_path / "r.control"
    control.write_text('-log_dir "/tmp/VICTIM DIR"\n-output_dir output\n', encoding="utf-8")
    values, parse_errors = staged_control_directories(str(control))
    assert values["log"] == "/tmp/VICTIM DIR"
    assert parse_errors == []
    assert classify_run_directory_value(values["log"]) == "escaping"


def test_preflight_rejects_quoted_escaping_path(tmp_path):
    """!
    @brief The quoted-path bypass must be refused at submission, as escaping.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, _ = preflight_staged_run_directories(
        staged(tmp_path, '-log_dir "/tmp/VICTIM DIR"\n-output_dir output\n')
    )
    assert any("escapes the run directory" in e for e in errors)


def test_preflight_treats_malformed_quoting_as_an_error(tmp_path):
    """!
    @brief An uninterpretable line must fail preflight rather than be skipped.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, _ = preflight_staged_run_directories(staged(tmp_path, '-log_dir "/tmp/unclosed\n'))
    assert len(errors) == 1 and "malformed quoting" in errors[0]


def test_preflight_honours_comment_lines(tmp_path):
    """!
    @brief A commented-out directory flag must not be read as configuration.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    errors, _ = preflight_staged_run_directories(
        staged(tmp_path, "# -log_dir /tmp/victim\n-log_dir logs\n-output_dir output\n")
    )
    assert errors == []


@pytest.mark.parametrize(
    "value,reason",
    [("my logs", "whitespace"), ('a"b', "double quote"), ("a'b", "single quote"),
     ("x#y", "comment marker"), ("a\tb", "whitespace")],
)
def test_unsafe_directory_characters_are_rejected(value, reason):
    """!
    @brief Directory values must be writable to an options line without quoting.
    @param[in] value Configured directory value.
    @param[in] reason Expected reason fragment.
    @return None.
    """
    assert reason in directory_value_charset_problem(value)
    errors, _ = validate_run_directory_containment(monitor(log=value), "m.yml")
    assert len(errors) == 1 and reason in errors[0]


@pytest.mark.parametrize("value", ["logs", "diagnostics/run1", "a_b-c.d"])
def test_ordinary_directory_names_remain_valid(value):
    """!
    @brief Plain relative names must not be caught by the charset contract.
    @param[in] value Configured directory value.
    @return None.
    """
    assert directory_value_charset_problem(value) == ""


def test_escaping_is_reported_before_formatting():
    """!
    @brief A dangerous absolute path must be reported as escaping, not as a format issue.
    @return None.
    """
    errors, _ = validate_run_directory_containment(monitor(log="/tmp/VICTIM DIR"), "m.yml")
    assert any("escapes the run directory" in e for e in errors)
    assert any("whitespace" in e for e in errors), "a non-waivable character problem must also report"


# --------------------------------------------------------------------------
# Bypasses found in the readiness review
# --------------------------------------------------------------------------

import os  # noqa: E402
import subprocess  # noqa: E402

from picurv_cli.core import (  # noqa: E402
    RESERVED_INDIRECTION_FLAGS,
    check_physical_containment,
    preflight_config_directories,
)


@pytest.mark.parametrize(
    "configured", [{"log": "."}, {"log": "config"}, {"log": "output"}, {"log": "a b"}]
)
def test_override_cannot_waive_self_destruction_or_malformed_values(configured):
    """!
    @brief The override covers deliberate external storage, never self-destruction.
    @param[in] configured Directory mapping under test.
    @return None.
    """
    from picurv_cli.core import effective_run_directories, evaluate_run_directories

    errors, _ = evaluate_run_directories(
        effective_run_directories(configured), override=True, explicit=set(configured)
    )
    assert errors, f"{configured} must not be waivable"
    assert any("cannot be overridden" in e for e in errors)


def test_override_still_waives_a_well_formed_external_path():
    """!
    @brief A deliberate external location remains permitted, loudly.
    @return None.
    """
    from picurv_cli.core import effective_run_directories, evaluate_run_directories

    errors, warnings = evaluate_run_directories(
        effective_run_directories({"log": "/central/picurv/logs"}), override=True, explicit={"log"}
    )
    assert errors == [] and len(warnings) == 1


def test_escaping_path_still_gets_character_validation():
    """!
    @brief Character rules apply to escaping values too; they were previously skipped.
    @return None.
    """
    from picurv_cli.core import effective_run_directories, evaluate_run_directories

    errors, _ = evaluate_run_directories(
        effective_run_directories({"log": '/tmp/a b"c'}), override=False, explicit={"log"}
    )
    assert any("double quote" in e or "whitespace" in e for e in errors)


def test_symlinked_directory_is_physically_detected(tmp_path):
    """!
    @brief A lexically contained name may still be a symlink out of the run tree.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "run"
    run.mkdir()
    victim = tmp_path / "victim"
    victim.mkdir()
    os.symlink(str(victim), str(run / "sneaky"))
    problems = check_physical_containment(str(run), {"log": "sneaky", "output": "output"})
    assert len(problems) == 1 and "outside the run directory" in problems[0]


def test_physical_containment_accepts_real_subdirectories(tmp_path):
    """!
    @brief An ordinary subdirectory must not be flagged.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "run"
    (run / "logs").mkdir(parents=True)
    assert check_physical_containment(str(run), {"log": "logs", "output": "output"}) == []


def test_study_members_are_preflighted(tmp_path):
    """!
    @brief Study controls live under cases/<member>/config, not at the study root.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    for member in ("a", "b"):
        config = tmp_path / "cases" / member / "config"
        config.mkdir(parents=True)
        (config / "r.control").write_text("-log_dir /tmp/victim\n-output_dir output\n", encoding="utf-8")
    assert len(preflight_config_directories(str(tmp_path))) == 2
    errors, _ = preflight_staged_run_directories(str(tmp_path))
    assert all(any(f"cases/{m}/config" in e for e in errors) for m in ("a", "b")), (
        "every study member must be checked"
    )


@pytest.mark.parametrize("flag", RESERVED_INDIRECTION_FLAGS)
def test_petsc_indirection_flags_are_rejected(flag):
    """!
    @brief An options file or alias could reintroduce a directory flag unchecked.
    @param[in] flag Indirection flag under test.
    @return None.
    """
    problems = validate_reserved_directory_flags(
        {"petsc_passthrough_options": {flag: "extra.opts"}}, "s.yml", "solver"
    )
    assert len(problems) == 1 and "indirection" in problems[0].lower() or "evaluates that" in problems[0]


def test_c_guard_rejects_unsafe_log_directories():
    """!
    @brief The C-side final guard must reject the same cases as the launcher.

    This is the last check before an irreversible recursive delete, and it is
    deliberately independent of anything Python validated.
    @return None.
    """
    source = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    assert "LogDirectoryIsSafeToWipe" in source
    # Structural, and only structural: the guard must stand between the decision and
    # the delete. What the guard decides is asserted behaviourally against the built
    # binary in tests/test_c_guard_behaviour.py, because a source scrape passes just as
    # happily on dead code.
    guard_at = source.index("LogDirectoryIsSafeToWipe(simCtx->log_dir")
    delete_at = source.index("PetscRMTree(simCtx->log_dir)")
    assert guard_at < delete_at, "the guard must run before the recursive delete"
    assert "Refusing to delete log directory" in source
    for helper in ("ClassifyLogDirectory", "ResolveDirectoryPhysically", "realpath", "getcwd"):
        assert helper in source, f"C guard does not use {helper}"


def test_validated_directory_flags_are_emitted_last_and_always(tmp_path):
    """!
    @brief The generated control must end with validated directories.

    PETSc takes the last occurrence, and emitting unconditionally removes the gap where
    an omitted flag let external configuration choose the directory.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    source = (REPO_ROOT / "picurv_cli" / "core.py").read_text(encoding="utf-8")
    marker = "# Canonical run-owned directories are fixed by the workspace contract."
    assert marker in source
    tail = source[source.index(marker) :]
    assert "-output_dir" in tail and "-log_dir" in tail


# --------------------------------------------------------------------------
# Python <-> C override reconciliation
# --------------------------------------------------------------------------

SIMULATOR = REPO_ROOT / "bin" / "simulator"


def c_guard_source() -> str:
    """!
    @brief The body of the C-side log-directory guard.
    @return Source text of the guard function.
    """
    source = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    start = source.index("static PetscBool LogDirectoryIsSafeToWipe")
    return source[start : source.index("\n}\n", start) + 3]


def test_authorization_flag_is_reserved_against_forgery():
    """!
    @brief Raw passthrough must not be able to forge the runtime authorization.
    @return None.
    """
    from picurv_cli.core import RESERVED_DIRECTORY_FLAGS

    assert "-allow_unsafe_log_dir" in RESERVED_DIRECTORY_FLAGS
    problems = validate_reserved_directory_flags(
        {"petsc_passthrough_options": {"-allow_unsafe_log_dir": "true"}}, "s.yml", "solver"
    )
    assert problems


def test_c_guard_waives_only_the_external_restriction():
    """!
    @brief The runtime authorization must not waive self-destructive or ambiguous values.

    Python and C must agree on exactly what the override covers, or a configuration can
    pass validation and then fail at runtime - or worse, pass both when it should not.
    @return None.
    """
    source = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    classifier = source[source.index("static DirectoryVerdict ClassifyLogDirectory"):]
    classifier = classifier[: classifier.index("\n}\n") + 3]
    # Classification is total and mentions no authorization at all: the waiver is
    # applied once, by the caller, to exactly one verdict. That is what makes "checked
    # before authorization is honoured" true by construction rather than by ordering.
    signature = classifier[: classifier.index(")")]
    assert "authorized" not in signature, (
        "ClassifyLogDirectory must not take the authorization; the caller applies it"
    )
    assert "if (authorized" not in classifier, (
        "ClassifyLogDirectory must not branch on the authorization"
    )
    for non_waivable in ("DirectoryValueIsWellFormed", "DirectoryHitsReservedName",
                         "DirectoriesOverlap", "DIR_VERDICT_RUN_ROOT",
                         "DIR_VERDICT_ANCESTOR", "DIR_VERDICT_RELATIVE_ESCAPE"):
        assert non_waivable in classifier, f"C guard lacks {non_waivable}"
    guard = c_guard_source()
    assert guard.count("DIR_VERDICT_EXTERNAL_ABSOLUTE && authorized") == 1, (
        "exactly one verdict may be waived by the authorization"
    )


def test_control_file_never_carries_obsolete_directory_authorization():
    """!
    @brief Canonical fixed paths make the former external-directory waiver unnecessary.
    @return None.
    """
    source = (REPO_ROOT / "picurv_cli" / "core.py").read_text(encoding="utf-8")
    assert "-allow_unsafe_log_dir true" not in source
    assert "-output_dir {CANONICAL_RUN_PATHS['output']}" in source
    assert "-log_dir {CANONICAL_RUN_PATHS['logs']}" in source


PICURV_CLI = REPO_ROOT / "picurv_cli" / "picurv"


@pytest.fixture(scope="module")
def staged_control(tmp_path_factory):
    """!
    @brief A real generated control file, staged once for the guard tests.

    @details The guard runs after the solver has loaded its configuration, so a
             hand-written stub fails earlier and never reaches it. Staging a genuine run
             is the only way to exercise the guard as a user would meet it.
    @param[in] tmp_path_factory Pytest temporary directory factory.
    @return Tuple of (working directory, control file text).
    """
    workspace = tmp_path_factory.mktemp("guard")
    case = workspace / "c"
    init = subprocess.run(
        [sys.executable, str(PICURV_CLI), "init", "flat_channel", "--dest", str(case)],
        cwd=str(workspace), capture_output=True, text=True, check=False, timeout=180,
    )
    if init.returncode != 0:
        pytest.skip("picurv init unavailable")
    staged = subprocess.run(
        [sys.executable, str(PICURV_CLI), "run", "--solve", "--no-submit",
         "--cluster", str(case / "slurm_cluster.yml"),
         "--case", str(case / "quickstart_flat_channel.yml"),
         "--solver", str(case / "Imp-MG-Standard.yml"),
         "--monitor", str(case / "quickstart_Standard_Output.yml"),
         "--post", str(case / "quickstart_standard_analysis.yml")],
        cwd=str(workspace), capture_output=True, text=True, check=False, timeout=300,
    )
    controls = sorted(workspace.glob("runs/*/config/*.control"))
    if staged.returncode != 0 or not controls:
        pytest.skip("could not stage a run for the guard tests")
    run_dir = controls[0].parent.parent
    return run_dir, controls[0].read_text(encoding="utf-8")


@pytest.mark.skipif(not SIMULATOR.is_file(), reason="solver binary not built")
@pytest.mark.parametrize(
    "log_dir,authorized,expect_refusal",
    [
        ("logs", False, False),
        ("/tmp/picurv-external-guard", False, True),
        ("/tmp/picurv-external-guard", True, False),
        ("config", True, True),
        (".", True, True),
    ],
)
def test_c_guard_end_to_end(staged_control, log_dir, authorized, expect_refusal):
    """!
    @brief Drive the real compiled guard for each documented case.

    @details Exercises the binary rather than reasoning about the source, so a guard
             weakened to make tests pass would be caught here.
    @param[in] staged_control Staged run directory and control text.
    @param[in] log_dir Configured log directory.
    @param[in] authorized Whether the runtime authorization is supplied.
    @param[in] expect_refusal Whether the guard must refuse.
    @return None.
    """
    run_dir, text = staged_control
    lines = [
        line for line in text.splitlines()
        if not line.startswith(("-log_dir", "-allow_unsafe_log_dir"))
    ]
    lines.append(f"-log_dir {log_dir}")
    if authorized:
        lines.append("-allow_unsafe_log_dir true")
    control = run_dir / "config" / "guard_probe.control"
    control.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = subprocess.run(
        [str(SIMULATOR), "-control_file", str(control)],
        cwd=str(run_dir), capture_output=True, text=True, check=False, timeout=180,
    )
    refused = "Refusing to delete log directory" in (result.stdout + result.stderr)
    assert refused is expect_refusal, (
        f"log_dir={log_dir!r} authorized={authorized}: expected refusal={expect_refusal}, got {refused}"
    )


# --------------------------------------------------------------------------
# Bypasses found by the independent containment review
# --------------------------------------------------------------------------

from picurv_cli.core import control_value  # noqa: E402


@pytest.mark.parametrize(
    "value",
    ["..", "../sibling", "../case_0000/config", "../case_0000/output",
     "./config", ".//config", "output/", "output/sub", "/abs/run/config", "a/b/config"],
)
def test_override_cannot_authorize_self_destruction_in_any_spelling(value):
    """!
    @brief The override must not waive a value that targets the run's own tree.

    Escaping values previously skipped the non-waivable checks entirely, so
    `log: ".."` with the override deleted the whole runs/ tree.
    @param[in] value Configured log directory.
    @return None.
    """
    from picurv_cli.core import effective_run_directories, evaluate_run_directories

    errors, _ = evaluate_run_directories(
        effective_run_directories({"log": value}), override=True, explicit={"log"}
    )
    assert errors, f"{value!r} must not be authorizable"


def test_relative_escape_is_never_authorizable():
    """!
    @brief A relative escape lands among sibling runs, so it is never permitted.
    @return None.
    """
    from picurv_cli.core import effective_run_directories, evaluate_run_directories

    errors, _ = evaluate_run_directories(
        effective_run_directories({"log": "../elsewhere"}), override=True, explicit={"log"}
    )
    assert any("relative traversal" in e for e in errors)


def test_absolute_external_path_remains_authorizable():
    """!
    @brief The authorized use case must keep working.
    @return None.
    """
    from picurv_cli.core import effective_run_directories, evaluate_run_directories

    errors, warnings = evaluate_run_directories(
        effective_run_directories({"log": "/central/picurv/logs"}), override=True, explicit={"log"}
    )
    assert errors == [] and warnings


def test_physical_containment_rejects_the_run_root(tmp_path):
    """!
    @brief A path resolving to the run root itself must be refused.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "run"
    run.mkdir()
    problems = check_physical_containment(str(run), {"log": ".", "output": "output"})
    assert any("run directory itself" in p for p in problems)


@pytest.mark.parametrize("flag", ["-options_file", "-options_file_yaml", "-alias"])
def test_staged_control_indirection_is_refused(tmp_path, flag):
    """!
    @brief A staged control must not smuggle directories through PETSc indirection.
    @param[in] tmp_path Pytest temporary directory.
    @param[in] flag Indirection flag.
    @return None.
    """
    control = tmp_path / "r.control"
    control.write_text(f"{flag} /tmp/evil.opts\n-log_dir logs\n", encoding="utf-8")
    _, parse_errors = staged_control_directories(str(control))
    assert parse_errors and flag in parse_errors[0]


@pytest.mark.parametrize("value", ["a\nb", "a\r\nb", "x\n-allow_unsafe_log_dir true"])
def test_newline_values_cannot_inject_option_lines(value):
    """!
    @brief A value with a newline would write extra option lines into the control file.
    @param[in] value Value under test.
    @return None.
    """
    with pytest.raises(SystemExit):
        control_value(value, "test.context")


def test_ordinary_control_values_pass_through():
    """!
    @brief The injection guard must not disturb ordinary values.
    @return None.
    """
    assert control_value("restart", "test.context") == "restart"
    assert control_value(1.5, "test.context") == 1.5


@pytest.mark.parametrize(
    "value,reserved",
    [("./config", True), (".//config", True), ("a/b/config", True), ("output/sub", True),
     ("logs", False), ("diagnostics/run1", False)],
)
def test_c_guard_reserved_walk_matches_every_segment(value, reserved):
    """!
    @brief The C reserved-name walk must inspect every segment, not just the first.

    A first-segment-only check let `./config` and `output/sub` delete run-owned
    directories with no authorization at all.
    @param[in] value Directory value.
    @param[in] reserved Whether it targets a reserved directory.
    @return None.
    """
    guard = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    start = guard.index("static PetscBool DirectoryHitsReservedName")
    body = guard[start : guard.index("\n}\n", start)]
    assert "while (*cursor)" in body, "reserved walk must iterate over all segments"
    assert "strcspn" in body


def test_c_guard_overlap_is_normalized():
    """!
    @brief Overlap must be compared on normalized segments, not with a bare strcmp.
    @return None.
    """
    guard = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    assert "DirectoriesOverlap" in guard
    start = guard.index("static PetscBool DirectoriesOverlap")
    body = guard[start : guard.index("\n}\n", start)]
    assert "strcspn" in body, "overlap must normalize before comparing"


def test_c_guard_rejects_unresolvable_parent_reference():
    """!
    @brief An unresolvable path containing `..` must not be assumed safe.
    @return None.
    """
    source = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    # A path that does not exist yet is resolved through its longest existing prefix,
    # so ".." is resolved rather than assumed safe. The behavioural counterpart lives
    # in tests/test_c_guard_behaviour.py.
    assert "ResolveDirectoryPhysically" in source
    assert "NormalizePathLexically" in source


# ---------------------------------------------------------------------------
# Typed physical verdicts.
#
# The earlier rule decided waivability by looking for the word "escapes" in the
# message. A finding reading "resolves to the run directory itself" contains no such
# word, so self-destruction was downgraded to a warning. These tests pin the rule to
# the verdict instead. None of them delete anything; the filesystem root is covered by
# classification only.
# ---------------------------------------------------------------------------

from picurv_cli.core import (  # noqa: E402
    PHYSICAL_VERDICT_ANCESTOR,
    PHYSICAL_VERDICT_EXTERNAL_ABSOLUTE,
    PHYSICAL_VERDICT_RELATIVE_ESCAPE,
    PHYSICAL_VERDICT_RUN_ROOT,
    WAIVABLE_PHYSICAL_VERDICTS,
    classify_physical_containment,
)


def _verdicts(run_dir, values):
    """!
    @brief Verdicts keyed by directory name.
    @param[in] run_dir Run directory.
    @param[in] values Effective directory mapping.
    @return Mapping of key to verdict.
    """
    return {key: verdict for key, verdict, _ in classify_physical_containment(run_dir, values)}


def test_only_an_external_absolute_location_is_ever_waivable():
    """!
    @brief Exactly one verdict may be waived by an authorization.
    @return None.
    """
    assert WAIVABLE_PHYSICAL_VERDICTS == {PHYSICAL_VERDICT_EXTERNAL_ABSOLUTE}


def test_absolute_path_equal_to_the_run_root_is_not_waivable(tmp_path):
    """!
    @brief An authorized absolute path naming the run itself must stay fatal.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "run"
    run.mkdir()
    verdict = _verdicts(str(run), {"log": str(run), "output": "output"})["log"]
    assert verdict == PHYSICAL_VERDICT_RUN_ROOT
    assert verdict not in WAIVABLE_PHYSICAL_VERDICTS


def test_absolute_ancestor_of_the_run_root_is_not_waivable(tmp_path):
    """!
    @brief An authorized absolute ancestor would destroy the run and its siblings.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "workspace" / "runs" / "case_0001"
    run.mkdir(parents=True)
    for ancestor in (run.parent, run.parent.parent, tmp_path):
        verdict = _verdicts(str(run), {"log": str(ancestor), "output": "output"})["log"]
        assert verdict == PHYSICAL_VERDICT_ANCESTOR, ancestor
        assert verdict not in WAIVABLE_PHYSICAL_VERDICTS


def test_filesystem_root_is_classified_as_an_ancestor(tmp_path):
    """!
    @brief The filesystem root contains every run. Classification only - nothing is run.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "run"
    run.mkdir()
    verdict = _verdicts(str(run), {"log": "/", "output": "output"})["log"]
    assert verdict == PHYSICAL_VERDICT_ANCESTOR
    assert verdict not in WAIVABLE_PHYSICAL_VERDICTS


def test_relative_symlink_escape_is_not_waivable(tmp_path):
    """!
    @brief A relative name resolving out through a symlink is never authorizable.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "run"
    run.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (run / "sneaky").symlink_to(outside)
    verdict = _verdicts(str(run), {"log": "sneaky", "output": "output"})["log"]
    assert verdict == PHYSICAL_VERDICT_RELATIVE_ESCAPE
    assert verdict not in WAIVABLE_PHYSICAL_VERDICTS


def test_absolute_symlink_to_an_external_location_stays_waivable(tmp_path):
    """!
    @brief A deliberate external absolute location remains an authorizable choice.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "workspace" / "runs" / "case"
    run.mkdir(parents=True)
    outside = tmp_path / "archive"
    outside.mkdir()
    verdict = _verdicts(str(run), {"log": str(outside), "output": "output"})["log"]
    assert verdict == PHYSICAL_VERDICT_EXTERNAL_ABSOLUTE
    assert verdict in WAIVABLE_PHYSICAL_VERDICTS


def test_a_contained_directory_produces_no_finding(tmp_path):
    """!
    @brief The ordinary case stays silent.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "run"
    (run / "logs").mkdir(parents=True)
    assert classify_physical_containment(str(run), {"log": "logs", "output": "output"}) == []


def test_no_message_matching_decides_waivability(tmp_path):
    """!
    @brief No non-waivable message may contain the token the old rule keyed on.

    @details The bug was structural: `"escapes" not in message` silently waived any
             finding phrased without that word. Pinning the phrasing would recreate it,
             so this asserts the opposite - the rule must not be recoverable from text.
    @param[in] tmp_path Temporary directory fixture.
    @return None.
    """
    run = tmp_path / "workspace" / "runs" / "case"
    run.mkdir(parents=True)
    outside = tmp_path / "outside"
    outside.mkdir()
    (run / "sneaky").symlink_to(outside)
    cases = [str(run), str(run.parent), "/", "sneaky"]
    for value in cases:
        for _, verdict, message in classify_physical_containment(
                str(run), {"log": value, "output": "output"}):
            if verdict in WAIVABLE_PHYSICAL_VERDICTS:
                continue
            assert "This cannot be overridden." in message, (value, message)


def test_staged_control_naming_the_run_root_is_refused_at_submission(tmp_path):
    """!
    @brief Submission preflight must refuse a staged control naming the run itself.

    @details The preflight previously decided waivability by looking for the word
             "escapes" in the message. This finding does not contain it, so an
             authorized staged run could pass preflight and then delete itself.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "runs" / "case_0001"
    (run / "config").mkdir(parents=True)
    control = run / "config" / "case_0001.control"
    control.write_text(
        f"-output_dir output\n-allow_unsafe_log_dir true\n-log_dir {run.resolve()}\n",
        encoding="utf-8",
    )
    errors, _ = preflight_staged_run_directories(str(run))
    assert errors, "an authorized staged control naming the run root must be refused"
    assert any("cannot be overridden" in message for message in errors)


def test_staged_control_naming_an_ancestor_is_refused_at_submission(tmp_path):
    """!
    @brief A staged control naming an ancestor of the run must be refused.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "workspace" / "runs" / "case_0001"
    (run / "config").mkdir(parents=True)
    control = run / "config" / "case_0001.control"
    control.write_text(
        f"-output_dir output\n-allow_unsafe_log_dir true\n"
        f"-log_dir {run.parent.parent.resolve()}\n",
        encoding="utf-8",
    )
    errors, _ = preflight_staged_run_directories(str(run))
    assert errors
    assert any("CONTAINS the run directory" in message for message in errors)


def test_staged_external_absolute_is_refused_at_submission(tmp_path):
    """!
    @brief The formerly waivable external absolute path is now refused at submission.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "runs" / "case_0001"
    (run / "config").mkdir(parents=True)
    external = tmp_path / "archive"
    external.mkdir()
    # Preflight reads the override from the staged monitor, the same source the
    # configuration layer used when the run was staged.
    (run / "config" / "monitor.yml").write_text(
        "io:\n  directories:\n    allow_unsafe_paths: true\n"
        f"    output: \"output\"\n    log: \"{external.resolve()}\"\n",
        encoding="utf-8",
    )
    control = run / "config" / "case_0001.control"
    control.write_text(
        f"-output_dir output\n-allow_unsafe_log_dir true\n-log_dir {external.resolve()}\n",
        encoding="utf-8",
    )
    errors, warnings = preflight_staged_run_directories(str(run))
    assert warnings == []
    assert any("the canonical value is 'logs'" in message for message in errors)
    assert any("outside the run directory" in message for message in errors)


def test_the_guard_does_not_spend_the_stack_on_path_buffers():
    """!
    @brief The guard's path buffers must be heap-allocated.

    @details They were stack arrays first: three PATH_MAX buffers in the classifier and
             five more in the physical resolution, about 32 KB inside one nested call.
             That overflowed far enough to corrupt the simulation context on rank zero,
             which diverged the ranks and failed `unit-newton-krylov` at two ranks with
             a collective mismatch - the guard broke the run it existed to protect. The
             behavioural regression is that suite; this pins the structural cause.
    @return None.
    """
    source = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    start = source.index("static DirectoryVerdict ClassifyLogDirectory")
    end = source.index("\n}\n", source.index("PetscFree(scratch)", start))
    classifier = source[start:end]
    assert "PetscMalloc1" in classifier and "PetscFree(scratch)" in classifier
    for helper in ("ClassifyLogDirectory", "ResolveDirectoryPhysically"):
        body_start = source.index(f"static {'DirectoryVerdict' if helper == 'ClassifyLogDirectory' else 'PetscBool'} {helper}")
        body = source[body_start: source.index("\n}\n", body_start)]
        declarations = [line for line in body.splitlines()
                        if "[PETSC_MAX_PATH_LEN]" in line and line.strip().startswith("char")]
        assert not declarations, (
            f"{helper} declares stack path buffers {declarations}; allocate them instead"
        )


# ---------------------------------------------------------------------------
# `~` is not an absolute path.
#
# Nothing expands it. The control file is read by PETSc, not by a shell, and every
# layer that uses the value resolves anything not starting with "/" relative to the
# run, so `~/logs` names a literal "~" directory inside the run tree.
#
# The defect was a disagreement, not a shell expansion: physical containment used to
# expand `~` and classify the result as an authorizable external absolute location,
# while nothing downstream expanded anything. The override then let it through. `~` is
# now its own non-waivable verdict at every layer, and no layer expands it.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", ["~", "~/logs", "~root/logs", "~/"])
def test_tilde_is_classified_as_its_own_verdict(value):
    """!
    @brief A `~` value is neither contained nor an ordinary escape.
    @param[in] value Configured directory value.
    @return None.
    """
    assert classify_run_directory_value(value) == "tilde"


@pytest.mark.parametrize("value", ["~/logs", "~", "~root/archive"])
def test_tilde_is_fatal_even_with_the_override(value):
    """!
    @brief The unsafe-paths override must not waive an unexpanded `~`.
    @param[in] value Configured directory value.
    @return None.
    """
    errors, warnings = evaluate_run_directories(
        {"log": value, "output": "output"}, True, {"log"}
    )
    assert errors, f"{value} was accepted with the override set"
    assert any("nothing expands" in message for message in errors)
    assert not any(value in message for message in warnings)


def test_tilde_is_not_expanded_by_the_physical_check(tmp_path):
    """!
    @brief Physical classification must not expand `~`.

    @details An earlier version did, which is what created the mismatch: this layer
             reasoned about an expanded home directory while every layer that used the
             value treated it literally.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "run"
    run.mkdir()
    findings = classify_physical_containment(str(run), {"log": "~/logs", "output": "output"})
    for _, verdict, _ in findings:
        assert verdict not in WAIVABLE_PHYSICAL_VERDICTS


def test_a_real_absolute_path_is_still_waivable(tmp_path):
    """!
    @brief Rejecting `~` must not reject the external locations that are legitimate.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    external = tmp_path / "archive"
    external.mkdir()
    errors, warnings = evaluate_run_directories(
        {"log": str(external), "output": "output"}, True, {"log"}
    )
    assert errors == []
    assert warnings


def test_staged_control_with_a_tilde_is_refused_at_submission(tmp_path):
    """!
    @brief Submission preflight must refuse a staged `~` log directory.
    @param[in] tmp_path Pytest temporary directory.
    @return None.
    """
    run = tmp_path / "runs" / "case_0001"
    (run / "config").mkdir(parents=True)
    (run / "config" / "monitor.yml").write_text(
        "io:\n  directories:\n    allow_unsafe_paths: true\n"
        '    output: "output"\n    log: "~/logs"\n', encoding="utf-8")
    (run / "config" / "case_0001.control").write_text(
        "-output_dir output\n-allow_unsafe_log_dir true\n-log_dir ~/logs\n", encoding="utf-8")
    errors, _ = preflight_staged_run_directories(str(run))
    assert errors
    assert any("nothing expands" in message for message in errors)


def test_no_layer_expands_a_tilde():
    """!
    @brief Neither the emitter nor the physical check may call expanduser on a value.

    @details The fix is that `~` is refused, not that one layer expands it. An
             expansion anywhere would recreate the mismatch in the opposite direction.
    @return None.
    """
    source = (REPO_ROOT / "picurv_cli" / "core.py").read_text(encoding="utf-8")
    start = source.index("def classify_physical_containment")
    body = source[start: source.index("\ndef ", start + 1)]
    assert "expanduser" not in body


def test_c_guard_gives_the_tilde_its_own_verdict():
    """!
    @brief The C classifier must refuse `~` before treating anything as absolute.
    @return None.
    """
    source = (REPO_ROOT / "src" / "setup.c").read_text(encoding="utf-8")
    assert "DIR_VERDICT_UNEXPANDED_TILDE" in source
    classifier = source[source.index("static DirectoryVerdict ClassifyLogDirectory"):]
    classifier = classifier[: classifier.index("\n}\n") + 3]
    tilde_at = classifier.index("log_dir[0] == '~'")
    absolute_at = classifier.index("log_dir[0] == '/'", tilde_at)
    assert tilde_at < absolute_at, "`~` must be refused before the absolute-path branch"
    guard = c_guard_source()
    assert "DIR_VERDICT_UNEXPANDED_TILDE" not in guard, (
        "the tilde verdict must not be waivable"
    )
