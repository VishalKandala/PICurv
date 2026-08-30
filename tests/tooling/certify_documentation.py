#!/usr/bin/env python3
"""Create a commit-scoped documentation certification after strict validation."""

from __future__ import annotations

import argparse
import json
import datetime as dt
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def run(command: list[str], label: str) -> None:
    """!
    @brief Run one certification gate and stream its output.
    @param[in] command Executable and arguments for the gate.
    @param[in] label Human-readable gate label.
    @return None.
    """

    print(f"==> {label}", flush=True)
    completed = subprocess.run(command, cwd=REPO_ROOT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed with exit code {completed.returncode}.")


def git_output(*args: str) -> str:
    """!
    @brief Return trimmed output from one Git command.
    @param[in] args Arguments passed after `git`.
    @return Decoded, trimmed Git stdout.
    """

    return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).strip()


def require_clean_revision() -> str:
    """!
    @brief Require that HEAD names a clean, commit-scoped source revision.
    @return Full SHA of the clean HEAD revision.
    """

    if git_output("status", "--porcelain"):
        raise RuntimeError(
            "Documentation certification requires a clean working tree so its result can be "
            "claimed for one commit. Commit or stash changes, then rerun."
        )
    return git_output("rev-parse", "HEAD")


def ensure_empty_doxygen_warning_log() -> None:
    """!
    @brief Fail certification when Doxygen reports a warning.
    @return None.
    """

    warning_log = REPO_ROOT / "logs" / "doxygen.warnings"
    warnings = warning_log.read_text(encoding="utf-8").strip() if warning_log.exists() else ""
    if warnings:
        raise RuntimeError(
            "Doxygen emitted warnings; documentation cannot be certified.\n" + warnings
        )


def invariant_contract_scope() -> str:
    """!
    @brief Describe how much of the invariant surface has an automated checker.

    Most invariant contracts are tracked rather than enforced, so the certificate must
    not imply that every guarantee is verified.
    @return Human-readable contract scope statement.
    """
    path = REPO_ROOT / "tests" / "tooling" / "contract_registry.json"
    try:
        contracts = json.loads(path.read_text(encoding="utf-8"))["contracts"]
    except (OSError, ValueError, KeyError):
        return "scope unknown (contract registry unreadable)"
    enforced = [c["id"] for c in contracts if c["status"] == "enforced"]
    other = [c for c in contracts if c["status"] != "enforced"]
    return (
        f"{len(enforced)} enforced and verified; "
        f"{len(other)} TRACKED OR PLANNED without an automated checker"
    )


def capability_coverage_scope() -> str:
    """!
    @brief Describe how much of the capability surface has enforced documentation coverage.

    The certificate must not imply that every capability is documented while families
    are still being backfilled, so it names what is enforced and what is advisory.
    @return Human-readable coverage scope statement.
    """
    registry_path = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
    try:
        families = json.loads(registry_path.read_text(encoding="utf-8"))["families"]
    except (OSError, ValueError, KeyError):
        return "scope unknown (capability registry unreadable)"
    enforced = sorted(f["id"] for f in families if f.get("coverage_enforced"))
    advisory = sorted(f["id"] for f in families if not f.get("coverage_enforced"))
    if not advisory:
        return f"enforced for all {len(enforced)} registered families"
    return (
        f"enforced for {', '.join(enforced) or 'no families'}; "
        f"ADVISORY ONLY for {', '.join(advisory)} (backfill pending)"
    )


def freshness_scope() -> str:
    """!
    @brief Describe the four freshness states separately.

    @details Mechanical verification, a current fingerprint, a soft suspicion, and a
             human review are four different claims. The certificate must not blur
             them: a matching digest records that a review happened, never that the
             prose it covers is correct.
    @return Human-readable freshness statement.
    """
    path = REPO_ROOT / "tests" / "tooling" / "freshness_manifest.json"
    try:
        surfaces = json.loads(path.read_text(encoding="utf-8"))["surfaces"]
    except (OSError, ValueError, KeyError):
        return "scope unknown (freshness manifest unreadable)"
    import hashlib

    def digest(paths: list) -> str:
        """!
        @brief Digest matching audit_freshness.digest_of.
        @param[in] paths Repository-relative paths.
        @return Hex digest with algorithm prefix.
        """
        accumulator = hashlib.sha256()
        for relative in sorted(paths):
            accumulator.update(relative.encode("utf-8"))
            accumulator.update(b"\0")
            accumulator.update((REPO_ROOT / relative).read_bytes())
            accumulator.update(b"\0")
        return f"sha256:{accumulator.hexdigest()}"

    hard_current = soft_current = suspect = never = 0
    for surface in surfaces:
        paths = ([surface["artifact"]] if surface["tier"] == "hard"
                 else list(surface.get("watched_paths", [])))
        attested = surface.get("attested_digest")
        if attested in (None, "unattested"):
            never += 1
        elif attested == digest(paths):
            if surface["tier"] == "hard":
                hard_current += 1
            else:
                soft_current += 1
        else:
            suspect += 1
    return (
        f"{hard_current} hard-current (normalized artifact unchanged since review); "
        f"{soft_current} soft-current (watched sources unchanged since review); "
        f"{suspect} SUSPECT; {never} NEVER ATTESTED against their sources"
    )


def write_certificate(sha: str, runtime_checked: bool) -> Path:
    """!
    @brief Write a reproducible certification record for one validated revision.
    @param[in] sha Full Git revision validated by the certification gates.
    @param[in] runtime_checked Whether the PETSc/MPI runtime suite was included.
    @return Path of the generated certificate.
    """

    certificate = REPO_ROOT / "logs" / f"documentation-certificate-{sha}.md"
    coverage_scope = capability_coverage_scope()
    contract_scope = invariant_contract_scope()
    freshness = freshness_scope()
    timestamp = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    runtime_gate = (
        "`make check-full` (including `unit-io` startup-banner and `unit-logging` contracts)"
        if runtime_checked
        else "not run"
    )
    scope_statement = (
        "This full certification is authoritative only for the commit named above. "
        "Any source, configuration, example, or documentation change requires a new certification.\n"
        if runtime_checked
        else "This structural certification is authoritative only for the commit named above and "
        "does not certify PETSc/MPI runtime behavior. A pre-push policy may reuse a recent full "
        "certificate only when it establishes that runtime-relevant files are unchanged.\n"
    )
    certificate.write_text(
        "# PICurv Documentation Certification\n\n"
        f"- Certified commit: `{sha}`\n"
        f"- Certified at (UTC): `{timestamp}`\n"
        "- Markdown links (all tracked and non-ignored Markdown): passed\n"
        "- Public-header and implementation comment audit: passed\n"
        "- User-facing C/Python reporting audit: passed\n"
        "- Starter template, example, and configuration audit: passed\n"
        "- Generic documentation-expansion debris audit: passed\n"
        f"- Capability parity (full source chain): passed\n"
        f"- Invariant contract registry: {contract_scope}\n"
        f"- Tier-2 documentation coverage: {coverage_scope}\n"
        "- Page-type coverage (every published page typed): passed\n"
        "- Subsystem lifecycle obligations: passed\n"
        f"- Documentation freshness: {freshness}\n"
        "- PETSc option-ingress manifest audit: passed\n"
        "- Example/template and configuration regression suite: passed\n"
        "- Doxygen HTML build: passed with zero warnings\n"
        "- Published-site URLs and navigation (validated against generated HTML): passed\n"
        f"- PETSc/MPI runtime validation: {runtime_gate}\n\n"
        + scope_statement
        + "\nThese lines record four different kinds of claim and must not be conflated. A "
          "gate that\npassed is **mechanically verified**. A hard-current surface means the "
          "normalized artifact a\npage describes has not changed since someone last compared "
          "them - not that the page is\ncorrect. A soft suspicion is advisory. A never-attested "
          "surface means no such comparison has\nhappened at all. **Scientific and visual review "
          "by the repository owner is a separate,\nhuman act that no gate here performs.**\n",
        encoding="utf-8",
    )
    return certificate


def parse_args(argv: list[str]) -> argparse.Namespace:
    """!
    @brief Parse documentation-certification command-line arguments.
    @param[in] argv Arguments excluding the executable name.
    @return Parsed command-line namespace.
    """

    parser = argparse.ArgumentParser(
        description="Validate and certify PICurv documentation for the current clean Git commit."
    )
    parser.add_argument(
        "--runtime",
        action="store_true",
        help=(
            "also run the PETSc/MPI `make check-full` runtime validation gate "
            "(including startup-banner and logging contracts)"
        ),
    )
    parser.add_argument(
        "--no-certificate",
        action="store_true",
        help="run the gates without writing logs/documentation-certificate-<sha>.md",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """!
    @brief Run the commit-scoped documentation certification workflow.
    @param[in] argv Optional command-line arguments excluding the executable name.
    @return Process status code.
    """

    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        sha = require_clean_revision()
        run([sys.executable, "tests/tooling/check_markdown_links.py"], "Markdown links")
        run([sys.executable, "tests/tooling/audit_function_docs.py"], "function documentation")
        run([sys.executable, "tests/tooling/audit_user_facing_reporting.py"], "user-facing reporting")
        run([sys.executable, "tests/tooling/audit_starter_content.py"], "starter templates and configuration")
        run([sys.executable, "tests/tooling/audit_generic_expansion.py"], "generic documentation-expansion debris")
        run([sys.executable, "tests/tooling/audit_path_literals.py"], "unmanaged run-path literals")
        run([sys.executable, "tests/tooling/audit_capability_coverage.py"], "capability parity and documentation coverage")
        run([sys.executable, "tests/tooling/audit_family_census.py"], "capability family census")
        run([sys.executable, "tests/tooling/generate_cli_reference.py", "--check"], "generated CLI reference")
        run([sys.executable, "tests/tooling/audit_ingress.py"], "option-ingress manifest")
        run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                "tests/test_repo_consistency.py",
                "tests/test_config_regressions.py",
                "tests/test_capability_tooling.py",
            ],
            "example/template, configuration, and capability-tooling validation",
        )
        run(["make", "--no-print-directory", "build-docs"], "Doxygen HTML build")
        ensure_empty_doxygen_warning_log()
        # Runs after the build: it validates canonical URLs and navigation tabs against the
        # generated HTML, so an excluded or renamed page cannot pass on a source declaration.
        run([sys.executable, "tests/tooling/audit_docs_site.py"], "published-site URLs and navigation")
        run([sys.executable, "tests/tooling/audit_page_types.py"], "page-type coverage")
        run([sys.executable, "tests/tooling/audit_inline_choices.py"],
            "named public choice sets")
        run([sys.executable, "tests/tooling/audit_subsystem_lifecycle.py"],
            "subsystem lifecycle obligations")
        run([sys.executable, "tests/tooling/audit_freshness.py"],
            "documentation freshness (hard blocking, soft advisory)")
        run([sys.executable, "tests/tooling/audit_contracts.py"], "invariant contract registry")
        if args.runtime:
            run(
                ["make", "--no-print-directory", "check-full"],
                "PETSc/MPI runtime validation (including startup-banner and logging contracts)",
            )
        if not args.no_certificate:
            certificate = write_certificate(sha, args.runtime)
            certification_kind = "full" if args.runtime else "structural"
            print(f"{certification_kind.capitalize()} documentation certification through commit {sha}: {certificate}")
        return 0
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Documentation certification failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
