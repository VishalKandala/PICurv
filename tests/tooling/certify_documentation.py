#!/usr/bin/env python3
"""Create a commit-scoped documentation certification after strict validation."""

from __future__ import annotations

import argparse
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


def write_certificate(sha: str, runtime_checked: bool) -> Path:
    """!
    @brief Write a reproducible certification record for one validated revision.
    @param[in] sha Full Git revision validated by the certification gates.
    @param[in] runtime_checked Whether the PETSc/MPI runtime suite was included.
    @return Path of the generated certificate.
    """

    certificate = REPO_ROOT / "logs" / f"documentation-certificate-{sha}.md"
    timestamp = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    runtime_gate = "`make check-full`" if runtime_checked else "not run"
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
        "- Markdown links: passed\n"
        "- Public-header and implementation comment audit: passed\n"
        "- PETSc option-ingress manifest audit: passed\n"
        "- Example/template and configuration regression suite: passed\n"
        "- Doxygen HTML build: passed with zero warnings\n"
        f"- PETSc/MPI runtime validation: {runtime_gate}\n\n"
        + scope_statement,
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
        help="also run the PETSc/MPI `make check-full` runtime validation gate",
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
        run([sys.executable, "tests/tooling/audit_ingress.py"], "option-ingress manifest")
        run(
            [sys.executable, "-m", "pytest", "-q", "tests/test_repo_consistency.py", "tests/test_config_regressions.py"],
            "example/template and configuration validation",
        )
        run(["make", "--no-print-directory", "build-docs"], "Doxygen HTML build")
        ensure_empty_doxygen_warning_log()
        if args.runtime:
            run(["make", "--no-print-directory", "check-full"], "PETSc/MPI runtime validation")
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
