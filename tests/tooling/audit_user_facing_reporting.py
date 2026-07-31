#!/usr/bin/env python3
"""Validate the repository's option-aware user-facing reporting contract."""

from __future__ import annotations

import re
import sys
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
ALLOWED_DIRECT_C_OUTPUT = {
    "io.c",             # Rank-zero effective-configuration banner.
    "logging.c",        # Logging and diagnostic rendering sinks.
    "setup.c",          # Bootstrap warnings before logging is configured.
    "runloop.c",        # Progress-bar line control.
    "postprocessor.c",  # Progress-bar line control.
    "momentumsolvers.c",  # Unconditional convergence warnings.
}
BANNER_REQUIRED = (
    "Momentum Equation Solver",
    "Initial Pseudo-CFL (Courant)",
    "Newton-Krylov PETSc Controls",
    "Pseudo-Time Controller      : NOT APPLICABLE",
    "Console Log Level",
    "Profiling Timestep Output",
    "Runtime Memory Log",
    "Solution Convergence Mode",
)
CONTRACT_PATH = REPO_ROOT / "tests" / "tooling" / "user_facing_reporting_contract.json"
CLI_REQUIRED = ("Momentum Solver: Dual Time Picard Jameson RK", "Momentum Solver: Newton Krylov",
                "Momentum Solver: Explicit RK (no pseudo-time controller)",
                "Dual-Time Pseudo-Time Controls", "Newton--Krylov Controls")


def strip_c_comments(text: str) -> str:
    """!
    @brief Replace C comments while preserving source line positions.
    @param[in] text C source text to normalize.
    @return Source text with comments replaced by whitespace/newlines.
    """

    text = re.sub(r"/\*.*?\*/", lambda match: "\n" * match.group(0).count("\n"), text, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", text)


def audit_direct_c_output() -> list[str]:
    """!
    @brief Return raw C console emissions outside their approved rendering modules.
    @return Human-readable findings for every disallowed emission.
    """

    findings: list[str] = []
    pattern = re.compile(r"\b(?:PetscPrintf|PetscSynchronizedPrintf|printf)\s*\(")
    for path in sorted(SRC_DIR.glob("*.c")):
        if path.name in ALLOWED_DIRECT_C_OUTPUT:
            continue
        stripped = strip_c_comments(path.read_text(encoding="utf-8", errors="ignore"))
        for match in pattern.finditer(stripped):
            line = stripped.count("\n", 0, match.start()) + 1
            findings.append(
                f"{path.relative_to(REPO_ROOT)}:{line}: raw console output must use the logging policy "
                "or be moved to an approved reporting surface."
            )
    return findings


def require_fragments(path: Path, fragments: tuple[str, ...], surface: str) -> list[str]:
    """!
    @brief Return missing reporting-contract fragments for one user-facing surface.
    @param[in] path Source file containing the surface.
    @param[in] fragments Required literal contract fragments.
    @param[in] surface Human-readable surface name for failures.
    @return Findings for missing fragments.
    """

    text = path.read_text(encoding="utf-8")
    return [f"{surface}: missing required option-aware reporting fragment {fragment!r}."
            for fragment in fragments if fragment not in text]


def audit_cli_contract(contract: dict) -> list[str]:
    """!
    @brief Verify that every declared CLI command has a parser, dispatch handler, and contextual message.
    @param[in] contract Parsed user-facing reporting contract mapping.
    @return Findings for missing CLI command-reporting contract elements.
    """

    cli_text = (REPO_ROOT / "picurv_cli" / "cli.py").read_text(encoding="utf-8")
    core_text = (REPO_ROOT / "picurv_cli" / "core.py").read_text(encoding="utf-8")
    findings: list[str] = []
    for command, entry in contract["cli_commands"].items():
        for field, text, surface in (
            ("parser", cli_text, "CLI parser"),
            ("handler", core_text, "CLI handler"),
            ("context", cli_text + core_text, "CLI context"),
        ):
            if entry[field] not in text:
                findings.append(f"{surface} for '{command}': missing contract fragment {entry[field]!r}.")
    return findings


def main() -> int:
    """!
    @brief Audit C banner, CLI summaries, and direct-output routing.
    @return Process status code.
    """

    findings = audit_direct_c_output()
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    findings.extend(require_fragments(
        REPO_ROOT / contract["banner"]["source"], tuple(contract["banner"]["required_fragments"]), "C startup banner"
    ))
    findings.extend(require_fragments(REPO_ROOT / "picurv_cli" / "core.py", CLI_REQUIRED, "Python CLI/summarize"))
    findings.extend(audit_cli_contract(contract))
    if findings:
        print("\n".join(findings), file=sys.stderr)
        return 1
    print("User-facing reporting audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
