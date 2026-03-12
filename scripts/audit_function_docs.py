#!/usr/bin/env python3
"""!
@file audit_function_docs.py
@brief Audits C and Python function documentation coverage across the repository.

This script enforces the repository's function-level documentation contract for:

- public C declarations in `include/`,
- C definitions in `src/` and `tests/c/`,
- Python functions in `scripts/`, `tests/`, and Python-backed executable scripts.

It is intentionally lightweight. The C side uses signature scanning instead of a
full parser, while the Python side uses `ast`.
"""

from __future__ import annotations

import ast
import re
import sys
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

C_HEADER_DIRS = (REPO_ROOT / "include",)
C_SOURCE_DIRS = (REPO_ROOT / "src", REPO_ROOT / "tests" / "c")
PYTHON_DIRS = (REPO_ROOT / "scripts", REPO_ROOT / "tests")
PYTHON_EXTRA_FILES = (
    REPO_ROOT / "scripts" / "picurv",
    REPO_ROOT / "scripts" / "grid.gen",
)

C_DECL_START_RE = re.compile(
    r"^\s*(?!typedef\b)(?!if\b)(?!for\b)(?!while\b)(?!switch\b)(?!return\b)(?!else\b)"
    r"(?:extern\s+)?(?:static\s+)?(?:inline\s+)?(?:const\s+)?(?:unsigned\s+|signed\s+)?"
    r"(?:[A-Za-z_][A-Za-z0-9_]*\s+)+(?:\*\s*)*"
    r"([A-Za-z_][A-Za-z0-9_]*)\s*\("
)
C_PARAM_RE = re.compile(r"@param(?:\[[^\]]+\])?\s+([A-Za-z_][A-Za-z0-9_]*)")


@dataclass(frozen=True)
class AuditFinding:
    """!
    @brief Represents one audit failure.
    @param[in] path Repository-relative path containing the failure.
    @param[in] line 1-based source line associated with the failure.
    @param[in] symbol Function symbol being audited.
    @param[in] message Human-readable failure description.
    """

    path: str
    line: int
    symbol: str
    message: str


def _iter_c_files(directories: tuple[Path, ...]) -> list[Path]:
    """!
    @brief Returns all C or header files below the configured directories.
    @param[in] directories Root directories to scan.
    @return Sorted list of matching file paths.
    """

    files: list[Path] = []
    for directory in directories:
        if not directory.exists():
            continue
        files.extend(sorted(path for path in directory.rglob("*") if path.suffix in {".c", ".h"}))
    return sorted(files)


def _iter_python_files() -> list[Path]:
    """!
    @brief Returns all Python source files covered by the audit.
    @return Sorted list of Python-backed source files.
    """

    files: set[Path] = set()
    for directory in PYTHON_DIRS:
        if not directory.exists():
            continue
        files.update(path for path in directory.rglob("*.py"))

    for path in PYTHON_EXTRA_FILES:
        if path.exists():
            files.add(path)

    return sorted(files)


def _read_lines(path: Path) -> list[str]:
    """!
    @brief Reads a text file into a list of lines.
    @param[in] path Path to read.
    @return File contents split into lines without trailing newline markers.
    """

    return path.read_text(encoding="utf-8", errors="ignore").splitlines()


def _relative_path(path: Path) -> str:
    """!
    @brief Returns a repository-relative path string.
    @param[in] path Absolute or repository-local path.
    @return POSIX-style repository-relative path.
    """

    return path.relative_to(REPO_ROOT).as_posix()


def _find_attached_doxygen_block(lines: list[str], start_line: int) -> tuple[int, int] | None:
    """!
    @brief Finds the Doxygen block immediately attached to a declaration or definition.
    @param[in] lines File content lines.
    @param[in] start_line 0-based line index where the symbol begins.
    @return `(start, end)` line indices for the attached block, or `None`.
    """

    probe = start_line - 1
    while probe >= 0 and lines[probe].strip() == "":
        probe -= 1

    if probe < 0 or "*/" not in lines[probe]:
        return None

    end = probe
    while probe >= 0 and "/**" not in lines[probe]:
        probe -= 1

    if probe < 0:
        return None

    return probe, end


def _split_c_parameters(signature: str) -> list[str]:
    """!
    @brief Splits a C signature parameter list into parameter names.
    @param[in] signature Full function signature text.
    @return Ordered list of parameter names excluding `void` and variadics.
    """

    start = signature.find("(")
    end = signature.rfind(")")
    if start < 0 or end < 0 or end <= start:
        return []

    raw = signature[start + 1:end]
    params: list[str] = []
    depth = 0
    current: list[str] = []
    for char in raw:
        if char == "," and depth == 0:
            params.append("".join(current).strip())
            current = []
            continue
        current.append(char)
        if char in "([{":
            depth += 1
        elif char in ")]}":
            depth -= 1
    if current:
        params.append("".join(current).strip())

    names: list[str] = []
    for param in params:
        if not param or param == "void" or param == "...":
            continue

        clean = re.sub(r"\b(const|volatile|restrict|extern|static|register|inline)\b", "", param)
        clean = clean.strip()
        match = re.search(r"([A-Za-z_][A-Za-z0-9_]*)\s*(?:\[[^\]]*\]\s*)*$", clean)
        if match:
            names.append(match.group(1))

    return names


def _c_return_type(signature: str, symbol: str) -> str:
    """!
    @brief Extracts the declared C return type prefix for one signature.
    @param[in] signature Full function signature text.
    @param[in] symbol Function name contained in the signature.
    @return Normalized return-type prefix.
    """

    prefix = signature.split(symbol, 1)[0]
    return " ".join(prefix.split())


def _return_tag_required(return_type: str) -> bool:
    """!
    @brief Reports whether a Doxygen `@return` tag is required for a C symbol.
    @param[in] return_type Normalized return-type prefix.
    @return `True` when the symbol does not return `void`.
    """

    stripped = return_type.replace("extern ", "").replace("static ", "").replace("inline ", "").strip()
    return not stripped.startswith("void")


def _collect_c_signatures(path: Path, require_terminator: str) -> list[tuple[int, str, str]]:
    """!
    @brief Collects C signatures from a header or source file.
    @param[in] path File to scan.
    @param[in] require_terminator Expected signature terminator, either `;` or `{`.
    @return List of `(start_line, symbol, signature_text)` tuples.
    """

    lines = _read_lines(path)
    signatures: list[tuple[int, str, str]] = []
    line_index = 0
    in_block_comment = False

    while line_index < len(lines):
        stripped = lines[line_index].lstrip()
        if in_block_comment:
            if "*/" in lines[line_index]:
                in_block_comment = False
            line_index += 1
            continue

        if "/*" in lines[line_index]:
            if "*/" not in lines[line_index]:
                in_block_comment = True
            line_index += 1
            continue

        if stripped.startswith(("#", "/*", "*", "//")) or "(" not in lines[line_index]:
            line_index += 1
            continue

        match = C_DECL_START_RE.match(lines[line_index])
        if not match:
            line_index += 1
            continue

        symbol = match.group(1)
        start_line = line_index
        signature = lines[line_index].rstrip()
        while line_index + 1 < len(lines) and require_terminator not in signature and ";" not in signature:
            line_index += 1
            signature += " " + lines[line_index].strip()

        if require_terminator == ";" and ";" in signature:
            signatures.append((start_line, symbol, signature))
        elif require_terminator == "{" and "{" in signature and ";" not in signature.split("{", 1)[0]:
            signatures.append((start_line, symbol, signature))

        line_index += 1

    return signatures


def _audit_c_header(path: Path) -> list[AuditFinding]:
    """!
    @brief Audits public C declarations in one header file.
    @param[in] path Header file to scan.
    @return Findings emitted for the header.
    """

    findings: list[AuditFinding] = []
    lines = _read_lines(path)
    for start_line, symbol, signature in _collect_c_signatures(path, ";"):
        block_range = _find_attached_doxygen_block(lines, start_line)
        if block_range is None:
            findings.append(AuditFinding(_relative_path(path), start_line + 1, symbol, "missing attached Doxygen block"))
            continue

        block = "\n".join(lines[block_range[0]:block_range[1] + 1])
        if "@brief" not in block:
            findings.append(AuditFinding(_relative_path(path), start_line + 1, symbol, "missing @brief tag"))

        declared_params = _split_c_parameters(signature)
        documented_params = set(C_PARAM_RE.findall(block))
        if set(declared_params) != documented_params:
            findings.append(
                AuditFinding(
                    _relative_path(path),
                    start_line + 1,
                    symbol,
                    f"documented @param names {sorted(documented_params)} do not match declaration {declared_params}",
                )
            )

        if _return_tag_required(_c_return_type(signature, symbol)) and "@return" not in block:
            findings.append(AuditFinding(_relative_path(path), start_line + 1, symbol, "missing @return tag"))

    return findings


def _audit_c_source(path: Path) -> list[AuditFinding]:
    """!
    @brief Audits function definitions in one C source file.
    @param[in] path Source file to scan.
    @return Findings emitted for the source file.
    """

    findings: list[AuditFinding] = []
    lines = _read_lines(path)
    for start_line, symbol, _signature in _collect_c_signatures(path, "{"):
        block_range = _find_attached_doxygen_block(lines, start_line)
        if block_range is None:
            findings.append(AuditFinding(_relative_path(path), start_line + 1, symbol, "missing attached Doxygen block"))
            continue

        block = "\n".join(lines[block_range[0]:block_range[1] + 1])
        if "@brief" not in block:
            findings.append(AuditFinding(_relative_path(path), start_line + 1, symbol, "missing @brief tag"))

    return findings


def _python_parameter_names(node: ast.FunctionDef | ast.AsyncFunctionDef) -> list[str]:
    """!
    @brief Returns the meaningful Python parameter names for one function node.
    @param[in] node Function AST node.
    @return Ordered list of parameters expected in `@param` tags.
    """

    names = [arg.arg for arg in node.args.posonlyargs + node.args.args + node.args.kwonlyargs]
    names = [name for name in names if name not in {"self", "cls"}]
    if node.args.vararg is not None:
        names.append(node.args.vararg.arg)
    if node.args.kwarg is not None:
        names.append(node.args.kwarg.arg)
    return names


def _python_requires_return(node: ast.FunctionDef | ast.AsyncFunctionDef) -> bool:
    """!
    @brief Reports whether one Python function should document a return value.
    @param[in] node Function AST node.
    @return `True` when the function returns a non-`None` value.
    """

    for child in ast.walk(node):
        if isinstance(child, ast.Return) and child.value is not None:
            if isinstance(child.value, ast.Constant) and child.value.value is None:
                continue
            return True
    return False


def _audit_python_file(path: Path) -> list[AuditFinding]:
    """!
    @brief Audits Python function docstrings in one file.
    @param[in] path Python source file to scan.
    @return Findings emitted for the Python file.
    """

    findings: list[AuditFinding] = []
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))

    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue

        docstring = ast.get_docstring(node)
        if docstring is None:
            findings.append(AuditFinding(_relative_path(path), node.lineno, node.name, "missing Python docstring"))
            continue

        if "@brief" not in docstring:
            findings.append(AuditFinding(_relative_path(path), node.lineno, node.name, "missing @brief tag"))

        declared_params = _python_parameter_names(node)
        documented_params = set(re.findall(r"@param(?:\[[^\]]+\])?\s+([A-Za-z_][A-Za-z0-9_]*)", docstring))
        if set(declared_params) != documented_params:
            findings.append(
                AuditFinding(
                    _relative_path(path),
                    node.lineno,
                    node.name,
                    f"documented @param names {sorted(documented_params)} do not match declaration {declared_params}",
                )
            )

        if _python_requires_return(node) and "@return" not in docstring:
            findings.append(AuditFinding(_relative_path(path), node.lineno, node.name, "missing @return tag"))

    return findings


def _collect_findings() -> list[AuditFinding]:
    """!
    @brief Runs the full repository documentation audit.
    @return Sorted list of all findings emitted by the audit.
    """

    findings: list[AuditFinding] = []
    for path in _iter_c_files(C_HEADER_DIRS):
        findings.extend(_audit_c_header(path))
    for path in _iter_c_files(C_SOURCE_DIRS):
        if path.suffix == ".c":
            findings.extend(_audit_c_source(path))
    for path in _iter_python_files():
        findings.extend(_audit_python_file(path))

    return sorted(findings, key=lambda item: (item.path, item.line, item.symbol, item.message))


def _print_findings(findings: list[AuditFinding]) -> None:
    """!
    @brief Prints findings in a grep-friendly format.
    @param[in] findings Findings to render.
    """

    for finding in findings:
        print(f"{finding.path}:{finding.line}: {finding.symbol}: {finding.message}")


def main() -> int:
    """!
    @brief Runs the repository function documentation audit from the command line.
    @return Process exit status.
    """

    findings = _collect_findings()
    if findings:
        _print_findings(findings)
        print(f"\nFound {len(findings)} documentation issue(s).", file=sys.stderr)
        return 1

    print("Function documentation audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
