#!/usr/bin/env python3
"""Distill Doxygen XML references into an optional, bounded review index."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


SCHEMA_VERSION = 1
SOURCE_ROOTS = ("include", "src", "picurv_cli", "generators")
SOURCE_SUFFIXES = {".c", ".cc", ".cpp", ".h", ".hpp", ".py", ".sh", ".flow"}


def source_files(repo_root: Path) -> list[Path]:
    """!
    @brief Discover source inputs whose dirty bytes determine cross-reference freshness.
    @param[in] repo_root Repository root.
    @return Sorted repository-relative paths.
    """

    found: list[Path] = []
    for root_name in SOURCE_ROOTS:
        root = repo_root / root_name
        if not root.is_dir():
            continue
        found.extend(
            path.relative_to(repo_root)
            for path in root.rglob("*")
            if path.is_file() and path.suffix in SOURCE_SUFFIXES and "__pycache__" not in path.parts
        )
    return sorted(set(found), key=lambda path: path.as_posix())


def digest_files(repo_root: Path, paths: list[Path]) -> str:
    """!
    @brief Hash path names and current file bytes deterministically.
    @param[in] repo_root Repository root.
    @param[in] paths Repository-relative paths to hash.
    @return SHA-256 digest with algorithm prefix.
    """

    accumulator = hashlib.sha256()
    for relative in sorted(paths, key=lambda path: path.as_posix()):
        accumulator.update(relative.as_posix().encode("utf-8"))
        accumulator.update(b"\0")
        accumulator.update((repo_root / relative).read_bytes())
        accumulator.update(b"\0")
    return f"sha256:{accumulator.hexdigest()}"


def file_digest(path: Path) -> str:
    """!
    @brief Hash one configuration file.
    @param[in] path File to hash.
    @return SHA-256 digest with algorithm prefix.
    """

    return f"sha256:{hashlib.sha256(path.read_bytes()).hexdigest()}"


def repository_path(repo_root: Path, raw: str) -> str:
    """!
    @brief Normalize a Doxygen location to a repository-relative POSIX path when possible.
    @param[in] repo_root Repository root.
    @param[in] raw Location value emitted by Doxygen.
    @return Normalized path string.
    """

    path = Path(raw)
    if path.is_absolute():
        try:
            path = path.resolve().relative_to(repo_root.resolve())
        except ValueError:
            return path.as_posix()
    return path.as_posix().lstrip("./")


def parse_xml(repo_root: Path, xml_dir: Path) -> dict[str, dict]:
    """!
    @brief Extract definition locations and Doxygen reference edges from XML output.
    @param[in] repo_root Repository root used to normalize locations.
    @param[in] xml_dir Doxygen XML directory.
    @return Mapping from Doxygen member id to distilled symbol record.
    """

    symbols: dict[str, dict] = {}
    for xml_path in sorted(xml_dir.glob("*.xml")):
        try:
            document = ET.parse(xml_path)
        except ET.ParseError as error:
            raise RuntimeError(f"cannot parse {xml_path}: {error}") from error
        for member in document.findall(".//memberdef"):
            refid = member.get("id")
            location = member.find("location")
            name = member.findtext("name")
            if not refid or location is None or not name or not location.get("file"):
                continue
            outgoing = sorted(
                {ref.get("refid") for ref in member.findall("references") if ref.get("refid")}
            )
            incoming = sorted(
                {ref.get("refid") for ref in member.findall("referencedby") if ref.get("refid")}
            )
            symbols[refid] = {
                "name": name,
                "qualified_name": member.findtext("qualifiedname") or name,
                "kind": member.get("kind", "unknown"),
                "definition": {
                    "path": repository_path(repo_root, location.get("file", "")),
                    "line": int(location.get("line", "0") or 0),
                },
                "outgoing": outgoing,
                "incoming": incoming,
            }
        document.getroot().clear()
    return distill_symbols(symbols)


def distill_symbols(symbols: dict[str, dict]) -> dict[str, dict]:
    """!
    @brief Keep production definitions and remap header references to source definitions.
    @param[in] symbols Raw Doxygen member records.
    @return Compact production-oriented symbol graph.
    """

    source_functions = {
        record["name"]: refid
        for refid, record in symbols.items()
        if record["kind"] == "function" and record["definition"]["path"].startswith("src/")
    }
    retained: set[str] = set()
    for refid, record in symbols.items():
        path = record["definition"]["path"]
        kind = record["kind"]
        if path.startswith(("src/", "picurv_cli/", "generators/")) and kind in {
            "function",
            "variable",
            "enum",
            "typedef",
        }:
            retained.add(refid)
        elif path.startswith("include/") and (
            kind in {"enum", "typedef"}
            or (kind == "function" and record["name"] not in source_functions)
        ):
            retained.add(refid)

    def canonical(refid: str) -> str | None:
        """!
        @brief Resolve a raw edge to a retained implementation-oriented member.
        @param[in] refid Raw Doxygen member id.
        @return Retained member id, or None when the edge is outside the compact index.
        """

        if refid in retained:
            return refid
        record = symbols.get(refid)
        if record and record["kind"] == "function":
            return source_functions.get(record["name"])
        return None

    distilled: dict[str, dict] = {}
    for refid in sorted(retained):
        record = dict(symbols[refid])
        record["incoming"] = sorted(
            {mapped for edge in record["incoming"] if (mapped := canonical(edge)) is not None}
        )
        record["outgoing"] = sorted(
            {mapped for edge in record["outgoing"] if (mapped := canonical(edge)) is not None}
        )
        distilled[refid] = record
    return distilled


def command_output(command: list[str], cwd: Path) -> str:
    """!
    @brief Return the first line from an orientation-only metadata command.
    @param[in] command Executable and arguments.
    @param[in] cwd Working directory.
    @return First stdout line, or `unavailable` when the command fails.
    """

    completed = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    if completed.returncode != 0 or not completed.stdout.strip():
        return "unavailable"
    return completed.stdout.strip().splitlines()[0]


def generate(repo_root: Path, xml_dir: Path, output: Path) -> None:
    """!
    @brief Generate the stamped JSON cross-reference index.
    @param[in] repo_root Repository root.
    @param[in] xml_dir Doxygen XML directory.
    @param[in] output Destination JSON path.
    @return None.
    @throws RuntimeError when XML or required configuration is unavailable.
    """

    if not xml_dir.is_dir() or not any(xml_dir.glob("*.xml")):
        raise RuntimeError(
            f"Doxygen XML is unavailable at {xml_dir}; run doxygen docs/Doxyfile first"
        )
    doxyfile = repo_root / "docs" / "Doxyfile"
    if not doxyfile.is_file():
        raise RuntimeError(f"missing Doxygen configuration: {doxyfile}")
    inputs = source_files(repo_root)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "source_digest": digest_files(repo_root, inputs),
        "source_files": [path.as_posix() for path in inputs],
        "doxyfile_digest": file_digest(doxyfile),
        "doxygen_version": command_output(["doxygen", "--version"], repo_root),
        "git_commit": command_output(["git", "rev-parse", "HEAD"], repo_root),
        "symbols": parse_xml(repo_root, xml_dir),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    """!
    @brief Parse generator command-line arguments.
    @param[in] argv Arguments excluding the executable name.
    @return Parsed command-line namespace.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--xml-dir", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """!
    @brief Generate the optional index and report actionable failures.
    @param[in] argv Optional arguments excluding the executable name.
    @return Zero on success and one on failure.
    """

    args = parse_args(sys.argv[1:] if argv is None else argv)
    repo_root = args.repo_root.resolve()
    xml_dir = (args.xml_dir or repo_root / "docs_build" / "xml").resolve()
    output = (args.output or repo_root / "docs_build" / "xref.json").resolve()
    try:
        generate(repo_root, xml_dir, output)
    except (OSError, RuntimeError) as error:
        print(f"Cross-reference generation failed: {error}", file=sys.stderr)
        return 1
    try:
        relative = output.relative_to(repo_root)
    except ValueError:
        relative = output
    print(f"Doxygen cross-reference index generated: {relative}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
