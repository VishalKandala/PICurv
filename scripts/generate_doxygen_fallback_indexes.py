#!/usr/bin/env python3
"""Generate robust Doxygen index pages and structured reference views."""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

HEADER_SUFFIXES = {".h", ".hpp"}
SOURCE_SUFFIXES = {".c", ".cc", ".cpp"}
SCRIPT_SUFFIXES = {".py", ".sh", ".flow"}
REPO_BLOB_URL = "https://github.com/VishalKandala/PICurv/blob/main/"
IGNORED_STRUCT_NAMES = {"Name"}

NAMED_STRUCT_RE = re.compile(r"\bstruct\s+([A-Za-z_]\w*)\s*\{")
TYPEDEF_START_RE = re.compile(r"^\s*typedef\s+struct(?:\s+([A-Za-z_]\w*))?")
TYPEDEF_END_RE = re.compile(r"^\s*}\s*([A-Za-z_]\w*)\s*;")


def doxygen_file_page(name: str) -> str:
    """Perform doxygen file page."""
    return name.replace("_", "__").replace(".", "_8") + ".html"


def doxygen_file_page_with_path(rel_path: str) -> str:
    """Perform doxygen file page with path."""
    return rel_path.replace("_", "__").replace("/", "_2").replace(".", "_8") + ".html"


def needs_files_fallback(path: Path) -> bool:
    """Perform needs files fallback."""
    if not path.exists():
        return True
    text = path.read_text(encoding="utf-8", errors="ignore")
    return "Detailed file index was not generated in this build." in text


def needs_structs_fallback(path: Path) -> bool:
    """Perform needs structs fallback."""
    if not path.exists():
        return True
    text = path.read_text(encoding="utf-8", errors="ignore")
    return "Detailed structure index was not generated in this build." in text


def resolve_doxygen_file_href(html_dir: Path, rel_path: str) -> str:
    """Resolve doxygen file href."""
    name = Path(rel_path).name
    candidate = doxygen_file_page(name)
    if (html_dir / candidate).exists():
        return candidate
    candidate_path = doxygen_file_page_with_path(rel_path)
    if (html_dir / candidate_path).exists():
        return candidate_path
    return ""


def make_repo_href(rel_path: str) -> str:
    """Perform make repo href."""
    return REPO_BLOB_URL + rel_path


def collect_file_rows(repo_root: Path, html_dir: Path, base_dir: str, suffixes: set[str]) -> list[tuple[str, str, str]]:
    """Collect file rows."""
    rows: list[tuple[str, str, str]] = []
    root = repo_root / base_dir
    if not root.exists():
        return rows
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if suffixes and path.suffix.lower() not in suffixes:
            continue
        rel = path.relative_to(repo_root).as_posix()
        href = resolve_doxygen_file_href(html_dir, rel) or make_repo_href(rel)
        rows.append((path.name, rel, href))
    return rows


def collect_all_source_like_files(repo_root: Path, html_dir: Path) -> list[tuple[str, str, str]]:
    """Collect all source like files."""
    rows: list[tuple[str, str, str]] = []
    rows.extend(collect_file_rows(repo_root, html_dir, "include", HEADER_SUFFIXES))
    rows.extend(collect_file_rows(repo_root, html_dir, "src", SOURCE_SUFFIXES))
    rows.extend(collect_file_rows(repo_root, html_dir, "scripts", SCRIPT_SUFFIXES))
    return rows


def collect_struct_rows(repo_root: Path, html_dir: Path) -> list[tuple[str, str, str]]:
    """Collect struct rows."""
    struct_to_header: dict[str, str] = {}
    include_dir = repo_root / "include"
    if not include_dir.exists():
        return []

    for header in sorted(include_dir.rglob("*.h")):
        text = header.read_text(encoding="utf-8", errors="ignore")
        header_rel = header.relative_to(repo_root).as_posix()
        names = extract_struct_names(text)
        for name in names:
            struct_to_header.setdefault(name, header_rel)

    rows: list[tuple[str, str, str]] = []
    for name, header_rel in sorted(struct_to_header.items(), key=lambda item: item[0].lower()):
        if name in IGNORED_STRUCT_NAMES:
            continue
        struct_page = f"struct{name}.html"
        if (html_dir / struct_page).exists():
            href = struct_page
        else:
            href = resolve_doxygen_file_href(html_dir, header_rel) or make_repo_href(header_rel)
        rows.append((name, header_rel, href))
    return rows


def extract_struct_names(text: str) -> set[str]:
    """Extract struct names."""
    names: set[str] = set()

    for match in NAMED_STRUCT_RE.finditer(text):
        names.add(match.group(1))

    in_typedef_struct = False
    brace_depth = 0
    typedef_tag_name: str | None = None
    for line in text.splitlines():
        if not in_typedef_struct:
            start = TYPEDEF_START_RE.search(line)
            if not start:
                continue
            in_typedef_struct = True
            typedef_tag_name = start.group(1)
            if typedef_tag_name:
                names.add(typedef_tag_name)
            brace_depth = line.count("{") - line.count("}")
            if brace_depth <= 0:
                in_typedef_struct = False
                typedef_tag_name = None
            continue

        brace_depth += line.count("{") - line.count("}")
        end = TYPEDEF_END_RE.search(line)
        if end:
            names.add(end.group(1))
        if brace_depth <= 0:
            in_typedef_struct = False
            typedef_tag_name = None

    return names


def categorize_struct(name: str) -> str:
    """Categorize struct."""
    if name.startswith("BC") or "Boundary" in name or name == "FlowWave":
        return "Boundary Condition System"
    if name.startswith("IBM") or name in {"FSInfo", "SurfElmtInfo", "Cstart"}:
        return "Immersed Boundary and FSI"
    if name.startswith("Particle") or name in {"MigrationInfo"}:
        return "Particle Transport and Statistics"
    if name in {"SimCtx", "UserCtx", "UserMG", "MGCtx", "ScalingCtx", "DualMonitorCtx", "ProfiledFunction"}:
        return "Runtime Control and Solver Orchestration"
    if name.startswith("VTK") or name == "PostProcessParams":
        return "I/O and Postprocessing"
    if name in {"BoundingBox", "Cell", "Cmpnts", "Cmpnts2", "Cpt2D", "RankCellInfo", "RankNeighbors"}:
        return "Grid and Geometry"
    return "Generic Containers and Utilities"


def render_link(label: str, href: str) -> str:
    """Render link."""
    label_esc = html.escape(label)
    href_esc = html.escape(href)
    if href.startswith("http"):
        return f"<a class='el' href='{href_esc}' target='_blank' rel='noopener'>{label_esc}</a>"
    return f"<a class='el' href='{href_esc}'>{label_esc}</a>"


def render_rows(rows: list[tuple[str, str, str]], empty_msg: str) -> str:
    """Render rows."""
    if not rows:
        return f"<tr><td colspan='2'>{html.escape(empty_msg)}</td></tr>"
    out: list[str] = []
    for name, rel, href in rows:
        out.append(
            "<tr>"
            f"<td class='indexkey'>{render_link(name, href)}</td>"
            f"<td class='indexvalue'><code>{html.escape(rel)}</code></td>"
            "</tr>"
        )
    return "\n".join(out)


def section_table(title: str, rows_html: str) -> str:
    """Perform section table."""
    return (
        f"<h2>{html.escape(title)}</h2>\n"
        "<table class='doxtable'>\n"
        "<thead><tr><th>Name</th><th>Location</th></tr></thead>\n"
        f"<tbody>\n{rows_html}\n</tbody>\n"
        "</table>\n"
    )


def render_page(title: str, intro: str, body_html: str) -> str:
    """Render page."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>PICurv: {html.escape(title)}</title>
  <link href="doxygen.css" rel="stylesheet" />
  <link href="custom.css" rel="stylesheet" />
  <script type="text/javascript" src="theme-sync.js"></script>
</head>
<body>
  <div class="header">
    <div class="headertitle"><div class="title">{html.escape(title)}</div></div>
  </div>
  <div class="contents">
    <p>{html.escape(intro)}</p>
{body_html}
    <p>See <a href="Documentation_Map.html">Documentation Map</a> for structural navigation.</p>
  </div>
</body>
</html>
"""


def write_structured_file_index(repo_root: Path, html_dir: Path) -> None:
    """Write structured file index."""
    headers = collect_file_rows(repo_root, html_dir, "include", HEADER_SUFFIXES)
    sources = collect_file_rows(repo_root, html_dir, "src", SOURCE_SUFFIXES)
    scripts = collect_file_rows(repo_root, html_dir, "scripts", SCRIPT_SUFFIXES)
    body = (
        section_table("Header Files", render_rows(headers, "No header files found."))
        + section_table("Source Files", render_rows(sources, "No source files found."))
        + section_table("Scripts", render_rows(scripts, "No script files found."))
    )
    out = html_dir / "files_structured.html"
    out.write_text(
        render_page(
            "File List (Structured)",
            "Organized by file role: headers, source files, and scripts.",
            body,
        ),
        encoding="utf-8",
    )
    print(f"[index] wrote {out}")


def write_structured_struct_index(repo_root: Path, html_dir: Path) -> None:
    """Write structured struct index."""
    rows = collect_struct_rows(repo_root, html_dir)
    grouped: dict[str, list[tuple[str, str, str]]] = {}
    for row in rows:
        grouped.setdefault(categorize_struct(row[0]), []).append(row)

    ordered_sections = [
        "Runtime Control and Solver Orchestration",
        "Grid and Geometry",
        "Boundary Condition System",
        "Particle Transport and Statistics",
        "Immersed Boundary and FSI",
        "I/O and Postprocessing",
        "Generic Containers and Utilities",
    ]
    body_parts: list[str] = []
    for section in ordered_sections:
        body_parts.append(
            section_table(
                section,
                render_rows(grouped.get(section, []), f"No structures found for section: {section}."),
            )
        )
    out = html_dir / "annotated_structured.html"
    out.write_text(
        render_page(
            "Data Structures (By Module)",
            "Grouped by major solver modules and responsibilities.",
            "\n".join(body_parts),
        ),
        encoding="utf-8",
    )
    print(f"[index] wrote {out}")


def write_fallback_files_page(repo_root: Path, html_dir: Path) -> None:
    """Write fallback files page."""
    rows = collect_all_source_like_files(repo_root, html_dir)
    body = section_table("Files", render_rows(rows, "No source-like files found."))
    out = html_dir / "files.html"
    out.write_text(
        render_page(
            "File List",
            "Fallback file list generated from include/src/scripts.",
            body,
        ),
        encoding="utf-8",
    )
    print(f"[fallback] wrote {out}")


def write_fallback_struct_page(repo_root: Path, html_dir: Path) -> None:
    """Write fallback struct page."""
    rows = collect_struct_rows(repo_root, html_dir)
    body = section_table("Data Structures", render_rows(rows, "No C struct declarations found."))
    out = html_dir / "annotated.html"
    out.write_text(
        render_page(
            "Data Structures",
            "Fallback structure list generated from C headers.",
            body,
        ),
        encoding="utf-8",
    )
    print(f"[fallback] wrote {out}")


def main() -> int:
    """Entry point for this script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--html-dir", required=True, type=Path)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    html_dir = args.html_dir.resolve()

    write_structured_file_index(repo_root, html_dir)
    write_structured_struct_index(repo_root, html_dir)

    files_page = html_dir / "files.html"
    structs_page = html_dir / "annotated.html"
    if needs_files_fallback(files_page):
        write_fallback_files_page(repo_root, html_dir)
    if needs_structs_fallback(structs_page):
        write_fallback_struct_page(repo_root, html_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
