#!/usr/bin/env python3
"""Generate non-empty fallback pages for Doxygen file/structure indices."""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

CODE_SUFFIXES = {".h", ".hpp", ".c", ".cc", ".cpp"}

TYPEDEF_STRUCT_RE = re.compile(
    r"typedef\s+struct(?:\s+([A-Za-z_]\w*))?\s*\{.*?\}\s*([A-Za-z_]\w*)\s*;",
    re.S,
)
NAMED_STRUCT_RE = re.compile(r"\bstruct\s+([A-Za-z_]\w*)\s*\{")


def doxygen_file_page(filename: str) -> str:
    """Map `grid.h` -> `grid_8h.html`, `BC_Handlers.h` -> `BC__Handlers_8h.html`."""
    return filename.replace("_", "__").replace(".", "_8") + ".html"


def needs_files_fallback(path: Path) -> bool:
    if not path.exists():
        return True
    text = path.read_text(encoding="utf-8", errors="ignore")
    return "Detailed file index was not generated in this build." in text


def needs_structs_fallback(path: Path) -> bool:
    if not path.exists():
        return True
    text = path.read_text(encoding="utf-8", errors="ignore")
    return "Detailed structure index was not generated in this build." in text


def collect_source_files(repo_root: Path, html_dir: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    for subdir in ("include", "src"):
        base = repo_root / subdir
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.suffix.lower() not in CODE_SUFFIXES:
                continue
            rel = path.relative_to(repo_root).as_posix()
            page = doxygen_file_page(path.name)
            href = page if (html_dir / page).exists() else ""
            rows.append((path.name, rel, href))
    return rows


def collect_structs(repo_root: Path, html_dir: Path) -> list[tuple[str, str, str]]:
    struct_to_header: dict[str, str] = {}
    include_dir = repo_root / "include"
    if not include_dir.exists():
        return []

    for header in sorted(include_dir.rglob("*.h")):
        text = header.read_text(encoding="utf-8", errors="ignore")
        header_rel = header.relative_to(repo_root).as_posix()

        names: set[str] = set()
        for match in NAMED_STRUCT_RE.finditer(text):
            names.add(match.group(1))
        for match in TYPEDEF_STRUCT_RE.finditer(text):
            tag_name = match.group(1)
            alias_name = match.group(2)
            if tag_name:
                names.add(tag_name)
            names.add(alias_name)

        for name in names:
            struct_to_header.setdefault(name, header_rel)

    rows: list[tuple[str, str, str]] = []
    for name, header_rel in sorted(struct_to_header.items(), key=lambda item: item[0].lower()):
        struct_page = f"struct{name}.html"
        if (html_dir / struct_page).exists():
            href = struct_page
        else:
            header_page = doxygen_file_page(Path(header_rel).name)
            href = header_page if (html_dir / header_page).exists() else ""
        rows.append((name, header_rel, href))
    return rows


def render_rows_file(rows: list[tuple[str, str, str]]) -> str:
    if not rows:
        return "<tr><td colspan='2'>No source files were detected in <code>include/</code> or <code>src/</code>.</td></tr>"
    rendered = []
    for name, rel, href in rows:
        name_html = html.escape(name)
        rel_html = html.escape(rel)
        if href:
            name_html = f"<a class='el' href='{html.escape(href)}'>{name_html}</a>"
        rendered.append(
            f"<tr><td class='indexkey'>{name_html}</td><td class='indexvalue'><code>{rel_html}</code></td></tr>"
        )
    return "\n".join(rendered)


def render_rows_struct(rows: list[tuple[str, str, str]]) -> str:
    if not rows:
        return "<tr><td colspan='2'>No C struct declarations were detected under <code>include/</code>.</td></tr>"
    rendered = []
    for name, rel, href in rows:
        name_html = html.escape(name)
        rel_html = html.escape(rel)
        if href:
            name_html = f"<a class='el' href='{html.escape(href)}'>{name_html}</a>"
        rendered.append(
            f"<tr><td class='indexkey'>{name_html}</td><td class='indexvalue'><code>{rel_html}</code></td></tr>"
        )
    return "\n".join(rendered)


def render_page(title: str, subtitle: str, rows_html: str) -> str:
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
    <p>{html.escape(subtitle)}</p>
    <table class="doxtable">
      <thead>
        <tr><th>Name</th><th>Location</th></tr>
      </thead>
      <tbody>
{rows_html}
      </tbody>
    </table>
    <p>See <a href="Documentation_Catalog.html">Documentation Catalog</a> for the structural index.</p>
  </div>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--html-dir", required=True, type=Path)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    html_dir = args.html_dir.resolve()
    files_page = html_dir / "files.html"
    structs_page = html_dir / "annotated.html"

    if needs_files_fallback(files_page):
        file_rows = collect_source_files(repo_root, html_dir)
        files_page.write_text(
            render_page(
                "File List",
                "Doxygen file index fallback generated from include/src source tree.",
                render_rows_file(file_rows),
            ),
            encoding="utf-8",
        )
        print(f"[fallback] wrote {files_page}")

    if needs_structs_fallback(structs_page):
        struct_rows = collect_structs(repo_root, html_dir)
        structs_page.write_text(
            render_page(
                "Data Structures",
                "Doxygen structure index fallback generated from C struct declarations.",
                render_rows_struct(struct_rows),
            ),
            encoding="utf-8",
        )
        print(f"[fallback] wrote {structs_page}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
