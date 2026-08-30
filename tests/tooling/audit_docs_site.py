#!/usr/bin/env python3
"""Verify published-site integrity against the generated HTML, not against source declarations."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from repo_files import enumerate_repository_files  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "tests" / "tooling" / "docs_site_contract.json"
SKIP_DIRS = {".git", "docs_build", "obj", "bin", "stubs", "runs", "studies", "__pycache__", ".pytest_cache"}
SCAN_SUFFIXES = {".md", ".js", ".xml", ".html", ".yml", ".yaml"}


def load_contract() -> dict:
    """!
    @brief Load the published-site integrity contract.
    @return Parsed contract mapping.
    """
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def git_tracked_files() -> list[Path] | None:
    """!
    @brief List tracked plus non-ignored untracked files, so the scan set is reproducible
           for a commit and generated scratch is excluded without excluding tracked docs.
    @return Sorted paths, or None when git enumeration is unavailable.
    """
    return enumerate_repository_files(REPO_ROOT)


def iter_scannable():
    """!
    @brief Yield repository files whose text may carry project-owned URLs.
    @return Generator of file paths.
    """
    candidates = git_tracked_files()
    if candidates is None:
        candidates = sorted(REPO_ROOT.rglob("*"))
    for path in candidates:
        if not path.is_file() or path.suffix not in SCAN_SUFFIXES:
            continue
        if SKIP_DIRS.intersection(path.relative_to(REPO_ROOT).parts):
            continue
        yield path


def require_built_html(contract: dict) -> Path:
    """!
    @brief Resolve the generated HTML tree, failing loudly when it has not been built.
    @param[in] contract Parsed contract mapping.
    @return Path to the generated HTML directory.
    @throws RuntimeError when the publication artifact is absent.
    """
    html_dir = REPO_ROOT / contract["html_dir"]
    if not (html_dir / "index.html").is_file():
        raise RuntimeError(
            f"{contract['html_dir']}/index.html is missing. This audit validates against the built "
            "publication artifact; run 'make build-docs' first."
        )
    return html_dir


def check_forbidden_urls(contract: dict) -> list[str]:
    """!
    @brief Reject project-owned URL forms that are known not to resolve.
    @param[in] contract Parsed contract mapping.
    @return Violation lines.
    """
    violations: list[str] = []
    for rule in contract["forbidden_url_patterns"]:
        pattern = re.compile(rule["pattern"])
        for path in iter_scannable():
            for match in pattern.findall(path.read_text(encoding="utf-8", errors="replace")):
                violations.append(
                    f"{path.relative_to(REPO_ROOT)}: dead project URL '{match}'\n"
                    f"      reason: {rule['reason']}\n"
                    f"      use:    {rule['replacement_hint']}"
                )
    return violations


def check_canonical_urls(contract: dict, html_dir: Path) -> list[str]:
    """!
    @brief Verify every canonical documentation URL names a page the build actually generates.
    @param[in] contract Parsed contract mapping.
    @param[in] html_dir Generated HTML directory.
    @return Violation lines.
    """
    base = re.escape(contract["canonical_docs_base"])
    pattern = re.compile(base + r"([A-Za-z0-9_.-]+\.html)")
    violations: list[str] = []
    for path in iter_scannable():
        for page in sorted(set(pattern.findall(path.read_text(encoding="utf-8", errors="replace")))):
            if not (html_dir / page).is_file():
                violations.append(
                    f"{path.relative_to(REPO_ROOT)}: canonical URL points at '{page}', "
                    f"which the build does not generate.\n"
                    f"      The page may be excluded from Doxyfile, renamed, or never to have existed."
                )
    return violations


def check_layout(contract: dict, html_dir: Path) -> list[str]:
    """!
    @brief Verify every Doxygen layout tab resolves to a generated page.
    @param[in] contract Parsed contract mapping.
    @param[in] html_dir Generated HTML directory.
    @return Violation lines.
    """
    external = set(contract["external_layout_urls"])
    violations: list[str] = []
    for tab in ET.parse(REPO_ROOT / contract["layout_file"]).getroot().iter("tab"):
        url = tab.get("url")
        if not url or url in external or url.startswith(("http://", "https://", "/")):
            continue
        if not (html_dir / url).is_file():
            violations.append(
                f"{contract['layout_file']}: tab '{tab.get('title')}' -> {url}\n"
                f"      no such file in the generated site; the published tab will 404"
            )
    return violations


def check_orphan_pages(contract: dict, html_dir: Path) -> list[str]:
    """!
    @brief Verify every page the build publishes is reachable from the navigation.

    A page is reachable when another page adopts it with a subpage directive, or when a layout tab
    points at it directly. A page that is generated but adopted by nothing renders in the
    site with no route to it except search.
    @param[in] contract Parsed contract mapping.
    @param[in] html_dir Generated HTML directory.
    @return Violation lines.
    """
    declared: dict[str, Path] = {}
    adopted: set[str] = set()
    for directory in ("docs/pages", "docs"):
        for markdown in sorted((REPO_ROOT / directory).glob("*.md")):
            text = markdown.read_text(encoding="utf-8", errors="replace")
            match = re.search(r"@page\s+([A-Za-z0-9_]+)", text)
            if match:
                declared.setdefault(match.group(1), markdown)
            adopted.update(re.findall(r"@subpage\s+([A-Za-z0-9_]+)", text))
    tabs = {
        (tab.get("url") or "")[: -len(".html")]
        for tab in ET.parse(REPO_ROOT / contract["layout_file"]).getroot().iter("tab")
        if (tab.get("url") or "").endswith(".html")
    }
    violations = []
    for page_id, source in sorted(declared.items()):
        if not (html_dir / f"{page_id}.html").is_file():
            continue  # excluded from the build; not published, so not an orphan
        if page_id in adopted or page_id in tabs:
            continue
        violations.append(
            f"{source.relative_to(REPO_ROOT)}: page '{page_id}' is published but orphaned\n"
            f"      no page adopts it with a subpage directive and no layout tab points at it"
        )
    return violations


def check_generated_fragment_links(html_dir: Path) -> list[str]:
    """!
    @brief Verify every link emitted by a generated fragment resolves in the rendered page.

    Generated tables link into Tier-2 entries. Doxygen does not validate raw HTML
    inserted through an HTML include, so a fragment pointing at an anchor that was
    never written renders as a dead link and passes every other gate.
    @param[in] html_dir Generated HTML directory.
    @return Violation lines.
    """
    generated = REPO_ROOT / "docs" / "generated"
    if not generated.is_dir():
        return []
    violations: list[str] = []
    for fragment in sorted(generated.glob("*.html")):
        text = fragment.read_text(encoding="utf-8")
        hosts = [
            page
            for page in sorted((REPO_ROOT / "docs" / "pages").glob("*.md"))
            if fragment.name in page.read_text(encoding="utf-8")
        ]
        if not hosts:
            violations.append(f"docs/generated/{fragment.name}: generated but no page includes it")
            continue
        for host in hosts:
            match = re.search(r"@page\s+([A-Za-z0-9_]+)", host.read_text(encoding="utf-8"))
            if not match:
                continue
            rendered = html_dir / f"{match.group(1)}.html"
            if not rendered.is_file():
                continue
            markup = rendered.read_text(encoding="utf-8")
            ids = set(re.findall(r'id="([A-Za-z0-9_]+)"', markup))
            for anchor in sorted(set(re.findall(r'href="#([A-Za-z0-9_]+)"', text))):
                if anchor not in ids:
                    violations.append(
                        f"docs/generated/{fragment.name}: links to #{anchor}, which does not "
                        f"exist in the rendered {rendered.name}"
                    )
            for page_ref, anchor in sorted(set(re.findall(r'href="([A-Za-z0-9_]+)\.html#([A-Za-z0-9_]+)"', text))):
                target = html_dir / f"{page_ref}.html"
                if not target.is_file():
                    violations.append(f"docs/generated/{fragment.name}: links to {page_ref}.html, which is not generated")
                elif f'id="{anchor}"' not in target.read_text(encoding="utf-8"):
                    violations.append(
                        f"docs/generated/{fragment.name}: links to {page_ref}.html#{anchor}, "
                        f"which does not exist in the rendered page"
                    )
    return violations


def strip_non_prose(text: str) -> str:
    """!
    @brief Blank out fenced code, inline code, and HTML comments.

    @details A link shown inside an example is documentation of a link, not a link.
             Replacing the spans with blanks rather than deleting them keeps line
             numbers meaningful for any future line-accurate reporting.
    @param[in] text Raw Markdown text.
    @return Text with non-prose spans blanked.
    """
    def blank(match: "re.Match") -> str:
        """!
        @brief Replace a matched span with spaces, preserving newlines.
        @param[in] match Regular-expression match to blank out.
        @return Whitespace of the same shape as the matched text.
        """
        return re.sub(r"[^\n]", " ", match.group(0))

    text = re.sub(r"<!--.*?-->", blank, text, flags=re.S)
    text = re.sub(r"^```.*?^```", blank, text, flags=re.S | re.M)
    text = re.sub(r"~~~.*?~~~", blank, text, flags=re.S)
    text = re.sub(r"`[^`\n]*`", blank, text)
    return text


def markdown_fragment_links(text: str) -> list:
    """!
    @brief Extract Markdown links carrying a fragment, from prose only.

    @details Covers `[t](#frag)`, `[t](#frag "title")`, and `[t](page.md#frag)`, and
             accepts every legal fragment character rather than only word characters -
             a hyphenated anchor is the common Markdown style and was previously
             invisible to this check.
    @param[in] text Raw Markdown text.
    @return List of (target_page_or_None, fragment) pairs.
    """
    prose = strip_non_prose(text)
    pattern = re.compile(
        r"\]\(\s*"
        r"(?P<path>[^)\s#]*)"
        r"#(?P<fragment>[^)\s\"']+)"
        r"(?:\s+[\"'][^)]*[\"'])?"
        r"\s*\)"
    )
    results = []
    for match in pattern.finditer(prose):
        target = match.group("path") or None
        results.append((target, match.group("fragment")))
    return results


def tracked_markdown() -> list:
    """!
    @brief Every Markdown file the current commit carries.

    @details Uses the same git-backed enumeration as the link checker, so fragment
             validation covers README, example documentation, and every `guide.md`
             rather than only the Doxygen page tree.
    @return Sorted Markdown paths.
    """
    found = enumerate_repository_files(REPO_ROOT, ".md", frozenset(SKIP_DIRS))
    if found is None:
        return sorted(path for path in REPO_ROOT.rglob("*.md") if path.is_file())
    return found


def rendered_ids(path: Path) -> set:
    """!
    @brief Anchor ids present in a rendered HTML page.
    @param[in] path Rendered page.
    @return Set of ids.
    """
    if not path.is_file():
        return set()
    markup = path.read_text(encoding="utf-8", errors="replace")
    ids = set(re.findall(r'id="([^"]+)"', markup))
    ids.update(re.findall(r'name="([^"]+)"', markup))
    return ids


def heading_anchor(text: str) -> str:
    """!
    @brief GitHub-style anchor slug for a Markdown heading.
    @param[in] text Heading text without its leading hashes.
    @return Anchor slug.
    """
    slug = text.strip().lower()
    slug = re.sub(r"[`*_~\[\]()]", "", slug)
    slug = re.sub(r"[^\w\s-]", "", slug)
    return re.sub(r"\s+", "-", slug).strip("-")


def markdown_anchors(path: Path) -> set:
    """!
    @brief Anchors a plain Markdown file offers.

    @details A file that Doxygen does not render still has targets: explicit HTML
             anchors, and the heading slugs a Markdown viewer generates.
    @param[in] path Markdown file.
    @return Set of available anchor names.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    anchors = set(re.findall(r'<a\s+(?:id|name)="([^"]+)"', text))
    anchors.update(re.findall(r"^@anchor\s+(\S+)\s*$", text, re.M))
    for heading in re.findall(r"^#{1,6}\s+(.+?)\s*$", strip_non_prose(text), re.M):
        anchors.add(heading_anchor(heading))
    return anchors


def check_page_cross_references(html_dir: Path) -> list[str]:
    """!
    @brief Verify hand-written Markdown fragment links resolve, repository-wide.

    @details Doxygen validates its own reference graph, but a Markdown fragment link is
             checked by nothing. Cross-page targets are resolved **relative to the
             linking file**, so the many `guide.md` files in this repository are not
             conflated by basename.
    @param[in] html_dir Generated HTML directory.
    @return Violation lines.
    """
    page_ids: dict = {}
    for markdown in tracked_markdown():
        match = re.search(r"^@page\s+([A-Za-z0-9_]+)", markdown.read_text(encoding="utf-8"), re.M)
        if match:
            page_ids[markdown.resolve()] = match.group(1)

    def anchors_for(path: Path) -> set:
        """!
        @brief Anchors available in one Markdown file, rendered or not.
        @param[in] path Markdown file.
        @return Set of anchor names.
        """
        page_id = page_ids.get(path.resolve())
        if page_id:
            return rendered_ids(html_dir / f"{page_id}.html")
        return markdown_anchors(path)

    violations: list[str] = []
    for markdown in tracked_markdown():
        text = markdown.read_text(encoding="utf-8", errors="replace")
        for target, fragment in markdown_fragment_links(text):
            if target is None:
                resolved = markdown
            else:
                resolved = (markdown.parent / target).resolve()
                if not resolved.is_file() or resolved.suffix != ".md":
                    continue  # the link checker owns missing-file reporting
            available = anchors_for(resolved)
            if not available:
                continue  # nothing to check against; not evidence of a broken link
            if fragment not in available:
                where = markdown.name if target is None else target
                violations.append(
                    f"{markdown.relative_to(REPO_ROOT)}: fragment link '#{fragment}' does not "
                    f"resolve in {where}"
                )
    return violations


def main() -> int:
    """!
    @brief Fail when a project-owned URL is dead or a navigation tab has no generated page.
    @return Process status code.
    """
    contract = load_contract()
    try:
        html_dir = require_built_html(contract)
    except RuntimeError as error:
        print(f"Docs-site audit could not run: {error}", file=sys.stderr)
        return 1
    violations = (
        check_forbidden_urls(contract)
        + check_canonical_urls(contract, html_dir)
        + check_layout(contract, html_dir)
        + check_orphan_pages(contract, html_dir)
        + check_generated_fragment_links(html_dir)
        + check_page_cross_references(html_dir)
    )
    if violations:
        print("Published-site integrity violations:", file=sys.stderr)
        for violation in violations:
            print(f"  {violation}", file=sys.stderr)
        return 1
    print(
        f"Docs-site audit passed: canonical URLs and navigation tabs all resolve to pages "
        f"generated under {contract['html_dir']}; no published page is orphaned; generated fragment links resolve."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
