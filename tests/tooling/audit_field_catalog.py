#!/usr/bin/env python3
"""Enforce that the documented field catalog matches the compiled one."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
EULERIAN_CATALOG = REPO_ROOT / "src" / "field_catalog.c"
PARTICLE_CATALOG = REPO_ROOT / "src" / "particle_field_catalog.c"
LAYOUT_HEADER = REPO_ROOT / "include" / "field_catalog.h"
PAGE = REPO_ROOT / "docs" / "pages" / "56_Field_Identity_and_Layout_Catalog.md"

# Section 6 groups the catalog by layout; each bullet label names one enum value.
GROUP_LAYOUTS = {
    "node-centered": "FIELD_LAYOUT_NODE_CENTERED",
    "shifted cell-centered": "FIELD_LAYOUT_CELL_CENTERED",
    "component-staggered": "FIELD_LAYOUT_COMPONENT_STAGGERED",
    "I-face family": "FIELD_LAYOUT_I_FACE",
    "J-face family": "FIELD_LAYOUT_J_FACE",
    "K-face family": "FIELD_LAYOUT_K_FACE",
}

ENTRY_RE = re.compile(
    r'FIELD_(?:COORDINATE_)?ENTRY\(\s*FIELD_ID_\w+\s*,\s*"([^"]+)"'
    r'[^)]*?(FIELD_LAYOUT_\w+)',
    re.S,
)
PARTICLE_ENTRY_RE = re.compile(r'PARTICLE_FIELD_ENTRY\(\s*PARTICLE_FIELD_ID_\w+\s*,\s*"([^"]+)"')
LAYOUT_ENUM_RE = re.compile(r"^\s*(FIELD_LAYOUT_\w+)", re.M)
BACKTICKED_RE = re.compile(r"`([A-Za-z_][A-Za-z0-9_]*)`")


def compiled_eulerian() -> dict:
    """!
    @brief Canonical name to layout for every compiled Eulerian catalog entry.
    @return Mapping of field name to `FIELD_LAYOUT_*` value.
    """
    return dict(ENTRY_RE.findall(EULERIAN_CATALOG.read_text(encoding="utf-8")))


def compiled_particle() -> set:
    """!
    @brief Canonical names of every compiled particle catalog entry.
    @return Set of particle field names.
    """
    return set(PARTICLE_ENTRY_RE.findall(PARTICLE_CATALOG.read_text(encoding="utf-8")))


def compiled_layouts() -> set:
    """!
    @brief Layout vocabulary the header defines.
    @return Set of `FIELD_LAYOUT_*` enumerators.
    """
    text = LAYOUT_HEADER.read_text(encoding="utf-8")
    body = text[text.index("FIELD_LAYOUT_NODE_CENTERED"):]
    return set(LAYOUT_ENUM_RE.findall(body[: body.index("} FieldLayout;")]))


def page_section(heading: str) -> str:
    """!
    @brief One `@section` body from the catalog page.
    @param[in] heading Section id to extract.
    @return Section text up to the next section.
    """
    text = PAGE.read_text(encoding="utf-8")
    start = text.index(heading)
    following = text.find("@section", start + len(heading))
    return text[start: following if following != -1 else len(text)]


def documented_eulerian(section: str) -> tuple:
    """!
    @brief Field-to-layout mapping the inventory section publishes.
    @param[in] section Inventory section body.
    @return Tuple of the mapping and any structural problems found.
    """
    problems: list = []
    block = re.search(r"((?:^- .*(?:\n(?!\n)(?!- ).*)*\n)+)", section, re.M)
    if not block:
        return {}, ["section 6 no longer contains a layout-grouped list"]
    documented: dict = {}
    for item in re.findall(r"^- (.*?)(?=^- |\Z)", block.group(1), re.S | re.M):
        label, _, body = item.partition(":")
        layout = GROUP_LAYOUTS.get(label.strip())
        if layout is None:
            problems.append(f"section 6 lists an unrecognized layout group '{label.strip()}'")
            continue
        for name in BACKTICKED_RE.findall(body):
            documented[name] = layout
    return documented, problems


def documented_particle(section: str) -> set:
    """!
    @brief Particle field names the inventory section publishes.
    @param[in] section Inventory section body.
    @return Set of documented particle field names.
    """
    match = re.search(r"The particle inventory is:(.*?)\.\s", section, re.S)
    return set(BACKTICKED_RE.findall(match.group(1))) if match else set()


def documented_layouts() -> set:
    """!
    @brief Layout vocabulary the layout section tabulates.
    @return Set of documented `FIELD_LAYOUT_*` enumerators.
    """
    section = page_section("@section p56_layout_sec")
    return set(re.findall(r"\|\s*`(FIELD_LAYOUT_\w+)`\s*\|", section))


def main() -> int:
    """!
    @brief Fail when the published catalog inventory and the compiled catalog disagree.

    @details Checks identity sets, layout assignment, and layout vocabulary. It does not
             check degrees of freedom, DM family, synchronization class, availability,
             capability flags, or whether the runtime actually allocates a field.
    @return Process status code.
    """
    code = compiled_eulerian()
    section = page_section("@section p56_inventory_sec")
    documented, problems = documented_eulerian(section)

    for name in sorted(set(code) - set(documented)):
        problems.append(f"{name}: in the compiled catalog, absent from section 6")
    for name in sorted(set(documented) - set(code)):
        problems.append(f"{name}: listed in section 6, absent from the compiled catalog")
    for name in sorted(set(code) & set(documented)):
        if code[name] != documented[name]:
            problems.append(
                f"{name}: compiled as {code[name]}, documented under {documented[name]}"
            )

    header_layouts, page_layouts = compiled_layouts(), documented_layouts()
    for layout in sorted(header_layouts - page_layouts):
        problems.append(f"{layout}: defined in the header, absent from the section 4 table")
    for layout in sorted(page_layouts - header_layouts):
        problems.append(f"{layout}: tabulated in section 4, absent from the header")

    particles_code, particles_page = compiled_particle(), documented_particle(section)
    for name in sorted(particles_code - particles_page):
        problems.append(f"{name}: in the compiled particle catalog, absent from section 6")
    for name in sorted(particles_page - particles_code):
        problems.append(f"{name}: listed as a particle field, absent from the compiled catalog")

    if problems:
        print("Field catalog documentation does not match the compiled catalog:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        print(
            "\nUpdate docs/pages/56_Field_Identity_and_Layout_Catalog.md section 6 to match\n"
            "src/field_catalog.c and src/particle_field_catalog.c. The compiled catalog is\n"
            "authoritative; the page is the published record of it.",
            file=sys.stderr,
        )
        return 1

    print(
        f"Field catalog audit passed: {len(code)} Eulerian and {len(particles_code)} particle "
        f"identities, {len(header_layouts)} layouts. Degrees of freedom, DM family, "
        "synchronization class, availability, capabilities, and runtime allocation are "
        "not checked."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
