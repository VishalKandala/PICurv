#!/usr/bin/env python3
"""Extract the public capability inventory from executable sources and render it for the docs."""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = REPO_ROOT / "tests" / "tooling" / "capability_families.json"
GENERATED_DIR = REPO_ROOT / "docs" / "generated"


def load_registry() -> dict:
    """!
    @brief Load the capability family registry.
    @return Parsed registry mapping.
    """
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


def literal(node: ast.AST):
    """!
    @brief Evaluate a literal AST node, additionally accepting the bare `set()` call
           that appears in the boundary-handler specs for empty parameter sets.
    @param[in] node Parsed AST node.
    @return The Python value the node denotes.
    """
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "set":
        return set(literal(node.args[0])) if node.args else set()
    if isinstance(node, ast.Dict):
        return {literal(k): literal(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.Set):
        return {literal(e) for e in node.elts}
    if isinstance(node, (ast.List, ast.Tuple)):
        return [literal(e) for e in node.elts]
    return ast.literal_eval(node)


def module_syntax(module: str) -> ast.Module:
    """!
    @brief Parse a dotted module into one syntax tree, whether file or package.

    @details A public surface may be a package: `picurv_cli.storage` exposes its
             constants through its `__init__`, but they are defined across its
             modules. Concatenating their bodies lets the readers below stay written
             against a single tree, which is what they mean by "the module".
    @param[in] module Dotted module name.
    @return Syntax tree covering every file the module resolves to.
    """
    relative = Path(*module.split("."))
    package = REPO_ROOT / relative
    if package.is_dir():
        combined = ast.Module(body=[], type_ignores=[])
        for child in sorted(package.glob("*.py")):
            if child.name == "__init__.py":
                continue
            combined.body.extend(ast.parse(child.read_text(encoding="utf-8")).body)
        return combined
    return ast.parse((REPO_ROOT / relative.with_suffix(".py")).read_text(encoding="utf-8"))


def python_dict_values(module: str, symbol: str) -> dict[str, dict]:
    """!
    @brief Read a public-surface dictionary from the CLI package without importing PETSc.
    @param[in] module Dotted module name.
    @param[in] symbol Module-level dictionary name.
    @return Mapping of selector value to its declared parameter contract.
    """
    tree = module_syntax(module)
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if symbol not in targets:
            continue
        raw = literal(node.value)
        result: dict[str, dict] = {}
        for key, value in raw.items():
            if isinstance(value, dict):
                result[key] = {
                    "types": sorted(value.get("types", [])),
                    "required_params": sorted(value.get("required_params", [])),
                    "optional_params": sorted(value.get("optional_params", [])),
                }
            else:
                result[key] = {"maps_to": value}
        return result
    raise RuntimeError(f"{symbol} not found in {module}")


def python_normalizer_values(module: str, symbol: str) -> dict[str, dict]:
    """!
    @brief Read the canonical selector strings accepted by a normalizer function.
    @param[in] module Dotted module name.
    @param[in] symbol Normalizer function name.
    @return Mapping of canonical value to the runtime token it maps to.
    """
    tree = module_syntax(module)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == symbol:
            for inner in ast.walk(node):
                if isinstance(inner, ast.Dict) and inner.keys:
                    try:
                        mapping = ast.literal_eval(inner)
                    except ValueError:
                        continue
                    if all(isinstance(k, str) for k in mapping):
                        return {k: {"maps_to": v} for k, v in mapping.items()}
    raise RuntimeError(f"no canonical mapping found in {module}.{symbol}")


def function_body(path: str, function: str) -> str:
    """!
    @brief Return the source text of one C function, so extraction never spans the file.
    @param[in] path Repository-relative C source path.
    @param[in] function Function name to isolate.
    @return Source text of that function body.
    @throws RuntimeError when the function cannot be located.
    """
    text = (REPO_ROOT / path).read_text(encoding="utf-8")
    match = re.search(rf"^[A-Za-z_][\w \*]*\b{re.escape(function)}\s*\(", text, re.M)
    if not match:
        raise RuntimeError(f"function {function} not found in {path}")
    index = text.index("{", match.end() - 1)
    depth = 0
    for offset in range(index, len(text)):
        if text[offset] == "{":
            depth += 1
        elif text[offset] == "}":
            depth -= 1
            if depth == 0:
                return text[index : offset + 1]
    raise RuntimeError(f"unbalanced braces while reading {function} in {path}")


def python_membership_values(module: str, symbol: str) -> dict:
    """!
    @brief Read the accepted values from a normalizer that validates by set membership.

    @details Some normalizers check `if value not in {...}` rather than mapping through a
             dict. The accepted set is still the public surface, so it is extracted the
             same way rather than being hand-listed.
    @param[in] module Dotted module name.
    @param[in] symbol Normalizer function name.
    @return Mapping of accepted value to the token it resolves to.
    """
    tree = module_syntax(module)
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef) and node.name == symbol):
            continue
        collected: dict = {}
        for inner in ast.walk(node):
            if isinstance(inner, ast.Set):
                try:
                    members = {literal(element) for element in inner.elts}
                except (ValueError, TypeError):
                    continue
                if all(isinstance(member, str) for member in members):
                    for member in members:
                        collected.setdefault(member, {"maps_to": member})
        if collected:
            return collected
    raise RuntimeError(f"no membership set found in {module}.{symbol}")


def python_equality_chain_values(module: str, symbol: str) -> dict:
    """!
    @brief Read accepted values from a normalizer that compares against string literals.

    @details A third normalizer shape: `if normalized == "ucat": ... elif ... == "ucont"`.
             The compared literals are the public surface, so they are extracted rather
             than hand-listed, keeping the inventory tied to the code.
    @param[in] module Dotted module name.
    @param[in] symbol Normalizer function name.
    @return Mapping of accepted value to itself.
    """
    tree = module_syntax(module)
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef) and node.name == symbol):
            continue
        collected: dict = {}
        for inner in ast.walk(node):
            if not isinstance(inner, ast.Compare):
                continue
            if not any(isinstance(op, (ast.Eq, ast.NotEq)) for op in inner.ops):
                continue
            for comparator in inner.comparators:
                if isinstance(comparator, ast.Constant) and isinstance(comparator.value, str):
                    if comparator.value:
                        collected.setdefault(comparator.value, {"maps_to": comparator.value})
        if collected:
            return collected
    raise RuntimeError(f"no equality chain found in {module}.{symbol}")


def c_string_map_values(path: str, function: str) -> dict[str, str]:
    """!
    @brief Extract the case-insensitive selector strings a C parser accepts.
    @param[in] path Repository-relative C source path.
    @param[in] function Parser function name.
    @return Mapping of accepted selector string to the enum constant it selects.
    """
    body = function_body(path, function)
    pairs = re.findall(
        r'strcasecmp\(\s*str\s*,\s*"([^"]+)"\s*\)\s*==\s*0\s*\)\s*\*handler_out\s*=\s*([A-Z][A-Z0-9_]+)',
        body,
    )
    return dict(pairs)


def c_token_map_values(path: str, function: str, variable: str, field: str) -> dict[str, str]:
    """!
    @brief Extract an exact-match token chain that assigns an enum to a context field.

    Handles the `strcmp(buf, "TOKEN") == 0 ... field = ENUM;` shape used for
    generated PETSc option tokens, including chains where several tokens share one
    assignment (an alias arm).
    @param[in] path Repository-relative C source path.
    @param[in] function Enclosing function name.
    @param[in] variable Name of the char buffer holding the option value.
    @param[in] field Assigned context field, for example `mom_solver_type`.
    @return Mapping of accepted token to the enum constant it selects.
    """
    body = function_body(path, function)
    mapping: dict[str, str] = {}
    pattern = re.compile(
        r'((?:strcmp\(\s*' + re.escape(variable) + r'\s*,\s*"[^"]+"\s*\)\s*==\s*0\s*\|?\|?\s*)+)'
        r"[^{]*\{[^}]*?\b" + re.escape(field) + r"\s*=\s*([A-Z][A-Z0-9_]+)\s*;",
        re.S,
    )
    for arm, enum in pattern.findall(body):
        for token in re.findall(r'"([^"]+)"', arm):
            mapping[token] = enum
    return mapping


def c_switch_values(path: str, function: str, prefix: str) -> set[str]:
    """!
    @brief Extract the enum constants a C factory switch dispatches on.
    @param[in] path Repository-relative C source path.
    @param[in] function Enclosing function name.
    @param[in] prefix Enum constant prefix.
    @return Set of dispatched enum constants.
    """
    body = function_body(path, function)
    return set(re.findall(rf"case\s+({re.escape(prefix)}[A-Z0-9_]+)\s*:", body))


def c_dispatch_values(path: str, function: str, prefix: str) -> set[str]:
    """!
    @brief Extract the enum constants an if/else dispatch chain compares against.
    @param[in] path Repository-relative C source path.
    @param[in] function Enclosing function name.
    @param[in] prefix Enum constant prefix.
    @return Set of enum constants the dispatch acts on.
    """
    body = function_body(path, function)
    return set(re.findall(rf"==\s*({re.escape(prefix)}[A-Z0-9_]+)", body))


def c_enum_values(path: str, symbol: str) -> set[str]:
    """!
    @brief Extract the members of a C enum.
    @param[in] path Repository-relative header path.
    @param[in] symbol Enum type name.
    @return Set of enum member names.
    """
    text = (REPO_ROOT / path).read_text(encoding="utf-8")
    match = re.search(r"typedef\s+enum\s*\{([^{}]*)\}\s*" + re.escape(symbol) + r"\s*;", text, re.S)
    if not match:
        raise RuntimeError(f"enum {symbol} not found in {path}")
    # Split on commas rather than requiring one member per line: a single-line enum
    # would otherwise be silently under-extracted, which is the failure class these
    # extractors exist to eliminate.
    members = set()
    body = re.sub(r"/\*.*?\*/", "", match.group(1), flags=re.S)
    body = re.sub(r"//[^\n]*", "", body)
    for chunk in body.split(","):
        name = re.match(r"\s*([A-Z][A-Z0-9_]*)\s*(?:=|$)", chunk)
        if name:
            members.add(name.group(1))
    return members


def python_constant_values(module: str, symbol: str) -> dict:
    """!
    @brief Read the accepted values from a named module-level choice set.

    @details The preferred shape. A choice set written as an inline literal at its point
             of use is invisible to the census, so the rule is that it must be a named
             module-level constant - a tuple, list, set, or dict of strings. A dict maps
             each accepted spelling to what it resolves to; a sequence maps each value
             to itself.
    @param[in] module Dotted module name.
    @param[in] symbol Constant name.
    @return Mapping of accepted value to the token it resolves to.
    """
    tree = module_syntax(module)
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == symbol
                   for target in node.targets):
            continue
        value = literal(node.value)
        if isinstance(value, dict):
            if not all(isinstance(k, str) for k in value):
                raise RuntimeError(f"{module}.{symbol} is not keyed by strings")
            # A dict of strings declares spellings and what each resolves to. A dict
            # whose values are anything else - a per-value specification, say - declares
            # the keys as the choice set and carries its own detail alongside.
            if all(isinstance(v, str) for v in value.values()):
                return {k: {"maps_to": v} for k, v in value.items()}
            return {k: {"maps_to": k} for k in value}
        if isinstance(value, (tuple, list, set, frozenset)):
            members = list(value)
            if not all(isinstance(member, str) for member in members):
                raise RuntimeError(f"{module}.{symbol} is not a sequence of strings")
            return {member: {"maps_to": member} for member in members}
        raise RuntimeError(f"{module}.{symbol} is not a choice set")
    raise RuntimeError(f"no module-level constant {symbol} in {module}")


def collect(family: dict) -> dict:
    """!
    @brief Build one family's inventory record from its declared sources.
    @param[in] family Family registry entry.
    @return Inventory record for the family.
    """
    surface = family["public_surface"]
    if surface["kind"] == "python_dict":
        values = python_dict_values(surface["module"], surface["symbol"])
    elif surface["kind"] == "python_normalizer":
        values = python_normalizer_values(surface["module"], surface["symbol"])
    elif surface["kind"] == "python_membership":
        values = python_membership_values(surface["module"], surface["symbol"])
    elif surface["kind"] == "python_equality_chain":
        values = python_equality_chain_values(surface["module"], surface["symbol"])
    elif surface["kind"] == "python_constant":
        values = python_constant_values(surface["module"], surface["symbol"])
    else:
        raise RuntimeError(f"unknown public_surface kind: {surface['kind']}")

    parity = []
    for source in family.get("parity_sources", []):
        kind = source["kind"]
        if kind == "c_string_map":
            found = c_string_map_values(source["path"], source["function"])
        elif kind == "c_token_map":
            found = c_token_map_values(
                source["path"], source["function"], source["variable"], source["field"]
            )
        elif kind == "c_switch":
            found = c_switch_values(source["path"], source["function"], source["prefix"])
        elif kind == "c_dispatch":
            found = c_dispatch_values(source["path"], source["function"], source["prefix"])
        elif kind == "c_enum":
            found = c_enum_values(source["path"], source["symbol"])
        else:
            raise RuntimeError(f"unknown parity source kind: {kind}")
        if isinstance(found, dict):
            parity.append({"source": source, "values": sorted(found), "mapping": found})
        else:
            parity.append({"source": source, "values": sorted(found)})

    return {
        "id": family["id"],
        "title": family["title"],
        "selector": family["selector"],
        "family_page": family["family_page"],
        "public_values": values,
        "parity": parity,
    }


def apply_metadata(inventory: list[dict], registry: dict) -> None:
    """!
    @brief Merge declared per-value metadata (status, alias target) into the inventory.
    @param[in,out] inventory Collected family records.
    @param[in] registry Parsed registry mapping.
    @return None.
    """
    entries = {entry["id"]: entry for entry in registry["families"]}
    for family in inventory:
        declared = entries[family["id"]].get("value_metadata", {})
        for name, spec in family["public_values"].items():
            meta = declared.get(name, {})
            if meta.get("alias_of"):
                spec["alias_of"] = meta["alias_of"]
            if meta.get("spelling_of"):
                spec["spelling_of"] = meta["spelling_of"]
            # A spelling has no status of its own: it inherits the canonical value's,
            # so a defective capability cannot look supported under another name. A
            # deprecated alias keeps its own status - being retired is a fact about the
            # alias, not about what it resolves to.
            spelling_target = meta.get("spelling_of")
            if spelling_target:
                spec["status"] = declared.get(spelling_target, {}).get("status", "supported")
            elif meta.get("status"):
                spec["status"] = meta["status"]


def classify(inventory: list[dict]) -> dict[str, int]:
    """!
    @brief Count selectable, alias, and latent values separately.

    A single total conflates three different things: what a user can choose, what is
    only kept readable for old configs, and what is declared but unreachable.
    @param[in] inventory Collected family records.
    @return Mapping of category to count.
    """
    counts = {"selectable": 0, "spelling": 0, "alias": 0, "latent": 0}
    for family in inventory:
        for spec in family["public_values"].values():
            if spec.get("reachability") == "latent":
                counts["latent"] += 1
            elif spec.get("alias_of"):
                counts["alias"] += 1
            elif spec.get("spelling_of"):
                counts["spelling"] += 1
            else:
                counts["selectable"] += 1
    return counts


def apply_reachability(inventory: list[dict], registry: dict) -> None:
    """!
    @brief Mark declared values that no other family can actually satisfy as latent.

    A boundary type is only selectable if some public handler accepts it. Listing a
    type no handler supports advertises a capability every complete configuration
    would be rejected for.
    @param[in,out] inventory Collected family records.
    @param[in] registry Parsed registry mapping.
    @return None.
    """
    by_id = {family["id"]: family for family in inventory}
    for entry in registry["families"]:
        spec = entry.get("reachable_from")
        if not spec:
            continue
        provider = by_id.get(spec["family"])
        target = by_id.get(entry["id"])
        if provider is None or target is None:
            raise RuntimeError(f"reachable_from names an unknown family: {spec['family']}")
        reachable: set[str] = set()
        for value in provider["public_values"].values():
            reachable.update(value.get(spec["field"], []))
        for name, value in target["public_values"].items():
            resolved = value.get("maps_to", name)
            # Reachability is an independent axis. Overwriting `status` here destroyed
            # the declared lifecycle, which silently disabled lifecycle enforcement for
            # every family using reachable_from.
            is_reachable = resolved in reachable
            value["reachable"] = is_reachable
            value["reachability"] = "reachable" if is_reachable else "latent"


def entry_anchor(registry_entry: dict, value: str) -> str:
    """!
    @brief Anchor name of the Tier-2 entry for one selector value.
    @param[in] registry_entry Registry entry carrying the anchor prefix.
    @param[in] value Public selector value.
    @return Anchor name.
    """
    return registry_entry["entry_anchor_prefix"] + re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def html_escape(text: str) -> str:
    """!
    @brief Escape text for inclusion in generated HTML.
    @param[in] text Raw text.
    @return Escaped text.
    """
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def documented_entries(family: dict, registry_entry: dict) -> set:
    """!
    @brief Values whose Tier-2 entry anchor is actually present on the family page.

    Generated tables must not link to an entry that does not exist: a deferred or
    latent value has no anchor, and a dead in-page link is worse than plain text.
    @param[in] family Collected family record.
    @param[in] registry_entry Registry entry naming the family page.
    @return Set of value names that have an entry.
    """
    page = REPO_ROOT / "docs" / "pages" / f"{registry_entry['family_page']}.md"
    if not page.is_file():
        return set()
    text = page.read_text(encoding="utf-8")
    anchors = set(re.findall(r"^@anchor\s+([A-Za-z0-9_]+)\s*$", text, re.M))
    return {
        name
        for name in family["public_values"]
        if entry_anchor(registry_entry, name) in anchors
    }


def render_family(family: dict, registry_entry: dict) -> str:
    """!
    @brief Render one family's value table as a Doxygen-includable HTML fragment.

    HTML rather than Markdown because Doxygen's plain include command inserts Markdown
    verbatim as a code block, while its HTML include command inserts real markup. Each value links to its
    Tier-2 entry so the inventory is a route into the documentation, not a dead list.
    @param[in] family Collected family record.
    @param[in] registry_entry Registry entry for the same family.
    @return HTML fragment text.
    """
    values = family["public_values"]
    documented_values = documented_entries(family, registry_entry)
    has_params = any("required_params" in v for v in values.values())
    latent_present = any(
        v.get("reachability") == "latent"
        or v.get("alias_of")
        or v.get("spelling_of")
        or v.get("status", "supported") != "supported"
        for v in values.values()
    )

    out = [
        f"<!-- GENERATED FILE - do not edit by hand. Family: {family['id']}.",
        "     Regenerate with: make docs-inventory -->",
        '<table class="markdownTable">',
        "<tr>",
    ]
    headers = ["Value", "Applies to", "Required parameters", "Optional parameters"] if has_params \
        else ["Value", "Maps to"]
    if latent_present:
        headers.append("Status")
    out += [f'<th class="markdownTableHeadNone">{h}</th>' for h in headers]
    out.append("</tr>")

    for name, spec in sorted(values.items()):
        cells = []
        label = f"<code>{html_escape(name)}</code>"
        # A spelling has no entry of its own: link it to the canonical value's entry.
        # A value whose entry is absent (latent, or deliberately deferred) must not
        # link at all rather than link into nothing.
        target = spec.get("spelling_of") or name
        documented = target in documented_values
        if documented and spec.get("reachability") != "latent":
            label = f'<a href="#{entry_anchor(registry_entry, target)}">{label}</a>'
        cells.append(label)
        if has_params:
            cells.append(", ".join(f"<code>{html_escape(x)}</code>" for x in spec.get("types", [])) or "-")
            cells.append(
                ", ".join(f"<code>{html_escape(x)}</code>" for x in spec.get("required_params", [])) or "none"
            )
            cells.append(
                ", ".join(f"<code>{html_escape(x)}</code>" for x in spec.get("optional_params", [])) or "none"
            )
        else:
            cells.append(f"<code>{html_escape(str(spec.get('maps_to', '-')))}</code>")
        if latent_present:
            if spec.get("reachability") == "latent":
                cells.append("<b>latent - not selectable</b>")
            elif spec.get("alias_of"):
                cells.append(
                    "<b>deprecated</b> - alias of <code>"
                    + html_escape(str(spec["alias_of"]))
                    + "</code>"
                )
            elif spec.get("spelling_of"):
                cells.append(
                    "accepted spelling of <code>"
                    + html_escape(str(spec["spelling_of"]))
                    + "</code>"
                )
            else:
                # Show the lifecycle status, not merely that the value can be typed. A
                # known-defective capability reading "selectable" is the failure this
                # column exists to prevent.
                status = spec.get("status", "supported")
                cells.append(
                    f"<b>{html_escape(status)}</b>" if status != "supported" else "supported"
                )
        out.append("<tr>" + "".join(f'<td class="markdownTableBodyNone">{c}</td>' for c in cells) + "</tr>")

    out.append("</table>")
    return "\n".join(out) + "\n"


def family_fragment_path(family_id: str) -> Path:
    """!
    @brief Path of the per-family includable fragment.
    @param[in] family_id Family identifier.
    @return Fragment path under the generated directory.
    """
    return GENERATED_DIR / f"capability_inventory_{family_id.replace('.', '_')}.html"


def render_evidence_matrix(inventory: list[dict], registry: dict) -> str:
    """!
    @brief Render the project-wide capability-by-evidence matrix as an HTML fragment.

    A scientist deciding whether a result is credible needs to see, in one place,
    what confidence the project claims for each capability. An empty row is a real
    answer - it says "implemented only".
    @param[in] inventory Collected family records.
    @param[in] registry Parsed registry mapping.
    @return HTML fragment text.
    """
    facets = registry["evidence_facets"]
    order = ["unit", "integration", "analytical", "benchmark", "reference", "production"]
    entries = {e["id"]: e for e in registry["families"]}
    out = [
        "<!-- GENERATED FILE - do not edit by hand. Regenerate with: make docs-inventory -->",
        '<table class="markdownTable">',
        "<tr>",
        '<th class="markdownTableHeadNone">Capability</th>',
        '<th class="markdownTableHeadNone">Family</th>',
    ]
    out += [f'<th class="markdownTableHeadNone">{html_escape(facets[f])}</th>' for f in order]
    out.append("</tr>")
    for family in inventory:
        meta = entries[family["id"]].get("value_metadata", {})
        for name, spec in sorted(family["public_values"].items()):
            if spec.get("spelling_of") or spec.get("reachability") == "latent":
                continue
            have = dict(meta.get(name, {}).get("evidence", {}) or {})
            label = f"<code>{html_escape(name)}</code>"
            if name in documented_entries(family, entries[family["id"]]):
                anchor = entry_anchor(entries[family["id"]], name)
                label = f'<a href="{family["family_page"]}.html#{anchor}">{label}</a>'
            cells = [label, f"<code>{html_escape(family['id'])}</code>"]
            cells += [
                (
                    '<span title="'
                    + html_escape("; ".join(have[f]))
                    + '">&#10003;</span>'
                )
                if f in have
                else "&#8212;"
                for f in order
            ]
            out.append("<tr>" + "".join(f'<td class="markdownTableBodyNone">{c}</td>' for c in cells) + "</tr>")
    out.append("</table>")
    return "\n".join(out) + "\n"


def render_markdown(inventory: list[dict]) -> str:
    """!
    @brief Render the inventory as a Doxygen-includable Markdown fragment.
    @param[in] inventory Collected family records.
    @return Markdown text.
    """
    lines = [
        "<!-- GENERATED FILE - do not edit by hand.",
        "     Regenerate with: make docs-inventory",
        "     Source of truth: the Python validation layer named per family below. -->",
        "",
    ]
    for family in inventory:
        lines.append(f"### {family['title']}")
        lines.append("")
        lines.append(f"Selector: `{family['selector']}`")
        lines.append("")
        has_params = any("required_params" in v for v in family["public_values"].values())
        if has_params:
            lines.append("| Value | Applies to | Required parameters | Optional parameters |")
            lines.append("|---|---|---|---|")
            for name, spec in sorted(family["public_values"].items()):
                types = ", ".join(f"`{t}`" for t in spec.get("types", [])) or "-"
                req = ", ".join(f"`{p}`" for p in spec.get("required_params", [])) or "none"
                opt = ", ".join(f"`{p}`" for p in spec.get("optional_params", [])) or "none"
                lines.append(f"| `{name}` | {types} | {req} | {opt} |")
        else:
            lines.append("| Value | Maps to |")
            lines.append("|---|---|")
            for name, spec in sorted(family["public_values"].items()):
                lines.append(f"| `{name}` | `{spec.get('maps_to', '-')}` |")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    """!
    @brief Generate the capability inventory artifacts.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(description="Generate the public capability inventory.")
    parser.add_argument("--check", action="store_true", help="Fail if generated output is stale.")
    args = parser.parse_args()

    registry = load_registry()
    entries = {f["id"]: f for f in registry["families"]}
    inventory = [collect(family) for family in registry["families"]]
    apply_metadata(inventory, registry)
    apply_reachability(inventory, registry)
    snapshot = json.dumps(inventory, indent=2, sort_keys=True) + "\n"
    markdown = render_markdown(inventory)

    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    json_path = GENERATED_DIR / "capability_inventory.json"
    md_path = GENERATED_DIR / "capability_inventory.md"

    managed = {json_path: snapshot, md_path: markdown}
    managed[GENERATED_DIR / "evidence_matrix.html"] = render_evidence_matrix(inventory, registry)
    for family in inventory:
        managed[family_fragment_path(family["id"])] = render_family(family, entries[family["id"]])

    # Orphan detection is scoped to the files THIS generator owns. Other generators
    # write into the same directory; claiming their output would report false orphans.
    owned_prefixes = ("capability_inventory", "evidence_matrix")
    existing = {
        path
        for path in GENERATED_DIR.iterdir()
        if path.is_file() and path.name.startswith(owned_prefixes)
    }
    orphans = sorted(existing - set(managed))

    if args.check:
        problems = [
            f"stale: {path.relative_to(REPO_ROOT)}"
            for path, content in managed.items()
            if not path.is_file() or path.read_text(encoding="utf-8") != content
        ]
        problems += [f"orphan: {path.relative_to(REPO_ROOT)} (no family produces it)" for path in orphans]
        if problems:
            print("Generated capability inventory is out of date:", file=sys.stderr)
            for problem in problems:
                print(f"  {problem}", file=sys.stderr)
            print("\nRegenerate with: make docs-inventory", file=sys.stderr)
            return 1
        print(f"Capability inventory is current ({len(inventory)} families, {len(managed)} managed files).")
        return 0

    for path, content in managed.items():
        path.write_text(content, encoding="utf-8")
    for path in orphans:
        path.unlink()
        print(f"Removed orphaned generated file: {path.relative_to(REPO_ROOT)}")
    counts = classify(inventory)
    print(
        f"Wrote capability inventory: {len(inventory)} families, "
        f"{counts['selectable']} canonical, {counts['spelling']} accepted spelling, "
        f"{counts['alias']} deprecated alias, {counts['latent']} latent; "
        f"{len(managed)} managed files."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
