#!/usr/bin/env python3
"""Generate the CLI reference by introspecting the fully assembled argparse parser."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
GENERATED_DIR = REPO_ROOT / "docs" / "generated"
HTML_PATH = GENERATED_DIR / "cli_reference.html"
JSON_PATH = GENERATED_DIR / "cli_reference.json"


def build_parser():
    """!
    @brief Build the real parser, including registrars delegated from other modules.

    @details Regex over `cli.py` would miss `add_storage_parser()`, which contributes a
             whole top-level command from `storage.py`. Only the assembled parser knows
             the true command set.
    @return The assembled argparse parser.
    """
    sys.path.insert(0, str(REPO_ROOT))
    from picurv_cli.cli import build_main_parser

    return build_main_parser()


def escape(text) -> str:
    """!
    @brief Escape text for generated HTML.
    @param[in] text Value to escape.
    @return Escaped string.
    """
    return (
        str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def describe_action(action) -> dict:
    """!
    @brief Capture the user-relevant contract of one argparse action.
    @param[in] action Parser action.
    @return Mapping describing the action.
    """
    return {
        "flags": list(action.option_strings),
        "dest": action.dest,
        "required": bool(getattr(action, "required", False)),
        "choices": sorted(str(choice) for choice in action.choices) if action.choices else [],
        "default": None if action.default in (None, argparse.SUPPRESS) else str(action.default),
        "nargs": None if action.nargs is None else str(action.nargs),
        "help": (action.help or "").strip(),
        "is_flag": isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)),
    }


#: argparse's own title for the group of options nobody explicitly grouped. Python
#: renamed it from "optional arguments" to "options" in 3.10 (bpo-9694); a group's
#: title is otherwise a literal string this repo wrote, so only this one default
#: needs normalizing. Left un-normalized, the generated reference would flip between
#: spellings depending on which Python happened to run the generator, making `--check`
#: fail on any interpreter that disagrees with whichever one last regenerated the file.
_DEFAULT_GROUP_TITLE_SPELLINGS = frozenset({"optional arguments", "options"})
_CANONICAL_DEFAULT_GROUP_TITLE = "options"


def _normalize_group_title(title) -> str:
    """!
    @brief Resolve one argparse group title to a Python-version-independent form.
    @param[in] title Raw `group.title` from argparse, or None for the positional group.
    @return Canonical title.
    """
    if title in _DEFAULT_GROUP_TITLE_SPELLINGS:
        return _CANONICAL_DEFAULT_GROUP_TITLE
    return title or "arguments"


def subcommands(parser) -> dict:
    """!
    @brief The subparser map of a parser, if it has one.
    @param[in] parser Parser to inspect.
    @return Mapping of command name to subparser.
    """
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    return {}


def describe_parser(name: str, parser) -> dict:
    """!
    @brief Recursively capture a command and its subcommands.
    @param[in] name Command name.
    @param[in] parser Parser for that command.
    @return Nested command description.
    """
    groups: dict = {}
    for group in parser._action_groups:
        actions = [
            describe_action(action)
            for action in group._group_actions
            if not isinstance(action, argparse._SubParsersAction)
            and action.dest != argparse.SUPPRESS
        ]
        if actions:
            groups.setdefault(_normalize_group_title(group.title), []).extend(actions)
    children = {
        child: describe_parser(child, child_parser)
        for child, child_parser in sorted(subcommands(parser).items())
    }
    return {
        "name": name,
        "help": (parser.description or "").strip().splitlines()[0] if parser.description else "",
        "groups": groups,
        "subcommands": children,
    }


def snapshot() -> dict:
    """!
    @brief The full normalized parser snapshot.
    @return Snapshot mapping.
    """
    parser = build_parser()
    return {
        "program": parser.prog,
        "commands": {
            name: describe_parser(name, sub)
            for name, sub in sorted(subcommands(parser).items())
        },
    }


def render_command(command: dict, depth: int = 0, prefix: str = "") -> list:
    """!
    @brief Render one command and its subcommands as HTML.

    @details A subcommand heading carries its full invocation path. Rendering `status`
             on its own would sit indistinguishably beside the top-level
             `status-source`, and a reader scanning headings could not tell which
             command a flag table belongs to.
    @param[in] command Command description.
    @param[in] depth Nesting depth.
    @param[in] prefix Invocation path of the parent command, empty at the top level.
    @return HTML lines.
    """
    heading = "h3" if depth == 0 else "h4"
    invocation = f"{prefix}{command['name']}"
    slug = invocation.replace("-", "_").replace(" ", "_")
    out = [f'<{heading} id="cli_{slug}"><code>picurv {escape(invocation)}</code></{heading}>']
    if command["help"]:
        out.append(f"<p>{escape(command['help'])}</p>")
    for title, actions in command["groups"].items():
        if not actions:
            continue
        out.append(f"<p><b>{escape(title)}</b></p>")
        out.append('<table class="markdownTable"><tr>'
                   '<th class="markdownTableHeadNone">Flag</th>'
                   '<th class="markdownTableHeadNone">Required</th>'
                   '<th class="markdownTableHeadNone">Choices</th>'
                   '<th class="markdownTableHeadNone">Default</th>'
                   '<th class="markdownTableHeadNone">Description</th></tr>')
        for action in actions:
            flags = ", ".join(f"<code>{escape(flag)}</code>" for flag in action["flags"]) or \
                    f"<code>{escape(action['dest'])}</code> (positional)"
            choices = ", ".join(f"<code>{escape(c)}</code>" for c in action["choices"]) or "&#8212;"
            default = f"<code>{escape(action['default'])}</code>" if action["default"] else "&#8212;"
            out.append(
                "<tr>"
                f'<td class="markdownTableBodyNone">{flags}</td>'
                f'<td class="markdownTableBodyNone">{"yes" if action["required"] else "no"}</td>'
                f'<td class="markdownTableBodyNone">{choices}</td>'
                f'<td class="markdownTableBodyNone">{default}</td>'
                f'<td class="markdownTableBodyNone">{escape(action["help"])}</td>'
                "</tr>"
            )
        out.append("</table>")
    for child in command["subcommands"].values():
        out += render_command(child, depth + 1, f"{invocation} ")
    return out


def render(data: dict) -> str:
    """!
    @brief Render the whole reference as an includable HTML fragment.
    @param[in] data Parser snapshot.
    @return HTML text.
    """
    out = [
        "<!-- GENERATED FILE - do not edit by hand.",
        "     Source of truth: the assembled argparse parser (picurv_cli.cli.build_main_parser),",
        "     including registrars delegated from other modules such as storage.py.",
        "     Regenerate with: make docs-cli-reference -->",
        f"<p>{len(data['commands'])} top-level commands.</p>",
    ]
    for command in data["commands"].values():
        out += render_command(command)
    return "\n".join(out) + "\n"


def main() -> int:
    """!
    @brief Write or verify the generated CLI reference.
    @return Process status code.
    """
    parser = argparse.ArgumentParser(description="Generate the CLI reference from the live parser.")
    parser.add_argument("--check", action="store_true", help="Fail if the generated output is stale.")
    args = parser.parse_args()

    data = snapshot()
    payload = json.dumps(data, indent=2, sort_keys=True) + "\n"
    html = render(data)

    if args.check:
        stale = [
            path.relative_to(REPO_ROOT)
            for path, content in ((JSON_PATH, payload), (HTML_PATH, html))
            if not path.is_file() or path.read_text(encoding="utf-8") != content
        ]
        if stale:
            print("Generated CLI reference is stale:", file=sys.stderr)
            for path in stale:
                print(f"  {path}", file=sys.stderr)
            print("\nThe parser changed. Regenerate with: make docs-cli-reference", file=sys.stderr)
            return 1
        print(f"CLI reference is current ({len(data['commands'])} commands).")
        return 0

    GENERATED_DIR.mkdir(parents=True, exist_ok=True)
    JSON_PATH.write_text(payload, encoding="utf-8")
    HTML_PATH.write_text(html, encoding="utf-8")
    total = sum(1 for _ in data["commands"])
    subs = sum(len(c["subcommands"]) for c in data["commands"].values())
    print(f"Wrote CLI reference: {total} commands, {subs} subcommands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
