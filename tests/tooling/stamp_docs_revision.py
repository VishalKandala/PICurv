#!/usr/bin/env python3
"""Stamp generated Doxygen HTML with the Git revision it documents."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def git_output(*args: str) -> str:
    """!
    @brief Return trimmed output from one Git command.
    @param[in] args Arguments passed after `git`.
    @return Decoded, trimmed Git stdout.
    """

    return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).strip()


def parse_args() -> argparse.Namespace:
    """!
    @brief Parse generated-documentation stamping arguments.
    @return Parsed command-line namespace.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html-dir", type=Path, required=True, help="generated Doxygen HTML directory")
    return parser.parse_args()


def main() -> int:
    """!
    @brief Write revision metadata and load it on every generated HTML page.
    @return Process status code.
    """

    args = parse_args()
    html_dir = args.html_dir.resolve()
    if not html_dir.is_dir():
        raise SystemExit(f"HTML directory does not exist: {html_dir}")
    sha = git_output("rev-parse", "HEAD")
    remote = git_output("remote", "get-url", "origin")
    if remote.startswith("git@github.com:"):
        remote = "https://github.com/" + remote[len("git@github.com:"):]
    if remote.endswith(".git"):
        remote = remote[:-4]
    revision = {
        "sha": sha,
        "short_sha": sha[:12],
        "commit_url": f"{remote}/commit/{sha}",
        "clean": not bool(git_output("status", "--porcelain")),
    }
    (html_dir / "picurv-docs-revision.js").write_text(
        "window.PICURV_DOCS_REVISION = " + json.dumps(revision, sort_keys=True) + ";\n",
        encoding="utf-8",
    )
    for page in html_dir.rglob("*.html"):
        relative = page.relative_to(html_dir)
        script_path = "../picurv-docs-revision.js" if relative.parts[0] == "search" else "picurv-docs-revision.js"
        content = page.read_text(encoding="utf-8")
        if "picurv-docs-revision.js" not in content:
            content = content.replace("</head>", f'<script src="{script_path}"></script>\n</head>')
            page.write_text(content, encoding="utf-8")
    print(f"Stamped generated documentation through commit {sha}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
