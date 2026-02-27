# Pages Guide

This directory contains authored Doxygen pages (`@page`) that populate the website navigation tabs.

## Authoring Rules

- Start each standalone page with a unique `@page <ID> <Title>` line.
- Use `@tableofcontents` for long pages.
- Prefer `@subpage` links for page-to-page navigation and `@ref` for API/function links.
- Keep IDs stable once published to avoid broken links.

## Navigation Wiring

Primary navigation ordering is controlled by:

- `docs/DoxygenLayout.xml`
- `docs/mainpage.md`

When adding a new page, wire it in both places if you want it discoverable from nav and landing page.

## Link Hygiene

Before pushing docs changes, run:

```bash
python3 scripts/check_markdown_links.py
```

This validates Markdown file links in `README.md` and `docs/**/*.md` (excluding generated docs output).
