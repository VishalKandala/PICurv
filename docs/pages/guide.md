# Pages Guide

`docs/pages/` contains authored Doxygen pages (`@page`) that define the main information architecture of the PICurv documentation website. These pages are intentionally workflow-centered: installation, case setup, run orchestration, numerical methods, runtime mapping, testing, and extension playbooks.

A good page in this directory should answer three questions for the reader:
- what problem this page solves,
- where this topic appears in the run/development lifecycle,
- which next page the reader should use after finishing this one.

## Authoring Rules

- Start each standalone page with a unique `@page <ID> <Title>` line.
- Use `@tableofcontents` for long pages to preserve scanability.
- Prefer `@subpage` for page-to-page navigation and `@ref` for API/function linking.
- Keep IDs stable once published to avoid broken inbound links.
- Write explanatory paragraphs around lists so users understand why a step matters.

## Structural Pattern For New Pages

Use this sequence unless a topic demands a different structure:

1. Problem context and scope.
2. Runtime/config/data flow mapping.
3. Procedure or method explanation.
4. Validation and troubleshooting signals.
5. Related pages for onward navigation.

This pattern helps CFD users move from concept to execution without switching between too many documents.

## Navigation Wiring

Primary navigation order is controlled by:

- `docs/DoxygenLayout.xml`
- `docs/mainpage.md`

When adding a new page, wire it in both places so it is discoverable from site navigation and the landing map.

## Link Hygiene and Validation

Before pushing docs changes, run:

```bash
python3 scripts/check_markdown_links.py
```

This validates Markdown file links in `README.md` and `docs/**/*.md` (excluding generated docs output).
