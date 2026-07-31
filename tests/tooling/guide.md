# Test Tooling Guide

This directory owns repository quality gates and audit utilities. These tools
are used by Make targets, CI, and focused developer checks; they are not part of
the user-facing PICurv runtime.

- `audit_function_docs.py`: C/Python function documentation audit.
- `audit_ingress.py`: PETSc option-ingress audit.
- `audit_ingress_manifest.json`: expected PETSc option manifest.
- `python_coverage_gate.py`: Python line-coverage gate.
- `c_coverage_gate.py`: C gcov summary and threshold gate.
- `check_markdown_links.py`: local Markdown-link checker.
- `certify_documentation.py`: clean-commit documentation certification; use
  `make certify-docs` for the full PETSc/MPI-backed certificate or
  `make certify-docs-fast` when only the structural/configuration gates are available.
- `stamp_docs_revision.py`: embeds the generated site's documented Git revision
  and commit link on every HTML page.
- `generate_doxygen_fallback_indexes.py`: Doxygen fallback index generator.
