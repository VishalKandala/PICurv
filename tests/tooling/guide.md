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
- `generate_doxygen_fallback_indexes.py`: Doxygen fallback index generator.
