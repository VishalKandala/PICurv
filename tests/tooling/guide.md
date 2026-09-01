# Test Tooling Guide

This directory owns repository quality gates and audit utilities. These tools
are used by Make targets, CI, and focused developer checks; they are not part of
the user-facing PICurv runtime.

- `audit_function_docs.py`: C/Python function documentation audit.
- `audit_ingress.py`: PETSc option-ingress audit.
- `audit_ingress_manifest.json`: expected PETSc option manifest.
- `python_coverage_gate.py`: Python line-coverage gate.
- `c_coverage_gate.py`: C gcov summary and threshold gate.
- `check_markdown_links.py`: local Markdown-link checker across every repository Markdown file.
- `audit_generic_expansion.py`: rejects generic documentation-expansion debris (contract: `generic_expansion_contract.json`).
- `audit_docs_site.py`: verifies project-owned documentation URLs and Doxygen navigation targets (contract: `docs_site_contract.json`).
- `generate_capability_inventory.py`: extracts the public capability inventory from the
  executable sources named in `capability_families.json` and writes the generated
  fragments under `docs/generated/`. Run via `make docs-inventory`. It owns an exact
  file set: files under `docs/generated/` that no family produces are deleted on
  regeneration and reported as orphans by `--check`.
- `capability_families.json`: the capability family registry. Per family it declares the
  authoritative public surface (the Python validation layer), the parity sources forming
  the rest of the chain, the family page and entry-anchor prefix, per-value metadata
  (status, canonical/alias), and whether coverage is enforced or advisory.
- `audit_capability_coverage.py`: verifies full-chain parity, value metadata, and Tier-2
  documentation coverage. Run via `make audit-capability`. An unrecognized parity-source
  kind is a violation, not a silent pass.
- `scaffold_documentation.py`: generates documentation skeletons that already satisfy the
  contracts - capability entries, deprecated-alias stubs, typed pages, and subsystem
  charters. Run via `make scaffold ARGS="..."`. Tests assert the output passes the audits.
- `audit_path_literals.py`: rejects unmanaged run-path literals in narrative prose
  (contract: `artifact_topology.json`); runnable command blocks and the changelog are exempt.
- `contract_registry.json` / `audit_contracts.py`: the invariant-contract registry and its
  generic runner. Validation fails closed on unknown vocabulary, missing fields, or bad
  checker paths; checkers declaring `requires_built_docs` are skipped with a note rather
  than run against a stale build.
- `artifact_topology.json` / `extract_artifact_topology.py`: logical run-artifact identities
  and the normalized snapshot extracted from the CLI's own dry-run planner.
- `audit_family_census.py` / `family_census.json`: classifies every public selector surface
  as covered by a capability family or deliberately out of scope. A new public surface that
  no family covers fails the census. Run via `make audit-family-census`.
- `generate_cli_reference.py`: introspects the fully assembled argparse parser
  (`picurv_cli.cli.build_main_parser`, including delegated registrars such as storage) and
  writes the generated CLI reference. `--check` fails when the artifact no longer matches
  the live parser. Never regex the source for this.
- `page_types.json` / `audit_page_types.py`: the document type of every published page -
  Tutorial, How-to, Reference, Explanation, Hub. Central rather than per-page so typing 70
  pages costs no visible chrome; a page's inline `@pagemeta` must agree with the registry.
  Run via `make audit-page-types`.
- `subsystem_records.json` / `audit_subsystem_lifecycle.py`: what each subsystem's claimed
  lifecycle status owes in documentation, and which transitions are legal. Obligations are
  cumulative up the ladder. A reasoned not-applicable is accepted; a bare `N/A` is not.
  Negative tests are in `tests/test_subsystem_lifecycle.py`. Run via `make audit-subsystems`.
- `freshness_manifest.json` / `audit_freshness.py`: tiered staleness suspicion. Hard
  surfaces fingerprint a normalized artifact and block; soft surfaces fingerprint raw
  sources and are advisory unless explicitly promoted with a stated reason; a malformed
  manifest or a vanished watched path is an integrity failure and always blocks. Re-attest
  with `make attest-freshness ARGS="<surface-id>"` **after** actually comparing the page
  with its sources. Attestation records that a review happened, not that it found nothing.
- `review_packet.py`: joins declared routing for a page, contract, capability family,
  subsystem, freshness surface, or the current changed set. It prints whether registry
  routing is complete, incomplete, or unavailable and keeps that claim separate from code
  behavior. Run via `make review-packet PAGE=44`, `CONTRACT=<id>`, `CAPABILITY=<id>`,
  `SUBSYSTEM=<id>`, `SURFACE=<id>`, or `CHANGED=working-tree`.
- `generate_xref_index.py`: distills Doxygen XML into the ignored
  `docs_build/xref.json` source-reference cache, stamped from current source bytes and the
  Doxyfile. `make docs-xref` builds it; review packets consume it only while current and
  never treat its edges as a semantic call graph.
- `audit_agent_setup.py`: verifies that `AGENTS.md` is the regular canonical instruction
  file, `CLAUDE.md` imports it, Claude's local settings have an exact repository ignore,
  and the three `.claude/skills/` trees are materialized byte-identical copies of
  `.agents/skills/`. Run `make sync-agent-skills` after an intentional canonical skill
  edit and `make audit-agent-setup` in verification.
- `capability_scope_records.json`: non-rendered records for documentation that is scoped
  but deliberately not yet written, with the conditions that force it to be published.
- `inject_theme_sync.py`: injects the theme-sync script into the generated HTML so the
  built tree is publishable without further mutation.
- `certify_documentation.py`: clean-commit documentation certification; use
  `make certify-docs` for the full PETSc/MPI-backed certificate or
  `make certify-docs-fast` when only the structural/configuration gates are available.
- `stamp_docs_revision.py`: embeds the generated site's documented Git revision
  and commit link on every HTML page.
- `generate_doxygen_fallback_indexes.py`: Doxygen fallback index generator.
