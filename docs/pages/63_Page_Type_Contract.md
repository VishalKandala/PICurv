@page 63_Page_Type_Contract Page Type Contract

@anchor _Page_Type_Contract

@pagemeta{Reference, Documentation authors, Enforced by make audit-page-types}

Every documentation page serves exactly one reader intent. Mixing intents is the
main reason a documentation set feels simultaneously bloated and unhelpful: a
tutorial that grows into a reference stops being followable, and a reference that
accumulates advice stops being scannable.

This page defines the four **content** types and one **structural** type, what each
owes the reader, and what each must refuse to do.

@tableofcontents

@section p63_types_sec 1. The Five Page Types

| Type | Reader intent | Shape |
|---|---|---|
| Tutorial | "Teach me by doing" | Ordered steps, one verified path, guaranteed to work end to end |
| How-to | "Help me accomplish X" | Task recipe, assumes competence, may branch on situation |
| Reference | "Tell me exactly what X is" | Exhaustive, structured, consultable in any order |
| Explanation | "Help me understand why" | Discursive, motivates design, discusses trade-offs |
| Hub | "Where do I go for X" | Structural only: owns the `@subpage` hierarchy and routes; carries no reference prose of its own |

The first four are content types and describe what a page teaches. **Hub** is a
fifth, *structural* type: it exists because the navigation tree is built from the
`@subpage` hierarchy, so some page must own that structure. A hub is judged by
whether its routes are complete and current, not by what it explains.

Every published page is assigned exactly one of these types in
`tests/tooling/page_types.json`. The registry is central rather than per-page so that
typing 71 pages costs no visible chrome on any of them; a page may still declare its
type inline with `@pagemeta`, and where both exist `make audit-page-types` requires
them to agree.

@section p63_owes_sec 2. What Each Type Owes

**Tutorial** owes a working outcome. Every command must run as written, in the
order given, from a stated starting state. It states its outcome up front and what
success looks like at the end. It links out for depth rather than explaining
inline. If a step can fail commonly, it names the failure and the fix.

**How-to** owes a completed task. It assumes the reader can already run PICurv and
wants a specific result. Unlike a tutorial it may branch, and it need not be
exhaustive about alternatives - it recommends one route.

**Reference** owes accuracy and completeness within its declared scope. It is
consultable, not readable: structure, tables, and predictable ordering matter more
than prose flow. Where a generated artifact exists, reference embeds or defers to
it rather than restating it by hand.

**Explanation** owes understanding. It may discuss history, alternatives
considered, and known limitations. It is the correct home for "why is it like
this" and for scientific formulation.

@section p63_refuses_sec 3. What Each Type Must Refuse

- A **tutorial** must not become exhaustive. If it lists every option, it is a
  reference wearing a tutorial's clothes.
- A **how-to** must not teach fundamentals. Link to the tutorial instead.
- A **reference** must not accumulate general advice. Operating guidance belongs
  in a how-to; rationale belongs in explanation.
- An **explanation** must not be the only place a fact is documented. If a reader
  needs it to run something, it belongs in reference or how-to too.

@warning The generic "CFD Reader Guidance and Practical Use" sections removed
during this overhaul were a textbook violation: identical advisory prose appended
to installation pages, API status pages, and method pages alike. Guidance that
would read the same on any page is not guidance. `make audit-docs-expansion`
now rejects it by exact signature.

@section p63_mixed_sec 4. Pages That Legitimately Mix

Pages 42 and 43 declare type `Hub` and are judged by the Hub rules above: they own
the `@subpage` structure the navigation tree renders, and they stay lean. A hub that
starts restating reference prose has become a duplicate of the pages it routes to.

Page 41 is deliberately **not** a hub. It declares `Tutorial`, because it now carries
a runnable first-run path rather than only routing. It also adopts subpages, which is
allowed: owning hierarchy is a hub's *only* job, not its exclusive privilege.

Where a page genuinely needs two types, prefer splitting. Where splitting would
scatter a single coherent contract across pages, keep it whole and say so.

@section p63_enforcement_sec 5. How the Contract Is Enforced

`make audit-page-types` fails when a published page has no assignment, when an
assignment names a type outside the five, when a page's inline `@pagemeta` disagrees
with the registry, or when the registry types a page that no longer exists. Coverage
is measured against the pages the build actually publishes, so a page added without a
type cannot certify. The audit runs in `make certify-docs` and in CI.

The registry types the type; it cannot judge the shape. Whether a page assigned
`Reference` actually reads as reference is a review question, and §3 is the standard
that review applies.

@section p63_related_sec 6. Related Documentation

- **@subpage 62_Capability_Status_Vocabulary** - status words and evidence facets
- **@subpage Documentation_Map** - the full page index
