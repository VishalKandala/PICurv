# Writing for the research log

The research log is the long-form record of PICurv's methods and studies, kept in LaTeX on
Overleaf (git remote `overleaf`, branch `main`; local clone `../picurv-research-log` beside the
repository). Read its README for the build, `notation.tex` (use the macros), `conventions.tex`
and `nomenclature.tex` before writing. Pull before editing and push after.

## Where an entry goes

- **Part I** (chapters 00-11d): one chapter per subsystem, each in the order overview, model,
  discretisation, implementation, limitations. A code change that alters a method updates the
  owning chapter so it describes the method as it now is. Code identifiers appear only in the
  "Implementation in PICurv" boxes and the nomenclature appendix.
- **Part II** (chapters 12-17b): studies, verification and campaigns, dated where the date
  matters. A campaign step or investigation adds or extends a study: question, set-up, runs
  (with record, product, run and archive ids), observations, interpretation, limitations.
- Figures are TikZ where possible, or built by scripts in `figures/` from product data, never
  pasted images. References go in `references.bib` only with a DOI or publisher record.

## How it reads

- Professional and impersonal, present tense; no conversational or self-referential text, no
  reference to chats, no status labels or reviewer markers. Open questions go to the owner,
  not into the document. Planned work does not appear in Part I.
- Prove, cite or cut every claim.
- Readers include a new student: define each term where it first appears and open each
  section with what it is for. Avoid jargon that needs outside context.
- A word must not mean one thing in ordinary use and another in the method. Before each
  update, sweep the edited text for: unknown, exact, closed, solution, error, converged,
  residual - and use "unclosed", "not expressible in resolved quantities", and similar where
  the ordinary reading would mislead.
- Build locally (pdfLaTeX + BibTeX) before pushing; fix every undefined reference and
  overfull box the edit introduces.
