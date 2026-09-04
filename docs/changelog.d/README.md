Changelog fragments

One file per push to `main`: `<short-git-sha>-<slug>.md`. Content is the
exact Markdown bullet(s) to insert under the `## Unreleased` section of
`../CHANGELOG.md` - sub-bullets are fine for one commit doing several
related things, the same way a single commit already reads as one
bullet there today. The pre-push hook refuses a push with no new
fragment here; there is no exemption for a commit judged "trivial" by
whoever is pushing it. If a push genuinely has nothing user-visible to
say - pure refactor, a typo, CI log noise - the fragment says exactly
that in one line, so the omission is a stated fact instead of a silent
judgment call repeated at every push.

At release time, every fragment here is folded into a new dated version
section in `../CHANGELOG.md`, this directory is emptied back to just this
file, `VERSION` is bumped, and the release is tagged `vX.Y.Z`. Nothing in
`../CHANGELOG.md` is a real release until it carries that tag - `##
Unreleased` accumulates indefinitely otherwise, exactly like this file's
directory does between releases.
