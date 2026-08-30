@page 62_Capability_Status_Vocabulary Capability Status Vocabulary

@anchor _Capability_Status_Vocabulary

This page defines the status words PICurv documentation uses to describe a
capability. The words are load-bearing: they tell a reader whether a feature can
be relied on for published work, and they gate what documentation a feature owes
before it may claim the next level.

Use these words exactly. Do not invent synonyms, and do not leave a public
capability unlabelled.

@note **Contract status: enforced.** Every canonical value in every registered
capability family carries one of these statuses, and `make audit-capability` fails on
a missing status, a status outside this list, a non-selectable status on a reachable
value, or a documentation part the claimed status owes. The same vocabulary applies to
subsystems, where `make audit-subsystems` enforces the obligations each status carries
(**@subpage 64_Documentation_Extension_Framework** §4).

@tableofcontents

@section p62_status_sec 1. Status Levels

Status describes what the project claims about a capability today.

| Status | Meaning | May appear in user reference? |
|---|---|---|
| `planned` | Designed or proposed; no usable implementation exists. | No - design/specification pages only |
| `internal` | Implemented, but not exposed as a user-selectable option. | No - developer documentation only |
| `experimental` | User-selectable, but not established for production or published work. | Yes, with explicit limitations |
| `supported` | User-selectable, documented, evidenced, and suitable for production use. | Yes |
| `deprecated` | Still functional, superseded; scheduled for removal. | Yes, with migration guidance |
| `removed` | No longer present; retained only as history and migration guidance. | Changelog/history only |
| `known-defective` | Reachable by a user but known to produce incorrect results. | Yes - disclosure is mandatory |

@section p62_defective_sec 2. The `known-defective` Rule

`known-defective` is the one status with a hard placement requirement:

@warning A known defect must be disclosed at **every decision point where a user
can select the affected capability** - the configuration reference, the capability
entry, the method page, and any example that offers it. Disclosure buried in a
subsystem page that the selecting user never reads does not satisfy this rule.

A capability may be `known-defective` while the runtime still accepts it. The
documentation states the defect; whether the runtime should refuse the option is a
separate engineering decision and is tracked outside the documentation program.

@section p62_evidence_sec 3. Status Is Not Evidence

Status and evidence are independent axes and must not be collapsed into one
ladder. Production use does not imply analytical verification, and a validation
comparison does not imply restart equivalence.

Evidence is recorded as facets, each linking to the test, example, report, or
dataset that establishes it:

- unit verified
- integration/regression verified
- analytical or manufactured-solution verified
- benchmark characterized
- externally/reference validated
- production exercised

A capability with no evidence facet is `implemented only`. That is a legitimate
state to record; it is not a legitimate state to hide.

@note A recorded facet is a **declared evidence source**: automation checks that the
named test, example, or report exists and that the capability entry cites it, but it
cannot check that the source establishes the property. Recording a facet is a human
assertion that remains reviewable.

@section p62_stale_sec 4. Observations That Have Gone Stale

An observation measured against code that has since changed is neither current
behavior nor a defect. Record it as:

```text
Previously observed; requires re-characterization at current HEAD.
```

State when it was measured and what changed underneath it. Do not silently delete
such an observation, and do not let it keep circulating as current fact. When it
is rerun, replace it with the new measured result.

@section p62_related_sec 5. Related Documentation

- **@subpage 12_Capabilities_Summary**
- **@subpage 29_Maintenance_Backlog**
- **@subpage 57_Future_Architecture_Specifications**
