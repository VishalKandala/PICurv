PICurv skill routing (project SessionStart hook). Load the matching skill as the FIRST step
of a task - before reading code to change it, editing, or staging a run - every time, even if
it was loaded for an earlier task in this session:

- code that adds or changes behaviour, an option, output or contract -> picurv-capability-change
- a loop over cells, faces or particles, a stencil, ghosts, ownership -> picurv-spatial-kernel-change
- a run, test or number that is wrong for an unknown reason          -> picurv-solver-debugging
- setting up, piloting, extending or submitting runs or campaigns     -> picurv-case-design
- writing, choosing or running tests                                  -> picurv-test-design
- any substantive explanation, plan or write-up (and the research log) -> picurv-technical-communication
- research data, products, records, papers, Drive                     -> research-data

Work is done when its documentation trail is complete (AGENTS.md, "Documentation trail").
