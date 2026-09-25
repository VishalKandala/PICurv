# Code and runtime explanations

Use this guidance when the question depends on implementation, control flow, distributed
state, or computational cost.

## Explain a live path, not an isolated symbol

Start from the user-visible ingress or owning caller and trace only far enough to establish
the behavior:

`configuration or CLI -> validation and normalization -> setup or dispatch -> producer -> synchronization or assembly -> consumer -> output`

Name the authoritative state and its owner. Distinguish process-lifetime state from
restart-persistent state, global from local PETSc vectors, owned physical entries from
ghost values, physical boundaries from decomposition halos, and canonical values from
temporary trial state.

When discussing distributed operations, state:

- what each rank may write;
- what representation supplies neighbor reads;
- which scatter, repair, assembly, or collective makes data current;
- what ordering the operation requires; and
- what result should be invariant under a different decomposition.

Do not describe a helper name as if it proves these semantics. Verify the caller, vectors,
DM or layout, and ordering when they affect the answer.

## Explain design and cost together

For a proposed change, identify the existing owner and reusable candidates inspected.
Explain why direct reuse or generalization is correct or unsuitable before proposing a new
helper, state field, or module. Mention allocation, copying, communication, collectives,
hot-loop indirection, and lifecycle effects when applicable.

For a bug, distinguish the demonstrated defect from nearby risks. Connect the symptom to
the causal execution path and the narrowest test or runtime check that can discriminate the
claim. A plausible code smell is not an observed failure.

Use file-and-line links for the few locations that establish the explanation. Avoid a file
tour or call-graph dump. Translate implementation terminology into the corresponding
numerical or user-visible consequence.
