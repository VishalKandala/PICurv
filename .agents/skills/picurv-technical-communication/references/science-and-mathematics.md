# Science and mathematics explanations

Use this guidance only when physical models, equations, discretization, or numerical
behavior materially affect the question.

## Establish the object being discussed

Define symbols near first use and state dimensions or nondimensionalization when relevant.
Identify whether a quantity is primitive, derived, modeled, imposed, or measured, and where
it lives: cell center, face, node, boundary, ghost region, block, or multigrid level.

Keep these layers distinct:

1. the continuous physical statement;
2. modeling assumptions and closures;
3. the discrete operator and boundary treatment;
4. the algebraic system and stopping criteria; and
5. the reported diagnostic or plotted quantity.

Do not use a continuous identity as proof of a property of the implemented discrete
operator. When writing an equation, connect each important term to its physical role and
then identify how PICurv represents or approximates it.

## Name the numerical consequence

Discuss only properties that bear on the question, choosing among:

- consistency and formal order;
- conservation and discrete balances;
- stability or timestep restrictions;
- conditioning and iterative convergence;
- grid, timestep, or sampling sensitivity;
- boundary and periodic treatment;
- curvilinear metrics and transformed units; and
- invariance under MPI decomposition.

Distinguish model error, discretization error, algebraic or iterative error, statistical
sampling uncertainty, and implementation defects. Similar-looking discrepancies can arise
from different members of this list.

For limiting cases or intuition, say which assumptions make the simplification valid and
which PICurv configurations satisfy them. Avoid implying that a pedagogical special case
describes the full solver.

## Match the reader

For a physics-first reader, introduce the numerical mechanism through its effect on
transport, pressure coupling, dissipation, conservation, or resolved scales. For a
numerics-first reader, begin from the operator, approximation, and error mechanism. For a
developer, connect both to array location, call order, and the verification oracle.

Prefer a small equation, dimensional check, or limiting-case example over a long generic
tutorial. State explicitly when a claim is general background rather than a behavior
verified in PICurv.
