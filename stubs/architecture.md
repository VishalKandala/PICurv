/** @page architecture Architecture
@tableofcontents

@section modules Modules
- **Grid & Metrics**: curvilinear metrics, index transforms.
- **Particles**: particle storage, migration, interpolation kernels.
- **Solvers**: Poisson/momentum RHS assembly, linear solves.
- **I/O & Config**: case files, logging, options parsing.

@section dataflow Data Flow
1. Initialize simulation context
2. Load grid/BCs → assemble operators
3. Advance particles ↔ scatter/gather with Eulerian fields
4. Output statistics & checkpoints

@section diagrams Diagrams
Enable `HAVE_DOT=YES` to render include graphs and call graphs.
*/
