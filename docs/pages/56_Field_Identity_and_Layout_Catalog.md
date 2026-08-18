@page 56_Field_Identity_and_Layout_Catalog Field Identity and Layout Catalog

@anchor _Field_Identity_and_Layout_Catalog

This page documents the first field-infrastructure phase: separate typed identity
and metadata catalogs for persistent Eulerian and solver-particle fields. The
catalogs remove repeated field-name dispatch from ghost updates, periodic
synchronization, field diagnostics, initialization, and Eulerian/particle
coupling without changing vector creation, numerical stencil/repair algorithms,
checkpoint formats, runtime configuration, or postprocessing recipe syntax.

@tableofcontents

@section p56_scope_sec 1. Implemented Scope

The catalog provides:

- `FieldId`, a compile-time identity for each field formerly recognized by
  @ref UpdateLocalGhosts;
- `FieldDescriptor`, immutable metadata describing name, aliases, degree of
  freedom, DM family, grid layout, ghost-repair class, availability conditions,
  capabilities, and the existing `UserCtx` vector binding;
- `FieldView`, a non-owning runtime view of an already-created DM/global/local
  vector tuple;
- @ref FieldIdFromName for resolving text at a true name-ingress boundary;
- @ref FieldGetDescriptor and @ref FieldGetView for shared consumers; and
- typed ghost-update, periodic-synchronization, initial-condition, and field
  diagnostic interfaces used by compiled callers.

The Eulerian catalog covers the 38 persistent fields supported by the previous
ghost dispatcher. The separate particle catalog covers the ten persistent
DMSwarm fields used by the solver, including PETSc-managed particle ID and rank
fields. Adding a field to a supported operation now requires metadata and a typed
call site rather than another string branch.

@section p56_lifecycle_sec 2. Identity and Storage Lifecycle

Field IDs are not assigned after PETSc vector creation. They are enum constants
known when the program is compiled. Storage is still created by the existing
setup lifecycle:

```text
program compiled
    -> FieldId values and immutable descriptors already exist
CreateSimulationContext
    -> runtime configuration is parsed as before
InitializeAllGridDMs
    -> da/fda and their actual periodic boundary topology are created as before
CreateAndInitializeAllVectors
    -> existing explicit UserCtx Vec creation remains unchanged
FieldGetView(user, id)
    -> resolves that existing DM and Vec pair; creates nothing
UpdateLocalGhosts(user, id)
    -> performs the existing scatter and periodic normal-face repair
FinalizeSimulation
    -> existing explicit destruction remains unchanged
```

Particle identities are likewise compile-time constants. During
@ref RegisterParticleFields, PICurv-owned descriptors drive the existing
DMSwarm registration helper; PETSc-owned `DMSwarm_pid` and `DMSwarm_rank`
descriptors are recorded but deliberately not re-registered. Particle IDs do
not depend on particle allocation order or on the number of local particles.

This separation is deliberate. `FieldId` identifies a concept, while
`FieldView` reports whether that concept has storage in a particular `UserCtx`.
Optional fields remain valid catalog entries even when their runtime conditions
are disabled; asking for their view then returns a wrong-state error instead of
silently dereferencing a null PETSc object.

@section p56_descriptor_sec 3. Descriptor Contract

Each @ref FieldDescriptor records:

| Member | Meaning |
| --- | --- |
| `id` | Typed identity used by compiled consumers |
| `canonical_name` | Stable printable/configuration-facing name |
| `alias_1`, `alias_2` | Existing accepted names resolved only at ingress |
| `dof` | PETSc degrees of freedom per DMDA entry: 1, 2, or 3 |
| `dm_kind` | `da`, `fda`, `fda2`, or PETSc coordinate storage |
| `layout` | node, shifted cell, I/J/K face, or component-staggered |
| `sync_class` | normal-face repair policy after global-to-local scatter |
| `availability` | setup conditions under which storage may exist |
| `capabilities` | operations for which the field is registered |
| vector offsets | binding to the existing global/local `UserCtx` members |

The availability flags describe requirements; they do not replace runtime state
validation. @ref FieldGetView always checks the actual DM and vector handles.

`ParticleFieldDescriptor` is separate because DMSwarm fields have a different
contract: canonical PETSc name, component count, `PetscDataType`, registration
owner, default initialization value, supported operations, and an optional
Eulerian scatter target. In particular, particle `Psi` maps explicitly to
Eulerian `FIELD_ID_PSI`; no name coincidence is used to infer that bridge.

The registered aliases are retained because current diagnostic and boundary
surfaces use them: `Eddy Viscosity` for `Nu_t`, `Cs` for `CS`,
`Center-Coordinates` for `Cent`, and the `X/Y/Z-Face-Centers` names for
`Centx/Centy/Centz`. Aliases do not create additional field identities.

@section p56_layout_sec 4. Layout, Boundaries, and Ghost Values

The catalog's topology is based on the solver's shifted/staggered architecture,
not just PETSc's owned-versus-ghosted partition.

| Catalog layout | Physical interpretation | Logical boundary convention |
| --- | --- | --- |
| `FIELD_LAYOUT_NODE_CENTERED` | grid nodes | direct indices; the extra high-side slot is not a physical node |
| `FIELD_LAYOUT_CELL_CENTERED` | cell-centered value | physical values occupy shifted indices `1..N-1`; indices `0` and `N` are boundary/dummy locations |
| `FIELD_LAYOUT_I_FACE` | one I-face family | node-like in I, cell-like in J and K |
| `FIELD_LAYOUT_J_FACE` | one J-face family | node-like in J, cell-like in I and K |
| `FIELD_LAYOUT_K_FACE` | one K-face family | node-like in K, cell-like in I and J |
| `FIELD_LAYOUT_COMPONENT_STAGGERED` | `Ucont`-style packed vector | x/y/z components live on I/J/K faces respectively |

Therefore, a local vector may contain two different categories that must not be
confused:

1. PETSc/MPI halo entries created by DMDA decomposition; and
2. solver-layout boundary, dummy, or unused indices that exist even on one MPI
   rank.

@ref UpdateLocalGhosts still calls the established
`DMGlobalToLocalBegin/End` pair. It then uses `sync_class`, instead of a field
name, to select the existing normal-face repair for I-, J-, K-face, or packed
component-staggered storage. Standard shifted-cell fields retain the same PETSc
scatter behavior.

@section p56_periodic_sec 5. Periodicity

Periodicity does not create a second identity for a field. It changes the
runtime DM topology and which physical boundary repairs are active:

- `InitializeAllGridDMs` still constructs the actual periodic/nonperiodic DMDA;
- @ref FieldGetView resolves that actual DM from the supplied `UserCtx`;
- the existing boundary-face state still decides which normal-face repairs run;
- geometric translation and periodic duplicate-plane algorithms in
  `Boundaries.c` are unchanged.

The periodic cell, single-face-family, and component-staggered synchronizers now
accept `FieldId` arrays. Catalog capabilities reject a field routed to an
incompatible synchronizer before the unchanged numerical copy loop runs.
Coordinate-face translation is selected through
`FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT`, not a field-name comparison. No
string comparison remains in the periodic field dispatcher.

@section p56_inventory_sec 6. Catalog Inventory

The initial catalog groups fields as follows:

- node-centered: `Coordinates`;
- shifted cell-centered: `Ucat`, `P`, `Nu_t`, `CS`, `Diffusivity`,
  `DiffusivityGradient`, `Nvert`, `Aj`, `Cent`, `GridSpace`, `Phi`, `Psi`,
  `Nvert_o`, `ParticleCount`, `K_Omega`, and `K_Omega_o`;
- component-staggered: `Ucont`, `Ucont_o`, and `Ucont_rm1`;
- I-face family: `Csi`, `Centx`, `ICsi`, `IEta`, `IZet`, and `IAj`;
- J-face family: `Eta`, `Centy`, `JCsi`, `JEta`, `JZet`, and `JAj`;
- K-face family: `Zet`, `Centz`, `KCsi`, `KEta`, `KZet`, and `KAj`.

`K_Omega` entries preserve the compiled RANS call surface, but the catalog does
not claim that the current setup path allocates their storage. The runtime view
check makes that state explicit.

The particle inventory is: `position`, `velocity`, `DMSwarm_CellID`, `weight`,
`Diffusivity`, `DiffusivityGradient`, `Psi`, `DMSwarm_location_status`,
`DMSwarm_pid`, and `DMSwarm_rank`. It is intentionally not mixed into
`FieldId`: particle fields are per-particle DMSwarm arrays, may be integer-valued,
and have no persistent `UserCtx` global/local Vec pair. Dynamic postprocessor
fields such as recipe-created `ske` or `disp` remain outside the persistent
solver-particle catalog.

@section p56_consumer_sec 7. Consumer Rules

Compiled code should pass IDs directly:

```c
PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));
PetscCall(InterpolateEulerFieldToSwarm(user, FIELD_ID_UCAT,
                                      PARTICLE_FIELD_ID_VELOCITY));
```

Code that genuinely receives text should resolve it once at its ingress:

```c
FieldId field_id;
PetscCall(FieldIdFromName(requested_name, &field_id));
PetscCall(UpdateLocalGhosts(user, field_id));
```

Consumers must not cache a @ref FieldView beyond the lifetime or re-creation of
the corresponding `UserCtx` PETSc objects. A view is non-owning and must never be
destroyed independently.

PETSc's DMSwarm API itself remains name-based. Code at that API boundary obtains
the canonical text with @ref ParticleFieldName; this does not require assigning
particle fields an Eulerian `FieldId` or wrapping every PETSc access function.
Generic gather/restart/output helpers also remain name-based for dynamic
postprocessor fields, but query the registered `PetscDataType` from DMSwarm
instead of inferring integer width from a field-name comparison.

@section p56_nongoals_sec 8. Explicit Non-Goals

The catalog deliberately does not:

- alter `CreateAndInitializeAllVectors`, setup ordering, or teardown ordering;
- infer or create fields from their IDs;
- change physical, dummy, or halo index ranges;
- rewrite periodic transfer or boundary-condition numerical algorithms;
- change function log filtering, runtime monitoring cadence/output, or profiling;
- change solver/postprocessor input or output formats;
- add statistics accumulation, moments, profiles, bins, or derived turbulence
  calculations, which are built on the catalog rather than in it (see
  @ref 58_Field_Statistics);
- assign IDs to functions; or
- serialize raw enum values into checkpoint files.

Runtime configuration strings and postprocessing recipe keywords are true text
ingress and remain textual. Function identity/profiling is intentionally deferred
to a later specification; it is not conflated with field identity.

These boundaries prevent a common-infrastructure change from silently becoming
a numerical, monitoring, I/O, or statistics change.

@section p56_validation_sec 9. Regression Coverage

`make unit-setup` contains direct catalog tests for:

- completeness and canonical-name round trips for every `FieldId`;
- case-insensitive canonical and alias lookup;
- cell, single-face, and component-staggered metadata;
- binding to existing scalar, vector, and PETSc coordinate storage;
- rejection of unknown names and unavailable optional storage;
- typed scalar, vector, scalar-face, and component-staggered ghost scatters;
- particle-catalog completeness, aliases, PETSc data types, registration owner,
  and the particle-`Psi`/Eulerian-`Psi` bridge;
- typed periodic synchronization across cell, face-family, and packed staggered
  storage; and
- catalog-driven field diagnostics and Eulerian/particle interpolation/scatter.

Existing `unit-periodic`, MPI, postprocessing, and runtime smoke tests remain the
behavioral regressions for the unchanged numerical paths.

@section p56_related_sec 10. Related Pages

- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@ref p54_geometric_periodic "Geometric Periodic Boundaries"**
- **@subpage 13_Code_Architecture**
- **@subpage 51_C_Test_Suite_Developer_Guide**
- **@subpage 57_Future_Architecture_Specifications**
- **@subpage 58_Field_Statistics**
