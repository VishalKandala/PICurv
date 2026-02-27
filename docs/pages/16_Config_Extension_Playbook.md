@page 16_Config_Extension_Playbook Configuration Extension Playbook

This playbook defines the standard workflow for adding new options/models without breaking YAML->C contract consistency.

@tableofcontents

@section workflow_sec 1. Standard Extension Workflow

1. Choose the schema home (`case.yml`, `solver.yml`, `monitor.yml`, or `post.yml`).
2. Add template documentation in `examples/master_template/*.yml`.
3. Add validation and normalization in `scripts/pic.flow`.
4. Emit the mapped control/post key into generated artifacts.
5. Add C-side ingestion default + parser wiring (`setup.c` or `io.c`).
6. Connect runtime consumer logic in domain module (`grid.c`, `poisson.c`, particle modules, etc.).
7. Update references:
   - `docs/pages/14_Config_Contract.md`
   - `docs/pages/15_Config_Ingestion_Map.md`
8. Update `scripts/audit_ingress_manifest.json` when PETSc ingress changes.
9. Verify with a config-generation run and static ingress audit.

@section design_sec 2. Design Rules

- Prefer structured schema first; use passthrough only for advanced/temporary gaps.
- Keep defaults in one place and preserve backward compatibility when keys are absent.
- Reject unsupported shapes early in Python validation (for example, list-vs-scalar contract mismatches).
- Keep C ingestion concentrated in setup/parsing layers (`setup.c`, `io.c`) unless data is truly file-local to a subsystem.

@section particle_sec 3. ParticlePhysics Extension Checklist

Use this checklist when adding a new particle model constant, source term, or initialization mode.

1. Schema:
   - Add structured keys under `models.physics.particles` (case) or `solver` (if solver-policy).
2. Generator:
   - Map to control flags in `scripts/pic.flow`.
   - Add validation for domain/range and mode compatibility.
3. C ingestion:
   - Add default in `CreateSimulationContext`.
   - Parse flag via `PetscOptionsGet*` in `setup.c`.
4. Runtime:
   - Wire into `ParticlePhysics.c`, `ParticleSwarm.c`, or `ParticleMotion.c` as needed.
   - Ensure post/statistics kernels are updated if new fields are emitted.
5. I/O:
   - If new field is restart/post relevant, add read/write hooks in `io.c` and post field registration where needed.
6. Docs and audit:
   - Update contract/map pages and audit manifest.

@section data_driven_sec 4. Data-Driven Closure Model Note

- Offline integration (recommended first):
  - Treat solver/post outputs as the canonical data contract and run ML training/inference in external Python workflows.
  - Keep C-side ingestion unchanged unless the model must feed back into timestep updates.
- Tightly coupled inference (next step):
  - Add a runtime-selectable closure model interface in `ParticlePhysics.c` (model type + config + evaluation hook).
  - Expose model selection and coefficients in YAML (`solver.yml` or `case.yml`), map through `pic.flow`, parse in `setup.c`, and wire runtime consumers.
  - Keep fallback behavior to current deterministic model when no ML model is selected.

@section verification_sec 5. Verification Checklist

- Scalar/list and type validation works as expected.
- Generated control/post artifacts contain intended keys.
- C parser consumes keys without unknown-key warnings.
- Runtime behavior changes only when new key is explicitly set.
- Existing templates/run paths still work with omitted new keys.

See also **@subpage 17_Workflow_Extensibility** for higher-level workflow expansion patterns.
