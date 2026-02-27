@page 28_IEM_and_Statistical_Averaging IEM Mixing and Statistical Averaging

This page summarizes two related model/analysis paths: particle scalar mixing via IEM and flow/statistics averaging utilities.

@tableofcontents

@section iem_sec 1. IEM Particle Mixing

PICurv includes an IEM-style particle scalar update path where particle scalar state relaxes toward a local/ensemble mean with model-controlled rate.

Code touchpoints:
- @ref UpdateParticleField
- @ref UpdateFieldForAllParticles
- @ref UpdateAllParticleFields

Config touchpoints:
- particle model controls through case/solver ingestion (current model constants are C-side managed, extensible via schema additions).

@section avg_sec 2. Averaging Paths

- Solver-side statistical averaging toggles (`simCtx->averaging`) write/read dedicated statistical fields.
- Postprocessing statistics pipeline supports reduction kernels such as MSD.

Code touchpoints:
- Post statistics kernel: @ref ComputeParticleMSD
- Particle-grid averaging helper path: @ref ScatterAllParticleFieldsToEulerFields

@section note_sec 3. Terminology Note

The codebase currently uses "averaging" for several contexts (time/statistical fields, particle-to-grid normalization, post reductions). Configure these explicitly by workflow stage (`solver.yml`/`monitor.yml`/`post.yml`).

@section refs_sec 4. References

- IEM model overview (combustion/mixing context): https://en.wikipedia.org/wiki/Probability_density_function_(turbulence)#Micromixing_models
