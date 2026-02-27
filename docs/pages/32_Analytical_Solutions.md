@page 32_Analytical_Solutions Analytical Solution Modes

This page documents analytical-mode options and where they are applied.

@tableofcontents

@section available_sec 1. Available Analytical Types

Configured via `-analytical_type` when `-euler_field_source analytical` is active.

Current types:
- `TGV3D`
- `ZERO_FLOW`

@section behavior_sec 2. Grid/Geometry Behavior

- Custom analytical geometry check: @ref AnalyticalTypeRequiresCustomGeometry
- Custom geometry builder: @ref SetAnalyticalGridInfo
- Runtime dispatcher: @ref AnalyticalSolutionEngine

Current policy:
- `TGV3D` requires custom geometry/decomposition.
- `ZERO_FLOW` uses standard programmatic-grid fallback path.

@section particle_sec 3. Particle Analytical Initialization

- Particle analytical setter entry point: @ref SetAnalyticalSolutionForParticles

This allows particle fields to stay consistent with analytical flow setups.

@section extension_sec 4. Adding a New Analytical Type

1. Add branch in @ref AnalyticalSolutionEngine.
2. Decide if custom geometry is required via @ref AnalyticalTypeRequiresCustomGeometry.
3. If needed, extend @ref SetAnalyticalGridInfo.
4. Add schema mapping (`solver.yml.operation_mode.analytical_type`) and docs update.
