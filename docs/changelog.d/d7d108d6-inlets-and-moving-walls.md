
- Inlets (experimental): a source on a face reads `r`, the distance from the face's centre; `constraints.flow_rate` targets a volumetric rate; `inlet_flow.csv` records each inlet per step.
- Inlets (experimental): a data source may hold Cartesian velocities (`components: cartesian`) and several frames, interpolated linearly in time and optionally repeated (`repeat: true`). `field_slice` records a plane through `field_files` and `frame_interval`.
- Sources: the case's top-level `constants:`, and `rho` and `nu` from `properties.fluid`, are available to every source.
- Walls (experimental): `handler: moving` is a wall moving along itself at a velocity source. Refused with a wall function.
- Internal: the option prefix `-inlet_source_<id>_*` is now `-face_source_<id>_*`.
