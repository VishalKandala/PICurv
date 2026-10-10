- The three inlet handlers share one face-application path, and the two inlet handlers
  that no boundary file could select (`interp_from_file`, `pulsatile_flux`) are removed.
  `constant_velocity` is documented as it behaves: only the face-normal component of
  `vx`/`vy`/`vz` is applied.
