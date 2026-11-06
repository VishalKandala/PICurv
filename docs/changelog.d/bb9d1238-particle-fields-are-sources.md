
- Breaking: `models.physics.particles.fields` takes a source, as any target does: a name under `sources:`, `{expression, constants}`, or `{sum: [...]}`. A bare number or expression, `params` and `value`, and `xn`, `yn`, `zn` are refused with the replacement named.
- Particles (experimental): a field's expression also reads `xi`, `eta`, `zeta`, `L_xi`, `L_eta`, `L_zeta`, and the case's `constants`, `rho` and `nu`.
- Sources: `{background, regions}` is part of every source's expression, so a velocity may use it too.
