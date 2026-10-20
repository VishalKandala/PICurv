- Experimental: particles are reflected at walls instead of being deleted.
  - Each boundary face has a `particle_policy` (`reflect` or `remove`). Walls reflect by
    default; inlets, outlets and periodic faces remove. A BC entry may override the
    default, except on a periodic face.
  - A reflected particle keeps its Psi, weight and ID, so the particle count and the
    carried scalar are conserved at walls. This changes the result of every existing case
    whose particles reach a wall; set `particle_policy: remove` to keep the old behaviour.
  - `search_metrics.csv` adds `reflections` and `reflection_fallbacks`, and the startup
    banner prints each face's policy.
