- Experimental: IEM micromixing now relaxes each particle toward the mean of the particles
  currently in its cell, not toward the mean scattered at the end of the previous step, before
  the particles moved. Runs with IEM change slightly; runs with `iem_constant: 0` are unchanged.
