- Experimental: `picurv run` checks a carried restart's statistics windows against the source checkpoint
  before launching or submitting, by running the simulator on one rank without building the grid; an
  incompatible window is refused at staging, with the changed property named, instead of failing in the queue.
