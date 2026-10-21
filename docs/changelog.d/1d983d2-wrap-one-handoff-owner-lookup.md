- Experimental: a particle crossing several periodic seams in one step (an edge or corner) is
  wrapped across all of them before it is handed off, so it migrates once instead of once per
  seam; on a fully periodic box split along every axis this removes the extra settlement passes.
  The owner rank of a cell is now found by a search along each axis, not a scan of every rank.
