- Experimental: particles sent to another rank or removed as lost are now dropped from the
  swarm in one pass instead of one at a time, and each rank's spare swarm capacity is 10% of
  its own particle count instead of the global count. Previously every removed particle
  re-zeroed a buffer the size of the whole run's particle count, which made migration
  dominate the step at large particle counts and gave each rank memory of order the global
  count. Results are bit-identical, including Brownian runs, because the survivors keep the
  same order. Cluster scaling evidence is pending.
