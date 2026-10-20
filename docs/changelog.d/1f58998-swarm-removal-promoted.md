- The one-pass swarm removal (71abdd9) is no longer experimental. On TAMU Grace it left checkpoints
  and search metrics byte-identical on 64 and 48 ranks, and cut the time per step from 29.1 to 6.4 s
  (64 ranks, 4.2M particles) and from 10.9 to 1.3 s (48 ranks, 1M particles). Migration fell to about
  1 ms per step, and peak per-rank memory dropped by up to 59%. Measurement record
  `swarm-one-pass-removal-2026-10-20`.
