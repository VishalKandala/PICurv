
- `picurv storage` is now `supported`, with its compression levels, offload policies, and the
  `checkpoints`, `logs`, `analysis`, `visualization`, and `inputs` retention components. A
  campaign against a Google Drive remote on 3.5 GiB turbulent-channel runs restored every
  file bit for bit, left exactly each policy's documented components, restarted from an
  offloaded-then-restored checkpoint, and measured each compression level (`fast` 70% of the
  raw size, `maximum` 65% at 2.3 times the time and 3.8 GB of memory). `raw-output` stays
  experimental. Page 61's compression guidance now reflects these measurements.
- A workspace restored to another directory no longer rewrites paths inside its asset
  objects, which left them failing their own checksums; checkpoint bundles are likewise left
  untouched, and rewritten files no longer change hard-linked copies.
- `storage protect` or `offload` with an explicit `--compression` no longer reuses an existing
  archive made at another level; without `--compression`, protect-then-offload still reuses it.
- `storage verify --workspace` verifies a workspace's own archive.
- `storage prune` reports an asset object it could not delete and exits non-zero, instead of
  counting it as removed.
- Restoring into a run's own directory skips components it already holds.
