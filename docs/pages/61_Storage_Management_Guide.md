@page 61_Storage_Management_Guide Storage Management Guide

@anchor _Storage_Management_Guide

PICurv can protect, offload, find, verify, and restore standalone runs and sweep-study data through an `rclone` remote. The storage commands preserve PICurv's existing directory model: runs remain under `runs/`, studies remain under `studies/`, and numbered study members remain under `studies/<study_id>/cases/case_####/`.

The storage catalog supplies the memorable names. Keep generated directory identities unchanged and attach a human-readable `--label` and repeatable `--tag` values when archiving.

@tableofcontents

@section p61_model_sec 1. The Storage Model

PICurv creates two kinds of top-level artifacts:

- a standalone run at `runs/<case-title>_<timestamp>/`, and
- a sweep at `studies/<study-title>_<timestamp>/`, containing numbered run directories at `cases/case_0000`, `cases/case_0001`, and so on.

Storage operations use those existing paths as selectors. No second naming scheme is introduced:

```bash
picurv storage status --run-dir runs/channel_20260824-143000
picurv storage status --study-dir studies/grid_study_20260824-150000
picurv storage status --study-dir studies/grid_study_20260824-150000 \
  --case-id case_0007
```

Three kinds of thing can be archived, and they are separate artifacts on purpose:

- a **run** or **study member**, which owns checkpoints and derived output;
- a **study**, which owns its members' shared scheduling and aggregate results;
- a **workspace**, which owns editable configuration, the input catalog, and the
  published asset store — but not its runs and studies, which carry their own
  archives. Duplicating them inside a workspace archive would store the same bytes
  twice and blur which object owns what.

```bash
picurv storage protect --workspace .                    # config, catalog, assets
picurv storage protect --workspace . --include-inputs    # and user-supplied files
```

User-supplied files under `<workspace>/inputs/` are excluded unless asked for, and no
offload policy prunes them or the editable configuration: protecting a workspace is a
backup, not custody of the user's own data.

There are two archival actions:

- `protect` uploads and verifies an immutable archive but keeps every local file;
- `offload` performs the same upload and verification, then removes the verified heavy payload while keeping a small, useful local skeleton.

What remains after `offload` is selected by a semantic policy. `metadata-only`
retains configuration, scheduler metadata, logs, and the storage marker;
`restart-ready` also retains run inputs and the newest committed checkpoint;
`analysis-ready` retains analysis and visualization products. An explicit
`--keep-latest-checkpoint` can add the newest checkpoint to another policy. Whole-study
offload applies the same component classification inside each numbered member.

The local state shown by `status` and `summarize` is:

- `LOCAL`: no PICurv archive marker exists;
- `PROTECTED`: a verified archive exists and all local files remain;
- `COLD`: heavy local payload was pruned;
- `PARTIAL`: selected checkpoints were restored, but the complete archive is not local;
- `BUSY`: a known local stage or Slurm job is active.

@section p61_setup_sec 2. One-Time Setup

Install `rclone`, configure credentials with `rclone config`, and test the remote using `rclone` itself. PICurv never stores cloud passwords, tokens, or access keys. Its workspace file contains only a remote path and storage policy.

From the workspace root:

```bash
picurv storage setup \
  --remote labstore:picurv-data \
  --profile archive \
  --compression auto \
  --chunk-size-gib 8 \
  --workers 16 \
  --offload-policy metadata-only
```

This verifies remote access, creates the remote `objects/` directory if necessary, and writes `.picurv-storage.yml`:

```yaml
default_profile: archive
profiles:
  archive:
    remote: labstore:picurv-data
    compression: auto
    chunk_size_gib: 8.0
    workers: 16
    offload_policy: metadata-only
    keep_latest_checkpoint: false
```

PICurv finds `.picurv-storage.yml` by searching upward from the current directory. Use `--storage-config /path/to/storage.yml` to select another file. Multiple profiles can share one file; select one with `--profile`. A fast local scratch filesystem can be selected during setup with `--staging-directory`.

Use a dry run to see what setup would write without touching the file or remote:

```bash
picurv storage setup --remote labstore:picurv-data --dry-run
```

@section p61_daily_sec 3. Day-to-Day Playbooks

@subsection p61_finished_run 3.1 A Finished Run Is Large and Can Leave the Cluster

First inspect the local state and the exact archive plan:

```bash
picurv storage status --run-dir runs/channel_20260824-143000
picurv storage plan --run-dir runs/channel_20260824-143000
```

Then offload it with a name and catalog fields that will still make sense months later:

```bash
picurv storage offload \
  --run-dir runs/channel_20260824-143000 \
  --policy restart-ready \
  --label "pilot 64, stable baseline" \
  --tag campaign=pilot-64 \
  --tag mesh=64 \
  --tag status=accepted
```

PICurv prints the archive ID. Save it in job notes if convenient, but the remote catalog also remains searchable. The original directory name does not need to be changed.

@subsection p61_keep_run 3.2 A Run Is Still Needed Locally

Protect it without pruning:

```bash
picurv storage protect \
  --run-dir runs/channel_20260824-143000 \
  --label "pilot 64, still processing"
```

This gives the run a verified remote copy while `--continue`, post-processing, and ordinary file access continue to use the full local tree.

@subsection p61_mixed_study 3.3 Only Some Cases in a Study Are Finished

Inspect all numbered members:

```bash
picurv storage status --study-dir studies/grid_study_20260824-150000
```

Offload selected members. Repeat `--case-id` to process several explicit cases:

```bash
picurv storage offload \
  --study-dir studies/grid_study_20260824-150000 \
  --case-id case_0003 \
  --case-id case_0007 \
  --label "completed grid members" \
  --tag campaign=grid-study
```

Each member receives its own archive ID. Other members remain untouched and can continue running or processing. PICurv refuses study continuation or reaggregation while required members are cold, so intentionally archived data cannot be mistaken for missing or incomplete output.

When the entire study is finished, archive it as one artifact by omitting `--case-id`:

```bash
picurv storage offload \
  --study-dir studies/grid_study_20260824-150000 \
  --label "complete grid convergence campaign"
```

@subsection p61_small_large 3.4 Small and Very Large Artifacts

Choose compression per operation when the profile default is not appropriate:

```bash
# Small/already-compressed data: package without compression.
picurv storage protect --run-dir runs/small_run --compression none

# Large data where transfer size matters more than compression time.
picurv storage offload --run-dir runs/large_run --compression maximum

# Lower CPU cost.
picurv storage offload --run-dir runs/medium_run --compression fast
```

The available compression policies are:

- `none`: uncompressed tar chunks;
- `fast`: gzip level 1;
- `balanced`: gzip level 6;
- `maximum`: xz level 9;
- `auto`: no compression below 256 MiB, balanced compression from 256 MiB through 20 GiB, and maximum compression at 20 GiB and above.

`--chunk-size-gib` is a transfer/staging target. PICurv places a file wholly in one chunk, so a single file larger than the target remains one larger chunk. Checkpoints are component-addressable, which enables selective restore.

Compression and restore are CPU-aware. `fast` and `balanced` use `pigz -p <workers>`
when installed; `maximum` uses `xz -T <workers>`. When native parallel compressors are
unavailable, PICurv uses Python codecs and parallelizes independent chunks. Restore
downloads, verifies, and extracts independent chunks concurrently. Set the reusable
profile default with `storage setup --workers N`, or override an operation with
`--workers N`. The plan reports the selected worker count, a broad compressed-size
range, retained local bytes, and bytes eligible for pruning.

@section p61_restore_sec 4. Finding and Restoring Old Data

@subsection p61_restore_skeleton 4.1 The Cold Local Skeleton Still Exists

The local marker knows the archive ID, so a full restore can use the familiar directory:

```bash
picurv storage restore --run-dir runs/channel_20260824-143000
```

For one individually archived study member:

```bash
picurv storage restore \
  --study-dir studies/grid_study_20260824-150000 \
  --case-id case_0007
```

A whole-study archive is restored with `--study-dir` and no case selector.

@subsection p61_restore_deleted 4.2 Every Local Directory Was Deleted

The catalog lives on the remote and does not depend on local marker files:

```bash
picurv storage list --search "pilot 64"
picurv storage show --archive-id 0123456789abcdef0123456789abcdef
picurv storage verify --archive-id 0123456789abcdef0123456789abcdef
```

Restore to the original absolute path recorded in the manifest:

```bash
picurv storage restore \
  --archive-id 0123456789abcdef0123456789abcdef
```

Or choose a new location:

```bash
picurv storage restore \
  --archive-id 0123456789abcdef0123456789abcdef \
  --to /cluster/new-workspace/runs/pilot-64-restored
```

When an individually archived study member is restored after the entire study was deleted, its archive also recreates the small parent study context, including `study.yml`, `study_manifest.json`, and `<run.scheduler>/case_index.tsv` when those files were present at archive time.

An alternate-location restore conservatively replaces the old run/study prefix in known generated text files (`.control`, `.run`, `.sbatch`, JSON, TSV, and YAML). Inspect regenerated scheduler scripts and any user-authored absolute paths before submitting on a different cluster. Remote archive bytes remain unchanged.

@subsection p61_restore_checkpoint 4.3 Only Particular Checkpoints Are Needed

Use repeatable checkpoint selectors to avoid downloading unrelated bulk output:

```bash
picurv storage restore \
  --archive-id 0123456789abcdef0123456789abcdef \
  --checkpoint 500 \
  --checkpoint 600
```

The result is `PARTIAL`. Restart and continuation proceed when their required checkpoint is present. Post-processing proceeds when all checkpoint steps in its requested start/end/interval window are present; otherwise PICurv prints the exact selective or full restore command to run.

`--force` permits a restore to merge into an existing directory that is not already the matching cold artifact. Because matching archived files can be replaced, use it only after inspecting the destination. Prefer a new `--to` path for recovery experiments.

@subsection p61_restore_component 4.4 Only One Semantic Component Is Needed

Use repeatable `--component` selectors when post-processing or inspection needs less
than the complete archive:

```bash
picurv storage restore --archive-id <id> --component analysis --component visualization
picurv storage restore --archive-id <id> --component inputs
```

Valid components are `inputs`, `raw-output`, `analysis`, `visualization`, `logs`, and
`assets`; checkpoint steps use `--checkpoint`. Metadata is always restored with a
partial selection so identity and config evidence remain available.

Archived run manifests also catalog the exact reusable assets in
`inputs/assets.lock.yml`. If the workspace's shared `assets/objects/` directory was
deleted, `picurv run --fetch-missing ...` searches remote manifests by provider-spec
hash, retrieves only the archived inputs component, verifies checksums, reconstructs
the immutable objects, and then binds the current case's asset set. If no matching
remote object exists, normal run behavior builds it; combine `--fetch-missing` with
`--require-precomputed` to refuse rebuilding.

@section p61_safety_sec 5. Verification and Failure Safety

PICurv uses this order for both `protect` and `offload`:

1. inventory the selected artifact without following symlinks;
2. reject known unsafe activity;
3. package bounded component chunks in a temporary staging directory, deleting each
   local chunk after its verified upload so staging does not grow to archive size;
4. upload each chunk with `rclone` and compare its remote SHA-256 with the local SHA-256;
5. upload a versioned manifest;
6. publish a `COMPLETE` marker bound to that manifest;
7. write the local storage marker;
8. for `offload` only, prune the verified heavy payload.

If packaging, transfer, or verification fails, local payload is not pruned. Incomplete remote objects have no valid completion marker and are ignored by `storage list`.

Archival/offload is refused when PICurv sees:

- an incomplete checkpoint directory;
- an active solver or post-processing lock;
- an active recorded Slurm job;
- recorded Slurm job IDs whose state cannot be checked because `squeue` is unavailable or failed; or
- another storage operation on the same artifact.

Use `storage plan` immediately before a real operation. Do not archive while an external process that PICurv cannot observe is modifying files.

Symlinks are stored as links and are never followed to pull unrelated data into an
archive. Explicit external input references are reported by `plan` and recorded in the
manifest, but are not copied automatically. Protect or restore those dependencies
separately.

Verify an archive again at any time:

```bash
picurv storage verify --archive-id 0123456789abcdef0123456789abcdef
```

@section p61_cli_integration_sec 6. Interaction With Existing Commands

The storage layer is additive and remains in the Python conductor. It does not change checkpoint bytes, numerical methods, or the C solver/postprocessor interfaces.

- `restart --from` behavior is represented by `picurv run --solve --restart-from ...`. A cold source is refused unless its requested start checkpoint has been restored.
- `run --solve --continue --run-dir ...` follows the same checkpoint rule.
- `run --post-process` accepts a partial restore only when the complete requested checkpoint window is local.
- `submit` refuses a cold run and refuses a study with cold members. Existing staged scripts and scheduler metadata are retained during offload.
- `sweep --continue` refuses cold members instead of interpreting pruned checkpoints as incomplete cases.
- `sweep --reaggregate` refuses cold members instead of replacing unavailable metrics with blank values.
- `summarize` remains usable from retained logs and displays the storage state, archive ID, and label. A requested summary that depends on pruned files can still report that those particular data are unavailable.
- `validate` checks YAML contracts and the optional workspace release requirement; it
  does not need bulk checkpoint payload.
- `cancel` continues to use retained scheduler metadata. Active Slurm jobs also prevent archival/offload in the first place.

Run creation snapshots initial YAML under `<run.config>`; continuations append revisions
under `<run.config.history>`; post recipes live under `<run.post_recipes>/<recipe-id>/`.
Storage classifies these as metadata, run-local materialized assets as inputs, and
canonical output subtrees by semantic component. Study creation uses the same fixed
member layout and stores aggregate analysis below `<run.analysis>`.

@subsection p61_prune_assets 4.5 Reclaiming a Shared Asset

Offloading a run does not remove the workspace assets it used, because another run may
still need them. Once nothing local does, and a workspace archive has been verified,
the local copies can go:

```bash
picurv storage prune --workspace . --assets --unused-locally --dry-run
```

Both selectors are required: `prune` removes nothing else, and it says what it decided
before it decides it.

```text
8df21abc0e4f1a92  grids
  referenced by remote runs       2
  referenced by active local runs 0
  remote protection               verified
  local removal                   safe
```

A run that is itself cold still *references* its assets — that is what keeps the remote
copy alive — but it no longer keeps the local one. An asset with any active local
reference, or without a verified remote copy, is reported as blocked and left alone.

@subsection p61_recover_workspace 4.6 Recovering From Complete Local Loss

Nothing needs to be remembered but the workspace name:

```bash
picurv storage setup --remote <remote>
picurv storage list --workspace-label pilot64
picurv storage restore --workspace-id pilot64 --to pilot64
```

That returns the editable configuration, the input catalog, and the asset store. Runs
and studies are restored individually from the same catalog, by id or by label.

@section p61_code_sec 7. Code and Configuration Reproducibility

Every remote manifest records:

- the generated run, study, and case identities;
- the human label and tags;
- original paths and creation time;
- source and stored sizes;
- compression and chunk checksums;
- checkpoint steps and checkpoint format version;
- SHA-256 values for canonical copied YAML files;
- the shared release version, Git commit, and dirty-worktree build identity;
- external paths and restart dependencies; and
- declared restore/continue/reprocess capabilities.

The archive preserves generated configs and data, but it deliberately does not bundle
the source tree, compiler/MPI libraries, or binaries. Use `picurv version` to inspect
the active build and `picurv versions activate <release>` to select a recorded tagged
release before processing old data. Exact binary reproduction still requires the
original build/container environment; the release and commit identify source, not
external libraries.

The storage schema is versioned. Unsupported future or older schema versions are rejected explicitly rather than guessed.

@section p61_remote_sec 8. Remote Layout and Catalog

Each archive is an immutable object below the configured remote:

```text
<remote>/objects/<archive-id>/
  chunks/
    00000_checkpoint-500.tar.xz
    00001_metadata.tar.xz
    ...
  manifest.json
  COMPLETE
```

The globally unique archive ID is the durable machine identity. `--label` and `--tag KEY=VALUE` are the human control surface. `storage list --search` searches archive IDs, labels, run/study/case identities, and tags case-insensitively. For scripts and audits, use JSON output:

```bash
picurv storage list --format json
picurv storage status --run-dir runs/channel_20260824-143000 --format json
picurv storage show --archive-id 0123456789abcdef0123456789abcdef
```

@section p61_reference_sec 9. Command Reference

```text
picurv storage setup    configure and verify a non-secret rclone profile
picurv storage status   inspect local lifecycle/storage state
picurv storage plan     preview chunks, compression, dependencies, and safety
picurv storage protect  upload and verify; retain all local files
picurv storage offload  upload and verify; prune verified heavy payload
picurv storage list     search completed remote archives
picurv storage show     print one complete remote manifest
picurv storage verify   checksum every chunk in one remote archive
picurv storage restore  restore a full archive or selected checkpoints
```

Use `picurv storage --help` and `picurv storage <action> --help` for the complete option set.

@section p61_cap_comp_sec 9.1 Compression Policy Entries

@htmlinclude generated/capability_inventory_storage_compression.html

@subsection p61_cap_comp_auto_sub auto

@anchor p61_cap_comp_auto

**Identity.** `picurv storage offload --compression auto`, or `compression: auto` in the storage profile. The default.

**What it does.** Chooses a concrete policy from the payload size: below 256 MiB it uses `none`, at or above 20 GiB it uses `maximum`, and everything between gets `balanced`.

**When to choose it.** Unless you have measured a reason not to. The thresholds encode the usual trade-off - small payloads are not worth the CPU, very large ones are worth the slowest compressor - and a fixed choice has to be re-justified as run sizes change.

**Parameters it owns.** None. The thresholds are constants (`AUTO_NO_COMPRESSION_BYTES`, `AUTO_MAXIMUM_COMPRESSION_BYTES`), not configuration.

**Interactions.** Resolves before the archive extension is chosen, so `auto` never appears in a filename. `picurv storage plan` reports the resolved policy without transferring anything.

**Diagnostics.** `picurv storage plan` prints the policy it resolved to and the payload size that decided it.

**Evidence.** Unit verified - `tests/test_storage.py` exercises the resolver and the archive path for this policy.**Limitations.** The thresholds are judgement, not measurement, and they take no account of the data's compressibility: an incompressible payload just above 256 MiB pays for `balanced` and gains nothing.

@subsection p61_cap_comp_none_sub none

@anchor p61_cap_comp_none

**Identity.** `--compression none` -> a `.tar` archive.

**What it does.** Archives without compressing. Bytes on disk equal bytes transferred.

**When to choose it.** For payloads that are already compressed - most binary checkpoint formats - and for anything where restore latency matters more than footprint. Also the right choice when the archive is a staging step before something that compresses anyway.

**Parameters it owns.** None.

**Interactions.** Produces the largest artefact and the fastest offload and restore of the four.

**Diagnostics.** The `.tar` suffix in the catalog listing identifies it after the fact.

**Evidence.** Unit verified - `tests/test_storage.py` exercises the resolver and the archive path for this policy.**Limitations.** Storage cost is the full payload size, which for a long campaign is the dominant term.

@subsection p61_cap_comp_fast_sub fast

@anchor p61_cap_comp_fast

**Identity.** `--compression fast` -> a `.tar.gz` archive.

**What it does.** Compresses with `pigz` level 1 across the requested CPU workers when
available, otherwise gzip level 1 across independent archive chunks.

**When to choose it.** When offload time is on the critical path - clearing scratch before a deadline, or archiving between queued jobs. It gives most of `balanced`'s size reduction for noticeably less CPU.

**Parameters it owns.** None.

**Interactions.** Produces the same `.tar.gz` extension as `balanced`, so the suffix
does not distinguish them; the catalog records the policy, worker count, and compressor.

**Diagnostics.** Catalog metadata records the policy used, which is the only way to tell `fast` and `balanced` archives apart afterwards.

**Evidence.** Implemented only.

**Limitations.** Shares an extension with `balanced`, so an archive inspected outside PICurv cannot be identified by filename alone.

@subsection p61_cap_comp_balanced_sub balanced

@anchor p61_cap_comp_balanced

**Identity.** `--compression balanced` -> a `.tar.gz` archive. What `auto` selects in the middle size band.

**What it does.** Compresses with `pigz` level 6 across requested workers when
available, otherwise gzip level 6 across independent chunks.

**When to choose it.** The sensible fixed choice when you want to pin a policy rather than let `auto` decide - typically so that archive sizes stay comparable across a campaign as run sizes drift across the `auto` thresholds.

**Parameters it owns.** None.

**Interactions.** Same extension as `fast`; the catalog record distinguishes them.

**Diagnostics.** Reported by `picurv storage plan` and recorded in the catalog entry.

**Evidence.** Unit verified - `tests/test_storage.py` exercises the resolver and the archive path for this policy.**Limitations.** Nothing establishes that gzip's default is the right point on the curve for this project's data; it is a conventional choice, not a measured one.

@subsection p61_cap_comp_maximum_sub maximum

@anchor p61_cap_comp_maximum

**Identity.** `--compression maximum` -> a `.tar.xz` archive. What `auto` selects at or above 20 GiB.

**What it does.** Compresses with `xz -T <workers> -9e`, which reaches the smallest
archives of the four at substantially higher CPU and memory cost.

**When to choose it.** For long-term archival of large campaigns that will rarely be restored - where storage is the recurring cost and decompression time is paid once, if ever.

**Parameters it owns.** None.

**Interactions.** Restore also parallelizes independent chunks. A single very large xz
chunk may still dominate latency, so chunk size and worker count should be planned
together.

**Diagnostics.** The `.tar.xz` suffix identifies it, and `picurv storage plan` reports it in advance.

**Evidence.** Unit verified — `tests/test_storage.py` creates and extracts a native xz
chunk while asserting the requested worker count.

**Limitations.** xz memory use scales with its window and worker count. Request the
compression on a compute node when cluster login-node policy or memory is restrictive.

@section p61_cap_policy_sec 9.2 Offload Policy Entries

@htmlinclude generated/capability_inventory_storage_offload_policy.html

@subsection p61_cap_policy_metadata_only_sub metadata-only

@anchor p61_cap_policy_metadata_only

**Identity.** `picurv storage offload --policy metadata-only`, or
`offload_policy: metadata-only` in a profile. The default.

**What it does.** Retains the local manifest, configuration, scheduler state, logs,
and storage marker after verified upload; other semantic payload is pruned.

**When to choose it.** Use it when the run is finished locally and the goal is maximum
scratch recovery while leaving enough identity to search, summarize logs, cancel, or
restore. Choose `restart-ready` for likely near-term continuation and `analysis-ready`
for local derived-data work.

**Parameters it owns.** `--keep-latest-checkpoint` optionally retains the newest
committed checkpoint in addition to the policy baseline.

**Interactions.** Post, restart, and continuation guards read the storage marker and
request the precise missing checkpoint or component instead of treating pruned data as
an incomplete run.

**Diagnostics.** `storage plan` reports retained and pruned byte estimates and names
the policy; `storage status` reports `COLD` after offload.

**Evidence.** Unit verified — `tests/test_storage.py` exercises the policy mapping,
pruning, state marker, and remote round trip.

**Limitations.** Logs and metadata can still be large in pathological runs; the policy
does not truncate them.

@subsection p61_cap_policy_restart_ready_sub restart-ready

@anchor p61_cap_policy_restart_ready

**Identity.** `picurv storage offload --policy restart-ready`.

**What it does.** Retains metadata, logs, all run-local inputs, and the latest committed
checkpoint while pruning older checkpoints and derived/bulk output.

**When to choose it.** Use it for completed segments that are likely to continue soon.
It consumes more local space than `metadata-only` but avoids downloading the common
restart path. Use `analysis-ready` instead when plots and post outputs matter locally.

**Parameters it owns.** Latest-checkpoint retention is intrinsic and cannot be disabled
without selecting another policy.

**Interactions.** `--continue` and `--restart-from` can proceed if they request the
retained step. Other steps remain cold and trigger a selective restore instruction.

**Diagnostics.** The plan names the exact retained checkpoint; the storage marker lists
`inputs` and `checkpoint:<step>` in `retained_components`.

**Evidence.** Unit verified — `tests/test_storage.py` exercises the policy mapping and
latest-checkpoint retention behavior.

**Limitations.** It retains only one checkpoint; branching from older states requires
a selective remote restore.

@subsection p61_cap_policy_analysis_ready_sub analysis-ready

@anchor p61_cap_policy_analysis_ready

**Identity.** `picurv storage offload --policy analysis-ready`.

**What it does.** Retains metadata, logs, `<run.analysis>`, and `<run.visualization>`
while pruning raw output, inputs, and checkpoints unless the latest checkpoint option
is added.

**When to choose it.** Use it after solver work is done but figures, statistics,
spectra, or comparisons remain active. Choose `metadata-only` when even derived output
can be cold, or `restart-ready` when continuation is more likely than analysis.

**Parameters it owns.** `--keep-latest-checkpoint` optionally adds restart readiness.

**Interactions.** Existing derived results remain readable; rerunning post may still
require restoring source checkpoints or raw output.

**Diagnostics.** `storage plan` separates retained analysis/visualization bytes from
pruned bytes, and the marker records both retained components.

**Evidence.** Unit verified — `tests/test_storage.py` exercises all named policy
mappings and semantic component classification.

**Limitations.** It cannot infer whether a particular future analysis needs raw fields;
restore those inputs explicitly when required.
