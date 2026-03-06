@page 39_Common_Fatal_Errors Common Fatal Errors and Fixes

@anchor _Common_Fatal_Errors

Use this page to map error output to concrete fix steps.

@tableofcontents

@section p39_structured_sec 1. Structured Validation Errors (ERROR &lt;CODE&gt;)

Validation and CLI-combination failures use:

`ERROR &lt;CODE&gt; | key=&lt;yaml.path or -&gt; | file=&lt;path[:line] or -&gt; | message=&lt;...&gt; | hint=&lt;...&gt;`

Code reference:

| Code | Meaning | Typical Fix |
|---|---|---|
| `CLI_USAGE_INVALID` | Invalid argument combination | Run `picurv &lt;command&gt; --help` and follow required flag groups. |
| `CFG_MISSING_SECTION` | Required YAML section missing | Add section from `examples/master_template/*.yml`. |
| `CFG_MISSING_KEY` | Required YAML key missing | Add required key at the reported path. |
| `CFG_INVALID_TYPE` | Wrong YAML value type | Change scalar/list/mapping type per contract docs. |
| `CFG_INVALID_VALUE` | Unsupported value/range | Use allowed enums/ranges from reference pages. |
| `CFG_FILE_NOT_FOUND` | Referenced file/path is missing | Fix path or create file. |
| `CFG_GRID_PARSE` | Grid payload/format error | Fix block count, dimensions, or coordinate rows. |
| `CFG_INCONSISTENT_COMBO` | Conflicting options/keys | Align related flags/keys (periodic pairs, scheduler, process counts). |

@section p39_legacy_sec 2. High-Frequency Fatal Messages (Runtime Paths)

| Message fragment | Likely cause | Fix command(s) |
|---|---|---|
| [FATAL] --post-process requires --run-dir when not used with --solve. | Post-only run missing input run directory | picurv run --post-process --run-dir runs/&lt;run_id&gt; --post post.yml |
| [FATAL] Could not automatically identify required config files in &lt;run&gt;/config | Missing/renamed case.yml, monitor.yml, or *.control | Re-run solve stage or restore expected files under runs/&lt;run_id&gt;/config/. |
| [FATAL] Source data directory for post-processing not found or empty | Post source data path points to missing/empty solver output | Verify monitor output directory and existing solver output files. |
| [FATAL] Unsupported scheduler '&lt;x&gt;'. Only Slurm is supported in v1. | cluster.yml has unsupported scheduler type | Set cluster.yml -> scheduler.type: slurm. |
| [FATAL] In cluster mode, --num-procs applies to the solver stage and must be 1 (auto) or exactly nodes*ntasks_per_node (...) | Solver MPI count mismatch against cluster resources | Set -n 1 or set solver -n to nodes * ntasks_per_node. |

@section p39_restart_sec 3. Restart Mistakes To Check First

| Mistake | Likely result | Fix |
|---|---|---|
| Set `start_step: 501` after a run that ended at `500` | Off-by-one mismatch against the saved restart state | Set `start_step: 500`; the first new step will be `501`. |
| Leave `eulerian_field_source` as `solve` for a restart | Fresh solve path instead of loading saved fields | Set `solver.yml -> operation_mode.eulerian_field_source: load`. |
| Omit or mis-set `particles.restart_mode` | Unexpected particle reseed/load behavior or warning-driven behavior | Set `restart_mode: load` or `restart_mode: init` explicitly. |
| Choose a `start_step` that was never written | Runtime file-not-found or missing-data restart failure | Verify the requested restart step exists in the saved output/restart files. |

@section p39_workflow_sec 4. Recommended Debug Workflow

1. Run config-only checks first:

```bash
./bin/picurv validate --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml
```

2. Preview launch plan:

```bash
./bin/picurv run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --dry-run
```

3. If cluster path is used, validate `cluster.yml` separately:

```bash
./bin/picurv validate --cluster cluster.yml
```

Related pages:
- **@subpage 02_Tutorial_Programmatic_Grid**
- **@subpage 14_Config_Contract**
- **@subpage 05_The_Conductor_Script**
- **@subpage 40_Testing_and_Quality_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Common Fatal Errors and Fixes** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/&lt;run_id&gt;/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
