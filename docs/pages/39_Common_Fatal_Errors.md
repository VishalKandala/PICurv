@page 39_Common_Fatal_Errors Common Fatal Errors and Fixes

Use this page to map error output to concrete fix steps.

@tableofcontents

@section structured_sec 1. Structured Validation Errors (ERROR <CODE>)

Validation and CLI-combination failures use:

`ERROR <CODE> | key=<yaml.path or -> | file=<path[:line] or -> | message=<...> | hint=<...>`

Code reference:

| Code | Meaning | Typical Fix |
|---|---|---|
| `CLI_USAGE_INVALID` | Invalid argument combination | Run `pic.flow <command> --help` and follow required flag groups. |
| `CFG_MISSING_SECTION` | Required YAML section missing | Add section from `examples/master_template/*.yml`. |
| `CFG_MISSING_KEY` | Required YAML key missing | Add required key at the reported path. |
| `CFG_INVALID_TYPE` | Wrong YAML value type | Change scalar/list/mapping type per contract docs. |
| `CFG_INVALID_VALUE` | Unsupported value/range | Use allowed enums/ranges from reference pages. |
| `CFG_FILE_NOT_FOUND` | Referenced file/path is missing | Fix path or create file. |
| `CFG_GRID_PARSE` | Grid payload/format error | Fix block count, dimensions, or coordinate rows. |
| `CFG_INCONSISTENT_COMBO` | Conflicting options/keys | Align related flags/keys (periodic pairs, scheduler, process counts). |

@section legacy_sec 2. High-Frequency Fatal Messages (Runtime Paths)

| Message fragment | Likely cause | Fix command(s) |
|---|---|---|
| `[FATAL] --post-process requires --run-dir when not used with --solve.` | Post-only run missing input run directory | `pic.flow run --post-process --run-dir runs/<run_id> --post post.yml` |
| `[FATAL] Could not automatically identify required config files in <run>/config` | Missing/renamed `case.yml`, `monitor.yml`, or `*.control` | Re-run solve stage or restore expected files under `runs/<run_id>/config/`. |
| `[FATAL] Source data directory for post-processing not found or empty` | Post source data path points to missing/empty solver output | Verify monitor output directory and existing solver output files. |
| `[FATAL] Unsupported scheduler '<x>'. Only Slurm is supported in v1.` | `cluster.yml` has unsupported scheduler type | Set `cluster.yml -> scheduler.type: slurm`. |
| `[FATAL] In cluster mode, --num-procs applies to the solver stage and must be 1 (auto) or exactly nodes*ntasks_per_node (...)` | Solver MPI count mismatch against cluster resources | Set `-n 1` or set solver `-n` to `nodes * ntasks_per_node`. |

@section restart_sec 3. Restart Mistakes To Check First

| Mistake | Likely result | Fix |
|---|---|---|
| Set `start_step: 501` after a run that ended at `500` | Off-by-one mismatch against the saved restart state | Set `start_step: 500`; the first new step will be `501`. |
| Leave `eulerian_field_source` as `solve` for a restart | Fresh solve path instead of loading saved fields | Set `solver.yml -> operation_mode.eulerian_field_source: load`. |
| Omit or mis-set `particles.restart_mode` | Unexpected particle reseed/load behavior or warning-driven behavior | Set `restart_mode: load` or `restart_mode: init` explicitly. |
| Choose a `start_step` that was never written | Runtime file-not-found or missing-data restart failure | Verify the requested restart step exists in the saved output/restart files. |

@section workflow_sec 4. Recommended Debug Workflow

1. Run config-only checks first:

```bash
./bin/pic.flow validate --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml
```

2. Preview launch plan:

```bash
./bin/pic.flow run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --dry-run
```

3. If cluster path is used, validate `cluster.yml` separately:

```bash
./bin/pic.flow validate --cluster cluster.yml
```

Related pages:
- **@subpage 02_Tutorial_Programmatic_Grid**
- **@subpage 14_Config_Contract**
- **@subpage 05_The_Conductor_Script**
- **@subpage 40_Testing_and_Quality_Guide**
