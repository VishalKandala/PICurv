# Geometric Periodic Hotfix: Cluster Test Plan

## 1) Scope of this patch
- Goal: validate geometric periodic correctness with modern helpers.
- Included in this patch:
  - `Ucont` is now part of `ApplyPeriodicBCs` orchestration.
  - `ApplyUcontPeriodicBCs` + `EnforceUcontPeriodicity` are called from the periodic pipeline.
  - Legacy periodic seam refresh in `VolumeFlux` is kept, with debug counters/logging.
- Not included in this patch:
  - driven periodic (`constant_flux`, `initial_flux`) behavior changes.

## 2) Preflight on cluster
1. Pull branch:
   - `git checkout hotfix/geometric-periodic-modern-helpers`
   - `git pull`
2. Build with your cluster PETSc toolchain.
3. Use your normal control flow (`case.control`, `case.bcs`, `whitelist.cfg`), no YAML required.

## 3) Required runtime options
In `case.control` (or command line), ensure:
- `-control_file <path/to/case.control>` is provided at launch (mandatory).
- `-bcs_files <path/to/case.bcs>` is set (or equivalent per block list).
- `-whitelist_config_file <path/to/whitelist.cfg>` is set.

Launch template:
```bash
export LOG_LEVEL=TRACE
mpiexec -n <NPROC> ./picsolver \
  -control_file /path/to/case.control \
  -bcs_files /path/to/case.bcs \
  -whitelist_config_file /path/to/whitelist.cfg
```

Notes:
- `DeterminePeriodicity` enforces periodic pairs (`-Xi/+Xi`, `-Eta/+Eta`, `-Zeta/+Zeta`) per block.
- Keep `NPROC >= 2` for at least one periodic test to validate MPI ghost/corner behavior.

## 4) Minimal whitelist.cfg for this hotfix
One function name per line, `#` comments allowed.

```txt
# Periodic pipeline
DeterminePeriodicity
SetupBoundaryConditions
ApplyBoundaryConditions
RefreshBoundaryGhostCells
ApplyPeriodicBCs
ApplyUcontPeriodicBCs
EnforceUcontPeriodicity
UpdateCornerNodes
UpdatePeriodicCornerNodes
UpdateLocalGhosts

# Solver path / overlap diagnostics
FlowSolver
VolumeFlux
```

## 5) Case matrix to run

### Case A (primary): mixed periodic/non-periodic
Configuration intent:
- 1 periodic pair (example: `-Xi/+Xi` periodic geometric)
- 1 inlet face, 1 outlet face
- remaining 2 faces walls

Example `case.bcs` pattern:
```txt
-Xi    PERIODIC  geometric
+Xi    PERIODIC  geometric
-Eta   WALL      noslip
+Eta   WALL      noslip
-Zeta  INLET     constant_velocity  vx=... vy=... vz=...
+Zeta  OUTLET    conservation
```

Pass criteria:
- run is stable (no BC parser or periodic consistency errors),
- no visible seam artifacts in periodic direction,
- continuity/divergence behavior remains acceptable for your baseline.

### Case B (requested): two periodic directions + inlet/outlet pair
Configuration intent:
- 4 periodic faces (`-Xi/+Xi` and `-Eta/+Eta`)
- inlet/outlet on remaining pair (`-Zeta/+Zeta`)

Example `case.bcs` pattern:
```txt
-Xi    PERIODIC  geometric
+Xi    PERIODIC  geometric
-Eta   PERIODIC  geometric
+Eta   PERIODIC  geometric
-Zeta  INLET     constant_velocity  vx=... vy=... vz=...
+Zeta  OUTLET    conservation
```

Pass criteria:
- stable run in both periodic directions,
- no corner/cross-edge mismatch symptoms,
- no abnormal jumps at periodic seams.

### Case C (regression): no periodic faces
Pass criteria:
- behavior matches pre-hotfix expectation,
- periodic path is effectively no-op (`ApplyPeriodicBCs` skip trace visible at `TRACE` level).

### Case D (MPI sanity): periodic case with multiple ranks
Run either Case A or B with `NPROC >= 2` (and your usual decomposition).

Pass criteria:
- no rank-dependent periodic mismatch,
- no MPI ghost/corner failures.

## 6) What logs to check
Look for these signatures:
- Periodic detection/setup:
  - `Global periodicity determined: I-periodic=..., J-periodic=..., K-periodic=...`
- Orchestrator activity:
  - `Applying periodic boundary conditions for all fields.`
  - `No periodic boundaries defined; skipping ApplyPeriodicBCs.` (Case C)
- Legacy overlap telemetry:
  - `VolumeFlux: entering legacy periodic Ucont seam refresh block.`
  - `VolumeFlux: legacy periodic Ucont seam refresh updated <count> values, max |delta|=<val>.`

Interpretation for `VolumeFlux` seam telemetry:
- This patch keeps that block intentionally.
- If `<count>` and `max |delta|` are near zero, modern periodic path is already doing most seam work.
- Nontrivial values indicate legacy block still contributes; keep note for next patch decision.

## 7) Report-back template (paste in chat)
```txt
Branch/commit tested:
MPI ranks used:
Case IDs run: A / B / C / D

Case A:
- Result: PASS/FAIL
- Stability notes:
- Seam artifact notes:
- Key logs (copy 2-5 lines):

Case B:
- Result: PASS/FAIL
- Stability notes:
- Corner/cross-edge notes:
- Key logs (copy 2-5 lines):

Case C:
- Result: PASS/FAIL
- No-op periodic path confirmed? (Y/N)
- Key logs:

Case D:
- Result: PASS/FAIL
- MPI-specific notes:

VolumeFlux seam telemetry summary:
- Representative `<count>`:
- Representative `max |delta|`:

Anything suspicious / candidate follow-up:
```

## 8) Expected next step after this validation
- If A/B/C/D pass cleanly: close geometric periodic hotfix PR.
- If issues remain: iterate on geometric periodic only.
- After geometric periodic is stable: start driven periodic `constant_flux` on top of this baseline.
