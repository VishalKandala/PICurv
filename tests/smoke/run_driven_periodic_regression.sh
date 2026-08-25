#!/usr/bin/env bash
#
# Driven-periodic and multigrid coarse-solve regression.
#
# Covers two invariants that used to be violated silently:
#
#  1. Poisson multigrid coarse solve (Part A). With a Krylov method at
#     level_0 the multigrid preconditioner is a nonlinear operator, so the
#     outer FGMRES tracked residual ("Unprecond Norm") stops describing the
#     true residual ("True Norm") and the convergence test passes on a number
#     that is orders of magnitude smaller than b-Ax. The failure was
#     rank-dependent: 6 and 10 ranks were clean, 4 and 8 were not. This runs
#     at 4 and 10 ranks and asserts the two norms agree, plus a divergence
#     bound the tracked-residual failure used to break.
#
#  2. Driven-periodic initial_flux (Part B). The handler measures the
#     volumetric flux of the initial condition and holds it. This asserts the
#     latched target matches the analytic flux of the initial condition, and
#     that a restart restores the latched target rather than re-measuring it
#     from a field that has since drifted.
#
# Verified discrimination on this fixture (33^3 cells, 2 multigrid levels):
#
#   level_0 {preonly, redundant}  ->  worst tracked-vs-true difference 0.0
#   level_0 {fgmres,  bjacobi }   ->  worst tracked-vs-true difference 1.7e-04, FAILS
#
# The separation is much larger on a production grid (see the worked table in
# docs/pages/25_Pressure_Poisson_GMRES_Multigrid.md, where the true residual was
# five orders of magnitude above the tracked one). Note also that the norms are
# parsed from formatted log output carrying about six significant digits, so
# roughly 1e-6 is the finest difference this check can resolve at all; the 1e-4
# threshold below is the "agree to at least four significant figures" criterion
# and sits comfortably above that floor.
#
# Usage: run_driven_periodic_regression.sh <simulator> <mpi-launcher> [ranks...]

set -euo pipefail

simulator_exe_input="${1:?missing simulator path}"
mpi_launcher="${2:?missing MPI launcher}"
shift 2
ranks=("$@")
if [[ ${#ranks[@]} -eq 0 ]]; then
  ranks=(4 10)
fi

simulator_exe="$(cd "$(dirname "${simulator_exe_input}")" && pwd)/$(basename "${simulator_exe_input}")"
picurv_exe="$(dirname "${simulator_exe}")/picurv"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
# Dedicated fixtures rather than the shipped example: the regression must run
# at 10 ranks, which constrains the grid and multigrid depth far more than the
# physics case does. See config/grids/plane_channel_regression.cfg.
case_dir="${repo_root}/tests/smoke/driven_periodic"
tmp_root="$(mktemp -d)"

cleanup() {
  if [[ "${KEEP_SMOKE_TMP:-0}" == "1" ]]; then
    echo "Driven-periodic debug: preserving workspace at '${tmp_root}'." >&2
    return
  fi
  rm -rf "${tmp_root}"
}
trap cleanup EXIT

die() { echo "Driven-periodic regression failure: $*" >&2; exit 1; }

[[ -x "${picurv_exe}" ]] || die "picurv is not executable at '${picurv_exe}'."
[[ -d "${case_dir}" ]]  || die "missing case directory '${case_dir}'."

# --- Stage a short variant of the shipped laminar case. -----------------------
# Only the step counts and output cadence change; the grid, the solver profile
# (including the level_0 coarse solve under test) and the boundary conditions
# are exactly what ships.
STEPS=4
stage_case() {
  local handler="$1" dest="$2" start_step="$3"
  mkdir -p "${dest}"
  cp "${case_dir}/solver.yml" "${dest}/solver.yml"
  cp "${case_dir}/monitor.yml" "${dest}/monitor.yml"
  sed -e "s|^  total_steps: .*|  total_steps: ${STEPS}|" \
      -e "s|^  start_step: .*|  start_step: ${start_step}|" \
      -e "s|\.\./\.\./\.\./config/grids/|${repo_root}/config/grids/|" \
      "${case_dir}/case_${handler}.yml" > "${dest}/case.yml"
}

run_case() {
  local dest="$1" nprocs="$2"; shift 2
  ( cd "${tmp_root}" && "${picurv_exe}" run --solve -n "${nprocs}" \
      --case "${dest}/case.yml" --solver "${dest}/solver.yml" \
      --monitor "${dest}/monitor.yml" "$@" ) > "${dest}/run.log" 2>&1 \
    || { cat "${dest}/run.log" >&2; die "solver run failed (see log above)."; }
  local run_dir
  run_dir="$(sed -n 's|^  Run directory  : ||p' "${dest}/run.log" | tail -1)"
  [[ -n "${run_dir}" ]] || die "could not determine the run directory from ${dest}/run.log."
  echo "${tmp_root}/${run_dir}"
}

# --- 1. Poisson tracked-vs-true residual agreement, per rank count. -----------
for nprocs in "${ranks[@]}"; do
  echo "--- Poisson coarse-solve check on ${nprocs} ranks ---"
  dest="${tmp_root}/poisson_n${nprocs}"
  stage_case constant_flux "${dest}" 0
  run_dir="$(run_case "${dest}" "${nprocs}")"

  poisson_log="${run_dir}/logs/Poisson_Solver_Convergence_History_Block_0.log"
  [[ -f "${poisson_log}" ]] || die "missing Poisson convergence log '${poisson_log}'."

  python3 - "${poisson_log}" "${nprocs}" <<'PY' || exit 1
import re, sys
path, nprocs = sys.argv[1], sys.argv[2]
row = re.compile(r"Unprecond Norm:\s*(\S+)\s*\|\s*True Norm:\s*(\S+)")
checked = worst = 0, 0.0
checked, worst, worst_line = 0, 0.0, ""
for line in open(path):
    m = row.search(line)
    if not m:
        continue
    tracked, true = float(m.group(1)), float(m.group(2))
    scale = max(abs(tracked), abs(true))
    if scale == 0.0:
        continue
    rel = abs(tracked - true) / scale
    checked += 1
    if rel > worst:
        worst, worst_line = rel, line.strip()
# 4 significant figures of agreement.
TOL = 1.0e-4
if checked == 0:
    sys.exit(f"[{nprocs} ranks] no Poisson residual rows found in {path}; "
             "is solver_monitoring.poisson.pic_true_residual enabled?")
# A solve that converges on its first iteration never exercises the recurrence,
# so agreement would be vacuous. Keep the fixture tolerances tight enough that
# the solve iterates, and fail loudly if that ever stops being true.
MIN_ROWS = 12
if checked < MIN_ROWS:
    sys.exit(f"[{nprocs} ranks] only {checked} Poisson residual rows (< {MIN_ROWS}); "
             "the solve converged too easily for this check to mean anything. "
             "Tighten poisson_solver tolerances in the regression fixture.")
if worst > TOL:
    sys.exit(f"[{nprocs} ranks] tracked and true Poisson residuals disagree by "
             f"{worst:.3e} (> {TOL:.0e}) over {checked} rows.\n  worst row: {worst_line}\n"
             "  A Krylov method configured at multigrid level_0 makes the\n"
             "  preconditioner nonlinear and decouples the two norms.")
print(f"[{nprocs} ranks] {checked} Poisson rows, worst tracked-vs-true "
      f"relative difference {worst:.3e} (<= {TOL:.0e}).")
PY

  continuity_log="${run_dir}/logs/Continuity_Metrics.log"
  [[ -f "${continuity_log}" ]] || die "missing continuity log '${continuity_log}'."
  python3 - "${continuity_log}" "${nprocs}" <<'PY' || exit 1
import sys
path, nprocs = sys.argv[1], sys.argv[2]
worst, rows = 0.0, 0
for line in open(path):
    parts = [p.strip() for p in line.split("|")]
    if len(parts) < 3:
        continue
    try:
        value = float(parts[2])
    except ValueError:
        continue
    worst = max(worst, abs(value))
    rows += 1
TOL = 1.0e-11
if rows == 0:
    sys.exit(f"[{nprocs} ranks] no divergence rows found in {path}.")
if worst > TOL:
    sys.exit(f"[{nprocs} ranks] max divergence {worst:.3e} exceeds {TOL:.0e}.")
print(f"[{nprocs} ranks] {rows} steps, max divergence {worst:.3e} (<= {TOL:.0e}).")
PY
done

# --- 2. The driven momentum source is constant across a physical timestep. ----
# Two pieces of state sit between the flux measurement and the applied force -
# bulkVelocityCorrection in the handler's PreStep, and a 50/50 EMA in
# ComputeDrivenChannelFlowSource - and BOTH run once per residual evaluation, not
# once per timestep. If either is ungated the applied force walks toward its
# target within the step (0.5, 0.75, 0.875 ... of the way), which is history
# dependence the matrix-free Newton residual forbids and which breaks the
# constant-forcing assumption behind the Picard shadow-Jacobian estimate.
echo "--- driven momentum source is frozen per timestep ---"
force_nprocs="${ranks[0]}"
dest="${tmp_root}/frozen_source"
stage_case constant_flux "${dest}" 0
# Raise verbosity and allow-list the body-force function so it reports the force
# it applies on every call.
sed -i -e "s|^  verbosity: .*|  verbosity: DEBUG|" \
       -e "s|^  enabled_functions: .*|  enabled_functions: [\"ComputeDrivenChannelFlowSource\"]|" \
       "${dest}/monitor.yml"
run_case "${dest}" "${force_nprocs}" > /dev/null

python3 - "${dest}/run.log" "${STEPS}" <<'PY' || exit 1
import re, sys, itertools
path, steps = sys.argv[1], int(sys.argv[2])
values = re.findall(r"New smoothed driving force:\s+(\S+)", open(path).read())
if not values:
    sys.exit("no driving-force records found; the DEBUG allow-list may not have taken effect.")
calls = len(values)
# The force may change ONLY at a physical step boundary, so the number of
# contiguous runs of identical values cannot exceed the number of steps. A step
# whose correction is negligible exits early and contributes no run, hence <=.
#
# Do NOT compare runs against the DISTINCT-value count instead: when the force
# changes on every call, every value is unique and each forms its own run, so
# runs == distinct holds trivially and the check would pass on the broken case.
groups = [(k, len(list(g))) for k, g in itertools.groupby(values)]
if len(groups) > steps:
    sys.exit(f"the applied driving force is NOT constant within a timestep: {calls} evaluations "
             f"over {steps} steps produced {len(groups)} contiguous runs (at most {steps} are "
             f"allowed).\n  run lengths: {[n for _, n in groups][:12]}\n"
             f"  first values: {values[:6]}\n"
             "  Both simCtx->bulkVelocityCorrection (src/BC_Handlers.c) and the smoothing EMA "
             "in ComputeDrivenChannelFlowSource (src/BodyForces.c) must be gated on "
             "simCtx->step; gating only one still lets the force walk toward its target "
             "across residual evaluations.")
if calls < 4 * len(groups):
    sys.exit(f"only {calls} force evaluations across {len(groups)} steps; this check needs several "
             "residual evaluations per step to be meaningful.")
print(f"applied force: {calls} evaluations over {steps} steps -> {len(groups)} contiguous runs "
      f"(lengths {[n for _, n in groups]}); constant within every timestep.")
PY

# --- 3. initial_flux latch, and its survival across a restart. ----------------
echo "--- initial_flux latch and restart persistence ---"
latch_nprocs="${ranks[0]}"
dest="${tmp_root}/initial_flux"
stage_case initial_flux "${dest}" 0
run_dir="$(run_case "${dest}" "${latch_nprocs}")"

read_target() {
  local meta="$1/checkpoint.meta"
  [[ -f "${meta}" ]] || die "missing checkpoint metadata '${meta}'."
  sed -n 's|^-checkpoint_driven_flux_target ||p' "${meta}"
}
read_latched() {
  sed -n 's|^-checkpoint_driven_flux_latched ||p' "$1/checkpoint.meta"
}

final_ckpt="$(ls -d "${run_dir}"/output/checkpoints/step_* | sort | tail -1)"
first_ckpt="$(ls -d "${run_dir}"/output/checkpoints/step_* | sort | head -1)"
latched_target="$(read_target "${final_ckpt}")"

# The shipped laminar case seeds a uniform streamwise velocity U_b = 1 over a
# cross-section of Lx * Ly = 1.0 * 2.0, so the flux the handler must latch is 2.0.
python3 - "${latched_target}" <<'PY' || exit 1
import sys
target, expected = float(sys.argv[1]), 2.0
rel = abs(target - expected) / expected
if rel > 1.0e-9:
    sys.exit(f"initial_flux latched {target!r}, expected the initial-condition "
             f"flux U_b*Lx*Ly = {expected} (relative error {rel:.3e}).")
print(f"initial_flux latched {target} == U_b*Lx*Ly = {expected}.")
PY

[[ "$(read_latched "${first_ckpt}")" == "false" ]] \
  || die "the step-0 checkpoint reports a latched target; the latch must not fire during setup, before the initial condition reaches lUcont."

# Restart from the halfway checkpoint and confirm the target is restored, not
# re-measured. Re-measuring is the subtle failure: by then the flux has drifted,
# so the controller would silently retarget.
restart_step=$(( STEPS / 2 ))
restart_ckpt="${run_dir}/output/checkpoints/$(printf 'step_%012d' "${restart_step}")"
[[ -d "${restart_ckpt}" ]] || die "missing restart checkpoint '${restart_ckpt}'."
restart_target="$(read_target "${restart_ckpt}")"

dest_r="${tmp_root}/initial_flux_restart"
stage_case initial_flux "${dest_r}" "${restart_step}"
restart_run_dir="$(run_case "${dest_r}" "${latch_nprocs}" --restart-from "${run_dir}")"
resumed_ckpt="$(ls -d "${restart_run_dir}"/output/checkpoints/step_* | sort | tail -1)"
resumed_target="$(read_target "${resumed_ckpt}")"

[[ "$(read_latched "${resumed_ckpt}")" == "true" ]] \
  || die "the restarted run reports no latched flux target."
[[ "${resumed_target}" == "${restart_target}" ]] \
  || die "restart re-measured the flux target: checkpoint held ${restart_target}, restarted run holds ${resumed_target}."
echo "initial_flux target ${resumed_target} survived the restart unchanged."

echo "Driven-periodic regression passed."
