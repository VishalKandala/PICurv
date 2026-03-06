#!/usr/bin/env bash

set -euo pipefail

simulator_exe_input="${1:?missing simulator path}"
postprocessor_exe_input="${2:?missing postprocessor path}"
simulator_exe="$(cd "$(dirname "${simulator_exe_input}")" && pwd)/$(basename "${simulator_exe_input}")"
postprocessor_exe="$(cd "$(dirname "${postprocessor_exe_input}")" && pwd)/$(basename "${postprocessor_exe_input}")"
mpi_launcher="${3:?missing MPI launcher}"
nprocs="${4:-1}"
picurv_exe="$(dirname "${simulator_exe}")/picurv"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
valid_fixtures_dir="${repo_root}/tests/fixtures/valid"
tmp_root="$(mktemp -d)"
LAST_RUN_DIR=""
LAST_RUN_LOG=""

cleanup() {
  if [[ "${KEEP_SMOKE_TMP:-0}" == "1" ]]; then
    echo "Smoke debug: preserving temporary workspace at '${tmp_root}'." >&2
    return
  fi
  rm -rf "${tmp_root}"
}
trap cleanup EXIT

die() {
  echo "Smoke failure: $*" >&2
  exit 1
}

require_executable() {
  local path="$1"
  local label="$2"
  if [[ ! -x "${path}" ]]; then
    die "${label} is not executable at '${path}'."
  fi
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "${path}" ]]; then
    die "missing ${label} at '${path}'."
  fi
}

require_dir() {
  local path="$1"
  local label="$2"
  if [[ ! -d "${path}" ]]; then
    die "missing ${label} at '${path}'."
  fi
}

require_count_ge() {
  local search_root="$1"
  local name_pattern="$2"
  local min_count="$3"
  local label="$4"
  local count
  count="$(find "${search_root}" -type f -name "${name_pattern}" | wc -l | tr -d '[:space:]')"
  if [[ "${count}" -lt "${min_count}" ]]; then
    die "expected at least ${min_count} ${label} in '${search_root}' (pattern: ${name_pattern}), found ${count}."
  fi
}

require_file_contains() {
  local file_path="$1"
  local pattern="$2"
  local label="$3"
  if ! grep -q -- "${pattern}" "${file_path}"; then
    die "expected '${file_path}' to contain '${pattern}' (${label})."
  fi
}

run_help_smoke() {
  local exe="$1"
  local banner="$2"
  local out_file="${tmp_root}/$(basename "${exe}").help.out"
  set +e
  "${mpi_launcher}" -n "${nprocs}" "${exe}" -help >"${out_file}" 2>&1
  local rc=$?
  set -e

  if ! grep -q "${banner}" "${out_file}"; then
    echo "Smoke failure: missing banner '${banner}' in ${exe} output." >&2
    sed -n '1,80p' "${out_file}" >&2
    exit 1
  fi

  if [[ ${rc} -ne 0 ]]; then
    echo "Smoke note: ${exe} -help exited ${rc}, but banner check passed."
  fi
}

run_case_init_smoke() {
  local case_dir="${tmp_root}/case-from-init"
  "${picurv_exe}" init flat_channel --dest "${case_dir}" >/dev/null

  require_executable "${case_dir}/picurv" "case-local picurv"
  require_executable "${case_dir}/simulator" "case-local simulator"
  require_executable "${case_dir}/postprocessor" "case-local postprocessor"
  require_file "${case_dir}/.picurv-origin.json" "case origin metadata"

  local status_json="${tmp_root}/status.json"
  "${case_dir}/picurv" status-source --case-dir "${case_dir}" --format json >"${status_json}"
  python3 - "${status_json}" <<'PY'
import json
import sys
path = sys.argv[1]
with open(path, "r", encoding="utf-8") as f:
    payload = json.load(f)
assert "source_repo_root" in payload
assert "binaries" in payload
PY
}

run_dry_run_plan_smoke() {
  local plan_json="${tmp_root}/dry-plan.json"
  "${picurv_exe}" run \
    --solve \
    --post-process \
    --case "${valid_fixtures_dir}/case.yml" \
    --solver "${valid_fixtures_dir}/solver.yml" \
    --monitor "${valid_fixtures_dir}/monitor.yml" \
    --post "${valid_fixtures_dir}/post.yml" \
    --dry-run \
    --format json >"${plan_json}"

  python3 - "${plan_json}" <<'PY'
import json
import sys
path = sys.argv[1]
with open(path, "r", encoding="utf-8") as f:
    payload = json.load(f)
assert payload.get("mode") == "dry-run"
stages = payload.get("stages", {})
assert "solve" in stages
assert "post-process" in stages
assert stages["solve"]["launch_command"][0].endswith("simulator")
assert stages["post-process"]["launch_command"][0].endswith("postprocessor")
PY
}

run_restart_resolution_smoke() {
  local restart_root="${tmp_root}/restart-smoke"
  local prior_run="${restart_root}/old_run"
  local prior_restart_dir="${prior_run}/prior_results"
  local case_cfg="${restart_root}/case.yml"
  local solver_cfg="${restart_root}/solver.yml"
  local monitor_cfg="${restart_root}/monitor.yml"
  local post_cfg="${restart_root}/post.yml"
  local plan_json="${tmp_root}/restart-plan.json"

  mkdir -p "${restart_root}" "${prior_run}/config" "${prior_restart_dir}"
  cp "${valid_fixtures_dir}/case.yml" "${case_cfg}"
  cp "${valid_fixtures_dir}/solver.yml" "${solver_cfg}"
  cp "${valid_fixtures_dir}/monitor.yml" "${monitor_cfg}"
  cp "${valid_fixtures_dir}/post.yml" "${post_cfg}"
  cat >"${prior_run}/config/monitor.yml" <<'YAML'
io:
  directories:
    restart: "prior_results"
YAML

  python3 - "${case_cfg}" "${solver_cfg}" <<'PY'
import sys
import yaml
case_path, solver_path = sys.argv[1], sys.argv[2]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = 5
case_cfg["run_control"]["restart_from_run_dir"] = "old_run"
solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "load"
with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
PY

  (
    cd "${restart_root}"
    "${picurv_exe}" run \
      --solve \
      --post-process \
      --case "${case_cfg}" \
      --solver "${solver_cfg}" \
      --monitor "${monitor_cfg}" \
      --post "${post_cfg}" \
      --dry-run \
      --format json >"${plan_json}"
  )

  python3 - "${plan_json}" "${prior_restart_dir}" <<'PY'
import json
import os
import sys
plan_path = sys.argv[1]
expected_restart = os.path.abspath(sys.argv[2])
with open(plan_path, "r", encoding="utf-8") as f:
    payload = json.load(f)
solve_stage = payload.get("stages", {}).get("solve", {})
assert solve_stage.get("restart_source_directory") == expected_restart
PY
}

prepare_flat_case_les() {
  local case_dir="$1"
  python3 - "${case_dir}/flat_channel.yml" "${case_dir}/Imp-MG-Standard.yml" "${case_dir}/Standard_Output.yml" "${case_dir}/standard_analysis.yml" <<'PY'
import sys
import yaml
case_path, solver_path, monitor_path, post_path = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
with open(monitor_path, "r", encoding="utf-8") as f:
    monitor_cfg = yaml.safe_load(f)
with open(post_path, "r", encoding="utf-8") as f:
    post_cfg = yaml.safe_load(f)

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = 0
case_cfg["run_control"]["total_steps"] = 3
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 8
grid["jm"] = 8
grid["km"] = 16
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}
case_cfg["models"]["physics"]["turbulence"] = {"les": False}

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "solve"

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 3
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/les_smoke"
post_cfg["io"]["output_filename_prefix"] = "Field"

with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
with open(monitor_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(monitor_cfg, f, sort_keys=False)
with open(post_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(post_cfg, f, sort_keys=False)
PY
}

prepare_flat_case_particles_base() {
  local case_dir="$1"
  python3 - "${case_dir}/flat_channel.yml" "${case_dir}/Imp-MG-Standard.yml" "${case_dir}/Standard_Output.yml" "${case_dir}/standard_analysis.yml" <<'PY'
import sys
import yaml
case_path, solver_path, monitor_path, post_path = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
with open(monitor_path, "r", encoding="utf-8") as f:
    monitor_cfg = yaml.safe_load(f)
with open(post_path, "r", encoding="utf-8") as f:
    post_cfg = yaml.safe_load(f)

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = 0
case_cfg["run_control"]["total_steps"] = 3
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 8
grid["jm"] = 8
grid["km"] = 16
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {
    "count": 32,
    "init_mode": "PointSource",
    "restart_mode": "init",
    "point_source": {"x": 0.5, "y": 0.5, "z": 0.5},
}
case_cfg["models"]["physics"]["turbulence"] = {"les": False}

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "solve"

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 3
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/particle_smoke"
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = True

with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
with open(monitor_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(monitor_cfg, f, sort_keys=False)
with open(post_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(post_cfg, f, sort_keys=False)
PY
}

prepare_flat_restart_variant() {
  local case_path="$1"
  local solver_path="$2"
  local post_path="$3"
  local restart_run_dir="$4"
  local particle_restart_mode="$5"
  local output_dir="$6"
  python3 - "${case_path}" "${solver_path}" "${post_path}" "${restart_run_dir}" "${particle_restart_mode}" "${output_dir}" <<'PY'
import sys
import yaml
case_path, solver_path, post_path, restart_dir, particle_mode, output_dir = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
with open(post_path, "r", encoding="utf-8") as f:
    post_cfg = yaml.safe_load(f)

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = 1
case_cfg["run_control"]["total_steps"] = 1
case_cfg["run_control"]["restart_from_run_dir"] = restart_dir
case_cfg.setdefault("models", {}).setdefault("physics", {}).setdefault("particles", {})
case_cfg["models"]["physics"]["particles"]["restart_mode"] = particle_mode

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "load"

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 1
post_cfg["run_control"]["end_step"] = 2
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = output_dir
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = True

with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
with open(post_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(post_cfg, f, sort_keys=False)
PY
}

prepare_brownian_case_analytical() {
  local case_dir="$1"
  python3 - "${case_dir}/brownian_motion.yml" "${case_dir}/Analytical-Zero.yml" "${case_dir}/Standard_Output.yml" "${case_dir}/brownian_analysis.yml" <<'PY'
import sys
import yaml
case_path, solver_path, monitor_path, post_path = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
with open(monitor_path, "r", encoding="utf-8") as f:
    monitor_cfg = yaml.safe_load(f)
with open(post_path, "r", encoding="utf-8") as f:
    post_cfg = yaml.safe_load(f)

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = 0
case_cfg["run_control"]["total_steps"] = 2
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 8
grid["jm"] = 8
grid["km"] = 8
case_cfg.setdefault("models", {}).setdefault("physics", {}).setdefault("particles", {})
case_cfg["models"]["physics"]["particles"]["count"] = 64
case_cfg["models"]["physics"]["particles"]["restart_mode"] = "init"

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "analytical"
solver_cfg["operation_mode"]["analytical_type"] = "ZERO_FLOW"

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 2
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/brownian_smoke"
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = True
post_cfg.setdefault("statistics_pipeline", {})
post_cfg["statistics_pipeline"]["output_prefix"] = "BrownianStats"
post_cfg["statistics_pipeline"]["tasks"] = [{"task": "msd"}]

with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
with open(monitor_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(monitor_cfg, f, sort_keys=False)
with open(post_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(post_cfg, f, sort_keys=False)
PY
}

run_case_workflow() {
  local case_dir="$1"
  local case_file="$2"
  local solver_file="$3"
  local monitor_file="$4"
  local post_file="$5"
  local label="$6"
  local output_log="${tmp_root}/${label}.log"
  local before_runs="${tmp_root}/${label}.runs.before"
  local after_runs="${tmp_root}/${label}.runs.after"
  local created_run=""
  mkdir -p "${case_dir}/runs"
  find "${case_dir}/runs" -mindepth 1 -maxdepth 1 -type d | sort >"${before_runs}"
  (
    cd "${case_dir}"
    ./picurv run \
      --solve \
      --post-process \
      -n "${nprocs}" \
      --case "${case_file}" \
      --solver "${solver_file}" \
      --monitor "${monitor_file}" \
      --post "${post_file}" >"${output_log}" 2>&1
  ) || {
    echo "Smoke failure: workflow '${label}' failed." >&2
    sed -n '1,220p' "${output_log}" >&2
    exit 1
  }

  find "${case_dir}/runs" -mindepth 1 -maxdepth 1 -type d | sort >"${after_runs}"
  created_run="$(comm -13 "${before_runs}" "${after_runs}" | tail -n 1)"
  if [[ -z "${created_run}" ]]; then
    created_run="$(tail -n 1 "${after_runs}")"
  fi
  if [[ -z "${created_run}" ]]; then
    die "workflow '${label}' completed but no run directory was found."
  fi
  LAST_RUN_DIR="${created_run}"
  LAST_RUN_LOG="${output_log}"
}

run_full_runtime_smoke() {
  local flat_les_case="${tmp_root}/flat-les"
  local flat_particles_case="${tmp_root}/flat-particles"
  local brownian_case="${tmp_root}/brownian"

  "${picurv_exe}" init flat_channel --dest "${flat_les_case}" >/dev/null
  prepare_flat_case_les "${flat_les_case}"
  run_case_workflow \
    "${flat_les_case}" \
    "${flat_les_case}/flat_channel.yml" \
    "${flat_les_case}/Imp-MG-Standard.yml" \
    "${flat_les_case}/Standard_Output.yml" \
    "${flat_les_case}/standard_analysis.yml" \
    "flat_les"

  local flat_les_run
  flat_les_run="${LAST_RUN_DIR}"
  require_dir "${flat_les_run}/results" "flat LES results directory"
  require_count_ge "${flat_les_run}/results" "*.dat" 1 "flat LES result data files"
  require_count_ge "${flat_les_run}/viz/les_smoke" "*.vts" 1 "flat LES post VTS files"

  "${picurv_exe}" init flat_channel --dest "${flat_particles_case}" >/dev/null
  prepare_flat_case_particles_base "${flat_particles_case}"
  run_case_workflow \
    "${flat_particles_case}" \
    "${flat_particles_case}/flat_channel.yml" \
    "${flat_particles_case}/Imp-MG-Standard.yml" \
    "${flat_particles_case}/Standard_Output.yml" \
    "${flat_particles_case}/standard_analysis.yml" \
    "flat_particles_base"

  local base_particles_run
  base_particles_run="${LAST_RUN_DIR}"
  require_count_ge "${base_particles_run}/results/particles" "*.dat" 1 "particle snapshot files"
  require_count_ge "${base_particles_run}/viz/particle_smoke" "*.vtp" 1 "particle VTP files"

  local case_restart_load="${flat_particles_case}/case_restart_load.yml"
  local solver_restart_load="${flat_particles_case}/solver_restart_load.yml"
  local post_restart_load="${flat_particles_case}/post_restart_load.yml"
  cp "${flat_particles_case}/flat_channel.yml" "${case_restart_load}"
  cp "${flat_particles_case}/Imp-MG-Standard.yml" "${solver_restart_load}"
  cp "${flat_particles_case}/standard_analysis.yml" "${post_restart_load}"
  prepare_flat_restart_variant "${case_restart_load}" "${solver_restart_load}" "${post_restart_load}" "${base_particles_run}" "load" "viz/restart_load"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_load}" \
    "${solver_restart_load}" \
    "${flat_particles_case}/Standard_Output.yml" \
    "${post_restart_load}" \
    "flat_particles_restart_load"

  local restart_load_run
  restart_load_run="${LAST_RUN_DIR}"
  require_count_ge "${restart_load_run}/results/particles" "*.dat" 1 "restart-load particle snapshots"
  require_file_contains "${LAST_RUN_LOG}" "Particle Restart Mode: load" "restart load branch"

  local case_restart_init="${flat_particles_case}/case_restart_init.yml"
  local solver_restart_init="${flat_particles_case}/solver_restart_init.yml"
  local post_restart_init="${flat_particles_case}/post_restart_init.yml"
  cp "${flat_particles_case}/flat_channel.yml" "${case_restart_init}"
  cp "${flat_particles_case}/Imp-MG-Standard.yml" "${solver_restart_init}"
  cp "${flat_particles_case}/standard_analysis.yml" "${post_restart_init}"
  prepare_flat_restart_variant "${case_restart_init}" "${solver_restart_init}" "${post_restart_init}" "${base_particles_run}" "init" "viz/restart_init"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_init}" \
    "${solver_restart_init}" \
    "${flat_particles_case}/Standard_Output.yml" \
    "${post_restart_init}" \
    "flat_particles_restart_init"

  local restart_init_run
  restart_init_run="${LAST_RUN_DIR}"
  require_count_ge "${restart_init_run}/results/particles" "*.dat" 1 "restart-init particle snapshots"
  require_file_contains "${LAST_RUN_LOG}" "Particle Restart Mode: init" "restart init branch"

  "${picurv_exe}" init brownian_motion --dest "${brownian_case}" >/dev/null
  prepare_brownian_case_analytical "${brownian_case}"
  run_case_workflow \
    "${brownian_case}" \
    "${brownian_case}/brownian_motion.yml" \
    "${brownian_case}/Analytical-Zero.yml" \
    "${brownian_case}/Standard_Output.yml" \
    "${brownian_case}/brownian_analysis.yml" \
    "brownian_analytical"

  local brownian_run
  brownian_run="${LAST_RUN_DIR}"
  require_count_ge "${brownian_run}/viz/brownian_smoke" "*.vts" 1 "brownian eulerian VTS files"
  require_count_ge "${brownian_run}/viz/brownian_smoke" "*.vtp" 1 "brownian particle VTP files"
  require_file "${brownian_run}/BrownianStats_msd.csv" "brownian MSD statistics CSV"
  require_file_contains "${LAST_RUN_LOG}" "Analytical Solution Type" "analytical runtime branch"
}

require_executable "${simulator_exe}" "simulator"
require_executable "${postprocessor_exe}" "postprocessor"
require_executable "${picurv_exe}" "picurv conductor"
python3 -c "import yaml" >/dev/null 2>&1 || die "python dependency 'pyyaml' is required for smoke profile mutation."

echo "==> PICurv smoke: simulator help"
run_help_smoke "${simulator_exe}" "PICurv Simulator"

echo "==> PICurv smoke: postprocessor help"
run_help_smoke "${postprocessor_exe}" "Unified Post-Processing Tool"

echo "==> PICurv smoke: case init and source metadata"
run_case_init_smoke

echo "==> PICurv smoke: dry-run execution plan"
run_dry_run_plan_smoke

echo "==> PICurv smoke: restart source resolution in dry-run plan"
run_restart_resolution_smoke

echo "==> PICurv smoke: full end-to-end runtime sequences"
run_full_runtime_smoke

echo "PICurv smoke completed successfully."
