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

cleanup() {
  rm -rf "${tmp_root}"
}
trap cleanup EXIT

require_executable() {
  local path="$1"
  local label="$2"
  if [[ ! -x "${path}" ]]; then
    echo "Smoke failure: ${label} is not executable at '${path}'." >&2
    exit 1
  fi
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "${path}" ]]; then
    echo "Smoke failure: missing ${label} at '${path}'." >&2
    exit 1
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

require_executable "${simulator_exe}" "simulator"
require_executable "${postprocessor_exe}" "postprocessor"
require_executable "${picurv_exe}" "picurv conductor"

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

echo "PICurv smoke completed successfully."
