#!/usr/bin/env bash

set -euo pipefail

simulator_exe_input="${1:?missing simulator path}"
postprocessor_exe_input="${2:?missing postprocessor path}"
simulator_exe="$(cd "$(dirname "${simulator_exe_input}")" && pwd)/$(basename "${simulator_exe_input}")"
postprocessor_exe="$(cd "$(dirname "${postprocessor_exe_input}")" && pwd)/$(basename "${postprocessor_exe_input}")"
mpi_launcher="${3:?missing MPI launcher}"
nprocs="${4:-1}"
smoke_mode="${5:-standard}"
picurv_exe="$(dirname "${simulator_exe}")/picurv"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
valid_fixtures_dir="${repo_root}/tests/fixtures/valid"
tmp_root="$(mktemp -d)"
LAST_RUN_DIR=""
LAST_RUN_LOG=""
LAST_SOLVER_LOG=""
declare -a mpi_launcher_cmd=()

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

require_file_not_contains() {
  local file_path="$1"
  local pattern="$2"
  local label="$3"
  if grep -q -- "${pattern}" "${file_path}"; then
    die "expected '${file_path}' to omit '${pattern}' (${label})."
  fi
}

parse_mpi_launcher() {
  local launcher="$1"
  local token=""

  mpi_launcher_cmd=()
  while IFS= read -r token; do
    mpi_launcher_cmd+=("${token}")
  done < <(
    python3 - "${launcher}" <<'PY'
import shlex
import sys

launcher = sys.argv[1]
parts = shlex.split(launcher)
if not parts:
    raise SystemExit(1)

for part in parts:
    print(part)
PY
  ) || die "failed to parse MPI launcher '${launcher}'."

  if [[ ${#mpi_launcher_cmd[@]} -eq 0 ]]; then
    die "parsed empty MPI launcher from '${launcher}'."
  fi
}

run_with_mpi_launcher() {
  "${mpi_launcher_cmd[@]}" -n "${nprocs}" "$@"
}

compare_continuity_max_divergence() {
  local first_log="$1"
  local second_log="$2"
  local target_step="$3"
  local tolerance="$4"
  local label="$5"
  python3 - "${first_log}" "${second_log}" "${target_step}" "${tolerance}" "${label}" <<'PY'
import math
import sys

first_log, second_log, target_step_raw, tolerance_raw, label = sys.argv[1:]
target_step = int(target_step_raw)
tolerance = float(tolerance_raw)


def extract_max_divergence(path: str, step: int) -> float:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            if "|" not in raw_line:
                continue
            columns = [token.strip() for token in raw_line.split("|")]
            if len(columns) < 3:
                continue
            if not columns[0].isdigit():
                continue
            if int(columns[0]) != step:
                continue
            return float(columns[2])
    raise RuntimeError(f"{label}: could not find step {step} in {path}")


first_value = extract_max_divergence(first_log, target_step)
second_value = extract_max_divergence(second_log, target_step)
delta = abs(first_value - second_value)

if not math.isfinite(delta):
    raise RuntimeError(f"{label}: non-finite divergence difference for step {target_step}")
if delta > tolerance:
    raise RuntimeError(
        f"{label}: step {target_step} max-divergence mismatch {delta:.6e} exceeds tolerance {tolerance:.6e} "
        f"({first_value:.6e} vs {second_value:.6e})"
    )
PY
}

run_help_smoke() {
  local exe="$1"
  local banner="$2"
  local out_file="${tmp_root}/$(basename "${exe}").help.out"
  set +e
  run_with_mpi_launcher "${exe}" -help >"${out_file}" 2>&1
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

  require_file "${case_dir}/.picurv-origin.json" "case origin metadata"
  require_file "${case_dir}/.picurv-execution.yml" "case runtime execution config"

  if [[ -e "${case_dir}/picurv" || -e "${case_dir}/simulator" || -e "${case_dir}/postprocessor" ]]; then
    echo "Smoke failure: init should not create case-local executables without --pin-binaries." >&2
    exit 1
  fi

  local status_json="${tmp_root}/status.json"
  "${picurv_exe}" status-source --case-dir "${case_dir}" --format json >"${status_json}"
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

run_template_matrix_smoke() {
  local matrix_root="${tmp_root}/template-matrix"
  mkdir -p "${matrix_root}"

  python3 - "${picurv_exe}" "${matrix_root}" <<'PY'
import glob
import json
import os
import subprocess
import sys
import yaml

picurv_exe = os.path.abspath(sys.argv[1])
matrix_root = os.path.abspath(sys.argv[2])
templates = ["flat_channel", "bent_channel", "brownian_motion"]


def fail(message: str) -> None:
    raise RuntimeError(message)


def discover_case_bundle(case_dir: str, template_name: str):
    ymls = sorted(os.path.basename(path) for path in glob.glob(os.path.join(case_dir, "*.yml")))
    case_file = f"{template_name}.yml"
    if case_file not in ymls:
        fail(f"{template_name}: expected case file '{case_file}' in {case_dir}, found {ymls}")

    monitor_file = "Standard_Output.yml" if "Standard_Output.yml" in ymls else None
    if not monitor_file:
        fail(f"{template_name}: expected monitor file 'Standard_Output.yml' in {case_dir}")

    post_candidates = [name for name in ymls if name.lower().endswith("_analysis.yml")]
    if not post_candidates:
        fail(f"{template_name}: expected *_analysis.yml post profile in {case_dir}")
    post_file = sorted(post_candidates)[0]

    excluded = {
        case_file,
        monitor_file,
        post_file,
        "slurm_cluster.yml",
        "grid_independence_study.yml",
        "timestep_sensitivity_study.yml",
    }

    def is_execution_example(name: str) -> bool:
        path = os.path.join(case_dir, name)
        with open(path, "r", encoding="utf-8") as f:
            payload = yaml.safe_load(f) or {}
        if not isinstance(payload, dict):
            return False
        return any(key in payload for key in ("default_execution", "local_execution", "cluster_execution"))

    solver_candidates = [
        name for name in ymls
        if (
            name not in excluded
            and "study" not in name.lower()
            and "cluster" not in name.lower()
            and not is_execution_example(name)
        )
    ]
    if len(solver_candidates) != 1:
        fail(
            f"{template_name}: expected exactly one solver profile, found {solver_candidates} "
            f"(all yml: {ymls})"
        )
    return case_file, solver_candidates[0], monitor_file, post_file


for template_name in templates:
    case_dir = os.path.join(matrix_root, template_name)
    init_cmd = [picurv_exe, "init", template_name, "--dest", case_dir]
    init_res = subprocess.run(init_cmd, text=True, capture_output=True)
    if init_res.returncode != 0:
        fail(f"{template_name}: init failed:\n{init_res.stdout}\n{init_res.stderr}")

    case_file, solver_file, monitor_file, post_file = discover_case_bundle(case_dir, template_name)

    validate_cmd = [
        picurv_exe,
        "validate",
        "--case", os.path.join(case_dir, case_file),
        "--solver", os.path.join(case_dir, solver_file),
        "--monitor", os.path.join(case_dir, monitor_file),
        "--post", os.path.join(case_dir, post_file),
    ]
    validate_res = subprocess.run(validate_cmd, cwd=case_dir, text=True, capture_output=True)
    if validate_res.returncode != 0:
        fail(f"{template_name}: validate failed:\n{validate_res.stdout}\n{validate_res.stderr}")

    dry_cmd = [
        picurv_exe,
        "run",
        "--solve",
        "--post-process",
        "--case", os.path.join(case_dir, case_file),
        "--solver", os.path.join(case_dir, solver_file),
        "--monitor", os.path.join(case_dir, monitor_file),
        "--post", os.path.join(case_dir, post_file),
        "--dry-run",
        "--format", "json",
    ]
    dry_res = subprocess.run(dry_cmd, cwd=case_dir, text=True, capture_output=True)
    if dry_res.returncode != 0:
        fail(f"{template_name}: dry-run failed:\n{dry_res.stdout}\n{dry_res.stderr}")

    payload = json.loads(dry_res.stdout)
    stages = payload.get("stages", {})
    if "solve" not in stages or "post-process" not in stages:
        fail(f"{template_name}: dry-run payload missing solve/post-process stages: {payload}")
    solve_launch = stages["solve"].get("launch_command", [])
    post_launch = stages["post-process"].get("launch_command", [])
    if not solve_launch or not any(str(token).endswith("simulator") for token in solve_launch):
        fail(f"{template_name}: solve launch command does not target simulator: {solve_launch}")
    if not post_launch or not any(str(token).endswith("postprocessor") for token in post_launch):
        fail(f"{template_name}: post launch command does not target postprocessor: {post_launch}")

print("template matrix smoke validated: " + ", ".join(templates))
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
assert any(str(token).endswith("simulator") for token in stages["solve"]["launch_command"])
assert any(str(token).endswith("postprocessor") for token in stages["post-process"]["launch_command"])
PY
}

run_restart_resolution_smoke() {
  local restart_root="${tmp_root}/restart-smoke"
  local prior_run="${restart_root}/old_run"
  local prior_output_dir="${prior_run}/output"
  local case_cfg="${restart_root}/case.yml"
  local solver_cfg="${restart_root}/solver.yml"
  local monitor_cfg="${restart_root}/monitor.yml"
  local plan_json="${tmp_root}/restart-plan.json"

  mkdir -p "${restart_root}" "${prior_run}/config" "${prior_output_dir}/eulerian"
  cp "${valid_fixtures_dir}/case.yml" "${case_cfg}"
  cp "${valid_fixtures_dir}/solver.yml" "${solver_cfg}"
  cp "${valid_fixtures_dir}/monitor.yml" "${monitor_cfg}"
  cat >"${prior_run}/config/monitor.yml" <<'YAML'
io:
  directories:
    output: "output"
    restart: "restart"
YAML

  # Create fake step files for load mode validation (steps 5-15)
  for step in $(seq 5 15); do
    printf -v stepfmt "%05d" "${step}"
    touch "${prior_output_dir}/eulerian/ufield${stepfmt}_0.dat"
  done

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
case_cfg["run_control"]["total_steps"] = 10
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
      --restart-from "${prior_run}" \
      --case "${case_cfg}" \
      --solver "${solver_cfg}" \
      --monitor "${monitor_cfg}" \
      --dry-run \
      --format json >"${plan_json}"
  )

  python3 - "${plan_json}" "${prior_output_dir}" <<'PY'
import json
import os
import sys
plan_path = sys.argv[1]
expected_restart = os.path.abspath(sys.argv[2])
with open(plan_path, "r", encoding="utf-8") as f:
    payload = json.load(f)
solve_stage = payload.get("stages", {}).get("solve", {})
assert solve_stage.get("restart_source_directory") == expected_restart, \
    f"Expected {expected_restart}, got {solve_stage.get('restart_source_directory')}"
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

prepare_flat_case_particles_corner_averaged() {
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
solver_cfg["interpolation"] = {"method": "CornerAveraged"}

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 3
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/particle_corner_averaged_smoke"
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

prepare_flat_restart_equivalence_case() {
  local case_dir="$1"
  local start_step="$2"
  local total_steps="$3"
  local eulerian_source="$4"
  local restart_run_dir="$5"
  local post_output_dir="$6"
  local post_start_step="$7"
  local post_end_step="$8"
  python3 - \
    "${case_dir}/flat_channel.yml" \
    "${case_dir}/Imp-MG-Standard.yml" \
    "${case_dir}/Standard_Output.yml" \
    "${case_dir}/standard_analysis.yml" \
    "${start_step}" \
    "${total_steps}" \
    "${eulerian_source}" \
    "${restart_run_dir}" \
    "${post_output_dir}" \
    "${post_start_step}" \
    "${post_end_step}" <<'PY'
import sys
import yaml

(
    case_path,
    solver_path,
    monitor_path,
    post_path,
    start_step,
    total_steps,
    eulerian_source,
    restart_run_dir,
    post_output_dir,
    post_start_step,
    post_end_step,
) = sys.argv[1:]

with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
with open(monitor_path, "r", encoding="utf-8") as f:
    monitor_cfg = yaml.safe_load(f)
with open(post_path, "r", encoding="utf-8") as f:
    post_cfg = yaml.safe_load(f)

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = int(start_step)
case_cfg["run_control"]["total_steps"] = int(total_steps)
case_cfg["run_control"]["dt_physical"] = 0.001
# restart_from_run_dir removed from case.yml; restart is now handled via
# --restart-from CLI flag. The restart_run_dir variable is passed to picurv
# via --restart-from on the command line instead.

case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 8
grid["jm"] = 8
grid["km"] = 16

case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}
case_cfg["models"]["physics"]["turbulence"] = {"les": False}

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = eulerian_source

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = int(post_start_step)
post_cfg["run_control"]["end_step"] = int(post_end_step)
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = post_output_dir
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = False

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

prepare_bent_case_tiny() {
  local case_dir="$1"
  python3 - "${case_dir}/bent_channel.yml" "${case_dir}/Imp-MG-Standard.yml" "${case_dir}/Standard_Output.yml" "${case_dir}/standard_analysis.yml" <<'PY'
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
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}

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
post_cfg["io"]["output_directory"] = "viz/bent_smoke"
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

prepare_flat_case_particles_stress() {
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
case_cfg["run_control"]["total_steps"] = 6
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 10
grid["jm"] = 10
grid["km"] = 20
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {
    "count": 96,
    "init_mode": "PointSource",
    "restart_mode": "init",
    "point_source": {"x": 0.5, "y": 0.5, "z": 0.5},
}
case_cfg["models"]["physics"]["turbulence"] = {"les": False}

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "solve"

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 2
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 2

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 6
post_cfg["run_control"]["step_interval"] = 2
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/particle_stress"
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

prepare_flat_case_parabolic_stress() {
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
case_cfg["run_control"]["total_steps"] = 4
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 8
grid["jm"] = 8
grid["km"] = 16
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}
case_cfg["boundary_conditions"] = [
    {"face": "-Xi", "type": "WALL", "handler": "noslip"},
    {"face": "+Xi", "type": "WALL", "handler": "noslip"},
    {"face": "-Eta", "type": "WALL", "handler": "noslip"},
    {"face": "+Eta", "type": "WALL", "handler": "noslip"},
    {"face": "-Zeta", "type": "INLET", "handler": "parabolic", "params": {"v_max": 1.5}},
    {"face": "+Zeta", "type": "OUTLET", "handler": "conservation"},
]

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "solve"

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 4
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/parabolic_stress"
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = False

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

prepare_flat_case_periodic_flux_stress() {
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
case_cfg["run_control"]["total_steps"] = 4
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("grid", {}).setdefault("programmatic_settings", {})
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"] = 8
grid["jm"] = 8
grid["km"] = 16
case_cfg.setdefault("models", {}).setdefault("domain", {})
case_cfg["models"]["domain"]["blocks"] = 1
case_cfg["models"]["domain"]["i_periodic"] = False
case_cfg["models"]["domain"]["j_periodic"] = False
case_cfg["models"]["domain"]["k_periodic"] = True
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}
case_cfg["boundary_conditions"] = [
    {"face": "-Xi", "type": "WALL", "handler": "noslip"},
    {"face": "+Xi", "type": "WALL", "handler": "noslip"},
    {"face": "-Eta", "type": "WALL", "handler": "noslip"},
    {"face": "+Eta", "type": "WALL", "handler": "noslip"},
    {"face": "-Zeta", "type": "PERIODIC", "handler": "constant_flux", "params": {"target_flux": 1.0, "apply_trim": True}},
    {"face": "+Zeta", "type": "PERIODIC", "handler": "constant_flux", "params": {"target_flux": 1.0, "apply_trim": True}},
]

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "solve"

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 4
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/periodic_flux_stress"
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = False

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
  local restart_from="${7:-}"  # Optional: --restart-from <run_dir>
  local output_log="${tmp_root}/${label}.log"
  local before_runs="${tmp_root}/${label}.runs.before"
  local after_runs="${tmp_root}/${label}.runs.after"
  local created_run=""
  local restart_args=()
  if [[ -n "${restart_from}" ]]; then
    restart_args=("--restart-from" "${restart_from}")
  fi
  mkdir -p "${case_dir}/runs"
  find "${case_dir}/runs" -mindepth 1 -maxdepth 1 -type d | sort >"${before_runs}"
  (
    cd "${case_dir}"
    "${picurv_exe}" run \
      --solve \
      --post-process \
      -n "${nprocs}" \
      --case "${case_file}" \
      --solver "${solver_file}" \
      --monitor "${monitor_file}" \
      --post "${post_file}" \
      "${restart_args[@]}" >"${output_log}" 2>&1
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
  LAST_SOLVER_LOG="$(find "${created_run}/scheduler" -maxdepth 1 -type f -name '*_solver.log' | sort | tail -n 1)"
  if [[ -z "${LAST_SOLVER_LOG}" ]]; then
    die "workflow '${label}' completed but no solver runtime log was found under '${created_run}/scheduler'."
  fi
}

run_restart_equivalence_smoke() {
  local continuous_case="${tmp_root}/restart-equivalence-continuous"
  local split_case="${tmp_root}/restart-equivalence-split"
  local continuous_run=""
  local split_base_run=""
  local split_restart_run=""
  local continuous_continuity_log=""
  local split_continuity_log=""

  "${picurv_exe}" init flat_channel --dest "${continuous_case}" >/dev/null
  prepare_flat_restart_equivalence_case "${continuous_case}" 0 4 "solve" "" "viz/restart_equivalence_continuous" 0 4
  run_case_workflow \
    "${continuous_case}" \
    "${continuous_case}/flat_channel.yml" \
    "${continuous_case}/Imp-MG-Standard.yml" \
    "${continuous_case}/Standard_Output.yml" \
    "${continuous_case}/standard_analysis.yml" \
    "restart_equiv_continuous"

  continuous_run="${LAST_RUN_DIR}"
  continuous_continuity_log="${continuous_run}/logs/Continuity_Metrics.log"
  require_file "${continuous_continuity_log}" "continuous continuity metrics log"

  "${picurv_exe}" init flat_channel --dest "${split_case}" >/dev/null
  prepare_flat_restart_equivalence_case "${split_case}" 0 3 "solve" "" "viz/restart_equivalence_split_base" 0 3
  run_case_workflow \
    "${split_case}" \
    "${split_case}/flat_channel.yml" \
    "${split_case}/Imp-MG-Standard.yml" \
    "${split_case}/Standard_Output.yml" \
    "${split_case}/standard_analysis.yml" \
    "restart_equiv_split_base"

  split_base_run="${LAST_RUN_DIR}"
  prepare_flat_restart_equivalence_case "${split_case}" 2 2 "solve" "${split_base_run}" "viz/restart_equivalence_split_restart" 2 4
  run_case_workflow \
    "${split_case}" \
    "${split_case}/flat_channel.yml" \
    "${split_case}/Imp-MG-Standard.yml" \
    "${split_case}/Standard_Output.yml" \
    "${split_case}/standard_analysis.yml" \
    "restart_equiv_split_restart" \
    "${split_base_run}"

  split_restart_run="${LAST_RUN_DIR}"
  split_continuity_log="${split_restart_run}/logs/Continuity_Metrics.log"
  require_file "${split_continuity_log}" "restart continuity metrics log"
  compare_continuity_max_divergence \
    "${continuous_continuity_log}" \
    "${split_continuity_log}" \
    "4" \
    "5.0e-10" \
    "restart-equivalence continuity max divergence"
}

run_full_runtime_smoke() {
  local flat_les_case="${tmp_root}/flat-les"
  local bent_case="${tmp_root}/bent"
  local flat_particles_case="${tmp_root}/flat-particles"
  local flat_particles_ca_case="${tmp_root}/flat-particles-ca"
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
  require_dir "${flat_les_run}/output/eulerian" "flat LES Eulerian output directory"
  require_count_ge "${flat_les_run}/output/eulerian" "*.dat" 1 "flat LES Eulerian data files"
  require_count_ge "${flat_les_run}/viz/les_smoke" "*.vts" 1 "flat LES post VTS files"
  require_file_contains "${LAST_SOLVER_LOG}" "Run Mode                   : Full Simulation" "runtime banner run mode"
  require_file_contains "${LAST_SOLVER_LOG}" "Field/Restart Cadence      : every 1 step(s)" "runtime banner field cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Immersed Boundary          : DISABLED" "runtime banner immersed-boundary state"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 0" "runtime banner particle count"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence" "runtime banner particle console cadence omission"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling" "runtime banner particle row sampling omission"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Initialization Mode" "runtime banner particle init omission"

  "${picurv_exe}" init bent_channel --dest "${bent_case}" >/dev/null
  prepare_bent_case_tiny "${bent_case}"
  run_case_workflow \
    "${bent_case}" \
    "${bent_case}/bent_channel.yml" \
    "${bent_case}/Imp-MG-Standard.yml" \
    "${bent_case}/Standard_Output.yml" \
    "${bent_case}/standard_analysis.yml" \
    "bent_les"

  local bent_run
  bent_run="${LAST_RUN_DIR}"
  require_dir "${bent_run}/output/eulerian" "bent Eulerian output directory"
  require_count_ge "${bent_run}/output/eulerian" "*.dat" 1 "bent Eulerian data files"
  require_count_ge "${bent_run}/viz/bent_smoke" "*.vts" 1 "bent post VTS files"

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
  require_count_ge "${base_particles_run}/output/particles" "*.dat" 1 "particle snapshot files"
  require_count_ge "${base_particles_run}/viz/particle_smoke" "*.vtp" 1 "particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 32" "particle runtime banner particle count"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence   : DISABLED" "particle runtime banner disabled particle console cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling  : every 1 particle(s)" "particle runtime banner row sampling"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Initialization Mode: Point Source" "particle runtime banner initialization mode"
  require_file_contains "${LAST_SOLVER_LOG}" "Interpolation Method       : Trilinear (direct cell-center)" "particle runtime banner default interpolation method"

  "${picurv_exe}" init flat_channel --dest "${flat_particles_ca_case}" >/dev/null
  prepare_flat_case_particles_corner_averaged "${flat_particles_ca_case}"
  run_case_workflow \
    "${flat_particles_ca_case}" \
    "${flat_particles_ca_case}/flat_channel.yml" \
    "${flat_particles_ca_case}/Imp-MG-Standard.yml" \
    "${flat_particles_ca_case}/Standard_Output.yml" \
    "${flat_particles_ca_case}/standard_analysis.yml" \
    "flat_particles_corner_averaged"

  require_count_ge "${LAST_RUN_DIR}/output/particles" "*.dat" 1 "corner-averaged particle snapshot files"
  require_count_ge "${LAST_RUN_DIR}/viz/particle_corner_averaged_smoke" "*.vtp" 1 "corner-averaged particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Interpolation Method       : CornerAveraged (legacy)" "corner-averaged runtime banner interpolation method"

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
    "flat_particles_restart_load" \
    "${base_particles_run}"

  local restart_load_run
  restart_load_run="${LAST_RUN_DIR}"
  require_count_ge "${restart_load_run}/output/particles" "*.dat" 1 "restart-load particle snapshots"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : load" "restart load branch"

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
    "flat_particles_restart_init" \
    "${base_particles_run}"

  local restart_init_run
  restart_init_run="${LAST_RUN_DIR}"
  require_count_ge "${restart_init_run}/output/particles" "*.dat" 1 "restart-init particle snapshots"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : init" "restart init branch"

  run_restart_equivalence_smoke

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
  require_file_contains "${LAST_SOLVER_LOG}" "Analytical Solution Type" "analytical runtime branch"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 64" "brownian runtime banner particle count"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence   : DISABLED" "brownian runtime banner disabled particle console cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling  : every 1 particle(s)" "brownian runtime banner row sampling"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Initialization Mode: Point Source" "brownian runtime banner initialization mode"
}

run_multi_rank_runtime_smoke() {
  local flat_case="${tmp_root}/flat-mpi"
  local bent_case="${tmp_root}/bent-mpi"
  local flat_particles_case="${tmp_root}/flat-particles-mpi"

  "${picurv_exe}" init flat_channel --dest "${flat_case}" >/dev/null
  prepare_flat_case_les "${flat_case}"
  run_case_workflow \
    "${flat_case}" \
    "${flat_case}/flat_channel.yml" \
    "${flat_case}/Imp-MG-Standard.yml" \
    "${flat_case}/Standard_Output.yml" \
    "${flat_case}/standard_analysis.yml" \
    "flat_mpi"

  local flat_run
  flat_run="${LAST_RUN_DIR}"
  require_count_ge "${flat_run}/output/eulerian" "*.dat" 1 "flat MPI Eulerian data files"
  require_count_ge "${flat_run}/viz/les_smoke" "*.vts" 1 "flat MPI post VTS files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "flat MPI rank count in runtime summary"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 0" "flat MPI runtime banner particle count"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence" "flat MPI runtime banner particle console cadence omission"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling" "flat MPI runtime banner particle row sampling omission"

  "${picurv_exe}" init bent_channel --dest "${bent_case}" >/dev/null
  prepare_bent_case_tiny "${bent_case}"
  run_case_workflow \
    "${bent_case}" \
    "${bent_case}/bent_channel.yml" \
    "${bent_case}/Imp-MG-Standard.yml" \
    "${bent_case}/Standard_Output.yml" \
    "${bent_case}/standard_analysis.yml" \
    "bent_mpi"

  local bent_run
  bent_run="${LAST_RUN_DIR}"
  require_count_ge "${bent_run}/output/eulerian" "*.dat" 1 "bent MPI Eulerian data files"
  require_count_ge "${bent_run}/viz/bent_smoke" "*.vts" 1 "bent MPI post VTS files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "bent MPI rank count in runtime summary"

  "${picurv_exe}" init flat_channel --dest "${flat_particles_case}" >/dev/null
  prepare_flat_case_particles_base "${flat_particles_case}"
  run_case_workflow \
    "${flat_particles_case}" \
    "${flat_particles_case}/flat_channel.yml" \
    "${flat_particles_case}/Imp-MG-Standard.yml" \
    "${flat_particles_case}/Standard_Output.yml" \
    "${flat_particles_case}/standard_analysis.yml" \
    "flat_particles_mpi_base"

  local base_particles_run
  base_particles_run="${LAST_RUN_DIR}"
  require_count_ge "${base_particles_run}/output/particles" "*.dat" 1 "MPI particle snapshot files"
  require_count_ge "${base_particles_run}/viz/particle_smoke" "*.vtp" 1 "MPI particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "MPI particle run rank count in runtime summary"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 32" "MPI particle runtime banner particle count"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence   : DISABLED" "MPI particle runtime banner disabled particle console cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling  : every 1 particle(s)" "MPI particle runtime banner row sampling"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Initialization Mode: Point Source" "MPI particle runtime banner initialization mode"

  local case_restart_load="${flat_particles_case}/case_restart_load.yml"
  local solver_restart_load="${flat_particles_case}/solver_restart_load.yml"
  local post_restart_load="${flat_particles_case}/post_restart_load.yml"
  cp "${flat_particles_case}/flat_channel.yml" "${case_restart_load}"
  cp "${flat_particles_case}/Imp-MG-Standard.yml" "${solver_restart_load}"
  cp "${flat_particles_case}/standard_analysis.yml" "${post_restart_load}"
  prepare_flat_restart_variant "${case_restart_load}" "${solver_restart_load}" "${post_restart_load}" "${base_particles_run}" "load" "viz/restart_load_mpi"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_load}" \
    "${solver_restart_load}" \
    "${flat_particles_case}/Standard_Output.yml" \
    "${post_restart_load}" \
    "flat_particles_mpi_restart_load" \
    "${base_particles_run}"
  require_count_ge "${LAST_RUN_DIR}/output/particles" "*.dat" 1 "MPI restart-load particle snapshots"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : load" "MPI restart load branch"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "MPI restart-load rank count in runtime summary"

  local case_restart_init="${flat_particles_case}/case_restart_init.yml"
  local solver_restart_init="${flat_particles_case}/solver_restart_init.yml"
  local post_restart_init="${flat_particles_case}/post_restart_init.yml"
  cp "${flat_particles_case}/flat_channel.yml" "${case_restart_init}"
  cp "${flat_particles_case}/Imp-MG-Standard.yml" "${solver_restart_init}"
  cp "${flat_particles_case}/standard_analysis.yml" "${post_restart_init}"
  prepare_flat_restart_variant "${case_restart_init}" "${solver_restart_init}" "${post_restart_init}" "${base_particles_run}" "init" "viz/restart_init_mpi"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_init}" \
    "${solver_restart_init}" \
    "${flat_particles_case}/Standard_Output.yml" \
    "${post_restart_init}" \
    "flat_particles_mpi_restart_init" \
    "${base_particles_run}"
  require_count_ge "${LAST_RUN_DIR}/output/particles" "*.dat" 1 "MPI restart-init particle snapshots"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : init" "MPI restart init branch"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "MPI restart-init rank count in runtime summary"
}

run_stress_smoke() {
  local particle_case="${tmp_root}/flat-particles-stress"
  local restart_case="${tmp_root}/restart-chain-stress"
  local parabolic_case="${tmp_root}/flat-parabolic-stress"
  local periodic_case="${tmp_root}/flat-periodic-flux-stress"
  local mpi_particle_case="${tmp_root}/flat-particles-stress-mpi"
  local saved_nprocs="${nprocs}"
  local stress_mpi_nprocs=2
  local periodic_nprocs=2

  "${picurv_exe}" init flat_channel --dest "${particle_case}" >/dev/null
  prepare_flat_case_particles_stress "${particle_case}"
  run_case_workflow \
    "${particle_case}" \
    "${particle_case}/flat_channel.yml" \
    "${particle_case}/Imp-MG-Standard.yml" \
    "${particle_case}/Standard_Output.yml" \
    "${particle_case}/standard_analysis.yml" \
    "flat_particles_stress"
  require_count_ge "${LAST_RUN_DIR}/output/particles" "*.dat" 2 "stress particle snapshots"
  require_count_ge "${LAST_RUN_DIR}/viz/particle_stress" "*.vtp" 1 "stress particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 96" "stress particle count in runtime summary"

  "${picurv_exe}" init flat_channel --dest "${restart_case}" >/dev/null
  prepare_flat_restart_equivalence_case "${restart_case}" 0 3 "solve" "" "viz/restart_chain_base" 0 3
  run_case_workflow \
    "${restart_case}" \
    "${restart_case}/flat_channel.yml" \
    "${restart_case}/Imp-MG-Standard.yml" \
    "${restart_case}/Standard_Output.yml" \
    "${restart_case}/standard_analysis.yml" \
    "restart_chain_base"
  local restart_base_run="${LAST_RUN_DIR}"
  prepare_flat_restart_equivalence_case "${restart_case}" 2 2 "solve" "${restart_base_run}" "viz/restart_chain_mid" 2 4
  run_case_workflow \
    "${restart_case}" \
    "${restart_case}/flat_channel.yml" \
    "${restart_case}/Imp-MG-Standard.yml" \
    "${restart_case}/Standard_Output.yml" \
    "${restart_case}/standard_analysis.yml" \
    "restart_chain_mid" \
    "${restart_base_run}"
  local restart_mid_run="${LAST_RUN_DIR}"
  prepare_flat_restart_equivalence_case "${restart_case}" 4 2 "solve" "${restart_mid_run}" "viz/restart_chain_final" 4 6
  run_case_workflow \
    "${restart_case}" \
    "${restart_case}/flat_channel.yml" \
    "${restart_case}/Imp-MG-Standard.yml" \
    "${restart_case}/Standard_Output.yml" \
    "${restart_case}/standard_analysis.yml" \
    "restart_chain_final" \
    "${restart_mid_run}"
  require_file "${LAST_RUN_DIR}/logs/Continuity_Metrics.log" "restart-chain continuity metrics log"

  "${picurv_exe}" init flat_channel --dest "${parabolic_case}" >/dev/null
  prepare_flat_case_parabolic_stress "${parabolic_case}"
  run_case_workflow \
    "${parabolic_case}" \
    "${parabolic_case}/flat_channel.yml" \
    "${parabolic_case}/Imp-MG-Standard.yml" \
    "${parabolic_case}/Standard_Output.yml" \
    "${parabolic_case}/standard_analysis.yml" \
    "flat_parabolic_stress"
  require_count_ge "${LAST_RUN_DIR}/viz/parabolic_stress" "*.vts" 1 "parabolic stress VTS files"

  if [[ "${saved_nprocs}" -gt 1 ]]; then
    stress_mpi_nprocs=$((saved_nprocs + 1))
    periodic_nprocs="${saved_nprocs}"
  fi

  nprocs="${periodic_nprocs}"
  "${picurv_exe}" init flat_channel --dest "${periodic_case}" >/dev/null
  prepare_flat_case_periodic_flux_stress "${periodic_case}"
  (
    cd "${periodic_case}"
    "${picurv_exe}" validate \
      --case "${periodic_case}/flat_channel.yml" \
      --solver "${periodic_case}/Imp-MG-Standard.yml" \
      --monitor "${periodic_case}/Standard_Output.yml" \
      --post "${periodic_case}/standard_analysis.yml" >/dev/null
    "${picurv_exe}" run \
      --solve \
      --post-process \
      -n "${periodic_nprocs}" \
      --case "${periodic_case}/flat_channel.yml" \
      --solver "${periodic_case}/Imp-MG-Standard.yml" \
      --monitor "${periodic_case}/Standard_Output.yml" \
      --post "${periodic_case}/standard_analysis.yml" \
      --dry-run \
      --format json >"${tmp_root}/flat_periodic_flux_stress_plan.json"
  )
  python3 - "${tmp_root}/flat_periodic_flux_stress_plan.json" <<'PY'
import json
import sys
path = sys.argv[1]
with open(path, "r", encoding="utf-8") as f:
    payload = json.load(f)
stages = payload.get("stages", {})
assert payload.get("mode") == "dry-run"
assert "solve" in stages and "post-process" in stages
solve_launch = [str(item) for item in stages["solve"].get("launch_command", [])]
post_launch = [str(item) for item in stages["post-process"].get("launch_command", [])]
assert any(item.endswith("simulator") for item in solve_launch)
assert any(item.endswith("postprocessor") for item in post_launch)
PY

  nprocs="${stress_mpi_nprocs}"
  "${picurv_exe}" init flat_channel --dest "${mpi_particle_case}" >/dev/null
  prepare_flat_case_particles_stress "${mpi_particle_case}"
  python3 - "${mpi_particle_case}/standard_analysis.yml" <<'PY'
import sys
import yaml
post_path = sys.argv[1]
with open(post_path, "r", encoding="utf-8") as f:
    post_cfg = yaml.safe_load(f)
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/particle_stress_mpi"
with open(post_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(post_cfg, f, sort_keys=False)
PY
  run_case_workflow \
    "${mpi_particle_case}" \
    "${mpi_particle_case}/flat_channel.yml" \
    "${mpi_particle_case}/Imp-MG-Standard.yml" \
    "${mpi_particle_case}/Standard_Output.yml" \
    "${mpi_particle_case}/standard_analysis.yml" \
    "flat_particles_stress_mpi"
  nprocs="${saved_nprocs}"
  require_count_ge "${LAST_RUN_DIR}/viz/particle_stress_mpi" "*.vtp" 1 "MPI particle stress VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${stress_mpi_nprocs}" "MPI stress rank count in runtime summary"
}

run_periodic_dev_smoke() {
  local periodic_case="${tmp_root}/flat-periodic-flux-dev"
  local saved_nprocs="${nprocs}"
  local periodic_nprocs="${nprocs}"

  if [[ "${periodic_nprocs}" -lt 2 ]]; then
    periodic_nprocs=2
  fi

  nprocs="${periodic_nprocs}"
  "${picurv_exe}" init flat_channel --dest "${periodic_case}" >/dev/null
  prepare_flat_case_periodic_flux_stress "${periodic_case}"
  run_case_workflow \
    "${periodic_case}" \
    "${periodic_case}/flat_channel.yml" \
    "${periodic_case}/Imp-MG-Standard.yml" \
    "${periodic_case}/Standard_Output.yml" \
    "${periodic_case}/standard_analysis.yml" \
    "flat_periodic_flux_dev"
  require_file "${LAST_RUN_DIR}/logs/Continuity_Metrics.log" "periodic dev continuity metrics log"
  require_count_ge "${LAST_RUN_DIR}/viz/periodic_flux_stress" "*.vts" 1 "periodic dev post VTS files"
  nprocs="${saved_nprocs}"
}

require_executable "${simulator_exe}" "simulator"
require_executable "${postprocessor_exe}" "postprocessor"
require_executable "${picurv_exe}" "picurv conductor"
parse_mpi_launcher "${mpi_launcher}"
# Export so that 'picurv run' also uses the same MPI launcher
# (picurv checks PICURV_MPI_LAUNCHER -> MPI_LAUNCHER env -> .picurv-execution.yml -> fallback mpiexec).
export PICURV_MPI_LAUNCHER="${mpi_launcher}"
python3 -c "import yaml" >/dev/null 2>&1 || die "python dependency 'pyyaml' is required for smoke profile mutation."
case "${smoke_mode}" in
  standard|stress|periodic-dev) ;;
  *) die "unknown smoke mode '${smoke_mode}' (expected 'standard', 'stress', or 'periodic-dev')" ;;
esac

echo "==> PICurv smoke: simulator help"
run_help_smoke "${simulator_exe}" "PICurv Simulator"

echo "==> PICurv smoke: postprocessor help"
run_help_smoke "${postprocessor_exe}" "Unified Post-Processing Tool"

if [[ "${smoke_mode}" == "periodic-dev" ]]; then
  echo "==> PICurv smoke: periodic development runtime harness"
  run_periodic_dev_smoke
  echo "PICurv periodic-dev smoke completed successfully."
  exit 0
fi

echo "==> PICurv smoke: case init and source metadata"
run_case_init_smoke

echo "==> PICurv smoke: template matrix init/validate/dry-run checks"
run_template_matrix_smoke

echo "==> PICurv smoke: dry-run execution plan"
run_dry_run_plan_smoke

echo "==> PICurv smoke: restart source resolution in dry-run plan"
run_restart_resolution_smoke

if [[ "${nprocs}" -gt 1 ]]; then
  echo "==> PICurv smoke: multi-rank runtime sequences (flat+bent)"
  run_multi_rank_runtime_smoke
else
  echo "==> PICurv smoke: full end-to-end runtime sequences"
  run_full_runtime_smoke
fi

if [[ "${smoke_mode}" == "stress" ]]; then
  echo "==> PICurv smoke: stress extensions"
  run_stress_smoke
fi

echo "PICurv smoke completed successfully."
