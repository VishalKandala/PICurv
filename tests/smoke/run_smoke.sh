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
LAST_POST_LOG=""
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

require_files_identical() {
  local left="$1"
  local right="$2"
  local label="$3"
  require_file "${left}" "${label} (left operand)"
  require_file "${right}" "${label} (right operand)"
  if ! cmp -s "${left}" "${right}"; then
    die "expected '${left}' and '${right}' to be byte-identical (${label})."
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

create_mock_committed_checkpoint_series() {
  local output_root="$1"
  local first_step="$2"
  local last_step="$3"

  python3 - "${output_root}" "${first_step}" "${last_step}" <<'PY'
import hashlib
import pathlib
import sys

output_root = pathlib.Path(sys.argv[1])
first_step = int(sys.argv[2])
last_step = int(sys.argv[3])
fields = ("Ucat", "Ucont", "Ucont_rm1", "P", "Nvert")

for step in range(first_step, last_step + 1):
    bundle = output_root / "checkpoints" / f"step_{step:012d}"
    block_dir = bundle / "eulerian" / "block_0000"
    block_dir.mkdir(parents=True, exist_ok=True)
    payloads = []
    for field in fields:
        relative = pathlib.Path("eulerian") / "block_0000" / f"{field}.dat"
        payload = bundle / relative
        payload.write_bytes(b"smoke-checkpoint-payload\n")
        payloads.append((relative.as_posix(), field, payload.stat().st_size))

    lines = [
        "-checkpoint_format picurv-checkpoint",
        "-checkpoint_version 1",
        f"-checkpoint_step {step}",
        f"-checkpoint_time {step * 0.1:.17g}",
        "-checkpoint_dt 0.1",
        "-checkpoint_reason smoke_fixture",
        "-checkpoint_geometry_sha256 " + "0" * 64,
        "-checkpoint_block_count 1",
        "-checkpoint_particles false",
        "-checkpoint_particle_count 0",
        "-checkpoint_les false",
        "-checkpoint_rans false",
        f"-checkpoint_payload_count {len(payloads)}",
        "-checkpoint_block_0_im 4",
        "-checkpoint_block_0_jm 4",
        "-checkpoint_block_0_km 4",
        "-checkpoint_periodic 0,0,0",
    ]
    for index, (relative, field, byte_count) in enumerate(payloads):
        lines.extend([
            f"-checkpoint_payload_{index}_path {relative}",
            f"-checkpoint_payload_{index}_kind eulerian",
            f"-checkpoint_payload_{index}_field {field}",
            f"-checkpoint_payload_{index}_block 0",
            f"-checkpoint_payload_{index}_layout smoke",
            f"-checkpoint_payload_{index}_components 1",
            f"-checkpoint_payload_{index}_logical_type PetscScalar",
            f"-checkpoint_payload_{index}_global_size 1",
            f"-checkpoint_payload_{index}_encoding smoke",
            f"-checkpoint_payload_{index}_bytes {byte_count}",
        ])
    metadata = bundle / "checkpoint.meta"
    metadata.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (bundle / "COMMITTED").write_text(
        hashlib.sha256(metadata.read_bytes()).hexdigest() + "\n",
        encoding="ascii",
    )
PY
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


# A workspace may ship more than one bundle. `quickstart_` is the documented
# short-run variant of the same case, and it is the bundle the getting-started page
# tells a new user to run, so it is smoke-tested rather than excluded.
BUNDLE_PREFIXES = ("", "quickstart_")

# Roles an initialized workspace is guaranteed to place at a canonical path. Before
# the workspace layout existed this had to be inferred from a flat pile of YAML -
# "exactly one file that is not the case, monitor, post, study, or cluster profile".
# The layout makes that guessing unnecessary, and a template whose roles did not land
# where init promises is exactly the failure this phase should report.
CANONICAL_ROLES = ("case", "solver", "monitor", "post")


def config_dir(case_dir: str) -> str:
    """The workspace configuration home init writes canonical roles into."""
    return os.path.join(case_dir, "config")


def bundle_prefixes(case_dir: str) -> list:
    """Which bundle prefixes this initialized workspace actually ships."""
    present = [os.path.basename(path)
               for path in glob.glob(os.path.join(config_dir(case_dir), "*.yml"))]
    return [p for p in BUNDLE_PREFIXES
            if p == "" or any(name.startswith(p) for name in present)]


def discover_case_bundle(case_dir: str, template_name: str, prefix: str = ""):
    """Resolve one bundle's four role files from the canonical workspace layout.

    A prefixed bundle overrides whichever roles it ships and inherits the rest from
    the canonical set, which is a legitimate way to author a short-run variant: the
    quickstart bundles reuse the default solver profile rather than duplicating it.
    """
    configs = config_dir(case_dir)
    label = f"{template_name}[{prefix}]" if prefix else template_name
    if not os.path.isdir(configs):
        fail(f"{label}: init produced no config/ directory in {case_dir}")

    present = {os.path.basename(path)
               for path in glob.glob(os.path.join(configs, "*.yml"))}
    resolved = []
    for role in CANONICAL_ROLES:
        canonical = f"{role}.yml"
        # The quickstart bundles keep the example's own names for the roles they
        # override, so match on the prefix rather than on a canonical name.
        override = sorted(
            name for name in present
            if prefix and name.startswith(prefix) and _role_of(configs, name) == role
        )
        chosen = override[0] if override else canonical
        if chosen not in present:
            fail(f"{label}: expected {role} profile '{chosen}' in {configs}, "
                 f"found {sorted(present)}")
        resolved.append(os.path.join("config", chosen))
    return tuple(resolved)


def _role_of(configs: str, name: str) -> str:
    """Classify one prefixed profile by the top-level keys it carries."""
    with open(os.path.join(configs, name), "r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    if not isinstance(payload, dict):
        return ""
    if "grid" in payload or "run_control" in payload:
        return "case"
    if "strategy" in payload or "operation_mode" in payload:
        return "solver"
    if "io" in payload or "logging" in payload:
        return "monitor"
    if "eulerian_pipeline" in payload or "lagrangian_pipeline" in payload:
        return "post"
    return ""


checked_bundles = []

for template_name in templates:
    case_dir = os.path.join(matrix_root, template_name)
    init_cmd = [picurv_exe, "init", template_name, "--dest", case_dir]
    init_res = subprocess.run(init_cmd, text=True, capture_output=True)
    if init_res.returncode != 0:
        fail(f"{template_name}: init failed:\n{init_res.stdout}\n{init_res.stderr}")

    for prefix in bundle_prefixes(case_dir):
        label = f"{template_name}[{prefix}]" if prefix else template_name
        checked_bundles.append(label)

        case_file, solver_file, monitor_file, post_file = discover_case_bundle(case_dir, template_name, prefix)

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
            fail(f"{label}: validate failed:\n{validate_res.stdout}\n{validate_res.stderr}")

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
            fail(f"{label}: dry-run failed:\n{dry_res.stdout}\n{dry_res.stderr}")

        payload = json.loads(dry_res.stdout)
        stages = payload.get("stages", {})
        if "solve" not in stages or "post-process" not in stages:
            fail(f"{label}: dry-run payload missing solve/post-process stages: {payload}")
        solve_launch = stages["solve"].get("launch_command", [])
        post_launch = stages["post-process"].get("launch_command", [])
        if not solve_launch or not any(str(token).endswith("simulator") for token in solve_launch):
            fail(f"{label}: solve launch command does not target simulator: {solve_launch}")
        if not post_launch or not any(str(token).endswith("postprocessor") for token in post_launch):
            fail(f"{label}: post launch command does not target postprocessor: {post_launch}")

print("template matrix smoke validated: " + ", ".join(checked_bundles))
PY
}

run_dry_run_plan_smoke() {
  local plan_json="${tmp_root}/dry-plan.json"
  local monitor_cfg="${tmp_root}/dry-run-monitor-diagnostics.yml"
  cp "${valid_fixtures_dir}/monitor.yml" "${monitor_cfg}"
  python3 - "${monitor_cfg}" <<'PY'
import sys
import yaml
monitor_path = sys.argv[1]
with open(monitor_path, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
cfg["diagnostics"] = {
    "petsc": {
        "malloc_debug": True,
        "log_view": True,
        "malloc_view": True,
    },
    "runtime_memory_log": {
        "enabled": True,
        "file": "Runtime_Memory.log",
    },
}
with open(monitor_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
PY
  "${picurv_exe}" run \
    --solve \
    --post-process \
    --case "${valid_fixtures_dir}/case.yml" \
    --solver "${valid_fixtures_dir}/solver.yml" \
    --monitor "${monitor_cfg}" \
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
solve_launch = [str(token) for token in stages["solve"]["launch_command"]]
assert "-malloc_debug" in solve_launch
assert "-log_view" in solve_launch
artifacts = [str(item) for item in payload.get("artifacts", [])]
assert any(item.endswith("Runtime_Memory.log") for item in artifacts)
assert any(item.endswith("PETSc_LogView_Solver.log") for item in artifacts)
PY
}

run_petsc_diagnostics_smoke() {
  local diag_case="${tmp_root}/petsc-diagnostics"
  local output_log="${tmp_root}/petsc-diagnostics.log"
  local before_runs="${tmp_root}/petsc-diagnostics.runs.before"
  local after_runs="${tmp_root}/petsc-diagnostics.runs.after"
  local created_run=""
  local solver_log=""
  local malloc_view_log=""

  "${picurv_exe}" init flat_channel --dest "${diag_case}" >/dev/null
  prepare_flat_case_les "${diag_case}"
  python3 - \
    "${diag_case}/config/case.yml" \
    "${diag_case}/config/monitor.yml" <<'PY'
import sys
import yaml

case_path, monitor_path = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
with open(monitor_path, "r", encoding="utf-8") as f:
    monitor_cfg = yaml.safe_load(f)

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["total_steps"] = 1
monitor_cfg["diagnostics"] = {
    "petsc": {
        "malloc_debug": True,
        "malloc_dump": True,
        "malloc_view": True,
        "log_view": True,
        "log_view_memory": True,
    },
    "runtime_memory_log": {
        "enabled": True,
        "file": "Runtime_Memory.log",
    },
}

with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(monitor_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(monitor_cfg, f, sort_keys=False)
PY

  mkdir -p "${diag_case}/runs"
  find "${diag_case}/runs" -mindepth 1 -maxdepth 1 -type d | sort >"${before_runs}"
  (
    cd "${diag_case}"
    "${picurv_exe}" run \
      --solve \
      -n 1 \
      --case "${diag_case}/config/case.yml" \
      --solver "${diag_case}/config/solver.yml" \
      --monitor "${diag_case}/config/monitor.yml" >"${output_log}" 2>&1
  ) || {
    echo "Smoke failure: PETSc diagnostics solve failed." >&2
    sed -n '1,220p' "${output_log}" >&2
    exit 1
  }

  find "${diag_case}/runs" -mindepth 1 -maxdepth 1 -type d | sort >"${after_runs}"
  created_run="$(comm -13 "${before_runs}" "${after_runs}" | tail -n 1)"
  if [[ -z "${created_run}" ]]; then
    created_run="$(tail -n 1 "${after_runs}")"
  fi
  if [[ -z "${created_run}" ]]; then
    die "PETSc diagnostics smoke completed but no run directory was found."
  fi

  solver_log="$(find "${created_run}/scheduler" -maxdepth 1 -type f -name '*_solver.log' | sort | tail -n 1)"
  if [[ -z "${solver_log}" ]]; then
    die "PETSc diagnostics smoke completed but no solver runtime log was found."
  fi

  require_file "${created_run}/logs/Runtime_Memory.log" "diagnostics runtime memory log"
  malloc_view_log="$(find "${created_run}/logs" -maxdepth 1 -type f -name 'PETSc_MallocView_Solver.log*' | sort | head -n 1)"
  if [[ -z "${malloc_view_log}" ]]; then
    die "missing PETSc malloc view log under '${created_run}/logs'."
  fi
  require_file "${created_run}/logs/PETSc_LogView_Solver.log" "PETSc log view log"
  require_file_contains "${malloc_view_log}" "Memory usage sorted by function" "PETSc malloc view summary"
  require_file_contains "${created_run}/logs/PETSc_LogView_Solver.log" "Event Stage" "PETSc log view event table"
  require_file_contains "${created_run}/logs/Runtime_Memory.log" "Process Current MB Max" "runtime memory log header"
  require_file_not_contains "${solver_log}" "Memory corruption" "PETSc malloc debug corruption report"
  require_file_not_contains "${solver_log}" "Bad location or corrupted memory" "PETSc malloc debug bad-location report"
  if grep -Eq '^\[[[:space:]]*[0-9]+\][[:space:]]+[0-9]+ bytes' "${solver_log}"; then
    echo "Smoke failure: PETSc malloc_dump reported unfreed allocations in '${solver_log}'." >&2
    grep -En '^\[[[:space:]]*[0-9]+\][[:space:]]+[0-9]+ bytes' "${solver_log}" | head -n 20 >&2
    exit 1
  fi
}

run_restart_resolution_smoke() {
  local restart_root="${tmp_root}/restart-smoke"
  local prior_run="${restart_root}/old_run"
  local prior_output_dir="${prior_run}/output"
  local case_cfg="${restart_root}/case.yml"
  local solver_cfg="${restart_root}/solver.yml"
  local monitor_cfg="${restart_root}/monitor.yml"
  local plan_json="${tmp_root}/restart-plan.json"

  mkdir -p "${restart_root}" "${prior_run}/config" "${prior_output_dir}"
  cp "${valid_fixtures_dir}/case.yml" "${case_cfg}"
  cp "${valid_fixtures_dir}/solver.yml" "${solver_cfg}"
  cp "${valid_fixtures_dir}/monitor.yml" "${monitor_cfg}"
  cat >"${prior_run}/config/monitor.yml" <<'YAML'
io:
  data_output_frequency: 1
YAML

  # Create complete fake bundles for load mode validation (steps 5-15).
  create_mock_committed_checkpoint_series "${prior_output_dir}" 5 15

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
      --format json >"${plan_json}" 2>"${plan_json}.err" || {
        echo "restart resolution dry-run failed:" >&2
        cat "${plan_json}.err" >&2
        exit 1
      }
  )

  python3 - "${plan_json}" "${restart_root}" <<'PY'
import json
import os
import sys
plan_path, restart_root = sys.argv[1], os.path.abspath(sys.argv[2])
with open(plan_path, "r", encoding="utf-8") as f:
    payload = json.load(f)
solve_stage = payload.get("stages", {}).get("solve", {})
resolved = solve_stage.get("restart_source_directory")
# Load mode materializes the consumed sequence into the new run rather than reading
# the source run in place, so the plan must name this run's own restart input.
assert resolved and resolved.startswith(os.path.join(restart_root, "runs") + os.sep), \
    f"restart source is not inside this workspace's runs tree: {resolved}"
assert resolved.endswith(os.path.join("inputs", "restart")), \
    f"restart source is not the canonical restart input: {resolved}"
# A dry run promises to write nothing, so the directory it names must not exist yet.
assert not os.path.exists(resolved), f"dry run materialized {resolved}"
PY
}

prepare_flat_case_les() {
  local case_dir="$1"
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
    "${case_dir}/config/case.yml" \
    "${case_dir}/config/solver.yml" \
    "${case_dir}/config/monitor.yml" \
    "${case_dir}/config/post.yml" \
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

prepare_flat_case_field_statistics() {
  local case_dir="$1"
  local statistics_enabled="$2"
  local console_frequency="$3"
  local verbosity="$4"
  local post_output_dir="$5"
  python3 - \
    "${case_dir}/config/case.yml" \
    "${case_dir}/config/solver.yml" \
    "${case_dir}/config/monitor.yml" \
    "${case_dir}/config/post.yml" \
    "${statistics_enabled}" \
    "${console_frequency}" \
    "${verbosity}" \
    "${post_output_dir}" <<'PY'
import sys
import yaml

(
    case_path,
    solver_path,
    monitor_path,
    post_path,
    statistics_enabled,
    console_frequency,
    verbosity,
    post_output_dir,
) = sys.argv[1:]

statistics_enabled = statistics_enabled == "1"

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
case_cfg["models"]["physics"]["turbulence"] = {"les": False}

# An analytically prescribed field is evaluated pointwise, so the same grid point
# carries the same value whatever rank owns it. That is what makes the accumulator
# payloads comparable byte for byte across a decomposition change; a solved field
# would differ in the last bits through the Krylov reductions alone.
solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "analytical"
solver_cfg["operation_mode"]["analytical_type"] = "UNIFORM_FLOW"
solver_cfg["operation_mode"]["uniform_flow"] = {"u": 0.05, "v": -0.02, "w": 0.1}

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 2
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1
monitor_cfg["io"]["statistics_console_output_frequency"] = int(console_frequency)
monitor_cfg.setdefault("logging", {})
monitor_cfg["logging"]["verbosity"] = verbosity

if statistics_enabled:
    monitor_cfg["field_statistics"] = {
        "enabled": True,
        "windows": [
            {
                "name": "equivalence",
                "start_time": 0.0,
                "weighting": "sample",
                "step_cadence": 1,
                "fields": [
                    {"field": "Ucat", "moments": ["first", "second"]},
                    {"field": "P", "moments": ["first"]},
                ],
                "covariances": [["Ucat", "P"]],
            }
        ],
    }
else:
    monitor_cfg.pop("field_statistics", None)

post_cfg.setdefault("run_control", {})
# Two processed steps, each reading the accumulator state committed at that step, so
# the CSV carries a real convergence history rather than one row repeated.
post_cfg["run_control"]["start_step"] = 2
post_cfg["run_control"]["end_step"] = 4
post_cfg["run_control"]["step_interval"] = 2
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = post_output_dir
post_cfg["io"]["output_filename_prefix"] = "Field"
post_cfg["io"]["output_particles"] = False

if statistics_enabled:
    post_cfg["field_statistics"] = {
        "windows": ["equivalence"],
        "outputs": ["mean", "reynolds_stress", "rms", "tke", "flux"],
        "formats": ["vtk", "csv"],
    }
else:
    post_cfg.pop("field_statistics", None)

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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
prepare_bent_case_analytical_uniform_tiny() {
  local case_dir="$1"
  # bent_channel.yml's own grid.mode is grid_gen, which validate_case refuses to combine
  # with a non-TGV3D analytical carrier (UNIFORM_FLOW here) - only 'file' or
  # 'programmatic_c' are allowed there. This variant exists specifically to exercise the
  # analytical branch on the bent (curvilinear) domain, so it overrides to 'file' and
  # points at the one bent-channel grid the examples still ship, rather than dropping to
  # a programmatic Cartesian box that would no longer exercise curvilinear metrics. The
  # workspace layer refuses an absolute source_file outside the case (it wants an
  # explicit import), so the grid is copied in and referenced relatively instead.
  # source_file resolves relative to the case root, not to config/case.yml itself.
  cp "${repo_root}/examples/search_robustness/bent_channel_coarse.picgrid" \
     "${case_dir}/bent_channel_coarse.picgrid"
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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

case_cfg["grid"] = {"mode": "file", "source_file": "bent_channel_coarse.picgrid"}

case_cfg.setdefault("run_control", {})
case_cfg["run_control"]["start_step"] = 0
case_cfg["run_control"]["total_steps"] = 3
case_cfg["run_control"]["dt_physical"] = 0.001
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}

solver_cfg.setdefault("operation_mode", {})
solver_cfg["operation_mode"]["eulerian_field_source"] = "analytical"
solver_cfg["operation_mode"]["analytical_type"] = "UNIFORM_FLOW"
solver_cfg["operation_mode"]["uniform_flow"] = {"u": 0.05, "v": -0.02, "w": 0.1}

monitor_cfg.setdefault("io", {})
monitor_cfg["io"]["data_output_frequency"] = 1
monitor_cfg["io"]["particle_console_output_frequency"] = 0
monitor_cfg["io"]["particle_log_interval"] = 1

post_cfg.setdefault("run_control", {})
post_cfg["run_control"]["start_step"] = 0
post_cfg["run_control"]["end_step"] = 3
post_cfg["run_control"]["step_interval"] = 1
post_cfg.setdefault("io", {})
post_cfg["io"]["output_directory"] = "viz/bent_analytical_smoke"
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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
  python3 - "${case_dir}/config/case.yml" "${case_dir}/config/solver.yml" "${case_dir}/config/monitor.yml" "${case_dir}/config/post.yml" <<'PY'
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
case_cfg.setdefault("models", {}).setdefault("physics", {})
case_cfg["models"]["physics"]["particles"] = {"count": 0}
case_cfg["boundary_conditions"] = [
    {"face": "-Xi", "type": "WALL", "handler": "noslip"},
    {"face": "+Xi", "type": "WALL", "handler": "noslip"},
    {"face": "-Eta", "type": "WALL", "handler": "noslip"},
    {"face": "+Eta", "type": "WALL", "handler": "noslip"},
    {"face": "-Zeta", "type": "PERIODIC", "handler": "constant_flux", "params": {"target_flux": 1.0, "enforce_seam_flux": True}},
    {"face": "+Zeta", "type": "PERIODIC", "handler": "constant_flux", "params": {"target_flux": 1.0, "enforce_seam_flux": True}},
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
  LAST_POST_LOG="$(
    python3 - "${created_run}" <<'PY'
import json
import os
import sys

run_dir = sys.argv[1]
with open(os.path.join(run_dir, "scheduler", "submission.json"), encoding="utf-8") as stream:
    submission = json.load(stream)
post_stage = submission.get("stages", {}).get("post-process")
if post_stage:
    print(os.path.join(run_dir, post_stage["log_file"]))
PY
  )"
  if [[ -n "${LAST_POST_LOG}" ]]; then
    require_file "${LAST_POST_LOG}" "postprocessor runtime log"
    require_file_contains "${LAST_POST_LOG}" "Postprocessor MPI processes: ${nprocs}" \
      "postprocessor MPI rank count in runtime summary"
  fi
  require_file "${created_run}/logs/Runtime_Memory.log" "runtime memory log"
  require_file_contains "${created_run}/logs/Runtime_Memory.log" "Process Current MB Max" "runtime memory log header"
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
    "${continuous_case}/config/case.yml" \
    "${continuous_case}/config/solver.yml" \
    "${continuous_case}/config/monitor.yml" \
    "${continuous_case}/config/post.yml" \
    "restart_equiv_continuous"

  continuous_run="${LAST_RUN_DIR}"
  continuous_continuity_log="${continuous_run}/logs/Continuity_Metrics.log"
  require_file "${continuous_continuity_log}" "continuous continuity metrics log"

  "${picurv_exe}" init flat_channel --dest "${split_case}" >/dev/null
  prepare_flat_restart_equivalence_case "${split_case}" 0 3 "solve" "" "viz/restart_equivalence_split_base" 0 3
  run_case_workflow \
    "${split_case}" \
    "${split_case}/config/case.yml" \
    "${split_case}/config/solver.yml" \
    "${split_case}/config/monitor.yml" \
    "${split_case}/config/post.yml" \
    "restart_equiv_split_base"

  split_base_run="${LAST_RUN_DIR}"
  prepare_flat_restart_equivalence_case "${split_case}" 2 2 "solve" "${split_base_run}" "viz/restart_equivalence_split_restart" 2 4
  run_case_workflow \
    "${split_case}" \
    "${split_case}/config/case.yml" \
    "${split_case}/config/solver.yml" \
    "${split_case}/config/monitor.yml" \
    "${split_case}/config/post.yml" \
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

# Exercise the complete CLI restart path for both Newton--Krylov PC models:
# output checkpoint -> --continue staging -> C field load/history seed -> SNES.
run_newton_krylov_continue_smoke() {
  local profile_name profile_source case_dir base_run continue_log summary_file saved_nprocs="${nprocs}"

  # The NK runtime fixture must use its supported all-geometric-periodic scope.
  nprocs=1

  for profile_name in pcnone frozen_point_block; do
    if [[ "${profile_name}" == "pcnone" ]]; then
      profile_source="${repo_root}/config/solvers/Newton-Krylov-Standard.yml"
    else
      profile_source="${repo_root}/config/solvers/Newton-Krylov-Frozen-Momentum-Point-Block.yml"
    fi

    case_dir="${tmp_root}/newton-krylov-restart-${profile_name}"
    cp -R "${repo_root}/examples/periodic_test/constant_uniform_flow" "${case_dir}"
    cp "${profile_source}" "${case_dir}/solver.yml"
    python3 - "${case_dir}/case.yml" "${case_dir}/solver.yml" <<'PY'
import sys
import yaml
case_path, solver_path = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
case_cfg["run_control"]["start_step"] = 0
case_cfg["run_control"]["total_steps"] = 2
with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
# The 16-node periodic smoke grid supports the two-level hierarchy used by
# the periodic fixture, not the four-level production default.
mg = solver_cfg["poisson_solver"]["multigrid"]
mg["levels"] = 2
for level in ("level_2", "level_3"):
    mg.get("level_solvers", {}).pop(level, None)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
PY
    run_case_workflow \
      "${case_dir}" \
      "${case_dir}/case.yml" \
      "${case_dir}/solver.yml" \
      "${case_dir}/monitor.yml" \
      "${case_dir}/post.yml" \
      "newton_krylov_${profile_name}_base"
    base_run="${LAST_RUN_DIR}"

    python3 - "${case_dir}/case.yml" <<'PY'
import sys
import yaml
with open(sys.argv[1], "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
case_cfg["run_control"]["start_step"] = 2
case_cfg["run_control"]["total_steps"] = 1
with open(sys.argv[1], "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
PY
    continue_log="${tmp_root}/newton_krylov_${profile_name}_continue.log"
    (
      cd "${case_dir}"
      "${picurv_exe}" run --solve --continue --run-dir "${base_run}" -n "${nprocs}" \
        --case "${case_dir}/case.yml" \
        --solver "${case_dir}/solver.yml" \
        --monitor "${case_dir}/monitor.yml" >"${continue_log}" 2>&1
    ) || {
      echo "Smoke failure: Newton--Krylov ${profile_name} --continue workflow failed." >&2
      sed -n '1,220p' "${continue_log}" >&2
      exit 1
    }

    summary_file="${base_run}/logs/Momentum_Solver_Newton_Krylov_Summary_Block_0.log"
    require_file "${summary_file}" "Newton--Krylov ${profile_name} summary"
    require_file_contains "${continue_log}" "Continue mode: logs will be appended" "Newton--Krylov ${profile_name} CLI continuation"
    require_file_contains "${continue_log}" "Restart source:" "Newton--Krylov ${profile_name} restart staging"
    require_file_contains "${summary_file}" "step: 3 | block: 0 | solver: Newton Krylov" "Newton--Krylov ${profile_name} continued solve"
    if [[ "${profile_name}" == "pcnone" ]]; then
      require_file_contains "${summary_file}" "Preconditioner: none" "Newton--Krylov PCNONE summary"
    else
      require_file_contains "${summary_file}" "Preconditioner: frozen_momentum_jacobian / point_block" "Newton--Krylov frozen point-block summary"
    fi
  done
  nprocs="${saved_nprocs}"
}

# Exercise the physical-wall/inlet/conservation-outlet startup path separately
# from the all-periodic restart fixture.  In particular, this guards the first
# BDF1 residual before any restart state exists.
run_newton_krylov_flat_channel_startup_smoke() {
  local profile_name profile_source case_dir summary_file

  for profile_name in pcnone frozen_point_block; do
    if [[ "${profile_name}" == "pcnone" ]]; then
      profile_source="${repo_root}/config/solvers/Newton-Krylov-Standard.yml"
    else
      profile_source="${repo_root}/config/solvers/Newton-Krylov-Frozen-Momentum-Point-Block.yml"
    fi
    case_dir="${tmp_root}/newton-krylov-flat-channel-${profile_name}"
    "${picurv_exe}" init flat_channel --dest "${case_dir}" >/dev/null
    cp "${profile_source}" "${case_dir}/solver.yml"
    python3 - "${case_dir}/config/case.yml" "${case_dir}/solver.yml" <<'PY'
import sys
import yaml

case_path, solver_path = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
case_cfg["run_control"]["start_step"] = 0
case_cfg["run_control"]["total_steps"] = 1
grid = case_cfg["grid"]["programmatic_settings"]
grid["im"], grid["jm"], grid["km"] = 8, 8, 16
with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)

with open(solver_path, "r", encoding="utf-8") as f:
    solver_cfg = yaml.safe_load(f)
mg = solver_cfg["poisson_solver"]["multigrid"]
mg["levels"] = 2
for level in ("level_2", "level_3"):
    mg.get("level_solvers", {}).pop(level, None)
with open(solver_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(solver_cfg, f, sort_keys=False)
PY
    run_case_workflow \
      "${case_dir}" \
      "${case_dir}/config/case.yml" \
      "${case_dir}/solver.yml" \
      "${case_dir}/config/monitor.yml" \
      "${case_dir}/config/post.yml" \
      "newton_krylov_flat_channel_${profile_name}"
    summary_file="${LAST_RUN_DIR}/logs/Momentum_Solver_Newton_Krylov_Summary_Block_0.log"
    require_file "${summary_file}" "Newton--Krylov ${profile_name} flat-channel summary"
    python3 - "${summary_file}" "${profile_name}" <<'PY'
import math
import re
import sys

path, profile = sys.argv[1:]
rows = [line.strip() for line in open(path, encoding="utf-8") if line.startswith("step:")]
assert len(rows) == 1, f"{profile}: expected one BDF1 Newton summary row, got {len(rows)}"
row = rows[0]
match = re.search(r"\| initial: ([^|]+) \| final:", row)
assert match, f"{profile}: missing initial residual in {row}"
initial = float(match.group(1))
assert math.isfinite(initial), f"{profile}: non-finite initial residual {initial}"
assert "state: committed" in row, f"{profile}: first BDF1 Newton solve did not commit"
PY
  done
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
    "${flat_les_case}/config/case.yml" \
    "${flat_les_case}/config/solver.yml" \
    "${flat_les_case}/config/monitor.yml" \
    "${flat_les_case}/config/post.yml" \
    "flat_les"

  local flat_les_run
  flat_les_run="${LAST_RUN_DIR}"
  require_dir "${flat_les_run}/output/checkpoints" "flat LES checkpoint directory"
  require_count_ge "${flat_les_run}/output/checkpoints" "COMMITTED" 1 "flat LES committed checkpoints"
  require_count_ge "${flat_les_run}/output/checkpoints" "Ucat.dat" 1 "flat LES Eulerian checkpoint payloads"
  require_count_ge "${flat_les_run}/output/visualization" "*.vts" 1 "flat LES post VTS files"
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
    "${bent_case}/config/case.yml" \
    "${bent_case}/config/solver.yml" \
    "${bent_case}/config/monitor.yml" \
    "${bent_case}/config/post.yml" \
    "bent_les"

  local bent_run
  bent_run="${LAST_RUN_DIR}"
  require_dir "${bent_run}/output/checkpoints" "bent checkpoint directory"
  require_count_ge "${bent_run}/output/checkpoints" "COMMITTED" 1 "bent committed checkpoints"
  require_count_ge "${bent_run}/output/checkpoints" "Ucat.dat" 1 "bent Eulerian checkpoint payloads"
  require_count_ge "${bent_run}/output/visualization" "*.vts" 1 "bent post VTS files"

  local bent_analytical_case="${tmp_root}/bent-analytical"
  "${picurv_exe}" init bent_channel --dest "${bent_analytical_case}" >/dev/null
  prepare_bent_case_analytical_uniform_tiny "${bent_analytical_case}"
  run_case_workflow \
    "${bent_analytical_case}" \
    "${bent_analytical_case}/config/case.yml" \
    "${bent_analytical_case}/config/solver.yml" \
    "${bent_analytical_case}/config/monitor.yml" \
    "${bent_analytical_case}/config/post.yml" \
    "bent_analytical_uniform"

  local bent_analytical_run
  bent_analytical_run="${LAST_RUN_DIR}"
  require_dir "${bent_analytical_run}/output/checkpoints" "bent analytical checkpoint directory"
  require_count_ge "${bent_analytical_run}/output/checkpoints" "COMMITTED" 1 "bent analytical committed checkpoints"
  require_count_ge "${bent_analytical_run}/output/checkpoints" "Ucat.dat" 1 "bent analytical Eulerian checkpoint payloads"
  require_count_ge "${bent_analytical_run}/output/visualization" "*.vts" 1 "bent analytical post VTS files"
  require_file_contains "${LAST_SOLVER_LOG}" "Analytical Solution Type : UNIFORM_FLOW" "bent analytical uniform-flow runtime branch"
  "${picurv_exe}" init flat_channel --dest "${flat_particles_case}" >/dev/null
  prepare_flat_case_particles_base "${flat_particles_case}"
  run_case_workflow \
    "${flat_particles_case}" \
    "${flat_particles_case}/config/case.yml" \
    "${flat_particles_case}/config/solver.yml" \
    "${flat_particles_case}/config/monitor.yml" \
    "${flat_particles_case}/config/post.yml" \
    "flat_particles_base"

  local base_particles_run
  base_particles_run="${LAST_RUN_DIR}"
  require_count_ge "${base_particles_run}/output/checkpoints" "position.dat" 1 "particle checkpoint payloads"
  require_count_ge "${base_particles_run}/output/visualization" "*.vtp" 1 "particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 32" "particle runtime banner particle count"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence   : DISABLED" "particle runtime banner disabled particle console cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling  : every 1 particle(s)" "particle runtime banner row sampling"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Initialization Mode: Point Source" "particle runtime banner initialization mode"
  require_file_contains "${LAST_SOLVER_LOG}" "Interpolation Method       : Trilinear (direct cell-center)" "particle runtime banner default interpolation method"

  "${picurv_exe}" init flat_channel --dest "${flat_particles_ca_case}" >/dev/null
  prepare_flat_case_particles_corner_averaged "${flat_particles_ca_case}"
  run_case_workflow \
    "${flat_particles_ca_case}" \
    "${flat_particles_ca_case}/config/case.yml" \
    "${flat_particles_ca_case}/config/solver.yml" \
    "${flat_particles_ca_case}/config/monitor.yml" \
    "${flat_particles_ca_case}/config/post.yml" \
    "flat_particles_corner_averaged"

  require_count_ge "${LAST_RUN_DIR}/output/checkpoints" "position.dat" 1 "corner-averaged particle checkpoint payloads"
  require_count_ge "${LAST_RUN_DIR}/output/visualization" "*.vtp" 1 "corner-averaged particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Interpolation Method       : CornerAveraged (legacy)" "corner-averaged runtime banner interpolation method"

  local case_restart_load="${flat_particles_case}/case_restart_load.yml"
  local solver_restart_load="${flat_particles_case}/solver_restart_load.yml"
  local post_restart_load="${flat_particles_case}/post_restart_load.yml"
  cp "${flat_particles_case}/config/case.yml" "${case_restart_load}"
  cp "${flat_particles_case}/config/solver.yml" "${solver_restart_load}"
  cp "${flat_particles_case}/config/post.yml" "${post_restart_load}"
  prepare_flat_restart_variant "${case_restart_load}" "${solver_restart_load}" "${post_restart_load}" "${base_particles_run}" "load" "viz/restart_load"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_load}" \
    "${solver_restart_load}" \
    "${flat_particles_case}/config/monitor.yml" \
    "${post_restart_load}" \
    "flat_particles_restart_load" \
    "${base_particles_run}"

  local restart_load_run
  restart_load_run="${LAST_RUN_DIR}"
  require_count_ge "${restart_load_run}/output/checkpoints" "position.dat" 1 "restart-load particle checkpoint payloads"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : load" "restart load branch"

  local case_restart_init="${flat_particles_case}/case_restart_init.yml"
  local solver_restart_init="${flat_particles_case}/solver_restart_init.yml"
  local post_restart_init="${flat_particles_case}/post_restart_init.yml"
  cp "${flat_particles_case}/config/case.yml" "${case_restart_init}"
  cp "${flat_particles_case}/config/solver.yml" "${solver_restart_init}"
  cp "${flat_particles_case}/config/post.yml" "${post_restart_init}"
  prepare_flat_restart_variant "${case_restart_init}" "${solver_restart_init}" "${post_restart_init}" "${base_particles_run}" "init" "viz/restart_init"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_init}" \
    "${solver_restart_init}" \
    "${flat_particles_case}/config/monitor.yml" \
    "${post_restart_init}" \
    "flat_particles_restart_init" \
    "${base_particles_run}"

  local restart_init_run
  restart_init_run="${LAST_RUN_DIR}"
  require_count_ge "${restart_init_run}/output/checkpoints" "position.dat" 1 "restart-init particle checkpoint payloads"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : init" "restart init branch"

  run_restart_equivalence_smoke
  run_newton_krylov_continue_smoke

  "${picurv_exe}" init brownian_motion --dest "${brownian_case}" >/dev/null
  prepare_brownian_case_analytical "${brownian_case}"
  run_case_workflow \
    "${brownian_case}" \
    "${brownian_case}/config/case.yml" \
    "${brownian_case}/config/solver.yml" \
    "${brownian_case}/config/monitor.yml" \
    "${brownian_case}/config/post.yml" \
    "brownian_analytical"

  local brownian_run
  brownian_run="${LAST_RUN_DIR}"
  require_count_ge "${brownian_run}/output/visualization" "*.vts" 1 "brownian eulerian VTS files"
  require_count_ge "${brownian_run}/output/visualization" "*.vtp" 1 "brownian particle VTP files"
  require_count_ge "${brownian_run}/output/analysis/statistics" "BrownianStats_msd.csv" 1 \
    "brownian MSD statistics CSV"
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
    "${flat_case}/config/case.yml" \
    "${flat_case}/config/solver.yml" \
    "${flat_case}/config/monitor.yml" \
    "${flat_case}/config/post.yml" \
    "flat_mpi"

  local flat_run
  flat_run="${LAST_RUN_DIR}"
  require_count_ge "${flat_run}/output/checkpoints" "Ucat.dat" 1 "flat MPI Eulerian checkpoint payloads"
  require_count_ge "${flat_run}/output/visualization" "*.vts" 1 "flat MPI post VTS files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "flat MPI rank count in runtime summary"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 0" "flat MPI runtime banner particle count"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence" "flat MPI runtime banner particle console cadence omission"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling" "flat MPI runtime banner particle row sampling omission"

  "${picurv_exe}" init bent_channel --dest "${bent_case}" >/dev/null
  prepare_bent_case_tiny "${bent_case}"
  run_case_workflow \
    "${bent_case}" \
    "${bent_case}/config/case.yml" \
    "${bent_case}/config/solver.yml" \
    "${bent_case}/config/monitor.yml" \
    "${bent_case}/config/post.yml" \
    "bent_mpi"

  local bent_run
  bent_run="${LAST_RUN_DIR}"
  require_count_ge "${bent_run}/output/checkpoints" "Ucat.dat" 1 "bent MPI Eulerian checkpoint payloads"
  require_count_ge "${bent_run}/output/visualization" "*.vts" 1 "bent MPI post VTS files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "bent MPI rank count in runtime summary"

  "${picurv_exe}" init flat_channel --dest "${flat_particles_case}" >/dev/null
  prepare_flat_case_particles_base "${flat_particles_case}"
  run_case_workflow \
    "${flat_particles_case}" \
    "${flat_particles_case}/config/case.yml" \
    "${flat_particles_case}/config/solver.yml" \
    "${flat_particles_case}/config/monitor.yml" \
    "${flat_particles_case}/config/post.yml" \
    "flat_particles_mpi_base"

  local base_particles_run
  base_particles_run="${LAST_RUN_DIR}"
  require_count_ge "${base_particles_run}/output/checkpoints" "position.dat" 1 "MPI particle checkpoint payloads"
  require_count_ge "${base_particles_run}/output/visualization" "*.vtp" 1 "MPI particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "MPI particle run rank count in runtime summary"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 32" "MPI particle runtime banner particle count"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Console Cadence   : DISABLED" "MPI particle runtime banner disabled particle console cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Log Row Sampling  : every 1 particle(s)" "MPI particle runtime banner row sampling"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Initialization Mode: Point Source" "MPI particle runtime banner initialization mode"

  local case_restart_load="${flat_particles_case}/case_restart_load.yml"
  local solver_restart_load="${flat_particles_case}/solver_restart_load.yml"
  local post_restart_load="${flat_particles_case}/post_restart_load.yml"
  cp "${flat_particles_case}/config/case.yml" "${case_restart_load}"
  cp "${flat_particles_case}/config/solver.yml" "${solver_restart_load}"
  cp "${flat_particles_case}/config/post.yml" "${post_restart_load}"
  prepare_flat_restart_variant "${case_restart_load}" "${solver_restart_load}" "${post_restart_load}" "${base_particles_run}" "load" "viz/restart_load_mpi"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_load}" \
    "${solver_restart_load}" \
    "${flat_particles_case}/config/monitor.yml" \
    "${post_restart_load}" \
    "flat_particles_mpi_restart_load" \
    "${base_particles_run}"
  require_count_ge "${LAST_RUN_DIR}/output/checkpoints" "position.dat" 1 "MPI restart-load particle checkpoint payloads"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : load" "MPI restart load branch"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "MPI restart-load rank count in runtime summary"

  local case_restart_init="${flat_particles_case}/case_restart_init.yml"
  local solver_restart_init="${flat_particles_case}/solver_restart_init.yml"
  local post_restart_init="${flat_particles_case}/post_restart_init.yml"
  cp "${flat_particles_case}/config/case.yml" "${case_restart_init}"
  cp "${flat_particles_case}/config/solver.yml" "${solver_restart_init}"
  cp "${flat_particles_case}/config/post.yml" "${post_restart_init}"
  prepare_flat_restart_variant "${case_restart_init}" "${solver_restart_init}" "${post_restart_init}" "${base_particles_run}" "init" "viz/restart_init_mpi"
  run_case_workflow \
    "${flat_particles_case}" \
    "${case_restart_init}" \
    "${solver_restart_init}" \
    "${flat_particles_case}/config/monitor.yml" \
    "${post_restart_init}" \
    "flat_particles_mpi_restart_init" \
    "${base_particles_run}"
  require_count_ge "${LAST_RUN_DIR}/output/checkpoints" "position.dat" 1 "MPI restart-init particle checkpoint payloads"
  require_file_contains "${LAST_SOLVER_LOG}" "Particle Restart Mode      : init" "MPI restart init branch"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${nprocs}" "MPI restart-init rank count in runtime summary"
}

# Restart a run under a different MPI rank count than the one that wrote the
# checkpoint. Committed Eulerian payloads are stored in DMDA natural ordering via
# GatherVectorToRankZero, so persisted state must be independent of the
# decomposition that wrote it.
#
# The assertion is split in two because a restart does not simply echo what it
# read: re-establishing a boundary-consistent state at branch_start recomputes
# Ucat and Ucont, which a same-rank restart does too. That recomputation is not a
# decomposition effect, so it is held constant rather than asserted away:
#
#   1. Payloads that are pure load-and-store (P, Nvert, Ucont_rm1) must come back
#      byte identical to the multi-rank base bundle after a single-rank restart.
#      This is the direct natural-ordering round-trip proof.
#   2. Restarting the one identical base bundle at two different rank counts must
#      produce byte-identical branch_start bundles for every payload. Same input
#      bytes, same code path, only the decomposition differs, so any difference
#      here would be a genuine rank dependence.
#
# Only meaningful when the harness rank count differs from one, so this runs in
# the multi-rank tier.
run_rank_change_restart_smoke() {
  local case_dir="${tmp_root}/rank-change-restart"
  local base_run="" restart_one_run="" restart_many_run=""
  local saved_nprocs="${nprocs}"
  local step_dir="output/checkpoints/step_000000000002"
  local block_dir="${step_dir}/eulerian/block_0000"
  local payload=""

  if [[ "${saved_nprocs}" -le 1 ]]; then
    echo "    (skipped: rank-change restart needs a multi-rank harness, have ${saved_nprocs})"
    return 0
  fi

  "${picurv_exe}" init flat_channel --dest "${case_dir}" >/dev/null
  prepare_flat_restart_equivalence_case "${case_dir}" 0 2 "solve" "" "viz/rank_change_base" 0 2
  run_case_workflow \
    "${case_dir}" \
    "${case_dir}/config/case.yml" \
    "${case_dir}/config/solver.yml" \
    "${case_dir}/config/monitor.yml" \
    "${case_dir}/config/post.yml" \
    "rank_change_base"
  base_run="${LAST_RUN_DIR}"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${saved_nprocs}" \
    "rank-change base run rank count"
  require_dir "${base_run}/${step_dir}" "rank-change base committed bundle at the resume step"

  # Restart the multi-rank bundle on a single rank.
  nprocs=1
  prepare_flat_restart_equivalence_case "${case_dir}" 2 2 "solve" "${base_run}" "viz/rank_change_restart_one" 2 4
  run_case_workflow \
    "${case_dir}" \
    "${case_dir}/config/case.yml" \
    "${case_dir}/config/solver.yml" \
    "${case_dir}/config/monitor.yml" \
    "${case_dir}/config/post.yml" \
    "rank_change_restart_one" \
    "${base_run}"
  restart_one_run="${LAST_RUN_DIR}"
  nprocs="${saved_nprocs}"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : 1" \
    "rank-change single-rank restart rank count"
  require_dir "${restart_one_run}/${step_dir}" "rank-change single-rank branch_start bundle"

  # (1) Load-and-store payloads must survive the decomposition change untouched.
  for payload in P.dat Nvert.dat Ucont_rm1.dat; do
    require_files_identical \
      "${base_run}/${block_dir}/${payload}" \
      "${restart_one_run}/${block_dir}/${payload}" \
      "rank-change natural-ordering round trip for ${payload} (${saved_nprocs} ranks -> 1 rank)"
  done

  # Restart the very same bundle again, this time at the original rank count.
  prepare_flat_restart_equivalence_case "${case_dir}" 2 2 "solve" "${base_run}" "viz/rank_change_restart_many" 2 4
  run_case_workflow \
    "${case_dir}" \
    "${case_dir}/config/case.yml" \
    "${case_dir}/config/solver.yml" \
    "${case_dir}/config/monitor.yml" \
    "${case_dir}/config/post.yml" \
    "rank_change_restart_many" \
    "${base_run}"
  restart_many_run="${LAST_RUN_DIR}"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${saved_nprocs}" \
    "rank-change multi-rank restart rank count"
  require_dir "${restart_many_run}/${step_dir}" "rank-change multi-rank branch_start bundle"

  # (2) Identical input bundle, differing only in rank count, must reload identically.
  require_count_ge "${base_run}/${block_dir}" "*.dat" 3 "rank-change base Eulerian payloads to compare"
  while IFS= read -r payload; do
    require_files_identical \
      "${restart_one_run}/${block_dir}/${payload}" \
      "${restart_many_run}/${block_dir}/${payload}" \
      "rank-change restart payload ${payload} (1 rank vs ${saved_nprocs} ranks from one bundle)"
  done < <(find "${base_run}/${block_dir}" -type f -name '*.dat' -printf '%f\n' | sort)

  # The restart must not merely reload; it must continue marching to the new horizon.
  require_dir "${restart_one_run}/output/checkpoints/step_000000000004" \
    "rank-change restart final bundle after resuming on a different rank count"
}

# Each variant gets its own case directory: run identifiers are stamped to the second,
# so two workflows finishing inside the same second would otherwise share a run tree.
run_field_statistics_variant() {
  local label="$1"
  local statistics_enabled="$2"
  local console_frequency="$3"
  local verbosity="$4"
  local case_dir="${tmp_root}/${label}"

  "${picurv_exe}" init flat_channel --dest "${case_dir}" >/dev/null
  prepare_flat_case_field_statistics "${case_dir}" "${statistics_enabled}" "${console_frequency}" \
    "${verbosity}" "viz/${label}"
  run_case_workflow \
    "${case_dir}" \
    "${case_dir}/config/case.yml" \
    "${case_dir}/config/solver.yml" \
    "${case_dir}/config/monitor.yml" \
    "${case_dir}/config/post.yml" \
    "${label}"
}

run_field_statistics_smoke() {
  local reported_run="" quiet_run="" disabled_run=""
  local stats_dir="output/checkpoints/step_000000000004/statistics/window_0000/block_0000"

  # (1) Statistics accumulating, console reporting live.
  run_field_statistics_variant "field_statistics_reported" 1 2 "INFO"
  reported_run="${LAST_RUN_DIR}"
  require_file_contains "${LAST_SOLVER_LOG}" "Statistics Console Cadence : every 2 step(s), 1 window(s)" \
    "startup banner reports the configured statistics console cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "Statistics windows at step" \
    "console snapshot is emitted at the configured cadence"
  require_file_contains "${LAST_SOLVER_LOG}" "equivalence .*samples=" \
    "console snapshot reports the named window's accepted sample count"
  require_dir "${reported_run}/${stats_dir}" "committed accumulator payloads"
  require_count_ge "${reported_run}/${stats_dir}" "*.dat" 5 \
    "accumulator payloads for two fields, a second moment, and a covariance"
  require_count_ge "${reported_run}/output/visualization" \
    "Field_statistics_equivalence_*.vts" 2 "derived statistics VTK output per processed step"
  require_count_ge "${reported_run}/output/analysis/statistics" \
    "Field_statistics_equivalence.csv" 1 "statistics convergence history CSV"
  # The derived VTK field and the convergence CSV come from paths that share no code,
  # so comparing them catches a layout boundary an interior-only producer left unwritten.
  if ! python3 "${repo_root}/tests/tooling/check_statistics_nodal_consistency.py" \
       "${reported_run}" "equivalence"; then
    die "derived statistics VTK output disagrees with the convergence history CSV."
  fi
  require_file_contains \
    "$(find "${reported_run}/output/analysis/statistics" -type f \
        -name Field_statistics_equivalence.csv | head -1)" \
    "step,state,samples" "statistics convergence history header"

  # (2) The same accumulation below LOG_INFO: the banner still records the cadence, so a
  #     quiet log stays interpretable, but no snapshot is emitted.
  run_field_statistics_variant "field_statistics_quiet" 1 2 "WARNING"
  quiet_run="${LAST_RUN_DIR}"
  require_file_contains "${LAST_SOLVER_LOG}" "Statistics Console Cadence : every 2 step(s), 1 window(s)" \
    "startup banner reports the cadence even below LOG_INFO"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Statistics windows at step" \
    "no console snapshot is emitted below LOG_INFO"
  require_dir "${quiet_run}/${stats_dir}" "accumulation continues while the console is quiet"

  # (3) A run that configured no window must be untouched by the subsystem on every path.
  run_field_statistics_variant "field_statistics_disabled" 0 0 "INFO"
  disabled_run="${LAST_RUN_DIR}"
  require_file_contains "${LAST_SOLVER_LOG}" "Statistics Console Cadence : DISABLED (no window configured)" \
    "startup banner records that no statistics window was configured"
  require_file_not_contains "${LAST_SOLVER_LOG}" "Statistics windows at step" \
    "a disabled run emits no statistics console output"
  if [[ -d "${disabled_run}/output/checkpoints/step_000000000004/statistics" ]]; then
    die "a disabled run wrote a statistics subtree into its checkpoint bundle."
  fi
}

run_field_statistics_rank_equivalence_smoke() {
  local serial_run="" parallel_run=""
  local saved_nprocs="${nprocs}"
  local stats_dir="output/checkpoints/step_000000000004/statistics/window_0000/block_0000"
  local payload=""

  if [[ "${saved_nprocs}" -le 1 ]]; then
    echo "    (skipped: rank equivalence needs a multi-rank harness, have ${saved_nprocs})"
    return 0
  fi

  nprocs=1
  run_field_statistics_variant "field_statistics_serial" 1 2 "INFO"
  serial_run="${LAST_RUN_DIR}"
  nprocs="${saved_nprocs}"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : 1" \
    "statistics serial reference run rank count"

  run_field_statistics_variant "field_statistics_parallel" 1 2 "INFO"
  parallel_run="${LAST_RUN_DIR}"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${saved_nprocs}" \
    "statistics parallel run rank count"

  # Accumulation is pointwise and the payloads are written in natural ordering, so a
  # decomposition change may not move a single byte. The per-point count and weight
  # payloads carry the sampled-point mask, which is what would shift first if a rank
  # boundary were mistaken for a domain boundary.
  require_count_ge "${serial_run}/${stats_dir}" "*.dat" 5 "statistics payloads to compare"
  while IFS= read -r payload; do
    require_files_identical \
      "${serial_run}/${stats_dir}/${payload}" \
      "${parallel_run}/${stats_dir}/${payload}" \
      "statistics payload ${payload} (1 rank vs ${saved_nprocs} ranks)"
  done < <(find "${serial_run}/${stats_dir}" -type f -name '*.dat' -printf '%f\n' | sort)
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
    "${particle_case}/config/case.yml" \
    "${particle_case}/config/solver.yml" \
    "${particle_case}/config/monitor.yml" \
    "${particle_case}/config/post.yml" \
    "flat_particles_stress"
  require_count_ge "${LAST_RUN_DIR}/output/checkpoints" "position.dat" 2 "stress particle checkpoint payloads"
  require_count_ge "${LAST_RUN_DIR}/output/visualization" "*.vtp" 1 "stress particle VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of Particles         : 96" "stress particle count in runtime summary"

  "${picurv_exe}" init flat_channel --dest "${restart_case}" >/dev/null
  prepare_flat_restart_equivalence_case "${restart_case}" 0 3 "solve" "" "viz/restart_chain_base" 0 3
  run_case_workflow \
    "${restart_case}" \
    "${restart_case}/config/case.yml" \
    "${restart_case}/config/solver.yml" \
    "${restart_case}/config/monitor.yml" \
    "${restart_case}/config/post.yml" \
    "restart_chain_base"
  local restart_base_run="${LAST_RUN_DIR}"
  prepare_flat_restart_equivalence_case "${restart_case}" 2 2 "solve" "${restart_base_run}" "viz/restart_chain_mid" 2 4
  run_case_workflow \
    "${restart_case}" \
    "${restart_case}/config/case.yml" \
    "${restart_case}/config/solver.yml" \
    "${restart_case}/config/monitor.yml" \
    "${restart_case}/config/post.yml" \
    "restart_chain_mid" \
    "${restart_base_run}"
  local restart_mid_run="${LAST_RUN_DIR}"
  prepare_flat_restart_equivalence_case "${restart_case}" 4 2 "solve" "${restart_mid_run}" "viz/restart_chain_final" 4 6
  run_case_workflow \
    "${restart_case}" \
    "${restart_case}/config/case.yml" \
    "${restart_case}/config/solver.yml" \
    "${restart_case}/config/monitor.yml" \
    "${restart_case}/config/post.yml" \
    "restart_chain_final" \
    "${restart_mid_run}"
  require_file "${LAST_RUN_DIR}/logs/Continuity_Metrics.log" "restart-chain continuity metrics log"

  "${picurv_exe}" init flat_channel --dest "${parabolic_case}" >/dev/null
  prepare_flat_case_parabolic_stress "${parabolic_case}"
  run_case_workflow \
    "${parabolic_case}" \
    "${parabolic_case}/config/case.yml" \
    "${parabolic_case}/config/solver.yml" \
    "${parabolic_case}/config/monitor.yml" \
    "${parabolic_case}/config/post.yml" \
    "flat_parabolic_stress"
  require_count_ge "${LAST_RUN_DIR}/output/visualization" "*.vts" 1 "parabolic stress VTS files"

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
      --case "${periodic_case}/config/case.yml" \
      --solver "${periodic_case}/config/solver.yml" \
      --monitor "${periodic_case}/config/monitor.yml" \
      --post "${periodic_case}/config/post.yml" >/dev/null
    "${picurv_exe}" run \
      --solve \
      --post-process \
      -n "${periodic_nprocs}" \
      --case "${periodic_case}/config/case.yml" \
      --solver "${periodic_case}/config/solver.yml" \
      --monitor "${periodic_case}/config/monitor.yml" \
      --post "${periodic_case}/config/post.yml" \
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
  python3 - "${mpi_particle_case}/config/post.yml" <<'PY'
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
    "${mpi_particle_case}/config/case.yml" \
    "${mpi_particle_case}/config/solver.yml" \
    "${mpi_particle_case}/config/monitor.yml" \
    "${mpi_particle_case}/config/post.yml" \
    "flat_particles_stress_mpi"
  nprocs="${saved_nprocs}"
  require_count_ge "${LAST_RUN_DIR}/output/visualization" "*.vtp" 1 "MPI particle stress VTP files"
  require_file_contains "${LAST_SOLVER_LOG}" "Number of MPI Processes     : ${stress_mpi_nprocs}" "MPI stress rank count in runtime summary"
}

run_geometric_periodic_smoke() {
  local periodic_case="${tmp_root}/constant-uniform-geometric-periodic"
  local streamwise_case="${tmp_root}/constant-streamwise-geometric-periodic"
  local file_ic_case="${tmp_root}/file-ic-geometric-periodic"
  local source_ufield=""
  local saved_nprocs="${nprocs}"

  cp -R "${repo_root}/examples/periodic_test/constant_uniform_flow" "${periodic_case}"
  nprocs=1
  run_case_workflow \
    "${periodic_case}" \
    "${periodic_case}/case.yml" \
    "${periodic_case}/solver.yml" \
    "${periodic_case}/monitor.yml" \
    "${periodic_case}/post.yml" \
    "constant_uniform_geometric_periodic"
  require_file "${LAST_RUN_DIR}/logs/Continuity_Metrics.log" "geometric-periodic continuity metrics log"
  require_file_contains "${LAST_SOLVER_LOG}" "Periodic Axes (BC-derived)  : I=YES, J=YES, K=YES" "BC-derived periodic axes"
  source_ufield="${LAST_RUN_DIR}/output/checkpoints/step_000000000000/eulerian/block_0000/Ucat.dat"
  require_file "${source_ufield}" "file-IC source Ucat"

  cp -R "${repo_root}/examples/periodic_test/constant_uniform_flow" "${file_ic_case}"
  python3 - "${file_ic_case}/case.yml" "${source_ufield}" <<'PY'
import sys
import yaml

case_path, source_ufield = sys.argv[1:]
with open(case_path, "r", encoding="utf-8") as f:
    case_cfg = yaml.safe_load(f)
case_cfg["run_control"]["total_steps"] = 1
case_cfg["properties"]["initial_conditions"] = {
    "mode": "file",
    "field": "Ucat",
    "source_file": source_ufield,
}
with open(case_path, "w", encoding="utf-8") as f:
    yaml.safe_dump(case_cfg, f, sort_keys=False)
PY
  run_case_workflow \
    "${file_ic_case}" \
    "${file_ic_case}/case.yml" \
    "${file_ic_case}/solver.yml" \
    "${file_ic_case}/monitor.yml" \
    "${file_ic_case}/post.yml" \
    "file_ic_geometric_periodic"
  require_file_contains "${LAST_SOLVER_LOG}" "Eulerian State Source       : initial condition (File)" "file-IC banner source"
  require_file_contains "${LAST_SOLVER_LOG}" "Initial Velocity File       : field=Ucat" "file-IC banner field"

  cp -R "${repo_root}/examples/periodic_test/constant_streamwise_flow" "${streamwise_case}"
  run_case_workflow \
    "${streamwise_case}" \
    "${streamwise_case}/case.yml" \
    "${streamwise_case}/solver.yml" \
    "${streamwise_case}/monitor.yml" \
    "${streamwise_case}/post.yml" \
    "constant_streamwise_geometric_periodic"
  require_file_contains "${LAST_SOLVER_LOG}" "initial condition (Streamwise Constant)" "streamwise-IC banner source"
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
  standard|stress|periodic) ;;
  *) die "unknown smoke mode '${smoke_mode}' (expected 'standard', 'stress', or 'periodic')" ;;
esac

echo "==> PICurv smoke: simulator help"
run_help_smoke "${simulator_exe}" "PICurv Simulator"

echo "==> PICurv smoke: postprocessor help"
run_help_smoke "${postprocessor_exe}" "Unified Post-Processing Tool"

if [[ "${smoke_mode}" == "periodic" ]]; then
  echo "==> PICurv smoke: geometric-periodic runtime harness"
  run_geometric_periodic_smoke
  echo "PICurv geometric-periodic smoke completed successfully."
  exit 0
fi

echo "==> PICurv smoke: case init and source metadata"
run_case_init_smoke

echo "==> PICurv smoke: template matrix init/validate/dry-run checks"
run_template_matrix_smoke

echo "==> PICurv smoke: dry-run execution plan"
run_dry_run_plan_smoke

echo "==> PICurv smoke: PETSc diagnostics runtime gate"
run_petsc_diagnostics_smoke

echo "==> PICurv smoke: restart source resolution in dry-run plan"
run_restart_resolution_smoke

echo "==> PICurv smoke: field statistics accumulation, monitoring, and post-processing"
run_field_statistics_smoke

if [[ "${nprocs}" -gt 1 ]]; then
  echo "==> PICurv smoke: multi-rank runtime sequences (flat+bent)"
  run_multi_rank_runtime_smoke
  echo "==> PICurv smoke: restart across a changed MPI rank count"
  run_rank_change_restart_smoke
  echo "==> PICurv smoke: field statistics across a changed MPI rank count"
  run_field_statistics_rank_equivalence_smoke
else
  echo "==> PICurv smoke: Newton--Krylov flat-channel BDF1 startup"
  run_newton_krylov_flat_channel_startup_smoke
  echo "==> PICurv smoke: full end-to-end runtime sequences"
  run_full_runtime_smoke
fi

if [[ "${smoke_mode}" == "stress" ]]; then
  echo "==> PICurv smoke: stress extensions"
  run_stress_smoke
fi

echo "PICurv smoke completed successfully."
