#!/usr/bin/env bash

set -euo pipefail

simulator_exe="${1:?missing simulator path}"
postprocessor_exe="${2:?missing postprocessor path}"
mpi_launcher="${3:?missing MPI launcher}"
nprocs="${4:-1}"

run_help_smoke() {
  local exe="$1"
  local banner="$2"
  local out_file
  out_file="$(mktemp)"
  set +e
  "${mpi_launcher}" -n "${nprocs}" "${exe}" -help >"${out_file}" 2>&1
  local rc=$?
  set -e

  if ! grep -q "${banner}" "${out_file}"; then
    echo "Smoke failure: missing banner '${banner}' in ${exe} output." >&2
    sed -n '1,80p' "${out_file}" >&2
    rm -f "${out_file}"
    exit 1
  fi

  if [[ ${rc} -ne 0 ]]; then
    echo "Smoke note: ${exe} -help exited ${rc}, but banner check passed."
  fi

  rm -f "${out_file}"
}

echo "==> PICurv smoke: simulator help"
run_help_smoke "${simulator_exe}" "PICurv Simulator"

echo "==> PICurv smoke: postprocessor help"
run_help_smoke "${postprocessor_exe}" "Unified Post-Processing Tool"

echo "PICurv smoke completed successfully."
