#!/usr/bin/env bash

set -euo pipefail

simulator_exe="${1:?missing simulator path}"
postprocessor_exe="${2:?missing postprocessor path}"
mpi_launcher="${3:?missing MPI launcher}"
nprocs="${4:-1}"

echo "==> PICurv smoke: simulator help"
"${mpi_launcher}" -n "${nprocs}" "${simulator_exe}" -help >/dev/null

echo "==> PICurv smoke: postprocessor help"
"${mpi_launcher}" -n "${nprocs}" "${postprocessor_exe}" -help >/dev/null

echo "PICurv smoke completed successfully."
