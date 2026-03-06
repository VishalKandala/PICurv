#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

INSTALL_PETSC=0
SKIP_SYSTEM_DEPS=0
PETSC_VERSION="3.20.3"
PETSC_PREFIX="${HOME}/software"
PETSC_ARCH="arch-linux-c-debug"
PYTHON_BIN=""

usage() {
  cat <<'EOF'
Usage: scripts/bootstrap_install.sh [options]

Automates local PICurv setup:
1) installs base system dependencies (Debian/Ubuntu),
2) installs required Python packages,
3) optionally installs PETSc with DMSwarm,
4) verifies PETSc + DMSwarm visibility,
5) builds PICurv binaries.

Options:
  --install-petsc            Build/install PETSc from source.
  --petsc-version <ver>      PETSc version tag (default: 3.20.3).
  --petsc-prefix <dir>       Parent directory for PETSc source/build.
  --petsc-arch <arch>        PETSc arch name (default: arch-linux-c-debug).
  --python-bin <path>        Python interpreter to use.
  --skip-system-deps         Skip apt package installation.
  -h, --help                 Show this help.

Examples:
  scripts/bootstrap_install.sh
  scripts/bootstrap_install.sh --install-petsc
EOF
}

log() {
  printf '[INFO] %s\n' "$*"
}

die() {
  printf '[ERROR] %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command '$1' is missing."
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --install-petsc) INSTALL_PETSC=1; shift ;;
    --skip-system-deps) SKIP_SYSTEM_DEPS=1; shift ;;
    --petsc-version) PETSC_VERSION="${2:?missing value for --petsc-version}"; shift 2 ;;
    --petsc-prefix) PETSC_PREFIX="${2:?missing value for --petsc-prefix}"; shift 2 ;;
    --petsc-arch) PETSC_ARCH="${2:?missing value for --petsc-arch}"; shift 2 ;;
    --python-bin) PYTHON_BIN="${2:?missing value for --python-bin}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

if [[ -z "${PYTHON_BIN}" ]]; then
  if command -v python3.10 >/dev/null 2>&1; then
    PYTHON_BIN="python3.10"
  else
    PYTHON_BIN="python3"
  fi
fi

require_cmd "${PYTHON_BIN}"
require_cmd git

if [[ "${SKIP_SYSTEM_DEPS}" -eq 0 ]]; then
  require_cmd apt-get
  log "Installing base system dependencies (Debian/Ubuntu)..."
  sudo apt-get update
  sudo apt-get install -y \
    build-essential gfortran mpich git make pkg-config \
    libx11-dev python3 python3-pip python3-venv
fi

log "Installing Python dependencies with ${PYTHON_BIN}..."
"${PYTHON_BIN}" -m pip install --upgrade pip
"${PYTHON_BIN}" -m pip install pyyaml numpy

if [[ "${INSTALL_PETSC}" -eq 1 ]]; then
  require_cmd mpicc
  PETSC_SRC="${PETSC_PREFIX}/petsc-${PETSC_VERSION}"
  mkdir -p "${PETSC_PREFIX}"

  if [[ ! -d "${PETSC_SRC}" ]]; then
    log "Cloning PETSc v${PETSC_VERSION} into ${PETSC_SRC}..."
    git clone -b "v${PETSC_VERSION}" https://gitlab.com/petsc/petsc.git "${PETSC_SRC}"
  else
    log "Using existing PETSc source at ${PETSC_SRC}."
  fi

  log "Configuring PETSc with DMSwarm support..."
  pushd "${PETSC_SRC}" >/dev/null
  ./configure --PETSC_ARCH="${PETSC_ARCH}" \
              --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
              --download-fblaslapack --download-metis --download-parmetis \
              --with-dmswarm=1 --with-debugging=1
  make all
  make check
  popd >/dev/null

  export PETSC_DIR="${PETSC_SRC}"
  export PETSC_ARCH="${PETSC_ARCH}"
fi

if [[ -z "${PETSC_DIR:-}" || -z "${PETSC_ARCH:-}" ]]; then
  die "PETSC_DIR and PETSC_ARCH must be set (or use --install-petsc)."
fi

PETSC_CONF="${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h"
PETSC_DMSWARM_HEADER="${PETSC_DIR}/include/petscdmswarm.h"
[[ -f "${PETSC_CONF}" ]] || die "PETSc config not found: ${PETSC_CONF}"
[[ -f "${PETSC_DMSWARM_HEADER}" ]] || die "DMSwarm header not found: ${PETSC_DMSWARM_HEADER}"

log "Verified PETSc config and DMSwarm header."

log "Building PICurv binaries..."
cd "${REPO_ROOT}"
"${PYTHON_BIN}" ./scripts/picurv build

[[ -x "${REPO_ROOT}/bin/simulator" ]] || die "Missing binary: bin/simulator"
[[ -x "${REPO_ROOT}/bin/postprocessor" ]] || die "Missing binary: bin/postprocessor"
[[ -x "${REPO_ROOT}/bin/picurv" ]] || die "Missing binary: bin/picurv"

log "Bootstrap complete."
log "Run: ./bin/picurv --help"
