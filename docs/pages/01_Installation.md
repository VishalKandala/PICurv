@page 01_Installation Installation Guide

@anchor _Installation

This guide covers dependency setup, PETSc configuration, and binary build verification for PICurv.

@tableofcontents

@section p01_prereqs_sec 1. Prerequisites

Install these first:

- C compiler (`gcc` or `clang`)
- MPI implementation (`mpich` or `openmpi`)
- GNU Make
- Python 3.10+ + pip for the managed PICurv CLI environment
- Git
- PETSc 3.20.3+ with `DMSwarm` support

@section p01_automated_sec 2. Automated Install (Recommended)

From the PICurv repo root:

```bash
export PETSC_DIR=/path/to/petsc
# Set PETSC_ARCH only for old-style in-tree PETSc builds:
# export PETSC_ARCH=arch-linux-c-debug
./scripts/bootstrap_install.sh --install-shell-hook
```

If PETSc is not installed yet, let the script build it:

```bash
./scripts/bootstrap_install.sh --install-petsc
```

The script installs system and Python dependencies, verifies PETSc/DMSwarm visibility, then builds:

- `bin/simulator` (compiled C solver)
- `bin/postprocessor` (compiled C post-processor)
- `bin/picurv` (launcher for `scripts/picurv`, the Python conductor)

By default, bootstrap creates `.picurv-venv/` under the repo and installs the
Python-side CLI dependencies there. PETSc, MPI, compilers, and scheduler tools
remain provided by your loaded system or cluster modules. Bootstrap also writes
`.picurv-python-env`, which records the seed Python runtime library path needed
to launch the managed venv after you switch to a different module stack.

On an existing HPC cluster where modules already provide compilers, MPI, and
PETSc, skip OS package installation:

```bash
module load <compiler-mpi-petsc-stack>
./scripts/bootstrap_install.sh --skip-system-deps --install-shell-hook
source ~/.bashrc
picurv --help
```

Useful variants:

```bash
./scripts/bootstrap_install.sh --with-plotting
./scripts/bootstrap_install.sh --venv-dir /path/to/picurv-venv
./scripts/bootstrap_install.sh --python-bin /path/to/python3.11
./scripts/bootstrap_install.sh --install-shell-hook
./scripts/bootstrap_install.sh --no-venv
```

Use `--no-venv` when your site requires Python packages to come from modules or
a centrally managed environment. If the only visible interpreter is Python 3.6,
load a newer Python module before the default bootstrap path; otherwise use
`--no-venv` and site-approved package versions.

@section p01_install_tools_sec 3. Install Base Toolchain (Manual Path)

Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran mpich git make pkg-config libx11-dev python3 python3-pip python3-venv
python3 -m venv .picurv-venv
.picurv-venv/bin/python -m pip install --upgrade pip
.picurv-venv/bin/python -m pip install pyyaml numpy
```

RHEL/CentOS/Fedora:

```bash
sudo yum groupinstall -y "Development Tools"
sudo yum install -y mpich-devel python3 python3-pip git
python3 -m venv .picurv-venv
.picurv-venv/bin/python -m pip install --upgrade pip
.picurv-venv/bin/python -m pip install pyyaml numpy
```

Optional (study plot generation):

```bash
.picurv-venv/bin/python -m pip install matplotlib
```

@section p01_petsc_sec 4. Install PETSc

Recommended source install:

```bash
git clone -b v3.20.3 https://gitlab.com/petsc/petsc.git
cd petsc
```

Debug build example:

```bash
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            --download-fblaslapack --download-metis --download-parmetis \
            --with-dmswarm=1 --with-debugging=1
make all
make check
```

Optimized build example:

```bash
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            --download-fblaslapack --download-metis --download-parmetis \
            --with-dmswarm=1 --with-debugging=0 \
            --COPTFLAGS='-O3' --CXXOPTFLAGS='-O3' --FOPTFLAGS='-O3'
make all
make check
```

Official references:

- https://petsc.org/release/install/
- https://petsc.org/release/docs/manual/

@section p01_env_sec 5. Configure Environment Variables

Add to your shell profile (`~/.bashrc` or equivalent):

```bash
export PETSC_DIR=/path/to/petsc
# PETSC_ARCH is optional for prefix installs from EasyBuild/Spack/system packages.
# export PETSC_ARCH=arch-linux-c-debug
source /path/to/PICurv/etc/picurv.sh
```

The `etc/picurv.sh` script sets `PICURV_DIR`, exports `PICURV_PYTHON` when a
managed venv or bootstrap-selected Python is available, adds `bin/` to your
`PATH` for compiled executables, and also exposes `scripts/` as a fallback so
`picurv` still resolves if `bin/picurv` is temporarily absent before a rebuild.
It is idempotent and safe to source multiple times.

If you want bootstrap to add this setup to `~/.bashrc`, pass
`--install-shell-hook`. The hook is written as a managed block, so rerunning the
installer updates it instead of appending duplicate source lines. Use
`--shell-rc <path>` to target another shell startup file.

Reload and verify:

```bash
source ~/.bashrc
echo "$PETSC_DIR"
echo "$PETSC_ARCH"
picurv --help
```

Verify PETSc has DMSwarm headers:

```bash
test -f "$PETSC_DIR/include/petscdmswarm.h" && echo "DMSwarm header found"
test -f "$PETSC_DIR/$PETSC_ARCH/include/petscconf.h" || test -f "$PETSC_DIR/include/petscconf.h"
```

@section p01_clone_sec 6. Clone PICurv

```bash
git clone https://github.com/VishalKandala/PICurv.git
cd PICurv
```

@section p01_build_sec 7. Build with picurv

```bash
./scripts/picurv build
```

Expected binaries:

- `bin/simulator` (compiled C solver)
- `bin/postprocessor` (compiled C post-processor)
- `bin/picurv` (launcher → `scripts/picurv`)

Useful variants:

```bash
./scripts/picurv build clean-project
./scripts/picurv build SYSTEM=cluster
make audit-build
```

After `source etc/picurv.sh`, use `picurv` directly from any directory.
`./scripts/picurv build` writes `logs/build.log` in the source repo. If you are
auditing compiler warnings from a direct Make invocation, use `make audit-build`
to generate both `logs/build.log` and `logs/build.warnings.log`.

@section p01_verify_sec 8. Verify Installation

```bash
make doctor
```

Recommended sequence after a successful build:

1. `make doctor`
2. `make smoke` (stronger executable-level sanity check)
3. `make check` (full local validation sweep)
4. `make check-full` (comprehensive MPI-inclusive validation sweep)

What these prove:

- `make doctor` proves the local PETSc installation can build and run a minimal PETSc-backed binary.
- `make smoke` proves the compiled PICurv executables still launch.
- `make check` runs Python regressions plus PETSc-backed validation.
- `make check-full` additionally validates dedicated MPI unit coverage, fixed-size multi-rank smoke, and rank-matrix smoke in one pass.

What `make doctor` does not prove:

- it does not prove a full solver case is numerically correct
- it does not replace case-specific validation or convergence testing

@section p01_common_sec 9. Common Installation Failures

- `PETSC_DIR` unset or pointing at the wrong PETSc installation; `PETSC_ARCH` is only required for old-style in-tree PETSc builds.
- MPI compiler wrappers unavailable in PATH.
- Python interpreter too old for default bootstrap (load Python 3.10+, pass `--python-bin`, or use `--no-venv`).
- Visualization modules leaking incompatible Python packages into `PYTHONPATH`; prefer the managed venv launcher from bootstrap for normal CLI use.
- Managed venv Python cannot find `libpython`; rerun bootstrap after pulling current source so `.picurv-python-env` records the seed Python library path.
- Old checkouts where `bin/picurv` is still a symlink; rerun `make -B conductor` after pulling current source.
- PETSc configured without required downloaded dependencies.
- missing X11 dev library at link time (`cannot find -lX11` on Linux; install `libx11-dev`).
- stale object files after toolchain changes (use `clean-project`).

For runtime-level failures after successful build, see **@subpage 39_Common_Fatal_Errors**.

For the full testing model after installation, see **@subpage 40_Testing_and_Quality_Guide**.

@section p01_next_steps_sec 10. Next Steps

- First simulation walkthrough: **@subpage 02_Tutorial_Programmatic_Grid**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Installation Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
