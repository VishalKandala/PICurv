@page 01_Installation Installation

@anchor _Installation

@pagemeta{How-to, New users, Local and HPC}

@htmlonly
<header class="pic-install-hero">
  <div class="pic-install-hero-copy">
    <span class="pic-eyebrow">Installation</span>
    <h1>Build the solver.<br />Prove the toolchain.</h1>
    <p>Set up PICurv, connect it to PETSc and MPI, and leave the installation with three verified command-line tools.</p>
    <nav class="pic-install-actions" aria-label="Installation shortcuts">
      <a class="pic-button pic-button-primary" href="#p01_automated_sec">Use the bootstrap installer</a>
      <a class="pic-text-link" href="#p01_install_tools_sec">Follow the manual path <span aria-hidden="true">→</span></a>
    </nav>
  </div>
  <div class="pic-install-outcome" aria-label="A completed PICurv installation">
    <span class="pic-install-outcome-label">A complete install</span>
    <div class="pic-install-stack" aria-hidden="true">
      <span>C compiler</span><b>+</b><span>MPI</span><b>+</b><span>PETSc</span>
    </div>
    <div class="pic-install-binaries">
      <code>simulator</code>
      <code>postprocessor</code>
      <code>picurv</code>
    </div>
    <p><strong>Compile.</strong> <strong>Launch.</strong> <strong>Validate.</strong></p>
  </div>
</header>

<section class="pic-install-route" aria-labelledby="pic-install-route-title">
  <header>
    <span class="pic-eyebrow">The route</span>
    <h2 id="pic-install-route-title">From dependencies to a verified build.</h2>
  </header>
  <nav aria-label="Installation stages">
    <a href="#p01_prereqs_sec"><b>01</b><span><strong>Prepare</strong><small>Check the toolchain</small></span></a>
    <a href="#p01_clone_sec"><b>02</b><span><strong>Clone</strong><small>Get the source</small></span></a>
    <a href="#p01_automated_sec"><b>03</b><span><strong>Install</strong><small>Bootstrap PICurv</small></span></a>
    <a href="#p01_verify_sec"><b>04</b><span><strong>Verify</strong><small>Exercise the build</small></span></a>
  </nav>
</section>

<section class="pic-install-paths" aria-labelledby="pic-install-paths-title">
  <header>
    <span class="pic-eyebrow">Choose your path</span>
    <h2 id="pic-install-paths-title">Match the install to your machine.</h2>
  </header>
  <div>
    <a class="pic-install-path pic-install-path--recommended" href="#p01_automated_sec">
      <span>Recommended</span><strong>Local workstation</strong>
      <small>Use bootstrap to create the managed Python environment and build the tools.</small><b aria-hidden="true">→</b>
    </a>
    <a class="pic-install-path" href="#p01_automated_sec">
      <span>Existing modules</span><strong>HPC cluster</strong>
      <small>Keep the site toolchain and ask bootstrap to skip operating-system packages.</small><b aria-hidden="true">→</b>
    </a>
    <a class="pic-install-path" href="#p01_install_tools_sec">
      <span>Full control</span><strong>Manual setup</strong>
      <small>Install the base tools and PETSc yourself, then build through the PICurv CLI.</small><b aria-hidden="true">→</b>
    </a>
  </div>
</section>
@endhtmlonly

@tableofcontents

@section p01_prereqs_sec 1. Prerequisites

PICurv sits on a compact scientific-computing toolchain. Have these components
available before starting the build:

@htmlonly
<div class="pic-dependency-grid" aria-label="PICurv installation prerequisites">
  <article><span>Compiler</span><strong>C compiler</strong><small><code>gcc</code> or <code>clang</code></small></article>
  <article><span>Parallel runtime</span><strong>MPI</strong><small><code>mpich</code> or <code>openmpi</code></small></article>
  <article><span>Build</span><strong>GNU Make</strong><small><code>make</code></small></article>
  <article><span>CLI runtime</span><strong>Python 3.10+</strong><small>with <code>pip</code></small></article>
  <article><span>Source control</span><strong>Git</strong><small><code>git</code></small></article>
  <article><span>Numerics</span><strong>PETSc 3.20.3+</strong><small>with <code>DMSwarm</code></small></article>
</div>
@endhtmlonly

If PETSc is the only missing component, bootstrap can build it during the
automated path.

@section p01_clone_sec 2. Clone PICurv

```bash
git clone https://github.com/VishalKandala/PICurv.git
cd PICurv
```

The remaining commands run from this repository root.

@section p01_automated_sec 3. Automated Install (Recommended)

@htmlonly
<p class="pic-section-lead">The bootstrap path is the shortest route to a repeatable local installation. It checks the native toolchain, isolates the Python CLI, and builds the three PICurv launch targets.</p>
@endhtmlonly

From the PICurv repo root:

```bash
export PETSC_DIR=/path/to/petsc
# Set PETSC_ARCH only for old-style in-tree PETSc builds:
# export PETSC_ARCH=arch-linux-c-debug
./bootstrap_install.sh --install-shell-hook
```

If PETSc is not installed yet, let the script build it:

```bash
./bootstrap_install.sh --install-petsc
```

The script installs system and Python dependencies, verifies PETSc/DMSwarm visibility, then builds:

- `bin/simulator` (compiled C solver)
- `bin/postprocessor` (compiled C post-processor)
- `bin/picurv` (launcher for `picurv_cli/picurv`, the Python conductor entrypoint)

By default, bootstrap creates `.picurv-venv/` under the repo and installs the
Python-side CLI dependencies there. PETSc, MPI, compilers, and scheduler tools
remain provided by your loaded system or cluster modules. Bootstrap also writes
`.picurv-python-env`, which records the seed Python runtime library path needed
to launch the managed venv after you switch to a different module stack.

On an existing HPC cluster where modules already provide compilers, MPI, and
PETSc, skip OS package installation:

```bash
module load <compiler-mpi-petsc-stack>
./bootstrap_install.sh --skip-system-deps --install-shell-hook
source ~/.bashrc
picurv --help
```

Useful variants:

```bash
./bootstrap_install.sh --venv-dir /path/to/picurv-venv
./bootstrap_install.sh --python-bin /path/to/python3.11
./bootstrap_install.sh --upgrade-pip
./bootstrap_install.sh --install-shell-hook
./bootstrap_install.sh --no-venv
```

Use `--no-venv` when your site requires Python packages to come from modules or
a centrally managed environment. If the only visible interpreter is Python 3.6,
load a newer Python module before the default bootstrap path; otherwise use
`--no-venv` and site-approved package versions.

Bootstrap does not upgrade pip during a normal install. This keeps routine
updates smaller and avoids unnecessary package churn on quota-constrained
cluster home directories. Pass `--upgrade-pip` when an explicit upgrade is
needed.

Bootstrap writes `.picurv-python` and `.picurv-python-env` atomically, then checks
that `picurv`, `simulator`, and `postprocessor` all report the release in the root
`VERSION` file. A partial interpreter record or mixed-version native build therefore
fails installation instead of becoming a latent cluster-launch error.

@section p01_install_tools_sec 4. Install Base Toolchain (Manual Path)

@htmlonly
<aside class="pic-install-manual-note"><span>Manual path</span><p>Use the remaining setup sections when system packages, cluster policy, or a custom PETSc build make bootstrap unsuitable.</p><a href="#p01_automated_sec">Back to the recommended path ↑</a></aside>
@endhtmlonly

Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran mpich git make pkg-config libx11-dev python3 python3-pip python3-venv
python3 -m venv .picurv-venv
.picurv-venv/bin/python -m pip install --upgrade pip
.picurv-venv/bin/python -m pip install pyyaml numpy packaging matplotlib
```

RHEL/CentOS/Fedora:

```bash
sudo yum groupinstall -y "Development Tools"
sudo yum install -y mpich-devel python3 python3-pip git
python3 -m venv .picurv-venv
.picurv-venv/bin/python -m pip install --upgrade pip
.picurv-venv/bin/python -m pip install pyyaml numpy packaging matplotlib
```

`matplotlib` is part of the standard Python dependency set because it powers
`summarize --plot` and study plot generation.

Bootstrap verifies `yaml`, `numpy`, `packaging`, and `matplotlib` imports with
`PYTHONPATH` and user-site packages disabled. This matches the isolated Python
environment used by the `picurv` launcher and catches dependencies that are
visible only through a loaded cluster module.

@section p01_petsc_sec 5. Install PETSc

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

@section p01_env_sec 6. Configure Environment Variables

Add to your shell profile (`~/.bashrc` or equivalent):

```bash
export PETSC_DIR=/path/to/petsc
# PETSC_ARCH is optional for prefix installs from EasyBuild/Spack/system packages.
# export PETSC_ARCH=arch-linux-c-debug
source /path/to/PICurv/etc/picurv.sh
```

The `etc/picurv.sh` script sets `PICURV_DIR`, exports `PICURV_PYTHON` when a
managed venv or bootstrap-selected Python is available, adds `bin/` to your
`PATH` for compiled executables, and also exposes `picurv_cli/` as a fallback so
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

@section p01_build_sec 7. Build with picurv

```bash
./picurv_cli/picurv build
```

Expected binaries:

- `bin/simulator` (compiled C solver)
- `bin/postprocessor` (compiled C post-processor)
- `bin/picurv` (launcher → `picurv_cli/picurv`)

Useful variants:

```bash
./picurv_cli/picurv build clean-project
./picurv_cli/picurv build SYSTEM=cluster
make audit-build
```

After `source etc/picurv.sh`, use `picurv` directly from any directory.
`./picurv_cli/picurv build` writes `<repo>/logs/build.log` in the source repo. If you are
auditing compiler warnings from a direct Make invocation, use `make audit-build`
to generate both `<repo>/logs/build.log` and `<repo>/logs/build.warnings.log`.

@section p01_verify_sec 8. Verify Installation

@htmlonly
<p class="pic-section-lead">Verification progresses from the installed PETSc toolchain to the complete MPI validation sweep. Start small, then choose the depth appropriate for your machine.</p>
@endhtmlonly

```bash
make doctor
```

Recommended sequence after a successful build:

```bash
picurv version
simulator --version
postprocessor --version
```

All three must report the same release. `picurv version --format json` also reports
the Git commit, dirty-tree state, active workspace, and optional workspace requirement.

@htmlonly
<ol class="pic-verification-ladder">
  <li><span>01</span><code>make doctor</code><p>Builds and runs a minimal PETSc-backed binary.</p></li>
  <li><span>02</span><code>make smoke</code><p>Confirms that the compiled PICurv executables launch.</p></li>
  <li><span>03</span><code>make check</code><p>Runs Python regressions and PETSc-backed validation.</p></li>
  <li><span>04</span><code>make check-full</code><p>Adds MPI units, fixed-size multi-rank smoke, and rank-matrix smoke.</p></li>
</ol>
@endhtmlonly

What `make doctor` does not prove:

- it does not prove a full solver case is numerically correct
- it does not replace case-specific validation or convergence testing

@section p01_common_sec 9. Common Installation Failures

@htmlonly
<div class="pic-failure-list">
  <ul>
    <li><code>PETSC_DIR</code> is unset or points at the wrong PETSc installation; <code>PETSC_ARCH</code> is only required for old-style in-tree PETSc builds.</li>
    <li>MPI compiler wrappers are unavailable in <code>PATH</code>.</li>
    <li>The Python interpreter is too old for default bootstrap; load Python 3.10+, pass <code>--python-bin</code>, or use <code>--no-venv</code>.</li>
    <li>Visualization modules leak incompatible Python packages into <code>PYTHONPATH</code>; prefer the managed venv launcher for normal CLI use.</li>
    <li>The managed venv cannot find <code>libpython</code>; rerun bootstrap so <code>.picurv-python-env</code> records the seed runtime library path.</li>
    <li>An old checkout still uses a <code>bin/picurv</code> symlink; rerun <code>make -B conductor</code> after pulling current source.</li>
    <li>PETSc was configured without the required downloaded dependencies.</li>
    <li>The X11 development library is missing at link time; install <code>libx11-dev</code> when Linux reports <code>cannot find -lX11</code>.</li>
    <li>Object files are stale after a toolchain change; use <code>clean-project</code>.</li>
  </ul>
</div>
@endhtmlonly

For runtime-level failures after successful build, see **@subpage 39_Common_Fatal_Errors**.

For the full testing model after installation, see **@subpage 40_Testing_and_Quality_Guide**.

@section p01_next_steps_sec 10. Next Steps

@htmlonly
<nav class="pic-next-cards pic-install-next-cards" aria-label="Recommended next documentation">
  <a href="41_Getting_Started_Index.html">
    <span>Run your first case</span>
    <strong>Create, validate, solve, and inspect the Quick Start simulation.</strong>
  </a>
  <a href="67_Troubleshooting.html">
    <span>Resolve a failure</span>
    <strong>Trace installation, configuration, MPI, and runtime symptoms.</strong>
  </a>
  <a href="40_Testing_and_Quality_Guide.html">
    <span>Understand the checks</span>
    <strong>See what each test layer verifies and where its evidence comes from.</strong>
  </a>
</nav>

<nav class="pic-page-pager" aria-label="Manual page navigation">
  <a class="pic-page-pager-next" href="41_Getting_Started_Index.html">
    <span>Next</span>
    <strong>Quick Start →</strong>
  </a>
</nav>
@endhtmlonly
