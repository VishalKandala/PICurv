@page 41_Getting_Started_Index Quick Start

@anchor _Getting_Started_Index
@anchor p41_quicklook_sec

@pagemeta{Tutorial, New users, Current workflow}

@htmlonly
<header class="pic-quickstart-intro">
  <div class="pic-quickstart-stage">
    <img src="paraview_flat_channel.png" alt="Cropped preview of the flat-channel velocity field produced by this Quick Start" />
    <div class="pic-quickstart-copy">
      <span class="pic-eyebrow">Quick Start</span>
      <h1>Your first simulation.</h1>
      <p class="pic-quickstart-lead">Create, validate, run, and inspect a complete
      PICurv case from a clean workspace.</p>
      <nav class="pic-quickstart-actions" aria-label="Quick Start shortcuts">
        <a class="pic-button pic-button-primary" href="#p41_prereq_sec">Start the walkthrough</a>
        <a class="pic-text-link" href="#p41_result_sec">Preview the result <span aria-hidden="true">↓</span></a>
      </nav>
    </div>
    <p class="pic-preview-label">Cropped preview · <strong>Ucat_nodal</strong> at step 20</p>
  </div>
  <div class="pic-quickstart-context" aria-label="Simulation overview">
    <article>
      <span class="pic-context-label">The case</span>
      <h2>Laminar flat-channel flow</h2>
      <p>Flow at <strong>Re = 200</strong> moves through a straight rectangular channel with no-slip walls, a constant-velocity inlet, and a conservative outlet.</p>
      <dl class="pic-case-facts">
        <div><dt>9 × 9 × 17</dt><dd>grid cells</dd></div>
        <div><dt>20</dt><dd>timesteps</dd></div>
        <div><dt>2</dt><dd>MPI ranks</dd></div>
      </dl>
    </article>
    <article>
      <span class="pic-context-label">The solver</span>
      <h2>Dual-time + multigrid</h2>
      <p>Dual-time Picard–Jameson RK advances momentum. FGMRES with a three-level geometric-multigrid preconditioner solves the pressure correction.</p>
      <a href="08_Solver_Reference.html">Explore the solver reference <span aria-hidden="true">→</span></a>
    </article>
  </div>
</header>

<section class="pic-quickstart-route" aria-labelledby="pic-quickstart-route-title">
  <header>
    <span class="pic-eyebrow">The route</span>
    <h2 id="pic-quickstart-route-title">One case. Four clear moves.</h2>
  </header>
  <nav aria-label="Quick Start steps">
    <a href="#p41_init_sec"><b>01</b><span><strong>Create</strong><small>Copy the starter case</small></span></a>
    <a href="#p41_validate_sec"><b>02</b><span><strong>Validate</strong><small>Check every profile</small></span></a>
    <a href="#p41_run_sec"><b>03</b><span><strong>Run</strong><small>Solve and post-process</small></span></a>
    <a href="#p41_result_sec"><b>04</b><span><strong>Inspect</strong><small>Open the VTK result</small></span></a>
  </nav>
</section>
@endhtmlonly

@tableofcontents

@section p41_prereq_sec Before you start

PICurv must be built with an MPI runtime available. Load the project environment:

@htmlonly
<div class="pic-command-block" data-label="Terminal">
@endhtmlonly
```bash
source <repo>/etc/picurv.sh
picurv --help
```
@htmlonly
</div>
@endhtmlonly

If `picurv --help` fails, complete the @ref 01_Installation "Installation guide".
Run the remaining commands from the workspace where the new case should live.

@section p41_init_sec 1. Create a case

Initialize the checked-in flat-channel template and enter the new case directory:

@htmlonly
<div class="pic-command-block" data-label="Terminal">
@endhtmlonly
```bash
picurv init flat_channel --dest my-first-run
cd my-first-run
```
@htmlonly
</div>
@endhtmlonly

A successful initialization ends with:

@htmlonly
<pre class="pic-terminal-output" aria-label="Successful initialization output"><code>[SUCCESS] Case directory is ready.
          Runtime binaries (simulator, postprocessor) are resolved from bin/ automatically.</code></pre>
@endhtmlonly

This run uses four profiles:

@htmlonly
<div class="pic-profile-grid" aria-label="The four configuration profiles used by this Quick Start">
  <article class="pic-profile-card">
    <span>Case</span>
    <code>quickstart_flat_channel.yml</code>
    <p>Grid and physics.</p>
  </article>
  <article class="pic-profile-card">
    <span>Solver</span>
    <code>Imp-MG-Standard.yml</code>
    <p>Numerical methods.</p>
  </article>
  <article class="pic-profile-card">
    <span>Monitor</span>
    <code>quickstart_Standard_Output.yml</code>
    <p>Logs and checkpoints.</p>
  </article>
  <article class="pic-profile-card">
    <span>Post-process</span>
    <code>quickstart_standard_analysis.yml</code>
    <p>VTK output.</p>
  </article>
</div>
@endhtmlonly

@section p41_validate_sec 2. Validate the configuration

Check the four profiles together:

@htmlonly
<div class="pic-command-block" data-label="Terminal">
@endhtmlonly
```bash
picurv validate \
  --case quickstart_flat_channel.yml \
  --solver Imp-MG-Standard.yml \
  --monitor quickstart_Standard_Output.yml \
  --post quickstart_standard_analysis.yml
```
@htmlonly
</div>
@endhtmlonly

The final line is:

@htmlonly
<pre class="pic-terminal-output" aria-label="Successful validation output"><code>[SUCCESS] Validation completed for 5 file(s).</code></pre>
@endhtmlonly

@section p41_run_sec 3. Run the solver and postprocessor

Launch the 20-step case on two MPI ranks and post-process its final checkpoint:

@htmlonly
<div class="pic-command-block" data-label="Terminal">
@endhtmlonly
```bash
picurv run --solve --post-process -n 2 \
  --case quickstart_flat_channel.yml \
  --solver Imp-MG-Standard.yml \
  --monitor quickstart_Standard_Output.yml \
  --post quickstart_standard_analysis.yml
```
@htmlonly
</div>
@endhtmlonly

The quickstart grid contains `9 x 9 x 17` cells. A completed run ends with:

@htmlonly
<pre class="pic-terminal-output" aria-label="Abbreviated successful run output"><code>[INFO] Created new self-contained run directory:
       runs/quickstart_flat_channel_&lt;timestamp&gt;

Progress: [==================================================] 100% (Step 20/20)
[SUCCESS] Execution finished successfully.

RUN SUMMARY
  Stages          : solve, post-process
  Solver MPI procs: 2
  Post MPI procs  : 2
  Steps run       : 20
  Post output     : runs/quickstart_flat_channel_&lt;timestamp&gt;/visualization/standard_analysis</code></pre>
@endhtmlonly

@section p41_result_sec 4. View the result

The postprocessor writes one VTK structured-grid file for the final step:

```text
runs/quickstart_flat_channel_<timestamp>/visualization/standard_analysis/eulerian_data_00020.vts
```

Open the file in ParaView, add a `Slice`, and color it by `Ucat_nodal`:

@htmlonly
<figure class="pic-result-figure">
  <img src="paraview_flat_channel.png" alt="Velocity-field slice from the PICurv flat-channel Quick Start in ParaView" />
  <figcaption>Velocity in the flat-channel Quick Start after nodal averaging at step 20.</figcaption>
</figure>
@endhtmlonly

Use the @ref 04_Visualization_Tutorial "Visualization tutorial" for the complete ParaView workflow and
other output fields.

@section p41_understand_sec The pattern behind every simulation

This Quick Start is deliberately small, but it does not use a special execution path.
Every PICurv simulation follows the same workflow: define its case and supporting profiles,
validate them together, materialize a self-contained run, execute the MPI solver, and
post-process checkpoints into analysis-ready output.

@htmlonly
<ol class="pic-run-pipeline" aria-label="PICurv simulation workflow">
  <li>Define profiles</li>
  <li>Validate controls</li>
  <li>Self-contained run</li>
  <li>MPI solve</li>
  <li>Write checkpoints</li>
  <li>Post-process outputs</li>
</ol>
@endhtmlonly

The grid, physics, solver, monitoring, and analysis choices can change from one case to
another; this simulation workflow stays the same.

For a failed command or incomplete run, use @ref 67_Troubleshooting "Troubleshooting".

@section p41_next_sec Where to go next

@htmlonly
<nav class="pic-next-cards" aria-label="Recommended next documentation">
  <a href="14_Config_Contract.html">
    <span>Understand the configuration</span>
    <strong>Learn how the four YAML roles fit together.</strong>
  </a>
  <a href="70_Case_Design_Guide.html">
    <span>Build your own case</span>
    <strong>Choose the grid, physics, boundaries, and solver.</strong>
  </a>
  <a href="03_Tutorial_File-Based_Grid.html">
    <span>Explore a curved domain</span>
    <strong>Run the file-based bent-channel tutorial.</strong>
  </a>
</nav>

<nav class="pic-page-pager" aria-label="Manual page navigation">
  <a class="pic-page-pager-next" href="06_Simulation_Anatomy.html">
    <span>Next</span>
    <strong>Simulation Anatomy →</strong>
  </a>
</nav>
@endhtmlonly
