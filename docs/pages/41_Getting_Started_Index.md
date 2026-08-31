@page 41_Getting_Started_Index Quick Start

@anchor _Getting_Started_Index
@anchor p41_quicklook_sec

@pagemeta{Tutorial, New users, Current workflow}

@htmlonly
<header class="pic-quickstart-intro">
  <span class="pic-eyebrow">Manual · Quick Start</span>
  <p class="pic-quickstart-lead">Run a complete PICurv simulation, post-process it,
  and inspect the resulting velocity field.</p>
  <ul class="pic-quickstart-facts" aria-label="Quick Start facts">
    <li>20 timesteps</li>
    <li>2 MPI ranks</li>
    <li>Programmatic grid</li>
  </ul>
</header>
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

@section p41_understand_sec How the run flows

PICurv validates the four profiles and materializes a self-contained run:

@htmlonly
<ol class="pic-run-pipeline" aria-label="Quick Start execution pipeline">
  <li>YAML profiles</li>
  <li>Validated controls</li>
  <li>MPI solver</li>
  <li>Checkpoint</li>
  <li>Postprocessor</li>
  <li>VTK output</li>
</ol>
@endhtmlonly

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
