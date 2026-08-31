@mainpage PICurv Solver Documentation

@htmlonly
<section class="pic-hero" aria-labelledby="pic-hero-title">
  <div class="pic-hero-copy">
    <div class="pic-hero-wordmark" aria-label="PICurv">
      <img data-pic-brand-asset src="picurv-wordmark.svg" alt="" />
    </div>
    <h1 id="pic-hero-title">
    Fast,Modern and robust C/C++ software for simulating incompressible fluid dynamics and particle transport in complex, curved domains.
    </h1>
  </div>
  <div class="pic-hero-logo" role="img" aria-label="PICurv graphical logo">
    <span class="pic-logo-fallback" aria-hidden="true">P</span>
    <img data-pic-brand-asset src="picurv-mark-dark.svg" alt="" />
  </div>
</section>

<section class="pic-highlights" aria-label="PICurv highlights">
  <article class="pic-highlight-card">
    <img class="pic-highlight-visual" src="pic-facecard.svg" alt="" aria-hidden="true" />
    <span class="pic-highlight-index">1</span>
    <div class="pic-highlight-copy">
      <h2>Coupled</h2>
      <p>Interpolate and project fields between curvilinear grids and tracking particles.</p>
    </div>
  </article>
  <article class="pic-highlight-card">
    <img class="pic-highlight-visual" src="config-facecard.svg" alt="" aria-hidden="true" />
    <span class="pic-highlight-index">2</span>
    <div class="pic-highlight-copy">
      <h2>Configurable</h2>
      <p>Compose and swap physics,solvers, data acquisition and analysis pipelines.</p>
    </div>
  </article>
  <article class="pic-highlight-card">
    <img class="pic-highlight-visual" src="hpc-facecard.svg" alt="" aria-hidden="true" />
    <span class="pic-highlight-index">3</span>
    <div class="pic-highlight-copy">
      <h2>Scalable</h2>
      <p>Scale and resolve simulations by deploying massive grids and millions of particles.</p>
    </div>
  </article>
</section>

<nav class="pic-home-actions" aria-label="PICurv documentation shortcuts">
  <a class="pic-button pic-button-primary" href="01_Installation.html">Install</a>
  <a class="pic-button pic-button-secondary" href="41_Getting_Started_Index.html#p41_quicklook_sec">Quick Start</a>
</nav>

<div class="pic-home-knowledge">
  <section class="pic-lower-section pic-doc-atlas" aria-labelledby="pic-knowledge-title">
    <header class="pic-section-heading">
      <span class="pic-section-kicker">Explore the documentation</span>
      <h2 id="pic-knowledge-title">The PICurv Knowledge Base</h2>
      <p>Deep, transparent documentation for building cases, understanding the numerics, and extending the solver.</p>
    </header>

    <div class="pic-atlas-grid">
      <a class="pic-atlas-card pic-atlas-card--wide" href="41_Getting_Started_Index.html">
        <span class="pic-atlas-icon" aria-hidden="true"><svg viewBox="0 0 24 24"><path d="M5 19V5h9l5 5v9H5Z M14 5v5h5 M8 14h8 M8 17h6" /></svg></span>
        <span class="pic-atlas-copy"><strong>Step-by-Step Tutorials</strong><span>Installation, structured-grid setup, curvilinear grid ingestion, execution, and visualization.</span></span>
        <span class="pic-card-arrow" aria-hidden="true">→</span>
      </a>
      <a class="pic-atlas-card" href="49_Workflow_Recipes_and_Config_Cookbook.html">
        <span class="pic-atlas-icon" aria-hidden="true"><svg viewBox="0 0 24 24"><path d="M6 4h12v16H6z M9 8h6 M9 12h6 M9 16h4" /></svg></span>
        <span class="pic-atlas-copy"><strong>YAML Recipes</strong><span>Reusable <code>case.yml</code>, <code>solver.yml</code>, <code>monitor.yml</code>, and <code>post.yml</code> profiles.</span></span>
        <span class="pic-card-arrow" aria-hidden="true">→</span>
      </a>
      <a class="pic-atlas-card" href="21_Methods_Overview.html">
        <span class="pic-atlas-icon" aria-hidden="true"><svg viewBox="0 0 24 24"><path d="M4 18 9 6l4 12 3-8 4 8 M4 21h16" /></svg></span>
        <span class="pic-atlas-copy"><strong>Methods and Numerics</strong><span>Eulerian-Lagrangian finite-volume methods and their runtime implementations.</span></span>
        <span class="pic-card-arrow" aria-hidden="true">→</span>
      </a>
      <a class="pic-atlas-card" href="23_Fractional_Step_Method.html">
        <span class="pic-atlas-icon" aria-hidden="true"><svg viewBox="0 0 24 24"><circle cx="7" cy="8" r="2"/><circle cx="17" cy="7" r="2"/><circle cx="14" cy="17" r="2"/><path d="m9 8 6-1 M8 10l5 5 M16 9l-2 6" /></svg></span>
        <span class="pic-atlas-copy"><strong>Equations and Particle Models</strong><span>Incompressible flow, projection, stochastic transport, and particle coupling.</span></span>
        <span class="pic-card-arrow" aria-hidden="true">→</span>
      </a>
      <a class="pic-atlas-card" href="annotated_structured.html">
        <span class="pic-atlas-icon pic-atlas-icon--code" aria-hidden="true">{ }</span>
        <span class="pic-atlas-copy"><strong>C/C++ API Reference</strong><span>Doxygen-generated structures, persistent fields, descriptors, and core routines.</span></span>
        <span class="pic-card-arrow" aria-hidden="true">→</span>
      </a>
    </div>
    <nav class="pic-related-links" aria-label="Related technical references">
      <span>Related:</span>
      <a href="34_Particle_Model_Overview.html">Particle coupling pipeline</a>
      <a href="56_Field_Identity_and_Layout_Catalog.html">Field and descriptor catalog</a>
    </nav>
  </section>

  <section class="pic-lower-section pic-open-development" aria-labelledby="pic-github-title">
    <header class="pic-section-heading">
      <span class="pic-section-kicker">Built in the open</span>
      <h2 id="pic-github-title">Open Development</h2>
      <p>Follow solver work, propose changes, and inspect the verification gates behind PICurv.</p>
    </header>
    <div class="pic-development-grid">
      <article class="pic-github-feature">
        <svg class="pic-github-mark" viewBox="0 0 24 24" aria-hidden="true"><path fill="currentColor" d="M12 .7A11.3 11.3 0 0 0 8.4 22.8c.6.1.8-.2.8-.6v-2.1c-3.3.7-4-1.4-4-1.4-.5-1.3-1.3-1.7-1.3-1.7-1.1-.7.1-.7.1-.7 1.2.1 1.8 1.2 1.8 1.2 1.1 1.8 2.8 1.3 3.5 1 .1-.8.4-1.3.8-1.6-2.7-.3-5.5-1.3-5.5-5.8 0-1.3.5-2.3 1.2-3.1-.1-.3-.5-1.6.1-3.1 0 0 1-.3 3.3 1.2a11.3 11.3 0 0 1 6 0c2.3-1.5 3.3-1.2 3.3-1.2.6 1.5.2 2.8.1 3.1.8.8 1.2 1.8 1.2 3.1 0 4.5-2.8 5.5-5.5 5.8.4.4.8 1.1.8 2.1v3.2c0 .4.2.7.8.6A11.3 11.3 0 0 0 12 .7Z"/></svg>
        <span class="pic-github-label">VishalKandala / PICurv</span>
        <h3>Solver development, documentation, and research roadmap</h3>
        <p>Track the code and the evidence used to support its documented behavior.</p>
        <a class="pic-inline-button" href="https://github.com/VishalKandala/PICurv">View repository <span aria-hidden="true">↗</span></a>
      </article>
      <nav class="pic-development-links" aria-label="PICurv development resources">
        <a href="https://github.com/VishalKandala/PICurv/issues"><span><strong>Issues and questions</strong><small>Bugs, numerical questions, missing documentation, custom grids, and Slurm workflows.</small></span><b aria-hidden="true">→</b></a>
        <a href="https://github.com/VishalKandala/PICurv/pulls"><span><strong>Pull requests</strong><small>Concrete solver, verification, and documentation contributions.</small></span><b aria-hidden="true">→</b></a>
        <a href="57_Future_Architecture_Specifications.html"><span><strong>Future architecture</strong><small>Explicitly labeled experimental and longer-term solver directions.</small></span><b aria-hidden="true">→</b></a>
        <a href="40_Testing_and_Quality_Guide.html"><span><strong>Testing and quality</strong><small>Automated gates used to verify proposed changes.</small></span><b aria-hidden="true">→</b></a>
      </nav>
    </div>
    <p class="pic-status-line"><span aria-hidden="true"></span>GitHub Discussions is not currently enabled; use Issues for scientific and workflow questions.</p>
  </section>

  <section class="pic-lower-section pic-contact-section" aria-labelledby="pic-contact-title">
    <header class="pic-section-heading">
      <span class="pic-section-kicker">Choose the right route</span>
      <h2 id="pic-contact-title">Getting in Touch</h2>
      <p>Ask about a first curvilinear simulation, report a boundary-condition defect, or start a research collaboration.</p>
    </header>
    <div class="pic-contact-grid">
      <article class="pic-contact-card">
        <span class="pic-contact-index">01</span><h3>Ask or report</h3>
        <p>For bugs, numerical questions, missing documentation, and reproducible feature requests.</p>
        <a href="https://github.com/VishalKandala/PICurv/issues">Open an issue <span aria-hidden="true">→</span></a>
      </article>
      <article class="pic-contact-card">
        <span class="pic-contact-index">02</span><h3>Contribute</h3>
        <p>For a concrete code, verification, workflow, or documentation contribution.</p>
        <a href="https://github.com/VishalKandala/PICurv/pulls">Open a pull request <span aria-hidden="true">→</span></a>
      </article>
      <article class="pic-contact-card">
        <span class="pic-contact-index">03</span><h3>Collaborate</h3>
        <p>Connect with the Scientific Computing and Biofluids Lab at Texas A&amp;M.</p>
        <div class="pic-contact-links"><a href="https://scbl.engr.tamu.edu/">Visit SCBL <span aria-hidden="true">↗</span></a><a href="mailto:vishalkandala@tamu.edu">Email directly</a></div>
      </article>
    </div>
  </section>

  <section class="pic-lower-section pic-citing-section" aria-labelledby="pic-citing-title">
    <header class="pic-section-heading">
      <span class="pic-section-kicker">Research record</span>
      <h2 id="pic-citing-title">Presentations and Related Research</h2>
      <p>Selected presentations describing the PICurv solver and related high-resolution biofluid applications.</p>
    </header>
    <div class="pic-citation-list">
      <a class="pic-citation-card" href="https://schedule.aps.org/dfd/2025/events/T21/4">
        <span class="pic-citation-year">2025</span>
        <span class="pic-citation-copy"><small>Primary solver presentation</small><strong>A Parallel Particle-In-Cell (PIC) Solver on Curvilinear Grids for Turbulent Flow Simulation</strong><span>V. I. Kandala and I. Borazjani</span><em>78th Annual Meeting of the APS Division of Fluid Dynamics · Houston, Texas</em></span>
        <span class="pic-citation-arrow" aria-hidden="true">↗</span>
      </a>
      <a class="pic-citation-card" href="https://meetings.aps.org/Meeting/DFD23/Session/L06.9">
        <span class="pic-citation-year">2023</span>
        <span class="pic-citation-copy"><small>Related biofluid application</small><strong>High Resolution Numerical Simulations of LVAD Outflow Graft Haemodynamics</strong><span>V. I. Kandala, M. Kassi, M. Moreno, and I. Borazjani</span><em>76th Annual Meeting of the APS Division of Fluid Dynamics · Washington, DC</em></span>
        <span class="pic-citation-arrow" aria-hidden="true">↗</span>
      </a>
    </div>
  </section>

  <section class="pic-project-note" aria-labelledby="pic-final-card-title">
    <div class="pic-project-note-copy">
      <span class="pic-section-kicker">Cite PICurv</span>
      <h2 id="pic-final-card-title">Using PICurv in your work?</h2>
      <p class="pic-project-note-intro">Please cite the primary solver presentation when PICurv supports your research, teaching, or modeling.</p>
      <div class="pic-primary-citation">
        <p id="pic-primary-citation-text">Kandala, V. I., and Borazjani, I. (2025). “A Parallel Particle-In-Cell (PIC) Solver on Curvilinear Grids for Turbulent Flow Simulation.” Presented at the 78th Annual Meeting of the APS Division of Fluid Dynamics, Houston, Texas.</p>
        <button type="button" data-pic-copy-target="pic-primary-citation-text">Copy citation</button>
      </div>
      <p class="pic-project-credit">PICurv is developed by <a href="https://vishalkandala.me">Vishal Kandala</a>.</p>
      <div class="pic-project-note-links"><a href="https://schedule.aps.org/dfd/2025/events/T21/4">View presentation <span aria-hidden="true">↗</span></a><a href="https://vishalkandala.me">Developer website <span aria-hidden="true">↗</span></a><a href="mailto:vishalkandala@tamu.edu">vishalkandala@tamu.edu</a></div>
    </div>
  </section>
</div>
@endhtmlonly
