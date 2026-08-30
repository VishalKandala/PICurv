// Sync Doxygen light/dark theme with the host site's stored preference.
(function () {
  function resolveTheme() {
    var stored = null;
    try {
      stored = window.localStorage ? window.localStorage.getItem("theme") : null;
    } catch (e) {
      stored = null;
    }
    if (stored === "dark" || stored === "light") return stored;
    if (window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches) {
      return "dark";
    }
    return "light";
  }

  function applyTheme(theme) {
    document.documentElement.setAttribute("data-theme", theme);
  }

  applyTheme(resolveTheme());

  window.addEventListener("storage", function (event) {
    if (event.key === "theme") applyTheme(resolveTheme());
  });

  // Make the project title act as "home" navigation without a dedicated navbar tab.
  function wireProjectTitleHomeLink() {
    var projectName = document.getElementById("projectname");
    if (!projectName) return;
    projectName.style.cursor = "pointer";
    projectName.addEventListener("click", function () {
      window.location.href = "index.html";
    });
  }

  function addRevisionBanner() {
    var revision = window.PICURV_DOCS_REVISION;
    if (!revision || !revision.sha || document.getElementById("picurv-docs-revision-banner")) return;
    var banner = document.createElement("div");
    banner.id = "picurv-docs-revision-banner";
    banner.className = "picurv-docs-revision-banner";
    // Provenance only: this banner reports what the build knows -- which commit the
    // HTML was generated from.  It deliberately makes no claim about whether the
    // prose was semantically reviewed; that is tracked separately per page.
    banner.appendChild(document.createTextNode(
      revision.clean ? "Built from commit " :
        "Built with uncommitted changes after commit "
    ));
    var link = document.createElement("a");
    link.href = revision.commit_url;
    link.textContent = revision.short_sha || revision.sha;
    link.title = revision.sha;
    banner.appendChild(link);
    banner.appendChild(document.createTextNode(revision.clean ? "." : "; the build is not reproducible from that commit alone."));
    // Doxygen renders its badge in `address.footer`.  Keep the revision claim
    // in that same footer so it reads as one coherent provenance row.
    var footer = document.querySelector("address.footer") || document.querySelector("footer") || document.body;
    footer.appendChild(banner);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () {
      wireProjectTitleHomeLink();
      addRevisionBanner();
    });
  } else {
    wireProjectTitleHomeLink();
    addRevisionBanner();
  }
})();
