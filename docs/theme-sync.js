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
    banner.appendChild(document.createTextNode(
      revision.clean ? "Documentation is up to date through commit " :
        "Documentation build includes uncommitted changes after commit "
    ));
    var link = document.createElement("a");
    link.href = revision.commit_url;
    link.textContent = revision.short_sha || revision.sha;
    link.title = revision.sha;
    banner.appendChild(link);
    banner.appendChild(document.createTextNode(revision.clean ? "." : "; it is not commit-certified."));
    var footer = document.querySelector("#nav-path ul") || document.querySelector("footer") || document.body;
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
