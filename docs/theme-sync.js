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

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", wireProjectTitleHomeLink);
  } else {
    wireProjectTitleHomeLink();
  }
})();
