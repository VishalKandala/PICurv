// Sync Doxygen light/dark theme with the host site's stored preference.
(function () {
  function resolveTheme() {
    try {
      var storedTheme = window.localStorage.getItem("theme");
      if (storedTheme === "light" || storedTheme === "dark") return storedTheme;
    } catch (error) {
      // Storage can be unavailable in privacy-restricted browser contexts.
    }
    if (window.matchMedia && window.matchMedia("(prefers-color-scheme: light)").matches) {
      return "light";
    }
    return "dark";
  }

  function syncThemeControls(theme) {
    var controls = document.querySelectorAll("[data-pic-theme-toggle]");
    var nextTheme = theme === "light" ? "dark" : "light";
    for (var i = 0; i < controls.length; i += 1) {
      controls[i].setAttribute("aria-label", "Switch to " + nextTheme + " mode");
      controls[i].title = "Switch to " + nextTheme + " mode";
    }
  }

  function applyTheme(theme) {
    document.documentElement.setAttribute("data-theme", theme);
    syncThemeControls(theme);
  }

  applyTheme(resolveTheme());

  window.addEventListener("storage", function (event) {
    if (event.key === "theme") applyTheme(resolveTheme());
  });

  function addBrandAsset(container, className, fallbackClass, fallbackText, source) {
    var frame = document.createElement("span");
    frame.className = className;

    var fallback = document.createElement("span");
    fallback.className = fallbackClass;
    fallback.setAttribute("aria-hidden", "true");
    fallback.textContent = fallbackText;
    frame.appendChild(fallback);

    var image = document.createElement("img");
    image.src = source;
    image.alt = "";
    image.setAttribute("data-pic-brand-asset", "");
    frame.appendChild(image);
    container.appendChild(frame);
  }

  function addMenuGroup(container, label, links) {
    var details = document.createElement("details");
    details.className = "pic-topbar-menu";

    var summary = document.createElement("summary");
    summary.textContent = label;
    details.appendChild(summary);

    var panel = document.createElement("div");
    panel.className = "pic-topbar-dropdown";
    for (var i = 0; i < links.length; i += 1) {
      var link = document.createElement("a");
      link.href = links[i][1];
      link.textContent = links[i][0];
      panel.appendChild(link);
    }
    details.appendChild(panel);
    container.appendChild(details);
    return details;
  }

  // Doxygen's own index is disabled. Build one small, independent top bar so a
  // navigation failure can never turn the former sidebar into page content.
  function addSiteTopbar() {
    if (document.getElementById("pic-site-topbar")) return;

    var header = document.createElement("header");
    header.id = "pic-site-topbar";
    header.className = "pic-site-topbar";

    var inner = document.createElement("nav");
    inner.className = "pic-site-topbar-inner";
    inner.setAttribute("aria-label", "Primary documentation navigation");

    var brand = document.createElement("a");
    brand.id = "pic-site-brand";
    brand.className = "pic-site-brand";
    brand.href = "index.html";
    brand.setAttribute("aria-label", "PICurv documentation home");
    addBrandAsset(brand, "pic-site-banner", "pic-brand-banner-fallback", "PICurv", "picurv-banner.svg");
    inner.appendChild(brand);

    var menus = document.createElement("div");
    menus.className = "pic-topbar-menus";
    var groups = [
      ["Manual", [
        ["Quick Start", "41_Getting_Started_Index.html#p41_quicklook_sec"],
        ["Installation", "01_Installation.html"],
        ["User Guide", "42_User_Guide_Index.html"],
        ["Configuration Reference", "14_Config_Contract.html"],
        ["Case Design", "70_Case_Design_Guide.html"]
      ]],
      ["Workflows", [
        ["Examples", "65_Example_Catalog.html"],
        ["Recipes", "49_Workflow_Recipes_and_Config_Cookbook.html"],
        ["Run Lifecycle", "52_Run_Artifact_Lifecycle_Contract.html"],
        ["Cluster Runs", "36_Cluster_Run_Guide.html"],
        ["Sweep Studies", "37_Sweep_Studies_Guide.html"],
        ["Visualization", "04_Visualization_Tutorial.html"]
      ]],
      ["Concepts", [
        ["Methods and Models", "21_Methods_Overview.html"],
        ["Simulation Anatomy", "06_Simulation_Anatomy.html"],
        ["Curvilinear Formulation", "22_CURVIB_Method.html"],
        ["Particle Coupling", "34_Particle_Model_Overview.html"],
        ["Boundary Conditions", "44_Boundary_Conditions_Guide.html"]
      ]],
      ["Resources", [
        ["Capabilities", "12_Capabilities_Summary.html"],
        ["Troubleshooting", "67_Troubleshooting.html"],
        ["Glossary", "68_Glossary.html"],
        ["Documentation Map", "Documentation_Map.html"],
        ["Evidence Matrix", "66_Evidence_Matrix.html"],
        ["C API", "annotated_structured.html"]
      ]],
      ["Developer", [
        ["Developer Portal", "43_Developer_Portal_Index.html"],
        ["Code Architecture", "13_Code_Architecture.html"],
        ["Runtime Execution Map", "46_C_Runtime_Execution_Map.html"],
        ["Extension Playbook", "16_Config_Extension_Playbook.html"],
        ["Testing and Quality", "40_Testing_and_Quality_Guide.html"]
      ]]
    ];
    var dropdowns = [];
    for (var i = 0; i < groups.length; i += 1) {
      dropdowns.push(addMenuGroup(menus, groups[i][0], groups[i][1]));
    }
    inner.appendChild(menus);

    var search = document.getElementById("MSearchBox");
    if (search) {
      var searchWrap = document.createElement("div");
      searchWrap.className = "pic-topbar-search";
      searchWrap.setAttribute("role", "search");
      var searchField = search.querySelector("#MSearchField");
      if (searchField) searchField.setAttribute("aria-label", "Search PICurv documentation");
      searchWrap.appendChild(search);
      inner.appendChild(searchWrap);
    }

    var themeToggle = document.createElement("button");
    themeToggle.className = "pic-theme-toggle";
    themeToggle.type = "button";
    themeToggle.setAttribute("data-pic-theme-toggle", "");
    var themeSvg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    themeSvg.setAttribute("viewBox", "0 0 24 24");
    themeSvg.setAttribute("aria-hidden", "true");
    var sunGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
    sunGroup.setAttribute("class", "pic-theme-sun");
    sunGroup.setAttribute("fill", "none");
    sunGroup.setAttribute("stroke", "currentColor");
    sunGroup.setAttribute("stroke-width", "1.8");
    sunGroup.setAttribute("stroke-linecap", "round");
    var sunCircle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    sunCircle.setAttribute("cx", "12");
    sunCircle.setAttribute("cy", "12");
    sunCircle.setAttribute("r", "3.5");
    var sunRays = document.createElementNS("http://www.w3.org/2000/svg", "path");
    sunRays.setAttribute("d", "M12 2v2 M12 20v2 M4.93 4.93l1.42 1.42 M17.65 17.65l1.42 1.42 M2 12h2 M20 12h2 M4.93 19.07l1.42-1.42 M17.65 6.35l1.42-1.42");
    sunGroup.appendChild(sunCircle);
    sunGroup.appendChild(sunRays);
    var moonPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
    moonPath.setAttribute("class", "pic-theme-moon");
    moonPath.setAttribute("fill", "none");
    moonPath.setAttribute("stroke", "currentColor");
    moonPath.setAttribute("stroke-width", "1.8");
    moonPath.setAttribute("stroke-linecap", "round");
    moonPath.setAttribute("stroke-linejoin", "round");
    moonPath.setAttribute("d", "M20.5 15.4A8.5 8.5 0 0 1 8.6 3.5 8.5 8.5 0 1 0 20.5 15.4Z");
    themeSvg.appendChild(sunGroup);
    themeSvg.appendChild(moonPath);
    themeToggle.appendChild(themeSvg);
    themeToggle.addEventListener("click", function () {
      var nextTheme = document.documentElement.getAttribute("data-theme") === "light" ? "dark" : "light";
      try {
        window.localStorage.setItem("theme", nextTheme);
      } catch (error) {
        // The switch remains useful for the current page even without storage.
      }
      applyTheme(nextTheme);
    });
    inner.appendChild(themeToggle);
    syncThemeControls(document.documentElement.getAttribute("data-theme") || resolveTheme());

    var github = document.createElement("a");
    github.className = "pic-topbar-github";
    github.href = "https://github.com/VishalKandala/PICurv";
    github.target = "_blank";
    github.rel = "noopener noreferrer";
    github.setAttribute("aria-label", "PICurv on GitHub");
    github.title = "PICurv on GitHub";
    var githubSvg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    githubSvg.setAttribute("viewBox", "0 0 24 24");
    githubSvg.setAttribute("aria-hidden", "true");
    var githubPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
    githubPath.setAttribute("fill", "currentColor");
    githubPath.setAttribute("d", "M12 .3a12 12 0 0 0-3.8 23.4c.6.1.8-.3.8-.6v-2.3c-3.3.7-4-1.4-4-1.4-.5-1.4-1.3-1.7-1.3-1.7-1.1-.7.1-.7.1-.7 1.2.1 1.8 1.2 1.8 1.2 1.1 1.8 2.8 1.3 3.5 1 .1-.8.4-1.3.8-1.6-2.7-.3-5.5-1.3-5.5-5.9 0-1.3.5-2.4 1.2-3.2-.1-.3-.5-1.6.1-3.2 0 0 1-.3 3.3 1.2a11.5 11.5 0 0 1 6 0C17.3 5 18.3 5.3 18.3 5.3c.6 1.6.2 2.9.1 3.2.8.8 1.2 1.9 1.2 3.2 0 4.6-2.8 5.6-5.5 5.9.4.4.8 1.1.8 2.2v3.3c0 .3.2.7.8.6A12 12 0 0 0 12 .3");
    githubSvg.appendChild(githubPath);
    github.appendChild(githubSvg);
    inner.appendChild(github);

    header.appendChild(inner);
    document.body.insertBefore(header, document.body.firstChild);

    for (var j = 0; j < dropdowns.length; j += 1) {
      (function (current) {
        function supportsMenuHover() {
          return window.matchMedia && window.matchMedia("(hover: hover) and (pointer: fine)").matches;
        }

        current.addEventListener("mouseenter", function () {
          if (supportsMenuHover()) current.open = true;
        });

        current.addEventListener("mouseleave", function () {
          if (supportsMenuHover() && !current.contains(document.activeElement)) current.open = false;
        });

        current.addEventListener("toggle", function () {
          if (!current.open) return;
          for (var k = 0; k < dropdowns.length; k += 1) {
            if (dropdowns[k] !== current) dropdowns[k].open = false;
          }
        });
      })(dropdowns[j]);
    }
    document.addEventListener("click", function (event) {
      if (header.contains(event.target)) return;
      for (var k = 0; k < dropdowns.length; k += 1) dropdowns[k].open = false;
    });
    document.addEventListener("keydown", function (event) {
      if (event.key !== "Escape") return;
      for (var k = 0; k < dropdowns.length; k += 1) dropdowns[k].open = false;
    });
  }

  function prepareBrandAssetFallbacks() {
    var assets = document.querySelectorAll("img[data-pic-brand-asset]");
    for (var i = 0; i < assets.length; i += 1) {
      (function (image) {
        var fallback = image.previousElementSibling;

        function showAsset() {
          image.style.display = "";
          if (fallback) fallback.style.display = "none";
        }

        function showFallback() {
          image.style.display = "none";
          if (fallback) fallback.style.display = "";
        }

        image.addEventListener("load", showAsset);
        image.addEventListener("error", showFallback);
        if (image.complete) {
          if (image.naturalWidth) showAsset();
          else showFallback();
        }
      })(assets[i]);
    }
  }

  function currentPageName() {
    return window.location.pathname.split("/").pop() || "index.html";
  }

  function markPageContext() {
    var page = currentPageName();
    if (!page || page === "index.html") document.body.classList.add("pic-home-page");
    if (page === "41_Getting_Started_Index.html") {
      document.body.classList.add("pic-manual-page", "pic-quick-start-page");
    }
    if (page === "01_Installation.html") {
      document.body.classList.add("pic-manual-page", "pic-installation-page");
    }
  }

  function buildLinkList(className, label, links) {
    var nav = document.createElement("nav");
    nav.className = className;
    nav.setAttribute("aria-label", label);
    for (var i = 0; i < links.length; i += 1) {
      var link = document.createElement("a");
      link.href = links[i][1];
      link.textContent = links[i][0];
      if (links[i][2]) {
        link.className = "is-current";
        link.setAttribute("aria-current", "page");
      }
      nav.appendChild(link);
    }
    return nav;
  }

  function buildCollapsibleNav(className, summaryText, label, links) {
    var details = document.createElement("details");
    details.className = className;
    var summary = document.createElement("summary");
    summary.textContent = summaryText;
    details.appendChild(summary);
    details.appendChild(buildLinkList(className + "-links", label, links));
    return details;
  }

  // The redesigned entry pages share one manual shell and navigation model.
  function decorateQuickStartSections(contents) {
    var headings = Array.prototype.slice.call(contents.querySelectorAll(".textblock > h1"));
    for (var i = 0; i < headings.length; i += 1) {
      var heading = headings[i];
      var nextHeading = headings[i + 1] || null;
      var anchor = heading.querySelector("a.anchor[id]");
      var title = heading.textContent.replace(/\s+/g, " ").trim();
      var stepMatch = title.match(/^([1-4])\.\s*(.+)$/);
      var section = document.createElement("section");
      section.className = "pic-quickstart-section";
      if (anchor) section.setAttribute("data-section", anchor.id);

      if (stepMatch) {
        section.classList.add("pic-quickstart-section--step");
        section.setAttribute("data-step", ("0" + stepMatch[1]).slice(-2));
        heading.setAttribute("data-step", ("0" + stepMatch[1]).slice(-2));
        heading.setAttribute("aria-label", "Step " + stepMatch[1] + ": " + stepMatch[2]);
        for (var j = 0; j < heading.childNodes.length; j += 1) {
          if (heading.childNodes[j].nodeType === 3) {
            heading.childNodes[j].textContent = heading.childNodes[j].textContent.replace(
              /^\s*[1-4]\.\s*/, ""
            );
          }
        }
      } else {
        section.classList.add("pic-quickstart-section--support");
      }

      heading.parentNode.insertBefore(section, heading);
      var node = heading;
      while (node && node !== nextHeading) {
        var nextNode = node.nextSibling;
        section.appendChild(node);
        node = nextNode;
      }
    }
  }

  function decorateInstallationSections(contents) {
    var headings = Array.prototype.slice.call(contents.querySelectorAll(".textblock > h1"));
    for (var i = 0; i < headings.length; i += 1) {
      var heading = headings[i];
      var nextHeading = headings[i + 1] || null;
      var anchor = heading.querySelector("a.anchor[id]");
      var title = heading.textContent.replace(/\s+/g, " ").trim();
      var stepMatch = title.match(/^(\d+)\.\s*(.+)$/);
      var section = document.createElement("section");
      section.className = "pic-install-section";
      if (anchor) section.setAttribute("data-section", anchor.id);

      if (stepMatch) {
        var step = ("0" + stepMatch[1]).slice(-2);
        section.setAttribute("data-step", step);
        heading.setAttribute("aria-label", "Step " + stepMatch[1] + ": " + stepMatch[2]);
        for (var j = 0; j < heading.childNodes.length; j += 1) {
          if (heading.childNodes[j].nodeType === 3) {
            heading.childNodes[j].textContent = heading.childNodes[j].textContent.replace(
              /^\s*\d+\.\s*/, ""
            );
          }
        }
      }

      if (anchor && anchor.id === "p01_automated_sec") {
        section.classList.add("pic-install-section--featured");
      }
      if (anchor && anchor.id === "p01_verify_sec") {
        section.classList.add("pic-install-section--verify");
      }
      if (anchor && anchor.id === "p01_common_sec") {
        section.classList.add("pic-install-section--failures");
      }
      if (anchor && anchor.id === "p01_next_steps_sec") {
        section.classList.add("pic-install-section--next");
      }

      heading.parentNode.insertBefore(section, heading);
      var node = heading;
      while (node && node !== nextHeading) {
        var nextNode = node.nextSibling;
        section.appendChild(node);
        node = nextNode;
      }
    }

    var fragments = contents.querySelectorAll(".pic-install-section .fragment");
    for (var k = 0; k < fragments.length; k += 1) {
      var fragment = fragments[k];
      if (fragment.parentNode.classList.contains("pic-command-block")) continue;
      var command = document.createElement("div");
      command.className = "pic-command-block";
      command.setAttribute("data-label", "Terminal");
      fragment.parentNode.insertBefore(command, fragment);
      command.appendChild(fragment);
    }
  }

  function trackQuickStartSections() {
    if (!("IntersectionObserver" in window)) return;
    var sections = document.querySelectorAll(
      ".pic-quickstart-section[data-section], .pic-install-section[data-section]"
    );
    var links = document.querySelectorAll(".pic-onpage-links a, .pic-onpage-mobile-links a");
    if (!sections.length || !links.length) return;

    function setCurrent(sectionId) {
      for (var i = 0; i < links.length; i += 1) {
        var current = links[i].getAttribute("href") === "#" + sectionId;
        links[i].classList.toggle("is-active", current);
        if (current) links[i].setAttribute("aria-current", "location");
        else links[i].removeAttribute("aria-current");
      }
    }

    var observer = new IntersectionObserver(function (entries) {
      var visible = entries.filter(function (entry) { return entry.isIntersecting; });
      if (!visible.length) return;
      visible.sort(function (a, b) { return a.boundingClientRect.top - b.boundingClientRect.top; });
      setCurrent(visible[0].target.getAttribute("data-section"));
    }, { rootMargin: "-18% 0px -68% 0px", threshold: 0 });

    for (var i = 0; i < sections.length; i += 1) observer.observe(sections[i]);
  }

  function addQuickStartShell() {
    var page = currentPageName();
    var isQuickStart = page === "41_Getting_Started_Index.html";
    var isInstallation = page === "01_Installation.html";
    if (!isQuickStart && !isInstallation) return;
    if (document.querySelector(".pic-doc-shell")) return;

    var pageDoc = document.querySelector(".PageDoc");
    if (!pageDoc) return;
    var pageHeader = pageDoc.querySelector(":scope > .header");
    var contents = pageDoc.querySelector(":scope > .contents");
    if (!pageHeader || !contents) return;

    var manualLinks = [
      ["Quick Start", "41_Getting_Started_Index.html", isQuickStart],
      ["Installation", "01_Installation.html", isInstallation],
      ["Simulation Anatomy", "06_Simulation_Anatomy.html"],
      ["Configuration Model", "14_Config_Contract.html"],
      ["Case Design", "70_Case_Design_Guide.html"],
      ["Grids and Domains", "48_Grid_Generator_Guide.html"],
      ["Fields and Initial Conditions", "33_Initial_Conditions.html"],
      ["Boundary Conditions", "44_Boundary_Conditions_Guide.html"],
      ["Particles", "34_Particle_Model_Overview.html"],
      ["Running PICurv", "05_The_Conductor_Script.html"],
      ["Simulation Tips", "11_User_How_To_Guides.html"]
    ];

    var sectionLinks = [];
    var headings = contents.querySelectorAll(".textblock h1");
    for (var i = 0; i < headings.length; i += 1) {
      var anchor = headings[i].querySelector("a.anchor[id]");
      if (!anchor) continue;
      var title = headings[i].textContent.replace(/\s+/g, " ").trim();
      sectionLinks.push([title, "#" + anchor.id]);
      if (isQuickStart && /^[1-4]\.\s/.test(title)) headings[i].classList.add("pic-numbered-step");
    }

    if (isQuickStart) decorateQuickStartSections(contents);
    if (isInstallation) decorateInstallationSections(contents);

    var shell = document.createElement("div");
    shell.className = "pic-doc-shell";

    var manualAside = document.createElement("aside");
    manualAside.className = "pic-manual-sidebar";
    var manualTitle = document.createElement("a");
    manualTitle.className = "pic-sidebar-title";
    manualTitle.href = "42_User_Guide_Index.html";
    manualTitle.textContent = "Manual";
    manualAside.appendChild(manualTitle);
    manualAside.appendChild(buildLinkList("pic-sidebar-links", "Manual pages", manualLinks));

    var article = document.createElement("main");
    article.className = "pic-doc-article";
    article.appendChild(pageHeader);
    article.appendChild(buildCollapsibleNav(
      "pic-manual-mobile", "Manual", "Manual pages", manualLinks
    ));
    if (sectionLinks.length) {
      article.appendChild(buildCollapsibleNav(
        "pic-onpage-mobile", "On this page", "Sections on this page", sectionLinks
      ));
    }
    article.appendChild(contents);

    var onPageAside = document.createElement("aside");
    onPageAside.className = "pic-onpage-sidebar";
    var onPageTitle = document.createElement("span");
    onPageTitle.className = "pic-sidebar-title";
    onPageTitle.textContent = "On this page";
    onPageAside.appendChild(onPageTitle);
    onPageAside.appendChild(buildLinkList("pic-onpage-links", "Sections on this page", sectionLinks));

    shell.appendChild(manualAside);
    shell.appendChild(article);
    shell.appendChild(onPageAside);
    pageDoc.appendChild(shell);
    trackQuickStartSections();
  }

  function copyText(value) {
    if (navigator.clipboard && navigator.clipboard.writeText) {
      return navigator.clipboard.writeText(value);
    }
    return new Promise(function (resolve, reject) {
      var field = document.createElement("textarea");
      field.value = value;
      field.setAttribute("readonly", "");
      field.style.position = "fixed";
      field.style.opacity = "0";
      document.body.appendChild(field);
      field.select();
      try {
        document.execCommand("copy");
        resolve();
      } catch (error) {
        reject(error);
      }
      document.body.removeChild(field);
    });
  }

  function addQuickStartCopyButtons() {
    var page = currentPageName();
    if (page !== "41_Getting_Started_Index.html" && page !== "01_Installation.html") return;
    var blocks = document.querySelectorAll(".pic-command-block .fragment");
    for (var i = 0; i < blocks.length; i += 1) {
      (function (block) {
        if (block.querySelector(".pic-copy-button")) return;
        block.classList.add("pic-copyable");
        var button = document.createElement("button");
        button.type = "button";
        button.className = "pic-copy-button";
        button.textContent = "Copy";
        button.setAttribute("aria-label", "Copy command to clipboard");
        button.addEventListener("click", function () {
          var lines = block.querySelectorAll(".line");
          var text = [];
          for (var j = 0; j < lines.length; j += 1) text.push(lines[j].textContent);
          copyText(text.join("\n")).then(function () {
            button.textContent = "Copied";
            window.setTimeout(function () { button.textContent = "Copy"; }, 1400);
          });
        });
        block.appendChild(button);
      })(blocks[i]);
    }
  }

  function addCitationCopyButtons() {
    var buttons = document.querySelectorAll("button[data-pic-copy-target]");
    for (var i = 0; i < buttons.length; i += 1) {
      (function (button) {
        if (button.getAttribute("data-pic-copy-ready") === "true") return;
        button.setAttribute("data-pic-copy-ready", "true");
        button.addEventListener("click", function () {
          var target = document.getElementById(button.getAttribute("data-pic-copy-target"));
          if (!target) return;
          copyText(target.textContent.replace(/\s+/g, " ").trim()).then(function () {
            button.textContent = "Copied";
            window.setTimeout(function () { button.textContent = "Copy citation"; }, 1400);
          });
        });
      })(buttons[i]);
    }
  }

  function addRevisionBanner() {
    var revision = window.PICURV_DOCS_REVISION;
    if (!revision || !revision.sha || document.getElementById("picurv-docs-revision-banner")) return;
    var banner = document.createElement("div");
    banner.id = "picurv-docs-revision-banner";
    banner.className = "picurv-docs-revision-banner";
    // Provenance only: this banner reports what the build knows -- which release and
    // which commit the HTML was generated from.  It deliberately makes no claim about
    // whether the prose was semantically reviewed; that is tracked separately per page.
    // The label leads with the plain release ("PICurv 0.2.0") only for a build that
    // actually is that release -- a clean checkout of its tag.  Anything else states
    // the full build id, so a development or dirty build is never shown as if it were
    // the release it is only heading toward.
    var releaseLabel = revision.released ? revision.release_version : revision.build_id;
    banner.appendChild(document.createTextNode(
      (releaseLabel ? "PICurv " + releaseLabel + " — built" : "Built") +
        (revision.clean ? " from commit " : " with uncommitted changes after commit ")
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
      markPageContext();
      addSiteTopbar();
      addQuickStartShell();
      addQuickStartCopyButtons();
      addCitationCopyButtons();
      prepareBrandAssetFallbacks();
      addRevisionBanner();
    });
  } else {
    markPageContext();
    addSiteTopbar();
    addQuickStartShell();
    addQuickStartCopyButtons();
    addCitationCopyButtons();
    prepareBrandAssetFallbacks();
    addRevisionBanner();
  }
})();
