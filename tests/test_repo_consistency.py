import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PICURV = REPO_ROOT / "scripts" / "picurv"
AUDIT_FUNCTION_DOCS = REPO_ROOT / "scripts" / "audit_function_docs.py"

EXAMPLE_BUNDLES = [
    {
        "case": "examples/flat_channel/flat_channel.yml",
        "solver": "examples/flat_channel/Imp-MG-Standard.yml",
        "monitor": "examples/flat_channel/Standard_Output.yml",
        "post": "examples/flat_channel/standard_analysis.yml",
    },
    {
        "case": "examples/bent_channel/bent_channel.yml",
        "solver": "examples/bent_channel/Imp-MG-Standard.yml",
        "monitor": "examples/bent_channel/Standard_Output.yml",
        "post": "examples/bent_channel/standard_analysis.yml",
    },
    {
        "case": "examples/brownian_motion/brownian_motion.yml",
        "solver": "examples/brownian_motion/Analytical-Zero.yml",
        "monitor": "examples/brownian_motion/Standard_Output.yml",
        "post": "examples/brownian_motion/brownian_analysis.yml",
    },
    {
        "case": "examples/master_template/master_case.yml",
        "solver": "examples/master_template/master_solver.yml",
        "monitor": "examples/master_template/master_monitor.yml",
        "post": "examples/master_template/master_postprocessor.yml",
    },
]

STUDY_BUNDLES = [
    {
        "cluster": "examples/flat_channel/slurm_cluster.yml",
        "study": "examples/flat_channel/grid_independence_study.yml",
    },
    {
        "cluster": "examples/bent_channel/slurm_cluster.yml",
        "study": "examples/bent_channel/timestep_sensitivity_study.yml",
    },
    {
        "cluster": "examples/master_template/master_cluster.yml",
        "study": "examples/master_template/master_study.yml",
    },
]


def run_picurv(args):
    """!
    @brief Run picurv.
    @param[in] args Command-line style argument list supplied to the function.
    @return Value returned by `run_picurv()`.
    """
    cmd = [sys.executable, str(PICURV)] + list(args)
    return subprocess.run(cmd, cwd=str(REPO_ROOT), text=True, capture_output=True, timeout=60, check=False)


def _read_text(path: Path) -> str:
    """!
    @brief Helper for read text.
    @param[in] path Filesystem path argument passed to `_read_text()`.
    @return Value returned by `_read_text()`.
    """
    return path.read_text(encoding="utf-8")


def test_all_example_bundles_validate():
    """!
    @brief Test that all example bundles validate.
    """
    for bundle in EXAMPLE_BUNDLES:
        args = ["validate"]
        for role in ("case", "solver", "monitor", "post"):
            args.extend([f"--{role}", str(REPO_ROOT / bundle[role])])

        result = run_picurv(args)
        assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_all_example_study_and_cluster_bundles_validate():
    """!
    @brief Test that all example study and cluster bundles validate.
    """
    for bundle in STUDY_BUNDLES:
        result = run_picurv(
            [
                "validate",
                "--cluster",
                str(REPO_ROOT / bundle["cluster"]),
                "--study",
                str(REPO_ROOT / bundle["study"]),
            ]
        )
        assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_docs_and_examples_do_not_use_legacy_post_run_control_yaml_keys():
    """!
    @brief Test that docs and examples do not use legacy post run control yaml keys.
    """
    legacy_key_pattern = re.compile(r"(?m)^\s*(startTime|endTime|timeStep):")

    scanned_roots = [
        REPO_ROOT / "docs",
        REPO_ROOT / "examples",
        REPO_ROOT / "tests",
    ]
    matched = []
    for root in scanned_roots:
        for path in root.rglob("*"):
            if path.suffix not in {".md", ".yml"}:
                continue
            if legacy_key_pattern.search(_read_text(path)):
                matched.append(str(path.relative_to(REPO_ROOT)))

    assert matched == [], f"Legacy post YAML keys found in: {matched}"


def test_docs_and_examples_do_not_contain_stale_post_contract_terms():
    """!
    @brief Test that docs and examples do not contain stale post contract terms.
    """
    forbidden_literals = [
        "post.cfg",
        "output/BrownianStats_msd.csv",
        "configured output location",
    ]

    scanned_roots = [
        REPO_ROOT / "docs",
        REPO_ROOT / "examples",
        REPO_ROOT / "tests",
    ]
    hits = []
    for root in scanned_roots:
        for path in root.rglob("*"):
            if path.suffix not in {".md", ".yml"}:
                continue
            text = _read_text(path)
            for literal in forbidden_literals:
                if literal in text:
                    hits.append(f"{path.relative_to(REPO_ROOT)}::{literal}")

    assert hits == [], f"Stale contract terms found: {hits}"


def test_yaml_templates_do_not_use_unquoted_off_for_mode_keys():
    """!
    @brief Test that yaml templates do not use unquoted off for mode keys.
    """
    bare_off_pattern = re.compile(r"(?m)^\s*mode:\s+off\s*(?:#.*)?$")

    scanned_roots = [
        REPO_ROOT / "examples",
        REPO_ROOT / "tests",
    ]
    matched = []
    for root in scanned_roots:
        for path in root.rglob("*.yml"):
            if bare_off_pattern.search(_read_text(path)):
                matched.append(str(path.relative_to(REPO_ROOT)))

    assert matched == [], f"Unquoted YAML 'mode: off' found in: {matched}"


def test_c_and_python_function_docs_pass_repository_audit():
    """!
    @brief Test that C and Python function docs pass the repository audit.
    """
    result = subprocess.run(
        [sys.executable, str(AUDIT_FUNCTION_DOCS)],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
