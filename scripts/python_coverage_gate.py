#!/usr/bin/env python3
"""Run pytest under stdlib trace and enforce a line-coverage threshold."""

from __future__ import annotations

import argparse
import os
import sys
import trace
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGETS = [
    "scripts/picurv",
]


def normalize_path(path: str | Path) -> str:
    """!
    @brief Normalize path.
    @param[in] path Filesystem path argument passed to `normalize_path()`.
    @return Value returned by `normalize_path()`.
    """
    return str(Path(path).resolve())


def parse_args() -> argparse.Namespace:
    """!
    @brief Parse args.
    @return Value returned by `parse_args()`.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/python_coverage_gate.py\n"
            "  python3 scripts/python_coverage_gate.py --target scripts/picurv --target scripts/grid.gen\n"
            "  python3 scripts/python_coverage_gate.py --pytest-args -- -q tests/test_cli_smoke.py\n"
        ),
    )
    parser.add_argument(
        "--min-line",
        type=float,
        default=70.0,
        help="Minimum required weighted line coverage percent (default: 70.0).",
    )
    parser.add_argument(
        "--target",
        action="append",
        default=[],
        help=(
            "Repository-relative file to include in coverage computation "
            "(repeatable). Defaults to core runtime scripts."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="coverage/python",
        help="Repository-relative output directory for coverage artifacts (default: coverage/python).",
    )
    parser.add_argument(
        "--pytest-args",
        nargs=argparse.REMAINDER,
        default=["-q"],
        help="Arguments passed to pytest (prefix with --pytest-args -- ...).",
    )
    return parser.parse_args()


def build_trace_ignoredirs() -> list[str]:
    """!
    @brief Build trace ignoredirs.
    @return Value returned by `build_trace_ignoredirs()`.
    """
    ignoredirs = {
        normalize_path(sys.prefix),
        normalize_path(sys.exec_prefix),
        normalize_path(Path(sys.prefix) / "lib"),
    }
    for raw in list(sys.path):
        if not raw:
            continue
        path = Path(raw)
        if not path.exists():
            continue
        resolved = normalize_path(path)
        if "/site-packages" in resolved or "/dist-packages" in resolved:
            ignoredirs.add(resolved)
    return sorted(ignoredirs)


def collect_counts(results: trace.CoverageResults) -> dict[str, dict[int, int]]:
    """!
    @brief Collect counts.
    @param[in] results Argument passed to `collect_counts()`.
    @return Value returned by `collect_counts()`.
    """
    counts_by_file: dict[str, dict[int, int]] = {}
    for (filename, lineno), count in results.counts.items():
        file_key = normalize_path(filename)
        counts_by_file.setdefault(file_key, {})[lineno] = count
    return counts_by_file


def compute_file_coverage(target: Path, counts_by_file: dict[str, dict[int, int]]) -> tuple[int, int, float]:
    """!
    @brief Compute file coverage.
    @param[in] target Argument passed to `compute_file_coverage()`.
    @param[in] counts_by_file Argument passed to `compute_file_coverage()`.
    @return Value returned by `compute_file_coverage()`.
    """
    finder = getattr(trace, "find_executable_linenos", None)
    if finder is None:
        finder = trace._find_executable_linenos  # type: ignore[attr-defined]
    executable = finder(str(target))
    executable_lines = set(executable.keys())
    total = len(executable_lines)
    if total == 0:
        return 0, 0, 100.0

    observed = counts_by_file.get(normalize_path(target), {})
    covered = sum(1 for line in executable_lines if observed.get(line, 0) > 0)
    percent = (100.0 * covered) / total
    return covered, total, percent


def main() -> int:
    """!
    @brief Entry point for this script.
    @return Value returned by `main()`.
    """
    args = parse_args()
    output_dir = (REPO_ROOT / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    targets_raw = args.target if args.target else DEFAULT_TARGETS
    targets = [Path(REPO_ROOT / rel).resolve() for rel in targets_raw]
    for target in targets:
        if not target.exists():
            raise SystemExit(f"[coverage-python] target not found: {target}")

    pytest_args = list(args.pytest_args)
    if pytest_args and pytest_args[0] == "--":
        pytest_args = pytest_args[1:]
    if not pytest_args:
        pytest_args = ["-q"]

    ignoredirs = build_trace_ignoredirs()

    import pytest

    tracer = trace.Trace(count=True, trace=False, ignoredirs=ignoredirs)
    exit_code = tracer.runfunc(pytest.main, pytest_args)
    results = tracer.results()

    counts_by_file = collect_counts(results)

    print("[coverage-python] per-file line coverage")
    print("[coverage-python] -----------------------------------------------")

    total_cov = 0
    total_exec = 0
    for target in targets:
        covered, executable, percent = compute_file_coverage(target, counts_by_file)
        total_cov += covered
        total_exec += executable
        rel = target.relative_to(REPO_ROOT)
        print(f"[coverage-python] {rel}: {covered}/{executable} ({percent:.2f}%)")

    overall = 100.0 if total_exec == 0 else (100.0 * total_cov) / total_exec
    print("[coverage-python] -----------------------------------------------")
    print(f"[coverage-python] weighted total: {total_cov}/{total_exec} ({overall:.2f}%)")
    print(f"[coverage-python] minimum required: {args.min_line:.2f}%")

    summary_path = output_dir / "summary.txt"
    summary_path.write_text(
        "\n".join(
            [
                f"weighted_total={overall:.4f}",
                f"covered_lines={total_cov}",
                f"executable_lines={total_exec}",
                f"minimum_required={args.min_line:.4f}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    if int(exit_code) != 0:
        print(f"[coverage-python] pytest failed with exit code {exit_code}.", file=sys.stderr)
        return int(exit_code)
    if overall < args.min_line:
        print(
            f"[coverage-python] FAIL: coverage {overall:.2f}% is below required {args.min_line:.2f}%.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
