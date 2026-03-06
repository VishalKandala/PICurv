#!/usr/bin/env python3
"""Generate gcov line-coverage summary for src/*.c and enforce a threshold."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


LINE_RE = re.compile(r"^\s*([^:]+):\s*([0-9]+):(.*)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src-dir", default="src", help="directory containing C source files")
    parser.add_argument("--obj-dir", default="obj", help="directory containing coverage objects (*.gcda/*.gcno)")
    parser.add_argument("--output-dir", default="coverage/c", help="directory for gcov artifacts")
    parser.add_argument("--min-line", type=float, default=55.0, help="minimum required weighted line coverage percent")
    return parser.parse_args()


def run_gcov(src_files: list[Path], obj_dir: Path, repo_root: Path, output_dir: Path) -> None:
    for src in src_files:
        cmd = ["gcov", "-o", str(obj_dir), str(src)]
        proc = subprocess.run(cmd, cwd=str(repo_root), text=True, capture_output=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(
                f"gcov failed for {src}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )
    for gcov_file in repo_root.glob("*.gcov"):
        gcov_file.replace(output_dir / gcov_file.name)


def parse_gcov_file(gcov_path: Path, repo_root: Path) -> tuple[Path | None, int, int]:
    source_path = None
    covered = 0
    total = 0

    for raw_line in gcov_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Source:" in raw_line and ":Source:" in raw_line:
            source_path = raw_line.split("Source:", 1)[1].strip()
            continue

        match = LINE_RE.match(raw_line)
        if not match:
            continue
        count_token = match.group(1).strip()
        if count_token == "-":
            continue

        if count_token.startswith("#####") or count_token.startswith("====="):
            total += 1
            continue

        numeric = "".join(ch for ch in count_token if ch.isdigit())
        if not numeric:
            continue

        total += 1
        if int(numeric) > 0:
            covered += 1

    if not source_path:
        return None, covered, total

    source = Path(source_path)
    if not source.is_absolute():
        source = (repo_root / source).resolve()
    else:
        source = source.resolve()
    return source, covered, total


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    src_dir = (repo_root / args.src_dir).resolve()
    obj_dir = (repo_root / args.obj_dir).resolve()
    output_dir = (repo_root / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not src_dir.is_dir():
        raise SystemExit(f"[coverage-c] src directory missing: {src_dir}")
    if not obj_dir.is_dir():
        raise SystemExit(f"[coverage-c] obj directory missing: {obj_dir}")

    for old in output_dir.glob("*.gcov"):
        old.unlink()

    src_files = sorted(src_dir.glob("*.c"))
    if not src_files:
        raise SystemExit(f"[coverage-c] no source files found in {src_dir}")

    run_gcov(src_files, obj_dir, repo_root, output_dir)

    by_source: dict[Path, tuple[int, int]] = {}
    for gcov_path in sorted(output_dir.glob("*.gcov")):
        source, covered, total = parse_gcov_file(gcov_path, repo_root)
        if source is None:
            continue
        if source.parent != src_dir:
            continue
        prev_cov, prev_total = by_source.get(source, (0, 0))
        by_source[source] = (prev_cov + covered, prev_total + total)

    print("[coverage-c] per-file line coverage")
    print("[coverage-c] -----------------------------------------------")

    total_cov = 0
    total_exec = 0
    missing = []
    for src in src_files:
        covered, executable = by_source.get(src, (0, 0))
        if executable == 0:
            missing.append(src)
            percent = 0.0
        else:
            percent = (100.0 * covered) / executable
        total_cov += covered
        total_exec += executable
        rel = src.relative_to(repo_root)
        print(f"[coverage-c] {rel}: {covered}/{executable} ({percent:.2f}%)")

    overall = 0.0 if total_exec == 0 else (100.0 * total_cov) / total_exec
    print("[coverage-c] -----------------------------------------------")
    print(f"[coverage-c] weighted total: {total_cov}/{total_exec} ({overall:.2f}%)")
    print(f"[coverage-c] minimum required: {args.min_line:.2f}%")

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

    if missing:
        print("[coverage-c] WARNING: missing or zero-coverage gcov data for:", file=sys.stderr)
        for src in missing:
            print(f"[coverage-c]   - {src.relative_to(repo_root)}", file=sys.stderr)

    if overall < args.min_line:
        print(
            f"[coverage-c] FAIL: coverage {overall:.2f}% is below required {args.min_line:.2f}%.",
            file=sys.stderr,
        )
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
