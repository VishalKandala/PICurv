#!/usr/bin/env python3
"""Enforce that every eligible C function in src/ is profiled and closes its timer on every return.

Rules (docs/pages/09_Monitor_Reference.md, profiling):
- every function defined in src/*.c that is not `static inline` (or `static PICURV_ALWAYS_INLINE`) and is not on the exclusion
  list below begins exactly one timer, PROFILE_FUNCTION_BEGIN or PROFILE_HELPER_BEGIN;
- every explicit `return` and `PetscFunctionReturn`/`PetscFunctionReturnVoid` after the timer's
  BEGIN is immediately preceded by PROFILE_FUNCTION_END or PROFILE_HELPER_END, and a function whose
  body can fall off its end ends with one; error returns inside SETERRQ/PetscCall/CHKERRQ are
  not explicit returns and are exempt;
- no file defines or reads `__FUNCT__`: timers, logging and PETSc's error traces all take the
  name from `__func__`.

A missed END leaves the function's PETSc event open, after which its calls are no longer recorded.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCES = sorted((REPO_ROOT / "src").glob("*.c"))

#: Functions deliberately untimed, with the reason. Keep reasons specific.
EXCLUDED = {
    "main": "process entry; runs before PetscInitialize and after PetscFinalize",
    "RuntimeShutdownSignalHandler": "asynchronous signal handler; PETSc calls are not async-signal-safe",
    "_ResolveLogLevel": "first read of the level by every logging macro; a timer would tax every log site",
    "is_function_allowed": "read by every logging macro that passes the level check",
    "set_allowed_functions": "configures logging before the profiler exists",
    "_ProfilingRegister": "the profiler itself",
    "_ProfilingEnsureStarted": "the profiler itself",
    "_ProfilingNameListed": "the profiler itself",
    "_ProfilingHelperEnabled": "the profiler itself",
    "_ProfilingIsPicurvClass": "the profiler itself",
    "_ProfilingEventTotals": "the profiler itself",
    "_ProfilingSnapshotStep": "the profiler itself",
    "_ProfilingReduceSection": "the profiler itself",
    "_CompareNames": "qsort comparator inside the profiler's report",
    "_CompareProfiledFunctionRanks": "qsort comparator inside the profiler's report",
    "ProfilingInitialize": "the profiler itself",
    "ProfilingResetTimestepCounters": "the profiler itself",
    "ProfilingLogTimestepSummary": "the profiler itself",
    "ProfilingFinalize": "the profiler itself; it reads every event's totals",
}

#: Arithmetic kernels left untimed because they run per particle or per cell: a timer, even off,
#: stops the compiler inlining them, which cost 10-15% of particle step time (Grace, optimized
#: PETSc, 6bd5b49 vs 151e81a). Rates are measured calls per step, the largest over a particle,
#: an implicit-momentum flow and a Newton-Krylov run with every helper enabled. Their time
#: appears in their callers' timers.
KERNELS = {
    "ApplyLESTestFilter": "29/cell per step",
    "ApplyLESTestFilterSymTensor": "4.0/cell per step",
    "ApplyVolumeWeightedBoxFilter": "29/cell per step",
    "BoundaryFaceHasWallShear": "211/cell per step",
    "CalculateBrownianDisplacement": "0.9/particle per step",
    "CalculateDistancesToCellFacePlanes": "2.3/particle per step",
    "ClassifyMomentumRow": "248/cell per step",
    "ClipModelCoefficient": "1.3/cell per step",
    "ComputeCellFilterWidth": "4.7/cell per step",
    "ComputeVectorFieldDerivatives": "4.9/cell per step",
    "EddyViscosityFromCoefficient": "1.7/cell per step",
    "EdgeSeparatingFace": "11/cell per step",
    "FillEdgeCell": "5.4/cell per step",
    "GenerateGaussianNoise": "1.7/particle per step",
    "GermanoModelTensor": "1.3/cell per step",
    "GetOwnedCellRange": "0.5/particle per step",
    "Gidx": "9/cell per step",
    "LESStrainSampleIndex": "4.7/cell per step",
    "LeonardStress": "1.3/cell per step",
    "MetricGatherCellFacePlanes": "2.3/particle per step",
    "MetricGetCellVertices": "1.3/particle per step",
    "MetricGetFacePlaneArrays": "1.0/particle per step",
    "MetricRestoreFacePlaneArrays": "1.0/particle per step",
    "Normalized": "0.4/particle per step",
    "PicurvMix64": "7.7/particle per step",
    "PicurvParticleDrawHash": "0.9/particle per step",
    "PicurvUnitFromHash": "3.4/particle per step",
    "PoissonMultigrid_InterpolationParent": "57/cell per step",
    "PoissonOperator_AddFaceFlux": "2.3/cell per step",
    "PoissonOperator_StencilSlot": "23/cell per step",
    "SpatialAverageSlot": "2.7/cell per step",
    "SpatialTargetPlanMaskAllows": "7.3/cell per step",
    "StrainRateFromGradients": "4.9/cell per step",
    "SubgridKineticEnergy": "1.7/cell per step",
    "SymTensorCombine": "40/cell per step",
    "SymTensorContract": "7.6/cell per step",
    "SymTensorDeviator": "1.3/cell per step",
    "SymTensorNormSq": "6.2/cell per step",
    "SymTensorSelfOuter": "37/cell per step",
    "SymTensorTrace": "1.3/cell per step",
    "TransformDerivativesToPhysical": "4.9/cell per step",
    "TrilinearInterpolation_Scalar": "1.0/particle per step",
    "TrilinearInterpolation_Vector": "2.0/particle per step",
    "UpdateParticleField": "0.9/particle per step",
    "WallNoSlipVelocity": "43/cell per step",
}

#: Whole files deliberately untimed: dependency-free leaf modules that include only petscsys.h.
EXCLUDED_FILES = {
    "src/checksum.c": "SHA-256 and digest byte helpers; a leaf module with no PICurv headers",
    "src/statistics_moments.c": "per-sample moment arithmetic called per cell per step; a leaf module with no PICurv headers",
}

BEGIN_RE = re.compile(r"\bPROFILE_(?:FUNCTION|HELPER)_BEGIN\b")
END_TAIL_RE = re.compile(r"PROFILE_(?:FUNCTION|HELPER)_END\s*;\s*$")
RETURN_RE = re.compile(r"\b(?:return\b|PetscFunctionReturn(?:Void)?\s*\()")
FUNCT_RE = re.compile(r"\b__FUNCT__\b")
CONTROL_WORDS = {"if", "for", "while", "switch", "return", "sizeof", "do", "else"}


def blank_comments_and_strings(text: str) -> str:
    """!
    @brief Replace comments, string and character literals with spaces, keeping offsets and newlines.
    @param[in] text C source.
    @return Text of the same length in which only code remains.
    """
    out = list(text)
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            for k in range(i, j):
                out[k] = " "
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            for k in range(i, j):
                if out[k] != "\n":
                    out[k] = " "
            i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            for k in range(i + 1, min(j, n)):
                if out[k] != "\n":
                    out[k] = " "
            i = j + 1
        elif c == "#" and (i == 0 or text[i - 1] == "\n"):
            j = i
            while True:  # Preprocessor line, with continuations
                e = text.find("\n", j)
                e = n if e < 0 else e
                if e > 0 and text[e - 1] == "\\":
                    j = e + 1
                    continue
                break
            for k in range(i, e):
                out[k] = " "
            i = e
        else:
            i += 1
    return "".join(out)


def find_functions(code: str):
    """!
    @brief Locate top-level function definitions in comment- and string-blanked C code.
    @param[in] code Output of blank_comments_and_strings().
    @return List of (name, signature_start, body_open, body_close, is_static_inline).
    """
    functions = []
    depth, i, n = 0, 0, len(code)
    stmt_start = 0
    while i < n:
        c = code[i]
        if c == "{":
            if depth == 0:
                head = code[stmt_start:i]
                m = re.search(r"\b([A-Za-z_]\w*)\s*\(([^;{}]*)\)\s*$", head, re.S)
                if m and m.group(1) not in CONTROL_WORDS and not re.search(r"[=]\s*$", head.split("(")[0]):
                    j, d = i + 1, 1
                    while j < n and d:
                        d += (code[j] == "{") - (code[j] == "}")
                        j += 1
                    sig = head[: m.start()]
                    words = set(re.findall(r"\b\w+\b", sig))
                    inline = "inline" in words or "PICURV_ALWAYS_INLINE" in words
                    functions.append((m.group(1), stmt_start, i, j - 1, "static" in words and inline))
                    i, stmt_start = j, j
                    continue
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                stmt_start = i + 1
        elif c == ";" and depth == 0:
            stmt_start = i + 1
        i += 1
    return functions


def check_function(code: str, body_open: int, body_close: int):
    """!
    @brief Check one function body against the timer rules.
    @param[in] code Blanked source.
    @param[in] body_open Offset of the body's opening brace.
    @param[in] body_close Offset of the body's closing brace.
    @return List of problem strings with source offsets.
    """
    body = code[body_open + 1 : body_close]
    problems = []
    begins = len(BEGIN_RE.findall(body))
    if begins != 1:
        problems.append((body_open, f"has {begins} timer BEGIN macros, expected 1"))
        return problems
    timer_start = BEGIN_RE.search(body).start()
    for m in RETURN_RE.finditer(body):
        if m.start() < timer_start:
            continue  # An early return before the timer starts closes nothing
        before = body[: m.start()].rstrip()
        if not END_TAIL_RE.search(before):
            problems.append((body_open + 1 + m.start(), "return not preceded by PROFILE_*_END"))
    if not ends_in_return(body) and not END_TAIL_RE.search(body.rstrip()):
        problems.append((body_close, "body can fall off its end without PROFILE_*_END"))
    return problems


def ends_in_return(body: str) -> bool:
    """!
    @brief Whether a function body's last top-level statement is a return.
    @details A body ending in a block (`}`) or in any other statement can fall off its end.
    @param[in] body Blanked function body, without its braces.
    @return True when the final statement returns.
    """
    tail = body.rstrip()
    if not tail.endswith(";"):
        return False
    depth, start = 0, 0
    for i, c in enumerate(tail[:-1]):
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
            if c == "}" and depth == 0:
                start = i + 1
        elif c == ";" and depth == 0:
            start = i + 1
    return bool(RETURN_RE.match(tail[start:].lstrip()))


def line_of(text: str, offset: int) -> int:
    """!
    @brief 1-based line number of an offset.
    @param[in] text Source text.
    @param[in] offset Character offset.
    @return Line number.
    """
    return text.count("\n", 0, offset) + 1


def main() -> int:
    """!
    @brief Run the audit over src/*.c, include/*.h and tests/c, print problems.
    @return 0 when clean, 1 otherwise.
    """
    problems = []
    seen = set()
    for path in SOURCES:
        text = path.read_text(encoding="utf-8")
        code = blank_comments_and_strings(text)
        rel = path.relative_to(REPO_ROOT)
        if str(rel) in EXCLUDED_FILES:
            continue
        for name, _, body_open, body_close, static_inline in find_functions(code):
            seen.add(name)
            if name in KERNELS and not static_inline and BEGIN_RE.search(code[body_open + 1 : body_close]):
                problems.append(f"{rel}:{line_of(text, body_open)}: {name}: is a declared per-particle/per-cell kernel "
                                "but begins a timer; remove the timer or drop it from KERNELS")
                continue
            if static_inline or name in EXCLUDED or name in KERNELS:
                continue
            for offset, msg in check_function(code, body_open, body_close):
                problems.append(f"{rel}:{line_of(text, offset)}: {name}: {msg}")
    for path in sorted((REPO_ROOT / "src").glob("*.c")) + sorted((REPO_ROOT / "include").glob("*.h")) + sorted((REPO_ROOT / "tests" / "c").glob("*.[ch]")):
        text = path.read_text(encoding="utf-8")
        for m in FUNCT_RE.finditer(blank_comments_and_strings(text)):
            problems.append(f"{path.relative_to(REPO_ROOT)}:{line_of(text, m.start())}: uses __FUNCT__; use __func__")
    stale = sorted((set(EXCLUDED) | set(KERNELS)) - seen)
    for name in stale:
        problems.append(f"tests/tooling/audit_profiling_coverage.py: excluded function '{name}' no longer exists")
    for p in problems:
        print(p)
    if problems:
        print(f"\nFound {len(problems)} profiling coverage issue(s).", file=sys.stderr)
        return 1
    print("Profiling coverage audit passed: every eligible src/ function times itself and closes its timer.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
