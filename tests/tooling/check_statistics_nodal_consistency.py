#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""!
@file check_statistics_nodal_consistency.py
@brief Cross-check derived statistics VTK output against the convergence-history CSV.

The two are produced by paths that share no code: the CSV mean comes from
`PicurvWindowSpatialMean`, which reduces the cell-centred accumulators over the
resolved spatial target, while the VTK field comes from `ComputeNodalAverage`
interpolating the derived staging buffer onto grid nodes. Agreement between them
is therefore a real check rather than a restatement.

It catches the specific failure this check was written for: an interior-only
producer leaving structural zeros on the layout boundary, which halves every
boundary node and biases a whole-domain mean without any other symptom.

@code
check_statistics_nodal_consistency.py RUN_DIR WINDOW_NAME
@endcode
"""

#: Usage line reported when the arguments do not parse.
USAGE = "usage: check_statistics_nodal_consistency.py RUN_DIR WINDOW_NAME"

import csv
import os
import re
import sys


def sole_recipe_dir(run_dir: str, canonical_relative: str) -> str:
    """!
    @brief Resolve the one per-recipe subdirectory beneath a canonical output home.
    @param[in] run_dir Run directory being inspected.
    @param[in] canonical_relative Run-relative canonical home, such as output/visualization.
    @return Absolute path to the single recipe subdirectory.
    @throws ValueError when the home is absent or holds anything but one recipe.
    """
    root = os.path.join(run_dir, canonical_relative)
    if not os.path.isdir(root):
        raise ValueError(f"{canonical_relative} does not exist under {run_dir}.")
    entries = sorted(name for name in os.listdir(root)
                     if os.path.isdir(os.path.join(root, name)))
    if len(entries) != 1:
        raise ValueError(
            f"expected exactly one recipe directory under {root}, found {entries}."
        )
    return os.path.join(root, entries[0])


def require_numpy():
    """!
    @brief Import NumPy, which this checker needs to read the binary VTK payload.
    @return Imported NumPy module.
    """
    import numpy
    return numpy


def read_vts_point_arrays(path):
    """!
    @brief Read appended raw Float64 point arrays from a PICurv `.vts`.
    @param[in] path VTK structured-grid file path.
    @return Tuple of the node extent per direction and a name-to-array mapping.
    """
    numpy = require_numpy()
    raw = open(path, "rb").read()
    head = raw[:8192].decode("utf-8", "replace")
    extent = re.search(r'WholeExtent="0 (\d+) 0 (\d+) 0 (\d+)"', head)
    if not extent:
        raise ValueError(f"{path}: no WholeExtent in header.")
    nodes = tuple(int(extent.group(i)) + 1 for i in (1, 2, 3))
    marker = raw.find(b"<AppendedData")
    if marker < 0:
        raise ValueError(f"{path}: no appended data section.")
    start = raw.find(b"_", marker) + 1
    arrays = {}
    for match in re.finditer(
        r'Name="([^"]+)" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"', head
    ):
        name, ncomp, offset = match.group(1), int(match.group(2)), int(match.group(3))
        base = start + offset
        nbytes = int(numpy.frombuffer(raw, dtype="<u4", count=1, offset=base)[0])
        values = numpy.frombuffer(raw, dtype="<f8", count=nbytes // 8, offset=base + 4)
        arrays[name] = values.reshape(-1, ncomp) if ncomp > 1 else values
    return nodes, arrays


def read_checkpoint_periodicity(run_dir):
    """!
    @brief Read the periodicity the solver recorded in any committed checkpoint.
    @param[in] run_dir Run directory to search.
    @return Tuple of three booleans for i, j, k, or None when no checkpoint is present.
    """
    for root, _dirs, files in os.walk(os.path.join(run_dir, "output")):
        if "checkpoint.meta" not in files:
            continue
        with open(os.path.join(root, "checkpoint.meta"), "r", encoding="utf-8") as stream:
            for line in stream:
                if line.startswith("-checkpoint_periodic "):
                    flags = line.split(None, 1)[1].strip().split(",")
                    return tuple(token.strip() == "1" for token in flags)
    return None


def newest_window_vts(viz_dir, window):
    """!
    @brief Locate the highest-step derived statistics file for one window.
    @param[in] viz_dir Directory holding the post-processing output.
    @param[in] window Window name.
    @return Path to the selected file.
    """
    pattern = re.compile(rf"_statistics_{re.escape(window)}_(\d+)\.vts$")
    best, best_step = None, -1
    for name in sorted(os.listdir(viz_dir)):
        match = pattern.search(name)
        if match and int(match.group(1)) > best_step:
            best, best_step = os.path.join(viz_dir, name), int(match.group(1))
    if best is None:
        raise ValueError(f"no derived statistics VTK output for window '{window}' in {viz_dir}.")
    return best


def check_periodic_wrap(arrays, nodes, periodic):
    """!
    @brief Verify every derived array wraps across each periodic layout boundary.

    @details Independent of whether the flow is turbulent, so a steady case still guards
             the layout boundary. An interior-only producer leaves the two boundary
             planes holding different partial averages, which this detects directly.

    @param[in] arrays Name-to-array mapping from `read_vts_point_arrays()`.
    @param[in] nodes Node extent per direction, in i, j, k order.
    @param[in] periodic Per-direction periodicity in i, j, k order.
    @return Zero when every array wraps, or one on the first disagreement.
    """
    numpy = require_numpy()
    checked = 0
    for name, values in sorted(arrays.items()):
        if name == "Position":
            continue
        ncomp = values.shape[1] if values.ndim > 1 else 1
        field = values.reshape(nodes[2], nodes[1], nodes[0], ncomp)
        # checkpoint_periodic is recorded in i, j, k order; the array is k, j, i.
        for flag, axis, label in zip(periodic, (2, 1, 0), ("i", "j", "k")):
            if not flag:
                continue
            low = numpy.take(field, 0, axis=axis)
            high = numpy.take(field, field.shape[axis] - 1, axis=axis)
            if not numpy.allclose(low, high):
                print(f"[FAIL] '{name}': the two layout boundary planes in the periodic "
                      f"{label} direction differ, so the wrap was not applied. See "
                      f"ExtendToLayoutBoundary().", file=sys.stderr)
                return 1
            checked += 1
    print(f"[INFO] periodic layout boundary wraps verified on {checked} array/direction pair(s).")
    return 0


def main(argv=None):
    """!
    @brief Run the consistency check.
    @param[in] argv Optional argument override for tests.
    @return Process exit code.
    """
    numpy = require_numpy()
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2:
        print(USAGE, file=sys.stderr)
        return 2
    run_dir, window = args
    # The two artifacts no longer share a directory: derived VTK lands in the run's
    # visualization home and the convergence history in its statistics home, each
    # under the producing recipe's own identity.
    viz_dir = sole_recipe_dir(run_dir, os.path.join("output", "visualization"))
    stats_dir = sole_recipe_dir(run_dir, os.path.join("output", "analysis", "statistics"))

    try:
        vts_path = newest_window_vts(viz_dir, window)
        nodes, arrays = read_vts_point_arrays(vts_path)

        periodic = read_checkpoint_periodicity(run_dir)
        if periodic and any(periodic):
            if check_periodic_wrap(arrays, nodes, periodic) != 0:
                return 1

        tke_name = next((n for n in arrays if n.endswith("_tke")), None)
        if tke_name is None:
            print(f"[SKIP] {os.path.basename(vts_path)} carries no turbulent kinetic energy field.")
            return 0
        field = arrays[tke_name].reshape(nodes[2], nodes[1], nodes[0])

        csv_path = next(
            (os.path.join(stats_dir, n) for n in sorted(os.listdir(stats_dir))
             if n.endswith(f"_statistics_{window}.csv")), None
        )
        if csv_path is None:
            raise ValueError(f"no convergence-history CSV for window '{window}' in {stats_dir}.")
        rows = list(csv.DictReader(open(csv_path, "r", encoding="utf-8")))
        if not rows or "mean_tke" not in rows[0]:
            print(f"[SKIP] {os.path.basename(csv_path)} records no mean_tke column.")
            return 0
        reference = float(rows[-1]["mean_tke"])
        if reference <= 0.0:
            print(f"[SKIP] window '{window}' has no accumulated energy to compare.")
            return 0

        # Interior nodes only: the outermost node layer is a periodic duplicate of the
        # opposite face, so counting it weights that plane twice.
        interior = field[1:-1, 1:-1, 1:-1].mean()
        deviation = abs(interior / reference - 1.0)
        tolerance = 0.05
        print(f"[INFO] {window}: nodal interior mean {interior:.6e} vs CSV mean_tke "
              f"{reference:.6e} ({100*deviation:+.2f}%)")
        if deviation > tolerance:
            print(
                f"[FAIL] derived statistics VTK and convergence CSV disagree by "
                f"{100*deviation:.2f}%, above the {100*tolerance:.0f}% interpolation "
                f"allowance. A layout boundary left unwritten by an interior-only "
                f"producer is the usual cause; see ExtendToLayoutBoundary().",
                file=sys.stderr,
            )
            return 1

        if periodic and all(periodic):
            # On a fully periodic layout every boundary node is defined, so the whole
            # domain must also agree. This is the assertion that catches an unwritten
            # layout boundary: the interior stays correct either way, so an
            # interior-only comparison would pass straight through the defect.
            whole = field.mean()
            whole_deviation = abs(whole / reference - 1.0)
            print(f"[INFO] {window}: whole-domain nodal mean {whole:.6e} "
                  f"({100*whole_deviation:+.2f}%)")
            if whole_deviation > tolerance:
                print(
                    f"[FAIL] every boundary node is defined on a fully periodic layout, "
                    f"but the whole-domain nodal mean is off by {100*whole_deviation:.2f}%. "
                    f"The layout boundary was left unwritten; see ExtendToLayoutBoundary().",
                    file=sys.stderr,
                )
                return 1

        print("[PASS] derived statistics VTK is consistent with the convergence CSV.")
        return 0
    except Exception as exc:  # noqa: BLE001 - a check reports rather than raises
        print(f"[ERROR] statistics nodal consistency check failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
