#!/usr/bin/env python3
"""!
@file wall_normal_profile.py
@brief Reduce a field-statistics window into wall-normal profiles for DNS comparison.

@details
`field_statistics` accumulates per-cell time moments and is BC-agnostic, so it
works unchanged under periodic boundaries. What the postprocessor does not have
is a *spatial* reduction: nothing averages a statistics field over homogeneous
directions. Comparing against a channel DNS needs exactly that, so this script
does the reduction outside the postprocessor, reading the window payloads
directly.

It averages the window's mean and second-moment fields over the two homogeneous
(periodic) directions, converts to wall units using the friction velocity implied
by the driven body force, and writes a CSV of

    y, y+, U+, u'+, v'+, w'+, -<u'v'>+

together with the log-law and viscous-sublayer reference curves.

The PICGRID and PETSc-binary readers are imported from `generators/spectra.gen`
rather than duplicated; that module already owns the DMDA interior-extraction
convention (a cell-centred payload is sized `(IM+1, JM+1, KM+1)` and the physical
interior is `[1:KM, 1:JM, 1:IM]`).

Usage:

@code
    wall_normal_profile.py --checkpoint CHECKPOINT_DIR --window stationary \\
        --grid GRID.picgrid --wall-axis Eta --stream-axis Zeta \\
        --viscosity 3.5714286e-04 --output profile.csv
@endcode
"""

import argparse
import importlib.machinery
import importlib.util
import math
import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
SPECTRA_GEN = os.path.join(REPO_ROOT, "generators", "spectra.gen")

# Axis token -> index into the (k, j, i) array order the readers produce.
AXIS_TO_KJI = {"Xi": 2, "Eta": 1, "Zeta": 0}
# Axis token -> Cartesian velocity component index.
AXIS_TO_COMPONENT = {"Xi": 0, "Eta": 1, "Zeta": 2}


def load_spectra_helpers():
    """!
    @brief Import the PICGRID and PETSc-Vec readers from generators/spectra.gen.
    @return The loaded module.
    """
    if not os.path.isfile(SPECTRA_GEN):
        raise SystemExit(f"cannot find {SPECTRA_GEN}; run this from a PICurv checkout.")
    loader = importlib.machinery.SourceFileLoader("picurv_spectra_gen", SPECTRA_GEN)
    spec = importlib.util.spec_from_loader("picurv_spectra_gen", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def find_payload(checkpoint_dir, window, field, moment, block):
    """!
    @brief Locate one statistics-window payload inside a checkpoint bundle.
    @param[in] checkpoint_dir Checkpoint directory holding the window payloads.
    @param[in] window Window name from monitor.yml.
    @param[in] field Field name, e.g. Ucat.
    @param[in] moment Moment name, e.g. first or second.
    @param[in] block Block index.
    @return Absolute path to the payload file.
    """
    wanted = (window.lower(), field.lower(), moment.lower())
    candidates = []
    for root, _dirs, files in os.walk(checkpoint_dir):
        for name in files:
            if not name.endswith(".dat"):
                continue
            haystack = os.path.join(root, name).lower()
            if all(token in haystack for token in wanted) and f"block_{block:04d}" in haystack:
                candidates.append(os.path.join(root, name))
    if not candidates:
        raise SystemExit(
            f"no payload found for window={window!r} field={field!r} moment={moment!r} "
            f"block={block} under {checkpoint_dir}.\n"
            "Check that monitor.yml requested this moment and that the window reached "
            "'active' or 'complete' before the checkpoint was written."
        )
    if len(candidates) > 1:
        raise SystemExit("ambiguous payloads:\n  " + "\n  ".join(sorted(candidates)))
    return candidates[0]


def homogeneous_average(numpy, field, wall_kji):
    """!
    @brief Average a cell array over both directions that are not the wall-normal one.
    @param[in] numpy The numpy module.
    @param[in] field Cell array shaped (nk, nj, ni) or (nk, nj, ni, 3).
    @param[in] wall_kji Index of the wall-normal axis in (k, j, i) order.
    @return Array indexed by the wall-normal coordinate, keeping any trailing component axis.
    """
    homogeneous = tuple(axis for axis in (0, 1, 2) if axis != wall_kji)
    return numpy.mean(field, axis=homogeneous)


def main(argv=None):
    """!
    @brief Entry point.
    @param[in] argv Command-line style argument list supplied to the function.
    @return Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--checkpoint", required=True,
                        help="Checkpoint directory containing the statistics window payloads.")
    parser.add_argument("--grid", required=True, help="Canonical PICGRID path for the run.")
    parser.add_argument("--window", default="stationary", help="Statistics window name.")
    parser.add_argument("--block", type=int, default=0, help="Block index.")
    parser.add_argument("--wall-axis", default="Eta", choices=sorted(AXIS_TO_KJI),
                        help="Wall-normal axis token.")
    parser.add_argument("--stream-axis", default="Zeta", choices=sorted(AXIS_TO_KJI),
                        help="Streamwise (driven) axis token.")
    parser.add_argument("--viscosity", type=float, required=True,
                        help="Kinematic viscosity in solver units (1/Re for length_ref=velocity_ref=1).")
    parser.add_argument("--half-height", type=float, default=1.0,
                        help="Channel half-height h in solver units.")
    parser.add_argument("--body-force", type=float,
                        help="Converged driving body force f. If given, u_tau = sqrt(f*h) is used "
                             "instead of the wall-shear estimate, which is the exact mean force balance.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    args = parser.parse_args(argv)

    helpers = load_spectra_helpers()
    numpy = helpers.require_numpy()

    blocks = helpers.read_picgrid_blocks(args.grid)
    if args.block >= len(blocks):
        raise SystemExit(f"block {args.block} is out of range; the grid holds {len(blocks)}.")
    node_dims = blocks[args.block]["dims"]
    nodes = blocks[args.block]["coords"]

    wall_kji = AXIS_TO_KJI[args.wall_axis]
    stream_c = AXIS_TO_COMPONENT[args.stream_axis]
    wall_c = AXIS_TO_COMPONENT[args.wall_axis]
    span_c = ({0, 1, 2} - {stream_c, wall_c}).pop()

    mean_path = find_payload(args.checkpoint, args.window, "Ucat", "first", args.block)
    second_path = find_payload(args.checkpoint, args.window, "Ucat", "second", args.block)
    mean = helpers.extract_interior_cells(helpers.read_petsc_vec_binary(mean_path), node_dims)
    second = helpers.extract_interior_cells(helpers.read_petsc_vec_binary(second_path), node_dims)

    # <u_i> and <u_i u_j> both reduce over the homogeneous directions; the
    # fluctuation intensities then follow as <u_i^2> - <u_i>^2. Averaging the
    # moments before differencing (rather than after) is what makes this a true
    # ensemble average over the homogeneous plane as well as over time.
    mean_p = homogeneous_average(numpy, mean, wall_kji)
    second_p = homogeneous_average(numpy, second, wall_kji)
    variance = numpy.maximum(second_p - mean_p ** 2, 0.0)

    # Wall-normal cell-centre coordinates, taken from the node coordinates of the
    # wall-normal axis. nodes is (KM, JM, IM, 3) in the same (k, j, i) order.
    axis_nodes = {0: nodes[:, 0, 0, :], 1: nodes[0, :, 0, :], 2: nodes[0, 0, :, :]}[wall_kji]
    axis_coord = axis_nodes[:, wall_c]
    y_centres = 0.5 * (axis_coord[:-1] + axis_coord[1:])
    if y_centres.size != mean_p.shape[0]:
        raise SystemExit(f"wall-normal cell count mismatch: grid gives {y_centres.size}, "
                         f"statistics give {mean_p.shape[0]}.")

    nu, h = args.viscosity, args.half_height
    if args.body_force is not None:
        u_tau = math.sqrt(args.body_force * h)
        source = f"sqrt(f*h) with f={args.body_force:.8e}"
    else:
        # Wall-shear estimate from the first cell: du/dy at the wall.
        dudy = abs(mean_p[0, stream_c]) / abs(y_centres[0] - axis_coord[0])
        u_tau = math.sqrt(nu * dudy)
        source = "sqrt(nu*du/dy) from the first cell (approximate; prefer --body-force)"
    if u_tau <= 0.0:
        raise SystemExit("computed a non-positive friction velocity; check the inputs.")

    wall_distance = numpy.minimum(y_centres - axis_coord[0], axis_coord[-1] - y_centres)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(f"# u_tau = {u_tau:.10e}  ({source})\n")
        handle.write(f"# Re_tau = u_tau*h/nu = {u_tau * h / nu:.6f}\n")
        handle.write(f"# nu = {nu:.10e}, h = {h:.10e}, window = {args.window}\n")
        handle.write("y,y_plus,U_plus,u_rms_plus,v_rms_plus,w_rms_plus,"
                     "uv_plus,U_plus_loglaw,U_plus_sublayer\n")
        for index in range(y_centres.size):
            y_plus = wall_distance[index] * u_tau / nu
            loglaw = (1.0 / 0.41) * math.log(y_plus) + 5.2 if y_plus > 0.0 else float("nan")
            handle.write(
                f"{y_centres[index]:.10e},{y_plus:.10e},"
                f"{mean_p[index, stream_c] / u_tau:.10e},"
                f"{math.sqrt(variance[index, stream_c]) / u_tau:.10e},"
                f"{math.sqrt(variance[index, wall_c]) / u_tau:.10e},"
                f"{math.sqrt(variance[index, span_c]) / u_tau:.10e},"
                f"{(second_p[index, stream_c] - mean_p[index, stream_c] ** 2) / (u_tau ** 2):.10e},"
                f"{loglaw:.10e},{y_plus:.10e}\n")

    print(f"u_tau = {u_tau:.10e}  ({source})")
    print(f"Re_tau = {u_tau * h / nu:.4f}")
    print(f"wrote {y_centres.size} wall-normal stations to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
