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
directly from one committed checkpoint bundle.

The window is looked up by name in the bundle's `checkpoint.meta`, which also
supplies its sample count and effective time bounds; its payloads are
`statistics/window_NNNN/block_NNNN/{Ucat_mean,Ucat_m2,weight}.dat`. `Ucat_m2`
holds the six centred, weighted sums M2 in the order (xx, xy, xz, yy, yz, zz),
so the per-cell temporal covariance is M2/W. The covariance over the
homogeneous plane adds the spatial covariance of the per-cell time means:

    <u_a' u_b'> = mean(M2_ab / W) + mean(U_a U_b) - mean(U_a) mean(U_b)

where mean() runs over the two directions that are not wall-normal. The two
walls are then folded onto one half-channel (the Reynolds shear stress changes
sign under the reflection), and the result is written as a CSV of

    y, y+, U+, u'+, v'+, w'+, -<u'v'>+

together with the log-law and viscous-sublayer reference curves.

The friction velocity comes from one of:

- `--body-force f`: u_tau = sqrt(f h), the exact mean force balance of a
  body-force-driven channel;
- `--wall-model-csv`: the run's `wall_model.csv`, as sqrt(<u_tau^2>) over the
  rows inside the window's effective bounds, which is the mean wall shear a
  wall-modelled run actually applies (a wall-modelled LES resolves no wall
  gradient to fit);
- `--u-tau`: a value supplied directly;
- otherwise sqrt(nu dU/dy) from the first cell, which is only meaningful when the
  first cell sits in the viscous sublayer.

The PICGRID and PETSc-binary readers are imported from `generators/spectra.gen`,
which owns the DMDA interior-extraction convention (a cell-centred payload is
sized `(IM+1, JM+1, KM+1)` and the physical interior is `[1:KM, 1:JM, 1:IM]`), and
the metadata parser from `picurv_cli`, rather than duplicated.

Usage:

@code
    wall_normal_profile.py \\
        --checkpoint RUN/output/checkpoints/step_000000050000 --window stationary \\
        --grid RUN/inputs/grid/grid.run --wall-axis Eta --stream-axis Zeta \\
        --viscosity 5.0e-05 --wall-model-csv RUN/output/analysis/metrics/wall_model.csv \\
        --output profile.csv
@endcode
"""

import argparse
import csv
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
# Component pairs of a dof-6 second-moment payload, in its stored order.
M2_PAIRS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


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


def read_window(checkpoint_dir, window, block):
    """!
    @brief Resolve a statistics window by name inside one committed checkpoint bundle.
    @param[in] checkpoint_dir Checkpoint bundle directory holding `checkpoint.meta`.
    @param[in] window Window name from monitor.yml.
    @param[in] block Block index.
    @return Tuple of (payload directory, window metadata dict, block node dims or None).
    """
    sys.path.insert(0, REPO_ROOT)
    from picurv_cli.core import _read_checkpoint_options

    metadata = os.path.join(checkpoint_dir, "checkpoint.meta")
    if not os.path.isfile(metadata):
        raise SystemExit(f"{checkpoint_dir} is not a committed checkpoint bundle: no checkpoint.meta.")
    options = _read_checkpoint_options(metadata)
    names = []
    for index in range(int(options.get("checkpoint_statistics_window_count", 0) or 0)):
        prefix = f"checkpoint_statistics_window_{index}_"
        names.append(options.get(prefix + "name"))
        if names[-1] != window:
            continue
        info = {key[len(prefix):]: value for key, value in options.items() if key.startswith(prefix)}
        payload_dir = os.path.join(checkpoint_dir, "statistics", f"window_{index:04d}", f"block_{block:04d}")
        if not os.path.isdir(payload_dir):
            raise SystemExit(f"window {window!r} has no payloads for block {block} under {payload_dir}.")
        dims = None
        if f"checkpoint_block_{block}_im" in options:
            dims = tuple(int(options[f"checkpoint_block_{block}_{axis}"]) for axis in ("im", "jm", "km"))
        return payload_dir, info, dims
    raise SystemExit(f"no statistics window named {window!r} in {metadata}; it holds {names}.")


def friction_velocity_from_wall_model(csv_path, start, end):
    """!
    @brief Mean wall shear the wall model applied over a time interval, as a velocity.
    @details Each wall_model.csv row carries the wall-face mean and standard deviation
             of u_tau, so the row's mean of u_tau^2 is mean^2 + rms^2.
    @param[in] csv_path Path to the run's wall_model.csv.
    @param[in] start Window effective start, in solver time.
    @param[in] end Window effective end, in solver time.
    @return Tuple of (sqrt(<u_tau^2>), number of rows used).
    """
    total, rows = 0.0, 0
    with open(csv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if start <= float(row["time"]) <= end:
                total += float(row["u_tau_mean"]) ** 2 + float(row["u_tau_rms"]) ** 2
                rows += 1
    if rows == 0:
        raise SystemExit(f"{csv_path} has no rows with time in the window [{start}, {end}].")
    return math.sqrt(total / rows), rows


def main(argv=None):
    """!
    @brief Entry point.
    @param[in] argv Command-line style argument list supplied to the function.
    @return Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--checkpoint", required=True,
                        help="Committed checkpoint bundle directory (output/checkpoints/step_NNN).")
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
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--body-force", type=float,
                        help="Converged driving body force f; u_tau = sqrt(f*h), the exact mean "
                             "force balance.")
    source.add_argument("--wall-model-csv",
                        help="The run's wall_model.csv; u_tau = sqrt(<u_tau^2>) over the window.")
    source.add_argument("--u-tau", type=float, help="Friction velocity supplied directly.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    args = parser.parse_args(argv)

    helpers = load_spectra_helpers()
    numpy = helpers.require_numpy()

    payload_dir, window, meta_dims = read_window(args.checkpoint, args.window, args.block)
    blocks = helpers.read_picgrid_blocks(args.grid)
    if args.block >= len(blocks):
        raise SystemExit(f"block {args.block} is out of range; the grid holds {len(blocks)}.")
    node_dims = blocks[args.block]["dims"]
    nodes = blocks[args.block]["coords"]
    if meta_dims is not None and tuple(meta_dims) != tuple(node_dims):
        raise SystemExit(f"grid dimensions {node_dims} do not match the checkpoint's {meta_dims}.")

    wall_kji = AXIS_TO_KJI[args.wall_axis]
    stream_c = AXIS_TO_COMPONENT[args.stream_axis]
    wall_c = AXIS_TO_COMPONENT[args.wall_axis]
    span_c = ({0, 1, 2} - {stream_c, wall_c}).pop()
    homogeneous = tuple(axis for axis in (0, 1, 2) if axis != wall_kji)

    def payload(name, components):
        values = helpers.read_petsc_vec_binary(os.path.join(payload_dir, f"{name}.dat"))
        return helpers.extract_interior_cells(values, node_dims, components)

    mean = payload("Ucat_mean", 3)
    m2 = payload("Ucat_m2", 6)
    weight = payload("weight", 1)[..., 0]
    if not numpy.all(weight > 0.0):
        raise SystemExit(f"window {args.window!r} has cells with no accepted weight.")

    plane_mean = numpy.mean(mean, axis=homogeneous)
    temporal = numpy.mean(m2 / weight[..., None], axis=homogeneous)

    def covariance(a, b):
        pair = M2_PAIRS.index(tuple(sorted((a, b))))
        spatial = numpy.mean(mean[..., a] * mean[..., b], axis=homogeneous) \
            - plane_mean[:, a] * plane_mean[:, b]
        return temporal[:, pair] + spatial

    # Wall-normal cell-centre coordinates, taken from the node coordinates of the
    # wall-normal axis. nodes is (KM, JM, IM, 3) in the same (k, j, i) order.
    axis_nodes = {0: nodes[:, 0, 0, :], 1: nodes[0, :, 0, :], 2: nodes[0, 0, :, :]}[wall_kji]
    axis_coord = axis_nodes[:, wall_c]
    y_centres = 0.5 * (axis_coord[:-1] + axis_coord[1:])
    if y_centres.size != plane_mean.shape[0]:
        raise SystemExit(f"wall-normal cell count mismatch: grid gives {y_centres.size}, "
                         f"statistics give {plane_mean.shape[0]}.")

    # Fold the upper half onto the lower one. Reflection reverses the wall-normal
    # velocity, so <u'v'> changes sign; the orientation term makes the reported
    # value the shear stress relative to the wall the point is measured from.
    orientation = 1.0 if axis_coord[-1] > axis_coord[0] else -1.0
    half = (y_centres.size + 1) // 2

    def fold(values, sign=1.0):
        return 0.5 * (values[:half] + sign * values[::-1][:half])

    wall_distance = fold(numpy.minimum(numpy.abs(y_centres - axis_coord[0]),
                                       numpy.abs(axis_coord[-1] - y_centres)))
    U = fold(plane_mean[:, stream_c])
    uu = fold(covariance(stream_c, stream_c))
    vv = fold(covariance(wall_c, wall_c))
    ww = fold(covariance(span_c, span_c))
    uv = orientation * fold(covariance(stream_c, wall_c), -1.0)

    nu, h = args.viscosity, args.half_height
    if args.body_force is not None:
        u_tau = math.sqrt(args.body_force * h)
        source_note = f"sqrt(f*h) with f={args.body_force:.8e}"
    elif args.wall_model_csv is not None:
        start, end = float(window["effective_start"]), float(window["effective_end"])
        u_tau, rows = friction_velocity_from_wall_model(args.wall_model_csv, start, end)
        source_note = f"sqrt(<u_tau^2>) over {rows} wall_model.csv rows, t in [{start:g}, {end:g}]"
    elif args.u_tau is not None:
        u_tau = args.u_tau
        source_note = "supplied with --u-tau"
    else:
        u_tau = math.sqrt(nu * abs(U[0]) / wall_distance[0])
        source_note = "sqrt(nu*U/y) from the first cell (approximate; needs a resolved sublayer)"
    if not u_tau > 0.0:
        raise SystemExit("computed a non-positive friction velocity; check the inputs.")

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(f"# u_tau = {u_tau:.10e}  ({source_note})\n")
        handle.write(f"# Re_tau = u_tau*h/nu = {u_tau * h / nu:.6f}\n")
        handle.write(f"# nu = {nu:.10e}, h = {h:.10e}, window = {args.window} "
                     f"({window.get('state')}, {window.get('sample_count')} samples, "
                     f"represented time {float(window.get('represented_time', 'nan')):g})\n")
        handle.write("y,y_plus,U_plus,u_rms_plus,v_rms_plus,w_rms_plus,"
                     "minus_uv_plus,U_plus_loglaw,U_plus_sublayer\n")
        for index in range(half):
            y_plus = wall_distance[index] * u_tau / nu
            loglaw = (1.0 / 0.41) * math.log(y_plus) + 5.2 if y_plus > 0.0 else float("nan")
            handle.write(
                f"{wall_distance[index]:.10e},{y_plus:.10e},"
                f"{U[index] / u_tau:.10e},"
                f"{math.sqrt(max(uu[index], 0.0)) / u_tau:.10e},"
                f"{math.sqrt(max(vv[index], 0.0)) / u_tau:.10e},"
                f"{math.sqrt(max(ww[index], 0.0)) / u_tau:.10e},"
                f"{-uv[index] / (u_tau ** 2):.10e},"
                f"{loglaw:.10e},{y_plus:.10e}\n")

    print(f"u_tau = {u_tau:.10e}  ({source_note})")
    print(f"Re_tau = {u_tau * h / nu:.4f}")
    print(f"wrote {half} wall-normal stations to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
