#!/usr/bin/env python3
"""!
@file wall_normal_profile.py
@brief Reduce a field-statistics window into wall-normal profiles for DNS comparison.

@details
The reduction over the homogeneous directions is the post-processor's. A recipe
deriving the window with `field_statistics.homogeneous` set to the streamwise and
spanwise directions and `outputs` including `mean` and `reynolds_stress` writes a
profile CSV whose rows hold, per wall-normal cell, the plane mean velocity and the
covariance over the plane's pooled samples, the per-cell temporal covariance plus
the spread of the cell means. This script reads that profile back, folds the two
walls onto one half-channel (the Reynolds shear stress changes sign under the
reflection), and writes a CSV of

    y, y+, U+, u'+, v'+, w'+, -<u'v'>+

together with the log-law and viscous-sublayer reference curves. The window's
effective time bounds come from the `checkpoint.meta` of the bundle the profile was
derived from, and the wall positions from the run's grid.

The friction velocity comes from one of:

- `--body-force f`: u_tau = sqrt(f h), the exact mean force balance of a
  body-force-driven channel;
- `--driven-flow-csv`: the same balance with f the driving acceleration from the run's
  `driven_flow.csv`, averaged over the window's effective bounds;
- `--wall-model-csv`: the run's `wall_model.csv`, as sqrt(<u_tau^2>) over the
  rows inside the window's effective bounds, which is the mean wall shear a
  wall-modelled run actually applies (a wall-modelled LES resolves no wall
  gradient to fit);
- `--u-tau`: a value supplied directly;
- otherwise sqrt(nu dU/dy) from the first cell, which is only meaningful when the
  first cell sits in the viscous sublayer.

The PICGRID reader is imported from `generators/spectra.gen` and the metadata parser
from `picurv_cli`, rather than duplicated. The profile must be in solver units: a
recipe with `dimensionalize` off, or a case with unit reference scales.

Usage:

@code
    wall_normal_profile.py \\
        --profile RUN/output/analysis/statistics/RECIPE/Field_profile_stationary_50000.csv \\
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
import json
import math
import os
import sys

def find_repo_root():
    """!
    @brief Locate the PICurv checkout from the shipped example or from a case `picurv init` made.
    @details In the repository this file sits four levels below the root. A copy inside an
             initialized case does not, so the case's `.picurv-origin.json`, found by
             walking up from this file, supplies the checkout it was made from.
    @return Absolute path to the checkout.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    shipped = os.path.abspath(os.path.join(here, "..", "..", "..", ".."))
    if os.path.isfile(os.path.join(shipped, "generators", "spectra.gen")):
        return shipped
    current = here
    while True:
        origin = os.path.join(current, ".picurv-origin.json")
        if os.path.isfile(origin):
            with open(origin, encoding="utf-8") as handle:
                return json.load(handle).get("source_repo_root", shipped)
        parent = os.path.dirname(current)
        if parent == current:
            return shipped
        current = parent


REPO_ROOT = find_repo_root()
SPECTRA_GEN = os.path.join(REPO_ROOT, "generators", "spectra.gen")

# Axis token -> index into the (k, j, i) array order the readers produce.
AXIS_TO_KJI = {"Xi": 2, "Eta": 1, "Zeta": 0}
# Axis token -> Cartesian velocity component index.
AXIS_TO_COMPONENT = {"Xi": 0, "Eta": 1, "Zeta": 2}
# Component pairs of the symmetric covariance, in the stored (xx, xy, xz, yy, yz, zz) order.
M2_PAIRS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def load_spectra_helpers():
    """!
    @brief Import the PICGRID and PETSc-Vec readers from generators/spectra.gen.
    @return The loaded module.
    """
    if not os.path.isfile(SPECTRA_GEN):
        raise SystemExit(f"cannot find {SPECTRA_GEN}; run this from a PICurv checkout "
                         "or a case initialized from one.")
    loader = importlib.machinery.SourceFileLoader("picurv_spectra_gen", SPECTRA_GEN)
    spec = importlib.util.spec_from_loader("picurv_spectra_gen", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def read_window(checkpoint_dir, window):
    """!
    @brief Read a statistics window's scalars from one committed checkpoint bundle.
    @details Only `checkpoint.meta` is read: the window's state, sample count and
             effective time bounds, which the friction-velocity sources average over.
    @param[in] checkpoint_dir Checkpoint bundle directory holding `checkpoint.meta`.
    @param[in] window Window name from monitor.yml.
    @return Dict of the window's checkpoint scalars.
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
        if names[-1] == window:
            return {key[len(prefix):]: value for key, value in options.items() if key.startswith(prefix)}
    raise SystemExit(f"no statistics window named {window!r} in {metadata}; it holds {names}.")


def profile_statistics(profile_csv, window_name, grid_path, block, retained):
    """!
    @brief Read a window's velocity statistics merged over its homogeneous directions.
    @details The reduction is the post-processor's: a recipe with
             `field_statistics.homogeneous` naming the homogeneous directions and
             `outputs` including `mean` and `reynolds_stress` writes, per retained
             cell, the plane mean and the covariance over the plane's pooled samples
             (the per-cell temporal covariance plus the spread of the cell means). This
             reads that CSV back into arrays indexed by the retained cells in (k, j, i)
             order, so the folding and wall units below work on them directly. The
             values are taken as they were written, so the recipe must leave them in
             solver units: `dimensionalize` off, or unit reference scales.
    @param[in] profile_csv The `<prefix>_profile_<window>_<step>.csv` the recipe wrote.
    @param[in] window_name Window the profile was written for.
    @param[in] grid_path Canonical PICGRID path for the run.
    @param[in] block Block index.
    @param[in] retained Axes the profile must retain, as indices in (k, j, i) order.
    @return Dict with `numpy`, `nodes` ((KM, JM, IM, 3) node coordinates), `mean`
            (retained axes + (3,)) and `covariance` (retained axes + (6,), in
            M2_PAIRS order).
    """
    helpers = load_spectra_helpers()
    numpy = helpers.require_numpy()
    blocks = helpers.read_picgrid_blocks(grid_path)
    if block >= len(blocks):
        raise SystemExit(f"block {block} is out of range; the grid holds {len(blocks)}.")
    node_dims = blocks[block]["dims"]
    cells_kji = (node_dims[2] - 1, node_dims[1] - 1, node_dims[0] - 1)
    index_names = ("k", "j", "i")
    with open(profile_csv, encoding="utf-8") as handle:
        rows = [line.strip().split(",") for line in handle if line.strip()]
    header, rows = rows[0], rows[1:]
    present = tuple(axis for axis in range(3) if index_names[axis] in header)
    if present != tuple(sorted(retained)):
        raise SystemExit(f"{profile_csv} retains {[index_names[a] for a in present]}, but this "
                         f"reduction needs {[index_names[a] for a in sorted(retained)]}; set the "
                         "recipe's field_statistics.homogeneous to the other directions.")
    labels = "xyz"
    columns = [f"{window_name}_Ucat_mean_{labels[c]}" for c in range(3)]
    columns += [f"{window_name}_Ucat_R_{labels[a]}{labels[b]}" for a, b in M2_PAIRS]
    missing = [name for name in columns if name not in header]
    if missing:
        raise SystemExit(f"{profile_csv} has no {missing}; the recipe's outputs must include "
                         "mean and reynolds_stress for this window.")
    shape = tuple(cells_kji[axis] for axis in present)
    values = numpy.full(shape + (9,), numpy.nan)
    for row in rows:
        cell = tuple(int(row[header.index(index_names[axis])]) for axis in present)
        values[cell] = [float(row[header.index(name)]) for name in columns]
    if numpy.isnan(values).any():
        raise SystemExit(f"{profile_csv} does not cover every retained cell; the window must "
                         "have sampled every cell.")
    return {"numpy": numpy, "nodes": blocks[block]["coords"],
            "mean": values[..., :3], "covariance": values[..., 3:]}


def read_metrics_rows(csv_path):
    """!
    @brief Rows of a runtime metrics CSV, one per time, as a continued run leaves them.
    @details A continued run writes a "# Continuation from step N" line into the file,
             and one continued from an earlier checkpoint than the step it stopped at
             writes the intervening steps a second time. Comment lines are skipped and,
             for a repeated time, the later row - the one the continued run wrote - is kept.
    @param[in] csv_path Path to the metrics CSV.
    @return List of row dictionaries in time order.
    """
    rows = {}
    with open(csv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(line for line in handle if not line.startswith("#")):
            rows[float(row["time"])] = row
    return [rows[t] for t in sorted(rows)]


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
    for row in read_metrics_rows(csv_path):
        if start <= float(row["time"]) <= end:
            total += float(row["u_tau_mean"]) ** 2 + float(row["u_tau_rms"]) ** 2
            rows += 1
    if rows == 0:
        raise SystemExit(f"{csv_path} has no rows with time in the window [{start}, {end}].")
    return math.sqrt(total / rows), rows


def mean_driving_acceleration(csv_path, start, end):
    """!
    @brief Time-averaged driving acceleration a driven-periodic controller applied.
    @details In a statistically stationary driven flow the mean force balances the mean
             wall stress: <a> h = u_tau^2 for a channel of half-height h, and <a> A = tau_w P
             for a duct of cross-section area A and wetted perimeter P.
    @param[in] csv_path Path to the run's driven_flow.csv.
    @param[in] start Window effective start, in solver time.
    @param[in] end Window effective end, in solver time.
    @return Tuple of (mean driving acceleration, number of rows used).
    """
    total, rows = 0.0, 0
    for row in read_metrics_rows(csv_path):
        if start <= float(row["time"]) <= end:
            total += float(row["driving_acceleration"])
            rows += 1
    if rows == 0:
        raise SystemExit(f"{csv_path} has no rows with time in the window [{start}, {end}].")
    return total / rows, rows


def main(argv=None):
    """!
    @brief Entry point.
    @param[in] argv Command-line style argument list supplied to the function.
    @return Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", required=True,
                        help="Profile CSV a post recipe wrote for the window with homogeneous set "
                             "to the two directions other than the wall-normal one.")
    parser.add_argument("--checkpoint", required=True,
                        help="Committed checkpoint bundle the profile was derived from; only its "
                             "checkpoint.meta is read, for the window's time bounds.")
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
    source.add_argument("--driven-flow-csv",
                        help="The run's driven_flow.csv; u_tau = sqrt(<a> h) from the driving "
                             "acceleration averaged over the window.")
    source.add_argument("--u-tau", type=float, help="Friction velocity supplied directly.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    args = parser.parse_args(argv)

    wall_kji = AXIS_TO_KJI[args.wall_axis]
    stream_c = AXIS_TO_COMPONENT[args.stream_axis]
    wall_c = AXIS_TO_COMPONENT[args.wall_axis]
    span_c = ({0, 1, 2} - {stream_c, wall_c}).pop()

    stats = profile_statistics(args.profile, args.window, args.grid, args.block, (wall_kji,))
    numpy, nodes = stats["numpy"], stats["nodes"]
    window = read_window(args.checkpoint, args.window)
    plane_mean = stats["mean"]

    def covariance(a, b):
        return stats["covariance"][:, M2_PAIRS.index(tuple(sorted((a, b))))]

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
    elif args.driven_flow_csv is not None:
        start, end = float(window["effective_start"]), float(window["effective_end"])
        accel, rows = mean_driving_acceleration(args.driven_flow_csv, start, end)
        u_tau = math.sqrt(accel * h)
        source_note = f"sqrt(<a>*h) with <a>={accel:.8e} over {rows} driven_flow.csv rows"
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
