#!/usr/bin/env python3
"""!
@file cross_section_profile.py
@brief Reduce a duct's field-statistics window to cross-section fields and bisector profiles.

@details
A square duct has one homogeneous direction, the streamwise one, so its statistics
reduce to a cross-section rather than a wall-normal line. The reduction itself, the
window lookup in `checkpoint.meta` and the covariance over the homogeneous set, is
the one `driven_channel/tools/wall_normal_profile.py` performs and is imported from
it; this script averages along the stream only and then reads off what the case's
acceptance criteria name (`driven_duct.md` section 5):

- `<prefix>_cross_section.csv`: every cell's mean velocity and the six covariance
  components, in solver units, for contouring the secondary flow;
- `<prefix>_wall_bisector.csv`: wall units along the wall bisector, with the
  wall-normal secondary velocity;
- `<prefix>_corner_bisector.csv`: wall units along the corner bisector, with the
  secondary velocity along it (square cross-sections only);
- `<prefix>_wall_shear.csv`: the local wall shear around the perimeter, normalized
  by its perimeter mean.

The cross-section is folded onto one octant by default: its four reflections and,
when the two cross-stream axes carry identical node distributions, the diagonal
swap. A reflection reverses the velocity component normal to its mirror line, so
that component's mean and its covariances with the other two change sign. Folding
averages eight images of one long average, which is what makes a 1-3% secondary
flow readable; `--no-fold` keeps the raw cross-section, whose asymmetry is a
measure of the remaining sampling error.

The friction velocity is the perimeter mean of the wall shear: from `--body-force`
as sqrt(f A / P), the exact mean force balance; from `--u-tau`; or by default from
the first cells' resolved gradient, nu U / d, which is meaningful only when the
first cells lie in the viscous sublayer.

Usage:

@code
    cross_section_profile.py \\
        --checkpoint RUN/output/checkpoints/step_000000200000 --window stationary \\
        --grid RUN/inputs/grid/grid.run --viscosity 4.5351474e-04 \\
        --output-prefix duct
@endcode
"""

import argparse
import importlib.machinery
import importlib.util
import math
import os
import sys

CHANNEL_TOOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                            "driven_channel", "tools", "wall_normal_profile.py")


def load_channel_tool():
    """!
    @brief Import the window reduction from the channel tool rather than duplicating it.
    @return The loaded module.
    """
    loader = importlib.machinery.SourceFileLoader("picurv_wall_normal_profile", CHANNEL_TOOL)
    spec = importlib.util.spec_from_loader("picurv_wall_normal_profile", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def transform(numpy, mean, covariance, pairs, permutation, sign):
    """!
    @brief Apply a component permutation with sign changes to a mean and its covariance.
    @param[in] numpy The numpy module.
    @param[in] mean Array ending in 3 velocity components.
    @param[in] covariance Array ending in 6 covariance components, in `pairs` order.
    @param[in] pairs Component pairs of the covariance layout.
    @param[in] permutation New component index of each old component.
    @param[in] sign Sign applied to each old component.
    @return Transformed (mean, covariance).
    """
    new_mean = numpy.empty_like(mean)
    new_cov = numpy.empty_like(covariance)
    for c in range(3):
        new_mean[..., permutation[c]] = sign[c] * mean[..., c]
    for index, (a, b) in enumerate(pairs):
        target = pairs.index(tuple(sorted((permutation[a], permutation[b]))))
        new_cov[..., target] = sign[a] * sign[b] * covariance[..., index]
    return new_mean, new_cov


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
    parser.add_argument("--stream-axis", default="Zeta", choices=("Xi", "Eta", "Zeta"),
                        help="Streamwise (driven, homogeneous) axis token.")
    parser.add_argument("--viscosity", type=float, required=True,
                        help="Kinematic viscosity in solver units.")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--body-force", type=float,
                        help="Converged driving body force f; u_tau = sqrt(f*A/P).")
    source.add_argument("--u-tau", type=float, help="Friction velocity supplied directly.")
    parser.add_argument("--no-fold", action="store_true",
                        help="Keep the raw cross-section instead of folding it onto one octant.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for the four output CSVs.")
    args = parser.parse_args(argv)

    channel = load_channel_tool()
    pairs = channel.M2_PAIRS
    stream_kji = channel.AXIS_TO_KJI[args.stream_axis]
    stream_c = channel.AXIS_TO_COMPONENT[args.stream_axis]
    stats = channel.homogeneous_statistics(args.checkpoint, args.window, args.grid, args.block,
                                           (stream_kji,))
    numpy, nodes, window = stats["numpy"], stats["nodes"], stats["window"]
    mean, covariance = stats["mean"], stats["covariance"]

    # The two remaining axes, in (k, j, i) order, and the velocity component normal
    # to each one's walls. The reduced arrays are indexed (first, second).
    cross_kji = [axis for axis in (0, 1, 2) if axis != stream_kji]
    kji_to_component = {0: 2, 1: 1, 2: 0}
    comp = [kji_to_component[axis] for axis in cross_kji]
    names = "xyz"

    def axis_nodes(kji, component):
        line = {0: nodes[:, 0, 0, :], 1: nodes[0, :, 0, :], 2: nodes[0, 0, :, :]}[kji]
        return line[:, component]

    edges = [axis_nodes(kji, c) for kji, c in zip(cross_kji, comp)]
    # Bisectors and wall shear need a rectilinear cross-section: each cross-stream
    # coordinate must depend on its own index alone.
    span = [abs(e[-1] - e[0]) for e in edges]
    for kji, c, e, length in zip(cross_kji, comp, edges, span):
        full = nodes[..., c]
        reference = {0: e[:, None, None], 1: e[None, :, None], 2: e[None, None, :]}[kji]
        if numpy.max(numpy.abs(full - reference)) > 1.0e-9 * length:
            raise SystemExit("the cross-section is not rectilinear; this reduction assumes it is.")
    centres = [0.5 * (e[:-1] + e[1:]) for e in edges]
    widths = [numpy.abs(numpy.diff(e)) for e in edges]
    n0, n1 = mean.shape[0], mean.shape[1]
    if (centres[0].size, centres[1].size) != (n0, n1):
        raise SystemExit(f"cross-section cell counts {(centres[0].size, centres[1].size)} "
                         f"do not match the statistics {(n0, n1)}.")

    def symmetric(c, e, length):
        return numpy.max(numpy.abs((c - e[0]) - (e[-1] - c[::-1]))) <= 1.0e-9 * length

    folded = "raw (--no-fold)"
    if not args.no_fold:
        if not all(symmetric(c, e, length) for c, e, length in zip(centres, edges, span)):
            raise SystemExit("the cross-section grid is not mirror-symmetric; rerun with --no-fold.")
        identity = (0, 1, 2)
        images = []
        for flip0 in (False, True):
            for flip1 in (False, True):
                m, v = mean, covariance
                sign = [1.0, 1.0, 1.0]
                if flip0:
                    m, v = m[::-1, :], v[::-1, :]
                    sign[comp[0]] = -1.0
                if flip1:
                    m, v = m[:, ::-1], v[:, ::-1]
                    sign[comp[1]] = -1.0
                images.append(transform(numpy, m, v, pairs, identity, sign))
        diagonal = n0 == n1 and numpy.max(
            numpy.abs((edges[0] - edges[0][0]) - (edges[1] - edges[1][0]))) <= 1.0e-9 * span[0]
        if diagonal:
            swap = list(identity)
            swap[comp[0]], swap[comp[1]] = comp[1], comp[0]
            images += [transform(numpy, numpy.swapaxes(m, 0, 1), numpy.swapaxes(v, 0, 1),
                                 pairs, swap, (1.0, 1.0, 1.0)) for m, v in images]
        mean = sum(m for m, _ in images) / len(images)
        covariance = sum(v for _, v in images) / len(images)
        folded = f"folded over {len(images)} symmetry images"

    nu = args.viscosity
    U = mean[..., stream_c]
    area = numpy.outer(widths[0], widths[1])
    bulk = float(numpy.sum(U * area) / numpy.sum(area))

    # Local wall shear from the first cell on each of the four walls: nu U / d.
    first = [abs(centres[0][0] - edges[0][0]), abs(edges[0][-1] - centres[0][-1]),
             abs(centres[1][0] - edges[1][0]), abs(edges[1][-1] - centres[1][-1])]
    walls = [
        (f"{names[comp[0]]}_min", centres[1], widths[1], nu * U[0, :] / first[0]),
        (f"{names[comp[0]]}_max", centres[1], widths[1], nu * U[-1, :] / first[1]),
        (f"{names[comp[1]]}_min", centres[0], widths[0], nu * U[:, 0] / first[2]),
        (f"{names[comp[1]]}_max", centres[0], widths[0], nu * U[:, -1] / first[3]),
    ]
    perimeter = sum(float(numpy.sum(w)) for _, _, w, _ in walls)
    tau_mean = sum(float(numpy.sum(t * w)) for _, _, w, t in walls) / perimeter

    if args.body_force is not None:
        cross_area = span[0] * span[1]
        u_tau = math.sqrt(args.body_force * cross_area / perimeter)
        source_note = f"sqrt(f*A/P) with f={args.body_force:.8e}"
    elif args.u_tau is not None:
        u_tau = args.u_tau
        source_note = "supplied with --u-tau"
    else:
        u_tau = math.sqrt(tau_mean)
        source_note = "sqrt of the perimeter-mean nu*U/d from the first cells (needs a resolved sublayer)"
    if not u_tau > 0.0:
        raise SystemExit("computed a non-positive friction velocity; check the inputs.")

    secondary = numpy.hypot(mean[..., comp[0]], mean[..., comp[1]])
    peak = numpy.unravel_index(numpy.argmax(secondary), secondary.shape)
    half_width = 0.5 * min(span)
    header = (f"# u_tau = {u_tau:.10e}  ({source_note})\n"
              f"# Re_tau = u_tau*a/nu = {u_tau * half_width / nu:.6f} with half-width a = {half_width:.10e}\n"
              f"# U_b = {bulk:.10e}, max secondary speed / U_b = {float(secondary[peak]) / bulk:.6e} "
              f"at ({centres[0][peak[0]]:.6e}, {centres[1][peak[1]]:.6e})\n"
              f"# nu = {nu:.10e}, window = {args.window} ({window.get('state')}, "
              f"{window.get('sample_count')} samples), {folded}\n")
    prefix = args.output_prefix
    pair_names = [names[a] + names[b] for a, b in pairs]

    with open(prefix + "_cross_section.csv", "w", encoding="utf-8") as handle:
        handle.write(header)
        handle.write(f"{names[comp[0]]},{names[comp[1]]},u_x,u_y,u_z,"
                     + ",".join("cov_" + p for p in pair_names) + "\n")
        for a in range(n0):
            for b in range(n1):
                handle.write(f"{centres[0][a]:.10e},{centres[1][b]:.10e},"
                             + ",".join(f"{value:.10e}" for value in mean[a, b])
                             + "," + ",".join(f"{value:.10e}" for value in covariance[a, b]) + "\n")

    def cov(a_index, b_index, c0, c1):
        return covariance[a_index, b_index, pairs.index(tuple(sorted((c0, c1))))]

    # Wall bisector: from the first axis's lower wall, at the centre of the second.
    middle = [n1 // 2] if n1 % 2 else [n1 // 2 - 1, n1 // 2]
    orientation = 1.0 if edges[0][-1] > edges[0][0] else -1.0
    with open(prefix + "_wall_bisector.csv", "w", encoding="utf-8") as handle:
        handle.write(header)
        handle.write("d,d_plus,U_plus,u_rms_plus,v_rms_plus,w_rms_plus,minus_uv_plus,V_over_Ub\n")
        for a in range((n0 + 1) // 2):
            pick = lambda f: sum(f(a, b) for b in middle) / len(middle)
            d = abs(centres[0][a] - edges[0][0])
            values = (
                pick(lambda i, j: U[i, j]) / u_tau,
                math.sqrt(max(pick(lambda i, j: cov(i, j, stream_c, stream_c)), 0.0)) / u_tau,
                math.sqrt(max(pick(lambda i, j: cov(i, j, comp[0], comp[0])), 0.0)) / u_tau,
                math.sqrt(max(pick(lambda i, j: cov(i, j, comp[1], comp[1])), 0.0)) / u_tau,
                -orientation * pick(lambda i, j: cov(i, j, stream_c, comp[0])) / u_tau ** 2,
                orientation * pick(lambda i, j: mean[i, j, comp[0]]) / bulk,
            )
            handle.write(f"{d:.10e},{d * u_tau / nu:.10e}," + ",".join(f"{v:.10e}" for v in values) + "\n")

    if n0 == n1:
        # Corner bisector: the diagonal cells from the (min, min) corner. Positive
        # along-diagonal velocity points away from the corner.
        sign0 = 1.0 if edges[0][-1] > edges[0][0] else -1.0
        sign1 = 1.0 if edges[1][-1] > edges[1][0] else -1.0
        with open(prefix + "_corner_bisector.csv", "w", encoding="utf-8") as handle:
            handle.write(header)
            handle.write("d,d_plus,U_plus,along_diagonal_over_Ub,k_plus\n")
            for a in range((n0 + 1) // 2):
                d = math.hypot(centres[0][a] - edges[0][0], centres[1][a] - edges[1][0])
                along = (sign0 * mean[a, a, comp[0]] + sign1 * mean[a, a, comp[1]]) / math.sqrt(2.0)
                k = 0.5 * sum(cov(a, a, c, c) for c in range(3))
                handle.write(f"{d:.10e},{d * u_tau / nu:.10e},{U[a, a] / u_tau:.10e},"
                             f"{along / bulk:.10e},{k / u_tau ** 2:.10e}\n")

    with open(prefix + "_wall_shear.csv", "w", encoding="utf-8") as handle:
        handle.write(header)
        handle.write("wall,s,tau_w,tau_w_over_mean\n")
        for name, coordinate, _width, tau in walls:
            for s_value, t_value in zip(coordinate, tau):
                handle.write(f"{name},{s_value:.10e},{t_value:.10e},{t_value / tau_mean:.10e}\n")

    print(header, end="")
    print(f"wrote {prefix}_cross_section.csv, {prefix}_wall_bisector.csv, "
          + (f"{prefix}_corner_bisector.csv, " if n0 == n1 else "")
          + f"{prefix}_wall_shear.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
