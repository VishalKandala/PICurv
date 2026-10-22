"""!
@file test_driven_periodic_profile_tools.py
@brief Pytest coverage for the wall-profile and cross-section reducers in `generators/reduce.gen`.

@details Both read the profile CSV a post recipe writes with
`field_statistics.homogeneous` and the run's PICGRID; the conductor supplies the
window's time bounds. They are exercised on a synthetic profile laid out as
`ComputeWindowStatisticsProfile` writes one; the merge over the homogeneous directions
itself is covered by the C statistics suites.
"""

import importlib.machinery
import importlib.util
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
REDUCE_GEN = REPO_ROOT / "generators" / "reduce.gen"
PAIRS = ("xx", "xy", "xz", "yy", "yz", "zz")


def load_tool(path, name):
    """!
    @brief Import one example tool as a module.
    @param[in] path Tool path supplied to the function.
    @param[in] name Module name supplied to the function.
    @return Loaded module.
    """
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def write_case(root, x, y, z, retained, mean, covariance):
    """!
    @brief Write a PICGRID and a window's profile CSV.
    @param[in] root Directory to write into.
    @param[in] x Node coordinates along i.
    @param[in] y Node coordinates along j.
    @param[in] z Node coordinates along k.
    @param[in] retained Index names the profile retains, in (k, j, i) order.
    @param[in] mean Slot means, retained cells + (3,).
    @param[in] covariance Slot covariances, retained cells + (6,) in PAIRS order.
    @return Tuple of (profile path, grid path).
    """
    dims = (x.size, y.size, z.size)
    grid = root / "grid.picgrid"
    rows = [f"{xi} {yj} {zk}" for zk in z for yj in y for xi in x]
    grid.write_text("PICGRID\n1\n" + " ".join(map(str, dims)) + "\n" + "\n".join(rows) + "\n",
                    encoding="utf-8")
    order = [name for name in ("i", "j", "k") if name in retained]
    columns = order + ["x", "y", "z"] + [f"stationary_Ucat_mean_{c}" for c in "xyz"]
    columns += [f"stationary_Ucat_R_{pair}" for pair in PAIRS]
    lines = [",".join(columns)]
    for cell in np.ndindex(mean.shape[:-1]):
        index = dict(zip(retained, cell))
        lines.append(",".join([str(index[name]) for name in order] + ["0", "0", "0"]
                              + [repr(float(v)) for v in mean[cell]] + [repr(float(v)) for v in covariance[cell]]))
    profile = root / "Field_profile_stationary_00100.csv"
    profile.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return profile, grid


def read_csv(path):
    """!
    @brief Read a tool CSV, skipping its comment header.
    @param[in] path CSV path supplied to the function.
    @return Tuple of (column names, data array).
    """
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    return lines[0].split(","), np.array([[float(v) for v in line.split(",")] for line in lines[1:]])


def test_channel_profile_folds_and_reads_the_wall_model(tmp_path):
    """!
    @brief Verify the channel profile read, the wall fold, the shear-stress sign, and wall-model u_tau.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    tool = load_tool(REDUCE_GEN, "reduce_gen_wall")
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 2, 7), np.linspace(0, 1, 5)
    yc = 0.5 * (y[1:] + y[:-1])
    mean = np.zeros((yc.size, 3))
    mean[:, 2] = yc * (2 - yc) + 0.01 * np.random.default_rng(3).standard_normal(yc.size)
    covariance = np.zeros((yc.size, 6))
    covariance[:, 5] = 0.04 + 0.01 * yc                    # zz
    covariance[:, 4] = -(1 - yc)                           # yz: antisymmetric about the centre
    profile, grid = write_case(tmp_path, x, y, z, ("j",), mean, covariance)
    csv_path = tmp_path / "wall_model.csv"
    csv_path.write_text("step,time,u_tau_mean,u_tau_rms\n1,0.5,9,0\n2,1.5,0.3,0.4\n3,2.5,0.3,0.4\n4,3.5,9,0\n",
                        encoding="utf-8")

    out = tmp_path / "profile.csv"
    common = ["wall-profile", "--profile", str(profile), "--grid", str(grid), "--viscosity", "0.01",
              "--window-start", "1", "--window-end", "3", "--output", str(out)]
    assert tool.main(common + ["--u-tau", "1"]) == 0
    names, data = read_csv(out)
    np.testing.assert_allclose(data[:, names.index("y")], yc[:3], rtol=1e-9)
    np.testing.assert_allclose(data[:, names.index("U_plus")], 0.5 * (mean[:3, 2] + mean[::-1, 2][:3]), rtol=1e-9)
    expected_uu = 0.5 * (covariance[:3, 5] + covariance[::-1, 5][:3])
    np.testing.assert_allclose(data[:, names.index("u_rms_plus")], np.sqrt(expected_uu), rtol=1e-8)
    # Positive near both walls once the upper half is reflected onto the lower.
    np.testing.assert_allclose(data[:, names.index("minus_uv_plus")], (1 - yc)[:3], rtol=1e-8)

    assert tool.main(common + ["--wall-model-csv", str(csv_path)]) == 0
    assert out.read_text(encoding="utf-8").startswith("# u_tau = 5.0000000000e-01")   # sqrt(0.3^2 + 0.4^2)

    driven = tmp_path / "driven_flow.csv"
    driven.write_text("step,time,direction,target_flux,measured_flux,cross_section_area,bulk_velocity,"
                      "bulk_velocity_correction,driving_acceleration,physical_time\n"
                      "1,0.5,Z,2,2,2,1,0,5.0,0.5\n2,1.5,Z,2,2,2,1,0,0.08,1.5\n3,2.5,Z,2,2,2,1,0,0.10,2.5\n",
                      encoding="utf-8")
    assert tool.main(common + ["--driven-flow-csv", str(driven)]) == 0
    assert out.read_text(encoding="utf-8").startswith("# u_tau = 3.0000000000e-01")   # sqrt(<0.08, 0.10> * 1)


def test_channel_profile_rejects_a_profile_over_the_wrong_directions(tmp_path):
    """!
    @brief A profile that keeps the streamwise direction instead of the wall-normal one is refused.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    tool = load_tool(REDUCE_GEN, "reduce_gen_wrong_axes")
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 2, 7), np.linspace(0, 1, 5)
    profile, grid = write_case(tmp_path, x, y, z, ("k",), np.zeros((4, 3)), np.zeros((4, 6)))
    try:
        tool.main(["wall-profile", "--profile", str(profile), "--grid", str(grid),
                   "--viscosity", "0.01", "--u-tau", "1", "--output", str(tmp_path / "out.csv")])
    except SystemExit as error:
        assert "homogeneous" in str(error)
    else:
        raise AssertionError("a profile retaining the wrong direction must be refused")


def test_duct_cross_section_folds_onto_one_octant(tmp_path):
    """!
    @brief Verify the duct profile read, symmetry fold, bulk velocity, and body-force u_tau.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    tool = load_tool(REDUCE_GEN, "reduce_gen_duct")
    x = y = 1.0 - np.cos(np.linspace(0, np.pi, 7))             # clustered, mirror-symmetric
    z = np.linspace(0, 3, 5)
    n = x.size - 1
    rng = np.random.default_rng(5)
    mean = rng.standard_normal((n, n, 3))                      # indexed [j, i]
    mean[..., 2] += 1.0
    covariance = rng.random((n, n, 6))
    profile, grid = write_case(tmp_path, x, y, z, ("j", "i"), mean, covariance)

    prefix = tmp_path / "duct"
    common = ["cross-section", "--profile", str(profile), "--grid", str(grid), "--viscosity", "0.01",
              "--window-start", "1", "--window-end", "3", "--output-prefix", str(prefix)]
    assert tool.main(common + ["--body-force", "0.02"]) == 0
    names, data = read_csv(Path(str(prefix) + "_cross_section.csv"))
    field = {name: data[:, names.index(name)].reshape(n, n) for name in names}   # [y, x]
    np.testing.assert_allclose(field["u_z"], field["u_z"][:, ::-1], atol=1e-12)
    np.testing.assert_allclose(field["u_z"], field["u_z"].T, atol=1e-12)
    np.testing.assert_allclose(field["u_x"], -field["u_x"][:, ::-1], atol=1e-12)
    np.testing.assert_allclose(field["u_x"], field["u_y"].T, atol=1e-12)
    np.testing.assert_allclose(field["cov_xz"], -field["cov_xz"][:, ::-1], atol=1e-12)
    np.testing.assert_allclose(field["cov_xz"], field["cov_yz"].T, atol=1e-12)

    header = Path(str(prefix) + "_wall_bisector.csv").read_text(encoding="utf-8")
    assert header.startswith("# u_tau = 1.0000000000e-01")                  # sqrt(0.02 * 4 / 8)
    dx = np.diff(x)
    area = np.outer(dx, dx)
    bulk = float((mean[..., 2] * area).sum() / area.sum())
    assert f"U_b = {bulk:.10e}" in header

    driven = tmp_path / "driven_flow.csv"
    driven.write_text("step,time,direction,target_flux,measured_flux,cross_section_area,bulk_velocity,"
                      "bulk_velocity_correction,driving_acceleration,physical_time\n"
                      "1,2.0,Z,4,4,4,1,0,0.02,2.0\n", encoding="utf-8")
    assert tool.main(common + ["--driven-flow-csv", str(driven)]) == 0
    assert Path(str(prefix) + "_wall_bisector.csv").read_text(encoding="utf-8").startswith(
        "# u_tau = 1.0000000000e-01")                                      # sqrt(0.02 * 4 / 8)

    assert tool.main(common + ["--no-fold"]) == 0
    names, data = read_csv(Path(str(prefix) + "_cross_section.csv"))
    np.testing.assert_allclose(data[:, names.index("u_z")], mean[..., 2].ravel(), rtol=1e-9)


def test_metrics_readers_handle_a_continued_run(tmp_path):
    """!
    @brief Test that the metrics readers skip continuation markers and repeated steps.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    tool = load_tool(REDUCE_GEN, "reduce_gen_metrics")
    driven = tmp_path / "driven_flow.csv"
    driven.write_text(
        "step,time,driving_acceleration\n"
        "1,1.0,2.0\n2,2.0,2.0\n3,3.0,9.0\n"
        "# Continuation from step 2\n"
        "3,3.0,4.0\n4,4.0,4.0\n", encoding="utf-8")
    accel, rows = tool.mean_driving_acceleration([str(driven)], 1.0, 4.0)
    assert rows == 4
    assert accel == (2.0 + 2.0 + 4.0 + 4.0) / 4

    wall = tmp_path / "wall_model.csv"
    wall.write_text(
        "step,time,u_tau_mean,u_tau_rms\n1,1.0,3.0,4.0\n"
        "# Continuation from step 1\n2,2.0,3.0,4.0\n", encoding="utf-8")
    u_tau, rows = tool.friction_velocity_from_wall_model([str(wall)], 0.0, 5.0)
    assert rows == 2
    assert u_tau == 5.0



def test_profile_series_writes_thickness_histories(tmp_path):
    """!
    @brief A tanh layer of growing width yields its vorticity and momentum thickness per step.
    @details For U = (1 + tanh(y / w)) / 2 the vorticity thickness is 2 w and the momentum
             thickness w / 2.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    tool = load_tool(REDUCE_GEN, "reduce_gen_series")
    for step, width in ((10, 1.0), (20, 2.0)):
        y = np.linspace(-20.0, 20.0, 4001)
        lines = ["j,x,y,z,layer_Ucat_mean_x"] + [
            f"{n},0.5,{value},0.5,{0.5 * (1.0 + np.tanh(value / width))}" for n, value in enumerate(y)]
        (tmp_path / f"Field_profile_layer_{step:05d}.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    out = tmp_path / "history.csv"
    assert tool.main(["profile-series", "--profile-prefix", str(tmp_path / "Field_profile_layer"),
                      "--dt", "0.5", "--output", str(out)]) == 0
    names, data = read_csv(out)
    np.testing.assert_allclose(data[:, names.index("time")], [5.0, 10.0])
    np.testing.assert_allclose(data[:, names.index("layer_Ucat_mean_x.max_gradient_thickness")], [2.0, 4.0], rtol=1e-4)
    np.testing.assert_allclose(data[:, names.index("layer_Ucat_mean_x.momentum_thickness")], [0.5, 1.0], rtol=1e-4)
    assert "max_gradient_thickness" not in tool.profile_functionals([0.0, 1.0], [1.0, 1.0])


def test_metrics_readers_follow_a_branched_lineage(tmp_path):
    """!
    @brief A parent's rows count only up to the step its child branched from.
    @details The parent ran on to step 4 after the checkpoint at step 2 the child
             branched from, so its rows past step 2 belong to a trajectory the child
             never saw and must not enter the child's window average.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    tool = load_tool(REDUCE_GEN, "reduce_gen_lineage")
    parent = tmp_path / "parent.csv"
    parent.write_text("step,time,driving_acceleration\n1,1.0,2.0\n2,2.0,2.0\n3,3.0,9.0\n4,4.0,9.0\n",
                      encoding="utf-8")
    child = tmp_path / "child.csv"
    child.write_text("step,time,driving_acceleration\n3,3.0,4.0\n", encoding="utf-8")
    accel, rows = tool.mean_driving_acceleration([f"{parent}@2", str(child)], 0.0, 10.0)
    assert rows == 3
    assert accel == (2.0 + 2.0 + 4.0) / 3
