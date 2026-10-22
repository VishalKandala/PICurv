"""!
@file test_driven_periodic_profile_tools.py
@brief Pytest coverage for the driven-periodic statistics reducers shipped with the examples.

@details Both tools read the profile CSV a post recipe writes with
`field_statistics.homogeneous`, the run's PICGRID, and the window scalars in the
checkpoint's `checkpoint.meta`. They are exercised on a synthetic profile laid out as
`ComputeWindowStatisticsProfile` writes one; the merge over the homogeneous directions
itself is covered by the C statistics suites.
"""

import importlib.machinery
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
CHANNEL_TOOL = REPO_ROOT / "examples" / "periodic_test" / "driven_channel" / "tools" / "wall_normal_profile.py"
DUCT_TOOL = REPO_ROOT / "examples" / "periodic_test" / "driven_duct" / "tools" / "cross_section_profile.py"
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
    @brief Write a PICGRID, a one-window checkpoint.meta, and the window's profile CSV.
    @param[in] root Directory to write into.
    @param[in] x Node coordinates along i.
    @param[in] y Node coordinates along j.
    @param[in] z Node coordinates along k.
    @param[in] retained Index names the profile retains, in (k, j, i) order.
    @param[in] mean Slot means, retained cells + (3,).
    @param[in] covariance Slot covariances, retained cells + (6,) in PAIRS order.
    @return Tuple of (profile path, bundle directory, grid path).
    """
    dims = (x.size, y.size, z.size)
    grid = root / "grid.picgrid"
    rows = [f"{xi} {yj} {zk}" for zk in z for yj in y for xi in x]
    grid.write_text("PICGRID\n1\n" + " ".join(map(str, dims)) + "\n" + "\n".join(rows) + "\n",
                    encoding="utf-8")
    bundle = root / "step_000000000100"
    bundle.mkdir()
    (bundle / "checkpoint.meta").write_text(
        "-checkpoint_statistics_window_count 1\n"
        "-checkpoint_statistics_window_0_name stationary\n"
        "-checkpoint_statistics_window_0_state complete\n"
        "-checkpoint_statistics_window_0_sample_count 10\n"
        "-checkpoint_statistics_window_0_represented_time 2\n"
        "-checkpoint_statistics_window_0_effective_start 1\n"
        "-checkpoint_statistics_window_0_effective_end 3\n",
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
    return profile, bundle, grid


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
    tool = load_tool(CHANNEL_TOOL, "wall_normal_profile_under_test")
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 2, 7), np.linspace(0, 1, 5)
    yc = 0.5 * (y[1:] + y[:-1])
    mean = np.zeros((yc.size, 3))
    mean[:, 2] = yc * (2 - yc) + 0.01 * np.random.default_rng(3).standard_normal(yc.size)
    covariance = np.zeros((yc.size, 6))
    covariance[:, 5] = 0.04 + 0.01 * yc                    # zz
    covariance[:, 4] = -(1 - yc)                           # yz: antisymmetric about the centre
    profile, bundle, grid = write_case(tmp_path, x, y, z, ("j",), mean, covariance)
    csv_path = tmp_path / "wall_model.csv"
    csv_path.write_text("step,time,u_tau_mean,u_tau_rms\n1,0.5,9,0\n2,1.5,0.3,0.4\n3,2.5,0.3,0.4\n4,3.5,9,0\n",
                        encoding="utf-8")

    out = tmp_path / "profile.csv"
    common = ["--profile", str(profile), "--checkpoint", str(bundle), "--grid", str(grid),
              "--viscosity", "0.01", "--output", str(out)]
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
    tool = load_tool(CHANNEL_TOOL, "wall_normal_profile_wrong_axes")
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 2, 7), np.linspace(0, 1, 5)
    profile, bundle, grid = write_case(tmp_path, x, y, z, ("k",), np.zeros((4, 3)), np.zeros((4, 6)))
    try:
        tool.main(["--profile", str(profile), "--checkpoint", str(bundle), "--grid", str(grid),
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
    tool = load_tool(DUCT_TOOL, "cross_section_profile_under_test")
    x = y = 1.0 - np.cos(np.linspace(0, np.pi, 7))             # clustered, mirror-symmetric
    z = np.linspace(0, 3, 5)
    n = x.size - 1
    rng = np.random.default_rng(5)
    mean = rng.standard_normal((n, n, 3))                      # indexed [j, i]
    mean[..., 2] += 1.0
    covariance = rng.random((n, n, 6))
    profile, bundle, grid = write_case(tmp_path, x, y, z, ("j", "i"), mean, covariance)

    prefix = tmp_path / "duct"
    common = ["--profile", str(profile), "--checkpoint", str(bundle), "--grid", str(grid),
              "--viscosity", "0.01", "--output-prefix", str(prefix)]
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


def test_tools_find_the_checkout_from_an_initialized_case(tmp_path):
    """!
    @brief Verify a copy outside the repository resolves the checkout through `.picurv-origin.json`.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    case = tmp_path / "deep" / "case"
    shutil.copytree(CHANNEL_TOOL.parent, case / "driven_channel" / "tools")
    (case / ".picurv-origin.json").write_text(json.dumps({"source_repo_root": str(REPO_ROOT)}),
                                              encoding="utf-8")
    copy = load_tool(case / "driven_channel" / "tools" / CHANNEL_TOOL.name, "wall_normal_profile_copy")
    assert Path(copy.REPO_ROOT) == REPO_ROOT


def test_metrics_readers_handle_a_continued_run(tmp_path):
    """!
    @brief Test that the metrics readers skip continuation markers and repeated steps.
    @param[in] tmp_path Pytest temporary-directory fixture supplied to the function.
    """
    tool = load_tool(CHANNEL_TOOL, "wall_normal_profile_continued")
    driven = tmp_path / "driven_flow.csv"
    driven.write_text(
        "step,time,driving_acceleration\n"
        "1,1.0,2.0\n2,2.0,2.0\n3,3.0,9.0\n"
        "# Continuation from step 2\n"
        "3,3.0,4.0\n4,4.0,4.0\n", encoding="utf-8")
    accel, rows = tool.mean_driving_acceleration(str(driven), 1.0, 4.0)
    assert rows == 4
    assert accel == (2.0 + 2.0 + 4.0 + 4.0) / 4

    wall = tmp_path / "wall_model.csv"
    wall.write_text(
        "step,time,u_tau_mean,u_tau_rms\n1,1.0,3.0,4.0\n"
        "# Continuation from step 1\n2,2.0,3.0,4.0\n", encoding="utf-8")
    u_tau, rows = tool.friction_velocity_from_wall_model(str(wall), 0.0, 5.0)
    assert rows == 2
    assert u_tau == 5.0

