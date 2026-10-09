"""!
@file test_driven_periodic_profile_tools.py
@brief Pytest coverage for the driven-periodic statistics reducers shipped with the examples.

@details Both tools read a committed checkpoint bundle directly, so they are exercised on
a synthetic bundle laid out as the solver writes one: `checkpoint.meta`, and the window's
`Ucat_mean` (dof 3), `Ucat_m2` (dof 6 centred weighted sums) and `weight` payloads as
PETSc binary vectors over the `(IM+1, JM+1, KM+1)` cell-centred DMDA.
"""

import importlib.machinery
import importlib.util
import json
import shutil
import struct
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
CHANNEL_TOOL = REPO_ROOT / "examples" / "periodic_test" / "driven_channel" / "tools" / "wall_normal_profile.py"
DUCT_TOOL = REPO_ROOT / "examples" / "periodic_test" / "driven_duct" / "tools" / "cross_section_profile.py"
PETSC_VEC_CLASSID = 1211214


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


def write_vec(path, interior, dims):
    """!
    @brief Write interior cell values as a PETSc binary Vec over the padded DMDA.
    @param[in] path Output path supplied to the function.
    @param[in] interior Array shaped (KM-1, JM-1, IM-1, dof).
    @param[in] dims PICGRID node dimensions (IM, JM, KM).
    """
    im, jm, km = dims
    full = np.zeros((km + 1, jm + 1, im + 1, interior.shape[-1]))
    full[1:km, 1:jm, 1:im] = interior
    values = full.ravel()
    path.write_bytes(struct.pack(">ii", PETSC_VEC_CLASSID, values.size) + values.astype(">f8").tobytes())


def write_bundle(root, x, y, z, mean, m2, weight):
    """!
    @brief Write a one-window checkpoint bundle and its PICGRID file.
    @param[in] root Directory to write into.
    @param[in] x Node coordinates along i.
    @param[in] y Node coordinates along j.
    @param[in] z Node coordinates along k.
    @param[in] mean Per-cell time means, (nk, nj, ni, 3).
    @param[in] m2 Per-cell centred weighted sums, (nk, nj, ni, 6).
    @param[in] weight Per-cell window weight, (nk, nj, ni).
    @return Tuple of (bundle directory, grid path).
    """
    dims = (x.size, y.size, z.size)
    grid = root / "grid.picgrid"
    rows = [f"{xi} {yj} {zk}" for zk in z for yj in y for xi in x]
    grid.write_text("PICGRID\n1\n" + " ".join(map(str, dims)) + "\n" + "\n".join(rows) + "\n",
                    encoding="utf-8")
    bundle = root / "step_000000000100"
    payloads = bundle / "statistics" / "window_0000" / "block_0000"
    payloads.mkdir(parents=True)
    write_vec(payloads / "Ucat_mean.dat", mean, dims)
    write_vec(payloads / "Ucat_m2.dat", m2, dims)
    write_vec(payloads / "weight.dat", weight[..., None], dims)
    (bundle / "checkpoint.meta").write_text(
        "-checkpoint_block_count 1\n"
        f"-checkpoint_block_0_im {dims[0]}\n-checkpoint_block_0_jm {dims[1]}\n-checkpoint_block_0_km {dims[2]}\n"
        "-checkpoint_statistics_window_count 1\n"
        "-checkpoint_statistics_window_0_name stationary\n"
        "-checkpoint_statistics_window_0_state complete\n"
        "-checkpoint_statistics_window_0_sample_count 10\n"
        "-checkpoint_statistics_window_0_represented_time 2\n"
        "-checkpoint_statistics_window_0_effective_start 1\n"
        "-checkpoint_statistics_window_0_effective_end 3\n",
        encoding="utf-8")
    return bundle, grid


def read_csv(path):
    """!
    @brief Read a tool CSV, skipping its comment header.
    @param[in] path CSV path supplied to the function.
    @return Tuple of (column names, data array).
    """
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    return lines[0].split(","), np.array([[float(v) for v in line.split(",")] for line in lines[1:]])


def test_channel_profile_reduces_folds_and_reads_the_wall_model(tmp_path):
    """!
    @brief Verify the channel reduction, the wall fold, the shear-stress sign, and wall-model u_tau.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    tool = load_tool(CHANNEL_TOOL, "wall_normal_profile_under_test")
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 2, 7), np.linspace(0, 1, 5)
    yc = 0.5 * (y[1:] + y[:-1])
    shape = (z.size - 1, y.size - 1, x.size - 1)
    rng = np.random.default_rng(3)
    mean = np.zeros(shape + (3,))
    mean[..., 2] = (yc * (2 - yc))[None, :, None] + 0.01 * rng.standard_normal(shape)
    weight = np.full(shape, 2.0)
    uv_true = (1 - yc)[None, :, None] * np.ones(shape)        # antisymmetric about the centre
    m2 = np.zeros(shape + (6,))
    m2[..., 5] = 0.04 * weight                                  # zz: <w'w'> = 0.04 per cell
    m2[..., 4] = -uv_true * weight                              # yz: <v'w'> = -(1 - y)
    bundle, grid = write_bundle(tmp_path, x, y, z, mean, m2, weight)
    csv_path = tmp_path / "wall_model.csv"
    csv_path.write_text("step,time,u_tau_mean,u_tau_rms\n1,0.5,9,0\n2,1.5,0.3,0.4\n3,2.5,0.3,0.4\n4,3.5,9,0\n",
                        encoding="utf-8")

    out = tmp_path / "profile.csv"
    assert tool.main(["--checkpoint", str(bundle), "--grid", str(grid), "--viscosity", "0.01",
                      "--u-tau", "1", "--output", str(out)]) == 0
    names, data = read_csv(out)
    U = mean[..., 2].mean(axis=(0, 2))
    spatial = (mean[..., 2] ** 2).mean(axis=(0, 2)) - U ** 2
    np.testing.assert_allclose(data[:, names.index("y")], yc[:3], rtol=1e-9)
    np.testing.assert_allclose(data[:, names.index("U_plus")], 0.5 * (U[:3] + U[::-1][:3]), rtol=1e-9)
    expected_uu = 0.04 + 0.5 * (spatial[:3] + spatial[::-1][:3])
    np.testing.assert_allclose(data[:, names.index("u_rms_plus")], np.sqrt(expected_uu), rtol=1e-8)
    # Positive near both walls once the upper half is reflected onto the lower.
    np.testing.assert_allclose(data[:, names.index("minus_uv_plus")], (1 - yc)[:3], rtol=1e-8)

    assert tool.main(["--checkpoint", str(bundle), "--grid", str(grid), "--viscosity", "0.01",
                      "--wall-model-csv", str(csv_path), "--output", str(out)]) == 0
    assert out.read_text(encoding="utf-8").startswith("# u_tau = 5.0000000000e-01")   # sqrt(0.3^2 + 0.4^2)

    driven = tmp_path / "driven_flow.csv"
    driven.write_text("step,time,direction,target_flux,measured_flux,cross_section_area,bulk_velocity,"
                      "bulk_velocity_correction,driving_acceleration,physical_time\n"
                      "1,0.5,Z,2,2,2,1,0,5.0,0.5\n2,1.5,Z,2,2,2,1,0,0.08,1.5\n3,2.5,Z,2,2,2,1,0,0.10,2.5\n",
                      encoding="utf-8")
    assert tool.main(["--checkpoint", str(bundle), "--grid", str(grid), "--viscosity", "0.01",
                      "--driven-flow-csv", str(driven), "--output", str(out)]) == 0
    assert out.read_text(encoding="utf-8").startswith("# u_tau = 3.0000000000e-01")   # sqrt(<0.08, 0.10> * 1)


def test_duct_cross_section_folds_onto_one_octant(tmp_path):
    """!
    @brief Verify the duct reduction's symmetry fold, bulk velocity, and body-force u_tau.
    @param[in] tmp_path Pytest temporary directory fixture supplied to the function.
    """
    tool = load_tool(DUCT_TOOL, "cross_section_profile_under_test")
    x = y = 1.0 - np.cos(np.linspace(0, np.pi, 7))             # clustered, mirror-symmetric
    z = np.linspace(0, 3, 5)
    shape = (z.size - 1, y.size - 1, x.size - 1)
    rng = np.random.default_rng(5)
    mean = rng.standard_normal(shape + (3,))
    mean[..., 2] += 1.0
    weight = np.full(shape, 1.0)
    m2 = rng.random(shape + (6,))
    bundle, grid = write_bundle(tmp_path, x, y, z, mean, m2, weight)

    prefix = tmp_path / "duct"
    assert tool.main(["--checkpoint", str(bundle), "--grid", str(grid), "--viscosity", "0.01",
                      "--body-force", "0.02", "--output-prefix", str(prefix)]) == 0
    names, data = read_csv(Path(str(prefix) + "_cross_section.csv"))
    n = x.size - 1
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
    bulk = float((mean[..., 2].mean(axis=0) * area).sum() / area.sum())
    assert f"U_b = {bulk:.10e}" in header

    driven = tmp_path / "driven_flow.csv"
    driven.write_text("step,time,direction,target_flux,measured_flux,cross_section_area,bulk_velocity,"
                      "bulk_velocity_correction,driving_acceleration,physical_time\n"
                      "1,2.0,Z,4,4,4,1,0,0.02,2.0\n", encoding="utf-8")
    assert tool.main(["--checkpoint", str(bundle), "--grid", str(grid), "--viscosity", "0.01",
                      "--driven-flow-csv", str(driven), "--output-prefix", str(prefix)]) == 0
    assert Path(str(prefix) + "_wall_bisector.csv").read_text(encoding="utf-8").startswith(
        "# u_tau = 1.0000000000e-01")                                      # sqrt(0.02 * 4 / 8)

    assert tool.main(["--checkpoint", str(bundle), "--grid", str(grid), "--viscosity", "0.01",
                      "--no-fold", "--output-prefix", str(prefix)]) == 0
    names, data = read_csv(Path(str(prefix) + "_cross_section.csv"))
    np.testing.assert_allclose(data[:, names.index("u_z")], mean[..., 2].mean(axis=0).ravel(), rtol=1e-9)


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
