"""Invariant tests for `resampled_velocity`, the import of an external periodic velocity field."""

import importlib.machinery
import importlib.util
import json
import os
import struct
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
LENGTHS = (2*np.pi, 4.0, 3.0)


def load(path, name):
    """! @brief Load a Python entry point. @param[in] path Path. @param[in] name Module name. @return Module. """
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


IC = load(ROOT / "generators" / "ic.gen", "picurv_resampled_velocity_ic_tests")
PS = load(ROOT / "generators" / "periodic_spectral.py", "picurv_resampled_velocity_spectral_tests")
CORE = load(ROOT / "picurv_cli" / "core.py", "picurv_resampled_velocity_core_tests")
PERIODIC = [{"type": "PERIODIC", "handler": "geometric"} for _ in range(6)]


def nodes(cells, lengths=LENGTHS):
    """! @brief Uniform Cartesian PICGRID nodes. @param[in] cells Cells along x, y, z. @param[in] lengths Periods. @return Nodes `[k, j, i, 3]`. """
    axes = [np.linspace(0.0, lengths[a], cells[a] + 1) for a in range(3)]
    z, y, x = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
    return np.stack((x, y, z), axis=-1)


def analytic(x, y, z, lengths=LENGTHS):
    """!
    @brief A band-limited, continuum-solenoidal velocity: a streamfunction in x-y plus w(x, y).
    @param[in] x Coordinates. @param[in] y Coordinates. @param[in] z Coordinates.
    @param[in] lengths Periods.
    @return Velocity `[..., 3]`.
    """
    kx, ky = 2*np.pi/lengths[0], 2*np.pi/lengths[1]
    u = -3*ky*np.sin(2*kx*x + 0.3)*np.sin(3*ky*y) + 0.5*ky*np.cos(kx*x)*np.cos(ky*y + 1.0)
    v = -(2*kx*np.cos(2*kx*x + 0.3)*np.cos(3*ky*y) - 0.5*kx*np.sin(kx*x)*np.sin(ky*y + 1.0))
    w = np.sin(kx*x + 2*ky*y + 0.7) + 0.0*z
    return np.stack((u, v, w), axis=-1)


def samples(cells, offset, lengths=LENGTHS):
    """! @brief Analytic samples at `(index + offset) * spacing`. @param[in] cells Counts x, y, z. @param[in] offset Fraction. @param[in] lengths Periods. @return `[k, j, i, 3]`. """
    axes = [(np.arange(cells[a]) + offset)*lengths[a]/cells[a] for a in range(3)]
    z, y, x = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
    return analytic(x, y, z, lengths)


def provider_params(path, **overrides):
    """! @brief Generator-side parameters for an `.npy` source. @param[in] path Source. @param[in] overrides Replacements. @return Mapping. """
    params = {"field": "Ucat", "source_file": str(path), "format": "npy", "box_length": list(LENGTHS),
              "sample_offset": [0.0, 0.0, 0.0], "layout": {"fastest_axis": "x", "components": "interleaved"},
              "filter": {"type": "none"}, "projection": "none", "remove_mean": False, "velocity_divisor": 1.0}
    params.update(overrides)
    return params


@pytest.fixture
def source(tmp_path):
    """! @brief A node-sampled 32-cubed analytic source saved as `.npy`. @param[in] tmp_path Temporary path. @return Path. """
    path = tmp_path / "source.npy"
    np.save(path, samples((32, 32, 32), 0.0))
    return path


@pytest.mark.parametrize("target", [(16, 16, 16), (48, 40, 24)])
def test_band_limited_field_is_resampled_exactly_to_cell_centres(source, target):
    """!
    @brief Coarsening and refining reproduce the analytic field at the target cell centres.
    @details The source is sampled on nodes, the target on cell centres, so this also checks
             the half-spacing phase shift; every mode both grids resolve is kept.
    @param[in] source Analytic source fixture.
    @param[in] target Target cell counts.
    """
    full, _summary = IC.generate_resampled_velocity(nodes(target), provider_params(source))
    assert np.abs(full[1:-1, 1:-1, 1:-1] - samples(target, 0.5)).max() < 1e-12


@pytest.mark.parametrize("kind,setting,expected", [
    ("gaussian", {"width": 0.6}, lambda k: np.exp(-k*k*0.36/24.0)),
    ("box", {"width": 0.6}, lambda k: np.sinc(0.6*k/(2*np.pi))),
    ("sharp", {"cutoff": 2.5}, lambda k: float(k <= 2.5)),
])
def test_filters_multiply_each_mode_by_their_transfer_function(kind, setting, expected):
    """!
    @brief A single mode comes through a filter scaled by exactly the kernel's transfer function.
    @param[in] kind Filter type. @param[in] setting Filter parameter. @param[in] expected Transfer of |k|.
    """
    cells = (16, 16, 16)
    for mode in (1, 2, 3):
        k = 2*np.pi*mode/LENGTHS[0]
        x = (np.arange(16) + 0.5)*LENGTHS[0]/16
        field = np.zeros((16, 16, 16, 3))
        field[..., 1] = np.cos(k*x)[None, None, :]
        first = [0.5*LENGTHS[a]/16 for a in range(3)]
        out = PS.fourier_resample(field, cells, LENGTHS, first, first, {"type": kind, **setting})
        assert out[..., 1] == pytest.approx(expected(k)*field[..., 1], abs=1e-12)


def test_ucat_projection_makes_the_staged_field_discretely_solenoidal(source):
    """!
    @brief `field: Ucat` with `picurv_discrete` has zero centred-difference divergence.
    @param[in] source Analytic source fixture.
    """
    full, summary = IC.generate_resampled_velocity(nodes((16, 16, 16)),
                                                   provider_params(source, projection="picurv_discrete"))
    assert summary["pre_projection_picurv_discrete_divergence_normalized"] > 1e-2
    assert summary["picurv_discrete_divergence_normalized"] < 1e-14
    assert 0.0 < summary["projection_removed_energy_fraction"] < 0.01
    assert summary["measurement_state"] == "staged_ucat"


def test_ucont_stages_face_centre_fluxes_in_picurv_layout(source):
    """!
    @brief `field: Ucont` writes velocity at each face centre times its area, where `generate_ucont` puts it.
    @param[in] source Analytic source fixture.
    """
    cells = (16, 12, 20)
    full, _summary = IC.generate_resampled_velocity(nodes(cells), provider_params(source, field="Ucont"))
    ni, nj, nk = cells
    d = [LENGTHS[a]/cells[a] for a in range(3)]
    centres = [(np.arange(cells[a]) + 0.5)*d[a] for a in range(3)]
    planes = [np.arange(cells[a] + 1)*d[a] for a in range(3)]
    for a, view in enumerate((full[1:nk+1, 1:nj+1, 0:ni+1, 0], full[1:nk+1, 0:nj+1, 1:ni+1, 1],
                              full[0:nk+1, 1:nj+1, 1:ni+1, 2])):
        axes = [planes[b] if b == a else centres[b] for b in range(3)]
        z, y, x = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
        area = d[(a + 1) % 3]*d[(a + 2) % 3]
        assert np.abs(view - analytic(x, y, z)[..., a]*area).max() < 1e-12


def test_ucont_projection_closes_every_cell_exactly(source):
    """!
    @brief Face-difference projection leaves every cell's net flux at round-off.
    @param[in] source Analytic source fixture.
    """
    cells = (16, 12, 20)
    full, summary = IC.generate_resampled_velocity(nodes(cells),
                                                   provider_params(source, field="Ucont", projection="picurv_discrete"))
    ni, nj, nk = cells
    fx = full[1:nk+1, 1:nj+1, 0:ni+1, 0]
    fy = full[1:nk+1, 0:nj+1, 1:ni+1, 1]
    fz = full[0:nk+1, 1:nj+1, 1:ni+1, 2]
    divergence = (fx[:, :, 1:] - fx[:, :, :-1]) + (fy[:, 1:, :] - fy[:, :-1, :]) + (fz[1:] - fz[:-1])
    assert np.abs(divergence).max() < 1e-13
    assert summary["flux_divergence_max_per_volume"] < 1e-12
    assert summary["measurement_state"] == "runtime_reconstructed_ucat_from_staged_fluxes"


def test_runtime_reconstruction_is_written_for_the_spectrum(source, tmp_path):
    """!
    @brief A `Ucont` result also writes the velocity the runtime reconstructs, when asked.
    @param[in] source Analytic source fixture. @param[in] tmp_path Temporary path.
    """
    target = tmp_path / "runtime_ucat.dat"
    full, summary = IC.generate_resampled_velocity(nodes((16, 16, 16)), provider_params(source, field="Ucont"),
                                                   {"runtime_ucat_output": str(target)})
    raw = target.read_bytes()
    count = struct.unpack(">ii", raw[:8])[1]
    values = np.frombuffer(raw[8:], dtype=">f8").reshape(18, 18, 18, 3)[1:-1, 1:-1, 1:-1]
    assert count == full.size
    assert 0.5*float(np.mean(np.sum(values*values, axis=-1))) == pytest.approx(summary["resolved_kinetic_energy"])


def test_every_container_and_layout_reads_the_same_field(source, tmp_path):
    """!
    @brief Raw (big-endian, header, z-fastest, separate), npy, and npz read identically.
    @param[in] source Analytic source fixture. @param[in] tmp_path Temporary path.
    """
    reference, _ = IC.generate_resampled_velocity(nodes((16, 16, 16)), provider_params(source))
    data = np.load(source)
    raw = tmp_path / "source.bin"
    with open(raw, "wb") as stream:
        stream.write(b"H"*100)
        stream.write(np.moveaxis(np.transpose(data, (2, 1, 0, 3)), -1, 0).astype(">f4").tobytes())
    raw_params = provider_params(raw, format="raw", cells=[32, 32, 32], layout={
        "fastest_axis": "z", "components": "separate", "dtype": "float32", "byte_order": "big", "header_bytes": 100})
    from_raw, _ = IC.generate_resampled_velocity(nodes((16, 16, 16)), raw_params)
    assert np.abs(from_raw - reference).max() < 1e-5
    npz = tmp_path / "source.npz"
    np.savez(npz, u=data[..., 0], v=data[..., 1], w=data[..., 2], uvw=data)
    for datasets in (["u", "v", "w"], ["uvw"]):
        layout = {"fastest_axis": "x", "components": "interleaved", "datasets": datasets}
        from_npz, _ = IC.generate_resampled_velocity(nodes((16, 16, 16)),
                                                     provider_params(npz, format="npz", layout=layout))
        assert np.array_equal(from_npz, reference)


def test_a_layout_that_does_not_fit_the_file_is_refused(source, tmp_path):
    """!
    @brief A raw file whose size disagrees with the declared layout, or an array with other cells, is refused.
    @param[in] source Analytic source fixture. @param[in] tmp_path Temporary path.
    """
    raw = tmp_path / "short.bin"
    raw.write_bytes(np.zeros(10, dtype="<f8").tobytes())
    with pytest.raises(ValueError, match="declared layout needs exactly"):
        IC.generate_resampled_velocity(nodes((16, 16, 16)), provider_params(
            raw, format="raw", cells=[32, 32, 32],
            layout={"fastest_axis": "x", "components": "interleaved", "dtype": "float64",
                    "byte_order": "little", "header_bytes": 0}))
    with pytest.raises(ValueError, match="params.cells declares"):
        IC.generate_resampled_velocity(nodes((16, 16, 16)), provider_params(source, cells=[16, 32, 32]))
    with pytest.raises(ValueError, match="resampling changes resolution, not the domain"):
        IC.generate_resampled_velocity(nodes((16, 16, 16)), provider_params(source, box_length=[1.0, 4.0, 3.0]))


def test_hdf5_without_h5py_names_the_fix(source, tmp_path, monkeypatch):
    """!
    @brief The optional HDF5 reader says how to install h5py rather than failing obscurely.
    @param[in] source Analytic source fixture. @param[in] tmp_path Temporary path. @param[in] monkeypatch Fixture.
    """
    monkeypatch.setitem(__import__("sys").modules, "h5py", None)
    path = tmp_path / "source.h5"
    path.write_bytes(b"")
    params = provider_params(path, format="hdf5",
                             layout={"fastest_axis": "x", "components": "interleaved", "datasets": ["u"]})
    with pytest.raises(ValueError, match="needs h5py"):
        IC.generate_resampled_velocity(nodes((16, 16, 16)), params)


def test_conductor_converts_units_and_the_provider_writes_solver_units(source):
    """!
    @brief Physical box lengths and filter widths reach the provider in solver units; velocities are divided.
    @param[in] source Analytic source fixture.
    """
    scales = {"length_ref": 2.0, "velocity_ref": 4.0}
    physical = [2.0*value for value in LENGTHS]
    resolved = CORE.resolve_initial_condition_config(
        {"mode": "generated", "generator": "resampled_velocity",
         "params": {"source_file": str(source), "format": "npy", "box_length": physical, "field": "Ucat",
                    "filter": {"type": "gaussian", "width": 1.0}, "projection": "none", "remove_mean": False}},
        [PERIODIC], scales)
    assert CORE.is_generated_ic_provider(resolved)
    params = resolved["params"]
    assert params["box_length"] == pytest.approx(list(LENGTHS))
    assert params["filter"]["width"] == pytest.approx(0.5)
    assert params["velocity_divisor"] == 4.0
    full, _ = IC.generate_resampled_velocity(nodes((16, 16, 16)), dict(params, filter={"type": "none"}))
    assert np.abs(full[1:-1, 1:-1, 1:-1] - samples((16, 16, 16), 0.5)/4.0).max() < 1e-12


def test_conductor_defaults_to_ucont_and_refuses_invalid_contracts(source):
    """!
    @brief Ucont is the default; walls, unknown keys, wrong-format layout keys, and bad values are refused.
    @param[in] source Analytic source fixture.
    @return None.
    """
    def resolve(params, blocks=None):
        """! @brief Resolve one configuration. @param[in] params Provider parameters. @param[in] blocks Boundaries. @return Contract. """
        return CORE.resolve_initial_condition_config(
            {"mode": "generated", "generator": "resampled_velocity", "params": params},
            blocks or [PERIODIC], {"length_ref": 1.0, "velocity_ref": 1.0})
    good = {"source_file": str(source), "format": "npy", "box_length": list(LENGTHS)}
    resolved = resolve(dict(good))
    assert resolved["field_code"] == 1 and resolved["params"]["field"] == "Ucont"
    assert resolved["params"]["projection"] == "picurv_discrete"
    walled = [dict(face) for face in PERIODIC]
    walled[2] = {"type": "WALL", "handler": "noslip"}
    cases = [
        (dict(good), "PERIODIC", [walled]),
        (dict(good, seed=1), "unsupported params", None),
        (dict(good, layout={"dtype": "float32"}), "do not apply to format: npy", None),
        (dict(good, format="raw", layout={"dtype": "float32"}), "requires params.cells", None),
        (dict(good, format="raw", cells=[32, 32, 32]), "layout.dtype", None),
        (dict(good, format="npz"), "layout.datasets", None),
        (dict(good, sample_offset=1.0), "in \\[0, 1\\)", None),
        (dict(good, filter={"type": "gaussian"}), "takes exactly", None),
        (dict(good, projection="spectral"), "projection must be", None),
        (dict(good, format="vtk"), "format must be one of", None),
    ]
    for params, message, blocks in cases:
        with pytest.raises(ValueError, match=message):
            resolve(params, *(blocks,) if blocks else ())


def test_mode_file_payload_for_another_grid_points_to_resampled_velocity(tmp_path):
    """!
    @brief A `mode: file` vector sized for another grid is refused at staging, naming the import provider.
    @param[in] tmp_path Temporary path.
    """
    grid = tmp_path / "grid.run"
    grid.write_text("PICGRID\n1\n5 5 5\n" + "0 0 0\n"*125, encoding="utf-8")
    good = {"path": "x.dat", "scalar_count": 6*6*6*3}
    CORE._check_file_ic_matches_grid(good, str(grid), {"field_name": "ufield"})
    with pytest.raises(ValueError, match="generator: resampled_velocity"):
        CORE._check_file_ic_matches_grid({"path": "x.dat", "scalar_count": 9*9*9*3}, str(grid),
                                         {"field_name": "ufield"})
    CORE._check_file_ic_matches_grid({"path": "x.dat", "scalar_count": 7}, str(tmp_path / "absent"), {})


def test_source_file_content_is_part_of_the_asset_identity(source, tmp_path):
    """!
    @brief Changing the source file's bytes changes the initial-condition provider's identity.
    @param[in] source Analytic source fixture. @param[in] tmp_path Temporary path.
    """
    case_path = tmp_path / "case.yml"
    case_path.write_text("", encoding="utf-8")
    spec = {"source_file": os.path.relpath(source, tmp_path)}
    before = CORE._provider_source_fingerprints(spec, str(case_path))
    np.save(source, samples((32, 32, 32), 0.0)*2.0)
    after = CORE._provider_source_fingerprints(spec, str(case_path))
    assert set(before) == {"source_file"}
    assert before["source_file"]["sha256"] != after["source_file"]["sha256"]


def test_generator_command_line_writes_vector_and_summary(source, tmp_path):
    """!
    @brief The `ic.gen --generator resampled_velocity` entry point writes a PETSc vector and a JSON summary.
    @param[in] source Analytic source fixture. @param[in] tmp_path Temporary path.
    """
    grid = tmp_path / "grid.run"
    grid_nodes = nodes((8, 8, 8))
    rows = "\n".join(" ".join(f"{value:.17g}" for value in node) for node in grid_nodes.reshape(-1, 3))
    grid.write_text(f"PICGRID\n1\n9 9 9\n{rows}\n", encoding="utf-8")
    output, summary_path = tmp_path / "ic.dat", tmp_path / "summary.json"
    code = IC.main(["--generator", "resampled_velocity", "--grid", str(grid), "--output", str(output),
                    "--params-json", json.dumps(provider_params(source, field="Ucont", projection="picurv_discrete")),
                    "--summary-json", str(summary_path)])
    assert code == 0
    assert struct.unpack(">ii", output.read_bytes()[:8])[1] == 10*10*10*3
    assert json.loads(summary_path.read_text(encoding="utf-8"))["field"] == "Ucont"
