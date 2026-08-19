"""Layout tests for the profile.gen field-slice extractor."""

import importlib.machinery
import importlib.util
import struct
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(path, name):
    """! @brief Load a Python entry point. @param[in] path Path. @param[in] name Module name. @return Module. """
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


PROFILE = load(ROOT / "generators" / "profile.gen", "picurv_profile_field_slice_tests")

NODES = (5, 6, 7)  # (IM, JM, KM)


def write_picgrid(path, dims, lengths=(1.0, 2.0, 3.0)):
    """! @brief Write a canonical Cartesian PICGRID. @param[in] path Path. @param[in] dims Node dims. @param[in] lengths Box lengths. @return Path. """
    im, jm, km = dims
    x = np.linspace(0.0, lengths[0], im)
    y = np.linspace(0.0, lengths[1], jm)
    z = np.linspace(0.0, lengths[2], km)
    with open(path, "w", encoding="utf-8") as fout:
        fout.write(f"PICGRID\n1\n{im} {jm} {km}\n")
        for kk in range(km):
            for jj in range(jm):
                for ii in range(im):
                    fout.write(f"{x[ii]:.16e} {y[jj]:.16e} {z[kk]:.16e}\n")
    return str(path)


def write_dmda_ucat(path, dims, field):
    """! @brief Write a cell-centered Ucat payload in DMDA order. @param[in] path Path. @param[in] dims Node dims. @param[in] field DMDA-shaped array. @return Path. """
    im, jm, km = dims
    assert field.shape == (km + 1, jm + 1, im + 1, 3)
    values = np.asarray(field, dtype=">f8").reshape(-1)
    with open(path, "wb") as fout:
        fout.write(struct.pack(">ii", 1211214, values.size))
        fout.write(values.tobytes())
    return str(path)


def uniform_dmda_field(dims, velocity):
    """! @brief Build a DMDA-shaped field holding one velocity everywhere. @param[in] dims Node dims. @param[in] velocity Velocity triple. @return Array. """
    im, jm, km = dims
    field = np.zeros((km + 1, jm + 1, im + 1, 3), dtype=float)
    field[..., :] = np.asarray(velocity, dtype=float)
    return field


def test_expected_length_follows_the_cell_centered_dmda(tmp_path):
    """! @brief A node-sized vector must be rejected and the DMDA size named. @param[in] tmp_path Temp dir. """
    im, jm, km = NODES
    grid = write_picgrid(tmp_path / "grid.run", NODES)
    values = np.zeros(im*jm*km*3, dtype=">f8")
    path = tmp_path / "node_sized.dat"
    with open(path, "wb") as fout:
        fout.write(struct.pack(">ii", 1211214, values.size))
        fout.write(values.tobytes())
    with pytest.raises(ValueError, match="field_file vector length mismatch"):
        PROFILE.generate_field_slice_profile(
            str(tmp_path / "out.picslice"), str(path), grid, grid, 0, 0,
            {"face": "+Zeta"}, "-Zeta", "opposite", 0.99, 1.0,
        )


def test_dummy_layers_are_never_sampled(tmp_path):
    """! @brief Poisoned dummy layers must not reach the extracted profile. @param[in] tmp_path Temp dir. """
    im, jm, km = NODES
    grid = write_picgrid(tmp_path / "grid.run", NODES)
    field = uniform_dmda_field(NODES, (0.0, 0.0, 1.0))
    # Dummy layers carry values that would produce a negative normal speed if read.
    field[0, :, :, :] = -50.0
    field[-1, :, :, :] = -50.0
    field[:, 0, :, :] = -50.0
    field[:, -1, :, :] = -50.0
    field[:, :, 0, :] = -50.0
    field[:, :, -1, :] = -50.0
    path = write_dmda_ucat(tmp_path / "ucat.dat", NODES, field)

    summary = PROFILE.generate_field_slice_profile(
        str(tmp_path / "out.picslice"), path, grid, grid, 0, 0,
        {"face": "+Zeta"}, "-Zeta", "opposite", 0.99, 1.0,
    )
    assert summary["min_speed"] == pytest.approx(1.0)
    assert summary["max_speed"] == pytest.approx(1.0)
    assert summary["dims"] == [jm - 1, im - 1]


@pytest.mark.parametrize(
    ("axis", "node_sample", "expected"),
    [
        ("Xi", (0, 0, 0), (1, 1, 1)),
        ("Xi", (3, 2, NODES[0] - 1), (4, 3, NODES[0] - 1)),
        ("Eta", (0, 0, 0), (1, 1, 1)),
        ("Eta", (3, NODES[1] - 1, 2), (4, NODES[1] - 1, 3)),
        ("Zeta", (0, 0, 0), (1, 1, 1)),
        ("Zeta", (NODES[2] - 1, 3, 2), (NODES[2] - 1, 4, 3)),
    ],
)
def test_cell_index_mapping_clamps_the_normal_axis(axis, node_sample, expected):
    """! @brief Tangential indices shift by one; the normal index clamps into the interior. @param[in] axis Axis. @param[in] node_sample Sample. @param[in] expected Expected index. """
    assert PROFILE._cell_index_for_sample(axis, node_sample, NODES) == expected


def test_every_face_reads_its_adjacent_interior_cell(tmp_path):
    """! @brief Each face must resolve to a distinct interior plane, not a dummy one. @param[in] tmp_path Temp dir. """
    im, jm, km = NODES
    for face, axis, limit in (("-Xi", "Xi", im), ("+Xi", "Xi", im),
                              ("-Eta", "Eta", jm), ("+Eta", "Eta", jm),
                              ("-Zeta", "Zeta", km), ("+Zeta", "Zeta", km)):
        spec = PROFILE._slice_specs_from_face(face, NODES)
        index = PROFILE._cell_index_for_sample(axis, (0, 0, 0) if spec["index"] == 0 else
                                               {"Xi": (0, 0, spec["index"]),
                                                "Eta": (0, spec["index"], 0),
                                                "Zeta": (spec["index"], 0, 0)}[axis], NODES)
        normal_component = {"Xi": index[2], "Eta": index[1], "Zeta": index[0]}[axis]
        assert 1 <= normal_component <= limit - 1
