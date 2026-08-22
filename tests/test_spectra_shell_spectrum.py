"""Invariant tests for the standalone shell-spectrum generator."""

import importlib.machinery
import importlib.util
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


IC = load(ROOT / "generators" / "ic.gen", "picurv_ic_shell_spectrum_tests")
SPECTRA = load(ROOT / "generators" / "spectra.gen", "picurv_spectra_shell_spectrum_tests")

BOX = (2*np.pi, 2*np.pi, 2*np.pi)


def cartesian_nodes(cells=(16, 16, 16), lengths=BOX):
    """! @brief Build Cartesian nodes. @param[in] cells Cell counts. @param[in] lengths Box lengths. @return Nodes. """
    x, y, z = (np.linspace(0.0, lengths[i], cells[i] + 1) for i in range(3))
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    return np.stack((xx, yy, zz), axis=-1)


def write_picgrid(path, nodes):
    """! @brief Write a canonical single-block PICGRID. @param[in] path Output path. @param[in] nodes Node coordinates. @return Path. """
    km, jm, im, _ = nodes.shape
    with open(path, "w", encoding="utf-8") as fout:
        fout.write("PICGRID\n1\n")
        fout.write(f"{im} {jm} {km}\n")
        for point in nodes.reshape(-1, 3):
            fout.write(f"{point[0]:.16e} {point[1]:.16e} {point[2]:.16e}\n")
    return str(path)


def params(seed=12345, k0=2.0, k_cut=3.0):
    """! @brief Build spectral-generator parameters. @param[in] seed Seed. @param[in] k0 Peak. @param[in] k_cut Cutoff. @return Parameters. """
    return {"field": "Ucat", "seed": seed,
            "random": {"distribution": "gaussian", "mean": [0.0, 0.0, 0.0]},
            "spectrum": {"type": "k4_exponential", "k0": k0, "k_cut": k_cut},
            "projection": {"type": "solenoidal", "operator": "picurv_discrete"},
            "normalization": {"type": "component_rms", "target": 1.0},
            "remove_mean": True}


@pytest.fixture(name="staged_case")
def fixture_staged_case(tmp_path):
    """! @brief Generate an IC and stage it as PICGRID plus PETSc binary. @param[in] tmp_path Temp dir. @return Staged paths and the IC summary. """
    nodes = cartesian_nodes()
    full, summary = IC.generate_spectral_random_velocity(nodes, params())
    grid_path = write_picgrid(tmp_path / "grid.run", nodes)
    field_path = str(tmp_path / "ufield00000_0.dat")
    IC.write_petsc_vec_binary(field_path, full)
    return {"grid": grid_path, "field": field_path, "summary": summary, "full": full}


def test_measurement_is_deterministic_for_a_fixed_seed(staged_case):
    """! @brief Repeated measurement of one field must be bit-identical. @param[in] staged_case Staged case. """
    first = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    second = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    assert first["shell_spectrum"] == second["shell_spectrum"]


def test_recovers_the_prescribed_envelope(staged_case):
    """! @brief The measured spectrum must peak at k0 and vanish above k_cut. @param[in] staged_case Staged case. """
    result = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    # params() requests k0 = 2.0 with a hard cutoff at k_cut = 3.0.
    assert result["spectrum_peak_k"] == pytest.approx(2.0)
    above = [row["energy"] for row in result["shell_spectrum"] if row["k"] > 3.0 + 1e-12]
    assert max(above) < 1e-24
    active = [row["energy"] for row in result["shell_spectrum"] if 0.0 < row["k"] <= 3.0]
    assert min(active) > 0.0


def test_interior_extraction_matches_dmda_convention(staged_case):
    """! @brief Interior extraction must equal the generator's own physical slice. @param[in] staged_case Staged case. """
    values = SPECTRA.read_petsc_vec_binary(staged_case["field"])
    blocks = SPECTRA.read_picgrid_blocks(staged_case["grid"])
    interior = SPECTRA.extract_interior_cells(values, blocks[0]["dims"])
    assert interior.shape == staged_case["full"][1:-1, 1:-1, 1:-1].shape
    assert np.array_equal(interior, staged_case["full"][1:-1, 1:-1, 1:-1])


@pytest.mark.parametrize("symbol", ["continuum", "discrete"])
def test_parseval_closes_for_both_symbols(staged_case, symbol):
    """! @brief Summed shell energy must equal the resolved kinetic energy. @param[in] staged_case Staged case. @param[in] symbol Binning abscissa. """
    result = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, symbol, "none"
    )
    assert result["parseval_residual"] < 1e-12
    assert result["spectrum_total_energy"] == pytest.approx(
        result["resolved_kinetic_energy"], rel=1e-12
    )


def test_discrete_symbol_bins_below_the_continuum_ceiling(staged_case):
    """! @brief The centered-difference symbol saturates, so its support is shorter. @param[in] staged_case Staged case. """
    continuum = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    discrete = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "discrete", "none"
    )
    assert discrete["max_resolved_wavenumber"] < continuum["max_resolved_wavenumber"]
    assert discrete["shell_count"] < continuum["shell_count"]


def test_subtract_mean_domain_removes_the_zero_mode(tmp_path):
    """! @brief Domain-mean removal must empty the zero wavenumber. @param[in] tmp_path Temp dir. """
    nodes = cartesian_nodes()
    offset = {**params(), "random": {"distribution": "gaussian", "mean": [2.0, -1.0, 0.5]},
              "remove_mean": False}
    full, _ = IC.generate_spectral_random_velocity(nodes, offset)
    grid_path = write_picgrid(tmp_path / "grid.run", nodes)
    field_path = str(tmp_path / "ufield00000_0.dat")
    IC.write_petsc_vec_binary(field_path, full)

    kept = SPECTRA.generate_shell_spectrum(field_path, grid_path, 0, "continuum", "none")
    removed = SPECTRA.generate_shell_spectrum(field_path, grid_path, 0, "continuum", "domain")
    assert kept["zero_mode_energy"] > 1e-6
    assert removed["zero_mode_energy"] < 1e-24
    assert removed["parseval_residual"] < 1e-12


def test_stretched_grid_is_rejected_before_the_field_is_read(tmp_path):
    """! @brief A non-uniform axis must fail on geometry, not on field contents. @param[in] tmp_path Temp dir. """
    nodes = cartesian_nodes()
    nodes[..., 0] = np.sinh(1.4*nodes[..., 0]/BOX[0])/np.sinh(1.4)*BOX[0]
    grid_path = write_picgrid(tmp_path / "stretched.run", nodes)
    with pytest.raises(ValueError, match="uniform positive x spacing"):
        SPECTRA.generate_shell_spectrum("/nonexistent-field.dat", grid_path, 0, "continuum", "none")


def test_curvilinear_bent_channel_grid_is_rejected():
    """! @brief The bent duct has no homogeneous direction and must be refused. """
    grid = ROOT / "examples" / "bent_channel" / "bent_channel_coarse.picgrid"
    with pytest.raises(ValueError, match="shell_spectrum requires"):
        SPECTRA.generate_shell_spectrum("/nonexistent-field.dat", str(grid), 0, "continuum", "none")


def test_vector_length_mismatch_is_reported_with_expected_layout(staged_case, tmp_path):
    """! @brief A field sized for the wrong grid must name the expected scalar count. @param[in] staged_case Staged case. @param[in] tmp_path Temp dir. """
    other = write_picgrid(tmp_path / "other.run", cartesian_nodes(cells=(6, 6, 6)))
    with pytest.raises(ValueError, match="field vector length mismatch"):
        SPECTRA.generate_shell_spectrum(staged_case["field"], other, 0, "continuum", "none")


def test_block_index_outside_range_is_rejected(staged_case):
    """! @brief An out-of-range block index must be refused. @param[in] staged_case Staged case. """
    with pytest.raises(ValueError, match="outside the PICGRID block range"):
        SPECTRA.generate_shell_spectrum(staged_case["field"], staged_case["grid"], 2, "continuum", "none")


def test_cli_writes_a_spectrum_csv_matching_the_ic_schema(staged_case, tmp_path):
    """! @brief The CSV must carry the same schema ic.gen emits today. @param[in] staged_case Staged case. @param[in] tmp_path Temp dir. """
    csv_path = tmp_path / "spectrum.csv"
    code = SPECTRA.main([
        "shell-spectrum",
        "--field-file", staged_case["field"],
        "--source-grid", staged_case["grid"],
        "--spectrum-csv", str(csv_path),
    ])
    assert code == 0
    expected = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    lines = csv_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "k,energy"
    assert len(lines) == 1 + expected["shell_count"]


def test_scalar_interior_extraction_matches_the_vector_convention(staged_case):
    """! @brief A single-component payload must use the same interior window. @param[in] staged_case Staged case. """
    im, jm, km = SPECTRA.read_picgrid_blocks(staged_case["grid"])[0]["dims"]
    values = np.arange((im + 1)*(jm + 1)*(km + 1), dtype=float)
    interior = SPECTRA.extract_interior_cells_scalar(values, (im, jm, km))
    assert interior.shape == (km - 1, jm - 1, im - 1)
    assert interior[0, 0, 0] == values.reshape((km + 1, jm + 1, im + 1))[1, 1, 1]


def test_scalar_payload_length_mismatch_is_rejected(staged_case):
    """! @brief A scalar payload sized for another grid must be refused. @param[in] staged_case Staged case. """
    im, jm, km = SPECTRA.read_picgrid_blocks(staged_case["grid"])[0]["dims"]
    with pytest.raises(ValueError, match="scalar payload length mismatch"):
        SPECTRA.extract_interior_cells_scalar(np.zeros(7), (im, jm, km))


def test_field_mean_subtraction_removes_the_supplied_field(staged_case, tmp_path):
    """! @brief Subtracting a field must leave exactly the residual. @param[in] staged_case Staged case. @param[in] tmp_path Temp dir. """
    mean_path = str(tmp_path / "mean.dat")
    IC.write_petsc_vec_binary(mean_path, staged_case["full"])
    result = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "field", mean_path
    )
    # Subtracting a snapshot from itself leaves no energy at any wavenumber.
    assert result["resolved_kinetic_energy"] == pytest.approx(0.0, abs=1e-28)
    assert result["spectrum_total_energy"] == pytest.approx(0.0, abs=1e-28)


def test_field_mode_without_a_mean_is_refused(staged_case):
    """! @brief The field mode must name what it is missing. @param[in] staged_case Staged case. """
    with pytest.raises(ValueError, match="requires a mean field"):
        SPECTRA.generate_shell_spectrum(
            staged_case["field"], staged_case["grid"], 0, "continuum", "field"
        )


def test_unsampled_window_mean_is_refused(staged_case, tmp_path):
    """! @brief A mean whose count is empty must not be subtracted. @param[in] staged_case Staged case. @param[in] tmp_path Temp dir. """
    im, jm, km = SPECTRA.read_picgrid_blocks(staged_case["grid"])[0]["dims"]
    mean_path = str(tmp_path / "mean.dat")
    IC.write_petsc_vec_binary(mean_path, np.zeros_like(staged_case["full"]))
    count_path = str(tmp_path / "count.dat")
    IC.write_petsc_vec_binary(count_path, np.zeros((km + 1, jm + 1, im + 1)))
    with pytest.raises(ValueError, match="carry no sample"):
        SPECTRA.generate_shell_spectrum(
            staged_case["field"], staged_case["grid"], 0, "continuum", "field",
            mean_path, count_path,
        )


def test_integral_length_scale_is_finite_under_the_discrete_symbol(staged_case):
    """! @brief The 1/k moment must not divide by the vanishing Nyquist symbol. @param[in] staged_case Staged case. """
    continuum = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    discrete = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "discrete", "none"
    )
    # sin(k*dx)/dx vanishes on the Nyquist planes, so a 1/k moment taken on the
    # discrete symbol diverges rather than merely shifting. Both abscissae must
    # report the same physical integral scale, and it must stay inside the box.
    assert discrete["integral_length_scale"] == pytest.approx(
        continuum["integral_length_scale"], rel=1e-12
    )
    assert 0.0 < discrete["integral_length_scale"] < min(BOX)
    # The k^2 moment keeps its abscissa: the discrete operator resolves less
    # gradient, so it must report strictly less dissipation.
    assert discrete["dissipation_over_viscosity"] < continuum["dissipation_over_viscosity"]


def test_spectral_integrals_are_reported_without_viscosity(staged_case):
    """! @brief Length scales and the dissipation factor must be present and positive. @param[in] staged_case Staged case. """
    result = SPECTRA.generate_shell_spectrum(
        staged_case["field"], staged_case["grid"], 0, "continuum", "none"
    )
    assert result["integral_length_scale"] > 0.0
    assert result["taylor_microscale"] > 0.0
    assert result["dissipation_over_viscosity"] == pytest.approx(2.0*result["spectral_moment_k2"])
