"""Invariant tests for configurable spectral random-velocity generation."""

import importlib.machinery
import importlib.util
import struct
from copy import deepcopy
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


IC = load(ROOT / "generators" / "ic.gen", "picurv_spectral_random_velocity_tests")
CORE = load(ROOT / "picurv_cli" / "core.py", "picurv_core_spectral_random_velocity_tests")


def cartesian_nodes(cells=(16, 16, 16), lengths=(2*np.pi, 2*np.pi, 2*np.pi)):
    """! @brief Build Cartesian nodes. @param[in] cells Cell counts. @param[in] lengths Box lengths. @return Nodes. """
    ni, nj, nk = cells
    x, y, z = (np.linspace(0.0, lengths[i], cells[i] + 1) for i in range(3))
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    return np.stack((xx, yy, zz), axis=-1)


def params(operator="picurv_discrete"):
    """! @brief Build provider options. @param[in] operator Projection operator. @return Parameters. """
    return {"field": "Ucat", "seed": 12345,
            "random": {"distribution": "gaussian", "mean": [0.25, -0.5, 1.0]},
            "spectrum": {"type": "k4_exponential", "k0": 2.0, "k_cut": 5.0},
            "projection": {"type": "solenoidal", "operator": operator},
            "normalization": {"type": "component_rms", "target": 0.75},
            "remove_mean": True}


def independent_symbols(shape, lengths=(2*np.pi, 2*np.pi, 2*np.pi)):
    """! @brief Build audit symbols independently. @param[in] shape Storage shape. @param[in] lengths Box lengths. @return Symbols. """
    nk, nj, ni = shape
    dx, dy, dz = lengths[0]/ni, lengths[1]/nj, lengths[2]/nk
    kz = 2*np.pi*np.fft.fftfreq(nk, d=dz)[:, None, None]
    ky = 2*np.pi*np.fft.fftfreq(nj, d=dy)[None, :, None]
    kx = 2*np.pi*np.fft.fftfreq(ni, d=dx)[None, None, :]
    return (kx, ky, kz), (np.sin(kx*dx)/dx, np.sin(ky*dy)/dy, np.sin(kz*dz)/dz)


def normalized_divergence(field, symbols):
    """! @brief Compute the certification divergence ratio. @param[in] field Velocity field. @param[in] symbols Operator symbols. @return Ratio. """
    coeff = np.fft.fftn(field, axes=(0, 1, 2))
    dot = sum(symbols[c]*coeff[..., c] for c in range(3))
    denominator = np.sqrt(np.sum(sum(symbol*symbol for symbol in symbols)*np.sum(np.abs(coeff)**2, axis=-1)))
    return float(np.sqrt(np.sum(np.abs(dot)**2))/denominator)


@pytest.mark.parametrize(
    ("spectrum", "projection", "normalization"),
    [
        ({"type": "white"}, {"type": "none"}, {"type": "none"}),
        ({"type": "white"}, {"type": "none"}, {"type": "component_rms", "target": 1.0}),
        ({"type": "white"}, {"type": "solenoidal", "operator": "continuum"}, {"type": "none"}),
        ({"type": "white"}, {"type": "solenoidal", "operator": "picurv_discrete"}, {"type": "none"}),
        ({"type": "k4_exponential", "k0": 2.0, "k_cut": 3.0},
         {"type": "solenoidal", "operator": "continuum"}, {"type": "component_rms", "target": 1.0}),
        ({"type": "k4_exponential", "k0": 2.0, "k_cut": 3.0},
         {"type": "solenoidal", "operator": "picurv_discrete"}, {"type": "component_rms", "target": 1.0}),
    ],
)
def test_public_option_matrix(spectrum, projection, normalization):
    """! @brief Exercise every supported complete selector composition. @param[in] spectrum Spectrum. @param[in] projection Projection. @param[in] normalization Normalization. """
    configured = {"field": "Ucat", "seed": 9,
                  "random": {"distribution": "gaussian", "mean": [0.0, 0.0, 0.0]},
                  "spectrum": spectrum, "projection": projection,
                  "normalization": normalization, "remove_mean": True}
    full, summary = IC.generate_spectral_random_velocity(cartesian_nodes((9, 9, 9)), configured)
    physical = full[1:-1, 1:-1, 1:-1]
    continuum, discrete = independent_symbols(physical.shape[:3])
    chi_c = normalized_divergence(physical, continuum)
    chi_d = normalized_divergence(physical, discrete)
    assert summary["continuum_divergence_normalized"] == pytest.approx(chi_c, abs=2e-16)
    assert summary["picurv_discrete_divergence_normalized"] == pytest.approx(chi_d, abs=2e-16)
    if projection.get("operator") == "continuum":
        assert chi_c < 5e-16
    elif projection.get("operator") == "picurv_discrete":
        assert chi_d < 5e-16
    else:
        assert chi_c > 1e-2 and chi_d > 1e-2
    if normalization["type"] == "component_rms":
        assert np.sqrt(np.mean(np.sum(physical**2, axis=-1))/3.0) == pytest.approx(normalization["target"])


def test_seed_mean_normalization_storage_and_selected_projection_invariants():
    """! @brief Verify reusable numerical and storage contracts. """
    first, summary = IC.generate_spectral_random_velocity(cartesian_nodes(), params(), {"kinematic_viscosity": 0.0125})
    repeat, _ = IC.generate_spectral_random_velocity(cartesian_nodes(), params(), {"kinematic_viscosity": 0.0125})
    changed_params = params(); changed_params["seed"] = 54321
    changed, _ = IC.generate_spectral_random_velocity(cartesian_nodes(), changed_params)
    assert np.array_equal(first, repeat)
    assert not np.array_equal(first, changed)
    assert first.shape == (18, 18, 18, 3)
    physical = first[1:-1, 1:-1, 1:-1]
    assert np.allclose(physical.mean(axis=(0, 1, 2)), [0.25, -0.5, 1.0], atol=2e-14)
    fluctuations = physical - physical.mean(axis=(0, 1, 2), keepdims=True)
    assert np.sqrt(np.mean(np.sum(fluctuations**2, axis=-1))/3.0) == pytest.approx(0.75, abs=2e-14)
    reconstructed = np.empty_like(physical)
    for component, axis in enumerate((2, 1, 0)):
        reconstructed[..., component] = 0.25*(
            np.roll(physical[..., component], 1, axis=axis)
            + 2.0*physical[..., component]
            + np.roll(physical[..., component], -1, axis=axis)
        )
    assert summary["picurv_reconstructed_component_rms"] == pytest.approx(
        np.sqrt(np.mean(reconstructed**2, axis=(0, 1, 2))), abs=2e-14
    )
    assert summary["picurv_reconstructed_kinetic_energy"] == pytest.approx(
        0.5*np.mean(np.sum(reconstructed**2, axis=-1)), abs=2e-14
    )
    assert summary["picurv_discrete_divergence_max"] < 2e-14
    assert summary["kinematic_viscosity"] == pytest.approx(0.0125)


def test_projection_choice_retains_two_polarizations():
    """! @brief Verify selectable projectors have rank two, not a common rank-one nullspace. """
    for symbol in (np.array([1.0, 2.0, 3.0]), np.sin(np.array([1.25, -0.75, 0.5]))):
        projector = np.eye(3) - np.outer(symbol, symbol)/np.dot(symbol, symbol)
        assert np.linalg.matrix_rank(projector, tol=1e-12) == 2
        assert np.linalg.eigvalsh(projector) == pytest.approx([0.0, 1.0, 1.0], abs=1e-12)
        symbols = tuple(np.full((1, 1, 1), value) for value in symbol)
        realized = np.stack([IC.project_coefficients(np.eye(3, dtype=complex)[index].reshape(1,1,1,3), symbols)[0,0,0]
                             for index in range(3)], axis=1)
        assert np.linalg.matrix_rank(realized, tol=1e-12) == 2
    field, summary = IC.generate_spectral_random_velocity(cartesian_nodes(), params("continuum"))
    assert summary["continuum_divergence_max"] < 2e-14
    assert summary["picurv_discrete_divergence_rms"] > 1e-5
    assert np.all(np.isfinite(field))


def test_white_unprojected_unnormalized_field_preserves_defined_sampling_semantics():
    """! @brief Cover explicit white, none-projection, none-normalization, and mean-removal choices. """
    raw_params = {"field": "Ucat", "seed": 77,
                  "random": {"distribution": "gaussian", "mean": [2.0, -1.0, 0.5]},
                  "spectrum": {"type": "white"}, "projection": {"type": "none"},
                  "normalization": {"type": "none"}, "remove_mean": False}
    full, summary = IC.generate_spectral_random_velocity(cartesian_nodes((9, 9, 9)), raw_params)
    physical = full[1:-1, 1:-1, 1:-1]
    expected = np.random.default_rng(77).standard_normal((9, 9, 9, 3)) + [2.0, -1.0, 0.5]
    assert np.allclose(physical, expected, rtol=0.0, atol=1e-14)
    assert summary["continuum_divergence_rms"] > 0.1
    assert summary["picurv_discrete_divergence_rms"] > 0.1


def test_seeded_solenoidal_field_has_balanced_component_energy_and_cutoff():
    """! @brief Check finite-sample isotropic component balance and spectral cutoff statistically. """
    statistical = params()
    statistical["random"]["mean"] = [0.0, 0.0, 0.0]
    statistical["spectrum"] = {"type": "k4_exponential", "k0": 4.0, "k_cut": 8.0}
    statistical["normalization"]["target"] = 1.0
    full, _ = IC.generate_spectral_random_velocity(cartesian_nodes((24, 24, 24)), statistical)
    physical = full[1:-1, 1:-1, 1:-1]
    component_rms = np.sqrt(np.mean(physical**2, axis=(0, 1, 2)))
    assert np.max(component_rms)/np.min(component_rms) < 1.08
    coeff = np.fft.fftn(physical, axes=(0, 1, 2))
    continuum, _, _ = IC.spectral_symbols(physical.shape[:3], (2*np.pi, 2*np.pi, 2*np.pi))
    kmag = np.sqrt(sum(symbol*symbol for symbol in continuum))
    assert np.max(np.abs(coeff[kmag > 8.0 + 1e-12])) < 2e-12


def test_k4_exponential_finite_shell_target_peaks_at_k0():
    """! @brief Verify finite mode-density correction recovers the requested mean shell spectrum. """
    accumulated = None
    for seed in range(8):
        configured = {"field": "Ucat", "seed": seed,
                      "random": {"distribution": "gaussian", "mean": [0.0, 0.0, 0.0]},
                      "spectrum": {"type": "k4_exponential", "k0": 4.0, "k_cut": 8.0},
                      "projection": {"type": "none"}, "normalization": {"type": "none"},
                      "remove_mean": True}
        _, summary = IC.generate_spectral_random_velocity(cartesian_nodes((24, 24, 24)), configured)
        energy = np.array([row["energy"] for row in summary["realized_shell_spectrum"]])
        accumulated = energy if accumulated is None else accumulated + energy
    mean_energy = accumulated/8.0
    assert int(np.argmax(mean_energy)) == 4
    target = np.arange(mean_energy.size, dtype=float)**4*np.exp(-2*(np.arange(mean_energy.size)/4.0)**2)
    active = np.arange(mean_energy.size) <= 8
    assert np.corrcoef(mean_energy[active], target[active])[0, 1] > 0.995


@pytest.mark.parametrize("operator", ["continuum", "picurv_discrete"])
def test_seed_ensemble_directional_isotropy_and_spectral_tensor(operator):
    """! @brief Bound directional and spectral-tensor statistics across seeds. @param[in] operator Projection operator. """
    directional_energy = np.zeros(3)
    directional_weight = np.zeros(3)
    component_energy = np.zeros(3)
    realized_tensor = np.zeros((3, 3), dtype=complex)
    expected_tensor = np.zeros((3, 3))
    for seed in range(8):
        configured = {"seed": seed, "spectrum": {"type": "white"},
                      "projection": {"type": "solenoidal", "operator": operator},
                      "normalization": {"type": "none"}, "remove_mean": True}
        full, _ = IC.generate_spectral_random_velocity(cartesian_nodes(), configured)
        physical = full[1:-1, 1:-1, 1:-1]
        component_energy += np.mean(physical**2, axis=(0, 1, 2))
        coeff = np.fft.fftn(physical, axes=(0, 1, 2))
        continuum, discrete = independent_symbols(physical.shape[:3])
        symbols = continuum if operator == "continuum" else discrete
        continuum_arrays = np.broadcast_arrays(*continuum)
        selected_arrays = np.broadcast_arrays(*symbols)
        kmag = np.sqrt(sum(value*value for value in continuum_arrays))
        active = (kmag >= 2.0) & (kmag <= 5.0)
        modal_energy = 0.5*np.sum(np.abs(coeff)**2, axis=-1)
        magnitude = np.stack([np.abs(value) for value in continuum_arrays], axis=-1)
        maximum = magnitude.max(axis=-1)
        ties = np.sum(np.isclose(magnitude, maximum[..., None]), axis=-1)
        for axis in range(3):
            weight = (np.isclose(magnitude[..., axis], maximum) & active)/ties
            directional_energy[axis] += np.sum(modal_energy*weight)
            directional_weight[axis] += np.sum(weight)
        vectors = np.column_stack([coeff[..., component][active] for component in range(3)])
        realized_tensor += vectors.conj().T@vectors
        selected = np.column_stack([value[active] for value in selected_arrays])
        for row in selected:
            expected_tensor += np.eye(3)-np.outer(row, row)/np.dot(row, row)
    component_ratio = component_energy/component_energy.mean()
    directional_ratio = (directional_energy/directional_weight)/(directional_energy/directional_weight).mean()
    scale = np.vdot(expected_tensor, realized_tensor.real).real/np.vdot(expected_tensor, expected_tensor).real
    tensor_error = np.linalg.norm(realized_tensor.real-scale*expected_tensor)/np.linalg.norm(scale*expected_tensor)
    assert np.max(np.abs(component_ratio-1.0)) < 0.04
    assert np.max(np.abs(directional_ratio-1.0)) < 0.04
    assert tensor_error < 0.05


def test_axis_mapping_nyquist_reality_and_petsc_layout(tmp_path):
    """! @brief Guard x/i mapping, real reconstruction, dummy planes, and PETSc layout. @param[in] tmp_path Temporary path. """
    configured = params("continuum")
    configured["random"]["mean"] = [0.0, 0.0, 0.0]
    full, _ = IC.generate_spectral_random_velocity(cartesian_nodes(), configured)
    physical = full[1:-1, 1:-1, 1:-1]
    continuum, _ = independent_symbols(physical.shape[:3])
    coeff = np.fft.fftn(physical, axes=(0, 1, 2))
    correct = sum(continuum[c]*coeff[..., c] for c in range(3))
    swapped = continuum[2]*coeff[..., 0] + continuum[1]*coeff[..., 1] + continuum[0]*coeff[..., 2]
    scale = np.sqrt(np.sum(sum(symbol*symbol for symbol in continuum)*np.sum(np.abs(coeff)**2, axis=-1)))
    assert np.sqrt(np.sum(np.abs(correct)**2))/scale < 5e-16
    assert np.sqrt(np.sum(np.abs(swapped)**2))/scale > 1e-2
    assert np.isrealobj(full) and np.all(np.isfinite(full))
    assert np.array_equal(full[0], full[-2]) and np.array_equal(full[-1], full[1])
    assert np.array_equal(full[:, 0], full[:, -2]) and np.array_equal(full[:, -1], full[:, 1])
    assert np.array_equal(full[:, :, 0], full[:, :, -2]) and np.array_equal(full[:, :, -1], full[:, :, 1])
    path = tmp_path/"field.dat"
    IC.write_petsc_vec_binary(str(path), full)
    raw = path.read_bytes()
    assert struct.unpack(">ii", raw[:8]) == (IC.PETSC_VEC_CLASS_ID, full.size)
    assert np.array_equal(np.frombuffer(raw[8:], dtype=">f8").reshape(full.shape), full)


def test_composable_schema_and_grid_safeguards_reject_invalid_contracts():
    """! @brief Ensure production validation rejects invalid reusable contracts. """
    invalid = [{**params(), "field": "Ucont"},
               {**params(), "random": {"distribution": "uniform"}},
               {**params(), "random": {"distribution": "gaussian", "mean": [0.0, 1.0]}},
               {**params(), "random": {"distribution": "gaussian", "mean": [0.0, np.nan, 0.0]}},
               {**params(), "spectrum": {"type": "unknown"}},
               {**params(), "spectrum": {"type": "white", "k_cut": 3.0}},
               {**params(), "spectrum": {"type": "k4_exponential", "k0": 0.0, "k_cut": 3.0}},
               {**params(), "spectrum": {"type": "k4_exponential", "k0": 2.0, "k_cut": np.inf}},
               {**params(), "projection": {"type": "solenoidal"}},
               {**params(), "projection": {"type": "none", "operator": "continuum"}},
               {**params(), "projection": {"type": "unknown"}},
               {**params(), "projection": {"type": "solenoidal", "operator": "unknown"}},
               {**params(), "normalization": {"type": "component_rms", "target": 0.0}},
               {**params(), "normalization": {"type": "component_rms"}},
               {**params(), "normalization": {"type": "unknown"}},
               {**params(), "remove_mean": 1},
               {**params(), "seed": True},
               {**params(), "viscosity": 0.1}]
    for value in invalid:
        with pytest.raises(ValueError):
            IC.generate_spectral_random_velocity(cartesian_nodes(), value)
    nonuniform = cartesian_nodes(); nonuniform[..., 0] **= 2
    with pytest.raises(ValueError, match="uniform"):
        IC.generate_spectral_random_velocity(nonuniform, params())
    beyond = params()
    beyond["spectrum"] = {"type": "k4_exponential", "k0": 2.0, "k_cut": 100.0}
    with pytest.raises(ValueError, match="two-thirds-resolution"):
        IC.generate_spectral_random_velocity(cartesian_nodes(), beyond)


def test_provider_registry_periodicity_and_shared_fluid_scaling_contract():
    """! @brief Verify generic dispatch capabilities and density-aware scaling. """
    case = {"properties": {"scaling": {"length_ref": 2.0, "velocity_ref": 4.0},
                            "fluid": {"density": 5.0, "viscosity": 0.2}}}
    scaling = CORE.resolve_fluid_scaling(case)
    assert scaling["physical_kinematic_viscosity"] == pytest.approx(0.04)
    assert scaling["nondimensional_kinematic_viscosity"] == pytest.approx(0.005)
    assert scaling["reynolds"] == pytest.approx(200.0)
    nonfinite_case = {"properties": {"scaling": {"length_ref": np.nan, "velocity_ref": 4.0},
                                     "fluid": {"density": 5.0, "viscosity": 0.2}}}
    with pytest.raises(ValueError, match="finite"):
        CORE.resolve_fluid_scaling(nonfinite_case)
    assert "spectral_random_velocity" in CORE.GENERATED_IC_PROVIDERS
    periodic = [{"type": "PERIODIC", "handler": "geometric"} for _ in range(6)]
    resolved = CORE.resolve_initial_condition_config(
        {"mode": "generated", "generator": "spectral_random_velocity", "params": params()},
        [periodic], U_ref=4.0, provider_context={"kinematic_viscosity": 0.005})
    assert CORE.is_generated_ic_provider(resolved)
    nonperiodic = [dict(face) for face in periodic]; nonperiodic[0] = {"type": "WALL", "handler": "noslip"}
    with pytest.raises(ValueError, match="PERIODIC"):
        CORE.resolve_initial_condition_config(
            {"mode": "generated", "generator": "spectral_random_velocity", "params": params()},
            [nonperiodic], U_ref=4.0)


def test_spectral_ic_is_subordinate_to_restart_eulerian_state(tmp_path):
    """!
    @brief Verify restart authority suppresses validation and staging of every IC provider.
    @param[in] tmp_path Pytest temporary-directory fixture.
    """
    example = ROOT / "examples" / "decaying_isotropic_turbulence"
    case = CORE.read_yaml_file(str(example / "case.yml"))
    solver = CORE.read_yaml_file(str(example / "solver.yml"))
    monitor = CORE.read_yaml_file(str(example / "monitor.yml"))
    restarted_case = deepcopy(case)
    restarted_case["run_control"]["start_step"] = 100

    # The C runtime loads the checkpoint when start_step > 0, and control-file
    # generation suppresses IC staging in that state.  Validation must follow
    # the same source-authority contract.
    CORE.validate_solver_configs(
        restarted_case, solver, monitor,
        "restart-case.yml", "solver.yml", "monitor.yml",
    )

    case_path = tmp_path / "case.yml"
    solver_path = example / "solver.yml"
    monitor_path = example / "monitor.yml"
    CORE.write_yaml_file(str(case_path), restarted_case)
    run_dir = tmp_path / "run"
    (run_dir / "config").mkdir(parents=True)
    monitor_files = CORE.prepare_monitor_files(
        str(run_dir), "restart", monitor,
        {"Case": str(case_path), "Solver": str(solver_path), "Monitor": str(monitor_path)},
    )
    control_path = CORE.generate_solver_control_file(
        str(run_dir), "restart",
        {
            "case": restarted_case, "case_path": str(case_path),
            "solver": solver, "solver_path": str(solver_path),
            "monitor": monitor, "monitor_path": str(monitor_path),
        },
        1, monitor_files,
        restart_source_dir=str(tmp_path / "restart"), continue_mode=True,
    )
    control = Path(control_path).read_text(encoding="utf-8")
    assert "-finit 0" in control
    assert "-ic_dir " not in control
    assert not (run_dir / "config" / "initial_condition.generated.dat").exists()

    for inactive_ic in (
        {"mode": "file", "field": "Ucat", "source_file": "does-not-exist.dat"},
        {"mode": "generated", "generator": "not-a-provider", "params": {}},
    ):
        inactive_case = deepcopy(restarted_case)
        inactive_case["properties"]["initial_conditions"] = inactive_ic
        CORE.validate_solver_configs(
            inactive_case, solver, monitor,
            "restart-case.yml", "solver.yml", "monitor.yml",
        )


def test_conductor_and_generator_defaults_and_finite_validation_agree():
    """! @brief Keep public defaults and finite-number rejection consistent across both dispatch layers. """
    periodic = [[{"type": "PERIODIC", "handler": "geometric"} for _ in range(6)]]
    resolved = CORE.resolve_initial_condition_config(
        {"mode": "generated", "generator": "spectral_random_velocity", "params": {}},
        periodic, U_ref=1.0)
    normalized = IC.validate_spectral_random_velocity_params({})
    assert resolved["params"] == normalized == {
        "field": "Ucat", "seed": 12345,
        "random": {"distribution": "gaussian", "mean": [0.0, 0.0, 0.0]},
        "spectrum": {"type": "white"}, "projection": {"type": "none"},
        "normalization": {"type": "none"}, "remove_mean": True,
    }
    for bad in (
        {"random": {"distribution": "gaussian", "mean": [0.0, np.nan, 0.0]}},
        {"spectrum": {"type": "white", "k_cut": 3.0}},
        {"spectrum": {"type": "k4_exponential", "k0": True, "k_cut": 3.0}},
        {"normalization": {"type": "component_rms", "target": np.inf}},
    ):
        with pytest.raises(ValueError):
            CORE.resolve_initial_condition_config(
                {"mode": "generated", "generator": "spectral_random_velocity", "params": bad},
                periodic, U_ref=1.0)
