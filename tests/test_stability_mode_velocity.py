"""Tests for the `stability_mode_velocity` initial-condition provider in `generators/ic.gen`."""

import importlib.machinery
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(path, name):
    """!
    @brief Load a module from a file path.
    @param[in] path Module source path.
    @param[in] name Module name.
    @return Loaded module.
    """
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


IC = load(ROOT / "generators" / "ic.gen", "picurv_stability_mode_tests")
CORE = load(ROOT / "picurv_cli" / "core.py", "picurv_core_stability_mode_tests")


def box_nodes(cells, lengths, origin, stretch_cross=False):
    """!
    @brief Axis-aligned PICGRID nodes in k/j/i/component order.
    @param[in] cells Cell counts in Cartesian (i, j, k) order.
    @param[in] lengths Box lengths in Cartesian order.
    @param[in] origin Box origin in Cartesian order.
    @param[in] stretch_cross Whether to cluster the j axis toward its middle.
    @return Node array.
    """
    axes = [np.linspace(origin[a], origin[a] + lengths[a], cells[a] + 1) for a in range(3)]
    if stretch_cross:
        t = np.linspace(-1.0, 1.0, cells[1] + 1)
        axes[1] = origin[1] + lengths[1] * (1.0 + np.sinh(1.5 * t) / np.sinh(1.5)) / 2.0
    k, j, i = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
    return np.stack([i, j, k], axis=-1)


def most_unstable(profile, alphas, extent=15.0, points=401):
    """!
    @brief Scan the Rayleigh growth rate over wavenumbers.
    @param[in] profile Base velocity as a function of y.
    @param[in] alphas Wavenumbers to scan.
    @param[in] extent Half-width of the cross-stream interval.
    @param[in] points Eigen-grid points.
    @return The wavenumber and growth rate at the parabolic peak of the scanned growth rates.
    """
    y = np.linspace(-extent, extent, points)
    rates = []
    for alpha in alphas:
        c, _ = IC.rayleigh_unstable_mode(y, profile(y), alpha)
        rates.append(0.0 if c is None else alpha * c.imag)
    best = int(np.clip(np.argmax(rates), 1, len(rates) - 2))
    a, b, c = np.polyfit(alphas[best - 1:best + 2], rates[best - 1:best + 2], 2)
    return -b / (2 * a), c - b * b / (4 * a)


def test_tanh_layer_reproduces_michalke():
    """!
    @brief U = (1 + tanh y)/2 grows fastest at alpha = 0.4449 with omega_i = 0.0949 (Michalke 1964).
    """
    alpha, rate = most_unstable(lambda y: 0.5 * (1.0 + np.tanh(y)), np.linspace(0.40, 0.49, 10))
    assert alpha == pytest.approx(0.4449, abs=0.003)
    assert rate == pytest.approx(0.0949, rel=0.01)


def test_neutral_wavenumbers_of_tanh_layer_and_bickley_jet():
    """!
    @brief The tanh layer is neutral at alpha = 1 and the Bickley jet's sinuous mode at alpha = 2.
    """
    y = np.linspace(-15.0, 15.0, 401)
    assert IC.rayleigh_unstable_mode(y, np.tanh(y), 1.05)[0] is None
    assert IC.rayleigh_unstable_mode(y, np.tanh(y), 0.9)[0].imag > 0.0
    jet = 1.0 / np.cosh(y) ** 2
    assert IC.rayleigh_unstable_mode(y, jet, 1.9)[0].imag > 0.0
    assert IC.rayleigh_unstable_mode(y, jet, 2.1)[0] is None


def layer_params(**changes):
    """!
    @brief Parameters for a tanh layer forced at the fundamental and subharmonic of two wavelengths.
    @param[in] changes Parameter overrides.
    @return Parameter mapping.
    """
    params = {"streamwise_axis": "k", "cross_axis": "j", "base_profile": "0.5*tanh(y)",
              "modes": [{"harmonic": 2, "amplitude": 0.01}, {"harmonic": 1, "amplitude": 0.01, "phase": 0.5}]}
    params.update(changes)
    return params


LAYER_LENGTH = 2.0 * 2.0 * np.pi / 0.4449


@pytest.mark.parametrize("stretch", [False, True])
def test_field_is_divergence_free_planar_and_normalized(stretch):
    """!
    @brief The field is discretely divergence-free, has no normal velocity or normal variation,
           and each mode's streamwise fluctuation peaks at amplitude * dU.
    @param[in] stretch Whether the cross-stream axis is stretched.
    """
    nodes = box_nodes((4, 64, 96), (0.4, 30.0, LAYER_LENGTH), (0.0, -15.0, 0.0), stretch)
    noise = {"rms": 0.02, "planar": True, "seed": 3, "spectrum": {"type": "white"}}
    full, summary = IC.generate_stability_mode_velocity(nodes, layer_params(noise=noise))
    interior = full[1:-1, 1:-1, 1:-1]
    assert summary["picurv_discrete_divergence_max"] < 1e-12
    assert np.max(np.abs(interior[..., 0])) == 0.0
    assert np.max(np.std(interior[..., 1:], axis=2)) < 1e-14
    assert summary["velocity_difference"] == pytest.approx(1.0, rel=1e-9)
    assert summary["fastest_growing_harmonic"]["harmonic"] == 2
    assert summary["modes"][0]["growth_rate"] == pytest.approx(0.0949, rel=0.01)
    single, _ = IC.generate_stability_mode_velocity(nodes, layer_params(modes=[{"harmonic": 2, "amplitude": 0.03}]))
    fluctuation = single[1:-1, 1:-1, 1:-1, 2] - single[1:-1, 1:-1, 1:-1, 2].mean(axis=0, keepdims=True)
    peak = np.max(np.sqrt(2.0 * np.mean(fluctuation ** 2, axis=(0, 2))))
    assert peak == pytest.approx(0.03, rel=1e-9)


def test_profile_is_physical_and_scaled_by_reference_velocity():
    """!
    @brief The base profile is evaluated in physical units and written in solver units.
    """
    nodes = box_nodes((4, 32, 48), (0.4, 30.0, LAYER_LENGTH), (0.0, -15.0, 0.0))
    params = layer_params(base_profile="5.0*tanh(y/2)", modes=[{"harmonic": 1}])
    full, summary = IC.generate_stability_mode_velocity(nodes, params, {"length_ref": 2.0, "velocity_ref": 10.0})
    # Physical y = 2 * grid y, so the solver-unit profile is 0.5 tanh(grid y).
    assert summary["velocity_difference"] == pytest.approx(1.0, rel=1e-9)
    assert full[1:-1, -2, 1:-1, 2].mean() == pytest.approx(0.5 * np.tanh(15.0 - 30.0 / 64), abs=2e-2)


@pytest.mark.parametrize("bad,message", [
    ({"modes": [{"harmonic": 40}]}, "two-thirds"),
    ({"base_profile": "y"}, "no unstable Rayleigh mode"),
    ({"base_profile": "tanh(z)"}, ""),
    ({"surprise": 1}, "accepts only"),
    ({"modes": [{"harmonic": 1}, {"harmonic": 1}]}, "only once"),
    ({"cross_axis": "k"}, "different axes"),
])
def test_invalid_parameters_are_refused(bad, message):
    """!
    @brief Unresolved harmonics, stable profiles, wrong variables and malformed modes are refused.
    @param[in] bad Parameter override.
    @param[in] message Expected message fragment.
    """
    nodes = box_nodes((4, 32, 96), (0.4, 30.0, LAYER_LENGTH), (0.0, -15.0, 0.0))
    with pytest.raises(ValueError, match=message):
        IC.generate_stability_mode_velocity(nodes, layer_params(**bad))


def prepared(cross_type="WALL"):
    """!
    @brief One block of boundary faces: periodic i and k, cross-stream j of the given type.
    @param[in] cross_type Mathematical type of the j faces.
    @return Prepared boundary blocks.
    """
    faces = []
    for face in ("-Xi", "+Xi", "-Zeta", "+Zeta"):
        faces.append({"face": face, "type": "PERIODIC", "handler": "geometric"})
    for face in ("-Eta", "+Eta"):
        faces.append({"face": face, "type": cross_type,
                      "handler": "geometric" if cross_type == "PERIODIC" else "noslip"})
    return [faces]


def test_conductor_routes_parameters_and_reference_scales():
    """!
    @brief The resolver accepts the layer, converts velocity_difference and hands the provider the scales.
    """
    ic = {"mode": "generated", "generator": "stability_mode_velocity",
          "params": layer_params(velocity_difference=4.0)}
    resolved = CORE.resolve_initial_condition_config(ic, prepared(), scales={"length_ref": 2.0, "velocity_ref": 2.0})
    assert resolved["kind"] == "stability_mode_velocity"
    assert resolved["params"]["velocity_difference"] == pytest.approx(2.0)
    assert resolved["provider_context"]["length_ref"] == 2.0
    assert resolved["provider_context"]["velocity_ref"] == 2.0
    with pytest.raises(ValueError, match="non-periodic cross-stream"):
        CORE.resolve_initial_condition_config(ic, prepared("PERIODIC"), scales={"length_ref": 1.0, "velocity_ref": 1.0})
