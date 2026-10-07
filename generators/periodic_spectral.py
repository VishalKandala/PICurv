"""!
@file periodic_spectral.py
@brief Fourier machinery for the uniform periodic axes of PICurv's Cartesian grids.
@details Shared by `ic.gen`, which synthesizes and resamples initial conditions, and
         `spectra.gen`, which measures spectra. Both need the same answer to three
         questions: which axes may be transformed, what the continuum and PICurv-discrete
         wavenumbers are, and how a periodic field is moved between resolutions. Arrays are
         in PICurv storage order, `[k, j, i, component]`; lengths and wavenumbers are in
         Cartesian `x, y, z` order. Every function imports NumPy through
         `require_numpy()`, so importing this module never imports NumPy itself.
"""

import sys


def _drop_imported_package(package_name: str):
    """!
    @brief Remove a failed/partial import package tree from sys.modules.
    @param[in] package_name Top-level package name.
    """
    prefix = package_name + "."
    for module_name in list(sys.modules):
        if module_name == package_name or module_name.startswith(prefix):
            sys.modules.pop(module_name, None)


def _prune_incompatible_python_site_paths(paths):
    """!
    @brief Remove site-package paths for a different Python major/minor version.
    @param[in] paths Candidate sys.path entries.
    @return Filtered path list.
    """
    current = (sys.version_info[0], sys.version_info[1])
    filtered = []
    for path in paths:
        text = str(path)
        marker = "python"
        idx = text.lower().find(marker)
        if idx >= 0 and ("site-packages" in text or "dist-packages" in text):
            version_text = text[idx + len(marker):idx + len(marker) + 4].strip("-/")
            parts = version_text.split(".")
            try:
                path_version = (int(parts[0]), int(parts[1]))
            except (IndexError, ValueError):
                filtered.append(path)
                continue
            if path_version != current:
                continue
        filtered.append(path)
    return filtered


def require_numpy():
    """!
    @brief Import NumPy with the same compatibility retry used by the conductor.
    @return Imported NumPy module.
    """
    try:
        import numpy
        return numpy
    except Exception as exc:
        first_error = exc
        original_path = list(sys.path)
        try:
            _drop_imported_package("numpy")
            sys.path = _prune_incompatible_python_site_paths(original_path)
            import numpy
            return numpy
        except Exception as retry_exc:
            raise RuntimeError(
                "NumPy is required for PICurv spectral generators, but no compatible "
                "NumPy could be imported. "
                f"First error: {first_error}. Retry error: {retry_exc}"
            ) from retry_exc
        finally:
            sys.path = original_path


def validate_spectral_grid(nodes, transform_axes=(0, 1, 2), subject="spectra require"):
    """!
    @brief Validate the uniform axis-aligned Cartesian grid an FFT on periodic axes requires.
    @param[in] nodes PICGRID node-coordinate array shaped `(KM, JM, IM, 3)`.
    @param[in] transform_axes Cartesian axes that are transformed and so must be uniform;
                              the remaining axes may stretch.
    @param[in] subject What needs the grid and its verb, opening each error message.
    @return Cell counts in storage order and physical lengths in Cartesian order.
    """
    numpy = require_numpy()
    km, jm, im, _ = nodes.shape
    cells = (km - 1, jm - 1, im - 1)
    if min(cells) < 4:
        raise ValueError(f"{subject} at least four cells per axis; found {cells}.")
    axes = (nodes[0, 0, :, 0], nodes[0, :, 0, 1], nodes[:, 0, 0, 2])
    for axis, (physical_axis, values) in enumerate(zip("xyz", axes)):
        delta = numpy.diff(values)
        if not numpy.all(numpy.isfinite(delta)) or not numpy.all(delta > 0) or (
                axis in transform_axes and not numpy.allclose(delta, delta[0], rtol=1e-10, atol=1e-12)):
            raise ValueError(f"{subject} uniform positive {physical_axis} spacing on every "
                             "transformed axis.")
    if not (numpy.allclose(nodes[..., 0], axes[0][None, None, :]) and
            numpy.allclose(nodes[..., 1], axes[1][None, :, None]) and
            numpy.allclose(nodes[..., 2], axes[2][:, None, None])):
        raise ValueError(f"{subject} an axis-aligned Cartesian grid (x->i, y->j, z->k).")
    return cells, tuple(float(values[-1] - values[0]) for values in axes)


def spectral_symbols(cells_kji, lengths_xyz):
    """!
    @brief Build continuum and centered-discrete Fourier symbols.
    @param[in] cells_kji Cell counts in storage order.
    @param[in] lengths_xyz Physical box lengths in Cartesian order.
    @return Continuum symbols, centered-discrete symbols, and Cartesian spacings.
    """
    numpy = require_numpy()
    nk, nj, ni = cells_kji
    lx, ly, lz = lengths_xyz
    dx, dy, dz = lx / ni, ly / nj, lz / nk
    kz = 2*numpy.pi*numpy.fft.fftfreq(nk, d=dz)[:, None, None]
    ky = 2*numpy.pi*numpy.fft.fftfreq(nj, d=dy)[None, :, None]
    kx = 2*numpy.pi*numpy.fft.fftfreq(ni, d=dx)[None, None, :]
    discrete = (numpy.sin(kx*dx)/dx, numpy.sin(ky*dy)/dy, numpy.sin(kz*dz)/dz)
    return (kx, ky, kz), discrete, (dx, dy, dz)


#: Filter kernels `fourier_resample` applies, by name.
RESAMPLE_FILTERS = ("none", "sharp", "gaussian", "box")


def _resample_axis(values, storage_axis, target_count, length, shift):
    """!
    @brief Move one periodic axis of an array of samples to another uniform resolution.
    @details The samples are read as a truncated Fourier series, which is re-evaluated at
             the target points. Modes are kept only where both resolutions represent them
             unambiguously, |m| < min(N, M)/2, so each Nyquist mode is dropped and the
             result stays real.
    @param[in] values Real samples; `storage_axis` is the axis being moved.
    @param[in] storage_axis Array axis to resample.
    @param[in] target_count Number of target samples along that axis.
    @param[in] length Physical period of the axis.
    @param[in] shift Target first-sample position minus source first-sample position.
    @return Real samples with `target_count` points along `storage_axis`.
    """
    numpy = require_numpy()
    source_count = values.shape[storage_axis]
    coeff = numpy.fft.fft(values, axis=storage_axis) / source_count
    modes = numpy.rint(numpy.fft.fftfreq(source_count) * source_count).astype(int)
    keep = numpy.abs(modes) < min(source_count, target_count) / 2.0
    phase = numpy.exp(1j * 2.0 * numpy.pi * modes[keep] * shift / length)
    shape = [1] * values.ndim
    shape[storage_axis] = int(keep.sum())
    kept = numpy.compress(keep, coeff, axis=storage_axis) * phase.reshape(shape)
    out_shape = list(values.shape)
    out_shape[storage_axis] = target_count
    target = numpy.zeros(out_shape, dtype=complex)
    index = [slice(None)] * values.ndim
    index[storage_axis] = modes[keep] % target_count
    target[tuple(index)] = kept
    return (numpy.fft.ifft(target, axis=storage_axis) * target_count).real


def filter_transfer(continuum_symbols, filter_spec):
    """!
    @brief Transfer function of a resampling filter on the target grid's wavenumbers.
    @param[in] continuum_symbols Continuum `(kx, ky, kz)` from `spectral_symbols()`.
    @param[in] filter_spec Mapping with `type` in `RESAMPLE_FILTERS`; `sharp` takes
                           `cutoff` (a wavenumber), `gaussian` and `box` take `width`
                           (the filter width Delta, a length).
    @return Broadcastable transfer array, or None for `none`.
    """
    numpy = require_numpy()
    kind = filter_spec.get("type", "none")
    kx, ky, kz = continuum_symbols
    if kind == "none":
        return None
    if kind == "sharp":
        return (kx*kx + ky*ky + kz*kz <= float(filter_spec["cutoff"])**2).astype(float)
    width = float(filter_spec["width"])
    if kind == "gaussian":
        return numpy.exp(-(kx*kx + ky*ky + kz*kz) * width * width / 24.0)
    if kind == "box":
        return numpy.sinc(kx*width/(2*numpy.pi)) * numpy.sinc(ky*width/(2*numpy.pi)) * \
            numpy.sinc(kz*width/(2*numpy.pi))
    raise ValueError(f"filter type must be one of {RESAMPLE_FILTERS}; got {kind!r}.")


def fourier_resample(values, target_cells_kji, lengths_xyz, source_first_xyz, target_first_xyz,
                     filter_spec=None):
    """!
    @brief Resample a triply periodic vector field onto another uniform grid, optionally filtered.
    @details Exact for band-limited data: every mode both grids resolve is carried over with
             its amplitude and phase, then the filter multiplies its transfer function. Each
             axis is moved separately, so memory stays at one array of the larger size.
    @param[in] values Real samples `[k, j, i, component]` of the source field.
    @param[in] target_cells_kji Target sample counts in storage order.
    @param[in] lengths_xyz Physical periods, shared by source and target, in Cartesian order.
    @param[in] source_first_xyz Position of the first source sample on each axis, relative
                                to the box origin.
    @param[in] target_first_xyz Position of the first target sample on each axis.
    @param[in] filter_spec Optional filter mapping for `filter_transfer()`.
    @return Real target samples `[k, j, i, component]`.
    """
    numpy = require_numpy()
    out = numpy.asarray(values, dtype=float)
    for axis in range(3):
        storage_axis = 2 - axis
        out = _resample_axis(out, storage_axis, int(target_cells_kji[storage_axis]),
                             float(lengths_xyz[axis]),
                             float(target_first_xyz[axis]) - float(source_first_xyz[axis]))
    if filter_spec and filter_spec.get("type", "none") != "none":
        continuum, _discrete, _spacing = spectral_symbols(tuple(target_cells_kji), lengths_xyz)
        transfer = filter_transfer(continuum, filter_spec)
        coeff = numpy.fft.fftn(out, axes=(0, 1, 2))
        out = numpy.fft.ifftn(coeff * transfer[..., None], axes=(0, 1, 2)).real
    return out
