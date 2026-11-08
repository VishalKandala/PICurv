"""!
@file picgrid.py
@brief The canonical PICGRID text grid format: its coordinate precision, header reader and
       coordinate reader, shared by the conductor and the generators.

Not a command. `core.py`, `grid.gen`, `ic.gen`, `profile.gen`, `reduce.gen` and `spectra.gen`
import it from `picurv_cli/`, as they import `petsc_binary.py`. A PICGRID file is the token
`PICGRID`, a block count, one `IM JM KM` node-count row per block, then every block's node
coordinates in block order, one `x y z` row per node with i fastest, then j, then k. Blank lines
and lines beginning with `#` are not part of the content.
"""

import math
import os

#: Coordinate text precision; 17 significant digits round-trip any binary64 value.
PICGRID_FLOAT_FORMAT = ".17e"


def format_picgrid_coordinate(value: float) -> str:
    """!
    @brief Format a coordinate with round-trip-safe binary64 precision.
    @param[in] value Coordinate value.
    @return Formatted coordinate.
    """
    return format(value, PICGRID_FLOAT_FORMAT)


def content_lines(stream):
    """!
    @brief Yield the content lines of a PICurv text payload, skipping blank and `#` lines.
    @param[in] stream Open text stream.
    @return Generator of `(line_number, stripped_line)` pairs, numbered from 1.
    """
    for lineno, raw in enumerate(stream, start=1):
        line = raw.strip()
        if line and not line.startswith("#"):
            yield lineno, line


def read_picgrid_header(lines, path: str, expected_nblk: int = None) -> list:
    """!
    @brief Consume a PICGRID header from a content-line iterator.
    @details Leaves the iterator at the first coordinate row.
    @param[in] lines Iterator from content_lines().
    @param[in] path File path, for errors.
    @param[in] expected_nblk Optional block count the caller requires.
    @return List of `(IM, JM, KM)` node-count tuples, one per block.
    @throws ValueError on a malformed header.
    """
    try:
        _, first_token = next(lines)
    except StopIteration:
        raise ValueError(f"Grid file '{path}' is empty.")
    if first_token != "PICGRID":
        raise ValueError(f"Grid file '{path}' must begin with the canonical PICGRID header token.")

    try:
        _, nblk_line = next(lines)
        nblk = int(nblk_line)
    except StopIteration:
        raise ValueError(f"Grid file '{path}' missing block count after PICGRID header.")
    except ValueError:
        raise ValueError(f"Invalid block count '{nblk_line}' in grid file '{path}'.")
    if nblk <= 0:
        raise ValueError(f"Grid file '{path}' has non-positive block count: {nblk}.")
    if expected_nblk is not None and nblk != expected_nblk:
        raise ValueError(f"Grid file block count mismatch: case expects {expected_nblk}, grid contains {nblk}.")

    dims = []
    for bi in range(nblk):
        try:
            lineno, dim_line = next(lines)
        except StopIteration:
            raise ValueError(f"Grid file '{path}' missing dimensions for block {bi}.")
        parts = dim_line.split()
        if len(parts) != 3:
            raise ValueError(f"Invalid dimensions line at {path}:{lineno}. Expected 3 integers, got: '{dim_line}'.")
        try:
            im, jm, km = (int(parts[0]), int(parts[1]), int(parts[2]))
        except ValueError:
            raise ValueError(f"Invalid dimensions line at {path}:{lineno}. Non-integer values: '{dim_line}'.")
        if im <= 0 or jm <= 0 or km <= 0:
            raise ValueError(f"Invalid block dimensions at {path}:{lineno}: ({im}, {jm}, {km}). Must be > 0.")
        dims.append((im, jm, km))
    return dims


def read_picgrid_dimensions(path: str, expected_nblk: int = None) -> list:
    """!
    @brief Read only the header of a PICGRID file.
    @param[in] path PICGRID path.
    @param[in] expected_nblk Optional block count the caller requires.
    @return List of `(IM, JM, KM)` node-count tuples, one per block.
    @throws ValueError on a missing file or malformed header.
    """
    if not os.path.isfile(path):
        raise ValueError(f"Grid file not found: {path}")
    with open(path, "r", encoding="utf-8") as stream:
        return read_picgrid_header(content_lines(stream), path, expected_nblk)


def parse_coordinate_row(lineno: int, line: str, path: str) -> tuple:
    """!
    @brief Parse one coordinate row into three finite floats.
    @param[in] lineno Line number, for errors.
    @param[in] line Stripped content line.
    @param[in] path File path, for errors.
    @return `(x, y, z)`.
    @throws ValueError when the row is not three finite numbers.
    """
    parts = line.split()
    if len(parts) != 3:
        raise ValueError(f"Invalid coordinate row at {path}:{lineno}. Expected 3 floats, got: '{line}'.")
    try:
        point = (float(parts[0]), float(parts[1]), float(parts[2]))
    except ValueError:
        raise ValueError(f"Invalid coordinate row at {path}:{lineno}. Non-numeric values: '{line}'.")
    if not all(math.isfinite(value) for value in point):
        raise ValueError(f"Invalid coordinate row at {path}:{lineno}. Coordinates must be finite: '{line}'.")
    return point


def read_picgrid_blocks(path: str, expected_nblk: int = None) -> list:
    """!
    @brief Read every block of a PICGRID file as node-coordinate arrays.
    @details For consumers that build cell geometry, so each axis needs at least two nodes.
    @param[in] path PICGRID path.
    @param[in] expected_nblk Optional block count the caller requires.
    @return List of `{"dims": (IM, JM, KM), "coords": array shaped (KM, JM, IM, 3)}`.
    @throws ValueError on a missing or malformed file.
    """
    from periodic_spectral import require_numpy

    numpy = require_numpy()
    if not os.path.isfile(path):
        raise ValueError(f"Grid file not found: {path}")
    with open(path, "r", encoding="utf-8") as stream:
        lines = content_lines(stream)
        dims = read_picgrid_header(lines, path, expected_nblk)
        blocks = []
        for block_idx, (im, jm, km) in enumerate(dims):
            if min(im, jm, km) < 2:
                raise ValueError(f"PICGRID block {block_idx} dimensions {(im, jm, km)} are too small.")
            coords = []
            for lineno, line in lines:
                coords.append(parse_coordinate_row(lineno, line, path))
                if len(coords) == im * jm * km:
                    break
            if len(coords) != im * jm * km:
                raise ValueError(f"Grid file '{path}' missing coordinates for block {block_idx}.")
            blocks.append({"dims": (im, jm, km),
                           "coords": numpy.asarray(coords, dtype=float).reshape((km, jm, im, 3))})
        if next(lines, None) is not None:
            raise ValueError(f"Grid file '{path}' has trailing coordinate rows.")
    return blocks
