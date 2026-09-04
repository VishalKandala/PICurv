"""Shared grid-generator stage tests: axis construction, placement, and transforms."""

import importlib.machinery
import importlib.util
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def load_grid_generator():
    """! @brief Load the canonical grid generator. @return Loaded module. """
    path = ROOT / "generators" / "grid.gen"
    loader = importlib.machinery.SourceFileLoader("picurv_grid_generator_tests", str(path))
    spec = importlib.util.spec_from_loader("picurv_grid_generator_tests", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


GRID = load_grid_generator()


def unit_box(ncells=4):
    """! @brief Build an unwarped Cartesian block. @param[in] ncells Cells per axis.
    @return Coordinates and realized stretch factors. """
    return GRID.generate_curvilinear_grid(
        ncells_i=ncells, ncells_j=ncells, ncells_k=ncells,
        x_min=0.0, x_max=1.0, y_min=0.0, y_max=2.0, z_min=0.0, z_max=3.0,
        A=0.0, B=0.0, C=0.0, origin=[0.0, 0.0, 0.0],
        stretch_i=0.0, stretch_j=0.0, stretch_k=0.0,
    )


def bounds(X, Y, Z):
    """! @brief Axis-aligned bounding box. @param[in] X x coordinates. @param[in] Y y
    coordinates. @param[in] Z z coordinates. @return Per-axis (min, max) pairs. """
    return [(float(a.min()), float(a.max())) for a in (X, Y, Z)]


def test_single_node_axis_returns_a_stretch_factor():
    """! @brief A one-node axis must not break the caller's unpack. """
    values, factor = GRID.create_stretched_linspace(1, 0.0, 2.0, 2.0)
    assert values.tolist() == [2.0]
    assert factor == 0.0


def test_default_placement_puts_the_bounding_box_corner_on_origin():
    """! @brief Placement keeps the historical `origin` meaning. """
    X, Y, Z, _ = GRID.generate_curvilinear_grid(
        ncells_i=3, ncells_j=3, ncells_k=3,
        x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, z_min=0.0, z_max=1.0,
        A=0.0, B=0.0, C=0.0, origin=[5.0, -2.0, 0.5],
        stretch_i=0.0, stretch_j=0.0, stretch_k=0.0,
    )
    assert bounds(X, Y, Z) == [(5.0, 6.0), (-2.0, -1.0), (0.5, 1.5)]


def test_geometries_report_the_stretch_factors_the_report_consumes():
    """! @brief Realized factors travel with the grid, not through module state. """
    _, _, _, factors = GRID.generate_curvilinear_grid(
        ncells_i=8, ncells_j=8, ncells_k=8,
        x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, z_min=0.0, z_max=1.0,
        A=0.0, B=0.0, C=0.0, origin=[0.0, 0.0, 0.0],
        stretch_i=0.0, stretch_j=2.5, stretch_k=0.0,
    )
    assert factors['j'] == pytest.approx(2.5)
    assert factors['i'] == 0.0 and factors['k'] == 0.0

    X, Y, Z, _ = unit_box()
    stats = GRID.analyze_grid_quality(X, Y, Z, factors)
    assert stats['computed_stretch_j'] == pytest.approx(2.5)


def test_transform_list_is_applied_in_the_order_given():
    """! @brief Transforms compose, and composition is not commutative. """
    X, Y, Z, _ = unit_box()
    first = GRID.apply_transforms(X, Y, Z, ["translate:dx=10", "scale:sx=2,about=origin"])
    second = GRID.apply_transforms(X, Y, Z, ["scale:sx=2,about=origin", "translate:dx=10"])
    assert bounds(*first)[0] == (20.0, 22.0)
    assert bounds(*second)[0] == (10.0, 12.0)


def test_transforms_compose_on_top_of_origin_placement():
    """! @brief A transform list refines placement rather than replacing it. """
    X, Y, Z, _ = GRID.generate_curvilinear_grid(
        ncells_i=3, ncells_j=3, ncells_k=3,
        x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, z_min=0.0, z_max=1.0,
        A=0.0, B=0.0, C=0.0, origin=[1.0, 1.0, 1.0],
        stretch_i=0.0, stretch_j=0.0, stretch_k=0.0,
        transforms=["translate:dx=4"],
    )
    assert bounds(X, Y, Z) == [(5.0, 6.0), (1.0, 2.0), (1.0, 2.0)]


@pytest.mark.parametrize("token, expected", [
    ("translate:dx=1,dy=2,dz=3", [(1.0, 2.0), (2.0, 4.0), (3.0, 6.0)]),
    ("anchor:mode=bbox_center,x=0,y=0,z=0", [(-0.5, 0.5), (-1.0, 1.0), (-1.5, 1.5)]),
    ("scale:sx=2,sy=2,sz=2,about=bbox_min", [(0.0, 2.0), (0.0, 4.0), (0.0, 6.0)]),
    ("permute:order=zyx", [(0.0, 3.0), (0.0, 2.0), (0.0, 1.0)]),
])
def test_each_transform_kind_moves_the_grid_as_declared(token, expected):
    """! @brief Every advertised transform is reachable and does what it says.
    @param[in] token Transform token under test. @param[in] expected Resulting bounds. """
    X, Y, Z, _ = unit_box()
    assert bounds(*GRID.apply_transforms(X, Y, Z, [token])) == expected


def test_rotation_about_the_box_centre_preserves_extent():
    """! @brief A quarter turn swaps the two in-plane extents. """
    X, Y, Z, _ = unit_box()
    rx, ry, rz = GRID.apply_transforms(X, Y, Z, ["rotate:axis=z,deg=90,about=bbox_center"])
    (x0, x1), (y0, y1), _ = bounds(rx, ry, rz)
    assert x1 - x0 == pytest.approx(2.0)
    assert y1 - y0 == pytest.approx(1.0)


@pytest.mark.parametrize("token, message", [
    ("bogus:x=1", "unknown transform 'bogus'"),
    ("translate:dx", "is not 'key=value'"),
    ("translate:dx=east", "must be a number"),
    ("rotate:axis=q,deg=1", "axis must be x, y or z"),
    ("rotate:deg=nope", "deg must be a number"),
    ("scale:sx=0", "factors must be non-zero"),
    ("mirror:plane=diagonal", "plane must be one of"),
    ("permute:order=xxy", "order must be one of"),
    ("anchor:mode=nowhere", "unknown reference point"),
])
def test_malformed_transforms_are_refused_with_their_reason(token, message):
    """! @brief A bad transform names what is wrong instead of silently doing nothing.
    @param[in] token Malformed token. @param[in] message Expected diagnostic fragment. """
    X, Y, Z, _ = unit_box()
    with pytest.raises(ValueError, match=message):
        GRID.apply_transforms(X, Y, Z, [token])


def test_handedness_is_reported_separately_from_sign_consistency():
    """! @brief A uniformly inverted grid is self-consistent, so consistency cannot see it. """
    X, Y, Z, factors = unit_box()
    upright = GRID.analyze_grid_quality(X, Y, Z, factors)
    assert upright['jacobian_consistent'] and upright['right_handed']

    flipped = GRID.analyze_grid_quality(
        *GRID.apply_transforms(X, Y, Z, ["permute:order=zyx"]), factors)
    assert flipped['jacobian_consistent'], "an inverted grid still agrees with itself"
    assert not flipped['right_handed']


def test_even_axis_permutation_preserves_handedness():
    """! @brief Only an odd permutation flips orientation. """
    X, Y, Z, factors = unit_box()
    cycled = GRID.analyze_grid_quality(
        *GRID.apply_transforms(X, Y, Z, ["permute:order=yzx"]), factors)
    assert cycled['right_handed']


def test_config_style_transform_tokens_parse_identically_to_cli_tokens():
    """! @brief The one grammar serves the config file and the command line. """
    kind, params = GRID.parse_transform_token("rotate:axis=z,deg=90,about=centroid")
    assert kind == "rotate"
    assert params == {"axis": "z", "deg": "90", "about": "centroid"}
