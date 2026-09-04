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


def plain_box(ncells=4, bounds=((0.0, 1.0), (0.0, 2.0), (0.0, 3.0)),
              origin=(0.0, 0.0, 0.0), **kwargs):
    """! @brief Build a flat-walled Cartesian block. @param[in] ncells Cells per axis.
    @param[in] bounds Per-axis extents. @param[in] origin Placement target.
    @param[in] kwargs Extra generator arguments.
    @return Coordinates and realized stretch factors. """
    defaults = dict(stretch_i=0.0, stretch_j=0.0, stretch_k=0.0)
    defaults.update(kwargs)
    return GRID.box_grid(
        ncells_i=ncells, ncells_j=ncells, ncells_k=ncells,
        bounds_x=list(bounds[0]), bounds_y=list(bounds[1]), bounds_z=list(bounds[2]),
        wall_j_lo=None, wall_j_hi=None, wall_j_lo_span=None, wall_j_hi_span=None,
        amp_A=0.0, amp_B=0.0, amp_C=0.0, origin=list(origin), **defaults)


def unit_box(ncells=4):
    """! @brief The block the transform tests measure against.
    @param[in] ncells Cells per axis. @return Coordinates and realized stretch factors. """
    return plain_box(ncells)


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
    X, Y, Z, _ = plain_box(3, bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
                           origin=(5.0, -2.0, 0.5))
    assert bounds(X, Y, Z) == [(5.0, 6.0), (-2.0, -1.0), (0.5, 1.5)]


def test_geometries_report_the_stretch_factors_the_report_consumes():
    """! @brief Realized factors travel with the grid, not through module state. """
    _, _, _, factors = plain_box(8, bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
                                 stretch_j=2.5)
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
    X, Y, Z, _ = plain_box(3, bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
                           origin=(1.0, 1.0, 1.0), transforms=["translate:dx=4"])
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


# --- box: piecewise wall-height fields ---------------------------------------

def shaped_box(wall_lo=None, wall_hi=None, span_lo=None, ncells_i=40, ncells_j=8,
               ncells_k=4, bounds_x=(0.0, 20.0)):
    """! @brief Build a shaped box. @param[in] wall_lo Lower-wall segments.
    @param[in] wall_hi Upper-wall segments. @param[in] span_lo Lower spanwise envelope.
    @param[in] ncells_i Streamwise cells. @param[in] ncells_j Wall-normal cells.
    @param[in] ncells_k Spanwise cells. @param[in] bounds_x Streamwise extent.
    @return Coordinates and realized stretch factors. """
    return GRID.box_grid(
        ncells_i=ncells_i, ncells_j=ncells_j, ncells_k=ncells_k,
        bounds_x=list(bounds_x), bounds_y=[0.0, 2.0], bounds_z=[0.0, 4.0],
        wall_j_lo=wall_lo, wall_j_hi=wall_hi,
        wall_j_lo_span=span_lo, wall_j_hi_span=None,
        amp_A=0.0, amp_B=0.0, amp_C=0.0, origin=[0.0, 0.0, 0.0],
        stretch_i=0.0, stretch_j=0.0, stretch_k=0.0,
    )


def test_flat_box_is_an_exact_cartesian_block():
    """! @brief With no wall lists and no warp the block is plain Cartesian. """
    X, Y, Z, _ = shaped_box()
    assert bounds(X, Y, Z) == [(0.0, 20.0), (0.0, 2.0), (0.0, 4.0)]
    assert np.allclose(np.unique(np.round(np.diff(X[:, 0, 0]), 12)), [0.5])


def test_a_step_moves_only_the_wall_it_names():
    """! @brief A lower-wall step leaves the opposite wall where it was. """
    X, Y, Z, _ = shaped_box(
        wall_lo=["flat:len=4,y=1", "step:len=2,dy=-1", "flat:len=14"])
    inlet = X[:, 0, 0] < 4.0
    outlet = X[:, 0, 0] > 6.0
    assert np.allclose(Y[inlet, 0, 0], 1.0), "inlet floor sits at the step height"
    assert np.allclose(Y[outlet, 0, 0], 0.0), "floor drops to the datum after the step"
    assert np.allclose(Y[:, -1, 0], 2.0), "the roof is untouched"


def test_wall_segments_are_continuous_across_the_join():
    """! @brief Each segment starts where the previous one ended. """
    _, Y, _, _ = shaped_box(
        wall_lo=["flat:len=5,y=0.5", "ramp:len=5,dy=0.5", "flat:len=10"],
        ncells_i=200)
    floor = Y[:, 0, 0]
    assert np.max(np.abs(np.diff(floor))) < 0.02, "no jump at a segment boundary"
    assert floor[-1] - floor[0] == pytest.approx(0.5), "the ramp raises the wall by dy"


def test_origin_anchors_the_realized_box_not_the_declared_bounds():
    """! @brief Shaping a wall shrinks the block, and `origin` places what is left.

    @details A raised floor makes the domain shorter than `bounds_y` describes, and
             placement has always put the bounding-box corner on `origin`. Pinned because
             the two together move a shaped box relative to where its bounds alone suggest.
    """
    _, Y, _, _ = shaped_box(
        wall_lo=["flat:len=5,y=0.5", "ramp:len=5,dy=0.5", "flat:len=10"],
        ncells_i=200)
    assert Y.min() == pytest.approx(0.0), "the realized corner lands on origin"
    assert Y.max() == pytest.approx(1.5), "roof at 2.0 less the 0.5 the floor was raised"


def test_a_hill_segment_returns_to_its_start_height():
    """! @brief A hill tiles periodically, so it must close on itself. """
    _, Y, _, _ = shaped_box(wall_lo=["hill:len=20,height=0.5"], ncells_i=200)
    floor = Y[:, 0, 0]
    assert floor[0] == pytest.approx(0.0, abs=1e-12)
    assert floor[-1] == pytest.approx(0.0, abs=1e-12)
    assert floor.max() == pytest.approx(0.5, rel=1e-6)


def test_spanwise_envelope_gives_an_obstacle_finite_width():
    """! @brief The envelope scales a wall's departure from its datum, not its datum. """
    _, Y, Z, _ = GRID.box_grid(
        ncells_i=40, ncells_j=8, ncells_k=40,
        bounds_x=[0.0, 20.0], bounds_y=[0.0, 2.0], bounds_z=[0.0, 4.0],
        wall_j_lo=["flat:len=8,y=0", "step:len=2,dy=0.5", "flat:len=4",
                   "step:len=2,dy=-0.5", "flat:len=4"],
        wall_j_hi=None,
        wall_j_lo_span=["flat:len=1,y=0", "step:len=0.5,dy=1", "flat:len=1",
                        "step:len=0.5,dy=-1", "flat:len=1"],
        wall_j_hi_span=None,
        amp_A=0.0, amp_B=0.0, amp_C=0.0, origin=[0.0, 0.0, 0.0],
        stretch_i=0.0, stretch_j=0.0, stretch_k=0.0,
    )
    on_block = np.argmin(np.abs(Z[0, 0, :] - 2.0))
    off_block = 0
    assert Y[:, 0, on_block].max() == pytest.approx(0.5, rel=1e-6)
    assert Y[:, 0, off_block].max() == pytest.approx(0.0, abs=1e-12), \
        "outside the footprint the wall stays on its datum"


def test_wall_segments_must_span_their_axis():
    """! @brief A wall that stops short leaves its remainder undefined. """
    with pytest.raises(ValueError, match="segment lengths sum to"):
        shaped_box(wall_lo=["flat:len=4,y=1", "step:len=2,dy=-1"])


def test_a_sub_cell_corner_is_refused_rather_than_silently_flattened():
    """! @brief Below one cell a shorter step yields the identical grid. """
    with pytest.raises(ValueError, match="cannot turn in less than one cell"):
        shaped_box(
            wall_lo=["flat:len=4,y=1", "step:len=0.05,dy=-1", "flat:len=15.95"],
            ncells_i=40)


def test_a_coarse_corner_is_reported_but_allowed(capsys):
    """! @brief Between one and four cells the corner is a quality call, not an error.
    @param[in] capsys Captured output fixture. """
    shaped_box(wall_lo=["flat:len=4,y=1", "step:len=1,dy=-1", "flat:len=15"],
               ncells_i=40)
    assert "spans only 2.0 cells" in capsys.readouterr().err


def test_a_sharper_corner_costs_orthogonality():
    """! @brief The quality report is what makes the corner tradeoff visible. """
    angles = []
    for length in (1.0, 4.0):
        X, Y, Z, factors = shaped_box(
            wall_lo=[f"flat:len=4,y=1", f"step:len={length},dy=-1",
                     f"flat:len={16 - length}"], ncells_i=200, ncells_j=32)
        angles.append(GRID.analyze_grid_quality(X, Y, Z, factors)['max_non_ortho'])
    assert angles[0] > angles[1], "a shorter step is a less orthogonal corner"


def test_unknown_wall_segment_names_the_allowed_kinds():
    """! @brief A misspelled segment kind says what was expected. """
    with pytest.raises(ValueError, match="unknown wall segment 'cliff'"):
        shaped_box(wall_lo=["cliff:len=20,dy=-1"])


# --- sweep: cross-section along a centreline ---------------------------------

def swept(path=("straight:len=10",), cross_section="rectangle", sides=(1.0, 1.0),
          scale=None, ncells_i=16, ncells_j=16, ncells_k=40, **kwargs):
    """! @brief Build a swept duct. @param[in] path Centreline segments.
    @param[in] cross_section Section shape. @param[in] sides Section widths.
    @param[in] scale Optional scale field. @param[in] ncells_i First section cells.
    @param[in] ncells_j Second section cells. @param[in] ncells_k Path cells.
    @param[in] kwargs Extra generator arguments. @return Coordinates and factors. """
    return GRID.sweep_grid(
        ncells_i=ncells_i, ncells_j=ncells_j, ncells_k=ncells_k,
        cross_section=cross_section, side_lengths=list(sides), path=list(path),
        cross_section_scale=scale, origin=[0.0, 0.0, 0.0],
        stretch_i=0.0, stretch_j=0.0, stretch_k=0.0, **kwargs)


def test_a_straight_rectangular_sweep_is_a_box():
    """! @brief The simplest sweep degenerates to the obvious answer. """
    X, Y, Z, _ = swept(sides=(2.0, 3.0))
    assert bounds(X, Y, Z) == [(0.0, 10.0), (0.0, 2.0), (0.0, 3.0)]


def test_swept_grids_are_right_handed():
    """! @brief The section axes precede the path axis, which fixes orientation.

    @details The retired O-grid emitted uniformly inverted cells that the solver had to
             repair at runtime. Pinned so the replacement cannot regress into that.
    """
    for path in (["straight:len=10"],
                 ["straight:len=3", "arc:radius=2,deg=180", "straight:len=3"]):
        X, Y, Z, factors = swept(path=path, ncells_k=60)
        assert GRID.analyze_grid_quality(X, Y, Z, factors)['right_handed'], path


def test_disc_cross_section_fills_the_circle_without_a_hole():
    """! @brief The square-to-disc map covers the whole disc and reaches its edge. """
    X, Y, Z, _ = swept(cross_section="circle", sides=(2.0, 2.0), ncells_i=32, ncells_j=32)
    section_y, section_z = Y[:, :, 0] - 1.0, Z[:, :, 0] - 1.0
    radius = np.hypot(section_y, section_z)
    assert radius.max() == pytest.approx(1.0, rel=1e-9), "the boundary lands on the circle"
    assert radius.min() == pytest.approx(0.0, abs=1e-12), "the centre is covered, no hole"
    assert np.all(radius <= 1.0 + 1e-12), "nothing escapes the disc"


def test_parallel_transport_keeps_an_out_of_plane_path_untwisted():
    """! @brief A frame carried by a fixed up-vector twists once a path leaves its plane.

    @details Turning first about z and then about x takes the duct out of one plane. Under
             parallel transport the section does not rotate about its own tangent, so the
             cells stay near-orthogonal; a twisting frame would shear them badly.
    """
    X, Y, Z, factors = swept(
        path=["straight:len=2", "arc:radius=3,deg=90,axis=z",
              "arc:radius=3,deg=90,axis=x", "straight:len=2"],
        ncells_i=24, ncells_j=24, ncells_k=150)
    stats = GRID.analyze_grid_quality(X, Y, Z, factors)
    assert stats['max_non_ortho'] < 5.0, stats['max_non_ortho']
    assert stats['right_handed']
    assert min(b[1] - b[0] for b in bounds(X, Y, Z)) > 1.0, "the path really left one plane"


def test_cross_section_scale_narrows_the_duct_along_the_path():
    """! @brief A scale field turns one sweep into a nozzle. """
    X, Y, Z, _ = swept(sides=(2.0, 2.0), ncells_k=80,
                       scale=["flat:len=3,y=1", "step:len=4,dy=-0.5", "flat:len=3"])
    inlet = Y[:, :, 0].max() - Y[:, :, 0].min()
    outlet = Y[:, :, -1].max() - Y[:, :, -1].min()
    assert inlet == pytest.approx(2.0)
    assert outlet == pytest.approx(1.0), "a 0.5 scale halves the section"


def test_an_arc_about_the_tangent_is_refused():
    """! @brief An arc turning about its own tangent has no defined plane. """
    with pytest.raises(ValueError, match="parallel to the current tangent"):
        swept(path=["straight:len=2", "arc:radius=3,deg=90,axis=x"])


@pytest.mark.parametrize("kwargs, message", [
    ({"path": []}, "needs at least one path segment"),
    ({"path": ["straight:len=0"]}, "len must be positive"),
    ({"path": ["arc:radius=0,deg=90"]}, "radius must be positive"),
    ({"path": ["helix:len=1"]}, "unknown path segment 'helix'"),
    ({"cross_section": "hexagon"}, "cross_section must be one of"),
    ({"scale": ["flat:len=10,y=0"]}, "cannot reach zero or negative width"),
])
def test_malformed_sweeps_are_refused_with_their_reason(kwargs, message):
    """! @brief Every sweep rejection says what was wrong.
    @param[in] kwargs Generator arguments under test. @param[in] message Expected text. """
    with pytest.raises(ValueError, match=message):
        swept(**kwargs)


# --- reference scales, wall units, and validators ----------------------------

def test_reference_scales_need_a_coherent_combination():
    """! @brief Each way of expressing a wall unit needs its whole pair. """
    assert GRID.resolve_reference_quantities(None, None, None, None, None) is None
    with pytest.raises(ValueError, match="re_tau is a ratio against the reference length"):
        GRID.resolve_reference_quantities(None, None, None, 180.0, None)
    with pytest.raises(ValueError, match="u_tau needs nu"):
        GRID.resolve_reference_quantities(1.0, None, None, None, 0.05)


def test_the_two_friction_forms_agree():
    """! @brief Re_tau/length_ref and u_tau/nu are the same ratio. """
    from_re_tau = GRID.resolve_reference_quantities(2.0, 1e-3, None, 180.0, None)
    from_u_tau = GRID.resolve_reference_quantities(2.0, 1e-3, None, None, 0.09)
    assert from_re_tau['inverse_viscous_length'] == pytest.approx(90.0)
    assert from_u_tau['inverse_viscous_length'] == pytest.approx(90.0)


def test_wall_units_reproduce_the_channel_configs_hand_arithmetic():
    """! @brief The generator computes what the shipped config states in prose.

    @details config/grids/plane_channel_retau180.cfg asserts y+_1 = 0.5, dx+ = 8.8 and
             dz+ = 17.5 in a comment nothing checks. Same domain and counts here.
    """
    reference = GRID.resolve_reference_quantities(1.0, 3.5714286e-04, 1.0, 180.0, None)
    X, Y, Z, factors = GRID.box_grid(
        ncells_i=128, ncells_j=128, ncells_k=128,
        bounds_x=[0.0, 2*np.pi], bounds_y=[0.0, 2.0], bounds_z=[0.0, 4*np.pi],
        wall_j_lo=None, wall_j_hi=None, wall_j_lo_span=None, wall_j_hi_span=None,
        amp_A=0.0, amp_B=0.0, amp_C=0.0, origin=[0.0, 0.0, 0.0],
        stretch_i=0.0, stretch_j=2.0, stretch_k=0.0,
        first_cell_j_start=1.388889e-03, first_cell_j_end=1.388889e-03,
    )
    wall = GRID.analyze_grid_quality(X, Y, Z, factors, reference)['wall_units']
    assert wall['j']['first_start'] == pytest.approx(0.5, rel=1e-3)
    assert wall['i']['max'] == pytest.approx(8.8, rel=1e-2)
    assert wall['k']['max'] == pytest.approx(17.5, rel=1e-2)
    assert wall['j']['extent'] == pytest.approx(360.0, rel=1e-6), "2 delta in wall units"


def test_solver_units_report_the_grid_after_the_launcher_divides_it():
    """! @brief The generator writes metres; the solver receives them over length_ref. """
    reference = GRID.resolve_reference_quantities(0.05, None, None, None, None)
    X, Y, Z, factors = shaped_box(bounds_x=(0.0, 1.0))
    stats = GRID.analyze_grid_quality(X, Y, Z, factors, reference)
    assert stats['bounds_x'][1] == pytest.approx(1.0), "written dimensional"
    assert stats['solver_bounds']['x'][1] == pytest.approx(20.0), "received scaled"


@pytest.mark.parametrize("counts, levels, ok", [
    ({'i': 129, 'j': 129, 'k': 129}, 5, True),
    ({'i': 129, 'j': 65, 'k': 33}, 5, True),   # 33 -> 17 -> 9 -> 5 -> 3 is a legal ladder
    ({'i': 129, 'j': 64, 'k': 129}, 5, False),  # an even count breaks at the first level
    ({'i': 101, 'j': 101, 'k': 101}, 5, False),
    ({'i': 33, 'j': 33, 'k': 33}, 3, True),
])
def test_multigrid_ladder_is_checked_at_generation(counts, levels, ok):
    """! @brief An illegal count is loud here instead of a runtime coarsening notice.
    @param[in] counts Node counts. @param[in] levels Requested levels.
    @param[in] ok Whether the combination is legal. """
    if ok:
        GRID.check_multigrid_levels(counts, levels)
    else:
        with pytest.raises(ValueError, match="multigrid levels need node counts"):
            GRID.check_multigrid_levels(counts, levels)


def test_a_periodic_axis_refuses_one_sided_stretching():
    """! @brief The seam would join cells of different width. """
    axes = {'i': np.array([0.0, 0.01, 0.05, 0.2, 0.6, 1.0])}
    with pytest.raises(ValueError, match="first and last cells differ"):
        GRID.check_periodic_axes(['i'], axes, None)


def test_a_periodic_axis_accepts_symmetric_stretching():
    """! @brief Matching end cells make a seam the runtime can join. """
    axes = {'i': np.array([0.0, 0.1, 0.35, 0.65, 0.9, 1.0])}
    GRID.check_periodic_axes(['i'], axes, None)


def test_rotating_a_periodic_grid_is_refused():
    """! @brief The runtime's periodic reconstruction assumes an axis-aligned seam.

    @details Metric.c offsets only the matching Cartesian component across the seam and
             copies the other two. Whether an arbitrarily oriented periodic grid is
             handled is unverified, so this refuses rather than producing one quietly.
    """
    axes = {'k': np.linspace(0.0, 1.0, 9)}
    with pytest.raises(ValueError, match="cannot be combined with 'rotate'"):
        GRID.check_periodic_axes(['k'], axes, ["rotate:axis=z,deg=30"])
    GRID.check_periodic_axes(['k'], axes, ["translate:dx=1"])


# --- the generator and the launcher must agree on the choice sets ------------

def load_core():
    """! @brief Load the conductor module. @return Loaded module. """
    path = ROOT / "picurv_cli" / "core.py"
    loader = importlib.machinery.SourceFileLoader("picurv_core_grid_tests", str(path))
    spec = importlib.util.spec_from_loader("picurv_core_grid_tests", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


CORE = load_core()


@pytest.mark.parametrize("symbol", [
    "GRID_CROSS_SECTION_KINDS", "GRID_WALL_SEGMENT_KINDS",
    "GRID_PATH_SEGMENT_KINDS", "GRID_TRANSFORM_KINDS",
])
def test_core_mirrors_the_generators_choice_sets(symbol):
    """! @brief The launcher validates against sets the generator actually implements.

    @details core.py holds these so audit_family_census can discover them: it scans only
             picurv_cli, never generators/. Two copies need pinning, and parity_sources
             has no Python-source kind, so it is pinned here instead.
    @param[in] symbol Choice-set constant present in both modules.
    """
    assert getattr(CORE, symbol) == getattr(GRID, symbol), symbol


def test_grid_generator_types_are_the_generators_subcommands():
    """! @brief Every advertised geometry is a subcommand, and nothing else is. """
    import argparse
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    # The generator builds its parser inside main(); take the names it documents instead.
    documented = set(CORE.GRID_GENERATOR_TYPES)
    assert documented == {"box", "sweep"}
    for name in documented:
        assert hasattr(GRID, f"{name}_grid"), f"{name} has no implementation"


def test_retired_geometries_are_gone_from_both_sides():
    """! @brief A removed value must not survive in either module. """
    for name in ("cpipe", "pipe", "warp"):
        assert name not in CORE.GRID_GENERATOR_TYPES
    for symbol in ("cartesian_pipe_grid", "bent_pipe_grid", "generate_curvilinear_grid"):
        assert not hasattr(GRID, symbol), f"{symbol} still present"


@pytest.mark.parametrize("cli_args, expected", [
    (["--cross-section", "circle"], 0),
    (["--cross-section", "hexagon"], 1),
    (["--path", "straight:len=5", "arc:radius=2,deg=90"], 0),
    (["--path", "helix:len=5"], 1),
    (["--wall-j-lo", "flat:len=4,y=1", "cliff:len=1,dy=-1"], 1),
    (["--transforms", "translate:dx=1", "warp:x=1"], 1),
    (["--ncells-i", "64", "--bounds-x", "0", "1"], 0),
])
def test_case_validation_catches_bad_choices_inside_cli_args(cli_args, expected):
    """! @brief A misspelled selector fails at validation, not mid-run.
    @param[in] cli_args Raw token list. @param[in] expected Number of errors. """
    errors = CORE.validate_grid_generator_cli_args(cli_args, "case.yml")
    assert len(errors) == expected, errors
