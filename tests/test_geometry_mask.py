"""objdel.geometry.to_mask: parts -> even-odd mask + legality, on the CUDA checker (to_masks of one object)."""
import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st
from oracles import raster_oracle as RO
from oracles import validity_oracle as VO

from objdel.geometry import GRID, canonicalize, to_mask

SQUARE = np.array([[0, 0], [1000, 0], [1000, 1000], [0, 1000]])  # clockwise on screen
HOLE = np.array([[250, 250], [250, 750], [750, 750], [750, 250]])  # counter-clockwise


def test_full_square_with_a_hole_on_a_4x4_crop():
    g = to_mask([(SQUARE, [HOLE])], 4, 4)
    # pixel centres at 125, 375, 625, 875: the hole (250..750) removes the inner 2x2 block
    expected = np.ones((4, 4), bool)
    expected[1:3, 1:3] = False
    assert g.legal and g.defects == {}
    assert g.mask.dtype == bool and np.array_equal(g.mask, expected)


def test_mask_shape_is_rows_by_columns():
    g = to_mask([(SQUARE, [])], 5, 3)
    assert g.mask.shape == (3, 5) and g.mask.all()


def test_no_parts_is_an_empty_legal_mask():
    g = to_mask([], 7, 5)
    assert g.mask.shape == (5, 7) and g.mask.dtype == bool and not g.mask.any()
    assert g.legal is True and g.defects == {}


def test_roles_hole_belongs_to_its_own_outer():
    """Two parts, each with a hole: holes are checked against their own outer (roles -1, 0, -1, 2)."""
    outer = np.array([[0, 0], [400, 0], [400, 400], [0, 400]])
    a = (outer, [np.array([[100, 100], [100, 300], [300, 300], [300, 100]])])
    b = (a[0] + 500, [a[1][0] + 500])
    assert to_mask([a, b], 10, 10).legal
    swapped = [(a[0], b[1]), (b[0], a[1])]  # each hole now lies outside its outer
    g = to_mask(swapped, 10, 10)
    assert not g.legal and g.defects.get("hole_outside_outer") == 2


def test_winding_is_part_of_legality():
    assert to_mask([(SQUARE[::-1], [])], 4, 4).defects == {"wrong_winding": 1}
    assert to_mask([(SQUARE, [HOLE[::-1]])], 4, 4).defects == {"wrong_winding": 1}
    # the mask itself is the even-odd fill and does not depend on direction
    assert np.array_equal(to_mask([(SQUARE[::-1], [HOLE[::-1]])], 4, 4).mask, to_mask([(SQUARE, [HOLE])], 4, 4).mask)


def test_self_crossing_ring_is_illegal_but_rasterized():
    bowtie = np.array([[0, 0], [1000, 1000], [1000, 0], [0, 1000]])
    g = to_mask([(bowtie, [])], 8, 8)
    assert not g.legal and g.defects.get("self_crossing", 0) >= 1 and g.mask.any()


coord = st.integers(0, GRID)
ring = st.lists(st.tuples(coord, coord), min_size=3, max_size=12).map(lambda p: np.array(p, np.int64))
parts = st.lists(st.tuples(ring, st.lists(ring, max_size=2)), min_size=1, max_size=3)


@settings(max_examples=200, deadline=None)
@given(parts, st.integers(1, 40), st.integers(1, 40))
def test_to_mask_equals_the_oracles(ps, W, H):
    rings, roles = [], []
    for o, hs in ps:
        roles += [-1] + [len(rings)] * len(hs)
        rings += [o.tolist(), *(h.tolist() for h in hs)]
    expected = {d: c for d, c in zip(VO.DEFECTS, VO.ref_totals(rings, roles, GRID), strict=True) if c}
    g = to_mask(ps, W, H)
    assert np.array_equal(g.mask, RO.winding_np(rings, W, H, GRID)[1] % 2 == 1)
    assert g.defects == expected and g.legal == (not expected)
    # canonicalization changes neither the mask nor the defect counts of a legal shape
    if g.legal:
        cg = to_mask(canonicalize({"x": ps})["x"], W, H)
        assert cg.legal and np.array_equal(cg.mask, g.mask)
