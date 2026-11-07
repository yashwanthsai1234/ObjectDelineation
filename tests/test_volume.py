"""objdel.checker.volume: slice masks stacked into a voxel volume and its Betti numbers (b0 parts, 26-connected; b1
tunnels; b2 cavities, 6-connected background), the 3-D form of the per-slice topology. Expected values are worked by
hand for each shape."""

import numpy as np
import pytest

from objdel.checker import volume as V


def box(shape, fill=True):
    return np.full(shape, fill, bool)


def shell(n=3):
    v = box((n, n, n))
    v[1:-1, 1:-1, 1:-1] = False
    return v


def ring(n=3, depth=1):
    v = box((depth, n, n))
    v[:, 1:-1, 1:-1] = False
    return v


@pytest.mark.parametrize(("vol", "expected"), [
    (box((1, 1, 1)), (1, 0, 0)),             # one voxel
    (box((3, 3, 3)), (1, 0, 0)),             # solid cube
    (shell(3), (1, 0, 1)),                   # hollow cube: one cavity
    (shell(5), (1, 0, 1)),
    (ring(3), (1, 1, 0)),                    # a square ring: one tunnel
    (ring(5, depth=3), (1, 1, 0)),           # a thick ring
    (box((0, 0, 0)), (0, 0, 0)),             # empty
])  # fmt: skip
def test_betti_numbers_of_worked_shapes(vol, expected):
    assert V.betti3d(vol) == expected


def test_voxels_touching_at_a_corner_are_one_part():
    """26-connected parts (closed unit cubes): a shared vertex joins two voxels, and makes no tunnel."""
    v = np.zeros((2, 2, 2), bool)
    v[0, 0, 0] = v[1, 1, 1] = True
    assert V.betti3d(v) == (1, 0, 0)
    v2 = np.zeros((3, 3, 3), bool)
    v2[0, 0, 0] = v2[2, 2, 2] = True  # no shared vertex
    assert V.betti3d(v2) == (2, 0, 0)


def test_two_rings_and_a_shell_add_up():
    """Betti numbers add over separate parts."""
    v = np.zeros((5, 5, 15), bool)
    v[1, 1:4, 0:3] = ring(3)[0]
    v[1, 1:4, 5:8] = ring(3)[0]
    v[:5, :5, 10:15] = shell(5)
    assert V.betti3d(v) == (3, 2, 1)


def test_stack_places_slices_by_index_and_leaves_gaps_empty():
    a, b = box((2, 3)), np.eye(2, 3, dtype=bool)
    z0, vol = V.stack({7: a, 9: b})
    assert z0 == 7 and vol.shape == (3, 2, 3)
    assert vol[0].all() and not vol[1].any() and (vol[2] == b).all()
    with pytest.raises(ValueError, match="same shape"):
        V.stack({0: a, 1: box((3, 3))})


def test_stacked_slices_of_a_hollow_nodule_have_a_cavity():
    """Slice-wise outlines (a disc, rings, a disc) stack into a closed shell: per slice one hole, in 3-D one cavity."""
    disc, ring2d = box((3, 3)), ring(3)[0]
    _, vol = V.stack({0: disc, 1: ring2d, 2: disc})
    assert V.betti3d(vol) == (1, 0, 1)


def test_topology_is_scored_against_the_set_of_reader_volumes():
    """Set-valued, as per slice (topo_any): the prediction matches when any reader's volume has its Betti numbers."""
    readers = [box((3, 3, 3)), shell(3)]
    assert V.topology_any(shell(3), readers) and V.topology_any(box((2, 2, 2)), readers)
    assert not V.topology_any(ring(3), readers)
    assert not V.topology_any(box((1, 1, 1)), [])
