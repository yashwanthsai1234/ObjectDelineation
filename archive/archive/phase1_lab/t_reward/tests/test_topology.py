"""Task 1: exact per-instance topology (cv2 Spaghetti) vs brute force."""
import os

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from objdel.checker.bitquad import bitquads
from objdel.checker.topology import betti
from rewardcore import fixtures, synth
from rewardcore.reference import instance_topology_bfs, instance_topology_bruteforce
from rewardcore.topology import instance_topology, instance_topology_euler, instance_topology_full

N_RANDOM = int(os.environ.get("TOPO_N_RANDOM", 5000))


def check(mask):
    m = np.asarray(mask, bool)
    full = instance_topology_full(m)
    ref = instance_topology_bruteforce(m)
    assert full.as_list() == ref
    assert instance_topology_euler(m).as_list() == ref  # independent method agrees
    b0, b1 = betti(m)
    assert full.betti == (b0, b1)
    assert bitquads(m).chi8 == sum(1 - h for _, h in ref)  # certificate: chi8 = sum(1 - holes_k)
    assert full.areas.sum() == m.sum()


@pytest.mark.parametrize(
    "mask, expected",
    [
        (np.zeros((5, 5), bool), []),
        (np.ones((4, 6), bool), [(24, 0)]),
        (np.eye(4, dtype=bool), [(4, 0)]),  # diagonal chain is one 8-component
        (np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], bool), [(4, 1)]),  # diamond: 4-connected hole
        (np.pad(np.pad(np.pad(np.ones((1, 1), bool), 1), 1, constant_values=True), 1), [(16, 1), (1, 0)]),
        (synth.nested_rings(13, 13, 1), [(48, 1), (32, 1), (16, 1), (1, 0)]),
        (synth.checkerboard(4, 4), [(8, 2)]),  # two interior 4-isolated background pixels
    ],
)
def test_known_shapes(mask, expected):
    assert instance_topology(mask) == expected
    check(mask)


def test_hole_touching_island_is_not_counted_for_enclosing_component():
    m = np.zeros((9, 9), bool)
    m[1:8, 1:8] = True
    m[2:7, 2:7] = False  # one hole of the outer ring
    m[3, 3:6] = m[5, 3:6] = m[3:6, 3] = m[3:6, 5] = True  # island ring inside the hole, with its own hole
    assert instance_topology(m) == [(24, 1), (8, 1)]
    check(m)


def test_random_masks_vs_bruteforce():
    rng = np.random.default_rng(12345)
    kinds = {}
    for _ in range(N_RANDOM):
        name, m = synth.random_mask(rng, 8, 256)
        check(m)
        kinds[name] = kinds.get(name, 0) + 1
    assert sum(kinds.values()) == N_RANDOM and len(kinds) == 4


@pytest.mark.parametrize("dataset", ["riga", "lidc"])
def test_fixture_masks_vs_bruteforce(dataset):
    for stack in fixtures.load(dataset).values():
        for m in stack:
            check(m)
            check(~m)  # complement: object with a hole, frame touching the border


@settings(max_examples=300, deadline=None)
@given(arrays(bool, st.tuples(st.integers(1, 12), st.integers(1, 12))))
def test_small_masks_vs_pure_python_bfs(m):
    assert instance_topology(m) == instance_topology_bfs(m)
    check(m)
