"""tests/oracles/crops_oracle.py: the phase 1 lab's pad_div crop, ported so no test runs lab code (CLEAN4)."""

import numpy as np

from tests.oracles.crops_oracle import mask_sha, pad_div


def test_small_mask_pads_to_the_next_divisor_of_1000():
    m = np.ones((3, 3), bool)  # divisors of 1000: 1 2 4 5 8 ...; the next at or above 3 is 4
    out = pad_div({"a": m})["a"]
    assert out.shape == (4, 4) and out.dtype == bool and out[:3, :3].all() and out.sum() == 9


def test_non_square_mask_is_centred_with_floor_offsets():
    m = np.ones((3, 6), bool)  # side 8, offsets (8-3)//2 = 2 and (8-6)//2 = 1
    out = pad_div({"a": m})["a"]
    assert out.shape == (8, 8) and out[2:5, 1:7].all() and out.sum() == 18


def test_empty_mask_stays_empty():
    out = pad_div({"a": np.zeros((5, 5), bool)})["a"]
    assert out.shape == (5, 5) and not out.any()


def test_large_mask_is_block_reduced_first():
    m = np.ones((600, 600), bool)  # factor 2 -> 300 x 300, then side 500 and offset 100
    out = pad_div({"a": m})["a"]
    assert out.shape == (500, 500) and out.sum() == 300 * 300 and out[100:400, 100:400].all()


def test_block_majority_ties_go_to_foreground():
    m = np.zeros((600, 600), bool)
    m[0, 0:2] = True  # a 2 x 2 block with 2 of 4 pixels set: a tie, kept
    m[2, 2] = True  # a block with 1 of 4: dropped
    out = pad_div({"a": m})["a"]
    assert out.sum() == 1 and out[100, 100]


def test_mask_sha_depends_on_shape_and_content():
    a, b = np.zeros((4, 4), bool), np.zeros((2, 8), bool)
    assert mask_sha({"x": a}) == mask_sha({"x": a.copy()}) != mask_sha({"x": b})
    c = a.copy()
    c[1, 1] = True
    assert mask_sha({"x": a}) != mask_sha({"x": c}) and len(mask_sha({"x": a})) == 64
