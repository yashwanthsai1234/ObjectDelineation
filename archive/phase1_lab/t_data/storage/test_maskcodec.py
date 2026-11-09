import numpy as np
import pytest
from scipy import ndimage

from maskcodec import poly_decode, poly_encode, rle_decode, rle_encode


def random_masks(n=300, seed=0):
    rng = np.random.default_rng(seed)
    for i in range(n):
        h, w = rng.integers(1, 80, 2)
        m = rng.random((h, w)) < rng.uniform(0.05, 0.9)
        if i % 3 == 0:
            m = ndimage.binary_opening(m) if m.any() else m
        yield m
    yield np.zeros((5, 7), bool)
    yield np.ones((5, 7), bool)
    yield np.eye(9, dtype=bool)  # diagonal-only links (saddles)
    yield np.eye(9, dtype=bool)[::-1] | np.eye(9, dtype=bool)


@pytest.mark.parametrize("enc,dec", [(poly_encode, poly_decode), (rle_encode, rle_decode)])
def test_roundtrip_bit_exact(enc, dec):
    for m in random_masks():
        assert np.array_equal(dec(enc(m)), m)
