"""The cache EDT on the engine (Phase D, D7; tests P08, U15, D-R06 and D-R07 of Implementation.md D.9, B15 of
D.spec.md).

`objdel.reward.cache.edt` calls `gc.edt`, the table build's integer EDT kernels on one map. The oracle is scipy's
exact transform cast to float32; equality is on bytes (so -0.0 and NaN would show). The stored caches stay valid only
if the fingerprint and every value `ExpertCache.build` computes are unchanged, so the data test rebuilds real caches
with both transforms.
"""
import os
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

import objdel.reward as R
from objdel.geom import engine as gc
from objdel.reward import ExpertCache
from objdel.reward import cache as rc

DATA = Path(os.environ.get("OBJDEL_DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
# the fingerprint before D7 (numpy 2.5.4, scipy 1.18.1)
FINGERPRINT = "671df12eab916e086c17cbe0dddab5c0c17d0a3224eae2f0b060a213b3cb54a3"


def scipy_edt(nz) -> np.ndarray:
    return ndimage.distance_transform_edt(nz).astype(np.float32)


def _same(nz):
    nz = np.ascontiguousarray(nz, np.uint8)
    got = gc.edt(nz)
    assert got.dtype == np.float32 and got.shape == nz.shape
    want = np.full(nz.shape, rc._BIG, np.float32) if nz.all() else scipy_edt(nz)  # scipy has no answer without a zero
    assert got.tobytes() == want.tobytes(), nz.shape


def test_edt_bits_equal_scipy_on_random_and_thin_maps():
    """P08: 300 random maps (sides 1..300, density 0.001..0.999), single zero and single nonzero pixels, zeros at the
    corners, thin lines, 1xN, Nx1 and 4096x2."""
    rng = np.random.default_rng(0)
    for _ in range(300):
        nz = rng.random((int(rng.integers(1, 301)), int(rng.integers(1, 301)))) < rng.uniform(0.001, 0.999)
        nz.flat[int(rng.integers(nz.size))] = False  # at least one zero pixel (no-zero maps: test below)
        _same(nz)
    for h, w in ((1, 1), (5, 7), (64, 33), (1, 300), (300, 1), (4096, 2), (2, 4096)):
        for y, x in ((0, 0), (0, w - 1), (h - 1, 0), (h - 1, w - 1), (h // 2, w // 2)):
            one_zero = np.ones((h, w), np.uint8)
            one_zero[y, x] = 0
            _same(one_zero)
            one_px = np.zeros((h, w), np.uint8)
            one_px[y, x] = 1
            _same(one_px)
    lines = np.ones((40, 50), np.uint8)
    lines[:, 7] = lines[13, :] = 0
    _same(lines)
    _same(lines.T)


def test_edt_special_cases():
    """U15: no zero pixel -> _BIG everywhere (1e6, pinned from the pre-D source); all-zero -> zeros; 0-size -> an
    empty (H, W) map; equidistant ties; one zero at a corner of 4096 x 4096 (squared distances above 2^24)."""
    assert rc._BIG == np.float32(1e6) and rc._BIG.dtype == np.float32
    for shape in ((7, 7), (1, 1), (7, 5)):
        got = rc.edt(np.ones(shape, bool))
        assert got.dtype == np.float32 and got.shape == shape and (got == rc._BIG).all()
    assert rc.edt(np.zeros((6, 9), bool)).tobytes() == np.zeros((6, 9), np.float32).tobytes()
    for shape in ((0, 5), (5, 0), (0, 0)):
        assert rc.edt(np.zeros(shape, bool)).shape == shape
    ties = np.ones((21, 21), np.uint8)
    for dy, dx in ((3, 4), (-3, 4), (3, -4), (-3, -4), (4, 3), (-4, -3), (5, 0), (-5, 0), (0, 5), (0, -5)):
        ties[10 + dy, 10 + dx] = 0
    _same(ties)
    assert rc.edt(ties)[10, 10] == np.float32(5.0)
    big = np.ones((4096, 4096), np.uint8)
    big[0, 0] = 0
    got = rc.edt(big)
    assert got.tobytes() == scipy_edt(big).tobytes()
    assert got[-1, -1] == np.float32(np.sqrt(2.0 * 4095 ** 2))


def test_edt_input_limits():
    """Sides up to 32767 (the cache build pads a cropped 4096 frame to 4098: P12's all-ones 4096 x 4096 mask)."""
    padded = np.ones((4098, 3), np.uint8)
    padded[0, 0] = 0
    _same(padded)
    for shape in ((32768, 1), (1, 32768)):
        with pytest.raises(ValueError):
            gc.edt(np.ones(shape, np.uint8))
    with pytest.raises(TypeError):  # 2-D uint8 only (cache.edt converts)
        gc.edt(np.ones((3, 3, 3), np.uint8))


def test_fingerprint_is_unchanged():
    """D-R06: the fingerprint the stored caches carry; D7 replaces scipy's EDT by the engine's without changing it."""
    assert R.fingerprint() == FINGERPRINT


def _caches(ds: str, limit: int) -> list[Path]:
    files = sorted((DATA / "prepared" / "caches" / ds).rglob("*.npz"))[:limit]
    if not files:
        pytest.skip(f"no {ds} caches under {DATA}")
    return files


@pytest.mark.data
@pytest.mark.slow
@pytest.mark.parametrize("ds, limit", [("lidc", 1000), ("riga", 10 ** 6)])
def test_cache_build_is_unchanged_on_real_caches(ds, limit, monkeypatch):
    """D-R07 (D-G10): the first 1,000 LIDC caches by sorted path (and every RIGA cache, Mac only) rebuilt with
    ExpertCache.build: tau, T5, dist_out_union and dist_in_core bytes equal those of the scipy EDT."""
    for f in _caches(ds, limit):
        with np.load(f, allow_pickle=False) as z:
            raters = z["raters"].astype(bool)
        new = ExpertCache.build(raters)
        maps = (new.dist_out_union.tobytes(), new.dist_in_core.tobytes())
        with monkeypatch.context() as m:
            m.setattr(rc, "edt", lambda nz: np.full(np.shape(nz), rc._BIG, np.float32) if np.all(nz)
                      else scipy_edt(np.ascontiguousarray(nz, np.uint8)))
            old = ExpertCache.build(raters)
            assert (new.tau.hex(), new.t5.hex()) == (old.tau.hex(), old.t5.hex()), f
            assert maps == (old.dist_out_union.tobytes(), old.dist_in_core.tobytes()), f
