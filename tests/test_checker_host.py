"""C2 (experiments/phase2/perf/implementation.md; issues CK-1 to CK-3): the checker returns the even-odd mask only,
bit-packed per row as uint32 words (LSB = leftmost pixel), and objdel.geometry unpacks it. The references here are
independent of the code under test: masks packed bit by bit in Python, and polycore's to_mask."""
import numpy as np
import pytest

from objdel import _gpucheck
from objdel.geometry import shape, to_mask, to_masks

WIDTHS = (1, 31, 32, 33, 63, 64, 65, 512)


def pack_reference(mask: np.ndarray) -> np.ndarray:
    """(H, W) bool -> (H * ceil(W / 32),) uint32, one bit at a time: pixel x of a row is bit x % 32 of word x // 32."""
    H, W = mask.shape
    nw = -(-W // 32)
    words = [0] * (H * nw)
    for y in range(H):
        for x in range(W):
            if mask[y, x]:
                words[y * nw + x // 32] |= 1 << (x % 32)
    return np.array(words, np.uint32)


@pytest.mark.parametrize("W", WIDTHS)
def test_unpack_matches_a_reference_mask(W):
    rng = np.random.default_rng(W)
    for H in (1, 3, 17):
        mask = rng.random((H, W)) < 0.5
        mask[0, 0], mask[-1, -1] = True, True  # the first and last pixel: bit 0 of the first word, the last bit
        words = pack_reference(mask)
        assert words.size == H * -(-W // 32)
        got = shape.unpack_rows(words, H, W)
        assert got.dtype == bool and got.shape == (H, W) and np.array_equal(got, mask)


def test_unpack_ignores_the_padding_bits():
    """Bits past W in a row's last word belong to no pixel."""
    assert not shape.unpack_rows(np.array([0xFFFFFFFE], np.uint32), 1, 1).any()
    assert shape.unpack_rows(np.array([0x80000001, 0xFFFFFFFF], np.uint32), 1, 33).tolist() == [
        [True] + [False] * 30 + [True, True]]


def test_to_masks_unpacks_the_checker_words(monkeypatch):
    """A fake checker returns known packed words for objects of every width: to_masks cuts and unpacks each one."""
    rng = np.random.default_rng(0)
    masks = [rng.random((H, W)) < 0.4 for W in WIDTHS for H in (1, 5)]
    words = [pack_reference(m) for m in masks]
    defects = np.zeros((len(masks), 13), np.int64)
    defects[3, 0] = 2

    def fake(xy, ring_start, shape_start, roles, W, H, G):
        assert W.tolist() == [m.shape[1] for m in masks] and H.tolist() == [m.shape[0] for m in masks]
        return {"word_start": np.cumsum([0] + [w.size for w in words]).astype(np.int64),
                "evenodd": np.concatenate(words), "defects": defects}

    monkeypatch.setattr(_gpucheck, "check_polygons", fake)
    out = to_masks([([], m.shape[1], m.shape[0]) for m in masks])
    assert [np.array_equal(g.mask, m) for g, m in zip(out, masks, strict=True)] == [True] * len(masks)
    assert [g.legal for g in out] == [i != 3 for i in range(len(masks))]


def _objects(seed: int, n: int):
    """Random legal and illegal objects (outers with holes, bow-ties, several parts) on crops of every WIDTH."""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        W, H = WIDTHS[i % len(WIDTHS)], int(rng.integers(1, 70))
        parts = []
        for _ in range(int(rng.integers(0, 3))):
            x0, y0 = rng.integers(0, 600, 2)
            x1, y1 = x0 + rng.integers(30, 400), y0 + rng.integers(30, 400)
            outer = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
            if rng.random() < 0.2:
                outer = outer[[0, 2, 1, 3]]  # a bow-tie: self-crossing
            holes = []
            if rng.random() < 0.4:
                cx, cy = (x0 + x1) // 2, (y0 + y1) // 2
                holes.append(np.array([[cx - 10, cy - 10], [cx - 10, cy + 10], [cx + 10, cy + 10], [cx + 10, cy - 10]]))
            parts.append((outer, holes))
        out.append((parts, W, H))
    return out


def test_to_masks_unchanged_on_fixtures():
    """The checker's masks and defects equal polycore's (to_mask) bit for bit, at every width around a word edge."""
    objs = _objects(1, 400)
    ref = [to_mask(p, W, H) for p, W, H in objs]
    got = to_masks(objs)
    bad = [i for i, (a, b) in enumerate(zip(ref, got, strict=True))
           if not (a.legal == b.legal and a.defects == b.defects and np.array_equal(a.mask, b.mask))]
    assert bad == []
    assert any(not g.legal for g in ref) and sum(g.mask.any() for g in ref) > 100
