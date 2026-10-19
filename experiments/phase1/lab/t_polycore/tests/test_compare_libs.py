"""Cross-checks: the numba prototype (exact agreement), skimage / rasterio (even-odd, away
from ties), GEOS containment in general position."""
import random

import numpy as np
import pytest

import oracle as O
import polycore as pc
from helpers import k1_case, k2_ring, on_boundary_centres


def _clean(ring):
    pts = [p for i, p in enumerate(ring) if i == 0 or p != ring[i - 1]]
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def test_numba_prototype_raster_agrees():
    proto = pytest.importorskip("numba_proto")
    n = 0
    for i in range(3000):
        _, rings, W, _ = k1_case(i, seed=5)
        S = min(W, 512)
        ours = pc.rasterize(rings, S, S, winding=True)["winding"]
        theirs = np.zeros((S, S), np.int32)
        for r in rings:
            if len(r):
                theirs += proto.scan(np.asarray(r, np.int64) * 2 * S, S)
        assert (ours == theirs).all(), i
        n += 1
    assert n == 3000


def test_numba_prototype_validity_agrees():
    proto = pytest.importorskip("numba_proto")
    for i in range(10000):
        _, ring = k2_ring(i, seed=8)
        pts = _clean(ring)
        if len(set(pts)) < 3:
            continue
        code = proto.seg_quadratic(np.asarray(pts, np.int64))
        row = pc.validate([pts])["ring"][0]
        D = {k: j for j, k in enumerate(pc.DEFECTS)}
        mine = 2 if row[D["self_crossing"]] else (1 if any(row[D[k]] for k in pc.SELF_TOUCHING) else 0)
        assert code == mine, (i, pts)


def _pixel_poly(ring, W, H):
    return [(x * W / 1000.0, y * H / 1000.0) for x, y in ring]


def test_skimage_rasterio_agree_off_boundary():
    skd = pytest.importorskip("skimage.draw")
    rio = pytest.importorskip("rasterio.features")
    from rasterio.transform import Affine
    import shapely
    off_sk = off_rio = on_px = tot = 0
    rng = random.Random(0)
    for i in range(600):
        _, rings, W, H = k1_case(rng.randrange(10**6), seed=17)
        rings = [r for r in rings if len(set(r)) >= 3][:1]
        if not rings:
            continue
        W, H = min(W, 333), min(H, 333)
        ours = pc.rasterize(rings, W, H)["evenodd"]
        B = on_boundary_centres(rings, W, H)
        p = np.array(_pixel_poly(rings[0], W, H))
        sk = np.zeros((H, W), bool)
        rr, cc = skd.polygon(p[:, 1] - 0.5, p[:, 0] - 0.5, (H, W))
        sk[rr, cc] = True
        g = shapely.Polygon(p)
        ri = rio.rasterize([(g, 1)], out_shape=(H, W), transform=Affine.identity(), fill=0, dtype="uint8").astype(bool)
        off_sk += int(((sk != ours) & ~B).sum())
        off_rio += int(((ri != ours) & ~B).sum())
        on_px += int(B.sum())
        tot += 1
    print(f"shapes={tot} on-edge centres={on_px} off-boundary mismatches skimage={off_sk} rasterio={off_rio}")
    assert off_sk == 0
    assert off_rio == 0


def test_containment_vs_geos_general_position():
    shapely = pytest.importorskip("shapely")
    rng = random.Random(2)
    checked = 0
    for _ in range(2000):
        a = O.random_simple_polygon(rng, rng.randint(3, 12))
        b = O.random_simple_polygon(rng, rng.randint(3, 10))[::-1]
        res = pc.validate([a, b], [-1, 0])
        if res["pair_cross"][0, 1] or res["pair_touch"][0, 1]:
            continue
        inside = shapely.Polygon(a).contains(shapely.Polygon(b))
        assert bool(res["ring"][1][pc.DEFECTS.index("hole_outside_outer")]) == (not inside)
        checked += 1
    assert checked > 1000
