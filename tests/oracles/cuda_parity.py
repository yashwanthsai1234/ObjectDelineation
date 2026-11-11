"""Pack cases for objdel._gpucheck, run it, and compare it with the references: for polygons, the results of the CPU
checker it replaced (polycore, archived in archive/cpu_checker), recorded per case in tests/goldens/polycore.npz
before the archive; for masks, the numpy polygonizer (tests/oracles/polygonize_oracle.py, never the engine under
test), the bit-quad oracle (tests/oracles/bitquad_oracle.py) and cv2's 8-connected labels (ring parts), run here.

Shared by tests/geom/test_cuda_checker*.py, tests/geom/test_checker_host.py, tests/verify/test_gpu_scoring.py and
tools/checker/eval.py. The packed formats are those of data/checker_eval
(archive/cpu_checker/experiments/prepare_data.py): polygons as (xy, ring_start, shape_start, roles, W, H), masks as
(bits, hw, start). The case generators are the repo's existing ones (tests/oracles, tests/geom/test_topology.py),
imported, not copied; the rollout, int32-extreme and at-the-caps families are defined here.
"""
import binascii
import functools
import hashlib
import random
import sys
from pathlib import Path

import cv2
import numpy as np

from objdel.geom import engine as _gpucheck
from objdel.geom.shape import pack_rings, unpack_rows
from tests.oracles.bitquad_oracle import bitquads
from tests.oracles.polygonize_oracle import polygonize

TESTS = Path(__file__).resolve().parents[2] / "tests"
if str(TESTS) not in sys.path:
    sys.path.append(str(TESTS))  # tests/oracles and the generator functions of the test modules

GOLDEN = TESTS / "goldens" / "polycore.npz"  # archive/cpu_checker/record_goldens.py
POLY_FIELDS = ("legal", "defects", "evenodd")  # compared in this order; the first difference is reported
MASK_FIELDS = ("quads", "rings", "b0", "b1", "betti", "parts", "ring_sha")


# ------------------------------------------------------------------ packing and running
def pack_polygons(cases) -> dict:
    """[(rings, roles or None, W, H)] -> the packed arrays of polygons.npz (roles default to all outer)."""
    xy, ring_start, shape_start, roles, W, H = [], [0], [0], [], [], []
    for rings, rl, w, h in cases:
        pxy, off = pack_rings(rings)
        xy.append(pxy)
        ring_start += (ring_start[-1] + off[1:]).tolist()
        shape_start.append(shape_start[-1] + len(rings))
        roles += [-1] * len(rings) if rl is None else list(rl)
        W.append(w), H.append(h)
    return dict(xy=np.concatenate(xy) if xy else np.zeros((0, 2), np.int32),
                ring_start=np.array(ring_start, np.int64), shape_start=np.array(shape_start, np.int64),
                roles=np.array(roles, np.int32), W=np.array(W, np.int32), H=np.array(H, np.int32))


def run_polygons(p: dict, G: int = 1000) -> dict:
    return _gpucheck.check_polygons(p["xy"], p["ring_start"], p["shape_start"], p["roles"], p["W"], p["H"], G)


def shape_of(p: dict, s: int):
    """Shape s of a packed set -> (rings as int32 arrays, roles list, W, H)."""
    r0, r1 = int(p["shape_start"][s]), int(p["shape_start"][s + 1])
    rs = p["ring_start"]
    return ([p["xy"][rs[r]:rs[r + 1]] for r in range(r0, r1)], p["roles"][r0:r1].tolist(), int(p["W"][s]),
            int(p["H"][s]))


def gpu_polygon(out: dict, p: dict, s: int) -> dict:
    """Shape s of a check_polygons result: its defects and its even-odd mask, unpacked from the packed rows."""
    a, b = out["word_start"][s], out["word_start"][s + 1]
    return dict(defects=out["defects"][s], evenodd=unpack_rows(out["evenodd"][a:b], int(p["H"][s]), int(p["W"][s])))


def pack_masks(masks) -> dict:
    ms = [np.ascontiguousarray(m != 0, np.uint8) for m in masks]
    start = np.zeros(len(ms) + 1, np.int64)
    start[1:] = np.cumsum([m.size for m in ms])
    return dict(bits=np.concatenate([m.ravel() for m in ms]) if ms else np.zeros(0, np.uint8),
                hw=np.array([m.shape for m in ms], np.int32).reshape(-1, 2), start=start)


def run_masks(p: dict) -> dict:
    return _gpucheck.check_masks(p["bits"], p["hw"], p["start"])


def mask_of(p: dict, m: int) -> np.ndarray:
    return p["bits"][p["start"][m]:p["start"][m + 1]].reshape(p["hw"][m]).astype(bool)


def gpu_rings(out: dict, m: int) -> list:
    """Rings of mask m as [(pts int64 (k, 2), signed area)], in polygonize's order."""
    rings = []
    for r in range(out["ring_start"][m], out["ring_start"][m + 1]):
        pts = out["pts"][out["pts_start"][r]:out["pts_start"][r + 1]].astype(np.int64)
        rings.append((pts, int(out["ring_area"][r])))
    return rings


# ------------------------------------------------------------------ references and comparison
def ring_sha(rings) -> bytes:
    """prepare_data.py's ring_sha: SHA-256 of every ring's vertices (int64) and signed area (int64 LE), in order."""
    return hashlib.sha256(b"".join(np.ascontiguousarray(pts, np.int64).tobytes()
                                   + int(area).to_bytes(8, "little", signed=True) for pts, area in rings)
                          ).hexdigest().encode()


def mask_sha(mask) -> bytes:
    """prepare_data.py's evenodd_sha."""
    return hashlib.sha256(np.packbits(np.asarray(mask, bool)).tobytes()).hexdigest().encode()


def object_case(parts) -> tuple[list, list]:
    """An object's parts [(outer, holes)] -> (rings, roles), as objdel.geom.shape.pack_objects packs them."""
    rings, roles = [], []
    for outer, holes in parts:
        roles += [-1] + [len(rings)] * len(holes)
        rings += [outer, *holes]
    return rings, roles


def case_key(rings, roles, W, H, G=1000) -> bytes:
    """Key of one polygon case in GOLDEN: the first 16 bytes of the SHA-256 of its grid and packed arrays."""
    p = pack_polygons([(rings, roles, W, H)])
    h = hashlib.sha256(int(G).to_bytes(8, "little"))
    for k in ("xy", "ring_start", "roles", "W", "H"):
        h.update(np.ascontiguousarray(p[k]).tobytes())
    return h.digest()[:16]


@functools.cache
def goldens() -> dict:
    z = np.load(GOLDEN, allow_pickle=False)
    return {k.tobytes(): (d, binascii.hexlify(e.tobytes())) for k, d, e in zip(z["key"], z["defects"], z["evenodd"],
                                                                              strict=True)}


def polygon_reference(rings, roles, W, H, G=1000) -> dict:
    """polycore's result for one case: its 13 defect counts and its even-odd mask's SHA-256 (mask_sha), recorded
    before the CPU checker was archived. roles: one per ring (-1 = outer). A case never recorded (the K1 "holes"
    generator rounds math.cos/sin, whose last bit differs between platforms: R3 met 4 such cases of ~7,600 on Linux)
    takes its reference from the exact integer oracles instead (tests/oracles: crossing parity, ref_totals), which
    the recorded results equal (tests/geom/test_checker_oracle.py)."""
    if (got := goldens().get(case_key(rings, roles, W, H, G))) is not None:
        defects, evenodd = got
        return dict(defects=defects, evenodd=evenodd)
    from oracles import raster_oracle as RO
    from oracles import validity_oracle as VO

    rings = [[(int(x), int(y)) for x, y in r] for r in rings]  # Python ints: the oracles' exact integer arithmetic
    roles = [int(r) for r in roles]
    return dict(defects=np.array(VO.ref_totals(rings, roles, G), np.int64),
                evenodd=mask_sha(RO.winding_and_parity(rings, W, H, G)[1] % 2 == 1))


def geometry_diff(g, ref: dict):
    """polygon_diff of an objdel.geom.shape.Geometry (to_mask / to_masks)."""
    return polygon_diff({"defects": np.array([g.defects.get(d, 0) for d in _gpucheck.DEFECTS]), "evenodd": g.mask}, ref)


def polygon_diff(gpu: dict, ref: dict):
    """First field of POLY_FIELDS where the GPU result differs from the reference (None if all agree). Masks may be
    arrays or SHA-256 digests (mask_sha)."""
    same = {
        "legal": (gpu["defects"].sum() == 0) == (np.asarray(ref["defects"]).sum() == 0),
        "defects": np.array_equal(gpu["defects"], ref["defects"]),
    }
    r = ref["evenodd"]
    same["evenodd"] = mask_sha(gpu["evenodd"]) == r if isinstance(r, bytes) else np.array_equal(gpu["evenodd"], r)
    return next((k for k in POLY_FIELDS if not same[k]), None)


def ring_parts(mask, rings) -> tuple:
    """Per ring, the raster rank (by first pixel) of the cv2 8-connected component of the foreground pixel right of the
    ring's first edge: an outer ring's own component, a hole's surrounding one (the engine's ring_part). cv2 numbers
    components in its own scan order (2x2 blocks), so its labels are ranked by first pixel here."""
    if not rings:
        return ()
    lab = cv2.connectedComponents(np.asarray(mask, np.uint8), connectivity=8)[1]
    labels, first = np.unique(lab.ravel(), return_index=True)
    fg = labels > 0
    rank = np.full(labels.max() + 1, -1, np.int64)
    rank[labels[fg][np.argsort(first[fg])]] = np.arange(fg.sum())
    out = []
    for pts, _ in rings:
        (x, y), (dx, dy) = pts[0], np.sign(pts[1] - pts[0])  # east: pixel (x, y); south (x-1, y); west, north
        out.append(int(rank[lab[y - (dx < 0 or dy < 0), x - (dx < 0 or dy > 0)]]))
    return tuple(out)


def mask_reference(mask) -> dict:
    q = bitquads(mask)
    rings = [(g.pts.astype(np.int64), g.area) for g in polygonize(mask)]
    return dict(quads=(q.n1, q.n2, q.nD, q.n3, q.n4), n_rings=len(rings), b0=sum(a > 0 for _, a in rings),
                b1=sum(a < 0 for _, a in rings), parts=ring_parts(mask, rings), ring_sha=ring_sha(rings))


def mask_diff(out: dict, m: int, ref: dict):
    """First field of MASK_FIELDS where the GPU result for mask m differs from ref (mask_reference, or the
    certify columns of masks.npz: b0, b1, area, corners and ring_sha). `betti` is the engine's own (b0, b1) count,
    `parts` its ring_part."""
    rings = gpu_rings(out, m)
    r0, r1 = out["ring_start"][m], out["ring_start"][m + 1]
    q = out["quads"][m]
    n1, n2, nd, n3, n4 = (int(v) for v in q)
    if "quads" in ref:
        quads_ok = tuple(int(v) for v in q) == tuple(ref["quads"])
    else:  # certify's columns: area, perimeter, corners, chi8 are functions of the quads
        quads_ok = ((n1 + 2 * n2 + 3 * n3 + 4 * n4 + 2 * nd) // 4, n1 + n2 + n3 + 2 * nd, n1 + n3 + 2 * nd,
                    (n1 - n3 - 2 * nd) // 4) == (ref["area"], ref["perimeter"], ref["corners"], ref["chi8"])
    same = {"quads": quads_ok,
            "rings": "n_rings" not in ref or len(rings) == ref["n_rings"],
            "b0": sum(a > 0 for _, a in rings) == ref["b0"],
            "b1": sum(a < 0 for _, a in rings) == ref["b1"],
            "betti": tuple(int(v) for v in out["betti"][m]) == (ref["b0"], ref["b1"]),
            "parts": "parts" not in ref or tuple(int(v) for v in out["ring_part"][r0:r1]) == ref["parts"],
            "ring_sha": ring_sha(rings) == ref["ring_sha"]}
    return next((k for k in MASK_FIELDS if not same[k]), None)


# ------------------------------------------------------------------ case generators
def huge_shape(rng: random.Random):
    """int32-extreme coordinates and grids (the __int128 paths) -> (rings, W, H, G)."""
    i32 = np.iinfo(np.int32)
    pick = [i32.min, i32.min + 1, -10**9, -1, 0, 1, 37, 500, 999, 1000, 10**9, i32.max]
    rings = [[(rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000),
               rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000))
              for _ in range(rng.randint(1, 7))] for _ in range(rng.randint(1, 2))]
    return rings, rng.randint(1, 9), rng.randint(1, 9), rng.choice([1000, 7, 2**30])


def caps_shape(rng: random.Random, size: int):
    """At the input caps: 10 rings of 1000 vertices (one outer and nine holes, or ten parts) on a size^2 crop."""
    t = np.linspace(0, 2 * np.pi, 1000, endpoint=False)
    def circle(cx, cy, r, ccw):
        pts = np.c_[cx + r * np.cos(t), cy + r * np.sin(t)].round().astype(np.int64)
        return pts[::-1] if ccw else pts
    if rng.random() < .5:
        rings = [circle(500, 500, 450, False)] + [circle(200 + 150 * (k % 3), 200 + 150 * (k // 3), 50, True)
                                                    for k in range(9)]
        return rings, [-1] + [0] * 9, size, size
    return [circle(100 + 80 * k, 500, 30, False) for k in range(10)], None, size, size


def blob(rng, n, cx, cy, r, ccw=False):
    """A smooth random ring of n vertices around (cx, cy), radius about r, on the 1000 grid (repeats dropped)."""
    k = np.arange(1, 5)
    a, ph = rng.normal(0, 0.25 / k), rng.uniform(0, 2 * np.pi, k.size)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rad = np.clip(r * (1 + (a[:, None] * np.cos(k[:, None] * t + ph[:, None])).sum(0)), 0.3 * r, None)
    pts = np.clip(np.rint(np.c_[cx + rad * np.cos(t), cy + rad * np.sin(t)]), 0, 1000).astype(np.int64)
    pts = pts[np.r_[True, np.any(np.diff(pts, axis=0) != 0, axis=1)]]
    return pts[::-1] if ccw else pts


def rollouts(n: int, seed: int):
    """GRPO-like rollouts -> [(rings, roles, S)]: 1-3 smooth rings (holes or second parts), 16-200 vertices, crops
    64-512 px (the generator of the archived test_polycore.py in archive/cpu_checker; its ground-truth blob is still
    drawn, so the cases are the same)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        S = int(rng.integers(64, 513))
        cx, cy = rng.integers(350, 651, 2)
        R = float(rng.uniform(150, 330))
        rings, rl = [blob(rng, int(rng.integers(16, 201)), cx, cy, R)], [-1]
        for _ in range(int(rng.integers(0, 3))):
            if rng.random() < 0.5:
                rings.append(blob(rng, int(rng.integers(16, 201)), cx, cy, R * rng.uniform(.15, .3), ccw=True))
                rl.append(0)
            else:
                rings.append(blob(rng, int(rng.integers(16, 201)), *rng.integers(100, 901, 2), R * .2))
                rl.append(-1)
        blob(rng, 64, cx + 9, cy - 7, R)  # the archived generator's ground truth
        out.append((rings, rl, S))
    return out


def fuzz_polygons(n_per_family: int, seed: int = 0):
    """Seeded fuzz cases -> [(family, rings, roles, W, H, G)], from the repo's generators: the K1 raster families
    (star, soup, grid, halfpixel, holes, multipart, tiny, degenerate, rect), the K2 validity rings, multi_ring
    (outers with holes and roles, mutated), GRPO-like rollouts, int32-extreme coordinates, and shapes at the caps."""
    from oracles import raster_oracle as RO
    from oracles import validity_oracle as VO
    rng = random.Random(seed)
    out = []
    for i in range(n_per_family * len(RO.K1_FAMILIES)):
        fam, rings, W, H = RO.k1_case(i, seed=12345 + seed)
        out.append((fam, rings, None, W, H, 1000))
    for i in range(n_per_family * len(VO.K2_FAMILIES)):
        fam, ring = VO.k2_ring(i, seed=777 + seed)
        out.append(("ring_" + fam, [ring], None, rng.choice([7, 16, 64]), rng.choice([7, 16, 64]), 1000))
    for _ in range(n_per_family):
        rings, roles = VO.multi_ring(rng)
        out.append(("roles", rings, roles, rng.choice([16, 64, 125]), rng.choice([16, 64, 125]), 1000))
    out += [("rollout", rings, roles, S, S, 1000) for rings, roles, S in rollouts(n_per_family, seed=7 + seed)]
    for _ in range(n_per_family):
        rings, W, H, G = huge_shape(rng)
        out.append(("huge", rings, None, W, H, G))
    for k in range(max(1, n_per_family // 50)):
        rings, roles, W, H = caps_shape(rng, (64, 512, 4096)[k % 3])
        out.append(("caps", rings, roles, W, H, 1000))
    return out


def fuzz_masks(n: int, seed: int = 0):
    """[(family, mask)]: tests/geom/test_topology.random_masks (ten families, 1-200 px) plus fixed edge cases."""
    from tests.geom.test_topology import random_masks
    edge = {"empty_0x0": np.zeros((0, 0), bool), "empty_0x5": np.zeros((0, 5), bool), "pixel": np.ones((1, 1), bool),
            "full": np.ones((9, 7), bool), "checker": (np.indices((12, 13)).sum(0) % 2) == 0,
            "saddle_6": np.array([[0, 1], [1, 0]], bool), "saddle_9": np.array([[1, 0], [0, 1]], bool),
            "frame": np.pad(np.zeros((3, 3), bool), 1, constant_values=True)}
    return list(edge.items()) + [("random", m) for m in random_masks(n, seed)]
