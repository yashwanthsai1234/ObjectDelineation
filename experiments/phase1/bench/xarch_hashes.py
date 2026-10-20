"""Cross-architecture determinism check (Phase 1, M9): SHA-256 of every hashed output of the library on a fixed
synthetic corpus. Run it on the Mac (arm64) and in the linux/amd64 container (x86.sh); the JSON outputs must be
equal section by section (`--compare a.json b.json`).

    .venv/bin/python -I -m experiments.phase1.bench.xarch_hashes --out build/x86/hashes-mac.json
    .venv/bin/python -I -m experiments.phase1.bench.xarch_hashes --compare build/x86/hashes-mac.json build/x86/hashes-x86.json

No data/ is read: the corpus is generated from integer draws of numpy's PCG64 only (integer generation is
platform-independent; no libm, no float sampling), so it is bit-identical on every machine, and the `corpus` section
checks that. RIGA-like items (6 raters, disc and cup, RGB) and LIDC-like items (2-4 raters, small nodules with holes,
specks, pinches and abstaining raters, int16 Hounsfield units) go through the production path:

    frame (crop; LIDC: lung window + bicubic upscale)  -> frame.masks, frame.image, png
    polycore rasterize / validate / score_batch (random rings, threads 1 and 3) -> polycore.*
    checker.topology on every frame mask and polygonize   -> topology, polygonize
    targets.build_target (reward size rule, and all rules) -> targets, target_set_hash
    reward.ExpertCache.build (tau, T5, digest)             -> caches, fingerprint
    reward.score_batch on 2 groups of 8 completions per item, advantages_v4, topology_reward -> reward.*
    eval.metrics on 16 completions per item                -> eval
    data.splits: lidc_split, riga_split, riga.descriptor + near_duplicates, riga.keypoints + keypoint_duplicates

Floats are hashed exactly (float.hex / raw bytes). Every section stores one short digest per record, so a mismatch
names the first differing record. `fingerprint` and the cache digests include the numpy / scipy / OpenCV version
strings (reward.spec.fingerprint): they match only when the versions match, which x86.sh ensures by installing the
Mac's exact pins.

References:
    numpy random (PCG64; "the same seed gives the same stream on every platform" for integer draws, with the
        Generator stream-compatibility policy): https://numpy.org/doc/stable/reference/random/compatibility.html
    IEEE 754 exactness of float.hex: https://docs.python.org/3/library/stdtypes.html#float.hex
    Determinism of floating-point results across ISAs (summation order, FMA contraction, SIMD dispatch):
        https://numpy.org/doc/stable/reference/simd/index.html ;
        https://docs.opencv.org/4.x/db/de0/group__core__utils.html (setUseOptimized: optimised SIMD code paths)
"""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import cv2
import numpy as np

from objdel.checker import betti, instance_topology, min_area, polycore, polygonize, signature, size_filter
from objdel.data import Item, splits
from objdel.data import riga as riga_mod
from objdel.eval import metrics
from objdel.geometry import SCHEMAS, parse, serialize
from objdel.reward import ExpertCache, advantages_v4, fingerprint, score_batch, topology_reward
from objdel.targets import BuildConfig, build_target, frame, target_set_hash

SEED = 20261019
# Accepted across architectures (coordinator, 2026-10-19: F4 skipped): raw split scores, whose decisions (ncc_pairs,
# keypoint_dups, the split assignment) are compared exactly; splits are computed once and pinned by the store
# manifest. --strict (x86 run-to-run) accepts nothing. Every other section has been identical on arm64 and x86-64
# since F1-F3 (frame.image / png: Pillow bicubic; caches, reward.*, eval: scipy EDT and frozen S_m), so a difference
# there is a regression; the diagnosis history is in experiments/phase1/impl/M9_x86.md.
CAUSES = {  # why each accepted section differs (M9_x86.md finding 4)
    "splits.descriptor": "np.linalg.norm: float32 BLAS sdot kernels differ per arch [F4, skipped]",
    "splits.ncc_values": "FFT rounding; the decisions are splits.ncc_pairs",
    "splits.keypoints": "ORB float paths (CLAHE, resize, Harris); check splits.keypoint_pairs decisions >= ORB_DUP",
    "splits.keypoint_pairs": "inlier counts (threshold 0) follow keypoints; duplicate decisions >= ORB_DUP matched",
}
ACCEPTED = tuple(CAUSES)


def _h(*parts) -> str:
    h = hashlib.sha256()
    for p in parts:
        if isinstance(p, np.ndarray):
            h.update(repr((p.dtype.str, p.shape)).encode())
            h.update(np.ascontiguousarray(p).tobytes())
        elif isinstance(p, float):
            h.update(p.hex().encode())
        else:
            h.update(repr(p).encode())
        h.update(b"\x00")
    return h.hexdigest()


def _exact(x):
    """JSON-like value with every float as its hex string (exact)."""
    if isinstance(x, float):
        return x.hex()
    if isinstance(x, dict):
        return {str(k): _exact(v) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    if isinstance(x, list | tuple):
        return [_exact(v) for v in x]
    if isinstance(x, np.ndarray):
        return _exact(x.tolist())
    if isinstance(x, np.generic):
        return _exact(x.item())
    return x


# ---- corpus (integer-only generation) ------------------------------------------------------------------------------

def _ellipse(H, W, cy, cx, ry, rx):
    y, x = np.ogrid[:H, :W]
    return (x - cx) ** 2 * ry * ry + (y - cy) ** 2 * rx * rx <= rx * rx * ry * ry


def _blob(rng, H, W, cy, cx, r):
    m = np.zeros((H, W), bool)
    for _ in range(int(rng.integers(1, 4))):
        dy, dx = (int(v) for v in rng.integers(-r // 3, r // 3 + 1, 2))
        m |= _ellipse(H, W, cy + dy, cx + dx, max(1, r + int(rng.integers(-r // 4, r // 4 + 1))),
                      max(1, r + int(rng.integers(-r // 4, r // 4 + 1))))
    return m


def _jitter(rng, m):
    """One rater's version: integer dilation / erosion and random edits along the boundary."""
    m8 = m.astype(np.uint8)
    op = int(rng.integers(0, 3))
    if op == 1:
        m8 = cv2.dilate(m8, np.ones((3, 3), np.uint8), iterations=int(rng.integers(1, 3)))
    elif op == 2:
        m8 = cv2.erode(m8, np.ones((3, 3), np.uint8), iterations=1)
    edge = np.argwhere(cv2.dilate(m8, np.ones((3, 3), np.uint8)) != cv2.erode(m8, np.ones((3, 3), np.uint8)))
    if len(edge):
        pick = edge[rng.integers(0, len(edge), len(edge) // 4)]
        m8[pick[:, 0], pick[:, 1]] ^= 1
    return m8.astype(bool)


def corpus(n_riga=6, n_lidc=12) -> list[Item]:
    rng = np.random.default_rng(SEED)
    items = []
    for i in range(n_riga):
        H, W = 240 + 8 * i, 320
        cy, cx, r = int(rng.integers(90, H - 90)), int(rng.integers(90, W - 90)), int(rng.integers(40, 70))
        disc0 = _blob(rng, H, W, cy, cx, r)
        discs, cups = [], []
        for _ in range(6):
            d = _jitter(rng, disc0)
            c = d & _blob(rng, H, W, cy + int(rng.integers(-4, 5)), cx + int(rng.integers(-4, 5)), r // 2)
            discs.append(d)
            cups.append(c)
        image = rng.integers(0, 256, (H, W, 3), dtype=np.uint8)
        items.append(Item(f"riga/syn{i}", "riga", f"riga:syn{i}", "train", image,
                          {"disc": np.stack(discs), "cup": np.stack(cups)}, [str(j + 1) for j in range(6)], {}))
    for i in range(n_lidc):
        H = W = 160
        k = int(rng.integers(2, 5))
        cy, cx, r = int(rng.integers(40, 120)), int(rng.integers(40, 120)), int(rng.integers(2, 16))
        base = _blob(rng, H, W, cy, cx, r)
        if i % 3 == 0 and r >= 8:  # a hole
            base &= ~_ellipse(H, W, cy, cx, r // 3, r // 3)
        if i % 4 == 1:  # a speck
            base[cy + r + 4, cx] = True
        if i % 5 == 2:  # a pinch: two squares touching at a corner
            base[cy - r - 9:cy - r - 5, cx:cx + 4] = True
            base[cy - r - 5:cy - r - 1, cx + 4:cx + 8] = True
        layers = [_jitter(rng, base) for _ in range(k)]
        if i % 6 == 5:  # an abstaining rater: one sub-filter speck
            layers[-1] = np.zeros((H, W), bool)
            layers[-1][cy, cx] = True
        image = rng.integers(-1100, 600, (H, W), dtype=np.int16)
        items.append(Item(f"lidc/syn{i}", "lidc", f"lidc:syn{i}", "train", image, {"nodule": np.stack(layers)},
                          [f"r{j}" for j in range(k)], {}))
    return items


# ---- completions ---------------------------------------------------------------------------------------------------

def _mutate(text, names, how, rng):
    p = parse(text, names)
    if not p.ok:
        return text
    objs = {}
    for n in names:
        parts = []
        for o, hs in p.objs[n]:
            rings = [o, *hs]
            if how == "shift":
                d = rng.integers(-12, 13, 2)
                rings = [np.clip(r + d, 0, 1000) for r in rings]
            elif how == "scale":
                num = int(rng.integers(8, 13))
                c = o.sum(0) // len(o)
                rings = [np.clip(c + (r - c) * num // 10, 0, 1000) for r in rings]
            elif how == "drop":
                rings = [r[::2] if len(r) >= 6 else r for r in rings]
            elif how == "bowtie" and len(o) >= 4:
                o = o.copy()
                o[[1, 2]] = o[[2, 1]]
                rings = [o, *hs]
            parts.append((rings[0], list(rings[1:])))
        objs[n] = parts
    return serialize(objs, names)


def completions(texts, names, rng, n):
    """n completions (text, truncated): rater targets and mutations, garbage, truncation, an empty answer."""
    texts = [t for t in texts if t] or ['{"' + '":[],"'.join(names) + '":[]}']
    hows = ["target", "shift", "scale", "drop", "bowtie", "garbage", "truncated", "empty", "target", "shift"]
    out = []
    for j in range(n):
        how, t = hows[j % len(hows)], texts[j % len(texts)]
        if how == "garbage":
            out.append(("not a polygon " + str(j), False))
        elif how == "truncated":
            out.append((t[:len(t) // 2], True))
        elif how == "empty":
            out.append((serialize({nm: [] for nm in names}, names), False))
        elif how == "target":
            out.append((t, False))
        else:
            out.append((_mutate(t, names, how, rng), False))
    return out


# ---- sections ------------------------------------------------------------------------------------------------------

def run() -> dict:
    rng = np.random.default_rng(SEED + 1)
    items = corpus()
    sec: dict[str, list[str]] = {k: [] for k in (
        "corpus", "frame.masks", "frame.image", "png", "polycore.rasterize", "polycore.validate",
        "polycore.score_batch", "polycore.threads_equal", "topology", "polygonize", "targets", "targets.all_rules",
        "target_set_hash", "caches", "fingerprint", "reward.terms", "reward.advantages", "eval",
        "splits.lidc", "splits.riga", "splits.descriptor", "splits.ncc_values", "splits.ncc_pairs",
        "splits.keypoints", "splits.keypoint_pairs", "splits.keypoint_dups")}

    for it in items:
        sec["corpus"].append(_h(it.item_id, it.image, *(it.masks[n] for n in sorted(it.masks))))

    # polycore on random integer rings (incl. off-grid vertices and self-crossings)
    shapes, Ws, Hs, gts = [], [], [], []
    for _ in range(300):
        W, H = (int(v) for v in rng.choice([50, 100, 125, 200, 250, 500], 2))
        rings = [rng.integers(-50, 1051, (int(rng.integers(3, 40)), 2)) for _ in range(int(rng.integers(1, 4)))]
        r = polycore.rasterize(rings, W, H)
        sec["polycore.rasterize"].append(_h(r["evenodd"], r["nonzero"], r["ring_wmin"], r["ring_wmax"],
                                            *(int(r[k]) for k in ("area_evenodd", "area_nonzero", "px_multi",
                                                                   "rings_mixed_sign"))))
        v = polycore.validate(rings)
        sec["polycore.validate"].append(_h(v["total"], v["per_ring"], v["area2"], v["n_clean"],
                                           v["pair_cross"], v["pair_touch"],
                                           [polycore.signed_area2(rr) for rr in rings]))
        shapes.append(rings)
        Ws.append(W)
        Hs.append(H)
        gts.append(rng.integers(0, 2, (H, W)).astype(bool))
    b1 = polycore.score_batch(shapes, np.array(Ws), np.array(Hs), gts=gts, threads=1)
    b3 = polycore.score_batch(shapes, np.array(Ws), np.array(Hs), gts=gts, threads=3)
    for k in sorted(b1):
        sec["polycore.score_batch"].append(_h(k, np.asarray(b1[k])))
        sec["polycore.threads_equal"].append(_h(k, bool(np.array_equal(b1[k], b3[k]))))

    names_of = {}
    rows = []
    for it in items:
        names = SCHEMAS[it.dataset]
        names_of[it.item_id] = names
        f = frame(it)
        sec["frame.masks"].append(_h(it.item_id, f.W, f.H, f.box, f.scale, *(f.masks[n] for n in names)))
        sec["frame.image"].append(_h(it.item_id, f.image))
        ok, png = cv2.imencode(".png", cv2.cvtColor(f.image, cv2.COLOR_RGB2BGR))  # as dataset.prepare writes it
        sec["png"].append(_h(it.item_id, bool(ok), png))
        stacks = dict(f.masks) | ({"rim": f.masks["disc"] & ~f.masks["cup"]} if it.dataset == "riga" else {})
        ma = {n: min_area(int(s.any(0).sum())) for n, s in stacks.items()}
        for n, s in stacks.items():
            for r, m in enumerate(s):
                t = instance_topology(m)
                sec["topology"].append(_h(it.item_id, n, r, t.areas, t.holes, t.hole_areas, betti(m),
                                          [signature(m, a) for a in (0, 5, 20, ma[n])], size_filter(m, ma[n])))
                sec["polygonize"].append(_h(it.item_id, n, r, [(rg.pts, rg.hole) for rg in polygonize(m)]))
        texts = []
        for r, rater in enumerate(it.raters):
            src = {n: f.masks[n][r] for n in names}
            for key, kw in (("targets", ma), ("targets.all_rules", None)):
                try:
                    tg = build_target(src, names, BuildConfig(), kw)
                    info = {k: v for k, v in tg.info.items() if k != "time"}
                    sec[key].append(_h(it.item_id, rater, tg.text, json.dumps(_exact(info), sort_keys=True)))
                    if kw is ma:
                        texts.append(tg.text)
                        rows.append((it.item_id, rater, tg.text))
                except ValueError as e:
                    sec[key].append(_h(it.item_id, rater, "error", str(e)))
                    if kw is ma:
                        texts.append(None)
        caches = {n: ExpertCache.build(f.masks[n]) for n in names}
        for n, c in caches.items():
            sec["caches"].append(_h(it.item_id, n, c.tau, c.t5, c.min_area, c.rater_sigs, c.voters, c.box,
                                    c.perimeter, c.dist_out_union, c.dist_in_core, c.digest))
        comps = completions(texts, names, rng, 16)
        jobs = [(t, tr, names, caches, f.W, f.H) for t, tr in comps]
        terms = score_batch(jobs)
        numeric = {k: v for k, v in terms.items() if isinstance(v, np.ndarray)}
        sec["reward.terms"].append(_h(it.item_id, *(numeric[k] for k in sorted(numeric)), terms["error"],
                                      json.dumps(terms["defects"], sort_keys=True)))
        sec["reward.advantages"].append(_h(it.item_id, advantages_v4(terms), topology_reward(terms)))
        m = metrics([t for t, _ in comps], {n: f.masks[n] for n in names}, f.W, f.H, names,
                    [tr for _, tr in comps], caches)
        sec["eval"].append(_h(it.item_id, json.dumps(_exact(m), sort_keys=True)))
    sec["target_set_hash"].append(target_set_hash(rows))
    sec["fingerprint"].append(fingerprint())

    # splits: hash-dealt assignments, FFT-based NCC near-duplicates, ORB + RANSAC keypoint duplicates
    pats = {f"P{i:04d}": splits.Patient(int(rng.integers(0, 60)), int(rng.integers(0, 3)) * int(rng.integers(0, 2)))
            for i in range(300)}
    sec["splits.lidc"].append(_h(sorted(splits.lidc_split(pats).items())))
    subsets = {f"img{i}": riga_mod.SUBSETS[i % len(riga_mod.SUBSETS)] for i in range(120)}
    groups = splits.group_ids(list(subsets), [(f"img{i}", f"img{i + 1}") for i in range(0, 119, 7)])
    sec["splits.riga"].append(_h(sorted(splits.riga_split(subsets, groups).items())))
    photos = []
    for i in range(8):  # fundus-like photographs: noisy field of view with vessel-like random lines; pairs 2j, 2j+1
        if i % 2 == 0:  # are the same eye reframed (shifted crop of one larger texture)
            big = np.zeros((560, 760, 3), np.uint8)
            big[..., 0] = 120
            big[..., 1] = rng.integers(30, 90, (560, 760), dtype=np.uint8)
            for _ in range(60):
                p0, p1 = (tuple(int(v) for v in rng.integers(0, 760, 2)) for _ in range(2))
                cv2.line(big, p0, p1, (150, 200, 40), int(rng.integers(1, 4)))
        dy, dx = (0, 0) if i % 2 == 0 else (int(rng.integers(3, 12)), int(rng.integers(3, 12)))
        img = np.ascontiguousarray(big[dy:dy + 480, dx:dx + 640])
        img[~_ellipse(480, 640, 240, 320, 230, 300)] = 0
        photos.append(img)
    desc = np.stack([riga_mod.descriptor(p) for p in photos])
    for d in desc:
        sec["splits.descriptor"].append(_h(d))
    ncc = splits.near_duplicates(desc, threshold=-1.0)  # every pair, to see raw values
    sec["splits.ncc_values"].append(_h([(i, j, float(c)) for i, j, c in ncc]))
    sec["splits.ncc_pairs"].append(_h([(i, j) for i, j, c in ncc if c >= splits.NCC_DUP]))
    feats = [riga_mod.keypoints(p) for p in photos]
    for xy, des in feats:
        sec["splits.keypoints"].append(_h(xy, des if des is not None else "none"))
    kp = splits.keypoint_duplicates(feats, threshold=0)  # every pair's inlier count
    sec["splits.keypoint_pairs"].append(_h(kp))
    sec["splits.keypoint_dups"].append(_h([(i, j) for i, j, n in kp if n >= splits.ORB_DUP]))  # the decisions

    return {"platform": {"machine": platform.machine(), "system": platform.system(), "python": sys.version.split()[0],
                         "numpy": np.__version__, "opencv": cv2.__version__,
                         "cpu_simd": _simd()},
            "sections": {k: {"sha256": hashlib.sha256("".join(v).encode()).hexdigest(), "n": len(v),
                             "records": [x[:16] for x in v]} for k, v in sec.items()}}


def _simd() -> dict:
    out = {}
    try:
        from numpy._core._multiarray_umath import __cpu_baseline__, __cpu_dispatch__, __cpu_features__
        out["numpy_baseline"] = list(__cpu_baseline__)
        out["numpy_dispatch_active"] = [f for f in __cpu_dispatch__ if __cpu_features__.get(f)]
    except ImportError:
        pass
    out["opencv_optimized"] = bool(cv2.useOptimized())
    out["opencv_ipp"] = bool(cv2.ipp.useIPP())
    return out


def compare(a_path, b_path, strict=False) -> int:
    a, b = (json.loads(Path(p).read_text()) for p in (a_path, b_path))
    print(f"A: {a['platform']['machine']} {a['platform']['system']}   B: {b['platform']['machine']} "
          f"{b['platform']['system']}")
    bad = acc = 0
    for k in a["sections"]:
        sa, sb = a["sections"][k], b["sections"].get(k)
        same = sb is not None and sa["sha256"] == sb["sha256"]
        first = "" if same or sb is None else next(
            (f" first differing record {i}" for i, (x, y) in enumerate(zip(sa["records"], sb["records"],
                                                                        strict=False)) if x != y), "")
        n_diff = 0 if same or sb is None else sum(x != y for x, y in zip(sa["records"], sb["records"], strict=False))
        ok = same or (not strict and k in ACCEPTED)
        print(f"{'OK  ' if same else 'ACPT' if ok else 'DIFF'} {k:24s} n={sa['n']:4d} {sa['sha256'][:16]}"
              + ("" if same else f" vs {sb['sha256'][:16] if sb else 'missing'} ({n_diff} records){first}"
                 + (f"\n       known cause: {CAUSES[k]}" if k in CAUSES else "")))
        bad += not ok
        acc += ok and not same
    print(f"{len(a['sections']) - bad - acc} / {len(a['sections'])} sections identical"
          + (f", {acc} accepted ({', '.join(ACCEPTED)})" if acc else "") + ("" if bad else ": PASS"))
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", help="write the JSON report here (default: stdout summary only)")
    ap.add_argument("--compare", nargs=2, metavar=("A", "B"),
                    help="compare two reports; exit 1 on a difference outside ACCEPTED (raw split scores)")
    ap.add_argument("--strict", action="store_true", help="with --compare: accept no difference (run-to-run checks)")
    ap.add_argument("--opencv", choices=("default", "no-ipp", "no-opt"), default="default",
                    help="diagnostic: run with OpenCV's Intel IPP HAL off (cv2.ipp.setUseIPP(False)) or with every "
                         "optimised path off (cv2.setUseOptimized(False))")
    args = ap.parse_args(argv)
    if args.compare:
        return compare(*args.compare, strict=args.strict)
    if args.opencv == "no-ipp":
        cv2.ipp.setUseIPP(False)
    elif args.opencv == "no-opt":
        cv2.setUseOptimized(False)
    rep = run()
    for k, v in rep["sections"].items():
        print(f"{v['sha256']}  {k} (n={v['n']})")
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
