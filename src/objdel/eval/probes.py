"""Synthetic precision (D19) and topology (D20) probes: build them as eval-view directories, score sampled answers.

    stats = lidc_stats("data/store")                                   # train diameters, hole and part records
    build_precision(out, patches("data/store", "val"), stats["diams"], n=1000, seed=2026)
    build_topology(out, patches("data/store", "val"), stats["holes"], stats["multis"], stats["diams"], seed=2027)
    score_precision(rows, out); score_topology(rows, out)             # rows: {"item_id", "samples": [text, ...]}

Promoted from archive/phase2_lab/B3 (probe_lib, precision_probe, probe_rule, topo_probe_v3 = generator v3e, hole_census,
power_c2) and archive/phase2_lab/A4/probe_proto; REPORT.md P2 and Addendum A hold the evidence.

Generator 2 (D19): star polygons r(t) = r0 (1 + E e cos(2t + p2) + A sum_{k=3..12} c_k cos(kt + pk)), E = 0.2,
A = 0.15 (calibrated: centroid-disc IoU 0.822 vs real 0.806), spicules with p = 0.3, rasterized on the 0-1000 grid by
objdel.geometry.to_mask; frame side = the targets.frame divisor rule with a margin drawn from [8 px, 60 %] (scale
jitter), object placed uniformly with >= 2 px clearance (translation jitter). Families of the precision probe:
fourier 60 % / lobulated 30 % / hole 5 % / two 5 %, diameters resampled from the LIDC train voting raters.
Topology probe (D20, v3e): every hole jointly resamples one train hole record (part diameter, hole area, centroid
offset, wall; accepted at area +- max(1 px, 10 %) inside the record's bin, wall +- 0.45 px, offset +- 0.15), every
two-part item one train multi-part record (diameter, part ratio, gap); every hole >= the size rule (signature (1,)
raw and at min_area(|GT|)). P5.2 fidelity failed narrowly (seed 7: hole-area KS p 0.004; seed 8 passes): a ~1 px
area under-draw of 5-px holes in the smallest bin.

Directory = the data/prepared layout (images/<item_id>.png, caches/, lidc/<split>.jsonl, manifest.json); each
record has one "rater" ("gt") whose certified target text is pixel-exact to the ground truth (checked at build), plus
"view" and "probe" (family, diam, bin, hole_area). Images: the object composited on a real LIDC lung patch
(partial-volume blur, solid or ground-glass HU, CT noise), lung window, Pillow bicubic upscale to 512 px.
Precision: lidc/val.jsonl = probe-test (1,000 items on patches of LIDC val patients, never seen in SFT training),
lidc/train.jsonl = probe-train (the 20k capacity set, built on demand, on patches of LIDC train patients).

Metrics (frame px): boundary = foreground pixels with a 4-adjacent background pixel; ASSD pooled over both
directions; HD95 = max of the directed 95th percentiles; centroid error between pixel-centre centroids.
score_precision (D19, §5.1 G-precision) on 6-16 px items, every sample: legal >= 95 %; medians over legal samples
ASSD <= 0.6 px, IoU >= 0.80, centroid <= 0.5 px; per-item mean IoU (illegal = 0) beats the equal-area centroid disc
on >= 60 % of items. score_topology (D20): per hole item the mean over its samples (all rows of an item pooled, so
seeds pool) of "has a hole", per bin, pooled = train-share-weighted bins with a within-bin bootstrap (B = 2,000);
gate (a) = simple controls answered (1, 0) >= 95 %; hole recall is reported against 0.30.
"""

import hashlib
import json
import time
from pathlib import Path

import cv2
import numpy as np
from PIL import Image
from scipy import ndimage, stats

from objdel.checker import signature, size_filter
from objdel.checker.topology import min_area
from objdel.geometry import GRID, canonicalize, parse, serialize, to_mask
from objdel.reward.cache import boundary  # foreground pixels with a 4-adjacent background pixel
from objdel.targets import build_target
from objdel.targets.frame import DIVISORS, LUNG_WINDOW

SHAPE = (0.2, 0.15)
MIX = (("fourier", 0.60), ("lobulated", 0.30), ("hole", 0.05), ("two", 0.05))
SIG = {"fourier": (0,), "lobulated": (0,), "hole": (1,), "two": (0, 0)}
BINS = (10, 30)  # hole area px: <= 10 | 11-30 | > 30
BIN_WEIGHTS = (0.4386, 0.3614, 0.2000)  # train hole shares (hole_census.json; 125 / 103 / 57 of 285)
PRECISION_RULE = {"legal": 0.95, "assd": 0.6, "iou": 0.80, "centroid": 0.5, "win": 0.60, "diam": (6.0, 16.0)}
TOPOLOGY_RULE = {"controls": 0.95, "recall_report": 0.30}
WALL_TOL = 0.45
C4, C8 = ndimage.generate_binary_structure(2, 1), ndimage.generate_binary_structure(2, 2)
NAMES = ("nodule",)


# ------------------------------------------------------------------------------------------------- generator
def _star(cx, cy, r0, shape, rng, n=96, spikes=True):
    E, A = shape
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rr = 1 + E * rng.uniform(0.5, 1) * np.cos(2 * th + rng.uniform(0, 2 * np.pi))
    for k in range(3, 13):
        rr += A * rng.uniform(0, 1) / k * np.cos(k * th + rng.uniform(0, 2 * np.pi))
    if spikes and rng.random() < 0.3:
        for _ in range(rng.integers(1, 5)):
            dth = np.angle(np.exp(1j * (th - rng.uniform(0, 2 * np.pi))))
            rr += rng.uniform(0.2, 0.6) * max(A, E) * np.exp(-0.5 * (dth / 0.08) ** 2)
    rr = np.maximum(rr, 0.35)
    return np.c_[cx + r0 * rr * np.cos(th), cy + r0 * rr * np.sin(th)]


def _ring(poly_px, S):
    g = np.clip(np.rint(np.asarray(poly_px) * (GRID / S)).astype(np.int64), 0, GRID)
    return g[np.any(g != np.roll(g, 1, 0), 1)]


def _compose(prims, S, dx, dy):
    """Union (+1) / carve (-1) of pixel-unit polygons, each rasterized by the checker; None if one is illegal."""
    m = np.zeros((S, S), bool)
    for sign, poly in prims:
        g = _ring(poly + [dx, dy], S)
        geo = to_mask(canonicalize({"x": [(g, [])]})["x"], S, S) if len(g) >= 3 else None
        if geo is None or not geo.legal:
            return None
        m = (m | geo.mask) if sign > 0 else (m & ~geo.mask)
    return m


def _place(polys, rng):
    """Frame side (divisor rule, margin in [8 px, 60 %]) and the placed offset; None when it does not fit."""
    pts = np.concatenate(polys)
    lo, hi = pts.min(0), pts.max(0)
    side = int(np.ceil((hi - lo).max())) + 1
    need = max(side + 2 * max(8, int(np.ceil(rng.uniform(0.0, 0.60) * side))), 50)  # one draw (not per divisor)
    S = next(x for x in DIVISORS if x >= need)
    if side + 6 > S or S > 500:
        return None
    x0, y0 = rng.uniform(2, S - 2 - (hi - lo)[0]), rng.uniform(2, S - 2 - (hi - lo)[1])
    return S, x0 - lo[0], y0 - lo[1]


def _prims(family, d, rng):
    (E, A), r0 = SHAPE, d / 2
    if family == "fourier":
        return [(1, _star(0, 0, r0, SHAPE, rng))]
    if family == "lobulated":
        out = []
        for _ in range(rng.integers(2, 5)):
            a, rr, off = rng.uniform(0, 2 * np.pi), r0 * rng.uniform(0.5, 0.85), r0 * rng.uniform(0.2, 0.55)
            out.append((1, _star(off * np.cos(a), off * np.sin(a), rr, (E * 0.5, A * 0.3), rng, spikes=False)))
        return out
    if family == "hole":
        hr = r0 * rng.uniform(0.25, 0.45)
        return [
            (1, _star(0, 0, r0, (E, A * 0.5), rng, spikes=False)),
            (
                -1,
                _star(
                    rng.uniform(-0.15, 0.15) * r0,
                    rng.uniform(-0.15, 0.15) * r0,
                    hr,
                    (0.1, 0.2),
                    rng,
                    n=48,
                    spikes=False,
                ),
            ),
        ]
    d1 = d * rng.uniform(0.55, 0.8)  # "two"
    d2 = np.sqrt(max(d * d - d1 * d1, 9.0))
    a = rng.uniform(0, 2 * np.pi)
    sep = (d1 + d2) / 2 * (1 + 0.5 * max(A, E)) + rng.uniform(2.5, 6)
    return [
        (1, _star(-sep / 2 * np.cos(a), -sep / 2 * np.sin(a), d1 / 2, SHAPE, rng, spikes=False)),
        (1, _star(sep / 2 * np.cos(a), sep / 2 * np.sin(a), d2 / 2, SHAPE, rng, spikes=False)),
    ]


def _finish(gt, family, S, **extra):
    """The item when its ground truth has the family's signature raw and at min_area(|GT|) and certifies."""
    if gt is None or not gt.any():
        return None
    ma = min_area(int(gt.sum()))
    if signature(gt) != SIG[family] or signature(gt, ma) != SIG[family]:
        return None
    try:
        text = build_target({"nodule": gt}, NAMES, min_area={"nodule": ma}).text
    except ValueError:
        return None
    return {"family": family, "S": S, "gt": gt, "text": text, "diam": float(2 * np.sqrt(gt.sum() / np.pi)), **extra}


def make_item(family, d, rng, tries=50):
    """Generator-2 item of a family and drawn diameter d (precision probe; simple controls)."""
    for _ in range(tries):
        prims = _prims(family, d, rng)
        p = _place([q for s, q in prims if s > 0], rng)
        if p is not None and (it := _finish(_compose(prims, p[0], p[1], p[2]), family, p[0])) is not None:
            return it
    return None


def make_hole(rec, rng, tries=100):
    """v3e: the hole at the record's offset, the outline dented toward it so the wall is the record's."""
    d, a, wall, off = rec["part_diam"], rec["area"], rec["wall"], rec["offset"]
    for _ in range(tries):
        outer = _star(0, 0, d / 2, SHAPE, rng)
        th, r = np.arctan2(outer[:, 1], outer[:, 0]), np.hypot(outer[:, 0], outer[:, 1])
        ang = rng.uniform(-np.pi, np.pi)
        hr = np.sqrt(a / np.pi) * rng.uniform(0.85, 1.1)
        o = off * d / 2 * rng.uniform(0.9, 1.1)
        target = o + hr + wall + rng.uniform(0.3, 1.0)
        w = max(hr + 1.5, 2.0) / max(o + hr, 2.0)
        r = r - (r - target) * np.exp(-0.5 * (np.angle(np.exp(1j * (th - ang))) / w) ** 2)
        outer = np.c_[r * np.cos(th), r * np.sin(th)]
        inner = _star(o * np.cos(ang), o * np.sin(ang), hr, (0.1, 0.1), rng, n=32, spikes=False)
        p = _place([outer], rng)
        if p is None or (gt := _compose([(1, outer), (-1, inner)], *p)) is None:
            continue
        g = geometry(gt)
        if g["n_holes"] != 1 or g["n_parts"] != 1:
            continue
        h = g["holes"][0]
        if (
            bin_of(h["area"]) != bin_of(a)
            or abs(h["area"] - a) > max(1.0, 0.10 * a)
            or abs(h["wall"] - wall) > WALL_TOL
            or abs(h["offset"] - off) > 0.15
        ):
            continue
        if (it := _finish(gt, "hole", p[0], hole_area=h["area"], bin=bin_of(h["area"]))) is not None:
            return it
    return None


def make_two(rec, rng, tries=60):
    D, ratio, gap = rec["diam"], rec["part_ratio"], rec["gap"]
    a1 = np.pi * (D / 2) ** 2 / (1 + ratio)
    d1, d2 = 2 * np.sqrt(a1 / np.pi), 2 * np.sqrt(ratio * a1 / np.pi)
    for _ in range(tries):
        ang = rng.uniform(0, 2 * np.pi)
        sep = d1 / 2 + d2 / 2 + max(gap, 1.0) * rng.uniform(0.9, 1.2) + 1
        p1 = _star(0, 0, d1 / 2, SHAPE, rng, spikes=False)
        p2 = _star(sep * np.cos(ang), sep * np.sin(ang), max(d2 / 2, 1.3), SHAPE, rng, spikes=False)
        p = _place([p1, p2], rng)
        if p is not None and (it := _finish(_compose([(1, p1), (1, p2)], *p), "two", p[0])) is not None:
            return it
    return None


def _draw(recs, n, maker, rng):
    out = []
    while len(out) < n:
        if (it := maker(recs[rng.integers(len(recs))], rng)) is not None:
            out.append(it)
    return out


def bin_of(area):
    return 0 if area <= BINS[0] else 1 if area <= BINS[1] else 2


def geometry(m):
    """Parts (8-conn.) and holes (4-conn. background not touching the border): area, fraction of the filled part,
    centroid offset / equivalent radius, wall (min distance to outside the filled part - 1), part ratio, gap."""
    parts, n_parts = ndimage.label(m, structure=C8)
    bg, n = ndimage.label(~m, structure=C4)
    border = set(np.unique(np.r_[bg[0], bg[-1], bg[:, 0], bg[:, -1]]).tolist())
    keep = [i for i in range(1, n + 1) if i not in border]
    rec = {
        "area": int(m.sum()),
        "diam": float(2 * np.sqrt(m.sum() / np.pi)),
        "n_parts": n_parts,
        "n_holes": len(keep),
        "holes": [],
    }
    for i in keep:
        hm = bg == i
        ring = ndimage.binary_dilation(hm, structure=C4) & ~hm
        filled = ndimage.binary_fill_holes(parts == int(np.bincount(parts[ring][parts[ring] > 0]).argmax()))
        (fy, fx), (hy, hx) = np.nonzero(filled), np.nonzero(hm)
        r_eq = np.sqrt(filled.sum() / np.pi)
        dist_out = ndimage.distance_transform_edt(np.pad(filled, 1))[1:-1, 1:-1]
        rec["holes"].append(
            {
                "area": int(hm.sum()),
                "fraction": float(hm.sum() / filled.sum()),
                "offset": float(np.hypot(hx.mean() - fx.mean(), hy.mean() - fy.mean()) / r_eq),
                "wall": float(dist_out[hm].min() - 1),
                "part_diam": float(2 * r_eq),
            }
        )
    if n_parts >= 2:
        sizes = np.bincount(parts.ravel())[1:]
        a, b = np.argsort(sizes)[::-1][:2] + 1
        rec.update(
            part_ratio=float(sizes[b - 1] / sizes[a - 1]),
            gap=float(ndimage.distance_transform_edt(parts != a)[parts == b].min() - 1),
        )
    return rec


def fidelity(holes, multis, seed, n_holes=300, n_two=200):
    """P5.2 on a natural draw: two-sample KS p of six geometry variables vs the records; hole-area bin shares."""
    rng = np.random.default_rng(seed)
    gh = [geometry(it["gt"])["holes"][0] for it in _draw(holes, n_holes, make_hole, rng)]
    g2 = [geometry(it["gt"]) for it in _draw(multis, n_two, make_two, rng)]
    ks = {
        f"hole_{k}": stats.ks_2samp([h[k] for h in gh], [h[k] for h in holes]).pvalue
        for k in ("area", "fraction", "offset", "wall")
    }
    ks |= {
        f"part_{k.split('_')[-1]}": stats.ks_2samp([g[k] for g in g2], [r[k] for r in multis]).pvalue
        for k in ("part_ratio", "gap")
    }
    share = np.bincount([bin_of(h["area"]) for h in gh], minlength=3) / len(gh)
    real = np.bincount([bin_of(h["area"]) for h in holes], minlength=3) / len(holes)
    return {
        "ks_p": {k: round(float(v), 4) for k, v in ks.items()},
        "gen_bin_shares": [round(float(s), 4) for s in share],
        "pass": bool(min(ks.values()) >= 0.01 and np.abs(share - real).max() <= 0.05),
    }


# ------------------------------------------------------------------------------------------------- real data
def lidc_stats(store_root):
    """Train diameters of the voting raters (>= 2.5 px), the median per frame side, and the hole / multi-part
    geometry records of the non-simple targets after the size rule
    (archive/phase2_lab/B3 census + hole_census, same order)."""
    from objdel.data import store
    from objdel.reward.cache import voters
    from objdel.targets import frame

    diams, by_w, holes, multis = [], {}, [], []
    for it in store.load("lidc", "train", store_root):
        f = frame(it)
        st = f.masks["nodule"]
        ma = min_area(int(st.any(0).sum()))
        sigs = [signature(m, ma) for m in st]
        vote = voters(sigs)
        for i in vote:
            if (d := 2 * np.sqrt(st[i].sum() / np.pi)) >= 2.5:
                diams.append(float(d))
                by_w.setdefault(f.W, []).append(float(d))
        if any(sigs[i] != (0,) for i in vote):
            for i in vote:
                g = geometry(size_filter(st[i], ma))
                holes += g["holes"]
                multis += [g] if g["n_parts"] >= 2 else []
    return {
        "diams": np.array(diams),
        "prior": {w: float(np.median(v)) for w, v in by_w.items()},
        "holes": holes,
        "multis": multis,
    }


def patches(store_root, split, n=200, side=125, per_patient=3, seed=11):
    """[(group_id, HU patch)]: lung-parenchyma windows away from the nodules
    (archive/phase2_lab/B3 background_patches)."""
    from objdel.data import store

    rng, out, count = np.random.default_rng(seed), [], {}
    for it in store.load("lidc", split, store_root):
        g = it.group_id
        if len(out) >= n or count.get(g, 0) >= per_patient:
            continue
        hu, busy = it.image, ndimage.binary_dilation(it.masks["nodule"].any(0), iterations=12)
        for _ in range(20):
            y, x = rng.integers(0, hu.shape[0] - side), rng.integers(0, hu.shape[1] - side)
            win = hu[y : y + side, x : x + side]
            if (
                not busy[y : y + side, x : x + side].any()
                and -920 < np.median(win) < -750
                and (win < -980).mean() < 0.02
                and (win > -300).mean() < 0.05
            ):
                out.append((g, win.copy()))
                count[g] = count.get(g, 0) + 1
                break
    return out


# ------------------------------------------------------------------------------------------------- writer
def render(gt, patch, rng, image_side=512):
    """uint8 RGB model image: the object on a lung patch (reflect-padded when S exceeds it), lung window, upscale."""
    S = gt.shape[0]
    if S > min(patch.shape):
        patch = np.pad(patch, [(0, max(S - patch.shape[0], 0)), (0, max(S - patch.shape[1], 0))], mode="reflect")
    y0, x0 = rng.integers(0, patch.shape[0] - S + 1), rng.integers(0, patch.shape[1] - S + 1)
    bg = patch[y0 : y0 + S, x0 : x0 + S].astype(float)
    hu = rng.uniform(-100, 100) if rng.random() < 0.7 else rng.uniform(-700, -300)
    alpha = ndimage.gaussian_filter(gt.astype(float), rng.uniform(0.5, 1.0))
    noisy = bg * (1 - alpha) + hu * alpha + rng.normal(0, rng.uniform(10, 30), bg.shape)
    lo, hi = LUNG_WINDOW
    u8 = np.clip(np.rint((noisy - lo) / (hi - lo) * 255), 0, 255).astype(np.uint8)
    img = Image.fromarray(np.repeat(u8[..., None], 3, -1)).resize((image_side,) * 2, Image.Resampling.BICUBIC)
    return np.asarray(img)


def _pixels_sha(image):  # as objdel.train.dataset._pixels_sha (not imported: that module needs `datasets`)
    return hashlib.sha256(repr(image.shape).encode() + image.tobytes()).hexdigest()


def _write(out, name, split, items, patch_list, seed, meta):
    out, rng, recs = Path(out), np.random.default_rng(seed + 1), []
    for k, it in enumerate(items):
        if not (to_mask(parse(it["text"], NAMES).objs["nodule"], it["S"], it["S"]).mask == it["gt"]).all():
            raise AssertionError("certified target not pixel-exact")  # never seen in 1,200 B3 items
        g, patch = patch_list[rng.integers(len(patch_list))]
        img = render(it["gt"], patch, rng)
        iid = f"lidc/{name}/{split}/{k:05d}"
        path = out / "images" / (iid + ".png")
        path.parent.mkdir(parents=True, exist_ok=True)
        cv2.imwrite(str(path), cv2.cvtColor(img, cv2.COLOR_RGB2BGR))
        S = it["S"]
        recs.append(
            {
                "item_id": iid,
                "group_id": g,
                "W": S,
                "H": S,
                "box": [0, 0, S, S],
                "scale": 1,
                "raters": ["gt"],
                "targets": [it["text"]],
                "rl": False,
                "disagreement": False,
                "image_sha256": _pixels_sha(img),
                "caches": {},
                "errors": [],
                "view": {"offset": [0, 0], "side": S},
                "probe": {k2: it[k2] for k2 in ("family", "diam", "bin", "hole_area") if k2 in it},
            }
        )
    (out / "caches").mkdir(exist_ok=True)
    (out / "lidc").mkdir(exist_ok=True)
    text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in recs)
    (out / "lidc" / f"{split}.jsonl").write_text(text)
    mpath = out / "manifest.json"
    man = json.loads(mpath.read_text()) if mpath.is_file() else {"files": {}, "sets": {}}
    man["files"][f"lidc/{split}.jsonl"] = {"sha256": hashlib.sha256(text.encode()).hexdigest()}
    fam = {f: sum(r["probe"]["family"] == f for r in recs) for f in sorted({r["probe"]["family"] for r in recs})}
    man["sets"][f"lidc/{split}"] = {
        "items": len(recs),
        "patients": len({r["group_id"] for r in recs}),
        "families": fam,
        "seed": seed,
        **meta,
    }
    mpath.write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
    return man


def build_precision(out, patch_list, diams, n=1000, seed=2026, split="val"):
    """Precision probe: probe-test (split "val", patches of LIDC val patients) or probe-train ("train", the capacity
    set, patches of LIDC train patients)."""
    t0, rng = time.time(), np.random.default_rng(seed)
    fams, p = [f for f, _ in MIX], np.array([w for _, w in MIX])
    items = []
    while len(items) < n:
        fam, d = fams[rng.choice(len(fams), p=p / p.sum())], float(rng.choice(diams))
        if not (fam in ("hole", "two") and d < 8.0) and (it := make_item(fam, d, rng)) is not None:
            items.append(it)
    return _write(
        out,
        "probe-precision",
        split,
        items,
        patch_list,
        seed,
        {"generator": {"E": SHAPE[0], "A": SHAPE[1], "mix": dict(MIX)}, "build_s": time.time() - t0},
    )


def build_topology(out, patch_list, holes, multis, diams, seed=2027, per_bin=120, n_two=160, n_simple=160):
    """Topology probe (D20): per_bin one-hole items per hole-area bin, n_two two-part items, n_simple controls."""
    t0, rng = time.time(), np.random.default_rng(seed)
    items = [
        it for b in range(3) for it in _draw([h for h in holes if bin_of(h["area"]) == b], per_bin, make_hole, rng)
    ]
    items += _draw(multis, n_two, make_two, rng)
    ctrl = []
    # generator-2 controls, fourier / lobulated alternating (archive/phase2_lab/B3 topo_probe_v3)
    while len(ctrl) < n_simple:
        if (it := make_item(("fourier", "lobulated")[len(ctrl) % 2], float(rng.choice(diams)), rng)) is not None:
            ctrl.append(it | {"family": "simple"})
    return _write(
        out,
        "probe-topology",
        "val",
        items + ctrl,
        patch_list,
        seed,
        {"generator": "v3e", "per_bin": per_bin, "bin_weights": BIN_WEIGHTS, "build_s": time.time() - t0},
    )


# ------------------------------------------------------------------------------------------------- scoring
def records(root, split="val"):
    with open(Path(root, "lidc", f"{split}.jsonl")) as fh:
        return [json.loads(line) for line in fh]


def gt_mask(rec):
    return to_mask(parse(rec["targets"][0], NAMES).objs["nodule"], rec["W"], rec["W"]).mask


def surface(gt, pred):
    """Pooled ASSD, HD95, NSD@1 between the boundaries; None when either mask is empty."""
    bg, bp = boundary(gt), boundary(pred)
    if not bg.any() or not bp.any():
        return None
    dg, dp = ndimage.distance_transform_edt(~bp)[bg], ndimage.distance_transform_edt(~bg)[bp]
    d = np.concatenate([dg, dp])
    return {
        "assd": float(d.mean()),
        "hd95": float(max(np.percentile(dg, 95), np.percentile(dp, 95))),
        "nsd1": float((d <= 1).mean()),
    }


def centroid(m):
    ys, xs = np.nonzero(m)
    return np.array([xs.mean() + 0.5, ys.mean() + 0.5])


def _pred(rec, text):
    p = parse(text, NAMES)
    g = to_mask(p.objs["nodule"], rec["W"], rec["W"]) if p.ok else None
    return g.mask if g is not None and g.legal else None


def score_answer(rec, text, gt=None):
    """legal, IoU, ASSD, HD95, centroid error of one answer (illegal: legal False only)."""
    gt, pr = gt_mask(rec) if gt is None else gt, _pred(rec, text)
    if pr is None:
        return {"legal": False}
    out = {"legal": True, "iou": np.count_nonzero(gt & pr) / np.count_nonzero(gt | pr)}
    if not pr.any():
        return out | {"assd": np.inf, "hd95": np.inf, "centroid": np.inf}
    return out | surface(gt, pr) | {"centroid": float(np.linalg.norm(centroid(pr) - centroid(gt)))}


def score_precision(rows, root="data/eval_views/probe-precision", split="val", rule=PRECISION_RULE, capacity=None):
    """G-precision on the zero-shot rows; `capacity` = the capacity run's rows on the same items, scored alike
    under "capacity" (gates.phase2_entry: pass, else capacity pass with a warning)."""
    if capacity is not None:
        return score_precision(rows, root, split, rule) | {"capacity": score_precision(capacity, root, split, rule)}
    recs = {r["item_id"]: r for r in records(root, split)}
    lo, hi = rule["diam"]
    samples, wins = [], []
    for row in rows:
        rec = recs[row["item_id"]]
        if not lo <= rec["probe"]["diam"] < hi:
            continue
        gt = gt_mask(rec)
        s = [score_answer(rec, t, gt) for t in row["samples"]]
        disc = score_answer(rec, model_centroid_disc(rec, None), gt).get("iou", 0.0)
        wins.append(np.mean([a.get("iou", 0.0) for a in s]) > disc)
        samples += s
    lg = [a for a in samples if a["legal"]]
    med = {
        k: float(np.median([a[k] for a in lg])) if lg else (0.0 if k == "iou" else np.inf)
        for k in ("assd", "iou", "centroid", "hd95")
    }
    m = {
        "n_items": len(wins),
        "n_samples": len(samples),
        "legal": float(np.mean([a["legal"] for a in samples])),
        **med,
        "win_share": float(np.mean(wins)),
    }
    m["checks"] = {
        "legal": m["legal"] >= rule["legal"],
        "assd": m["assd"] <= rule["assd"],
        "iou": m["iou"] >= rule["iou"],
        "centroid": m["centroid"] <= rule["centroid"],
        "win": m["win_share"] >= rule["win"],
    }
    return m | {"pass": all(m["checks"].values())}


def stratified_ci(per_bin, weights, rng, B=2000):
    """Weighted mean of the bin means and its within-bin percentile bootstrap 95 % CI
    (archive/phase2_lab/B3 power_c2.pooled; also c2_analyze's pooled CI).
    Rounded to 12 decimals: an exact 288 / 960 hits sums to 0.29999999999999993 in floats and would miss the
    C2 floor 0.30 (BK2-4; the same tolerance as objdel.eval.signal.DIGITS)."""
    per_bin = [np.asarray(x, float) for x in per_bin]
    est = sum(w * x[rng.integers(0, len(x), (B, len(x)))].mean(1) for w, x in zip(weights, per_bin, strict=True))
    lo, hi = np.percentile(est, [2.5, 97.5])
    point = sum(w * x.mean() for w, x in zip(weights, per_bin, strict=True))
    return round(float(point), 12), round(float(lo), 12), round(float(hi), 12)


def score_topology(rows, root="data/eval_views/probe-topology", split="val", B=2000, seed=0, rule=TOPOLOGY_RULE):
    recs = {r["item_id"]: r for r in records(root, split)}
    texts = {}
    for row in rows:  # rows of one item (seeds) pool
        texts.setdefault(row["item_id"], []).extend(row["samples"])
    mean, anyk, part, ok, false_ns, legal = [[], [], []], [[], [], []], [], [], [], []
    for iid in (i for i in recs if i in texts):  # record order, so per-seed "bins" align item by item
        rec, ts = recs[iid], texts[iid]
        ma = min_area(int(gt_mask(rec).sum()))
        sigs = [signature(p, ma) if (p := _pred(rec, t)) is not None else None for t in ts]
        legal += [s is not None for s in sigs]
        fam = rec["probe"]["family"]
        if fam == "hole":
            h = [bool(s) and sum(s) > 0 for s in sigs]
            mean[rec["probe"]["bin"]].append(np.mean(h))
            anyk[rec["probe"]["bin"]].append(float(any(h)))
        elif fam == "two":
            part.append(np.mean([s is not None and len(s) > 1 for s in sigs]))
        else:
            ok.append(np.mean([s == (0,) for s in sigs]))
            false_ns.append(np.mean([s is not None and s != (0,) for s in sigs]))
    pooled = stratified_ci(mean, BIN_WEIGHTS, np.random.default_rng(seed), B)
    res = {
        "bins": [[float(v) for v in x] for x in mean],  # per-item mean hit (c2_analyze input)
        "hole_recall": {
            "bins": [float(np.mean(x)) for x in mean],
            "n_bins": [len(x) for x in mean],
            "pooled": list(pooled),
            "any_of_k": list(stratified_ci(anyk, BIN_WEIGHTS, np.random.default_rng(seed), B)),
            "vs_report_line": pooled[0] >= rule["recall_report"],
        },
        "part_recall": float(np.mean(part)),
        "controls_correct": float(np.mean(ok)),
        "false_nonsimple": float(np.mean(false_ns)),
        "legal": float(np.mean(legal)),
    }
    return res | {
        "gate": {"controls": res["controls_correct"] >= rule["controls"]},
        "pass": res["controls_correct"] >= rule["controls"],  # rule (a) only; hole recall is report-only (C2)
    }


# ------------------------------------------------------------------------------------------------- synthetic models
def _text(parts):
    return serialize(canonicalize({"nodule": parts}), NAMES)


def _disc_text(cx, cy, r, S, n=48):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return _text([(_ring(np.c_[cx + r * np.cos(th), cy + r * np.sin(th)], S), [])])


def model_centroid_disc(rec, rng):
    gt = gt_mask(rec)
    return _disc_text(*centroid(gt), np.sqrt(gt.sum() / np.pi), rec["W"])
