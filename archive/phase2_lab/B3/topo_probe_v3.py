"""Synthetic topology probe v3: hole and part geometry drawn from the real LIDC train targets (P5.2-P5.4; judge
MUST-FIX 1). Outlines are generator 2 (probe_lib, E = 0.2, A = 0.15, calibrated in precision_probe run 2).

    .venv/bin/python -I -m experiments.phase2.lab.B3.topo_probe_v3 --records <scratch>/hole_records.jsonl \
        --census <scratch>/census.jsonl --out <lab>/topo_probe_v3.json

- hole item (v3c): resample one TRAIN hole record (its part's filled diameter, hole area, centroid offset, wall);
  outline = generator-2 star of that diameter, dented (or bulged) toward the hole with a Gaussian in angle so the
  wall is the record's; hole = a small star of that area at the record's offset; accepted when the ground truth
  has exactly one hole (signature (1,) raw and at min_area(|GT|), so every hole is >= the size rule), area within
  max(1 px, 10 %) and in the record's bin (v3e), wall within 0.45 px (v3d), offset within 0.15. Earlier versions failed
  fidelity: v3a (offset placement: walls too thick, areas too large; topo_probe_v3a.json), v3b (wall placement:
  offsets too large; topo_probe_v3b.json), v3c (wall tolerance 0.8 px: diagonal-touch walls lost;
  topo_probe_v3c.json), v3d (area +-25 % + 1 px: +1 px area bias, failed the fresh-seed confirmation;
  topo_probe_v3d_seed5/6.json). v3e is confirmed on two fresh seeds (7, 8) because v3a-v3d were tuned on seed 5;
- two-part item: resample one TRAIN multi-part record (object diameter, part-area ratio, gap); two generator-2 stars
  with the record's areas, centres separated so the gap is about the record's; accepted when the signature is
  (0, 0);
- fidelity (P5.2): a natural (unstratified) draw of 300 hole and 200 two-part items, measured with
  hole_census.geometry, KS against the train records; bin shares;
- probe (P5.3): 120 one-hole items per hole-area bin (<= 10, 11-30, > 30 px; resampling records within the bin),
  160 two-part items, 160 simple controls (generator 2, fourier / lobulated);
- sanity (P5.4): exact, A0 (holes emptied), centroid disc, sigma = 0.5 px (K = 8 draws per item: single draw,
  mean over 8, any of 8), per bin, pooled = train-share-weighted mean (post-stratified).
"""
import argparse
import json

import numpy as np
from scipy import stats

from objdel.checker import signature, topology
from objdel.geometry import parse, to_mask
from objdel.targets import build_target

from . import prereg
from . import probe_lib as P
from .hole_census import geometry
from .precision_probe import train_diams

SHAPE = (0.2, 0.15)
WALL_TOL = 0.45  # v3d/v3e; 0.8 in v3c let 0.414 px (diagonal-touch) walls drift to 1.0 px


def bin_of(area):
    lo, hi = prereg.P5_BINS
    return 0 if area <= lo else 1 if area <= hi else 2


def _finish(gt, family, extra):
    ma = topology.min_area(int(gt.sum()))
    want = (1,) if family == "hole" else (0, 0)
    if signature(gt) != want or signature(gt, ma) != want:
        return None
    try:
        t = build_target({"nodule": gt}, ("nodule",), min_area={"nodule": ma})
    except ValueError:
        return None
    return {"family": family, "gt": gt, "min_area": ma, "text": t.text, "objs": t.objs["nodule"],
            "diam": float(2 * np.sqrt(gt.sum() / np.pi)), **extra}


def make_hole(rec, rng, tries=100):
    """v3c: the hole sits at the record's centroid offset and the outline is dented (or bulged) toward it so the
    wall is the record's: real holed outlines come close to their holes (median wall 1.2 px at median offset
    0.46 r). v3a placed by offset only (walls about 2 px too thick), v3b by wall only (offsets too large)."""
    d, a, wall, off = rec["part_diam"], rec["area"], rec["wall"], rec["offset"]
    for _ in range(tries):
        outer = P._star(0, 0, d / 2, SHAPE, rng)
        th = np.arctan2(outer[:, 1], outer[:, 0])
        r = np.hypot(outer[:, 0], outer[:, 1])
        ang = rng.uniform(-np.pi, np.pi)
        hr = np.sqrt(a / np.pi) * rng.uniform(0.85, 1.1)
        o = off * d / 2 * rng.uniform(0.9, 1.1)
        target = o + hr + wall + rng.uniform(0.3, 1.0)
        w = max(hr + 1.5, 2.0) / max(o + hr, 2.0)
        dth = np.angle(np.exp(1j * (th - ang)))
        r = r - (r - target) * np.exp(-0.5 * (dth / w) ** 2)
        outer = np.c_[r * np.cos(th), r * np.sin(th)]
        inner = P._star(o * np.cos(ang), o * np.sin(ang), hr, (0.1, 0.1), rng, n=32, spikes=False)
        lo, hi = outer.min(0), outer.max(0)
        side = int(np.ceil((hi - lo).max())) + 1
        S = P._frame_side(side, rng)
        if side + 6 > S or S > 500:
            continue
        x0, y0 = rng.uniform(2, S - 2 - (hi - lo)[0]), rng.uniform(2, S - 2 - (hi - lo)[1])
        gt = P._compose([(1, outer), (-1, inner)], S, x0 - lo[0], y0 - lo[1])
        if gt is None:
            continue
        g = geometry(gt)
        if g["n_holes"] != 1 or g["n_parts"] != 1:
            continue
        h = g["holes"][0]
        if (bin_of(h["area"]) != bin_of(a) or abs(h["area"] - a) > max(1.0, 0.10 * a) or abs(h["wall"] - wall) > WALL_TOL
                or abs(h["offset"] - off) > 0.15):
            continue
        it = _finish(gt, "hole", {"S": S, "hole_area": h["area"], "bin": bin_of(h["area"]), "rec_area": a})
        if it is not None:
            return it
    return None


def make_two(rec, rng, tries=60):
    D, ratio, gap = rec["diam"], rec["part_ratio"], rec["gap"]
    total = np.pi * (D / 2) ** 2
    a1 = total / (1 + ratio)
    d1, d2 = 2 * np.sqrt(a1 / np.pi), 2 * np.sqrt(ratio * a1 / np.pi)
    for _ in range(tries):
        ang = rng.uniform(0, 2 * np.pi)
        sep = d1 / 2 + d2 / 2 + max(gap, 1.0) * rng.uniform(0.9, 1.2) + 1
        p1 = P._star(0, 0, d1 / 2, SHAPE, rng, spikes=False)
        p2 = P._star(sep * np.cos(ang), sep * np.sin(ang), max(d2 / 2, 1.3), SHAPE, rng, spikes=False)
        pts = np.concatenate([p1, p2])
        lo, hi = pts.min(0), pts.max(0)
        side = int(np.ceil((hi - lo).max())) + 1
        S = P._frame_side(side, rng)
        if side + 6 > S or S > 500:
            continue
        x0, y0 = rng.uniform(2, S - 2 - (hi - lo)[0]), rng.uniform(2, S - 2 - (hi - lo)[1])
        gt = P._compose([(1, p1), (1, p2)], S, x0 - lo[0], y0 - lo[1])
        if gt is None:
            continue
        it = _finish(gt, "two", {"S": S})
        if it is not None:
            return it
    return None


def draw(recs, n, maker, rng):
    out, fails = [], 0
    while len(out) < n:
        it = maker(recs[rng.integers(len(recs))], rng)
        if it is None:
            fails += 1
            continue
        out.append(it)
    return out, fails


def pred_sig(item, text):
    p = parse(text, ("nodule",))
    if not p.ok:
        return None
    g = to_mask(p.objs["nodule"], item["S"], item["S"])
    return signature(g.mask, item["min_area"]) if g.legal else None


def ci(x, rng, B=4000):
    x = np.asarray(x, float)
    e = x[rng.integers(0, len(x), (B, len(x)))].mean(1)
    return [round(float(x.mean()), 4), round(float(np.percentile(e, 2.5)), 4), round(float(np.percentile(e, 97.5)), 4)]


def pooled_ci(per_bin, w, rng, B=4000):
    """Post-stratified mean sum_b w_b mean_b, stratified bootstrap."""
    est = np.zeros(B)
    point = 0.0
    for x, wb in zip(per_bin, w, strict=True):
        x = np.asarray(x, float)
        est += wb * x[rng.integers(0, len(x), (B, len(x)))].mean(1)
        point += wb * x.mean()
    return [round(float(point), 4), round(float(np.percentile(est, 2.5)), 4), round(float(np.percentile(est, 97.5)), 4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=5)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    with open(a.records) as fh:
        recs = [json.loads(line) for line in fh]
    train = [r for r in recs if r["split"] == "train"]
    holes = [h for r in train for h in r["holes"]]
    multis = [r for r in train if r["n_parts"] >= 2]
    shares = np.bincount([bin_of(h["area"]) for h in holes], minlength=3) / len(holes)
    res = {"train_holes": len(holes), "train_multipart_targets": len(multis),
           "train_bin_shares": [round(float(s), 4) for s in shares]}

    # P5.2 fidelity: natural draw
    nat_h, fh_ = draw(holes, 300, make_hole, rng)
    nat_t, ft_ = draw(multis, 200, make_two, rng)
    gh = [geometry(it["gt"])["holes"][0] for it in nat_h]
    gt2 = [geometry(it["gt"]) for it in nat_t]
    ks = {
        "hole_area": stats.ks_2samp([h["area"] for h in gh], [h["area"] for h in holes]).pvalue,
        "hole_fraction": stats.ks_2samp([h["fraction"] for h in gh], [h["fraction"] for h in holes]).pvalue,
        "hole_offset": stats.ks_2samp([h["offset"] for h in gh], [h["offset"] for h in holes]).pvalue,
        "hole_wall": stats.ks_2samp([h["wall"] for h in gh], [h["wall"] for h in holes]).pvalue,
        "part_ratio": stats.ks_2samp([g["part_ratio"] for g in gt2], [r["part_ratio"] for r in multis]).pvalue,
        "part_gap": stats.ks_2samp([g["gap"] for g in gt2], [r["gap"] for r in multis]).pvalue,
    }
    gen_shares = np.bincount([bin_of(h["area"]) for h in gh], minlength=3) / len(gh)
    q = (10, 50, 90)
    res["fidelity"] = {
        "ks_p": {k: round(float(v), 4) for k, v in ks.items()},
        "gen_bin_shares": [round(float(s), 4) for s in gen_shares],
        "generated": {"hole_area": [float(np.percentile([h["area"] for h in gh], p)) for p in q],
                      "hole_fraction": [round(float(np.percentile([h["fraction"] for h in gh], p)), 3) for p in q],
                      "hole_offset": [round(float(np.percentile([h["offset"] for h in gh], p)), 3) for p in q],
                      "hole_wall": [round(float(np.percentile([h["wall"] for h in gh], p)), 3) for p in q],
                      "part_ratio": [round(float(np.percentile([g["part_ratio"] for g in gt2], p)), 3) for p in q],
                      "part_gap": [round(float(np.percentile([g["gap"] for g in gt2], p)), 3) for p in q]},
        "real_train": {"hole_area": [float(np.percentile([h["area"] for h in holes], p)) for p in q],
                       "hole_fraction": [round(float(np.percentile([h["fraction"] for h in holes], p)), 3) for p in q],
                       "hole_offset": [round(float(np.percentile([h["offset"] for h in holes], p)), 3) for p in q],
                       "hole_wall": [round(float(np.percentile([h["wall"] for h in holes], p)), 3) for p in q],
                       "part_ratio": [round(float(np.percentile([r["part_ratio"] for r in multis], p)), 3) for p in q],
                       "part_gap": [round(float(np.percentile([r["gap"] for r in multis], p)), 3) for p in q]},
        "rejections": {"hole": fh_, "two": ft_}}
    res["fidelity"]["pass"] = bool(min(ks.values()) >= prereg.P5_KS_P
                                   and np.abs(gen_shares - shares).max() <= prereg.P5_SHARE_TOL)
    print("fidelity", json.dumps(res["fidelity"]), flush=True)

    # P5.3 probe: stratified hole items, two-part, controls
    probe_holes = []
    for b in range(3):
        rb = [h for h in holes if bin_of(h["area"]) == b]
        its, _ = draw(rb, prereg.P5_PER_BIN, make_hole, rng)
        probe_holes += its
    two, _ = draw(multis, 160, make_two, rng)
    diams, _ = train_diams(a.census)
    ctrl = []
    while len(ctrl) < 160:
        fam = ("fourier", "lobulated")[len(ctrl) % 2]
        it = P.make_item(fam, float(rng.choice(diams)), SHAPE, rng)
        if it is not None:
            ctrl.append(it)
    res["probe"] = {"hole_items_per_bin": [sum(it["bin"] == b for it in probe_holes) for b in range(3)],
                    "hole_area_per_bin_median": [float(np.median([it["hole_area"] for it in probe_holes
                                                                  if it["bin"] == b])) for b in range(3)],
                    "holed_object_diam_median": float(np.median([it["diam"] for it in probe_holes])),
                    "two_part": len(two), "controls": len(ctrl)}

    # P5.4 sanity models
    K = 8
    models = {"exact": (P.model_exact, 1), "A0_holes_emptied": (P.model_a0, 1),
              "centroid_disc": (P.model_centroid_disc, 1), "sigma_0.5": (P.model_noise(0.5), K)}
    res["models"] = {}
    for name, (m, k) in models.items():
        mrng = np.random.default_rng(17)
        hit = {"single": [[], [], []], "mean8": [[], [], []], "any8": [[], [], []]}
        for it in probe_holes:
            h = [bool((s := pred_sig(it, m(it, mrng))) and sum(s) > 0) for _ in range(k)]
            hit["single"][it["bin"]].append(h[0])
            hit["mean8"][it["bin"]].append(np.mean(h))
            hit["any8"][it["bin"]].append(any(h))
        e = {}
        for metric, per_bin in hit.items():
            if k == 1 and metric != "single":
                continue
            e[f"hole_recall_{metric}"] = {"bins": [ci(x, mrng) for x in per_bin],
                                          "pooled_train_weighted": pooled_ci(per_bin, shares, mrng)}
        parts = [bool((s := pred_sig(it, m(it, mrng))) and len(s) > 1) for it in two]
        fp = [bool((s := pred_sig(it, m(it, mrng))) is not None and s != (0,)) for it in ctrl]
        e["part_recall_descriptive"] = ci(parts, mrng)
        e["false_nonsimple_controls"] = ci(fp, mrng)
        res["models"][name] = e
        print(name, json.dumps(e), flush=True)
    s = res["models"]
    res["sanity"] = {
        "exact_all_1": all(b[0] == 1.0 for b in s["exact"]["hole_recall_single"]["bins"])
        and s["exact"]["false_nonsimple_controls"][0] == 0.0,
        "A0_all_0": all(b[0] == 0.0 for b in s["A0_holes_emptied"]["hole_recall_single"]["bins"]),
        "disc_all_0": all(b[0] == 0.0 for b in s["centroid_disc"]["hole_recall_single"]["bins"])}
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res["sanity"]), json.dumps(res["probe"]))


if __name__ == "__main__":
    main()
