"""Precision probe (A4 section 6) and synthetic topology probe (H-A1-7): calibration, exactness, sensitivity,
discrimination. Pass criteria: prereg.py (P2).

    .venv/bin/python -I -m experiments.phase2.lab.B3.precision_probe --census <scratch>/census.jsonl \
        --out <lab>/precision_probe.json --scratch <scratch>

Steps:
1. real reference (A4 S3): equal-area disc at each LIDC reader's centroid, val every 3rd item (A4's sample) and train
   every 7th item, scored with probe_lib's metrics;
2. calibration grid over (E, A) (H-A4-7; probe_lib docstring) on 500 items each (seed 100), family mix fourier
   60 % / lobulated 30 % / hole 5 % / two 5 %; the chosen point minimises the distance to the real disc IoU and the
   real targets' median vertex count. The vertex-noise legality curve is NOT fitted: it is the out-of-sample check;
3. a fresh probe set at A* (N items, new seed): P2.1 exactness, P2.2 brute force, P2.3 1-px shifts, P2.4 noise curve,
   P2.5 discrimination by synthetic models, plus the diameter-stratified table;
4. rendering on real LIDC val lung patches (sanity: the noise-free composite thresholded half-way recovers the
   ground truth); a few example images to the scratch directory;
5. the synthetic topology probe: 1/3 one-hole, 1/3 two-part, 1/3 simple items; metrics of exact, A0 (holes
   emptied), centroid-disc and sigma = 0.5 answers with 95 % bootstrap CIs.
Reads untrusted data (store): run with python -I.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image

from objdel.checker import signature, topology
from objdel.data import store
from objdel.geometry import parse, to_mask
from objdel.targets import frame

from . import prereg
from . import probe_lib as P

MIX = (("fourier", 0.60), ("lobulated", 0.30), ("hole", 0.05), ("two", 0.05))
BINS = (0, 6, 10, 16, 25, 1e9)


def train_diams(census):
    d, by_w = [], {}
    with open(census) as fh:
        for line in fh:
            r = json.loads(line)
            if r["split"] != "train":
                continue
            for i in r["voters"]:
                if r["rater_diam"][i] >= 2.5:
                    d.append(r["rater_diam"][i])
                    by_w.setdefault(r["W"], []).append(r["rater_diam"][i])
    return np.array(d), {w: float(np.median(v)) for w, v in by_w.items()}


def real_reference(split, every, root="data/store"):
    rng = np.random.default_rng(0)
    ious = []
    for i, it in enumerate(store.load("lidc", split, root)):
        if i % every:
            continue
        f = frame(it)
        st = f.masks["nodule"]
        for m in st:
            if m.sum() < 5:
                continue
            item = {"gt": m, "S": f.W, "min_area": topology.min_area(int(st.any(0).sum()))}
            ious.append(P.score(item, P.model_centroid_disc(item, rng)).get("iou", 0.0))
    return {"n": len(ious), "median": float(np.median(ious)), "p10": float(np.percentile(ious, 10))}


def draw_items(n, shape, diams, rng, mix=MIX, min_d_topo=8.0):
    fams = [f for f, _ in mix]
    p = np.array([w for _, w in mix])
    items = []
    while len(items) < n:
        fam = fams[rng.choice(len(fams), p=p / p.sum())]
        d = float(rng.choice(diams))
        if fam in ("hole", "two") and d < min_d_topo:
            continue
        it = P.make_item(fam, d, shape, rng)
        if it is not None:
            items.append(it)
    return items


def run_model(items, model, rng):
    rows = []
    for it in items:
        s = P.score(it, model(it, rng))
        rows.append({**s, "diam": it["diam"], "family": it["family"]})
    return rows


def summary(rows):
    lg = [r for r in rows if r["legal"]]
    q = lambda k: float(np.median([r[k] for r in lg])) if lg else None  # noqa: E731
    return {"n": len(rows), "legal": float(np.mean([r["legal"] for r in rows])), "iou": q("iou"), "assd": q("assd"),
            "hd95": q("hd95"), "nsd1": q("nsd1"), "centroid": q("centroid")}


def by_bins(rows):
    out = {}
    for lo, hi in zip(BINS[:-1], BINS[1:], strict=True):
        rs = [r for r in rows if lo <= r["diam"] < hi]
        if rs:
            out[f"{lo}-{hi if hi < 1e9 else 'inf'}"] = summary(rs)
    return out


def boot_ci(x, rng, B=2000):
    x = np.asarray(x, float)
    if not len(x):
        return None
    est = x[rng.integers(0, len(x), (B, len(x)))].mean(1)
    return [round(float(x.mean()), 4), round(float(np.percentile(est, 2.5)), 4),
            round(float(np.percentile(est, 97.5)), 4)]


def background_patches(n, rng, root="data/store", side=125):
    from scipy import ndimage
    out = []
    for it in store.load("lidc", "val", root):
        if len(out) >= n:
            break
        hu = it.image
        busy = ndimage.binary_dilation(it.masks["nodule"].any(0), iterations=12)
        for _ in range(20):
            y, x = rng.integers(0, hu.shape[0] - side), rng.integers(0, hu.shape[1] - side)
            win = hu[y:y + side, x:x + side]
            if busy[y:y + side, x:x + side].any():
                continue
            # lung parenchyma only: not outside-body air (about -1000 HU), table, chest wall or mediastinum
            if -920 < np.median(win) < -750 and (win < -980).mean() < 0.02 and (win > -300).mean() < 0.05:
                out.append(win.copy())
                break
    return out


def topo_probe(diams, shape, rng, n_each=160):
    items = []
    for fam, k in (("hole", n_each), ("two", n_each), ("fourier", n_each // 2), ("lobulated", n_each - n_each // 2)):
        items += draw_items(k, shape, diams, rng, mix=((fam, 1.0),))
    models = {"exact": P.model_exact, "A0_holes_emptied": P.model_a0, "centroid_disc": P.model_centroid_disc,
              "sigma_0.5": P.model_noise(0.5)}
    out = {"items": {f: sum(it["family"] == f for it in items) for f in ("hole", "two", "fourier", "lobulated")}}
    for name, m in models.items():
        hit, hole_rec, part_rec, fp, hole_pred_ok, part_pred_ok = [], [], [], [], [], []
        for it in items:
            text = m(it, rng)
            p = parse(text, ("nodule",))
            pred_sig = None
            if p.ok:
                g = to_mask(p.objs["nodule"], it["S"], it["S"])
                pred_sig = signature(g.mask, it["min_area"]) if g.legal else None
            nonsimple = it["family"] in ("hole", "two")
            has_hole = bool(pred_sig) and sum(pred_sig) > 0
            multi = bool(pred_sig) and len(pred_sig) > 1
            if nonsimple:
                hit.append(pred_sig == P.SIG[it["family"]])
            else:
                fp.append(pred_sig is not None and pred_sig != (0,))
            if it["family"] == "hole":
                hole_rec.append(has_hole)
            if it["family"] == "two":
                part_rec.append(multi)
            if has_hole:
                hole_pred_ok.append(it["family"] == "hole")
            if multi:
                part_pred_ok.append(it["family"] == "two")
        out[name] = {"nonsimple_hit": boot_ci(hit, rng), "hole_recall": boot_ci(hole_rec, rng),
                     "hole_precision": boot_ci(hole_pred_ok, rng), "part_recall": boot_ci(part_rec, rng),
                     "part_precision": boot_ci(part_pred_ok, rng), "simple_false_nonsimple": boot_ci(fp, rng)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--n", type=int, default=1200)
    ap.add_argument("--grid-e", default="0.1,0.2,0.3,0.4,0.5")
    ap.add_argument("--grid-a", default="0.0,0.15,0.3,0.45,0.6")
    a = ap.parse_args()
    t0 = time.time()
    res = {}
    diams, prior = train_diams(a.census)
    res["real_diam_px"] = {q: round(float(np.percentile(diams, q)), 2) for q in (10, 50, 90)}
    res["prior_diam_by_S"] = prior
    res["real_reference"] = {"val_every3": real_reference("val", 3), "train_every7": real_reference("train", 7)}
    target = res["real_reference"]["val_every3"]["median"]
    print("real", res["real_reference"], f"{time.time() - t0:.0f}s", flush=True)

    # 2. calibration on (E, A): disc IoU and the certified targets' median vertex count (real LIDC train targets)
    verts = []
    with open("data/targets/lidc/train.jsonl") as fh:
        for line in fh:
            o = json.loads(json.loads(line)["text"])["nodule"]
            verts.append(sum(len(q["outer"]) + sum(len(h) for h in q["holes"]) for q in o))
    v_target = float(np.median(verts))
    res["real_target_vertices"] = {q: float(np.percentile(verts, q)) for q in (10, 50, 90)}
    cal = {}
    for E in [float(x) for x in a.grid_e.split(",")]:
        for A in [float(x) for x in a.grid_a.split(",")]:
            rng = np.random.default_rng(100)
            its = draw_items(500, (E, A), diams, rng)
            rows = run_model(its, P.model_centroid_disc, rng)
            nv = [sum(len(r) for o, hs in it["objs"] for r in (o, *hs)) for it in its]
            c = {"disc_iou_median": float(np.median([r.get("iou", 0) for r in rows])),
                 "verts_median": float(np.median(nv)), "verts_p90": float(np.percentile(nv, 90))}
            c["loss"] = ((c["disc_iou_median"] - target) / prereg.P2_DISC_IOU[1]) ** 2 + \
                ((c["verts_median"] - v_target) / 2.0) ** 2
            cal[f"{E},{A}"] = c
            print("E,A", E, A, c, flush=True)
    key = min(cal, key=lambda k: cal[k]["loss"])
    shape = tuple(float(x) for x in key.split(","))
    res["calibration"] = {"grid": cal, "shape_star": shape, "target_disc": target, "target_vertices": v_target,
                          "loss": "((disc - target) / 0.03)^2 + ((vertices - target) / 2)^2"}

    # 3. fresh probe set
    rng = np.random.default_rng(2026)
    items = draw_items(a.n, shape, diams, rng)
    res["probe_set"] = {"n": len(items), "families": {f: sum(it["family"] == f for it in items) for f, _ in MIX},
                        "diam_px": {q: round(float(np.percentile([it["diam"] for it in items], q)), 2)
                                    for q in (10, 50, 90)},
                        "S": {str(s): sum(it["S"] == s for it in items) for s in sorted({it["S"] for it in items})},
                        "target_vertices": {q: float(np.percentile([sum(len(r) for o, hs in it["objs"]
                                                                        for r in (o, *hs)) for it in items], q))
                                            for q in (50, 90)},
                        "gen_s": round(time.time() - t0, 1)}
    # P2.1 exactness
    ex_raster = [P.surface(it["gt"], it["gt"]) for it in items]
    ex_text = run_model(items, P.model_exact, rng)
    pixel_exact = [bool((to_mask(it["objs"], it["S"], it["S"]).mask == it["gt"]).all()) for it in items]
    res["P2_1_exact"] = {
        "raster_all_zero": all(e["assd"] == 0 and e["hd95"] == 0 and e["nsd1"] == 1 for e in ex_raster),
        "text_pixel_exact_share": float(np.mean(pixel_exact)),
        "text": summary(ex_text),
        "text_all_perfect": all(r["legal"] and r["iou"] == 1 and r["assd"] == 0 and r["hd95"] == 0
                                and r["centroid"] == 0 and r["topology"] for r in ex_text)}
    # P2.2 brute force
    maxdiff, npairs = 0.0, 0
    for it in items:
        for m in (P.model_noise(1.0), P.model_shift, P.model_centroid_disc):
            g = to_mask(parse(m(it, rng), ("nodule",)).objs["nodule"], it["S"], it["S"])
            if not g.mask.any():
                continue
            s1, s2 = P.surface(it["gt"], g.mask), P.surface_brute(it["gt"], g.mask)
            maxdiff = max(maxdiff, abs(s1["assd"] - s2["assd"]), abs(s1["hd95"] - s2["hd95"]))
            npairs += 1
    res["P2_2_brute"] = {"pairs": npairs, "max_abs_diff": maxdiff, "pass": npairs >= 1000 and maxdiff <= prereg.P2_BRUTE_TOL}
    # P2.3 1-px shifts in 4 directions
    sh = [P.score(it, P.shift_text(it, dx, dy)) for it in items for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))]
    res["P2_3_shift"] = {"n": len(sh), "assd_pos_share": float(np.mean([s["assd"] > 0 for s in sh])),
                         "centroid_eq_1_share": float(np.mean([abs(s["centroid"] - 1) < 1e-9 for s in sh])),
                         "assd_median": float(np.median([s["assd"] for s in sh])),
                         "iou_median": float(np.median([s["iou"] for s in sh]))}
    # P2.4 / P2.5 models
    models = {"exact": P.model_exact, "sigma_0.25": P.model_noise(0.25), "sigma_0.5": P.model_noise(0.5),
              "sigma_1": P.model_noise(1.0), "sigma_2": P.model_noise(2.0), "centroid_disc": P.model_centroid_disc,
              "centred_prior": P.model_centred_prior(prior), "shift_1px": P.model_shift}
    res["models"] = {}
    for name, m in models.items():
        rows = run_model(items, m, np.random.default_rng(7))
        res["models"][name] = {"all": summary(rows), "bins": by_bins(rows),
                               "pass_rule_6_16": P.pass_rule(rows, prereg.PASS_RULE, prereg.PASS_BINS)}
        print(name, res["models"][name]["all"], res["models"][name]["pass_rule_6_16"], flush=True)
    cal_ok = {}
    for sig, (iou_ref, legal_ref) in prereg.P2_NOISE_REF.items():
        s = res["models"][f"sigma_{sig:g}"]["all"]
        cal_ok[f"sigma_{sig:g}"] = {"iou": s["iou"], "legal": s["legal"],
                                    "ok": abs(s["iou"] - iou_ref) <= prereg.P2_NOISE_TOL[0]
                                    and abs(s["legal"] - legal_ref) <= prereg.P2_NOISE_TOL[1]}
    disc = res["models"]["centroid_disc"]["all"]["iou"]
    res["P2_4_calibration"] = {"disc_iou": disc, "target": target,
                               "disc_ok": abs(disc - prereg.P2_DISC_IOU[0]) <= prereg.P2_DISC_IOU[1],
                               "noise": cal_ok}
    res["P2_5_discrimination"] = {
        "should_pass": {n: res["models"][n]["pass_rule_6_16"]["pass"] for n in prereg.P2_SHOULD_PASS},
        "should_fail": {n: not res["models"][n]["pass_rule_6_16"]["pass"] for n in prereg.P2_SHOULD_FAIL}}

    # 4. rendering
    rng = np.random.default_rng(11)
    bgs = background_patches(200, rng)
    thr, out_dir = [], Path(a.scratch) / "probe_examples"
    out_dir.mkdir(parents=True, exist_ok=True)
    for k, it in enumerate(items):
        img, clean, hu, bg_med = P.render(it, bgs[k % len(bgs)], rng)
        if hu > -300:  # solid: the half-way iso-level recovers the boundary (inside the object's box + 2 px,
            # so that background vessels elsewhere in the frame do not count)
            ys, xs = np.nonzero(it["gt"])
            win = (slice(max(ys.min() - 2, 0), ys.max() + 3), slice(max(xs.min() - 2, 0), xs.max() + 3))
            rec, g = clean[win] > (hu + bg_med) / 2, it["gt"][win]
            thr.append(np.count_nonzero(rec & g) / np.count_nonzero(rec | g))
        if k < 12:
            Image.fromarray(img).save(out_dir / f"item{k:02d}_{it['family']}_S{it['S']}.png")
    res["render"] = {"backgrounds": len(bgs), "solid_items": len(thr),
                     "threshold_iou_median": float(np.median(thr)), "threshold_iou_p10": float(np.percentile(thr, 10))}
    print("render", res["render"], flush=True)

    # 5. synthetic topology probe
    res["topology_probe"] = topo_probe(diams, shape, np.random.default_rng(31))
    print("topology", json.dumps(res["topology_probe"]), flush=True)
    res["rendered_examples"] = str(out_dir)
    res["wall_s"] = round(time.time() - t0, 1)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("P2_1_exact", "P2_2_brute", "P2_3_shift", "P2_4_calibration",
                                         "P2_5_discrimination")}, indent=1, default=float))


if __name__ == "__main__":
    main()
