"""R1: single-spec round trip over every fixture target (all raters, hybrid and lossless modes).

    .venv/bin/python -I r1_roundtrip.py [nproc=3] [limit_items=0]

target builder -> production serialize -> strict parse -> polycore validate + rasterize (crop size) ->
reward-v2 terms vs the item's rater set. Also: polycore raster == fastgeom raster bit for bit, and
canonicalization invariance (rotate ring start; reverse then re-fix winding, new start).
"""
import json
import os
import sys
import time
import zlib
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402

import prod  # noqa: E402
from targets.builder import RIGA_DERIVED, BuildConfig, build_target, rasterize_parts  # noqa: E402
from targets.crops import apply_crop  # noqa: E402

FIX = prod.HERE.parent / "fixtures"
MODES = {"hybrid": BuildConfig(), "lossless": BuildConfig("minseg_exact", {}, iou_min=1.0)}
CROP = "pad_div"  # chosen crop policy (t_targets H6)


def items():
    out = []
    z = np.load(FIX / "riga_experts.npz")
    for it in sorted({k.rsplit("__", 1)[0] for k in z.files}):
        out.append(("riga", it))
    z = np.load(FIX / "lidc_experts.npz")
    for k in sorted(z.files):
        out.append(("lidc", k.rsplit("__", 1)[0]))
    return out


_Z = {}


def load(ds, it):
    if ds not in _Z:
        _Z[ds] = np.load(FIX / f"{ds}_experts.npz")
    z = _Z[ds]
    if ds == "riga":
        return {"disc": z[f"{it}__disc"].astype(bool), "cup": z[f"{it}__cup"].astype(bool)}
    return {"nodule": z[f"{it}__nodule"].astype(bool)}


def perturb(parts, rng, how):
    out = []
    for o, hs in parts:
        rs = []
        for r, hole in [(o, False), *[(h, True) for h in hs]]:
            r = np.asarray(r)
            if how == "rotate":
                r = np.roll(r, -int(rng.integers(1, len(r))), axis=0)
            else:  # reverse, then re-fix winding (outer clockwise on screen = positive area2), then new start
                r = r[::-1]
                if (prod.polycore.signed_area2(r) > 0) == hole:
                    r = r[::-1]
                r = np.roll(r, -int(rng.integers(0, len(r))), axis=0)
            rs.append(r)
        out.append((rs[0], rs[1:]))
    return out


def eval_text(text, names, caches, W, H):
    p = prod.parse(text, names)
    assert p.ok, p.error
    res, masks = {}, {}
    for n in names:
        m, legal, defects = prod.geometry(p.objs[n], W, H)
        t = prod.terms(m, caches[n])
        masks[n] = m
        res[n] = {"legal": legal, "defects": defects, **t, "sig": tuple(t["sig"])}
    return p, res, masks


def run(job):
    ds, it = job
    stacks = load(ds, it)
    names = list(stacks)
    k = next(iter(stacks.values())).shape[0]
    # crop every rater with the same policy; caches per structure (and the derived RIGA rim)
    per_rater = [apply_crop({n: stacks[n][r] for n in names}, CROP) for r in range(k)]
    H, W = per_rater[0][1]
    cstack = {n: np.stack([pr[0][n] for pr in per_rater]) for n in names}
    caches = {n: prod.build_cache(cstack[n]) for n in names}
    rim_cache = prod.build_cache(cstack["disc"] & ~cstack["cup"]) if ds == "riga" else None
    rng = np.random.default_rng(zlib.crc32(it.encode()))
    rows = []
    for mode, cfg in MODES.items():
        for r in range(k):
            masks = per_rater[r][0]
            t0 = time.process_time()
            _, objs, info = build_target(masks, {n: "parts" for n in names}, cfg, RIGA_DERIVED if ds == "riga" else ())
            row = {"ds": ds, "item": it, "rater": r, "mode": mode, "W": W, "H": H,
                   "status": {n: objs[n].status for n in names}, "build_s": time.process_time() - t0}
            parts = {n: objs[n].parts for n in names}
            text = prod.serialize(parts, names)
            row["chars"] = len(text)
            p, res, pm = eval_text(text, names, caches, W, H)
            # exact round trip of the parts
            row["roundtrip_exact"] = all(
                len(p.objs[n]) == len(parts[n]) and all(
                    np.array_equal(a[0], b[0]) and len(a[1]) == len(b[1]) and all(
                        np.array_equal(x, y) for x, y in zip(a[1], b[1])) for a, b in zip(p.objs[n], parts[n]))
                for n in names)
            row["has_holes_key"] = text.count('"outer"') == text.count('"holes"')
            for n in names:
                fg = rasterize_parts(parts[n], W, H) if parts[n] else np.zeros((H, W), bool)
                row[f"{n}_raster_eq_fastgeom"] = bool(np.array_equal(np.asarray(fg, bool), pm[n]))
                row[f"{n}_legal"] = res[n]["legal"]
                row[f"{n}_defects"] = res[n]["defects"]
                row[f"{n}_sig_eq"] = res[n]["sig"] == caches[n].rater_sigs[r]
                row[f"{n}_sig"] = str(res[n]["sig"])
                row[f"{n}_src_sig"] = str(caches[n].rater_sigs[r])
                row[f"{n}_iou_tc"] = res[n]["iou_tc"]
                row[f"{n}_f5"] = res[n]["f5"]
                row[f"{n}_iou_src"] = float((pm[n] & cstack[n][r]).sum() / max(1, (pm[n] | cstack[n][r]).sum()))
            if ds == "riga":
                rim = pm["disc"] & ~pm["cup"]
                tr = prod.terms(rim, rim_cache)
                row["rim_sig_eq"] = tuple(tr["sig"]) == rim_cache.rater_sigs[r]
                row["rim_iou_tc"], row["rim_f5"] = tr["iou_tc"], tr["f5"]
            # canonicalization invariance
            inv = True
            for how in ("rotate", "reverse_fix"):
                text2 = prod.serialize({n: perturb(parts[n], rng, how) for n in names}, names)
                _, res2, pm2 = eval_text(text2, names, caches, W, H)
                for n in names:
                    a, b = res[n], res2[n]
                    same = (a["legal"] == b["legal"] and a["sig"] == b["sig"] and a["iou_tc"] == b["iou_tc"]
                            and a["f5"] == b["f5"] and a["iou_vote"] == b["iou_vote"]
                            and a["topo_freq"] == b["topo_freq"] and np.array_equal(pm[n], pm2[n]))
                    inv &= bool(same)
            row["canon_invariant"] = inv
            # pure reversal (wrong winding, not re-fixed) must be flagged illegal by polycore
            text3 = prod.serialize({n: [(o[::-1], [h[::-1] for h in hs]) for o, hs in parts[n]] for n in names}, names)
            _, res3, _ = eval_text(text3, names, caches, W, H)
            row["reversed_flagged"] = all((not res3[n]["legal"]) or not parts[n] for n in names)
            rows.append(row)
    return rows


def main(nproc=3, limit=0):
    js = items()
    if int(limit):
        js = [j for j in js if j[0] == "riga"][: int(limit)] + [j for j in js if j[0] == "lidc"][: int(limit)]
    t0 = time.time()
    rows = []
    with Pool(int(nproc)) as pool:
        for i, rs in enumerate(pool.imap_unordered(run, js, chunksize=4)):
            rows.extend(rs)
            if i % 100 == 0:
                print(i, "/", len(js), f"{time.time() - t0:.0f}s", flush=True)
    out = prod.HERE / "results"
    out.mkdir(exist_ok=True)
    with open(out / "r1_rows.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r, default=str) + "\n")
    summary = {"targets": len(rows), "wall_s": round(time.time() - t0, 1)}
    for ds in ("riga", "lidc"):
        for mode in MODES:
            R = [r for r in rows if r["ds"] == ds and r["mode"] == mode]
            if not R:
                continue
            names = ["disc", "cup"] if ds == "riga" else ["nodule"]
            objs = [(r, n) for r in R for n in names]
            ious = np.array([r[f"{n}_iou_tc"] for r, n in objs])
            f5 = np.array([r[f"{n}_f5"] for r, n in objs])
            s = {
                "targets": len(R), "objects": len(objs),
                "builder_status": dict(Counter(r["status"][n] for r, n in objs)),
                "roundtrip_exact": float(np.mean([r["roundtrip_exact"] for r in R])),
                "holes_key_always": float(np.mean([r["has_holes_key"] for r in R])),
                "legal": float(np.mean([r[f"{n}_legal"] for r, n in objs])),
                "sig_eq_source": float(np.mean([r[f"{n}_sig_eq"] for r, n in objs])),
                "iou_tc_min": float(ious.min()), "iou_tc_median": float(np.median(ious)),
                "iou_tc_ge_0.99": float(np.mean(ious >= 0.99)), "iou_tc_eq_1": float(np.mean(ious == 1.0)),
                "iou_src_min": float(min(r[f"{n}_iou_src"] for r, n in objs)),
                "f5_min": float(f5.min()), "f5_ge_0.95": float(np.mean(f5 >= 0.95)),
                "raster_eq_fastgeom": float(np.mean([r[f"{n}_raster_eq_fastgeom"] for r, n in objs])),
                "canon_invariant": float(np.mean([r["canon_invariant"] for r in R])),
                "reversed_flagged_illegal": float(np.mean([r["reversed_flagged"] for r in R])),
                "defects_seen": dict(Counter(d for r, n in objs for d in r[f"{n}_defects"])),
                "build_s_median": float(np.median([r["build_s"] for r in R])),
            }
            if ds == "riga":
                s["rim_sig_eq"] = float(np.mean([r["rim_sig_eq"] for r in R]))
                s["rim_iou_tc_min"] = float(min(r["rim_iou_tc"] for r in R))
                s["rim_f5_ge_0.95"] = float(np.mean([r["rim_f5"] >= 0.95 for r in R]))
            iou_req = 1.0 if mode == "lossless" else 0.99
            s["PASS"] = bool(s["legal"] == 1 and s["sig_eq_source"] == 1 and ious.min() >= iou_req
                             and s["f5_ge_0.95"] >= 0.99 and s["raster_eq_fastgeom"] == 1
                             and s["canon_invariant"] == 1 and s["roundtrip_exact"] == 1)
            summary[f"{ds}_{mode}"] = s
            print(ds, mode, json.dumps(s), flush=True)
    (out / "r1_summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:])
