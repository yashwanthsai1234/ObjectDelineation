"""A4 sandbox: cost and correctness of dihedral augmentation with rebuilt certified targets.

For val items of RIGA and LIDC: frame -> per-rater original target; for each of the 7 non-identity dihedral maps
(D4): (a) naive = transform the original target's vertices, rasterize, compare with the transformed rater mask;
(b) rebuilt = build_target on the transformed masks (certified). Also checks the image transform commutes with the
LIDC Pillow upscale (flip/rot of the upscaled image == upscale of the flipped crop).
Run: .venv/bin/python -I aug_cost.py --riga 20 --lidc 200 --out aug_cost.json
"""
import argparse
import json
import time

import numpy as np
from PIL import Image

from objdel.checker import signature
from objdel.data import store
from objdel.geometry import GRID, canonicalize, to_mask
from objdel.targets import FRAMES, FrameConfig, build_target, frame, min_areas

# (name, mask transform on (..., H, W), vertex transform on (n, 2) grid coords (x, y))
OPS = {
    "hflip": (lambda m: m[..., :, ::-1], lambda v: np.c_[GRID - v[:, 0], v[:, 1]]),
    "vflip": (lambda m: m[..., ::-1, :], lambda v: np.c_[v[:, 0], GRID - v[:, 1]]),
    "rot90": (lambda m: np.rot90(m, 1, (-2, -1)), lambda v: np.c_[v[:, 1], GRID - v[:, 0]]),
    "rot180": (lambda m: np.rot90(m, 2, (-2, -1)), lambda v: np.c_[GRID - v[:, 0], GRID - v[:, 1]]),
    "rot270": (lambda m: np.rot90(m, 3, (-2, -1)), lambda v: np.c_[GRID - v[:, 1], v[:, 0]]),
    "transpose": (lambda m: np.swapaxes(m, -1, -2), lambda v: np.c_[v[:, 1], v[:, 0]]),
    "antitranspose": (lambda m: np.rot90(np.swapaxes(m, -1, -2), 2, (-2, -1)),
                      lambda v: np.c_[GRID - v[:, 1], GRID - v[:, 0]]),
}


def iou(a, b):
    u = np.count_nonzero(a | b)
    return 1.0 if u == 0 else np.count_nonzero(a & b) / u


def run(dataset, n_items, every):
    rows, img_checks = [], []
    t_orig = t_aug = 0.0
    n_t = 0
    for i, it in enumerate(store.load(dataset, "val", "data/store")):
        if i % every:
            continue
        if len(rows) and len({r["item"] for r in rows}) >= n_items:
            break
        if it.image is None:
            continue
        fr = frame(it)
        names = list(fr.masks)
        ma = min_areas(fr.masks)
        k = next(iter(fr.masks.values())).shape[0]
        # image commutation (LIDC upscale): crop at frame res, transform, upscale vs upscale, transform
        if FRAMES[dataset].image_side:
            small = frame(it, FrameConfig(image_side=None)).image  # the frame-resolution crop before the upscale
            for op, (fm, _) in OPS.items():
                a = np.ascontiguousarray(fm(np.moveaxis(fr.image, -1, 0)))
                b = np.array(Image.fromarray(np.ascontiguousarray(np.moveaxis(fm(np.moveaxis(small, -1, 0)), 0, -1)))
                             .resize((fr.image.shape[1], fr.image.shape[0]), Image.Resampling.BICUBIC))
                img_checks.append((op, int(np.abs(np.moveaxis(a, 0, -1).astype(int) - b.astype(int)).max()), float(np.abs(np.moveaxis(a, 0, -1).astype(int) - b.astype(int)).mean())))
        for r in range(k):
            masks = {n: fr.masks[n][r] for n in names}
            if not all(m.any() for m in masks.values()) and dataset == "riga":
                continue
            t0 = time.perf_counter()
            try:
                tgt = build_target(masks, names, min_area=ma)
            except ValueError:
                continue
            t_orig += time.perf_counter() - t0
            n_t += 1
            for op, (fm, fv) in OPS.items():
                tm = {n: np.ascontiguousarray(fm(m)) for n, m in masks.items()}
                ma2 = min_areas({n: np.ascontiguousarray(fm(fr.masks[n])) for n in names})
                assert ma2 == ma
                rec = {"item": it.item_id, "rater": r, "op": op, "S": fr.W}
                for n in names:
                    parts = [(fv(o), [fv(h) for h in hs]) for o, hs in tgt.objs[n]]
                    parts = canonicalize({n: parts})[n]  # re-orient (a reflection reverses every ring) and restart
                    g = to_mask(parts, fr.W, fr.H)
                    diff = int(np.count_nonzero(g.mask ^ tm[n]))
                    rec[n] = {"naive_iou": iou(g.mask, tm[n]), "naive_diff_px": diff,
                              "naive_legal": bool(g.legal),
                              "naive_sig_ok": signature(g.mask, ma[n]) == signature(tm[n], ma[n]),
                              "orig_iou": tgt.info[n]["iou"], "orig_vertices": tgt.info[n]["vertices"],
                              "orig_tau": tgt.info[n]["tau"],
                              "naive_meets_cert": bool(g.legal) and iou(g.mask, tm[n]) >= 0.99
                              and signature(g.mask, ma[n]) == signature(tm[n], ma[n])}
                t0 = time.perf_counter()
                try:
                    t2 = build_target(tm, names, min_area=ma)
                    ok = True
                except ValueError:
                    ok = False
                t_aug += time.perf_counter() - t0
                rec["rebuilt_ok"] = ok
                if ok:
                    for n in names:
                        rec[n].update(rebuilt_iou=t2.info[n]["iou"], rebuilt_vertices=t2.info[n]["vertices"],
                                      rebuilt_tau=t2.info[n]["tau"])
                    rec["len_orig"], rec["len_rebuilt"] = len(tgt.text), len(t2.text)
                rows.append(rec)
    return rows, img_checks, t_orig / max(n_t, 1), t_aug / max(len(rows), 1), n_t


def summarize(rows, names):
    out = {"targets_x_ops": len(rows), "rebuilt_ok": float(np.mean([r["rebuilt_ok"] for r in rows]))}
    for n in names:
        rs = [r[n] for r in rows]
        out[n] = {
            "naive_meets_cert": float(np.mean([x["naive_meets_cert"] for x in rs])),
            "naive_exact": float(np.mean([x["naive_diff_px"] == 0 for x in rs])),
            "naive_iou_min": float(np.min([x["naive_iou"] for x in rs])),
            "naive_iou_median": float(np.median([x["naive_iou"] for x in rs])),
            "naive_diff_px_median": float(np.median([x["naive_diff_px"] for x in rs])),
            "naive_sig_ok": float(np.mean([x["naive_sig_ok"] for x in rs])),
            "naive_legal": float(np.mean([x["naive_legal"] for x in rs])),
            "rebuilt_iou_min": float(np.min([x.get("rebuilt_iou", np.nan) for x in rs])),
            "vertices_ratio_median": float(np.median([x["rebuilt_vertices"] / x["orig_vertices"]
                                                      for x in rs if "rebuilt_vertices" in x])),
            "vertices_ratio_p95": float(np.percentile([x["rebuilt_vertices"] / x["orig_vertices"]
                                                       for x in rs if "rebuilt_vertices" in x], 95)),
        }
    lr = [r["len_rebuilt"] / r["len_orig"] for r in rows if "len_rebuilt" in r]
    out["text_len_ratio_median"], out["text_len_ratio_p95"] = float(np.median(lr)), float(np.percentile(lr, 95))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--riga", type=int, default=20)
    ap.add_argument("--lidc", type=int, default=200)
    ap.add_argument("--out", default="aug_cost.json")
    a = ap.parse_args()
    res = {}
    for ds, n, every, names in (("riga", a.riga, 3, ["disc", "cup"]), ("lidc", a.lidc, 5, ["nodule"])):
        t0 = time.perf_counter()
        rows, img, t_o, t_a, n_t = run(ds, n, every)
        res[ds] = {"summary": summarize(rows, names), "targets": n_t,
                   "build_s_per_target_orig": t_o, "build_s_per_target_aug": t_a,
                   "img_commute": ({op: {"maxabs": max(x[1] for x in img if x[0] == op), "meanabs": float(np.mean([x[2] for x in img if x[0] == op]))} for op in OPS} if img else None),
                   "wall_s": time.perf_counter() - t0}
        print(ds, json.dumps(res[ds], indent=1), flush=True)
    json.dump(res, open(a.out, "w"), indent=1)
