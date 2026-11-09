"""Offline dihedral (D4) augmentation of every train target, rebuilt and certified from transformed masks
(H-A6-6, H-A4-5; issue I7: a flip must rebuild targets from flipped masks). Pass criteria: prereg.py (P3).

    .venv/bin/python -I -m experiments.phase2.lab.B3.dihedral --out <scratch>/d4_run1 --procs 2 --checks
    .venv/bin/python -I -m experiments.phase2.lab.B3.dihedral --out <scratch>/d4_run2 --procs 1

Per train item of LIDC and RIGA (data/store, the same items and raters as data/targets): objdel.targets.frame, the
item's size rule min_areas (D4-invariant: union areas do not change), then for each of the 8 ops and each rater
build_target(op(masks), names, BuildConfig(), min_areas). Writes <out>/<ds>/<op>.jsonl (item_id, rater, text) and
<out>/result.json with the per-variant target_set_hash and the build CPU time (Target.info["time"], process time).
With --checks, also (P3.1-2, P3.6):
- text-only re-check: parse -> polycore raster (legal) -> size-filtered signature == op(mask)'s at min_areas,
  IoU >= 0.99; RIGA rim (disc & ~cup) signature too, as experiments/phase1/build_targets;
- inverse op of the variant raster vs the original rater mask: IoU, signature, pixel identity with the identity
  variant's raster;
- F5 and IoU_tc of the inverse-transformed raster against the original item's stored reward cache (data/prepared),
  for voting raters of prepared items (R1 as in Phase 1). D4 maps the pixel lattice onto itself and the exact EDT,
  boundaries, connectivity and IoU are D4-invariant, so this equals scoring the variant against a transformed cache;
- answer tokens per tokenizer (HF_HUB_OFFLINE=1; add_special_tokens=False as build_targets);
- PNG bytes of the frame-resolution model crop (LIDC: before the 512 px upscale, A4 S1) for the identity op, and
  for all 8 ops on every 10th item (disk estimate).
Reads untrusted data: run with python -I.
"""
import argparse
import hashlib
import io
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image

from objdel.checker import signature
from objdel.data import store
from objdel.geometry import SCHEMAS, parse, to_mask
from objdel.reward import ExpertCache
from objdel.reward.score import score_mask
from objdel.targets import BuildConfig, FrameConfig, build_target, frame, min_areas, target_set_hash

from . import prereg

OPS = {  # op on (..., H, W) arrays, and its inverse's name
    "identity": (lambda m: m, "identity"),
    "hflip": (lambda m: m[..., :, ::-1], "hflip"),
    "vflip": (lambda m: m[..., ::-1, :], "vflip"),
    "rot90": (lambda m: np.rot90(m, 1, (-2, -1)), "rot270"),
    "rot180": (lambda m: np.rot90(m, 2, (-2, -1)), "rot180"),
    "rot270": (lambda m: np.rot90(m, 3, (-2, -1)), "rot90"),
    "transpose": (lambda m: np.swapaxes(m, -1, -2), "transpose"),
    "antitranspose": (lambda m: np.rot90(np.swapaxes(m, -1, -2), 2, (-2, -1)), "antitranspose"),
}
TOKENIZERS = {"smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct", "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
              "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct", "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B"}
_TOK = {}
PREPARED = Path("data/prepared")


def _op(name, a):
    return np.ascontiguousarray(OPS[name][0](a))


def _tokens(text):
    if not _TOK:
        from transformers import AutoTokenizer
        _TOK.update({k: AutoTokenizer.from_pretrained(v) for k, v in TOKENIZERS.items()})
    return {k: len(t(text, add_special_tokens=False)["input_ids"]) for k, t in _TOK.items()}


def _png_bytes(img):
    a = img[..., 0] if img.ndim == 3 and (img[..., 0] == img[..., 2]).all() and (img[..., 0] == img[..., 1]).all() \
        else img
    b = io.BytesIO()
    Image.fromarray(a).save(b, format="PNG", optimize=False, compress_level=6)
    return b.tell()


def _iou(a, b):
    u = np.count_nonzero(a | b)
    return 1.0 if u == 0 else np.count_nonzero(a & b) / u


def _job(job):
    item_id, ds, raters, masks, image, checks, idx, rec = job
    names = SCHEMAS[ds]
    ma = min_areas(masks)
    out, ref_raster = [], {}
    caches = None
    if checks and rec is not None:
        caches = {n: ExpertCache.load(PREPARED / "caches" / item_id / f"{n}.npz", d) for n, d in rec["caches"].items()}
    for op in OPS:
        tm = {n: _op(op, masks[n]) for n in names}
        assert min_areas(tm) == ma
        for r, rater in enumerate(raters):
            src = {n: tm[n][r] for n in names}
            row = {"item_id": item_id, "rater": rater, "op": op}
            try:
                t = build_target(src, names, BuildConfig(), ma)
            except ValueError as e:
                row["error"] = str(e)
                out.append(row)
                continue
            row.update(text=t.text, cpu=t.info["time"])
            if checks:
                t0 = time.process_time()
                p = parse(t.text, names)
                ok, inv_ok, pix_same, f5s, ious_tc = p.ok, True, True, [], []
                ras = {}
                for n in names:
                    g = to_mask(p.objs[n], src[n].shape[1], src[n].shape[0])
                    ras[n] = g.mask
                    ok &= g.legal and _iou(g.mask, src[n]) >= prereg.P3_IOU and \
                        signature(g.mask, ma[n]) == signature(src[n], ma[n])
                    back = _op(OPS[op][1], g.mask)
                    orig = masks[n][r]
                    inv_ok &= _iou(back, orig) >= prereg.P3_IOU and signature(back, ma[n]) == signature(orig, ma[n])
                    if op == "identity":
                        ref_raster[(r, n)] = g.mask
                    pix_same &= bool((back == ref_raster[(r, n)]).all())
                    if caches is not None and r in caches[n].voters:
                        terms = score_mask(back, caches[n])
                        f5s.append(terms.f5)
                        ious_tc.append(terms.iou_tc)
                if "disc" in names:
                    rim, rim_src = ras["disc"] & ~ras["cup"], src["disc"] & ~src["cup"]
                    ok &= signature(rim, ma["rim"]) == signature(rim_src, ma["rim"])
                row.update(ok=bool(ok), inv_ok=bool(inv_ok), pixel_same_as_identity=bool(pix_same),
                           f5=min(f5s) if f5s else None, iou_tc=min(ious_tc) if ious_tc else None,
                           tokens=_tokens(t.text), check_cpu=time.process_time() - t0)
            out.append(row)
    png = None
    if checks and image is not None:
        png = {"identity": _png_bytes(image)}
        if idx % 10 == 0:
            png.update({op: _png_bytes(_op_img(op, image)) for op in OPS if op != "identity"})
    return out, png


def _op_img(op, img):
    return np.ascontiguousarray(np.moveaxis(OPS[op][0](np.moveaxis(img, -1, 0)), 0, -1))


def _jobs(root, ds, checks, limit=None):
    recs = {}
    if checks:
        with open(PREPARED / ds / "train.jsonl") as fh:
            for line in fh:
                r = json.loads(line)
                recs[r["item_id"]] = r
    for idx, it in enumerate(store.load(ds, "train", root)):
        if limit is not None and idx >= limit:
            return
        f = frame(it)
        img = None
        if checks and it.image is not None:
            img = frame(it, FrameConfig(image_side=None)).image  # frame-resolution crop (before any upscale)
        yield it.item_id, ds, it.raters, f.masks, img, checks, idx, recs.get(it.item_id)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--checks", action="store_true")
    ap.add_argument("--datasets", default="lidc,riga")
    ap.add_argument("--limit", type=int, default=None, help="first n items per dataset (smoke test)")
    a = ap.parse_args()
    if a.checks:
        os.environ.setdefault("HF_HUB_OFFLINE", "1")
    out = Path(a.out)
    res = {"procs": a.procs, "checks": a.checks, "datasets": {}}
    t_all = time.time()
    with Pool(a.procs) as pool:
        for ds in a.datasets.split(","):
            t0 = time.time()
            rows, pngs = [], []
            for rs, png in pool.imap(_job, _jobs(a.store, ds, a.checks, a.limit), chunksize=4):
                rows += rs
                if png:
                    pngs.append(png)
            (out / ds).mkdir(parents=True, exist_ok=True)
            d = {"wall_s": round(time.time() - t0, 1), "targets_per_op": {}, "hash": {}, "errors": [],
                 "build_cpu_s": round(sum(r.get("cpu", 0) for r in rows), 1)}
            for op in OPS:
                rs = [r for r in rows if r["op"] == op]
                d["errors"] += [(r["item_id"], r["rater"], op, r["error"]) for r in rs if "error" in r]
                ok = [r for r in rs if "text" in r]
                with open(out / ds / f"{op}.jsonl", "w") as fh:
                    for r in ok:
                        fh.write(json.dumps({"item_id": r["item_id"], "rater": r["rater"], "text": r["text"]}) + "\n")
                d["targets_per_op"][op] = len(ok)
                d["hash"][op] = target_set_hash((r["item_id"], r["rater"], r["text"]) for r in ok)
                d.setdefault("text_bytes", {})[op] = (out / ds / f"{op}.jsonl").stat().st_size
            if a.checks:
                d["checks"] = summarize_checks(rows)
                d["png"] = summarize_png(pngs)
            d["identity_matches_data_targets"] = d["hash"]["identity"] == prereg.IDENTITY_HASH[f"{ds}/train"]
            res["datasets"][ds] = d
            print(ds, json.dumps({k: v for k, v in d.items() if k != "checks"}), flush=True)
            if a.checks:
                print(ds, json.dumps(d["checks"]), flush=True)
    res["wall_s"] = round(time.time() - t_all, 1)
    res["build_cpu_min"] = round(sum(d["build_cpu_s"] for d in res["datasets"].values()) / 60, 2)
    res["all_hash"] = hashlib.sha256(json.dumps({ds: d["hash"] for ds, d in res["datasets"].items()},
                                                sort_keys=True).encode()).hexdigest()
    with open(out / "result.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "datasets"}))


def summarize_checks(rows):
    ok = [r for r in rows if "text" in r]
    out = {"n": len(ok), "ok": float(np.mean([r["ok"] for r in ok])), "inv_ok": float(np.mean([r["inv_ok"] for r in ok])),
           "check_cpu_s": round(sum(r["check_cpu"] for r in ok), 1), "per_op": {}}
    f5 = [r["f5"] for r in ok if r["f5"] is not None]
    itc = [r["iou_tc"] for r in ok if r["iou_tc"] is not None]
    out.update(f5_n=len(f5), f5_min=float(min(f5)) if f5 else None,
               f5_ge=float(np.mean(np.array(f5) >= prereg.P3_F5)) if f5 else None,
               iou_tc_min=float(min(itc)) if itc else None,
               iou_tc_ge=float(np.mean(np.array(itc) >= prereg.P3_IOU)) if itc else None)
    ident = {(r["item_id"], r["rater"]): r for r in ok if r["op"] == "identity"}
    for op in OPS:
        rs = [r for r in ok if r["op"] == op]
        e = {"n": len(rs), "ok": float(np.mean([r["ok"] for r in rs])),
             "pixel_same_as_identity": float(np.mean([r["pixel_same_as_identity"] for r in rs])), "tokens": {}}
        for k in TOKENIZERS:
            rel = np.array([r["tokens"][k] / ident[(r["item_id"], r["rater"])]["tokens"][k] - 1 for r in rs])
            tok = np.array([r["tokens"][k] for r in rs])
            idt = np.array([ident[(r["item_id"], r["rater"])]["tokens"][k] for r in rs])
            budget = {"smolvlm-256m": 1536, "lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536, "qwen3.5-0.8b": 1536}[k]
            e["tokens"][k] = {"median": float(np.median(tok)), "p99": float(np.percentile(tok, 99)),
                              "max": int(tok.max()), "rel_change_median": float(np.median(rel)),
                              "rel_change_p5_p95": [float(np.percentile(rel, 5)), float(np.percentile(rel, 95))],
                              "over_budget": int((tok > budget).sum()),
                              "over_budget_identity_ok": int(((tok > budget) & (idt <= budget)).sum())}
        out["per_op"][op] = e
    return out


def summarize_png(pngs):
    ident = np.array([p["identity"] for p in pngs])
    full = [p for p in pngs if len(p) == 8]
    ratio = {op: float(np.median([p[op] / p["identity"] for p in full])) for op in OPS} if full else {}
    return {"items": len(pngs), "identity_total_mb": round(ident.sum() / 2 ** 20, 1),
            "identity_median_bytes": float(np.median(ident)), "op_size_ratio_median": ratio,
            "all8_estimate_mb": round(ident.sum() * sum(ratio.values()) / 2 ** 20, 1) if ratio else None}


if __name__ == "__main__":
    main()
