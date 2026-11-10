"""A6 micro-benchmark for H1/H2: per-object cost of score_mask's parts, interleaved A/B per object so machine load
drifts cancel (the Mac is shared with other agents). Objects: the masks of the RIGA / LIDC K3-style batches.

Variants per object (same mask, same cache):
- lib: objdel.reward.score_mask (numpy + OpenCV on the ROI);
- roi_kernel: cv2.boundingRect + the per-pixel-code kernel (proto_fused.fused_code, numba nogil) on the ROI only;
- signature: objdel.checker.topology.signature on the ROI (2 OpenCV CCL; the part a C++ port would also own);
- parse, to_mask: the per-completion front end (Python json/regex; polycore through its Python wrapper).

    .venv/bin/python -I -m experiments.phase2.lab.A6.micro_score_mask --out <dir>
"""
import argparse
import time

import cv2
import numpy as np

from objdel.checker.topology import signature
from objdel.geometry import parse, to_mask
from objdel.reward import score_mask

from .common import dump, env
from .prof_reward import batches, load_prompts
from .proto_fused import _code, fused_code


def roi_kernel(m, c):
    x, y, w, h = cv2.boundingRect(m.view(np.uint8))
    b = c.box
    if w:
        y0, y1, x0, x1 = (min(y, b[0]), max(y + h, b[1]), min(x, b[2]), max(x + w, b[3])) if b[1] > b[0] else \
            (y, y + h, x, x + w)
    else:
        y0, y1, x0, x1 = b
    code = _code(c)
    sl = (slice(y0, y1), slice(x0, x1))
    out = fused_code(m[sl], code[sl], c.dist_out_union[sl], c.dist_in_core[sl], c.k)
    return out, sl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    prompts, _ = load_prompts(64, 200)
    res = {"env": env()}
    for label, share in (("riga", 1.0), ("lidc", 0.0)):
        jobs = batches(prompts, share, 1, seed=5)[0]
        objs, texts = [], []
        for text, tr, names, cs, W, H in jobs:
            p = parse(text, names, tr)
            texts.append((text, names, W, H))
            if p.ok:
                for n in names:
                    objs.append((to_mask(p.objs[n], W, H).mask, cs[n], p.objs[n], W, H))
        for m, c, *_ in objs:  # warm up (JIT, code images)
            roi_kernel(m, c)
        T = {k: [] for k in ("lib", "roi_kernel", "signature", "to_mask", "lib_minus_parts")}
        mism = 0
        for rep in range(5):
            for m, c, parts, W, H in objs:
                t0 = time.perf_counter_ns()
                ref = score_mask(m, c)
                t1 = time.perf_counter_ns()
                (area, inter, iv, d, *_), sl = roi_kernel(m, c)
                t2 = time.perf_counter_ns()
                sig = signature(m[sl], c.min_area)
                t3 = time.perf_counter_ns()
                to_mask(parts, W, H)
                t4 = time.perf_counter_ns()
                if rep == 0:
                    f5 = 1.0 / (1.0 + (2.0 * d / c.perimeter) / (c.t5 * c.t5))
                    mism += (f5 != ref.f5) + (area != ref.area) + (sig != ref.sig)
                T["lib"].append((t1 - t0) / 1e3)
                T["roi_kernel"].append((t2 - t1) / 1e3)
                T["signature"].append((t3 - t2) / 1e3)
                T["to_mask"].append((t4 - t3) / 1e3)
        P = []
        for _ in range(5):
            for text, names, _W, _H in texts:
                t0 = time.perf_counter_ns()
                parse(text, names)
                P.append((time.perf_counter_ns() - t0) / 1e3)
        med = {k: float(np.median(v)) for k, v in T.items() if v}
        res[label] = {"objects": len(objs), "completions": len(texts), "mismatches": int(mism),
                      "median_us_per_object": med, "mean_us_per_object": {k: float(np.mean(v)) for k, v in T.items() if v},
                      "parse_median_us": float(np.median(P)),
                      "batch_estimate_ms": {
                          "lib score_mask": float(np.mean(T["lib"])) * len(objs) / 5 / 1e3 * 5,
                          "kernel+signature": (float(np.mean(T["roi_kernel"])) + float(np.mean(T["signature"])))
                          * len(objs) / 1e3,
                          "to_mask": float(np.mean(T["to_mask"])) * len(objs) / 1e3,
                          "parse": float(np.mean(P)) * len(texts) / 1e3}}
        print(label, res[label], flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/micro_score_mask.json")


if __name__ == "__main__":
    main()
