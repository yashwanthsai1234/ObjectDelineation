"""A6 prototype H1: fuse score_mask's pixel work (bounding box, area, intersections with the raters / vote, egg-yolk
band depth) into one GIL-free compiled pass, as a stand-in for a C++ port inside polycore. numba njit(nogil=True)
compiles through LLVM like clang, so its speed approximates a C++ kernel; Python threads then run it in parallel.

Exactness: every count is an integer; the band sums add float32 distances >= 1 (frames <= 500 px: < 2^18 terms,
each < 2^10, all multiples of 2^-23), so the float64 sum is exact in any order and equals numpy's pairwise sum
bit for bit. The check below compares every NUMERIC term of every completion with score_completion.

    .venv/bin/python -I -m experiments.phase2.lab.A6.proto_fused --out <dir>
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

import cv2
import numpy as np
from numba import njit

from objdel.checker.topology import signature
from objdel.geometry import parse, to_mask
from objdel.reward import NUMERIC, advantages_v4, rater_answers, score_batch, score_completion
from objdel.reward.score import _FAILED, area_guard, sig_id

from .common import dump, env, stats_ms
from .prof_reward import batches, load_prompts


@njit(nogil=True, cache=True)
def _bbox(m):
    H, W = m.shape
    y0, y1, x0, x1 = H, 0, W, 0
    for y in range(H):
        for x in range(W):
            if m[y, x]:
                if y < y0:
                    y0 = y
                if y >= y1:
                    y1 = y + 1
                if x < x0:
                    x0 = x
                if x >= x1:
                    x1 = x + 1
    return y0, y1, x0, x1


@njit(nogil=True, cache=True)
def fused(m, box, raters, vote, union, core, d_out, d_in):
    """Over the ROI = bbox(m) | box: area, |m & rater_i| for every rater, |m & vote|, and the band depth D."""
    y0, y1, x0, x1 = _bbox(m)
    if box[1] > box[0]:
        if y1 <= y0:
            y0, y1, x0, x1 = box[0], box[1], box[2], box[3]
        else:
            y0, y1, x0, x1 = min(y0, box[0]), max(y1, box[1]), min(x0, box[2]), max(x1, box[3])
    k = raters.shape[0]
    inter = np.zeros(k, np.int64)
    area = 0
    iv = 0
    spill = 0.0
    nspill = 0
    miss = 0.0
    nmiss = 0
    for y in range(y0, y1):
        for x in range(x0, x1):
            p = m[y, x]
            if p:
                area += 1
                if vote[y, x]:
                    iv += 1
                for i in range(k):
                    if raters[i, y, x]:
                        inter[i] += 1
                if not union[y, x]:
                    spill += d_out[y, x]
                    nspill += 1
            elif core[y, x]:
                miss += d_in[y, x]
                nmiss += 1
    d = (spill - 0.5 * nspill) + (miss - 0.5 * nmiss)
    return area, inter, iv, d, y0, y1, x0, x1


@njit(nogil=True, cache=True)
def fused_code(m, code, d_out, d_in, k):
    """One pass over the frame with a per-pixel code (bit i = rater i, bit k = vote, k+1 = union, k+2 = core):
    histogram of the codes under m, band sums, and m's bounding box."""
    H, W = m.shape
    hist = np.zeros(1 << (k + 3), np.int64)
    ub, cb = 1 << (k + 1), 1 << (k + 2)
    spill = 0.0
    nspill = 0
    miss = 0.0
    nmiss = 0
    y0, y1, x0, x1 = H, 0, W, 0
    for y in range(H):
        for x in range(W):
            c = code[y, x]
            if m[y, x]:
                hist[c] += 1
                if not c & ub:
                    spill += d_out[y, x]
                    nspill += 1
                if y < y0:
                    y0 = y
                y1 = y + 1
                if x < x0:
                    x0 = x
                if x >= x1:
                    x1 = x + 1
            elif c & cb:
                miss += d_in[y, x]
                nmiss += 1
    area = hist.sum()
    inter = np.zeros(k + 1, np.int64)
    for c in range(hist.size):
        if hist[c]:
            for i in range(k + 1):
                if c >> i & 1:
                    inter[i] += hist[c]
    return area, inter[:k], inter[k], (spill - 0.5 * nspill) + (miss - 0.5 * nmiss), y0, y1, x0, x1


_CODES = {}


def _code(c):
    x = _CODES.get(id(c))
    if x is None:
        code = np.zeros(c.union.shape, np.uint16)
        for i, r in enumerate(c.raters):
            code |= r.astype(np.uint16) << i
        code |= c.vote.astype(np.uint16) << c.k
        code |= c.union.astype(np.uint16) << (c.k + 1)
        code |= c.core.astype(np.uint16) << (c.k + 2)
        x = _CODES[id(c)] = (c, code)
    return x[1]


MODE = {"kernel": "code"}


def score_mask_fast(m, c):
    if MODE["kernel"] == "code":
        area, inter, iv, d, y0, y1, x0, x1 = fused_code(m, _code(c), c.dist_out_union, c.dist_in_core, c.k)
        if y1 <= y0:
            y0 = y1 = x0 = x1 = 0
    else:
        area, inter, iv, d, y0, y1, x0, x1 = fused(m, np.array(c.box, np.int64), c.raters, c.vote, c.union, c.core,
                                                   c.dist_out_union, c.dist_in_core)
    sig = signature(m[y0:y1, x0:x1], c.min_area)
    same = c.same_sig.get(sig, ())
    if same:
        ious = [inter[i] / u if (u := area + c.rater_areas[i] - inter[i]) else 1.0 for i in same]
        iou_tc = max(ious)
    else:
        u = area + c.n_vote - iv
        iou_tc = iv / u if u else 1.0
    f5 = 1.0 / (1.0 + (2.0 * d / c.perimeter) / (c.t5 * c.t5))
    return iou_tc, f5, area_guard(area, c), sig, len(same) / len(c.voters)


def score_fast(text, truncated, names, caches, W, H):
    answers = rater_answers(caches, names)
    p = parse(text, names, truncated)
    if not p.ok:
        return _FAILED | {"error": p.error, "defects": {}}
    legal, defects, ts = True, {}, []
    for n in names:
        g = to_mask(p.objs[n], W, H)
        legal &= g.legal
        for dd, cnt in g.defects.items():
            defects[dd] = defects.get(dd, 0) + cnt
        ts.append(score_mask_fast(g.mask, caches[n]))
    sigs = tuple(t[3] for t in ts)
    joint = len(answers.get(sigs, ())) / sum(map(len, answers.values()))
    return {"format_ok": 1.0, "legal": float(legal), "iou_tc": float(np.mean([t[0] for t in ts])),
            "iou_min": float(min(t[0] for t in ts)), "f5": float(np.mean([t[1] for t in ts])),
            "g_area": float(np.mean([t[2] for t in ts])), "topo_any": float(joint > 0), "topo_freq": joint,
            "sig_id": sig_id(sigs), "error": "", "defects": defects}


def batch_fast(jobs, pool=None, chunk=8):
    if pool is None:
        rows = [score_fast(*j) for j in jobs]
    else:
        rows = [r for rs in pool.map(lambda js: [score_fast(*j) for j in js],
                                     [jobs[i:i + chunk] for i in range(0, len(jobs), chunk)]) for r in rs]
    return {k: np.array([r[k] for r in rows], np.float64) for k in NUMERIC}


def wall(fn, bl, reps=3):
    out = []
    for jobs in bl:
        for rep in range(reps):
            t0 = time.perf_counter()
            advantages_v4(fn(jobs))
            if rep:
                out.append(1e3 * (time.perf_counter() - t0))
    return stats_ms(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    prompts, _ = load_prompts(64, 200)
    res = {"env": env()}
    for label, share in (("riga", 1.0), ("mixed", 0.4)):
        bl = batches(prompts, share, 3, seed=11)
        mism = 0
        for jobs in bl:
            ref = [score_completion(*j) for j in jobs]
            new = [score_fast(*j) for j in jobs]
            mism += sum(r[k] != x[k] for r, x in zip(ref, new, strict=True) for k in NUMERIC)
        r = {"mismatched_terms": int(mism)}
        r["library_serial"] = wall(lambda j: score_batch(j), bl)
        for kern in ("roi", "code"):
            MODE["kernel"] = kern
            mism = sum(score_completion(*j)[k] != score_fast(*j)[k] for jobs in bl for j in jobs for k in NUMERIC)
            r[f"{kern}_mismatched_terms"] = int(mism)
            r[f"{kern}_serial"] = wall(lambda j: batch_fast(j), bl)
            for nt in (2, 4, 8):
                with ThreadPoolExecutor(nt) as ex:
                    r[f"{kern}_threads{nt}"] = wall(lambda j, ex=ex: batch_fast(j, ex), bl)
        for nt in (2, 4, 8):
            with ThreadPoolExecutor(nt) as ex:
                r[f"library_threads{nt}"] = wall(lambda j, ex=ex: score_batch(j, ex), bl)
        res[label] = r
        print(label, {k: (round(v["median"], 1) if isinstance(v, dict) else v) for k, v in r.items()}, flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/proto_fused.json")


if __name__ == "__main__":
    main()
