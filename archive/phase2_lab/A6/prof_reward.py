"""A6 (d): reward scoring of 256 rollouts (32 prompts x 8) on the prepared caches: wall time serial and threaded,
per-stage breakdown, cProfile hot spots, cache-load cost, and the C++ floor (polycore.score_batch on the same shapes).

Workload: experiments.phase1.reward_gates.make_jobs (the K3 generator: rater masks jittered / shifted / grown /
exact, traced with findContours + approxPolyDP 1 px), but the caches come from data/prepared (ExpertCache.load,
the production path) instead of the fixtures.

    .venv/bin/python -I -m experiments.phase2.lab.A6.prof_reward --out <dir>
"""
import argparse
import time
from concurrent.futures import ThreadPoolExecutor

import cv2
import numpy as np

from experiments.phase1.reward_gates import make_jobs
from objdel.checker import polycore
from objdel.checker.topology import instance_topology, signature, size_filter
from objdel.geometry import GRID, parse, to_mask
from objdel.reward import advantages_v4, rater_answers, score_batch, score_completion
from objdel.reward.cache import band_depth

from .common import Timer, by_module, caches, dump, env, profile, records, stats_ms


def load_prompts(n_riga: int, n_lidc: int):
    t = Timer()
    out = {}
    for ds, n in (("riga", n_riga), ("lidc", n_lidc)):
        recs = [r for r in records(ds, "train") if r["rl"]][:n]
        ps = []
        for r in recs:
            t0 = time.perf_counter()
            c = caches(r)
            t.add(f"load_{ds}", time.perf_counter() - t0)
            names = tuple(c)
            ps.append((("disc", "cup") if ds == "riga" else ("nodule",), c))
            assert set(names) == set(ps[-1][0])
        out[ds] = ps
    return out, {k: t.t[k] / t.n[k] for k in t.t}


def batches(prompts, share_riga: float, n: int, seed: int):
    rng = np.random.default_rng(seed)
    out = []
    for b in range(n):
        n_r = int(round(32 * share_riga))
        pr = [prompts["riga"][i] for i in rng.choice(len(prompts["riga"]), n_r, replace=False)]
        pl = [prompts["lidc"][i] for i in rng.choice(len(prompts["lidc"]), 32 - n_r, replace=False)] if n_r < 32 else []
        out.append(make_jobs(pr + pl, 1000 + b))
    return out


def staged(job, t: Timer):
    """score_objects step by step (same calls), timing each stage; returns the library's result for a check."""
    text, truncated, names, cs, W, H = job
    t0 = time.perf_counter()
    rater_answers(cs, names)
    t1 = time.perf_counter()
    t.add("rater_answers", t1 - t0)
    p = parse(text, names, truncated)
    t2 = time.perf_counter()
    t.add("parse", t2 - t1)
    if not p.ok:
        return
    for n in names:
        rings, roles = [], []
        for outer, holes in p.objs[n]:
            roles += [-1] + [len(rings)] * len(holes)
            rings += [outer, *holes]
        a = time.perf_counter()
        shape = polycore.pack(rings)
        b = time.perf_counter()
        polycore.validate(shape, roles=np.array(roles, np.int32), G=GRID)
        c = time.perf_counter()
        m = np.asarray(polycore.rasterize(shape, W, H, G=GRID)["evenodd"], bool)
        d = time.perf_counter()
        t.add("polycore.pack", b - a)
        t.add("polycore.validate", c - b)
        t.add("polycore.rasterize", d - c)
        cache = cs[n]
        # score_mask, stage by stage
        e = time.perf_counter()
        x, y, w, h = cv2.boundingRect(m.view(np.uint8))
        boxes = [bb for bb in (cache.box, (y, y + h, x, x + w)) if bb[1] > bb[0]]
        roi = (slice(min(bb[0] for bb in boxes), max(bb[1] for bb in boxes)),
               slice(min(bb[2] for bb in boxes), max(bb[3] for bb in boxes)))
        mm = m[roi]
        area = int(np.count_nonzero(mm))
        f = time.perf_counter()
        sig = signature(mm, cache.min_area)
        g = time.perf_counter()
        same = cache.same_sig.get(sig, ())
        if same:
            for i in same:
                inter = int(np.count_nonzero(mm & cache.raters[i][roi]))
                _ = inter / max(1, area + cache.rater_areas[i] - inter)
        else:
            inter = int(np.count_nonzero(mm & cache.vote[roi]))
        hh = time.perf_counter()
        band_depth(mm, cache.union[roi], cache.core[roi], cache.dist_out_union[roi], cache.dist_in_core[roi])
        i_ = time.perf_counter()
        t.add("score_mask.roi+area", f - e)
        t.add("score_mask.signature", g - f)
        t.add("score_mask.iou", hh - g)
        t.add("score_mask.band_depth", i_ - hh)


def signature_parts(jobs, t: Timer):
    """Split signature(): instance_topology (2 CCL) vs size_filter + second instance_topology."""
    for text, truncated, names, cs, W, H in jobs:
        p = parse(text, names, truncated)
        if not p.ok:
            continue
        for n in names:
            m = to_mask(p.objs[n], W, H).mask
            a = time.perf_counter()
            tp = instance_topology(m)
            b = time.perf_counter()
            t.add("instance_topology", b - a)
            ma = cs[n].min_area
            if (tp.areas < ma).any() or (tp.hole_areas < ma).any():
                c = time.perf_counter()
                instance_topology(size_filter(m, ma))
                t.add("size_filter+instance_topology (slow path)", time.perf_counter() - c)


def wall(batch_list, pool, reps=2):
    walls = []
    for jobs in batch_list:
        for rep in range(reps):
            t0 = time.perf_counter()
            advantages_v4(score_batch(jobs, pool))
            if rep:
                walls.append(1e3 * (time.perf_counter() - t0))
    return stats_ms(walls)


def cxx_floor(jobs, threads):
    """polycore.score_batch (raster stats + validity + intersection with the vote) on the jobs' shapes: what the
    C++ part of the reward costs when it is called once per batch."""
    shapes, Ws, Hs, roles, gts = [], [], [], [], []
    for text, truncated, names, cs, W, H in jobs:
        p = parse(text, names, truncated)
        if not p.ok:
            continue
        for n in names:
            rings, rl = [], []
            for outer, holes in p.objs[n]:
                rl += [-1] + [len(rings)] * len(holes)
                rings += [outer, *holes]
            shapes.append(polycore.pack(rings))
            Ws.append(W)
            Hs.append(H)
            roles.append(np.array(rl, np.int32))
            gts.append(cs[n].vote)
    out = []
    for _ in range(5):
        t0 = time.perf_counter()
        polycore.score_batch(shapes, Ws, Hs, roles=roles, gts=gts, threads=threads)
        out.append(1e3 * (time.perf_counter() - t0))
    return float(np.median(out)), len(shapes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--batches", type=int, default=4)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    res = {"env": env()}
    prompts, load_ms = load_prompts(64, 200)
    res["cache_load_ms_per_item"] = load_ms
    print("cache load ms/item", load_ms, flush=True)
    for label, share in (("mixed", 0.4), ("riga", 1.0), ("lidc", 0.0)):
        bl = batches(prompts, share, a.batches, seed=7)
        r = {"chars_per_completion": float(np.mean([len(j[0]) for b in bl for j in b]))}
        lat = []
        for jobs in bl:
            for j in jobs:
                t0 = time.perf_counter()
                score_completion(*j)
                lat.append(1e3 * (time.perf_counter() - t0))
        r["per_completion_ms"] = stats_ms(lat)
        r["serial"] = wall(bl, None)
        for nt in (2, 4, 8):
            with ThreadPoolExecutor(nt) as ex:
                r[f"threads{nt}"] = wall(bl, ex)
        t = Timer()
        for jobs in bl:
            for j in jobs:
                staged(j, t)
        nb = len(bl)
        r["stages_ms_per_batch"] = {k: v / nb for k, v in sorted(t.t.items(), key=lambda kv: -kv[1])}
        t2 = Timer()
        for jobs in bl:
            signature_parts(jobs, t2)
        r["signature_split_ms_per_batch"] = {k: v / nb for k, v in t2.t.items()}
        r["signature_slow_path_share"] = t2.n["size_filter+instance_topology (slow path)"] / max(1, t2.n["instance_topology"])
        _, rows, total = profile(lambda: [score_batch(j) for j in bl], top=200)
        r["cprofile_total_ms_per_batch"] = total / nb
        r["cprofile_top"] = [(round(tt / nb, 3), round(ct / nb, 3), nc // nb, name) for tt, ct, nc, name in rows[:22]]
        r["cprofile_by_file_ms_per_batch"] = {k: round(v / nb, 3) for k, v in list(by_module(rows).items())[:15]}
        r["cxx_floor_ms"] = {f"threads{th}": cxx_floor(bl[0], th) for th in (1, 8)}
        res[label] = r
        print(label, "serial", r["serial"]["median"], "t8", r["threads8"]["median"], flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/reward.json")


if __name__ == "__main__":
    main()
