"""K3/K4 benchmark: score 256 realistic rollouts (rasterize eo+nz, IoU counts vs GT, full
validity) with polycore (inline, internal pool, Python ThreadPoolExecutor), the numba
prototype, and the current pure-Python repo path (PIL + O(n^2) ring_problems).

    python bench/bench_k3k4.py --out results/bench.json [--skip-repo]
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "tests"))

import numpy as np  # noqa: E402

import polycore as pc  # noqa: E402
from rollouts import blob, make_rollouts  # noqa: E402


def pct(xs, q):
    return float(np.percentile(np.asarray(xs), q))


def lat_stats(ts):
    ts = [t * 1e3 for t in ts]
    return dict(p50_ms=round(pct(ts, 50), 4), p99_ms=round(pct(ts, 99), 4), max_ms=round(max(ts), 4),
                mean_ms=round(statistics.fmean(ts), 4))


def per_rollout(fn, rolls, reps):
    """Median-of-reps wall and thread-CPU time per rollout (single thread)."""
    wall, cpu = [], []
    for r in rolls:
        w, c = [], []
        for _ in range(reps):
            t0, c0 = time.perf_counter(), time.thread_time()
            fn(r)
            w.append(time.perf_counter() - t0)
            c.append(time.thread_time() - c0)
        wall.append(statistics.median(w))
        cpu.append(statistics.median(c))
    return dict(wall=lat_stats(wall), thread_cpu=lat_stats(cpu))


def batch_wall(fn, reps):
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return dict(median_ms=round(statistics.median(ts) * 1e3, 3), best_ms=round(min(ts) * 1e3, 3),
                p90_ms=round(pct(ts, 90) * 1e3, 3))


def loadavg():
    try:
        return [round(x, 1) for x in os.getloadavg()]
    except OSError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--reps", type=int, default=15)
    ap.add_argument("--skip-repo", action="store_true")
    ap.add_argument("--skip-numba", action="store_true")
    a = ap.parse_args()
    rolls = make_rollouts(256, seed=0)
    R = {"machine": dict(platform=platform.platform(), processor=platform.processor(), cpus=os.cpu_count(),
                         python=sys.version.split()[0]),
         "load_before": loadavg(),
         "workload": dict(rollouts=len(rolls), rings=[len(r["rings"]) for r in rolls].count(1),
                          vertices_total=int(sum(r["xy"].shape[0] for r in rolls)),
                          S_min=min(r["S"] for r in rolls), S_max=max(r["S"] for r in rolls))}
    R["workload"]["rings_hist"] = {k: [len(r["rings"]) for r in rolls].count(k) for k in (1, 2, 3)}
    del R["workload"]["rings"]

    args = [(r["xy"], r["off"], r["S"], r["S"], r["roles"], r["gt"]) for r in rolls]
    xys, offs, Ws, Hs, roles, gts = map(list, zip(*args))

    # ---------- polycore
    score = pc._score
    one = lambda r: score(r["xy"], r["off"], r["S"], r["S"], r["roles"], r["gt"], 1000)  # noqa: E731
    for r in rolls[:8]:
        one(r)
    R["polycore_latency_single_thread"] = per_rollout(one, rolls, a.reps)
    R["polycore_raster_only"] = per_rollout(lambda r: pc._rasterize(r["xy"], r["off"], r["S"], r["S"], 1000, False), rolls, a.reps)
    R["polycore_validate_only"] = per_rollout(lambda r: pc._validate(r["xy"], r["off"], r["roles"], 1000, "sweep"), rolls, a.reps)
    R["polycore_validate_only_brute"] = per_rollout(lambda r: pc._validate(r["xy"], r["off"], r["roles"], 1000, "brute"), rolls, a.reps)
    sb = pc._score_batch
    R["polycore_batch"] = {}
    for t in (1, 2, 4, 8):
        R["polycore_batch"][f"internal_pool_{t}"] = batch_wall(lambda: sb(xys, offs, Ws, Hs, roles, gts, 1000, t), a.reps * 2)
    for t in (1, 8):
        ex = ThreadPoolExecutor(t)
        list(ex.map(one, rolls))
        R["polycore_batch"][f"py_threadpool_{t}"] = batch_wall(lambda: list(ex.map(one, rolls)), a.reps * 2)
        ex.shutdown()

    # ---------- worst case: 3 rings x 200 vertices at S = 512
    rng = np.random.default_rng(1)
    worst = []
    for _ in range(64):
        rings = [blob(rng, 200, 500, 500, 330), blob(rng, 200, 500, 500, 80, ccw=True), blob(rng, 200, 200, 200, 90)]
        xy, off = pc.pack(rings)
        worst.append(dict(xy=xy, off=off, roles=np.array([-1, 0, -1], np.int32), S=512,
                          gt=np.zeros((512, 512), np.uint8)))
    R["polycore_worst_3x200_S512"] = per_rollout(one, worst, a.reps)

    # ---------- validity scaling: sweep vs brute (one blob ring)
    R["validate_scaling_us"] = {}
    for n in (16, 64, 200, 500, 1000, 2000):
        xy, off = pc.pack([blob(rng, n, 500, 500, 300)])
        row = {}
        for m in ("sweep", "brute"):
            pc._validate(xy, off, None, 1000, m)
            ts = []
            for _ in range(max(5, a.reps)):
                t0 = time.perf_counter()
                pc._validate(xy, off, None, 1000, m)
                ts.append(time.perf_counter() - t0)
            row[m] = round(statistics.median(ts) * 1e6, 2)
        R["validate_scaling_us"][str(n)] = row

    # ---------- numba prototype
    if not a.skip_numba:
        import numba_proto as nbp
        prepared = [(np.asarray(r["xy"], np.int64), np.asarray(r["off"], np.int64), r["S"], r["gt"].astype(np.bool_))
                    for r in rolls]
        f = lambda p: nbp.reward(p[0], p[1], p[2], p[3])  # noqa: E731
        t0 = time.perf_counter()
        f(prepared[0])
        R["numba_jit_warmup_s"] = round(time.perf_counter() - t0, 2)
        R["numba_latency_single_thread"] = per_rollout(f, prepared, a.reps)
        R["numba_batch"] = {}
        for t in (1, 8):
            ex = ThreadPoolExecutor(t)
            list(ex.map(f, prepared))
            R["numba_batch"][f"py_threadpool_{t}"] = batch_wall(lambda: list(ex.map(f, prepared)), a.reps * 2)
            ex.shutdown()

    # ---------- current repo pure-Python path
    if not a.skip_repo:
        from objdel.checker.polygon_io import rasterize as pil_rasterize
        from objdel.checker.polygon_io import ring_problems

        def repo(r):
            S = r["S"]
            m = np.zeros((S, S), bool)
            labels = []
            for ring in r["rings"]:
                pts = np.asarray(ring, float)
                labels.append(ring_problems(pts, bounds=(0, 1000)))
                m ^= pil_rasterize(pts * S / 1000.0, S, S)
            g = r["gt"].astype(bool)
            return (m & g).sum() / max((m | g).sum(), 1), labels

        R["repo_latency_single_thread"] = per_rollout(repo, rolls, 1)
        R["repo_batch"] = {"serial": batch_wall(lambda: [repo(r) for r in rolls], 1)}
    R["load_after"] = loadavg()
    print(json.dumps(R, indent=1))
    if a.out:
        Path(a.out).write_text(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
