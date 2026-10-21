"""K3/K4 reward-speed gate on a training box (plan §4.7 row 6, D14), with no RIGA data on the host.

RIGA-like prompts are synthetic (below; caches built at start, ~0.3 s per item); LIDC prompts are the first --lidc
RL train items of a prepared root, after the production preload of every LIDC cache (serve.preload, timed). Three
workloads (RIGA-only, mixed 40 % RIGA, LIDC-only), --batches batches of 32 prompts x 8 completions
(reward_gates.make_jobs). Per workload: the warm wall (second of two runs of advantages_v4(score_batch)) serially,
on --threads threads and on a fork pool of each --workers size, serial re-timed per batch beside the pool (paired
speed-up); per-rollout latency (serial); every batch's hash (reward_gates.digest), which every executor must
reproduce. Gate: K3 = best of serial / pool median <= 250 ms on RIGA-only and mixed; K4 = per-rollout p99 <= 5 ms;
hashes equal; D14's C++ trigger = pool-4 RIGA-only median > 125 ms or paired speed-up < 3.

    python -I -m experiments.phase2.k3_harness --out k3.json [--root data/prepared --workers 2,4,8 --batches 5]

Synthetic items (lab B4 synth.py, H-K3SYN): 500 x 500 frames, 6 raters' disc and cup, smooth outlines
r(t) = R (1 + sum of low harmonics) with per-rater perturbations, the cup inside the disc, parameters from the
5-95 % quantiles of aggregate RIGA statistics only (no RIGA pixel, mask or id). Python floats and math.cos/sin feed
the 1/16 px vertices, so the corpus is identical on arm64 and x86-64 (digest of corpus(64, 0): e31d8f1e...).
"""

import argparse
import hashlib
import json
import math
import os
import platform
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import cv2
import numpy as np

from experiments.phase1.reward_gates import digest, make_jobs
from objdel.geometry import SCHEMAS
from objdel.reward import ExpertCache, advantages_v4, score_batch, score_completion
from objdel.train import serve

SIDE, K, N_VERT = 500, 6, 720
DISC_HF, CUP_HF = 0.0055, 0.013  # high-harmonic roughness: sets the traced vertex count (completion length)
WORKLOADS = {"riga": 1.0, "mixed": 0.4, "lidc": 0.0}


def _outline(cy, cx, r, asp, harm) -> np.ndarray:
    pts = []
    for i in range(N_VERT):
        t = 2.0 * math.pi * i / N_VERT
        rr = r * (1.0 + sum(a * math.cos(k * t + p) for k, a, p in harm))
        pts.append((round(16.0 * (cx + rr * math.cos(t))), round(16.0 * (cy + asp * rr * math.sin(t)))))
    m = np.zeros((SIDE, SIDE), np.uint8)
    cv2.fillPoly(m, [np.array(pts, np.int32)], 1, lineType=cv2.LINE_8, shift=4)
    return m.astype(bool)


def _harm(rng, ks, amp) -> list[tuple[int, float, float]]:
    return [(k, float(rng.normal(0.0, amp)), float(rng.uniform(0.0, 2.0 * math.pi))) for k in ks]


def synth_item(rng) -> dict[str, np.ndarray]:
    """{"disc": (K, SIDE, SIDE) bool, "cup": ...}: one RIGA-like item."""
    R, asp = float(rng.uniform(87.0, 116.0)), float(rng.uniform(0.99, 1.20))
    cy, cx = SIDE / 2 + float(rng.uniform(-4.0, 4.0)), SIDE / 2 + float(rng.uniform(-4.0, 4.0))
    disc_shape = _harm(rng, (2, 3, 4), 0.02)
    rho = float(rng.uniform(0.29, 0.64))
    cdy, cdx = float(rng.uniform(-0.07, 0.07)) * R, float(rng.uniform(-0.13, 0.11)) * R
    casp = float(rng.uniform(0.77, 1.18))
    cup_shape = _harm(rng, (2, 3), 0.04)
    discs, cups = [], []
    for _ in range(K):
        rj = R * (1.0 + float(rng.normal(0.0, 0.018)))
        oy, ox = (float(v) for v in rng.normal(0.0, 1.5, 2))
        dh = [(k, a + float(rng.normal(0.0, 0.025)), p) for k, a, p in disc_shape] + _harm(rng, range(5, 21), DISC_HF)
        discs.append(_outline(cy + oy, cx + ox, rj, asp, dh))
        rc = rho * R * (1.0 + float(rng.normal(0.0, 0.13)))
        ch = [(k, a + float(rng.normal(0.0, 0.05)), p) for k, a, p in cup_shape] + _harm(rng, range(4, 15), CUP_HF)
        cups.append(_outline(cy + cdy + 2 * oy, cx + cdx + 2 * ox, rc, casp, ch) & discs[-1])
    return {"disc": np.stack(discs), "cup": np.stack(cups)}


def corpus(n: int, seed: int = 0) -> list[dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    return [synth_item(rng) for _ in range(n)]


def corpus_digest(items) -> str:
    h = hashlib.sha256()
    for it in items:
        for name in ("disc", "cup"):
            h.update(np.packbits(it[name]).tobytes())
    return h.hexdigest()


def batch(riga, lidc, share: float, seed: int) -> list[tuple]:
    """round(32 share) RIGA prompts, the rest LIDC (with replacement only when a set is smaller), x 8 completions."""
    rng = np.random.default_rng(seed)
    n_r = round(32 * share)
    pick = lambda ps, n: [ps[i] for i in rng.choice(len(ps), n, replace=len(ps) < n)] if n else []  # noqa: E731
    return make_jobs(pick(riga, n_r) + pick(lidc, 32 - n_r), 1000 + seed)


def warm_wall(jobs, pool=None) -> tuple[float, str]:
    """(ms of the second of two runs of advantages_v4(score_batch(jobs, pool)), that run's hash)."""
    for _ in range(2):
        t0 = time.perf_counter()
        T = score_batch(jobs, pool)
        advantages_v4(T)
        ms = 1e3 * (time.perf_counter() - t0)
    return ms, digest(T)


def stats(xs) -> dict:
    a = np.asarray(xs, np.float64)
    return {"n": int(a.size), "median": float(np.median(a)), "p99": float(np.percentile(a, 99)), "max": float(a.max())}


def per_rollout_ms(jobs) -> list[float]:
    out = []
    for j in jobs:
        t0 = time.perf_counter()
        score_completion(*j)
        out.append(1e3 * (time.perf_counter() - t0))
    return out


def env() -> dict:
    cpu = platform.processor()
    if Path("/proc/cpuinfo").is_file():
        cpu = next((ln.split(":", 1)[1].strip() for ln in open("/proc/cpuinfo") if ln.startswith("model name")), cpu)
    return {
        "python": platform.python_version(),
        "machine": platform.machine(),
        "cpu": cpu,
        "nproc": os.cpu_count(),
        "loadavg_1m": os.getloadavg()[0],
        "time": time.strftime("%Y-%m-%d %H:%M:%S"),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="JSON report path")
    ap.add_argument("--root", default="data/prepared", help="prepared root with the LIDC records and caches")
    ap.add_argument("--synth", type=int, default=64, help="synthetic RIGA-like items")
    ap.add_argument("--lidc", type=int, default=200, help="LIDC prompts (first RL train items)")
    ap.add_argument("--batches", type=int, default=5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workers", default="2,4,8", help="fork pool sizes")
    a = ap.parse_args(argv)
    cv2.setNumThreads(0)
    res = {"env": env()}
    t0 = time.perf_counter()
    res["lidc_items_preloaded"] = serve.preload(a.root, ("lidc",))
    res["preload_s"] = time.perf_counter() - t0
    root = str(Path(a.root).resolve())
    lidc = [(SCHEMAS["lidc"], serve.CacheRef(k, root)) for k in list(serve.STORE[root])[: a.lidc]]  # train first
    t0 = time.perf_counter()
    items = corpus(a.synth)
    serve.STORE["synthetic"] = {
        f"synth/{i}": {n: ExpertCache.build(it[n]) for n in SCHEMAS["riga"]} for i, it in enumerate(items)
    }
    riga = [(SCHEMAS["riga"], serve.CacheRef(f"synth/{i}", "synthetic")) for i in range(a.synth)]
    res |= {"synthetic_digest": corpus_digest(items), "synthetic_build_s": time.perf_counter() - t0}
    batches = {w: [batch(riga, lidc, s, b) for b in range(a.batches)] for w, s in WORKLOADS.items()}
    with ThreadPoolExecutor(a.threads) as ex:
        for w, bl in batches.items():
            walls, hashes = zip(*map(warm_wall, bl), strict=True)
            tw, th = zip(*[warm_wall(jobs, ex) for jobs in bl], strict=True)
            res[w] = {
                "serial": stats(walls),
                "hashes": list(hashes),
                "per_rollout": stats([x for jobs in bl for x in per_rollout_ms(jobs)]),
                f"threads{a.threads}": {**stats(tw), "hash_equal": th == hashes},
            }
            print(w, "serial", round(np.median(walls), 1), "threads", round(np.median(tw), 1), flush=True)
    for n in map(int, a.workers.split(",")):
        pool = serve.fork_pool(n)
        for w, bl in batches.items():
            sw, pw, ph = [], [], []
            for jobs in bl:  # serial and pool interleaved per batch: the paired speed-up resists load swings
                sw.append(warm_wall(jobs)[0])
                ms, h = warm_wall(jobs, pool)
                pw.append(ms)
                ph.append(h)
            res[w][f"fork{n}"] = {
                **stats(pw),
                "hash_equal": ph == res[w]["hashes"],
                "paired_speedup": float(np.median(np.array(sw) / np.array(pw))),
            }
            print(
                w, f"fork{n}", round(np.median(pw), 1), "x", round(res[w][f"fork{n}"]["paired_speedup"], 2), flush=True
            )
        pool.shutdown()
    best = {
        w: min(v["median"] for k, v in res[w].items() if k == "serial" or k.startswith("fork"))
        for w in ("riga", "mixed")
    }
    f4 = res["riga"].get("fork4")
    res["gate"] = {
        "K3_best_ms": best,
        "K3_pass": max(best.values()) <= 250.0,
        "K4_p99_ms": max(res[w]["per_rollout"]["p99"] for w in WORKLOADS),
        "all_hashes_equal": all(
            v["hash_equal"] for w in WORKLOADS for v in res[w].values() if isinstance(v, dict) and "hash_equal" in v
        ),
        "d14_cpp_trigger": None if f4 is None else f4["median"] > 125.0 or f4["paired_speedup"] < 3.0,
    }
    res["gate"]["K4_pass"] = res["gate"]["K4_p99_ms"] <= 5.0
    res["env_end"] = env()
    print(res["gate"], flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
