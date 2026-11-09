"""K3/K4 reward-speed harness for a Vast host (H-A6-1 with H-K3SYN): needs no RIGA data.

RIGA prompts are synthetic RIGA-like items (synth.py, 64 items, caches built at start-up, ~20 s); LIDC prompts are
the first 200 RL train items of data/prepared (LIDC may go to any host, A3 §5). With --real-riga (datacenter hosts
only, A3 governance), the first 64 RIGA RL train items are used instead.

Reports, per workload (RIGA-only, mixed 40 % RIGA, LIDC-only; 10 batches of 32 prompts x 8, make_jobs):
- K3: warm batch wall ms (median, max), serial, 8 threads, and a fork pool (created after the preload) of 2 / 4 / 8
  workers; K4: per-rollout p50 / p99 / max ms (serial);
- the hash of every batch (reward_gates.digest) and whether every executor reproduced serial's;
- CPU model, nproc, load average; the gate verdict (K3: best of serial / pool <= 250 ms on RIGA-only and mixed;
  K4: p99 <= 5 ms).

    python -I -m experiments.phase2.lab.B4.k3_harness --out <dir> [--real-riga] [--workers 2,4,8] [--batches 10]
"""
import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor

import cv2
import numpy as np

from experiments.phase1.reward_gates import make_jobs
from experiments.phase2.lab.A6.common import dump, stats_ms
from objdel.reward import ExpertCache, score_completion

from . import serve
from .common import ROOT, env, first_keys, prompt_sets, warm_wall
from .synth import corpus, corpus_digest

NAMES = ("disc", "cup")
WORKLOADS = (("riga", 1.0), ("mixed", 0.4), ("lidc", 0.0))


def synthetic_prompts(n: int = 64, seed: int = 0) -> tuple[list, str]:
    """[(names, CacheRef)] of n synthetic items registered in serve.STORE as synth/<i>; and the corpus digest."""
    items = corpus(n, seed)
    for i, it in enumerate(items):
        serve.STORE[f"synth/{i}"] = {nm: ExpertCache.build(it[nm]) for nm in NAMES}
    return [(NAMES, serve.CacheRef(f"synth/{i}")) for i in range(n)], corpus_digest(items)


def batch(riga, lidc, share_riga: float, seed: int) -> list[tuple]:
    """One K3 batch: round(32 share) RIGA prompts and the rest LIDC, drawn with rng(seed); make_jobs(seed)."""
    rng = np.random.default_rng(seed)
    n_r = int(round(32 * share_riga))
    pr = [riga[i] for i in rng.choice(len(riga), n_r, replace=False)] if n_r else []
    pl = [lidc[i] for i in rng.choice(len(lidc), 32 - n_r, replace=False)] if n_r < 32 else []
    return make_jobs(pr + pl, 1000 + seed)


def per_rollout_ms(jobs) -> list[float]:
    out = []
    for j in jobs:
        t0 = time.perf_counter()
        score_completion(*j)
        out.append(1e3 * (time.perf_counter() - t0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--real-riga", action="store_true")
    ap.add_argument("--workers", default="2,4,8")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--batches", type=int, default=10)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    res = {"env": env(), "nproc": os.cpu_count(), "riga_source": "real" if a.real_riga else "synthetic"}
    t0 = time.perf_counter()
    keys = first_keys(64 if a.real_riga else 0, 200)
    serve.preload(ROOT, keys)
    lidc = prompt_sets(keys)["lidc"]
    if a.real_riga:
        riga = prompt_sets(keys)["riga"]
    else:
        riga, res["synthetic_digest"] = synthetic_prompts()
    res["setup_s"] = time.perf_counter() - t0
    batches = {w: [batch(riga, lidc, s, b) for b in range(a.batches)] for w, s in WORKLOADS}
    for w, bl in batches.items():
        lat = [x for jobs in bl for x in per_rollout_ms(jobs)]
        walls, hashes = zip(*[warm_wall(jobs) for jobs in bl], strict=True)
        res[w] = {"per_rollout": stats_ms(lat), "serial": stats_ms(walls), "hashes": list(hashes),
                  "chars": float(np.mean([len(j[0]) for jobs in bl for j in jobs]))}
        with ThreadPoolExecutor(a.threads) as ex:
            tw, th = zip(*[warm_wall(jobs, ex) for jobs in bl], strict=True)
        res[w][f"threads{a.threads}"] = {**stats_ms(tw), "hash_equal": list(th) == list(hashes)}
        print(w, "serial", round(np.median(walls), 1), "threads", round(np.median(tw), 1), flush=True)
    for n in map(int, a.workers.split(",")):
        pool = serve.fork_pool(n)
        for w, bl in batches.items():  # serial re-timed per batch, interleaved with the pool (paired speed-up)
            sw = []
            pw, ph = [], []
            for jobs in bl:
                sw.append(warm_wall(jobs)[0])
                p, h = warm_wall(jobs, pool)
                pw.append(p)
                ph.append(h)
            res[w][f"fork{n}"] = {**stats_ms(pw), "hash_equal": ph == res[w]["hashes"],
                                  "paired_speedup_median": float(np.median(np.array(sw) / np.array(pw)))}
            print(w, f"fork{n}", round(np.median(pw), 1), "paired x", round(res[w][f"fork{n}"]["paired_speedup_median"], 2),
                  flush=True)
        pool.shutdown()
    best = {w: min(v["median"] for k, v in res[w].items() if k == "serial" or k.startswith("fork")) for w in ("riga", "mixed")}
    res["gate"] = {"K3_best_ms": best, "K3_pass": all(v <= 250.0 for v in best.values()),
                   "K4_p99_ms": max(res[w]["per_rollout"]["p99"] for w, _ in WORKLOADS),
                   "all_hashes_equal": all(v.get("hash_equal", True) for w, _ in WORKLOADS for v in res[w].values()
                                           if isinstance(v, dict))}
    res["gate"]["K4_pass"] = res["gate"]["K4_p99_ms"] <= 5.0
    res["env_end"] = env()
    print(res["gate"], flush=True)
    dump(res, f"{a.out}/k3_harness_{res['riga_source']}_{env()['machine']}.json")


if __name__ == "__main__":
    main()
