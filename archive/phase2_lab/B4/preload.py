"""B4 H-A6-4: preload every prepared RL train cache at start-up, then replay TRL-like steps through TRLReward with
the production LRU loader (cache_loader, maxsize 64) and with the preloaded store; optionally a fork pool over the
full preload (copy-on-write memory at full size).

Steps: a seeded permutation of the RL train items of one dataset (TRL's RepeatSampler permutes the rows and repeats
each G times in a row), 32 prompts per step, 8 make_jobs completions per prompt, passed to the 8 TRLReward functions
with TRL's kwargs. Per step: the LRU path's call (cold: new prompts), the preloaded path's call and plain
advantages_v4(score_batch(...)) on the same jobs (warm walls), the LRU misses and the store's loads.

    .venv/bin/python -I -m experiments.phase2.lab.B4.preload --out <dir> [--steps 10 --threads 2 --pool 4]

Memory: holds every RIGA + LIDC train cache (~4.8 GB) in one process; with --pool, plus the workers' private pages.
"""
import argparse
import os
import resource

import cv2
import numpy as np

from experiments.phase1.reward_gates import make_jobs
from experiments.phase2.lab.A6.common import dump, stats_ms
from objdel.geometry import SCHEMAS
from objdel.reward import advantages_v4, score_batch
from objdel.train.dataset import cache_loader
from objdel.train.rewards import AGGREGATOR, TRLReward

from . import serve
from .common import EOS, ROOT, env, trl_call, warm_wall
from .prereg import verdicts


def steps(ds: str, keys: list[str], n: int, seed: int) -> list[list[tuple]]:
    perm = np.random.default_rng(seed).permutation(len(keys))
    return [make_jobs([(SCHEMAS[ds], serve.CacheRef(keys[i])) for i in perm[32 * s:32 * s + 32]], 2000 + s)
            for s in range(n)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--pool", type=int, default=0)
    ap.add_argument("--max-riga", type=int, default=0, help="preload only the first N RIGA items (container memory)")
    a = ap.parse_args()
    os.environ.setdefault("B4_SCRATCH", a.out)
    cv2.setNumThreads(0)
    res = {"env": env(), "rss_start_mb": serve.mem_mb(os.getpid())["rss"]}
    dig = serve.rl_digests(ROOT, ("riga", "lidc"), ("train",))
    if a.max_riga:
        drop = [k for k in dig if k.startswith("riga")][a.max_riga:]
        dig = {k: v for k, v in dig.items() if k not in set(drop)}
    for ds in ("riga", "lidc"):
        sub = {k: v for k, v in dig.items() if k.startswith(ds)}
        res[f"preload_{ds}"] = {"items": len(sub), "caches": sum(map(len, sub.values())),
                                "s": serve.preload(ROOT, sub, a.threads)}
        res[f"rss_after_{ds}_mb"] = serve.mem_mb(os.getpid())["rss"]
        print(ds, res[f"preload_{ds}"], res[f"rss_after_{ds}_mb"], flush=True)
    res["preload_s"] = res["preload_riga"]["s"] + res["preload_lidc"]["s"]
    res["maxrss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20 if os.uname().sysname == "Darwin" else 2**10)
    loads0 = serve.LOADS["n"]
    for ds in ("riga", "lidc"):
        keys = sorted(k for k in dig if k.startswith(ds))
        lru = TRLReward(cache_loader(ROOT, maxsize=64), {EOS})
        pre = TRLReward(serve.ref_loader, {EOS})
        rows = []
        for s, jobs in enumerate(steps(ds, keys, a.steps, seed=42)):
            m0 = lru.caches.cache_info().misses
            t_lru, o_lru = trl_call(lru, jobs)
            miss = lru.caches.cache_info().misses - m0
            trl_call(pre, jobs)  # warm, as warm_wall
            t_pre, o_pre = trl_call(pre, jobs)
            plain, h = warm_wall(jobs)
            rows.append({"step": s, "lru_s": t_lru, "pre_s": t_pre, "plain_s": plain / 1e3, "misses": miss,
                         "same": o_lru == o_pre, "hash": h, "adv_equal_plain": o_pre[AGGREGATOR] == advantages_v4(score_batch(jobs)).tolist()})
            print(ds, rows[-1], flush=True)
        res[ds] = {"rows": rows, "lru_s": stats_ms([r["lru_s"] for r in rows]),
                   "pre_s": stats_ms([r["pre_s"] for r in rows]), "plain_s": stats_ms([r["plain_s"] for r in rows]),
                   "misses_mean": float(np.mean([r["misses"] for r in rows])),
                   "delta_median_s": float(np.median([r["lru_s"] - r["pre_s"] for r in rows])),
                   "pre_over_plain_median": float(np.median([r["pre_s"] / r["plain_s"] for r in rows]))}
    res["loads_during_steps"] = serve.LOADS["n"] - loads0
    res["verdicts"] = verdicts({
        "P4-1 preload wall s (Mac, 2 threads)": res["preload_s"],
        "P4-2 RSS after preload GB": res["rss_after_lidc_mb"] / 1024,
        "P4-3a preloaded loads during 10 steps": res["loads_during_steps"],
        "P4-3b LRU misses per step (mean)": res["riga"]["misses_mean"],
        "P4-4 LRU minus preloaded step reward call, median s": res["riga"]["delta_median_s"],
        "P4-5 preloaded TRLReward / plain score_batch, median": res["riga"]["pre_over_plain_median"],
        "P4-6 step hashes identical LRU vs preloaded": all(r["same"] for ds in ("riga", "lidc") for r in res[ds]["rows"])})
    print(res["verdicts"], flush=True)
    dump(res, f"{a.out}/preload_{env()['machine']}.json")
    if a.pool:
        keys = sorted(k for k in dig if k.startswith("riga"))
        bl = steps("riga", keys, min(10, len(keys) // 32), seed=7)
        ser = [warm_wall(j) for j in bl]
        pool = serve.fork_pool(a.pool)
        m0 = serve.workers_mem(pool, a.pool)
        par = [warm_wall(j, pool) for j in bl]
        m1 = serve.workers_mem(pool, a.pool)
        res["full_pool"] = {"workers": a.pool, "serial_ms": stats_ms([s for s, _ in ser]),
                            "pool_ms": stats_ms([p for p, _ in par]),
                            "hash_equal": [h for _, h in ser] == [h for _, h in par], "mem_start": m0, "mem_end": m1,
                            "parent": serve.mem_mb(os.getpid()),
                            "uss_growth_max": max(m1[p]["uss"] - m0[p]["uss"] for p in m1 if p in m0)}
        pool.shutdown()
        print(res["full_pool"], flush=True)
        res["env_end"] = env()
        dump(res, f"{a.out}/preload_{env()['machine']}.json")


if __name__ == "__main__":
    main()
