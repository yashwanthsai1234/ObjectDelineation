"""B4 H-A6-2: 256-rollout reward batches serial vs a process pool forked after the preload (fork / forkserver),
2 / 4 / 8 workers, interleaved per batch with serial; hashes, worker memory, pickling and pipe overhead; a 4-thread
pool for reference.

    .venv/bin/python -I -m experiments.phase2.lab.B4.pool_bench --out <dir> [--riga 64 --lidc 200 --batches 10]

Peak memory: the parent holds the subset's caches (64 RIGA items ~0.5 GB + 200 LIDC) and the workers share them.
"""
import argparse
import os
import pickle
import time
from concurrent.futures import ThreadPoolExecutor

import cv2
import numpy as np

from experiments.phase1.reward_gates import digest
from experiments.phase2.lab.A6.common import dump, stats_ms
from objdel.reward import score_batch
from objdel.reward.score import REWARD_V4

from . import serve
from .common import ROOT, env, first_keys, k3_batches, prompt_sets, timed, warm_wall

WORKLOADS = (("riga", 1.0), ("mixed", 0.4), ("lidc", 0.0))


def overhead_ms(jobs, pool) -> float:
    """Warm wall of the same chunked map with a no-op scorer: pickling both ways + pipe + scheduling."""
    chunks = [jobs[i:i + REWARD_V4.group] for i in range(0, len(jobs), REWARD_V4.group)]
    out = []
    for _ in range(3):
        t0 = time.perf_counter()
        list(pool.map(serve._noop_chunk, chunks))
        out.append(1e3 * (time.perf_counter() - t0))
    return float(np.median(out[1:]))


def pickle_cost(jobs) -> dict:
    t0 = time.perf_counter()
    b = pickle.dumps(jobs, protocol=pickle.HIGHEST_PROTOCOL)
    t1 = time.perf_counter()
    pickle.loads(b)
    t2 = time.perf_counter()
    res = score_batch(jobs)
    rows = [{k: res[k][i] for k in res} for i in range(len(jobs))]
    rb = pickle.dumps(rows, protocol=pickle.HIGHEST_PROTOCOL)
    return {"jobs_bytes": len(b), "dumps_ms": 1e3 * (t1 - t0), "loads_ms": 1e3 * (t2 - t1), "results_bytes": len(rb)}


def run_pool(label, make_pool, n, batches, parent_pid) -> dict:
    t0 = time.perf_counter()
    pool = make_pool(n)
    start_s = time.perf_counter() - t0
    pids = serve.worker_pids(pool)
    mem0 = serve.workers_mem(pool, n)
    r = {"start_s": start_s, "parent_mem_mb": serve.mem_mb(parent_pid)}
    for w, bl in batches.items():
        ser, par, same, cs, cp = [], [], [], [], []
        for jobs in bl:
            s, c1, hs = timed(jobs)
            p, c2, hp = timed(jobs, pool)
            ser.append(s)
            par.append(p)
            cs.append(c1)
            cp.append(c2)
            same.append(hs == hp)
        r[w] = {"serial": stats_ms(ser), "pool": stats_ms(par), "speedup_median": float(np.median(ser) / np.median(par)),
                "serial_cpu": stats_ms(cs), "pool_cpu": stats_ms(cp),
                "cpu_ratio_median": float(np.median(np.array(cp) / np.array(cs))),
                "paired_ratio_median": float(np.median(np.array(ser) / np.array(par))), "hash_equal": all(same),
                "overhead_ms": overhead_ms(bl[0], pool)}
        print(label, n, w, round(np.median(ser), 1), round(np.median(par), 1), all(same), "cpu", round(np.median(cs), 1),
              round(np.median(cp), 1), "load", round(os.getloadavg()[0], 1), flush=True)
    mem1 = serve.workers_mem(pool, n)
    r["worker_mem_mb"] = {"start": mem0, "end": mem1, "uss_growth_max": max(
        mem1[p].get("uss", np.nan) - mem0[p].get("uss", np.nan) for p in pids if p in mem0 and p in mem1)}
    if label == "forkserver":
        from multiprocessing import forkserver
        r["forkserver_mem_mb"] = serve.mem_mb(forkserver._forkserver._forkserver_pid)
    r["pids_stable"] = serve.worker_pids(pool) == pids
    pool.shutdown()
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--riga", type=int, default=64)
    ap.add_argument("--lidc", type=int, default=200)
    ap.add_argument("--batches", type=int, default=10)
    ap.add_argument("--workers", default="2,4,8")
    ap.add_argument("--methods", default="fork,forkserver")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.environ.setdefault("B4_SCRATCH", a.out)
    cv2.setNumThreads(0)
    res = {"env": env()}
    keys = first_keys(a.riga, a.lidc)
    res["preload_s"] = serve.preload(ROOT, keys)
    res["parent_mem_after_preload_mb"] = serve.mem_mb(os.getpid())
    prompts = prompt_sets(keys)
    batches = {w: k3_batches(prompts, share, a.batches, seed=7) for w, share in WORKLOADS}
    # CacheRef jobs score like dict jobs (same objects behind them)
    j = batches["mixed"][0]
    res["cacheref_equals_dict"] = digest(score_batch(j)) == digest(score_batch([(*x[:3], dict(x[3]), *x[4:]) for x in j]))
    res["pickle"] = {w: pickle_cost(bl[0]) for w, bl in batches.items()}
    res["pickle_dict_jobs_riga_bytes"] = serve.pickled_bytes([(*x[:3], dict(x[3]), *x[4:]) for x in batches["riga"][0]])
    with ThreadPoolExecutor(4) as ex:  # reference: A6's thread result
        th = [warm_wall(jobs, ex)[0] for jobs in batches["riga"]]
        se = [warm_wall(jobs)[0] for jobs in batches["riga"]]
    res["threads4_riga"] = {"threads": stats_ms(th), "serial": stats_ms(se)}
    makers = {"fork": serve.fork_pool, "forkserver": lambda n: serve.forkserver_pool(n, ROOT, keys)}
    for m in a.methods.split(","):
        for n in map(int, a.workers.split(",")):
            res[f"{m}{n}"] = run_pool(m, makers[m], n, batches, os.getpid())
            dump(res, f"{a.out}/pool_bench_{env()['machine']}{a.tag}.json")
    res["env_end"] = env()
    dump(res, f"{a.out}/pool_bench_{env()['machine']}{a.tag}.json")


if __name__ == "__main__":
    main()
