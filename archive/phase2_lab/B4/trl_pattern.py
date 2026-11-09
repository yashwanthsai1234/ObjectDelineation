"""B4 H-A6-2 (P2-6) and pool robustness: TRLReward on a forked pool under TRL's reward-function call pattern.

1. pattern: 10 generation batches (mixed workload, 256 completions) and one eval-sized batch (64); per batch the 8
   TRLReward functions are called in TRL's order with TRL's kwargs (preload.trl_call). Counted: pool dispatches per
   batch (must be 1), outputs == the serial TRLReward's, worker pids unchanged. For "fork" and "forkserver".
2. threads at fork: OS thread count and Python's fork DeprecationWarning when the pool is forked (a) right after the
   preload and (b) after torch has started its intra-op thread pool (a CPU matmul; no CUDA here); then a batch is
   scored with a timeout (a deadlock would show as a timeout).
3. failures: a task that raises (unknown key) -> exception in the parent, pool still usable; a worker SIGKILLed
   mid-batch -> ProcessPoolExecutor raises BrokenProcessPool (time to raise), multiprocessing.Pool.map_async does
   not return (timeout).

    .venv/bin/python -I -m experiments.phase2.lab.B4.trl_pattern --out <dir>
"""
import argparse
import multiprocessing
import os
import signal
import time
import warnings
from concurrent.futures import TimeoutError as FutTimeout
from concurrent.futures.process import BrokenProcessPool

import cv2
import psutil

from experiments.phase2.lab.A6.common import dump
from objdel.reward.score import _score_chunk
from objdel.train.rewards import TRLReward

from . import serve
from .common import EOS, ROOT, env, first_keys, k3_batches, prompt_sets, trl_call
from .prereg import verdicts


class Counting:
    """Executor proxy that counts map() dispatches."""

    def __init__(self, ex):
        self.ex, self.n = ex, 0

    def map(self, *a, **k):
        self.n += 1
        return self.ex.map(*a, **k)


def pattern(label, pool, batches, loader) -> dict:
    cp = Counting(pool)
    par, ser = TRLReward(loader, {EOS}, pool=cp), TRLReward(serve.ref_loader, {EOS})
    pids0 = serve.worker_pids(pool)
    rows = []
    for jobs in batches:
        n0 = cp.n
        tp, op = trl_call(par, jobs)
        ts, os_ = trl_call(ser, jobs)
        rows.append({"n": len(jobs), "dispatches": cp.n - n0, "equal": op == os_, "pool_s": tp, "serial_s": ts})
    ok = all(r["dispatches"] == 1 and r["equal"] for r in rows) and serve.worker_pids(pool) == pids0
    print(label, ok, [(r["dispatches"], r["equal"], round(r["pool_s"], 3), round(r["serial_s"], 3)) for r in rows],
          flush=True)
    return {"rows": rows, "pids_stable": serve.worker_pids(pool) == pids0, "ok": ok}


def fork_with_warnings(make) -> tuple[object, list[str], int]:
    nthreads = psutil.Process().num_threads()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        pool = make()
    return pool, [f"{x.category.__name__}: {str(x.message)[:120]}" for x in w], nthreads


def score_with_timeout(pool, jobs, timeout=120.0) -> str:
    try:
        futs = [pool.submit(_score_chunk, jobs[i:i + 8]) for i in range(0, len(jobs), 8)]
        [f.result(timeout=timeout) for f in futs]
        return "ok"
    except FutTimeout:
        return "timeout"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.environ.setdefault("B4_SCRATCH", a.out)
    cv2.setNumThreads(0)
    res = {"env": env()}
    keys = first_keys(64, 200)
    serve.preload(ROOT, keys)
    prompts = prompt_sets(keys)
    batches = k3_batches(prompts, 0.4, 10, seed=11)
    batches.append(batches[0][:64])  # an eval-sized call (8 prompts)

    # 1 + 2a: fork right after the preload
    pool, w, nt = fork_with_warnings(lambda: serve.fork_pool(4))
    res["fork_after_preload"] = {"os_threads_before_fork": nt, "warnings": w}
    res["pattern_fork"] = pattern("fork", pool, batches, serve.ref_loader)

    # 3: failures on the fork pool
    t0 = time.perf_counter()
    try:
        list(pool.map(_score_chunk, [[(*batches[0][0][:3], serve.CacheRef("nope"), *batches[0][0][4:])]]))
        res["task_exception"] = "no exception"
    except KeyError as e:
        res["task_exception"] = f"KeyError raised in parent after {time.perf_counter() - t0:.3f} s: {e}"
    res["usable_after_exception"] = score_with_timeout(pool, batches[1])
    victim = serve.worker_pids(pool)[0]
    futs = [pool.submit(_score_chunk, batches[2][i:i + 8]) for i in range(0, 256, 8)]
    time.sleep(0.05)
    t0 = time.perf_counter()
    os.kill(victim, signal.SIGKILL)
    try:
        [f.result(timeout=60) for f in futs]
        res["ppe_worker_killed"] = "completed (victim was idle)"
    except BrokenProcessPool as e:
        res["ppe_worker_killed"] = f"BrokenProcessPool after {time.perf_counter() - t0:.3f} s: {str(e)[:100]}"
    except FutTimeout:
        res["ppe_worker_killed"] = "timeout"
    pool.shutdown(wait=False, cancel_futures=True)
    mp = multiprocessing.get_context("fork").Pool(2)
    r = mp.map_async(_score_chunk, [batches[3][i:i + 8] for i in range(0, 256, 8)])
    time.sleep(0.05)
    os.kill(mp._pool[0].pid, signal.SIGKILL)
    try:
        r.get(timeout=30)
        res["mp_pool_worker_killed"] = "completed (map finished or the worker was replaced in time)"
    except multiprocessing.TimeoutError:
        res["mp_pool_worker_killed"] = "no result after 30 s (map_async never returns)"
    mp.terminate()
    print({k: res[k] for k in ("task_exception", "usable_after_exception", "ppe_worker_killed",
                               "mp_pool_worker_killed")}, flush=True)

    # 1: forkserver
    pool = serve.forkserver_pool(4, ROOT, keys)
    from multiprocessing import forkserver
    try:
        res["forkserver_os_threads"] = psutil.Process(forkserver._forkserver._forkserver_pid).num_threads()
    except psutil.AccessDenied:
        res["forkserver_os_threads"] = "access denied"
    res["pattern_forkserver"] = pattern("forkserver", pool, batches, serve.CacheRef)
    pool.shutdown()

    # 2b: fork after torch's intra-op pool has started
    import torch
    torch.set_num_threads(4)
    x = torch.randn(1024, 1024)
    for _ in range(5):
        x = (x @ x).tanh()
    pool, w, nt = fork_with_warnings(lambda: serve.fork_pool(2))
    res["fork_after_torch_threads"] = {"os_threads_before_fork": nt, "warnings": w,
                                       "torch_threads": torch.get_num_threads(),
                                       "score": score_with_timeout(pool, batches[4])}
    pool.shutdown()
    pool, w, nt = fork_with_warnings(lambda: serve.forkserver_pool(2, ROOT, keys))
    res["forkserver_after_torch_threads"] = {"os_threads_before": nt, "warnings": w,
                                             "score": score_with_timeout(pool, batches[5])}
    pool.shutdown()
    print({k: res[k] for k in ("fork_after_preload", "fork_after_torch_threads", "forkserver_after_torch_threads")},
          flush=True)
    res["verdicts"] = verdicts({"P2-6 TRL pattern: one dispatch per batch, outputs ==, 10 steps":
                                res["pattern_fork"]["ok"] and res["pattern_forkserver"]["ok"]})
    res["env_end"] = env()
    print(res["verdicts"], flush=True)
    dump(res, f"{a.out}/trl_pattern_{env()['machine']}.json")


if __name__ == "__main__":
    main()
