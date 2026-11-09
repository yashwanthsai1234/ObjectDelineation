"""Shared helpers of the B4 lab: prompt sets as CacheRefs, K3 batches, warm wall times, environment record.

Run every B4 script as a module with the isolated interpreter (data/ is untrusted):
    .venv/bin/python -I -m experiments.phase2.lab.B4.<script> --out <scratch>/phase2/B4
"""
import os
import platform
import time
from pathlib import Path

import numpy as np

from experiments.phase1.reward_gates import digest, make_jobs
from objdel.geometry import SCHEMAS
from objdel.reward import advantages_v4, score_batch

from . import serve

ROOT = Path("data/prepared").resolve()
EOS = 0  # stop id of the simulated completions


def trl_call(reward, jobs) -> tuple[float, dict]:
    """The 8 reward functions in TRL's order with TRL's kwargs (objdel.train.rewards); (s, {name: list})."""
    funcs, _ = reward.funcs()
    kw = {"prompts": [[{"role": "user", "content": "p"}]] * len(jobs),
          "completions": [[{"role": "assistant", "content": j[0]}] for j in jobs],
          "completion_ids": [[7, EOS] for _ in jobs], "key": [j[3].key for j in jobs], "W": [j[4] for j in jobs],
          "H": [j[5] for j in jobs], "names": [list(j[2]) for j in jobs], "trainer_state": None,
          "log_metric": lambda *_: None, "log_extra": lambda *_: None}
    t0 = time.perf_counter()
    out = {f.__name__: f(**kw) for f in funcs}
    return time.perf_counter() - t0, out


def first_keys(n_riga: int, n_lidc: int, root=ROOT) -> dict[str, dict[str, str]]:
    """The A6/M6 K3 prompt sets: the first n RL train items of each dataset, with their cache digests."""
    out = {}
    for ds, n in (("riga", n_riga), ("lidc", n_lidc)):
        out |= dict(list(serve.rl_digests(root, (ds,)).items())[:n])
    return out


def prompt_sets(keys) -> dict[str, list]:
    """{"riga": [(names, CacheRef)], "lidc": [...]} for keys already in serve.STORE."""
    out = {"riga": [], "lidc": []}
    for k in keys:
        ds = "riga" if k.startswith("riga") else "lidc"
        out[ds].append((SCHEMAS[ds], serve.CacheRef(k)))
    return out


def k3_batches(prompts, share_riga: float, n: int, seed: int) -> list[list[tuple]]:
    """A6 prof_reward.batches: n batches of 32 prompts (round(32 share) RIGA) x 8 completions (make_jobs)."""
    rng = np.random.default_rng(seed)
    out = []
    for b in range(n):
        n_r = int(round(32 * share_riga))
        pr = [prompts["riga"][i] for i in rng.choice(len(prompts["riga"]), n_r, replace=False)] if n_r else []
        pl = [prompts["lidc"][i] for i in rng.choice(len(prompts["lidc"]), 32 - n_r, replace=False)] if n_r < 32 else []
        out.append(make_jobs(pr + pl, 1000 + b))
    return out


def warm_wall(jobs, pool=None) -> tuple[float, str]:
    """(ms of the second of two runs of advantages_v4(score_batch(jobs, pool)), hash of that run)."""
    for _ in range(2):
        t0 = time.perf_counter()
        T = score_batch(jobs, pool)
        advantages_v4(T)
        dt = 1e3 * (time.perf_counter() - t0)
    return dt, digest(T)


def timed(jobs, pool=None) -> tuple[float, float, str]:
    """warm_wall plus CPU: (wall ms, CPU ms of this process + the pool's workers, hash) of the second of two runs."""
    pids = serve.worker_pids(pool) if pool is not None else []
    for _ in range(2):
        c0, w0 = time.process_time(), serve.cpu_s(pids)
        t0 = time.perf_counter()
        T = score_batch(jobs, pool)
        advantages_v4(T)
        dt = 1e3 * (time.perf_counter() - t0)
        cpu = 1e3 * (time.process_time() - c0 + serve.cpu_s(pids) - w0)
    return dt, cpu, digest(T)


def env() -> dict:
    return {"python": platform.python_version(), "machine": platform.machine(), "system": platform.system(),
            "cpus": os.cpu_count(), "loadavg_1m": os.getloadavg()[0], "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "cpuinfo": _cpu()}


def _cpu() -> str:
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()
