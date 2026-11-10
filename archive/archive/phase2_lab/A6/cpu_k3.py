"""A6: K3 workload timed in process CPU time (time.process_time) as well as wall time, serial, on the prepared and
the fixture caches: CPU time is far less sensitive than wall time to other processes on a shared machine.

    .venv/bin/python -I -m experiments.phase2.lab.A6.cpu_k3 --out <dir>
"""
import argparse
import time

import cv2
import numpy as np

from objdel.reward import advantages_v4, score_batch

from .common import dump, env
from .control_k3 import fixture_prompts
from .prof_reward import batches, load_prompts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    res = {"env": env()}
    sets = {"prepared": load_prompts(64, 200)[0], "fixtures": fixture_prompts()}
    for name, prompts in sets.items():
        for label, share in (("mixed", 0.4), ("riga", 1.0), ("lidc", 0.0)):
            cpu, wall = [], []
            for jobs in batches(prompts, share, 6, seed=21):
                score_batch(jobs)  # warm
                for _ in range(2):
                    c0, w0 = time.process_time(), time.perf_counter()
                    advantages_v4(score_batch(jobs))
                    cpu.append(1e3 * (time.process_time() - c0))
                    wall.append(1e3 * (time.perf_counter() - w0))
            res[f"{name}/{label}"] = {"cpu_ms_median": float(np.median(cpu)), "cpu_ms_max": float(np.max(cpu)),
                                      "wall_ms_median": float(np.median(wall))}
            print(name, label, res[f"{name}/{label}"], flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/cpu_k3.json")


if __name__ == "__main__":
    main()
