"""Task 5, KPIs K3/K4: wall time to score a GRPO batch of 256 rollouts, p99 per-rollout latency.

Batches: 32 prompts x 8 rollouts, prompts drawn uniformly from the fixtures ("mixed": ~40 % RIGA
crops of ~470x430 px, ~60 % LIDC crops of ~27x28 px) and RIGA-only (worst case). Per-prompt
caches are built once (timed separately: they can be precomputed offline for the whole dataset).

    ../.venv/bin/python scripts/k3_throughput.py [n_batches]
"""
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cv2  # noqa: E402
import numpy as np  # noqa: E402

from rewardcore.batch import BatchScorer, advantages_v4  # noqa: E402
from rewardcore.reward import ExpertCache, RewardConfig  # noqa: E402
from rewardcore.workload import make_batch  # noqa: E402

cv2.setNumThreads(1)
CFG = RewardConfig(band_formula="F5")
N_WORKERS = int(__import__("os").environ.get("K3_WORKERS", 8))  # iteration 3: run with 4 (shared machine)
BACKENDS = (("serial", 1), ("thread", N_WORKERS), ("process", N_WORKERS))


def main(n_batches=10):
    print(f"load average before: {os.getloadavg()} (8-core M3; other jobs share the machine)")
    for dataset in (None, "riga"):
        name = "mixed" if dataset is None else "riga-only"
        batches = [make_batch(1000 + b, dataset=dataset) for b in range(n_batches)]
        shapes = np.array([m.shape for m in batches[0][2]])
        print(f"\n[{name}] crop px per rollout: median {np.median(shapes.prod(1)):.0f}, "
              f"max {shapes.prod(1).max()}, mean {shapes.prod(1).mean():.0f}")
        t0 = time.perf_counter()
        build_ms, build_cpu = [], []
        caches_per_batch = []
        for st_list, _, _ in batches:
            cs = []
            for st in st_list:
                t1, c1 = time.perf_counter(), time.thread_time()
                cs.append(ExpertCache.build(st))
                build_ms.append((time.perf_counter() - t1) * 1e3)
                build_cpu.append((time.thread_time() - c1) * 1e3)
            caches_per_batch.append(cs)
        print(f"  cache build per prompt (single thread): wall median {np.median(build_ms):.1f} ms, "
              f"p99 {np.percentile(build_ms, 99):.1f} ms | CPU median {np.median(build_cpu):.1f} ms, "
              f"p99 {np.percentile(build_cpu, 99):.1f} ms; CPU per 32 prompts {np.sum(build_cpu) / n_batches:.0f} ms"
              f" (wall {1e3 * (time.perf_counter() - t0) / n_batches:.0f} ms)")
        for backend, workers in BACKENDS:
            walls, lats, cpus, batch_cpu = [], [], [], []
            for b, (_, idx, masks) in enumerate(batches):
                with BatchScorer(caches_per_batch[b], CFG, workers, "thread" if backend == "serial" else backend) as sc:
                    sc(idx, masks)  # warm-up
                    t1 = time.perf_counter()
                    terms, lat, cpu = sc(idx, masks)
                    advantages_v4(terms)  # reward v4 group advantages (included in the batch wall time)
                    walls.append((time.perf_counter() - t1) * 1e3)
                    lats.append(lat / 1e6)
                    cpus.append(cpu / 1e6)
                    batch_cpu.append(cpu.sum() / 1e6)
            lats, cpus = np.concatenate(lats), np.concatenate(cpus)
            k3 = "PASS" if np.median(walls) <= 250 else "FAIL"
            k4 = "PASS" if np.percentile(lats, 99) <= 5 else "FAIL"
            print(f"  {backend:7s} x{workers}: batch wall median {np.median(walls):6.1f} ms, max {np.max(walls):6.1f} ms [K3 {k3}]"
                  f" | per-rollout wall p50 {np.median(lats):.2f} ms, p99 {np.percentile(lats, 99):.2f} ms [K4 {k4}]"
                  f" | CPU per rollout p50 {np.median(cpus):.2f} p99 {np.percentile(cpus, 99):.2f} max {cpus.max():.2f} ms,"
                  f" CPU per batch {np.median(batch_cpu):.0f} ms")
    print(f"\nload average after: {os.getloadavg()}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
