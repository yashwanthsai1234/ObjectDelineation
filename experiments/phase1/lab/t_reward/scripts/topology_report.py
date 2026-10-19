"""Task 1: verify instance_topology against brute force and benchmark it.

    ../.venv/bin/python scripts/topology_report.py [n_random]
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cv2  # noqa: E402
import numpy as np  # noqa: E402

from objdel.checker.bitquad import bitquads  # noqa: E402
from objdel.checker.topology import betti  # noqa: E402
from rewardcore import fixtures, synth  # noqa: E402
from rewardcore.reference import instance_topology_bruteforce  # noqa: E402
from rewardcore.topology import instance_topology, instance_topology_euler  # noqa: E402

cv2.setNumThreads(1)


def verify(n_random):
    rng = np.random.default_rng(2026)
    stats = {}
    for _ in range(n_random):
        name, m = synth.random_mask(rng, 8, 256)
        ref = instance_topology_bruteforce(m)
        s = stats.setdefault(name, [0, 0, 0, 0])
        s[0] += 1
        s[1] += instance_topology(m) != ref
        s[2] += instance_topology_euler(m).as_list() != ref
        s[3] += len(ref)
    print(f"random masks (8-256 px), seed 2026: {n_random}")
    for name, (n, bad, bad_e, comps) in sorted(stats.items()):
        print(f"  {name:12s} n={n:5d} components={comps:8d} mismatches spaghetti={bad} euler={bad_e}")
    for ds in ("riga", "lidc"):
        n = bad = 0
        for stack in fixtures.load(ds).values():
            for m in stack:
                for mm in (m, ~m):
                    n += 1
                    bad += instance_topology(mm) != instance_topology_bruteforce(mm)
        print(f"fixture {ds}: {n} masks (rater masks and complements), mismatches={bad}")


def bench_one(fn, m, reps):
    fn(m)
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter_ns()
        fn(m)
        ts.append(time.perf_counter_ns() - t0)
    return np.median(ts) / 1e3


def current_path(m):
    return betti(m), bitquads(m).chi8


def bench():
    rng = np.random.default_rng(7)
    print("\nsingle thread, median us per mask (20 masks x reps)")
    print(f"{'size':>6} {'kind':>6} | {'objdel betti+bitquad':>21} {'brute per-inst':>15} {'spaghetti per-inst':>19} {'euler per-inst':>15} | speedup vs current")
    for size in (64, 256, 1024):
        for kind in ("blob", "noisy"):
            masks = [synth.GENERATORS[kind](rng, size, size) for _ in range(20)]
            reps = {64: 200, 256: 50, 1024: 5}[size]
            cur = np.median([bench_one(current_path, m, reps) for m in masks])
            brute = np.median([bench_one(instance_topology_bruteforce, m, max(1, reps // 10)) for m in masks[:5]])
            new = np.median([bench_one(instance_topology, m, reps) for m in masks])
            eul = np.median([bench_one(lambda x: instance_topology_euler(x).as_list(), m, reps) for m in masks])
            print(f"{size:>5}² {kind:>6} | {cur:21.1f} {brute:15.1f} {new:19.1f} {eul:15.1f} | {cur / new:5.1f}x")
    sizes = []
    for ds in ("riga", "lidc"):
        for stack in fixtures.load(ds).values():
            sizes.append((ds, stack[0]))
    for ds in ("riga", "lidc"):
        ms = [m for d, m in sizes if d == ds][:100]
        cur = np.median([bench_one(current_path, m, 20) for m in ms])
        new = np.median([bench_one(instance_topology, m, 20) for m in ms])
        print(f"fixture {ds} rater masks (median shape {tuple(int(v) for v in np.median([m.shape for m in ms], 0))}): current {cur:.1f} us, spaghetti {new:.1f} us, {cur / new:.1f}x")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
    if n:
        verify(n)
    bench()
