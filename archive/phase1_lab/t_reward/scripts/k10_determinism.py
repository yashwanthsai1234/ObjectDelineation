"""Task 6, KPI K10: score a 256-rollout batch under several execution settings, print SHA-256 digests.

    ../.venv/bin/python scripts/k10_determinism.py
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cv2  # noqa: E402
import numpy as np  # noqa: E402

from rewardcore.batch import BatchScorer, advantages_v4, digest  # noqa: E402
from rewardcore.reward import ExpertCache, RewardConfig  # noqa: E402
from rewardcore.workload import make_batch  # noqa: E402

CFG = RewardConfig(band_formula="F5")
SNIPPET = (
    "import sys; sys.path.insert(0, %r)\n"
    "import numpy as np\n"
    "from rewardcore.batch import BatchScorer, advantages_v4, digest\n"
    "from rewardcore.reward import ExpertCache, RewardConfig\n"
    "from rewardcore.workload import make_batch\n"
    "stacks, idx, masks = make_batch(seed=%d)\n"
    "with BatchScorer([ExpertCache.build(s) for s in stacks], RewardConfig(band_formula='F5'), %d, %r) as sc:\n"
    "    t = sc(idx, masks)[0]\n"
    "    print(digest(np.concatenate([t.ravel(), advantages_v4(t)])))\n"
)


def main():
    results = {}
    for seed in (11, 12):
        stacks, idx, masks = make_batch(seed=seed)
        caches = [ExpertCache.build(s) for s in stacks]
        nw = int(__import__("os").environ.get("K10_WORKERS", 8))  # iteration 3: 4 (shared machine)
        for label, workers, backend in (("serial run 1", 1, "thread"), ("serial run 2", 1, "thread"),
                                        (f"{nw} threads", nw, "thread"), ("1 process", 1, "process"),
                                        (f"{nw} processes", nw, "process")):
            with BatchScorer(caches, CFG, workers, backend) as sc:
                t = sc(idx, masks)[0]
                results[(seed, label)] = digest(np.concatenate([t.ravel(), advantages_v4(t)]))
        cv2.setNumThreads(nw)
        with BatchScorer(caches, CFG, 1, "thread") as sc:  # BatchScorer pins OpenCV to 1 thread ...
            cv2.setNumThreads(nw)  # ... override it for this run
            t = sc(idx, masks)[0]
            results[(seed, f"serial, OpenCV {nw} threads")] = digest(np.concatenate([t.ravel(), advantages_v4(t)]))
        cv2.setNumThreads(1)
        for workers, backend in ((1, "thread"), (nw, "process")):
            out = subprocess.run([sys.executable, "-c", SNIPPET % (str(ROOT), seed, workers, backend)],
                                 capture_output=True, text=True, check=True)
            results[(seed, f"fresh interpreter, {workers} {backend}")] = out.stdout.strip()
    for seed in (11, 12):
        ds = {lab: d for (s, lab), d in results.items() if s == seed}
        print(f"batch seed {seed} (256 rollouts, 32 prompts):")
        for lab, d in ds.items():
            print(f"  {lab:32s} {d[:16]}...")
        print(f"  identical across {len(ds)} settings: {len(set(ds.values())) == 1}")
    print("K10:", "PASS" if all(len({d for (s, _), d in results.items() if s == seed}) == 1 for seed in (11, 12)) else "FAIL")


if __name__ == "__main__":
    main()
