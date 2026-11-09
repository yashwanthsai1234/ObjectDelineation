"""Record the phase 2 lab outputs on seeded synthetic inputs into tests/goldens/p2_lab_refs.json (no data is read).

Archived with the lab. Its imports (`experiments.phase2.lab.B2`, `.B3`) work only in a checkout of the tag
pre-clean-2026-11-08, where it ran as `python -I -m experiments.phase2.lab.record_p2_lab_refs` (needs torch)."""

import json
from pathlib import Path

import numpy as np

from experiments.phase2.lab.B2 import analyze as A
from experiments.phase2.lab.B2 import posthoc as PH
from experiments.phase2.lab.B2.score import score_view
from experiments.phase2.lab.B3 import power_c2
from objdel.eval import metrics
from objdel.eval import signal as S
from objdel.reward import ExpertCache
from objdel.targets import build_target

OUT = Path(__file__).resolve().parents[3] / "tests" / "goldens" / "p2_lab_refs.json"


def rand_groups(rng, n=40, k=8):  # same draws as tests/test_p2_signal.py::rand_groups
    out = []
    for i in range(n):
        legal = (rng.random(k) < rng.random()).astype(float)
        R = rng.random(k) * legal
        out.append({"item_id": f"lidc/1.2.{i // 3}/n000/z{i:03d}", "R": R.tolist(), "legal": legal.tolist(),
                    "iou_tc": (R * 1.5).clip(0, 1).tolist(), "advantage": (R - R.mean()).tolist(),
                    "sigma": float(R.std() > 0.05), "dup": 0.0, "distinct": 1.0})  # fmt: skip
    return out


def main():
    g = rand_groups(np.random.default_rng(3))
    v = {k: [x[k] for x in g] for k in g[0]}
    group = {k: np.asarray(x).tolist() for k, x in (A.group_metrics(v) | PH.parts(v)).items()}

    yy, xx = np.ogrid[:50, :50]
    disk = lambda cy, cx, r: (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r  # noqa: E731
    rat = np.stack([disk(25, 25, 12), disk(25, 27, 14), disk(24, 25, 11)])
    cache = ExpertCache.build(rat)
    texts = [build_target({"nodule": disk(25, 25, r)}, ("nodule",), min_area={"nodule": cache.min_area}).text
             for r in (9, 11, 12, 13, 15)]  # fmt: skip
    samples = texts + ["garbage", texts[0], "{"]
    item = {"item_id": "lidc/x/n000/z001", "caches": {"nodule": cache}, "W": 50, "H": 50, "names": ("nodule",)}
    raw = score_view([{"item_id": item["item_id"], "samples": samples, "truncated": [False] * 8}], [item])
    cols = ("R", "legal", "iou_tc", "advantage", "sigma", "dup", "distinct")
    # the live scorer's row 0 on the same item must equal the lab's before it is frozen
    m = metrics(samples, {"nodule": rat}, 50, 50, ("nodule",), [False] * 8, {"nodule": cache})
    mine = S.group(item["item_id"], samples, m)
    assert all(mine[k] == raw[k][0] for k in cols)

    rng = np.random.default_rng(11)
    bins = [rng.random(120) for _ in range(3)]
    rng = np.random.default_rng(1)
    per_bin = [rng.random(120), rng.random(120) * 0.5, (rng.random(120) > 0.3).astype(float)]
    gold = {
        "seeds": {"rand_groups": 3, "c2_bins": [11, 0], "probes_bins": [1, 5]},
        "group": group,
        "score_view_row0": {k: np.asarray(raw[k][0]).tolist() for k in cols},
        "c2_pooled": [float(x) for x in power_c2.pooled(bins, np.random.default_rng(0))],
        "probes_pooled": [float(x) for x in power_c2.pooled(per_bin, np.random.default_rng(5))],
        "W": power_c2.W.tolist(),
    }
    OUT.write_text(json.dumps(gold, indent=1, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
