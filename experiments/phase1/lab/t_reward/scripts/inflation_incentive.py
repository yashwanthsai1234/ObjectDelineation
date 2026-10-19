"""Iteration 4: is there a learnable inflation (or deflation) incentive?

A policy cannot output the rater union without knowing the raters, but it CAN learn a systematic
size bias. Test: in honest groups (8 realistic rollouts: scripts/honest.py 'real' mode), replace
rollout 0 by itself dilated / eroded by d in {1 px, 0.25 tau, 0.5 tau} (others unchanged) and
measure the change of ITS advantage. A positive mean change for dilation = the reward pays for
inflating. 6 groups per item, all fixture items; 95 % CI over items.

    ../.venv/bin/python scripts/inflation_incentive.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from honest import honest_pools  # noqa: E402
from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.advantages import dm_loo_topology, reward_v3, reward_v4  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402

W, F = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])
SCHEMES = (("v3", None, None), ("v4 gate zero", "gate", 7), ("v4 term zero", "term", 7), ("v4 gate tauP", "gate", 8))
MOVES = (("dilate 1px", 1.0, None), ("dilate 0.25tau", None, 0.25), ("dilate 0.5tau", None, 0.5),
         ("erode 1px", -1.0, None), ("erode 0.25tau", None, -0.25), ("erode 0.5tau", None, -0.5))


def pack(t):
    return (t.iou, t.iou_tc, t.band["F5"], t.topo_any, t.topo_freq, signature_id(t.sig), t.area, t.g_area, t.g_area_tauP)


def adv(T, scheme):
    _, guard, gcol = scheme
    T = np.asarray(T, float)
    iou_tc, band, topo, freq, sid = T[..., 1], T[..., 2], T[..., 3], T[..., 4], T[..., 5]
    legal = np.ones_like(band)
    if guard is None:
        r = dm_loo_topology(sid, freq, (topo > 0) & (iou_tc >= 0.5))
        return reward_v3(np.stack([iou_tc, band, legal], -1), W, F, ("iou", "band", "legal"), r, c_topo=2.0)[0]
    g = T[..., gcol]
    r = dm_loo_topology(sid, freq, (topo > 0) & (iou_tc * (g if guard == "gate" else 1) >= 0.5))
    return reward_v4(np.stack([iou_tc, band, legal, g], -1), W, F, ("iou", "band", "legal", "g_area"), r, c_topo=2.0, guard=guard)[0]


def main():
    rng = np.random.default_rng(0)
    res = {(s[0], m[0]): {} for s in SCHEMES for m in MOVES}
    for _, key, st in fixtures.load_all():
        cache = ExpertCache.build(st)
        pools = honest_pools(st, list(cache.rater_sigs), cache, rng, "real", 6, lambda x: x, lambda x, c: x)
        for _ in range(6):
            masks = [pools[r][rng.integers(6)] for r in rng.integers(0, len(st), 8)]
            base = [pack(score(m, cache)) for m in masks]
            for mname, px, frac in MOVES:
                d = px if px is not None else frac * cache.tau
                moved = perturb.grow(masks[0], d)
                T1 = [[pack(score(moved, cache))] + base[1:]]
                for s in SCHEMES:
                    da = adv(T1, s)[0, 0] - adv([base], s)[0, 0]
                    res[(s[0], mname)].setdefault(key, []).append(da)
    print("mean change of the moved rollout's advantage (others fixed), mean over items [95% CI]; > 0 = incentive")
    print(f"  {'move':16s} " + " ".join(f"{s[0]:>26s}" for s in SCHEMES))
    for mname, _, _ in MOVES:
        cells = []
        for s in SCHEMES:
            per_item = np.array([np.mean(v) for v in res[(s[0], mname)].values()])
            boots = [np.mean(rng.choice(per_item, len(per_item))) for _ in range(2000)]
            lo, hi = np.percentile(boots, [2.5, 97.5])
            cells.append(f"{per_item.mean():+.3f} [{lo:+.3f},{hi:+.3f}]")
        print(f"  {mname:16s} " + " ".join(f"{c:>26s}" for c in cells))


if __name__ == "__main__":
    main()
