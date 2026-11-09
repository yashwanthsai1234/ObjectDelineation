"""Iteration 2: K9 under reward v2 (group-level: the DM-loo topology term depends on the group).

Same matched pairs as k9_discrimination.py (|dIoU vs vote| <= 0.01, random correct partner).
Each pair is placed in a group of G = 8: [correct, wrong, 6 other topology-correct variants of the
same item]; groups are batched 256 at a time. A pair is ranked correctly if
advantage(correct) > advantage(wrong). Reported for v1 (gdpo_floor_gated on IoU-vote, F5,
topo_any, legal) and v2 (IoU_topocond + F5 z-scored with floors, DM-loo topology mean-centred).

    ../.venv/bin/python scripts/k9_v2.py
"""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from k9_discrimination import MAX_DIOU, composite, correct_candidates, wrong_variants  # noqa: E402
from rewardcore import fixtures  # noqa: E402
from rewardcore.advantages import dm_loo_topology, gdpo_floor_gated, reward_v2  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402
from rewardcore.topology import instance_topology_full  # noqa: E402

G = 8
W1, F1 = np.array([0.3, 0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25, 0.25])
W2, F2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])
FIELDS = ("iou", "iou_tc", "band", "topo_any", "topo_freq", "sig")


def pack(t):
    return (t.iou, t.iou_tc, t.band["F5"], t.topo_any, t.topo_freq, t.sig)


def item_groups(args):
    name, subset, stack, two_part, seed = args
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack)
    out = []
    for base in (cache.vote, stack[int(rng.integers(len(stack)))]):
        if instance_topology_full(base).signature() not in cache.signatures:
            continue
        cands = [t for _, m in correct_candidates(base, rng) if (t := score(m, cache)).topo_any == 1.0]
        if len(cands) < G:
            continue
        for wname, m in wrong_variants(base, cache, rng, two_part):
            tw = score(m, cache)
            if tw.topo_any == 1.0:
                continue
            near = [j for j, c in enumerate(cands) if abs(c.iou - tw.iou) <= MAX_DIOU]
            if not near:
                continue
            jc = near[int(rng.integers(len(near)))]
            fill = [cands[j] for j in rng.choice([j for j in range(len(cands)) if j != jc], G - 2, replace=False)]
            group = [cands[jc], tw] + fill
            out.append(dict(subset=subset, defect=wname.rstrip("0123456789"), hi=min(cands[jc].iou, tw.iou) >= 0.5,
                            terms=[pack(t) for t in group]))
    return out


def advantages(groups, c_topo=1.0):
    T = [g["terms"] for g in groups]
    get = lambda i: np.array([[t[i] for t in row] for row in T], dtype=float)  # noqa: E731
    iou, iou_tc, band, topo, freq = (get(i) for i in range(5))
    legal = np.ones_like(iou)
    sig_ids = np.array([[signature_id(t[5]) for t in row] for row in T])
    a1 = gdpo_floor_gated(np.stack([iou, band, topo, legal], -1), W1, floors=F1, names=("iou", "band", "topo", "legal"))[0]
    r_topo = dm_loo_topology(sig_ids, freq, (topo > 0) & (iou_tc >= 0.5), lam=0.5)
    a2 = reward_v2(np.stack([iou_tc, band, legal], -1), W2, F2, ("iou", "band", "legal"), r_topo, c_topo=c_topo)[0]
    return a1, a2


def main():
    items = fixtures.load_all()
    args = [(k, fixtures.subset(ds, k), st, False, j) for j, (ds, k, st) in enumerate(items)]
    for sub in ("lidc/nodule", "riga/cup"):
        grp = [(k, st) for ds, k, st in items if fixtures.subset(ds, k) == sub]
        for j in range(0, len(grp) - 1, 2):
            (ka, a), (kb, b) = grp[j], grp[j + 1]
            args.append((f"{ka}+{kb}", sub + "+2part", composite(a, b, gap=3 + j % 6), True, 10_000 + j))
    with Pool(4) as pool:  # moderate load
        groups = [g for gs in pool.map(item_groups, args, chunksize=4) for g in gs]
    rng = np.random.default_rng(0)
    order = rng.permutation(len(groups))
    hi = np.array([g["hi"] for g in groups])
    dtc = np.array([g["terms"][0][1] - g["terms"][1][1] for g in groups])  # IoU_tc(correct) - IoU_tc(wrong)
    print(f"IoU_tc(correct) - IoU_tc(wrong) on vote-IoU-matched pairs: median {np.median(dtc):+.3f}, "
          f"p5 {np.percentile(dtc, 5):+.3f}, share |d| <= 0.01: {np.mean(np.abs(dtc) <= 0.01):.1%}")
    tc_matched = np.abs(dtc) <= MAX_DIOU
    for c_topo in (1.0, 2.0, 3.0):
        ok1, ok2 = np.zeros(len(groups), bool), np.zeros(len(groups), bool)
        for s in range(0, len(groups), 256):
            idx = order[s:s + 256]
            a1, a2 = advantages([groups[i] for i in idx], c_topo)
            ok1[idx], ok2[idx] = a1[:, 0] > a1[:, 1], a2[:, 0] > a2[:, 1]
        print(f"c_topo={c_topo}: v2 all {ok2.mean():.2%} | IoU>=0.5 {ok2[hi].mean():.2%} | matched on IoU_tc "
              f"(n={tc_matched.sum()}) {ok2[tc_matched].mean():.2%}, high-IoU {ok2[tc_matched & hi].mean():.2%}"
              f" | v1 on IoU_tc-matched {ok1[tc_matched].mean():.2%}")
    ok1, ok2 = np.zeros(len(groups), bool), np.zeros(len(groups), bool)
    for s in range(0, len(groups), 256):
        idx = order[s:s + 256]
        a1, a2 = advantages([groups[i] for i in idx])
        ok1[idx], ok2[idx] = a1[:, 0] > a1[:, 1], a2[:, 0] > a2[:, 1]
    print("\ndetail at c_topo = 1 (vote-IoU-matched pairs):")
    iou_ok = np.array([g["terms"][0][0] > g["terms"][1][0] for g in groups])
    print(f"groups (pairs): {len(groups)}; high-IoU (both >= 0.5): {hi.sum()}; IoU-vote alone correct {iou_ok.mean():.1%}")
    for name, ok in (("v1 gdpo_floor_gated", ok1), ("v2 (IoU_tc + F5 + DM-loo)", ok2)):
        cells = " ".join(f"{d} {ok[[g['defect'] == d for g in groups]].mean():.1%}" for d in ("pinhole", "speck", "bridge"))
        print(f"  {name:28s} all {ok.mean():.2%} | IoU>=0.5 {ok[hi].mean():.2%} | IoU<0.5 {ok[~hi].mean():.1%} (n={(~hi).sum()})"
              f" | hard (IoU favours wrong, high-IoU) {ok[hi & ~iou_ok].mean():.2%} | {cells}")
    for sub in sorted({g["subset"] for g in groups}):
        m = np.array([g["subset"] == sub for g in groups])
        print(f"    {sub:20s} v1 {ok1[m].mean():.1%}  v2 {ok2[m].mean():.1%}  (n={m.sum()})")


if __name__ == "__main__":
    main()
