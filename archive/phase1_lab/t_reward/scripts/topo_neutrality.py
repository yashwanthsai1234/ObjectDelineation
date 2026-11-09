"""Iteration 2, task 2c: is any topology class systematically favoured when the policy equals the raters?

For every fixture item whose raters disagree on per-instance topology (14 LIDC items), groups of
G = 8 are drawn i.i.d. from the rater set:
  "exact"  : the raters' own masks;
  "jitter" : a random topology-preserving 1-px smooth jitter of the drawn rater's mask (a model
             sample that follows that rater).
One group per item forms a batch (batch-level normalization), 400 batches. Class of a sample:
"mode" = the most frequent rater topology of the item, "minority" = any other rater topology.
Reported: mean advantage per class (items weighted equally) with a 95 % bootstrap CI over items,
and the mean of the per-item (minority - mode) difference.

Schemes:
  v1            gdpo_floor_gated on [IoU vs vote, band F5, topo_any, legal]       (iteration 1)
  v2 c_topo=c   reward_v2: geo [IoU_topocond, F5, legal] + c * DM-loo topology, mean-centred
  v2 geo only   reward_v2 without the topology term (c_topo = 0)
  v2 vote-IoU   reward_v2 with IoU vs vote instead of IoU_topocond
  v2 IoU_tc only / band only   decomposition of the geometric part

    ../.venv/bin/python scripts/topo_neutrality.py
"""
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.advantages import dm_loo_topology, gdpo_floor_gated, reward_v2  # noqa: E402
from rewardcore.reward import ExpertCache, score  # noqa: E402
from rewardcore.topology import instance_topology_full  # noqa: E402

G, N_BATCH, N_JITTER = 8, 400, 12
W1, F1 = np.array([0.3, 0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25, 0.25])
W2, F2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])


def item_pools(rng):
    items = []
    for _, key, st in fixtures.load_all():
        sigs = [instance_topology_full(m).signature() for m in st]
        if len(set(sigs)) == 1:
            continue
        cache = ExpertCache.build(st)
        mode = Counter(sigs).most_common(1)[0][0]
        pools = {"exact": [], "jitter": []}
        for j, m in enumerate(st):
            pools["exact"].append([(score(m, cache), sigs[j] != mode)])
            jit = []
            for _ in range(N_JITTER):
                for _try in range(30):
                    x = perturb.jitter(m, 1.0, rng, corr=3.0)
                    if instance_topology_full(x).signature() == sigs[j]:
                        jit.append((score(x, cache), sigs[j] != mode))
                        break
            pools["jitter"].append(jit or pools["exact"][-1])
        items.append(dict(key=key, k=len(st), pools=pools, sig_index={s: i for i, s in enumerate(sorted(set(sigs)))}))
    return items


def advantages(terms, sig_index):
    """terms: (B, G) list of Terms. Returns dict scheme -> (B, G) advantages."""
    B = len(terms)
    get = lambda f: np.array([[f(t) for t in row] for row in terms])  # noqa: E731
    iou_v, iou_tc, band = get(lambda t: t.iou), get(lambda t: t.iou_tc), get(lambda t: t.band["F5"])
    topo, legal = get(lambda t: t.topo_any), np.ones((B, G))
    freq = get(lambda t: t.topo_freq)
    sig_ids = np.array([[sig_index[b].get(t.sig, -1) for t in row] for b, row in enumerate(terms)])
    out = {}
    out["v1"] = gdpo_floor_gated(np.stack([iou_v, band, topo, legal], -1), W1, floors=F1,
                                 names=("iou", "band", "topo", "legal"))[0]
    names2 = ("iou", "band", "legal")
    for name, iou, c, w in (("v2 c_topo=1", iou_tc, 1.0, W2), ("v2 c_topo=2", iou_tc, 2.0, W2),
                            ("v2 c_topo=3", iou_tc, 3.0, W2), ("v2 geo only", iou_tc, 0.0, W2),
                            ("v2 IoU_tc only", iou_tc, 0.0, np.array([0.3, 0.0, 0.1])),
                            ("v2 band only", iou_tc, 0.0, np.array([0.0, 0.3, 0.1])),
                            ("v2 vote-IoU", iou_v, 1.0, W2)):
        r_topo = dm_loo_topology(sig_ids, freq, (topo > 0) & (iou >= 0.5), lam=0.5)
        out[name] = reward_v2(np.stack([iou, band, legal], -1), w, F2, names2, r_topo, c_topo=c)[0]
    return out


def main():
    rng = np.random.default_rng(0)
    items = item_pools(rng)
    print(f"items: {len(items)} (LIDC; RIGA has no topology disagreement); G={G}, {N_BATCH} batches of one group per item")
    for mode in ("exact", "jitter"):
        acc = {}  # scheme -> item -> class -> list of advantages
        for _ in range(N_BATCH):
            terms, minority, sig_index = [], [], []
            for it in items:
                draws = rng.integers(0, it["k"], G)
                row = [it["pools"][mode][d][rng.integers(len(it["pools"][mode][d]))] for d in draws]
                terms.append([t for t, _ in row])
                minority.append([mi for _, mi in row])
                sig_index.append(it["sig_index"])
            for scheme, a in advantages(terms, sig_index).items():
                for b, it in enumerate(items):
                    for g in range(G):
                        acc.setdefault(scheme, {}).setdefault(it["key"], {0: [], 1: []})[int(minority[b][g])].append(a[b, g])
        print(f"\n[{mode} rater masks] mean advantage per class (items equally weighted), 95% CI over items")
        print(f"  {'scheme':14s} {'mode class':>24s} {'minority class':>24s} {'minority - mode':>26s}")
        for scheme, per_item in acc.items():
            m0 = np.array([np.mean(v[0]) for v in per_item.values()])
            m1 = np.array([np.mean(v[1]) for v in per_item.values()])
            d = m1 - m0
            ci = lambda x: np.percentile([np.mean(rng.choice(x, len(x))) for _ in range(4000)], [2.5, 97.5])  # noqa: E731
            c0, c1, cd = ci(m0), ci(m1), ci(d)
            print(f"  {scheme:14s} {m0.mean():+.3f} [{c0[0]:+.3f},{c0[1]:+.3f}]   {m1.mean():+.3f} [{c1[0]:+.3f},{c1[1]:+.3f}]"
                  f"   {d.mean():+.3f} [{cd[0]:+.3f},{cd[1]:+.3f}]  (items with minority < mode: {np.sum(d < 0)}/{len(d)})")


if __name__ == "__main__":
    main()
