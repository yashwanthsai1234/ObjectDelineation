"""Iteration 3 (R4): disagreement-scale re-validation on every full-LIDC TRAIN item whose raters
disagree on the size-filtered per-instance topology (251 items).

1. topo_gap (production emulation): new samples = 1-px smooth jitter of each rater's mask (size-
   filtered topology preserved, 10 per rater; scripts/honest.py, amplitude matched across raters)
   scored against references from all raters. Gap per
   item = mean term of samples whose topology matches the vote's - mean of the others; 95 % CI over
   items. Terms: IoU vs vote, IoU_tc (reward), band F5 (topology-conditional, reward).
2. topo_neutrality: groups of G = 8 drawn i.i.d. from the raters ("exact": their masks; "jitter":
   a random pool sample of the drawn rater); 300 batches of one group per item. Mode class = most
   frequent rater topology (ties: the vote's topology, else the first in sorted order); minority =
   the others. Per item: mean advantage(minority) - mean advantage(mode); mean over items with a
   95 % bootstrap CI. Acceptance (jitter): CI half-width <= 0.10 and the CI contains 0.
   Grid: reward v3 with c_topo in {1, 1.5, 2, 3} x lambda in {0.25, 0.5, 1}; v1 for reference.

    ../.venv/bin/python scripts/r4_scale.py
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from rewardcore import fixtures  # noqa: E402
from rewardcore.advantages import dm_loo_topology, gdpo_floor_gated, reward_v3, reward_v4  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402
from rewardcore.topology import filtered_signature  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from honest import honest_pools  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
G, N_BATCH, N_JIT = 8, 300, 10
W1, F1 = np.array([0.3, 0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25, 0.25])
W2, F2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])
C_TOPO, LAMS = (1.0, 1.5, 2.0, 3.0), (0.25, 0.5, 1.0)


def pack(t, t_global):
    """(IoU vote, IoU_tc, F5 topology-conditional band, topo_any, f, signature id, F5 global band)."""
    return (t.iou, t.iou_tc, t.band["F5"], t.topo_any, t.topo_freq, signature_id(t.sig), t_global.band["F5"],
            t_global.g_area, t_global.g_area_tauP, t_global.iou_tc_mean)


def build_items(rng):
    items = []
    for key, st in fixtures.load_lidc_full("train", disagree_only=True).items():
        cache = ExpertCache.build(st, topology_conditional_band=True)
        cache_g = ExpertCache.build(st, topology_conditional_band=False)
        sigs = list(cache.rater_sigs)
        vsig = filtered_signature(cache.vote, cache.min_area)
        cnt = Counter(sigs)
        top = max(cnt.values())
        tied = sorted(s for s, c in cnt.items() if c == top)
        mode = vsig if vsig in tied else tied[0]
        def score_pair(x, _cache, cache=cache, cache_g=cache_g):
            """Packed terms under the topology-conditional band (reference) and the global band."""
            return pack(score(x, cache), score(x, cache_g))

        def identity(t):
            return t

        exact = [[score_pair(m, None)] for m in st]
        jit = honest_pools(st, sigs, cache, rng, "jit1px", N_JIT, identity, score_pair)
        real = honest_pools(st, sigs, cache, rng, "real", N_JIT, identity, score_pair)
        items.append(dict(key=key, k=len(st), minority=[s != mode for s in sigs], match_vote=[s == vsig for s in sigs],
                          exact=exact, jitter=jit, real=real))
    return items


def topo_gap(items, rng):
    print("1. per-topology gap, production emulation (1-px jittered rater masks vs all raters);"
          " gap = matches vote topology - not, mean over items [95% CI]")
    for i, name in ((0, "IoU vs vote (v1)"), (1, "IoU_tc (reward)"), (2, "band F5, topology-conditional (reward)")):
        gaps = []
        for it in items:
            m1 = [t[i] for j, pool in enumerate(it["jitter"]) if it["match_vote"][j] for t in pool]
            m0 = [t[i] for j, pool in enumerate(it["jitter"]) if not it["match_vote"][j] for t in pool]
            if m1 and m0:
                gaps.append(np.mean(m1) - np.mean(m0))
        g = np.array(gaps)
        lo, hi = np.percentile([np.mean(rng.choice(g, len(g))) for _ in range(4000)], [2.5, 97.5])
        print(f"   {name:42s} items {len(g)}  mean {g.mean():+.4f} [{lo:+.4f}, {hi:+.4f}]  median {np.median(g):+.4f}"
              f"  |gap| median {np.median(np.abs(g)):.4f}")


def advantages(T, scheme, c_topo=2.0, lam=0.5, debias=True, global_band=False):
    T = np.asarray(T, float)
    iou, iou_tc, band, topo, freq, sid = (T[..., i] for i in range(6))
    if global_band:
        band = T[..., 6]
    if scheme.startswith("m"):  # m3 / m4 ...: mean-over-class IoU instead of max (iteration 4)
        iou_tc = T[..., 9]
        scheme = "v" + scheme[1:]
    if scheme.startswith("v4"):
        _, guard, tol = scheme.split()
        g = T[..., 7] if tol == "zero" else T[..., 8]
        gate_iou = iou_tc * (g if guard == "gate" else 1.0)
        r_topo = dm_loo_topology(sid, freq, (topo > 0) & (gate_iou >= 0.5), lam=lam, debias=debias)
        return reward_v4(np.stack([iou_tc, band, np.ones_like(iou), g], -1), W2, F2, ("iou", "band", "legal", "g_area"),
                         r_topo, c_topo=c_topo, guard=guard)[0]
    legal = np.ones_like(iou)
    if scheme == "v1":
        return gdpo_floor_gated(np.stack([iou, band, topo, legal], -1), W1, floors=F1, names=("iou", "band", "topo", "legal"))[0]
    r_topo = dm_loo_topology(sid, freq, (topo > 0) & (iou_tc >= 0.5), lam=lam, debias=debias)
    return reward_v3(np.stack([iou_tc, band, legal], -1), W2, F2, ("iou", "band", "legal"), r_topo, c_topo=c_topo)[0]


def neutrality(items, rng):
    # frozen design: global band (topology_conditional_band=False); class band reported for reference
    schemes = ([("v1", 0, 0, True, True)] + [("v3", c, lam, True, True) for c in C_TOPO for lam in LAMS]
               + [(f"v4 {gd} {tl}", 2.0, 0.5, True, True) for gd in ("gate", "term") for tl in ("zero", "tauP")]
               + [("m3", 2.0, 0.5, True, True), ("m4 gate zero", 2.0, 0.5, True, True), ("m4 term zero", 2.0, 0.5, True, True)]
               + [("v3 no-debias", 2.0, 0.5, False, True), ("v3 class band", 2.0, 0.5, True, False),
                  ("geo only", 0.0, 0.5, True, True), ("geo only, class band", 0.0, 0.5, True, False)])
    res = {}
    for mode in ("exact", "jitter", "real"):
        acc = {s: [[[], []] for _ in items] for s in schemes}
        for _ in range(N_BATCH):
            T, cls = [], []
            for it in items:
                d = rng.integers(0, it["k"], G)
                T.append([it[mode][r][rng.integers(len(it[mode][r]))] for r in d])
                cls.append([it["minority"][r] for r in d])
            T = np.array(T, float)
            for s in schemes:
                a = advantages(T, s[0] if s[0][:2] in ("v4", "m3", "m4") else ("v1" if s[0] == "v1" else "v3"),
                               s[1], s[2], s[3], s[4])
                for b in range(len(items)):
                    for g in range(G):
                        acc[s][b][int(cls[b][g])].append(a[b, g])
        print(f"\n2. neutrality [{mode} rater masks], minority - mode advantage, mean over items [95% CI], half-width")
        for s in schemes:
            d = np.array([np.mean(x[1]) - np.mean(x[0]) for x in acc[s] if x[0] and x[1]])
            boots = [np.mean(rng.choice(d, len(d))) for _ in range(4000)]
            lo, hi = np.percentile(boots, [2.5, 97.5])
            ok = (hi - lo) / 2 <= 0.10 and lo <= 0 <= hi
            res[(mode,) + s] = (d.mean(), lo, hi)
            name = s[0] if s[0] == "v1" else f"{s[0]} c={s[1]:g} lam={s[2]:g}"
            print(f"   {name:34s} n={len(d)}  {d.mean():+.3f} [{lo:+.3f}, {hi:+.3f}]  hw {(hi - lo) / 2:.3f}"
                  f"  minority loses on {np.mean(d < 0):.0%} of items  {'ACCEPT' if ok else 'reject'}")
    return res


def main():
    rng = np.random.default_rng(0)
    items = build_items(rng)
    print(f"items: {len(items)} full-LIDC train items with size-filtered rater topology disagreement;"
          f" k: {dict(Counter(it['k'] for it in items))}; minority share of raters {np.mean([np.mean(it['minority']) for it in items]):.2f}\n")
    topo_gap(items, rng)
    res = neutrality(items, rng)
    (OUT / "iter3_r4_neutrality.json").write_text(json.dumps({"|".join(map(str, k)): v for k, v in res.items()}))


if __name__ == "__main__":
    main()
