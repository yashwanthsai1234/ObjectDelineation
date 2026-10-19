"""Iteration 3 (R3): K9 with the IoU_tc vote fallback, size-filtered topology, mixed groups, item CIs.

Pairs as in k9_discrimination.py (defects at/above the size filter; correct partner drawn at random
among topology-preserving variants with |IoU vs vote difference| <= 0.01). Each pair is placed in a
group of G = 8: [correct, wrong] + 6 fillers, of which U{0..3} are other wrong-topology variants of the
same item (so the wrong sample is not always the DM-loo minimum) and the rest topology-correct
variants. 256 groups per batch. Correct if advantage(correct) > advantage(wrong). 95 % CIs: cluster
bootstrap over items (2000 resamples).

    ../.venv/bin/python scripts/k9_v3.py
"""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from k9_discrimination import MAX_DIOU, composite, correct_candidates, wrong_variants  # noqa: E402
from rewardcore import fixtures  # noqa: E402
from rewardcore.advantages import dm_loo_topology, gdpo_floor_gated, reward_v3, reward_v4  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402

G = 8
W1, F1 = np.array([0.3, 0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25, 0.25])
W2, F2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])


def pack(t):
    return (t.iou, t.iou_tc, t.band["F5"], t.topo_any, t.topo_freq, signature_id(t.sig), t.g_area, t.g_area_tauP)


V4 = (("v4 gate zero", "gate", 6), ("v4 gate tauP", "gate", 7), ("v4 term zero", "term", 6), ("v4 term tauP", "term", 7))


def advantages_v4(groups, guard, gcol, c_topo=2.0):
    T = np.array([g["terms"] for g in groups], dtype=float)
    iou_tc, band, topo, freq, sid, g = T[..., 1], T[..., 2], T[..., 3], T[..., 4], T[..., 5], T[..., gcol]
    legal = np.ones_like(band)
    gate_iou = iou_tc * (g if guard == "gate" else 1.0)
    r_topo = dm_loo_topology(sid, freq, (topo > 0) & (gate_iou >= 0.5), lam=0.5)
    return reward_v4(np.stack([iou_tc, band, legal, g], -1), W2, F2, ("iou", "band", "legal", "g_area"), r_topo,
                     c_topo=c_topo, guard=guard)[0]


def item_groups(args):
    name, subset, stack, two_part, seed = args
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack)
    out = []
    for base in (cache.vote, stack[int(rng.integers(len(stack)))]):
        if score(base, cache).topo_any != 1.0:
            continue
        cands = [t for _, m in correct_candidates(base, rng) if (t := score(m, cache)).topo_any == 1.0]
        wrongs = [(w, t) for w, m in wrong_variants(base, cache, rng, two_part) if (t := score(m, cache)).topo_any == 0.0]
        if len(cands) < G:
            continue
        for wi, (wname, tw) in enumerate(wrongs):
            near = [j for j, c in enumerate(cands) if abs(c.iou - tw.iou) <= MAX_DIOU]
            if not near:
                continue
            jc = near[int(rng.integers(len(near)))]
            others_w = [t for j, (_, t) in enumerate(wrongs) if j != wi]
            n_w = min(len(others_w), int(rng.integers(0, 4)))
            fill_w = [others_w[j] for j in rng.choice(len(others_w), n_w, replace=False)] if n_w else []
            pool = [j for j in range(len(cands)) if j != jc]
            fill_c = [cands[j] for j in rng.choice(pool, G - 2 - n_w, replace=False)]
            out.append(dict(item=name, subset=subset, defect=wname.rstrip("0123456789"),
                            hi=min(cands[jc].iou, tw.iou) >= 0.5, terms=[pack(t) for t in [cands[jc], tw] + fill_w + fill_c]))
    return out


def advantages(groups, c_topo, lam=0.5):
    """v1 (iteration 1) and the frozen v3: reward_v3 + debiased, clipped DM-loo, global F5 band."""
    T = np.array([g["terms"] for g in groups], dtype=float)  # (B, G, 8)
    iou, iou_tc, band, topo, freq, sid = (T[..., i] for i in range(6))
    legal = np.ones_like(iou)
    a1 = gdpo_floor_gated(np.stack([iou, band, topo, legal], -1), W1, floors=F1, names=("iou", "band", "topo", "legal"))[0]
    r_topo = dm_loo_topology(sid, freq, (topo > 0) & (iou_tc >= 0.5), lam=lam)
    a2 = reward_v3(np.stack([iou_tc, band, legal], -1), W2, F2, ("iou", "band", "legal"), r_topo, c_topo=c_topo)[0]
    return a1, a2


def cluster_ci(ok, items, rng, n_boot=2000):
    keys = sorted(set(items))
    idx = {k: [] for k in keys}
    for j, k in enumerate(items):
        idx[k].append(j)
    per = [(ok[idx[k]].sum(), len(idx[k])) for k in keys]
    s = np.array(per, float)
    boots = []
    for _ in range(n_boot):
        b = s[rng.integers(0, len(s), len(s))]
        boots.append(b[:, 0].sum() / b[:, 1].sum())
    return np.percentile(boots, [2.5, 97.5])


def build_groups(n_proc=3, full_lidc=False):
    if full_lidc:  # R4: every full-LIDC train item with size-filtered rater topology disagreement
        args = [(k, "lidc-full/train-disagree", st, False, j)
                for j, (k, st) in enumerate(fixtures.load_lidc_full("train", disagree_only=True).items())]
        with Pool(n_proc) as pool:
            return [g for gs in pool.map(item_groups, args, chunksize=4) for g in gs]
    items = fixtures.load_all()
    args = [(k, fixtures.subset(ds, k), st, False, j) for j, (ds, k, st) in enumerate(items)]
    for sub in ("lidc/nodule", "riga/cup"):
        grp = [(k, st) for ds, k, st in items if fixtures.subset(ds, k) == sub]
        for j in range(0, len(grp) - 1, 2):
            (ka, a), (kb, b) = grp[j], grp[j + 1]
            args.append((f"{ka}+{kb}", sub + "+2part", composite(a, b, gap=3 + j % 6), True, 10_000 + j))
    with Pool(n_proc) as pool:
        return [g for gs in pool.map(item_groups, args, chunksize=4) for g in gs]


def main():
    full = "--full-lidc" in sys.argv
    groups = build_groups(full_lidc=full)
    print("item set:", "full-LIDC train disagreement items" if full else "fixtures + 2-part composites")
    rng = np.random.default_rng(0)
    order = rng.permutation(len(groups))
    hi = np.array([g["hi"] for g in groups])
    items = [g["item"] for g in groups]
    print(f"pairs: {len(groups)} from {len(set(items))} items; high-IoU {hi.sum()}; "
          f"groups with extra wrong fillers {np.mean([sum(t[3] == 0 for t in g['terms'][2:]) > 0 for g in groups]):.0%}")
    iou_ok = np.array([g["terms"][0][0] > g["terms"][1][0] for g in groups])
    lo, hi_ci = cluster_ci(iou_ok, items, rng)
    print(f"IoU (vote) alone correct: {iou_ok.mean():.1%} [{lo:.1%}, {hi_ci:.1%}]")
    for name, guard, gcol in V4:
        ok = np.zeros(len(groups), bool)
        for s in range(0, len(groups), 256):
            idx = order[s:s + 256]
            a = advantages_v4([groups[i] for i in idx], guard, gcol)
            ok[idx] = a[:, 0] > a[:, 1]
        lo, up = cluster_ci(ok, items, rng)
        hlo, hup = cluster_ci(ok[hi], [i for i, h in zip(items, hi) if h], rng)
        print(f"  {name:22s} all {ok.mean():.2%} [{lo:.2%}, {up:.2%}] | IoU>=0.5 {ok[hi].mean():.2%} [{hlo:.2%}, {hup:.2%}]"
              f" | IoU<0.5 {ok[~hi].mean():.1%}")
    for c_topo in (2.0,) if "--v4-only" in sys.argv else (1.0, 1.5, 2.0, 3.0):
        ok1, ok2 = np.zeros(len(groups), bool), np.zeros(len(groups), bool)
        for s in range(0, len(groups), 256):
            idx = order[s:s + 256]
            a1, a2 = advantages([groups[i] for i in idx], c_topo)
            ok1[idx], ok2[idx] = a1[:, 0] > a1[:, 1], a2[:, 0] > a2[:, 1]
        for name, ok in ((f"v3 c_topo={c_topo}", ok2),) + ((("v1 gdpo_floor_gated", ok1),) if c_topo == 1.0 else ()):
            lo, up = cluster_ci(ok, items, rng)
            hlo, hup = cluster_ci(ok[hi], [i for i, h in zip(items, hi) if h], rng)
            cells = " ".join(f"{d} {ok[[g['defect'] == d for g in groups]].mean():.1%}" for d in ("pinhole", "speck", "bridge"))
            print(f"  {name:22s} all {ok.mean():.2%} [{lo:.2%}, {up:.2%}] | IoU>=0.5 {ok[hi].mean():.2%} [{hlo:.2%}, {hup:.2%}]"
                  f" | IoU<0.5 {ok[~hi].mean():.1%} (n={(~hi).sum()}) | {cells}")
        subs = sorted({g["subset"] for g in groups})
        print("    per subset (v3): " + ", ".join(f"{s} {ok2[[g['subset'] == s for g in groups]].mean():.1%}" for s in subs))


if __name__ == "__main__":
    main()
