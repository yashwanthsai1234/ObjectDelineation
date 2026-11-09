"""Iteration 4: adversary suite with realistic honest rollouts and the area guard.

Honest rollouts (7 per group, raters drawn uniformly; scripts/honest.py: one amplitude per pool entry
shared by all raters of the item, so topology classes get the same perturbation quality):
  jit1px : rater mask + smooth 1-px boundary jitter (iteration 3; unrealistically poor on tiny nodules)
  real   : 25 % exact rater masks; otherwise smooth boundary perturbation with amplitude ~ U(0, tau)
           and correlation length max(3, tau) px (tau = the item's P95 rater-pair contour distance)
  (size-filtered topology preserved in both; 12 pool entries per item)
Adversaries: dilated_union (2 tau), hull, union, core (rater intersection: the deflation hedge),
band_speck, holes_stamped, minority_hole.
Schemes: v3 (no guard); v4 with the area guard g = min(1, A_hi/A) * min(1, A/A_lo) computed from ALL
raters, tolerance 'zero' (A_hi = A_max, A_lo = A_min) or 'tauP' (+- tau * P_U), applied as a 'gate'
(multiplies IoU_tc and F5) or as a separate z-scored 'term' (w 0.3, floor 0.05).
KPI (inflation / deflation / speck): adversary advantage < median honest advantage in >= 95 % of
groups, 95 % CI lower bound >= 95 % (cluster bootstrap over items).
Re-registered hole KPI: win-rate of the adversary against the median honest rollout of the SAME
topology class <= 60 %, CI upper bound <= 65 % (no systematic advantage). The original KPI (below
the overall median) is reported too.

    ../.venv/bin/python scripts/iter4_adversaries.py
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from adversaries import band_speck, ci, holes_of, hull, load_items  # noqa: E402
from honest import honest_pools  # noqa: E402
from rewardcore import perturb  # noqa: E402
from rewardcore.advantages import dm_loo_topology, reward_v3, reward_v4  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402
from rewardcore.topology import filtered_signature  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
G, N_GROUPS, N_POOL = 8, 20, 12
W, F = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])
NAMES = ("iou", "band", "legal", "g_area")
INFLATION = ("dilated_union", "hull", "union", "core", "band_speck")
HOLES = ("holes_stamped", "minority_hole")
SCHEMES = (("v3", None, None), ("v4 gate zero", "gate", "zero"), ("v4 term zero", "term", "zero"),
           ("v3 meanIoU", None, None), ("meanIoU gate zero", "gate", "zero"), ("meanIoU term zero", "term", "zero"))


def pack(t):
    return (t.iou, t.iou_tc, t.band["F5"], t.topo_any, t.topo_freq, signature_id(t.sig), t.area, t.g_area, t.g_area_tauP,
            t.iou_tc_mean)


def item_data(args):
    ds, key, stack, seed = args
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack)
    sigs = list(cache.rater_sigs)
    pools = {mode: honest_pools(stack, sigs, cache, rng, mode, N_POOL, pack, score) for mode in ("jit1px", "real")}
    adv = {"dilated_union": perturb.dilate(cache.union, 2 * cache.tau), "hull": hull(cache.union),
           "union": cache.union, "core": cache.core}
    if not cache.core.any():
        del adv["core"]
    sp = band_speck(cache, rng)
    if sp is not None:
        adv["band_speck"] = sp
    holes = np.zeros_like(cache.vote)
    for m in stack:
        holes |= holes_of(m, cache.min_area)
    if (cache.vote & holes).any():
        adv["holes_stamped"] = cache.vote & ~holes
    vote_sig = filtered_signature(cache.vote, cache.min_area)
    for j, m in enumerate(stack):
        if sigs[j] != vote_sig and sum(sigs[j]) > sum(vote_sig):
            h = holes_of(m, cache.min_area)
            cand = cache.vote & ~h
            if (cache.vote & h).any() and filtered_signature(cand, cache.min_area) == sigs[j]:
                adv["minority_hole"] = cand
                break
    packed = {a: pack(score(x, cache)) for a, x in adv.items()}
    adv_sig = {a: filtered_signature(x, cache.min_area) for a, x in adv.items()}
    same_class = {a: [j for j in range(len(stack)) if sigs[j] == adv_sig[a]] for a in adv}
    return dict(key=key, ds=ds, k=len(stack), pools=pools, adv=packed, same_class=same_class)


def advantages(T, scheme):
    name, guard, tol = scheme
    T = np.asarray(T, float)
    iou_tc, band, topo, freq, sid = T[..., 1], T[..., 2], T[..., 3], T[..., 4], T[..., 5]
    if "meanIoU" in name:
        iou_tc = T[..., 9]
    g = T[..., 7] if tol == "zero" else T[..., 8]
    legal = np.ones_like(band)
    gate_iou = iou_tc * (g if guard == "gate" else 1.0)
    r_topo = dm_loo_topology(sid, freq, (topo > 0) & (gate_iou >= 0.5), lam=0.5)
    if guard is None:
        return reward_v3(np.stack([iou_tc, band, legal], -1), W, F, ("iou", "band", "legal"), r_topo, c_topo=2.0)[0]
    return reward_v4(np.stack([iou_tc, band, legal, g], -1), W, F, NAMES, r_topo, c_topo=2.0, guard=guard)[0]


def main():
    rng = np.random.default_rng(0)
    items = load_items()
    with Pool(3) as pool:
        data = pool.map(item_data, [(ds, k, st, j) for j, (ds, k, st) in enumerate(items)], chunksize=8)
    print(f"items: {len(items)} (fixtures + full-LIDC train disagreement items); KPI cells: share [95% CI]")
    summary = {}
    for mode in ("jit1px", "real"):
        print(f"\n[honest rollouts: {mode}]  inflation/deflation/speck: adversary < median honest (PASS = CI low >= 95 %)")
        print(f"  {'adversary':14s} {'items':>5s} | " + " | ".join(f"{s[0]:>24s}" for s in SCHEMES)
              + " | adv (IoU_tc, F5, g) vs honest")
        for name in INFLATION + HOLES:
            T, it_of, same_idx = [], [], []
            for d in data:
                if name not in d["adv"]:
                    continue
                if d["ds"] == "lidc" and d["key"].startswith("full:") and name in INFLATION:
                    continue  # fixtures suffice for these
                sc = d["same_class"][name]
                if name in HOLES and not sc:
                    continue  # adversary topology no rater drew: invalid, trivially loses
                for _ in range(N_GROUPS):
                    draws = rng.integers(0, d["k"], G - 1)
                    if name in HOLES:
                        draws[0] = sc[int(rng.integers(len(sc)))]
                    T.append([d["adv"][name]] + [d["pools"][mode][r][rng.integers(N_POOL)] for r in draws])
                    it_of.append(d["key"])
                    same_idx.append([1 + g for g, r in enumerate(draws) if r in sc])
            if not T:
                continue
            T = np.array(T, float)
            cells = []
            for s in SCHEMES:
                a = np.concatenate([advantages(T[b:b + 256], s) for b in range(0, len(T), 256)])
                below = a[:, 0] < np.median(a[:, 1:], axis=1)
                lo, hi = ci(below, it_of, rng)
                if name in HOLES:
                    win = np.array([a[i, 0] > np.median(a[i, idx]) for i, idx in enumerate(same_idx)])
                    wlo, whi = ci(win, it_of, rng)
                    ok = win.mean() <= 0.60 and whi <= 0.65
                    cells.append(f"win {win.mean():4.0%} [{wlo:3.0%},{whi:3.0%}]{'P' if ok else 'F'}")
                    summary[(mode, name, s[0])] = dict(win=win.mean(), win_ci=(wlo, whi), below_median=below.mean())
                else:
                    cells.append(f"{below.mean():5.1%} [{lo:5.1%},{hi:5.1%}]{'P' if lo >= 0.95 else 'F'}")
                    summary[(mode, name, s[0])] = dict(below_median=below.mean(), ci=(lo, hi))
            honest = T[:, 1:]
            print(f"  {name:14s} {len(set(it_of)):5d} | " + " | ".join(f"{c:>24s}" for c in cells)
                  + f" | ({T[:, 0, 1].mean():.3f}, {T[:, 0, 2].mean():.3f}, {T[:, 0, 7].mean():.3f})"
                  f" vs ({honest[..., 1].mean():.3f}, {honest[..., 2].mean():.3f}, {honest[..., 7].mean():.3f})")
        print("  (hole rows: re-registered KPI = win-rate vs the median honest rollout of the same topology class,"
              " PASS if <= 60 % and CI upper <= 65 %)")
        for name in HOLES:
            cells = [f"{s[0]}: {summary[(mode, name, s[0])]['below_median']:.1%}" for s in SCHEMES if (mode, name, s[0]) in summary]
            print(f"  original KPI (adversary < overall median honest), {name}: " + ", ".join(cells))
    (OUT / "iter4_adversaries.json").write_text(json.dumps({"|".join(k): v for k, v in summary.items()}, default=float))


if __name__ == "__main__":
    main()
