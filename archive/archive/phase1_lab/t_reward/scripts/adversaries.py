"""Iteration 3 (R3): set-specific adversaries scored inside honest GRPO groups.

Honest rollouts: rater masks (rater drawn uniformly) with smooth boundary jitter of amplitude
`amp` (1 px as specified; 0.5 tau as a stress setting), size-filtered topology preserved.
A group = 1 adversary + 7 honest rollouts; 256 groups per batch; 20 groups per item and adversary.
KPI: adversary advantage < median honest advantage of its group in >= 95 % of groups
(95 % CI: cluster bootstrap over items).

Adversaries
  dilated_union   rater union dilated by 2 tau
  hull            convex hull of the rater union
  holes_stamped   vote minus every rater's holes (only items where that differs from the vote)
  band_speck      vote + a detached compact speck of min_area px inside the band (U minus vote)
  minority_hole   vote with the minority rater's holes punched (its topology matches that rater);
                  must also score below that rater's own jittered mask (placed in the group)
  reversed ring   skipped: a mask has no orientation; winding is a polygon-legality matter (t_polycore)

    ../.venv/bin/python scripts/adversaries.py
"""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import cv2  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.advantages import dm_loo_topology, gdpo_floor_gated, reward_v2, reward_v3  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402
from rewardcore.topology import filtered_signature  # noqa: E402

G, N_GROUPS, N_POOL = 8, 20, 12
W1, F1 = np.array([0.3, 0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25, 0.25])
W2, F2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])
ADVERSARIES = ("dilated_union", "hull", "holes_stamped", "band_speck", "minority_hole")


def pack(t):
    return (t.iou, t.iou_tc, t.band["F5"], t.topo_any, t.topo_freq, signature_id(t.sig), t.area)


def hull(m):
    pts = cv2.findNonZero(m.astype(np.uint8))
    out = np.zeros(m.shape, np.uint8)
    cv2.fillConvexPoly(out, cv2.convexHull(pts), 1)
    return out.astype(bool) | m


def holes_of(m, min_area):
    h = ndimage.binary_fill_holes(m) & ~m
    lab, n = ndimage.label(h)
    if not n:
        return h
    sizes = np.bincount(lab.ravel())
    return (sizes >= min_area)[lab] & h


def band_speck(cache, rng):
    A = cache.min_area
    allowed = cache.union & ~ndimage.binary_dilation(cache.vote, np.ones((3, 3), bool), iterations=2)
    off = perturb._disk_offsets(A)
    cand = np.argwhere(allowed)
    rng.shuffle(cand)
    for c in cand[:400]:
        p = c + off
        if (p.min() >= 0 and p[:, 0].max() < allowed.shape[0] and p[:, 1].max() < allowed.shape[1]
                and allowed[p[:, 0], p[:, 1]].all()):
            out = cache.vote.copy()
            out[p[:, 0], p[:, 1]] = True
            return out
    return None


def item_data(args):
    ds, key, stack, amp_mode, seed = args[:5]
    cache_kw = args[5] if len(args) > 5 else {}
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack, **cache_kw)
    amp = 1.0 if amp_mode == "1px" else max(1.0, 0.5 * cache.tau)
    pools = []
    for j, m in enumerate(stack):
        pool = []
        for _ in range(N_POOL):
            for _try in range(20):
                x = perturb.jitter(m, amp, rng, corr=3.0 if amp_mode == "1px" else max(3.0, cache.tau))
                if filtered_signature(x, cache.min_area) == cache.rater_sigs[j]:
                    pool.append(pack(score(x, cache)))
                    break
        pools.append(pool or [pack(score(m, cache))])
    adv = {}
    adv["dilated_union"] = perturb.dilate(cache.union, 2 * cache.tau)
    adv["hull"] = hull(cache.union)
    holes = np.zeros_like(cache.vote)
    for m in stack:
        holes |= holes_of(m, cache.min_area)
    if (cache.vote & holes).any():
        adv["holes_stamped"] = cache.vote & ~holes
    sp = band_speck(cache, rng)
    if sp is not None:
        adv["band_speck"] = sp
    vote_sig = filtered_signature(cache.vote, cache.min_area)
    minority_rater = None
    for j, m in enumerate(stack):
        if cache.rater_sigs[j] != vote_sig and sum(cache.rater_sigs[j]) > sum(vote_sig):
            h = holes_of(m, cache.min_area)
            cand = cache.vote & ~h
            if (cache.vote & h).any() and filtered_signature(cand, cache.min_area) == cache.rater_sigs[j]:
                adv["minority_hole"] = cand
                minority_rater = j
                break
    return dict(key=key, subset=fixtures.subset(ds, key), k=len(stack), pools=pools,
                adv={a: pack(score(m, cache)) for a, m in adv.items()}, minority_rater=minority_rater,
                union_area=cache.n_union, max_rater_area=max(cache.rater_areas), gscale=cache.tol5 * cache.perimeter)


def advantages(T, scheme, c_topo=2.0, lam=0.5):
    T = np.asarray(T, float)
    iou, iou_tc, band, topo, freq, sid = (T[..., i] for i in range(6))
    legal = np.ones_like(iou)
    if scheme == "v1":
        return gdpo_floor_gated(np.stack([iou, band, topo, legal], -1), W1, floors=F1, names=("iou", "band", "topo", "legal"))[0]
    r_topo = dm_loo_topology(sid, freq, (topo > 0) & (iou_tc >= 0.5), lam=lam)
    fn = reward_v2 if scheme == "v2" else reward_v3  # v3+guard: T already guarded
    return fn(np.stack([iou_tc, band, legal], -1), W2, F2, ("iou", "band", "legal"), r_topo, c_topo=c_topo)[0]


SCHEMES = (("v1", 0), ("v3", 1.0), ("v3", 2.0), ("v3+guard", 2.0))
GUARD_KAPPA = 0.3


def guard(T, amax, scale):
    """Area-inflation guard (experiment): g = 1 / (1 + (x / kappa)^2), x = max(0, area - max rater
    area) / (T5 * P_U); multiplies IoU_tc and F5."""
    x = np.maximum(0.0, T[..., 6] - amax[:, None]) / scale[:, None]
    g = 1.0 / (1.0 + (x / GUARD_KAPPA) ** 2)
    T = T.copy()
    T[..., 1] *= g
    T[..., 2] *= g
    return T


def ci(ok, items, rng, n_boot=2000):
    keys = sorted(set(items))
    pos = {k: i for i, k in enumerate(keys)}
    s = np.zeros((len(keys), 2))
    for o, k in zip(ok, items):
        s[pos[k]] += (o, 1)
    b = [s[i, 0].sum() / s[i, 1].sum() for i in (rng.integers(0, len(s), len(s)) for _ in range(n_boot))]
    return np.percentile(b, [2.5, 97.5])


def load_items():
    """Fixture items + full-LIDC TRAIN items whose raters disagree under the size filter (for the
    hole adversaries, which need a minority hole)."""
    items = [(ds, k, st) for ds, k, st in fixtures.load_all()]
    for k, st in fixtures.load_lidc_full(split="train", disagree_only=True).items():
        items.append(("lidc", "full:" + k, st))
    return items


def run(amp_mode, rng, items):
    with Pool(3) as pool:
        data = pool.map(item_data, [(ds, k, st, amp_mode, j) for j, (ds, k, st) in enumerate(items)], chunksize=8)
    print(f"\n[honest rollouts: rater masks + {amp_mode} smooth jitter]")
    print(f"  {'adversary':14s} {'items':>5s} | " + " | ".join(f"{a + ('' if a == 'v1' else f' c={c:g}'):>26s}" for a, c in SCHEMES)
          + " | adversary mean terms (IoU_tc, F5, topo) vs honest")
    for name in ADVERSARIES:
        groups, items_of, subs, amax, gscale = [], [], [], [], []
        for d in data:
            if name not in d["adv"] or (name == "minority_hole" and d["minority_rater"] is None):
                continue
            if name in ("dilated_union", "hull", "band_speck") and d["key"].startswith("full:"):
                continue  # fixtures suffice for these
            for _ in range(N_GROUPS):
                draws = rng.integers(0, d["k"], G - 1)
                if name == "minority_hole":
                    draws[0] = d["minority_rater"]
                honest = [d["pools"][r][rng.integers(len(d["pools"][r]))] for r in draws]
                groups.append([d["adv"][name]] + honest)
                items_of.append(d["key"])
                subs.append(d["subset"])
                amax.append(d["max_rater_area"])
                gscale.append(d["gscale"])
        if not groups:
            print(f"  {name:14s} {0:5d} | n/a")
            continue
        T = np.array(groups, float)
        Tg = guard(T, np.array(amax, float), np.array(gscale, float))
        cells, per_sub = [], {}
        for scheme, c in SCHEMES:
            X = Tg if scheme == "v3+guard" else T
            a = np.concatenate([advantages(X[s:s + 256], scheme, c) for s in range(0, len(X), 256)])
            ok = a[:, 0] < np.median(a[:, 1:], axis=1)
            if name == "minority_hole":
                ok = ok & (a[:, 0] < a[:, 1])  # also below the minority rater's own jittered mask
            lo, hi = ci(ok, items_of, rng)
            cells.append(f"{ok.mean():5.1%} [{lo:5.1%},{hi:5.1%}]{'P' if lo >= 0.95 else ('p' if ok.mean() >= 0.95 else 'F')}")
            per_sub[(scheme, c)] = {sb: ok[[x == sb for x in subs]].mean() for sb in sorted(set(subs))}
        print(f"  {name:14s} {len(set(items_of)):5d} | " + " | ".join(f"{c:>26s}" for c in cells)
              + f" | adv ({T[:, 0, 1].mean():.3f}, {T[:, 0, 2].mean():.3f}, {T[:, 0, 3].mean():.2f})"
              f" honest ({T[:, 1:, 1].mean():.3f}, {T[:, 1:, 2].mean():.3f}, {T[:, 1:, 3].mean():.2f})")
        print("      per subset, v3 c=2: " + ", ".join(f"{sb} {v:.1%}" for sb, v in per_sub[("v3", 2.0)].items())
              + " | v3+guard: " + ", ".join(f"{sb} {v:.1%}" for sb, v in per_sub[("v3+guard", 2.0)].items()))
    return data


def main():
    rng = np.random.default_rng(0)
    print("KPI: share of groups with adversary advantage < median honest advantage, 95 % CI over items;"
          " PASS = lower CI >= 95 %, pass = point >= 95 %")
    print("reversed-winding ring: skipped (masks carry no orientation; covered by polygon legality in t_polycore)")
    items = load_items()
    print(f"items: {len(items)} (fixtures + full-LIDC train items with size-filtered rater disagreement)")
    for amp_mode in ("1px", "0.5tau"):
        run(amp_mode, rng, items)


if __name__ == "__main__":
    main()
