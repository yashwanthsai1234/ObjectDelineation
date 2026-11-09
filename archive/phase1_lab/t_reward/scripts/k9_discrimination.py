"""Task 3, KPI K9: does the full reward rank a topology-correct mask above a topology-wrong one of
(almost) the same IoU?

For every fixture item (and composite two-part items made by placing two fixture items side by
side), from a base mask (the strict-majority vote, or one rater's mask):
  wrong   = base + one topology defect: pinhole (1/4/9 px), detached speck (1/4/9 px, 2-12 px away),
            bridge between the two parts (composites, width 1/2 px). Must not match any rater's
            per-instance topology.
  correct = a topology-preserving perturbation (Euclidean grow/shrink by 0.5-3 px, shifts up to
            2 px, smooth boundary jitter) drawn at random among those with |IoU - IoU(wrong)| <= 0.01;
            its per-instance topology must match some rater. No such candidate -> no pair.
KPI: share of pairs with reward(correct) > reward(wrong) >= 95 %; IoU alone is reported (~50 %).

    ../.venv/bin/python scripts/k9_discrimination.py
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.reward import ExpertCache, RewardConfig, reward, score  # noqa: E402
from rewardcore.topology import bbox, instance_topology_full  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
MAX_DIOU = 0.01
CONFIGS = {
    "default (0.3 IoU + 0.3 band F5 + 0.4 topo_any)": RewardConfig(band_formula="F5"),
    "w_topo 0.2": RewardConfig(w_iou=0.4, w_band=0.4, w_topo=0.2, band_formula="F5"),
    "w_topo 0.1": RewardConfig(w_iou=0.45, w_band=0.45, w_topo=0.1, band_formula="F5"),
    "w_topo 0.05": RewardConfig(w_iou=0.475, w_band=0.475, w_topo=0.05, band_formula="F5"),
    "topo_near (0.4)": RewardConfig(band_formula="F5", topo_term="near"),
    "topo_freq (0.4)": RewardConfig(band_formula="F5", topo_term="freq"),
    "band F0 current, 0.4 topo_any": RewardConfig(band_formula="F0"),
    "repo default (F0 band, topo_near)": RewardConfig(band_formula="F0", topo_term="near"),
}


def composite(a: np.ndarray, b: np.ndarray, gap: int) -> np.ndarray:
    """Raters paired by index: part a on the left, part b on the right, `gap` px between union boxes."""
    k = min(len(a), len(b))
    ya0, ya1, xa0, xa1 = bbox(a.any(0))
    yb0, yb1, xb0, xb1 = bbox(b.any(0))
    a, b = a[:k, ya0:ya1, xa0:xa1], b[:k, yb0:yb1, xb0:xb1]
    pad = 12
    h = max(a.shape[1], b.shape[1]) + 2 * pad
    w = a.shape[2] + gap + b.shape[2] + 2 * pad
    out = np.zeros((k, h, w), bool)
    out[:, pad:pad + a.shape[1], pad:pad + a.shape[2]] = a
    out[:, pad:pad + b.shape[1], pad + a.shape[2] + gap:pad + a.shape[2] + gap + b.shape[2]] = b
    return out


def correct_candidates(base, rng):
    for d in (-3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3):
        yield f"grow{d:+}", perturb.grow(base, d)
    for dy in (-2, -1, 0, 1, 2):
        for dx in (-2, -1, 0, 1, 2):
            if dy or dx:
                yield f"shift{dy},{dx}", perturb.shift(base, dy, dx)
    for amp in (0.5, 1.0, 1.5, 2.0, 3.0):
        for _ in range(2):
            yield f"jitter{amp}", perturb.jitter(base, amp, rng)


def wrong_variants(base, cache, rng, two_part):
    """Topology defects at and above the signature size filter A = cache.min_area (A and 2A px);
    smaller defects are invisible to the topology term by design (R6)."""
    A = cache.min_area
    for a in (A, 2 * A):
        m = perturb.punch_hole(base, a, rng)
        if m is not None:
            yield f"pinhole{a}", m
    for a in (A, 2 * A):
        reach = int(np.ceil(np.sqrt(a / np.pi))) + 2
        m = perturb.add_speck(base, a, float(rng.uniform(reach + 2, reach + 12)), rng)
        if m is not None:
            yield f"speck{a}", m
    if two_part:
        top = instance_topology_full(base)
        lab = top.full_labels(base.shape)
        if top.n == 2:
            for w in (1, 2):
                yield f"bridge{w}", perturb.bridge(base, lab == 1, lab == 2, w)


def item_pairs(args):
    name, subset, stack, two_part, seed = args
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack)
    rows = []
    for base_name, base in (("vote", cache.vote), ("rater", stack[int(rng.integers(len(stack)))])):
        if instance_topology_full(base).signature() not in cache.signatures:
            continue
        cands = []
        for cname, m in correct_candidates(base, rng):
            t = score(m, cache)
            if t.topo_any == 1.0:
                cands.append((cname, t))
        if not cands:
            continue
        for wname, m in wrong_variants(base, cache, rng, two_part):
            tw = score(m, cache)
            if tw.topo_any == 1.0:
                continue  # defect happens to match some rater: not a topology error
            near = [c for c in cands if abs(c[1].iou - tw.iou) <= MAX_DIOU]
            if not near:
                continue
            cname, tc = near[int(rng.integers(len(near)))]
            rows.append(dict(item=name, subset=subset, base=base_name, defect=wname.rstrip("0123456789"),
                             defect_full=wname, correct=cname, iou_c=tc.iou, iou_w=tw.iou,
                             **{f"R_c|{k}": reward(tc, c) for k, c in CONFIGS.items()},
                             **{f"R_w|{k}": reward(tw, c) for k, c in CONFIGS.items()}))
    return rows


def main():
    items = fixtures.load_all()
    args = [(k, fixtures.subset(ds, k), st, False, j) for j, (ds, k, st) in enumerate(items)]
    # composites: consecutive pairs within each subset, gap 3-8 px
    for sub in ("lidc/nodule", "riga/cup"):
        group = [(k, st) for ds, k, st in items if fixtures.subset(ds, k) == sub]
        for j in range(0, len(group) - 1, 2):
            (ka, a), (kb, b) = group[j], group[j + 1]
            st = composite(a, b, gap=3 + j % 6)
            args.append((f"{ka}+{kb}", sub + "+2part", st, True, 10_000 + j))
    with Pool(8) as pool:
        rows = [r for rs in pool.map(item_pairs, args, chunksize=4) for r in rs]
    (OUT / "k9_pairs.json").write_text(json.dumps(rows))
    n = len(rows)
    print(f"matched pairs: {n} (|dIoU| <= {MAX_DIOU}); median IoU of wrong {np.median([r['iou_w'] for r in rows]):.3f}")
    iou_ok = np.mean([r["iou_c"] > r["iou_w"] for r in rows])
    print(f"IoU alone ranks correctly: {iou_ok:.1%} (ties count as wrong: {np.mean([r['iou_c'] == r['iou_w'] for r in rows]):.1%} ties)")
    defects = sorted({r["defect"] for r in rows})
    subsets = sorted({r["subset"] for r in rows})
    print(f"\n{'config':48s} {'ALL':>7} " + " ".join(f"{d:>8}" for d in defects) + " | " + " ".join(f"{s:>18}" for s in subsets))
    for k in CONFIGS:
        ok = np.array([r[f"R_c|{k}"] > r[f"R_w|{k}"] for r in rows])
        cells = [f"{np.mean([o for o, r in zip(ok, rows) if r['defect'] == d]):8.1%}" for d in defects]
        scells = [f"{np.mean([o for o, r in zip(ok, rows) if r['subset'] == s]):18.1%}" for s in subsets]
        print(f"{k:48s} {ok.mean():7.1%} " + " ".join(cells) + " | " + " ".join(scells))
    hi = [r for r in rows if min(r["iou_c"], r["iou_w"]) >= 0.5]
    hard = [r for r in hi if r["iou_c"] < r["iou_w"]]
    k0 = next(iter(CONFIGS))
    for name, sub in (("high-IoU pairs (both IoU >= 0.5, the topology gate)", hi), ("hard: high-IoU and IoU favours wrong", hard)):
        print(f"\n{name}: n={len(sub)}; IoU alone {np.mean([r['iou_c'] > r['iou_w'] for r in sub]):.1%}; "
              + "; ".join(f"{k.split(' (')[0]} {np.mean([r[f'R_c|{k}'] > r[f'R_w|{k}'] for r in sub]):.2%}" for k in CONFIGS))
    low = [r for r in rows if min(r["iou_c"], r["iou_w"]) < 0.5]
    print(f"\nlow-IoU pairs (< 0.5): n={len(low)}, {k0}: {np.mean([r[f'R_c|{k0}'] > r[f'R_w|{k0}'] for r in low]):.1%}"
          " (no topology credit below the gate by design)")
    print("\npairs per defect:", {d: sum(r["defect"] == d for r in rows) for d in defects})
    print("pairs per subset:", {s: sum(r["subset"] == s for r in rows) for s in subsets})


if __name__ == "__main__":
    main()
