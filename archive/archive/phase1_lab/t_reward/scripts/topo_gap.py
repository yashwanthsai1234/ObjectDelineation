"""Iteration 2, task 1: the real per-topology IoU gap on the fixtures.

Items whose raters disagree on per-instance topology. Each rater's mask is a "sample". It is
scored against references built from the OTHER raters (leave-own-rater-out, as for a new
expert) and, for the vote, also from all raters (the iteration-1 definition):
  (c)  IoU vs the strict-majority vote               [all raters | other raters]
  (a)  IoU vs the strict-majority consensus of the other raters with the same topology
       (fallback: vote of the other raters)
  (b)  max IoU over the other raters with the same topology (fallback: vote of the others)
Gap per item = mean IoU of samples whose topology matches the vote's topology minus the mean
of samples whose topology does not (items where both groups exist).

Every minority topology here has a single rater, so leave-own-rater-out forces (a)/(b) back to
the vote for minority samples. Part 2 therefore emulates production: new samples are jittered
copies (smooth boundary noise, 1 px, topology preserved, 10 per rater) of each rater's mask,
scored against references built from ALL raters (the sample is not one of them).

    ../.venv/bin/python scripts/topo_gap.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.topology import instance_topology_full  # noqa: E402


def iou(a, b):
    u = np.count_nonzero(a | b)
    return np.count_nonzero(a & b) / u if u else 1.0


def vote(ms):
    return 2 * np.sum(ms, 0) > len(ms)


def main():
    rows = []
    n_items = {}
    for ds, key, st in fixtures.load_all():
        sigs = [instance_topology_full(m).signature() for m in st]
        n_items.setdefault(ds, [0, 0])
        n_items[ds][0] += 1
        if len(set(sigs)) == 1:
            continue
        n_items[ds][1] += 1
        v_all = vote(st)
        vsig = instance_topology_full(v_all).signature()
        for j, (m, s) in enumerate(zip(st, sigs)):
            others = [i for i in range(len(st)) if i != j]
            v_oth = vote(st[others])
            same = [i for i in others if sigs[i] == s]
            cons = vote(st[same]) if same else v_oth
            rows.append(dict(key=key, k=len(st), sig=s, match=s == vsig, freq=sigs.count(s) / len(st),
                             c_all=iou(m, v_all), c_loo=iou(m, v_oth), a=iou(m, cons),
                             b=max((iou(m, st[i]) for i in same), default=iou(m, v_oth)), has_same=bool(same)))
    print("items with rater topology disagreement:", {d: f"{n[1]}/{n[0]}" for d, n in n_items.items()})
    print(f"rater-samples: {len(rows)}; vote topology matched by {sum(r['match'] for r in rows)}; "
          f"samples with another rater of the same topology: {sum(r['has_same'] for r in rows)}")
    keys = sorted({r["key"] for r in rows})
    print("\nper-item gap = mean IoU(matches vote topology) - mean IoU(does not)")
    print(f"{'reference':52s} {'n items':>7} {'mean gap':>9} {'median':>7} {'[min, max]':>16} {'mean IoU match / other':>24}")
    names = {"c_all": "(c) vote, all raters (iteration-1 IoU)", "c_loo": "(c') vote, other raters",
             "a": "(a) consensus of same-topology others (fb: vote)", "b": "(b) max over same-topology others (fb: vote)"}
    rng = np.random.default_rng(0)
    for f, name in names.items():
        gaps, mm, mo = [], [], []
        for k in keys:
            g = [r for r in rows if r["key"] == k]
            m1 = [r[f] for r in g if r["match"]]
            m0 = [r[f] for r in g if not r["match"]]
            if m1 and m0:
                gaps.append(np.mean(m1) - np.mean(m0))
                mm += m1
                mo += m0
        gaps = np.array(gaps)
        boots = [np.mean(rng.choice(gaps, len(gaps))) for _ in range(5000)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        print(f"{name:52s} {len(gaps):7d} {gaps.mean():+9.3f} {np.median(gaps):+7.3f} [{gaps.min():+.3f}, {gaps.max():+.3f}]"
              f" {np.mean(mm):11.3f} / {np.mean(mo):.3f}   95% CI of mean gap [{lo:+.3f}, {hi:+.3f}]")
    print("\nper item (k, rater topologies, gap c_all / a / b):")
    for k in keys:
        g = [r for r in rows if r["key"] == k]
        if not (any(r["match"] for r in g) and not all(r["match"] for r in g)):
            print(f"  {k:28s} k={g[0]['k']} sigs={[r['sig'] for r in g]} (vote topology matches no rater / all)")
            continue
        def gap(f, g=g):
            return np.mean([r[f] for r in g if r["match"]]) - np.mean([r[f] for r in g if not r["match"]])
        print(f"  {k:28s} k={g[0]['k']} sigs={[r['sig'] for r in g]} gap {gap('c_all'):+.3f} / {gap('a'):+.3f} / {gap('b'):+.3f}")


def production_emulation(n_draws=10):
    rng = np.random.default_rng(1)
    rows = []
    for _, key, st in fixtures.load_all():
        sigs = [instance_topology_full(m).signature() for m in st]
        if len(set(sigs)) == 1:
            continue
        v_all = vote(st)
        vsig = instance_topology_full(v_all).signature()
        cons = {s: vote(st[[i for i in range(len(st)) if sigs[i] == s]]) for s in set(sigs)}
        for m, s in zip(st, sigs):
            same = [i for i in range(len(st)) if sigs[i] == s]
            for _ in range(n_draws):
                for _try in range(20):
                    x = perturb.jitter(m, 1.0, rng, corr=3.0)
                    if instance_topology_full(x).signature() == s:
                        break
                else:
                    continue
                rows.append(dict(key=key, match=s == vsig, c=iou(x, v_all), a=iou(x, cons[s]),
                                 b=max(iou(x, st[i]) for i in same)))
    keys = sorted({r["key"] for r in rows})
    print(f"\nproduction emulation: {len(rows)} jittered samples (1 px smooth jitter, topology kept)")
    rng = np.random.default_rng(0)
    for f, name in (("c", "(c) vote, all raters"), ("a", "(a) consensus of same-topology raters"),
                    ("b", "(b) max over same-topology raters")):
        gaps = []
        for k in keys:
            g = [r for r in rows if r["key"] == k]
            m1, m0 = [r[f] for r in g if r["match"]], [r[f] for r in g if not r["match"]]
            if m1 and m0:
                gaps.append(np.mean(m1) - np.mean(m0))
        gaps = np.array(gaps)
        lo, hi = np.percentile([np.mean(rng.choice(gaps, len(gaps))) for _ in range(5000)], [2.5, 97.5])
        print(f"  {name:40s} items {len(gaps)}  mean gap {gaps.mean():+.3f} [{lo:+.3f}, {hi:+.3f}]  median {np.median(gaps):+.3f}"
              f"  |gap| median {np.median(np.abs(gaps)):.3f}  range [{gaps.min():+.3f}, {gaps.max():+.3f}]")


if __name__ == "__main__":
    main()
    production_emulation()
