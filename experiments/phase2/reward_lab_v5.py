"""B3.7 lab round (offline, $0): re-score saved eval answers (K = 8 per item, final_eval's answers.jsonl) on the
CPU checker and compare reward v4 and v5 advantages on val-shape groups (Implementation.md B3.7).

    python -I -m experiments.phase2.reward_lab_v5 --answers runs/a3/eval_val --view data/eval_views/lidc-val-shifted \
        --out runs/b3/lab_v5.json

Val-shape items only (splits.val_half), so val-select stays clean for Gate B3. A group is one item's K answers; a
missed answer has gated IoU_tc < miss_iou. Report keys (rule numbers of B3.7):
  1 loc_signal: per all-miss group with >= 3 distinct centroid distances d, Spearman(A, -d / r) (r = sqrt(radius2));
    groups with constant A count 0; mean and a patient (group_id) bootstrap 95 % CI for v4, v5 and v5 - v4
  2 drop_share: share of all-miss groups with A = 0 for the whole group, v4 and v5
  3 variance_share: share of the variance of v5 A in live all-miss groups from each weighted term of R
    (sum of A_term x A over sum of A^2, A_term = the term's RLOO advantage; the shares sum to 1)
  4 miss_beats_hit: share of groups where a missed answer has a higher v5 R than a hit
  5 distance: median and quartiles of d / r over missed answers
  6 adv_size: median |A| on non-dropped groups for v4 and v5, r = their ratio, beta = 1e-3 / r, adam_epsilon = 1e-8 / r
  7 hit_agreement: mean Spearman(A_v4, A_v5) over groups with every IoU_tc >= miss_iou
  8 topology: jump = c_topo x median r_topo of valid answers; hit_share = the c_topo r_topo share of the variance of
    v5 A in live all-hit groups (Gate B3: <= 0.5)
  9 runaway: share of truncated answers and their mean v5 R
 10 allmiss_grad_share: S = sum over all-miss groups of sum |A| / sum over all groups, v4 and v5, and S_v5 / S_v4
Rule 8's c_topo lower bound (q95 of IoU vs vote minus the best same-signature IoU for undrawn signatures) is not
computed here: the scorer returns no per-signature IoU.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from objdel.reward import REWARD_V5, advantages_v4, advantages_v5, score_many
from objdel.reward.advantages import _groups, _topology, reward_parts_v5, topology_valid


def spearman(a, b) -> float:
    """Spearman correlation with average ranks for ties; 0 when either side is constant."""
    ra, rb = rankdata(a), rankdata(b)
    if np.ptp(ra) == 0 or np.ptp(rb) == 0:
        return 0.0
    return float(np.corrcoef(ra, rb)[0, 1])


def bootstrap(values, clusters, n=2000, seed=0) -> list[float]:
    """Patient bootstrap 95 % CI of the mean of values (one per group), resampling clusters."""
    values, clusters = np.asarray(values, float), np.asarray(clusters)
    ids = np.unique(clusters)
    if values.size == 0:
        return [float("nan")] * 2
    by = [values[clusters == c] for c in ids]
    rng = np.random.default_rng(seed)
    means = [np.concatenate([by[i] for i in rng.integers(0, len(by), len(by))]).mean() for _ in range(n)]
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def rloo(x: np.ndarray) -> np.ndarray:
    """Leave-one-out advantage of each row's values (B, G), as advantages_v5 takes it."""
    return x - (x.sum(1, keepdims=True) - x) / (x.shape[1] - 1)


def variance_share(parts: dict, a: np.ndarray, rows: np.ndarray) -> dict:
    """Share of sum(A^2) over the rows from each part: sum(rloo(part) A) / sum(A^2); nan without variance."""
    tot = (a[rows] ** 2).sum()
    return {k: float((rloo(v)[rows] * a[rows]).sum() / tot) if tot else float("nan") for k, v in parts.items()}


def report(T: dict, d_over_r: np.ndarray, patients: list[str], truncated: np.ndarray, spec=REWARD_V5) -> dict:
    """The B3.7 report of score terms T (flat, K = spec.group per item), centroid distances d / r (nan for an empty
    answer), the group's patient id per group and the truncated flag per answer."""
    g = spec.group
    a4, a5 = advantages_v4(T, g).reshape(-1, g), advantages_v5(T, g, spec).reshape(-1, g)
    G = _groups(T, g, ("loc",))
    r_topo = _topology(G, g, spec)
    parts = reward_parts_v5(G, r_topo, spec)
    R = sum(parts.values())
    iou = G["iou_tc"] * (G["legal"] >= 0.5)
    miss = iou < spec.miss_iou
    allmiss = miss.all(1)
    dr = np.asarray(d_over_r, float).reshape(-1, g)
    sig = {"v4": [], "v5": []}
    keep = []
    for i in np.flatnonzero(allmiss):
        if len(np.unique(dr[i][np.isfinite(dr[i])])) >= 3:
            keep.append(i)
            x = -np.nan_to_num(dr[i], nan=np.inf)
            for name, a in (("v4", a4), ("v5", a5)):
                sig[name].append(spearman(a[i], x) if np.ptp(a[i]) > 0 else 0.0)
    pid = np.asarray(patients)[keep]
    diff = np.subtract(sig["v5"], sig["v4"])
    loc_signal = {k: {"mean": float(np.mean(v)) if v else float("nan"), "ci": bootstrap(v, pid)}
                  for k, v in sig.items()} | {"v5_minus_v4": {"mean": float(diff.mean()) if diff.size else
                                                               float("nan"), "ci": bootstrap(diff, pid)},
                                              "groups": len(keep)}
    hits = ~miss
    beats = [(R[i][miss[i]].max() > R[i][hits[i]].min()) for i in range(len(R)) if miss[i].any() and hits[i].any()]
    live4, live5 = np.ptp(a4, 1) > 0, np.ptp(a5, 1) > 0
    m4, m5 = float(np.median(np.abs(a4[live4]))), float(np.median(np.abs(a5[live5])))
    r = m4 / m5 if m5 else float("nan")
    allhit = hits.all(1)
    s4, s5 = (float(np.abs(a[allmiss]).sum() / np.abs(a).sum()) if np.abs(a).sum() else 0.0 for a in (a4, a5))
    d_miss = dr[miss & np.isfinite(dr)]
    valid = topology_valid(G, spec)
    cut = np.asarray(truncated, bool).reshape(-1, g)
    return {
        "groups": len(R), "allmiss_groups": int(allmiss.sum()),
        "loc_signal": loc_signal,
        "drop_share": {k: float((np.ptp(a[allmiss], 1) == 0).mean()) if allmiss.any() else float("nan")
                       for k, a in (("v4", a4), ("v5", a5))},
        "variance_share": {"allmiss": variance_share(parts, a5, allmiss & live5)},
        "topology": {"jump": spec.c_topo * float(np.median(r_topo[valid])) if valid.any() else float("nan"),
                     "hit_share": variance_share(parts, a5, allhit & live5)["topo"]},
        "miss_beats_hit": float(np.mean(beats)) if beats else 0.0,
        "distance": {"median": float(np.median(d_miss)) if d_miss.size else float("nan"),
                     "q25_q75": [float(np.percentile(d_miss, q)) for q in (25, 75)] if d_miss.size else []},
        "adv_size": {"median_v4": m4, "median_v5": m5, "r": r, "beta": 1e-3 / r, "adam_epsilon": 1e-8 / r},
        "hit_agreement": float(np.mean([spearman(a4[i], a5[i]) for i in np.flatnonzero(allhit)]))
        if allhit.any() else float("nan"),
        "runaway": {"share": float(cut.mean()), "mean_v5_R": float(R[cut].mean()) if cut.any() else float("nan")},
        "allmiss_grad_share": {"v4": s4, "v5": s5, "ratio": s5 / s4 if s4 else float("nan")},
    }


def distance_over_r(terms: dict, caches: dict, names) -> float:
    """Mean over objects of the centroid distance over r; nan when any object is empty or unparsed."""
    out = []
    for n in names:
        t = terms.get(n)
        if t is None or t.area == 0:
            return float("nan")
        rx, ry, k = caches[n].ref_sum
        out.append(np.hypot(t.sum_x / t.area - rx / k, t.sum_y / t.area - ry / k) / np.sqrt(caches[n].radius2))
    return float(np.mean(out))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--answers", type=Path, required=True, help="a final_eval output folder with answers.jsonl")
    ap.add_argument("--view", type=Path, required=True, help="the eval view the answers were sampled on")
    ap.add_argument("--half", choices=("shape",), default="shape")
    ap.add_argument("--out", type=Path, required=True, help="report JSON")
    a = ap.parse_args(argv)
    from objdel.data.splits import val_half
    from objdel.train import dataset, serve

    rows = {r["key"]: r for r in dataset.prompts(a.view, "lidc", "val")}
    serve.preload(a.view, ["lidc"])
    caches = serve.loader(a.view)
    answers = [json.loads(line) for line in (a.answers / "answers.jsonl").open()]
    answers = [d for d in answers if val_half(d["group_id"]) == a.half and len(d["samples"]) == REWARD_V5.group]
    jobs, meta = [], []
    for d in answers:
        row, c = rows[d["item_id"]], caches(d["item_id"])
        names = tuple(row["names"])
        for text, cut in zip(d["samples"], d["truncated"], strict=True):
            jobs.append((text, cut, names, c, row["W"], row["H"]))
            meta.append((names, c, cut))
    scored = score_many(jobs)
    T = {k: np.array([t[k] for t, _ in scored], np.float64) for k in scored[0][0] if k not in ("error", "defects")}
    dr = np.array([distance_over_r(o, c, n) for (_, o), (n, c, _) in zip(scored, meta, strict=True)])
    out = report(T, dr, [d["group_id"] for d in answers], np.array([m[2] for m in meta]))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out | {"answers": str(a.answers), "view": str(a.view), "items": len(answers)},
                                indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
