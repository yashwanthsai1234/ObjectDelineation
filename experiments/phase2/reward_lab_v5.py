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
  4 miss_beats_hit: share of groups where a missed answer has a higher v5 R than a hit; area_gate = that share > 1 %
  5 distance: median and quartiles of d / r over missed answers; loc_width = the kernel width when the median is <= 3,
    else the width that gives the median missed answer loc = 0.1
  6 adv_size: median |A| on non-dropped groups for v4 and v5, r = their ratio, beta = 1e-3 / r, adam_epsilon = 1e-8 / r
    (b3_grpo); with --spec both go into the B3 run spec's train stage (beta 0.0 with --beta-zero, owner decision
    B.8.3 item 5)
  7 hit_agreement: mean Spearman(A_v4, A_v5) over groups with every IoU_tc >= miss_iou
  8 topology: jump = c_topo x median r_topo of valid answers; hit_share = the c_topo r_topo share of the variance of
    v5 A in live all-hit groups (Gate B3: <= 0.5); fallback_margin_q95 = q95 of (IoU_tc of a legal answer no doctor
    drew) - (best IoU_tc of a valid answer of its group), over groups with a valid answer; c_topo_min =
    max(0, q95) / eps, the lower bound of c_topo (nan without a pair)
  9 runaway: share of truncated answers and their mean v5 R
 10 allmiss_grad_share: S = sum over all-miss groups of sum |A| / sum over all groups, v4 and v5, and S_v5 / S_v4;
    pass = S_v5 >= 0.5 S_v4
`gate` holds the pre-run Gate B3 lab rules (B.6.2, lab_gate): loc v5 mean > 0.3 with CI lower bound > 0 and
(v5 - v4) mean > 0; drop share v5 <= v4; miss beats hit <= 1 %; rule 10 pass; topology hit share <= 0.5.
Rule 8's margin compares an undrawn-signature answer (IoU_tc = IoU vs the vote) with the best valid answer of its
own group, the terms the scorer already returns, so no scorer change is needed.

Test surface: kept by the CLEAN16 decision (2026-11-10); only tests import it, so a change here is checked there.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import yaml
from scipy.stats import rankdata

from objdel.verify.reward.advantages import (
    _groups,
    _topology,
    advantages_v4,
    advantages_v5,
    reward_parts_v5,
    topology_valid,
)
from objdel.verify.reward.score import score_many
from objdel.verify.reward.spec import REWARD_V5

AREA_GATE_SHARE = 0.01  # rule 4
LOC_FAR, LOC_MIN = 3.0, 0.1  # rule 5: a median missed d / r above LOC_FAR gets loc >= LOC_MIN
GRAD_SHARE_MIN = 0.5  # rule 10
TOPO_SHARE_MAX = 0.5  # rule 8
LOC_MEAN_MIN = 0.3  # rule 1


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


def fallback_margins(G: dict, valid: np.ndarray) -> np.ndarray:
    """IoU_tc of every legal answer no doctor drew minus the best IoU_tc of a valid answer in its group, for the groups
    that have a valid answer (B3.7 rule 8; the H1 hack of H2_rl_math.md s3.4)."""
    best = np.where(valid, G["iou_tc"], -np.inf).max(1, keepdims=True)
    fallback = (G["legal"] >= 0.5) & (G["topo_any"] == 0) & np.isfinite(best)
    return (G["iou_tc"] - best)[fallback]


def allmiss_grad_share(a4: np.ndarray, a5: np.ndarray, allmiss: np.ndarray) -> dict:
    """Rule 10 on advantages (groups, G): S per version (0 without any advantage), S_v5 / S_v4 and the pass flag."""
    s4, s5 = (float(np.abs(a[allmiss]).sum() / np.abs(a).sum()) if np.abs(a).sum() else 0.0 for a in (a4, a5))
    return {"v4": s4, "v5": s5, "ratio": s5 / s4 if s4 else float("nan"), "pass": s5 >= GRAD_SHARE_MIN * s4}


def lab_gate(rep: dict) -> dict[str, bool]:
    """The pre-run Gate B3 lab rules (module docstring) of a report."""
    loc = rep["loc_signal"]
    return {"loc_signal": loc["v5"]["mean"] > LOC_MEAN_MIN and loc["v5"]["ci"][0] > 0
            and loc["v5_minus_v4"]["mean"] > 0,
            "drop_share": rep["drop_share"]["v5"] <= rep["drop_share"]["v4"],
            "miss_beats_hit": rep["miss_beats_hit"] <= AREA_GATE_SHARE,
            "allmiss_grad_share": rep["allmiss_grad_share"]["pass"],
            "topology": rep["topology"]["hit_share"] <= TOPO_SHARE_MAX}


def b3_grpo(r: float, beta_zero: bool = False) -> dict:
    """Rule 6: the B3 run's beta (0.0 by owner decision) and adam_epsilon, from the same r (ADR-0014)."""
    return {"beta": 0.0 if beta_zero else 1e-3 / r, "adam_epsilon": 1e-8 / r}


def write_grpo(spec: Path, values: dict) -> None:
    """Set values in the grpo mapping of the run spec's one `grpo:` line (the B3 train stage), keeping its other keys;
    the spec must still load."""
    from objdel.core.spec import load_spec

    text = spec.read_text()
    found = re.findall(r"^( +grpo: )(\{.*\})$", text, re.M)
    if len(found) != 1:
        raise ValueError(f"{spec}: expected one grpo line, found {len(found)}")
    head, old = found[0]
    new = yaml.safe_dump(yaml.safe_load(old) | values, default_flow_style=True, sort_keys=False).strip()
    spec.write_text(text.replace(head + old, head + new))
    load_spec(spec)


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
    d_miss = dr[miss & np.isfinite(dr)]
    valid = topology_valid(G, spec)
    margin = fallback_margins(G, valid)
    q95 = float(np.percentile(margin, 95)) if margin.size else float("nan")
    cut = np.asarray(truncated, bool).reshape(-1, g)
    med = float(np.median(d_miss)) if d_miss.size else float("nan")
    miss_beats = float(np.mean(beats)) if beats else 0.0
    return {
        "groups": len(R), "allmiss_groups": int(allmiss.sum()),
        "loc_signal": loc_signal,
        "drop_share": {k: float((np.ptp(a[allmiss], 1) == 0).mean()) if allmiss.any() else float("nan")
                       for k, a in (("v4", a4), ("v5", a5))},
        "variance_share": {"allmiss": variance_share(parts, a5, allmiss & live5)},
        "topology": {"jump": spec.c_topo * float(np.median(r_topo[valid])) if valid.any() else float("nan"),
                     "hit_share": variance_share(parts, a5, allhit & live5)["topo"],
                     "fallback_margin_q95": q95, "c_topo_min": (max(0.0, q95) if margin.size else q95) / spec.eps},
        "miss_beats_hit": miss_beats, "area_gate": miss_beats > AREA_GATE_SHARE,
        "distance": {"median": med,
                     "q25_q75": [float(np.percentile(d_miss, q)) for q in (25, 75)] if d_miss.size else [],
                     "loc_width": med / np.sqrt(-2 * np.log(LOC_MIN)) if med > LOC_FAR else spec.loc_width},
        "adv_size": {"median_v4": m4, "median_v5": m5, "r": r} | b3_grpo(r),
        "hit_agreement": float(np.mean([spearman(a4[i], a5[i]) for i in np.flatnonzero(allhit)]))
        if allhit.any() else float("nan"),
        "runaway": {"share": float(cut.mean()), "mean_v5_R": float(R[cut].mean()) if cut.any() else float("nan")},
        "allmiss_grad_share": allmiss_grad_share(a4, a5, allmiss),
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
    ap.add_argument("--spec", type=Path, help="B3 run spec (specs/b3_v5.yaml): write rule 6's beta and adam_epsilon")
    ap.add_argument("--beta-zero", action="store_true", help="owner decision B.8.3 item 5: beta 0.0 for B3")
    a = ap.parse_args(argv)
    from objdel.data.splits import val_half
    from objdel.policy import rows as dataset
    from objdel.verify import store as serve

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
    out["gate"] = lab_gate(out)
    if a.spec:
        write_grpo(a.spec, b3_grpo(out["adv_size"]["r"], a.beta_zero))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out | {"answers": str(a.answers), "view": str(a.view), "items": len(answers)},
                                indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
