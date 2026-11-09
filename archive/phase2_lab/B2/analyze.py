"""B2 analysis: GRPO-signal metrics per group of G = 8 with nodule-cluster bootstrap CIs, paired depth / seed /
interpolation differences, and the pre-registered verdicts (prereg.py).

    .venv/bin/python -I -m experiments.phase2.lab.B2.analyze [--out <scratch>/phase2/B2/analysis.json]

Reads <scratch>/phase2/B2/scored/<run>.json (score.py). Paired differences use the same resampled nodules for both
runs (items aligned by item_id). Percentile CIs, B = 2,000.
"""

import argparse
import json

import numpy as np

from . import prereg as P
from .score import B2

INTERP = {a: f"lfm_a{a:g}" for a in P.ALPHAS}


def group_metrics(v: dict) -> dict[str, np.ndarray]:
    R, legal, iou = (np.asarray(v[k], float) for k in ("R", "legal", "iou_tc"))
    adv = np.abs(np.asarray(v["advantage"], float))
    L = legal > 0
    return {
        "r_std": R.std(1),
        "sigma": np.asarray(v["sigma"], float),
        "mixed": (L.any(1) & ~L.all(1)).astype(float),
        "headroom": R.max(1) - R.mean(1),
        "distinct": np.asarray(v["distinct"], float),
        "raster_distinct": 1 - np.asarray(v["dup"], float),
        "legal_wrong": (L & (iou < P.IOU_OK)).mean(1),
        "active": (adv.max(1) > 1e-12).astype(float),
        "abs_adv": adv.mean(1),
        "legal": legal.mean(1),
        "r_mean": R.mean(1),
        "iou_mean": iou.mean(1),
    }


def load(run: str) -> dict:
    s = json.loads((B2 / "scored" / f"{run}.json").read_text())
    return {view: {"ids": v["item_id"], "m": group_metrics(v)} for view, v in s.items()}


def draws(ids: list[str], seed: int = 0) -> list[np.ndarray]:
    """B nodule-cluster bootstrap index draws (nodule = item id up to the slice)."""
    nod = [i.rsplit("/", 1)[0] for i in ids]
    groups = [np.array([j for j, x in enumerate(nod) if x == u]) for u in dict.fromkeys(nod)]
    rng = np.random.default_rng(seed)
    return [np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))]) for _ in range(P.BOOT)]


def ci(x: np.ndarray, D: list[np.ndarray]) -> list[float]:
    reps = np.array([x[d].mean() for d in D])
    q = (1 - P.CI) / 2
    return [round(float(x.mean()), 4), round(float(np.quantile(reps, q)), 4), round(float(np.quantile(reps, 1 - q)), 4)]


def table(run: dict) -> dict:
    out = {}
    for view, v in run.items():
        D = draws(v["ids"])
        out[view] = {k: ci(v["m"][k], D) for k in P.METRICS} | {"items": len(v["ids"])}
    return out


def paired(a: dict, b: dict) -> dict:
    """a - b per view and metric, on the items both have."""
    out = {}
    for view in a:
        ids = [i for i in a[view]["ids"] if i in set(b[view]["ids"])]
        ia = [a[view]["ids"].index(i) for i in ids]
        ib = [b[view]["ids"].index(i) for i in ids]
        D = draws(ids)
        out[view] = {k: ci(a[view]["m"][k][ia] - b[view]["m"][k][ib], D) for k in P.METRICS} | {"items": len(ids)}
    return out


def lo_gt0(d):
    return d[1] > 0


def hi_lt0(d):
    return d[2] < 0


def verdicts(depth: dict, seed: dict, interp: dict, tables: dict) -> dict:
    v = {}
    view = P.PRIMARY_VIEW
    noise = {k: abs(seed["lfm2.5-vl-450m"][view][k][0]) for k in P.METRICS}
    per_arm = {}
    for arm, d in depth.items():
        dd = d[view]
        per_arm[arm] = {
            "falls": {k: hi_lt0(dd[k]) and abs(dd[k][0]) > noise[k] for k in P.PRIMARY},
            "sig_fall_any": any(hi_lt0(dd[k]) for k in P.PRIMARY),
            "sig_rise_any": any(lo_gt0(dd[k]) for k in P.PRIMARY),
            "distinct_drop": dd["distinct"][0] <= -P.H1B_DISTINCT_DROP and hi_lt0(dd["distinct"]),
        }
    n_acc = sum(all(x["falls"].values()) for x in per_arm.values())
    if n_acc >= P.H1_ARMS_NEEDED:
        h1 = "ACCEPT"
    elif (
        not any(x["sig_fall_any"] for x in per_arm.values())
        and sum(x["sig_rise_any"] for x in per_arm.values()) >= P.H1_ARMS_NEEDED
    ):
        h1 = "REJECT"
    else:
        h1 = "INCONCLUSIVE"
    v["H-B2-1"] = {"verdict": h1, "per_arm": per_arm, "lfm_seed_noise": {k: noise[k] for k in P.PRIMARY}}
    v["H-B2-1b"] = sum(x["distinct_drop"] for x in per_arm.values()) >= P.H1_ARMS_NEEDED
    if interp:
        early = tables["lfm2.5-vl-450m_s0"][view]["sigma"][0]
        per_alpha = {}
        for alpha, d in interp.items():
            vs1, vs36 = d["vs_alpha1"], d["vs_36"]
            c1 = all(vs1[vw]["legal"][0] >= -P.INTERP_LEGAL_TOL for vw in vs1)
            c2 = any(lo_gt0(vs1[view][k]) for k in P.PRIMARY)
            c3 = tables[INTERP[alpha]][view]["sigma"][0] >= early - P.INTERP_SIGMA_TOL
            loses = any(hi_lt0(vs1[vw]["legal"]) for vw in vs1)
            lower = all(hi_lt0(vs1[view][k]) for k in P.PRIMARY)
            dominates36 = all(lo_gt0(vs36[vw]["legal"]) for vw in vs36) and not any(
                hi_lt0(vs36[view][k]) for k in P.PRIMARY
            )
            per_alpha[alpha] = {
                "i_legal_kept": c1,
                "ii_signal_up": c2,
                "iii_sigma_vs_early": c3,
                "loses_legal": loses,
                "signal_lower": lower,
                "dominates_36": dominates36,
            }
        if any(x["i_legal_kept"] and x["ii_signal_up"] and x["iii_sigma_vs_early"] for x in per_alpha.values()):
            h2 = "ACCEPT"
        elif all(x["loses_legal"] or x["signal_lower"] for x in per_alpha.values()):
            h2 = "REJECT"
        else:
            h2 = "INCONCLUSIVE"
        v["H-B2-2"] = {"verdict": h2, "per_alpha": per_alpha}
    v["D12_rule"] = "R-A (signal-first lexicographic)" if h1 == "ACCEPT" else "R-B (deepest eligible, signal veto)"
    v["D12_repair"] = "interpolation" if v.get("H-B2-2", {}).get("verdict") == "ACCEPT" else "next-earlier branch"
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(B2 / "analysis.json"))
    a = ap.parse_args()
    runs = sorted({r for p in (*P.DEPTH_RUNS.values(), *P.SEED_RUNS.values()) for r in p})
    runs += [r for r in INTERP.values() if (B2 / "scored" / f"{r}.json").exists()]
    data = {r: load(r) for r in runs}
    tables = {r: table(data[r]) for r in runs}
    depth = {arm: paired(data[deep], data[shallow]) for arm, (shallow, deep) in P.DEPTH_RUNS.items()}
    seed = {arm: paired(data[s1], data[s0]) for arm, (s0, s1) in P.SEED_RUNS.items()}
    interp = {
        alpha: {"vs_alpha1": paired(data[r], data[P.INTERP_RUN]), "vs_36": paired(data[r], data["lfm2.5-vl-450m_s0"])}
        for alpha, r in INTERP.items()
        if r in data
    }
    res = {
        "tables": tables,
        "depth_150_minus_36": depth,
        "seed_s1_minus_s0": seed,
        "interp": interp,
        "verdicts": verdicts(depth, seed, interp, tables),
    }
    open(a.out, "w").write(json.dumps(res, indent=1, default=str))
    print(json.dumps(res["verdicts"], indent=1, default=str))


if __name__ == "__main__":
    main()
