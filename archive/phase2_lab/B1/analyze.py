"""B1 analysis: tables with item-bootstrap CIs and the pre-registered verdicts (criteria fixed below before any run).

    .venv/bin/python -I -m experiments.phase2.lab.B1.analyze [--out <scratch>/analysis.json]

Reads <scratch>/zeroshot/<arm>.json (zeroshot.py) and <scratch>/runs/<arm>_s<seed>/{result,per_item}.json
(train_eval.py). Paired comparisons resample the same items for both arms (B = 2,000, percentile 95 % CI).
"""
import argparse
import json
from itertools import combinations

import numpy as np

from . import common as C
from .prereg import (
    FORMAT_LEGAL,
    FORMAT_PARSE,
    H1_FAIL_IOU,
    H1_IOU,
    H1_PARSE,
    H8_DUP,
    H8_LEGAL,
    LEADING_ZERO_MATERIAL,
    SMOL_GAP,
    STEPS_SECONDARY,
)

ARMS = list(C.ARMS)


def load_runs(seed=0, steps=36, suffix=""):
    out = {}
    for arm in ARMS:
        d = C.SCRATCH / "runs" / (f"{arm}_s{seed}" + (f"_n{steps}" if steps != 36 else "") + suffix)
        if (d / "per_item.json").exists():
            out[arm] = {"res": json.loads((d / "result.json").read_text()),
                        "items": json.loads((d / "per_item.json").read_text())}
    return out


def arrays(p):
    a = {k: np.asarray(p[k], float) for k in ("iou_tc", "format", "legal", "tokens", "zero_led")}
    a["dup"], a["signal"], a["prior"] = (np.asarray(p[k], float) for k in ("dup", "signal", "prior"))
    a["lead0"] = np.array([[e == "leading_zero" for e in row] for row in p["error"]], float)
    return a


STATS = {
    "format": lambda a, i: a["format"][i].mean(), "legal": lambda a, i: a["legal"][i].mean(),
    "iou_tc_median": lambda a, i: np.median(a["iou_tc"][i].ravel()),
    "iou_tc_mean": lambda a, i: a["iou_tc"][i].mean(),
    "iou_ge_0.5": lambda a, i: (a["iou_tc"][i] >= 0.5).mean(),
    "dup": lambda a, i: a["dup"][i].mean(), "signal": lambda a, i: a["signal"][i].mean(),
    "tokens": lambda a, i: a["tokens"][i].mean(),
    "minus_prior": lambda a, i: (a["iou_tc"][i].mean(1) - a["prior"][i]).mean(),
}


def arm_table(runs):
    out = {}
    for arm, r in runs.items():
        out[arm] = {}
        for view, p in r["items"].items():
            a = arrays(p)
            n = len(a["prior"])
            row = {k: C.ci(lambda i, f=f: f(a, i), n) for k, f in STATS.items()}
            row["prior_median"] = round(float(np.median(a["prior"])), 4)
            row["win_share_vs_prior"] = round(float((a["iou_tc"].mean(1) > a["prior"]).mean()), 4)
            row["zero_led_samples"] = round(float((a["zero_led"] > 0).mean()), 4)
            row["leading_zero_fail"] = round(float(a["lead0"].mean()), 4)
            row["truncated"] = r["res"]["eval"][view]["truncated"]
            out[arm][view] = row
        out[arm]["train"] = {k: r["res"]["train"].get(k) for k in
                             ("s_per_step", "s_per_step_mean", "peak_mps_gb", "wall_s", "steps", "rows",
                              "target_tokens_max", "trainable_params")}
        out[arm]["train"]["final_loss"] = r["res"]["train"]["loss"][-1] if r["res"]["train"].get("loss") else None
        out[arm]["eval_peak_mps_gb"] = r["res"]["eval"]["peak_mps_gb"]
        out[arm]["gen_s_per_item"] = {v: r["res"]["eval"][v]["gen_s_per_item"] for v in r["items"]}
    return out


def paired(runs, view, key="iou_tc_median"):
    out = {}
    for x, y in combinations(runs, 2):
        ax, ay = arrays(runs[x]["items"][view]), arrays(runs[y]["items"][view])
        assert runs[x]["items"][view]["item_id"] == runs[y]["items"][view]["item_id"]
        f = STATS[key]
        out[f"{x} - {y}"] = C.ci(lambda i: f(ax, i) - f(ay, i), len(ax["prior"]))
    return out


def sig(d):  # +1 / -1 when the CI excludes 0, else 0
    return 1 if d[1] > 0 else -1 if d[2] < 0 else 0


def verdicts(tab, pairs, zs):
    v = {}
    q = zs.get("qwen3.5-0.8b", {})
    sh, ce = q.get("riga/shifted/native"), q.get("riga/centred/native")
    if sh:
        ok = sh["parse"][0] >= H1_PARSE and sh["box_iou_median"][0] >= H1_IOU and sh["minus_prior_mean"][1] > 0
        fail = all(q.get(f"{ds}/shifted/native", {}).get("box_iou_median", [1])[0] < H1_FAIL_IOU
                   for ds in ("riga", "lidc"))
        v["H-A1-1"] = "accept" if ok else "reject" if fail else "inconclusive"
        v["H-A1-1 centred (A1 original)"] = ce["parse"][0] >= H1_PARSE and ce["box_iou_median"][0] >= H1_IOU
    view = "shifted"
    learn = {arm: t[view]["format"][0] >= FORMAT_PARSE and t[view]["legal"][0] >= FORMAT_LEGAL for arm, t in tab.items()}
    v["learns_format_shifted"] = learn
    v["learns_format_centred"] = {arm: t["centred"]["format"][0] >= FORMAT_PARSE
                                  and t["centred"]["legal"][0] >= FORMAT_LEGAL for arm, t in tab.items()}
    v["uses_image_shifted"] = {arm: t[view]["minus_prior"][1] > 0 for arm, t in tab.items()}
    L, Q, S = "lfm2.5-vl-450m", "qwen3.5-0.8b", "smolvlm-256m"
    p = pairs[view]

    def d(x, y):
        return sig(p[f"{x} - {y}"]) if f"{x} - {y}" in p else -sig(p[f"{y} - {x}"])

    if all(a in tab for a in (L, Q, S)):
        lq, qs = d(L, Q), d(Q, S)
        smol_fails = not learn[S] and learn[Q]
        if lq >= 0 and (qs > 0 or smol_fails):
            v["H-A1-2"] = "accept"
        elif lq < 0 or qs < 0:
            v["H-A1-2"] = "reject"
        else:
            v["H-A1-2"] = "inconclusive"
        best = max(tab[a][view]["iou_tc_median"][0] for a in tab)
        v["drop_smolvlm (A1 rule)"] = best - tab[S][view]["iou_tc_median"][0] > SMOL_GAP
    v["H-A1-8"] = {arm: {vw: t[vw]["legal"][0] >= H8_LEGAL and t[vw]["dup"][0] <= H8_DUP for vw in ("centred", "shifted")}
                   for arm, t in tab.items()}
    if L in tab:
        v["lfm_leading_zero_material"] = max(tab[L][vw]["leading_zero_fail"] for vw in ("centred", "shifted")) \
            >= LEADING_ZERO_MATERIAL
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(C.SCRATCH / "analysis.json"))
    a = ap.parse_args()
    zs = {arm: json.loads((C.SCRATCH / "zeroshot" / f"{arm}.json").read_text()) for arm in ARMS
          if (C.SCRATCH / "zeroshot" / f"{arm}.json").exists()}
    res = {"zeroshot": zs}
    for seed, steps, suffix in ((0, 36, ""), (1, 36, ""), (0, STEPS_SECONDARY, ""), (0, 36, "_lm")):  # _lm: AMEND 5
        runs = load_runs(seed, steps, suffix)
        if not runs:
            continue
        key = f"seed{seed}" + (f"_n{steps}" if steps != 36 else "") + suffix
        tab = arm_table(runs)
        pairs = {view: {k: paired(runs, view, k) for k in ("iou_tc_median", "iou_tc_mean", "legal")}
                 for view in ("centred", "shifted")}
        res[key] = {"table": tab, "paired": pairs,
                              "verdicts": verdicts(tab, {vw: pairs[vw]["iou_tc_median"] for vw in pairs}, zs)}
    if "seed1" in res:  # seed variation for arms with both seeds
        res["seed_delta"] = {arm: {vw: {k: round(res["seed1"]["table"][arm][vw][k][0]
                                             - res["seed0"]["table"][arm][vw][k][0], 4)
                                        for k in ("format", "legal", "iou_tc_median", "dup", "signal")}
                                   for vw in ("centred", "shifted")}
                             for arm in res["seed1"]["table"]}
    open(a.out, "w").write(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: v.get("verdicts") for k, v in res.items() if k.startswith("seed")}, indent=1, default=str))


if __name__ == "__main__":
    main()
