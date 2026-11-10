"""Question A experiments: does GRPO with a set-valued topology reward collapse the topology
distribution, and which fix preserves it? Writes results/a_*.json, results/a_tables.md, figs/a_*.png."""
import json
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import grpo_toy as gt
from plotstyle import PALETTE, style

S = 20  # seeds per configuration
T = 300  # GRPO visits per prompt
MIXED = [n for n in gt.SUITE if n != "unanimous"]
OUT = {}
MD = []


def ci(x):
    x = np.asarray(x, float)
    m = x.mean()
    h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0
    return m, h


def summarize(o, prompts=None):
    names = o["names"]
    idx = [names.index(n) for n in (prompts or MIXED) if n in names]
    hf = -(o["freq"] * np.log(np.where(o["freq"] > 0, o["freq"], 1))).sum(-1)
    res = {}
    for k in ["invalid", "tv", "js", "ent", "minret"]:
        v = o[k][-1][:, idx]
        if k == "ent":
            v = v - hf[idx]
        res[k] = v.mean(1)  # per seed
    inv = o["invalid"][:, :, idx].mean(2)  # n_log,S
    hit = inv < 0.01
    first = np.where(hit.any(0), o["steps"][hit.argmax(0)], np.nan)
    res["t01"] = first
    res["tv_curve"] = o["tv"][:, :, idx].mean(2)  # n_log,S
    res["inv_curve"] = inv
    res["ent_curve"] = (o["ent"][:, :, idx] - hf[idx]).mean(2)
    return res


def fmt(x, d=3):
    m, h = ci(x[~np.isnan(x)]) if np.isnan(x).any() else ci(x)
    s = f"{m:.{d}f} ± {h:.{d}f}"
    if np.isnan(x).any():
        s += f" ({np.isnan(x).sum()}/{len(x)} never)"
    return s


def table(title, rows, cols=("invalid", "tv", "js", "ent", "minret", "t01")):
    hdr = {"invalid": "invalid mass", "tv": "TV to raters", "js": "JS (bits)", "ent": "H(π)−H(raters) (nats)",
           "minret": "min-mode retention", "t01": "visits to invalid<1%"}
    MD.append(f"\n### {title}\n")
    MD.append("| variant | " + " | ".join(hdr[c] for c in cols) + " |")
    MD.append("|---" * (len(cols) + 1) + "|")
    for name, r in rows:
        cells = []
        for c in cols:
            d = 0 if c == "t01" else (4 if c == "invalid" else 3)
            cells.append(fmt(r[c], d))
        MD.append(f"| {name} | " + " | ".join(cells) + " |")


def save(key, rows):
    OUT[key] = {n: {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()
                    if not k.endswith("curve")} for n, r in rows}


VARIANTS = [
    ("(1) set (any rater)", dict(reward="set")),
    ("(1b) near (repo default topo_near)", dict(reward="near")),
    ("(2) freq (share of raters)", dict(reward="freq")),
    ("(3) majority only", dict(reward="majority")),
    ("(4) set + coverage λ=0.25", dict(reward="cov", lam=0.25)),
    ("(4) set + coverage λ=0.5", dict(reward="cov", lam=0.5)),
    ("(4) set + coverage λ=1", dict(reward="cov", lam=1.0)),
    ("(4b) set + unique-answer bonus λ=0.5", dict(reward="uniq", lam=0.5)),
    ("(5) set + KL(SFT) β=0.01", dict(reward="set", beta=0.01)),
    ("(5) set + KL(SFT) β=0.04", dict(reward="set", beta=0.04)),
    ("(5) set + KL(SFT) β=0.1", dict(reward="set", beta=0.1)),
    ("(5) set + KL(SFT) β=0.3", dict(reward="set", beta=0.3)),
    ("(6a) log-freq reward", dict(reward="logfreq")),
    ("(6b) set + share λ=0.25", dict(reward="share", lam=0.25)),
    ("(6b) set + share λ=0.5", dict(reward="share", lam=0.5)),
    ("(6b) set + share λ=1", dict(reward="share", lam=1.0)),
    ("(6c) DM: log f − log π (FlowRL-style)", dict(reward="dm")),
    ("(6d) set + DM λ=0.1", dict(reward="set+dm", lam=0.1)),
    ("(6d) set + DM λ=0.2", dict(reward="set+dm", lam=0.2)),
    ("(6d) set + DM λ=0.5", dict(reward="set+dm", lam=0.5)),
    ("(6e) set + DM-emp λ=0.2 (π̂ = n/G)", dict(reward="set+dm_emp", lam=0.2)),
    ("(6e) set + DM-emp λ=0.5 (π̂ = n/G)", dict(reward="set+dm_emp", lam=0.5)),
    ("(6f) set + DM-loo λ=0.2 (π̂ = (n−1)/(G−1))", dict(reward="set+dm_loo", lam=0.2)),
    ("(6f) set + DM-loo λ=0.5 (π̂ = (n−1)/(G−1))", dict(reward="set+dm_loo", lam=0.5)),
]
KEY = ["(1) set (any rater)", "(2) freq (share of raters)", "(4) set + coverage λ=0.5",
       "(5) set + KL(SFT) β=0.04", "(6b) set + share λ=0.5", "(6d) set + DM λ=0.2",
       "(6e) set + DM-emp λ=0.5 (π̂ = n/G)"]


def a1():
    rows, curves, pis = [], {}, {}
    for name, kw in VARIANTS:
        o = gt.run(gt.Cfg(**kw), S=S, T=T, seed=1)
        r = summarize(o)
        rows.append((name, r))
        curves[name] = (o["steps"], r)
        pis[name] = o["pi"]
    table(f"A1. Topology-only reward, {S} seeds, {T} visits/prompt, mixed prompts (mean ± 95% CI)", rows)
    save("A1", rows)
    # per-prompt breakdown for key variants
    MD.append("\n#### A1 per-prompt final distribution on the valid answers (mean over seeds)\n")
    MD.append("| prompt (rater freq) | " + " | ".join(KEY) + " |")
    MD.append("|---" * (len(KEY) + 1) + "|")
    names = list(gt.SUITE)
    for p, n in enumerate(names):
        f = gt.hist(gt.SUITE[n])
        va = np.where(f > 0)[0]
        lab = ", ".join(f"{gt.ANSWERS[k][0]}{gt.ANSWERS[k][1]}:{f[k]:.2f}" for k in va)
        cells = []
        for v in KEY:
            pm = pis[v][:, p].mean(0)
            cells.append(" / ".join(f"{pm[k]:.2f}" for k in va))
        MD.append(f"| {n} ({lab}) | " + " | ".join(cells) + " |")
    MD.append("\n#### A1 per-prompt final TV to raters (mean ± 95% CI over seeds) — exposes per-seed collapse "
              "that seed-averaged distributions hide on symmetric prompts\n")
    MD.append("| prompt | " + " | ".join(KEY) + " |")
    MD.append("|---" * (len(KEY) + 1) + "|")
    for p, n in enumerate(names):
        f = gt.hist(gt.SUITE[n])
        cells = [fmt(0.5 * np.abs(pis[v][:, p] - f).sum(-1), 2) for v in KEY]
        MD.append(f"| {n} | " + " | ".join(cells) + " |")
    return curves, pis


def plot_a1(curves, pis):
    style()
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.0))
    for i, v in enumerate(KEY):
        steps, r = curves[v]
        for j, (k, lab) in enumerate([("tv_curve", "TV(π, raters)"), ("inv_curve", "invalid mass"),
                                       ("ent_curve", "H(π) − H(raters) [nats]")]):
            m = r[k].mean(1)
            h = stats.t.ppf(0.975, S - 1) * r[k].std(1, ddof=1) / np.sqrt(S)
            ax[j].plot(steps, m, color=PALETTE[i], lw=2, label=v)
            ax[j].fill_between(steps, m - h, m + h, color=PALETTE[i], alpha=0.15, lw=0)
            ax[j].set_title(lab)
            ax[j].set_xlabel("GRPO visits per prompt")
    ax[1].set_yscale("log")
    ax[2].axhline(0, color="#888", lw=1)
    h_, l_ = ax[0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=8.5, frameon=False, ncol=4, loc="lower center")
    fig.suptitle("A1 · topology-only reward, mixed-rater prompts, mean ± 95% CI over 20 seeds")
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig("figs/a1_curves.png", dpi=130)
    plt.close(fig)

    names = list(gt.SUITE)
    show = ["hole 1/4", "3-way 2/1/1", "hole 2/4"]
    fig, ax = plt.subplots(1, len(show), figsize=(15, 4.2))
    labs = ["raters"] + [k.split(" ", 1)[1][:22] for k in KEY]
    for j, n in enumerate(show):
        p = names.index(n)
        f = gt.hist(gt.SUITE[n])
        va = list(np.where(f > 0)[0])
        mats = [f] + [pis[v][:, p].mean(0) for v in KEY]
        bottom = np.zeros(len(mats))
        for c, k in enumerate(va + ["inv"]):
            vals = np.array([m[k] if k != "inv" else m[f == 0].sum() for m in mats])
            lab = "invalid" if k == "inv" else f"(β0,β1)={gt.ANSWERS[k]}"
            ax[j].bar(range(len(mats)), vals, bottom=bottom, color=PALETTE[c] if k != "inv" else "#bbb",
                      edgecolor="white", lw=1.5, label=lab)
            bottom += vals
        ax[j].set_xticks(range(len(mats)))
        ax[j].set_xticklabels(labs, rotation=35, ha="right", fontsize=8)
        ax[j].set_title(f"prompt '{n}'")
        ax[j].legend(fontsize=8, frameon=False, loc="lower right")
    fig.suptitle(f"A1 · final topology distribution after {T} visits (mean of 20 seeds)")
    fig.tight_layout()
    fig.savefig("figs/a1_final_dists.png", dpi=130)
    plt.close(fig)


def a2_long():
    rows = []
    for v in KEY:
        kw = dict(VARIANTS)[v]
        o = gt.run(gt.Cfg(**kw), S=S, T=1500, seed=2, log_every=25)
        rows.append((v + " @1500", summarize(o)))
    table(f"A2. Long horizon: 1500 visits/prompt ({S} seeds)", rows)
    save("A2", rows)


COMBS = ["sum", "sum-mean", "gdpo"]
SDS = [0.02, 0.05]
GAPS = [0.0, 0.03, 0.1]
COLS11 = ["#2a78d6", "#7aa9e6", "#eb6834", "#1baf7a", "#eda100", "#c98500", "#e87ba4", "#008300", "#4a3aa7",
          "#e34948", "#555555"]
IOU_FIXES = [("(1) set", dict(reward="set")), ("(1b) near (repo default)", dict(reward="near")),
             ("(2) freq", dict(reward="freq")), ("(4) set + coverage λ=0.5", dict(reward="cov", lam=0.5)),
             ("(5) set + KL β=0.04", dict(reward="set", beta=0.04)),
             ("(5) set + KL β=0.1", dict(reward="set", beta=0.1)),
             ("(6b) set + share λ=0.5", dict(reward="share", lam=0.5)),
             ("(6d) set + DM λ=0.2", dict(reward="set+dm", lam=0.2)),
             ("(6d) set + DM λ=0.5", dict(reward="set+dm", lam=0.5)),
             ("(6e) set + DM-emp λ=0.5", dict(reward="set+dm_emp", lam=0.5)),
             ("(6f) set + DM-loo λ=0.5", dict(reward="set+dm_loo", lam=0.5))]


def iou_table(gap, base=0.85, inv_drop=0.02):
    """mean IoU per (prompt, answer): minority / hole answers are `gap` harder to draw."""
    mu = np.full((len(gt.SUITE), gt.K), base - inv_drop)
    for p, n in enumerate(gt.SUITE):
        f = gt.hist(gt.SUITE[n])
        for k in np.where(f > 0)[0]:
            harder = gt.ANSWERS[k] != (1, 0)
            mu[p, k] = base - (gap if harder else 0)
    return mu


def a3_iou():
    grid = []
    rows_all = []
    for combine in COMBS:
        for sd in SDS:
            for gap in GAPS:
                for name, kw in IOU_FIXES:
                    kw = dict(kw)
                    if combine == "sum-mean":
                        kw.update(adv="mean", lr=1.0)
                    cfg = gt.Cfg(combine=combine.split("-")[0], iou_mu=iou_table(gap), iou_sd=sd, **kw)
                    o = gt.run(cfg, S=S, T=T, seed=3)
                    r = summarize(o)
                    grid.append(dict(combine=combine, sd=sd, gap=gap, variant=name,
                                     tv=ci(r["tv"]), inv=ci(r["invalid"]), minret=ci(r["minret"])))
                    rows_all.append((f"{combine} σ={sd} gap={gap} · {name}", r))
    OUT["A3"] = grid
    MD.append(f"\n### A3. Topology term combined with a noisy continuous IoU/band term ({S} seeds, {T} visits)\n")
    MD.append("`sum`: r = 0.5·IoU + 0.5·topo then GRPO z-score. `gdpo`: the repo's group_advantages "
              "(z(IoU)·0.3 + z(band)·0.3 + z(topo·[IoU≥0.5])·0.4, re-z-scored). IoU ~ N(μ_answer, σ), "
              "μ = 0.85 for valid answers, 0.83 for invalid, minus `gap` for non-(1,0) valid answers.\n")
    MD.append("Cells: final TV to raters (mean ± 95% CI) / invalid mass.\n")
    cols = [(c, s, g) for c in COMBS for s in SDS for g in GAPS]
    for c in COMBS:
        sub = [x for x in cols if x[0] == c]
        MD.append(f"\n**combine = {c}**\n")
        MD.append("| variant | " + " | ".join(f"σ={s} gap={g}" for _, s, g in sub) + " |")
        MD.append("|---" * (len(sub) + 1) + "|")
        for name, _ in IOU_FIXES:
            cells = []
            for c_, s, g in sub:
                e = next(x for x in grid if x["combine"] == c_ and x["sd"] == s and x["gap"] == g and x["variant"] == name)
                cells.append(f"{e['tv'][0]:.2f}±{e['tv'][1]:.2f} / {e['inv'][0]:.3f}")
            MD.append(f"| {name} | " + " | ".join(cells) + " |")
    # plot
    style()
    fig, ax = plt.subplots(len(COMBS), 1, figsize=(13, 4.0 * len(COMBS)), sharey=True)
    for j, c in enumerate(COMBS):
        conds = [(s, g) for s in SDS for g in GAPS]
        w = 0.8 / len(IOU_FIXES)
        for i, (name, _) in enumerate(IOU_FIXES):
            m = [next(x for x in grid if x["combine"] == c and x["sd"] == s and x["gap"] == g and x["variant"] == name)["tv"]
                 for s, g in conds]
            ax[j].bar(np.arange(len(conds)) + (i - len(IOU_FIXES) / 2 + 0.5) * w, [x[0] for x in m], w * 0.9,
                      yerr=[x[1] for x in m], color=COLS11[i],
                      label=name, error_kw=dict(lw=0.8))
        ax[j].set_xticks(range(len(conds)))
        ax[j].set_xticklabels([f"σ={s}\ngap={g}" for s, g in conds])
        ax[j].set_title(f"combine = {c}")
        ax[j].set_ylabel("final TV(π, raters)")
    ax[0].legend(fontsize=7.5, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.suptitle("A3 · topology diversity when combined with a noisy IoU term (mean ± 95% CI, 20 seeds)")
    fig.tight_layout()
    fig.savefig("figs/a3_iou.png", dpi=130)
    plt.close(fig)


def a4_shared():
    """Population with a shared bias parameter: 75% of prompts unanimous (1,0)."""
    suite = {}
    for i in range(24):
        suite[f"u{i}"] = {(1, 0): 4}
    for i in range(2):
        for n in MIXED:
            suite[f"{n}#{i}"] = gt.SUITE[n]
    mixed = [n for n in suite if not n.startswith("u")]
    rows = []
    for name, kw in [("(1) set", dict(reward="set")), ("(4) set + coverage λ=0.5", dict(reward="cov", lam=0.5)),
                     ("(5) set + KL β=0.04", dict(reward="set", beta=0.04)),
                     ("(6b) set + share λ=0.5", dict(reward="share", lam=0.5)),
                     ("(6d) set + DM λ=0.2", dict(reward="set+dm", lam=0.2)),
                     ("(6e) set + DM-emp λ=0.5", dict(reward="set+dm_emp", lam=0.5))]:
        for slr in [0.0, 0.3, 1.0]:
            cfg = gt.Cfg(shared_lr=slr, **kw)
            o = gt.run(cfg, suite=suite, S=S, T=T, seed=4)
            rows.append((f"{name}, shared-bias lr={slr}", summarize(o, mixed)))
    table("A4. Shared parameter (generalisation): 24 unanimous + 10 mixed prompts, bias vector shared "
          "across prompts", rows)
    save("A4", rows)


def a5_sensitivity():
    rows = []
    base = dict(reward="set")
    fixes = {"set": dict(reward="set"), "set+DM λ=0.2": dict(reward="set+dm", lam=0.2),
             "set+share λ=0.5": dict(reward="share", lam=0.5), "set+KL β=0.04": dict(reward="set", beta=0.04),
             "set+cov λ=0.5": dict(reward="cov", lam=0.5),
             "set+DM-emp λ=0.5": dict(reward="set+dm_emp", lam=0.5),
             "set+DM-loo λ=0.5": dict(reward="set+dm_loo", lam=0.5)}
    sweeps = [("lr", [0.1, 1.0]), ("G", [4, 16]), ("mu", [1, 4]), ("adv", ["mean"]),
              ("sharpen", [2.0]), ("invalid_mass", [0.05, 0.3]), ("kl", ["exact"])]
    for fname, fkw in fixes.items():
        rows.append((f"{fname} · default (lr 0.5, G 8, μ 2, std-norm, calibrated SFT, 15% invalid)",
                     summarize(gt.run(gt.Cfg(**fkw), S=S, T=T, seed=5))))
        for p, vals in sweeps:
            for v in vals:
                if p == "kl" and "KL" not in fname:
                    continue
                kw = dict(fkw)
                runkw = {}
                if p in ("sharpen", "invalid_mass"):
                    runkw[p] = v
                else:
                    kw[p] = v
                if p == "adv" and v == "mean":
                    kw["lr"] = 0.5 * 2.0  # rescale: mean-only advantages are ~half the size for binary rewards
                o = gt.run(gt.Cfg(**kw), S=S, T=T, seed=5, **runkw)
                rows.append((f"{fname} · {p}={v}", summarize(o)))
    table(f"A5. Sensitivity ({S} seeds, {T} visits). 'sharpen=2' = SFT under-produces minority modes "
          "(p ∝ f²); KL reference is always the SFT policy", rows)
    save("A5", rows)


if __name__ == "__main__":
    t0 = time.time()
    curves, pis = a1()
    plot_a1(curves, pis)
    a2_long()
    a3_iou()
    a4_shared()
    a5_sensitivity()
    with open("results/a_results.json", "w") as fh:
        json.dump(OUT, fh, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else str(x))
    with open("results/a_tables.md", "w") as fh:
        fh.write("\n".join(MD) + "\n")
    print(f"done in {time.time() - t0:.0f}s")
