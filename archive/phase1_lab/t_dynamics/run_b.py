"""Question B experiments: typed-error curriculum vs uniform / difficulty filtering / learnability,
plus the train-reward-up & val-typed-error-up alarm. Writes results/b_*.json|md and figs/b_*.png."""
import json
import sys
import time
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import curriculum_sim as cs
from plotstyle import PALETTE, style

SEEDS = 20
SAMPLERS = list(cs.SAMPLERS)
MD = []


def job(args):
    scen, sampler, seed, alarm = args
    w = cs.world_for(scen)
    r = cs.run(w, sampler, seed, alarm=alarm)
    return (scen, sampler, seed, alarm), dict(
        step=r["step"], test_fail=r["test_fail"], test_worst=r["test_worst"], train_fail=r["train_fail"],
        test_typed=r["test_typed"], alarm_at=r["alarm_at"], best=r["best_val_ckpt"], gini=r["visits_gini"])


def ci(x):
    x = np.asarray(x, float)
    h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return x.mean(), h


def steps_to(curve, steps, target):
    hit = np.where(curve <= target)[0]
    return steps[hit[0]] if len(hit) else np.nan


def main():
    scens = list(cs.SCENARIOS)
    jobs = [(sc, sm, s, al) for sc in scens for sm in SAMPLERS for s in range(SEEDS) for al in (False, True)]
    t0 = time.time()
    with Pool(8) as pool:
        res = dict(pool.map(job, jobs, chunksize=4))
    print(f"{len(jobs)} runs in {time.time() - t0:.0f}s", file=sys.stderr)

    summary = {}
    for sc in scens:
        steps = res[(sc, "uniform", 0, False)]["step"]
        uni_best = np.mean([res[(sc, "uniform", s, False)]["test_fail"].min() for s in range(SEEDS)])
        uni_final = np.mean([res[(sc, "uniform", s, False)]["test_fail"][-1] for s in range(SEEDS)])
        target = min(uni_best, uni_final)
        for sm in SAMPLERS:
            R = [res[(sc, sm, s, False)] for s in range(SEEDS)]
            RA = [res[(sc, sm, s, True)] for s in range(SEEDS)]
            summary[(sc, sm)] = dict(
                final=np.array([r["test_fail"][-1] for r in R]),
                auc=np.array([r["test_fail"].mean() for r in R]),
                worst=np.array([r["test_worst"][-1] for r in R]),
                gap=np.array([r["test_fail"][-1] - r["train_fail"][-1] for r in R]),
                t2t=np.array([steps_to(r["test_fail"], steps, target) for r in R]),
                best_ckpt=np.array([r["best"][0] for r in R]),
                best_step=np.array([r["best"][1] for r in R]),
                onset=np.array([steps[np.argmin(r["test_fail"])] for r in R]),
                gini=np.array([r["gini"] for r in R]),
                alarm_rate=np.mean([r["alarm_at"] is not None for r in RA]),
                alarm_at=np.array([r["alarm_at"] if r["alarm_at"] is not None else np.nan for r in RA]),
                final_alarm=np.array([r["test_fail"][-1] for r in RA]),
                curve=np.array([r["test_fail"] for r in R]),
                curve_worst=np.array([r["test_worst"] for r in R]),
                curve_alarm=np.array([r["test_fail"] for r in RA]),
                typed=np.array([r["test_typed"][-1] for r in R]),
            )
        summary[(sc, "_target")] = target
        summary[(sc, "_steps")] = steps

    def cell(x, d=3):
        m, h = ci(x)
        return f"{m:.{d}f} ± {h:.{d}f}"

    def paired(sc, a, b, key):
        d = summary[(sc, a)][key] - summary[(sc, b)][key]
        m, h = ci(d)
        sig = "*" if abs(m) > h else ""
        return f"{m:+.3f} ± {h:.3f}{sig}"

    short = {"uniform": "uniform", "difficulty-filter": "diff-filter", "learnability": "learnability",
             "typed-flat": "typed-flat", "typed-balance": "typed-balance", "typed-stall": "typed-stall",
             "typed-lp": "typed-LP"}
    for key, title, d in [("final", "B1. Final held-out fail rate P(any typed error) at T=600 (mean ± 95% CI, 20 seeds)", 3),
                          ("auc", "B2. Mean held-out fail rate over training (AUC/T; lower = faster)", 3),
                          ("worst", "B3. Final held-out worst-type error rate (max over types of per-type error on relevant items)", 3)]:
        MD.append(f"\n### {title}\n")
        MD.append("| scenario | " + " | ".join(short[s] for s in SAMPLERS) + " |")
        MD.append("|---" * (len(SAMPLERS) + 1) + "|")
        for sc in scens:
            vals = [summary[(sc, s)][key].mean() for s in SAMPLERS]
            best = int(np.argmin(vals))
            cells = []
            for i, s in enumerate(SAMPLERS):
                c = cell(summary[(sc, s)][key], d)
                cells.append(f"**{c}**" if i == best else c)
            MD.append(f"| {sc} | " + " | ".join(cells) + " |")

    MD.append("\n### B4. Steps to reach the uniform sampler's best held-out fail rate (median over seeds; "
              "'n/m' = seeds that never reached it within 600 steps)\n")
    MD.append("| scenario | target | " + " | ".join(short[s] for s in SAMPLERS) + " |")
    MD.append("|---" * (len(SAMPLERS) + 2) + "|")
    for sc in scens:
        cells = []
        for s in SAMPLERS:
            x = summary[(sc, s)]["t2t"]
            miss = np.isnan(x).sum()
            med = np.nanmedian(x) if miss < len(x) else np.nan
            cells.append(f"{med:.0f}" + (f" ({miss}/{len(x)} never)" if miss else ""))
        MD.append(f"| {sc} | {summary[(sc, '_target')]:.3f} | " + " | ".join(cells) + " |")

    MD.append("\n### B5. Paired differences (same seed = same item population): typed variants minus scalar "
              "learnability p(1−p). Negative = typed better; * = 95% CI excludes 0\n")
    MD.append("| scenario | metric | " + " | ".join(short[s] for s in SAMPLERS if s.startswith("typed")) +
              " | uniform − learnability | diff-filter − learnability |")
    MD.append("|---" * (len([s for s in SAMPLERS if s.startswith('typed')]) + 4) + "|")
    for sc in scens:
        for key in ["final", "auc", "worst"]:
            cells = [paired(sc, s, "learnability", key) for s in SAMPLERS if s.startswith("typed")]
            cells += [paired(sc, "uniform", "learnability", key), paired(sc, "difficulty-filter", "learnability", key)]
            MD.append(f"| {sc} | {key} | " + " | ".join(cells) + " |")

    MD.append("\n### B6. Overfitting / hacking alarm (train-reward EMA ↑ while some val typed error's 5-eval mean "
              "rises above its best earlier 5-eval mean by > max(0.03, 4·SE), for 3 consecutive evals). Action on fire: add the missing guard (EP weight 0.02 → 0.2) "
              "and unlearn the hack. Same seeds as B1.\n")
    MD.append("| scenario | sampler | fire rate | median fire step | median step of held-out optimum (no alarm) | "
              "final fail w/o alarm | final fail with alarm | best-val-checkpoint (early stop) |")
    MD.append("|---|---|---|---|---|---|---|---|")
    for sc in scens:
        for s in ["uniform", "learnability", "typed-balance"]:
            S_ = summary[(sc, s)]
            fs = np.nanmedian(S_["alarm_at"]) if not np.isnan(S_["alarm_at"]).all() else np.nan
            MD.append(f"| {sc} | {short[s]} | {S_['alarm_rate']:.2f} | {fs:.0f} | {np.median(S_['onset']):.0f} | "
                      f"{cell(S_['final'])} | "
                      f"{cell(S_['final_alarm'])} | {cell(S_['best_ckpt'])} |")

    MD.append("\n### B7. Generalisation gap (held-out fail − train fail at T) and visit concentration (Gini of "
              "per-item visit counts)\n")
    MD.append("| scenario | " + " | ".join(f"{short[s]} gap / gini" for s in SAMPLERS) + " |")
    MD.append("|---" * (len(SAMPLERS) + 1) + "|")
    for sc in scens:
        MD.append(f"| {sc} | " + " | ".join(
            f"{summary[(sc, s)]['gap'].mean():.3f} / {summary[(sc, s)]['gini'].mean():.2f}" for s in SAMPLERS) + " |")

    MD.append("\n### B8. Final held-out per-type error (relevant items), S1 default and S2 unlearnable EP\n")
    for sc in ["S1 default (heterogeneous)", "S2 + unlearnable EP (floor 0.3)", "S8 rare type also slow (MH α=0.5, q=0.05)",
               "S10 ambiguous type (EP q=0.35, floor 0.45)"]:
        MD.append(f"\n**{sc}**\n")
        MD.append("| sampler | " + " | ".join(cs.TYPES) + " |")
        MD.append("|---" * (cs.K + 1) + "|")
        for s in SAMPLERS:
            t = summary[(sc, s)]["typed"].mean(0)
            MD.append(f"| {short[s]} | " + " | ".join(f"{v:.3f}" for v in t) + " |")

    # --------------- plots
    style()
    show = ["S0 homogeneous types", "S1 default (heterogeneous)", "S2 + unlearnable EP (floor 0.3)",
            "S8 rare type also slow (MH α=0.5, q=0.05)", "S10 ambiguous type (EP q=0.35, floor 0.45)",
            "S11 slow reward hack (α_h=0.6)"]
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))
    for ax, sc in zip(axs.ravel(), show):
        steps = summary[(sc, "_steps")]
        for i, s in enumerate(SAMPLERS):
            c = summary[(sc, s)]["curve"]
            m = c.mean(0)
            h = stats.t.ppf(0.975, SEEDS - 1) * c.std(0, ddof=1) / np.sqrt(SEEDS)
            ax.plot(steps, m, color=PALETTE[i], lw=1.8, label=short[s])
            ax.fill_between(steps, m - h, m + h, color=PALETTE[i], alpha=0.12, lw=0)
        ax.set_yscale("log")
        ax.set_title(sc, fontsize=10)
        ax.set_xlabel("GRPO steps")
        ax.set_ylabel("held-out fail rate")
    axs[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle("B · held-out fail rate by sampler (mean ± 95% CI, 20 seeds)")
    fig.tight_layout()
    fig.savefig("figs/b1_curves.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 4.8))
    w = 0.8 / len(SAMPLERS)
    for i, s in enumerate(SAMPLERS):
        vals = [ci(summary[(sc, s)]["final"] - summary[(sc, "learnability")]["final"]) for sc in scens]
        ax.bar(np.arange(len(scens)) + (i - len(SAMPLERS) / 2 + 0.5) * w, [v[0] for v in vals], w * 0.9,
               yerr=[v[1] for v in vals], color=PALETTE[i], label=short[s], error_kw=dict(lw=0.8))
    ax.axhline(0, color="#888", lw=1)
    ax.set_xticks(range(len(scens)))
    ax.set_xticklabels([sc.split(" ", 1)[0] for sc in scens], fontsize=9)
    ax.set_xlabel("scenario (see table B1 for definitions)")
    ax.set_ylabel("final held-out fail − learnability (paired)")
    ax.legend(fontsize=8, frameon=False, ncol=4)
    ax.set_title("B · sensitivity: each sampler vs scalar learnability p(1−p), paired by seed (mean ± 95% CI)")
    fig.tight_layout()
    fig.savefig("figs/b2_sensitivity.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
    sc = "S5 + reward hack"
    steps = summary[(sc, "_steps")]
    for i, s in enumerate(["uniform", "learnability", "typed-balance"]):
        for j, (k, ls) in enumerate([("curve", "--"), ("curve_alarm", "-")]):
            m = summary[(sc, s)][k].mean(0)
            ax[0].plot(steps, m, color=PALETTE[i], ls=ls, lw=1.8, label=f"{short[s]}{' + alarm' if j else ''}")
    ax[0].set_title("S5 reward hack: held-out fail rate, with/without alarm")
    ax[0].set_xlabel("GRPO steps")
    ax[0].legend(fontsize=8, frameon=False)
    fires = [np.mean([summary[(sc_, s_)]["alarm_rate"] for s_ in SAMPLERS]) for sc_ in scens]
    ax[1].barh(range(len(scens)), fires, color=[PALETTE[7] if "hack" in sc_ else PALETTE[0] for sc_ in scens])
    ax[1].set_yticks(range(len(scens)))
    ax[1].set_yticklabels(scens, fontsize=8)
    ax[1].set_xlabel("fraction of runs where the alarm fired (all 7 samplers × 20 seeds)")
    ax[1].set_title("alarm fire rate: red = true hack, blue = no hack (false alarms)")
    fig.tight_layout()
    fig.savefig("figs/b3_alarm.png", dpi=130)
    plt.close(fig)

    with open("results/b_tables.md", "w") as fh:
        fh.write("\n".join(MD) + "\n")
    js = {f"{k[0]} | {k[1]}": ({kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv) for kk, vv in v.items()
                                 if not kk.startswith("curve")} if isinstance(v, dict) else
                                (v.tolist() if isinstance(v, np.ndarray) else v)) for k, v in summary.items()}
    with open("results/b_results.json", "w") as fh:
        json.dump(js, fh)


if __name__ == "__main__":
    main()
