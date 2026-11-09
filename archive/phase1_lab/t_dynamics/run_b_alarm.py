"""Alarm threshold sensitivity: detection rate / delay on true hacks vs fire rate without a hack."""
from multiprocessing import Pool

import numpy as np

import curriculum_sim as cs

SEEDS = 10
SETTINGS = {
    "lenient (win 3, z 2, rise 0.02, persist 1)": dict(win=3, z=2.0, min_rise=0.02, persist=1),
    "medium (win 4, z 3, rise 0.02, persist 2)": dict(win=4, z=3.0, min_rise=0.02, persist=2),
    "default (win 5, z 4, rise 0.03, persist 3)": dict(win=5, z=4.0, min_rise=0.03, persist=3),
    "strict (win 6, z 5, rise 0.05, persist 4)": dict(win=6, z=5.0, min_rise=0.05, persist=4),
}
TRUE = ["S5 + reward hack", "S11 slow reward hack (α_h=0.6)"]


def job(a):
    sc, sm, seed, st = a
    r = cs.run(cs.world_for(sc), sm, seed, alarm=True, alarm_kw=SETTINGS[st])
    upto = r["step"] <= (r["alarm_at"] if r["alarm_at"] is not None else r["step"][-1])
    onset = r["step"][upto][np.argmin(r["test_fail"][upto])]  # held-out optimum before any intervention
    return a, (r["alarm_at"], onset, r["test_fail"][-1])


if __name__ == "__main__":
    jobs = [(sc, sm, s, st) for sc in cs.SCENARIOS for sm in ["uniform", "learnability"] for s in range(SEEDS)
            for st in SETTINGS]
    with Pool(8) as p:
        res = dict(p.map(job, jobs, chunksize=4))
    lines = ["\n### B9. Alarm threshold sensitivity (uniform + learnability samplers, 10 seeds each)\n",
             "| setting | fire rate, fast hack S5 | median delay after held-out optimum (S5) | fire rate, slow hack S11 | "
             "final fail S11 with alarm | fire rate, no-hack scenarios (S0–S4, S6–S10) | worst no-hack scenario |",
             "|---|---|---|---|---|---|---|"]
    for st in SETTINGS:
        def rows(sc):
            return [res[(sc, sm, s, st)] for sm in ["uniform", "learnability"] for s in range(SEEDS)]
        s5 = rows(TRUE[0])
        s11 = rows(TRUE[1])
        fr5 = np.mean([a is not None for a, _, _ in s5])
        delay = np.median([a - o for a, o, _ in s5 if a is not None]) if fr5 > 0 else np.nan
        fr11 = np.mean([a is not None for a, _, _ in s11])
        f11 = np.mean([f for _, _, f in s11])
        fa = {sc: np.mean([a is not None for a, _, _ in rows(sc)]) for sc in cs.SCENARIOS if sc not in TRUE}
        worst = max(fa, key=fa.get)
        lines.append(f"| {st} | {fr5:.2f} | {delay:.0f} | {fr11:.2f} | {f11:.3f} | {np.mean(list(fa.values())):.3f} | "
                     f"{worst.split(' ', 1)[0]} ({fa[worst]:.2f}) |")
    with open("results/b_tables.md", "a") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
