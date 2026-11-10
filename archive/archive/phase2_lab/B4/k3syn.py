"""B4 H-K3SYN: does the synthetic RIGA-like K3 workload time like the real RIGA one?

64 real RIGA K3 prompts (first 64 RL train items, prepared caches) vs 64 synthetic items (synth.corpus, caches
built). Per repetition r (20): the same 32 prompt indices and the same make_jobs seed for both sets; warm serial
batch wall of each, the order alternating with r. RIGA-only and mixed (40 % RIGA + the same 60 % LIDC prompts).
Statistic: median over repetitions of synthetic / real, with a paired percentile bootstrap CI (10,000 resamples of
repetitions); per-rollout p99 over all rollouts of each set. Verdicts against prereg.py (PS-1..PS-3; PS-4 is
PS-1's ratio when this runs in the x86 container; PS-5 compares the printed corpus digest across machines).

    .venv/bin/python -I -m experiments.phase2.lab.B4.k3syn --out <dir> [--reps 20]
"""
import argparse
import time

import cv2
import numpy as np

from experiments.phase2.lab.A6.common import dump, stats_ms

from . import serve
from .common import ROOT, env, first_keys, prompt_sets, warm_wall
from .k3_harness import batch, per_rollout_ms, synthetic_prompts
from .prereg import verdicts


def boot_ci(x, rng, n=10_000) -> tuple[float, float]:
    x = np.asarray(x)
    meds = np.median(x[rng.integers(0, len(x), (n, len(x)))], axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--tag", default="", help="output suffix (a replication run)")
    a = ap.parse_args()
    cv2.setNumThreads(0)
    res = {"env": env()}
    keys = first_keys(64, 200)
    serve.preload(ROOT, keys)
    ps = prompt_sets(keys)
    sets = {"real": ps["riga"]}
    sets["synthetic"], res["synthetic_digest"] = synthetic_prompts()
    for w, share in (("riga", 1.0), ("mixed", 0.4)):
        walls = {"real": [], "synthetic": []}
        cpu = {"real": [], "synthetic": []}
        lat = {"real": [], "synthetic": []}
        chars = {"real": [], "synthetic": []}
        for r in range(a.reps):
            order = ("real", "synthetic") if r % 2 == 0 else ("synthetic", "real")
            for s in order:
                jobs = batch(sets[s], ps["lidc"], share, 500 + r)
                c0 = time.process_time()
                walls[s].append(warm_wall(jobs)[0])
                cpu[s].append(1e3 * (time.process_time() - c0) / 2)  # mean CPU ms of the two runs (secondary)
                if w == "riga":
                    lat[s] += per_rollout_ms(jobs)
                chars[s] += [len(j[0]) for j in jobs]
        ratio = np.array(walls["synthetic"]) / np.array(walls["real"])
        res[w] = {"walls": {s: stats_ms(v) for s, v in walls.items()}, "ratio_median": float(np.median(ratio)),
                  "ratio_ci95": boot_ci(ratio, np.random.default_rng(0)), "ratios": ratio.tolist(),
                  "chars_mean": {s: float(np.mean(v)) for s, v in chars.items()},
                  "cpu_ms": {s: stats_ms(v) for s, v in cpu.items()},
                  "cpu_ratio_median": float(np.median(np.array(cpu["synthetic"]) / np.array(cpu["real"]))),
                  "cpu_ratio_ci95": boot_ci(np.array(cpu["synthetic"]) / np.array(cpu["real"]),
                                            np.random.default_rng(1))}
        if w == "riga":
            res[w]["per_rollout"] = {s: stats_ms(v) for s, v in lat.items()}
        print(w, res[w]["ratio_median"], res[w]["ratio_ci95"], "cpu", res[w]["cpu_ratio_median"], res[w]["cpu_ratio_ci95"], {s: round(v["median"], 1) for s, v in res[w]["walls"].items()},
              flush=True)
    vals = {"PS-1a RIGA-only synthetic/real median ratio (Mac)": res["riga"]["ratio_median"],
            "PS-1b its 95 % bootstrap CI": tuple(res["riga"]["ratio_ci95"]),
            "PS-2 mixed synthetic/real median ratio (Mac)": res["mixed"]["ratio_median"],
            "PS-3 |p99 synthetic - p99 real| per rollout ms": abs(res["riga"]["per_rollout"]["synthetic"]["p99"]
                                                                  - res["riga"]["per_rollout"]["real"]["p99"])}
    if env()["machine"] == "x86_64":
        vals = {"PS-4 RIGA-only ratio in Docker": res["riga"]["ratio_median"]}
    res["verdicts"] = verdicts(vals)
    res["env_end"] = env()
    print(res["verdicts"], res["synthetic_digest"], flush=True)
    dump(res, f"{a.out}/k3syn_{env()['machine']}{a.tag}.json")


if __name__ == "__main__":
    main()
