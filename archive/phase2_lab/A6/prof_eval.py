"""A6 (c): objdel.eval.metrics on K = 8 / 64 samples per item (val items, prepared caches), with the caches passed in
(loaded once) and without (eval_run's path: metrics builds ExpertCache per item); cProfile hot spots.

Samples: the K3 generator (rater masks jittered / shifted / grown / exact, traced to text), K per item, so most
samples are distinct, as at temperature 1.

    .venv/bin/python -I -m experiments.phase2.lab.A6.prof_eval --out <dir>
"""
import argparse
import time

import cv2
import numpy as np

from experiments.phase1.reward_gates import to_parts
from objdel.eval import metrics
from objdel.geometry import SCHEMAS, canonicalize, serialize
from tests.oracles import reward_harness as H

from .common import by_module, caches, dump, env, profile, records, stats_ms


def samples(cs, names, K, rng):
    Hh, W = cs[names[0]].union.shape
    out = []
    for _ in range(K):
        r = int(rng.integers(cs[names[0]].k))
        kind = int(rng.integers(0, 4))
        objs = {}
        for n in names:
            m = cs[n].raters[r]
            m = (H.jitter(m, float(rng.uniform(0.5, 3)), rng, 6.0) if kind == 0 else
                 H.shift(m, *map(int, rng.integers(-3, 4, 2))) if kind == 1 else
                 H.grow(m, float(rng.uniform(-2, 2))) if kind == 2 else m)
            objs[n] = to_parts(m, W, Hh)
        out.append(serialize(canonicalize(objs), names))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--riga", type=int, default=10)
    ap.add_argument("--lidc", type=int, default=60)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    res = {"env": env()}
    rng = np.random.default_rng(3)
    for ds, n in (("riga", a.riga), ("lidc", a.lidc)):
        names = SCHEMAS[ds]
        recs = [r for r in records(ds, "val") if r["rl"]]
        recs = recs[:: max(1, len(recs) // n)][:n]
        items = [(r, caches(r)) for r in recs]
        res[ds] = {"items": len(items)}
        for K in (8, 64):
            texts = [samples(cs, names, K, rng) for _, cs in items]
            res[ds][f"K{K}_distinct_share"] = float(np.mean([len(set(t)) / K for t in texts]))
            for mode in ("caches", "build"):
                ts = []
                for (r, cs), tx in zip(items, texts, strict=True):
                    raters = {nm: cs[nm].raters for nm in names}
                    t0 = time.perf_counter()
                    metrics(tx, raters, r["W"], r["H"], names, caches=cs if mode == "caches" else None)
                    ts.append(1e3 * (time.perf_counter() - t0))
                res[ds][f"K{K}_{mode}_ms_per_item"] = stats_ms(ts)
            _, rows, tot = profile(lambda: [metrics(tx, {nm: cs[nm].raters for nm in names}, r["W"], r["H"], names,
                                                    caches=cs) for (r, cs), tx in zip(items, texts, strict=True)],
                                   top=300)
            res[ds][f"K{K}_cprofile_top"] = [(round(tt / len(items), 3), round(ct / len(items), 3), nc, nm)
                                             for tt, ct, nc, nm in rows[:16]]
            res[ds][f"K{K}_cprofile_by_file"] = {k: round(v / len(items), 3)
                                                 for k, v in list(by_module(rows).items())[:10]}
            print(ds, K, {m: round(res[ds][f"K{K}_{m}_ms_per_item"]["median"], 1) for m in ("caches", "build")},
                  flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/eval.json")


if __name__ == "__main__":
    main()
