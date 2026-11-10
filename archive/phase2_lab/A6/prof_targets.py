"""A6 (b): on-the-fly augmentation cost: rebuild certified SFT targets (objdel.targets.build_target) from
flipped / rotated rater masks, per item, plus the GRPO cache rebuild (ExpertCache.build) the same augmentation
would need; cProfile hot spots inside build_target.

Masks: the frame-size rater stacks stored in data/prepared caches (ExpertCache.raters, the same arrays prepare()
built targets from). Check: the identity transform reproduces every stored target text.

    .venv/bin/python -I -m experiments.phase2.lab.A6.prof_targets --out <dir>
"""
import argparse
import time

import numpy as np

from objdel.geometry import SCHEMAS
from objdel.reward import ExpertCache
from objdel.targets import build_target, min_areas

from .common import by_module, caches, dump, env, profile, records, stats_ms

AUGS = {
    "identity": lambda m: m,
    "hflip": lambda m: m[..., :, ::-1],
    "vflip": lambda m: m[..., ::-1, :],
    "rot90": lambda m: np.rot90(m, 1, axes=(-2, -1)),
    "rot180": lambda m: np.rot90(m, 2, axes=(-2, -1)),
    "transpose": lambda m: np.swapaxes(m, -1, -2),
}


def item_targets(stacks, names, aug):
    st = {n: np.ascontiguousarray(AUGS[aug](s)) for n, s in stacks.items()}
    ma = min_areas(st)
    k = len(st[names[0]])
    out = []
    for r in range(k):
        try:
            out.append(build_target({n: st[n][r] for n in names}, names, None, ma))
        except ValueError:
            out.append(None)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--riga", type=int, default=24)
    ap.add_argument("--lidc", type=int, default=300)
    a = ap.parse_args()
    res = {"env": env()}
    for ds, n in (("riga", a.riga), ("lidc", a.lidc)):
        names = SCHEMAS[ds]
        recs = [r for r in records(ds, "train") if r["rl"]]
        recs = recs[:: max(1, len(recs) // n)][:n]
        items = []
        for r in recs:
            c = caches(r)
            items.append((r, {nm: c[nm].raters for nm in names}))
        same = total = 0
        per = {aug: [] for aug in AUGS}
        per_target = {aug: [] for aug in AUGS}
        taus = []
        for r, stacks in items:
            for aug in AUGS:
                t0 = time.perf_counter()
                ts = item_targets(stacks, names, aug)
                per[aug].append(1e3 * (time.perf_counter() - t0))
                per_target[aug] += [1e3 * t.info["time"] for t in ts if t]
                if aug == "identity":
                    for t, ref in zip(ts, r["targets"], strict=True):
                        total += 1
                        same += t is not None and t.text == ref
                else:
                    taus += [t.info[nm]["tau"] for t in ts if t for nm in names]
        res[ds] = {"items": len(items), "raters_per_item": float(np.mean([len(s[names[0]]) for _, s in items])),
                   "identity_reproduces_stored_targets": f"{same}/{total}",
                   "item_ms": {aug: stats_ms(v) for aug, v in per.items()},
                   "target_cpu_ms": {aug: stats_ms(v) for aug, v in per_target.items()},
                   "aug_tau_hist": {str(k): int(v) for k, v in zip(*np.unique([str(x) for x in taus],
                                                                                return_counts=True), strict=True)}}
        _, rows, tot = profile(lambda: [item_targets(s, names, "hflip") for _, s in items], top=300)
        res[ds]["cprofile_hflip_total_ms_per_item"] = tot / len(items)
        res[ds]["cprofile_top"] = [(round(tt / len(items), 3), round(ct / len(items), 3), nc, nm)
                                   for tt, ct, nc, nm in rows[:22]]
        res[ds]["cprofile_by_file_ms_per_item"] = {k: round(v / len(items), 3)
                                                   for k, v in list(by_module(rows).items())[:14]}
        # GRPO would also need the caches of the augmented masks
        cb = []
        for _, stacks in items[: (8 if ds == "riga" else 100)]:
            st = {nm: np.ascontiguousarray(AUGS["hflip"](s)) for nm, s in stacks.items()}
            t0 = time.perf_counter()
            for nm in names:
                ExpertCache.build(st[nm])
            cb.append(1e3 * (time.perf_counter() - t0))
        res[ds]["cache_build_ms_per_item"] = stats_ms(cb)
        print(ds, res[ds]["identity_reproduces_stored_targets"],
              {k: round(v["median"], 1) for k, v in res[ds]["item_ms"].items()},
              "cache build", round(res[ds]["cache_build_ms_per_item"]["median"], 1), flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/targets.json")


if __name__ == "__main__":
    main()
