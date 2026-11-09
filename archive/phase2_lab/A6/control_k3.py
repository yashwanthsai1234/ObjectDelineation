"""A6 control: the K3 workload on the fixture caches (M6's K3 set: first 64 RIGA disc+cup pairs, first 200 LIDC
fixture items, ExpertCache.build) vs the prepared caches (prof_reward's set), same generator, same process, so the
gap between M6's K3 numbers and prof_reward's is attributed to data (frame sizes) or to code / machine load.

    .venv/bin/python -I -m experiments.phase2.lab.A6.control_k3 --out <dir>
"""
import argparse
from pathlib import Path

import cv2
import numpy as np

from objdel.reward import ExpertCache
from tests.oracles import reward_harness as H

from .common import dump, env
from .prof_reward import batches, load_prompts, wall


def fixture_prompts():
    items = H.load_fixture_items(Path("data/fixtures"))
    riga = {}
    for ds, k, st in items:
        if ds == "riga":
            riga.setdefault(k.rsplit("__", 1)[0], {})[k.rsplit("__", 1)[1]] = st
    pr = [(("disc", "cup"), {n: ExpertCache.build(v[n]) for n in ("disc", "cup")}) for v in list(riga.values())[:64]]
    pl = [(("nodule",), {"nodule": ExpertCache.build(st)}) for ds, _, st in items if ds == "lidc"][:200]
    return {"riga": pr, "lidc": pl}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    res = {"env": env()}
    sets = {"fixtures": fixture_prompts(), "prepared": load_prompts(64, 200)[0]}
    for name, prompts in sets.items():
        px = [float(np.prod(c[n[0]].union.shape)) for n, c in prompts["riga"]]
        res[name] = {"riga_frame_px_mean": float(np.mean(px))}
        for label, share in (("mixed", 0.4), ("riga", 1.0)):
            res[name][label] = wall(batches(prompts, share, 3, seed=7), None)["median"]
    for rep in range(2):  # interleave again to expose load drift
        for name, prompts in sets.items():
            res[name][f"riga_rep{rep}"] = wall(batches(prompts, 1.0, 2, seed=9), None)["median"]
    res["env_end"] = env()
    print(res, flush=True)
    dump(res, f"{a.out}/control_k3.json")


if __name__ == "__main__":
    main()
