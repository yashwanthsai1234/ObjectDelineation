"""Build the precision / topology probe eval views (objdel.data.probes) and score sampled completions on them.

    .venv/bin/python -I -m experiments.phase2.probes build [--store data/store] [--out data/eval_views] \
        [--n-test 1000] [--n-train 0] [--only precision|topology]
    .venv/bin/python -I -m experiments.phase2.probes score --probe precision|topology --completions C.jsonl \
        [--root data/eval_views/probe-<probe>] [--out score.json]

build writes data/eval_views/probe-precision (lidc/val.jsonl = probe-test on the "val" half of the LIDC train
patients; lidc/train.jsonl = the probe-train capacity set when --n-train > 0) and data/eval_views/probe-topology
(D20: 120 one-hole items per bin, 160 two-part, 160 simple, on LIDC val patches), and records the v3e fidelity
check (P5.2) of the build seed in the topology manifest. Completions: {"item_id", "samples": [text, ...]} per line
(objdel.policy.sampler.import_completions); several lines of one item (seeds) are pooled. Reads untrusted data: -I.
"""

import argparse
import json
import time
from pathlib import Path

from objdel.core.io import read_jsonl, write_json
from objdel.data import probes as P


def build(a):
    t0 = time.time()
    stats = P.lidc_stats(a.store)
    print(
        f"stats: {len(stats['diams'])} diameters, {len(stats['holes'])} holes, {len(stats['multis'])} multi-part "
        f"({time.time() - t0:.0f} s)",
        flush=True,
    )
    out = Path(a.out)
    if a.only in (None, "precision"):
        for split, n, seed in (("val", a.n_test, 2026), ("train", a.n_train, 2025)):
            if n:
                pt = P.patches(a.store, split, n=200 if split == "val" else 400)
                P.build_precision(out / "probe-precision", pt, stats["diams"], n=n, seed=seed, split=split)
                print(f"precision {split}: {n} items, {len(pt)} patches ({time.time() - t0:.0f} s)", flush=True)
    if a.only in (None, "topology"):
        pt = P.patches(a.store, "val", n=200)
        man = P.build_topology(out / "probe-topology", pt, stats["holes"], stats["multis"], stats["diams"])
        man["sets"]["lidc/val"]["fidelity_seed_2027"] = P.fidelity(stats["holes"], stats["multis"], 2027)
        write_json(out / "probe-topology" / "manifest.json", man)
        print(f"topology: {man['sets']['lidc/val']} ({time.time() - t0:.0f} s)", flush=True)


def score(a):
    rows = read_jsonl(a.completions)
    root = a.root or f"data/eval_views/probe-{a.probe}"
    res = (P.score_precision if a.probe == "precision" else P.score_topology)(rows, root)
    text = json.dumps(res, indent=1, default=float)
    print(text)
    if a.out:
        Path(a.out).write_text(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--store", default="data/store")
    b.add_argument("--out", default="data/eval_views")
    b.add_argument("--n-test", type=int, default=1000)
    b.add_argument("--n-train", type=int, default=0, help="probe-train items (20000 for the capacity run)")
    b.add_argument("--only", choices=("precision", "topology"))
    s = sub.add_parser("score")
    s.add_argument("--probe", choices=("precision", "topology"), required=True)
    s.add_argument("--completions", required=True)
    s.add_argument("--root")
    s.add_argument("--out")
    a = ap.parse_args()
    (build if a.cmd == "build" else score)(a)


if __name__ == "__main__":
    main()
