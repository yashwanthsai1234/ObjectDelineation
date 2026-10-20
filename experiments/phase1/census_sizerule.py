"""Phase 0 topology census re-run under the Phase 1 size rule (plan §4.2, milestone M3).

An item is one slice (LIDC) or one image structure (RIGA) drawn by k >= 2 raters. The raters disagree on
topology under a rule when their per-rater keys differ. Rules:

    p0_raw     global (b0, b1), no filter                        (Phase 0 "raw")
    p0_fixed   global (b0, b1) after a fixed filter: 5 px LIDC, 100 px RIGA (Phase 0 "filtered")
    sig_raw    per-instance signature (holes per component), no filter
    sig_rule   per-instance signature, filter min_area(|union of the raters|)   (Phase 1)

Data (untrusted, hence -I): data/fixtures/{lidc_full, lidc_experts, riga_experts}.npz and the Phase 0b RIGA
masks in data/riga_true (0 background, 127 disc, 254 cup; one PNG per expert).

    .venv/bin/python -I -m experiments.phase1.census_sizerule [--data data]

References:
- Plan §4.2 (size rule, D3): experiments/phase1/PLAN.md; implementation note experiments/phase1/impl/M3_topology.md.
- Phase 0 censuses this re-runs: experiments/phase0/lidc_census.py, riga_census.py and RESULTS.md / RESULTS_0b.md.
- L. Vincent, "Morphological area openings and closings for grey-scale images", Shape in Picture, NATO ASI F 126,
  1994 (the area opening and closing of size_filter). https://www.doi.org/10.1007/978-3-662-03039-4_13
- T. Y. Kong, A. Rosenfeld, "Digital topology: introduction and survey", CVGIP 48(3):357-393, 1989 (b0, b1 under
  8-connected foreground, 4-connected background). https://doi.org/10.1016/0734-189X(89)90147-3
"""
import argparse
import csv
import re
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image

from objdel.checker import betti, min_area, signature, size_filter

RULES = ("p0_raw", "p0_fixed", "sig_raw", "sig_rule")
NAME = re.compile(r"^(.+)-(\d)\.png$")


def item(stack, fixed):
    """Per-rule disagreement flags, masks changed by the rule, and the rule's min_area."""
    a = min_area(int(stack.any(0).sum()))
    keys = {"p0_raw": [betti(m) for m in stack], "p0_fixed": [betti(size_filter(m, fixed)) for m in stack],
            "sig_raw": [signature(m) for m in stack], "sig_rule": [signature(m, a) for m in stack]}
    changed = sum(r != f for r, f in zip(keys["sig_raw"], keys["sig_rule"], strict=True))
    return {r: len(set(keys[r])) > 1 for r in RULES}, changed, a


def riga_image(files):
    planes = [np.asarray(Image.open(f).convert("L")) for f in files]
    disc = np.stack([p >= 100 for p in planes])
    cup = np.stack([p >= 200 for p in planes])
    return {s: item(st, 100) for s, st in (("disc", disc), ("cup", cup), ("rim", disc & ~cup))}


def summarize(label, results):
    n = len(results)
    dis = {r: sum(x[0][r] for x in results) for r in RULES}
    changed = sum(x[1] for x in results)
    areas = [x[2] for x in results]
    return (f"| {label} | {n:,} | " + " | ".join(f"{dis[r]:,} ({100 * dis[r] / max(n, 1):.1f} %)" for r in RULES)
            + f" | {changed:,} | {int(np.median(areas)) if areas else '-'} |")


def main(data="data"):
    data = Path(data)
    rows = []

    meta = {r["key"]: r for r in csv.DictReader(open(data / "fixtures" / "lidc_full_meta.csv"))}
    lidc = {}
    with np.load(data / "fixtures" / "lidc_full.npz") as z:
        for k in z.files:
            lidc[k] = item(z[k].astype(bool), 5)
    rows.append(summarize("LIDC full, all slices", list(lidc.values())))
    clean = {k: v for k, v in lidc.items() if meta[k]["dup_reader_in_cluster"] == "False"}
    for split in ("train", "val", "test"):
        rows.append(summarize(f"LIDC full, {split} (no duplicated reader)",
                              [v for k, v in clean.items() if meta[k]["split"] == split]))
    nodules = defaultdict(list)
    for k, v in lidc.items():
        nodules[(meta[k]["series"], meta[k]["cluster"])].append(v)
    rows.append(summarize("LIDC full, nodules (any slice)", [
        ({r: any(x[0][r] for x in vs) for r in RULES}, sum(x[1] for x in vs), min(x[2] for x in vs))
        for vs in nodules.values()]))

    for name, fixed in (("lidc_experts", 5), ("riga_experts", 100)):
        with np.load(data / "fixtures" / f"{name}.npz") as z:
            rows.append(summarize(f"fixture {name}", [item(z[k].astype(bool), fixed) for k in z.files]))

    groups = defaultdict(list)
    for f in sorted((data / "riga_true").glob("*.png")):
        m = NAME.match(f.name)
        if m:
            groups[m.group(1)].append(f)
    jobs = [fs for fs in groups.values() if len(fs) >= 2]
    with Pool(3) as pool:
        riga = pool.map(riga_image, jobs, chunksize=8)
    for s in ("disc", "cup", "rim"):
        rows.append(summarize(f"RIGA true (Phase 0b), {s}", [r[s] for r in riga]))

    print("| items | n | " + " | ".join(RULES) + " | masks changed by the rule | median min_area |")
    print("|---|---|" + "---|" * len(RULES) + "---|---|")
    print("\n".join(rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", default="data", help="data root holding fixtures/ and riga_true/ (default: data)")
    main(ap.parse_args().data)
