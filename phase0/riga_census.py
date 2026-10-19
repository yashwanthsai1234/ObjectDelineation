"""Phase 0 census on Segmentation-Ready-RIGA masks (up to 6 experts per image).

Mask values: 0 background, 127 disc, 254 cup (cup drawn over disc). Masks are
filled convex hulls of the experts' outlines, so disc and cup are convex by
construction; see phase0/PLAN.md.

    .venv/bin/python -m phase0.riga_census data/raw/riga/RIGA_preprocessed phase0/results/riga
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image

from checker import AnnotatorSet, certify

MIN_SIZES = (0, 100)
STRUCTURES = ("disc", "cup", "rim")
NAME = re.compile(r"image(\d+)-(\d+)\.png$")


def structures(path):
    a = np.asarray(Image.open(path))[..., 0]
    disc = a >= 100
    cup = a >= 200
    return {"disc": disc, "cup": cup, "rim": disc & ~cup}


def vertical_extent(m):
    rows = np.nonzero(m.any(axis=1))[0]
    return int(rows[-1] - rows[0] + 1) if rows.size else 0


def per_image(job):
    key, files = job
    loaded = {e: structures(f) for e, f in files}
    experts = sorted(loaded)
    out_masks, out_image = [], {"image": key, "experts": len(experts)}
    for e in experts:
        s = loaded[e]
        row = {"image": key, "expert": e,
               "vcdr": vertical_extent(s["cup"]) / max(1, vertical_extent(s["disc"]))}
        for name in STRUCTURES:
            for ms in MIN_SIZES:
                r = certify(s[name], min_size=ms)
                row[f"{name}_b0_{ms}"], row[f"{name}_b1_{ms}"] = r.b0, r.b1
                row[f"{name}_cert_{ms}"] = r.certified
                if not r.certified:
                    row.setdefault("failures", []).extend(r.failures)
            row[f"{name}_area"] = int(s[name].sum())
        out_masks.append(row)
    vc = [r["vcdr"] for r in out_masks]
    out_image["vcdr_min"], out_image["vcdr_max"] = min(vc), max(vc)
    for name in STRUCTURES:
        for ms in MIN_SIZES:
            aset = AnnotatorSet([loaded[e][name] for e in experts], min_size=ms)
            out_image[f"{name}_topos_{ms}"] = ";".join(f"{b0}/{b1}x{n}" for (b0, b1), n in sorted(aset.hist.items()))
            out_image[f"{name}_agree_{ms}"] = aset.topology_agrees
            if ms == 0:
                out_image[f"{name}_band"] = aset.band_width
        ious = []
        for a, b in combinations(experts, 2):
            ma, mb = loaded[a][name], loaded[b][name]
            u = np.logical_or(ma, mb).sum()
            ious.append(np.logical_and(ma, mb).sum() / u if u else 1.0)
        out_image[f"{name}_iou_mean"] = float(np.mean(ious)) if ious else 1.0
        out_image[f"{name}_iou_min"] = float(np.min(ious)) if ious else 1.0
    return out_masks, out_image


def collect(root):
    groups = defaultdict(list)
    for f in sorted(Path(root).rglob("*.png")):
        m = NAME.search(f.name)
        if m:
            subset = str(f.parent.relative_to(root))
            groups[f"{subset}/image{m.group(1)}"].append((int(m.group(2)), f))
    return sorted(groups.items())


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (";".join(v) if isinstance(v, list) else v) for k, v in r.items()})


def summarize(masks, images):
    s = {"images": len(images), "masks": len(masks),
         "experts_per_image": dict(Counter(i["experts"] for i in images)),
         "certificate_failures": sum(1 for m in masks for n in STRUCTURES for ms in MIN_SIZES
                                     if not m[f"{n}_cert_{ms}"])}
    for n in STRUCTURES:
        for ms in MIN_SIZES:
            s[f"{n}_mask_topology_{ms}"] = dict(Counter(f"{m[f'{n}_b0_{ms}']}/{m[f'{n}_b1_{ms}']}" for m in masks))
            s[f"{n}_images_disagree_{ms}"] = sum(1 for i in images if not i[f"{n}_agree_{ms}"])
        bw = np.array([i[f"{n}_band"] for i in images])
        iou = np.array([i[f"{n}_iou_mean"] for i in images])
        s[f"{n}_band_width"] = {"median": float(np.median(bw)), "p10": float(np.percentile(bw, 10)),
                                "p90": float(np.percentile(bw, 90))}
        s[f"{n}_pairwise_iou_mean"] = {"median": float(np.median(iou)), "p10": float(np.percentile(iou, 10))}
    spread = np.array([i["vcdr_max"] - i["vcdr_min"] for i in images])
    s["vcdr_range"] = {"median": float(np.median(spread)), "p90": float(np.percentile(spread, 90)),
                       "share_ge_0.2": float((spread >= 0.2).mean())}
    return s


def main(root, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    jobs = collect(root)
    with Pool() as pool:
        results = pool.map(per_image, jobs, chunksize=4)
    masks = [m for ms, _ in results for m in ms]
    images = [i for _, i in results]
    write_csv(out / "masks.csv", masks)
    write_csv(out / "images.csv", images)
    summary = summarize(masks, images)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
