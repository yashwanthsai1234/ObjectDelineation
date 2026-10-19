"""Phase 0b baseline A: box-prompted MedSAM, scored by the checker against each annotator set.

Boxes come from the union of the annotators' masks plus a margin, so the model is
told where the object is; topology errors under this prompt are a lower bound.

    .venv/bin/python -I -m experiments.phase0.baseline_medsam riga data/raw/riga data/riga_true experiments/phase0/results/medsam_riga.csv
    .venv/bin/python -I -m experiments.phase0.baseline_medsam lidc data/raw/lidc_xml data/lidc_slices experiments/phase0/results/medsam_lidc.csv
"""
import csv
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image

from experiments.phase0.lidc_census import NS, nodule_clusters
from experiments.phase0.riga_outlines import jobs as riga_jobs
from objdel.models.medsam import MedSAM, box_of
from objdel.reward import AnnotatorSet

TRUE = re.compile(r"^(.+)_image(\d+)-(\d)\.png$")
LUNG_WINDOW = (-1000, 400)


def write_rows(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def score_row(aset, pred):
    s = aset.score(pred)
    return {"iou": round(s.iou, 4), "band": round(s.band, 4), "b0": s.b0, "b1": s.b1,
            "topo_freq": round(s.topo_freq, 4), "expert_topos": ";".join(f"{a}/{b}x{n}" for (a, b), n in sorted(aset.hist.items())),
            "experts_agree": aset.topology_agrees, **s.errors}


def riga(raw_root, true_root, out_csv):
    primes = {(sub.replace("/", "_"), idx): p for sub, idx, p, _, _ in riga_jobs(raw_root, None)}
    experts = defaultdict(dict)
    for f in Path(true_root).glob("*.png"):
        if m := TRUE.match(f.name):
            experts[(m.group(1), int(m.group(2)))][int(m.group(3))] = f
    model = MedSAM()
    rows = []
    for n, key in enumerate(sorted(experts), 1):
        if key not in primes or len(experts[key]) < 2:
            continue
        ex = [np.asarray(Image.open(f)) for _, f in sorted(experts[key].items())]
        disc, cup = [e >= 100 for e in ex], [e >= 200 for e in ex]
        ud, uc = np.any(disc, 0), np.any(cup, 0)
        img = np.asarray(Image.open(primes[key]).convert("RGB"))
        # Crop to the disc box grown by half its size on each side, for resolution.
        x0, y0, x1, y1 = box_of(ud, 0)
        w, h = x1 - x0, y1 - y0
        cy0, cx0 = max(0, y0 - h // 2), max(0, x0 - w // 2)
        cy1, cx1 = min(img.shape[0], y1 + h // 2 + 1), min(img.shape[1], x1 + w // 2 + 1)
        crop = lambda a: a[cy0:cy1, cx0:cx1]
        cd, cc = [crop(d) for d in disc], [crop(c) for c in cup]
        bd = box_of(crop(ud), max(2, int(0.05 * w)))
        bc = box_of(crop(uc), max(2, int(0.05 * (np.ptp(np.nonzero(uc)[1]) + 1))))
        md, mc = model.segment(crop(img), [bd, bc])
        for name, pred, gt in (("disc", md, cd), ("cup", mc, cc),
                               ("rim", md & ~mc, [d & ~c for d, c in zip(cd, cc)])):
            rows.append({"subset": key[0], "image": key[1], "structure": name, "experts": len(gt),
                         **score_row(AnnotatorSet(gt, min_size=100), pred)})
        if n % 50 == 0:
            print(f"riga {n}/{len(experts)}", flush=True)
            write_rows(out_csv, rows)
    write_rows(out_csv, rows)


def lidc(xml_root, slices_root, out_csv, shard="0/1"):
    series_xml = {}
    for p in Path(xml_root).rglob("*.xml"):
        root = ET.parse(p).getroot()
        if root.tag.endswith("LidcReadMessage"):
            uid = root.find("n:ResponseHeader/n:SeriesInstanceUid", NS).text.strip()
            if uid not in series_xml or "correction" in p.name:
                series_xml[uid] = p
    model = MedSAM()
    out_csv = Path(out_csv)
    index, count = map(int, shard.split("/"))
    if count > 1:
        out_csv = out_csv.with_name(f"{out_csv.stem}.{index}{out_csv.suffix}")
    rows = load_rows(out_csv)
    # Series finished by any earlier run or shard are skipped.
    done = {r["series"] for p in out_csv.parent.glob(f"{Path(out_csv).stem.split('.')[0]}*.csv") for r in load_rows(p)}
    slices_root = Path(slices_root)
    files = [f for i, f in enumerate(sorted(slices_root.glob("*.npz"))) if i % count == index and f.stem not in done]
    for n, f in enumerate(files, 1):
        rows += lidc_series(model, f, series_xml[f.stem])
        print(f"lidc shard {index}/{count} series {n}/{len(files)} rows {len(rows)}", flush=True)
        write_rows(out_csv, rows)
    write_rows(out_csv, rows)


def reader_iou(masks, k):
    """Median IoU of each reader against the strict majority of the other readers."""
    out = []
    for i in range(k):
        rest = [m for j, m in enumerate(masks) if j != i]
        vote = 2 * np.sum(rest, axis=0) > len(rest)
        u = np.logical_or(masks[i], vote).sum()
        out.append(np.logical_and(masks[i], vote).sum() / u if u else 1.0)
    return float(np.median(out))


def zoom_window(box, shape, scale, min_side):
    """Square window centred on box [x0, y0, x1, y1], side = max(min_side, scale x box), clipped."""
    x0, y0, x1, y1 = box
    side = max(min_side, scale * max(x1 - x0 + 1, y1 - y0 + 1))
    cx, cy = (x0 + x1) // 2, (y0 + y1) // 2
    h, w = shape
    ys, xs = max(0, cy - side // 2), max(0, cx - side // 2)
    return slice(ys, min(h, ys + side)), slice(xs, min(w, xs + side))


def window(mask, margin):
    """Slices of the mask's bounding box grown by margin, clipped to the image."""
    ys, xs = np.nonzero(mask)
    return (slice(max(0, ys.min() - margin), ys.max() + margin + 1),
            slice(max(0, xs.min() - margin), xs.max() + margin + 1))


def lidc_series(model, f, xml_path):
    series = f.stem
    got = nodule_clusters(xml_path)
    if got is None:
        return []
    _, sessions, zi, clusters, _ = got
    z_of = {k: z for z, k in zi.items()}
    hu = np.load(f)
    rows = []
    for ci, cl in enumerate(clusters):
        if len(sessions) != 4 or len(cl["readers"]) < 2:
            continue
        for k in sorted({k for sl in cl["per_reader"].values() for k in sl}):
            masks = [sl[k] for sl in cl["per_reader"].values() if k in sl]
            key = f"z{z_of[k]:.2f}"
            if len(masks) < 2 or key not in hu:
                continue
            raw = np.any(masks, axis=0)
            if not raw.any():
                continue
            # Score inside small windows: same topology, IoU and band, far less work than 512 x 512.
            w1 = window(raw, 8)
            first = AnnotatorSet([m[w1] for m in masks], min_size=5)
            if not first.union.any():
                continue
            bx = box_of(first.union, 5)
            box = [bx[0] + w1[1].start, bx[1] + w1[0].start, bx[2] + w1[1].start, bx[3] + w1[0].start]
            lo, hi = LUNG_WINDOW
            g = ((np.clip(hu[key].astype(np.float32), lo, hi) - lo) / (hi - lo) * 255).astype(np.uint8)
            # Zoom: MedSAM sees a window 3x the box (at least 64 px), not the whole 512 slice,
            # so small nodules are not lost at its quarter-resolution mask output.
            zy, zx = zoom_window(box, g.shape, scale=3, min_side=64)
            zbox = [box[0] - zx.start, box[1] - zy.start, box[2] - zx.start, box[3] - zy.start]
            pred = np.zeros(g.shape, bool)
            pred[zy, zx] = model.segment(np.stack([g[zy, zx]] * 3, axis=2), [zbox])[0]
            w2 = window(raw | pred, 2)
            aset = AnnotatorSet([m[w2] for m in masks], min_size=5)
            pred = pred[w2]
            rows.append({"series": series, "cluster": ci, "z": z_of[k], "readers_slice": len(masks),
                         "readers_nodule": len(cl["readers"]), "union_px": int(aset.union.sum()),
                         "reader_iou": round(reader_iou(aset.masks, aset.k), 4), **score_row(aset, pred)})
    return rows


def load_rows(path):
    return list(csv.DictReader(open(path))) if Path(path).exists() else []


if __name__ == "__main__":
    {"riga": riga, "lidc": lidc}[sys.argv[1]](*sys.argv[2:6])
