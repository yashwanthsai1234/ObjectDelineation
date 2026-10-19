"""Phase 0b step 1: RIGA experts' raw outlines -> true disc/cup masks, no convex hull.

An expert's drawing is the difference between the annotated copy and the
original image. A clean drawing is two nested closed curves. The checker finds
closed curves as enclosed background regions; masks are filled from the curves.

Each drawing is also filled twice: as drawn, and with gaps of up to 2 px bridged
(binary closing), to separate pen gaps from JPEG/extraction noise.

    .venv/bin/python -m phase0.riga_outlines data/raw/riga data/riga_true phase0/results/riga_outlines
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

from checker import AnnotatorSet, betti, certify
from checker.topology import BG4, FG8

SUBSETS = ["MESSIDOR", "Magrabia/MagrabiFemale", "Magrabia/MagrabiaMale",
           "BinRushed/BinRushed2", "BinRushed/BinRushed3", "BinRushed/BinRushed4"]
PRIME = re.compile(r"^image(\d+)prime\.(tif|jpg)$", re.I)
ANNO = re.compile(r"^image(\d+)-([1-6])\.(tif|jpg)$", re.I)
MIN_STROKE = 15
MIN_SIZE = 100


def drawing(prime, anno, jpeg):
    d = np.abs(anno.astype(np.int16) - prime.astype(np.int16)).sum(axis=2)
    o = d > (60 if jpeg else 0)
    lab, n = ndimage.label(o, structure=FG8)
    if n:
        sizes = np.bincount(lab.ravel())
        keep = sizes >= MIN_STROKE
        keep[0] = False
        # Not pen strokes: regions touching the image border, and solid blobs
        # (a pen curve covers a small share of the area it encloses).
        border = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
        keep[border] = False
        for k, sl in enumerate(ndimage.find_objects(lab), 1):
            if keep[k]:
                stroke = lab[sl] == k
                if stroke.sum() > 0.5 * ndimage.binary_fill_holes(stroke).sum():
                    keep[k] = False
        o = keep[lab]
    return o


def curves(o):
    """Closed curves of a drawing: (filled mask, stroke pixels) per stroke that encloses something."""
    lab, n = ndimage.label(o, structure=FG8)
    out, open_strokes = [], 0
    for k, sl in enumerate(ndimage.find_objects(lab), 1):
        stroke = lab[sl] == k
        filled = ndimage.binary_fill_holes(stroke, structure=BG4)
        if filled.sum() > stroke.sum():
            full = np.zeros(o.shape, bool)
            full[sl] = filled
            out.append(full)
        else:
            open_strokes += 1
    out.sort(key=lambda m: -m.sum())
    return out, open_strokes, n


def interpret(o):
    closed, open_strokes, strokes = curves(o)
    status = "clean"
    disc = cup = None
    if closed:
        disc = closed[0]
        inner = [c for c in closed[1:] if np.all(disc[c])]
        if inner:
            cup = inner[0]
    if len(closed) == 1:
        # One stroke, but disc and cup may touch: enclosed regions of the drawing tell.
        _, nreg = ndimage.label(~np.pad(o, 1), structure=BG4)
        status = "touching" if nreg - 1 >= 2 else "curve_open"
    elif not closed:
        status = "curve_open"
    elif cup is None:
        status = "curve_open"
    if status == "clean" and (open_strokes or len(closed) > 2):
        status = "extra_strokes"
    return disc, cup, status, strokes, open_strokes, len(closed)


def per_image(job):
    subset, idx, prime_path, annos, out_dir = job
    jpeg = prime_path.suffix.lower() == ".jpg"
    prime = np.asarray(Image.open(prime_path).convert("RGB"))
    rows, masks, drawings = [], defaultdict(dict), {}
    for expert, path in sorted(annos):
        anno = np.asarray(Image.open(path).convert("RGB"))
        if anno.shape != prime.shape:
            rows.append({"subset": subset, "image": idx, "expert": expert, "status": "shape_mismatch"})
            continue
        drawings[expert] = drawing(prime, anno, jpeg)
    # Work inside the box around all large strokes; topology is unchanged by the crop.
    full_shape = prime.shape[:2]
    y0, x0, y1, x1 = 0, 0, full_shape[0], full_shape[1]
    if drawings:
        union = np.any(list(drawings.values()), axis=0)
        lab, n = ndimage.label(union, structure=FG8)
        if n:
            big = np.bincount(lab.ravel()) >= 200
            big[0] = False
            ys, xs = np.nonzero(big[lab])
            if ys.size:
                y0, y1 = max(0, ys.min() - 30), min(full_shape[0], ys.max() + 31)
                x0, x1 = max(0, xs.min() - 30), min(full_shape[1], xs.max() + 31)
    for expert, o_full in drawings.items():
        for variant in ("drawn", "bridged"):
            o = o_full[y0:y1, x0:x1]
            if variant == "bridged":
                o = ndimage.binary_closing(o, structure=np.ones((3, 3), bool), iterations=1)
            disc, cup, status, strokes, open_strokes, nclosed = interpret(o)
            row = {"subset": subset, "image": idx, "expert": expert, "variant": variant, "status": status,
                   "strokes": strokes, "open_strokes": open_strokes, "closed_curves": nclosed}
            if disc is not None and cup is not None:
                rim = disc & ~cup
                for name, m in (("disc", disc), ("cup", cup), ("rim", rim)):
                    r = certify(m, min_size=MIN_SIZE)
                    row[f"{name}_b0"], row[f"{name}_b1"], row[f"{name}_cert"] = r.b0, r.b1, r.certified
                    row[f"{name}_area"] = int(m.sum())
                masks[variant][expert] = (disc, cup)
            rows.append(row)
    image = {"subset": subset, "image": idx, "crop": f"{y0},{x0},{y1},{x1}"}
    for variant, ms in masks.items():
        if len(ms) < 2:
            continue
        for k, name in enumerate(("disc", "cup")):
            aset = AnnotatorSet([v[k] for v in ms.values()], min_size=MIN_SIZE)
            image[f"{variant}_{name}_experts"] = aset.k
            image[f"{variant}_{name}_agree"] = aset.topology_agrees
            image[f"{variant}_{name}_band"] = aset.band_width
        aset = AnnotatorSet([d & ~c for d, c in ms.values()], min_size=MIN_SIZE)
        image[f"{variant}_rim_agree"] = aset.topology_agrees
        image[f"{variant}_rim_topos"] = ";".join(f"{a}/{b}x{n}" for (a, b), n in sorted(aset.hist.items()))
        image[f"{variant}_rim_band"] = aset.band_width
    # Save the "drawn" masks for scoring models later: 0 bg, 127 disc, 254 cup.
    for expert, (disc, cup) in masks.get("drawn", {}).items():
        v = np.zeros(full_shape, np.uint8)
        v[y0:y1, x0:x1][disc] = 127
        v[y0:y1, x0:x1][cup] = 254
        Image.fromarray(v).save(out_dir / f"{subset.replace('/', '_')}_image{idx}-{expert}.png", optimize=True)
    return rows, image


def jobs(root, out_dir):
    out = []
    for subset in SUBSETS:
        d = Path(root) / subset
        primes, annos = {}, defaultdict(list)
        for f in d.iterdir():
            if m := PRIME.match(f.name):
                primes[int(m.group(1))] = f
            elif m := ANNO.match(f.name):
                annos[int(m.group(1))].append((int(m.group(2)), f))
        for idx in sorted(primes):
            if annos[idx]:
                out.append((subset, idx, primes[idx], annos[idx], out_dir))
    return out


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main(root, masks_out, out, per_subset=None):
    masks_out, out = Path(masks_out), Path(out)
    masks_out.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    js = jobs(root, masks_out)
    if per_subset:
        js = [j for sub in SUBSETS for j in [j for j in js if j[0] == sub][:int(per_subset)]]
    with Pool(6) as pool:
        res = pool.map(per_image, js, chunksize=2)
    rows = [r for rs, _ in res for r in rs]
    images = [i for _, i in res]
    write_csv(out / "annotations.csv", rows)
    write_csv(out / "images.csv", images)
    s = {"images": len(images), "annotations": sum(1 for r in rows if r.get("variant") == "drawn")}
    for variant in ("drawn", "bridged"):
        vr = [r for r in rows if r.get("variant") == variant]
        s[f"{variant}_status"] = dict(Counter(r["status"] for r in vr))
        s[f"{variant}_status_by_subset"] = {sub: dict(Counter(r["status"] for r in vr if r["subset"] == sub)) for sub in SUBSETS}
        s[f"{variant}_cert_failures"] = sum(1 for r in vr for n in ("disc", "cup", "rim") if r.get(f"{n}_cert") is False)
        for n in ("disc", "cup", "rim"):
            s[f"{variant}_{n}_topology"] = dict(Counter(f"{r[f'{n}_b0']}/{r[f'{n}_b1']}" for r in vr if f"{n}_b0" in r))
        im = [i for i in images if f"{variant}_rim_agree" in i]
        for n in ("disc", "cup", "rim"):
            s[f"{variant}_{n}_images_disagree"] = {"n": sum(1 for i in im if not i[f"{variant}_{n}_agree"]), "of": len(im)}
            bw = [i[f"{variant}_{n}_band"] for i in im]
            s[f"{variant}_{n}_band_median"] = float(np.median(bw)) if bw else None
    (out / "summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:5])
