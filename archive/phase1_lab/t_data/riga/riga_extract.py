"""Fixed RIGA expert-outline extractor (prototype; reuses experiments.phase0.riga_outlines).

Fixes over the Phase 0b extractor:
  F1  BinRushed1 has no prime (unannotated) image: the reference is the per-pixel
      median of the 6 rater copies (outside strokes the copies are byte-identical).
  F2  Magrabia/MagrabiFemale image48's prime is named `image48.tif` (no "prime").
  F3  Low-contrast pen on dark fundus (BinRushed): hysteresis threshold (25/60)
      instead of a single >60 cut, so faint strokes stay connected.
  F4  Open strokes (pen gaps) were removed by the solid-blob filter before the
      bridged variant could close them: bridge (binary closing) BEFORE filtering,
      and apply the blob test only to closed strokes (open blobs: by thickness).
  F5  Disc and cup drawn touching (one stroke, two enclosed regions): split the
      regions; the cup is the region with least contact with the disc exterior.
  F6  A rater file identical to the reference image carries no drawing: flagged
      `blank` (genuinely missing annotation), not silently dropped.

Stage A is the repo algorithm verbatim, so every mask it accepted is reproduced
bit-exactly; later stages only run for raters that stage A could not resolve.

    python -I riga_extract.py <repo>/data/raw/riga <out_dir> [--force-median SUBSET ...] [--limit N]
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO), str(REPO / "src")]

import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402
from scipy import ndimage  # noqa: E402
from skimage.filters import apply_hysteresis_threshold  # noqa: E402
from skimage.morphology import convex_hull_image  # noqa: E402

from experiments.phase0.riga_outlines import MIN_STROKE, drawing, interpret  # noqa: E402
from objdel.checker import certify  # noqa: E402
from objdel.checker.topology import BG4, FG8  # noqa: E402

SUBSETS = ["MESSIDOR", "Magrabia/MagrabiFemale", "Magrabia/MagrabiaMale", "BinRushed/BinRushed1",
           "BinRushed/BinRushed2", "BinRushed/BinRushed3", "BinRushed/BinRushed4"]
PRE = {"MESSIDOR": "MESSIDOR", "Magrabia/MagrabiFemale": "Magrabia/MagrabiaFemale",
       "Magrabia/MagrabiaMale": "Magrabia/MagrabiaMale", "BinRushed/BinRushed1": "BinRushed/BinRushed1-Corrected",
       "BinRushed/BinRushed2": "BinRushed/BinRushed2", "BinRushed/BinRushed3": "BinRushed/BinRushed3",
       "BinRushed/BinRushed4": "BinRushed/BinRushed4"}
PRIME = re.compile(r"^image(\d+)(?:prime)?\.(tif|jpg)$", re.I)  # F2
ANNO = re.compile(r"^image(\d+)-([1-6])\.(tif|jpg)$", re.I)
LOW, HIGH = 25, 60
MIN_SIZE = 100
PRIME_INK_MAX = 1500  # stroke pixels a clean prime may differ from the 6-copy reference
BRIDGES = (0, 1, 2, 3, 4, 6, 8)  # closing iterations: gaps up to ~2*k px


def load(p):
    return np.asarray(Image.open(p).convert("RGB"))


def discover(root, subsets=SUBSETS):
    out = []
    for subset in subsets:
        d = Path(root) / subset
        primes, annos = {}, defaultdict(dict)
        for f in sorted(d.iterdir()):
            if m := ANNO.match(f.name):
                annos[int(m.group(1))][int(m.group(2))] = f
            elif m := PRIME.match(f.name):
                primes[int(m.group(1))] = f
        for idx in sorted(set(primes) | set(annos)):
            out.append({"subset": subset, "image": idx, "prime": primes.get(idx), "annos": dict(sorted(annos[idx].items()))})
    return out


def median_reference(paths):
    """Reference without a prime image (F1).

    The pen replaces R and G by about (45, 49) and keeps B, and fundus pixels are
    red-dominant, so at each pixel the copy with the largest R+G is one that no
    rater drew on (if any). A plain median fails where >= 3 raters' disc
    outlines coincide, which is common.
    """
    stack = np.stack([load(p) for p in paths])
    k = stack[..., :2].astype(np.int16).sum(-1).argmax(0)
    return np.take_along_axis(stack, k[None, ..., None], 0)[0]


def pen_like(img):
    """The BinRushed pen: R and G replaced by about (45, 49), B kept."""
    r, g = img[..., 0].astype(np.int16), img[..., 1].astype(np.int16)
    return (np.abs(r - 45) <= 15) & (np.abs(g - 49) <= 15)


def ink(a, b):
    """Pen strokes of a over b: |diff| > 60 in components of >= 200 px."""
    d = np.abs(a.astype(np.int16) - b.astype(np.int16)).sum(axis=2) > HIGH
    lab, n = ndimage.label(d, structure=FG8)
    if not n:
        return d
    keep = np.bincount(lab.ravel()) >= 200
    keep[0] = False
    return keep[lab]


def drawing2(ref, anno, jpeg, bridge=0):
    """F3 + F4: hysteresis threshold, bridge first, then filter strokes."""
    d = np.abs(anno.astype(np.int16) - ref.astype(np.int16)).sum(axis=2)
    o = apply_hysteresis_threshold(d, LOW, HIGH) if jpeg else d > 0
    # F1b: where every copy is inked (all raters' outlines coincide) the reference
    # is pen-coloured too and the difference vanishes. Near the rater's own
    # stroke, a pixel that is pen-coloured in both counts as stroke.
    both = pen_like(anno) & pen_like(ref)
    if both.any():
        o |= both & ndimage.binary_dilation(o, iterations=15)
    if bridge:
        o = ndimage.binary_closing(o, structure=np.ones((3, 3), bool), iterations=bridge)
    lab, n = ndimage.label(o, structure=FG8)
    if not n:
        return o
    keep = np.bincount(lab.ravel()) >= MIN_STROKE
    keep[0] = False
    border = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
    keep[border] = False
    for k, sl in enumerate(ndimage.find_objects(lab), 1):
        if not keep[k]:
            continue
        stroke = lab[sl] == k
        filled = ndimage.binary_fill_holes(stroke, structure=BG4)
        # A curve encloses at least its own area (original solid-blob test); an
        # open or hairy stroke (tiny ringing holes) is kept if it is pen-thin.
        keep[k] = stroke.sum() <= 0.5 * filled.sum() or \
            ndimage.distance_transform_edt(np.pad(stroke, 1)).max() <= 10
    return keep[lab]


def split_touching(o):
    """F5: one stroke enclosing two regions -> (disc, cup) or None."""
    filled = ndimage.binary_fill_holes(o, structure=BG4)
    bg, n = ndimage.label(~np.pad(o, 1), structure=BG4)
    bg = bg[1:-1, 1:-1]
    outer = set(np.unique(np.concatenate([bg[0], bg[-1], bg[:, 0], bg[:, -1]]))) - {0}
    sizes = {k: int((bg == k).sum()) for k in range(1, n + 1) if k not in outer}
    regions = sorted(sizes, key=lambda k: -sizes[k])
    # Two real regions; pockets where the curves run side by side are ignored.
    if len(regions) < 2 or (len(regions) > 2 and sizes[regions[2]] >= 0.05 * sizes[regions[1]]):
        return None
    regions = regions[:2]
    width = max(2, int(round(2 * ndimage.distance_transform_edt(np.pad(o, 1)).max())))
    exterior = ~filled
    best = None
    for k in regions:
        r = bg == k
        grown = ndimage.binary_dilation(r, iterations=width + 2)
        contact = (grown & exterior).sum() / max(1, (ndimage.binary_dilation(r) & ~r).sum())
        if best is None or contact < best[0]:
            best = (contact, r)
    cup = ndimage.binary_dilation(best[1], iterations=width) & filled
    cup = ndimage.binary_fill_holes(cup, structure=BG4)
    return filled, cup


def validate(disc, cup):
    """Each structure one clean closed curve (b0=1, b1=0, certified), cup strictly inside disc."""
    out = {}
    for name, m in (("disc", disc), ("cup", cup), ("rim", disc & ~cup)):
        # Rim at the repo's min_size: a cup with 1-px diagonal links splits the
        # 4-connected rim hole into slivers that are not drawn features.
        r = certify(m, min_size=MIN_SIZE if name == "rim" else 0)
        out[f"{name}_b0"], out[f"{name}_b1"], out[f"{name}_cert"] = r.b0, r.b1, r.certified
        out[f"{name}_area"] = int(m.sum())
    out["cup_in_disc"] = bool(np.all(disc[cup]))
    out["nested_ok"] = (out["disc_b0"], out["disc_b1"], out["cup_b0"], out["cup_b1"]) == (1, 0, 1, 0) \
        and out["cup_in_disc"] and out["disc_cert"] and out["cup_cert"] and out["rim_cert"]
    out["strict_nesting"] = out["nested_ok"] and (out["rim_b0"], out["rim_b1"]) == (1, 1)
    out["cup_disc_ratio"] = out["cup_area"] / max(1, out["disc_area"])
    return out


def crop_box(drawings, shape, min_comp=200, margin=30):
    y0, x0, y1, x1 = 0, 0, shape[0], shape[1]
    if drawings:
        union = np.any(list(drawings.values()), axis=0)
        lab, n = ndimage.label(union, structure=FG8)
        if n:
            big = np.bincount(lab.ravel()) >= min_comp
            big[0] = False
            ys, xs = np.nonzero(big[lab])
            if ys.size:
                y0, y1 = max(0, ys.min() - margin), min(shape[0], ys.max() + margin + 1)
                x0, x1 = max(0, xs.min() - margin), min(shape[1], xs.max() + margin + 1)
    return y0, x0, y1, x1


def pre_masks(root, subset, idx, expert):
    p = Path(root) / "RIGA_preprocessed" / PRE[subset] / f"image{idx}-{expert}.png"
    if not p.exists():
        return None
    a = np.asarray(Image.open(p))
    a = a[..., 0] if a.ndim == 3 else a
    return a >= 100, a >= 200


def iou(a, b):
    u = np.logical_or(a, b).sum()
    return float(np.logical_and(a, b).sum() / u) if u else 1.0


def per_image(job):
    item, root, out_dir, force_median = job
    subset, idx = item["subset"], item["image"]
    annos = item["annos"]
    use_median = item["prime"] is None or subset in force_median
    ref_kind = "maxrg6" if use_median else "prime"
    alt_ref = median_reference(list(annos.values()))
    prime_ink = None
    if not use_median:
        ref = load(item["prime"])
        # F7: a "prime" that carries a drawing (Magrabia F48: image48.tif is
        # rater 1's annotated copy) must not be the reference. Trigger: a rater
        # file pixel-identical to the prime, confirmed by ink against the
        # 6-copy reference (the max-RG reference alone is not trusted as an
        # audit: it misfires on some Magrabia images).
        twin = [e for e, p in annos.items() if np.array_equal(load(p), ref)]
        if twin:
            prime_ink = int(ink(ref, alt_ref).sum())
            if prime_ink > PRIME_INK_MAX:
                use_median, ref_kind = True, f"maxrg6:prime_is_rater{twin[0]}_copy"
    if use_median:
        ref = alt_ref
    any_jpeg = use_median or item["prime"].suffix.lower() == ".jpg"
    rows, masks = [], {}
    full_shape = ref.shape[:2]
    imgs, stage_a = {}, {}
    for e, p in annos.items():
        a = load(p)
        row = {"subset": subset, "image": idx, "expert": e, "file": p.name, "ref": ref_kind}
        if a.shape != ref.shape:
            row["status"] = "shape_mismatch"
            rows.append(row)
            continue
        if not use_median and np.array_equal(a, ref):  # F6
            row["status"] = "blank"
            rows.append(row)
            continue
        imgs[e] = a
        jpeg = any_jpeg if use_median else item["prime"].suffix.lower() == ".jpg"
        stage_a[e] = drawing(ref, a, jpeg)
        rows.append(row)
    box_a = crop_box(stage_a, full_shape)
    for row in rows:
        e = row["expert"]
        if e not in imgs:
            continue
        box = box_a
        y0, x0, y1, x1 = box
        disc, cup, status, *_ = interpret(stage_a[e][y0:y1, x0:x1])
        stage = "A"
        attempts = [status]
        if disc is None or cup is None:
            refs = [(ref, ref_kind)]
            if not use_median:
                refs.append((alt_ref, "maxrg6"))
            done = False
            for rr, rk in refs:
                for bridge in BRIDGES:
                    jpeg = True if rk.startswith("maxrg6") else any_jpeg
                    full = drawing2(rr, imgs[e], jpeg, bridge)
                    # The stage-A box can miss strokes stage A did not see: own box.
                    box = crop_box({0: full}, full_shape, min_comp=MIN_STROKE, margin=40)
                    y0, x0, y1, x1 = box
                    ob = full[y0:y1, x0:x1]
                    disc, cup, status, *_ = interpret(ob)
                    attempts.append(f"{rk}/b{bridge}:{status}")
                    if disc is not None and cup is not None:
                        stage, done = f"B:{rk}/b{bridge}", True
                        break
                    if status == "touching":
                        got = split_touching(ob)
                        if got is not None:
                            disc, cup = got
                            if validate(disc, cup)["nested_ok"]:
                                stage, status, done = f"C:{rk}/b{bridge}", "touching_split", True
                                break
                        disc = cup = None
                if done:
                    break
        row["stage"], row["status"], row["attempts"] = stage, status, "|".join(attempts)
        if disc is None or cup is None:
            continue
        row.update(validate(disc, cup))
        y0, x0, y1, x1 = box
        pm = pre_masks(root, subset, idx, e)
        if pm is not None:
            pd, pc = pm[0][y0:y1, x0:x1], pm[1][y0:y1, x0:x1]
            row["hull_iou_disc"] = iou(convex_hull_image(disc), pd)
            row["hull_iou_cup"] = iou(convex_hull_image(cup), pc)
            row["pre_outside_crop"] = int(pm[0].sum() - pd.sum())
        masks[e] = (disc, cup, box)
    for e, (disc, cup, (y0, x0, y1, x1)) in masks.items():
        v = np.zeros(full_shape, np.uint8)
        v[y0:y1, x0:x1][disc] = 127
        v[y0:y1, x0:x1][cup] = 254
        Image.fromarray(v).save(Path(out_dir) / "masks" / f"{subset.replace('/', '_')}_image{idx}-{e}.png", optimize=True)
    # Reference image (prime or median) for hashing / storage.
    Image.fromarray(ref).save(Path(out_dir) / "refs" / f"{subset.replace('/', '_')}_image{idx}.png")
    image = {"subset": subset, "image": idx, "ref": ref_kind, "prime_ink": prime_ink, "crop": ",".join(map(str, box_a)),
             "rater_files": len(annos), "masks": len(masks), "h": full_shape[0], "w": full_shape[1]}
    return rows, image


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main(argv):
    root, out = Path(argv[0]), Path(argv[1])
    force_median, limit, subsets = set(), None, list(SUBSETS)
    rest = argv[2:]
    while rest:
        flag = rest.pop(0)
        if flag == "--force-median":
            force_median.add(rest.pop(0))
        elif flag == "--limit":
            limit = int(rest.pop(0))
        elif flag == "--subsets":
            subsets = rest.pop(0).split(",")
    (out / "masks").mkdir(parents=True, exist_ok=True)
    (out / "refs").mkdir(parents=True, exist_ok=True)
    items = discover(root, subsets)
    if limit:
        items = [i for s in subsets for i in [i for i in items if i["subset"] == s][:limit]]
    with Pool(7) as pool:
        res = pool.map(per_image, [(i, root, out, force_median) for i in items], chunksize=1)
    rows = [r for rs, _ in res for r in rs]
    images = [i for _, i in res]
    write_csv(out / "annotations.csv", rows)
    write_csv(out / "images.csv", images)
    s = {"images": len(images), "rater_files": len(rows),
         "masks": sum(i["masks"] for i in images),
         "images_with_6": sum(i["masks"] == 6 for i in images),
         "masks_hist": dict(sorted(Counter(i["masks"] for i in images).items())),
         "status": dict(Counter(r.get("status") for r in rows)),
         "stage": dict(Counter(r.get("stage") for r in rows if r.get("disc_area"))),
         "nested_ok": sum(1 for r in rows if r.get("nested_ok")),
         "strict_nesting": sum(1 for r in rows if r.get("strict_nesting")),
         "by_subset": {}}
    for sub in subsets:
        im = [i for i in images if i["subset"] == sub]
        rs = [r for r in rows if r["subset"] == sub]
        s["by_subset"][sub] = {"images": len(im), "masks": sum(i["masks"] for i in im),
                               "images_lt6": [(i["image"], i["masks"]) for i in im if i["masks"] < 6],
                               "stage": dict(Counter(r.get("stage") for r in rs if r.get("disc_area"))),
                               "not_recovered": [(r["image"], r["expert"], r.get("status"), r.get("attempts")) for r in rs if not r.get("disc_area")]}
    hd = [r["hull_iou_disc"] for r in rows if "hull_iou_disc" in r]
    hc = [r["hull_iou_cup"] for r in rows if "hull_iou_cup" in r]
    if hd:
        s["hull_iou_disc"] = {"n": len(hd), "min": min(hd), "p1": float(np.percentile(hd, 1)), "median": float(np.median(hd))}
        s["hull_iou_cup"] = {"n": len(hc), "min": min(hc), "p1": float(np.percentile(hc, 1)), "median": float(np.median(hc))}
    (out / "summary.json").write_text(json.dumps(s, indent=1, default=str))
    print(json.dumps(s, indent=1, default=str))


if __name__ == "__main__":
    main(sys.argv[1:])
