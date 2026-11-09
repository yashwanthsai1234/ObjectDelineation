"""RIGA expert masks from the raters' pen drawings (staged extractor F1–F7) and near-duplicate features.

RIGA (Almazroa et al. 2018) has 750 fundus photographs, each outlined (optic disc and cup) by six
ophthalmologists who drew on their own copy of the photograph. A rater's drawing is the difference between their
copy and the unannotated reference; a clean drawing is two nested closed curves, and a mask is the curve filled
(stroke included).

Stage A is the Phase 0 extractor unchanged (single threshold, solid-blob filter), so every mask it accepted is
reproduced bit for bit. Later stages run only for raters stage A could not resolve:

- F1  BinRushed1 has no reference photograph: rebuild it from the six copies by per-pixel max R+G.
- F2  Magrabia Female image48's reference is named ``image48.tif`` (no "prime").
- F3  Faint pen on dark fundus: hysteresis threshold (25 / 60) instead of one cut at 60.
- F4  Pen gaps: bridge (binary closing) before the stroke filter; open strokes are kept when pen-thin.
- F5  Disc and cup drawn touching (one stroke, two regions): split the regions.
- F6  A rater copy identical to a clean reference carries no drawing: that rater is missing.
- F7  A "reference" that is itself a rater's annotated copy is replaced by the F1 reconstruction.

Every returned mask pair is certified: disc and cup each one component without holes (bit-quad, Betti and ring
certificates agree), cup inside disc, rim certified.

Model image (B3-1): the BinRushed1 reconstruction (F1) extracts the masks but is not a model image. Where every
copy is inked it keeps the pen colour (18 / 50 images have >= 30 near-black pixels within 2 px of a rater edge), and
elsewhere the per-pixel switch between JPEG copies leaves a seam on the rater edges (green-channel Laplacian on the
2-px edge band over a band 8-14 px away: 1.8-6.0 on 50 / 50 BinRushed1 images, <= 1.5 on the 700 others); OpenCV
inpainting removes the dark pixels but not the seam (12 / 12 tried). BinRushed1 items keep their masks with
``image=None``. The F7 reconstruction (Magrabia Female image48, lossless TIFF copies) measures clean and stays.

References:
- Almazroa et al., "Retinal fundus images for glaucoma analysis: the RIGA dataset", Proc. SPIE 10579 (2018),
  https://doi.org/10.1117/12.2293584 ; data: https://deepblue.lib.umich.edu/data/concern/data_sets/3b591905z
- OpenCV inpainting (Telea fast marching, Navier-Stokes), tried and rejected for the model image:
  https://docs.opencv.org/4.x/df/d3d/tutorial_py_inpainting.html
- ORB (Rublee et al., ICCV 2011) and CLAHE in OpenCV: https://docs.opencv.org/4.x/d1/d89/tutorial_py_orb.html ,
  https://docs.opencv.org/4.x/d5/daf/tutorial_py_histogram_equalization.html
- Hysteresis thresholding as in scikit-image ``apply_hysteresis_threshold`` (strict ``>``, default 4-connected
  ``ndi.label``): https://github.com/scikit-image/scikit-image/blob/v0.25.2/skimage/filters/thresholding.py
- Lab evidence: archive/phase1_lab/t_data/REPORT.md (4,500 / 4,500 masks; 4,172 / 4,172 equal to Phase 0).
"""
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

from objdel.checker import certify_many

from .store import Item

SUBSETS = ("BinRushed/BinRushed1", "BinRushed/BinRushed2", "BinRushed/BinRushed3", "BinRushed/BinRushed4",
           "MESSIDOR", "Magrabia/MagrabiFemale", "Magrabia/MagrabiaMale")
PRIME = re.compile(r"^image(\d+)(?:prime)?\.(tif|jpg)$", re.I)  # F2: "prime" is optional
RATER = re.compile(r"^image(\d+)-([1-6])\.(tif|jpg)$", re.I)
FG8 = np.ones((3, 3), bool)  # 8-connected foreground
BG4 = ndimage.generate_binary_structure(2, 1)  # 4-connected background
LOW, HIGH = 25, 60  # pen-stroke thresholds on the summed |RGB| difference (F3)
MIN_STROKE = 15  # smallest stroke component, pixels
RIM_MIN_SIZE = 100  # rim certified after dropping slivers below this size (cup touching the disc by 1-px links)
PRIME_INK_MAX = 1500  # stroke pixels a clean reference may differ from the F1 reconstruction (F7)
BRIDGES = (0, 1, 2, 3, 4, 6, 8)  # closing iterations tried by stage B: gaps up to about 2k px
ORB_FEATURES, ORB_WIDTH, FIELD_INSET = 500, 800, 20  # keypoints(): features, thumbnail width px, field inset px


@dataclass(frozen=True)
class Photo:
    """One RIGA photograph: its reference file (None for BinRushed1) and the raters' annotated copies."""

    subset: str
    index: int
    prime: Path | None
    raters: dict[int, Path]

    @property
    def image_id(self) -> str:
        """Image id `<subset>/image<index:03d>`."""
        return f"{self.subset}/image{self.index:03d}"


def discover(root) -> list[Photo]:
    """All photographs under the RIGA root, ordered by image id."""
    photos = []
    for subset in SUBSETS:
        primes, raters = {}, {}
        for f in sorted((Path(root) / subset).iterdir()):
            if m := RATER.match(f.name):
                raters.setdefault(int(m.group(1)), {})[int(m.group(2))] = f
            elif m := PRIME.match(f.name):
                primes[int(m.group(1))] = f
        for idx in sorted(set(primes) | set(raters)):
            photos.append(Photo(subset, idx, primes.get(idx), dict(sorted(raters.get(idx, {}).items()))))
    return sorted(photos, key=lambda p: p.image_id)


def load(path) -> np.ndarray:
    """The image file as a uint8 (H, W, 3) RGB array."""
    return np.asarray(Image.open(path).convert("RGB"))


def max_rg(copies) -> np.ndarray:
    """F1: reference without a reference file.

    The pen replaces R and G by about (45, 49) and keeps B, and fundus pixels are red-dominant, so at each pixel the
    copy with the largest R+G is one no rater drew on (if any). A plain median fails where three or more raters'
    outlines coincide, which is common.
    """
    stack = np.stack(copies)
    k = stack[..., :2].astype(np.int16).sum(-1).argmax(0)
    return np.take_along_axis(stack, k[None, ..., None], 0)[0]


def pen_like(img) -> np.ndarray:
    """The BinRushed pen colour: R and G near (45, 49)."""
    r, g = img[..., 0].astype(np.int16), img[..., 1].astype(np.int16)
    return (np.abs(r - 45) <= 15) & (np.abs(g - 49) <= 15)


def hysteresis(d, low, high) -> np.ndarray:
    """Pixels above `low` that are 4-connected to a pixel above `high` (scikit-image semantics)."""
    lab, n = ndimage.label(d > low)
    seeded = np.zeros(n + 1, bool)
    seeded[lab[d > high]] = True
    seeded[0] = False
    return seeded[lab]


def diff(a, b) -> np.ndarray:
    """Summed absolute RGB difference of two images, an (H, W) int array (computed in int16)."""
    return np.abs(a.astype(np.int16) - b.astype(np.int16)).sum(axis=2)


def components(o, min_size) -> tuple[np.ndarray, np.ndarray]:
    """(labels, keep): 8-connected components of `o`; keep[label] is True for those of at least `min_size` px."""
    lab, _ = ndimage.label(o, structure=FG8)
    keep = np.bincount(lab.ravel()) >= min_size
    keep[0] = False
    return lab, keep


def ink(a, b) -> np.ndarray:
    """F7 audit: strokes of `a` over `b` (difference > 60 in components of at least 200 px)."""
    lab, keep = components(diff(a, b) > HIGH, 200)
    return keep[lab]


def drawing(ref, copy, jpeg, bridge=None) -> np.ndarray:
    """A rater's pen strokes: the difference of `copy` over `ref`, filtered to pen-like components.

    bridge=None is stage A (the Phase 0 rule: one threshold, closed strokes must enclose more than their own area).
    An integer is stage B: hysteresis (F3), coincident-outline recovery, bridging by `bridge` closing iterations
    before the filter (F4), and open strokes kept when pen-thin.
    """
    d = diff(copy, ref)
    if bridge is None:
        o = d > (HIGH if jpeg else 0)
    else:
        o = hysteresis(d, LOW, HIGH) if jpeg else d > 0
        # Where every copy is inked the F1 reference is pen-coloured too and the difference vanishes: near the
        # rater's own stroke, a pixel pen-coloured in both counts as stroke.
        both = pen_like(copy) & pen_like(ref)
        if both.any():
            o |= both & ndimage.binary_dilation(o, iterations=15)
        if bridge:
            o = ndimage.binary_closing(o, structure=FG8, iterations=bridge)
    lab, keep = components(o, MIN_STROKE)
    keep[np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))] = False  # not strokes: border regions
    for k, sl in enumerate(ndimage.find_objects(lab), 1):
        if keep[k]:
            stroke = lab[sl] == k
            keep[k] = stroke.sum() <= 0.5 * ndimage.binary_fill_holes(stroke, structure=BG4).sum() or (
                bridge is not None and ndimage.distance_transform_edt(np.pad(stroke, 1)).max() <= 10)
    return keep[lab]


def curves(o) -> list[np.ndarray]:
    """Filled closed curves of a drawing (strokes that enclose something), largest first."""
    lab, _ = ndimage.label(o, structure=FG8)
    out = []
    for k, sl in enumerate(ndimage.find_objects(lab), 1):
        stroke = lab[sl] == k
        filled = ndimage.binary_fill_holes(stroke, structure=BG4)
        if filled.sum() > stroke.sum():
            full = np.zeros(o.shape, bool)
            full[sl] = filled
            out.append(full)
    out.sort(key=lambda m: -m.sum())
    return out


def interpret(o):
    """(disc, cup, touching): the largest closed curve, the first closed curve inside it, and whether a single
    stroke encloses two or more regions (disc and cup drawn touching)."""
    closed = curves(o)
    if not closed:
        return None, None, False
    disc = closed[0]
    cup = next((c for c in closed[1:] if np.all(disc[c])), None)
    touching = len(closed) == 1 and ndimage.label(~np.pad(o, 1), structure=BG4)[1] - 1 >= 2
    return disc, cup, touching


def split_touching(o):
    """F5: one stroke enclosing two regions -> (disc, cup), or None. The cup is the region with least contact with
    the disc exterior."""
    filled = ndimage.binary_fill_holes(o, structure=BG4)
    bg, n = ndimage.label(~np.pad(o, 1), structure=BG4)
    bg = bg[1:-1, 1:-1]
    outer = set(np.unique(np.concatenate([bg[0], bg[-1], bg[:, 0], bg[:, -1]]))) - {0}
    sizes = {k: int((bg == k).sum()) for k in range(1, n + 1) if k not in outer}
    regions = sorted(sizes, key=lambda k: -sizes[k])
    # Two real regions; pockets where the curves run side by side are ignored.
    if len(regions) < 2 or (len(regions) > 2 and sizes[regions[2]] >= 0.05 * sizes[regions[1]]):
        return None
    width = max(2, int(round(2 * ndimage.distance_transform_edt(np.pad(o, 1)).max())))
    best = None
    for k in regions[:2]:
        r = bg == k
        contact = (ndimage.binary_dilation(r, iterations=width + 2) & ~filled).sum() / max(
            1, (ndimage.binary_dilation(r) & ~r).sum())
        if best is None or contact < best[0]:
            best = (contact, r)
    cup = ndimage.binary_fill_holes(ndimage.binary_dilation(best[1], iterations=width) & filled, structure=BG4)
    return filled, cup


def certified(disc, cup) -> bool:
    """Disc and cup each one hole-free component, cup inside disc, all three (with the rim) certified."""
    reports = certify_many([disc, cup, disc & ~cup], [0, 0, RIM_MIN_SIZE])
    return all(r.certified for r in reports) and (reports[0].b0, reports[0].b1, reports[1].b0, reports[1].b1) == (
        1, 0, 1, 0) and bool(np.all(disc[cup]))


def crop_box(drawings, shape, min_comp=200, margin=30):
    """(y0, x0, y1, x1) around the components of at least `min_comp` px of the union of drawings, plus a margin."""
    if drawings:
        lab, keep = components(np.any(drawings, axis=0), min_comp)
        ys, xs = np.nonzero(keep[lab])
        if ys.size:
            return (max(0, ys.min() - margin), max(0, xs.min() - margin),
                    min(shape[0], ys.max() + margin + 1), min(shape[1], xs.max() + margin + 1))
    return 0, 0, shape[0], shape[1]


def reference(photo, copies=None):
    """(reference image, kind, F1 reconstruction): the photograph the raters drew on. kind is "prime", "maxrg6" (F1)
    or "maxrg6:prime_is_rater<e>_copy" (F7)."""
    copies = copies if copies is not None else {e: load(p) for e, p in photo.raters.items()}
    rebuilt = max_rg(list(copies.values()))
    if photo.prime is None:
        return rebuilt, "maxrg6", rebuilt
    prime = load(photo.prime)
    twin = [e for e, a in copies.items() if np.array_equal(a, prime)]
    if twin and int(ink(prime, rebuilt).sum()) > PRIME_INK_MAX:
        return rebuilt, f"maxrg6:prime_is_rater{twin[0]}_copy", rebuilt
    return prime, "prime", rebuilt


def extract(photo):
    """(reference, kind, {rater: (disc, cup, stage)}) with full-frame bool masks; unresolved or blank raters are
    absent. stage is "A", "B:<ref>/b<k>" (bridged k times) or "C:<ref>/b<k>" (touching curves split)."""
    copies = {e: load(p) for e, p in photo.raters.items()}
    ref, kind, rebuilt = reference(photo, copies)
    clean_prime = kind == "prime"
    if clean_prime:  # F6: a copy identical to a clean reference has no drawing
        copies = {e: a for e, a in copies.items() if not np.array_equal(a, ref)}
    jpeg = not clean_prime or photo.prime.suffix.lower() == ".jpg"
    fallbacks = [(ref, kind), (rebuilt, "maxrg6")] if clean_prime else [(ref, kind)]
    shape = ref.shape[:2]
    stage_a = {e: drawing(ref, a, jpeg) for e, a in copies.items()}
    y0, x0, y1, x1 = box_a = crop_box(list(stage_a.values()), shape)
    out = {}
    for e, a in copies.items():
        disc, cup, _ = interpret(stage_a[e][y0:y1, x0:x1])
        stage, box = "A", box_a
        if disc is None or cup is None:
            disc, cup, stage, box = stage_bc(a, fallbacks, jpeg, shape)
        if disc is None or not certified(disc, cup):
            continue
        full = np.zeros((2, *shape), bool)
        full[:, box[0]:box[2], box[1]:box[3]] = disc, cup
        out[e] = (full[0], full[1], stage)
    return ref, kind, out


def stage_bc(copy, fallbacks, jpeg, shape):
    """Stages B and C over the references and bridge levels, first success wins: (disc, cup, stage, box)."""
    for rr, rk in fallbacks:
        for bridge in BRIDGES:
            full = drawing(rr, copy, jpeg or rk.startswith("maxrg6"), bridge)
            box = crop_box([full], shape, min_comp=MIN_STROKE, margin=40)  # own box: stage A may have missed strokes
            ob = full[box[0]:box[2], box[1]:box[3]]
            disc, cup, touching = interpret(ob)
            if disc is not None and cup is not None:
                return disc, cup, f"B:{rk}/b{bridge}", box
            if touching and (got := split_touching(ob)) is not None and certified(*got):
                return *got, f"C:{rk}/b{bridge}", box
    return None, None, "", None


def descriptor(ref) -> np.ndarray:
    """Vessel-texture descriptor for near-duplicate search: high-pass central crop of a 576 x 384 grey thumbnail,
    zero mean, unit norm. Global hashes are dominated by the shared fundus layout; vessel texture identifies the eye
    (lab: same-eye pairs NCC >= 0.49, unrelated pairs <= 0.33)."""
    g = np.asarray(Image.fromarray(ref).convert("L").resize((576, 384), Image.BILINEAR))
    c = g[96:288, 144:432].astype(np.float32).reshape(96, 2, 144, 2).mean((1, 3))
    hp = c - ndimage.gaussian_filter(c, 2)
    hp -= hp.mean()
    return hp / (np.linalg.norm(hp) + 1e-6)


def keypoints(ref) -> tuple[np.ndarray, np.ndarray | None]:
    """ORB keypoints for the rotation- and scale-tolerant near-duplicate check (B3-8): (xy float32 (k, 2),
    descriptors uint8 (k, 32) or None) on the CLAHE-equalised green channel resized to ORB_WIDTH px (green carries
    the vessel contrast), detected only inside the camera's field of view shrunk by FIELD_INSET px. The field edge
    and its notch are the same in every photograph of one camera; unmasked, they gave unrelated pairs up to 27 RANSAC
    inliers at an identity transform, masked <= 6."""
    import cv2
    size = (ORB_WIDTH, int(round(ORB_WIDTH * ref.shape[0] / ref.shape[1])))
    g = cv2.resize(np.ascontiguousarray(ref[..., 1]), size, interpolation=cv2.INTER_AREA)
    field = cv2.resize(np.ascontiguousarray(ref[..., 0]), size, interpolation=cv2.INTER_AREA) > 20
    lab, k = ndimage.label(field, structure=FG8)
    if k:
        field = lab == np.bincount(lab.ravel())[1:].argmax() + 1
    field = cv2.erode(field.astype(np.uint8) * 255, np.ones((3, 3), np.uint8), iterations=FIELD_INSET)
    g = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(g)
    kp, des = cv2.ORB_create(nfeatures=ORB_FEATURES).detectAndCompute(g, field)
    return np.float32([k.pt for k in kp]).reshape(-1, 2), des


def item(photo: Photo, ref, kind, masks, split: str, group: str, duplicates: list[str]) -> Item:
    """The store item of one photograph: masks {"disc", "cup"} of shape (k, H, W) in rater order; no image for the
    BinRushed1 reconstruction (kind "maxrg6", B3-1)."""
    raters = sorted(masks)
    return Item(item_id=f"riga/{photo.image_id}", dataset="riga", group_id=f"riga:{group}", split=split,
                image=None if kind == "maxrg6" else ref,
                masks={name: np.stack([masks[e][j] for e in raters]) for j, name in enumerate(("disc", "cup"))},
                raters=[str(e) for e in raters],
                meta={"subset": photo.subset, "source": photo.subset.split("/")[0], "index": photo.index,
                      "reference": kind, "stages": {str(e): masks[e][2] for e in raters}, "duplicates": duplicates})
