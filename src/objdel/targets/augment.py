"""Phase 2 training views (D9, D10; plan §6.0): exact dihedral ops, view draws, context crops and the renderer.

OPS: the dihedral group D4 on the last two axes (y, x) of an array: "id"; "r90", "r180", "r270" (np.rot90 with
k = 1, 2, 3, counter-clockwise as displayed); "fh" (x mirrored), "fv" (y mirrored), "tr" (transpose), "at"
(anti-transpose). LIDC uses all eight, RIGA ("id", "fh"): fundus orientation is physiological. The ops map the
pixel lattice onto itself, so masks are transformed exactly and targets are rebuilt from them (build_target: the
top-left tie rule is not mirror-symmetric, issue I7); vertices are never transformed.

draw(item, rng, side=None): one D10 window: side S uniform over `sides(item)`, the policy side and its LADDER
neighbours (scale 1), i.e. every side the shifted and rescaled eval views ask for (plan §5.0), then an offset uniform
over the windows that keep the union inside with CLEARANCE px and meet no other outlined nodule (rejection sampling,
so uniform over the admissible sides and offsets). neighbour_side: the eval "rescaled" view's S (plan §5.0); edge:
the eval "edge" view's offset, uniform on the boundary of the admissible rectangle (B1.6).
The side set is the eval views' on purpose: a margin draw ([8 px, 60 % of the union's extent], the R1 bundle) pinned
every union under 23 px at S = 50 (min_side), so 27 % of the rescaled val's nodules fill < 10 % of the frame against
0.5 % of those training rows.

Context crops (C-VIEW): per item the image region every admissible window of every side the scale draw can give
(and the policy window) can show, clipped to the image: int16 Hounsfield units for LIDC (window
jitter and noise act on raw HU, D10), uint8 RGB for RIGA. render(spec, context) = frame(item, offset, side, window)
on the context with the spec's HU noise added, the op applied at frame resolution, then the Pillow upscale; with
op "id", sigma 0 and the policy window and side it is the prepared PNG byte for byte (tests/test_p2_views.py).

Row seeds: row_rng(*keys) seeds numpy from SHA-256 of the keys, so every row's draw depends only on (seed, epoch,
item_id, rater_index, rep), not on the order or subset of items built.
"""
import hashlib
import zipfile
from pathlib import Path

import numpy as np

from objdel.data import Item

from .frame import CLEARANCE, DIVISORS, FRAMES, LUNG_WINDOW, centred, crowded, pixels, policy_side, rgb, union_box

OPS = ("id", "r90", "r180", "r270", "fh", "fv", "tr", "at")
DATASET_OPS = {"lidc": OPS, "riga": ("id", "fh")}
INVERSE = {"id": "id", "r90": "r270", "r180": "r180", "r270": "r90", "fh": "fh", "fv": "fv", "tr": "tr", "at": "at"}
LADDER = tuple(d for d in DIVISORS if 50 <= d <= 500)  # 50, 100, 125, 200, 250, 500
WINDOW_JITTER = ((-300, 60), (1400, 80))  # CT level and width, centre and half-range in HU (D10)
NOISE_MAX = 10.0  # HU, Gaussian sigma upper bound (D10)


_OP = {"id": lambda a: a, "r90": lambda a: np.rot90(a, 1, (-2, -1)), "r180": lambda a: np.rot90(a, 2, (-2, -1)),
       "r270": lambda a: np.rot90(a, 3, (-2, -1)), "fh": lambda a: a[..., :, ::-1], "fv": lambda a: a[..., ::-1, :],
       "tr": lambda a: np.swapaxes(a, -1, -2), "at": lambda a: np.rot90(np.swapaxes(a, -1, -2), 2, (-2, -1))}


def apply(a: np.ndarray, op: str) -> np.ndarray:
    """The D4 op on the last two axes (contiguous copy)."""
    return np.ascontiguousarray(_OP[op](a))


def row_rng(*keys) -> np.random.Generator:
    return np.random.default_rng(int.from_bytes(hashlib.sha256("|".join(map(str, keys)).encode()).digest()[:8]))


def _admissible(item: Item, side: int, box):
    """(centred corner, lo, hi): the offsets in [lo, hi] (inclusive) keep the union inside with CLEARANCE px."""
    x0, y0, x1, y1 = box
    cx, cy = centred(item, side, box)
    lo, hi = (x1 + CLEARANCE - side - cx, y1 + CLEARANCE - side - cy), (x0 - CLEARANCE - cx, y0 - CLEARANCE - cy)
    return (cx, cy), lo, hi


def shift(item: Item, side: int, rng, tries: int = 200, box=None) -> tuple[int, int] | None:
    """A uniform offset of a window of `side` from its centred corner that keeps the union inside with CLEARANCE px
    and meets no other outlined nodule (rejection sampling), or None."""
    (cx, cy), lo, hi = _admissible(item, side, box := box or union_box(item))
    if lo[0] > hi[0] or lo[1] > hi[1]:
        return None
    for _ in range(tries):
        dx, dy = (int(v) for v in rng.integers(lo, np.add(hi, 1)))
        if not crowded(item, (cx + dx, cy + dy, cx + dx + side, cy + dy + side)):
            return dx, dy
    return None


def edge(item: Item, side: int, rng, tries: int = 64) -> tuple[int, int] | None:
    """The guard view's offset (B1.6): a lattice point drawn uniformly on the boundary of shift's admissible
    rectangle (one coordinate at its lo or hi: the largest legal shift that way), the first of `tries` draws that
    meets no other outlined nodule, or None."""
    (cx, cy), lo, hi = _admissible(item, side, union_box(item))
    w, h = hi[0] - lo[0], hi[1] - lo[1]
    if w < 0 or h < 0:
        return None
    for _ in range(tries):  # t walks the perimeter: bottom, right, top, left (a segment or a point: each once)
        t = int(rng.integers(2 * (w + h) if w and h else w + h + 1))
        dx, dy = ((lo[0] + t, lo[1]) if t < w else (hi[0], lo[1] + t - w) if t < w + h
                  else (hi[0] - (t - w - h), hi[1]) if t < 2 * w + h else (lo[0], hi[1] - (t - 2 * w - h)))
        if not crowded(item, (cx + dx, cy + dy, cx + dx + side, cy + dy + side)):
            return dx, dy
    return None


def _neighbours(side: int) -> list[int]:
    i = LADDER.index(side)
    return [LADDER[j] for j in (i - 1, i + 1) if 0 <= j < len(LADDER)]


def sides(item: Item, box=None) -> list[int]:
    """The sides of the scale draw (D10): the policy side and its LADDER neighbours up to max_side, ascending."""
    cfg = FRAMES[item.dataset]
    S0 = policy_side(box or union_box(item), cfg)[1]
    return sorted(s for s in (S0, *_neighbours(S0)) if s <= cfg.max_side)


def draw(item: Item, rng, side: int | None = None, tries: int = 20) -> tuple[tuple[int, int], int] | None:
    """((dx, dy), S) of one D10 window (module docstring; a fixed `side`, as in the eval views, draws the shift
    only), or None when no admissible window was found."""
    box = union_box(item)
    choices = [side] if side else sides(item, box)
    for _ in range(tries):
        S = side or choices[int(rng.integers(len(choices)))]
        if (off := shift(item, S, rng, box=box)) is not None:
            return off, S
    return None


def neighbour_side(item: Item, side: int, rng) -> tuple[tuple[int, int], int] | None:
    """The rescaled view (plan §5.0): S moved to a neighbouring LADDER divisor where an admissible shifted window
    exists (one drawn uniformly when both do), with a uniform shift; None when neither neighbour is admissible."""
    for s in rng.permutation(_neighbours(side)).tolist():
        if s <= FRAMES[item.dataset].max_side and (off := shift(item, s, rng)) is not None:
            return off, s
    return None


def jitter(rng) -> dict:
    """CT window and noise of one LIDC row (D10)."""
    (lc, lr), (wc, wr) = WINDOW_JITTER
    level, width = lc + int(rng.integers(-lr, lr + 1)), wc + int(rng.integers(-wr, wr + 1))
    lo = level - width // 2
    return {"window": [lo, lo + width], "noise_sigma": round(float(rng.uniform(0, NOISE_MAX)), 3),
            "noise_seed": int(rng.integers(2 ** 31))}


def context_box(item: Item, policy_box) -> tuple[int, int, int, int]:
    """The region (clipped to the image) that any admissible window of a `sides` side, or the policy window, shows."""
    x0, y0, x1, y1 = box = union_box(item)
    S = max(sides(item, box))  # the largest side the scale draw can give
    H, W = next(iter(item.masks.values())).shape[-2:]
    bx0, by0 = min(x1 + CLEARANCE - S, policy_box[0]), min(y1 + CLEARANCE - S, policy_box[1])
    bx1, by1 = max(x0 - CLEARANCE + S, policy_box[2]), max(y0 - CLEARANCE + S, policy_box[3])
    return max(bx0, 0), max(by0, 0), min(bx1, W), min(by1, H)


def write_context(path: Path, item: Item, policy_box) -> None:
    """context/<item_id>.npz: pixels (int16 HU or uint8 RGB), x0, y0 and the image shape. Written with a fixed zip
    timestamp so two builds are byte-identical (np.savez stamps the current time)."""
    x0, y0, x1, y1 = context_box(item, policy_box)
    arrays = {"pixels": np.ascontiguousarray(item.image[y0:y1, x0:x1]), "x0": np.int64(x0), "y0": np.int64(y0),
              "shape": np.array(item.image.shape[:2], np.int64)}
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, a in arrays.items():
            with zf.open(zipfile.ZipInfo(f"{name}.npy", (1980, 1, 1, 0, 0, 0)), "w") as fh:
                np.lib.format.write_array(fh, np.asarray(a), allow_pickle=False)


def load_context(view_dir, item_id: str) -> dict:
    with np.load(Path(view_dir, "context", item_id + ".npz"), allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def render(spec: dict, context: dict) -> np.ndarray:
    """uint8 (side_img, side_img, 3) model image of one row spec (module docstring). The context holds every image
    pixel of the window, so outside it is outside the image: frame.pixels pads and windows it exactly as frame."""
    S, (x0, y0) = int(spec["side"]), spec["box"][:2]
    px, cx, cy = context["pixels"], int(context["x0"]), int(context["y0"])
    H, W = (int(v) for v in context["shape"])
    ix0, iy0, ix1, iy1 = max(x0, 0), max(y0, 0), min(x0 + S, W), min(y0 + S, H)  # window inside the image
    if ix0 < cx or iy0 < cy or ix1 > cx + px.shape[1] or iy1 > cy + px.shape[0]:
        raise ValueError(f"{spec.get('row_id')}: window {spec['box']} leaves the context crop")
    sigma = spec.get("noise_sigma")
    noise = np.rint(row_rng("noise", spec["noise_seed"]).normal(0, sigma, (S, S))).astype(np.int64) if sigma else None
    chw = pixels(px, x0 - cx, y0 - cy, S, tuple(spec.get("window") or LUNG_WINDOW), noise)
    return rgb(apply(chw, spec.get("op", "id")), FRAMES[spec["dataset"]].image_side)
