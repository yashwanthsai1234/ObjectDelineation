"""SFT targets: rater masks of one item -> certified, canonical polygon text (plan D11, gates K5-K7, R1).

Per object (structure) of one rater:
1. exact crack-edge rings (objdel.checker.polygonize), each ring's holes grouped with its 8-connected component;
2. every ring simplified by `simplify` at tolerance tau (min-# polygon, lossless or L-inf tau shortcuts);
3. corners mapped to the 0-1000 grid; every vertex the object uses twice (a pinch: two pixels touching at a corner)
   split into two points one grid unit along the pinch's incoming and outgoing crack steps ("bevel", legal and
   pixel-neutral for crops <= 500 px); repeated and straight-through vertices dropped;
4. certified by the CUDA checker through objdel.geometry.to_mask: legal, at most MAX_VERTICES per ring, the same
   per-instance signature as the mask under the size rule (plan §4.2: the rule applies to targets, reward and eval),
   and IoU >= iou_min. The rule is the reward's min_area (of the raters' union) when the caller passes it; without
   it the union is unknown, so the signature must agree at every rule a union on this crop can have, from the
   mask's own min_area to min_area(H * W) (break-it B1-3: agreement at one rule does not imply it at a larger one,
   and the signature is not even monotone in the rule). `_kept` checks only the rules where a signature can change:
   one above the size of a raw component (the opening drops it) or of a hole of the filtered mask (the closing
   fills it). For RIGA (disc, cup) the rim disc & ~cup must keep its signature too. The raw signature is
   not required: on jagged masks (RIGA BinRushed) the tau = 3 polygon rasterizes with sub-filter pinholes or specks,
   and demanding the raw signature pushed 1-2 % of RIGA objects to near-lossless answers of up to 6,000 tokens.
   The ring cap MAX_RINGS applies after the size rule (break-it B1-4): an object with more raw rings is traced from
   its size-filtered mask instead, and still certified against the raw mask.
   LIDC nodules (BuildConfig.exact) must also rasterize to the mask pixel for pixel (D9). tau is in grid units and
   a pixel is 1000 / S of them, so at S >= 167 tau = 3 is >= 1/2 px and a shortcut can cross a pixel centre
   (break-it BK1-1: 20 % of S = 200 view rows were certified but not exact). The check steps those down the
   ladder; every identity (policy-frame) LIDC target is already exact, so the Phase 1 target-set hashes do not
   change. If no rung is exact (rings traced from the size-filtered mask, or the vertex cap), the object falls
   back to the certificate alone. RIGA is exempt: its jagged BinRushed masks would be pushed to near-lossless
   answers, as the raw signature did (above).
5. If an object fails, it is rebuilt one rung lower on the tolerance ladder tau = 3, 2.75, ..., 0.25, then lossless
   only, then the exact rings, which pass by construction for crops <= 500 px (the lab's stress test, S = 125-501).
The answer is the canonical text of all objects (objdel.geometry.canonicalize / serialize); it must parse.

Repair: the lab localized defects (pinning points of offending segments, halving the tolerance per ring, then the
better of both). On the 3,169 fixture targets the per-object ladder with 0.25 steps needs 0.81x (disc) and 0.77x
(cup) the lab's vertices on RIGA, and the same LIDC texts (99.85 %; the rest differ only in the bevel), so the
localization was dropped (experiments/phase1/impl/M5_targets.md).

References:
    P. Selinger, "Potrace: a polygon-based tracing algorithm", 2003: https://potrace.sourceforge.net/potrace.pdf
    D. Estkowski, J. S. B. Mitchell, "Simplifying a polygonal subdivision while keeping it simple" (min-vertex
        simplification without crossings is hard; simplify per ring, then verify and repair):
        https://www.uni-trier.de/fileadmin/fb4/prof/INF/DEA/Seminar0708/Estkowski.pdf
    experiments/phase1/lab/t_targets/REPORT.md (minseg_hybrid, tau = 3, IoU >= 0.99, bevel, crop policy)
"""
import hashlib
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from typing import NamedTuple

import cv2
import numpy as np

from objdel.checker import as_mask, instance_topology, polygonize, signature, size_filter, topology
from objdel.geometry import GRID, MAX_RINGS, MAX_VERTICES, Objects, canonicalize, parse, serialize, to_mask

from .simplify import STEPS, lattice_path, simplify, to_grid

MAX_SIDE = 500  # crop policy; beyond it the exact-ring fallback is no longer guaranteed (lab stress test, S = 640)
_RIGHT = np.array([[0, 0], [-1, 0], [-1, -1], [0, -1]])  # pixel right of a step E S W N, relative to its start


@dataclass(frozen=True)
class BuildConfig:
    tau: float = 3.0  # largest L-inf tolerance, grid units
    step: float = 0.25  # ladder step down to lossless
    iou_min: float = 0.99
    exact: tuple = ("nodule",)  # objects whose target must rasterize to the mask pixel for pixel (step 4)

    def ladder(self) -> list:
        k = round(self.tau / self.step)
        return [self.step * i for i in range(k, -1, -1)] + [None]  # ..., 0 = lossless, None = exact rings


class Target(NamedTuple):
    text: str
    objs: Objects  # canonical, as parsed back from text
    info: dict  # per object: tau (None = exact rings), vertices, iou; "time" (CPU s)


def _rings(mask):
    """(lattice path, step directions, component label, is hole) of every ring, traced on the mask cut to its
    bounding box (crop before trace, Phase 2 B3 P4: every target text and hash unchanged, RIGA -23 % build time).
    The texts cannot change: boundary_edges pads by one pixel itself, the ring start key orders vertices by (y, x, d)
    for any width, so a translation keeps every ring's start and order, and connectedComponents labels in raster
    order, which a crop keeps."""
    mask = np.ascontiguousarray(mask, dtype=bool)
    x0, y0, w, h = cv2.boundingRect(mask.view(np.uint8))
    if not w:
        return []
    sub = np.ascontiguousarray(mask[y0:y0 + h, x0:x0 + w])
    _, label = cv2.connectedComponents(sub.view(np.uint8), connectivity=8)
    out = []
    for r in polygonize(sub, offset=(x0, y0)):
        path, dirs = lattice_path(r.pts)
        x, y = path[0] + _RIGHT[dirs[0]]  # the foreground pixel right of the first step
        out.append((path, dirs, int(label[y - y0, x - x0]), r.hole))
    return out


def _clean(v):
    """Drop repeated and straight-through vertices of a closed ring."""
    while len(v) >= 3:
        a, b = v - np.roll(v, 1, 0), np.roll(v, -1, 0) - v
        drop = np.all(a == 0, 1) | ((a[:, 0] * b[:, 1] == a[:, 1] * b[:, 0]) & ((a * b).sum(1) > 0))
        if not drop.any():
            break
        v = v[~drop]
    return v


def _bevel(rings, steps_in, steps_out):
    """Split each vertex the object uses more than once (a pinch) into two points one grid unit back along the crack
    step arriving there and one unit along the step leaving it. The cut stays inside the pinch's background pixel
    corner, so no pixel centre changes side for crops <= 500 px. (The lab stepped along the simplified edges; on a
    slanted edge that rounds to a pixel centre at S = 500 and flips it.)"""
    pts, counts = np.unique(np.concatenate(rings), axis=0, return_counts=True)
    twice = {tuple(p) for p in pts[counts > 1].tolist()}
    if not twice:
        return rings
    out = []
    for v, din, dout in zip(rings, steps_in, steps_out, strict=True):
        pts = []
        for p, a, b in zip(v, din, dout, strict=True):
            pts += [p - a, p + b] if tuple(p.tolist()) in twice else [p]
        out.append(np.array(pts, np.int64))
    return out


def _parts(rings, tau, W, H):
    vs, steps_in, steps_out = [], [], []
    for path, dirs, _, _ in rings:
        keep = np.flatnonzero(dirs != np.roll(dirs, 1)) if tau is None else simplify(path, dirs, W, H, tau)
        vs.append(to_grid(path[keep], W, H, GRID))
        steps_in.append(STEPS[dirs[keep - 1]])  # dirs[k] is the step from path[k] to path[k + 1]
        steps_out.append(STEPS[dirs[keep]])
    parts = {}
    for v, (_, _, comp, hole) in zip(_bevel(vs, steps_in, steps_out) if vs else vs, rings, strict=True):
        part = parts.setdefault(comp, [None, []])
        if hole:
            part[1].append(_clean(v))
        else:
            part[0] = _clean(v)
    return [(o, hs) for o, hs in parts.values()]


def _kept(a, b, lo: int, hi: int) -> bool:
    """signature(a, t) == signature(b, t) for every rule t in [lo, hi]. Between consecutive rules where either
    filtered mask can change (a raw component of t' - 1 px is removed at t'; a hole of the filtered mask of t' - 1
    px is filled at t'), the opening and the set of filled holes are constant, so one check per interval suffices."""
    if lo == hi:  # the reward's rule is known
        return signature(a, lo) == signature(b, lo)
    raw = np.concatenate([instance_topology(a).areas, instance_topology(b).areas])
    t = lo
    while t <= hi:
        ta, tb = instance_topology(size_filter(a, t)), instance_topology(size_filter(b, t))
        if sorted(ta.holes.tolist()) != sorted(tb.holes.tolist()):
            return False
        sizes = np.concatenate([raw, ta.hole_areas, tb.hole_areas])
        sizes = sizes[sizes >= t]
        if not sizes.size:
            return True
        t = int(sizes.min()) + 1
    return True


def min_areas(stacks: Mapping[str, np.ndarray]) -> dict[str, int]:
    """The reward's size rule per structure of one item's rater stacks {name: (k, H, W)}: min_area of the union of
    all raters, as ExpertCache computes it, and for RIGA (disc, cup) also of the rim union (disc & ~cup per rater).
    This is build_target's min_area argument for every rater of the item."""
    if {"disc", "cup"} <= set(stacks):
        stacks = {**stacks, "rim": stacks["disc"] & ~stacks["cup"]}
    return {n: topology.min_area(int(s.any(0).sum())) for n, s in stacks.items()}


def build_target(masks: Mapping[str, np.ndarray], names, cfg: BuildConfig | None = None,
                 min_area: Mapping[str, int] | None = None) -> Target:
    """One rater's masks {name: (H, W), any nonzero value is foreground} on one crop -> certified target.
    min_area: the reward's size filter per object name (and "rim" for RIGA), min_area of the union of all raters, as
    ExpertCache computes it (`min_areas`); the target keeps the size-filtered signature at that rule. A name without
    an entry is certified at every rule a union on this crop can have (module docstring), which is safe for any union
    but can need more vertices. Objects in cfg.exact must also be pixel-exact; when no rung is (rings traced from
    the size-filtered mask, or the vertex cap), the target is built again without that demand (module docstring).
    Raises ValueError when no target can be certified."""
    cfg = cfg or BuildConfig()
    try:
        return _build(masks, names, cfg, min_area or {})
    except ValueError:
        if not set(cfg.exact) & set(names):
            raise
        return _build(masks, names, replace(cfg, exact=()), min_area or {})


def _build(masks, names, cfg: BuildConfig, min_area: Mapping[str, int]) -> Target:
    t0 = time.process_time()
    src = {n: as_mask(masks[n]) for n in names}
    shapes = {m.shape for m in src.values()}
    if len(shapes) != 1 or not 0 < max(shape := shapes.pop()) <= MAX_SIDE:
        raise ValueError(f"masks must share one (H, W) shape with sides <= {MAX_SIDE}, got {shapes or shape}")
    H, W = shape
    rim = {"disc", "cup"} <= set(names)
    ref = {**src, "rim": src["disc"] & ~src["cup"]} if rim else src
    rules = {n: (min_area[n], min_area[n]) if n in min_area else
             (topology.min_area(np.count_nonzero(m)), topology.min_area(H * W)) for n, m in ref.items()}
    rings = {}
    for n in names:
        rings[n] = _rings(src[n])
        if len(rings[n]) > MAX_RINGS:  # the cap applies after the size rule (module docstring)
            rings[n] = _rings(size_filter(src[n], rules[n][0]))
        if len(rings[n]) > MAX_RINGS:  # the text would not parse (cap_rings)
            raise ValueError(f"{n}: {len(rings[n])} rings > MAX_RINGS = {MAX_RINGS} after the size rule")

    def kept(n, mask):  # same signature under the size rule
        return _kept(mask, ref[n], *rules[n])

    ladder, level = cfg.ladder(), dict.fromkeys(names, 0)
    while True:
        objs, info, ras = {}, {}, {}
        for n in names:
            while True:
                parts = _parts(rings[n], ladder[level[n]], W, H)
                if all(len(r) <= MAX_VERTICES for o, hs in parts for r in (o, *hs)):  # else the text would not parse
                    g = to_mask(parts, W, H)
                    u = np.count_nonzero(g.mask | src[n])
                    iou = 1.0 if u == 0 else np.count_nonzero(g.mask & src[n]) / u
                    if g.legal and iou >= cfg.iou_min and kept(n, g.mask) and (
                            n not in cfg.exact or np.array_equal(g.mask, src[n])):
                        break
                if level[n] == len(ladder) - 1:
                    raise ValueError(f"{n}: no certified target")
                level[n] += 1
            objs[n], ras[n] = parts, g.mask
            info[n] = {"tau": ladder[level[n]], "vertices": sum(len(r) for o, hs in parts for r in (o, *hs)),
                       "iou": iou}
        if not rim or kept("rim", ras["disc"] & ~ras["cup"]):
            break
        if all(level[n] == len(ladder) - 1 for n in names):
            raise ValueError("rim: no certified target")
        level = {n: min(level[n] + 1, len(ladder) - 1) for n in names}
    text = serialize(canonicalize(objs), names)
    p = parse(text, names)
    if not p.ok:
        raise ValueError(f"target text does not parse: {p.error}")
    info["time"] = time.process_time() - t0
    return Target(text, p.objs, info)


def target_set_hash(rows: Iterable[tuple[str, str, str]]) -> str:
    """SHA-256 of a target set: rows (item_id, rater, text), order-independent (sorted lines)."""
    h = hashlib.sha256()
    for item_id, rater, text in sorted(rows):
        h.update(f"{item_id}\t{rater}\t{text}\n".encode())
    return h.hexdigest()
