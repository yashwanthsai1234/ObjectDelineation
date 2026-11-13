"""SFT targets: rater masks of one item -> certified, canonical polygon text .

Per object (structure) of one rater:
1. exact crack-edge rings (objdel.geom.rings.trace), each ring's holes grouped with its 8-connected component;
2. every ring simplified by `simplify` at tolerance tau (min-# polygon, lossless or L-inf tau shortcuts);
3. corners mapped to the 0-1000 grid; every vertex the object uses twice (a pinch: two pixels touching at a corner)
   split into two points one grid unit along the pinch's incoming and outgoing crack steps ("bevel", legal and
   pixel-neutral for crops <= 500 px); repeated and straight-through vertices dropped;
4. certified on the CUDA engine (objdel.geom.engine.certify_polygons: raster, legality, IoU and the size-filtered
   signatures of the rung and the mask, only numbers leave the device): legal, at most MAX_VERTICES per ring, the
   same per-instance signature as the mask under the size rule (the rule applies to targets, reward and eval), and
   IoU >= iou_min. The rule is the reward's min_area of the raters' union, which every caller passes (`min_areas`).
   The raw signature is not required: on jagged masks the tau = 3 polygon rasterizes with sub-filter
   pinholes or specks, and demanding the raw signature pushed objects to near-lossless answers.
   The ring cap MAX_RINGS applies after the size rule: an object with more raw rings is traced from
   its size-filtered mask instead, and still certified against the raw mask.
   LIDC nodules (BuildConfig.exact) must also rasterize to the mask pixel for pixel. tau is in grid units and
   a pixel is 1000 / S of them, so at S >= 167 tau = 3 is >= 1/2 px and a shortcut can cross a pixel centre
   (20 % of S = 200 view rows were certified but not exact). The check steps those down the
   ladder; every identity (policy-frame) LIDC target is already exact, so the Phase 1 target-set hashes do not
   change. If no rung is exact (rings traced from the size-filtered mask, or the vertex cap), the object falls
   back to the certificate alone.
5. If an object fails, it is rebuilt one rung lower on the tolerance ladder tau = 3, 2.75, ..., 0.25, then lossless
   only, then the exact rings, which pass by construction for crops <= 500 px (stress test, S = 125-501).
   Rungs are certified in engine calls of 1, 2, 4, ... rungs from the current one (one call per start rung); the
   first certified rung wins, so the text is the one rung-by-rung certification gives. Simplifying every remaining
   rung up front would cost 3x the whole build on LIDC, where 99 % pass the first rung.
The answer is the canonical text of all objects (objdel.geom.shape.canonicalize / serialize); it must parse.

References:
    P. Selinger, "Potrace: a polygon-based tracing algorithm", 2003: https://potrace.sourceforge.net/potrace.pdf
    D. Estkowski, J. S. B. Mitchell, "Simplifying a polygonal subdivision while keeping it simple" (min-vertex
        simplification without crossings is hard; simplify per ring, then verify and repair):
        https://www.uni-trier.de/fileadmin/fb4/prof/INF/DEA/Seminar0708/Estkowski.pdf
"""
import hashlib
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from typing import NamedTuple

import cv2
import numpy as np

from objdel.data.targets.simplify import STEPS, lattice_path, simplify, to_grid
from objdel.geom import engine as _gpucheck
from objdel.geom import topology
from objdel.geom.rings import trace
from objdel.geom.shape import canonicalize, pack_objects
from objdel.geom.text import GRID, MAX_RINGS, MAX_VERTICES, Objects, parse, serialize
from objdel.geom.topology import as_mask, pack_masks, size_filter

MAX_SIDE = 500  # crop policy; beyond it the exact-ring fallback is no longer guaranteed (stress test, S = 640)


@dataclass(frozen=True)
class BuildConfig:
    """Target build settings: largest tolerance and ladder step (grid units, 0-1000), minimum IoU and the objects that
    must be exact.
    """
    tau: float = 3.0  # largest L-inf tolerance, grid units
    step: float = 0.25  # ladder step down to lossless
    iou_min: float = 0.99
    exact: tuple = ("nodule",)  # objects whose target must rasterize to the mask pixel for pixel (step 4)

    def ladder(self) -> list:
        """Tolerances to try, from `tau` down to 0 in steps of `step`, then None (exact rings)."""
        k = round(self.tau / self.step)
        return [self.step * i for i in range(k, -1, -1)] + [None]  # ..., 0 = lossless, None = exact rings


class Target(NamedTuple):
    """A built target: the text, its objects parsed back, and per-object build info."""
    text: str
    objs: Objects  # canonical, as parsed back from text
    info: dict  # per object: tau (None = exact rings), vertices, iou; "time" (CPU s)


def _rings(mask):
    """(lattice path, step directions, component, is hole) of every ring, traced on the mask cut to its bounding box
    (crop before trace: every target text and hash unchanged). The component is
    the ring's part (the 8-connected component it bounds; it replaced cv2.connectedComponents labels, which group
    the rings the same way). The texts cannot change: the ring start key orders vertices by (y, x, d) for any width,
    so a translation keeps every ring's start and order."""
    mask = np.ascontiguousarray(mask, dtype=bool)
    x0, y0, w, h = cv2.boundingRect(mask.view(np.uint8))
    if not w:
        return []
    return [(*lattice_path(r.pts), r.part, r.hole) for r in trace([mask[y0:y0 + h, x0:x0 + w]], [(x0, y0)])[0].rings]


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
    corner, so no pixel centre changes side for crops <= 500 px. (Stepping along the simplified edges instead
    would, on a slanted edge at S = 500, round to a pixel centre and flip it.)"""
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


def _certify(objs, refs, ref_of, min_area) -> dict:
    """certify_polygons of objects [(parts, W, H)] against as_mask reference masks: one engine call."""
    return _gpucheck.certify_polygons(
        *pack_objects(objs), GRID, *pack_masks(refs), np.array(ref_of, np.int32), np.array(min_area, np.int64))


def _sig(out: dict, key: str, i: int) -> list:
    s = out[key + "_start"]
    return out[key][s[i]:s[i + 1]].tolist()


def _climb(name, rings, ladder, level, ref, min_area, exact, iou_min, W, H):
    """(rung, parts, IoU) of the first certified rung from `level` on (module docstring, steps 4-5)."""
    size = 1
    while level < len(ladder):
        rungs = []
        for k in range(level, min(level + size, len(ladder))):
            parts = _parts(rings, ladder[k], W, H)
            if all(len(r) <= MAX_VERTICES for o, hs in parts for r in (o, *hs)):  # else the text would not parse
                rungs.append((k, parts))
        if rungs:
            c = _certify([(p, W, H) for _, p in rungs], [ref], [0] * len(rungs), [min_area] * len(rungs))
            for i, (k, parts) in enumerate(rungs):
                iou = topology.iou(int(c["inter"][i]), int(c["union"][i]))
                if c["legal"][i] and iou >= iou_min and _sig(c, "sig", i) == _sig(c, "ref_sig", i) and (
                        not exact or c["equal"][i]):
                    return k, parts, iou
        level, size = level + size, 2 * size
    raise ValueError(f"{name}: no certified target")


def min_areas(stacks: Mapping[str, np.ndarray]) -> dict[str, int]:
    """The reward's size rule per structure of one item's rater stacks {name: (k, H, W)}: min_area of the union of
    all raters, as ExpertCache computes it. This is build_target's min_area argument for every rater of the item."""
    return {n: topology.union_min_area(s) for n, s in stacks.items()}


def rater_targets(stacks: Mapping[str, np.ndarray], names, k: int, cfg: BuildConfig | None = None
                  ) -> tuple[dict[str, int], list["Target | ValueError"]]:
    """(min_areas(stacks), one entry per rater r < k): build_target of rater r's masks of one item's frame stacks
    {name: (k, H, W)} at that size rule, or the ValueError it raised (each caller handles a failed build its way)."""
    ma, out = min_areas(stacks), []
    for r in range(k):
        try:
            out.append(build_target({n: stacks[n][r] for n in names}, names, cfg, ma))
        except ValueError as e:
            out.append(e.with_traceback(None))  # the frames of the failed build are not kept alive
    return ma, out


def build_target(masks: Mapping[str, np.ndarray], names, cfg: BuildConfig | None = None,
                 min_area: Mapping[str, int] | None = None) -> Target:
    """One rater's masks {name: (H, W), any nonzero value is foreground} on one crop -> certified target.
    min_area (required): the reward's size filter per object name, min_area of the union of all raters, as
    ExpertCache computes it (`min_areas`); the target keeps the size-filtered signature at that rule. A missing
    entry raises ValueError. Objects in cfg.exact must also be pixel-exact; when no rung is (rings traced from the
    size-filtered mask, or the vertex cap), the target is built again without that demand (module docstring).
    Raises ValueError when no target can be certified."""
    cfg, min_area = cfg or BuildConfig(), min_area or {}
    for n in names:
        if n not in min_area:
            raise ValueError(f"min_area has no entry for {n!r}")
    try:
        return _build(masks, names, cfg, min_area)
    except ValueError:
        if not set(cfg.exact) & set(names):
            raise
        return _build(masks, names, replace(cfg, exact=()), min_area)


def _build(masks, names, cfg: BuildConfig, min_area: Mapping[str, int]) -> Target:
    t0 = time.process_time()
    src = {n: as_mask(masks[n]) for n in names}
    shapes = {m.shape for m in src.values()}
    if len(shapes) != 1 or not 0 < max(shape := shapes.pop()) <= MAX_SIDE:
        raise ValueError(f"masks must share one (H, W) shape with sides <= {MAX_SIDE}, got {shapes or shape}")
    H, W = shape
    rings = {}
    for n in names:
        rings[n] = _rings(src[n])
        if len(rings[n]) > MAX_RINGS:  # the cap applies after the size rule (module docstring)
            rings[n] = _rings(size_filter(src[n], min_area[n]))
        if len(rings[n]) > MAX_RINGS:  # the text would not parse (cap_rings)
            raise ValueError(f"{n}: {len(rings[n])} rings > MAX_RINGS = {MAX_RINGS} after the size rule")
    ladder = cfg.ladder()
    objs, info = {}, {}
    for n in names:
        level, objs[n], iou = _climb(n, rings[n], ladder, 0, src[n], min_area[n], n in cfg.exact, cfg.iou_min, W, H)
        info[n] = {"tau": ladder[level], "vertices": sum(len(r) for o, hs in objs[n] for r in (o, *hs)), "iou": iou}
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
