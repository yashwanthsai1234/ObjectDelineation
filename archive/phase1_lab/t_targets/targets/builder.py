"""SFT target builder: expert mask(s) -> legal, canonical, compact polygon JSON.

Pipeline per object (one structure mask of one rater, crop of W x H pixels):
  exact crack-edge rings (objdel.checker.polygonize)
  -> per-ring simplification on the pixel-corner lattice (targets.simplify)
  -> map to the integer 0..G grid (round half up, x by W, y by H)
  -> pinch bevel (vertices shared within an object are split by one grid unit)
  -> certify: legality rule L, exact re-rasterization (oracle rule), (b0, b1) and IoU vs the source
  -> local repair: split the offending simplified segments at their farthest skipped lattice point;
     iterate; fallback = exact rings (+ bevel). Terminates: every round adds >= 1 kept point or falls back.

Legality rule L (per object): every ring passes oracle.ring_problems == [] (>= 3 distinct points,
non-zero area, no self-crossing, no self-touching); no two rings of the object touch or cross;
outer rings clockwise on screen (positive y-down shoelace), holes counter-clockwise; every hole lies
inside its part's outer ring; coordinates are integers in [0, G].
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from objdel.checker.polygonize import polygonize
from objdel.checker.bitquad import bitquads
from objdel.checker.certify import crop as _crop

from . import fastgeom as F
from .simplify import METHODS, corners_only, lattice_path, to_grid

G_DEFAULT = 1000
FG8 = np.ones((3, 3), bool)
# FG pixel (col, row) on the right of a unit step starting at lattice point p with direction d
_RIGHT_PIX = {0: (0, 0), 1: (-1, 0), 2: (-1, -1), 3: (0, -1)}


@dataclass
class RingState:
    path: np.ndarray
    dirs: np.ndarray
    K: np.ndarray  # sorted kept lattice indices
    hole: bool
    comp: int
    level: int = 0  # repair level: 0 = configured tolerance, halved per level; then floor, then corners
    pins: np.ndarray = None  # lattice indices the simplified ring must keep ("pin" repair)


@dataclass
class ObjectTarget:
    name: str
    parts: list  # [(outer (k,2) grid, [holes...])] canonical
    status: str  # "ok" | "repaired" | "fallback" | "fail"
    iou: float
    betti_src: tuple
    betti_out: tuple
    vertices: int
    repair_rounds: int
    legal: bool = True
    problems: list = field(default_factory=list)


# ------------------------------------------------------------------ helpers

def betti(mask):
    """(b0, b1) with 8-connected foreground / 4-connected background (= objdel.checker.betti):
    b0 by labelling the bounding-box crop, b1 = b0 - Euler characteristic from bit-quads."""
    m, _ = _crop(mask)
    if m.size == 0:
        return (0, 0)
    _, b0 = ndimage.label(m, structure=FG8)
    return (int(b0), int(b0 - bitquads(m).chi8))


def _components(mask):
    lab, n = ndimage.label(mask, structure=FG8)
    return lab, n


def prepare(mask, method, params, G=G_DEFAULT):
    m = np.asarray(mask, bool)
    H, W = m.shape
    lab, _ = _components(m)
    states = []
    simp = METHODS[method]
    for r in polygonize(m):
        path, dirs = lattice_path(r.pts)
        ox, oy = _RIGHT_PIX[int(dirs[0])]
        comp = int(lab[path[0, 1] + oy, path[0, 0] + ox])
        assert comp > 0
        K = np.asarray(simp(path, dirs, W, H, G, **params), np.int64)
        states.append(RingState(path, dirs, np.unique(K), r.hole, comp))
    return states


def _remove_collinear(v, seg):
    """Drop consecutive duplicates and straight-through collinear vertices (cyclic)."""
    changed = True
    while changed and len(v) >= 3:
        changed = False
        keep = np.r_[True, np.any(v[1:] != v[:-1], axis=1)]
        if len(v) > 1 and np.all(v[0] == v[-1]):
            keep[-1] = False
        if not keep.all():
            v, seg = v[keep], seg[keep]
            changed = True
            continue
        p, n = np.roll(v, 1, axis=0), np.roll(v, -1, axis=0)
        cr = (v[:, 0] - p[:, 0]) * (n[:, 1] - p[:, 1]) - (v[:, 1] - p[:, 1]) * (n[:, 0] - p[:, 0])
        dot = (v[:, 0] - p[:, 0]) * (n[:, 0] - v[:, 0]) + (v[:, 1] - p[:, 1]) * (n[:, 1] - v[:, 1])
        straight = (cr == 0) & (dot > 0)
        if straight.any():
            i = int(np.argmax(straight))  # remove one at a time: deterministic and safe
            v, seg = np.delete(v, i, axis=0), np.delete(seg, i)
            changed = True
    return v, seg


def _step_toward(v, w):
    u = w - v
    m = int(np.abs(u).max())
    if m == 0:
        return v.copy()
    # one grid unit (L-inf) along the edge direction, rounded half away from zero
    return v + np.sign(u) * ((2 * np.abs(u) + m) // (2 * m))


def bevel_pinches(rings, segs):
    """Split every grid vertex that occurs more than once within an object by one grid unit
    along each of its two edges. Returns new rings/segs and the number of bevels."""
    from collections import Counter
    cnt = Counter(tuple(p) for r in rings for p in r.tolist())
    dup = {p for p, c in cnt.items() if c > 1}
    if not dup:
        return rings, segs, 0
    out_r, out_s, nb = [], [], 0
    for v, s in zip(rings, segs):
        nv, ns = [], []
        n = len(v)
        for i in range(n):
            if tuple(v[i].tolist()) in dup:
                a = _step_toward(v[i], v[i - 1])
                b = _step_toward(v[i], v[(i + 1) % n])
                nv += [a, b]
                ns += [s[i - 1], s[i]]
                nb += 1
            else:
                nv.append(v[i])
                ns.append(s[i])
        out_r.append(np.array(nv, np.int64).reshape(-1, 2))
        out_s.append(np.array(ns, np.int64))
    return out_r, out_s, nb


def emit(states, W, H, G=G_DEFAULT):
    rings, segs = [], []
    for st in states:
        v = to_grid(st.path[st.K], W, H, G)
        v, seg = _remove_collinear(v, np.arange(len(st.K)))
        rings.append(v)
        segs.append(seg)
    rings, segs, nb = bevel_pinches(rings, segs)
    out_r, out_s = [], []
    for v, s in zip(rings, segs):
        v, s = _remove_collinear(v, s)
        out_r.append(v)
        out_s.append(s)
    return out_r, out_s, nb


def legality_rule(rings, states):
    """Rule L. Returns (ok, problems list, offending (ring, edge) set)."""
    probs, bad = [], set()
    leg = F.legality(rings)
    for k, p in enumerate(leg["ring_problems"]):
        if p:
            probs.append((k, tuple(p)))
        if {"fewer_than_3_points", "zero_area"} & set(p):
            bad.update((k, e) for e in range(len(rings[k])))
    for ra, ea, rb, eb, kind in leg["pairs"]:
        bad.add((ra, ea))
        bad.add((rb, eb))
    if leg["ring_touch"] or leg["ring_cross"]:
        probs.append(("rings", leg["ring_touch"], leg["ring_cross"]))
    for k, (v, st) in enumerate(zip(rings, states)):
        if len(v) < 3:
            continue
        a2 = F.twice_area(v)
        if (a2 > 0) == st.hole or a2 == 0:
            probs.append((k, "orientation"))
            bad.update((k, e) for e in range(len(v)))
    outer_of = {st.comp: k for k, st in enumerate(states) if not st.hole}
    for k, (v, st) in enumerate(zip(rings, states)):
        if st.hole and len(v) >= 3:
            o = outer_of.get(st.comp)
            if o is None or len(rings[o]) < 3 or not F.point_in_ring(v[0], rings[o]):
                probs.append((k, "hole_outside"))
                bad.update((k, e) for e in range(len(v)))
    return not probs, probs, bad


def _simple_lut():
    """256-entry table: is the centre of a 3x3 neighbourhood a simple point for 8-FG / 4-BG?
    Bit k of the code = neighbour k in order NW N NE E SE S SW W (ring order)."""
    ring = [(-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1)]
    lut = np.zeros(256, bool)
    for code in range(256):
        fg = {ring[k] for k in range(8) if code >> k & 1}
        bg = set(ring) - fg

        def comps(cells, adj, seeds):
            seen, n = set(), 0
            for c in cells:
                if c in seen or c not in seeds:
                    continue
                n += 1
                stack = [c]
                while stack:
                    a = stack.pop()
                    if a in seen:
                        continue
                    seen.add(a)
                    for b in cells:
                        if b not in seen and adj(a, b):
                            stack.append(b)
            return n

        a8 = lambda a, b: max(abs(a[0] - b[0]), abs(a[1] - b[1])) == 1  # noqa: E731
        a4 = lambda a, b: abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1  # noqa: E731
        t8 = comps(sorted(fg), a8, fg)
        t4 = comps(sorted(bg), a4, {c for c in bg if abs(c[0]) + abs(c[1]) == 1})
        lut[code] = t8 == 1 and t4 == 1
    return lut


_SIMPLE = _simple_lut()
_RING = [(-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1)]


def nonsimple(mask):
    """Pixels whose flip changes the 8/4 topology of mask (non-simple points), vectorized."""
    p = np.pad(np.asarray(mask, bool), 1)
    H, W = mask.shape
    code = np.zeros((H, W), np.int64)
    for k, (dy, dx) in enumerate(_RING):
        code |= p[1 + dy:1 + dy + H, 1 + dx:1 + dx + W].astype(np.int64) << k
    return ~_SIMPLE[code]


def topo_defects(out, src):
    """Mismatching pixels that are non-simple in the output or the source (where topology broke);
    all mismatches if none is."""
    diff = out != src
    d = diff & (nonsimple(out) | nonsimple(src))
    return d if d.any() else diff


def _seg_index(K, idx):
    """Segment m with K[m] <= idx < K[m+1] (cyclic)."""
    m = int(np.searchsorted(K, idx, side="right")) - 1
    return m % len(K)


def farthest_point(st, m):
    """Farthest skipped lattice index of segment m (middle one among ties), or None for a unit step."""
    K, n = st.K, len(st.path)
    a, b = int(K[m]), int(K[(m + 1) % len(K)])
    L = (b - a) % n if len(K) > 1 else n
    if L < 2:
        return None
    idx = (a + np.arange(1, L)) % n
    P = st.path[idx].astype(np.float64)
    A, B = st.path[a].astype(np.float64), st.path[b].astype(np.float64)
    d = B - A
    nd = np.hypot(*d)
    dist = np.hypot(*(P - A).T) if nd == 0 else np.abs(d[0] * (P[:, 1] - A[1]) - d[1] * (P[:, 0] - A[0])) / nd
    cand = np.nonzero(dist >= dist.max() - 1e-12)[0]
    return int(idx[cand[len(cand) // 2]])


def split_segment(st, m):
    """Insert the farthest skipped lattice point of segment m. Returns True if a point was added."""
    K, n = st.K, len(st.path)
    a, b = int(K[m]), int(K[(m + 1) % len(K)])
    L = (b - a) % n
    if len(K) == 1:
        L = n
    if L < 2:
        return False
    idx = (a + np.arange(1, L)) % n
    P = st.path[idx].astype(np.float64)
    A, B = st.path[a].astype(np.float64), st.path[b].astype(np.float64)
    d = B - A
    nd = np.hypot(*d)
    if nd == 0:
        dist = np.hypot(*(P - A).T)
    else:
        dist = np.abs(d[0] * (P[:, 1] - A[1]) - d[1] * (P[:, 0] - A[0])) / nd
    best = dist.max()
    cand = np.nonzero(dist >= best - 1e-12)[0]
    pick = idx[cand[len(cand) // 2]]
    st.K = np.unique(np.r_[K, pick])
    return True


# ------------------------------------------------------------------ certification

def certify_object(rings, states, src, W, H, G=G_DEFAULT):
    ok_leg, probs, bad_edges = legality_rule(rings, states)
    ras = F.raster(rings, W, H, G) if all(len(r) >= 3 for r in rings) else np.zeros_like(src)
    inter = np.logical_and(ras, src).sum()
    union = np.logical_or(ras, src).sum()
    iou = 1.0 if union == 0 else inter / union
    return ok_leg, probs, bad_edges, ras, float(iou)


@dataclass
class BuildConfig:
    method: str = "minseg_hybrid"  # chosen: hybrid min-# (exact DSS or L-inf tau), tau = 3 grid units
    params: dict = field(default_factory=lambda: {"tau": 3})
    G: int = G_DEFAULT
    iou_min: float = 0.99
    # "auto": topology/legality defects -> pin the farthest skipped point of the offending segments and
    #         re-solve (local); IoU-only failures -> halve the offending ring's tolerance.
    # "retol": always halve the tolerance of offending rings. "split": insert vertices (any simplifier).
    # "pin": pins for everything. "best": auto, and if a repair was needed also retol; keep fewer vertices.
    repair: str = "best"
    max_level: int = 3  # retol: tolerance halvings before the lossless floor
    floor: str = "exact"  # retol floor before exact corners: "exact" (lossless min-#) or "corners"
    max_rounds: int = 60  # split: rounds before fallback to exact rings


_TOL_SCALE = {"tau": 0.5, "t_px": 0.5, "eps": 0.5, "area": 0.25}


def resimplify(st, cfg, W, H):
    """Recompute st.K for the ring's repair level. Returns False if already at exact corners."""
    floor_level = cfg.max_level + 1
    if st.level > floor_level + (cfg.floor == "exact"):
        return False
    if st.level == 0 or (st.level <= cfg.max_level and any(k in cfg.params for k in _TOL_SCALE)):
        p = {k: (v * _TOL_SCALE[k] ** st.level if k in _TOL_SCALE else v) for k, v in cfg.params.items()}
        if cfg.method.startswith("minseg") and st.pins is not None:
            p["pins"] = st.pins
        K = METHODS[cfg.method](st.path, st.dirs, W, H, cfg.G, **p)
        if st.pins is not None and not cfg.method.startswith("minseg"):
            K = np.r_[K, st.pins]
    elif cfg.floor == "exact" and st.level <= floor_level:
        K = METHODS["minseg_exact"](st.path, st.dirs, W, H, cfg.G, pins=st.pins)
    else:
        K = corners_only(st.path, st.dirs)
    st.K = np.unique(np.asarray(K, np.int64))
    return True


def _bad_pixels_to_splits(diff, states_by_obj):
    """Map mismatching pixel centres to (object, ring, segment) by the nearest lattice path point."""
    pts, owner = [], []
    for oi, states in states_by_obj:
        for ri, st in enumerate(states):
            pts.append(st.path.astype(np.float64))
            owner.append(np.stack([np.full(len(st.path), oi), np.full(len(st.path), ri), np.arange(len(st.path))], 1))
    if not pts:
        return set()
    tree = cKDTree(np.concatenate(pts))
    owner = np.concatenate(owner)
    rr, cc = np.nonzero(diff)
    q = np.stack([cc + 0.5, rr + 0.5], 1)
    _, nn = tree.query(q, k=1)
    out = set()
    lookup = dict(states_by_obj)
    for oi, ri, idx in owner[nn]:
        st = lookup[int(oi)][int(ri)]
        m = _seg_index(st.K, idx)
        out.add((int(oi), int(ri), m))
        if idx in st.K:  # pixel nearest to a kept vertex: both adjacent segments are suspects
            out.add((int(oi), int(ri), (m - 1) % len(st.K)))
    return out


def build_item(masks: dict, cfg: BuildConfig, derived=()):
    """masks: {name: bool (H, W)} objects of one item (same crop). derived: [(name, fn(dict of masks) -> mask)]
    checked for topology only (e.g. RIGA rim = disc & ~cup). Returns ({name: ObjectTarget}, info).

    repair="best": run "auto"; if any repair was needed also run "retol" and keep the passing result with
    fewer vertices (ties -> auto). Deterministic."""
    if cfg.repair == "best":
        from dataclasses import replace
        t0 = time.process_time()
        a = _build_item(masks, replace(cfg, repair="auto"), derived)
        if a[1]["rounds"]:
            b = _build_item(masks, replace(cfg, repair="retol"), derived)
            ok = lambda r: all(o.status != "fail" for o in r[0].values())  # noqa: E731
            nv = lambda r: sum(o.vertices for o in r[0].values())  # noqa: E731
            if (ok(b) and not ok(a)) or (ok(b) == ok(a) and nv(b) < nv(a)):
                a = b
        a[1]["time"] = time.process_time() - t0
        return a
    return _build_item(masks, cfg, derived)


def _build_item(masks: dict, cfg: BuildConfig, derived=()):
    t0 = time.process_time()
    names = list(masks)
    H, W = next(iter(masks.values())).shape
    G = cfg.G
    srcs = {k: np.asarray(v, bool) for k, v in masks.items()}
    states = {k: prepare(srcs[k], cfg.method, cfg.params, G) for k in names}
    b_src = {k: betti(srcs[k]) for k in names}
    d_src = {dn: betti(fn(srcs)) for dn, fn in derived}
    fallback = set()
    touched_split = {}
    rounds = 0
    while True:
        emitted, ras, splits, fail, res, tol = {}, {}, set(), set(), {}, set()
        for oi, k in enumerate(names):
            rings, segs, nb = emit(states[k], W, H, G)
            ok_leg, probs, bad_edges, r, iou = certify_object(rings, states[k], srcs[k], W, H, G)
            b_out = betti(r)
            emitted[k], ras[k], res[k] = (rings, segs, nb), r, (ok_leg, probs, iou, b_out)
            if not ok_leg or b_out != b_src[k] or iou < cfg.iou_min:
                fail.add(k)
                for ri, e in bad_edges:
                    s, nK = segs[ri], len(states[k][ri].K)
                    m0, m1 = int(s[e]), int(s[(e + 1) % len(s)])
                    m = m0
                    while True:
                        splits.add((oi, ri, m))
                        m = (m + 1) % nK
                        if m == m1 or m == m0:
                            break
                if b_out != b_src[k]:
                    splits |= _bad_pixels_to_splits(topo_defects(r, srcs[k]), [(oi, states[k])])
                elif iou < cfg.iou_min:
                    tol |= _bad_pixels_to_splits(r != srcs[k], [(oi, states[k])])
        for dn, fn in derived:
            r = fn(ras)
            if betti(r) != d_src[dn]:
                fail.add(dn)
                splits |= _bad_pixels_to_splits(topo_defects(r, fn(srcs)),
                                                [(oi, states[k]) for oi, k in enumerate(names)])
        if rounds == 0:
            raw = {k: {"legal": res[k][0], "topo": res[k][3] == b_src[k], "iou": res[k][2]} for k in names}
            raw.update({dn: {"topo": betti(fn(ras)) == d_src[dn]} for dn, fn in derived})
        if not fail:
            break
        rounds += 1
        progressed = False
        if cfg.repair == "auto":
            # topology / legality defects: pin the farthest skipped point of the offending segments and
            # re-solve the ring (local, keeps the rest optimal); IoU-only failures: halve the ring tolerance
            stuck = set()
            if rounds <= cfg.max_rounds:
                new = {}
                for oi, ri, m in sorted(splits):
                    q = farthest_point(states[names[oi]][ri], m)
                    if q is not None:
                        new.setdefault((oi, ri), set()).add(q)
                stuck = {(oi, ri) for oi, ri, _ in splits}
                for (oi, ri), qs in sorted(new.items()):
                    st = states[names[oi]][ri]
                    old = set() if st.pins is None else set(st.pins.tolist())
                    if qs - old:
                        st.pins = np.array(sorted(old | qs), np.int64)
                        progressed |= resimplify(st, cfg, W, H)
                        touched_split[names[oi]] = True
                        stuck.discard((oi, ri))
            else:
                stuck = {(oi, ri) for oi, ri, _ in splits}
            for oi, ri in sorted(stuck | {(oi, ri) for oi, ri, _ in tol}):
                st = states[names[oi]][ri]
                st.level += 1
                progressed |= resimplify(st, cfg, W, H)
        elif cfg.repair == "pin" and rounds <= cfg.max_rounds:
            # local repair that keeps the rest optimal: pin the farthest skipped point of every offending
            # segment and re-solve the ring; rings that cannot get new pins drop one tolerance level
            new = {}
            for oi, ri, m in sorted(splits | tol):
                q = farthest_point(states[names[oi]][ri], m)
                if q is not None:
                    new.setdefault((oi, ri), set()).add(q)
            stuck = {(oi, ri) for oi, ri, _ in splits | tol}
            for (oi, ri), qs in sorted(new.items()):
                st = states[names[oi]][ri]
                old = set() if st.pins is None else set(st.pins.tolist())
                if qs - old:
                    st.pins = np.array(sorted(old | qs), np.int64)
                    progressed |= resimplify(st, cfg, W, H)
                    touched_split[names[oi]] = True
                    stuck.discard((oi, ri))
            for oi, ri in sorted(stuck):
                st = states[names[oi]][ri]
                st.level += 1
                progressed |= resimplify(st, cfg, W, H)
        elif cfg.repair == "split" and rounds <= cfg.max_rounds:
            splits |= tol
            for oi, ri, m in sorted(splits):
                if names[oi] not in fallback and split_segment(states[names[oi]][ri], m):
                    progressed = True
                    touched_split[names[oi]] = True
        elif cfg.repair == "retol":
            for oi, ri in sorted({(oi, ri) for oi, ri, _ in splits | tol}):
                st = states[names[oi]][ri]
                st.level += 1
                progressed |= resimplify(st, cfg, W, H)
        if not progressed:
            # fallback: exact rings for every failing object (a derived failure -> all objects)
            targets = set(names) if any(k not in srcs for k in fail) else set(fail)
            newly = targets - fallback
            if not newly:
                break  # already exact: report failure
            for k in newly:
                for st in states[k]:
                    st.K = corners_only(st.path, st.dirs)
                    st.level = 99
                fallback.add(k)
    out = {}
    for k in names:
        rings, segs, nb = emitted[k]
        ok_leg, probs, iou, b_out = res[k]
        parts = canonical_parts(rings, [st.hole for st in states[k]], [st.comp for st in states[k]])
        bad = not ok_leg or b_out != b_src[k] or iou < cfg.iou_min
        touched = any(st.level > 0 for st in states[k]) or touched_split.get(k, False)
        status = "fail" if bad else ("fallback" if k in fallback else ("repaired" if touched else "ok"))
        out[k] = ObjectTarget(k, parts, status, iou, b_src[k], b_out, sum(len(r) for r in rings), rounds,
                              ok_leg, probs)
    info = {"time": time.process_time() - t0, "rounds": rounds, "fallback": sorted(fallback),
            "derived_ok": {dn: betti(fn(ras)) == d_src[dn] for dn, fn in derived},
            "bevels": sum(emitted[k][2] for k in names), "raw": raw,
            "levels": {k: max([st.level for st in states[k]] + [0]) for k in names}}
    return out, info


# ------------------------------------------------------------------ canonical form + serialization

def canon_ring(v, hole):
    v = np.asarray(v, np.int64).reshape(-1, 2)
    a2 = F.twice_area(v)
    if (a2 < 0) != hole:  # outer clockwise on screen (y down, positive), hole counter-clockwise
        v = v[::-1]
    i = int(np.lexsort((v[:, 0], v[:, 1]))[0])  # top-most, then left-most
    return np.roll(v, -i, axis=0)


def canonical_parts(rings, holes, comps):
    parts = {}
    for v, h, c in zip(rings, holes, comps):
        parts.setdefault(c, [None, []])
        if h:
            parts[c][1].append(canon_ring(v, True))
        else:
            parts[c][0] = canon_ring(v, False)
    out = []
    for c, (o, hs) in parts.items():
        hs = sorted(hs, key=lambda r: (int(r[0, 1]), int(r[0, 0])))
        out.append((o, hs))
    out.sort(key=lambda p: (int(p[0][0, 1]), int(p[0][0, 0])))
    return out


def _pts(v):
    return [[int(x), int(y)] for x, y in v]


def serialize(objs: dict, schema: dict):
    """objs: {name: parts}; schema: {name: "ring" | "parts"}. Compact, deterministic JSON.
    "ring": the structure must be one part without holes -> [[x,y],...].
    "parts": [{"outer":[[x,y],...],"holes":[[[x,y],...],...]},...]."""
    d = {}
    for k, parts in objs.items():
        if schema[k] == "ring":
            if len(parts) != 1 or parts[0][1]:
                raise ValueError(f"{k}: ring schema needs exactly one part without holes")
            d[k] = _pts(parts[0][0])
        else:
            d[k] = [{"outer": _pts(o), "holes": [_pts(h) for h in hs]} for o, hs in parts]
    return json.dumps(d, separators=(",", ":"))


def parse(text, schema: dict):
    """Inverse of serialize -> {name: [(outer, [holes])]} with (k,2) int64 arrays."""
    d = json.loads(text)
    out = {}
    for k, kind in schema.items():
        v = d[k]
        if kind == "ring":
            out[k] = [(np.array(v, np.int64).reshape(-1, 2), [])]
        else:
            out[k] = [(np.array(p["outer"], np.int64).reshape(-1, 2),
                       [np.array(h, np.int64).reshape(-1, 2) for h in p["holes"]]) for p in v]
    return out


def rasterize_parts(parts, W, H, G=G_DEFAULT):
    rings = [r for o, hs in parts for r in [o, *hs]]
    return F.raster(rings, W, H, G)


def is_canonical(parts):
    for o, hs in parts:
        for v, h in [(o, False), *[(x, True) for x in hs]]:
            if not np.array_equal(v, canon_ring(v, h)):
                return False
        if [tuple(x[0][::-1]) for x in hs] != sorted(tuple(x[0][::-1]) for x in hs):
            return False
    starts = [tuple(o[0][::-1]) for o, _ in parts]
    return starts == sorted(starts)


# ------------------------------------------------------------------ dataset front-ends

RIGA_SCHEMA = {"disc": "ring", "cup": "ring"}
LIDC_SCHEMA = {"nodule": "parts"}
RIGA_DERIVED = (("rim", lambda m: m["disc"] & ~m["cup"]),)


def build_target(masks: dict, schema: dict, cfg: BuildConfig = None, derived=()):
    """Masks of one item (same crop) -> (compact JSON text, {name: ObjectTarget}, info)."""
    cfg = cfg or BuildConfig()
    objs, info = build_item(masks, cfg, derived)
    sch = {k: (v if v == "parts" or (len(objs[k].parts) == 1 and not objs[k].parts[0][1]) else "parts")
           for k, v in schema.items()}
    info["schema"] = sch
    return serialize({k: objs[k].parts for k in schema}, sch), objs, info
