"""Exact per-instance topology of 2D binary masks, and the signature size rule.

Convention (digital topology): 8-connected foreground, 4-connected background. With this pair the digital
Jordan curve theorem holds, so a component's holes are well defined [Rosenfeld 1970; Kong-Rosenfeld 1989].
A hole is a bounded 4-connected background component; the holes of component k are those k has when it is
the only foreground (an island inside a hole of k does not change k's count; the island's holes are its own).

Algorithm, on the CUDA engine (objdel.geom.engine.topology: the reward scorer's pixel label chain, nvcc build on the
box,
CPU build of the same source on the Mac), O(pixels), one call for many masks:
1. Crop to the foreground bounding box (holes cannot leave it) and pad by one background pixel, so the whole
   outside is a single background component.
2. Label the foreground with 8-connectivity and the background with 4-connectivity (union-find, root = the first
   raster pixel; the background root at pixel 0 is the outside, every other background root is a hole), apply the
   size filter, label again.
3. The first raster pixel of a hole has a foreground pixel directly above it (a background pixel there would be
   4-adjacent, hence in the hole, and earlier in raster order), and that pixel belongs to the enclosing
   component: an island inside the hole lies strictly below the hole's top row.
The engine returns per-component numbers only (areas, hole counts, hole areas, signatures), and the filtered pixels
only when asked (size_filter). The OpenCV versions it replaced are tests/oracles/topology_oracle.py (cv2_*). The
bit-quad Euler number (the engine's quads, checker.rings.trace) is an independent certificate:
chi8 = b0 - b1 = sum_k (1 - holes_k).

Input: any 2-D array-like; every nonzero value is foreground (NaN included). `as_mask` is the one normalisation,
used here and by objdel.data.targets. A numpy bool array may hold bytes other than 0 and 1 (np.asarray(PIL "1" image,
bool) stores 255; views of other dtypes store anything), and bool -> bool copies (asarray, astype, ascontiguousarray,
pad) keep those bytes, while every cast to or from bool computes x != 0 [numpy casts]; so a bool array is cast
to uint8 (5 us at 500 x 500, 10x cheaper than a != 0 comparison) and viewed back (issue I1, break-it B1-5, B2-10).

Size rule (plan §4.2): `min_area(|U|) = max(5, ceil(0.5 % of the rater union))`. `size_filter` is an area
opening (drop 8-connected components with < min_area pixels) followed by an area closing (fill holes with
< min_area pixels) [Vincent 1994]. A hole's size is its own background pixel count, measured after the
opening, so removed specks inside it count toward it and an island it encloses does not; filling a small hole
merges the islands inside it into the enclosing component. The result is idempotent: every surviving component
and hole has >= min_area pixels.

References:
- A. Rosenfeld, "Connectivity in digital pictures", J. ACM 17(1):146-160, 1970.
  https://projects.csail.mit.edu/jacm/Authors/rosenfeldazriel.html
- T. Y. Kong, A. Rosenfeld, "Digital topology: introduction and survey", CVGIP 48(3):357-393, 1989.
  https://doi.org/10.1016/0734-189X(89)90147-3
- O. Kalentev, A. Rai, S. Kemnitz, R. Schneider, "Connected component labeling on a 2D grid using CUDA", JPDC
  71(4):615-620, 2011 (union-find labeling on the GPU). https://doi.org/10.1016/j.jpdc.2010.10.012
- L. Vincent, "Morphological area openings and closings for grey-scale images", Shape in Picture, NATO ASI F 126,
  1994 (binary area opening and its dual closing). https://www.doi.org/10.1007/978-3-662-03039-4_13
- numpy casts: _CONVERT_FN is (npy_bool)(x != 0) whenever the source or the destination is bool, in
  https://github.com/numpy/numpy/blob/main/numpy/_core/src/multiarray/lowlevel_strided_loops.c.src ; JAX's
  Array.view warns that viewing as bool gives non-canonical values and recommends comparing with 0:
  https://docs.jax.dev/en/latest/_autosummary/jax.Array.view.html
"""
from collections.abc import Sequence
from dataclasses import dataclass

import cv2
import numpy as np

from objdel.geom import engine as _gpucheck

_EMPTY = np.zeros(0, np.int64)
_EMPTY.flags.writeable = False  # shared by every empty result


@dataclass(frozen=True)
class InstanceTopology:
    """Canonical order: components by (area, holes) descending; hole areas descending."""

    areas: np.ndarray  # (n,) int64, pixels per 8-connected component
    holes: np.ndarray  # (n,) int64, holes per component
    hole_areas: np.ndarray  # (holes.sum(),) int64, background pixels per hole


def as_mask(mask, ndim: int = 2) -> np.ndarray:
    """New C-contiguous bool array of `ndim` dimensions whose bytes are 0 or 1; every nonzero value (any nonzero
    byte of a bool array) is foreground. Raises ValueError for any other number of dimensions."""
    m = np.asarray(mask)
    if m.ndim != ndim:
        raise ValueError(f"expected a {ndim}-D mask, got shape {m.shape}")
    if m.dtype == bool:  # bool -> uint8 casts to x != 0 (module docstring); viewing back is free
        return np.ascontiguousarray(m, np.uint8).view(bool)
    return np.ascontiguousarray(m, bool)  # x != 0 by value


def bbox(m: np.ndarray) -> tuple[int, int, int, int]:
    """Foreground bounding box (x, y, w, h) of an as_mask mask; (0, 0, 0, 0) if it has no foreground."""
    return cv2.boundingRect(m.view(np.uint8))  # box of the non-zero pixels


def _engine(masks: Sequence, min_area, want_masks: bool = False) -> tuple[list[np.ndarray], dict]:
    """One gc.topology call: the as_mask masks and the engine's result (min_area: one value or one per mask)."""
    ms = [as_mask(m) for m in masks]
    start = np.zeros(len(ms) + 1, np.int64)
    start[1:] = np.cumsum([m.size for m in ms])
    bits = np.concatenate([m.view(np.uint8).ravel() for m in ms]) if ms else np.zeros(0, np.uint8)
    hw = np.array([m.shape for m in ms], np.int32).reshape(-1, 2)
    rules = np.array(np.broadcast_to(np.asarray(min_area, np.int64), (len(ms),)))
    return ms, _gpucheck.topology(bits, hw, start, rules, want_masks)


def topologies(masks: Sequence) -> list[InstanceTopology]:
    """instance_topology of each mask, in one engine call."""
    _, out = _engine(masks, 0)
    cs, hs, res = out["comp_start"], out["hole_start"], []
    for i in range(len(cs) - 1):
        a, b = cs[i], cs[i + 1]
        if a == b:
            res.append(InstanceTopology(_EMPTY, _EMPTY, _EMPTY))
            continue
        areas, holes = out["comp_area"][a:b], out["comp_holes"][a:b].astype(np.int64)
        ha = out["hole_area"][hs[a]:hs[b]]
        order = np.lexsort((-holes, -areas))  # canonical, independent of the engine's raster order
        res.append(InstanceTopology(areas[order], holes[order], np.sort(ha)[::-1] if ha.size else _EMPTY))
    return res


def instance_topology(mask) -> InstanceTopology:
    """Area and hole count of every 8-connected component, and the area of every hole."""
    return topologies([mask])[0]


def betti(mask) -> tuple[int, int]:
    """(b0, b1): number of components and of holes."""
    t = instance_topology(mask)
    return int(t.areas.size), int(t.holes.sum())


def min_area(union_area: int) -> int:
    """The size filter of every topology signature: max(5 px, ceil(0.5 % of the rater union's area))."""
    return max(5, -(-int(union_area) // 200))  # integer ceil(union_area / 200)


def filtered(masks: Sequence, min_area) -> list[tuple[np.ndarray, int, int]]:
    """(size_filter(mask, min_area), b0, b1 of it) of each mask, in one engine call (min_area: one value or one per
    mask)."""
    ms, out = _engine(masks, min_area, want_masks=True)
    flat, cs, holes = out["masks"].view(bool), out["comp_start"], out["comp_holes"]
    res, at = [], 0
    for i, m in enumerate(ms):
        res.append((flat[at:at + m.size].reshape(m.shape), int(cs[i + 1] - cs[i]), int(holes[cs[i]:cs[i + 1]].sum())))
        at += m.size
    return res


def size_filter(mask, min_area: int) -> np.ndarray:
    """Remove components, then fill holes, with fewer than min_area pixels (see module docstring); a new array."""
    return filtered([mask], min_area)[0][0]


def signatures(masks: Sequence, min_area: int) -> list[tuple[int, ...]]:
    """signature of each mask at one rule, in one engine call."""
    _, out = _engine(masks, min_area)
    s, sig = out["sig_start"], out["sig"].tolist()
    return [tuple(sig[s[i]:s[i + 1]]) for i in range(len(s) - 1)]


def signature(mask, min_area: int = 0) -> tuple[int, ...]:
    """Holes per component, descending, after `size_filter(mask, min_area)`: (0,) = one simply connected part,
    (0, 0) = two parts, (1,) = one part with one hole, () = empty."""
    return signatures([mask], min_area)[0]
