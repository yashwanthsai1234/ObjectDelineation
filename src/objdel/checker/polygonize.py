"""Exact mask -> rectilinear rings, linked with array operations only.

This mirrors the planned GPU linker: emit unit boundary edges, link each edge to
its successor, find rings by min-label pointer jumping, order each ring by list
ranking from a canonical start. No sequential tracing loop.

Edges run clockwise on screen (y down) around foreground, so outer rings have
positive shoelace area and holes negative. At a saddle vertex (codes 6 and 9)
the walk turns left, which joins diagonal foreground pixels: 8-connected
foreground, 4-connected background, matching checker.topology.

The edges are the cracks between pixels (1-cells of the cell complex of the pixel grid) [Kovalevsky 1989]. A saddle
vertex is the ambiguous case of marching squares; turning left resolves it toward joining the foreground, which is
the (8, 4) pairing of digital topology (the ring then touches itself at that vertex) [Kong-Rosenfeld 1989].
ring_labels is min-label pointer jumping on the successor permutation, as in parallel connectivity
[Shiloach-Vishkin 1982; Hillis-Steele 1986]; list_rank is Wyllie's list ranking [Wyllie 1979]. Both take
ceil(log2 n) + 1 rounds. The signed area is the shoelace sum over the edges.

References:
- V. A. Kovalevsky, "Finite topology as applied to image analysis", CVGIP 46(2):141-161, 1989.
  https://doi.org/10.1016/0734-189X(89)90165-5
- T. Y. Kong, A. Rosenfeld, "Digital topology: introduction and survey", CVGIP 48(3):357-393, 1989.
  https://doi.org/10.1016/0734-189X(89)90147-3
- Marching squares, saddle disambiguation: https://en.wikipedia.org/wiki/Marching_squares
- J. C. Wyllie, "The complexity of parallel computations", Cornell TR79-387, 1979, ch. 4 (list ranking).
  https://hdl.handle.net/1813/7502
- Y. Shiloach, U. Vishkin, "An O(log n) parallel connectivity algorithm", J. Algorithms 3(1):57-67, 1982.
  https://doi.org/10.1016/0196-6774(82)90008-6
- W. D. Hillis, G. L. Steele, "Data parallel algorithms", CACM 29(12):1170-1183, 1986 (pointer doubling).
  https://doi.org/10.1145/7902.7903
- Pointer jumping (list ranking, root finding): https://en.wikipedia.org/wiki/Pointer_jumping
- Shoelace formula: https://en.wikipedia.org/wiki/Shoelace_formula
"""
from dataclasses import dataclass

import numpy as np

# Directions: 0 east, 1 south, 2 west, 3 north.
DX = np.array([1, 0, -1, 0])
DY = np.array([0, 1, 0, -1])


@dataclass
class Ring:
    pts: np.ndarray  # (k, 2) int64 corner vertices (x, y), canonical start first
    area: int  # signed: > 0 outer, < 0 hole

    @property
    def hole(self):
        return self.area < 0


def boundary_edges(mask):
    """Start vertex (x, y) and direction of every unit boundary edge."""
    m = np.asarray(mask, dtype=bool)
    p = np.pad(m, 1)
    up, down = ~p[:-2, 1:-1], ~p[2:, 1:-1]
    left, right = ~p[1:-1, :-2], ~p[1:-1, 2:]
    xs, ys, ds = [], [], []
    for side, (ox, oy), d in ((up, (0, 0), 0), (right, (1, 0), 1), (down, (1, 1), 2), (left, (0, 1), 3)):
        y, x = np.nonzero(m & side)
        xs.append(x + ox)
        ys.append(y + oy)
        ds.append(np.full(x.size, d))
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(ds)


def _jumps(n):
    return max(1, int(np.ceil(np.log2(max(n, 2)))) + 1)


def link(x, y, d, width):
    """Successor of every edge. At saddles take the left turn."""
    stride = width + 1
    start = y * stride + x
    end = (y + DY[d]) * stride + (x + DX[d])
    order = np.argsort(start, kind="stable")
    sorted_start = start[order]
    first = np.searchsorted(sorted_start, end, side="left")
    count = np.searchsorted(sorted_start, end, side="right") - first
    assert count.min() >= 1 and count.max() <= 2, "boundary vertex with 0 or >2 outgoing edges"
    succ = order[first]
    two = np.nonzero(count == 2)[0]
    if two.size:
        a, b = order[first[two]], order[first[two] + 1]
        want = (d[two] + 3) % 4
        assert np.all((d[a] == want) | (d[b] == want)), "no left turn at saddle"
        succ[two] = np.where(d[a] == want, a, b)
    return succ


def ring_labels(succ):
    """Min edge index on each cycle, by pointer jumping."""
    lab = np.arange(succ.size)
    nxt = succ.copy()
    for _ in range(_jumps(succ.size)):
        lab = np.minimum(lab, lab[nxt])
        nxt = nxt[nxt]
    return lab


def list_rank(succ, start_mask):
    """Steps from the start edge of each cycle, by Wyllie pointer jumping."""
    n = succ.size
    pred = np.empty(n, dtype=np.int64)
    pred[succ] = np.arange(n)
    nxt = succ.copy()
    dist = np.ones(n, dtype=np.int64)
    term = pred[start_mask]
    nxt[term] = term
    dist[term] = 0
    for _ in range(_jumps(n)):
        dist = dist + dist[nxt]
        nxt = nxt[nxt]
    return dist, pred


def polygonize(mask, offset=(0, 0)):
    """Rings of a mask. offset (x, y) is added to vertices, e.g. after cropping."""
    m = np.asarray(mask, dtype=bool)
    H, W = m.shape
    x, y, d = boundary_edges(m)
    if x.size == 0:
        return []
    succ = link(x, y, d, W)
    lab = ring_labels(succ)
    ring_ids, ring = np.unique(lab, return_inverse=True)
    R = ring_ids.size

    key = (y * (W + 1) + x) * 4 + d
    minkey = np.full(R, np.iinfo(np.int64).max)
    np.minimum.at(minkey, ring, key)
    is_start = key == minkey[ring]

    dist_to_end, pred = list_rank(succ, is_start)
    length = np.bincount(ring, minlength=R)
    rank = length[ring] - 1 - dist_to_end

    x1, y1 = x + DX[d], y + DY[d]
    area2 = np.bincount(ring, weights=x * y1 - x1 * y, minlength=R).astype(np.int64)

    corner = d != d[pred]
    order = np.lexsort((rank, ring))
    order = order[corner[order]]
    rings_sorted = ring[order]
    cuts = np.searchsorted(rings_sorted, np.arange(R + 1))
    out = []
    for r in range(R):
        idx = order[cuts[r]:cuts[r + 1]]
        pts = np.stack([x[idx] + offset[0], y[idx] + offset[1]], axis=1)
        out.append(Ring(pts, int(area2[r] // 2)))
    return out
