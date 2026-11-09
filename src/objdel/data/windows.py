"""Train windows for scheduled RL rows (ADR-0015; plan B1.1, B1.2): a geometry-only pool of windows per RL item and
the cells of the coverage target.

Cells: (size bin, offset bin) over SIZE_EDGES px and OFFSET_EDGES (SD3's cells). The size is the apparent size
`extent_px * s0 / side`, where `extent_px` is the larger side of the doctors' union box in slice pixels and s0 the
policy side; for centred and shifted windows side == s0. The offset is the distance of the union-box centre from
the window centre over the window side. A value on an edge goes to the upper bin; values are cast to float64 first.

pool(items, draws, seed, exclude): per RL item (two or more doctors) one centred window (the policy frame) and, for
each draw d, one shifted window (augment.draw at the policy side S0) and one rescaled window (augment.neighbour_side:
a neighbouring LADDER side), seeded per item by augment.row_rng (draw 0 with build_eval_views' keys, so draw 0 is
the eval views' window). Every window lies inside the slice, holds every doctor pixel with CLEARANCE px and meets no
other outlined nodule; a draw that breaks a rule is left out and counted, never replaced by another kind or side.
Geometry only: no pixel is read or written. `disagreement` is per slice: the dataset._prepare rule (the doctors'
joint topology answers differ) on the policy frame's masks.
"""
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from operator import attrgetter

import numpy as np

from objdel.geometry import SCHEMAS
from objdel.reward import ExpertCache, rater_answers, rl_eligible
from objdel.targets import augment, frame
from objdel.targets.frame import CLEARANCE, centred, crowded, union_box

SIZE_EDGES = (10, 20, 40)
OFFSET_EDGES = (0.15, 0.2, 0.3)
KINDS = ("centred", "shifted", "rescaled")
DROPS = ("excluded", "not_rl", "shifted_none", "rescaled_none", "padded", "centred_none")


def cell(apparent_px: float, offset: float) -> tuple[int, int]:
    return (int(np.searchsorted(SIZE_EDGES, np.float64(apparent_px), side="right")),
            int(np.searchsorted(OFFSET_EDGES, np.float64(offset), side="right")))


def offset_of(union_box, window_box) -> float:
    """Distance of the two half-open box centres over the window side (the window box's width)."""
    (x0, y0, x1, y1), (wx0, wy0, wx1, wy1) = union_box, window_box
    return float(np.hypot((x0 + x1 - wx0 - wx1) / 2, (y0 + y1 - wy0 - wy1) / 2) / (wx1 - wx0))


@dataclass(frozen=True)
class Window:
    """One train window of a source item: the frame of `side` px moved by (dx, dy) from where the policy centres it."""
    item_id: str
    group_id: str
    nodule: str
    kind: str
    draw: int
    dx: int
    dy: int
    side: int
    s0: int
    offset: float
    extent_px: int
    disagreement: bool

    @property
    def window_id(self) -> str:  # "~", not "/": the source item's nodule and series parse unchanged
        return f"{self.item_id}~{self.kind}{self.draw}"

    @property
    def apparent_px(self) -> float:
        return self.extent_px * self.s0 / self.side

    @property
    def cell(self) -> tuple[int, int]:
        return cell(self.apparent_px, self.offset)


def cell_shares(windows) -> np.ndarray:
    """(4, 4) float64 share of the windows in each cell; ValueError for no windows."""
    if not windows:
        raise ValueError("cell_shares of no windows")
    out = np.zeros((4, 4))
    for w in windows:
        out[w.cell] += 1
    return out / len(windows)


@dataclass(frozen=True)
class WindowPool:
    windows: tuple[Window, ...]  # sorted by window_id
    dropped: Mapping[str, int]  # every DROPS key: items excluded or not RL, draws without an admissible window


def _disagreement(masks, dataset: str) -> bool:
    """dataset._prepare's rule; tau and T5 do not enter the doctors' topology answers, so they are not computed."""
    names = SCHEMAS[dataset]
    return len(rater_answers({n: ExpertCache(masks[n], 1.0, 1.0) for n in names}, names)) > 1


def _windows(item, draws: int, seed: int, dropped: dict) -> list[Window]:
    bare = replace(item, image=None)
    f, ub = frame(bare), union_box(bare)
    H, W = next(iter(item.masks.values())).shape[-2:]
    S0, extent, dis = f.W, max(ub[2] - ub[0], ub[3] - ub[1]), _disagreement(f.masks, item.dataset)
    found = [("centred", 0, (0, 0), S0, f.box)]
    for d in range(draws):
        for kind in ("shifted", "rescaled"):
            rng = augment.row_rng(seed, kind, item.item_id) if d == 0 else augment.row_rng(seed, kind, d, item.item_id)
            got = augment.draw(bare, rng, side=S0) if kind == "shifted" else augment.neighbour_side(bare, S0, rng)
            if got is None:
                dropped[f"{kind}_none"] += 1
                continue
            (dx, dy), S = got
            x, y = centred(bare, S, ub)
            found.append((kind, d, (dx, dy), S, (x + dx, y + dy, x + dx + S, y + dy + S)))
    out = []
    for kind, d, (dx, dy), S, box in found:
        if box[0] < 0 or box[1] < 0 or box[2] > W or box[3] > H:
            dropped["padded"] += 1
        elif kind == "centred" and (crowded(bare, box) or min(ub[0] - box[0], ub[1] - box[1], box[2] - ub[2],
                                                               box[3] - ub[3]) < CLEARANCE):
            dropped["centred_none"] += 1  # the policy frame, moved inside the slice or holding another nodule
        else:
            out.append(Window(item.item_id, item.group_id, item.item_id.rsplit("/", 1)[0], kind, d, dx, dy, S, S0,
                              offset_of(ub, box), extent, dis))
    return out


def pool(items: Iterable, draws: int, seed: int = 0, exclude: frozenset[str] = frozenset()) -> WindowPool:
    """The windows of every RL item not in `exclude` (module docstring), sorted by window_id."""
    dropped, out = dict.fromkeys(DROPS, 0), []
    for it in items:
        if it.item_id in exclude:
            dropped["excluded"] += 1
        elif not rl_eligible(it.raters):
            dropped["not_rl"] += 1
        else:
            out += _windows(it, draws, seed, dropped)
    return WindowPool(tuple(sorted(out, key=attrgetter("window_id"))), dropped)
