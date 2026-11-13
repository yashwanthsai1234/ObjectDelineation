"""Train windows for scheduled RL rows (ADR-0015): a geometry-only pool of windows per RL item and
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
Geometry only: no pixel is read or written. `disagreement` is per slice: the prepare._prepare rule (the doctors'
joint topology answers differ) on the policy frame's masks.

schedule(windows, target, steps, ...): `steps * per_step` windows in step-major order, built in blocks of
`block_steps` steps, so every prefix of whole blocks holds exact counts:
1. kind counts per block = largest-remainder rounding (ties to the lower index) of shares x block rows, dealt over
   the block's steps by the same rounding of what is left (each step gets 4 or 5 centred rows at 32 per step);
2. per off-centre kind with n block slots, cell counts = the rounding of target x n; centred slots have no cell;
3. per (kind, cell) slot set of m slots, m_dis = the rounding of q x m disagreement slots, q = w p / (1 + (w - 1) p),
   p the disagreement share of that (kind, cell) in the pool (w = 1: q = p). Cell and disagreement labels go to the
   kind's slot positions by one shuffle with default_rng((seed, block));
4. each slot takes the next window of a seeded permutation of its (kind, cell, disagreement) windows whose item is
   not yet in the step. A window is used once in the run, except that a disagreement slot set that runs out starts
   another permutation of its windows (at most MAX_USES uses of one window, as R5 repeated each disagreement row 4
   times). ValueError names the kind, the cell and the part when a slot cannot be filled; no other kind or cell
   stands in. Nodule weighting is off.
"""
from bisect import bisect_right
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from operator import attrgetter

import numpy as np

from objdel.core import datasets
from objdel.data.targets import augment
from objdel.data.targets.frame import centred, clear, crowded, extent, frame, padded, union_box
from objdel.verify.reward.cache import ExpertCache, rl_eligible
from objdel.verify.reward.score import rater_answers

SIZE_EDGES = (10, 20, 40)
OFFSET_EDGES = (0.15, 0.2, 0.3)
KINDS = ("centred", "shifted", "rescaled")
CELLS = tuple((a, b) for a in range(4) for b in range(4))  # row-major, as target.reshape(-1)
SHARES: Mapping[str, float] = {"centred": 0.15, "shifted": 0.75, "rescaled": 0.10}  # config values [estimate]
MAX_USES = 4  # uses of one disagreement window in a run
DROPS = ("excluded", "not_rl", "shifted_none", "rescaled_none", "padded", "centred_none")


def cell(apparent_px: float, offset: float) -> tuple[int, int]:
    """The (size bin, offset bin) cell of an apparent size in pixels and a window offset."""
    if not (np.isfinite(apparent_px) and np.isfinite(offset)):
        raise ValueError(f"cell: apparent_px and offset must be finite, got {apparent_px}, {offset}")
    return bisect_right(SIZE_EDGES, float(apparent_px)), bisect_right(OFFSET_EDGES, float(offset))


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
        """Id of the window: the source item id, `~`, the kind and the draw number."""
        return f"{self.item_id}~{self.kind}{self.draw}"

    @property
    def apparent_px(self) -> float:
        """Apparent size of the union in pixels: `extent_px * s0 / side`."""
        return self.extent_px * self.s0 / self.side

    @property
    def cell(self) -> tuple[int, int]:
        """The (size bin, offset bin) cell of this window."""
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
    """The train window pool: windows sorted by id and the count of everything dropped, by reason."""
    windows: tuple[Window, ...]  # sorted by window_id
    dropped: Mapping[str, int]  # every DROPS key: items excluded or not RL, draws without an admissible window


def _disagreement(masks, dataset: str) -> bool:
    """prepare._prepare's rule; tau and T5 do not enter the doctors' topology answers, so they are not computed."""
    names = datasets.get(dataset).names
    return len(rater_answers({n: ExpertCache(masks[n], 1.0, 1.0) for n in names}, names)) > 1


def _windows(item, draws: int, seed: int, dropped: dict) -> list[Window]:
    bare = replace(item, image=None)
    f, ub = frame(bare), union_box(bare)
    H, W = next(iter(item.masks.values())).shape[-2:]
    S0, ext, dis = f.W, extent(ub), _disagreement(f.masks, item.dataset)
    found = [("centred", 0, (0, 0), S0, f.box)]
    for d in range(draws):
        for kind in ("shifted", "rescaled"):
            rng = augment.row_rng(seed, kind, item.item_id) if d == 0 else augment.row_rng(seed, kind, d, item.item_id)
            got = (augment.draw(bare, rng, side=S0) if kind == "shifted"
                   else augment.neighbour_side(bare, S0, rng, box=ub))
            if got is None:
                dropped[f"{kind}_none"] += 1
                continue
            (dx, dy), S = got
            x, y = centred(bare, S, ub)
            found.append((kind, d, (dx, dy), S, (x + dx, y + dy, x + dx + S, y + dy + S)))
    out = []
    for kind, d, (dx, dy), S, box in found:
        if padded(box, W, H):
            dropped["padded"] += 1
        elif kind == "centred" and (crowded(bare, box) or not clear(ub, box)):
            dropped["centred_none"] += 1  # the policy frame, moved inside the slice or holding another nodule
        else:
            out.append(Window(item.item_id, item.group_id, item.item_id.rsplit("/", 1)[0], kind, d, dx, dy, S, S0,
                              offset_of(ub, box), ext, dis))
    return out


def pool(items: Iterable, draws: int, seed: int = 0, exclude: frozenset[str] = frozenset()) -> WindowPool:
    """The windows of every RL item not in `exclude` (module docstring), sorted by window_id."""
    dropped, out, seen = dict.fromkeys(DROPS, 0), [], set()
    for it in items:
        if it.item_id in exclude or it.item_id in seen:  # a repeated item counts as excluded
            dropped["excluded"] += 1
        elif not rl_eligible(it.raters):
            dropped["not_rl"] += 1
        else:
            out += _windows(it, draws, seed, dropped)
        seen.add(it.item_id)
    return WindowPool(tuple(sorted(out, key=attrgetter("window_id"))), dropped)


def _lrr(x, n: int) -> np.ndarray:
    """Largest-remainder rounding of x (summing to n) to integers summing to n; ties go to the lower index."""
    x = np.asarray(x, np.float64)
    out = np.floor(x).astype(np.int64)
    out[np.argsort(out - x, kind="stable")[:n - out.sum()]] += 1
    return out


def _fractions(name, values, n=None) -> np.ndarray:
    v = np.asarray(values, np.float64)
    if (n and v.shape != n) or not np.isfinite(v).all() or (v < 0).any() or abs(v.sum() - 1) > 1e-9:
        raise ValueError(f"{name} must be {n or ''} non-negative finite shares summing to 1, got {values}")
    return v


def schedule(windows: Sequence[Window], target: np.ndarray, steps: int, per_step: int = 32,
             shares: Mapping[str, float] = SHARES, seed: int = 0, disagreement_weight: float = 4.0,
             block_steps: int = 20) -> list[Window]:
    """The run's rows (module docstring); `target` is the (4, 4) cell law of the off-centre kinds."""
    t = _fractions("target", target, (4, 4)).reshape(-1)
    if set(shares) != set(KINDS):
        raise ValueError(f"shares need exactly the kinds {KINDS}, got {sorted(shares)}")
    sh = _fractions("shares", [shares[k] for k in KINDS])
    w = float(disagreement_weight)
    if min(steps, per_step, block_steps) <= 0 or steps % block_steps or not w > 0 or not np.isfinite(w):
        raise ValueError(f"need steps a positive multiple of block_steps ({steps}, {block_steps}), per_step > 0 "
                         f"({per_step}) and a finite disagreement_weight > 0 ({disagreement_weight})")
    groups: dict[tuple, list[Window]] = {}
    for x in sorted(windows, key=attrgetter("window_id")):
        groups.setdefault((x.kind, None if x.kind == "centred" else x.cell, x.disagreement), []).append(x)
    pick = np.random.default_rng(seed)
    queues: dict[tuple, list[Window]] = {}
    rounds: dict[tuple, int] = {}

    def take(key, step_items):
        q = queues.setdefault(key, [])
        while True:
            for j, x in enumerate(q):
                if x.item_id not in step_items:
                    step_items.add(x.item_id)
                    return q.pop(j)
            if rounds.get(key, 0) >= (MAX_USES if key[2] else 1) or key not in groups:
                kind, c, dis = key
                raise ValueError(f"no free {kind} window in cell {c if c else 'any'} for a "
                                 f"{'disagreement' if dis else 'other'} slot (pool too small: raise draws)")
            rounds[key] = rounds.get(key, 0) + 1
            q += [groups[key][i] for i in pick.permutation(len(groups[key]))]

    nb, out = block_steps * per_step, []
    for b in range(steps // block_steps):
        left = _lrr(sh * nb, nb)
        kinds = []
        for s in range(block_steps):
            now = _lrr(left / (block_steps - s), per_step)
            left -= now
            kinds += [k for k, c in zip(KINDS, now, strict=True) for _ in range(c)]
        labels, rng = [None] * nb, np.random.default_rng((seed, b))
        for kind in KINDS:
            at = [i for i, k in enumerate(kinds) if k == kind]
            counts = {None: len(at)} if kind == "centred" else dict(zip(CELLS, _lrr(t * len(at), len(at)), strict=True))
            tags = []
            for c, m in ((c, int(m)) for c, m in counts.items() if m):
                n_dis, n_other = (len(groups.get((kind, c, d), ())) for d in (True, False))
                if not n_dis + n_other:
                    raise ValueError(f"no {kind} window in cell {c} for its {m} slots (pool too small: raise draws)")
                p = n_dis / (n_dis + n_other)
                q = w * p / (1 + (w - 1) * p)
                k_dis = int(_lrr(np.array([q, 1 - q]) * m, m)[0])
                tags += [(kind, c, True)] * k_dis + [(kind, c, False)] * (m - k_dis)
            for i, j in zip(at, rng.permutation(len(at)), strict=True):
                labels[i] = tags[j]
        for s in range(block_steps):
            step_items: set[str] = set()
            out += [take(labels[i], step_items) for i in range(s * per_step, (s + 1) * per_step)]
    return out
