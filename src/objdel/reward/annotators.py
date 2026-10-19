"""Scoring a mask against a set of annotator masks: band, topology set, typed errors."""
from collections import Counter
from dataclasses import dataclass

import numpy as np

from ..checker.topology import betti, size_filter


@dataclass
class Score:
    iou: float  # against the strict-majority mask
    band: float  # 1 - (missed core + spill outside union) / |union|
    topo_freq: float  # share of annotators whose (b0, b1) equals the sample's
    topo_near: float  # max over annotators of 1 / (1 + |db0| + |db1|)
    b0: int
    b1: int
    errors: dict


class AnnotatorSet:
    def __init__(self, masks, min_size=0):
        ms = [np.asarray(m, dtype=bool) for m in masks]
        if min_size > 1:
            ms = [size_filter(m, min_size) for m in ms]
        self.masks = ms
        self.k = len(ms)
        self.min_size = min_size
        s = np.sum(ms, axis=0)
        self.core = s == self.k
        self.union = s > 0
        self.vote = 2 * s > self.k
        self.topologies = [betti(m) for m in ms]
        self.hist = Counter(self.topologies)

    @property
    def band_width(self):
        """Share of the union that is ambiguous (not in every annotator's mask)."""
        u = self.union.sum()
        return float((u - self.core.sum()) / u) if u else 0.0

    @property
    def topology_agrees(self):
        return len(self.hist) == 1

    def score(self, mask):
        m = np.asarray(mask, dtype=bool)
        if self.min_size > 1:
            m = size_filter(m, self.min_size)
        inter, uni = np.logical_and(m, self.vote).sum(), np.logical_or(m, self.vote).sum()
        iou = float(inter / uni) if uni else 1.0
        miss_core = int(np.logical_and(self.core, ~m).sum())
        spill = int(np.logical_and(m, ~self.union).sum())
        un = int(self.union.sum())
        band = max(0.0, 1 - (miss_core + spill) / un) if un else float(spill == 0)
        b0, b1 = betti(m)
        freq = self.hist.get((b0, b1), 0) / self.k
        near = max(1 / (1 + abs(b0 - t0) + abs(b1 - t1)) for t0, t1 in self.topologies)
        mode0, mode1 = self.hist.most_common(1)[0][0]
        errors = {
            "extra_objects": max(0, b0 - mode0),
            "missing_objects": max(0, mode0 - b0),
            "extra_holes": max(0, b1 - mode1),
            "missing_holes": max(0, mode1 - b1),
            "missed_core_px": miss_core,
            "spill_px": spill,
        }
        return Score(iou, band, freq, near, b0, b1, errors)


def standardize(v):
    v = np.asarray(v, dtype=float)
    sd = v.std()
    return (v - v.mean()) / (sd if sd > 0 else 1.0)


def group_advantages(scores, w_iou=0.3, w_band=0.3, w_topo=0.4, topo="near", gate=0.5):
    """GRPO-style advantages for one group: standardize each term, mix, re-standardize.

    The topology term is zeroed for samples with IoU below `gate`.
    """
    t = np.array([getattr(s, "topo_" + topo) if s.iou >= gate else 0.0 for s in scores])
    mix = (w_iou * standardize([s.iou for s in scores])
           + w_band * standardize([s.band for s in scores])
           + w_topo * standardize(t))
    return standardize(mix)
