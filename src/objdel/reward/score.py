"""Reward-v4 terms of one rollout: completion text -> parse -> the CUDA checker scores it against the expert caches.

All checking and scoring runs in the CUDA checker (objdel._gpucheck.score_polygons, cuda/gpucheck.cpp): polygon ->
even-odd mask -> ring defects -> the terms below, against each cache's tables resident on the device
(ExpertCache.device). Masks never leave the device; only per-object numbers come back. One checker call per
score_many / score_batch. Per object, against its ExpertCache:
- IoU_tc: max IoU over the raters whose size-filtered signature equals the prediction's; IoU vs the vote
  when no rater drew that signature (lab ITERATION3, R3);
- F5 = 1 / (1 + (2 D / P) / T5^2): D = depth beyond the egg-yolk band, P = contour length of the union;
- g_area = min(1, A_max / A) * min(1, A / A_min), A_min / A_max over all raters, zero tolerance (ITERATION4);
  an empty prediction gets 0;
- the signature and its rater frequency;
- loc_raw (reward v5, Implementation.md B3.2): the location kernel of the mask's centroid (the checker's exact
  sum_x, sum_y) against the vote's centroid, scaled by the union's radius (ExpertCache.ref_sum, radius2).
The terms are bit-identical to the host scorer this replaced (tests/oracles/score_numpy.py; the exactness argument is
in docs/research/notes/gpu_scoring.md).

A completion with several objects (RIGA: disc, cup) gets the mean of the per-object IoU_tc, F5 and g_area; it
is legal when every object is legal; its topology answer is joint: topo_freq is the share of the voting raters
(`rater_answers`; ExpertCache docstring on abstaining raters) whose signatures equal the prediction's for every
object (rater r of one object's cache is rater r of the others; caches with different rater counts are refused),
and sig_id identifies the tuple of per-object signatures. This is the lab's rule (t_stack/prod.py score_text);
for one object it reduces to the per-object terms. iou_min, the smallest per-object IoU_tc, feeds the topology
validity gate (advantages_v4), so a good disc cannot carry a missed cup past it (B2-7; the eval pass criterion).
The completion's loc is the mean over objects of gated_loc: non-zero only for a missed answer (IoU_tc < miss_iou)
of doctor-like size (g_area) and a topology some doctor drew (topo_freq > 0). These combinations of the returned
per-object numbers are the only host arithmetic.

An unparseable or truncated completion has format_ok = legal = 0, every term 0 and sig_id -1. A parsed but
illegal one keeps its terms for logging; advantages_v4 zeroes them through the legality gate. An empty object
list is a legal empty mask, scored like any other mask.

References:
- Reward v4 definition: archive/phase1_lab/t_reward/results/ITERATION4.md; plan §3.2, §4.3.
- VLM-R1's count guard s_ovd = min(1, L_gt / L_pred) multiplying the detection reward against over-prediction
  (§4.3, Eq. 7-8; the hack it removes: §5.3.1), the analogue of g_area: https://arxiv.org/abs/2504.07615 (full text
  https://arxiv.org/html/2504.07615)
"""
import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from objdel import _gpucheck
from objdel.geometry import GRID, parse
from objdel.geometry.shape import pack_objects

from .cache import ExpertCache, device_many
from .spec import REWARD_V5

NUMERIC = ("format_ok", "legal", "iou_tc", "iou_min", "f5", "g_area", "topo_any", "topo_freq", "sig_id", "loc")
_FAILED = dict.fromkeys(NUMERIC, 0.0) | {"sig_id": -1.0}


@dataclass(frozen=True)
class ObjectTerms:
    """Per-object reward terms from the device checker: IoU_tc, F5, area term, signature, its voter share and area in
    pixels.
    """
    iou_tc: float
    f5: float
    g_area: float
    sig: tuple[int, ...]  # size-filtered signature: holes per component, descending
    topo_freq: float  # share of the voting raters who drew `sig`
    area: int
    sum_x: int  # sum of the column indices of the mask's pixels
    sum_y: int  # sum of the row indices
    loc_raw: float  # location kernel of the centroid (`location`)


def location(sum_x: int, sum_y: int, area: int, ref_sum: tuple[int, int, int], radius2: float,
             width: float = REWARD_V5.loc_width) -> float:
    """exp(-d2 / (2 width^2 radius2)), d2 the squared distance of the mask's centroid to the reference centroid
    (pixel indices on both sides); 0 for an empty mask or an empty reference."""
    rx, ry, n = ref_sum
    if area == 0 or n == 0:
        return 0.0
    d2 = (sum_x / area - rx / n) ** 2 + (sum_y / area - ry / n) ** 2
    return float(np.exp(-d2 / (2 * width**2 * radius2)))


def gated_loc(iou_tc: float, g_area: float, topo_freq: float, loc_raw: float,
              miss_iou: float = REWARD_V5.miss_iou) -> float:
    """An object's location term: loc_raw on a missed answer, ramped to 0 at IoU_tc = miss_iou (R stays continuous
    there), times g_area (no specks, no huge outlines) and only for a topology a doctor drew (no blob pairs around
    the centroid)."""
    return (iou_tc < miss_iou) * max(0.0, 1 - iou_tc / miss_iou) * g_area * (topo_freq > 0) * loc_raw


def _object(r: dict, i: int, cache: ExpertCache) -> ObjectTerms:
    """Object i of a device result, against its cache."""
    sig = tuple(r["sig"][r["sig_start"][i]:r["sig_start"][i + 1]].tolist())
    sx, sy, area = int(r["sum_x"][i]), int(r["sum_y"][i]), int(r["area"][i])
    return ObjectTerms(float(r["iou_tc"][i]), float(r["f5"][i]), float(r["g_area"][i]), sig,
                       float(r["topo_freq"][i]), area, sx, sy, location(sx, sy, area, cache.ref_sum, cache.radius2))


def score_mask(mask: np.ndarray, cache: ExpertCache) -> ObjectTerms:
    """Terms of one object's mask (H, W) against its cache, on the device (for masks that do not come from a
    completion: baselines, reward gates)."""
    m = np.ascontiguousarray(mask, dtype=bool)
    if m.shape != cache.union.shape:
        raise ValueError(f"mask shape {m.shape} != cache shape {cache.union.shape}")
    H, W = m.shape
    r = _gpucheck.score_masks(m.view(np.uint8).reshape(-1), np.array([[H, W]], np.int32),
                              np.array([0, H * W], np.int64), np.array([cache.device()], np.int64))
    return _object(r, 0, cache)


def sig_id(sigs: tuple) -> float:
    """Deterministic id of a (joint) signature: CRC32, independent of PYTHONHASHSEED; exact in float64."""
    return float(zlib.crc32(repr(sigs).encode()))


def rater_answers(caches: Mapping[str, ExpertCache], names: Sequence[str]) -> dict[tuple, list[int]]:
    """The raters' joint topology answers {(signature per object in `names` order): [raters]}, over the raters who
    vote in every object (all raters if none does). ValueError if the caches have different rater counts."""
    ks = {n: caches[n].k for n in names}
    if len(set(ks.values())) != 1:
        raise ValueError(f"the caches of one item need the same raters, got rater counts {ks}")
    k = ks[names[0]]
    vote = [r for r in range(k) if all(r in caches[n].voters for n in names)] or range(k)
    out: dict[tuple, list[int]] = {}
    for r in vote:
        out.setdefault(tuple(caches[n].rater_sigs[r] for n in names), []).append(r)
    return out


def score_completion(text: str, truncated: bool, names: Sequence[str], caches: Mapping[str, ExpertCache],
                     W: int, H: int) -> dict:
    """One completion -> {format_ok, legal, iou_tc, iou_min, f5, g_area, topo_any, topo_freq, sig_id, loc, error,
    defects}."""
    return score_many([(text, truncated, names, caches, W, H)])[0][0]


def score_many(jobs: Sequence[tuple]) -> list[tuple[dict, dict[str, ObjectTerms]]]:
    """(score_completion's terms, {object: ObjectTerms}, empty when the completion does not parse) of every job
    (tuples of score_completion's arguments), in order: one engine call that builds the jobs' cold caches
    (device_many), and one CUDA checker call for all parsed objects."""
    device_many([c for _, _, names, caches, _, _ in jobs for c in (caches[n] for n in names)])
    parsed = [parse(text, names, truncated) for text, truncated, names, _, _, _ in jobs]
    objs, tables = [], []
    for p, (_, _, names, caches, W, H) in zip(parsed, jobs, strict=True):
        if p.ok:
            objs += [(p.objs[n], W, H) for n in names]
            tables += [caches[n].device() for n in names]
    r = _gpucheck.score_polygons(*pack_objects(objs), GRID, np.array(tables, np.int64))
    out, i = [], 0
    for p, (_, _, names, caches, _, _) in zip(parsed, jobs, strict=True):
        a = rater_answers(caches, names)  # checks the caches of every job, parsed or not
        if not p.ok:
            out.append((_FAILED | {"error": p.error, "defects": {}}, {}))
            continue
        out.append(_combine(names, caches, a, r, range(i, i + len(names))))
        i += len(names)
    return out


def _combine(names, caches, answers, r: dict, idx) -> tuple[dict, dict[str, ObjectTerms]]:
    """A parsed completion's terms from its objects' rows idx of the device result r."""
    legal, defects, terms = True, {}, {}
    for n, i in zip(names, idx, strict=True):
        for d, c in zip(_gpucheck.DEFECTS, r["defects"][i].tolist(), strict=True):
            if c:
                legal, defects[d] = False, defects.get(d, 0) + c
        terms[n] = _object(r, i, caches[n])
    ts = list(terms.values())
    sigs = tuple(t.sig for t in ts)
    joint = len(answers.get(sigs, ())) / sum(map(len, answers.values()))
    return ({"format_ok": 1.0, "legal": float(legal), "iou_tc": float(np.mean([t.iou_tc for t in ts])),
             "iou_min": float(min(t.iou_tc for t in ts)), "f5": float(np.mean([t.f5 for t in ts])),
             "g_area": float(np.mean([t.g_area for t in ts])), "topo_any": float(joint > 0), "topo_freq": joint,
             "sig_id": sig_id(sigs),
             "loc": float(np.mean([gated_loc(t.iou_tc, t.g_area, t.topo_freq, t.loc_raw) for t in ts])),
             "error": "", "defects": defects}, terms)


def score_batch(jobs: Sequence[tuple]) -> dict:
    """score_completion over jobs (tuples of its arguments) in input order, in one CUDA checker call. Returns float64
    arrays for the numeric keys and lists for "error" and "defects"."""
    rows = [t for t, _ in score_many(list(jobs))]
    out = {k: np.array([r[k] for r in rows], np.float64) for k in NUMERIC}
    return out | {"error": [r["error"] for r in rows], "defects": [r["defects"] for r in rows]}
