"""Reward-v4 terms of one rollout: completion text -> parse -> polycore mask -> terms against the expert caches.

Per object, against its ExpertCache (all work on the bounding box of prediction | rater union):
- IoU_tc: max IoU over the raters whose size-filtered signature equals the prediction's; IoU vs the vote
  when no rater drew that signature (lab ITERATION3, R3);
- F5 = 1 / (1 + (2 D / P) / T5^2): D = depth beyond the egg-yolk band, P = contour length of the union;
- g_area = min(1, A_max / A) * min(1, A / A_min), A_min / A_max over all raters, zero tolerance (ITERATION4);
  an empty prediction gets 0;
- the signature and its rater frequency.

A completion with several objects (RIGA: disc, cup) gets the mean of the per-object IoU_tc, F5 and g_area; it
is legal when every object is legal; its topology answer is joint: topo_freq is the share of the voting raters
(`rater_answers`; ExpertCache docstring on abstaining raters) whose signatures equal the prediction's for every
object (rater r of one object's cache is rater r of the others; caches with different rater counts are refused),
and sig_id identifies the tuple of per-object signatures. This is the lab's rule (t_stack/prod.py score_text);
for one object it reduces to the per-object terms. iou_min, the smallest per-object IoU_tc, feeds the topology
validity gate (advantages_v4), so a good disc cannot carry a missed cup past it (B2-7; the eval pass criterion).

An unparseable or truncated completion has format_ok = legal = 0, every term 0 and sig_id -1. A parsed but
illegal one keeps its terms for logging; advantages_v4 zeroes them through the legality gate. An empty object
list is a legal empty mask (objdel.geometry.to_mask), scored like any other mask.

References:
- Reward v4 definition: experiments/phase1/lab/t_reward/results/ITERATION4.md; plan §3.2, §4.3.
- Python threads + GIL-releasing kernels (polycore, OpenCV): https://docs.python.org/3/library/concurrent.futures.html
- OpenCV releases the GIL in every wrapped call (ERRWRAP2 / PyAllowThreads):
  https://github.com/opencv/opencv/blob/4.x/modules/python/src2/cv2_util.hpp ; setNumThreads semantics (GCD backend
  accepts only values <= 0): https://github.com/opencv/opencv/blob/4.x/modules/core/include/opencv2/core/utility.hpp
- VLM-R1's count guard s_ovd = min(1, L_gt / L_pred) multiplying the detection reward against over-prediction
  (§4.3, Eq. 7-8; the hack it removes: §5.3.1), the analogue of g_area: https://arxiv.org/abs/2504.07615 (full text
  https://arxiv.org/html/2504.07615)
"""
import zlib
from collections.abc import Mapping, Sequence
from concurrent.futures import Executor
from dataclasses import dataclass

import cv2
import numpy as np

from objdel.checker.topology import signature
from objdel.geometry import parse, to_masks

from .cache import ExpertCache, band_depth
from .spec import REWARD_V4

NUMERIC = ("format_ok", "legal", "iou_tc", "iou_min", "f5", "g_area", "topo_any", "topo_freq", "sig_id")
_FAILED = dict.fromkeys(NUMERIC, 0.0) | {"sig_id": -1.0}


@dataclass(frozen=True)
class ObjectTerms:
    iou_tc: float
    f5: float
    g_area: float
    sig: tuple[int, ...]  # size-filtered signature: holes per component, descending
    topo_freq: float  # share of the voting raters who drew `sig`
    area: int


def _iou(a: np.ndarray, area_a: int, b: np.ndarray, area_b: int) -> float:
    inter = int(np.count_nonzero(a & b))
    union = area_a + area_b - inter
    return inter / union if union else 1.0


def area_guard(area: int, cache: ExpertCache) -> float:
    """min(1, A_max / A) * min(1, A / A_min) over all raters; 0 for an empty prediction."""
    if area <= 0:
        return 0.0
    a_min, a_max = min(cache.rater_areas), max(cache.rater_areas)
    return min(1.0, a_max / area) * (min(1.0, area / a_min) if a_min > 0 else 1.0)


def score_mask(mask: np.ndarray, cache: ExpertCache) -> ObjectTerms:
    """Terms of one object's mask (H, W) against its cache."""
    m = np.ascontiguousarray(mask, dtype=bool)
    if m.shape != cache.union.shape:
        raise ValueError(f"mask shape {m.shape} != cache shape {cache.union.shape}")
    x, y, w, h = cv2.boundingRect(m.view(np.uint8))
    boxes = [b for b in (cache.box, (y, y + h, x, x + w)) if b[1] > b[0]]  # every term lives in pred | union
    roi = (slice(min(b[0] for b in boxes), max(b[1] for b in boxes)),
           slice(min(b[2] for b in boxes), max(b[3] for b in boxes))) if boxes else (slice(0, 0), slice(0, 0))
    m = m[roi]
    area = int(np.count_nonzero(m))
    sig = signature(m, cache.min_area)
    same = cache.same_sig.get(sig, ())
    if same:
        iou_tc = max(_iou(m, area, cache.raters[i][roi], cache.rater_areas[i]) for i in same)
    else:
        iou_tc = _iou(m, area, cache.vote[roi], cache.n_vote)
    d = band_depth(m, cache.union[roi], cache.core[roi], cache.dist_out_union[roi], cache.dist_in_core[roi])
    f5 = 1.0 / (1.0 + (2.0 * d / cache.perimeter) / (cache.t5 * cache.t5))
    return ObjectTerms(iou_tc, f5, area_guard(area, cache), sig, len(same) / len(cache.voters), area)


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
    """One completion -> {format_ok, legal, iou_tc, iou_min, f5, g_area, topo_any, topo_freq, sig_id, error,
    defects}."""
    return score_many([(text, truncated, names, caches, W, H)])[0][0]


def score_objects(text: str, truncated: bool, names: Sequence[str], caches: Mapping[str, ExpertCache],
                  W: int, H: int) -> tuple[dict, dict[str, np.ndarray] | None, dict[str, ObjectTerms]]:
    """score_completion's terms plus, when the completion parses, its mask and ObjectTerms per object (else None, {}):
    one pass for callers that need both (objdel.eval)."""
    return score_many([(text, truncated, names, caches, W, H)])[0]


def score_many(jobs: Sequence[tuple]) -> list[tuple[dict, dict[str, np.ndarray] | None, dict[str, ObjectTerms]]]:
    """score_objects of every job (tuples of its arguments), in order. The masks of all parsed objects come from one
    objdel.geometry.to_masks call: one CUDA checker call per score_many."""
    return [_score_one(*w) for w in _checked(jobs)]


def _checked(jobs: Sequence[tuple]) -> list[tuple]:
    """Parse every job and make the masks of all parsed objects in one to_masks call (one CUDA checker call), in the
    calling process. -> per job (names, caches, parser error, [Geometry per name] or None): what _score_checked
    needs, and all a pool job carries."""
    parsed = [parse(text, names, truncated) for text, truncated, names, _, _, _ in jobs]
    objs = [(p.objs[n], W, H) for p, (_, _, names, _, W, H) in zip(parsed, jobs, strict=True) if p.ok for n in names]
    geoms = iter(to_masks(objs))
    return [(names, caches, p.error, [next(geoms) for _ in names] if p.ok else None)
            for p, (_, _, names, caches, _, _) in zip(parsed, jobs, strict=True)]


def _score_one(names, caches, error, geoms) -> tuple[dict, dict[str, np.ndarray] | None, dict[str, ObjectTerms]]:
    """Terms of one checked job (CPU only: score_mask against the caches, never the checker)."""
    a = rater_answers(caches, names)
    if geoms is None:
        return _FAILED | {"error": error, "defects": {}}, None, {}
    legal, defects, masks, terms = True, {}, {}, {}
    for n, g in zip(names, geoms, strict=True):
        legal &= g.legal
        for d, c in g.defects.items():
            defects[d] = defects.get(d, 0) + c
        masks[n], terms[n] = g.mask, score_mask(g.mask, caches[n])
    ts = list(terms.values())
    sigs = tuple(t.sig for t in ts)
    joint = len(a.get(sigs, ())) / sum(map(len, a.values()))
    return ({"format_ok": 1.0, "legal": float(legal), "iou_tc": float(np.mean([t.iou_tc for t in ts])),
             "iou_min": float(min(t.iou_tc for t in ts)), "f5": float(np.mean([t.f5 for t in ts])),
             "g_area": float(np.mean([t.g_area for t in ts])), "topo_any": float(joint > 0), "topo_freq": joint,
             "sig_id": sig_id(sigs), "error": "", "defects": defects}, masks, terms)


def score_batch(jobs: Sequence[tuple], pool: Executor | None = None, chunk: int = REWARD_V4.group) -> dict:
    """score_completion over jobs (tuples of its arguments) in input order, serially or on an executor. The calling
    process parses every job and makes all masks in one CUDA checker call (_checked); the executor only scores the
    masks on the CPU (_score_checked), so its workers never touch CUDA (a process forked after CUDA started cannot
    use it, and one context per worker would share the GPU with the trainer: RL-5, CK-5).

    Masks go to the executor in chunks of `chunk` consecutive jobs (one prompt's group), so a process pool pickles
    each prompt's caches once per chunk (ExpertCache pickles as its rater masks, tau and T5; serve.CacheRef as its
    key). With a thread pool, call cv2.setNumThreads(0) first: OpenCV then runs each call sequentially (1 is a no-op
    on the macOS GCD backend); results do not depend on it. Returns float64 arrays for the numeric keys and lists for
    "error" and "defects"."""
    work = _checked(list(jobs))
    if pool is None:
        rows = _score_checked(work)
    else:
        chunks = [work[i:i + chunk] for i in range(0, len(work), chunk)]
        rows = [r for rs in pool.map(_score_checked, chunks) for r in rs]
    out = {k: np.array([r[k] for r in rows], np.float64) for k in NUMERIC}
    return out | {"error": [r["error"] for r in rows], "defects": [r["defects"] for r in rows]}


def _score_checked(work: Sequence[tuple]) -> list[dict]:
    """The pool job: terms of checked jobs (_checked), on the CPU; never imports or calls objdel._gpucheck."""
    return [_score_one(*w)[0] for w in work]
