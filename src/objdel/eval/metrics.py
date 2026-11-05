"""Evaluation metrics of K sampled completions per item against its raters (plan §4.7, judge R10).

Everything goes through the training path: objdel.reward.score_many (score_completion's parse, CUDA checker raster,
size-filtered signatures and per-object score_mask, with the masks and per-object terms kept) and advantages_v4 on
the item's ExpertCaches. Per item (`metrics`):

- per sample: format_ok, legal, IoU_tc, F5, g_area (legality-gated), topology-any (the joint size-filtered signature
  is some rater's), the v4 advantage in its group, `passed` = legal and topology-any and IoU_tc >= 0.5 for every object;
- per object and sample: IoU_tc, F5, g_area, IoU and Dice vs the strict-majority vote, D_max (best IoU to any single
  rater), area A, A / A_max, A / A_min, |M| / |E| (A over the raters' union, the egg of the egg-yolk band),
  size-filtered components and their excess over the most fragmented rater;
- per object and sample: NSD@tau (surface Dice) against the best rater, tau = NSD_TAU per structure;
- per object: soft Dice (MRNet protocol), the squared generalized energy distance (GED) and HM-IoU;
- per item: TV and JS (base 2) between the K-sample topology histogram and the rater histogram, minimum-mode
  retention, mode coverage, pass@k for k in 1, 2, 4, ..., <= min(K, 64), and per GRPO group of G = 8 samples (one
  group of K when K is not a multiple of 8) the duplicate share and whether some geometric term's std is above its
  sigma floor.

Definitions and choices:
- A sample's topology answer is its joint size-filtered signature (score_completion's sig_id); unparseable or
  illegal samples share one "fail" answer that no rater gives. The rater histogram is over the voting raters
  (objdel.reward.rater_answers: a rater whose whole mask is below the size filter abstains), as in the reward.
  Histograms are per image; summaries average the per-image TV / JS, never histograms pooled over images (lab
  t_dynamics REPORT).
- Minimum-mode retention = min over rater answers a of pi_v(a) / f(a), with pi_v the sample histogram restricted to
  rater answers and renormalised (lab t_dynamics grpo_toy.py); 1 iff pi_v = f, 0 when some rater answer is never
  sampled or no sample gives a rater answer.
- Coverage = share of distinct rater answers drawn by at least one sample (Multi-Answer RL's coverage).
- Pass criterion: the reward's own validity gate (objdel.reward.topology_valid: legal, a rater drew the topology,
  IoU_tc >= REWARD_V4.iou_gate for every object, score_completion's iou_min), so a RIGA answer with a good disc and a
  missed cup fails. The Phase 2 geometry gate (summarize's iou_tc_pass) is its geometric part, legal and IoU_tc >= 0.5
  for every object; topology is gated by the synthetic topology probe (objdel.eval.probes; plan D20), and
  topo_any_unanimous is reported only (val has 0 unanimous non-simple items).
- Duplicates (GRPO groups) are outputs the reward cannot tell apart: equal format and legality flags and equal
  rasters of every object (SHA-256 of the packed masks), so re-spelled texts of one mask (sub-pixel vertex nudges,
  whitespace) count as duplicates, as the reward's low-spread drop treats them (B2-6).
- Failed (unparseable or truncated) samples get 0 for every score, an empty mask for the probability map and GED,
  and areas of an empty mask (the audit only reads parsed samples). Illegal samples (parsed, but some ring defective)
  are failures too, by the reward's legality gate (objdel.reward.legal_gate, which zeroes every term of an illegal
  rollout in advantages_v4): 0 for IoU_tc, F5, g_area and the vote / D_max scores, and an empty mask for the probability
  map and GED, so they never count toward the geometry gates (REVIEW #2). Their areas and components stay those of
  their raster, which the audit reads for every parsed sample.
- Soft Dice (Ji et al. 2021, MRNet code): the soft ground truth is the mean of the rater masks (uniform expertness
  1/k), the prediction the mean of the K sample masks; both are binarized with >= t for t in 0.1, 0.3, 0.5, 0.7, 0.9,
  Dice is computed per threshold (two empty masks: 1) and averaged. QUBIQ's BDice is the same with 9 thresholds.
- GED (Kohl et al. 2018, App. B): 2 E d(S, Y) - E d(S, S') - E d(Y, Y') over all sample / rater pairs (V-statistic),
  d = 1 - IoU, d(empty, empty) = 0.
- pass@k: the unbiased estimator 1 - C(n - c, k) / C(n, k), in Chen et al.'s numerically stable product form.
- HM-IoU (HPU-Net, CAR; "HM-IoU16" at K = 16): Hungarian matching of the K samples to the k raters, each rater
  duplicated so that it is matched floor(K / k) or ceil(K / k) times (exactly K / k when k divides K, the published
  protocol), maximising the total IoU; the mean IoU of the matched pairs (two empty masks: 1). Solved as one
  rectangular assignment (scipy linear_sum_assignment) on floor(K / k) mandatory copies of each rater, forced by a
  bonus, plus one optional copy per rater (impl/M4.md: the ceil-only tiling lets a rater go unmatched).
- NSD@tau (Nikolov et al. 2018; Metrics Reloaded's default boundary metric), as google-deepmind/surface-distance
  (compute_surface_distances + compute_surface_dice_at_tolerance, spacing 1 px): surfels on the corner grid of 2 x 2
  neighbourhood codes, weighted by their contour length; (|S_P within tau of S_G| + |S_G within tau of S_P|) /
  (|S_P| + |S_G|), distances by exact EDT; two empty masks: 1 (the reference gives 0 / 0), one empty: 0. Per sample
  the maximum over the raters (like D_max), so a rater's own target scores 1. tau at frame resolution (D25, A4 S2):
  2 px LIDC, 9 px disc, 16 px cup. 4-connected boundary pixels (A4's prototype) differ from the reference by up to
  0.14 on LIDC rater pairs (impl/M4.md).
- Token statistics: per sample completion token counts when the caller passes them (train.sample / sample_vllm);
  summarize reports the truncated share, mean and p95 tokens.

References:
- W. Ji et al., "Learning calibrated medical image segmentation via multi-rater agreement modeling", CVPR 2021
  (soft Dice at thresholds 0.1-0.9, RIGA BinRushed + MESSIDOR train / Magrabia test):
  https://mlanthology.org/cvpr/2021/ji2021cvpr-learning ; evaluation code (thresholds, >=, uniform 1/6 rater
  weights): https://github.com/jiwei0921/MRNet/blob/main/MRNet_Code/trainer/trainer_cvpr.py
- QUBIQ challenge, binarized Dice at 0.1, ..., 0.9 against the mean of the raters: https://arxiv.org/abs/2405.18435 ,
  https://qubiq.grand-challenge.org/participation/
- J. Orlando et al., REFUGE challenge (single-reference Dice of disc and cup): https://arxiv.org/abs/1910.03667
- S. Kohl et al., "A probabilistic U-Net for segmentation of ambiguous images", NeurIPS 2018 (GED, App. B):
  https://arxiv.org/abs/1806.05034
- M. Chen et al., "Evaluating large language models trained on code", 2021 (unbiased pass@k, Fig. 3):
  https://arxiv.org/abs/2107.03374
- Multi-Answer RL (coverage of a set of valid answers): https://arxiv.org/abs/2603.24844
- HM-IoU: C. Kohl et al., "A hierarchical probabilistic U-Net", 2019, https://arxiv.org/abs/1905.13077 ; E. Kassapis
  et al., CAR (HM-IoU16), https://arxiv.org/abs/2006.13144 ; scipy.optimize.linear_sum_assignment:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html
- NSD: S. Nikolov et al., https://arxiv.org/abs/1809.04430 ; reference code:
  https://github.com/google-deepmind/surface-distance (metrics.py, lookup_tables.py) ; Metrics Reloaded:
  https://arxiv.org/abs/2206.01653
- Per-image TV / JS and minimum-mode retention: experiments/phase1/lab/t_dynamics/REPORT.md, grpo_toy.py
- Duplicate hash on packed rasters (numpy.packbits):
  https://numpy.org/doc/stable/reference/generated/numpy.packbits.html ; identical rewards give a zero-signal group
  (DAPO dynamic sampling): https://arxiv.org/abs/2503.14476
"""
import hashlib
from collections import Counter
from collections.abc import Mapping, Sequence

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

from objdel.reward import (
    NUMERIC,
    REWARD_V4,
    TERMS,
    ExpertCache,
    advantages_v4,
    legal_gate,
    rater_answers,
    score_many,
    score_mask,
    sig_id,
    topology_valid,
)

THRESHOLDS = (0.1, 0.3, 0.5, 0.7, 0.9)  # MRNet soft Dice
PASS_KS = (1, 2, 4, 8, 16, 32, 64)
FAIL = -1.0  # topology answer of unparseable or illegal samples (sig_ids are CRC32 values >= 0)
_ZEROED = ("iou_tc", "f5", "g_area", "iou_vote", "dice_vote", "d_max")  # scores of failed and illegal samples
NSD_TAU = {"nodule": 2.0, "disc": 9.0, "cup": 16.0}  # px at frame resolution (D25)
_D = 0.5 * np.sqrt(2.0)  # surfel length per 2 x 2 neighbourhood code (surface-distance lookup table, 1 px spacing)
_SURFEL = np.array([0, _D, _D, 1, _D, 1, 2 * _D, _D, _D, 2 * _D, 1, _D, 1, _D, _D, 0])
_GEOMETRIC = [(t, f) for t, f in zip(TERMS, REWARD_V4.floors, strict=True) if t != "legal"]  # (term, sigma floor)


def pass_at_k(n: int, c: int, k: int) -> float:
    """Unbiased pass@k from n samples with c passes (Chen et al. 2021)."""
    if n - c < k:
        return 1.0
    return 1.0 - float(np.prod(1.0 - k / np.arange(n - c + 1, n + 1)))


def _unit(x) -> float:
    """x clamped to [0, 1]: float sums of shares can land 2.2e-16 outside (S1-15: NSD -2.2e-16 at tau = 0; TV, JS and
    retention 1 + 2.2e-16)."""
    return min(1.0, max(0.0, float(x)))


def tv_js(p: Mapping, q: Mapping) -> tuple[float, float]:
    """Total variation and Jensen-Shannon divergence (base 2, in [0, 1]) of two histograms {answer: probability}."""
    keys = sorted(set(p) | set(q))
    a, b = np.array([p.get(x, 0.0) for x in keys]), np.array([q.get(x, 0.0) for x in keys])
    m = (a + b) / 2

    def kl(x):
        nz = x > 0
        return float(np.sum(x[nz] * np.log2(x[nz] / m[nz])))

    return _unit(0.5 * float(np.abs(a - b).sum())), _unit(0.5 * kl(a) + 0.5 * kl(b))


def min_mode_retention(p: Mapping, f: Mapping) -> float:
    """min over rater answers a of p_v(a) / f(a), p_v = p restricted to f's support, renormalised; 0 if empty."""
    valid = sum(p.get(a, 0.0) for a in f)
    return _unit(min(p.get(a, 0.0) / valid / fa for a, fa in f.items())) if valid > 0 else 0.0


def _dice(a: np.ndarray, b: np.ndarray) -> float:
    s = np.count_nonzero(a) + np.count_nonzero(b)
    return 2.0 * np.count_nonzero(a & b) / s if s else 1.0


def soft_dice(samples: np.ndarray, raters: np.ndarray) -> float:
    """Mean over THRESHOLDS t of Dice(mean(samples) >= t, mean(raters) >= t); masks (n, H, W) and (k, H, W)."""
    p, r = samples.mean(0), raters.mean(0)
    return float(np.mean([_dice(p >= t, r >= t) for t in THRESHOLDS]))


def _joint_box(masks: np.ndarray) -> tuple[slice, slice] | None:
    """Bounding box (row slice, column slice) of the foreground of all masks (n, H, W); None if there is none."""
    rows, cols = np.flatnonzero(masks.any((0, 2))), np.flatnonzero(masks.any((0, 1)))
    return (slice(rows[0], rows[-1] + 1), slice(cols[0], cols[-1] + 1)) if rows.size else None


def _pairwise_iou(masks: np.ndarray) -> np.ndarray:
    """IoU of every pair of masks (n, H, W); two empty masks: 1. One matrix product over the distinct masks on their
    joint bounding box (duplicates, e.g. a deterministic model's K samples, are computed once)."""
    box = _joint_box(masks)
    if box is None:
        return np.ones((len(masks), len(masks)))
    crop = masks[(slice(None), *box)].reshape(len(masks), -1)
    _, first, inv = np.unique(np.packbits(crop, axis=1), axis=0, return_index=True, return_inverse=True)
    x = crop[first].astype(np.float64)
    inter = x @ x.T  # exact: integer counts far below 2**53
    area = np.diag(inter)
    union = area[:, None] + area[None, :] - inter
    iou = np.where(union > 0, inter / np.where(union > 0, union, 1.0), 1.0)
    inv = inv.reshape(-1)
    return iou[inv][:, inv]


def hm_iou(iou: np.ndarray) -> float:
    """Hungarian-matched IoU from the (K samples, k raters) IoU matrix (module docstring)."""
    K, k = iou.shape
    base, extra = divmod(K, k)
    cols = np.concatenate([np.tile(np.arange(k), base), np.arange(k) if extra else []]).astype(int)
    bonus = np.where(np.arange(len(cols)) < base * k, 2.0 * K + 2, 0.0)  # every mandatory copy is matched
    r, c = linear_sum_assignment(iou[:, cols] + bonus, maximize=True)
    return float(iou[r, cols[c]].mean())


def _surface(m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(surfel length on the (H + 1, W + 1) corner grid, 0 off the surface; distance of every corner to the surface,
    inf for an empty mask)."""
    p = np.pad(m, 1).astype(np.uint8)
    code = 8 * p[:-1, :-1] + 4 * p[:-1, 1:] + 2 * p[1:, :-1] + p[1:, 1:]
    border = (code != 0) & (code != 15)
    dist = ndimage.distance_transform_edt(~border) if border.any() else np.full(border.shape, np.inf)
    return _SURFEL[code], dist


def _nsd(a: tuple, b: tuple, tau: float) -> float:
    (la, da), (lb, db) = a, b
    total = la.sum() + lb.sum()
    return _unit(1.0 - (la[db > tau].sum() + lb[da > tau].sum()) / total) if total else 1.0  # exactly 1 if all in


def nsd(a: np.ndarray, b: np.ndarray, tau: float) -> float:
    """Normalized surface Dice of two masks at tolerance tau px (module docstring)."""
    return _nsd(_surface(a), _surface(b), tau)


def _nsd_best(drawn: np.ndarray, raters: np.ndarray, tau: float) -> np.ndarray:
    """Per sample, the maximum NSD over the raters, on the joint bounding box (distances are unchanged by the crop);
    each distinct mask once."""
    both = np.concatenate([drawn, raters])
    box = _joint_box(both)
    if box is None:
        return np.ones(len(drawn))
    crop = both[(slice(None), *box)]
    rs = [_surface(r) for r in crop[len(drawn):]]
    seen = {}
    for m in crop[:len(drawn)]:
        key = np.packbits(m).tobytes()
        if key not in seen:
            s = _surface(m)
            seen[key] = max(_nsd(s, r, tau) for r in rs)
    return np.array([seen[np.packbits(m).tobytes()] for m in crop[:len(drawn)]])


def _ged(iou: np.ndarray, n: int) -> float:
    """GED from the pairwise IoU of n samples followed by the raters."""
    d = 1.0 - iou
    return max(0.0, float(2 * d[:n, n:].mean() - d[:n, :n].mean() - d[n:, n:].mean()))  # BK2-9: -1.1e-16 at 0


def ged(samples: np.ndarray, raters: np.ndarray) -> float:
    """Squared generalized energy distance with d = 1 - IoU (Kohl et al. 2018), all pairs."""
    return _ged(_pairwise_iou(np.concatenate([samples, raters])), len(samples))


def _samples(keys: Sequence[tuple[str, bool]], names, caches, W: int, H: int) -> list[dict]:
    """Training-path terms of distinct completions (text, truncated), from one score_many call: per completion its
    terms, per-object masks and terms (empty masks for a failed sample), and its duplicate hash."""
    out = []
    for terms, masks, objects in score_many([(t, tr, names, caches, W, H) for t, tr in keys]):
        if masks is None:
            masks = {n: np.zeros((H, W), bool) for n in names}
            objects = {n: score_mask(masks[n], caches[n]) for n in names}
        h = hashlib.sha256(bytes([int(terms["format_ok"]), int(terms["legal"])]))
        for n in names:
            h.update(np.packbits(masks[n]).tobytes())
        out.append({"terms": terms, "masks": masks, "objects": objects, "hash": h.hexdigest()})
    return out


def _native(x) -> list:
    return [float(v) for v in np.asarray(x, np.float64)]


def metrics(samples: Sequence[str], raters: Mapping[str, np.ndarray], W: int, H: int,
            names: Sequence[str] | None = None, truncated: Sequence[bool] | None = None,
            caches: Mapping[str, ExpertCache] | None = None, tokens: Sequence[int] | None = None,
            nsd_tau: Mapping[str, float] = NSD_TAU) -> dict:
    """Metrics of one item: K completions (texts) against the rater masks {name: bool (k, H, W)} on a W x H frame;
    tokens: completion token counts (optional); nsd_tau: NSD tolerance per structure, px. The masks come from one
    CUDA checker call per item (objdel.reward.score_many).

    Returns JSON-native values: per-sample lists under "samples" and "objects"[name], item scalars otherwise
    (module docstring). Raises ValueError for fewer than 2 raters or samples, or mismatched lengths."""
    names = tuple(names or raters)
    truncated = [False] * len(samples) if truncated is None else [bool(t) for t in truncated]
    if len(samples) < 2 or len(truncated) != len(samples):
        raise ValueError("need >= 2 samples and one truncated flag per sample")
    if tokens is not None and len(tokens) != len(samples):
        raise ValueError("need one token count per sample")
    caches = caches or {n: ExpertCache.build(raters[n]) for n in names}  # ValueError for < 2 raters
    K, k = len(samples), caches[names[0]].k
    keys = [(str(t), tr) for t, tr in zip(samples, truncated, strict=True)]
    uniq = list(dict.fromkeys(keys))  # each distinct completion once
    distinct = dict(zip(uniq, _samples(uniq, names, caches, W, H), strict=True))
    rows = [distinct[x] for x in keys]
    T = {t: np.array([r["terms"][t] for r in rows]) for t in NUMERIC}  # failed samples: every term 0 (reward rule)
    legal = legal_gate(T)  # the reward's legality gate; legal implies parsed
    gated = {t: T[t] * legal for t, _ in _GEOMETRIC}  # illegal: 0, as advantages_v4 sees the terms
    passed = topology_valid(T)
    objects = {}
    for n in names:
        c, ot = caches[n], [r["objects"][n] for r in rows]
        masks = np.stack([r["masks"][n] for r in rows])
        drawn = masks & legal[:, None, None]  # what is scored: illegal rasters count as empty
        iou = _pairwise_iou(np.concatenate([drawn, c.raters, c.vote[None]]))
        area, comps = masks.sum((1, 2)), np.array([len(t.sig) for t in ot])
        a_min, a_max = min(c.rater_areas), max(c.rater_areas)
        o = {"iou_tc": [t.iou_tc for t in ot], "f5": [t.f5 for t in ot], "g_area": [t.g_area for t in ot],
             "iou_vote": iou[:K, -1], "dice_vote": 2 * iou[:K, -1] / (1 + iou[:K, -1]),  # Dice = 2 IoU / (1 + IoU)
             "d_max": iou[:K, K:K + k].max(1),
             "area": area, "a_amax": area / max(a_max, 1), "a_amin": area / a_min if a_min else np.ones(K),
             "m_e": area / max(int(np.count_nonzero(c.union)), 1),
             "components": comps, "extra_components": comps - max(len(s) for s in c.rater_sigs),
             "nsd": _nsd_best(drawn, c.raters, nsd_tau[n])}
        o = {key: np.asarray(v, np.float64) for key, v in o.items()}
        for key in _ZEROED:
            o[key][~legal] = 0.0
        objects[n] = {**{key: _native(v) for key, v in o.items()},
                      "soft_dice": soft_dice(drawn, c.raters), "ged": _ged(iou[:-1, :-1], K),
                      "hm_iou": hm_iou(iou[:K, K:K + k])}

    # topology answers: the joint size-filtered signature, as score_completion's sig_id; failures share FAIL
    answers = rater_answers(caches, names)
    votes = sum(map(len, answers.values()))
    f = {sig_id(a): len(rs) / votes for a, rs in answers.items()}
    p = {a: c / K for a, c in Counter(np.where(legal, T["sig_id"], FAIL).tolist()).items()}
    tv, js = tv_js(p, f)
    group = REWARD_V4.group if K % REWARD_V4.group == 0 else K
    starts = range(0, K, group)
    adv = np.concatenate([advantages_v4({t: v[i:i + group] for t, v in T.items()}, group) for i in starts])
    hashes = [r["hash"] for r in rows]
    return {
        "k": k, "K": K, "group": group, "modes": len(f),
        "samples": {"format_ok": _native(T["format_ok"]), "legal": _native(T["legal"]),
                    "iou_tc": _native(gated["iou_tc"]), "f5": _native(gated["f5"]), "g_area": _native(gated["g_area"]),
                    "topo_any": _native(T["topo_any"]),
                    "passed": _native(passed), "advantage": _native(adv), "error": [r["terms"]["error"] for r in rows],
                    "truncated": _native(truncated),
                    **({} if tokens is None else {"tokens": [int(t) for t in tokens]})},
        "objects": objects,
        "tv": tv, "js": js, "min_mode": min_mode_retention(p, f), "coverage": len(f.keys() & p.keys()) / len(f),
        "pass_at_k": {str(x): pass_at_k(K, int(passed.sum()), x) for x in PASS_KS if x <= K},
        "groups": {"dup": [1.0 - len(set(hashes[i:i + group])) / group for i in starts],
                   "signal": [any(float(gated[t][i:i + group].std()) > fl for t, fl in _GEOMETRIC)
                              for i in starts]},  # std of the legality-gated term, as advantages_v4 computes it
    }


def _mean(xs) -> float | None:
    return float(np.mean(xs)) if len(xs) else None


def summarize(items: Sequence[dict]) -> dict:
    """Per dataset (item["dataset"]): rates over all samples, per-structure medians and means over all samples
    (failed and illegal samples count as 0), per-image metrics averaged over items, and the same on disagreement items
    (> 1 rater topology answer) for the topology metrics; the truncated share and, when every item has token counts,
    their mean and p95 ("tokens": None otherwise)."""
    out = {}
    for ds in sorted({it["dataset"] for it in items}):
        its = [it for it in items if it["dataset"] == ds]
        S = {key: np.concatenate([it["samples"][key] for it in its])
             for key in ("format_ok", "legal", "topo_any", "passed", "truncated")}
        tokens = [it["samples"].get("tokens") for it in its]
        tokens = np.concatenate(tokens) if all(t is not None for t in tokens) else None
        unanimous = [x for it in its if it["modes"] == 1 for x in it["samples"]["topo_any"]]
        dis = [it for it in its if it["modes"] > 1]
        names = list(its[0]["objects"])
        objs = {n: {key: np.concatenate([it["objects"][n][key] for it in its]) for key in (*_ZEROED, "nsd")}
                for n in names}
        k_min = min(it["K"] for it in its)
        groups = {it["group"] for it in its}
        out[ds] = {
            "items": len(its), "samples": int(S["format_ok"].size),
            "format": float(S["format_ok"].mean()), "legal": float(S["legal"].mean()),
            "errors": dict(Counter(e for it in its for e in it["samples"]["error"] if e)),
            "topo_any": float(S["topo_any"].mean()), "topo_any_unanimous": _mean(unanimous),
            "pass_rate": float(S["passed"].mean()), "truncated": float(S["truncated"].mean()),
            "tokens": None if tokens is None else {"mean": float(tokens.mean()),
                                                   "p95": float(np.percentile(tokens, 95))},
            "iou_tc_pass": float(np.mean(np.min([objs[n]["iou_tc"] for n in names], 0) >= REWARD_V4.iou_gate)),
            "pass_at_k": {str(x): _mean([it["pass_at_k"][str(x)] for it in its]) for x in PASS_KS if x <= k_min},
            **{key: _mean([it[key] for it in its]) for key in ("tv", "js", "coverage")},
            "disagreement": {"items": len(dis), **{key: _mean([it[key] for it in dis])
                                                   for key in ("tv", "js", "min_mode", "coverage")}},
            "group": groups.pop() if len(groups) == 1 else None,
            "dup_rate": _mean([d for it in its for d in it["groups"]["dup"]]),
            "signal_share": _mean([s for it in its for s in it["groups"]["signal"]]),
            "objects": {n: {"iou_tc_median": float(np.median(objs[n]["iou_tc"])),
                            **{f"{key}_mean": float(objs[n][key].mean()) for key in _ZEROED},
                            "f5_median": float(np.median(objs[n]["f5"])),
                            "soft_dice": _mean([it["objects"][n]["soft_dice"] for it in its]),
                            "ged": _mean([it["objects"][n]["ged"] for it in its]),
                            "hm_iou": _mean([it["objects"][n]["hm_iou"] for it in its]),
                            "nsd_mean": float(objs[n]["nsd"].mean()), "nsd_median": float(np.median(objs[n]["nsd"]))}
                        for n in names},
        }
    return {"datasets": out}
