"""Fair human ceiling (Implementation.md A2): each doctor and the model are scored against the same k - 1 doctors.

Per item with k doctors (`cache.raters`) and voting doctors V (`cache.voters`; a doctor whose size-filtered signature
is empty has no vote, ExpertCache docstring):
- for each r in V the leave-one-out cache is built from the other k - 1 doctors (a non-voter stays in every stack).
  At k = 2 that is one doctor, which ExpertCache refuses, so the stack is [other, other]: IoU_tc against it is the
  plain IoU with the other doctor (its band is the mask itself; tau falls to its floor, so F5 there is meaningless);
- doctor = mean over r of IoU_tc(doctor r; LOO r);
- model = mean over r and over the K answers of IoU_tc(answer; LOO r), illegal or failed answers 0 (the reward's
  legality gate); model_all_k = the same against all k doctors (final_eval's iou_tc); diff = model - doctor;
- legacy_dmax = per r in V, the max plain IoU of doctor r with any other doctor (the old ceiling, 0.723 val / 0.732
  test), pooled per item (`legacy_loo_dmax_item`) or over all (item, r) pairs (`legacy_loo_dmax_pooled`).
Both sides see the same references, so the gap from the number of references cancels (notes/research/math-stats.md
§3). `model_all_k >= model` is not guaranteed item by item: IoU_tc falls back to the vote when no doctor drew the
answer's signature.

summarize: strata all, k>=3, k=2 (k = all doctors) and, with answers, found / miss (the K-mean rule on model_all_k,
signal.MISS). Each metric is [mean, lo, hi] under patient draws (the gate CI, `groups` = {item_id: group_id}); the
nodule CI is the second line under "nodule". `diff` is paired: the CI of the per-item differences under the same draws.
"""

from collections.abc import Mapping, Sequence
from typing import NamedTuple

import numpy as np

from objdel.reward import score
from objdel.reward.cache import ExpertCache

from . import signal

METRICS = ("doctor", "model", "model_all_k", "diff", "legacy_loo_dmax_item")
MODEL = ("model", "model_all_k", "diff")


class LeaveOut(NamedTuple):
    index: int  # the left-out voting doctor
    cache: ExpertCache  # the other k - 1 doctors ([other, other] at k = 2)
    single: bool  # k = 2: one reference doctor


def loo_caches(cache: ExpertCache) -> list[LeaveOut]:
    """One entry per voting doctor r, in increasing r: the cache of the other doctors (module docstring). k = 1:
    ValueError from ExpertCache.build (an empty stack)."""
    out = []
    for r in cache.voters:
        rest = np.delete(cache.raters, r, 0)
        single = len(rest) == 1
        out.append(LeaveOut(r, ExpertCache.build(np.concatenate([rest, rest]) if single else rest), single))
    return out


def _iou(a: np.ndarray, b: np.ndarray) -> float:
    u = np.count_nonzero(a | b)
    return np.count_nonzero(a & b) / u if u else 1.0


def item(cache: ExpertCache, texts: Sequence[str] | None, truncated: Sequence[bool] | None, W: int, H: int,
         name: str = "nodule") -> dict:
    """One item's ceiling row (module docstring); texts = None gives the doctor side only (model keys null)."""
    loos = loo_caches(cache)
    R = cache.raters
    dmax = [max(_iou(R[lo.index], R[j]) for j in range(cache.k) if j != lo.index) for lo in loos]
    row = {"k": cache.k, "n_voters": len(loos), "single": cache.k == 2,
           "doctor": float(np.mean([score.score_mask(R[lo.index], lo.cache).iou_tc for lo in loos])),
           "legacy_dmax": dmax, "legacy_loo_dmax_item": float(np.mean(dmax))}
    if texts is None:
        return row | dict.fromkeys(MODEL)
    truncated = [False] * len(texts) if truncated is None else list(truncated)
    if len(texts) != len(truncated):
        raise ValueError(f"{len(texts)} texts but {len(truncated)} truncated flags")
    refs = [lo.cache for lo in loos] + [cache]  # the LOO caches, then all k doctors (model_all_k)
    jobs = [(t, tr, (name,), {name: c}, W, H) for c in refs for t, tr in zip(texts, truncated, strict=True)]
    terms = [t for t, _ in score.score_many(jobs)]
    iou = np.array([t["iou_tc"] if t["legal"] >= 0.5 else 0.0 for t in terms])  # not *: NaN * 0 is NaN
    iou = iou.reshape(len(refs), len(texts))
    model = float(iou[:-1].mean())
    return row | {"model": model, "model_all_k": float(iou[-1].mean()), "diff": model - row["doctor"]}


def _stratum(rows: Sequence[dict], groups: Mapping[str, str] | None, boot: int, seed: int) -> dict:
    """n, patients, items with a non-voting doctor, each metric's [mean, lo, hi] (patient draws; [mean, None, None]
    without groups), the pooled legacy value and the nodule-draw line; null metrics when empty or without answers."""
    ids = [r["item_id"] for r in rows]
    out = {"n": len(rows), "n_patients": None if groups is None else len({groups[i] for i in ids}),
           "n_voters_lt_raters": sum(r["n_voters"] < r["k"] for r in rows)}
    if not rows:
        return out | dict.fromkeys(METRICS) | {"legacy_loo_dmax_pooled": None, "nodule": dict.fromkeys(METRICS)}

    def line(D):
        res = {}
        for m in METRICS:
            x = None if rows[0][m] is None else np.array([r[m] for r in rows], float)
            res[m] = None if x is None else [float(x.mean()), None, None] if D is None else signal.ci(x, D)
        return res

    patient = None if groups is None else signal.draws(ids, boot, seed, key=groups.__getitem__)
    pooled = float(np.mean([v for r in rows for v in r["legacy_dmax"]]))
    return out | line(patient) | {"legacy_loo_dmax_pooled": pooled, "nodule": line(signal.draws(ids, boot, seed))}


def summarize(rows: Sequence[dict], groups: Mapping[str, str] | None = None, boot: int = 2000, seed: int = 0) -> dict:
    """Strata all, k>=3, k=2 and, with answers, found / miss (module docstring); rows are `item` rows with item_id.
    ValueError when answers are given without groups, or only some rows have answers."""
    answered = [r["model"] is not None for r in rows]
    if any(answered) and not all(answered):
        raise ValueError(f"summarize: item {rows[answered.index(False)]['item_id']} has no answers, others have")
    if any(answered) and groups is None:
        raise ValueError("summarize: the answers' CIs need groups ({item_id: group_id}, patient draws)")
    strata = {"all": rows, "k>=3": [r for r in rows if r["k"] >= 3], "k=2": [r for r in rows if r["k"] == 2]}
    if any(answered):
        strata |= {"found": [r for r in rows if r["model_all_k"] >= signal.MISS],
                   "miss": [r for r in rows if r["model_all_k"] < signal.MISS]}
    out = {s: _stratum(v, groups, boot, seed) for s, v in strata.items()}
    out["k=2"]["single_reference"] = True  # one reference doctor: IoU_tc is the plain IoU, F5 is meaningless
    return out
