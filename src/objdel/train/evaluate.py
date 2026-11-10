"""Evaluation without a trainer (Implementation.md E4): rows of an eval view, one item's checker metrics, and the
per-item and summary statistics of final_eval (policy/evaluate.py after the E5 move).

- eval_rows(view, items, rl_only): the prompt rows of a view's lidc/val records (a seeded subset with `items`).
- score_items(rows, done, root): per-item metrics (objdel.eval.metrics on the preloaded caches, in this process: one
  CUDA context, D7) and the signal groups.
- per_item(result) and summarize(items, groups): final_eval's K-sample means, answer choices (A4) and the schema-2
  summary with item, nodule and patient bootstrap intervals (final_eval's module docstring names every key).
"""

from pathlib import Path

import numpy as np

from objdel.eval import metrics, select, signal
from objdel.train import serve

SUBSET_SEED = 0  # view_rows' `items`: one fixed subset for every branch and rung

# name -> (source, higher_is_better): "samples" and "objects" lists are averaged over the K answers, "item" scalars
# (one per object and item) are copied
METRICS = {"legal": ("samples", True), "iou_tc": ("samples", True), "f5": ("samples", True),
           "precision_vote": ("objects", True), "recall_vote": ("objects", True), "dice_vote": ("objects", True),
           "nsd": ("objects", True), "boundary_iou": ("objects", True), "boundary_iou_vote": ("objects", True),
           "soft_dice": ("item", True), "hm_iou": ("item", True), "ged": ("item", False)}
SELECTION = ("sel_mean_of_k", "sel_legal_mean", "sel_medoid", "sel_medoid_legal", "sel_oracle", "sel_vote")


def _score(job) -> dict:
    """One item's metrics on the preloaded caches."""
    item_id, texts, truncated, tokens, W, H, names, root = job
    m = metrics(texts, {}, W, H, names, truncated, caches=serve.CacheRef(item_id, root), tokens=tokens)
    return {"item_id": item_id, "dataset": item_id.split("/")[0], **m}


def view_rows(view: Path, items: int | None = None, rl_only: bool = True) -> list[dict]:
    """Prompt rows of an eval view's lidc/val records (images checked against their records): the scored items (rl,
    with caches) or, with rl_only=False, every record (probe records are rl: false); `items` = a seeded subset."""
    from objdel.train import dataset as D

    return D.prompts(view, "lidc", "val", rl_only, subset=items and (items, SUBSET_SEED))


def score_items(rows: list[dict], done: dict, root: Path) -> tuple[list[dict], list[dict]]:
    """_score in this process for every row (done: {key: completions row with samples, truncated and tokens}) ->
    (per-item metrics, signal groups); two CUDA checker calls per item (score_many, eval_sets), no forked worker."""
    jobs = [(r["key"], *(done[r["key"]][k] for k in ("samples", "truncated", "tokens")), r["W"], r["H"], r["names"],
             str(root)) for r in rows]  # fmt: skip
    results = [_score(j) for j in jobs]
    return results, [signal.group(r["item_id"], done[r["item_id"]]["samples"], r) for r in results]


def per_item(result: dict, obj: str = "nodule") -> dict:
    """K-sample means of one item's scores (metrics' output; METRICS), `miss` (the K-mean IoU_tc below signal.MISS),
    and the answer choices (A4): the mean of the K IoU_tc,
    the mean over legal answers (0 if none), the medoid's (with its legality), the best one (it reads the labels: an
    upper bound, never deployed) and the pixel vote's (analysis only: a mask is not an answer; absent when the
    result has no `pixel_vote_iou_tc`)."""
    o = result["objects"][obj]
    src = {"samples": result["samples"], "objects": o, "item": o}
    iou, legal = np.array(result["samples"]["iou_tc"], float), np.array(result["samples"]["legal"]) > 0
    i = select.medoid(select.sample_iou(result["objects"]), legal)
    vote = o.get("pixel_vote_iou_tc")
    means = {m: float(np.mean(src[s][m])) for m, (s, _) in METRICS.items()}
    return {"item_id": result["item_id"], **means, "miss": float(means["iou_tc"] < signal.MISS),
            "sel_mean_of_k": float(iou.mean()), "sel_legal_mean": float(iou[legal].mean()) if legal.any() else 0.0,
            "sel_medoid": float(iou[i]), "sel_medoid_legal": float(legal[i]), "sel_oracle": float(iou.max()),
            **({} if vote is None else {"sel_vote": float(vote)})}


def summarize(items: list[dict], groups: dict[str, str] | None = None, draws: int = 2000, seed: int = 0) -> dict:
    """Schema 2: per METRICS key the items hold, the mean and its item, nodule and patient (with `groups`,
    {item_id: group_id}) 95 % bootstrap intervals; the median of iou_tc; the miss rate. One draw set per unit."""
    ids = [it["item_id"] for it in items]
    if not ids:
        raise ValueError("summarize: no items")
    if len(set(ids)) != len(ids):
        raise ValueError("summarize: duplicate item_id")
    def unit_draws(sub, key=signal.cluster):
        return signal.point_or_draws(sub, draws, seed, key)

    D = {"ci95_nodule": unit_draws(ids), "ci95_patient": None if groups is None else unit_draws(ids,
        groups.__getitem__)}
    idx = np.random.default_rng(seed).integers(0, len(items), (draws, len(items)))

    def cis(v: np.ndarray) -> dict:
        lo, hi = np.percentile(v[idx].mean(1), [2.5, 97.5])
        return {"mean": float(v.mean()), **{k: None if d is None else signal.ci(v, d)[1:] for k, d in D.items()},
                "ci95_item": [float(lo), float(hi)]}

    out = {"schema": 2, "items": len(items), "nodules": len(set(map(signal.cluster, ids))),
           "patients": None if groups is None else len({groups[i] for i in ids}), "ci_level": 0.95, "boot": draws,
           "gate_ci": "ci95_patient"}
    for m, higher in [(m, h) for m, (_, h) in METRICS.items()] + [(m, True) for m in SELECTION]:
        held = [m in it for it in items]
        if not any(held):
            continue
        if not all(held):
            raise ValueError(f"summarize: item {ids[held.index(False)]} has no {m}")
        v = np.array([it[m] for it in items], float)
        bad = ~np.isfinite(v)
        if bad.any():
            raise ValueError(f"summarize: item {ids[bad.argmax()]} has {m} = {v[bad.argmax()]}")
        out[m] = cis(v) | {"higher_is_better": higher}
    iou = np.array([it["iou_tc"] for it in items], float)
    out["iou_tc"]["median"] = float(np.median(iou))
    out["miss_rate"] = cis((iou < signal.MISS).astype(float))
    found = iou >= signal.MISS
    fids = [i for i, f in zip(ids, found, strict=True) if f]
    out["iou_found"] = None if not fids else {
        "mean": float(iou[found].mean()), "items": len(fids),
        "ci95_patient": None if groups is None else signal.ci(iou[found], unit_draws(fids, groups.__getitem__))[1:],
        "ci95_nodule": signal.ci(iou[found], unit_draws(fids))[1:]}

    def paired(arm: str) -> dict | None:  # the per-item difference arm - sel_mean_of_k under the same draws
        if arm not in out:
            return None
        d = np.array([it[arm] - it["sel_mean_of_k"] for it in items])
        return {"mean": float(d.mean()), "ci95_patient": None if D["ci95_patient"] is None
                else signal.ci(d, D["ci95_patient"])[1:], "ci95_nodule": signal.ci(d, D["ci95_nodule"])[1:]}

    out["paired"] = {"medoid_minus_mean": paired("sel_medoid"), "vote_minus_mean": paired("sel_vote")}
    return out
