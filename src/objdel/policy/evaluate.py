"""Scoring of sampled answers per item, shared by checkpoint selection, the G-probe and the final exam (E4 adds
`evaluate`): eval-view prompt rows and the in-process metrics of each item's samples on the preloaded caches."""

from pathlib import Path

from objdel.eval import metrics
from objdel.train import serve

SUBSET_SEED = 0  # --items: one fixed subset for every branch and rung


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
    from objdel.eval import signal

    jobs = [(r["key"], *(done[r["key"]][k] for k in ("samples", "truncated", "tokens")), r["W"], r["H"], r["names"],
             str(root)) for r in rows]  # fmt: skip
    results = [_score(j) for j in jobs]
    return results, [signal.group(r["item_id"], done[r["item_id"]]["samples"], r) for r in results]
