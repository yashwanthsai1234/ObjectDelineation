"""Final evaluation of a tuned model on a whole split: K samples per item, scored by the checker.

    python -I -m experiments.phase2.final_eval --ckpt DIR --view data/eval_views/lidc-val-shifted --split val \
        --train data/prepared --out runs/eval/val [--k 8] [--seed 11] [--per-call 64]
    python -I -m experiments.phase2.final_eval --resummarize runs/r5/eval/eval_test --out runs/r5/eval_a1/eval_test \
        [--store data/store]

Every RL-eligible item of the view's <split> records (the G-probe uses a 128-item subset of the val view; this uses
all of them, and runs on the test split's view built by build_eval_views --split test). Samples with vLLM through
g_probe.make_trainer and objdel.train.sample (T = 1, as the G-probe), scores with select_ckpt.score_items. Writes
<out>/final_eval.json (summary, schema 2) and <out>/items.jsonl (per-item metrics). Per item the K-sample mean; then
the mean over items with three 95 % percentile bootstrap intervals (B = 2,000, seed 0): `ci95_patient` resamples
patients (`group_id`; the gate CI, named by `gate_ci`), `ci95_nodule` resamples nodules, `ci95_item` single items
(the schema-1 `ci95`). Illegal samples score 0 on every geometric term (the reward's rule).
- iou_tc: IoU against the doctors' band (the reward's term); f5: the boundary score (calibrated so a new doctor
  scores 0.95); precision / recall / dice: pixels against the doctors' strict-majority outline;
- miss_rate: the share of items whose K-mean IoU_tc is below 0.1.
--resummarize DIR: the schema-2 summary of a finished run's items.jsonl (no model, no scoring), with group_id joined
from the LIDC store of the run's split; written to --out, never into DIR; exit 2 when an item has no group_id.
Reads untrusted data: run with python -I.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from objdel.eval import signal

METRICS = {"legal": ("samples", True), "iou_tc": ("samples", True), "f5": ("samples", True),
           "precision_vote": ("objects", True), "recall_vote": ("objects", True), "dice_vote": ("objects", True)}
RUN_KEYS = ("ckpt", "view", "split", "K", "seed")


def per_item(result: dict, obj: str = "nodule") -> dict:
    """K-sample means of one item's scores (metrics' output)."""
    src = {"samples": result["samples"], "objects": result["objects"][obj]}
    return {"item_id": result["item_id"], **{m: float(np.mean(src[s][m])) for m, (s, _) in METRICS.items()}}


def summarize(items: list[dict], groups: dict[str, str] | None = None, draws: int = 2000, seed: int = 0) -> dict:
    """Schema 2: per METRICS key the items hold, the mean and its item, nodule and patient (with `groups`,
    {item_id: group_id}) 95 % bootstrap intervals; the median of iou_tc; the miss rate. One draw set per unit."""
    ids = [it["item_id"] for it in items]
    D = {"ci95_nodule": signal.draws(ids, draws, seed),
         "ci95_patient": None if groups is None else signal.draws(ids, draws, seed, key=groups.__getitem__)}
    idx = np.random.default_rng(seed).integers(0, len(items), (draws, len(items)))

    def cis(v: np.ndarray) -> dict:
        lo, hi = np.percentile(v[idx].mean(1), [2.5, 97.5])
        return {"mean": float(v.mean()), **{k: None if d is None else signal.ci(v, d)[1:] for k, d in D.items()},
                "ci95_item": [float(lo), float(hi)]}

    out = {"schema": 2, "items": len(items), "nodules": len(set(map(signal.cluster, ids))),
           "patients": None if groups is None else len({groups[i] for i in ids}), "ci_level": 0.95, "boot": draws,
           "gate_ci": "ci95_patient"}
    for m, (_, higher) in METRICS.items():
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
    return out


def view_groups(view: Path, split: str) -> dict[str, str]:
    """{item_id: group_id} of a view's LIDC records (never from the prompt rows: TRL forwards their columns)."""
    from objdel.train import dataset

    return {r["item_id"]: r["group_id"] for r in dataset.records(view, "lidc", split)}


def resummarize(run: Path, out: Path, store: Path) -> None:
    """Schema-2 final_eval.json of a finished run (module docstring); exit 2 when an item has no group_id."""
    import pyarrow.parquet as pq

    doc = json.loads((run / "final_eval.json").read_text())
    with open(run / "items.jsonl") as f:
        items = [json.loads(line) for line in f]
    table = pq.read_table(store / "lidc" / f"{doc['split']}.parquet", columns=["item_id", "group_id"])
    found = {i: g for i, g in zip(*table.to_pydict().values(), strict=True) if g}
    lost = [it["item_id"] for it in items if it["item_id"] not in found]
    if lost:
        print(f"resummarize: {len(lost)} items have no group_id in {store}, first {lost[0]}", file=sys.stderr)
        sys.exit(2)
    s = summarize(items, found)
    res = {k: doc[k] for k in RUN_KEYS} | s | {"metrics_absent": [m for m in METRICS if m not in s],
                                               "resummarized_from": str(run)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "final_eval.json").write_text(json.dumps(res, indent=1))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ckpt", type=Path)
    ap.add_argument("--view", type=Path, help="an eval view directory (data/prepared layout)")
    ap.add_argument("--split", default="val")
    ap.add_argument("--train", type=Path, help="the prepared training data (make_trainer's dataset)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--resummarize", type=Path, metavar="DIR", help="a finished run: summarize its items.jsonl again")
    ap.add_argument("--store", type=Path, default=Path("data/store"), help="the item store (group_id per item)")
    ap.add_argument("--model", default="qwen3-vl-2b")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--per-call", type=int, default=64, help="items per vLLM generation call")
    a = ap.parse_args(argv)
    if a.resummarize:
        if a.out.resolve() == a.resummarize.resolve():
            ap.error("--out must not be the --resummarize folder")
        return resummarize(a.resummarize, a.out, a.store)
    if not (a.ckpt and a.view and a.train):
        ap.error("the arguments --ckpt, --view and --train are required (unless --resummarize)")

    import torch
    from PIL import Image

    from experiments.phase2 import g_probe
    from experiments.phase2.select_ckpt import score_items
    from objdel.train import dataset, sample, serve

    t0 = time.time()
    for root in (a.train, a.view):
        serve.preload(root, ["lidc"])
    pool = serve.fork_pool(serve.pool_size())
    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=dataset.DISAGREEMENT_OVERSAMPLE)
    rows = [r | {"image": Image.open(r["image"]).convert("RGB")} for r in dataset.prompts(a.view, "lidc", a.split)]
    args = argparse.Namespace(model=a.model, out=a.out, grpo="{}", train=a.train)
    trainer, _, _ = g_probe.make_trainer(args, a.ckpt, train, 1)
    torch.manual_seed(a.seed)
    done = {d["item_id"]: d for d in sample.sample(trainer, rows, a.k, items=a.per_call)}
    t_gen = time.time() - t0
    results, _ = score_items(pool, rows, done, a.view.resolve())
    pool.shutdown()
    items = [per_item(r) for r in results]
    doc = {"ckpt": str(a.ckpt), "view": str(a.view), "split": a.split, "K": a.k, "seed": a.seed,
           **summarize(items, view_groups(a.view, a.split)),
           "seconds": {"load_and_sample": t_gen, "total": time.time() - t0}}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "items.jsonl").write_text("".join(json.dumps(it) + "\n" for it in items))
    (a.out / "final_eval.json").write_text(json.dumps(doc, indent=1))
    print("FINAL_EVAL", a.split, json.dumps({m: round(doc[m]["mean"], 4) for m in METRICS}),
          "items", doc["items"], f"s {doc['seconds']['total']:.0f}")  # fmt: skip


if __name__ == "__main__":
    main()
