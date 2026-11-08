"""Final evaluation of a tuned model on a whole split: K samples per item, scored by the checker.

    python -I -m experiments.phase2.final_eval --ckpt DIR --view data/eval_views/lidc-val-shifted --split val \
        --train data/prepared --out runs/eval/val [--k 8] [--seed 11] [--per-call 64]

Every RL-eligible item of the view's <split> records (the G-probe uses a 128-item subset of the val view; this uses
all of them, and runs on the test split's view built by build_eval_views --split test). Samples with vLLM through
g_probe.make_trainer and objdel.train.sample (T = 1, as the G-probe), scores with select_ckpt.score_items. Writes
<out>/final_eval.json (summary) and <out>/items.jsonl (per-item metrics). Per item the K-sample mean; then the mean
over items with a 95 % item bootstrap interval. Illegal samples score 0 on every geometric term (the reward's rule).
- iou_tc: IoU against the doctors' band (the reward's term); f5: the boundary score (calibrated so a new doctor
  scores 0.95); precision / recall / dice: pixels against the doctors' strict-majority outline.
Reads untrusted data: run with python -I.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

TERMS = ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote")


def per_item(result: dict, obj: str = "nodule") -> dict:
    """K-sample means of one item's scores (metrics' output)."""
    s, o = result["samples"], result["objects"][obj]
    return {"item_id": result["item_id"], "legal": float(np.mean(s["legal"])), "iou_tc": float(np.mean(s["iou_tc"])),
            "f5": float(np.mean(s["f5"])), **{t: float(np.mean(o[t])) for t in TERMS[2:]}}  # fmt: skip


def summarize(items: list[dict], draws: int = 2000, seed: int = 0) -> dict:
    """Mean over items of each per-item score, its 95 % item bootstrap interval, and the median of iou_tc."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(items), (draws, len(items)))
    out = {"items": len(items)}
    for t in ("legal", *TERMS):
        v = np.array([it[t] for it in items])
        lo, hi = np.percentile(v[idx].mean(1), [2.5, 97.5])
        out[t] = {"mean": float(v.mean()), "ci95": [float(lo), float(hi)]}
    out["iou_tc"]["median"] = float(np.median([it["iou_tc"] for it in items]))
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--view", type=Path, required=True, help="an eval view directory (data/prepared layout)")
    ap.add_argument("--split", default="val")
    ap.add_argument("--train", type=Path, required=True, help="the prepared training data (make_trainer's dataset)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--model", default="qwen3-vl-2b")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--per-call", type=int, default=64, help="items per vLLM generation call")
    a = ap.parse_args(argv)

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
           **summarize(items), "seconds": {"load_and_sample": t_gen, "total": time.time() - t0}}  # fmt: skip
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "items.jsonl").write_text("".join(json.dumps(it) + "\n" for it in items))
    (a.out / "final_eval.json").write_text(json.dumps(doc, indent=1))
    print("FINAL_EVAL", a.split, json.dumps({t: round(doc[t]["mean"], 4) for t in ("legal", *TERMS)}),
          "items", doc["items"], f"s {doc['seconds']['total']:.0f}")  # fmt: skip


if __name__ == "__main__":
    main()
