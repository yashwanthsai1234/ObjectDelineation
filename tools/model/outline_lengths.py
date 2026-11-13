"""True outline-length count for the completion limit (Implementation.md A8): the doctors' target texts on every
root, and the model's answers.

Outline length = the token ids of a target text without special tokens + 1 (the end token; the sft_rows rule), so
the LIDC train maximum is 1,182. A limit L loses the outlines longer than L. With --answers (A3's answers.jsonl and
the metrics.jsonl beside it) the answers' `tokens` already count the end token when the model stopped; 1,536-token
answers are the truncated ones. The limit choice stays with the owner.

.venv/bin/python -I -m tools.model.outline_lengths --train data/prepared \
    --views data/views/lidc-train-shifted data/views/lidc-train-rescaled --model qwen3-vl-2b \
    --datasets lidc --splits train val [--answers runs/a3/eval_val/answers.jsonl] --out runs/a8

Test surface: kept by the CLEAN16 decision (2026-11-10); only tests import it, so a change here is checked there.
"""

import argparse
import json
import os
from collections.abc import Callable, Iterable
from pathlib import Path

import numpy as np

from objdel.core.hashing import file_sha256
from objdel.core.io import read_jsonl
from objdel.policy import models as configs
from objdel.policy import rows as dataset
from objdel.verify import store as verify_store

LIMITS = (256, 512, 800, 1000, 1200, 1536)
LONGEST = 20
GOOD = 0.5  # n_legal_good: legal and IoU_tc >= this


def _above(lengths: list[int]) -> dict[int, int]:
    a = np.asarray(lengths)
    return {L: int((a > L).sum()) for L in LIMITS}


def count(records: Iterable[dict], tokens_of: Callable[[str], int],
          voting: Callable[[dict], set[int]] | None = None) -> dict:
    """Lengths of every non-None target (tokens_of as given); voting=None counts every target as voting."""
    n, lengths, voted, longest = 0, [], [], []
    for r in records:
        n += 1
        vote = None if voting is None else voting(r)
        for i, text in enumerate(r["targets"]):
            if text is None:
                continue
            t = tokens_of(text)
            lengths.append(t)
            longest.append({"item_id": r["item_id"], "rater": i, "tokens": t})
            if vote is None or i in vote:
                voted.append(t)
    a = np.asarray(lengths, float)
    stats = {k: None for k in ("mean", "p50", "p90", "p99", "p99_9", "max")} if not lengths else {
        "mean": float(a.mean()), "max": int(a.max()),
        **{k: float(np.percentile(a, q)) for k, q in (("p50", 50), ("p90", 90), ("p99", 99), ("p99_9", 99.9))}}
    return {"n_records": n, "n_targets": len(lengths), "n_voting_targets": len(voted), **stats,
            "above": _above(lengths), "above_voting": _above(voted),
            "longest": sorted(longest, key=lambda d: -d["tokens"])[:LONGEST]}


def answers(answer_rows: Iterable[dict], metric_rows: Iterable[dict]) -> dict:
    """Per limit L: the share of all answers with tokens > L, and over those answers the legal share, the mean IoU_tc
    and the count legal with IoU_tc >= 0.5 (None when no answer is above L). Joined by item_id; KeyError for an
    answer with no metric row, ValueError for a repeated item."""
    met = {m["item_id"]: m["samples"] for m in metric_rows}
    tok, legal, iou, seen = [], [], [], set()
    for a in answer_rows:
        i = a["item_id"]
        if i in seen:
            raise ValueError(f"answers: item {i} is repeated")
        seen.add(i)
        s = met[i]
        if not len(a["tokens"]) == len(s["legal"]) == len(s["iou_tc"]):
            raise ValueError(f"answers: item {i} has {len(a['tokens'])} token counts, {len(s['legal'])} scores")
        tok += a["tokens"]
        legal += s["legal"]
        iou += s["iou_tc"]
    t, lg, io = np.asarray(tok), np.asarray(legal, float), np.asarray(iou, float)
    out = {"n_answers": len(t)}
    for L in LIMITS:
        m = t > L
        out[L] = {"share_above": float(m.mean()), "legal_share": float(lg[m].mean()) if m.any() else None,
                  "mean_iou_tc": float(io[m].mean()) if m.any() else None,
                  "n_legal_good": int(((lg > 0) & (io >= GOOD) & m).sum())}
    return out


def tokens_of(tokenizer) -> Callable[[str], int]:
    """Outline length of a text: target tokens without special tokens + 1 end token."""
    return lambda text: dataset.target_tokens(tokenizer, text) + 1


def report(train: Path, views: list[Path], datasets: list[str], splits: list[str], tok: Callable[[str], int],
           answer_files: list[Path]) -> dict:
    """{"roots": {root: {"manifest_sha256", dataset: {split: count}}}, "model_answers": {file: answers}}; a
    (dataset, split) with no file under a root is left out."""
    roots = {}
    for root in [train, *views]:
        row = {"manifest_sha256": ((p := root / "manifest.json").is_file() and file_sha256(p) or None)}
        for ds in datasets:
            got = {sp: count(dataset.records(root, ds, sp), tok,
                             lambda r, ds=ds, root=root: verify_store.voting(root, r, ds))
                   for sp in splits if (root / ds / f"{sp}.jsonl").is_file()}
            if got:
                row[ds] = got
        roots[str(root)] = row
    return {"roots": roots, "model_answers": {str(f): answers(read_jsonl(f), read_jsonl(f.parent / "metrics.jsonl"))
                                              for f in answer_files}}


def table(doc: dict) -> str:
    """The owner's decision lines: what each limit loses on each root and cuts from the model's answers."""
    lines = []
    for root, row in doc["roots"].items():
        for ds, by in row.items():
            if ds == "manifest_sha256":
                continue
            for sp, c in by.items():
                for L in LIMITS:
                    n = c["above"][L]
                    lines.append(f"limit {L} loses {n} doctor outlines ({100 * n / max(c['n_targets'], 1):.3f} %) "
                                 f"on {root} {ds}/{sp} (voting {c['above_voting'][L]}; max {c['max']})")
    for f, a in doc["model_answers"].items():
        for L in LIMITS:
            m = round(a[L]["share_above"] * a["n_answers"])
            lines.append(f"limit {L} cuts {m} model answers of {a['n_answers']} in {f}, of which "
                         f"{a[L]['n_legal_good']} were legal with IoU_tc >= {GOOD}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", type=Path, required=True)
    ap.add_argument("--views", type=Path, nargs="*", default=[])
    ap.add_argument("--model", required=True, choices=sorted(configs.MODELS))
    ap.add_argument("--datasets", nargs="+", default=["lidc"])
    ap.add_argument("--splits", nargs="+", default=["train", "val"])
    ap.add_argument("--answers", type=Path, nargs="*", default=[])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    from transformers import AutoProcessor

    # the tokenizer does not depend on the dataset's image size: one load (export_prompts' call)
    tok = AutoProcessor.from_pretrained(configs.get(a.model).hf_id, **configs.processor_kwargs(a.model)).tokenizer
    doc = {"model": a.model} | report(a.train, a.views, a.datasets, a.splits, tokens_of(tok), a.answers)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "outline_lengths.json").write_text(json.dumps(doc, indent=1))
    print(table(doc))


if __name__ == "__main__":
    main()
