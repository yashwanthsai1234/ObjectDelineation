"""C11 step 1: size the per-row answer cap on the Mac ($0), on val only.

For each val image cap_i = min(L, ceil(factor x T_i) + margin) (objdel.rl.runaway). Reports, per factor, the
share of found answers (IoU_tc >= 0.1) longer than cap_i (gate G-RC1: < 0.5 %) and the decode seconds per step the cap
would remove (4.3 ms per token of the longest answer of a block of 32 images, replayed on the logged lengths).

.venv/bin/python -I -m experiments.phase2.runaway_cut --train data/prepared --answers runs/a3/eval_val/answers.jsonl \
    --model qwen3-vl-2b --out runs/c11
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np

from objdel.policy import models as configs
from objdel.policy import rows as dataset
from objdel.rl import runaway

FOUND = 0.1  # a found answer: IoU_tc at or above this (the BRIEF counts IoU < 0.1 as a miss)
SECONDS_PER_TOKEN = 0.0043  # B_gpu_audit.md:40


def _jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines()]


def refuse_test(answers: Path, split: str) -> None:
    """The cap is chosen on val only: refuse split test and a test directory."""
    if split == "test" or "test" in answers.parent.name:
        raise ValueError(f"{answers}, split {split}: the test split is not used to choose the cap")


def size(records, answer_rows, metric_rows, tok, limit: int, factors, block: int = 32) -> dict:
    """{factor: {cap_min, found, found_above_cap, share, seconds_removed_per_step}}; tok(text) counts no end token."""
    outline = {r["item_id"]: runaway.outline_tokens(tok, r["targets"]) for r in records}
    met = {m["item_id"]: m["samples"] for m in metric_rows}
    ans = [(a["item_id"], np.asarray(a["tokens"]), np.asarray(met[a["item_id"]]["iou_tc"])) for a in answer_rows]
    out = {}
    for f in factors:
        caps = dict(zip(outline, runaway.runaway_caps(outline.values(), limit, f, runaway.MARGIN), strict=True))
        assert all(outline[i] <= caps[i] for i in outline), f"factor {f}: a doctor outline is longer than its cap"
        found = above = 0
        longest, capped, saved = [], [], []
        for k, (i, t, iou) in enumerate(ans):
            m = iou >= FOUND
            found += int(m.sum())
            above += int((m & (t > caps[i])).sum())
            longest.append(t.max())
            capped.append(np.minimum(t, caps[i]).max())
            if (k + 1) % block == 0 or k + 1 == len(ans):
                saved.append(SECONDS_PER_TOKEN * (max(longest) - max(capped)))
                longest, capped = [], []
        out[f] = {"cap_min": min(caps.values()), "found": found, "found_above_cap": above,
                  "share": above / max(found, 1), "seconds_removed_per_step": float(np.mean(saved))}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", type=Path, required=True)
    ap.add_argument("--answers", type=Path, required=True)
    ap.add_argument("--model", required=True, choices=sorted(configs.MODELS))
    ap.add_argument("--split", default="val")
    ap.add_argument("--factors", type=float, nargs="+", default=[1.5, 2.0])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    refuse_test(a.answers, a.split)
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    from transformers import AutoProcessor

    tok = AutoProcessor.from_pretrained(configs.MODELS[a.model], **configs.processor_kwargs(a.model, "lidc")).tokenizer
    res = size(dataset.records(a.train, "lidc", a.split), _jsonl(a.answers), _jsonl(a.answers.parent / "metrics.jsonl"),
               lambda text: dataset.target_tokens(tok, text), configs.MAX_COMPLETION[a.model], a.factors)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "runaway_cut.json").write_text(json.dumps({str(f): r for f, r in res.items()}, indent=1))
    for f, r in res.items():
        print(f"factor {f}: {100 * r['share']:.3f} % of {r['found']} found answers above cap "
              f"(G-RC1 {'pass' if r['share'] < 0.005 else 'FAIL'}); {r['seconds_removed_per_step']:.2f} s per step")


if __name__ == "__main__":
    main()
