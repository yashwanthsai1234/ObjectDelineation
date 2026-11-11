"""Prepare the training set once: model frames (PNG), SFT targets and expert caches with hashes (objdel.data.prepare).

    .venv/bin/python -I -m experiments.phase1.prepare_rl [--store data/store] [--out data/prepared] [--procs 3]

Writes <out>/images, <out>/caches, <out>/<dataset>/<split>.jsonl and <out>/manifest.json (git-ignored), then prints
the per-split counts and the answer-token evidence for issue I5 (targets over each model's max_completion_length).
Reads untrusted data: run with python -I.

References:
- objdel.data.prepare.prepare (frames, targets, ExpertCache files, manifest) and impl/M7_train.md (I5 table)
- transformers AutoTokenizer (add_special_tokens=False counts the answer alone):
  https://huggingface.co/docs/transformers/main_classes/tokenizer
- TRL GRPOConfig max_completion_length (generation stops there, the completion is truncated):
  https://huggingface.co/docs/trl/main/en/grpo_trainer
"""
import argparse
import json
import time
from pathlib import Path

from objdel.core.io import read_json
from objdel.data import prepare as data_prepare
from objdel.policy import prompt as policy_prompt
from objdel.policy import rows as policy_rows


def token_budget(out):
    """Issue I5: per model, train targets whose answer (+1 end-of-turn token) exceeds max_completion_length."""
    from transformers import AutoTokenizer

    from objdel.policy.models import MODELS
    res = {}
    for m, entry in MODELS.items():
        tok = AutoTokenizer.from_pretrained(entry.hf_id)
        for ds in sorted({k.split("/")[0] for k in read_json(Path(out, "manifest.json"))["sets"]}):
            n = [len(tok(t, add_special_tokens=False)["input_ids"]) + 1
                 for r in policy_rows.records(out, ds, "train") for t in r["targets"] if t]
            over = sorted(x for x in n if x > entry.max_completion)
            res[f"{m}/{ds}"] = {"targets": len(n), "budget": entry.max_completion, "over": len(over), "max": max(n),
                                "over_lengths": over}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--out", default="data/prepared")
    ap.add_argument("--splits", default="train,val")
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args()
    t0 = time.time()
    m = data_prepare.prepare(a.store, a.out, a.splits.split(","), a.procs, prompts=policy_prompt.PROMPTS)
    print(json.dumps(m["sets"], indent=1), f"\n{time.time() - t0:.0f}s", flush=True)
    print(json.dumps(m["files"], indent=1))
    print(json.dumps(token_budget(a.out), indent=1))


if __name__ == "__main__":
    main()
