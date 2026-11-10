"""R0 step 7 / G-parity (PLAN §4.7, §5.1, D13): vLLM vs HF on a model's base weights, a few LIDC val items.

    python -I -m experiments.phase2.r0_parity --model lfm2.5-vl-450m --view <eval_views>/lidc-val-shifted \
        --out <run_dir>/parity-lfm [--items 32] [--logprob-items 8] [--k 8]

Base weights: snapshot_download(configs.MODELS[model]), one local directory that vLLM and transformers both load. The seeded --items subset (objdel.train.evaluate.eval_rows) is
sampled at n = --k by `sample_vllm.py --logprobs` through select_ckpt.sample_view (CUDA untouched, timeout, retries,
completeness; files in <out>/samples/). Then parity.compare in bfloat16: (a) H-A5-1 prompt ids equal on every item;
(b) H-A5-2 HF teacher forcing of the first --logprob-items x n (PLAN: 8 x 8) within parity.MEDIAN["bfloat16"]
(0.02); parity.SEQ_P99 (0.5) is reported only (D13, 2026-11-05). Writes <out>/parity.json (compare's report plus model, weights, view, items,
logprob_items, k); exit 1 on a fail (-> the HF sampler is that model's gate path, D13).
"""

import argparse
import json
import sys
from pathlib import Path

from experiments.phase2 import parity
from experiments.phase2 import select_ckpt as SC


def weights(model: str) -> Path:
    from huggingface_hub import snapshot_download

    from objdel.train import configs

    return Path(snapshot_download(configs.MODELS[model]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, help="configs.MODELS key")
    ap.add_argument("--view", type=Path, required=True, help="eval view with lidc/val (e.g. lidc-val-shifted)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--items", type=int, default=32, help="seeded items for the prompt-id check")
    ap.add_argument("--logprob-items", type=int, default=8, help="the first of them teacher-forced (x n = --k)")
    ap.add_argument("--k", type=int, default=8, help="samples per item")
    ap.add_argument("--seed", type=int, default=0, help="sampler seed")
    ap.add_argument("--timeout_s", type=float, default=2.0, help="sampler child timeout per sample, plus 600 s")
    ap.set_defaults(sampler="vllm")  # what SC.sample_view reads
    a = ap.parse_args(argv)
    w, rows, samples = weights(a.model), SC.eval_rows(a.view, a.items), a.out / "samples"
    SC.sample_view(a, w, rows, samples, a.k, extra=("--logprobs",))  # before anything touches CUDA (D13)
    teach = {r["key"] for r in rows[: a.logprob_items]}
    doc = parity.compare(w, samples / "prompts.jsonl", samples / "completions.jsonl", "bfloat16", teach) | {
        "model": a.model, "weights": str(w), "view": str(a.view), "items": len(rows), "k": a.k,
        "logprob_items": len(teach)}  # fmt: skip
    (a.out / "parity.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps({k: doc[k] for k in ("model", "pass", "median_abs", "seq_p99_abs", "prompt_ids")}))
    return 0 if doc["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
