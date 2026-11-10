"""R0 (PLAN §4.7 step 8, §4.9; D13 HF path): HF generate's sampling rate on a model's base weights, for project_r1.

    python -I -m experiments.phase2.hf_rate --model lfm2.5-vl-450m --view <eval_views>/lidc-val-shifted \
        --out <run_dir>/hf-rate-lfm [--items 32] [--k 8]

The seeded --items subset (policy.evaluate.view_rows) at n = --k through sft.select.sample_view with the HF sampler, as
R1's select_ckpt runs it (hf_sample at its default --batch_items; CUDA untouched here). hf_sample --timing ->
<out>/hf_rate.json: s_per_sample (generate seconds per item-sample), load_s (one child's model load), tokens_mean
and length_share (base weights answer unlike a trained branch, and a batch runs as long as its longest sample, so
the rate holds for the base model's answer lengths), plus model, weights and view.
"""

import argparse
import json
import sys
from pathlib import Path

from experiments.phase2 import r0_parity
from objdel.policy import evaluate
from objdel.sft import select as SC


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, help="configs.MODELS key")
    ap.add_argument("--view", type=Path, required=True, help="eval view with lidc/val (lidc-val-shifted)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--items", type=int, default=32)
    ap.add_argument("--k", type=int, default=8, help="samples per item")
    ap.add_argument("--seed", type=int, default=0, help="sampler seed")
    ap.add_argument("--timeout_s", type=float, default=SC.TIMEOUT_S, help="sampler timeout per sample, plus 600 s")
    ap.set_defaults(sampler="hf")  # what SC.sample_view reads
    a = ap.parse_args(argv)
    w, timing = r0_parity.weights(a.model), a.out / "timing.json"
    a.out.mkdir(parents=True, exist_ok=True)
    SC.sample_view(a, w, evaluate.view_rows(a.view, a.items), a.out / "samples", a.k, extra=("--timing", str(timing)))
    doc = json.loads(timing.read_text()) | {"model": a.model, "weights": str(w), "view": str(a.view)}
    (a.out / "hf_rate.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(doc))
    return 0


if __name__ == "__main__":
    sys.exit(main())
