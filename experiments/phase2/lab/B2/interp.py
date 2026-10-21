"""H-B2-2: the LFM2.5-VL-450M 150-step B1 adapter evaluated with its LoRA delta scaled by alpha toward the base
(every lora_B weight times alpha, then merged: W = W0 + alpha * s * B A), through B1's evaluation path unchanged
(train_eval.evaluate: 64 LIDC val items x K = 8 at T = 1, centred and shifted views, common.sample_k, seed 0).

    HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B2.interp --alpha 0.75 [--eval-items 64]
    HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B2.interp --alpha 1 --check 4

--check N: alpha = 1 on the first N centred items only, compared with B1's stored samples (path sanity and timing; not
an eval). Writes <scratch>/phase2/B2/runs/lfm_a<alpha>/{result,per_item,samples}.json in B1's formats.
Reads untrusted data: run with python -I.

References:
- PEFT LoRA (scaling = lora_alpha / r, merge_and_unload): https://huggingface.co/docs/peft/main/en/package_reference/lora
- WiSE-FT (weight-space interpolation of fine-tuned and initial weights): https://arxiv.org/abs/2109.01903 ; for
  reasoning SFT -> RL starts: https://arxiv.org/abs/2504.10478
"""

import argparse
import json
import time
from argparse import Namespace

import torch
from peft import PeftModel
from transformers import AutoModelForImageTextToText

from ..B1 import common as C
from ..B1.train_eval import evaluate
from .prereg import INTERP_RUN

OUT = C.SCRATCH.parent / "B2" / "runs"
ARM = "lfm2.5-vl-450m"


def scaled_model(alpha: float):
    base = AutoModelForImageTextToText.from_pretrained(C.ARMS[ARM]["hf"], dtype=torch.float32)
    peft = PeftModel.from_pretrained(base, str(C.SCRATCH / "runs" / INTERP_RUN / "adapter"))
    n = 0
    with torch.no_grad():
        for name, p in peft.named_parameters():
            if "lora_B" in name:
                p.mul_(alpha)
                n += 1
    print("scaled lora_B tensors:", n, "alpha", alpha, flush=True)
    return peft.merge_and_unload()


def check(model, proc, n: int):
    """alpha = 1 on the first n centred items vs B1's stored samples (same seed, same order)."""
    from objdel.train import dataset as D

    recs = C.eval_records("lidc", 64)[:n]
    V = C.views("lidc", recs)
    trainer = C.eval_trainer(model, proc, ARM, str(OUT / "check_gen"), 1024)
    rows = [{"prompt": [D._user("lidc")], "image": x["image"], "key": x["item_id"]} for x in V["centred"]]
    torch.manual_seed(0)
    t = time.time()
    gen = C.sample_k(trainer, rows, 8, C.zero_led_ids(proc.tokenizer))
    s = (time.time() - t) / n
    ref = json.loads((C.SCRATCH / "runs" / INTERP_RUN / "samples.json").read_text())["centred"][:n]
    same = [sum(a == b for a, b in zip(g["samples"], r["samples"], strict=True)) for g, r in zip(gen, ref, strict=True)]
    print(
        "CHECK",
        json.dumps(
            {
                "s_per_item": round(s, 2),
                "identical_samples_per_item": same,
                "projected_min_128_items": round(s * 128 / 60, 1),
            }
        ),
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--eval-items", type=int, default=64)
    ap.add_argument("--check", type=int, default=0)
    a = ap.parse_args()
    C.cap_mps()
    proc = C.processor(ARM, "lidc")
    t0 = time.time()
    model = scaled_model(a.alpha).to("mps")
    if a.check:
        check(model, proc, a.check)
        return
    out = OUT / f"lfm_a{a.alpha:g}"
    out.mkdir(parents=True, exist_ok=True)
    args = Namespace(arm=ARM, eval_items=a.eval_items, k=8, max_new=1024, seed=0)
    res, per_item, raw = evaluate(args, model, proc, out)
    (out / "per_item.json").write_text(json.dumps(per_item))
    (out / "samples.json").write_text(json.dumps(raw))
    (out / "result.json").write_text(
        json.dumps({"args": vars(a), "eval": res, "wall_s": round(time.time() - t0, 1)}, indent=1, default=str)
    )
    print("DONE", out.name, round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
