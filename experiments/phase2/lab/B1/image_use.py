"""POST-HOC diagnostic (CPU): does a trained checkpoint's next-token distribution depend on the image? For each of n
eval items, the base model + saved LoRA adapter scores the first answer tokens of the generation prompt with the
centred image, the shifted image and an all-black image; reports the mean total-variation distance between the
next-token distributions at each of the first T teacher-forced positions of the item's first rater target, and the
same through TRL's GRPOTrainer generation (greedy) to check the image reaches generate.
Reads untrusted data: run with python -I.

    .venv/bin/python -I -m experiments.phase2.lab.B1.image_use --run smolvlm-256m_s0_n150 --arm smolvlm-256m
"""
import argparse
import json

import numpy as np
import torch
from peft import PeftModel
from PIL import Image
from transformers import AutoModelForImageTextToText

from objdel.train import dataset as D

from . import common as C


def dists(model, proc, image, target: str, T: int) -> np.ndarray:
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": D.PROMPTS["lidc"]}]}]
    prompt = proc.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
    ids = proc(text=[prompt], images=[[image]], return_tensors="pt")
    tgt = proc.tokenizer(target, add_special_tokens=False, return_tensors="pt")["input_ids"][:, :T]
    batch = {k: v for k, v in ids.items()}
    n = batch["input_ids"].shape[1]
    batch["input_ids"] = torch.cat([batch["input_ids"], tgt], 1)
    batch["attention_mask"] = torch.ones_like(batch["input_ids"])
    if "mm_token_type_ids" in batch:
        batch["mm_token_type_ids"] = torch.cat([batch["mm_token_type_ids"], torch.zeros_like(tgt)], 1)
    with torch.no_grad():
        logits = model(**batch).logits[0, n - 1:n - 1 + tgt.shape[1]]
    return torch.softmax(logits.float(), -1).numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run directory name, or base (no adapter)")
    ap.add_argument("--arm", required=True, choices=list(C.ARMS))
    ap.add_argument("--items", type=int, default=8)
    ap.add_argument("--T", type=int, default=24)
    a = ap.parse_args()
    torch.set_num_threads(2)
    proc = C.processor(a.arm, "lidc")
    base = AutoModelForImageTextToText.from_pretrained(C.ARMS[a.arm]["hf"], dtype=torch.float32)
    model = base.eval() if a.run == "base" else \
        PeftModel.from_pretrained(base, str(C.SCRATCH / "runs" / a.run / "adapter")).merge_and_unload().eval()
    recs = C.eval_records("lidc", 64)[:a.items]
    V = C.views("lidc", recs)
    tv = {"shifted": [], "black": [], "other_item": []}
    for i, (rec, c, s) in enumerate(zip(recs, V["centred"], V["shifted"], strict=True)):
        target = next(t for t in rec["targets"] if t)
        p = dists(model, proc, c["image"], target, a.T)
        other = V["centred"][(i + 1) % len(recs)]["image"]
        black = Image.new("RGB", c["image"].size)
        for name, im in (("shifted", s["image"]), ("black", black), ("other_item", other)):
            q = dists(model, proc, im, target, a.T)
            tv[name].append(float(0.5 * np.abs(p - q).sum(-1).mean()))
    out = {"run": a.run, "items": a.items, "T": a.T,
           **{f"tv_centred_vs_{k}": round(float(np.mean(v)), 4) for k, v in tv.items()}}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
