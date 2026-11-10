"""Teacher-force vLLM's samples through HF SmolVLM-256M (fp32): per-token logprob difference vLLM - HF.
Unconstrained samples should agree (sampler parity); grammar-constrained ones differ by -log Z_t >= 0 (the mass the
grammar removed), which is what TRL's IS correction would see if rollouts were constrained."""
import json
import os
import sys

import numpy as np
import torch
from transformers import AutoModelForImageTextToText, AutoProcessor

sys.path.insert(0, os.path.dirname(__file__))
from vllm_check import REPO, image  # noqa: E402

R = json.load(open(sys.argv[1]))
edge = R["edge"]
proc = AutoProcessor.from_pretrained(REPO, size={"longest_edge": edge})
model = AutoModelForImageTextToText.from_pretrained(REPO, dtype=torch.float32).eval()
enc = proc(text=R["prompt_text"], images=[image()], return_tensors="pt")
P = enc["input_ids"].shape[1]
out = {}
for name in ("free", "regex"):
    d = []
    for s in R[f"so_{name}"]["samples"]:
        ids = torch.tensor([enc["input_ids"][0].tolist() + s["ids"]])
        with torch.no_grad():
            logits = model(input_ids=ids, attention_mask=torch.ones_like(ids), pixel_values=enc["pixel_values"],
                           pixel_attention_mask=enc.get("pixel_attention_mask")).logits[0, P - 1:-1]
        lp = torch.log_softmax(logits.float(), -1).gather(1, torch.tensor(s["ids"])[:, None])[:, 0].numpy()
        d.append(np.asarray(s["lp"]) - lp)
    allv = np.concatenate(d)
    out[name] = {"tokens": int(allv.size), "mean_diff": float(allv.mean()), "max_abs": float(np.abs(allv).max()),
                 "p50": float(np.median(allv)), "min": float(allv.min()),
                 "seq_logratio_vllm_minus_hf": [round(float(x.sum()), 3) for x in d]}
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[1].replace(".json", "_hf.json"), "w"))
