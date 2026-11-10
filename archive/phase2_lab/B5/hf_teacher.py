"""B5 / H-A5-2(a), step 3 (training env, repo .venv: torch 2.14.1, transformers 5.20.1): teacher-force vLLM's samples.

For every sample written by vllm_parity.py: HF inputs from the processor exactly as TRL GRPO builds them
(apply_chat_template(add_generation_prompt=True, tokenize=True, return_dict=True), processor_kwargs), the sampled ids
appended, one forward pass in the same dtype, log_softmax in float32, the sampled token's logprob. Reports per-token
Δ = vLLM - HF and the per-sequence sum. Pass (pre-registered, REPORT.md): fp32 median |Δ| <= 1e-3, bf16 <= 0.02,
p99 per-sequence |ΣΔ| <= 0.5 nats.
    .venv/bin/python -I -m experiments.phase2.lab.B5.hf_teacher <vllm_out.json> <prompts.json> <out.json>
"""
import json
import sys

import numpy as np
import torch
from PIL import Image
from transformers import AutoModelForImageTextToText, AutoProcessor
from trl.data_utils import prepare_multimodal_messages

from objdel.train.dataset import _user


def main(vpath: str, ppath: str, out: str) -> None:
    V = json.load(open(vpath))
    P = json.load(open(ppath))["models"][V["model"]]
    dtype = {"float32": torch.float32, "bfloat16": torch.bfloat16}[V["dtype"]]
    model = AutoModelForImageTextToText.from_pretrained(V["repo"], dtype=dtype).eval()
    torch.set_num_threads(2)
    res = {"model": V["model"], "dtype": V["dtype"], "datasets": {}}
    all_tok, all_seq = [], []
    for dataset, D in V["datasets"].items():
        items = {x["item_id"]: x for x in P["datasets"][dataset]["items"]}
        proc = AutoProcessor.from_pretrained(V["repo"], **P["datasets"][dataset]["processor_kwargs"])
        tok_d, seq_d = [], []
        for it in D["samples"]:
            x = items[it["item_id"]]
            img = Image.open(x["image"]).convert("RGB")
            prompt = prepare_multimodal_messages([_user(dataset)], images=[img])
            enc = proc.apply_chat_template(conversation=[prompt], add_generation_prompt=True, tokenize=True,
                                           return_dict=True, return_tensors="pt")
            n_prompt = enc["input_ids"].shape[1]
            assert enc["input_ids"][0].tolist() == x["gen_ids"]
            for s in it["samples"]:
                comp = torch.tensor([s["ids"]])
                kw = {}
                for k, v in enc.items():
                    if isinstance(v, torch.Tensor) and v.dim() == 2 and v.shape == (1, n_prompt):
                        pad = torch.ones_like(comp) if k == "attention_mask" else (
                            comp if k == "input_ids" else torch.zeros_like(comp))
                        kw[k] = torch.cat([v, pad.to(v.dtype)], 1)
                    else:
                        kw[k] = v.to(dtype) if isinstance(v, torch.Tensor) and v.is_floating_point() else v
                with torch.no_grad():
                    logits = model(**kw).logits[0, n_prompt - 1:-1]
                lp = torch.log_softmax(logits.float(), -1).gather(1, comp[0][:, None])[:, 0].numpy()
                d = np.asarray(s["lp"], dtype=np.float64) - lp
                tok_d.append(d)
                seq_d.append(float(d.sum()))
        t = np.abs(np.concatenate(tok_d))
        res["datasets"][dataset] = {"tokens": int(t.size), "sequences": len(seq_d), "median_abs": float(np.median(t)),
                                    "p99_abs": float(np.quantile(t, 0.99)), "max_abs": float(t.max()),
                                    "mean_signed": float(np.concatenate(tok_d).mean()),
                                    "seq_sum": [round(v, 5) for v in seq_d]}
        all_tok += tok_d
        all_seq += seq_d
        print(V["model"], dataset, {k: v for k, v in res["datasets"][dataset].items() if k != "seq_sum"}, flush=True)
    t = np.abs(np.concatenate(all_tok))
    s = np.abs(np.asarray(all_seq))
    res["all"] = {"tokens": int(t.size), "sequences": int(s.size), "median_abs": float(np.median(t)),
                  "p99_abs": float(np.quantile(t, 0.99)), "max_abs": float(t.max()),
                  "seq_p99_abs": float(np.quantile(s, 0.99)), "seq_max_abs": float(s.max())}
    print(V["model"], "all", res["all"], flush=True)
    json.dump(res, open(out, "w"))


if __name__ == "__main__":
    main(*sys.argv[1:4])
