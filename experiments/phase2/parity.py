"""H-A5-2(a) logprob parity on the sampled checkpoint (promoted from lab/B5/hf_teacher.py, S1-16; R0 runs it).

    .venv/bin/python -I -m experiments.phase2.parity --prompts P.jsonl --completions C.jsonl --model DIR --out R.json

Inputs: the bridge files of one sampling run (C-BRIDGE): the prompts JSONL of objdel.train.sample.export_prompts and
the completions of `sample_vllm.py --logprobs` (each line's "logprobs": vLLM's processed logprob of every sampled
token). Every completion is teacher-forced through HF in the training env: the processor applied to the exported
prompt text and image (the prompt ids TRL GRPO builds, B5: equal to vLLM's), the sampled ids appended, one forward
pass in --dtype, log_softmax in float32, the sampled token's logprob. Δ = vLLM - HF per token and its per-sequence
sum. Pass (pre-registered in lab/B5/REPORT.md): median |Δ| <= 1e-3 in float32 (0.02 in bfloat16) and the p99 of
the per-sequence |ΣΔ| <= 0.5 nats.
"""
import argparse
import json
from pathlib import Path

import numpy as np

MEDIAN = {"float32": 1e-3, "bfloat16": 0.02}
SEQ_P99 = 0.5


def teacher(model, processor, prompt: str, image, ids: list[int]) -> np.ndarray:
    """HF logprobs (float32 log_softmax) of the sampled ids after the prompt."""
    import torch
    enc = processor(text=[prompt], images=[[image]], add_special_tokens=False, return_tensors="pt")
    n, comp = enc["input_ids"].shape[1], torch.tensor([ids])
    kw = {}
    for k, v in enc.items():
        if isinstance(v, torch.Tensor) and v.shape == (1, n):  # token-aligned: extend over the completion
            tail = comp if k == "input_ids" else torch.ones_like(comp) if k == "attention_mask" else 0 * comp
            v = torch.cat([v, tail.to(v.dtype)], 1)
        elif isinstance(v, torch.Tensor) and v.is_floating_point():
            v = v.to(model.dtype)
        kw[k] = v
    with torch.no_grad():
        logits = model(**kw).logits[0, n - 1:-1].float()
    return torch.log_softmax(logits, -1).gather(1, comp[0][:, None])[:, 0].numpy()


def report(deltas: list[np.ndarray], dtype: str) -> dict:
    t, s = np.abs(np.concatenate(deltas)), np.abs([d.sum() for d in deltas])
    r = {"dtype": dtype, "tokens": int(t.size), "sequences": len(deltas), "median_abs": float(np.median(t)),
         "p99_abs": float(np.quantile(t, 0.99)), "max_abs": float(t.max()), "seq_p99_abs": float(np.quantile(s, 0.99))}
    return r | {"pass": r["median_abs"] <= MEDIAN[dtype] and r["seq_p99_abs"] <= SEQ_P99}


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prompts", required=True, help="JSONL from objdel.train.sample.export_prompts")
    ap.add_argument("--completions", required=True, help="JSONL from sample_vllm.py --logprobs")
    ap.add_argument("--model", required=True, help="the checkpoint directory vLLM sampled (with its processor)")
    ap.add_argument("--dtype", choices=sorted(MEDIAN), default="float32")
    ap.add_argument("--out", required=True, help="report JSON")
    a = ap.parse_args(argv)
    import torch
    from PIL import Image
    from transformers import AutoModelForImageTextToText, AutoProcessor

    prompts = {x["id"]: x for x in map(json.loads, Path(a.prompts).read_text().splitlines())}
    kwargs = {json.dumps(x["processor_kwargs"], sort_keys=True) for x in prompts.values()}
    if len(kwargs) != 1:
        raise SystemExit("parity: one processor configuration per prompts file")
    processor = AutoProcessor.from_pretrained(a.model, **json.loads(kwargs.pop()))
    model = AutoModelForImageTextToText.from_pretrained(a.model, dtype=getattr(torch, a.dtype)).eval()
    deltas = []
    for c in map(json.loads, Path(a.completions).read_text().splitlines()):
        if "logprobs" not in c:
            raise SystemExit("parity: completions without logprobs (run sample_vllm.py --logprobs)")
        x = prompts[c["id"]]
        hf = teacher(model, processor, x["prompt"], Image.open(x["image"]).convert("RGB"), c["token_ids"])
        deltas.append(np.asarray(c["logprobs"], np.float64) - hf)
    r = report(deltas, a.dtype) | {"model": a.model}
    Path(a.out).write_text(json.dumps(r, indent=1) + "\n")
    print(json.dumps(r))
    return r


if __name__ == "__main__":
    main()
