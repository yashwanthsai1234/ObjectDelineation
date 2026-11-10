"""H-A5-2(a) logprob parity on a sampled checkpoint (from archive/phase2_lab/B5/hf_teacher.py, S1-16; R0 called compare through
r0_parity, archived in CLEAN8).

    .venv/bin/python -I -m experiments.phase2.parity --prompts P.jsonl --completions C.jsonl --model DIR --out R.json

Inputs (C-BRIDGE): export_prompts' prompts JSONL and `sample_vllm.py --logprobs` completions (vLLM's processed
logprob of every sampled token). Each completion is teacher-forced through HF: the processor on the exported prompt
text and image (TRL GRPO's prompt ids, B5: equal to vLLM's), the sampled ids appended, one forward in --dtype,
float32 log_softmax. Δ = vLLM - HF per token and its per-sequence sum. Pass: median |Δ| <= 1e-3 in float32 (0.02
in bfloat16) and (H-A5-1) every item's prompt ids as vLLM echoed them equal the processor's. The p99 per-sequence
|ΣΔ| (pre-registered limit 0.5 nats, archive/phase2_lab/B5/REPORT.md) is report-only since the owner's decision of 2026-11-05
(PLAN D13): bf16 rounding adds up over ~500-token answers (Qwen3-VL-2B R0: token median 4.6e-6, sequence p99 1.79),
and the question it stood for, whether vLLM's answers score like HF's, is measured directly on the trained model
(select_ckpt's HF check, D13).
"""
import argparse
import json
from pathlib import Path

import numpy as np

MEDIAN = {"float32": 1e-3, "bfloat16": 0.02}
SEQ_P99 = 0.5  # report-only (D13, 2026-11-05)


def encode(processor, x: dict):
    """The processor on an exported prompt's text and image: the prompt ids TRL GRPO builds (B5: equal to vLLM's)."""
    from PIL import Image

    image = Image.open(x["image"]).convert("RGB")
    return processor(text=[x["prompt"]], images=[[image]], add_special_tokens=False, return_tensors="pt")


def teacher(model, processor, x: dict, ids: list[int]) -> np.ndarray:
    """HF logprobs (float32 log_softmax) of the sampled ids after prompt `x`."""
    import torch
    enc = encode(processor, x)
    n, comp = enc["input_ids"].shape[1], torch.tensor([ids])
    kw = {}
    for k, v in enc.items():
        if isinstance(v, torch.Tensor) and v.shape == (1, n):  # token-aligned: extend over the completion
            tail = comp if k == "input_ids" else torch.ones_like(comp) if k == "attention_mask" else 0 * comp
            v = torch.cat([v, tail.to(v.dtype)], 1)
        elif isinstance(v, torch.Tensor) and v.is_floating_point():
            v = v.to(model.dtype)
        kw[k] = v.to(model.device) if isinstance(v, torch.Tensor) else v
    with torch.no_grad():
        logits = model(**kw).logits[0, n - 1:-1].float().cpu()
    return torch.log_softmax(logits, -1).gather(1, comp[0][:, None])[:, 0].numpy()


def load(model_dir, prompts: dict, dtype: str):
    """(model, processor) of model_dir in dtype, on CUDA when there is one (call it after the sampler child has
    exited), the processor with the prompts' one processor configuration."""
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor

    kwargs = {json.dumps(x["processor_kwargs"], sort_keys=True) for x in prompts.values()}
    if len(kwargs) != 1:
        raise SystemExit("parity: one processor configuration per prompts file")
    processor = AutoProcessor.from_pretrained(model_dir, **json.loads(kwargs.pop()))
    model = AutoModelForImageTextToText.from_pretrained(model_dir, dtype=getattr(torch, dtype))
    return model.to("cuda" if torch.cuda.is_available() else "cpu").eval(), processor


def compare(model_dir, prompts_path, comps_path, dtype: str, teach=None) -> dict:
    """H-A5-1 on every item (vLLM's echoed prompt_token_ids equal the processor's) and H-A5-2 on the completions of
    the `teach` ids (default: all); "pass" needs both."""
    prompts = {x["id"]: x for x in map(json.loads, Path(prompts_path).read_text().splitlines())}
    comps = [json.loads(x) for x in Path(comps_path).read_text().splitlines()]
    if any("logprobs" not in c for c in comps):
        raise SystemExit("parity: completions without logprobs (run sample_vllm.py --logprobs)")
    model, processor = load(model_dir, prompts, dtype)
    echoed = {c["id"]: c.get("prompt_token_ids") for c in comps}
    bad = [i for i, ids in echoed.items() if encode(processor, prompts[i])["input_ids"][0].tolist() != ids]
    ids = {"items": len(echoed), "equal": len(echoed) - len(bad), "mismatched": bad[:5]}
    ids["pass"] = bool(echoed) and not bad
    deltas = [np.asarray(c["logprobs"], np.float64) - teacher(model, processor, prompts[c["id"]], c["token_ids"])
              for c in comps if teach is None or c["id"] in teach]  # fmt: skip
    t, s = np.abs(np.concatenate(deltas)), np.abs([d.sum() for d in deltas])
    r = {"dtype": dtype, "tokens": int(t.size), "sequences": len(deltas), "median_abs": float(np.median(t)),
         "p99_abs": float(np.quantile(t, 0.99)), "max_abs": float(t.max()), "seq_p99_abs": float(np.quantile(s, 0.99))}
    return r | {"thresholds": {"median_abs": MEDIAN[dtype]}, "report_only": {"seq_p99_abs": SEQ_P99}, "prompt_ids": ids,
                "pass": r["median_abs"] <= MEDIAN[dtype] and ids["pass"]}


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prompts", required=True, help="JSONL from objdel.train.sample.export_prompts")
    ap.add_argument("--completions", required=True, help="JSONL from sample_vllm.py --logprobs")
    ap.add_argument("--model", required=True, help="the checkpoint directory vLLM sampled (with its processor)")
    ap.add_argument("--dtype", choices=sorted(MEDIAN), default="float32")
    ap.add_argument("--out", required=True, help="report JSON")
    a = ap.parse_args(argv)
    r = compare(a.model, a.prompts, a.completions, a.dtype) | {"model": a.model}
    Path(a.out).write_text(json.dumps(r, indent=1) + "\n")
    print(json.dumps(r))
    return r


if __name__ == "__main__":
    main()
