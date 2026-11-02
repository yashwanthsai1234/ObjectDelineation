"""The HF generate sampler (PLAN D13 fallback, Qwen3.5 HF-only path, CPU smoke): sample_vllm's JSONL contract.

    python -I -m experiments.phase2.hf_sample --prompts P.jsonl --model DIR --out C.jsonl

A child of select_ckpt in its own process (select_ckpt never touches CUDA). Every D13 sampling field is passed
explicitly (configs.as_names(sampling, "hf"), eos_token_id = the stop ids) over a fresh GenerationConfig, so nothing
comes from the checkpoint's (LFM's card: min_p 0.15, repetition 1.05); HF has no presence / frequency penalty
(as_names refuses one). Prompt i is seeded seed + i. Output per sample: {"id", "k", "token_ids" (cut after the
first stop id), "finish_reason", "text"}.

Reference: https://huggingface.co/docs/transformers/main_classes/text_generation
"""

import argparse
import json
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--prompts", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    import torch
    from PIL import Image
    from transformers import AutoModelForImageTextToText, AutoProcessor, GenerationConfig

    from objdel.train import configs

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = AutoModelForImageTextToText.from_pretrained(a.model).to(dev).eval()
    model.generation_config = GenerationConfig()  # library defaults, never the card's
    proc = AutoProcessor.from_pretrained(a.model)
    with open(a.out, "w") as f:
        for i, line in enumerate(a.prompts.read_text().splitlines()):
            p = json.loads(line)
            s, kw = p["sampling"], {k: v for k, v in p["processor_kwargs"].items() if k != "padding_side"}
            gen = configs.as_names(s, "hf")
            x = proc(text=[p["prompt"]], images=[[Image.open(p["image"]).convert("RGB")]], return_tensors="pt", **kw)
            torch.manual_seed(s["seed"] + i)
            with torch.no_grad():
                ids = model.generate(**x.to(dev), **gen, num_return_sequences=p["n"],
                                     pad_token_id=proc.tokenizer.pad_token_id)  # fmt: skip
            for k, row in enumerate(ids[:, x["input_ids"].shape[1] :].tolist()):
                cut = next((j + 1 for j, t in enumerate(row) if t in gen["eos_token_id"]), None)
                toks = row[:cut] if cut else row
                out = {"id": p["id"], "k": k, "token_ids": toks, "finish_reason": "stop" if cut else "length"}
                f.write(json.dumps(out | {"text": proc.tokenizer.decode(toks, skip_special_tokens=True)}) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
