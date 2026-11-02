"""The HF generate sampler (PLAN D13 fallback, Qwen3.5 HF-only path, CPU smoke): sample_vllm's JSONL contract.

    python -I -m experiments.phase2.hf_sample --prompts P.jsonl --model DIR --out C.jsonl [--batch_items 16]
        [--timing T.json]

A child of select_ckpt in its own process (select_ckpt never touches CUDA). Every D13 sampling field is passed
explicitly (configs.as_names(sampling, "hf"), eos_token_id = the stop ids) over a fresh GenerationConfig, so nothing
comes from the checkpoint's (LFM's card: min_p 0.15, repetition 1.05); HF has no presence / frequency penalty
(as_names refuses one). Output per sample: {"id", "k", "token_ids" (cut after the row's own first stop id),
"finish_reason", "text"}.

Batches: --batch_items prompts per generate call (each expanded to its n samples, num_return_sequences), padded on
the left (the processor is loaded with padding_side="left"; each item's text carries its own image tokens and the
models mask the padding, LFM2's short convolutions and Qwen3.5's linear attention included), so a 450M model is not
overhead-bound at one item per call. All prompts must share sampling, n and processor_kwargs.

Seeding and reproducibility: batch b (prompts b*B .. b*B+B-1 of the file, B = --batch_items) is seeded seed + b*B.
The same prompts file, --batch_items, model and hardware / library versions give the same samples; another
--batch_items gives other (equally valid) draws. Resume: --out is kept: lines of whole finished batches stay, a
half-written batch or a cut last line is dropped, and the remaining batches are sampled with their own seeds, so a
resumed file equals an uninterrupted one. --timing writes generate seconds (model load apart) per item-sample.

Reference: https://huggingface.co/docs/transformers/main_classes/text_generation
"""

import argparse
import json
import sys
import time
from pathlib import Path


def load(path):
    """(model, processor, device): library generation defaults, never the card's; left padding for batches. The GPU
    samples in bf16 whatever the checkpoint holds (full FT saves fp32; from_pretrained's "auto" would keep it), the
    dtype of D13's vLLM path and of R0's parity and rate; CPU smoke runs use fp32."""
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor, GenerationConfig

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = torch.bfloat16 if dev == "cuda" else torch.float32
    model = AutoModelForImageTextToText.from_pretrained(path, dtype=dtype).to(dev).eval()
    model.generation_config = GenerationConfig()
    proc = AutoProcessor.from_pretrained(path, padding_side="left")
    if proc.tokenizer.padding_side != "left":
        raise SystemExit(f"hf_sample: {path}: the tokenizer pads on the {proc.tokenizer.padding_side}")
    return model, proc, dev


def records(rows: list[list[int]], prompts: list[dict], stop, decode) -> list[dict]:
    """Generated ids (n rows per prompt, in order) -> output lines, each row cut after its own first stop id."""
    out, n = [], prompts[0]["n"]
    for j, row in enumerate(rows):
        cut = next((i + 1 for i, t in enumerate(row) if t in stop), None)
        toks = row[:cut] if cut else row
        out.append({"id": prompts[j // n]["id"], "k": j % n, "token_ids": toks,
                    "finish_reason": "stop" if cut else "length", "text": decode(toks)})  # fmt: skip
    return out


def generate(model, proc, dev, batch: list[tuple[int, dict]]) -> list[dict]:
    """One generate call on a batch of (file index, prompt), seeded by its first index."""
    import torch
    from PIL import Image

    from objdel.train import configs

    ps = [p for _, p in batch]
    s, kw = ps[0]["sampling"], {k: v for k, v in ps[0]["processor_kwargs"].items() if k != "padding_side"}
    gen = configs.as_names(s, "hf")
    x = proc(text=[p["prompt"] for p in ps], images=[[Image.open(p["image"]).convert("RGB")] for p in ps],
             padding=True, return_tensors="pt", **kw)  # fmt: skip
    torch.manual_seed(s["seed"] + batch[0][0])
    with torch.no_grad():
        ids = model.generate(**x.to(dev), **gen, num_return_sequences=ps[0]["n"],
                             pad_token_id=proc.tokenizer.pad_token_id)  # fmt: skip
    rows = ids[:, x["input_ids"].shape[1] :].tolist()
    return records(rows, ps, set(gen["eos_token_id"]), lambda t: proc.tokenizer.decode(t, skip_special_tokens=True))


def todo(out: Path, batches: list[list[tuple[int, dict]]]) -> list:
    """Resume: rewrite `out` with the lines of every whole finished batch and return the other batches."""
    have = {}
    for x in out.read_text().splitlines() if out.exists() else []:
        try:
            c = json.loads(x)
            have.setdefault((c["id"], c["k"]), x)
        except (ValueError, KeyError, TypeError):
            continue  # a line cut by a kill
    keys = [[(p["id"], k) for _, p in b for k in range(p["n"])] for b in batches]
    done = [all(key in have for key in ks) for ks in keys]
    out.write_text("".join(have[key] + "\n" for ks, d in zip(keys, done, strict=True) if d for key in ks))
    return [b for b, d in zip(batches, done, strict=True) if not d]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--prompts", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch_items", type=int, default=16, help="prompts per generate call (x n samples each)")
    ap.add_argument("--timing", type=Path, help="write generate seconds per item-sample here (R0 rate, hf_rate)")
    a = ap.parse_args(argv)
    ps = [json.loads(x) for x in a.prompts.read_text().splitlines()]
    if len({json.dumps([p["sampling"], p["n"], p["processor_kwargs"]], sort_keys=True) for p in ps}) > 1:
        raise SystemExit("hf_sample: prompts differ in sampling, n or processor_kwargs; a batch needs one setting")
    indexed = list(enumerate(ps))
    left = todo(a.out, [indexed[i : i + a.batch_items] for i in range(0, len(ps), a.batch_items)])
    t0 = time.perf_counter()
    model, proc, dev = load(a.model) if left else (None, None, None)
    t1 = time.perf_counter()
    lines = []
    with open(a.out, "a") as f:
        for b in left:
            recs = generate(model, proc, dev, b)
            f.write("".join(json.dumps(c) + "\n" for c in recs))  # whole batches: a kill cuts at most one line
            f.flush()
            lines += recs
    if a.timing:
        gen_s, n = time.perf_counter() - t1, len(lines)
        a.timing.write_text(json.dumps({
            "items": sum(map(len, left)), "samples": n, "batch_items": a.batch_items, "load_s": t1 - t0,
            "generate_s": gen_s, "s_per_sample": gen_s / n if n else None,
            "tokens_mean": sum(len(c["token_ids"]) for c in lines) / n if n else None,
            "length_share": sum(c["finish_reason"] == "length" for c in lines) / n if n else None}) + "\n")  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
