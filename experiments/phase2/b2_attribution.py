"""B2.3b attribution of the trainer-vLLM gap, Mac part ($0, CPU, no training):

    python -I -m experiments.phase2.b2_attribution --out RD/b2_attribution.json

Scores the same answers four ways under torch.no_grad(), each LoRA adapter on the run's start checkpoint:
- (a) the FastGRPO default: packed shared-prefix forward, VisionCache, RopeCache, bf16 autocast;
- (b) TRL's plain forward (no FastGRPO patch), bf16 autocast;
- (c) TRL's plain forward of a model holding fast.lora_merged_params' values (W + delta rounded in W's dtype, what
  vLLM receives), no adapter, bf16 autocast;
- (d) TRL's plain forward in fp32 (base cast to fp32, no autocast): the reference.
Answers: per prompt, its first K saved answers with token ids (A3's when they exist; by default the e2 completions
of the R1 selection), for the first N prompts in file order.
Prints, per adapter and pair, the signed mean logprob gap, the mean |gap| and both per 64-token answer position bin.
Readings (Implementation.md B2.3b): (a)-(b) is the packed forward's part, (c)-(b) merge rounding's part, (b)-(d)
the trainer's bf16 part. The box part (vLLM lines (e), (f)) is not here.
"""

import argparse
import contextlib
import gc
import json
import sys
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import torch
from trl import GRPOTrainer

from objdel.train import fast

PAIRS = (("a", "d"), ("b", "d"), ("c", "d"), ("a", "b"), ("c", "b"))
BIN = 64
RUN = Path("runs/rl/train_rl/trainer")
ANSWERS = Path("runs/r1/r1-qwen3vl2b/select/e2/shifted-k16")


def forward_obj(cls, bf16: bool, **attrs):
    """A trainer of class `cls` with only what its logprob forward reads: TRL's autocast (bf16 as Accelerate runs a
    bf16 run, else none), the metrics dict and the model's mode (packed_chunk logs a token ratio)."""
    t = object.__new__(cls)
    cast = (lambda: torch.autocast("cpu", dtype=torch.bfloat16)) if bf16 else contextlib.nullcontext
    t.accelerator = SimpleNamespace(autocast=cast, is_main_process=True)  # TRL's profiling decorator reads it
    t._metrics = defaultdict(lambda: defaultdict(list))
    t.model = SimpleNamespace(training=False)
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


@torch.no_grad()
def score(t, model, batches, clear=()) -> list[torch.Tensor]:
    """Per-token answer logprobs (float32) of each batch (ids, mask, T, kw); `clear` caches are emptied per batch."""
    out = []
    for ids, mask, T, kw in batches:
        for c in clear:
            c.clear()
        out.append(t._get_per_token_logps_and_entropies(model, ids, mask, T, batch_size=len(ids), **dict(kw))[0].float())
    return out


def score_ways(peft, batches, load_base) -> dict:
    """{way: [logprobs per batch]}. peft: a PEFT model (base in its checkpoint dtype, fp32 adapter) with TRL's fused
    LM head; load_base(): a fresh base model with the fused head, for (c). peft's weights are as before afterwards;
    (a)'s caches stay installed on it."""
    plain = forward_obj(GRPOTrainer, bf16=True)
    lp = {"b": score(plain, peft, batches)}
    merged = load_base()
    missing, unexpected = merged.load_state_dict(dict(fast.lora_merged_params(peft)), strict=False)
    if unexpected or any(not k.endswith("lm_head.weight") for k in missing):  # a tied head is not a named parameter
        raise RuntimeError(f"merged weights do not fit the base: missing {missing[:3]}, unexpected {unexpected[:3]}")
    lp["c"] = score(plain, merged.eval(), batches)
    del merged
    gc.collect()
    low = {n: p.dtype for n, p in peft.named_parameters() if p.dtype != torch.float32}
    for p in peft.parameters():
        p.data = p.data.float()
    lp["d"] = score(forward_obj(GRPOTrainer, bf16=False), peft, batches)
    for n, p in peft.named_parameters():
        if n in low:
            p.data = p.data.to(low[n])  # exact: the values came from that dtype
    caches = [c for c in (fast.VisionCache.install(peft), fast.RopeCache.install(peft)) if c is not None]
    t = forward_obj(fast.FastGRPO, bf16=True, shared_prefix=True, loss_inputs=None, rope=caches[-1] if caches else None)
    lp["a"] = score(t, peft, batches, clear=caches)
    return lp


def report(lp: dict, masks: list, pairs=PAIRS) -> dict:
    """{"x-y": {signed, abs, bins: [{signed, abs, n} per BIN answer positions]}} of lp[x] - lp[y] over the real
    answer tokens (mask: each batch's answer mask, or its full attention mask, cut to the answer width)."""
    out = {}
    for x, y in pairs:
        s = defaultdict(float)
        for a, b, m in zip(lp[x], lp[y], masks, strict=True):
            m = m[:, -a.size(1):].bool()
            g = (a - b).double()
            pos = torch.arange(a.size(1)).expand_as(a)[m] // BIN
            for k, v in (("signed", g[m]), ("abs", g[m].abs()), ("n", torch.ones_like(g[m]))):
                for i, t in enumerate(torch.zeros(int(pos.max()) + 1 if len(pos) else 0, dtype=torch.float64)
                                      .index_add_(0, pos, v).tolist()):  # fmt: skip
                    s[(k, i)] += t
        n_bins = 1 + max((i for _, i in s), default=-1)
        bins = [{"signed": s[("signed", i)] / s[("n", i)], "abs": s[("abs", i)] / s[("n", i)], "n": int(s[("n", i)])}
                for i in range(n_bins) if s[("n", i)]]  # fmt: skip
        n = sum(b["n"] for b in bins)
        out[f"{x}-{y}"] = {"signed": sum(b["signed"] * b["n"] for b in bins) / n,
                           "abs": sum(b["abs"] * b["n"] for b in bins) / n, "bins": bins}  # fmt: skip
    return out


def pick(answers: Path, n_prompts: int, k: int) -> list[tuple[dict, list[list[int]]]]:
    """(prompt row, its first k answers' token ids) for the first n_prompts prompts (file order) with k answers."""
    prompts = {r["id"]: r for r in map(json.loads, (answers / "prompts.jsonl").read_text().splitlines())}
    ids = defaultdict(dict)
    for r in map(json.loads, (answers / "completions.jsonl").read_text().splitlines()):
        if r["k"] < k:
            ids[r["id"]][r["k"]] = r["token_ids"]
    full = [(prompts[i], [a[j] for j in range(k)]) for i, a in ids.items() if len(a) == k]
    return full[:n_prompts]


def batch(proc, prompt: dict, answers: list[list[int]], views: Path, pad: int):
    """TRL's layout for one group: the processed prompt (left-padded, one image) repeated per answer, answers
    right-padded with the pad id; (ids, mask, T, kw)."""
    from PIL import Image

    image = views / prompt["image"].split("eval_views/", 1)[1]  # the box path -> this machine's eval view
    p = proc(text=[prompt["prompt"]] * len(answers), images=[Image.open(image).convert("RGB")] * len(answers),
             return_tensors="pt", padding=True)  # fmt: skip
    T = max(map(len, answers))
    comp = torch.full((len(answers), T), pad)
    cmask = torch.zeros(len(answers), T, dtype=torch.long)
    for r, a in enumerate(answers):
        comp[r, : len(a)], cmask[r, : len(a)] = torch.tensor(a), 1
    kw = {"pixel_values": p["pixel_values"], "image_grid_thw": p["image_grid_thw"], "num_images": [1] * len(answers),
          "mm_token_type_ids": torch.cat([p["mm_token_type_ids"], torch.zeros_like(comp)], 1)}  # fmt: skip
    return torch.cat([p["input_ids"], comp], 1), torch.cat([p["attention_mask"], cmask], 1), T, kw


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--base", type=Path, default=Path("runs/r1/r1-qwen3vl2b/branches/e2"),
                    help="the run's start checkpoint (the full run started from the selected R1 branch e2)")
    ap.add_argument("--adapters", type=Path, nargs="+", default=[RUN / "checkpoint-100", RUN / "checkpoint-770"])
    ap.add_argument("--answers", type=Path, default=ANSWERS, help="prompts.jsonl + completions.jsonl with token_ids")
    ap.add_argument("--views", type=Path, default=Path("data/eval_views"), help="eval views root (images)")
    ap.add_argument("--prompts", type=int, default=32)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    from peft import PeftModel
    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl.trainer.utils import add_fused_lm_head

    from objdel.train import configs

    def load_base():
        m = AutoModelForImageTextToText.from_pretrained(a.base).eval()
        add_fused_lm_head(m, temperature=1.0, outputs=("log_probs", "entropy"))  # as GRPOTrainer (T = 1)
        return m

    proc = AutoProcessor.from_pretrained(a.base, **configs.processor_kwargs("qwen3-vl-2b", "lidc"))
    rows = pick(a.answers, a.prompts, a.k)
    batches = [batch(proc, p, ans, a.views, proc.tokenizer.pad_token_id) for p, ans in rows]
    masks = [m for _, m, _, _ in batches]
    doc = {"base": str(a.base), "answers": str(a.answers), "sequences": len(rows) * a.k,
           "tokens": int(sum(m[:, -T:].sum() for _, m, T, _ in batches)), "adapters": {}}  # fmt: skip
    for path in a.adapters:
        peft = PeftModel.from_pretrained(load_base(), path).eval()
        doc["adapters"][str(path)] = report(score_ways(peft, batches, load_base), masks)
        del peft
        gc.collect()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"{doc['sequences']} answers, {doc['tokens']} tokens; gap = way minus way (nats per token)")
    for path, pairs in doc["adapters"].items():
        print(path)
        for name, r in pairs.items():
            bins = " ".join(f"{b['abs']:.4f}" for b in r["bins"])
            print(f"  {name}  signed {r['signed']:+.5f}  |gap| {r['abs']:.5f}  |gap| per {BIN} tokens: {bins}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
