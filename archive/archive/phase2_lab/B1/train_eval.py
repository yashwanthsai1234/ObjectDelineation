"""H-A1-2 / H-A1-8: equal-budget LoRA SFT of one arm on a fixed LIDC train subset, then K samples at temperature 1
on fixed LIDC val items, centred and shifted crops, scored with objdel.eval (the Phase 2 gate definitions).

    HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.train_eval --arm lfm2.5-vl-450m --steps 200 \
        [--seed 0] [--pilot 5]

Training (identical across arms except the model): TRL SFTTrainer through objdel.train.configs.sft_config, rows from
objdel.train.dataset.sft_rows (every voting rater's target) restricted to --items LIDC train items (a seeded subset,
seed 0 for every arm and seed), the configured LoRA adapter (configs.lora: r16, alpha 32, all linear layers incl. the
vision tower), lr 5e-4 constant after 10 warmup steps, micro-batch 1 x 8 accumulation (8 rows per step), fp32 on MPS
(as experiments/phase1/smoke_grpo.py), MPS capped at 7.0 GB. --seed sets the trainer seed (adapter init and data
order).
Evaluation: the adapter is merged; a GRPOTrainer (configs.grpo_config, HF generation: temperature 1, top_p 1,
top_k 0, no repetition penalty) draws K = 8 completions per item (common.sample_k = objdel.train.sample.sample's
loop plus token ids), max_completion_length --max-new (1024: above every LIDC train target, max 881 tokens); scored
per item by objdel.eval.metrics on the item's view (centred or shifted crop, rater masks cut the same way) and
summarised by objdel.eval.summarize. The image-blind prior is scored on the same views. Per-item arrays are written
for the cross-arm bootstrap (analyze.py).
--pilot N: time N training steps only (no eval, nothing kept).
Reads untrusted data: run with python -I.

References:
- TRL 1.16 SFTTrainer (assistant_only_loss, training chat templates): https://huggingface.co/docs/trl/main/en/sft_trainer
- PEFT LoRA merge_and_unload: https://huggingface.co/docs/peft/main/en/package_reference/lora
- torch.mps memory API: https://docs.pytorch.org/docs/stable/mps.html
"""
import argparse
import gc
import json
import time

import numpy as np
import torch
from transformers import AutoModelForImageTextToText, TrainerCallback
from trl import SFTTrainer

from objdel.eval import metrics, summarize
from objdel.train import configs
from objdel.train import dataset as D

from . import common as C


class Clock(TrainerCallback):
    def __init__(self):
        self.t, self.steps, self.peak = None, [], 0.0

    def on_step_begin(self, args, state, control, **kw):
        self.t = time.time()

    def on_step_end(self, args, state, control, **kw):
        self.steps.append(time.time() - self.t)
        self.peak = max(self.peak, C.mps_gb())


def train(a, proc, out):
    recs = D.records(C.PREPARED, "lidc", "train")
    keys = {recs[i]["item_id"] for i in np.random.default_rng(0).permutation(len(recs))[:a.items]}
    rows = D.sft_rows(C.PREPARED, "lidc", "train", proc.tokenizer, configs.MAX_COMPLETION.get(C.ARMS[a.arm]["cfg"]))
    rows = rows.filter(lambda k: k in keys, input_columns="key")
    lens = [len(proc.tokenizer(m[1]["content"][0]["text"], add_special_tokens=False)["input_ids"])
            for m in rows["messages"]]
    model = AutoModelForImageTextToText.from_pretrained(C.ARMS[a.arm]["hf"], dtype=getattr(torch, a.base_dtype))
    steps = a.pilot or a.steps
    args = configs.sft_config(C.ARMS[a.arm]["cfg"], str(out / "sft"), max_steps=steps, per_device_train_batch_size=a.mb,
                              gradient_accumulation_steps=8 // a.mb, learning_rate=5e-4, warmup_steps=10, bf16=bool(a.bf16),
                              logging_steps=10, save_strategy="no", seed=a.seed, data_seed=0,
                              gradient_checkpointing=a.gc)
    clock = Clock()
    trainer = SFTTrainer(model=model, args=args, train_dataset=rows, processing_class=proc, peft_config=C.lora(a.lora_scope),
                         callbacks=[clock])
    t = time.time()
    trainer.train()
    s = clock.steps[1:] or clock.steps
    info = {"rows": len(rows), "items": len(keys), "target_tokens_median": float(np.median(lens)),
            "target_tokens_max": int(max(lens)), "steps": steps, "wall_s": round(time.time() - t, 1),
            "s_per_step": round(float(np.median(s)), 2), "s_per_step_mean": round(float(np.mean(s)), 2),
            "peak_mps_gb": round(clock.peak, 2),
            "loss": [[h["step"], round(h["loss"], 4)] for h in trainer.state.log_history if "loss" in h],
            "trainable_params": sum(p.numel() for p in trainer.model.parameters() if p.requires_grad)}
    print("TRAIN", json.dumps({k: v for k, v in info.items() if k != "loss"}), flush=True)
    if a.pilot:
        return None, info
    trainer.model.save_pretrained(out / "adapter")
    merged = trainer.model.merge_and_unload()
    del trainer
    gc.collect()
    torch.mps.empty_cache()
    return merged, info


def evaluate(a, model, proc, out):
    recs = C.eval_records("lidc", a.eval_items)
    V = C.views("lidc", recs)
    zero = C.zero_led_ids(proc.tokenizer)
    trainer = C.eval_trainer(model, proc, a.arm, str(out / "gen"), a.max_new)
    res, per_item, raw = {}, {}, {}
    for view, items in V.items():
        rows = [{"prompt": [D._user("lidc")], "image": x["image"], "key": x["item_id"]} for x in items]
        t = time.time()
        torch.manual_seed(a.seed)
        gen = C.sample_k(trainer, rows, a.k, zero)
        gen_s = time.time() - t
        its = []
        for x, g in zip(items, gen, strict=True):
            its.append({"item_id": x["item_id"], "dataset": "lidc",
                        **metrics(g["samples"], {n: c.raters for n, c in x["caches"].items()}, x["W"], x["H"],
                                  x["names"], g["truncated"], x["caches"])})
        s = summarize(its)["datasets"]["lidc"]
        tok = np.concatenate([g["tokens"] for g in gen])
        zl = np.concatenate([g["zero_led"] for g in gen])
        res[view] = {"summary": s, "gen_s": round(gen_s, 1), "gen_s_per_item": round(gen_s / len(items), 2),
                     "tokens_mean": round(float(tok.mean()), 1), "tokens_median": float(np.median(tok)),
                     "truncated": float(np.mean([t for g in gen for t in g["truncated"]])),
                     "zero_led_samples": float((zl > 0).mean()), "zero_led_tokens": int(zl.sum())}
        per_item[view] = {
            "item_id": [x["item_id"] for x in items],
            "iou_tc": [it["objects"]["nodule"]["iou_tc"] for it in its],
            "format": [it["samples"]["format_ok"] for it in its], "legal": [it["samples"]["legal"] for it in its],
            "dup": [it["groups"]["dup"][0] for it in its], "signal": [float(it["groups"]["signal"][0]) for it in its],
            "tokens": [g["tokens"] for g in gen], "zero_led": [g["zero_led"] for g in gen],
            "error": [it["samples"]["error"] for it in its],
            "prior": [C.prior_iou("lidc", x)["nodule"] for x in items]}
        raw[view] = gen
        print("EVAL", view, json.dumps({"format": s["format"], "legal": s["legal"],
                                        "iou_tc_median": s["objects"]["nodule"]["iou_tc_median"],
                                        "dup": s["dup_rate"], "signal": s["signal_share"],
                                        **{k: v for k, v in res[view].items() if k != "summary"}}), flush=True)
    res["peak_mps_gb"] = round(C.mps_gb(), 2)
    return res, per_item, raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=list(C.ARMS))
    ap.add_argument("--steps", type=int, default=36)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--items", type=int, default=1000, help="LIDC train items (all voting raters' targets)")
    ap.add_argument("--eval-items", type=int, default=64)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--max-new", type=int, default=1024)
    ap.add_argument("--pilot", type=int, default=0)
    ap.add_argument("--gc", type=int, default=1, help="gradient checkpointing (1/0)")
    ap.add_argument("--mb", type=int, default=1, help="micro-batch (8 rows per step in total)")
    ap.add_argument("--bf16", type=int, default=0, help="bf16 autocast in training (1/0)")
    ap.add_argument("--base-dtype", default="float32", help="frozen base weights (LoRA adapters stay fp32)")
    ap.add_argument("--lora-scope", default="all", choices=("all", "lm"))
    ap.add_argument("--base", action="store_true", help="evaluate the untrained base model (no SFT)")
    a = ap.parse_args()
    C.cap_mps()
    run = f"{a.arm}_s{a.seed}" + (f"_n{a.steps}" if a.steps != 36 else "") + ("_lm" if a.lora_scope == "lm" else "") + (f"_pilot_mb{a.mb}_bf{a.bf16}" if a.pilot else "") + ("_base" if a.base else "")
    out = C.SCRATCH / "runs" / run
    out.mkdir(parents=True, exist_ok=True)
    proc = C.processor(a.arm, "lidc")
    t0 = time.time()
    if a.base:
        model, info = AutoModelForImageTextToText.from_pretrained(C.ARMS[a.arm]["hf"], dtype=torch.float32), {}
    else:
        model, info = train(a, proc, out)
    res = {"args": vars(a), "train": info}
    if model is not None:
        model.to("mps")
        res["eval"], per_item, raw = evaluate(a, model, proc, out)
        (out / "per_item.json").write_text(json.dumps(per_item))
        (out / "samples.json").write_text(json.dumps(raw))
    res["wall_s"] = round(time.time() - t0, 1)
    (out / "result.json").write_text(json.dumps(res, indent=1, default=str))
    print("DONE", run, res["wall_s"], flush=True)


if __name__ == "__main__":
    main()
