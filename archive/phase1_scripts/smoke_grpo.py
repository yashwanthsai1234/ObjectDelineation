"""K12 with reward v4 on the Mac: a short SFT of SmolVLM-256M, then a few TRL GRPO steps with the production reward.

    .venv/bin/python -I -m experiments.phase1.smoke_grpo [--prepared data/prepared] [--out data/runs/smoke_grpo]

1. SFT (TRL SFTTrainer, objdel.train.configs.sft_config) on --sft-items LIDC train items, every rater's target, for
   --sft-steps steps: deliberately short, because an over-trained SFT collapses sampling diversity and gives GRPO no
   signal (lab t_stack K12). LoRA (configs.lora's adapter, r16 on all linear layers; merged afterwards) instead of
   the configured full fine-tune, to stay within the Mac's memory budget.
2. GRPO (TRL GRPOTrainer, grpo_config: dr_grpo, beta 1e-3, G = 8, max_completion_length 1,536, HF generation) on
   --grpo-items other LIDC train items (disagreement items oversampled), --grpo-steps steps of 2 prompts x 8
   completions (micro-batch 1), a LoRA adapter (the reference model is the adapter-disabled SFT model), TRLReward on
   the prepared caches. No held-out evaluation: TRL computes evaluation log-probabilities for a whole group (8
   completions) at once, which alone exceeds the memory budget here; that path (eval_rewards/*, Controller) is covered
   by tests/test_train_trl.py and was run once uncapped (experiments/phase1/impl/M7_train.md).
MPS memory is capped with --mps-fraction (default 0.35 of recommended_max_memory, 4.4 GB on a 16 GB M3).
Asserts: every per-term log present, no NaN, non-zero advantages in at least one group. Records timings, memory,
per-term means and, per generation batch, the plan §7 Phase 2 GRPO-signal gate numbers with the gates' own
definitions (objdel.eval.metrics / summarize, read by eval.gates grpo_duplicates and grpo_sigma): duplicate_rate
(equal rasters within a group of 8) and sigma_above_floor_any_term (share of groups with a legality-gated
geometric sigma above its floor), in <out>/smoke.json.
Reads untrusted data: run with python -I.

References:
- TRL (version pinned in pyproject.toml) SFTTrainer / GRPOTrainer (peft_config, reference model = adapter-disabled policy under PEFT):
  https://huggingface.co/docs/trl/main/en/sft_trainer , https://huggingface.co/docs/trl/main/en/grpo_trainer ,
  https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py
- PEFT LoRA (target_modules="all-linear", merge_and_unload): https://huggingface.co/docs/peft/main/en/package_reference/lora
- torch.mps.set_per_process_memory_fraction (fraction of recommendedMaxWorkingSetSize):
  https://docs.pytorch.org/docs/stable/generated/torch.mps.set_per_process_memory_fraction.html
- DAPO (groups with identical rewards carry no gradient; dynamic sampling): https://arxiv.org/abs/2503.14476
- Plan §7 (Phase 2 GRPO-signal gate) and lab t_stack K12: experiments/phase1/PLAN.md,
  experiments/phase1/lab/t_stack/REPORT.md
"""
import argparse
import gc
import json
import math
import resource
import time
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModelForImageTextToText, AutoProcessor
from trl import GRPOTrainer, SFTTrainer

from objdel.eval import metrics, summarize
from objdel.reward import REWARD_V4
from objdel.train import configs
from objdel.train import dataset as D
from objdel.train.rewards import AGGREGATOR, LOGGED, TRLReward, stop_ids

MODEL, DATASET, G = "smolvlm-256m", "lidc", REWARD_V4.group
LORA = configs.lora("qwen3-vl-2b")  # the configured adapter (r16, all linear layers), here for SmolVLM


def memory():
    mps = torch.mps.driver_allocated_memory() / 1e9 if torch.backends.mps.is_available() else 0.0
    return {"peak_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2), "mps_driver_gb": round(mps, 2)}


def group_stats(records, caches, adv):
    """One generation batch: the Phase 2 GRPO-signal gate numbers as objdel.eval computes them (eval.metrics on each
    group of G, eval.summarize): duplicate_rate = dup_rate (equal rasters within a group), sigma_above_floor_any_term
    = signal_share (groups with a legality-gated geometric sigma above its floor); and the groups with a non-zero
    advantage. records: TRLReward.recorded, (key, text, truncated, W, H, names) per completion, in order."""
    items = []
    for i in range(0, len(records), G):
        block = records[i:i + G]
        key, _, _, W, H, names = block[0]
        c = caches(key)
        items.append({"item_id": key, "dataset": DATASET, **metrics([r[1] for r in block], {n: c[n].raters for n in names},
                                                                    W, H, names, [r[2] for r in block], c)})
    s = summarize(items)["datasets"][DATASET]
    a = np.asarray(adv).reshape(-1, G)
    return {"duplicate_rate": s["dup_rate"], "groups": len(a), "sigma_above_floor_any_term": s["signal_share"],
            "nonzero_adv_groups": int((np.abs(a) > 0).any(1).sum()), "adv_abs_max": float(np.abs(a).max())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepared", default="data/prepared")
    ap.add_argument("--out", default="data/runs/smoke_grpo")
    ap.add_argument("--sft-items", type=int, default=48)
    ap.add_argument("--sft-steps", type=int, default=40)
    ap.add_argument("--grpo-items", type=int, default=8)
    ap.add_argument("--grpo-steps", type=int, default=3)
    ap.add_argument("--mps-fraction", type=float, default=0.35, help="MPS cap, share of recommended_max_memory")
    a = ap.parse_args()
    if torch.backends.mps.is_available():
        torch.mps.set_per_process_memory_fraction(a.mps_fraction)  # of recommended_max_memory (12.7 GB on a 16 GB M3)
    out, t0, res = Path(a.out), time.time(), {"args": vars(a)}
    mid = configs.MODELS[MODEL]
    proc = AutoProcessor.from_pretrained(mid, **configs.processor_kwargs(MODEL, DATASET))
    recs = [r for r in D.records(a.prepared, DATASET, "train") if r["rl"]]
    rng = np.random.default_rng(0)
    order = rng.permutation(len(recs))
    sft_keys = {recs[i]["item_id"] for i in order[:a.sft_items]}
    rest = [recs[i] for i in order[a.sft_items:]]
    grpo_keys = {r["item_id"] for r in sorted(rest, key=lambda r: not r["disagreement"])[:a.grpo_items]}

    # 1. SFT
    sft = D.sft_rows(a.prepared, DATASET, "train", proc.tokenizer, configs.MAX_COMPLETION[MODEL])
    sft = sft.filter(lambda k: k in sft_keys, input_columns="key")
    model = AutoModelForImageTextToText.from_pretrained(mid, dtype=torch.float32)
    args = configs.sft_config(MODEL, str(out / "sft"), max_steps=a.sft_steps, per_device_train_batch_size=1,
                              gradient_accumulation_steps=4,
                              learning_rate=5e-4, bf16=False, logging_steps=5, save_strategy="no")
    trainer = SFTTrainer(model=model, args=args, train_dataset=sft, processing_class=proc, peft_config=LORA)
    t = time.time()
    trainer.train()
    res["sft"] = {"rows": len(sft), "s": round(time.time() - t, 1), "s_per_step": round((time.time() - t) / a.sft_steps, 2),
                  "loss": [round(h["loss"], 4) for h in trainer.state.log_history if "loss" in h], **memory()}
    trainer.model.merge_and_unload().save_pretrained(out / "sft")
    print("SFT", json.dumps(res["sft"]), flush=True)
    del trainer, model
    gc.collect()  # the SFT trainer holds reference cycles: free its tensors before the GRPO model loads
    if torch.backends.mps.is_available():
        torch.mps.empty_cache()

    # 2. GRPO
    train = D.grpo_rows(a.prepared, DATASET, "train", D.DISAGREEMENT_OVERSAMPLE).filter(lambda k: k in grpo_keys, input_columns="key")
    model = AutoModelForImageTextToText.from_pretrained(out / "sft", dtype=torch.float32)
    caches = D.cache_loader(a.prepared)
    reward = TRLReward(caches, stop_ids(model, proc))
    reward.recorded = []  # the completions of each reward call, for group_stats
    funcs, weights = reward.funcs()
    args = configs.grpo_config(MODEL, str(out / "grpo"), max_steps=a.grpo_steps, per_device_train_batch_size=1,
                               gradient_accumulation_steps=2 * G, learning_rate=1e-5, bf16=False, logging_steps=1, save_strategy="no",
                               gradient_checkpointing=True)
    trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=args, train_dataset=train,
                          processing_class=proc, peft_config=LORA)
    batches, orig = [], trainer._generate_and_score_completions

    def hook(inputs):
        tg = time.time()
        o = orig(inputs)
        n = len(inputs)
        records, reward.recorded[:] = list(reward.recorded), []
        terms = {k: list(trainer._logs["rewards"][k])[-n:] for k in LOGGED}
        batches.append({"generate_and_score_s": round(time.time() - tg, 1),
                        **group_stats(records, caches, o["advantages"].cpu().numpy()),
                        **{f"mean_{k}": float(np.mean(v)) for k, v in terms.items()}, **memory()})
        print("BATCH", json.dumps(batches[-1]), flush=True)
        return o

    trainer._generate_and_score_completions = hook
    t = time.time()
    trainer.train()
    logs = [h for h in trainer.state.log_history if "loss" in h]
    res["grpo"] = {"rows": len(train), "steps": len(logs), "s": round(time.time() - t, 1),
                   "s_per_step": round((time.time() - t) / max(1, len(logs)), 1), "batches": batches,
                   "logs": logs, **memory()}
    res["wall_s"] = round(time.time() - t0, 1)
    keys = {f"rewards/{n}/mean" for n in (*LOGGED, AGGREGATOR)}
    res["checks"] = {
        "per_term_logs": all(keys <= set(h) for h in logs),
        "no_nan": not any(isinstance(v, float) and math.isnan(v) for h in trainer.state.log_history for v in h.values()),
        "nonzero_advantage_groups": sum(b["nonzero_adv_groups"] for b in batches)}
    (out / "smoke.json").write_text(json.dumps(res, indent=1, default=str))
    print("RESULT", json.dumps({k: res[k] for k in ("wall_s", "checks")}), flush=True)
    assert res["checks"]["per_term_logs"] and res["checks"]["no_nan"] and res["checks"]["nonzero_advantage_groups"] > 0


if __name__ == "__main__":
    main()
