"""K12 redo: real GRPO steps with the production stack (strict parser + polycore + reward v2 + TRL aggregator).

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I k12_grpo.py --temperature 1.0 --tag t10

SmolVLM-256M from runs/sft_smolvlm-256m_full (trained on {"disc":[{"outer":..,"holes":[]}]} = production grammar
with names ("disc",)), the same 4 RIGA crops (448 px, jittered), rater set = all experts' disc masks at the crop.
No vLLM. max_completion_length = 1536 (D10). G = 8, 2 prompts per generation (2 groups -> batch z over 16).
"""
import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import torch  # noqa: E402
from datasets import Dataset  # noqa: E402

import prod  # noqa: E402
from common import HERE, PROMPT_POLY, model_path, mps_mem_gb, peak_rss_gb, riga_samples  # noqa: E402

NAMES = ("disc",)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", default=str(HERE / "runs" / "sft_smolvlm-256m_full"))
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--temperature", type=float, default=1.0)
    ap.add_argument("--max-completion", type=int, default=1536)
    ap.add_argument("--tag", default="t10")
    ap.add_argument("--peft", type=int, default=1)
    a = ap.parse_args()

    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl import GRPOConfig, GRPOTrainer

    proc = AutoProcessor.from_pretrained(model_path("smolvlm-256m"), padding_side="left", size={"longest_edge": 512})
    model = AutoModelForImageTextToText.from_pretrained(a.init, dtype=torch.float32)
    samples = riga_samples(4, seed=1, side=448)
    caches, sizes = {}, {}
    rows = []
    for s in samples:
        caches[s["key"]] = {"disc": prod.build_cache(np.stack(s["disc_masks"]))}
        sizes[s["key"]] = s["image"].size
        rows.append({"prompt": [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPT_POLY}]}],
                     "image": s["image"], "key": s["key"]})
    ds = Dataset.from_list(rows)
    tok = proc.tokenizer
    eos_ids = {tok.eos_token_id, tok.convert_tokens_to_ids("<end_of_utterance>")}
    errors = []

    def provider(prompts, completions, completion_ids, key, **kw):
        out = {k: [] for k in ("format_ok", "legal", "iou_tc", "f5", "topo_any", "topo_freq", "sig_id")}
        for c, ids, k in zip(completions, completion_ids, key):
            text = c[0]["content"] if isinstance(c, list) else c
            trunc = len(ids) == 0 or ids[-1] not in eos_ids
            W, H = sizes[k]
            r = prod.score_text(text, trunc, NAMES, caches[k], W, H)
            errors.append(r.pop("error"))
            for n in out:
                out[n].append(r[n])
        return out

    agg = prod.TRLReward(provider, group=8, drop_spread=1e-2)
    funcs, weights = agg.funcs()
    peft_config = None
    if a.peft:
        from peft import LoraConfig
        peft_config = LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", task_type="CAUSAL_LM")
    cfg = GRPOConfig(
        output_dir=str(HERE / "runs" / f"k12_{a.tag}"), max_steps=a.steps, per_device_train_batch_size=8,
        gradient_accumulation_steps=2, num_generations=8, max_completion_length=a.max_completion,
        temperature=a.temperature, learning_rate=1e-5, lr_scheduler_type="constant", warmup_steps=0,
        loss_type="dr_grpo", epsilon=0.2, epsilon_high=0.28, beta=1e-3,
        multi_objective_aggregation="sum_then_normalize", scale_rewards="none", reward_weights=weights,
        mask_truncated_completions=True, use_vllm=False, bf16=False, gradient_checkpointing=True,
        logging_steps=1, report_to="none", save_strategy="no", seed=0, remove_unused_columns=False)
    trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=cfg, train_dataset=ds, processing_class=proc,
                          peft_config=peft_config)
    captured = []
    orig = trainer._generate_and_score_completions

    def hook(*x, **k):
        out = orig(*x, **k)
        captured.append(out["advantages"].detach().cpu().double().numpy().copy())
        return out
    trainer._generate_and_score_completions = hook
    t0 = time.time()
    status, err = "pass", None
    try:
        trainer.train()
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        status, err = "fail", f"{type(e).__name__}: {e}"[:600]
    wall = time.time() - t0
    logs = [h for h in trainer.state.log_history if "loss" in h]
    nan = any(isinstance(v, float) and math.isnan(v) for l in logs for v in l.values())
    # equality: trainer advantages vs the aggregator's own reference, per generation
    diffs = [float(np.abs(c - h["advantage"]).max()) for c, h in zip(captured, agg.history)]
    batches = []
    for h in agg.history:
        g = h["iou_tc"].reshape(-1, 8)
        batches.append({
            "format_ok": float(h["format_ok"].mean()), "legal": float(h["legal"].mean()),
            "iou_tc_mean": float(h["iou_tc"].mean()), "iou_tc_group_std": g.std(1).round(4).tolist(),
            "iou_tc_range": [float(h["iou_tc"].min()), float(h["iou_tc"].max())],
            "f5_mean": float(h["f5"].mean()), "f5_range": [float(h["f5"].min()), float(h["f5"].max())],
            "topo_any": float(h["topo_any"].mean()), "distinct_sigs": len(set(h["sig_id"].tolist())),
            "groups_dropped": int(sum(np.all(b == 0) for b in h["advantage"].reshape(-1, 8))),
            "adv_abs_max": float(np.abs(h["advantage"]).max())})
    keys = sorted({k for l in logs for k in l if k.startswith("rewards/")})
    res = {"args": vars(a), "status": status, "error": err, "wall_s": round(wall, 1),
           "s_per_step": round(wall / max(1, len(logs)), 1), "peak_rss_gb": round(peak_rss_gb(), 2),
           "mps_driver_gb": round(mps_mem_gb(), 2), "reward_log_keys": keys, "nan_in_logs": nan,
           "trl_vs_aggregator_adv_maxdiff": diffs, "batches": batches,
           "parse_errors": {e: errors.count(e) for e in set(errors)},
           "logs": [{k: (round(v, 5) if isinstance(v, float) else v) for k, v in l.items()
                     if k.startswith("rewards/") or k in ("step", "loss", "kl", "reward", "frac_reward_zero_std",
                                                          "completions/mean_length", "completions/clipped_ratio",
                                                          "step_time", "entropy")} for l in logs]}
    res["non_degenerate"] = bool(any(b["groups_dropped"] < len(b["iou_tc_group_std"]) for b in batches))
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / f"k12_{a.tag}.json").write_text(json.dumps(res, indent=1, default=str))
    print("RESULT", json.dumps({k: v for k, v in res.items() if k != "logs"}, default=str))
    for l in res["logs"]:
        print("LOG", json.dumps(l))


if __name__ == "__main__":
    main()
