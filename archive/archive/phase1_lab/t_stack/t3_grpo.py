"""Test 3: TRL GRPOTrainer smoke on a small VLM (HF generate, no vLLM) with objdel rewards + a controller callback.

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I t3_grpo.py --model smolvlm-256m --tag base
    ... --agg normalize_then_sum --beta 0.001 --ctrl-step 2 --beta-to 0.5 --weights-to 0,0,0,5 --tag ctrl

Writes results/t3_<model>_<tag>.json with every logged step, timing, memory and the controller trace.
"""
import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import torch  # noqa: E402
from datasets import Dataset  # noqa: E402
from transformers import TrainerCallback  # noqa: E402

import rewards as R  # noqa: E402
from common import HERE, PROMPT_POLY, model_path, mps_mem_gb, peak_rss_gb, riga_samples  # noqa: E402


class Controller(TrainerCallback):
    """Level-2 controller stand-in: at `step` it rewrites trainer.beta and trainer.reward_weights in place."""

    def __init__(self, step, beta_to=None, weights_to=None):
        self.step, self.beta_to, self.weights_to = step, beta_to, weights_to
        self.trainer = None
        self.trace = []

    def on_step_begin(self, args, state, control, **kw):
        t = self.trainer
        if state.global_step == self.step:
            if self.beta_to is not None:
                t.beta = self.beta_to
            if self.weights_to is not None:
                t.reward_weights = torch.tensor(self.weights_to, dtype=torch.float32)
        self.trace.append({"step_begin": state.global_step, "beta": t.beta, "weights": t.reward_weights.tolist()})


class LogSink(TrainerCallback):
    def __init__(self):
        self.logs = []

    def on_log(self, args, state, control, logs=None, **kw):
        self.logs.append({"step": state.global_step, **{k: v for k, v in (logs or {}).items()}})


def build_dataset(n, side):
    samples = riga_samples(n, seed=1, side=side)
    rows = []
    for s in samples:
        R.register(s["key"], s["disc_masks"])
        rows.append({"prompt": [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPT_POLY}]}],
                     "image": s["image"], "key": s["key"]})
    return Dataset.from_list(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="smolvlm-256m")
    ap.add_argument("--tag", default="base")
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--G", type=int, default=4)
    ap.add_argument("--max-completion", type=int, default=96)
    ap.add_argument("--agg", default="sum_then_normalize")
    ap.add_argument("--loss-type", default="dr_grpo")
    ap.add_argument("--scale-rewards", default="none")
    ap.add_argument("--beta", type=float, default=0.001)
    ap.add_argument("--epsilon-high", type=float, default=0.28)
    ap.add_argument("--ctrl-step", type=int, default=-1)
    ap.add_argument("--beta-to", type=float, default=None)
    ap.add_argument("--weights-to", default=None)
    ap.add_argument("--peft", action="store_true")
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--side", type=int, default=448)
    ap.add_argument("--grad-ckpt", type=int, default=1)
    ap.add_argument("--lr", type=float, default=1e-5)
    ap.add_argument("--init", default=None, help="local checkpoint dir (e.g. SFT output) to start from")
    a = ap.parse_args()

    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl import GRPOConfig, GRPOTrainer

    t0 = time.time()
    path = model_path(a.model)
    pkw = {"size": {"longest_edge": 512}} if a.model.startswith("smolvlm") else {}
    proc = AutoProcessor.from_pretrained(path, padding_side="left", **pkw)
    dtype = torch.bfloat16 if a.model in ("qwen3-vl-2b", "qwen3.5-0.8b") else torch.float32
    model = AutoModelForImageTextToText.from_pretrained(a.init or path, dtype=dtype)
    ds = build_dataset(4, a.side)

    peft_config = None
    if a.peft:
        from peft import LoraConfig
        peft_config = LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", task_type="CAUSAL_LM")

    funcs = [R.format_reward, R.iou_reward, R.topo_reward, R.legality_const]
    cfg = GRPOConfig(
        output_dir=str(HERE / "runs" / f"grpo_{a.model}_{a.tag}"),
        max_steps=a.steps, per_device_train_batch_size=a.G, num_generations=a.G, gradient_accumulation_steps=1,
        max_completion_length=a.max_completion, temperature=1.0, learning_rate=a.lr, lr_scheduler_type="constant", warmup_steps=0,
        loss_type=a.loss_type, scale_rewards=a.scale_rewards, multi_objective_aggregation=a.agg,
        beta=a.beta, epsilon=0.2, epsilon_high=a.epsilon_high, mask_truncated_completions=False,
        reward_weights=[1.0, 1.0, 1.0, 1.0],
        use_vllm=False, bf16=False, use_cpu=a.cpu, gradient_checkpointing=bool(a.grad_ckpt),
        logging_steps=1, report_to="none", save_strategy="no", log_completions=False, seed=0,
        dataloader_num_workers=0, remove_unused_columns=False,
    )
    ctrl = Controller(a.ctrl_step,
                      a.beta_to, [float(x) for x in a.weights_to.split(",")] if a.weights_to else None)
    sink = LogSink()
    trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=cfg, train_dataset=ds, processing_class=proc,
                          peft_config=peft_config, callbacks=[ctrl, sink])
    ctrl.trainer = trainer
    t_init = time.time() - t0
    print(f"init {t_init:.1f}s device={trainer.accelerator.device} ref_model={trainer.ref_model is not None}", flush=True)
    status, err = "pass", None
    t1 = time.time()
    try:
        trainer.train()
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        status, err = "fail", f"{type(e).__name__}: {e}"[:800]
    t_train = time.time() - t1

    steps = [l for l in sink.logs if "loss" in l]
    keys = sorted({k for l in steps for k in l if k.startswith("rewards/")})
    nan = any(isinstance(v, float) and math.isnan(v) for l in steps for v in l.values())
    # does the logged `reward` equal sum_k w_k * rewards/k/mean with the weights active at that step?
    names = [f.__name__ for f in funcs]
    checks = []
    for l in steps:
        tr = next((t for t in ctrl.trace if t["step_begin"] == l["step"] - 1), None)
        if tr and all(f"rewards/{n}/mean" in l for n in names):
            pred = sum(w * l[f"rewards/{n}/mean"] for w, n in zip(tr["weights"], names))
            checks.append({"step": l["step"], "weights": tr["weights"], "beta": tr["beta"], "logged_reward": l.get("reward"),
                           "sum_w_mean": pred, "kl": l.get("kl"), "loss": l.get("loss")})
    adv = trainer._logs["advantages"] if hasattr(trainer, "_logs") else []
    out = {
        "args": vars(a), "status": status, "error": err, "device": str(trainer.accelerator.device),
        "ref_model_loaded": trainer.ref_model is not None,
        "init_s": round(t_init, 1), "train_s": round(t_train, 1), "s_per_step": round(t_train / max(1, len(steps)), 1),
        "peak_rss_gb": round(peak_rss_gb(), 2), "mps_driver_gb": round(mps_mem_gb(), 2),
        "reward_log_keys": keys, "any_nan_in_logs": nan, "weight_checks": checks,
        "advantages_logged": [round(x, 4) for x in list(adv)][:64], "controller_trace": ctrl.trace, "logs": sink.logs,
        "rewards_logged": {k: list(v)[:64] for k, v in trainer._logs["rewards"].items()} if hasattr(trainer, "_logs") else {},
        "sample_completions": [str(c)[:300] for c in list(trainer._logs["completion"])[:4]] if hasattr(trainer, "_logs") else [],
    }
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / f"t3_{a.model}_{a.tag}.json").write_text(json.dumps(out, indent=1, default=str))
    print("RESULT", json.dumps({k: out[k] for k in ("status", "error", "device", "ref_model_loaded", "init_s", "train_s",
                                                   "s_per_step", "peak_rss_gb", "mps_driver_gb", "reward_log_keys",
                                                   "any_nan_in_logs")}, default=str))
    for c in checks:
        print("CHECK", json.dumps(c))


if __name__ == "__main__":
    main()
