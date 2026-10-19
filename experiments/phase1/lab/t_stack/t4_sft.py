"""Test 4: TRL SFTTrainer overfit of a small VLM on 4 (fundus crop, polygon JSON) pairs, assistant-only loss.

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I t4_sft.py --model smolvlm-256m --steps 50

Then greedy-generates for the same 4 images and scores the polygons with objdel. Saves the model to runs/sft_<model>.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import torch  # noqa: E402
from datasets import Dataset  # noqa: E402
from transformers import TrainerCallback  # noqa: E402

import rewards as R  # noqa: E402
from common import HERE, PROMPT_POLY, mask_to_polygon, model_path, mps_mem_gb, peak_rss_gb, polygon_json, riga_samples  # noqa: E402


class LogSink(TrainerCallback):
    def __init__(self):
        self.logs = []

    def on_log(self, args, state, control, logs=None, **kw):
        self.logs.append({"step": state.global_step, **(logs or {})})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="smolvlm-256m")
    ap.add_argument("--steps", type=int, default=50)
    ap.add_argument("--lr", type=float, default=5e-5)
    ap.add_argument("--vertices", type=int, default=24)
    ap.add_argument("--assistant-only", type=int, default=1)
    ap.add_argument("--peft", action="store_true")
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl import SFTConfig, SFTTrainer

    t0 = time.time()
    path = model_path(a.model)
    pkw = {"size": {"longest_edge": 512}} if a.model.startswith("smolvlm") else {}
    proc = AutoProcessor.from_pretrained(path, **pkw)
    dtype = torch.bfloat16 if a.model in ("qwen3-vl-2b", "qwen3.5-0.8b") else torch.float32
    model = AutoModelForImageTextToText.from_pretrained(path, dtype=dtype)

    samples = riga_samples(4, seed=1)  # same 4 crops as the GRPO smoke
    rows, targets = [], []
    for s in samples:
        R.register(s["key"], s["disc_masks"])
        tgt = polygon_json(mask_to_polygon(R.GT[s["key"]].vote, n=a.vertices))
        targets.append(tgt)
        rows.append({"messages": [
            {"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPT_POLY}]},
            {"role": "assistant", "content": [{"type": "text", "text": tgt}]}], "image": s["image"], "key": s["key"]})
    # the targets themselves should score ~ IoU of a 24-gon vs the majority mask
    target_iou = [R.iou_reward([t], [s["key"]])[0] for t, s in zip(targets, samples)]
    ds = Dataset.from_list(rows)

    peft_config = None
    if a.peft:
        from peft import LoraConfig
        peft_config = LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", task_type="CAUSAL_LM")
    tag = a.tag or ("lora" if a.peft else "full")
    out_dir = HERE / "runs" / f"sft_{a.model}_{tag}"
    cfg = SFTConfig(
        output_dir=str(out_dir), max_steps=a.steps, per_device_train_batch_size=4, gradient_accumulation_steps=1,
        learning_rate=a.lr, lr_scheduler_type="constant", warmup_steps=0, max_length=None,
        assistant_only_loss=bool(a.assistant_only), bf16=False, use_cpu=a.cpu, gradient_checkpointing=True,
        logging_steps=5, report_to="none", save_strategy="no", seed=0, dataloader_num_workers=0,
    )
    sink = LogSink()
    status, err = "pass", None
    try:
        trainer = SFTTrainer(model=model, args=cfg, train_dataset=ds, processing_class=proc,
                             peft_config=peft_config, callbacks=[sink])
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        out = {"args": vars(a), "status": "fail_init", "error": f"{type(e).__name__}: {e}"[:800]}
        (HERE / "results").mkdir(exist_ok=True)
        (HERE / "results" / f"t4_{a.model}_{tag}.json").write_text(json.dumps(out, indent=1))
        print("RESULT", json.dumps(out))
        return
    # inspect the label mask of one collated batch: how many label tokens, and do they decode to the target?
    batch = trainer.data_collator([trainer.train_dataset[0]])
    lab = batch["labels"][0]
    kept = lab[lab != -100]
    label_text = proc.tokenizer.decode(kept) if hasattr(proc, "tokenizer") else ""
    print(f"labels kept={len(kept)} of {len(lab)}: {label_text[:120]!r} ... {label_text[-40:]!r}", flush=True)
    t_init = time.time() - t0
    t1 = time.time()
    try:
        trainer.train()
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        status, err = "fail_train", f"{type(e).__name__}: {e}"[:800]
    t_train = time.time() - t1
    losses = [(l["step"], round(l["loss"], 4)) for l in sink.logs if "loss" in l]
    print("losses", losses, flush=True)

    # greedy generation on the training images
    gens = []
    model = trainer.model.eval()
    dev = trainer.accelerator.device
    t2 = time.time()
    for s in samples:
        msgs = [{"role": "user", "content": [{"type": "image", "image": s["image"]}, {"type": "text", "text": PROMPT_POLY}]}]
        inp = proc.apply_chat_template(msgs, tokenize=True, add_generation_prompt=True, return_dict=True,
                                       return_tensors="pt").to(dev)
        with torch.no_grad():
            o = model.generate(**inp, max_new_tokens=400, do_sample=False)
        txt = proc.batch_decode(o[:, inp["input_ids"].shape[1]:], skip_special_tokens=True)[0]
        pts, probs = R.parse(txt)
        iou = R.iou_reward([txt], [s["key"]])[0]
        gens.append({"key": s["key"], "parsed": pts is not None, "problems": ";".join(probs) or "ok",
                     "n_vertices": None if pts is None else len(pts), "iou": round(iou, 3), "text": txt[:300]})
        print(json.dumps(gens[-1]), flush=True)
    t_gen = time.time() - t2
    if status == "pass":
        trainer.save_model(str(out_dir))
        proc.save_pretrained(str(out_dir))
    out = {
        "args": vars(a), "status": status, "error": err, "device": str(dev), "init_s": round(t_init, 1),
        "train_s": round(t_train, 1), "s_per_step": round(t_train / max(1, a.steps), 2), "gen_s": round(t_gen, 1),
        "peak_rss_gb": round(peak_rss_gb(), 2), "mps_driver_gb": round(mps_mem_gb(), 2),
        "label_tokens": int(len(kept)), "seq_tokens": int(len(lab)), "label_text_head": label_text[:200],
        "losses": losses, "target_iou": [round(x, 3) for x in target_iou],
        "parse_rate": sum(g["parsed"] for g in gens) / len(gens),
        "valid_rate": sum(g["problems"] == "ok" for g in gens) / len(gens),
        "mean_iou": round(sum(g["iou"] for g in gens) / len(gens), 3), "generations": gens,
        "saved_to": str(out_dir) if status == "pass" else None,
    }
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / f"t4_{a.model}_{tag}.json").write_text(json.dumps(out, indent=1, default=str))
    print("RESULT", json.dumps({k: v for k, v in out.items() if k not in ("generations", "losses")}, default=str))


if __name__ == "__main__":
    main()
