"""G-probe (PLAN §5.1, §3.2 ④): can GRPO improve this checkpoint? 20 GRPO steps at the pre-registered probe LR.

    python -I -m experiments.phase2.g_probe --model lfm2.5-vl-450m (--ckpt runs/r1/branches/e1 | --select runs/r1/select)
        --train data/prepared --probe data/eval_views/lidc-val-shifted --out runs/r1/g_probe
    python -I -m experiments.phase2.g_probe --spec specs/r1_qwen3vl2b.yaml --stage probe [--rd RD] [--data DATA]

The flags are an objdel.core.spec.ProbeStage (cli_stage), with the trainer flags train_rl has (TrainerFlags: --lr,
--grpo JSON, --shared-prefix, --compile, --plain, --trace, the C flags).

configs.grpo_config with prereg.G_PROBE: LR 5e-6 (full FT) or 5e-5 (LoRA), constant, 20 steps, G = 8, 32 LIDC train
prompts per step (disagreement items x4), vLLM colocate generation (owner 2026-11-06, perf F1: HF generate of 256
sequences needs 20-51 GiB of KV cache; needs the Phase 3 single-venv image, D15), sampling from
configs.sampling(model, "grpo") (D13);
fast.FastGRPO (per-chunk trimmed forwards, a per-step vision cache: perf C1, C6), TRLReward on the
preloaded caches; the probe metrics run in this process (D7). Probe set: a seeded subset of `items`
(128) shifted-val records (select_ckpt.view_rows), K = 8 through the trainer's own generation (sample.sample), twice
at step 0 (the null, h0) and once at step 20; per item the raw composite R. Verdict: prereg.g_probe, step-0 legal =
the mean of both draws. Also: Spearman of step-0 sd_geom / r_std / mixed with the gain (D12 veto, §6.2), the
training-log form (mean logged R, steps 1-5 vs 11-20) and the step timing.

--select (select_ckpt's dir): probe the selected branch, then (not PASS) the D12(c) repair: the next-earlier
eligible branch, once. <out>/g_probe.json = the last probe's result + "probed" (names) + "earlier" (results before);
its "ckpt" received the final verdict (C2 depth); NaN (e.g. Spearman of constant gains) is written as null. When
R1's §5.1 gates failed (selected.json "gates_pass" not true) nothing is probed: {"verdict": null, "skipped": "§5.1
gates failed: <gate names>", "probed": [], "gates_branch"}, exit 0, so the chain (train_sft && select_ckpt &&
g_probe) ends cleanly with the reason on disk.
"""
import gc
import json
import math
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from experiments.phase2.prereg import G_PROBE, g_probe
from experiments.phase2.select_ckpt import score_items, view_rows
from objdel.eval import signal

SEEDS = (11, 12, 13)  # step-0 draw a, step-0 draw b, step 20


def score(trainer, rows, probe_root: Path, seed: int) -> dict:
    """K = 8 samples per probe row -> per-item panel metrics."""
    import torch

    from objdel.train import sample

    torch.manual_seed(seed)
    _, groups = score_items(rows, {d["item_id"]: d for d in sample.sample(trainer, rows, 8)}, probe_root)
    m = signal.per_item(groups)
    ids = [g["item_id"] for g in groups]
    return {"ids": ids, "per_item": m, "legal": float(np.mean(m["legal"])), "sigma_share": float(np.mean(m["sigma"]))}


def verdict(s0a: dict, s0b: dict, s20: dict) -> dict:
    """prereg.g_probe: gain = R20 - mean of the two step-0 draws; step-0 legal = mean of both draws."""
    r0a, r0b, r20 = (s["per_item"]["r_mean"] for s in (s0a, s0b, s20))
    return g_probe(r20, r0a, r0b, s20["ids"], (s0a["legal"] + s0b["legal"]) / 2, s20["legal"], s20["sigma_share"])


def train_log(log: list[dict]) -> dict:
    """§5.1 training-log form: mean logged raw composite R, steps 1-5 vs 11-20 (TRL logs rewards/<term>/mean)."""
    R = [sum(w * h.get(f"rewards/{t}/mean", np.nan) for t, w in signal.WEIGHTS.items()) for h in log
         if "rewards/iou_tc/mean" in h]  # fmt: skip
    return {"steps_1_5": float(np.mean(R[:5])), "steps_11_20": float(np.mean(R[10:20]))}


def candidates(sel: dict) -> list[str]:
    """D12(c): the selected branch, then the next-earlier (shallower) eligible branch."""
    b, depth = sel["branches"], sel["branches"][sel["selected"]]["depth"]
    earlier = sorted((x["depth"], n) for n, x in b.items() if x["eligible"] and x["depth"] < depth)
    return [sel["selected"]] + [n for _, n in earlier[-1:]]


def make_trainer(a, ckpt: Path, train, steps: int, eval_dataset=None, check=None, calls=(), reward=None, **over):
    """(trainer, processor, lr): FastGRPO (or, with a.plain, TRL's GRPOTrainer) on ckpt with configs.grpo_config,
    the model's GRPO sampling, the probe LR (LoRA or full FT; a.lr overrides when set), logging every step, then
    `over` and a.grpo's overrides. `a` is a stage with the objdel.core.spec.TrainerFlags fields (RlStage, ProbeStage). check(config) runs before the model loads (train_rl's sampler guards); calls:
    more callbacks; reward: the reward spec (objdel.reward.REWARDS; None: reward v4); eval_dataset: rows of the view
    a.heldout (its caches resolve first). Shared by the G-probe and the full run (train_rl)."""
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl import GRPOTrainer

    from objdel.reward import REWARD_V4
    from objdel.train import callbacks, configs, fast, rewards, serve, vllm_sync

    peft = configs.lora(a.model)
    lr = getattr(a, "lr", None) or configs.rl_lr(a.model)
    over = {"learning_rate": lr, "max_steps": steps, "save_strategy": "no"} | over
    over |= {"logging_steps": 1, "log_completions": False} | configs.sampling(a.model, "grpo") | dict(a.grpo)
    reward = reward or REWARD_V4
    cfg = configs.grpo_config(a.model, str(Path(a.out) / "trainer"), reward=reward, **over)
    if check:
        check(cfg)
    model = AutoModelForImageTextToText.from_pretrained(ckpt, **({"dtype": torch.float16} if cfg.fp16 else {}))
    proc = AutoProcessor.from_pretrained(ckpt, **configs.processor_kwargs(a.model, "lidc"))
    roots = [a.train] if eval_dataset is None else [a.heldout, a.train]
    funcs, _ = rewards.TRLReward(serve.loader(*roots), rewards.stop_ids(model, proc), spec=reward).funcs()
    cls = GRPOTrainer if getattr(a, "plain", False) else fast.FastGRPO
    win = getattr(a, "profile_window", None)  # --profile-window S,E: NVTX steps and the nsys capture window
    calls = [*calls] + ([callbacks.NvtxSteps(tuple(win))] if win else [])  # raises on S < 2 (TraceStep's wait < 0)
    if getattr(a, "trace", None):  # torch.profiler on step S (the nsys fallback), else on step 3
        calls.append(callbacks.TraceStep(a.trace, wait=win[0] - 2 if win else 1))
    if getattr(a, "no_ckpt_layers", None):  # C9 arm (a)
        calls.append(callbacks.NoCkptLayers(a.no_ckpt_layers))
    trainer = cls(model=model, reward_funcs=funcs, args=cfg, train_dataset=train, eval_dataset=eval_dataset,
                  processing_class=proc, peft_config=peft, callbacks=calls)  # fmt: skip
    # FastGRPO: pack each group's prompt once; off when compiled (packed widths are not bucketed: a recompile each)
    trainer.shared_prefix = bool(getattr(a, "shared_prefix", False)) and not getattr(a, "compile", None)
    trainer.log_by_prompt = bool(getattr(a, "log_by_prompt", False))
    if cfg.off_policy_mask_threshold is not None and cls is fast.FastGRPO:  # bf16-tis-opsm: B2.4 rule 4's live delta
        trainer.opsm_watch = callbacks.OpsmDelta(trainer, cfg.output_dir)
        trainer.add_callback(trainer.opsm_watch)
    trainer.group_requests = bool(getattr(a, "group_requests", False))  # C6: one n = G vLLM request per prompt
    trainer.runaway_cap = bool(getattr(a, "runaway_cap", 0))  # C11: each prompt's max_tokens from the rows' `cap`
    if cfg.vllm_lora and cls is fast.FastGRPO:  # option F: vLLM serves the unmerged adapter from a tmpfs folder
        shm = "/dev/shm" if Path("/dev/shm").is_dir() else None
        trainer.lora_dir = Path(tempfile.mkdtemp(prefix="objdel-lora-", dir=shm))
        if not vllm_sync.install(trainer, "lora-seqmask", trainer.lora_dir):  # C5's sync_lora each sync
            raise ValueError("--numerics lora-seqmask needs a colocate vLLM on one process with a LoRA model")
    elif getattr(a, "fast_sync", False) and not vllm_sync.install(trainer):  # C5: one load_weights call per sync
        raise ValueError("--fast-sync needs a colocate vLLM on one process with a LoRA model")
    trainer.count_syncs = bool(getattr(a, "count_syncs", False))  # C0
    trainer.host_meta = bool(getattr(a, "host_meta", False))  # C4
    trainer.check_meta_every = getattr(a, "check_meta_every", 50)
    trainer.check_rope_every = getattr(a, "check_rope_every", 50)  # C2
    trainer.checks_from_step = getattr(a, "checks_from_step", 1)  # C0: the check flags and recording act from step S
    rec = getattr(a, "record_batches", None)  # C0, opened (and checked empty) by train_rl before the load
    trainer.record_dir = Path(rec) if rec else None
    trainer.prompt_memo = bool(getattr(a, "prompt_memo", False))  # C2
    if getattr(a, "encoder_once", False) and not trainer.use_encoder_once():  # C3, after host_meta
        raise ValueError("--encoder-once needs a Qwen-VL vision tower, base and image processor")
    if getattr(a, "gpu_pixels", "off") == "auto":  # C1: 8 training images for the start-up self-check
        imgs = [train[i]["image"] for i in range(min(8, len(train)))]
        print(f"gpu_pixels/backend={trainer.use_gpu_pixels(imgs, getattr(a, 'allow_cpu_pixels', False))}")
    if getattr(a, "compile", None):
        fast.compile_text_layers(trainer.model, a.compile)
        trainer.buckets = fast.BUCKETS  # few shapes for the compiled layers
    return trainer, proc, lr


def probe(a, ckpt: Path, train, rows) -> dict:
    """One G-probe run on one checkpoint."""
    trainer, _, lr = make_trainer(a, ckpt, train, G_PROBE["steps"])
    cfg = trainer.args
    t0 = time.time()
    s0a, s0b = (score(trainer, rows, Path(a.probe).resolve(), s) for s in SEEDS[:2])
    t1 = time.time()
    trainer.train()
    t2 = time.time()
    s20 = score(trainer, rows, Path(a.probe).resolve(), SEEDS[2])
    r0a, r0b, r20 = (s["per_item"]["r_mean"] for s in (s0a, s0b, s20))
    step0 = {k: (s0a["per_item"][k] + s0b["per_item"][k]) / 2 for k in ("sd_geom", "r_std", "mixed")}
    log = trainer.state.log_history
    return verdict(s0a, s0b, s20) | {
        "lr": lr, "steps": G_PROBE["steps"], "items": len(s20["ids"]), "ckpt": str(ckpt), "seeds": SEEDS,
        "legal": {"step0": [s0a["legal"], s0b["legal"]], "step20": s20["legal"]},
        "spearman_gain": {k: signal.spearman_ci(v, r20 - (r0a + r0b) / 2, s20["ids"]) for k, v in step0.items()},
        "per_item": {"ids": s20["ids"], "r0a": r0a.tolist(), "r0b": r0b.tolist(), "r20": r20.tolist()},
        "timing_s": {"step0_sampling": t1 - t0, "train": t2 - t1, "per_step": (t2 - t1) / cfg.max_steps},
        "train_log_R": train_log(log), "log": log,
    }  # fmt: skip


def main(argv=None) -> int:
    from objdel.core.spec import ProbeStage, cli_stage
    from objdel.train import serve

    a = cli_stage(ProbeStage, sys.argv[1:] if argv is None else argv)
    out = Path(a.out)
    ckpts = [Path(a.ckpt)] if a.ckpt else []
    if a.select:
        sel = json.loads((Path(a.select) / "selected.json").read_text())
        if not (sel["gates_pass"] and sel["selected"]):  # R1 failed §5.1: a probe would not change the outcome
            why = ", ".join(sel.get("gates_failed") or ["no eligible branch (D12)"])
            return save(out, {"verdict": None, "skipped": f"§5.1 gates failed: {why}", "probed": [],
                            "gates_branch": sel["gates_branch"]})
        ckpts = [Path(sel["run"]) / "branches" / n for n in candidates(sel)]
    for root in (a.train, a.probe):
        serve.preload(root, ["lidc"])
    from PIL import Image

    from objdel.train import dataset

    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=dataset.DISAGREEMENT_OVERSAMPLE)
    rows = view_rows(a.probe, a.items or G_PROBE["items"])
    rows = [r | {"image": Image.open(r["image"]).convert("RGB")} for r in rows]
    docs = []
    for ckpt in ckpts:  # at most the selected branch and one repair branch
        docs.append(probe(a, ckpt, train, rows))
        gc.collect()  # the trainer's reference cycles hold the model
        if docs[-1]["verdict"] == "PASS":
            break
    return save(out, docs[-1] | {"probed": [c.name for c in ckpts[: len(docs)]], "earlier": docs[:-1]})


def finite(x):
    """NaN / inf -> None, recursively (Spearman of constant gains is NaN; strict JSON has no NaN)."""
    if isinstance(x, dict):
        return {k: finite(v) for k, v in x.items()}
    if isinstance(x, list):
        return [finite(v) for v in x]
    return None if isinstance(x, float) and not math.isfinite(x) else x


def save_json(path: Path, doc: dict, keys: tuple[str, ...]) -> int:
    """Strict JSON (numpy scalars as floats, NaN / inf as null) to path; prints doc's `keys` on one line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = finite(json.loads(json.dumps(doc, default=float)))  # numpy scalars -> floats first
    path.write_text(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: doc[k] for k in keys if k in doc}, ensure_ascii=False))
    return 0


def save(out: Path, doc: dict) -> int:
    return save_json(out / "g_probe.json", doc, ("verdict", "skipped", "h0", "gain", "probed"))


if __name__ == "__main__":
    sys.exit(main())
