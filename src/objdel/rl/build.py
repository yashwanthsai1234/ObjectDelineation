"""The GRPO trainer of a run: one builder for the G-probe, the full run (objdel.rl.run), the final exam, the C checks.

    built = build_trainer(model, ckpt, rows, out, roots, steps, flags, **over)   # Built(trainer, processor, lr)

The arguments are typed: `flags` is an objdel.core.spec.TrainerFlags (an RlStage or a ProbeStage passes itself), so no
caller builds a fake stage. GRPOConfig overrides merge in one fixed order, later wins:
{learning_rate, max_steps, save_strategy} | over | {logging_steps 1, log_completions False} | the model's GRPO
sampling | flags.grpo (golden 6 pins the result).
"""

import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NamedTuple

from objdel.core.spec import TrainerFlags


class Built(NamedTuple):
    """What build_trainer returns: the trainer, its processor and the learning rate it uses."""

    trainer: Any
    processor: Any
    lr: float


def build_trainer(model: str, ckpt: Path, rows, out: Path, roots: Sequence[Path], steps: int, flags: TrainerFlags, *,
                  eval_dataset=None, check=None, calls=(), reward=None, record_dir: Path | None = None,
                  **over) -> Built:
    """FastGRPO (or, with flags.plain, TRL's GRPOTrainer) on ckpt with configs.grpo_config, the model's GRPO sampling
    and the probe LR (configs.rl_lr; flags.lr overrides when set). rows: the training rows; out: the run folder
    (the trainer writes out/trainer); roots: the prepared roots whose caches the reward reads (the held-out view
    first when eval_dataset is given). check(config) runs before the model loads (train_rl's sampler guards); calls:
    more callbacks; reward: the reward spec (objdel.reward.REWARDS; None: reward v4); record_dir: C0's batch folder
    (train_rl opens it, checked empty, before the load)."""
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl import GRPOTrainer

    from objdel.reward import REWARD_V4
    from objdel.train import callbacks, configs, fast, rewards, serve, vllm_sync

    peft = configs.lora(model)
    lr = flags.lr or configs.rl_lr(model)
    over = {"learning_rate": lr, "max_steps": steps, "save_strategy": "no"} | over
    over |= {"logging_steps": 1, "log_completions": False} | configs.sampling(model, "grpo") | dict(flags.grpo)
    reward = reward or REWARD_V4
    cfg = configs.grpo_config(model, str(Path(out) / "trainer"), reward=reward, **over)
    if check:
        check(cfg)
    net = AutoModelForImageTextToText.from_pretrained(ckpt, **({"dtype": torch.float16} if cfg.fp16 else {}))
    proc = AutoProcessor.from_pretrained(ckpt, **configs.processor_kwargs(model, "lidc"))
    funcs, _ = rewards.TRLReward(serve.loader(*roots), rewards.stop_ids(net, proc), spec=reward).funcs()
    cls = GRPOTrainer if flags.plain else fast.FastGRPO
    win = flags.profile_window  # --profile-window S,E: NVTX steps and the nsys capture window
    calls = [*calls] + ([callbacks.NvtxSteps(tuple(win))] if win else [])
    if flags.trace:  # torch.profiler on step S (the nsys fallback; raises on S < 2), else on step 3
        calls.append(callbacks.TraceStep(flags.trace, wait=win[0] - 2 if win else 1))
    if flags.no_ckpt_layers:  # C9 arm (a)
        calls.append(callbacks.NoCkptLayers(flags.no_ckpt_layers))
    trainer = cls(model=net, reward_funcs=funcs, args=cfg, train_dataset=rows, eval_dataset=eval_dataset,
                  processing_class=proc, peft_config=peft, callbacks=calls)  # fmt: skip
    # FastGRPO: pack each group's prompt once; off when compiled (packed widths are not bucketed: a recompile each)
    trainer.shared_prefix = flags.shared_prefix and not flags.compile
    trainer.log_by_prompt = flags.log_by_prompt
    if cfg.off_policy_mask_threshold is not None and cls is fast.FastGRPO:  # bf16-tis-opsm: B2.4 rule 4's live delta
        trainer.opsm_watch = callbacks.OpsmDelta(trainer, cfg.output_dir)
        trainer.add_callback(trainer.opsm_watch)
    trainer.group_requests = flags.group_requests  # C6: one n = G vLLM request per prompt
    trainer.runaway_cap = bool(flags.runaway_cap)  # C11: each prompt's max_tokens from the rows' `cap`
    if cfg.vllm_lora and cls is fast.FastGRPO:  # option F: vLLM serves the unmerged adapter from a tmpfs folder
        shm = "/dev/shm" if Path("/dev/shm").is_dir() else None
        trainer.lora_dir = Path(tempfile.mkdtemp(prefix="objdel-lora-", dir=shm))
        if not vllm_sync.install(trainer, "lora-seqmask", trainer.lora_dir):  # C5's sync_lora each sync
            raise ValueError("--numerics lora-seqmask needs a colocate vLLM on one process with a LoRA model")
    elif flags.fast_sync and not vllm_sync.install(trainer):  # C5: one load_weights call per sync
        raise ValueError("--fast-sync needs a colocate vLLM on one process with a LoRA model")
    trainer.count_syncs = flags.count_syncs  # C0
    trainer.host_meta = flags.host_meta  # C4
    trainer.check_meta_every = flags.check_meta_every
    trainer.check_rope_every = flags.check_rope_every  # C2
    trainer.checks_from_step = flags.checks_from_step  # C0: the check flags and recording act from step S
    trainer.record_dir = Path(record_dir) if record_dir else None
    trainer.prompt_memo = flags.prompt_memo  # C2
    if flags.encoder_once and not trainer.use_encoder_once():  # C3, after host_meta
        raise ValueError("--encoder-once needs a Qwen-VL vision tower, base and image processor")
    if flags.gpu_pixels == "auto":  # C1: 8 training images for the start-up self-check
        imgs = [rows[i]["image"] for i in range(min(8, len(rows)))]
        print(f"gpu_pixels/backend={trainer.use_gpu_pixels(imgs, flags.allow_cpu_pixels)}")
    if flags.compile:
        fast.compile_text_layers(trainer.model, flags.compile)
        trainer.buckets = fast.BUCKETS  # few shapes for the compiled layers
    return Built(trainer, proc, lr)
