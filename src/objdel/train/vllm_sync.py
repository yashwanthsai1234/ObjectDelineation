"""C5: the vLLM weight sync in one `load_weights` call per step (train_rl --fast-sync), and the option-F LoRA sync.

TRL's colocate `sync_weights` calls `load_weights([(name, param)])` once per parameter (about 700 calls; vLLM's
AutoWeightsLoader rebuilds its maps on each). `push_weights` hands vLLM one generator of the same (name, value) pairs:
the same values reach the same vLLM tensors, with one call. `sync_lora` is the option-F sync (preset lora-seqmask,
B2.2): vLLM keeps the adapter apart from the base, as PEFT does in the trainer.

Sleep state: TRL's flag `vg._llm_weights_sleeping` is the only one read. Both syncs wake the weights only when it is
True and set it False right after the wake, in TRL's order (empty_cache, wake_up(tags=["weights"])). Level-2 sleep
discards the weights and the LoRA slots, so every sync after a sleep re-pushes both.
"""

import time
from typing import NamedTuple

from objdel.train.fast import lora_merged_params

LORA_NAME = "policy"  # the vLLM name of the trained adapter ("default" in PEFT); generate's LoRARequest uses it too


class PushStats(NamedTuple):
    """What one weight push to vLLM moved: tensor count, bytes and the load seconds."""
    n_tensors: int
    bytes: int
    seconds: float  # the load_weights / add_lora calls alone (the wake has its own timer, time/wake_s)


def _wake(vg) -> None:
    if vg._llm_weights_sleeping:
        from trl.generation import vllm_generation

        vllm_generation.empty_cache()
        vg.llm.wake_up(tags=["weights"])
        vg._llm_weights_sleeping = False


def _load(vg, pairs) -> PushStats:
    """One load_weights call on the in-process engine's model, over a generator (peak: one merged tensor)."""
    n = size = 0

    def counted():
        nonlocal n, size
        for name, t in pairs:
            n, size = n + 1, size + t.numel() * t.element_size()
            yield name, t

    t0 = time.perf_counter()
    vg.llm.llm_engine.model_executor.driver_worker.model_runner.model.load_weights(counted())
    return PushStats(n, size, time.perf_counter() - t0)


def push_weights(vg, model, values: str = "merged") -> PushStats:
    """Every weight in one call. values "merged": W + delta (fast.lora_merged_params, today's values); "base": W."""
    if values not in ("merged", "base"):
        raise ValueError(f"values {values!r}: 'merged' or 'base'")
    _wake(vg)
    stats = _load(vg, lora_merged_params(model, with_delta=values == "merged"))
    vg.llm.reset_prefix_cache()
    return stats


def write_adapter(model, adapter_dir, adapter_name: str = "default") -> tuple[int, int]:
    """The PEFT adapter folder vLLM loads: adapter_model.safetensors (bf16, keys without `.default.`) and
    adapter_config.json. Returns (tensors, bytes)."""
    from pathlib import Path

    import torch
    from peft import get_peft_model_state_dict
    from safetensors.torch import save_file

    d = Path(adapter_dir)
    d.mkdir(parents=True, exist_ok=True)
    sd = {k: v.detach().to(torch.bfloat16).contiguous() for k, v in
          get_peft_model_state_dict(model, adapter_name=adapter_name).items()}  # fmt: skip
    save_file(sd, str(d / "adapter_model.safetensors"))
    model.peft_config[adapter_name].save_pretrained(str(d))
    return len(sd), sum(v.numel() * 2 for v in sd.values())


def sync_lora(vg, model, adapter_dir, adapter_name: str = "default", slot: int = 1) -> PushStats:
    """The option-F sync: wake (if asleep); in mode B (sleep mode on) the base W pushed in one call; the adapter
    written to adapter_dir and loaded once with load_inplace=True (it also refills the slots level-2 sleep
    discards); the prefix cache reset (its key holds the LoRA name, not its weights)."""
    from vllm.lora.request import LoRARequest

    _wake(vg)
    base = _load(vg, lora_merged_params(model, with_delta=False)) if vg.enable_sleep_mode else PushStats(0, 0, 0.0)
    n, size = write_adapter(model, adapter_dir, adapter_name)
    t0 = time.perf_counter()
    name = LORA_NAME if adapter_name == "default" else adapter_name
    vg.llm.llm_engine.add_lora(LoRARequest(name, slot, str(adapter_dir), load_inplace=True))
    vg.llm.reset_prefix_cache()
    return PushStats(base.n_tensors + n, base.bytes + size, base.seconds + time.perf_counter() - t0)


def install(trainer, numerics: str = "bf16-seqmask", adapter_dir=None) -> bool:
    """Replaces the trainer's vLLM sync_weights (on the instance) with push_weights (merged presets) or sync_lora
    (lora-seqmask, into adapter_dir); False and nothing changed without a colocate vLLM on one process with PEFT."""
    vg = getattr(trainer, "vllm_generation", None)
    model = getattr(vg, "model", None)
    if vg is None or not hasattr(model, "peft_config") or getattr(vg, "mode", "colocate") != "colocate":
        return False
    if vg._dist.is_fsdp or trainer.accelerator.num_processes > 1:
        return False
    if numerics == "lora-seqmask":
        if adapter_dir is None:
            raise ValueError("lora-seqmask needs an adapter_dir (a tmpfs folder)")
        vg.sync_weights = lambda: sync_lora(vg, model, adapter_dir)  # generate's LoRARequest: FastGRPO.lora_dir
    else:
        vg.sync_weights = lambda: push_weights(vg, model)
    return True
