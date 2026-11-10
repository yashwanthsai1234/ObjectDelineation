"""B2.2 F step 5: the F check, run by `train_rl --f-check` inside the first run that uses `lora-seqmask` or
`fp16-seqmask` (one process, the run's own vLLM engine). The callback counts NaN/inf and the vLLM LoRA loads per RL
step over steps 1-15, runs the preset's probes, and writes `<out>/f_check.json`: the preset, each check's value, and
`passed`, the AND of every pass of that preset (`restart_preset` reads only `preset` and `passed`). On a failure it
raises after writing the file, so the plan step fails and the run restarts under B2.4 rule 0's fallback.

Fields under `fp16-seqmask`: nan, engine_dtype_ok. Under `lora-seqmask`: nan, comparator, base_hash_ok, lora_equal,
load_count, lora_over_base (plus prompt_expansion false when vLLM expanded an image differently). The probes are
callables; `run_comparator` measures the comparator, the restore base hash and the timing ratio before step 1 (B2.2 F
steps 4 and 4b) and they are passed in as values; `AfterSyncs` runs the packed-layout equality (`lora_equal`) right
after sync 1 and after the first sleep/wake cycle, while the engine is awake.

Restore base hash (B2.2 F step 4b.5, C5 step 5(b) method): the reference is vLLM's own tensors right after arm (ii)'s
sync, hashed on the host per name. In mode B that sync pushed `p.data` (no delta); in mode A they are the start-up
load of the same --ckpt. After the restore push every tensor is hashed again and compared per name. A direct compare
with `p.data` would need a map from PEFT names to vLLM's packed tensors (qkv_proj, gate_up_proj), unchecked here.

Output-level arms only on the comparator's model (B2.2 F step 4b.1): the --ckpt whose hf_revision.json (written by
hf_pull) names rl-full-3ep. On another model run_comparator gets `skip` (the reason): arms (0), (i), (ii) and their
scoring are skipped, `comparator` records {"skipped": reason, "informative": false, "ok": true}, and the tensor, NaN,
load-count, base-hash and timing checks decide.
"""
import json
import math
import time
from pathlib import Path

import numpy as np
from transformers import TrainerCallback

from objdel.io import read_json

STEPS = 15  # the check covers optimizer steps 1-15
PRESETS = ("fp16-seqmask", "lora-seqmask")
HALF = 0.5  # comparator: the LoRA gap above the floor is at most this share of the merged gap
FLOOR = 0.002  # below this merged gap above the floor the comparator cannot tell F from no F
LORA_OVER_BASE_MAX = 1.15  # timing: median of the LoRA / base generate-time ratios
LORA_KEYS = ("comparator", "base_hash_ok", "lora_equal", "load_count", "lora_over_base")
FIXED_LEN = 512  # tokens every row decodes in the timed arms (min_tokens = max_tokens, ignore_eos)
SEEDS = (0, 1, 2)  # timing seeds; the arms (0), (i), (ii) use seed 0
B1_SLOT = 2  # the comparator's vLLM LoRA slot; the run's own adapter is slot 1
COMPARATOR_CKPT = "rl-full-3ep"  # the HF revision the B1 adapter was trained on: the comparator's base
EXPANSION = "vLLM expanded the image"  # fast.vllm_generate_checked's error
# PEFT module -> (vLLM packed module, slice): vLLM packs q, k, v into qkv_proj and gate, up into gate_up_proj
PACKED = {"q_proj": ("qkv_proj", 0), "k_proj": ("qkv_proj", 1), "v_proj": ("qkv_proj", 2),
          "gate_proj": ("gate_up_proj", 0), "up_proj": ("gate_up_proj", 1)}


def comparator(f_merged: float, f_lora: float, f0: float) -> dict:
    """Output-level check: mean |logp gap| of the merged push, of the adapter push, and of the base at d = 0 (the
    floor). ratio = (f_lora - f0) / (f_merged - f0) passes at <= 0.5; a merged gap above the floor under 0.002 is not
    informative (ok, no ratio)."""
    gap = f_merged - f0
    out = {"f_merged": f_merged, "f_lora": f_lora, "f0": f0}
    if gap < FLOOR:
        return out | {"ratio": None, "informative": False, "ok": True}
    ratio = (f_lora - f0) / gap
    return out | {"ratio": ratio, "informative": True, "ok": ratio <= HALF}


def lora_over_base(t_lora: list[float], t_base: list[float]) -> float:
    """The timing value: the median of the t_A / t_B ratios (A = slot 2, B = `lora_request=None`), one per A B pair."""
    return float(np.median([a / b for a, b in zip(t_lora, t_base, strict=True)]))


def passes(rec: dict) -> bool:
    """The AND of every pass of the record's preset."""
    if rec["preset"] == "fp16-seqmask":
        return rec["nan"] == 0 and rec["engine_dtype_ok"]
    return (rec.get("prompt_expansion", True) is True and all(k in rec for k in LORA_KEYS) and rec["nan"] == 0
            and rec["comparator"]["ok"] and rec["base_hash_ok"] and rec["lora_equal"] and rec["load_count"] == 1
            and rec["lora_over_base"] is not None and rec["lora_over_base"] <= LORA_OVER_BASE_MAX)


def check_record(rec: dict) -> None:
    """Assert a written record: `passed` is the AND of the flags (never trusted alone) and it is True."""
    assert rec["passed"] == passes(rec), "passed differs from the AND of the checks"
    assert rec["passed"], {k: v for k, v in rec.items() if k != "passed"}


def count_calls(fn):
    """fn with a `calls` counter (the LoRA load count per RL step wraps `llm_engine.add_lora`)."""
    def wrapped(*a, **k):
        wrapped.calls += 1
        return fn(*a, **k)
    wrapped.calls = 0
    return wrapped


def comparator_skip(ckpt) -> str | None:
    """None when --ckpt is the comparator's model (its hf_revision.json names rl-full-3ep), else the reason the
    output-level check is skipped (B2.2 F step 4b.1)."""
    marker = Path(ckpt) / "hf_revision.json"
    if not marker.is_file():
        return f"{marker} missing: --ckpt is not a fetch_hf pull of {COMPARATOR_CKPT}"
    rev = read_json(marker).get("revision")
    return None if rev == COMPARATOR_CKPT else f"--ckpt is HF revision {rev!r}, not {COMPARATOR_CKPT}"


def host_hashes(vllm_model) -> dict:
    """{name: SHA-256 of the int16 view} of every vLLM parameter, hashed on the host one tensor at a time (C5 step
    5(b): no device copy of the model). Never call it while the engine sleeps (level 2 freed the weights)."""
    import hashlib

    import torch

    return {n: hashlib.sha256(p.detach().view(torch.int16).cpu().numpy().tobytes()).hexdigest()
            for n, p in vllm_model.named_parameters()}  # fmt: skip


def lora_slot_index(vg, lora_id: int) -> int:
    """The GPU slot index vLLM's adapter manager gave LoRA id `lora_id` (its stacked tensors are indexed by it)."""
    runner = vg.llm.llm_engine.model_executor.driver_worker.model_runner
    return list(runner.lora_manager._adapter_manager.lora_index_to_id).index(lora_id)


def lora_equal(vg, model, lora_id: int = 1, adapter_name: str = "default") -> bool:
    """B2.2 F step 4, packed-layout equality: for every PEFT LoRA module, the vLLM slice of the packed `lora_a`
    equals A.to(bf16) and the slice of `lora_b` equals scale x B.to(bf16) (vLLM folds alpha / r into lora_b), with
    torch.equal. Modules are matched on the name tail from the last `layers.`; a PEFT module with no vLLM twin, or
    none at all, fails."""
    import torch
    from peft.tuners.lora import LoraLayer

    from objdel.rl.vllm_sync import engine_model

    def tail(name):
        return name[name.rindex("layers."):] if "layers." in name else None

    served = {tail(n): m for n, m in engine_model(vg).named_modules() if hasattr(m, "lora_a_stacked") and tail(n)}
    idx, n = lora_slot_index(vg, lora_id), 0
    for name, m in model.named_modules():
        if not (isinstance(m, LoraLayer) and adapter_name in m.lora_A):
            continue
        parent, leaf = tail(name).rsplit(".", 1)
        packed, k = PACKED.get(leaf, (leaf, 0))
        v = served.get(f"{parent}.{packed}")
        if v is None:
            return False
        a, b = m.lora_A[adapter_name].weight, m.lora_B[adapter_name].weight
        va, vb = v.lora_a_stacked[k][idx], v.lora_b_stacked[k][idx]
        want_a = a.detach().to(torch.bfloat16).to(va.device)
        want_b = (b.detach().to(torch.bfloat16) * m.scaling[adapter_name]).to(vb.device)
        if not (torch.equal(va.reshape(want_a.shape), want_a) and torch.equal(vb.reshape(want_b.shape), want_b)):
            return False
        n += 1
    return n > 0


class AfterSyncs:
    """Wraps vg.sync_weights: check() right after syncs 1..n (mode B: n = 2, the second sync follows the first
    sleep/wake cycle; mode A: n = 3 syncs on 3 adapter states), while the engine is awake. ok(): every check ran
    and passed (None-safe: False before the n-th sync)."""

    def __init__(self, vg, check, n: int):
        self.check, self.n, self.results, self.sync = check, n, [], vg.sync_weights
        vg.sync_weights = self

    def __call__(self, *a, **k):
        """The wrapped sync, then the check while fewer than n ran."""
        out = self.sync(*a, **k)
        if len(self.results) < self.n:
            self.results.append(bool(self.check()))
        return out

    def ok(self) -> bool:
        """Every one of the n checks ran and passed."""
        return len(self.results) == self.n and all(self.results)


def sampling_params(vg, seed: int, fixed: bool = False):
    """The run's own SamplingParams (TRL's colocate fields, vllm_generation.py 676-697) with `seed`; fixed: every
    row decodes exactly FIXED_LEN tokens (B2.4 rule 0 timing)."""
    from vllm import SamplingParams

    kw = {"n": 1, "repetition_penalty": vg.repetition_penalty, "temperature": vg.temperature, "top_p": vg.top_p,
          "top_k": vg.top_k, "min_p": 0.0 if vg.min_p is None else vg.min_p, "max_tokens": vg.max_completion_length,
          "logprobs": vg.logprobs} | dict(vg.generation_kwargs or {}) | {"seed": seed}  # fmt: skip
    if fixed:
        kw |= {"max_tokens": FIXED_LEN, "min_tokens": FIXED_LEN, "ignore_eos": True}
    return SamplingParams(**kw)


def run_comparator(vg, model, prompts: list[dict], pad: int, b1_dir, tmp_dir, gap, sampling=None,
                   clock=time.perf_counter, skip: str | None = None) -> dict:
    """B2.2 F steps 4 and 4b before step 1, on the run's own engine. prompts: the 32 x 8 batch in TRL's vLLM form
    (processor-expanded ids, multi_modal_data), sent only through fast.vllm_generate_checked. gap(outputs, adapter):
    mean |vLLM logp - trainer logp| of an answer set, scored by way (b) with "b1" active (adapter True) or the
    adapter disabled. sampling(seed, fixed): the SamplingParams (default sampling_params on vg).
    Order: load "b1" beside "default"; (ii) sync_lora of b1 into slot 2; the reference hashes; (0) the floor arm on
    base W; the A B A B timing (A = slot 2, B = no adapter; one warm-up each; seeds 0-2; 512 tokens per row); (i)
    the merged push last; scoring; then the restore: base push, host hash of every base tensor against the
    reference, remove_lora(2), set_adapter("default"), delete_adapter("b1"), reset_prefix_cache, and in mode B
    sleep(level=2) last with TRL's flag set. Returns the FCheck fixed values (prompt_expansion False on an
    expansion error, then no other value). skip: the reason the output-level arms do not run (comparator_skip):
    then only the timing and the restore run, and `comparator` records the reason."""
    from vllm.lora.request import LoRARequest

    from objdel.rl.trl_adapter import vllm_generate_checked
    from objdel.rl.vllm_sync import engine_model, push_weights, sync_lora

    sampling = sampling or (lambda seed, fixed=False: sampling_params(vg, seed, fixed))
    llm, slot = vg.llm, LoRARequest("b1", B1_SLOT, str(tmp_dir))

    def gen(params, lora):
        return vllm_generate_checked(llm.generate, prompts, pad, sampling_params=params, use_tqdm=False,
                                     lora_request=lora)  # fmt: skip

    def timed(params, lora):
        t0 = clock()
        out = gen(params, lora)
        return clock() - t0, all(len(o.outputs[0].token_ids) == FIXED_LEN for o in out)

    model.load_adapter(str(b1_dir), adapter_name="b1")
    model.set_adapter("b1")
    sync_lora(vg, model, tmp_dir, adapter_name="b1", slot=B1_SLOT)
    if vg.enable_sleep_mode:
        llm.wake_up(tags=["kv_cache"])  # TRL wakes the KV cache only inside its own generate
    ref = host_hashes(engine_model(vg))
    try:
        if skip is None:
            f_lora = gen(sampling(0), slot)
            f_base = gen(sampling(0), None)
        t_a, t_b, full = [], [], True
        for lora in (slot, None):
            gen(sampling(0, True), lora)  # warm-up: fills the image and prefix caches of both arms alike
        for seed in SEEDS:
            for _ in range(2):
                for lora, times in ((slot, t_a), (None, t_b)):
                    t, ok = timed(sampling(seed, True), lora)
                    times.append(t)
                    full &= ok
        if skip is None:
            push_weights(vg, model, values="merged")  # W + d_B1 into vLLM's base tensors: every base-W arm ran before
            f_merged = gen(sampling(0), None)
    except RuntimeError as e:
        if EXPANSION not in str(e):
            raise
        return {"prompt_expansion": False, "error": str(e)}
    out = ({"skipped": skip, "informative": False, "ok": True} if skip is not None
           else comparator(gap(f_merged, True), gap(f_lora, True), gap(f_base, False)))
    fixed = {"prompt_expansion": True, "comparator": out,
             "lora_over_base": lora_over_base(t_a, t_b) if full else None, "fixed_len_ok": full}  # fmt: skip
    push_weights(vg, model, values="base")
    fixed["base_hash_ok"] = host_hashes(engine_model(vg)) == ref
    llm.llm_engine.remove_lora(B1_SLOT)
    model.set_adapter("default")
    model.delete_adapter("b1")
    llm.reset_prefix_cache()
    if vg.enable_sleep_mode:
        llm.sleep(level=2)  # last: level 2 frees the weights, so the hash above read live memory
        vg._llm_weights_sleeping = True
    return fixed


class FCheck(TrainerCallback):
    """Counts NaN/inf in the logs and the loads per step over steps 1-15, then writes `<out>/f_check.json`.
    probes: {name: callable} of the preset's checks (engine_dtype_ok; lora_equal, load_count's counter is `loads`);
    fixed: values measured before step 1 (comparator, base_hash_ok, lora_over_base)."""

    def __init__(self, out, preset: str, probes: dict, loads=None, fixed: dict | None = None):
        if preset not in PRESETS:
            raise ValueError(f"--f-check: no check for preset {preset!r}")
        self.out, self.preset, self.probes, self.loads, self.fixed = Path(out), preset, probes, loads, fixed or {}
        self.nan, self.counts, self.seen, self.done = 0, [], 0, False

    def on_log(self, args, state, control, logs=None, **kw):
        """NaN/inf count over the logged values of steps 1-15."""
        if state.global_step <= STEPS:
            self.nan += sum(isinstance(v, float) and not math.isfinite(v) for v in (logs or {}).values())

    def on_step_end(self, args, state, control, **kw):
        """The loads of this step; the record at step 15 or at the first failure."""
        if self.loads is not None:
            self.counts.append(self.loads.calls - self.seen)
            self.seen = self.loads.calls
        if not self.done and (state.global_step >= STEPS or self.nan or (self.counts and self.counts[-1] != 1)):
            self.finish()

    def finish(self) -> None:
        """Write f_check.json; raise when a check failed."""
        self.done = True
        rec = {"preset": self.preset, "nan": self.nan} | {k: p() for k, p in self.probes.items()} | self.fixed
        if self.loads is not None:
            rec["load_count"] = next((c for c in self.counts if c != 1), 1)
        rec["passed"] = bool(passes(rec))
        self.out.mkdir(parents=True, exist_ok=True)
        (self.out / "f_check.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
        if not rec["passed"]:
            raise RuntimeError(f"F check failed ({self.out / 'f_check.json'}): {rec}")
