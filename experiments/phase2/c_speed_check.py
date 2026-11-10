"""V2 box checks for C (module-run). One process, one trainer; the modes run in this order.

    .venv/bin/python -I -m experiments.phase2.c_speed_check --pixels --sync --model M --ckpt C --train ROWS --out RUN

--pixels is C1's G-EQ1: it first makes the five batches (5 x 32 seeded train rows of the B1 rows root, 8 answers each
from the trainer's own generation path, saved as RUN/c_check/batches/b0..b4.pt in the --record-batches format plus
seeded group-centred advantages, image keys and image SHA-256 values), then compares gpu_pixels.patchify with the
HF CPU processor on 2,000 views and on all 32 images of each batch, and counts the host syncs of one pass-2 image
call (0 expected). Result: RUN/c_check/pixels.json. --sync is C5's G-SY1 (push equals TRL's), G-SY2 (speed), G-SY3
(no leak): RUN/c_check/sync.json; under option F (--numerics lora-seqmask) it checks the packed LoRA slices instead
of the weight hashes. Other flags are train_rl's (the trainer is built as train_rl builds it, with --fast-sync forced
on for --sync and --gpu-pixels auto for --pixels). A run with `--grpo '{"vllm_enable_sleep_mode": true}'` is mode B:
every push follows a level-2 sleep.
"""

import argparse
import contextlib
import hashlib
import json
import re
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import torch

from objdel.train import fast, vllm_sync

MAX_PUSH_S, MAX_SPREAD = 0.2, 2.0  # G-SY2 median push seconds; G-SY3 max / median seconds
MAX_RISE_B, MAX_LEAK_B = 100 << 20, 1 << 20  # G-SY1 peak memory rise; G-SY3 allocated memory drift
DELTA_GAIN_S = 0.05  # C5 decision: a delta push is worth re-adding only if the touched-tensor push is this much faster
PUSHES, LEAK_PUSHES, PERTURBED, REPEATS = 20, 50, 3, 5
PROMPTS, GEN_TOKENS, LORA_ID = 16, 32, 1  # G-SY1 greedy generation; sync_lora's default slot
PEFT_KEY = re.compile(r"layers\.(\d+)\.(?:self_attn|mlp)\.(\w+)\.lora_([AB])\.weight$")
VLLM_KEY = re.compile(r"layers\.(\d+)\.(?:self_attn|mlp)\.(qkv_proj|o_proj|gate_up_proj|down_proj)$")
PACKED = {"q_proj": ("qkv_proj", 0), "k_proj": ("qkv_proj", 1), "v_proj": ("qkv_proj", 2), "o_proj": ("o_proj", 0),
          "gate_proj": ("gate_up_proj", 0), "up_proj": ("gate_up_proj", 1), "down_proj": ("down_proj", 0)}  # fmt: skip
BATCHES, PER_BATCH, G, VIEWS, CHUNK = 5, 32, 8, 2000, 64  # five batches of 32 prompts x G answers; G-EQ1's views


def parse(argv=None):
    """(--pixels and --sync, the train_rl stage of the other flags)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    ap.add_argument("--pixels", action="store_true", help="C1: the five batches, then G-EQ1")
    ap.add_argument("--sync", action="store_true", help="C5: G-SY1, G-SY2, G-SY3")
    flag, rest = ap.parse_known_args(argv)
    from experiments.phase2 import train_rl

    return flag, train_rl.parse(rest)


def refuse_val(root) -> None:
    """The five batches and the pixel views are train rows: a val view or eval folder raises ValueError."""
    for part in Path(root).parts:
        if re.fullmatch(r"(.*[-_])?val([-_].*)?|eval.*", part):
            raise ValueError(f"{root}: a val view ({part}); C reads train rows only")


def draw_batches(n_rows: int, seed: int = 0) -> np.ndarray:
    """(BATCHES, PER_BATCH) distinct row indices, seeded."""
    if n_rows < BATCHES * PER_BATCH:
        raise ValueError(f"{n_rows} rows: the five batches need {BATCHES * PER_BATCH} distinct rows (160)")
    return np.random.default_rng(seed).choice(n_rows, BATCHES * PER_BATCH, replace=False).reshape(BATCHES, PER_BATCH)


def centred_advantages(batches: int = BATCHES):
    """(batches, 256) float tensor of seeded normal values (seed 0), centred per group of G: G-EQ4 compares two code
    paths on the same loss, so the values need not be rewards."""
    import torch

    x = torch.randn(batches, PER_BATCH, G, generator=torch.Generator().manual_seed(0))
    return (x - x.mean(-1, keepdim=True)).reshape(batches, PER_BATCH * G)


def make_batches(generate, rows, n_rows: int, out) -> list[Path]:
    """The five check batches b0..b4.pt in `out`. rows: dataset.OpenOnce over the train rows (key, image); each
    batch is 32 drawn rows x G identical copies sharing one image object, as TRL's sampler hands them.
    generate(inputs, step) -> the batch's record (record.batch_record); this adds the seeded advantages and, per
    distinct image in first-appearance order, its row key and pixel SHA-256 (hashed on the host)."""
    import torch

    from objdel.hashing import pixels_sha256
    from objdel.train import record

    out, adv = record.open_dir(out), centred_advantages()
    paths = []
    for b, idx in enumerate(draw_batches(n_rows)):
        inputs = rows.__getitems__([int(i) for i in idx for _ in range(G)])
        rec = generate(inputs, b)
        distinct = [r for r in inputs[::G]]
        rec |= {"advantages": adv[b], "image_key": [r["key"] for r in distinct],
                "image_sha256": [pixels_sha256(np.asarray(r["image"].convert("RGB"))) for r in distinct]}
        paths.append(out / f"b{b}.pt")
        torch.save(rec, paths[-1])
    return paths


def compare_pixels(images, ip, device, backend: str, chunk: int = CHUNK) -> dict:
    """G-EQ1: gpu_pixels.patchify(backend) against the processor's own output, image by image, in chunks:
    {backend, views, mismatches, first (index of the first differing image or None), median_ms (per chunk call)}."""
    import torch

    from objdel.train import gpu_pixels

    bad, ms = [], []
    for start in range(0, len(images), chunk):
        part = images[start:start + chunk]
        want = ip(part, return_tensors="pt")
        t0 = time.perf_counter()
        pv, grid = gpu_pixels.patchify(part, ip, device=device, backend=backend)
        if torch.device(device).type == "cuda":
            torch.cuda.synchronize()
        ms.append((time.perf_counter() - t0) * 1e3)
        sizes = want["image_grid_thw"].prod(-1).tolist()
        bad += [start + i for i, (a, b) in enumerate(zip(pv.cpu().split(sizes), want["pixel_values"].split(sizes),
                                                         strict=True)) if not torch.equal(a, b)]
        if not torch.equal(grid, want["image_grid_thw"]):
            bad += [start + i for i in range(len(part))]
    return {"backend": backend, "views": len(images), "mismatches": len(bad), "first": min(bad, default=None),
            "median_ms": statistics.median(ms)}


def hash_params(pairs) -> dict[str, str]:
    """name -> SHA-256 of the raw bytes, hashed on the host one tensor at a time."""
    return {
        n: hashlib.sha256(t.detach().contiguous().view(-1).view(torch.uint8).cpu().numpy().tobytes()).hexdigest()
        for n, t in pairs
    }


def sync_verdict(seconds, mem, equal, rises) -> dict:
    """G-SY1: every `equal` check holds and max_memory_allocated rose by < 100 MB in every window (`rises`); G-SY2:
    median push seconds; G-SY3: allocated memory flat over the pushes and seconds stable."""
    med, rise = statistics.median(seconds), max(rises)
    ok = {
        "G-SY1": all(equal.values()) and rise < MAX_RISE_B,
        "G-SY2": med <= MAX_PUSH_S,
        "G-SY3": max(abs(m - mem[0]) for m in mem) <= MAX_LEAK_B and max(seconds) <= MAX_SPREAD * med,
    }
    return {
        "median_s": med,
        "max_s": max(seconds),
        "mem_drift_b": max(abs(m - mem[0]) for m in mem),
        "peak_rise_b": rise,
        "checks": equal,
        **ok,
        "passed": all(ok.values()),
    }


def _engine_model(vg):
    return vg.llm.llm_engine.model_executor.driver_worker.model_runner.model


def _engine_hashes(vg) -> dict[str, str]:
    return hash_params(_engine_model(vg).named_parameters())


def prompts(trainer) -> list[dict]:
    tok = trainer.processing_class.tokenizer
    return [{"prompt_token_ids": tok.encode(f"Count from {i} to {i + 20}:")} for i in range(PROMPTS)]


def greedy(vg, rows) -> list[list[int]]:
    from vllm import SamplingParams

    if vg.enable_sleep_mode:
        vg.llm.wake_up(tags=["kv_cache"])
    outs = vg.llm.generate(rows, sampling_params=SamplingParams(temperature=0.0, max_tokens=GEN_TOKENS), use_tqdm=False)
    return [list(o.outputs[0].token_ids) for o in outs]


def to_sleep(vg) -> None:
    """Mode B: level-2 sleep discards the weights, the KV cache and the LoRA slots, so the next push must restore
    them. TRL's flag is the one sleep state (an engine that sleeps already is not slept again)."""
    if vg.enable_sleep_mode and not vg._llm_weights_sleeping:
        vg.llm.sleep(level=2)
        vg._llm_weights_sleeping = True


@contextlib.contextmanager
def measured(rises):
    """Appends how far max_memory_allocated rose above the memory allocated at entry."""
    torch.cuda.reset_peak_memory_stats()
    base = torch.cuda.memory_allocated()
    yield
    rises.append(torch.cuda.max_memory_allocated() - base)


@contextlib.contextmanager
def counting_add_lora(engine, counts):
    """engine.add_lora calls are added to counts[-1]."""
    real = engine.add_lora

    def counted(*args, **kwargs):
        counts[-1] += 1
        return real(*args, **kwargs)

    engine.add_lora = counted
    try:
        yield
    finally:
        engine.add_lora = real


def lora_slices(vg) -> dict[tuple[int, str, int], tuple]:
    """(layer, vLLM module, slice) -> (lora_a, lora_b) of the adapter in slot LORA_ID, as vLLM's LoRA manager holds it
    (packed qkv_proj and gate_up_proj keep one slice per projection). Read from vLLM's private LoRA manager; unchecked
    on the Mac."""
    mgr = vg.llm.llm_engine.model_executor.driver_worker.model_runner.lora_manager._adapter_manager
    slot = mgr.lora_index_to_id.index(LORA_ID)
    out = {}
    for name, layer in mgr.modules.items():
        m = VLLM_KEY.search(name)
        if m:
            for k, (a, b) in enumerate(zip(layer.lora_a_stacked, layer.lora_b_stacked, strict=True)):
                out[(int(m[1]), m[2], k)] = (a[slot, 0], b[slot, 0])
    return out


def packed_mismatches(model, slices) -> list[str]:
    """The adapter tensors whose vLLM slice is not exactly A.to(bf16) (lora_a) or alpha / r * B.to(bf16) (lora_b: vLLM
    casts the adapter to bf16 and folds the scaling into lora_b). A slice vLLM lacks counts; so does reading none."""
    from peft import get_peft_model_state_dict

    cfg = model.peft_config["default"]
    bad, n = [], 0
    for key, t in get_peft_model_state_dict(model).items():
        m = PEFT_KEY.search(key)
        if m is None:
            continue
        n += 1
        a, b = slices.get((int(m[1]), *PACKED[m[2]]), (None, None))
        got, want = (a, t.to(torch.bfloat16)) if m[3] == "A" else (b, cfg.lora_alpha / cfg.r * t.to(torch.bfloat16))
        if got is None or not torch.equal(got, want):
            bad.append(key)
    return bad or ([] if n else ["no LoRA tensor in the model"])


def check_merged(vg, model, rows, rises, i) -> dict[str, bool]:
    """One perturbed state: our push, hash, greedy tokens; then TRL's per-parameter push of the same state, hash,
    tokens. Ours first, so a push that changes nothing differs from TRL's."""
    to_sleep(vg)
    with measured(rises):
        vg.sync_weights()
        got = _engine_hashes(vg)
    tokens = greedy(vg, rows)
    to_sleep(vg)
    with measured(rises):
        type(vg).sync_weights(vg)
        want = _engine_hashes(vg)
    return {f"state {i} hashes": got == want, f"state {i} tokens": tokens == greedy(vg, rows)}


def check_lora(vg, model, rises, i, detail) -> dict[str, bool]:
    """One perturbed state under F: the packed slices; in mode B also the base W the sync pushed, against a
    per-parameter push of p.data."""
    to_sleep(vg)
    with measured(rises):
        vg.sync_weights()
        bad = packed_mismatches(model, lora_slices(vg))
        out = {f"state {i} packed": not bad}
        detail.extend(bad[:5])
        if vg.enable_sleep_mode:
            got = _engine_hashes(vg)
            for pair in fast.lora_merged_params(model, with_delta=False):
                _engine_model(vg).load_weights([pair])
            out[f"state {i} base"] = got == _engine_hashes(vg)
    return out


def touched_names(model) -> set[str]:
    """The vLLM names of the weights a LoRA delta moves (196 for the 2B model)."""
    return {n.removeprefix("base_model.model.").replace(".base_layer", "") for n, _ in model.named_parameters()
            if n.endswith(".base_layer.weight")}  # fmt: skip


def median_seconds(fn, n) -> float:
    return statistics.median(fn() for _ in range(n))


def _timed(fn, *args) -> float:
    t0 = time.perf_counter()
    fn(*args)
    return time.perf_counter() - t0


def run_sync(trainer, lora: bool) -> dict:
    """The vLLM-side checks around `vg.sync_weights`, which make_trainer installed (push_weights, or sync_lora under
    --numerics lora-seqmask). lora: option F, so TRL's per-parameter loop is no reference (it would write W + d into the base the
    adapter rides on) and the checks are the packed slices and the adapter load count."""
    vg, model = trainer.vllm_generation, trainer.vllm_generation.model
    peft = [p for p in model.parameters() if p.requires_grad]
    gen = torch.Generator(device=peft[0].device).manual_seed(0)
    rows, equal, rises, adds, detail = prompts(trainer), {}, [], [], []
    with counting_add_lora(vg.llm.llm_engine, adds):
        for i in range(PERTURBED):  # G-SY1: PERTURBED adapter states
            with torch.no_grad():
                for p in peft:
                    p.add_(torch.randn(p.shape, generator=gen, device=p.device, dtype=p.dtype) * 1e-3)
            adds.append(0)
            equal |= check_lora(vg, model, rises, i, detail) if lora else check_merged(vg, model, rows, rises, i)
    if lora:
        equal["one add_lora per sync"] = adds == [1] * PERTURBED
    seconds = [vg.sync_weights().seconds for _ in range(PUSHES)]  # G-SY2
    trl_s = touched_s = None
    if not lora:  # the one-call push of the touched tensors alone and TRL's loop, REPEATS times each
        touched = touched_names(model)
        touched_s = median_seconds(lambda: vllm_sync._load(
            vg, ((n, t) for n, t in fast.lora_merged_params(model) if n in touched)).seconds, REPEATS)  # fmt: skip
        trl_s = median_seconds(lambda: _timed(type(vg).sync_weights, vg), REPEATS)
    mem = []
    for _ in range(LEAK_PUSHES):  # G-SY3
        vg.sync_weights()
        mem.append(torch.cuda.memory_allocated())
    gain = None if touched_s is None else statistics.median(seconds) - touched_s
    return sync_verdict(seconds, mem, equal, rises) | {
        "touched_tensors": None if lora else len(touched),
        "touched_s": touched_s,
        "trl_loop_s": trl_s,
        "readd_delta_push": None if gain is None else gain >= DELTA_GAIN_S,
        "add_lora_per_sync": adds if lora else None,
        "packed_mismatches": detail,
    }


def trainer_generate(trainer, scratch):
    """generate(inputs, step) of make_batches on a real trainer: one generation batch through its own path
    (_generate_and_score_completions in train mode), recorded by FastGRPO as --record-batches does and read back."""
    import torch

    from objdel.train import record

    trainer.record_dir, trainer.checks_from_step = record.open_dir(scratch), 1

    def generate(inputs, step):
        trainer.model.train()
        trainer.state.global_step = step
        trainer._generate_and_score_completions(inputs)
        f = Path(scratch) / record.file_name(step + 1)
        rec = torch.load(f, weights_only=True)
        f.unlink()
        return rec

    return generate


def run_pixels(trainer, rows, n_rows: int, out: Path) -> dict:
    """C1 on the box: the five batches, then G-EQ1 for each gpu_pixels backend on VIEWS drawn views and on the 32
    images of every batch, and the host syncs of one pass-2 image call (G-SP4: 0)."""
    import shutil

    import torch

    from objdel.train import gpu_pixels, nvtx

    check = out / "c_check"
    paths = make_batches(trainer_generate(trainer, check / "scratch"), rows, n_rows, check / "batches")
    shutil.rmtree(check / "scratch")
    ip = trainer.images.inner
    pool = np.random.default_rng(1).choice(n_rows, min(VIEWS, n_rows), replace=False)
    views = [r["image"] for r in rows.__getitems__([int(i) for i in pool])]
    batches = [[r["image"] for r in rows.__getitems__([int(i) for i in idx])] for idx in draw_batches(n_rows)]
    backends = {b: {"views": compare_pixels(views, ip, "cuda", b),
                    "batches": [compare_pixels(x, ip, "cuda", b) for x in batches]} for b in gpu_pixels.BACKENDS}
    ok = [b for b, r in backends.items() if not r["views"]["mismatches"]
          and not any(x["mismatches"] for x in r["batches"])]  # fmt: skip
    memo, rows256 = trainer.images, rows.__getitems__([int(i) for i in draw_batches(n_rows)[0] for _ in range(G)])
    memo.clear()
    memo.prefill([r["image"] for r in rows256], return_tensors="pt")
    with nvtx.count_syncs() as n:
        memo([r["image"] for r in rows256], return_tensors="pt")
    torch.cuda.synchronize()  # outside the count: it is itself a sync
    return {"backends": backends, "equal_backends": ok, "pass2_syncs": n.syncs, "batches": [p.name for p in paths],
            "passed": bool(ok) and n.syncs == 0}


def main(argv=None) -> int:
    flag, a = parse(argv)
    if not (flag.pixels or flag.sync):
        print("nothing to check: pass --pixels and/or --sync", file=sys.stderr)
        return 2
    if flag.pixels:
        refuse_val(a.train)
    from experiments.phase2.g_probe import make_trainer, save_json
    from objdel.train import dataset

    a = a.model_copy(update={"fast_sync": flag.sync or a.fast_sync, "gpu_pixels": "auto" if flag.pixels else a.gpu_pixels})
    out = Path(a.out)
    train = dataset.grpo_rows(a.train, "lidc", "train")
    trainer, _, _ = make_trainer(a, Path(a.ckpt), train, 1)
    outs = {}
    if flag.pixels:
        outs["pixels"] = run_pixels(trainer, dataset.OpenOnce(train), len(train), out)
    if flag.sync:
        outs["sync"] = run_sync(trainer, trainer.lora_dir is not None)  # option F: make_trainer set the adapter dir
    for name, res in outs.items():
        save_json(out / "c_check" / f"{name}.json", res)
        print(json.dumps({name: res}))
    return 0 if all(o["passed"] for o in outs.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
