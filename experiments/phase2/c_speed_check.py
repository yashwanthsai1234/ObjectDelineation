"""V2 box checks for C (module-run). One process, one trainer; the modes run in this order.

    .venv/bin/python -I -m experiments.phase2.c_speed_check --pixels --sync --model M --ckpt C --train ROWS --out RUN

--pixels is C1's G-EQ1: it first makes the five batches (5 x 32 seeded train rows of the B1 rows root, 8 answers each
from the trainer's own generation path, saved as RUN/c_check/batches/b0..b4.pt in the --record-batches format plus
seeded group-centred advantages, image keys and image SHA-256 values), then compares gpu_pixels.patchify with the
HF CPU processor on 2,000 views and on all 32 images of each batch, and counts the host syncs of one pass-2 image
call (0 expected). Result: RUN/c_check/pixels.json. --sync is C5's G-SY1 (push equals TRL's), G-SY2 (speed), G-SY3
(no leak): RUN/c_check/sync.json. Other flags are train_rl's (the trainer is built as train_rl builds it, with
--fast-sync forced on for --sync and --gpu-pixels auto for --pixels).
"""

import argparse
import hashlib
import json
import re
import statistics
import sys
import time
from pathlib import Path

import numpy as np

MAX_PUSH_S, MAX_LEAK_B, MAX_SPREAD = 0.2, 1 << 20, 2.0  # G-SY2 median; G-SY3 memory; G-SY3 max / median seconds
PUSHES, LEAK_PUSHES, PERTURBED, TRL_REPEATS = 20, 50, 3, 5
BATCHES, PER_BATCH, G, VIEWS, CHUNK = 5, 32, 8, 2000, 64  # five batches of 32 prompts x G answers; G-EQ1's views


def parse(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    ap.add_argument("--pixels", action="store_true", help="C1: the five batches, then G-EQ1")
    ap.add_argument("--sync", action="store_true", help="C5: G-SY1, G-SY2, G-SY3")
    flag, rest = ap.parse_known_args(argv)
    from experiments.phase2 import train_rl

    a = train_rl.parse(rest)
    a.pixels, a.sync = flag.pixels, flag.sync
    return a


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
    import torch

    return {
        n: hashlib.sha256(t.detach().contiguous().view(-1).view(torch.uint8).cpu().numpy().tobytes()).hexdigest()
        for n, t in pairs
    }


def sync_verdict(seconds, mem, equal) -> dict:
    med = statistics.median(seconds)
    ok = {
        "G-SY1": all(equal),
        "G-SY2": med <= MAX_PUSH_S,
        "G-SY3": max(abs(m - mem[0]) for m in mem) <= MAX_LEAK_B and max(seconds) <= MAX_SPREAD * med,
    }
    return {
        "median_s": med,
        "max_s": max(seconds),
        "mem_drift_b": max(abs(m - mem[0]) for m in mem),
        **ok,
        "passed": all(ok.values()),
    }


def _engine_hashes(vg) -> dict[str, str]:
    return hash_params(vg.llm.llm_engine.model_executor.driver_worker.model_runner.model.named_parameters())


def run_sync(trainer) -> dict:
    import torch

    from objdel.train import vllm_sync

    vg, model = trainer.vllm_generation, trainer.vllm_generation.model
    trl_sync = lambda: type(vg).sync_weights(vg)  # noqa: E731  TRL's per-parameter loop
    peft = [p for p in model.parameters() if p.requires_grad]
    gen = torch.Generator(device=peft[0].device).manual_seed(0)
    equal = []
    for _ in range(PERTURBED):  # G-SY1: three perturbed adapter states, TRL's push then ours
        with torch.no_grad():
            for p in peft:
                p.add_(torch.randn(p.shape, generator=gen, device=p.device, dtype=p.dtype) * 1e-3)
        trl_sync()
        want = _engine_hashes(vg)
        vllm_sync.push_weights(vg, model)
        equal.append(_engine_hashes(vg) == want)
    seconds = [vllm_sync.push_weights(vg, model).seconds for _ in range(PUSHES)]
    mem = []
    for _ in range(LEAK_PUSHES):  # G-SY3
        vllm_sync.push_weights(vg, model)
        mem.append(torch.cuda.memory_allocated())
    t0 = time.perf_counter()
    for _ in range(TRL_REPEATS):
        trl_sync()
    return sync_verdict(seconds, mem, equal) | {"trl_loop_s": (time.perf_counter() - t0) / TRL_REPEATS}


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
    a = parse(argv)
    if not (a.pixels or a.sync):
        print("nothing to check: pass --pixels and/or --sync", file=sys.stderr)
        return 2
    if a.pixels:
        refuse_val(a.train)
    from experiments.phase2.g_probe import make_trainer, save_json
    from objdel.train import dataset

    a.fast_sync, a.gpu_pixels = a.sync or a.fast_sync, "auto" if a.pixels else a.gpu_pixels
    train = dataset.grpo_rows(a.train, "lidc", "train")
    trainer, _, _ = make_trainer(a, a.ckpt, train, 1)
    outs = {}
    if a.pixels:
        outs["pixels"] = run_pixels(trainer, dataset.OpenOnce(train), len(train), a.out)
    if a.sync:
        outs["sync"] = run_sync(trainer)
    for name, out in outs.items():
        save_json(a.out / "c_check" / f"{name}.json", out)
        print(json.dumps({name: out}))
    return 0 if all(o["passed"] for o in outs.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
