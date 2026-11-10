"""V2 box checks for C (module-run): --sync is C5's G-SY1 (push equals TRL's), G-SY2 (speed), G-SY3 (no leak).

    .venv/bin/python -I -m experiments.phase2.c_speed_check --sync --model M --ckpt C --train ROWS --out RUN

Other flags are train_rl's (the trainer is built as train_rl builds it, with --fast-sync forced on).
"""

import argparse
import hashlib
import json
import statistics
import sys
import time
from pathlib import Path

MAX_PUSH_S, MAX_LEAK_B, MAX_SPREAD = 0.2, 1 << 20, 2.0  # G-SY2 median; G-SY3 memory; G-SY3 max / median seconds
PUSHES, LEAK_PUSHES, PERTURBED, TRL_REPEATS = 20, 50, 3, 5


def parse(argv=None):
    """(--sync, the train_rl stage of the other flags)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    ap.add_argument("--sync", action="store_true", help="C5: G-SY1, G-SY2, G-SY3")
    flag, rest = ap.parse_known_args(argv)
    from experiments.phase2 import train_rl

    return flag.sync, train_rl.parse(rest)


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


def main(argv=None) -> int:
    sync, a = parse(argv)
    if not sync:
        print("nothing to check: pass --sync", file=sys.stderr)
        return 2
    from experiments.phase2.g_probe import make_trainer, save_json
    from objdel.train import dataset

    a = a.model_copy(update={"fast_sync": True})
    trainer, _, _ = make_trainer(a, Path(a.ckpt), dataset.grpo_rows(a.train, "lidc", "train"), 1)
    out = run_sync(trainer)
    save_json(Path(a.out) / "c_check" / "sync.json", out)
    print(json.dumps(out))
    return 0 if out["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
