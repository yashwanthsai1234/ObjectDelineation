"""GPU profile of the RL step (perf D1, docs/archive/experiments/phase2/perf/implementation.md): short g_probe runs on R1's
selected branch, one process per arm, nvidia-smi sampled alongside. Replaces profile_sft (no SFT runs: owner).

    python -I -m tools.profiling.profile --ckpt CKPT --prepared PREP --probe VIEW --out DIR [--arms fast mb32]

Arms (each: STEPS optimizer steps of 32 prompts x G = 8 answers, ITEMS probe items, step 1 left out of s/step):
- trl: TRL 1.17 as shipped (g_probe --plain), the reference for the gains below;
- fast: objdel.rl.trl_adapter.FastGRPO (per-chunk trimmed forwards, vision cache, rewards overlapped), 16 x 16;
- mb32: fast at micro-batch 32 x 8 (fewer, larger launches);
- no_ckpt: fast without gradient checkpointing (may run out of memory: recorded);
- compile / graphs: fast with the text decoder layers compiled, without and with CUDA Graphs (Stage F: adopted
  only when >= 10 % faster at equal loss);
- trace: fast with torch.profiler on optimizer step 3 (GPU busy share, kernels, syncs); last, and not timed (the
  profiler slows the steps it records).
Writes DIR/<arm>/{g_probe.json, smi.csv, log.txt[, trace]} and DIR/profile.json: per arm the exit code, s/step, the
split (generate / rewards / old+ref logprobs / update), mean answer length, the minimum importance-sampling ratio,
nvidia-smi peak memory, mean utilization and power, and for the traced arm the update step's GPU busy share.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from objdel.core.io import read_json
from objdel.core.provenance import run_dir

ROOT = Path(__file__).resolve().parents[2]
STEPS, ITEMS = 4, 8
ARMS = {  # most important first: profile.json is written after each arm, so a run cut by its cap keeps these
    "fast": {},
    "trl": {"flags": ["--plain"]},
    "graphs": {"flags": ["--compile", "reduce-overhead"]},
    "mb32": {"grpo": {"per_device_train_batch_size": 32, "gradient_accumulation_steps": 8}},
    "compile": {"flags": ["--compile", "default"]},
    "no_ckpt": {"grpo": {"gradient_checkpointing": False}},
    "trace": {"flags": ["--trace"]},
}
TIMED = [a for a in ARMS if a != "trace"]
GPU_CATS = ("kernel", "gpu_memcpy", "gpu_memset")
SYNCS = ("aten::_local_scalar_dense", "aten::nonzero", "cudaStreamSynchronize", "cudaDeviceSynchronize")  # host waits
SMI = ["nvidia-smi", "--query-gpu=timestamp,utilization.gpu,power.draw,memory.used", "--format=csv,noheader,nounits",
       "-lms", "500"]


def arm_cmd(arm: str, ckpt, prepared, probe, out: Path) -> list[str]:
    a = ARMS[arm]
    grpo = {"max_steps": STEPS, "logging_steps": 1} | a.get("grpo", {})
    flags = [x for f in a.get("flags", []) for x in ([f, str(out / arm / "trace")] if f == "--trace" else [f])]
    return ["-I", "-m", "objdel.cli.train.probe", "--model", "qwen3-vl-2b", "--ckpt", str(ckpt), "--train",
            str(prepared), "--probe", str(probe), "--out", str(out / arm), "--items", str(ITEMS),
            "--grpo", json.dumps(grpo), *flags]  # fmt: skip


def summary(log: list[dict]) -> dict:
    """Per-step g_probe log -> s/step and (FastGRPO's timers; TRL as shipped has none) its split, step 1 left out
    (vLLM, allocator and compile warm-up)."""
    steps = [h for h in log if "step_time" in h]
    rest = steps[1:] or steps

    def mean(k):
        return float(np.mean([h[k] for h in rest]))

    out = {"steps": len(steps), "s_per_step": mean("step_time"),
           "mean_length": float(np.mean([h["completions/mean_length"] for h in steps])),
           "is_ratio_min": min(h["sampling/importance_sampling_ratio/min"] for h in steps)}  # fmt: skip
    if all("time/batch_s" in h for h in rest):
        gen, rew, batch = (mean(k) for k in ("time/generate_s", "time/rewards_s", "time/batch_s"))
        out["split_s"] = {"generate": gen, "rewards": rew, "logprobs": batch - gen - rew,
                          "update": out["s_per_step"] - batch}  # fmt: skip
    return out


def covered(intervals, t0: float, t1: float) -> float:
    """The length of the union of the (start, end) intervals clipped to [t0, t1] (0.0 when t1 <= t0; an interval
    with end <= start counts nothing). busy here and nsys_baseline share it."""
    total, end = 0.0, t0
    for a, b in sorted((max(a, t0), min(b, t1)) for a, b in intervals):
        a = max(a, end)
        if b > a:
            total += b - a
            end = b
    return total


def busy(trace: dict) -> dict:
    """The traced optimizer step (ProfilerStep#N span): wall ms, GPU busy share (union of GPU intervals in the span),
    kernel launches, CUDA runtime API ms, host-sync calls by name, top kernels by total ms."""
    ev = [e for e in trace["traceEvents"] if e.get("ph") == "X"]
    steps = [e for e in ev if e.get("name", "").startswith("ProfilerStep#")]
    if not steps:
        raise ValueError("no ProfilerStep span in the trace")
    t0 = min(e["ts"] for e in steps)
    t1 = max(e["ts"] + e["dur"] for e in steps)
    gpu = covered(((e["ts"], e["ts"] + e["dur"]) for e in ev if e.get("cat") in GPU_CATS), t0, t1)
    kernels = [e for e in ev if e.get("cat") == "kernel"]
    top: dict[str, float] = {}
    for e in kernels:
        top[e["name"]] = top.get(e["name"], 0.0) + e["dur"] / 1e3
    api = sum(e["dur"] for e in ev if e.get("cat") == "cuda_runtime") / 1e3
    syncs = {n: sum(e.get("name") == n for e in ev) for n in SYNCS}
    return {"step_ms": (t1 - t0) / 1e3, "gpu_busy": gpu / (t1 - t0), "kernels": len(kernels),
            "cuda_api_ms": api, "syncs": syncs, "top_kernels_ms": dict(sorted(top.items(), key=lambda x: -x[1])[:10])}


def smi_mean(path: Path) -> dict:
    rows = [r.split(", ") for r in path.read_text().splitlines() if r.count(",") == 3]
    vals = [(float(u), float(p), float(m)) for _, u, p, m in rows if u.strip().isdigit()]
    if not vals:
        return {}
    u, p, m = zip(*vals, strict=True)
    return {"util_mean": float(np.mean(u)), "power_w_mean": float(np.mean(p)), "mem_max_mib": max(m)}


def run_arm(arm: str, a) -> dict:
    d = a.out / arm
    d.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, *arm_cmd(arm, a.ckpt, a.prepared, a.probe, a.out)]
    t = time.time()
    with open(d / "smi.csv", "w") as smi_f, open(d / "log.txt", "w") as log:
        smi = subprocess.Popen(SMI, stdout=smi_f, stderr=subprocess.DEVNULL)
        try:
            env = os.environ | {"TORCH_LOGS": "recompiles"}  # recompiles in the measured steps show in log.txt
            code = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
        finally:
            smi.terminate()
            smi.wait()
    res = {"exit": code, "wall_s": time.time() - t, "arm": ARMS[arm], "smi": smi_mean(d / "smi.csv")}
    try:  # a malformed output is recorded; the next arm still runs
        if (d / "g_probe.json").exists() and arm in TIMED:
            res |= summary(read_json(d / "g_probe.json")["log"])
        if (d / "trace" / "trace.json.gz").exists():
            with gzip.open(d / "trace" / "trace.json.gz", "rt") as f:
                res["traced_step"] = busy(json.load(f))
    except (KeyError, ValueError, OSError) as e:
        res["parse_error"] = repr(e)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--prepared", type=Path, required=True)
    ap.add_argument("--probe", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=run_dir() / "profile")
    ap.add_argument("--arms", nargs="+", choices=list(ARMS), default=list(ARMS))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    report = {}
    for arm in a.arms:
        report[arm] = run_arm(arm, a)
        (a.out / "profile.json").write_text(json.dumps(report, indent=1))
        print(arm, json.dumps({k: report[arm].get(k) for k in ("exit", "s_per_step", "split_s")}), flush=True)
    return 0 if all(r["exit"] == 0 for r in report.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
